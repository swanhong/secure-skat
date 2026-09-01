import argparse
import csv
import gzip
import os
import struct
import subprocess
import time

BILLING_PROJECT = (
    os.environ.get("GOOGLE_CLOUD_PROJECT")
    or os.environ.get("GOOGLE_PROJECT")
)

DEFAULT_VAT_PATH = (
    "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/"
    "snpindel/aux/vat/vat_complete.bgz.tsv.gz"
)

# Raw VAT column indexes (zero-based).
COL_CONTIG = 2
COL_POSITION = 3
COL_REF = 4
COL_ALT = 5
COL_GVS_AF = 8
COL_GENE_SYMBOL = 43
COL_CONSEQUENCE = 46
COL_GENE_ID = 53
COL_CANONICAL = 55
COL_GNOMAD_AF = 56
COL_MANE_SELECT = 105
COL_LOF = 109

MIN_COLUMNS = 110
LINEAR_SHIFT = 14
PROGRESS_EVERY = 2_000_000

MISSING_VALUES = {"", "."}

OUTPUT_COLUMNS = (
    "variant_key",
    "gene_id",
    "gene_symbol",
    "LoF",
    "consequence",
    "gnomad_af",
    "gvs_af",
)

TRANSCRIPT_PRIORITY = (
    ("mane_select_name", COL_MANE_SELECT, "not_empty", None),
    (
        "is_canonical_transcript",
        COL_CANONICAL,
        "equals_ignore_case",
        "true",
    ),
    ("LoF", COL_LOF, "equals", "HC"),
)


def transcript_priority(fields: list[str]) -> int:
    for priority, (name, column, condition, expected) in enumerate(
        TRANSCRIPT_PRIORITY
    ):
        value = fields[column].strip()
        if condition == "not_empty":
            matches = value not in MISSING_VALUES
        elif condition == "equals_ignore_case":
            matches = value.casefold() == expected.casefold()
        elif condition == "equals":
            matches = value == expected
        else:
            raise ValueError(
                f"unsupported condition for {name}: {condition}"
            )

        if matches:
            return priority

    return len(TRANSCRIPT_PRIORITY)


def load_pvar_variants(pvar_path: str):
    with open(pvar_path) as pvar_file:
        for line in pvar_file:
            if line.startswith("#CHROM"):
                columns = line.lstrip("#").rstrip().split("\t")
                break
        else:
            raise ValueError(f"PVAR header not found: {pvar_path}")

        position_column = columns.index("POS")
        ref_column = columns.index("REF")
        alt_column = columns.index("ALT")
        filter_column = (
            columns.index("FILTER") if "FILTER" in columns else None
        )

        variants = set()
        minimum_position = None
        maximum_position = 0

        for line in pvar_file:
            fields = line.rstrip().split("\t")
            if (
                filter_column is not None
                and fields[filter_column] not in {"PASS", "."}
            ):
                continue

            alt = fields[alt_column]
            # remove if ALT contains a comma (multi-allelic)
            if "," in alt:
                continue

            position = int(fields[position_column])
            ref = fields[ref_column]
            variants.add(f"{position}:{ref}:{alt}")

            if minimum_position is None or position < minimum_position:
                minimum_position = position
            maximum_position = max(maximum_position, position)

    if not variants:
        raise ValueError(f"no PASS biallelic variants found in {pvar_path}")

    return variants, minimum_position, maximum_position


def annotation_candidate(fields: list[str], source_order: int):
    '''
        remove
            - lack of len(column)
            - both LoF & consequence are empty
            - gene_id is empty

        variant key = f"{contig}:{position}:{ref}:{alt}"
        gene symbol = fields[COL_GENE_SYMBOL].strip()
        group_key  = (variant key, gene_id)
        rank = (transcript priority, VAT source_order)
        output_row = final row that will be written to the output file
    '''
    if len(fields) < MIN_COLUMNS:
        return None

    lof = fields[COL_LOF].strip()
    consequence = fields[COL_CONSEQUENCE].strip()
    if lof in MISSING_VALUES and consequence in MISSING_VALUES:
        return None

    gene_id = fields[COL_GENE_ID].strip()
    if gene_id in MISSING_VALUES:
        return None

    variant_key = ":".join([
        fields[COL_CONTIG].strip(),
        fields[COL_POSITION].strip(),
        fields[COL_REF].strip(),
        fields[COL_ALT].strip(),
    ])

    group_key = (variant_key, gene_id)
    rank = (transcript_priority(fields), source_order)
    output_row = {
        "variant_key": variant_key,
        "gene_id": gene_id,
        "gene_symbol": fields[COL_GENE_SYMBOL].strip(),
        "LoF": lof,
        "consequence": consequence,
        "gnomad_af": fields[COL_GNOMAD_AF].strip(),
        "gvs_af": fields[COL_GVS_AF].strip(),
    }
    return group_key, rank, output_row


def select_best_rows(
    lines,
    chromosome: str,
    wanted_variants: set[str],
    maximum_position: int,
):
    best_rows = {}
    chromosome_started = False
    source_order = 0
    matched_rows = 0
    started_at = time.time()

    for line in lines:
        head = line.split("\t", 6)

        if len(head) < 7:
            continue

        row_chromosome = head[COL_CONTIG]
        if row_chromosome != chromosome:
            if chromosome_started:
                break
            continue

        chromosome_started = True
        source_order += 1

        if source_order % PROGRESS_EVERY == 0:
            elapsed = time.time() - started_at
            print(
                f"... {source_order:,} VAT rows, "
                f"{matched_rows:,} matched transcript rows, "
                f"{len(best_rows):,} variant-gene rows, "
                f"{source_order / elapsed:,.0f} rows/s",
                flush=True,
            )

        position = int(head[COL_POSITION])
        if position > maximum_position:
            break

        variant_key = f"{position}:{head[COL_REF]}:{head[COL_ALT]}"
        if variant_key not in wanted_variants:
            continue

        fields = line.rstrip("\n").split("\t")
        candidate = annotation_candidate(fields, source_order)
        if candidate is None:
            continue

        matched_rows += 1
        group_key, rank, output_row = candidate
        current = best_rows.get(group_key)
        if current is None or rank < current[0]:
            best_rows[group_key] = (rank, output_row)

    elapsed = time.time() - started_at
    print(
        f"scanned {source_order:,} {chromosome} VAT rows; "
        f"matched {matched_rows:,} transcript rows and selected "
        f"{len(best_rows):,} variant-gene rows in {elapsed:,.0f}s",
        flush=True,
    )
    return best_rows


def chromosome_offset(
    vat_path: str,
    chromosome: str,
    genomic_position: int,
) -> int:
    if not BILLING_PROJECT:
        raise ValueError(
            "GOOGLE_CLOUD_PROJECT or GOOGLE_PROJECT is required"
        )
    command = [
        "gsutil",
        "-u",
        BILLING_PROJECT,
        "cat",
        f"{vat_path}.tbi",
    ]
    compressed_index = subprocess.run(
        command, check=True, capture_output=True,
    ).stdout

    index = gzip.decompress(compressed_index)

    if index[:4] != b"TBI\x01":
        raise ValueError("invalid tabix index file")

    # read chromosome names from the index and find the target chromosome
    # byte 4-7: #chr count
    # byte 32-35: names length
    # byte 36-: null-terminated chromosome names
    reference_count = struct.unpack_from("<i", index, 4)[0]
    names_length = struct.unpack_from("<i", index, 32)[0]
    names = [
        name.decode()
        for name in index[36:36 + names_length].split(b"\x00")
        if name
    ]

    if chromosome not in names:
        raise ValueError(f"chromosome not found in VAT index: {chromosome}")

    target_reference = names.index(chromosome)
    position = 36 + names_length

    for reference in range(reference_count):
        bin_count = struct.unpack_from("<i", index, position)[0]
        position += 4

        for _ in range(bin_count):
            chunk_count = struct.unpack_from("<i", index, position + 4)[0]
            position += 8 + 16 * chunk_count

        interval_count = struct.unpack_from("<i", index, position)[0]
        position += 4

        intervals = struct.unpack_from(
            f"<{interval_count}Q",
            index,
            position,
        )
        position += 8 * interval_count

        if reference == target_reference:
            interval = min(
                genomic_position >> LINEAR_SHIFT,
                interval_count - 1,
            )
            while interval >= 0 and intervals[interval] == 0:
                interval -= 1
            if interval >= 0:
                return intervals[interval] >> 16

            raise ValueError(
                f"no data offset found for {chromosome}:"
                f"{genomic_position}"
            )

    raise ValueError(f"chromosome not found: {chromosome}")


def stream_vat_lines(vat_path: str, start: int):
    gsutil_process = subprocess.Popen(
        [
            "gsutil",
            "-u",
            BILLING_PROJECT,
            "cat",
            "-r",
            f"{start}-",
            vat_path,
        ],
        stdout=subprocess.PIPE,
    )

    gzip_process = subprocess.Popen(
        ["gzip", "-dc"],
        stdin=gsutil_process.stdout, stdout=subprocess.PIPE, text=True,
    )

    gsutil_process.stdout.close()

    try:
        for line in gzip_process.stdout:
            yield line
    finally:
        gzip_process.stdout.close()

        if gzip_process.poll() is None:
            gzip_process.terminate()
        if gsutil_process.poll() is None:
            gsutil_process.terminate()

        gzip_process.wait()
        gsutil_process.wait()


def write_annotations(output_path: str, best_rows) -> int:
    output_directory = os.path.dirname(output_path)
    if output_directory:
        os.makedirs(output_directory, exist_ok=True)

    with open(output_path, "w", newline="") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for _, output_row in best_rows.values():
            writer.writerow(output_row)
    return len(best_rows)


def main():
    parser = argparse.ArgumentParser(
        description="create a simplified AoU VAT annotation file."
    )
    parser.add_argument(
        "--chromosome",
        required=True,
        help="chromosome number"
    )
    parser.add_argument(
        "--vat-path",
        default=DEFAULT_VAT_PATH,
        help="path to the AoU VAT file (default: %(default)s)"
    )
    parser.add_argument(
        "--pvar",
        required=True,
        help="PVAR containing the PASS biallelic variants to annotate",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="path to the output TSV file"
    )
    args = parser.parse_args()

    if not BILLING_PROJECT:
        parser.error("GOOGLE CLOUD_PROJECT or GOOGLE_PROJECT environment variable is required")

    chromosome = args.chromosome
    if not chromosome.startswith("chr"):
        chromosome = f"chr{chromosome}"

    output_path = args.output
    if not output_path:
        output_path = os.path.expanduser(
            f"~/fed_prep_out/{chromosome}_annotation.tsv"
        )

    print(f"chromosome: {chromosome}")
    print(f"vat_path: {args.vat_path}")
    print(f"pvar_path: {args.pvar}")
    print(f"output_path: {output_path}")

    wanted_variants, minimum_position, maximum_position = load_pvar_variants(
        args.pvar
    )
    print(
        f"PVAR variants: {len(wanted_variants):,}; "
        f"range: {minimum_position:,}-{maximum_position:,}"
    )

    start = chromosome_offset(
        args.vat_path,
        chromosome,
        minimum_position,
    )
    print(f"VAT start byte: {start:,}")

    lines = stream_vat_lines(args.vat_path, start)
    try:
        best_rows = select_best_rows(
            lines,
            chromosome,
            wanted_variants,
            maximum_position,
        )
    finally:
        lines.close()

    if not best_rows:
        raise ValueError(f"No annotations found for chromosome {chromosome}")

    written = write_annotations(output_path, best_rows)
    print(f"Wrote {written:,} annotations to {output_path}")

if __name__ == "__main__":
    main()
