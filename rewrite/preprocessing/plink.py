import csv, subprocess, tempfile
from pathlib import Path
import numpy as np

# add items if input .raw file has other values than these, but PLINK2's --export A include-alt should not emit them
HARD_CALL = {"0": 0, "1": 1, "2": 2, "NA": 0}


def parse_raw_column(column: str) -> tuple[str, str]:
    """Split a .raw genotype column into its variant key and counted allele.

    `--export A include-alt` names each column `<variant_id>_<counted>(/<other>)`.
    """
    name = column.split("(", 1)[0]
    key, _, counted_allele = name.rpartition("_")
    if not key:
        raise ValueError(f"unexpected .raw genotype column: {column}")
    return key, counted_allele


def plink_extract(
        pgen_prefix: Path,
        sample_ids: tuple[str, ...],
        variant_keys: tuple[str, ...],
        plink2_bin: str = "plink2",
) -> tuple[
    np.ndarray,
    tuple[str, ...],
]:
    """
        Export ALT dosages and align rows to sample_ids.

    Args:
        pgen_prefix: The prefix of the PLINK pgen dataset (without the .pgen extension).
        sample_ids: A tuple of sample IDs to keep.
        variant_keys: A tuple of variant keys to extract.

    Returns:
        A tuple containing:
            - int8 dosages, shape (len(sample_ids), len(emitted_keys)), missing as 0
            - the variant keys PLINK2 emitted, in its own output order
    """
    empty = (np.zeros((len(sample_ids), 0), dtype=np.int8), ())
    if not variant_keys:
        return empty

    with tempfile.TemporaryDirectory() as tmpdir:
        work_dir = Path(tmpdir)
        keep_path = work_dir / "samples.keep"
        extract_path = work_dir / "variants.extract"
        allele_path = work_dir / "alleles.txt"
        out_prefix = work_dir / "genotypes"

        keep_path.write_text(
            "#IID\n" + "".join(f"{sample_id}\n" for sample_id in sample_ids)
        )
        extract_path.write_text(
            "".join(f"{variant_key}\n" for variant_key in variant_keys)
        )
        allele_path.write_text(
            "".join(f"{variant_key}\t{variant_key.rsplit(':', 1)[1]}\n" for variant_key in variant_keys)
        )

        completed = subprocess.run([
                plink2_bin,
                "--pfile", str(pgen_prefix),
                "--keep", str(keep_path),
                "--extract", str(extract_path),
                "--export", "A", "include-alt",
                "--export-allele", str(allele_path),
                "--max-alleles", "2",
                "--out", str(out_prefix),
            ],
            capture_output=True,
            text=True,
        )

        if completed.returncode != 0:
            if "No variants remaining" in completed.stdout + completed.stderr:
                return empty
            raise RuntimeError(
                f"plink2 failed ({completed.returncode}):\n"
                + (completed.stdout + completed.stderr)[-2000:]
            )

        with Path(f"{out_prefix}.raw").open(newline="") as raw_file:
            reader = csv.reader(raw_file, delimiter="\t")
            header = next(reader)
            iid_column = header.index("IID")
            genotype_start = header.index("PHENOTYPE") + 1

            emitted_keys = []
            for column in header[genotype_start:]:
                key, counted_allele = parse_raw_column(column)
                expected_allele = key.rsplit(":", 1)[1]
                if counted_allele != expected_allele:
                    raise ValueError(
                        f"{key}: counted allele {counted_allele!r} "
                        f"is not the key's ALT {expected_allele!r}"
                    )
                emitted_keys.append(key)

            row_of_sample = {
                sample_id: row for row, sample_id in enumerate(sample_ids)
            }
            matrix = np.zeros(
                (len(sample_ids), len(emitted_keys)),
                dtype=np.int8,
            )
            filled = np.zeros(len(sample_ids), dtype=bool)

            for fields in reader:
                row = row_of_sample.get(fields[iid_column])
                if row is None or filled[row]:
                    raise ValueError("sample set mismatch")
                filled[row] = True
                try:
                    matrix[row] = [
                        HARD_CALL[value]
                        for value in fields[genotype_start:]
                    ]
                except KeyError as missing:
                    raise ValueError(
                        f"{fields[iid_column]}: dosage {missing.args[0]!r} is "
                        "not a hard call; fractional dosages are unsupported"
                    ) from None

    if not filled.all():
        raise ValueError("sample set mismatch")

    return matrix, tuple(emitted_keys)
