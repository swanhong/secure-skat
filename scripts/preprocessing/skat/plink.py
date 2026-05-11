from __future__ import annotations

import argparse
from pathlib import Path

from cloud import cloud_cp, cloud_ls
from utils import basename_from_path_text, copy_local_file, is_gcs_uri, pfile_path, resolve_path, run


def materialize_input_file(path_text: str, raw_dir: Path, args: argparse.Namespace, default_name: str) -> Path:
    raw_dir.mkdir(parents=True, exist_ok=True)
    dst = raw_dir / basename_from_path_text(path_text, default_name)
    if is_gcs_uri(path_text):
        if not dst.exists() or args.force:
            cloud_cp(path_text, dst, args)
    else:
        copy_local_file(resolve_path(path_text), dst)
    return dst


def materialize_optional_index_files(path_text: str, dst_file: Path, args: argparse.Namespace) -> None:
    for suffix in (".tbi", ".csi"):
        index_text = path_text + suffix
        dst = Path(str(dst_file) + suffix)
        if is_gcs_uri(path_text):
            if cloud_ls(index_text, args):
                cloud_cp(index_text, dst, args)
            continue

        index_src = resolve_path(index_text)
        if index_src.exists():
            copy_local_file(index_src, dst)


def materialize_raw_pgen_prefix(args: argparse.Namespace, raw_dir: Path) -> Path:
    prefix_text = args.pgen_prefix
    prefix_name = Path(prefix_text.rstrip("/")).name
    if not prefix_name:
        prefix_name = f"pgen.chr{args.chromosome}"
    raw_prefix = raw_dir / prefix_name

    raw_dir.mkdir(parents=True, exist_ok=True)
    required_exts = [".pgen", ".psam"]

    if is_gcs_uri(prefix_text):
        for ext in required_exts:
            dst = pfile_path(raw_prefix, ext)
            if not dst.exists() or args.force:
                cloud_cp(prefix_text + ext, dst, args)

        pvar_dst = pfile_path(raw_prefix, ".pvar")
        pvar_zst_dst = pfile_path(raw_prefix, ".pvar.zst")
        if args.force or (not pvar_dst.exists() and not pvar_zst_dst.exists()):
            if cloud_ls(prefix_text + ".pvar", args):
                pvar_zst_dst.unlink(missing_ok=True)
                cloud_cp(prefix_text + ".pvar", pvar_dst, args)
            elif cloud_ls(prefix_text + ".pvar.zst", args):
                pvar_dst.unlink(missing_ok=True)
                cloud_cp(prefix_text + ".pvar.zst", pvar_zst_dst, args)
            else:
                raise FileNotFoundError(f"Could not find {prefix_text}.pvar or {prefix_text}.pvar.zst")
    else:
        src_prefix = resolve_path(prefix_text)
        for ext in required_exts:
            copy_local_file(pfile_path(src_prefix, ext), pfile_path(raw_prefix, ext))

        pvar_src = pfile_path(src_prefix, ".pvar")
        pvar_zst_src = pfile_path(src_prefix, ".pvar.zst")
        if pvar_src.exists():
            pfile_path(raw_prefix, ".pvar.zst").unlink(missing_ok=True)
            copy_local_file(pvar_src, pfile_path(raw_prefix, ".pvar"))
        elif pvar_zst_src.exists():
            pfile_path(raw_prefix, ".pvar").unlink(missing_ok=True)
            copy_local_file(pvar_zst_src, pfile_path(raw_prefix, ".pvar.zst"))
        else:
            raise FileNotFoundError(f"Could not find {pvar_src} or {pvar_zst_src}")

    print(f"Raw PGEN prefix ready: {raw_prefix}")
    return raw_prefix


def materialize_raw_vcf_as_pgen(args: argparse.Namespace, raw_dir: Path) -> Path:
    vcf_path = materialize_input_file(args.vcf, raw_dir, args, f"input.chr{args.chromosome}.vcf.gz")
    materialize_optional_index_files(args.vcf, vcf_path, args)

    raw_prefix = raw_dir / f"vcf_chr{args.chromosome}"
    keep_path = None
    if args.vcf_keep:
        keep_path = materialize_input_file(args.vcf_keep, raw_dir, args, f"vcf_chr{args.chromosome}.keep.txt")

    return run_vcf_to_pgen(
        vcf_path=vcf_path,
        out_prefix=raw_prefix,
        chromosome=args.chromosome,
        plink2=args.plink2,
        new_id_max_allele_len=args.new_id_max_allele_len,
        vcf_no_double_id=args.vcf_no_double_id,
        keep_path=keep_path,
        force=args.force,
    )


def run_vcf_to_pgen(
    *,
    vcf_path: Path,
    out_prefix: Path,
    chromosome: int,
    plink2: str,
    new_id_max_allele_len: int,
    vcf_no_double_id: bool = False,
    keep_path: Path | None = None,
    force: bool = False,
) -> Path:
    raw_prefix = out_prefix
    if pfile_path(raw_prefix, ".pgen").exists() and not force:
        print(f"Reusing existing VCF-derived PGEN prefix: {raw_prefix}")
        return raw_prefix

    if force:
        for pattern in (raw_prefix.name + ".*", raw_prefix.name + "-temporary.*"):
            for path in raw_prefix.parent.glob(pattern):
                path.unlink()

    raw_prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [plink2, "--vcf", str(vcf_path)]
    if not vcf_no_double_id:
        cmd.append("--double-id")
    if keep_path:
        cmd.extend(["--keep", str(keep_path)])
    cmd.extend(
        [
            "--chr",
            str(chromosome),
            "--max-alleles",
            "2",
            "--set-all-var-ids",
            "@:#:$r:$a",
            "--new-id-max-allele-len",
            str(new_id_max_allele_len),
            "--make-pgen",
            "vzs",
            "--out",
            str(raw_prefix),
        ]
    )
    run(cmd)

    print(f"Raw VCF converted to PGEN prefix: {raw_prefix}")
    return raw_prefix


def materialize_raw_pgen(args: argparse.Namespace, raw_dir: Path) -> Path:
    if args.pgen_prefix:
        return materialize_raw_pgen_prefix(args, raw_dir)
    return materialize_raw_vcf_as_pgen(args, raw_dir)


def detect_pfile_vzs(prefix: Path) -> bool:
    return pfile_path(prefix, ".pvar.zst").exists()


def run_plink_source(args: argparse.Namespace, raw_prefix: Path, work_dir: Path) -> Path:
    source_prefix = work_dir / f"source_chr{args.chromosome}_biallelic"
    if pfile_path(source_prefix, ".pgen").exists() and not args.force:
        print(f"Reusing existing PLINK source prefix: {source_prefix}")
        return source_prefix

    if args.force:
        for path in source_prefix.parent.glob(source_prefix.name + ".*"):
            path.unlink()

    cmd = [args.plink2, "--pfile", str(raw_prefix)]
    if detect_pfile_vzs(raw_prefix):
        cmd.append("vzs")
    cmd.extend(
        [
            "--keep",
            str(work_dir / "sample_keep_all.txt"),
            "--chr",
            str(args.chromosome),
            "--max-alleles",
            "2",
            "--set-all-var-ids",
            "@:#:$r:$a",
            "--new-id-max-allele-len",
            str(args.new_id_max_allele_len),
            "--make-pgen",
            "--out",
            str(source_prefix),
        ]
    )
    run(cmd)

    freq_cmd = [
        args.plink2,
        "--pfile",
        str(source_prefix),
        "--freq",
        "--out",
        str(source_prefix),
    ]
    run(freq_cmd)
    return source_prefix


def materialize_party_blocks(args: argparse.Namespace, source_prefix: Path, out_dataset: Path, windows: list[dict[str, object]]) -> None:
    for party_name in ("party1", "party2"):
        party_dir = out_dataset / party_name
        for row in windows:
            block_id = int(row["block_id"])
            out_prefix = party_dir / "geno" / f"chr{block_id}"
            cmd = [
                args.plink2,
                "--pfile",
                str(source_prefix),
                "--keep",
                str(party_dir / "sample_keep.txt"),
                "--extract",
                str(row["extract_path"]),
                "--indiv-sort",
                "none",
                "--make-pgen",
                "--out",
                str(out_prefix),
            ]
            run(cmd)
