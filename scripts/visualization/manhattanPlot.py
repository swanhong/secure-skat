#!/usr/bin/env python3

import argparse
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11 fallback
    import tomli as tomllib

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent

DEFAULT_INPUT_FILE = REPO_ROOT / "out/party1/assoc.txt"
DEFAULT_POS_FILE = REPO_ROOT / "example_data/party1/snp_pos.txt"
DEFAULT_GKEEP_FILE = REPO_ROOT / "cache/party1/gkeep.txt"
DEFAULT_OUTPUT_FILE = SCRIPT_DIR / "sfgwas_results.jpg"
DEFAULT_CONFIG_GLOBAL = REPO_ROOT / "config/configGlobal.toml"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Postprocess SF-GWAS association output and render a Manhattan plot."
    )
    parser.add_argument("--input-file", type=Path, default=DEFAULT_INPUT_FILE)
    parser.add_argument("--pos-file", type=Path, default=DEFAULT_POS_FILE)
    parser.add_argument("--gkeep-file", type=Path, default=DEFAULT_GKEEP_FILE)
    parser.add_argument("--output-file", type=Path, default=DEFAULT_OUTPUT_FILE)
    parser.add_argument(
        "--config-global",
        type=Path,
        default=DEFAULT_CONFIG_GLOBAL,
        help="Global TOML config used to infer total sample count and covariate count.",
    )
    parser.add_argument(
        "--num-inds",
        type=int,
        default=None,
        help="Override total number of individuals instead of reading it from config.",
    )
    parser.add_argument(
        "--num-cov",
        type=int,
        default=None,
        help="Override total covariate count (including intercept) instead of reading it from config.",
    )
    return parser.parse_args()


def load_analysis_dimensions(config_global: Path) -> tuple[int, int]:
    with config_global.open("rb") as handle:
        config = tomllib.load(handle)

    num_inds = config.get("num_inds")
    if not isinstance(num_inds, list) or len(num_inds) < 2:
        raise ValueError(f"Invalid or missing num_inds in {config_global}")

    num_ind_total = sum(int(x) for x in num_inds[1:])

    if "num_covs" not in config:
        raise ValueError(f"Missing num_covs in {config_global}")

    # SKAT prepends an intercept internally, so the downstream test uses
    # num_covs + 1 total covariates.
    num_cov = int(config["num_covs"]) + 1
    return num_ind_total, num_cov

def postprocess_assoc(
    new_assoc_file: str,
    assoc_file: str,
    pos_file: str,
    gkeep_file: str,
    num_ind_total: int,
    num_cov: int,
) -> None:
    import numpy as np
    from scipy.stats import chi2

    # new_assoc_file: Name of new assoc file (processed)
    # assoc_file: Name of original assoc file
    # pos_file: Path to pos.txt
    # gkeep_file: Path to gkeep.txt
    # num_ind_total: Total number of individuals
    # num_cov: Number of covariates

    # Load SNP filter
    gkeep = np.loadtxt(gkeep_file, dtype=bool)

    # Load and check dimension of output association stats
    assoc = np.loadtxt(assoc_file)
    assert len(assoc) == gkeep.sum()

    # Calculate p-values
    t2 = (assoc**2) * (num_ind_total - num_cov) / (1 - assoc**2 + 1e-10)
    log10p = np.log10(chi2.sf(t2, df=1))

    # Append SNP position information and write to a new file
    lineno = 0
    assoc_idx = 0

    with open(new_assoc_file, "w") as out:
        out.write("\t".join(["#CHROM", "POS", "R", "LOG10P"]) + "\n")

        for line in open(pos_file):
            pos = line.strip().split()

            if gkeep[lineno]:
                out.write(pos[0] + "\t" + pos[1] + "\t" + str(assoc[assoc_idx]) + "\t" + str(log10p[assoc_idx]) + "\n")
                assoc_idx += 1

            lineno += 1


def plot_assoc(plot_file: str, new_assoc_file: str) -> None:
    import matplotlib.pyplot as plt
    import pandas as pd
    from qmplot import manhattanplot

    # Load postprocessed assoc file and convert p-values
    tab = pd.read_table(new_assoc_file)
    tab["P"] = 10 ** tab["LOG10P"]

    # Create a Manhattan plot
    plt.figure()
    manhattanplot(
        data=tab,
        suggestiveline=None,  # type: ignore
        genomewideline=None,  # type: ignore
        marker=".",
        xticklabel_kws={"rotation": "vertical"},  # set vertical or any other degrees as you like.
    )
    plt.tight_layout()
    plt.savefig(plot_file)

def main():
    args = parse_args()
    num_inds_cfg, num_cov_cfg = load_analysis_dimensions(args.config_global)
    num_inds = args.num_inds if args.num_inds is not None else num_inds_cfg
    num_cov = args.num_cov if args.num_cov is not None else num_cov_cfg

    print("Plotting script called...")
    print(f"SF-GWAS output file: {args.input_file}")
    print(f"SNP position file: {args.pos_file}")
    print(f"QC filter file: {args.gkeep_file}")
    print(f"Global config file: {args.config_global}")
    print(f"Number of individuals: {num_inds}")
    print(f"Number of covariates: {num_cov}")
    
    processed_input = args.input_file.with_name(args.input_file.name + ".processed")
    postprocess_assoc(processed_input, args.input_file, args.pos_file, args.gkeep_file, num_inds, num_cov)
    plot_assoc(args.output_file, processed_input)

    print(f"Plot saved to {args.output_file}")

if __name__ == "__main__":
    main()
