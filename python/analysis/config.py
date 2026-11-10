from pathlib import Path
import tomllib


def normalize_chromosomes(values: list[int | str]) -> list[int]:
    if not isinstance(values, list) or not values:
        raise ValueError('chromosomes must be a nonempty integer list or ["all"]')
    if values == ["all"]:
        return list(range(1, 23))
    if any(type(value) is not int or not 1 <= value <= 22 for value in values):
        raise ValueError('chromosomes must contain integers from 1 to 22 or exactly ["all"]')
    if len(set(values)) != len(values):
        raise ValueError("chromosomes must not contain duplicates")
    return values.copy()


def analysis_groups(config: dict) -> list[tuple[str, float | None, str]]:
    """(ancestry, global_max_maf, output subdirectory) per analysis, matching Go's analysisLabel."""
    thresholds = config.get("global_max_maf") or [None]
    return [
        (ancestry, threshold, ancestry if threshold is None else f"{ancestry}/maf{threshold!r}")
        for ancestry in config["ancestries"]
        for threshold in thresholds
    ]


def load_party_config(directory: Path, party_id: int = 1) -> dict:
    config = {}
    local_name = "configLocal.toml"
    if not (directory / local_name).exists():
        local_name = f"configLocal.Party{party_id}.toml"
    for name in (
        "configGlobal.toml",
        local_name,
        "configPrepare.toml",
    ):
        with (directory / name).open("rb") as config_file:
            config.update(tomllib.load(config_file))
    config["chromosomes"] = normalize_chromosomes(config["chromosomes"])
    config["ancestries"] = [ancestry.strip().upper() for ancestry in config["ancestries"]]
    return config
