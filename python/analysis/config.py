from pathlib import Path
import tomllib


def load_party_config(directory: Path, party_id: int = 1) -> dict:
    config = {}
    for name in (
        "configGlobal.toml",
        f"configLocal.Party{party_id}.toml",
    ):
        with (directory / name).open("rb") as config_file:
            config.update(tomllib.load(config_file))
    return config
