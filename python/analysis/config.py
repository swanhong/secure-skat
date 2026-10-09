from pathlib import Path
import tomllib


def load_party_config(directory: Path, party_id: int = 1) -> dict:
    config = {}
    local_name = "configLocal.toml"
    if not (directory / local_name).exists():
        local_name = f"configLocal.Party{party_id}.toml"
    for name in (
        "configGlobal.toml",
        local_name,
    ):
        with (directory / name).open("rb") as config_file:
            config.update(tomllib.load(config_file))
    return config
