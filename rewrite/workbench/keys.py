import secrets
from pathlib import Path


KEY_BYTES = 32
KEY_FILES = (
    "shared_key_global.bin",
    "shared_key_0_1.bin",
    "shared_key_0_2.bin",
    "shared_key_1_2.bin",
)


def write_shared_keys(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=False)
    for name in KEY_FILES:
        (output_dir / name).write_bytes(secrets.token_bytes(KEY_BYTES))
