"""Shared metadata and CSV handling; no Hail or R dependency."""

import csv
import json
import math
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = Path.home() / "allxall-comparison" / "axa_test"


def settings(path=HERE / "settings.json"):
    return json.loads(Path(path).read_text())


def gene_id(value):
    return re.sub(r"^(ENSG\d+)\.\d+$", r"\1", value.strip())


def read_rows(path, delimiter=","):
    with Path(path).open(newline="") as source:
        return list(csv.DictReader(source, delimiter=delimiter))


def write_rows(path, rows, fields=None, delimiter=","):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields if fields is not None else rows[0], delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def index_rows(rows, key):
    indexed = {}
    for row in rows:
        value = key(row)
        if value in indexed:
            raise ValueError(f"Duplicate key: {value}")
        indexed[value] = row
    return indexed


def valid_p(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) and 0 <= number <= 1 else None


def p_text(value):
    # Keep the original representation, including p=0; no plotting floor here.
    return str(value) if valid_p(value) is not None else "NA"
