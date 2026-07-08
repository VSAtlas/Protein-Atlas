"""Report writers for the protein-prep benchmark."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Mapping, Sequence

BenchmarkRow = Mapping[str, object]


def write_outputs(rows: Sequence[BenchmarkRow], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(list(rows), handle, indent=2, sort_keys=True)
    if not rows:
        return
    with (out_dir / "summary.csv").open("w", encoding="utf-8", newline="") as handle:
        fieldnames = list(rows[0].keys())
        for row in rows[1:]:
            for key in row:
                if key not in fieldnames:
                    fieldnames.append(key)
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


__all__ = ["write_outputs"]
