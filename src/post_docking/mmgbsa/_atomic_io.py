"""Atomic file writers shared by MMGBSA helpers."""

from __future__ import annotations

import csv
import json
import os
import uuid
from pathlib import Path
from typing import Any, Mapping, Sequence


def tmp_path(final_path: Path) -> Path:
    return final_path.with_name(f"{final_path.name}.tmp.{uuid.uuid4().hex}")


def write_text_atomic(
    path: Path,
    text: str,
    *,
    mkdir: bool = True,
    require_nonempty: bool = False,
) -> None:
    if mkdir:
        path.parent.mkdir(parents=True, exist_ok=True)
    staged_path = tmp_path(path)
    staged_path.write_text(text, encoding="utf-8")
    if require_nonempty and (not staged_path.exists() or staged_path.stat().st_size == 0):
        if staged_path.exists():
            staged_path.unlink()
        raise RuntimeError(f"output not created: {staged_path}")
    os.replace(staged_path, path)


def write_json_atomic(
    path: Path,
    payload: Mapping[str, Any],
    *,
    mkdir: bool = True,
    require_nonempty: bool = False,
) -> None:
    write_text_atomic(
        path,
        json.dumps(dict(payload), indent=2, sort_keys=True) + "\n",
        mkdir=mkdir,
        require_nonempty=require_nonempty,
    )


def write_csv_rows_atomic(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    *,
    mkdir: bool = True,
) -> None:
    if mkdir:
        path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    staged_path = tmp_path(path)
    with staged_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or ["empty"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    os.replace(staged_path, path)
