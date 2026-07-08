from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Mapping

import pandas as pd


def read_source_table(path: str | Path, sheet_name: str | int | None = None) -> pd.DataFrame:
    table_path = Path(path)
    suffix = table_path.suffix.lower()
    if suffix in {".xlsx", ".xlsm", ".xls"}:
        return pd.read_excel(table_path, sheet_name=sheet_name or 0)
    if suffix == ".parquet":
        return pd.read_parquet(table_path)
    sep = "\t" if suffix in {".tsv", ".tab"} else ","
    comment = "#" if suffix in {".tsv", ".tab"} else None
    if suffix in {".tsv", ".tab"}:
        try:
            first_line = table_path.open("r", encoding="utf-8", errors="replace").readline()
        except OSError:
            first_line = ""
        if "GtoPdb Version" in first_line:
            return pd.read_csv(table_path, sep=sep, skiprows=1)
    return pd.read_csv(table_path, sep=sep, comment=comment)


def write_source_manifest(path: str | Path, payload: Mapping[str, Any]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(payload), indent=2, sort_keys=True), encoding="utf-8")


def download_to_cache(
    url: str,
    out_path: str | Path,
    *,
    retries: int = 3,
    sleep_sec: float = 2.0,
    overwrite: bool = False,
) -> Path:
    out = Path(out_path)
    if out.exists() and not overwrite:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "AtlasAnalysis/1.0"})
    last_error: Exception | None = None
    for attempt in range(max(1, retries)):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                out.write_bytes(response.read())
            return out
        except (urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
            if attempt + 1 < retries:
                time.sleep(max(0.0, sleep_sec) * float(attempt + 1))
    raise RuntimeError(f"failed to download {url} -> {out}: {last_error}")
