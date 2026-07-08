from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from analysis._common import read_csv_rows, write_csv_rows


def load_config(path: str | Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    config_path = Path(path)
    if not config_path.exists():
        return {}
    text = config_path.read_text(encoding="utf-8")
    if config_path.suffix.lower() in {".yaml", ".yml"}:
        try:
            import yaml  # type: ignore[import-untyped]
        except Exception as exc:  # pragma: no cover - depends on optional environment
            raise RuntimeError("YAML config requested but PyYAML is not installed") from exc
        loaded = yaml.safe_load(text)
    else:
        loaded = json.loads(text)
    return loaded if isinstance(loaded, dict) else {}


def read_table(path: str | Path) -> list[dict[str, str]]:
    return read_csv_rows(Path(path))


def write_table(path: str | Path, rows: list[Mapping[str, Any]], fieldnames: list[str] | None = None) -> None:
    write_csv_rows(Path(path), rows, fieldnames)


def write_json(path: str | Path, payload: Any) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def output_dirs(config: Mapping[str, Any]) -> tuple[Path, Path, Path]:
    raw_outputs = config.get("outputs")
    outputs: Mapping[str, Any] = raw_outputs if isinstance(raw_outputs, Mapping) else {}
    analysis_dir = Path(str(outputs.get("analysis_dir", "outputs/analysis")))
    figures_dir = Path(str(outputs.get("figures_dir", "outputs/figures")))
    models_dir = Path(str(outputs.get("models_dir", "outputs/models")))
    return analysis_dir, figures_dir, models_dir
