from __future__ import annotations

import json
import logging
import os
import re
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

_CHUNK_ID_SANITIZER = re.compile(r"[^A-Za-z0-9_.-]+")


def sanitize_chunk_id(raw: str) -> str:
    token = _CHUNK_ID_SANITIZER.sub("_", str(raw or "").strip())
    token = token.strip("._-")
    if not token:
        return "chunk"
    return token[:160]


def completion_marker_path(
    stage_dir: Path,
    *,
    engine: str,
    chunk_id: Optional[str] = None,
) -> Path:
    if chunk_id:
        safe_chunk = sanitize_chunk_id(str(chunk_id))
        return stage_dir / f"completion_{engine}__chunk_{safe_chunk}.json"
    return stage_dir / f"completion_{engine}.json"


def write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=path.name + ".",
            suffix=".tmp",
            dir=path.parent,
            delete=False,
        ) as handle:
            tmp_path = Path(handle.name)
            handle.write(json.dumps(payload, indent=2, sort_keys=False))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_path, path)
    finally:
        if tmp_path is not None and tmp_path.exists():
            try:
                tmp_path.unlink()
            except Exception:
                pass


def _read_marker(path: Path) -> Optional[dict[str, Any]]:
    try:
        loaded = json.loads(path.read_text(encoding="utf-8")) or {}
    except Exception:
        return None
    if not isinstance(loaded, dict):
        return None
    return loaded


def _to_int(value: Any) -> int:
    try:
        return int(value)
    except Exception:
        return 0


def aggregate_stage_chunk_markers(
    stage_dir: Path,
    *,
    engine: str,
    logger: Optional[logging.Logger] = None,
) -> Optional[dict[str, Any]]:
    if logger is None:
        logger = logging.getLogger("completion.aggregate")
    chunk_markers = sorted(stage_dir.glob(f"completion_{engine}__chunk_*.json"))
    if not chunk_markers:
        return None

    expected_ligands: set[str] = set()
    missing_ligands: set[str] = set()
    failure_markers: dict[str, str] = {}
    expected_count_sum = 0
    missing_before_sum = 0
    missing_after_sum = 0
    missing_after_strict_sum = 0
    run_id = ""
    pdb_id = ""
    variant = ""
    ph_label = ""
    stage_name = stage_dir.name
    valid_sources = 0

    for marker_path in chunk_markers:
        payload = _read_marker(marker_path)
        if payload is None:
            continue
        valid_sources += 1
        if not run_id:
            run_id = str(payload.get("run_id") or "").strip()
        if not pdb_id:
            pdb_id = str(payload.get("pdb_id") or "").strip()
        if not variant:
            variant = str(payload.get("variant") or "").strip()
        if not ph_label:
            ph_label = str(payload.get("ph_label") or "").strip()
        if not stage_name:
            stage_name = str(payload.get("stage") or "").strip()

        expected_count_sum += _to_int(payload.get("expected_count"))
        missing_before_sum += _to_int(payload.get("missing_count_before"))
        missing_after_sum += _to_int(payload.get("missing_count_after"))
        missing_after_strict_sum += _to_int(payload.get("missing_count_after_strict"))

        expected_raw = payload.get("expected_ligands")
        if isinstance(expected_raw, list):
            for item in expected_raw:
                token = str(item).strip()
                if token:
                    expected_ligands.add(token)

        missing_raw = payload.get("missing_ligands_after")
        if isinstance(missing_raw, list):
            for item in missing_raw:
                token = str(item).strip()
                if token:
                    missing_ligands.add(token)

        failed_raw = payload.get("failure_markers")
        if isinstance(failed_raw, dict):
            for key, value in failed_raw.items():
                k = str(key).strip()
                v = str(value).strip()
                if k and v:
                    failure_markers[k] = v

    if valid_sources <= 0:
        return None

    expected_count = len(expected_ligands) if expected_ligands else expected_count_sum
    missing_count = len(missing_ligands)
    aggregate_payload: dict[str, Any] = {
        "engine": str(engine),
        "stage": str(stage_name),
        "pdb_id": str(pdb_id),
        "variant": str(variant),
        "ph_label": str(ph_label),
        "run_id": str(run_id),
        "expected_count": int(expected_count),
        "missing_count_before": int(missing_before_sum),
        "missing_count_after": int(missing_count),
        "missing_count_after_strict": int(min(missing_after_strict_sum, missing_count)),
        "missing_ligands_after": sorted(missing_ligands),
        "expected_ligands": sorted(expected_ligands),
        "failure_markers": failure_markers,
        "chunk_marker_count": int(valid_sources),
        "sources": [p.name for p in chunk_markers],
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "success": bool(missing_count == 0),
        "marker_scope": "aggregate",
    }
    out_path = completion_marker_path(stage_dir, engine=engine, chunk_id=None)
    write_json_atomic(out_path, aggregate_payload)
    logger.info(
        "[completion.aggregate] stage=%s engine=%s chunk_markers=%d expected=%d missing=%d out=%s",
        stage_dir.name,
        str(engine),
        int(valid_sources),
        int(expected_count),
        int(missing_count),
        out_path,
    )
    return aggregate_payload


def aggregate_combo_chunk_markers(
    combo_dir: Path,
    *,
    engines: Sequence[str] = ("vina", "gnina", "ledock"),
    logger: Optional[logging.Logger] = None,
) -> int:
    if logger is None:
        logger = logging.getLogger("completion.aggregate")
    if not combo_dir.exists():
        return 0
    aggregated = 0
    for stage_dir in sorted(combo_dir.iterdir()):
        if not stage_dir.is_dir():
            continue
        for engine in engines:
            try:
                payload = aggregate_stage_chunk_markers(
                    stage_dir,
                    engine=str(engine),
                    logger=logger,
                )
            except Exception:
                logger.warning(
                    "[completion.aggregate] action=skip stage=%s engine=%s reason=exception",
                    stage_dir,
                    str(engine),
                    exc_info=True,
                )
                continue
            if payload is not None:
                aggregated += 1
    return int(aggregated)
