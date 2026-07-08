from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
import shutil
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Set, Tuple

from config.output_paths import runtime_root
from post_docking.rescoring.scorch_selection import pose_base_from_path
from post_docking.rescoring.scorch_types import StageSpec

HIDDEN_CACHE_DIRNAME = ".scorch_provisional_cache"
CACHE_SCHEMA_VERSION = 2


def _truthy(value: Any, *, default: bool = False) -> bool:
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return bool(value)
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def final_reuse_enabled(cfg: Mapping[str, Any]) -> bool:
    return _truthy(cfg.get("SCORCH_PROVISIONAL_REUSE"), default=True)


def cache_write_enabled(cfg: Mapping[str, Any]) -> bool:
    return _truthy(cfg.get("SCORCH_PROVISIONAL_CACHE_WRITE"), default=False)


def provisional_run_enabled(cfg: Mapping[str, Any]) -> bool:
    return _truthy(cfg.get("SCORCH_PROVISIONAL_ENABLE"), default=False)


def cache_root(cfg: Mapping[str, Any], run_id: str) -> Path:
    data_root = runtime_root(cfg, "DATA_DIR", "data")
    return Path(data_root).expanduser() / str(run_id) / HIDDEN_CACHE_DIRNAME


def artifact_post_root(cfg: Mapping[str, Any], run_id: str) -> Path:
    return cache_root(cfg, run_id) / "artifacts" / "post_docked"


def path_is_hidden_cache(path: Path) -> bool:
    return HIDDEN_CACHE_DIRNAME in {part for part in Path(path).parts}


def _safe_token(value: str) -> str:
    out = []
    for ch in str(value or ""):
        if ch.isalnum() or ch in {"_", "-", "."}:
            out.append(ch)
        else:
            out.append("_")
    token = "".join(out).strip("._")
    return token or "unknown"


def _stable_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _float_token(value: Any) -> str:
    if value is None:
        return ""
    try:
        return f"{float(value):.8g}"
    except Exception:
        return str(value).strip()


def _file_fingerprint(path: Optional[Path], *, hash_contents: bool) -> Dict[str, Any]:
    if path is None:
        return {"path": "", "exists": False}
    p = Path(path)
    try:
        stat = p.stat()
    except Exception:
        return {"path": str(p), "exists": False}
    payload: Dict[str, Any] = {
        "path": str(p),
        "exists": True,
        "size": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }
    if hash_contents:
        h = hashlib.sha256()
        try:
            with p.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    h.update(chunk)
            payload["sha256"] = h.hexdigest()
        except Exception:
            payload["sha256"] = ""
    return payload


def _fingerprint_matches(recorded: Any, current: Mapping[str, Any]) -> bool:
    if not isinstance(recorded, Mapping):
        return False
    if bool(recorded.get("exists")) != bool(current.get("exists")):
        return False
    if str(recorded.get("path", "")) and str(current.get("path", "")):
        if str(recorded.get("path", "")) != str(current.get("path", "")):
            return False
    if "sha256" in recorded or "sha256" in current:
        for key in ("size", "sha256"):
            if str(recorded.get(key, "")) != str(current.get(key, "")):
                return False
        return True
    for key in ("size", "mtime_ns"):
        if key in recorded or key in current:
            if str(recorded.get(key, "")) != str(current.get(key, "")):
                return False
    return True


def _scorch_identity(cfg: Mapping[str, Any]) -> Dict[str, str]:
    keys = (
        "SCORCH",
        "SCORCH_SCRIPT",
        "SCORCH_ENV",
        "SCORCH_ENV_PREFIX",
        "SCORCH_DEVICE_EFFECTIVE",
        "SCORCH_GPU_BACKEND_EFFECTIVE",
        "SCORCH_GPU_IDS_EFFECTIVE",
        "SCORCH_DEVICE_REASON",
    )
    return {key: str(cfg.get(key, "") or "") for key in keys if cfg.get(key)}


def _ligand_base_from_row(row: Mapping[str, str], *, decoy_prefix: str) -> str:
    for key in (
        "ligand_file",
        "Ligand_ID",
        "ligand",
        "Ligand",
        "ligand_id",
        "filename",
        "name",
    ):
        raw = str(row.get(key, "") or "").strip()
        if raw:
            return pose_base_from_path(Path(Path(raw).name), decoy_prefix=decoy_prefix)
    return ""


def _pose_fingerprints(
    ligands: Sequence[Path],
    *,
    decoy_prefix: str,
) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    for ligand_path in ligands:
        base = pose_base_from_path(Path(Path(ligand_path).name), decoy_prefix=decoy_prefix)
        if not base:
            continue
        out[base] = _file_fingerprint(Path(ligand_path), hash_contents=True)
    return out


def _row_key(
    *,
    run_id: str,
    combo: Tuple[str, str, str],
    spec: StageSpec,
    run_mode: str,
    decoy_prefix: str,
    base: str,
    selected_stage: str,
    selected_score: Any,
    pose_fingerprint: Mapping[str, Any],
    receptor_fingerprint: Mapping[str, Any],
    scorch_identity: Mapping[str, str],
) -> str:
    payload = {
        "schema_version": CACHE_SCHEMA_VERSION,
        "run_id": str(run_id),
        "pdb": str(combo[0]).upper(),
        "variant": str(combo[1]).upper(),
        "ph": str(combo[2]),
        "engine": str(spec.source),
        "run_mode": str(run_mode),
        "stage_dir": str(spec.stage_dir),
        "decoy_prefix": str(decoy_prefix),
        "ligand": str(base),
        "selected_stage": str(selected_stage or ""),
        "selected_score": _float_token(selected_score),
        "pose_sha256": str(pose_fingerprint.get("sha256", "")),
        "pose_size": int(pose_fingerprint.get("size", 0) or 0),
        "receptor_sha256": str(receptor_fingerprint.get("sha256", "")),
        "receptor_size": int(receptor_fingerprint.get("size", 0) or 0),
        "scorch_identity_sha256": _sha256_text(_stable_json(scorch_identity)),
    }
    return _sha256_text(_stable_json(payload))


def _task_match(
    manifest: Mapping[str, Any],
    *,
    run_id: str,
    combo: Tuple[str, str, str],
    spec: StageSpec,
    run_mode: str,
    decoy_prefix: str,
) -> bool:
    return (
        str(manifest.get("run_id", "")) == str(run_id)
        and str(manifest.get("pdb", "")).upper() == str(combo[0]).upper()
        and str(manifest.get("variant", "")).upper() == str(combo[1]).upper()
        and str(manifest.get("ph", "")) == str(combo[2])
        and str(manifest.get("engine", "")) == str(spec.source)
        and str(manifest.get("stage_dir", "")) == str(spec.stage_dir)
        and str(manifest.get("run_mode", "")) == str(run_mode)
        and str(manifest.get("decoy_prefix", "")) == str(decoy_prefix)
    )


def _read_csv_rows(path: Path) -> Tuple[list[str], list[Dict[str, str]]]:
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        rows = [dict(row) for row in reader]
    return fields, rows


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
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
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_path, path)
    finally:
        if tmp_path is not None and tmp_path.exists():
            try:
                tmp_path.unlink()
            except Exception:
                pass


def write_rows_csv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, str]]) -> None:
    ordered_fields: list[str] = []
    for field in fields:
        if field and field not in ordered_fields:
            ordered_fields.append(str(field))
    for row in rows:
        for field in row.keys():
            if field not in ordered_fields:
                ordered_fields.append(str(field))
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered_fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in ordered_fields})


def record_task_output(
    *,
    cfg: Mapping[str, Any],
    run_id: str,
    combo: Tuple[str, str, str],
    spec: StageSpec,
    run_mode: str,
    decoy_prefix: str,
    receptor: Path,
    score_csv: Optional[Path],
    allowed_bases: Optional[Set[str]],
    selected_stage_by_base: Optional[Mapping[str, str]],
    selected_score_by_base: Optional[Mapping[str, Any]],
    ligands: Sequence[Path],
    output_csv: Path,
    phase: str,
    logger: logging.Logger,
) -> Optional[Path]:
    if not output_csv.exists() or output_csv.stat().st_size <= 0:
        return None
    try:
        fields, rows = _read_csv_rows(output_csv)
    except Exception:
        logger.warning(
            "[scorch-provisional.cache] action=record status=skip reason=read_failed path=%s",
            output_csv,
            exc_info=True,
        )
        return None
    if not rows:
        return None

    selected_stage_by_base = selected_stage_by_base or {}
    selected_score_by_base = selected_score_by_base or {}
    pose_by_base = _pose_fingerprints(ligands, decoy_prefix=decoy_prefix)
    receptor_fp = _file_fingerprint(receptor, hash_contents=True)
    score_csv_fp = _file_fingerprint(score_csv, hash_contents=False)
    scorch_identity = _scorch_identity(cfg)
    row_keys: Dict[str, str] = {}
    for row in rows:
        base = _ligand_base_from_row(row, decoy_prefix=decoy_prefix)
        if not base:
            continue
        pose_fp = pose_by_base.get(base, {})
        if not pose_fp:
            continue
        row_keys[base] = _row_key(
            run_id=run_id,
            combo=combo,
            spec=spec,
            run_mode=run_mode,
            decoy_prefix=decoy_prefix,
            base=base,
            selected_stage=str(selected_stage_by_base.get(base, "")),
            selected_score=selected_score_by_base.get(base),
            pose_fingerprint=pose_fp,
            receptor_fingerprint=receptor_fp,
            scorch_identity=scorch_identity,
        )

    if not row_keys:
        return None

    identity = {
        "schema_version": CACHE_SCHEMA_VERSION,
        "run_id": str(run_id),
        "pdb": str(combo[0]).upper(),
        "variant": str(combo[1]).upper(),
        "ph": str(combo[2]),
        "engine": str(spec.source),
        "stage_dir": str(spec.stage_dir),
        "output_name": str(spec.output_name),
        "run_mode": str(run_mode),
        "decoy_prefix": str(decoy_prefix),
        "phase": str(phase or "provisional"),
        "allowed_sha256": _sha256_text(
            "\n".join(sorted(str(x) for x in (allowed_bases or set())))
        ),
        "row_key_sha256": _sha256_text("\n".join(sorted(row_keys.values()))),
    }
    digest = _sha256_text(_stable_json(identity))[:24]
    root = cache_root(cfg, run_id)
    task_dir = (
        root
        / _safe_token(str(phase or "provisional"))
        / _safe_token(str(spec.source))
        / digest
    )
    task_dir.mkdir(parents=True, exist_ok=True)
    cached_csv = task_dir / "scores.csv"
    tmp_csv: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=cached_csv.name + ".",
            suffix=".tmp",
            dir=task_dir,
            delete=False,
        ) as handle:
            tmp_csv = Path(handle.name)
            with Path(output_csv).open("rb") as src:
                shutil.copyfileobj(src, handle)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_csv, cached_csv)
    finally:
        if tmp_csv is not None and tmp_csv.exists():
            try:
                tmp_csv.unlink()
            except Exception:
                pass
    manifest: Dict[str, Any] = dict(identity)
    manifest.update(
        {
            "created_at": time.time(),
            "source_csv": str(output_csv),
            "cache_csv": str(cached_csv),
            "row_count": len(rows),
            "fields": fields,
            "allowed_bases": sorted(str(x) for x in (allowed_bases or set())),
            "selected_stage_by_base": dict(selected_stage_by_base),
            "selected_score_by_base": {
                str(k): _float_token(v) for k, v in selected_score_by_base.items()
            },
            "row_keys_by_base": row_keys,
            "pose_fingerprints_by_base": pose_by_base,
            "receptor": receptor_fp,
            "score_csv": score_csv_fp,
            "scorch_identity": scorch_identity,
        }
    )
    _atomic_write_json(task_dir / "manifest.json", manifest)
    logger.info(
        "[scorch-provisional.cache] action=record status=ok phase=%s rows=%d path=%s",
        str(phase or "provisional"),
        len(rows),
        cached_csv,
    )
    return cached_csv


def reusable_rows_for_task(
    *,
    cfg: Mapping[str, Any],
    run_id: str,
    combo: Tuple[str, str, str],
    spec: StageSpec,
    run_mode: str,
    decoy_prefix: str,
    allowed_bases: Optional[Set[str]],
    selected_stage_by_base: Optional[Mapping[str, str]],
    selected_score_by_base: Optional[Mapping[str, Any]],
    ligands: Sequence[Path],
    receptor: Path,
    score_csv: Optional[Path] = None,
    logger: logging.Logger,
) -> Tuple[list[str], list[Dict[str, str]], Set[str]]:
    wanted_bases = {str(x) for x in (allowed_bases or set()) if str(x).strip()}
    if not wanted_bases:
        return [], [], set()
    selected_stage_by_base = selected_stage_by_base or {}
    selected_score_by_base = selected_score_by_base or {}
    pose_by_base = _pose_fingerprints(ligands, decoy_prefix=decoy_prefix)
    receptor_fp = _file_fingerprint(receptor, hash_contents=True)
    score_csv_fp = _file_fingerprint(score_csv, hash_contents=False)
    scorch_identity = _scorch_identity(cfg)
    expected_keys: Dict[str, str] = {}
    for base in wanted_bases:
        pose_fp = pose_by_base.get(base)
        if not pose_fp:
            continue
        expected_keys[base] = _row_key(
            run_id=run_id,
            combo=combo,
            spec=spec,
            run_mode=run_mode,
            decoy_prefix=decoy_prefix,
            base=base,
            selected_stage=str(selected_stage_by_base.get(base, "")),
            selected_score=selected_score_by_base.get(base),
            pose_fingerprint=pose_fp,
            receptor_fingerprint=receptor_fp,
            scorch_identity=scorch_identity,
        )
    if not expected_keys:
        return [], [], set()

    root = cache_root(cfg, run_id)
    if not root.exists():
        return [], [], set()

    fields: list[str] = []
    rows_by_base: Dict[str, Dict[str, str]] = {}
    wanted_by_key = {row_key: base for base, row_key in expected_keys.items()}
    manifests: list[tuple[int, float, str, Path, Mapping[str, Any]]] = []
    for manifest_path in sorted(root.glob("*/*/*/manifest.json")):
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        if not isinstance(manifest, Mapping):
            continue
        if int(manifest.get("schema_version", 0) or 0) != CACHE_SCHEMA_VERSION:
            continue
        if not _task_match(
            manifest,
            run_id=run_id,
            combo=combo,
            spec=spec,
            run_mode=run_mode,
            decoy_prefix=decoy_prefix,
        ):
            continue
        if not _fingerprint_matches(manifest.get("receptor"), receptor_fp):
            continue
        recorded_score = manifest.get("score_csv")
        if score_csv is not None and not _fingerprint_matches(recorded_score, score_csv_fp):
            continue
        if dict(manifest.get("scorch_identity") or {}) != scorch_identity:
            continue
        phase = str(manifest.get("phase", "provisional")).strip().lower()
        phase_rank = 0 if phase == "final" else 1
        try:
            created_at = float(manifest.get("created_at") or 0.0)
        except Exception:
            created_at = 0.0
        manifests.append((phase_rank, -created_at, str(manifest_path), manifest_path, manifest))

    for _phase_rank, _created_rank, _path_key, manifest_path, manifest in sorted(manifests):
        manifest_keys = manifest.get("row_keys_by_base")
        if not isinstance(manifest_keys, Mapping):
            continue
        matching_bases = {
            str(base)
            for base, row_key in manifest_keys.items()
            if str(row_key) in wanted_by_key and str(base) == wanted_by_key[str(row_key)]
        }
        matching_bases -= set(rows_by_base.keys())
        if not matching_bases:
            continue
        csv_path = Path(str(manifest.get("cache_csv") or manifest_path.parent / "scores.csv"))
        if not csv_path.exists() or csv_path.stat().st_size <= 0:
            continue
        try:
            cached_fields, cached_rows = _read_csv_rows(csv_path)
        except Exception:
            continue
        for field in cached_fields:
            if field not in fields:
                fields.append(field)
        for row in cached_rows:
            base = _ligand_base_from_row(row, decoy_prefix=decoy_prefix)
            if base in matching_bases and base not in rows_by_base:
                rows_by_base[base] = row
        if set(rows_by_base) >= set(expected_keys):
            break

    covered = set(rows_by_base.keys())
    if covered:
        logger.info(
            "[scorch-provisional.cache] action=reuse status=hit source=%s run_mode=%s covered=%d requested=%d",
            spec.source,
            run_mode,
            len(covered),
            len(wanted_bases),
        )
    return fields, [rows_by_base[base] for base in sorted(rows_by_base)], covered


def remaining_ligands_after_reuse(
    ligands: Sequence[Path],
    covered_bases: Iterable[str],
    *,
    decoy_prefix: str,
) -> list[Path]:
    covered = {str(x) for x in covered_bases}
    out: list[Path] = []
    for ligand_path in ligands:
        base = pose_base_from_path(Path(Path(ligand_path).name), decoy_prefix=decoy_prefix)
        if base and base in covered:
            continue
        out.append(Path(ligand_path))
    return out
