from __future__ import annotations

import csv
import fcntl
import hashlib
import json
import logging
import math
import os
import re
import shutil
import tempfile
import threading
import time
from concurrent.futures import FIRST_COMPLETED, Future, ThreadPoolExecutor, wait
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

from config.output_paths import runtime_root
from post_docking.rescoring.scorch_orchestration_types import ComboKey, OrchestrationDeps
from post_docking.rescoring.scorch_types import ScorchTask, StageSpec

SCORCH_SHARD_PLAN_POLICY = "score_based_consensus_selection_v2"
_STAGE_RESULTS_KEY = "_ATLAS_SCORCH_STAGE_RESULTS"
_TIMEOUT_OVERRIDES_KEY = "_ATLAS_SCORCH_SHARD_TIMEOUT_OVERRIDES"


@dataclass
class ShardedExecutionResult:
    failed_jobs: int
    completed: int
    total: int
    all_complete: bool
    combo_failed_local: Dict[ComboKey, bool]


def _truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return bool(value)
    return str(value or "").strip().lower() in {"1", "true", "yes", "on", "y"}


def _falsey(value: Any) -> bool:
    if isinstance(value, bool):
        return not bool(value)
    return str(value or "").strip().lower() in {
        "0",
        "false",
        "no",
        "off",
        "n",
        "disabled",
        "none",
    }


def _positive_int_token(value: Any) -> Optional[int]:
    try:
        parsed = int(str(value or "").strip())
    except Exception:
        return None
    return parsed if parsed > 0 else None


def _float_cfg_or_env(
    cfg: Mapping[str, Any],
    *names: str,
    default: float,
) -> float:
    for name in names:
        raw = os.environ.get(name)
        if raw is None:
            raw = cfg.get(name)
        if raw is None:
            continue
        try:
            value = float(str(raw).strip())
        except Exception:
            continue
        if value >= 0.0:
            return float(value)
    return float(default)


def _int_cfg_or_env(
    cfg: Mapping[str, Any],
    *names: str,
    default: int,
) -> int:
    for name in names:
        raw = os.environ.get(name)
        if raw is None:
            raw = cfg.get(name)
        if raw is None:
            continue
        try:
            value = int(str(raw).strip())
        except Exception:
            continue
        if value >= 0:
            return int(value)
    return int(default)


def _adaptive_scorch_timeout_sec(ligand_count: int) -> float:
    count = max(1, int(ligand_count or 0))
    if count <= 1:
        return 300.0
    if count <= 16:
        return 600.0
    return float(min(2400.0, max(600.0, 120.0 + (12.0 * count))))


def _stage_result_for(cfg: Mapping[str, Any], shard_id: str) -> Dict[str, Any]:
    raw = cfg.get(_STAGE_RESULTS_KEY)
    if not isinstance(raw, Mapping):
        return {}
    payload = raw.get(str(shard_id))
    return dict(payload) if isinstance(payload, Mapping) else {}


def _timeout_for_attempt(
    cfg: Mapping[str, Any],
    record: Mapping[str, Any],
    *,
    attempt: int,
) -> float:
    configured_timeout_sec: Optional[float] = None
    for name in (
        "ATLAS_SCORCH_SUBPROCESS_TIMEOUT_SEC",
        "SCORCH_SUBPROCESS_TIMEOUT_SEC",
        "ATLAS_SCORCH_SHARD_TIMEOUT_SEC",
        "SCORCH_SHARD_TIMEOUT_SEC",
    ):
        raw = os.environ.get(name)
        if raw is None:
            raw = cfg.get(name)
        if raw is None:
            continue
        try:
            configured_timeout_sec = float(str(raw).strip())
        except Exception:
            continue
        if configured_timeout_sec > 0.0:
            break
        configured_timeout_sec = None
    try:
        allowed_count = int(record.get("allowed_count") or 0)
    except Exception:
        allowed_count = 0
    base_timeout = (
        float(configured_timeout_sec)
        if configured_timeout_sec is not None
        else _adaptive_scorch_timeout_sec(allowed_count)
    )
    if int(attempt) <= 1:
        return float(base_timeout)
    multiplier = max(
        1.0,
        _float_cfg_or_env(
            cfg,
            "ATLAS_SCORCH_SHARD_RETRY_TIMEOUT_MULTIPLIER",
            "SCORCH_SHARD_RETRY_TIMEOUT_MULTIPLIER",
            default=2.0,
        ),
    )
    max_timeout = max(
        base_timeout,
        _float_cfg_or_env(
            cfg,
            "ATLAS_SCORCH_SHARD_RETRY_TIMEOUT_MAX_SEC",
            "SCORCH_SHARD_RETRY_TIMEOUT_MAX_SEC",
            default=1800.0,
        ),
    )
    return float(min(max_timeout, base_timeout * (multiplier ** max(0, int(attempt) - 1))))


def _timeout_split_min(cfg: Mapping[str, Any]) -> int:
    return max(
        2,
        _int_cfg_or_env(
            cfg,
            "ATLAS_SCORCH_SHARD_TIMEOUT_SPLIT_MIN",
            "SCORCH_SHARD_TIMEOUT_SPLIT_MIN",
            default=2,
        ),
    )


def _timeout_ligand_fanout_max(cfg: Mapping[str, Any]) -> int:
    return max(
        0,
        _int_cfg_or_env(
            cfg,
            "ATLAS_SCORCH_SHARD_TIMEOUT_LIGAND_FANOUT_MAX",
            "SCORCH_SHARD_TIMEOUT_LIGAND_FANOUT_MAX",
            default=32,
        ),
    )


def _timeout_split_candidate_count(record: Mapping[str, Any]) -> int:
    expected_bases = _expected_output_bases(record)
    if expected_bases:
        return int(len(expected_bases))
    try:
        return max(0, int(record.get("allowed_count") or 0))
    except Exception:
        return 0


def _children_by_parent(
    records: Sequence[Mapping[str, Any]],
) -> Dict[str, List[Mapping[str, Any]]]:
    children: Dict[str, List[Mapping[str, Any]]] = {}
    for record in records:
        parent = str(record.get("parent_shard_id", "") or "").strip()
        if not parent:
            continue
        children.setdefault(parent, []).append(record)
    return children


def shard_mode_enabled(cfg: Mapping[str, Any]) -> bool:
    for name in ("ATLAS_SCORCH_SHARDS_ENABLE", "SCORCH_SHARDS_ENABLE"):
        if os.environ.get(name) is not None:
            return _truthy(os.environ.get(name))
        if cfg.get(name) is not None:
            return _truthy(cfg.get(name))
    mode = str(
        os.environ.get("ATLAS_DISTRIBUTED_MODE")
        or cfg.get("ATLAS_DISTRIBUTED_MODE")
        or cfg.get("DISTRIBUTED_MODE")
        or ""
    ).strip().lower().replace("-", "_")
    if mode in {"slurm", "array", "slurm_array", "distributed", "distributed_combo"}:
        return True
    if mode and not _falsey(mode):
        return True
    for name in ("SLURM_ARRAY_TASK_COUNT", "ATLAS_DIST_TASK_COUNT"):
        value = _positive_int_token(os.environ.get(name))
        if value is not None and value > 1:
            return True
    value = _positive_int_token(cfg.get("ATLAS_DIST_TASK_COUNT"))
    if value is not None and value > 1:
        return True
    return False


def _distributed_root(cfg: Mapping[str, Any], run_id: str) -> Path:
    return runtime_root(cfg, "MANIFESTS_DIR", "manifests") / str(run_id) / "distributed"


def _plans_dir(cfg: Mapping[str, Any], run_id: str) -> Path:
    return _distributed_root(cfg, run_id) / "scorch_shard_plans"


def _claims_dir(cfg: Mapping[str, Any], run_id: str) -> Path:
    return _distributed_root(cfg, run_id) / "scorch_shard_claims"


def _results_dir(cfg: Mapping[str, Any], run_id: str) -> Path:
    return _distributed_root(cfg, run_id) / "scorch_shard_results"


def _ligand_timeout_dir(cfg: Mapping[str, Any], run_id: str) -> Path:
    return _distributed_root(cfg, run_id) / "scorch_ligand_timeouts"


def _summary_path(cfg: Mapping[str, Any], run_id: str) -> Path:
    return _distributed_root(cfg, run_id) / "scorch_shard_summary.json"


def _safe_token(value: Any) -> str:
    token = str(value or "base").strip().replace("/", "_").replace("\\", "_")
    return "".join(ch if ch.isalnum() or ch in {"_", "-", "."} else "_" for ch in token) or "base"


def _plan_path(cfg: Mapping[str, Any], run_id: str, combo: ComboKey, decoy_prefix: str) -> Path:
    name = "__".join([_safe_token(part) for part in (*combo, decoy_prefix)])
    return _plans_dir(cfg, run_id) / f"{name}.json"


def _plan_lock_path(path: Path) -> Path:
    return path.with_name(f"{path.name}.lock")


def _claim_path(cfg: Mapping[str, Any], run_id: str, shard_id: str) -> Path:
    return _claims_dir(cfg, run_id) / f"{str(shard_id)}.claim.json"


def _result_path(cfg: Mapping[str, Any], run_id: str, shard_id: str) -> Path:
    return _results_dir(cfg, run_id) / f"{str(shard_id)}.json"


def _ligand_timeout_path(
    cfg: Mapping[str, Any],
    run_id: str,
    *,
    run_mode: str,
    ligand_base: str,
) -> Path:
    key = f"{str(run_mode or '')}|{str(ligand_base or '')}"
    digest = hashlib.sha1(key.encode("utf-8")).hexdigest()[:16]
    name = "__".join([_safe_token(run_mode or "base"), _safe_token(ligand_base), digest])
    return _ligand_timeout_dir(cfg, run_id) / f"{name}.json"


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=path.name + ".",
            suffix=".tmp",
            dir=str(path.parent),
            delete=False,
        ) as handle:
            tmp_path = Path(handle.name)
            json.dump(payload, handle, sort_keys=True)
            handle.flush()
            os.fsync(handle.fileno())
        tmp_path.replace(path)
        return path
    finally:
        try:
            if tmp_path is not None and tmp_path.exists():
                tmp_path.unlink()
        except Exception:
            pass


def _read_json(path: Path) -> Optional[Dict[str, Any]]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
    return payload if isinstance(payload, dict) else None


def _csv_data_rows(path: Path) -> int:
    try:
        with path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            return sum(1 for _ in reader)
    except Exception:
        return 0


_RUN_MODE_STAGE_SUFFIX_RE = re.compile(
    r"_(?:dud|decoy|decoys|fda|prod|production)_(?:stage|pose)\d+$",
    re.IGNORECASE,
)
_STAGE_SUFFIX_RE = re.compile(r"_(?:stage|pose)\d+$", re.IGNORECASE)


def _ligand_key(value: Any) -> str:
    token = Path(str(value or "").strip()).name
    if not token:
        return ""
    for suffix in (".pdbqt", ".pdb", ".sdf", ".mol2", ".csv"):
        if token.lower().endswith(suffix):
            token = token[: -len(suffix)]
            break
    token = token.replace(".sanitized", "")
    token = _RUN_MODE_STAGE_SUFFIX_RE.sub("", token)
    return _STAGE_SUFFIX_RE.sub("", token).strip().lower()


def _row_ligand_key(row: Mapping[str, Any]) -> str:
    for field in ("Ligand_ID", "ligand", "ligand_file", "ligand_id"):
        key = _ligand_key(row.get(field))
        if key:
            return key
    return ""


def _csv_ligand_bases(path: Path) -> Set[str]:
    bases: Set[str] = set()
    try:
        with path.open("r", newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                key = _row_ligand_key(row)
                if key:
                    bases.add(key)
    except Exception:
        return set()
    return bases


def _task_output_path(post_root: Path, task: ScorchTask, decoy_prefix: str) -> Path:
    spec, combo, _receptor, _allowed, _ctrl, run_mode, _dirs, _score, _stage, _sel_score, chunk_tag = task
    pdb_id, variant, ph = combo
    output_name = spec.output_name if run_mode == "fda" else f"{decoy_prefix}_{spec.output_name}"
    if chunk_tag:
        output_name = output_name.replace(".csv", f".{chunk_tag}.csv")
    return post_root / pdb_id / variant / ph / output_name


def _stable_hash(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha1(encoded).hexdigest()[:24]


def task_to_shard_record(
    *,
    run_id: str,
    task: ScorchTask,
    post_root: Path,
    decoy_prefix: str,
    top_fraction: float,
) -> Dict[str, Any]:
    (
        spec,
        combo,
        receptor,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
    ) = task
    identity = {
        "run_id": str(run_id),
        "combo": list(combo),
        "source": spec.source,
        "stage_dir": spec.stage_dir,
        "run_mode": str(run_mode),
        "chunk_tag": str(chunk_tag or "all"),
        "decoy_prefix": str(decoy_prefix),
        "top_fraction": float(top_fraction),
        "allowed_hash": _stable_hash({"allowed": sorted(str(x) for x in allowed_bases)}),
        "receptor": str(receptor),
    }
    shard_id = f"scorch_{_stable_hash(identity)}"
    out_path = _task_output_path(post_root, task, decoy_prefix)
    return {
        "schema": "atlas.scorch_shard.v1",
        "shard_id": shard_id,
        "identity": identity,
        "spec": {
            "source": spec.source,
            "stage_dir": spec.stage_dir,
            "output_name": spec.output_name,
        },
        "combo": list(combo),
        "receptor": str(receptor),
        "allowed_bases": sorted(str(x) for x in allowed_bases),
        "control_bases": sorted(str(x) for x in control_bases),
        "run_mode": str(run_mode),
        "stage_dirs_override": (
            None
            if stage_dirs_override is None
            else [str(x) for x in stage_dirs_override]
        ),
        "score_csv": None if score_csv is None else str(score_csv),
        "selected_stage_by_base": dict(selected_stage_by_base or {}),
        "selected_score_by_base": {
            str(k): float(v) for k, v in dict(selected_score_by_base or {}).items()
        },
        "chunk_tag": None if chunk_tag is None else str(chunk_tag),
        "output_csv": str(out_path),
        "allowed_count": int(len(allowed_bases)),
        "expected_output_rows": int(len(set(allowed_bases) - set(control_bases))),
        "selection_phase": "final",
        "provisional": False,
        "premature": False,
    }


def shard_record_to_task(record: Mapping[str, Any]) -> ScorchTask:
    spec_raw = record.get("spec") if isinstance(record.get("spec"), dict) else {}
    combo_raw = record.get("combo") if isinstance(record.get("combo"), list) else []
    combo = (
        str(combo_raw[0] if len(combo_raw) > 0 else ""),
        str(combo_raw[1] if len(combo_raw) > 1 else ""),
        str(combo_raw[2] if len(combo_raw) > 2 else ""),
    )
    score_csv_raw = record.get("score_csv")
    stage_dirs_raw = record.get("stage_dirs_override")
    return (
        StageSpec(
            source=str(spec_raw.get("source", "")),
            stage_dir=str(spec_raw.get("stage_dir", "")),
            output_name=str(spec_raw.get("output_name", "")),
        ),
        combo,
        Path(str(record.get("receptor", ""))),
        set(str(x) for x in record.get("allowed_bases", []) or []),
        set(str(x) for x in record.get("control_bases", []) or []),
        str(record.get("run_mode", "fda")),
        None if stage_dirs_raw is None else [str(x) for x in stage_dirs_raw],
        None if score_csv_raw is None else Path(str(score_csv_raw)),
        {str(k): str(v) for k, v in dict(record.get("selected_stage_by_base") or {}).items()},
        {
            str(k): float(v)
            for k, v in dict(record.get("selected_score_by_base") or {}).items()
        },
        None if record.get("chunk_tag") is None else str(record.get("chunk_tag")),
    )


def _write_plan_if_absent(path: Path, payload: Mapping[str, Any]) -> Dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload_policy = str(payload.get("selection_policy", "") or "")
    payload_signature = str(payload.get("plan_signature", "") or "")
    if path.exists():
        existing = _read_json(path)
        if existing is not None:
            existing_policy = str(existing.get("selection_policy", "") or "")
            existing_signature = str(existing.get("plan_signature", "") or "")
            if (
                existing_policy == payload_policy
                and existing_signature == payload_signature
            ):
                return existing
    tmp_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=path.name + ".",
            suffix=".tmp",
            dir=str(path.parent),
            delete=False,
        ) as handle:
            tmp_path = Path(handle.name)
            json.dump(payload, handle, sort_keys=True)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(str(tmp_path), str(path))
        except FileExistsError:
            existing = _read_json(path)
            if existing is not None:
                existing_policy = str(existing.get("selection_policy", "") or "")
                existing_signature = str(existing.get("plan_signature", "") or "")
                if (
                    existing_policy == payload_policy
                    and existing_signature == payload_signature
                ):
                    return existing
            tmp_path.replace(path)
            return dict(payload)
        return dict(payload)
    finally:
        try:
            if tmp_path is not None and tmp_path.exists():
                tmp_path.unlink()
        except Exception:
            pass


def ensure_shard_plan(
    *,
    cfg: Mapping[str, Any],
    run_id: str,
    combo: ComboKey,
    tasks: Sequence[ScorchTask],
    post_root: Path,
    decoy_prefix: str,
    top_fraction: float,
) -> Dict[str, Any]:
    records = [
        task_to_shard_record(
            run_id=run_id,
            task=task,
            post_root=post_root,
            decoy_prefix=decoy_prefix,
            top_fraction=top_fraction,
        )
        for task in tasks
    ]
    plan_signature = _stable_hash(
        {
            "selection_policy": SCORCH_SHARD_PLAN_POLICY,
            "run_id": str(run_id),
            "combo": list(combo),
            "decoy_prefix": str(decoy_prefix),
            "top_fraction": float(top_fraction),
            "shards": [
                {
                    "spec": record.get("spec"),
                    "run_mode": record.get("run_mode"),
                    "chunk_tag": record.get("chunk_tag"),
                    "output_csv": record.get("output_csv"),
                    "allowed_bases": record.get("allowed_bases"),
                    "control_bases": record.get("control_bases"),
                }
                for record in records
            ],
        }
    )
    payload = {
        "schema": "atlas.scorch_shard_plan.v1",
        "run_id": str(run_id),
        "combo": list(combo),
        "decoy_prefix": str(decoy_prefix),
        "top_fraction": float(top_fraction),
        "selection_policy": SCORCH_SHARD_PLAN_POLICY,
        "plan_signature": plan_signature,
        "created_at": float(time.time()),
        "shards": records,
    }
    return _write_plan_if_absent(_plan_path(cfg, run_id, combo, decoy_prefix), payload)


def read_all_shard_plans(cfg: Mapping[str, Any], run_id: str) -> List[Dict[str, Any]]:
    plans: List[Dict[str, Any]] = []
    plans_dir = _plans_dir(cfg, run_id)
    if not plans_dir.is_dir():
        return plans
    for path in sorted(plans_dir.glob("*.json")):
        payload = _read_json(path)
        if payload is None:
            continue
        plan = dict(payload)
        plan["plan_path"] = str(path)
        plans.append(plan)
    return plans


def read_all_shard_results(cfg: Mapping[str, Any], run_id: str) -> Dict[str, Dict[str, Any]]:
    results: Dict[str, Dict[str, Any]] = {}
    results_dir = _results_dir(cfg, run_id)
    if not results_dir.is_dir():
        return results
    for path in sorted(results_dir.glob("*.json")):
        payload = _read_json(path)
        if payload is None:
            continue
        shard_id = str(payload.get("shard_id") or path.stem).strip()
        if not shard_id:
            continue
        result = dict(payload)
        result["result_path"] = str(path)
        results[shard_id] = result
    return results


def read_all_shard_records(cfg: Mapping[str, Any], run_id: str) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    seen: Set[str] = set()
    for plan in read_all_shard_plans(cfg, run_id):
        shards = plan.get("shards")
        if not isinstance(shards, list):
            continue
        for idx, raw_record in enumerate(shards):
            if not isinstance(raw_record, Mapping):
                continue
            record = dict(raw_record)
            shard_id = str(record.get("shard_id", "") or "").strip()
            if shard_id and shard_id in seen:
                continue
            if shard_id:
                seen.add(shard_id)
            record["plan_path"] = str(plan.get("plan_path", "") or "")
            record["plan_signature"] = str(plan.get("plan_signature", "") or "")
            record["selection_policy"] = str(plan.get("selection_policy", "") or "")
            record["plan_created_at"] = plan.get("created_at")
            record["plan_record_index"] = int(idx)
            identity_raw = record.get("identity")
            identity = identity_raw if isinstance(identity_raw, Mapping) else {}
            record["decoy_prefix"] = str(
                record.get("decoy_prefix")
                or identity.get("decoy_prefix")
                or plan.get("decoy_prefix")
                or ""
            )
            records.append(record)
    return records


def _record_combo(record: Mapping[str, Any]) -> ComboKey:
    combo_raw = record.get("combo") if isinstance(record.get("combo"), list) else []
    return (
        str(combo_raw[0] if len(combo_raw) > 0 else ""),
        str(combo_raw[1] if len(combo_raw) > 1 else ""),
        str(combo_raw[2] if len(combo_raw) > 2 else ""),
    )


def _record_decoy_prefix(record: Mapping[str, Any]) -> str:
    identity_raw = record.get("identity")
    identity = identity_raw if isinstance(identity_raw, Mapping) else {}
    return str(record.get("decoy_prefix") or identity.get("decoy_prefix") or "")


def _record_scope_key(record: Mapping[str, Any]) -> Tuple[str, str, str, str, str, str, str, str]:
    combo = _record_combo(record)
    spec_raw = record.get("spec") if isinstance(record.get("spec"), Mapping) else {}
    return (
        combo[0],
        combo[1],
        combo[2],
        str(record.get("run_mode", "") or ""),
        str(spec_raw.get("source", "") or ""),
        str(spec_raw.get("stage_dir", "") or ""),
        str(spec_raw.get("output_name", "") or ""),
        _record_decoy_prefix(record),
    )


def _records_for_existing_scope(
    cfg: Mapping[str, Any],
    run_id: str,
    records: Sequence[Mapping[str, Any]],
) -> List[Mapping[str, Any]]:
    scope = {_record_scope_key(record) for record in records if isinstance(record, Mapping)}
    if not scope:
        return list(records)
    return [
        record
        for record in read_all_shard_records(cfg, run_id)
        if _record_scope_key(record) in scope
    ]


def _split_timeout_record(
    *,
    cfg: Optional[Mapping[str, Any]] = None,
    run_id: str,
    record: Mapping[str, Any],
    post_root: Path,
    decoy_prefix: str,
    top_fraction: float,
    attempt: int,
) -> List[Dict[str, Any]]:
    task = shard_record_to_task(record)
    (
        spec,
        combo,
        receptor,
        allowed_bases,
        control_bases,
        run_mode,
        stage_dirs_override,
        score_csv,
        selected_stage_by_base,
        selected_score_by_base,
        chunk_tag,
    ) = task
    ordered = sorted(str(base) for base in allowed_bases)
    if len(ordered) <= 1:
        return []
    fanout_max = _timeout_ligand_fanout_max(cfg or {})
    if fanout_max > 0 and len(ordered) <= fanout_max:
        chunks = [{base} for base in ordered]
    else:
        parts = 4 if len(ordered) >= 256 else 2
        chunk_size = max(1, int(math.ceil(len(ordered) / float(parts))))
        chunks = [
            set(ordered[idx : idx + chunk_size])
            for idx in range(0, len(ordered), chunk_size)
        ]
    chunks = [chunk for chunk in chunks if chunk]
    if len(chunks) <= 1:
        return []
    root_tag = str(chunk_tag or "part000")
    child_records: List[Dict[str, Any]] = []
    for idx, bases in enumerate(chunks):
        child_task: ScorchTask = (
            spec,
            combo,
            receptor,
            set(bases),
            set(control_bases),
            run_mode,
            stage_dirs_override,
            score_csv,
            {
                base: stage
                for base, stage in dict(selected_stage_by_base or {}).items()
                if base in bases
            },
            {
                base: score
                for base, score in dict(selected_score_by_base or {}).items()
                if base in bases
            },
            f"{root_tag}t{max(1, int(attempt)):02d}r{idx:02d}",
        )
        child = task_to_shard_record(
            run_id=run_id,
            task=child_task,
            post_root=post_root,
            decoy_prefix=decoy_prefix,
            top_fraction=top_fraction,
        )
        child["parent_shard_id"] = str(record.get("shard_id", ""))
        child["split_reason"] = "timeout"
        child["split_attempt"] = int(attempt)
        child_records.append(child)
    return child_records


def _append_split_children_to_plan(
    cfg: Mapping[str, Any],
    run_id: str,
    parent_record: Mapping[str, Any],
    child_records: Sequence[Mapping[str, Any]],
    *,
    decoy_prefix: str,
) -> int:
    if not child_records:
        return 0
    combo = _record_combo(parent_record)
    path = _plan_path(cfg, run_id, combo, decoy_prefix)
    lock_path = _plan_lock_path(path)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            plan = _read_json(path)
            if plan is None:
                return 0
            shards = plan.get("shards")
            if not isinstance(shards, list):
                shards = []
                plan["shards"] = shards
            existing = {
                str(record.get("shard_id", ""))
                for record in shards
                if isinstance(record, Mapping)
            }
            appended = 0
            for child in child_records:
                child_id = str(child.get("shard_id", ""))
                if not child_id or child_id in existing:
                    continue
                shards.append(dict(child))
                existing.add(child_id)
                appended += 1
            if appended <= 0:
                return 0
            runtime_splits = plan.get("runtime_splits")
            if not isinstance(runtime_splits, list):
                runtime_splits = []
                plan["runtime_splits"] = runtime_splits
            runtime_splits.append(
                {
                    "parent_shard_id": str(parent_record.get("shard_id", "")),
                    "child_shard_ids": [
                        str(child.get("shard_id", "")) for child in child_records
                    ],
                    "reason": "timeout",
                    "created_at": float(time.time()),
                }
            )
            _write_json_atomic(path, plan)
            return int(appended)
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def _record_top_fraction(record: Mapping[str, Any], default: float = 0.10) -> float:
    identity_raw = record.get("identity")
    identity = identity_raw if isinstance(identity_raw, Mapping) else {}
    for raw in (record.get("top_fraction"), identity.get("top_fraction")):
        if raw is None:
            continue
        try:
            value = float(raw)
        except Exception:
            continue
        if value > 0.0:
            return float(value)
    return float(default)


def _post_root_from_record(record: Mapping[str, Any]) -> Optional[Path]:
    output_raw = str(record.get("output_csv", "") or "").strip()
    if not output_raw:
        return None
    path = Path(output_raw).parent
    for part in _record_combo(record):
        if str(part or "").strip():
            path = path.parent
    return path


def repair_missing_split_child_plan_records(
    cfg: Mapping[str, Any],
    run_id: str,
    *,
    records: Optional[Sequence[Mapping[str, Any]]] = None,
    logger: Optional[logging.Logger] = None,
) -> int:
    active_logger = logger if logger is not None else logging.getLogger(__name__)
    appended_total = 0
    for _pass in range(8):
        current_records = (
            list(records)
            if records is not None and _pass == 0
            else read_all_shard_records(cfg, run_id)
        )
        by_id = {
            str(record.get("shard_id", "")): record
            for record in current_records
            if isinstance(record, Mapping) and str(record.get("shard_id", ""))
        }
        results = read_all_shard_results(cfg, run_id)
        appended_this_pass = 0
        for parent_id, parent_record in list(by_id.items()):
            result = results.get(parent_id)
            if not isinstance(result, Mapping):
                continue
            if str(result.get("status", "") or "").strip().lower() != "split":
                continue
            child_ids_raw = result.get("child_shard_ids")
            if not isinstance(child_ids_raw, list):
                continue
            child_ids = {str(child_id) for child_id in child_ids_raw if str(child_id)}
            missing_child_ids = sorted(child_ids - set(by_id))
            if not missing_child_ids:
                continue
            post_root = _post_root_from_record(parent_record)
            if post_root is None:
                active_logger.warning(
                    "[scorch-shard.repair] run_id=%s parent=%s status=skip reason=missing_post_root missing_children=%d",
                    run_id,
                    parent_id,
                    len(missing_child_ids),
                )
                continue
            try:
                parent_attempt = int(result.get("attempt") or 1)
            except Exception:
                parent_attempt = 1
            generated = _split_timeout_record(
                cfg=cfg,
                run_id=run_id,
                record=parent_record,
                post_root=post_root,
                decoy_prefix=_record_decoy_prefix(parent_record),
                top_fraction=_record_top_fraction(parent_record),
                attempt=max(1, parent_attempt + 1),
            )
            generated_by_id = {
                str(child.get("shard_id", "")): child
                for child in generated
                if str(child.get("shard_id", ""))
            }
            missing_generated = [
                child_id
                for child_id in missing_child_ids
                if child_id not in generated_by_id
            ]
            if missing_generated:
                active_logger.warning(
                    "[scorch-shard.repair] run_id=%s parent=%s status=partial reason=child_id_mismatch missing=%s",
                    run_id,
                    parent_id,
                    ",".join(missing_generated[:5]),
                )
            children_to_append = [
                generated_by_id[child_id]
                for child_id in missing_child_ids
                if child_id in generated_by_id
            ]
            appended = _append_split_children_to_plan(
                cfg,
                run_id,
                parent_record,
                children_to_append,
                decoy_prefix=_record_decoy_prefix(parent_record),
            )
            if appended > 0:
                appended_this_pass += int(appended)
                active_logger.info(
                    "[scorch-shard.repair] run_id=%s parent=%s action=append_missing_split_children count=%d",
                    run_id,
                    parent_id,
                    int(appended),
                )
        appended_total += appended_this_pass
        if appended_this_pass <= 0:
            break
        records = None
    return int(appended_total)


def planned_completed_output_csvs(
    cfg: Mapping[str, Any],
    run_id: str,
    combo: ComboKey,
    specs: Sequence[StageSpec],
    run_mode: str,
    decoy_prefix: str,
) -> Optional[List[Path]]:
    wanted_specs = {
        (str(spec.source), str(spec.stage_dir), str(spec.output_name)) for spec in specs
    }
    repair_missing_split_child_plan_records(
        cfg,
        run_id,
        logger=logging.getLogger(__name__),
    )
    prefix = str(decoy_prefix or "").strip()
    records: List[Mapping[str, Any]] = []
    for record in read_all_shard_records(cfg, run_id):
        combo_raw = record.get("combo") if isinstance(record.get("combo"), list) else []
        record_combo: ComboKey = (
            str(combo_raw[0] if len(combo_raw) > 0 else ""),
            str(combo_raw[1] if len(combo_raw) > 1 else ""),
            str(combo_raw[2] if len(combo_raw) > 2 else ""),
        )
        if record_combo != combo:
            continue
        record_prefix = str(record.get("decoy_prefix") or "").strip()
        if prefix and record_prefix and record_prefix != prefix:
            continue
        if str(record.get("run_mode", "") or "") != str(run_mode):
            continue
        spec_raw = record.get("spec") if isinstance(record.get("spec"), dict) else {}
        spec_key = (
            str(spec_raw.get("source", "") or ""),
            str(spec_raw.get("stage_dir", "") or ""),
            str(spec_raw.get("output_name", "") or ""),
        )
        if spec_key not in wanted_specs:
            continue
        records.append(record)
    if not records:
        return None
    repair_missing_completed_shard_results(
        cfg,
        run_id,
        records=records,
        logger=logging.getLogger(__name__),
    )
    out_paths: List[Path] = []
    seen: Set[Path] = set()
    incomplete = False
    children_by_parent = _children_by_parent(records)
    for record in records:
        shard_id = str(record.get("shard_id", ""))
        if shard_id in children_by_parent:
            continue
        result = _read_result(cfg, run_id, str(record.get("shard_id", ""))) or {}
        status = str(result.get("status", "") or "").strip().lower()
        if status == "split":
            continue
        if not _result_completed_and_valid(cfg, run_id, record):
            incomplete = True
            continue
        if bool(result.get("no_output", False)):
            continue
        output_raw = str(record.get("output_csv", "") or "").strip()
        if not output_raw:
            continue
        output_path = Path(output_raw)
        try:
            resolved = output_path.resolve()
        except Exception:
            resolved = output_path
        if resolved in seen:
            continue
        seen.add(resolved)
        out_paths.append(output_path)
    if incomplete:
        return []
    return out_paths


def planned_completed_output_bases(
    cfg: Mapping[str, Any],
    run_id: str,
    combo: ComboKey,
    specs: Sequence[StageSpec],
    run_mode: str,
    decoy_prefix: str,
) -> Optional[Set[str]]:
    wanted_specs = {
        (str(spec.source), str(spec.stage_dir), str(spec.output_name)) for spec in specs
    }
    prefix = str(decoy_prefix or "").strip()
    records: List[Mapping[str, Any]] = []
    for record in read_all_shard_records(cfg, run_id):
        combo_raw = record.get("combo") if isinstance(record.get("combo"), list) else []
        record_combo: ComboKey = (
            str(combo_raw[0] if len(combo_raw) > 0 else ""),
            str(combo_raw[1] if len(combo_raw) > 1 else ""),
            str(combo_raw[2] if len(combo_raw) > 2 else ""),
        )
        if record_combo != combo:
            continue
        record_prefix = str(record.get("decoy_prefix") or "").strip()
        if prefix and record_prefix and record_prefix != prefix:
            continue
        if str(record.get("run_mode", "") or "") != str(run_mode):
            continue
        spec_raw = record.get("spec") if isinstance(record.get("spec"), dict) else {}
        spec_key = (
            str(spec_raw.get("source", "") or ""),
            str(spec_raw.get("stage_dir", "") or ""),
            str(spec_raw.get("output_name", "") or ""),
        )
        if spec_key not in wanted_specs:
            continue
        records.append(record)
    if not records:
        return None
    children_by_parent = _children_by_parent(records)
    bases: Set[str] = set()
    for record in records:
        shard_id = str(record.get("shard_id", ""))
        if shard_id in children_by_parent:
            continue
        result = _read_result(cfg, run_id, shard_id) or {}
        status = str(result.get("status", "") or "").strip().lower()
        if status == "split":
            continue
        if not _result_completed_and_valid(cfg, run_id, record):
            continue
        bases.update(_expected_output_bases(record))
    return bases


def _expected_output_rows(record: Mapping[str, Any]) -> int:
    explicit = record.get("expected_output_rows")
    if explicit is not None:
        try:
            return max(0, int(explicit or 0))
        except Exception:
            pass
    allowed_raw = record.get("allowed_bases")
    control_raw = record.get("control_bases")
    if isinstance(allowed_raw, list):
        allowed = {str(x) for x in allowed_raw if str(x).strip()}
        controls = (
            {str(x) for x in control_raw if str(x).strip()}
            if isinstance(control_raw, list)
            else set()
        )
        return max(0, len(allowed - controls))
    try:
        return max(0, int(record.get("allowed_count") or 0))
    except Exception:
        return 0


def _expected_output_bases(record: Mapping[str, Any]) -> Set[str]:
    allowed_raw = record.get("allowed_bases")
    if not isinstance(allowed_raw, list):
        return set()
    control_raw = record.get("control_bases")
    allowed = {_ligand_key(x) for x in allowed_raw if _ligand_key(x)}
    controls = (
        {_ligand_key(x) for x in control_raw if _ligand_key(x)}
        if isinstance(control_raw, list)
        else set()
    )
    return allowed - controls


def _allowed_output_bases(record: Mapping[str, Any]) -> Set[str]:
    allowed_raw = record.get("allowed_bases")
    if not isinstance(allowed_raw, list):
        return set()
    return {_ligand_key(x) for x in allowed_raw if _ligand_key(x)}


def _output_valid(record: Mapping[str, Any]) -> Tuple[bool, int]:
    output_raw = str(record.get("output_csv", "") or "").strip()
    if not output_raw:
        return _expected_output_rows(record) <= 0, 0
    output_path = Path(output_raw)
    if not output_path.exists() or output_path.stat().st_size <= 0:
        return False, 0
    row_count = _csv_data_rows(output_path)
    expected_rows = _expected_output_rows(record)
    if expected_rows > 0 and row_count < expected_rows:
        return False, int(row_count)
    expected_bases = _expected_output_bases(record)
    if expected_bases:
        observed_bases = _csv_ligand_bases(output_path)
        if not expected_bases.issubset(observed_bases):
            return False, int(row_count)
        allowed_bases = _allowed_output_bases(record)
        unexpected_bases = observed_bases - (allowed_bases or expected_bases)
        if unexpected_bases:
            return False, int(row_count)
    return row_count > 0, int(row_count)


def _write_quarantine_output(
    record: Mapping[str, Any],
    *,
    reason: str,
    timeout_sec: float,
    returncode: Optional[int],
) -> Tuple[Optional[Path], int]:
    output_raw = str(record.get("output_csv", "") or "").strip()
    if not output_raw:
        return None, 0
    bases = sorted(_expected_output_bases(record))
    if not bases:
        return None, 0
    output_path = Path(output_raw)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combo_raw = record.get("combo") if isinstance(record.get("combo"), list) else []
    spec_raw = record.get("spec") if isinstance(record.get("spec"), Mapping) else {}
    selected_stage = {
        _ligand_key(k): str(v)
        for k, v in dict(record.get("selected_stage_by_base") or {}).items()
        if _ligand_key(k)
    }
    selected_score = {
        _ligand_key(k): v
        for k, v in dict(record.get("selected_score_by_base") or {}).items()
        if _ligand_key(k)
    }
    receptor_name = Path(str(record.get("receptor", "") or "")).name
    fields = [
        "pdb_id",
        "variant",
        "ph",
        "source",
        "stage_dir",
        "run_mode",
        "selected_stage",
        "rescored_stage",
        "selected_docking_score",
        "stage_match_flag",
        "stage_fallback_reason",
        "ligand_file",
        "Receptor",
        "Ligand_ID",
        "Pose_Number",
        "SCORCH_score",
        "SCORCH_certainty",
        "scorch_unscorable_flag",
        "scorch_status",
        "scorch_failure_reason",
        "scorch_timeout_sec",
        "scorch_returncode",
    ]
    rows: List[Dict[str, str]] = []
    for base in bases:
        score_raw = selected_score.get(base)
        try:
            score_text = "" if score_raw is None else f"{float(score_raw):.6g}"
        except Exception:
            score_text = ""
        rows.append(
            {
                "pdb_id": str(combo_raw[0] if len(combo_raw) > 0 else ""),
                "variant": str(combo_raw[1] if len(combo_raw) > 1 else ""),
                "ph": str(combo_raw[2] if len(combo_raw) > 2 else ""),
                "source": str(spec_raw.get("source", "") or ""),
                "stage_dir": str(spec_raw.get("stage_dir", "") or ""),
                "run_mode": str(record.get("run_mode", "") or ""),
                "selected_stage": selected_stage.get(base, ""),
                "rescored_stage": "quarantined",
                "selected_docking_score": score_text,
                "stage_match_flag": "0",
                "stage_fallback_reason": "scorch_quarantined",
                "ligand_file": base,
                "Receptor": receptor_name,
                "Ligand_ID": base,
                "Pose_Number": "1",
                "SCORCH_score": "0",
                "SCORCH_certainty": "0",
                "scorch_unscorable_flag": "1",
                "scorch_status": "quarantined",
                "scorch_failure_reason": str(reason or "unknown"),
                "scorch_timeout_sec": f"{float(timeout_sec):.1f}",
                "scorch_returncode": "" if returncode is None else str(returncode),
            }
        )
    tmp_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            newline="",
            encoding="utf-8",
            prefix=output_path.name + ".",
            suffix=".tmp",
            dir=str(output_path.parent),
            delete=False,
        ) as handle:
            tmp_path = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        tmp_path.replace(output_path)
        return output_path, len(rows)
    finally:
        try:
            if tmp_path is not None and tmp_path.exists():
                tmp_path.unlink()
        except Exception:
            pass


def _single_expected_ligand_base(record: Mapping[str, Any]) -> str:
    bases = sorted(_expected_output_bases(record))
    if len(bases) != 1:
        return ""
    return str(bases[0] or "").strip()


def _global_ligand_timeout_threshold(cfg: Mapping[str, Any]) -> int:
    return max(
        1,
        _int_cfg_or_env(
            cfg,
            "ATLAS_SCORCH_LIGAND_TIMEOUT_QUARANTINE_THRESHOLD",
            "SCORCH_LIGAND_TIMEOUT_QUARANTINE_THRESHOLD",
            default=2,
        ),
    )


def _explicit_ligand_quarantine_bases(cfg: Mapping[str, Any]) -> Set[str]:
    raw = os.environ.get("ATLAS_SCORCH_LIGAND_QUARANTINE_BASES")
    if raw is None:
        raw = os.environ.get("SCORCH_LIGAND_QUARANTINE_BASES")
    if raw is None:
        raw = cfg.get("ATLAS_SCORCH_LIGAND_QUARANTINE_BASES")
    if raw is None:
        raw = cfg.get("SCORCH_LIGAND_QUARANTINE_BASES")
    tokens = re.split(r"[\s,;]+", str(raw or ""))
    return {_ligand_key(token) for token in tokens if _ligand_key(token)}


def _record_ligand_timeout(
    cfg: Mapping[str, Any],
    run_id: str,
    record: Mapping[str, Any],
    *,
    attempt: int,
    timeout_sec: float,
    returncode: Optional[int],
    reason: str,
) -> Dict[str, Any]:
    ligand_base = _single_expected_ligand_base(record)
    if not ligand_base:
        return {}
    run_mode = str(record.get("run_mode", "") or "")
    path = _ligand_timeout_path(
        cfg,
        run_id,
        run_mode=run_mode,
        ligand_base=ligand_base,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = path.with_name(f"{path.name}.lock")
    threshold = _global_ligand_timeout_threshold(cfg)
    combo_raw = record.get("combo") if isinstance(record.get("combo"), list) else []
    combo = [
        str(combo_raw[0] if len(combo_raw) > 0 else ""),
        str(combo_raw[1] if len(combo_raw) > 1 else ""),
        str(combo_raw[2] if len(combo_raw) > 2 else ""),
    ]
    event = {
        "recorded_at": float(time.time()),
        "run_id": str(run_id),
        "ligand_base": ligand_base,
        "run_mode": run_mode,
        "combo": combo,
        "shard_id": str(record.get("shard_id", "") or ""),
        "chunk_tag": str(record.get("chunk_tag", "") or ""),
        "attempt": int(max(1, attempt)),
        "timeout_sec": float(timeout_sec),
        "returncode": None if returncode is None else int(returncode),
        "reason": str(reason or "timeout"),
    }
    with lock_path.open("a+", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            payload = _read_json(path) or {}
            events = payload.get("timeout_events")
            if not isinstance(events, list):
                events = []
            event_key = (
                str(event["combo"]),
                str(event["shard_id"]),
                int(event["attempt"]),
            )
            existing_keys = {
                (
                    str(item.get("combo", "")),
                    str(item.get("shard_id", "")),
                    int(item.get("attempt") or 0),
                )
                for item in events
                if isinstance(item, Mapping)
            }
            if event_key not in existing_keys:
                events.append(event)
            timed_out_combos = {
                tuple(item.get("combo", []) or [])
                for item in events
                if isinstance(item, Mapping)
            }
            timeout_count = len(timed_out_combos)
            payload = {
                "schema": "atlas.scorch_ligand_timeout.v1",
                "run_id": str(run_id),
                "ligand_base": ligand_base,
                "run_mode": run_mode,
                "timeout_count": int(timeout_count),
                "threshold": int(threshold),
                "quarantined": bool(timeout_count >= threshold),
                "last_timeout_sec": float(timeout_sec),
                "last_returncode": None if returncode is None else int(returncode),
                "last_reason": str(reason or "timeout"),
                "updated_at": float(time.time()),
                "timeout_events": events,
            }
            _write_json_atomic(path, payload)
            return dict(payload)
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def _global_ligand_quarantine_info(
    cfg: Mapping[str, Any],
    run_id: str,
    record: Mapping[str, Any],
) -> Dict[str, Any]:
    ligand_base = _single_expected_ligand_base(record)
    if not ligand_base:
        return {}
    if ligand_base in _explicit_ligand_quarantine_bases(cfg):
        return {
            "ligand_base": ligand_base,
            "quarantined": True,
            "reason": "explicit_ligand_quarantine",
            "timeout_count": 0,
            "threshold": _global_ligand_timeout_threshold(cfg),
        }
    run_mode = str(record.get("run_mode", "") or "")
    path = _ligand_timeout_path(
        cfg,
        run_id,
        run_mode=run_mode,
        ligand_base=ligand_base,
    )
    payload = _read_json(path) if path.exists() else None
    if not isinstance(payload, Mapping):
        return {}
    try:
        count = int(payload.get("timeout_count") or 0)
    except Exception:
        count = 0
    threshold = _global_ligand_timeout_threshold(cfg)
    if bool(payload.get("quarantined")) or count >= threshold:
        info = dict(payload)
        info["reason"] = str(info.get("reason") or "global_ligand_timeout")
        info["threshold"] = int(threshold)
        return info
    return {}


def _write_global_ligand_quarantine_result(
    cfg: Mapping[str, Any],
    run_id: str,
    record: Mapping[str, Any],
    *,
    reason: str,
    timeout_sec: float,
    returncode: Optional[int],
    logger: logging.Logger,
) -> bool:
    if not _single_expected_ligand_base(record):
        return False
    _remove_materialized_inputs(record)
    quarantine_path, quarantine_rows = _write_quarantine_output(
        record,
        reason=reason,
        timeout_sec=timeout_sec,
        returncode=returncode,
    )
    if quarantine_path is None or quarantine_rows <= 0:
        return False
    valid, _rows = _output_valid(record)
    if not valid:
        return False
    write_shard_result(
        cfg,
        run_id,
        record,
        status="quarantined",
        payload={
            "attempt": 0,
            "output_csv": str(quarantine_path),
            "row_count": int(quarantine_rows),
            "expected_row_count": int(_expected_output_rows(record)),
            "row_coverage_ok": True,
            "quarantine_output": True,
            "quarantined_ligand_count": int(quarantine_rows),
            "failure_reason": reason,
            "timed_out": True,
            "timeout_sec": float(timeout_sec),
            "returncode": returncode,
            "global_ligand_quarantine": True,
        },
    )
    logger.warning(
        "[scorch-shard.result] run_id=%s shard_id=%s status=quarantined reason=%s ligand=%s output=%s",
        run_id,
        str(record.get("shard_id", "") or ""),
        reason,
        _single_expected_ligand_base(record),
        str(quarantine_path),
    )
    return True


def _remove_invalid_output(record: Mapping[str, Any]) -> bool:
    valid, _rows = _output_valid(record)
    if valid:
        return False
    output_raw = str(record.get("output_csv", "") or "").strip()
    if not output_raw:
        return False
    output_path = Path(output_raw)
    if not output_path.exists():
        return False
    try:
        output_path.unlink()
        return True
    except Exception:
        return False


def _remove_materialized_inputs(record: Mapping[str, Any]) -> bool:
    output_raw = str(record.get("output_csv", "") or "").strip()
    if not output_raw:
        return False
    spec_raw = record.get("spec") if isinstance(record.get("spec"), dict) else {}
    stage_dir = str(spec_raw.get("stage_dir", "") or "").strip()
    run_mode = str(record.get("run_mode", "") or "").strip()
    chunk_tag = str(record.get("chunk_tag", "") or "").strip()
    if not stage_dir or not run_mode:
        return False
    input_dir = Path(output_raw).parent / ".scorch_inputs" / run_mode / stage_dir
    if chunk_tag:
        input_dir = input_dir / chunk_tag
    if not input_dir.exists():
        return False
    try:
        shutil.rmtree(input_dir)
        return True
    except Exception:
        return False


def _read_result(cfg: Mapping[str, Any], run_id: str, shard_id: str) -> Optional[Dict[str, Any]]:
    return _read_json(_result_path(cfg, run_id, shard_id))


def _result_completed_and_valid(
    cfg: Mapping[str, Any],
    run_id: str,
    record: Mapping[str, Any],
) -> bool:
    shard_id = str(record.get("shard_id", ""))
    result = _read_result(cfg, run_id, shard_id)
    status = str((result or {}).get("status", "")).lower()
    if not result or status not in {"completed", "quarantined"}:
        return False
    if bool(result.get("no_output", False)):
        return _expected_output_rows(record) <= 0
    if status == "quarantined" and not bool(result.get("quarantine_output", False)):
        return False
    valid, _rows = _output_valid(record)
    return bool(valid)


def _repair_completed_output(
    cfg: Mapping[str, Any],
    run_id: str,
    record: Mapping[str, Any],
    *,
    logger: logging.Logger,
) -> bool:
    shard_id = str(record.get("shard_id", ""))
    if _read_result(cfg, run_id, shard_id) is not None:
        return False
    valid, row_count = _output_valid(record)
    if not valid:
        return False
    write_shard_result(
        cfg,
        run_id,
        record,
        status="completed",
        payload={
            "output_csv": str(record.get("output_csv", "")),
            "row_count": int(row_count),
            "expected_row_count": int(_expected_output_rows(record)),
            "repaired_from_output": True,
            "repair_reason": "missing_completed_result_json",
            "result_kind": "final",
            "final": True,
            "provisional": False,
            "premature": False,
        },
    )
    logger.info(
        "[scorch-shard.repair] run_id=%s shard_id=%s status=completed output=%s rows=%d",
        run_id,
        shard_id,
        str(record.get("output_csv", "")),
        int(row_count),
    )
    return True


def repair_missing_completed_shard_results(
    cfg: Mapping[str, Any],
    run_id: str,
    *,
    records: Optional[Sequence[Mapping[str, Any]]] = None,
    logger: Optional[logging.Logger] = None,
) -> int:
    active_logger = logger if logger is not None else logging.getLogger(__name__)
    repaired = 0
    for record in list(records) if records is not None else read_all_shard_records(cfg, run_id):
        if not isinstance(record, Mapping):
            continue
        if _repair_completed_output(cfg, run_id, record, logger=active_logger):
            repaired += 1
    return int(repaired)


def try_claim_shard(
    cfg: Mapping[str, Any],
    run_id: str,
    record: Mapping[str, Any],
    *,
    lease_sec: float,
    max_attempts: int,
) -> Tuple[bool, int]:
    shard_id = str(record.get("shard_id", ""))
    result_path = _result_path(cfg, run_id, shard_id)
    result = _read_json(result_path) if result_path.exists() else None
    if result is not None:
        status = str(result.get("status", "")).strip().lower()
        try:
            prev_attempt = int(result.get("attempt") or 1)
        except Exception:
            prev_attempt = 1
        if status == "completed":
            if _result_completed_and_valid(cfg, run_id, record):
                return False, prev_attempt
            try:
                result_path.unlink()
            except Exception:
                return False, prev_attempt
        elif status == "split":
            return False, prev_attempt
        elif status == "quarantined":
            if not _retry_quarantined_on_resume(cfg):
                return False, prev_attempt
            _remove_materialized_inputs(record)
            try:
                output_raw = str(record.get("output_csv", "") or "").strip()
                if output_raw:
                    Path(output_raw).unlink(missing_ok=True)
            except Exception:
                pass
            try:
                result_path.unlink()
            except Exception:
                return False, prev_attempt
        elif status == "failed" and prev_attempt < max(1, int(max_attempts)):
            try:
                result_path.unlink()
            except Exception:
                return False, prev_attempt
        else:
            return False, prev_attempt
    claim_path = _claim_path(cfg, run_id, shard_id)
    claim_path.parent.mkdir(parents=True, exist_ok=True)
    now = time.time()
    attempt = 1
    if result is not None:
        try:
            attempt = int(result.get("attempt") or 0) + 1
        except Exception:
            attempt = 1
    payload = {
        "run_id": str(run_id),
        "shard_id": shard_id,
        "pid": int(os.getpid()),
        "hostname": os.uname().nodename if hasattr(os, "uname") else "",
        "slurm_job_id": str(os.environ.get("SLURM_JOB_ID") or ""),
        "slurm_array_job_id": str(os.environ.get("SLURM_ARRAY_JOB_ID") or ""),
        "slurm_array_task_id": str(os.environ.get("SLURM_ARRAY_TASK_ID") or ""),
        "slurm_array_task_count": str(os.environ.get("SLURM_ARRAY_TASK_COUNT") or ""),
        "atlas_dist_task_id": str(os.environ.get("ATLAS_DIST_TASK_ID") or ""),
        "atlas_dist_task_count": str(os.environ.get("ATLAS_DIST_TASK_COUNT") or ""),
        "attempt": int(max(1, attempt)),
        "claimed_at": float(now),
        "updated_at": float(now),
        "lease_sec": float(max(1.0, lease_sec)),
    }

    def _attempt_create() -> bool:
        try:
            fd = os.open(str(claim_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
        except FileExistsError:
            return False
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, sort_keys=True)
                handle.flush()
                os.fsync(handle.fileno())
            return True
        except Exception:
            try:
                claim_path.unlink()
            except Exception:
                pass
            return False

    if _attempt_create():
        return True, int(payload["attempt"])
    try:
        if claim_path.exists():
            claim = _read_json(claim_path) or {}
            active, _age, _reason = _claim_liveness(
                cfg,
                claim_path,
                claim,
                now=now,
                lease_sec=float(max(1.0, lease_sec)),
            )
            if not active:
                try:
                    claim_path.unlink()
                except Exception:
                    pass
                if _attempt_create():
                    return True, int(payload["attempt"])
    except Exception:
        pass
    return False, int(payload["attempt"])


def renew_shard_claim(
    cfg: Mapping[str, Any],
    run_id: str,
    shard_id: str,
    *,
    lease_sec: float,
) -> bool:
    claim_path = _claim_path(cfg, run_id, shard_id)
    if not claim_path.exists():
        return False
    now = time.time()
    try:
        payload = _read_json(claim_path) or {}
        payload.update(
            {
                "updated_at": float(now),
                "lease_sec": float(max(1.0, lease_sec)),
                "pid": int(os.getpid()),
            }
        )
        _write_json_atomic(claim_path, payload)
        os.utime(claim_path, (now, now))
        return True
    except Exception:
        return False


def write_shard_result(
    cfg: Mapping[str, Any],
    run_id: str,
    record: Mapping[str, Any],
    *,
    status: str,
    payload: Optional[Mapping[str, Any]] = None,
) -> Path:
    shard_id = str(record.get("shard_id", ""))
    claim_payload = _read_json(_claim_path(cfg, run_id, shard_id)) or {}
    status_value = str(status)
    status_lower = status_value.strip().lower()
    terminal_status = status_lower in {"completed", "split", "quarantined"}
    result_payload: Dict[str, Any] = {
        "schema": "atlas.scorch_shard_result.v1",
        "run_id": str(run_id),
        "shard_id": shard_id,
        "status": status_value,
        "completed_at": float(time.time()),
        "pid": int(os.getpid()),
        "hostname": os.uname().nodename if hasattr(os, "uname") else "",
        "slurm_job_id": str(os.environ.get("SLURM_JOB_ID") or ""),
        "slurm_array_job_id": str(os.environ.get("SLURM_ARRAY_JOB_ID") or ""),
        "slurm_array_task_id": str(os.environ.get("SLURM_ARRAY_TASK_ID") or ""),
        "slurm_array_task_count": str(os.environ.get("SLURM_ARRAY_TASK_COUNT") or ""),
        "atlas_dist_task_id": str(os.environ.get("ATLAS_DIST_TASK_ID") or ""),
        "atlas_dist_task_count": str(os.environ.get("ATLAS_DIST_TASK_COUNT") or ""),
        "claim_pid": claim_payload.get("pid"),
        "claim_hostname": claim_payload.get("hostname"),
        "claim_slurm_job_id": claim_payload.get("slurm_job_id"),
        "claim_slurm_array_task_id": claim_payload.get("slurm_array_task_id"),
        "claim_atlas_dist_task_id": claim_payload.get("atlas_dist_task_id"),
        "combo": list(record.get("combo", []) or []),
        "source": str((record.get("spec") or {}).get("source", "")) if isinstance(record.get("spec"), dict) else "",
        "run_mode": str(record.get("run_mode", "")),
        "chunk_tag": record.get("chunk_tag"),
        "result_kind": "final" if terminal_status else "provisional",
        "terminal": bool(terminal_status),
        "final": bool(status_lower == "completed"),
        "provisional": bool(not terminal_status),
        "premature": False,
    }
    if payload:
        result_payload.update({str(k): v for k, v in payload.items()})
    path = _write_json_atomic(_result_path(cfg, run_id, shard_id), result_payload)
    try:
        claim_path = _claim_path(cfg, run_id, shard_id)
        if claim_path.exists():
            claim_path.unlink()
    except Exception:
        pass
    return path


def _max_attempts_for_cfg(cfg: Mapping[str, Any]) -> int:
    return max(
        1,
        _int_cfg_or_env(
            cfg,
            "ATLAS_SCORCH_SHARD_MAX_ATTEMPTS",
            "SCORCH_SHARD_MAX_ATTEMPTS",
            default=2,
        ),
    )


def _retry_quarantined_on_resume(cfg: Mapping[str, Any]) -> bool:
    for name in (
        "ATLAS_SCORCH_RETRY_QUARANTINED_ON_RESUME",
        "SCORCH_RETRY_QUARANTINED_ON_RESUME",
    ):
        raw = os.environ.get(name)
        if raw is None:
            raw = cfg.get(name)
        if raw is None:
            continue
        return _truthy(raw)
    return False


def _current_hostname() -> str:
    try:
        return os.uname().nodename if hasattr(os, "uname") else ""
    except Exception:
        return ""


def _pid_is_alive(pid: Any) -> bool:
    try:
        parsed = int(pid)
    except Exception:
        return False
    if parsed <= 0:
        return False
    try:
        os.kill(parsed, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False


def _claim_orphan_timeout_sec(cfg: Mapping[str, Any], lease_sec: float) -> float:
    configured = _float_cfg_or_env(
        cfg,
        "ATLAS_SCORCH_SHARD_ORPHAN_CLAIM_SEC",
        "SCORCH_SHARD_ORPHAN_CLAIM_SEC",
        default=300.0,
    )
    return max(30.0, min(max(30.0, float(lease_sec)), float(configured)))


def _claim_liveness(
    cfg: Mapping[str, Any],
    claim_path: Path,
    claim: Mapping[str, Any],
    *,
    now: float,
    lease_sec: float,
) -> Tuple[bool, float, str]:
    try:
        age = max(0.0, now - float(claim_path.stat().st_mtime))
    except Exception:
        age = 0.0
    if age > max(30.0, float(lease_sec)):
        return False, float(age), "lease_expired"

    claim_job = str(claim.get("slurm_job_id") or "").strip()
    current_job = str(os.environ.get("SLURM_JOB_ID") or "").strip()
    if claim_job and claim_job != current_job:
        orphan_sec = _claim_orphan_timeout_sec(cfg, lease_sec)
        if age > orphan_sec:
            return False, float(age), "foreign_slurm_job_orphan"

    claim_host = str(claim.get("hostname") or "").strip()
    if claim_host and claim_host == _current_hostname():
        claim_pid = claim.get("pid")
        if claim_pid and not _pid_is_alive(claim_pid):
            return False, float(age), "dead_local_pid"

    return True, float(age), "active"


def _active_claim_info(
    cfg: Mapping[str, Any],
    run_id: str,
    shard_id: str,
    *,
    now: float,
) -> Tuple[bool, Optional[Dict[str, Any]]]:
    claim_path = _claim_path(cfg, run_id, shard_id)
    if not claim_path.exists():
        return False, None
    claim = _read_json(claim_path) or {}
    try:
        lease_sec = float(claim.get("lease_sec") or 3600.0)
    except Exception:
        lease_sec = 3600.0
    active, age, reason = _claim_liveness(
        cfg,
        claim_path,
        claim,
        now=now,
        lease_sec=lease_sec,
    )
    claim["claim_path"] = str(claim_path)
    claim["claim_age_sec"] = float(age)
    claim["claim_liveness"] = reason
    return bool(active), claim


def _shard_ledger_entry(
    cfg: Mapping[str, Any],
    run_id: str,
    record: Mapping[str, Any],
    result: Optional[Mapping[str, Any]],
    *,
    max_attempts: int,
    now: float,
) -> Dict[str, Any]:
    shard_id = str(record.get("shard_id", "") or "")
    status = str((result or {}).get("status", "") or "").strip().lower()
    valid, observed_rows = _output_valid(record)
    expected_rows = _expected_output_rows(record)
    no_output = bool((result or {}).get("no_output", False))
    completed_valid = bool(
        status == "completed"
        and (
            (no_output and expected_rows <= 0)
            or (valid and (observed_rows > 0 or expected_rows <= 0))
        )
    )
    split_final = bool(status == "split")
    quarantined = bool(status == "quarantined")
    try:
        attempt = int((result or {}).get("attempt") or 0)
    except Exception:
        attempt = 0
    failed_final = bool(status == "failed" and attempt >= max(1, int(max_attempts)))
    active_claim, claim = _active_claim_info(cfg, run_id, shard_id, now=now)
    premature = bool(status == "completed" and not completed_valid)
    terminal = bool(completed_valid or split_final or quarantined or failed_final)
    if completed_valid or split_final:
        state = "completed"
    elif quarantined:
        state = "quarantined"
    elif failed_final:
        state = "failed"
    elif premature:
        state = "premature"
    elif active_claim:
        state = "inflight"
    elif result is not None:
        state = "provisional"
    else:
        state = "pending"
    combo_raw = record.get("combo") if isinstance(record.get("combo"), list) else []
    spec_raw = record.get("spec") if isinstance(record.get("spec"), Mapping) else {}
    output_raw = str(record.get("output_csv", "") or "")
    return {
        "shard_id": shard_id,
        "combo": [
            str(combo_raw[0] if len(combo_raw) > 0 else ""),
            str(combo_raw[1] if len(combo_raw) > 1 else ""),
            str(combo_raw[2] if len(combo_raw) > 2 else ""),
        ],
        "source": str(spec_raw.get("source", "") or ""),
        "stage_dir": str(spec_raw.get("stage_dir", "") or ""),
        "output_name": str(spec_raw.get("output_name", "") or ""),
        "run_mode": str(record.get("run_mode", "") or ""),
        "selection_phase": str(record.get("selection_phase") or "final"),
        "chunk_tag": record.get("chunk_tag"),
        "parent_shard_id": str(record.get("parent_shard_id", "") or ""),
        "decoy_prefix": str(record.get("decoy_prefix", "") or ""),
        "plan_path": str(record.get("plan_path", "") or ""),
        "result_path": str((result or {}).get("result_path", "") or ""),
        "output_csv": output_raw,
        "status": status or "missing",
        "ledger_state": state,
        "result_kind": "final" if terminal else "provisional",
        "terminal": bool(terminal),
        "final": bool(completed_valid or split_final),
        "accounted": bool(completed_valid or split_final or quarantined),
        "provisional": bool(not terminal and not premature),
        "premature": bool(premature),
        "inflight": bool(active_claim),
        "output_valid": bool(valid),
        "row_count": int(observed_rows),
        "expected_row_count": int(expected_rows),
        "attempt": int(attempt),
        "allowed_count": int(record.get("allowed_count") or 0),
        "claim_path": "" if claim is None else str(claim.get("claim_path", "")),
        "claim_age_sec": None if claim is None else float(claim.get("claim_age_sec") or 0.0),
        "claim_liveness": "" if claim is None else str(claim.get("claim_liveness", "") or ""),
    }


def summarize_shard_ledger(
    cfg: Mapping[str, Any],
    run_id: str,
    *,
    records: Optional[Sequence[Mapping[str, Any]]] = None,
    max_attempts: Optional[int] = None,
    workers: int = 0,
    threads: int = 1,
    write: bool = False,
) -> Dict[str, Any]:
    if records is None:
        repair_missing_split_child_plan_records(
            cfg,
            run_id,
            logger=logging.getLogger(__name__),
        )
        records = read_all_shard_records(cfg, run_id)
    else:
        records = list(records)
    results = read_all_shard_results(cfg, run_id)
    attempts = _max_attempts_for_cfg(cfg) if max_attempts is None else max(1, int(max_attempts))
    now = time.time()
    entries = [
        _shard_ledger_entry(
            cfg,
            run_id,
            record,
            results.get(str(record.get("shard_id", "") or "")),
            max_attempts=attempts,
            now=now,
        )
        for record in records
    ]
    children = _children_by_parent(records)
    for entry in entries:
        shard_id = str(entry.get("shard_id", "") or "")
        child_count = len(children.get(shard_id, []))
        if child_count <= 0:
            entry["superseded_by_children"] = False
            entry["child_shard_count"] = 0
            continue
        entry["superseded_by_children"] = True
        entry["child_shard_count"] = int(child_count)
        if str(entry.get("status", "") or "") != "split":
            entry["ledger_state"] = "superseded"
            entry["result_kind"] = "final"
            entry["terminal"] = True
            entry["final"] = True
            entry["accounted"] = True
            entry["provisional"] = False
            entry["premature"] = False
            entry["inflight"] = False
    scored_completed = sum(1 for entry in entries if entry["ledger_state"] == "completed")
    quarantined = sum(1 for entry in entries if entry["ledger_state"] == "quarantined")
    superseded = sum(1 for entry in entries if entry["ledger_state"] == "superseded")
    completed = int(scored_completed + quarantined + superseded)
    failed = sum(1 for entry in entries if entry["ledger_state"] == "failed")
    inflight = sum(1 for entry in entries if entry["ledger_state"] == "inflight")
    pending = sum(1 for entry in entries if entry["ledger_state"] in {"pending", "premature", "provisional"})
    productive = sum(
        1
        for entry in entries
        if entry["ledger_state"] == "completed"
        and entry["status"] == "completed"
        and int(entry["row_count"]) > 0
    )
    zero_output = sum(
        1
        for entry in entries
        if entry["ledger_state"] == "completed"
        and entry["status"] == "completed"
        and int(entry["row_count"]) <= 0
    )
    row_total = sum(
        max(0, int(entry["row_count"]))
        for entry in entries
        if entry["ledger_state"] == "completed"
        and entry["status"] == "completed"
    )
    quarantined_rows = sum(
        max(0, int(entry["row_count"]))
        for entry in entries
        if entry["ledger_state"] == "quarantined" and entry["status"] == "quarantined"
    )
    elapsed_total = 0.0
    thread_seconds = 0.0
    for entry in entries:
        result = results.get(str(entry.get("shard_id", "") or "")) or {}
        if entry["ledger_state"] != "completed" or entry["status"] != "completed":
            continue
        try:
            elapsed_sec = float(result.get("elapsed_sec") or 0.0)
        except Exception:
            elapsed_sec = 0.0
        try:
            thread_count = int(result.get("threads") or threads or 1)
        except Exception:
            thread_count = max(1, int(threads))
        elapsed_total += max(0.0, elapsed_sec)
        thread_seconds += max(0.0, elapsed_sec) * max(1, int(thread_count))
    summary_payload: Dict[str, Any] = {
        "schema": "atlas.scorch_shard_summary.v1",
        "run_id": str(run_id),
        "total": int(len(entries)),
        "completed": int(completed),
        "scored_completed": int(scored_completed),
        "accounted": int(completed),
        "failed": int(failed),
        "inflight": int(inflight),
        "pending": int(pending),
        "all_complete": bool(len(entries) > 0 and completed >= len(entries)),
        "all_scored": bool(len(entries) > 0 and scored_completed >= len(entries)),
        "workers": int(workers),
        "threads": max(1, int(threads)),
        "productive": int(productive),
        "zero_output": int(zero_output),
        "split": int(sum(1 for entry in entries if entry["status"] == "split")),
        "superseded": int(superseded),
        "quarantined": int(quarantined),
        "quarantined_rows": int(quarantined_rows),
        "premature": int(sum(1 for entry in entries if entry["premature"])),
        "provisional": int(sum(1 for entry in entries if entry["provisional"])),
        "row_total": int(row_total),
        "elapsed_sec_total": float(elapsed_total),
        "thread_seconds_total": float(thread_seconds),
        "written_at": float(now),
        "shards": entries,
    }
    if write:
        _write_json_atomic(_summary_path(cfg, run_id), summary_payload)
    return summary_payload


def _summarize_records(
    cfg: Mapping[str, Any],
    run_id: str,
    records: Sequence[Mapping[str, Any]],
    *,
    max_attempts: int,
) -> Tuple[int, int, int, int]:
    completed = 0
    failed = 0
    inflight = 0
    pending = 0
    now = time.time()
    children = _children_by_parent(records)
    for record in records:
        shard_id = str(record.get("shard_id", ""))
        if children.get(shard_id):
            completed += 1
            continue
        result = _read_result(cfg, run_id, shard_id)
        if result is not None:
            status = str(result.get("status", "")).strip().lower()
            if status == "completed" and _result_completed_and_valid(cfg, run_id, record):
                completed += 1
                continue
            if status == "split":
                completed += 1
                continue
            if status == "quarantined":
                completed += 1
                continue
            try:
                attempt = int(result.get("attempt") or 1)
            except Exception:
                attempt = 1
            if status == "failed" and attempt >= max(1, int(max_attempts)):
                failed += 1
                continue
        claim_path = _claim_path(cfg, run_id, shard_id)
        if claim_path.exists():
            try:
                claim = _read_json(claim_path) or {}
                lease_sec = float(claim.get("lease_sec") or 3600.0)
                active, _age, _reason = _claim_liveness(
                    cfg,
                    claim_path,
                    claim,
                    now=now,
                    lease_sec=lease_sec,
                )
            except Exception:
                active = True
            if active:
                inflight += 1
                continue
        pending += 1
    return completed, failed, inflight, pending


def _claim_next_record(
    cfg: Mapping[str, Any],
    run_id: str,
    records: Sequence[Mapping[str, Any]],
    *,
    lease_sec: float,
    max_attempts: int,
    overwrite: bool,
    logger: logging.Logger,
    skip_shard_ids: Optional[Set[str]] = None,
) -> Optional[Tuple[Mapping[str, Any], int]]:
    skip_ids = set(skip_shard_ids or set())
    ordered = sorted(
        records,
        key=lambda rec: int(rec.get("allowed_count", 0) or 0),
        reverse=True,
    )
    children = _children_by_parent(records)
    for record in ordered:
        shard_id = str(record.get("shard_id", "") or "")
        if shard_id in skip_ids:
            continue
        if children.get(shard_id):
            continue
        if not overwrite and _result_completed_and_valid(cfg, run_id, record):
            continue
        if not overwrite:
            _repair_completed_output(cfg, run_id, record, logger=logger)
        if not overwrite and _result_completed_and_valid(cfg, run_id, record):
            continue
        if not overwrite:
            quarantine_info = _global_ligand_quarantine_info(cfg, run_id, record)
            if quarantine_info:
                reason = str(quarantine_info.get("reason") or "global_ligand_timeout")
                timeout_sec = 0.0
                try:
                    timeout_sec = float(quarantine_info.get("last_timeout_sec") or 0.0)
                except Exception:
                    timeout_sec = 0.0
                if _write_global_ligand_quarantine_result(
                    cfg,
                    run_id,
                    record,
                    reason=reason,
                    timeout_sec=timeout_sec,
                    returncode=None,
                    logger=logger,
                ):
                    continue
        claimed, attempt = try_claim_shard(
            cfg,
            run_id,
            record,
            lease_sec=lease_sec,
            max_attempts=max_attempts,
        )
        if claimed:
            return record, attempt
    return None


def execute_sharded_tasks(
    *,
    cfg: Dict[str, object],
    tasks: Sequence[ScorchTask],
    run_root: Path,
    post_root: Path,
    jobs: int,
    threads: int,
    overwrite: bool,
    component: str,
    logger: logging.Logger,
    deps: OrchestrationDeps,
    combo_failed_local: Dict[ComboKey, bool],
    decoy_prefix: str,
    top_fraction: float,
) -> ShardedExecutionResult:
    run_id = run_root.name
    if not tasks:
        return ShardedExecutionResult(0, 0, 0, True, combo_failed_local)
    combos = sorted({task[1] for task in tasks})
    records: List[Mapping[str, Any]] = []
    for combo in combos:
        combo_tasks = [task for task in tasks if task[1] == combo]
        plan = ensure_shard_plan(
            cfg=cfg,
            run_id=run_id,
            combo=combo,
            tasks=combo_tasks,
            post_root=post_root,
            decoy_prefix=decoy_prefix,
            top_fraction=top_fraction,
        )
        records.extend(list(plan.get("shards", []) or []))
    repaired_split_children = repair_missing_split_child_plan_records(
        cfg,
        run_id,
        records=records,
        logger=logger,
    )
    if repaired_split_children > 0:
        records = _records_for_existing_scope(cfg, run_id, records)
        logger.info(
            "[scorch-shard.repair] run_id=%s action=refresh_records reason=missing_split_children count=%d scoped_total=%d",
            run_id,
            int(repaired_split_children),
            len(records),
        )
    if overwrite:
        for record in records:
            try:
                result_path = _result_path(cfg, run_id, str(record.get("shard_id", "")))
                if result_path.exists():
                    result = _read_json(result_path) or {}
                    status = str(result.get("status", "") or "").strip().lower()
                    if status == "split":
                        continue
                    if status == "quarantined" and not _retry_quarantined_on_resume(cfg):
                        continue
                    result_path.unlink()
            except Exception:
                pass

    max_attempts = _max_attempts_for_cfg(cfg)
    lease_sec = _float_cfg_or_env(
        cfg,
        "ATLAS_SCORCH_SHARD_LEASE_SEC",
        "SCORCH_SHARD_LEASE_SEC",
        default=7200.0,
    )
    heartbeat_sec = _float_cfg_or_env(
        cfg,
        "ATLAS_SCORCH_SHARD_HEARTBEAT_SEC",
        "SCORCH_SHARD_HEARTBEAT_SEC",
        default=30.0,
    )
    drain_timeout_sec = _float_cfg_or_env(
        cfg,
        "ATLAS_SCORCH_SHARD_DRAIN_TIMEOUT_SEC",
        "SCORCH_SHARD_DRAIN_TIMEOUT_SEC",
        default=1800.0,
    )
    records_lock = threading.Lock()
    active_condition = threading.Condition()
    active_local = 0
    active_shard_ids: Set[str] = set()

    def _records_snapshot() -> List[Mapping[str, Any]]:
        with records_lock:
            return list(records)

    def _run_one(record: Mapping[str, Any], attempt: int) -> Tuple[bool, ComboKey]:
        task = shard_record_to_task(record)
        shard_id = str(record.get("shard_id", ""))
        combo = task[1]
        stop_event = threading.Event()

        def _heartbeat_loop() -> None:
            interval = min(max(1.0, float(heartbeat_sec)), max(1.0, float(lease_sec) / 3.0))
            while not stop_event.wait(interval):
                if not renew_shard_claim(cfg, run_id, shard_id, lease_sec=lease_sec):
                    logger.warning(
                        "[scorch-shard.heartbeat] run_id=%s shard_id=%s status=lost_claim",
                        run_id,
                        shard_id,
                    )
                    return

        hb = threading.Thread(
            target=_heartbeat_loop,
            name=f"scorch-shard-hb-{shard_id[-8:]}",
            daemon=True,
        )
        hb.start()
        started = time.time()
        try:
            (
                spec,
                combo,
                receptor,
                allowed_bases,
                control_bases,
                run_mode,
                stage_dirs_override,
                score_csv,
                selected_stage_by_base,
                selected_score_by_base,
                chunk_tag,
            ) = task
            removed_partial = _remove_invalid_output(record)
            if removed_partial:
                removed_inputs = _remove_materialized_inputs(record)
                logger.warning(
                    "[scorch-shard.output] run_id=%s shard_id=%s action=delete reason=partial_or_invalid_output output=%s inputs_deleted=%s",
                    run_id,
                    shard_id,
                    str(record.get("output_csv", "")),
                    bool(removed_inputs),
                )
            effective_overwrite = bool(overwrite or removed_partial)
            timeout_sec = _timeout_for_attempt(cfg, record, attempt=int(attempt))
            overrides_raw = cfg.get(_TIMEOUT_OVERRIDES_KEY)
            if not isinstance(overrides_raw, dict):
                overrides_raw = {}
                cfg[_TIMEOUT_OVERRIDES_KEY] = overrides_raw
            overrides_raw[shard_id] = float(timeout_sec)
            ok, out_path = deps.score_stage(
                cfg,
                spec,
                combo,
                run_root,
                post_root,
                receptor,
                max(1, int(threads)),
                effective_overwrite,
                logger,
                allowed_bases,
                control_bases,
                run_mode,
                stage_dirs_override,
                score_csv,
                selected_stage_by_base,
                selected_score_by_base,
                chunk_tag,
                shard_id,
                "primary",
                0,
            )
            stage_result = _stage_result_for(cfg, shard_id)
            row_count = _csv_data_rows(out_path) if out_path is not None else 0
            expected_rows = _expected_output_rows(record)
            output_valid = _output_valid(record)[0] if out_path is not None else False
            row_coverage_ok = (
                bool(out_path is not None)
                and int(row_count) > 0
                and bool(output_valid)
            )
            no_output_ok = bool(out_path is None and expected_rows <= 0)
            shard_ok = bool(ok and (row_coverage_ok or no_output_ok))
            timed_out = bool(stage_result.get("timed_out", False))
            failure_reason = str(stage_result.get("failure_reason") or "").strip()
            if not shard_ok and timed_out:
                failure_reason = "timeout"
            try:
                returncode = (
                    None
                    if stage_result.get("returncode") is None
                    else int(stage_result.get("returncode"))
                )
            except Exception:
                returncode = None
            timeout_sec_observed = float(stage_result.get("timeout_sec") or timeout_sec)
            ligand_timeout_info: Dict[str, Any] = {}
            if not shard_ok and timed_out:
                ligand_timeout_info = _record_ligand_timeout(
                    cfg,
                    run_id,
                    record,
                    attempt=int(attempt),
                    timeout_sec=timeout_sec_observed,
                    returncode=returncode,
                    reason=failure_reason or "timeout",
                )
            global_ligand_quarantine_now = bool(
                ligand_timeout_info.get("quarantined")
                and _single_expected_ligand_base(record)
            )
            split_children: List[Dict[str, Any]] = []
            split_candidate_count = _timeout_split_candidate_count(record)
            split_min = _timeout_split_min(cfg)
            if (
                not shard_ok
                and timed_out
                and not global_ligand_quarantine_now
                and int(attempt) < max_attempts
                and split_candidate_count >= split_min
            ):
                split_children = _split_timeout_record(
                    cfg=cfg,
                    run_id=run_id,
                    record=record,
                    post_root=post_root,
                    decoy_prefix=decoy_prefix,
                    top_fraction=top_fraction,
                    attempt=int(attempt) + 1,
                )
                appended = _append_split_children_to_plan(
                    cfg,
                    run_id,
                    record,
                    split_children,
                    decoy_prefix=decoy_prefix,
                )
                if appended > 0:
                    with records_lock:
                        existing_ids = {str(item.get("shard_id", "")) for item in records}
                        for child in split_children:
                            if str(child.get("shard_id", "")) not in existing_ids:
                                records.append(child)
                                existing_ids.add(str(child.get("shard_id", "")))
                    with active_condition:
                        active_condition.notify_all()
                    write_shard_result(
                        cfg,
                        run_id,
                        record,
                        status="split",
                        payload={
                            "attempt": int(attempt),
                            "allowed_count": int(len(allowed_bases)),
                            "expected_row_count": int(expected_rows),
                            "elapsed_sec": max(0.0, time.time() - started),
                            "failure_reason": "timeout_split",
                            "timed_out": True,
                            "returncode": returncode,
                            "timeout_sec": timeout_sec_observed,
                            "split_candidate_count": int(split_candidate_count),
                            "split_min": int(split_min),
                            "child_shard_ids": [
                                str(child.get("shard_id", "")) for child in split_children
                            ],
                            "child_count": int(appended),
                            "threads": max(1, int(threads)),
                        },
                    )
                    logger.warning(
                        "[scorch-shard.result] run_id=%s shard_id=%s status=split reason=timeout attempt=%d allowed=%d children=%d timeout_s=%.1f",
                        run_id,
                        shard_id,
                        int(attempt),
                        int(len(allowed_bases)),
                        int(appended),
                        timeout_sec_observed,
                    )
                    return False, combo
            status = "completed" if shard_ok else "failed"
            if (
                not shard_ok
                and timed_out
                and (int(attempt) >= max_attempts or global_ligand_quarantine_now)
                and split_candidate_count < split_min
            ):
                status = "quarantined"
                if global_ligand_quarantine_now:
                    failure_reason = "global_ligand_timeout"
                quarantine_path, quarantine_rows = _write_quarantine_output(
                    record,
                    reason=failure_reason or "timeout",
                    timeout_sec=timeout_sec_observed,
                    returncode=returncode,
                )
                if quarantine_path is not None and quarantine_rows > 0:
                    out_path = quarantine_path
                    row_count = int(quarantine_rows)
                    output_valid = _output_valid(record)[0]
                    row_coverage_ok = bool(output_valid)
            write_shard_result(
                cfg,
                run_id,
                record,
                status=status,
                payload={
                    "attempt": int(attempt),
                    "allowed_count": int(len(allowed_bases)),
                    "expected_row_count": int(expected_rows),
                    "elapsed_sec": max(0.0, time.time() - started),
                    "output_csv": "" if out_path is None else str(out_path),
                    "row_count": int(row_count),
                    "no_output": bool(out_path is None),
                    "productive": bool(status == "completed" and out_path is not None and row_count > 0),
                    "row_coverage_ok": bool(row_coverage_ok or no_output_ok),
                    "quarantine_output": bool(status == "quarantined" and out_path is not None and row_count > 0),
                    "quarantined_ligand_count": int(row_count if status == "quarantined" else 0),
                    "threads": max(1, int(threads)),
                    "returncode": returncode,
                    "failure_reason": failure_reason,
                    "timed_out": bool(timed_out),
                    "timeout_sec": timeout_sec_observed,
                    "stderr_tail": str(stage_result.get("stderr_tail") or "")[-2000:],
                    "global_ligand_quarantine": bool(global_ligand_quarantine_now),
                    "ligand_timeout_count": int(ligand_timeout_info.get("timeout_count") or 0),
                },
            )
            logger.info(
                "[scorch-shard.result] run_id=%s shard_id=%s status=%s attempt=%d rows=%d expected=%d output=%s reason=%s returncode=%s timeout_s=%.1f",
                run_id,
                shard_id,
                status,
                int(attempt),
                int(row_count),
                int(expected_rows),
                "" if out_path is None else str(out_path),
                failure_reason or "none",
                "" if returncode is None else str(returncode),
                timeout_sec_observed,
            )
            return bool(shard_ok), combo
        except Exception as exc:
            write_shard_result(
                cfg,
                run_id,
                record,
                status="failed",
                payload={
                    "attempt": int(attempt),
                    "elapsed_sec": max(0.0, time.time() - started),
                    "error": str(exc),
                },
            )
            logger.error(
                "[scorch-shard.result] run_id=%s shard_id=%s status=failed reason=exception error=%s",
                run_id,
                shard_id,
                exc,
            )
            return False, combo
        finally:
            stop_event.set()
            hb.join(timeout=1.0)

    def _worker_loop() -> None:
        nonlocal active_local
        while True:
            with records_lock:
                claimed = _claim_next_record(
                    cfg,
                    run_id,
                    list(records),
                    lease_sec=lease_sec,
                    max_attempts=max_attempts,
                    overwrite=overwrite,
                    logger=logger,
                    skip_shard_ids=set(active_shard_ids),
                )
                if claimed is not None:
                    active_shard_ids.add(str(claimed[0].get("shard_id", "") or ""))
            if claimed is None:
                with active_condition:
                    if active_local > 0:
                        active_condition.wait(
                            timeout=min(5.0, max(1.0, float(heartbeat_sec)))
                        )
                        continue
                return
            record, attempt = claimed
            shard_id = str(record.get("shard_id", ""))
            logger.info(
                "[scorch-shard.claim] run_id=%s shard_id=%s attempt=%d allowed=%s output=%s",
                run_id,
                shard_id,
                int(attempt),
                str(record.get("allowed_count", "")),
                str(record.get("output_csv", "")),
            )
            with active_condition:
                active_local += 1
            try:
                _run_one(record, int(attempt))
            finally:
                with records_lock:
                    active_shard_ids.discard(shard_id)
                with active_condition:
                    active_local = max(0, active_local - 1)
                    active_condition.notify_all()

    worker_count = max(1, int(jobs))
    with ThreadPoolExecutor(max_workers=worker_count, thread_name_prefix="scorch-shard") as pool:
        futures: Set[Future[Any]] = {pool.submit(_worker_loop) for _ in range(worker_count)}
        while futures:
            done, futures = wait(futures, timeout=5.0, return_when=FIRST_COMPLETED)
            for fut in done:
                fut.result()

    wait_started = time.time()
    while True:
        record_snapshot = _records_snapshot()
        completed, failed, inflight, pending = _summarize_records(
            cfg,
            run_id,
            record_snapshot,
            max_attempts=max_attempts,
        )
        logger.info(
            "[scorch-shard.summary] run_id=%s total=%d completed=%d failed=%d inflight=%d pending=%d workers=%d",
            run_id,
            len(record_snapshot),
            completed,
            failed,
            inflight,
            pending,
            worker_count,
        )
        if completed >= len(record_snapshot):
            break
        if failed > 0 and pending == 0 and inflight == 0:
            break
        if pending > 0 and inflight == 0:
            _worker_loop()
            continue
        if drain_timeout_sec > 0 and time.time() - wait_started >= drain_timeout_sec:
            logger.warning(
                "[scorch-shard.summary] run_id=%s status=timeout total=%d completed=%d inflight=%d pending=%d timeout_s=%.1f",
                run_id,
                len(record_snapshot),
                completed,
                inflight,
                pending,
                float(drain_timeout_sec),
            )
            break
        time.sleep(min(10.0, max(1.0, float(heartbeat_sec))))

    record_snapshot = _records_snapshot()
    repaired_split_children = repair_missing_split_child_plan_records(
        cfg,
        run_id,
        records=record_snapshot,
        logger=logger,
    )
    if repaired_split_children > 0:
        record_snapshot = _records_for_existing_scope(cfg, run_id, record_snapshot)
        with records_lock:
            records = list(record_snapshot)
        logger.info(
            "[scorch-shard.repair] run_id=%s action=refresh_records_after_workers reason=missing_split_children count=%d scoped_total=%d",
            run_id,
            int(repaired_split_children),
            len(record_snapshot),
        )
    repair_missing_completed_shard_results(cfg, run_id, records=record_snapshot, logger=logger)
    summary_payload = summarize_shard_ledger(
        cfg,
        run_id,
        records=record_snapshot,
        max_attempts=max_attempts,
        workers=worker_count,
        threads=max(1, int(threads)),
        write=False,
    )
    summarize_shard_ledger(
        cfg,
        run_id,
        max_attempts=max_attempts,
        workers=worker_count,
        threads=max(1, int(threads)),
        write=True,
    )
    failed = int(summary_payload.get("failed") or 0)
    completed = int(summary_payload.get("completed") or 0)
    total = int(summary_payload.get("total") or 0)
    all_complete = bool(summary_payload.get("all_complete", False))
    if failed > 0:
        for entry in summary_payload.get("shards", []) or []:
            if not isinstance(entry, Mapping):
                continue
            if str(entry.get("ledger_state", "") or "") != "failed":
                continue
            combo_raw = entry.get("combo") if isinstance(entry.get("combo"), list) else []
            combo = (
                str(combo_raw[0] if len(combo_raw) > 0 else ""),
                str(combo_raw[1] if len(combo_raw) > 1 else ""),
                str(combo_raw[2] if len(combo_raw) > 2 else ""),
            )
            combo_failed_local[combo] = True
    logger.info(
        "[scorch-shard.utilization] run_id=%s total=%d productive=%d zero_output=%d split=%d quarantined=%d rows=%d workers=%d threads=%d thread_seconds=%.2f",
        run_id,
        total,
        int(summary_payload.get("productive") or 0),
        int(summary_payload.get("zero_output") or 0),
        int(summary_payload.get("split") or 0),
        int(summary_payload.get("quarantined") or 0),
        int(summary_payload.get("row_total") or 0),
        worker_count,
        max(1, int(threads)),
        float(summary_payload.get("thread_seconds_total") or 0.0),
    )
    return ShardedExecutionResult(
        failed_jobs=int(failed + (0 if all_complete else 1)),
        completed=int(completed),
        total=int(total),
        all_complete=bool(all_complete),
        combo_failed_local=combo_failed_local,
    )
