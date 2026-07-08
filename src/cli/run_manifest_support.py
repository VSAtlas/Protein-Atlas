from __future__ import annotations

import datetime
import hashlib
import json
import logging
import os
import platform
import re
import subprocess
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterator, Mapping, MutableMapping, Optional, Sequence

from cli.cli_utils import _norm_pdb_id

try:
    import fcntl
except Exception:  # pragma: no cover
    fcntl = None  # type: ignore[assignment]

try:  # PyYAML is optional; we degrade gracefully if unavailable
    import yaml  # type: ignore[import-untyped]
except Exception:  # pragma: no cover - import guard
    yaml = None  # type: ignore[assignment]

STAGE_KEYS = ("prep", "pocket_detection", "docking", "postprocessing")
_INT_RE = re.compile(r"(\d+)")
_missing_yaml_logged = False


@dataclass
class PocketDetectionEvent:
    run_id: str
    pdb_id: str
    variant_label: Optional[str]
    ph_tag: Optional[str]
    method: Optional[str]
    center: Optional[Any]
    box_size: Optional[Any]


def utc_now_iso() -> str:
    return datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"


def require_yaml() -> Optional[Any]:
    global _missing_yaml_logged
    if yaml is None:
        if not _missing_yaml_logged:
            logging.warning("[run-manifest] PyYAML not available; manifest writes disabled")
            _missing_yaml_logged = True
        return None
    return yaml


def default_stage_entry() -> Dict[str, Any]:
    return {
        "status": "pending",
        "error": None,
        "timing": {"started_at": None, "finished_at": None, "wall_time_sec": None},
        "details": {},
    }


def default_protein_entry() -> Dict[str, Any]:
    return {
        "pdb_id": None,
        "library": None,
        "variant": None,
        "ph": None,
        "status": "pending",
        "error": None,
        "timing": {"started_at": None, "finished_at": None, "wall_time_sec": None},
        "stages": {k: default_stage_entry() for k in STAGE_KEYS},
    }


def normalize_pdb_id_token(token: Any) -> Optional[str]:
    if token is None:
        return None
    try:
        return _norm_pdb_id(str(token).strip())
    except Exception:
        return None


def normalize_string_token(value: Any) -> Optional[str]:
    if value is None:
        return None
    try:
        text = str(value).strip()
    except Exception:
        return None
    if len(text) >= 2 and text[0] == text[-1] and text[0] in {'"', "'"}:
        text = text[1:-1].strip()
    return text or None


def coerce_bool_token(value: Any) -> Optional[bool]:
    if isinstance(value, bool):
        return value
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return bool(value)
    try:
        token = str(value).strip().lower()
    except Exception:
        return None
    if not token:
        return None
    if token in {"1", "true", "yes", "on", "y", "t"}:
        return True
    if token in {"0", "false", "no", "off", "n", "f"}:
        return False
    return None


def refresh_summary(manifest: MutableMapping[str, Any]) -> None:
    summary = manifest.get("summary")
    if not isinstance(summary, MutableMapping):
        summary = {}
        manifest["summary"] = summary
    existing_protein_list = summary.get("total_protein_list")

    proteins = manifest.get("proteins")
    if not isinstance(proteins, MutableMapping):
        manifest["proteins"] = {}
        if "total_proteins_scheduled" not in summary:
            summary["total_proteins_scheduled"] = 0
        if "total_protein_list" not in summary:
            summary["total_protein_list"] = []
        summary["total_proteins_completed"] = 0
        summary["total_proteins_failed"] = 0
        return
    if not proteins:
        if isinstance(existing_protein_list, list):
            summary["total_protein_list"] = existing_protein_list
            summary["total_proteins_scheduled"] = len(existing_protein_list)
        else:
            summary.setdefault("total_protein_list", [])
            summary.setdefault("total_proteins_scheduled", 0)
        summary["total_proteins_completed"] = 0
        summary["total_proteins_failed"] = 0
        return

    try:
        ph_pairs: set[tuple[str, str]] = set()
        base_entries: dict[tuple[str, str], str] = {}
        by_pair: dict[tuple[str, str], list[str]] = {}
        for key, entry in proteins.items():
            raw_key = str(key)
            parts = raw_key.split("|", 2)
            if len(parts) != 3:
                continue
            pdb_part, variant_part, ph_token = parts
            pair = (pdb_part, variant_part)
            by_pair.setdefault(pair, []).append(raw_key)
            if ph_token and ph_token != "base":
                ph_pairs.add(pair)
            elif ph_token == "base":
                base_entries[pair] = raw_key

        for pdb_part, variant_part in ph_pairs:
            base_key = base_entries.get((pdb_part, variant_part))
            if not base_key:
                continue
            base_entry = proteins.get(base_key)
            if not isinstance(base_entry, MutableMapping):
                continue
            base_stages = base_entry.get("stages")
            base_details = None
            if isinstance(base_stages, MutableMapping):
                base_pocket = base_stages.get("pocket_detection")
                if isinstance(base_pocket, MutableMapping):
                    maybe_details = base_pocket.get("details")
                    if isinstance(maybe_details, MutableMapping) and maybe_details:
                        base_details = maybe_details

            if base_details is not None:
                for raw_key in by_pair.get((pdb_part, variant_part), []):
                    if raw_key == base_key:
                        continue
                    entry = proteins.get(raw_key)
                    if not isinstance(entry, MutableMapping):
                        continue
                    stages = entry.setdefault("stages", {})
                    if not isinstance(stages, MutableMapping):
                        entry["stages"] = {}
                        stages = entry["stages"]
                    pocket = stages.get("pocket_detection")
                    if not isinstance(pocket, MutableMapping):
                        pocket = default_stage_entry()
                        stages["pocket_detection"] = pocket
                    details = pocket.get("details")
                    if not isinstance(details, MutableMapping):
                        details = {}
                    merged = dict(base_details)
                    merged.update(details)
                    pocket["details"] = merged
            proteins.pop(base_key, None)

        grouped_status: dict[tuple[str, str], list[tuple[str, str]]] = {}
        for key, entry in proteins.items():
            if not isinstance(entry, Mapping):
                continue
            parts = str(key).split("|", 2)
            if len(parts) != 3:
                continue
            pdb_part, variant_part, ph_token = parts
            bucket_key = (
                str(entry.get("pdb_id") or pdb_part).upper(),
                str(entry.get("variant") or variant_part).upper(),
            )
            status = str(entry.get("status") or "").lower()
            grouped_status.setdefault(bucket_key, []).append((ph_token, status))

        scheduled_keys: set[tuple[str, str]] = set()
        completed = 0
        failed = 0
        for bucket_key, statuses in grouped_status.items():
            scheduled_keys.add(bucket_key)
            ph_statuses = [status for ph_token, status in statuses if ph_token != "base"]
            base_statuses = [status for ph_token, status in statuses if ph_token == "base"]
            effective_statuses = ph_statuses or base_statuses
            if not effective_statuses:
                continue
            if any(status == "failed" for status in effective_statuses):
                failed += 1
            elif all(status == "completed" for status in effective_statuses):
                completed += 1

        summary["total_proteins_scheduled"] = len(scheduled_keys)
        if isinstance(existing_protein_list, list) and existing_protein_list:
            summary["total_protein_list"] = existing_protein_list
        else:
            summary["total_protein_list"] = sorted(
                pdb_id for pdb_id, _variant in scheduled_keys
            )
        summary["total_proteins_completed"] = completed
        summary["total_proteins_failed"] = failed
    except Exception:
        logging.warning("[run-manifest] Failed to refresh manifest summary", exc_info=True)


def parse_int(raw: Any) -> Optional[int]:
    if raw is None:
        return None
    match = _INT_RE.search(str(raw))
    if not match:
        return None
    try:
        value = int(match.group(1))
    except Exception:
        return None
    return value if value >= 0 else None


def distributed_mode(cfg: Mapping[str, Any]) -> str:
    env_mode = os.environ.get("ATLAS_DISTRIBUTED_MODE")
    raw = env_mode if env_mode is not None else cfg.get("DISTRIBUTED_MODE")
    token = str(raw or "").strip().lower().replace("-", "_")
    explicit_env_off = env_mode is not None and token in {
        "",
        "off",
        "none",
        "false",
        "0",
        "disabled",
    }
    if token in {"", "off", "none", "false", "0", "disabled"}:
        slurm_task_count = (
            parse_int(os.environ.get("SLURM_ARRAY_TASK_COUNT"))
            or parse_int(os.environ.get("ATLAS_DIST_TASK_COUNT"))
            or 1
        )
        has_task_id = any(
            str(os.environ.get(name) or "").strip()
            for name in ("SLURM_ARRAY_TASK_ID", "ATLAS_DIST_TASK_ID")
        )
        if not explicit_env_off and int(slurm_task_count) > 1 and has_task_id:
            return "slurm_array"
        return "off"
    if token in {"slurm", "array", "slurm_array"}:
        return "slurm_array"
    return token


def distributed_task_count() -> int:
    return (
        parse_int(os.environ.get("SLURM_ARRAY_TASK_COUNT"))
        or parse_int(os.environ.get("ATLAS_DIST_TASK_COUNT"))
        or 1
    )


def distributed_task_id() -> int:
    return (
        parse_int(os.environ.get("SLURM_ARRAY_TASK_ID"))
        or parse_int(os.environ.get("ATLAS_DIST_TASK_ID"))
        or 0
    )


def distributed_enabled(cfg: Mapping[str, Any]) -> bool:
    return distributed_mode(cfg) == "slurm_array" and distributed_task_count() > 1


def suppress_distributed_base_write(cfg: Mapping[str, Any], *, ph_tag: Optional[str]) -> bool:
    if ph_tag is not None:
        return False
    if not distributed_enabled(cfg):
        return False
    return bool(coerce_bool_token(cfg.get("PH_ENSEMBLE")))


def distributed_protein_state_relpath(protein_key: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9._-]+", "_", protein_key).strip("_")
    if not sanitized:
        sanitized = "protein"
    digest = hashlib.sha1(protein_key.encode("utf-8")).hexdigest()[:16]
    return f"{sanitized[:64]}__{digest}.json"


def manifest_has_proteins(manifest: Optional[Mapping[str, Any]]) -> bool:
    if not isinstance(manifest, Mapping):
        return False
    proteins = manifest.get("proteins")
    return isinstance(proteins, Mapping) and len(proteins) > 0


def latest_mtime(paths: Sequence[Path]) -> float:
    latest = 0.0
    for path in paths:
        try:
            mtime = float(path.stat().st_mtime)
        except Exception:
            continue
        if mtime > latest:
            latest = mtime
    return latest


def file_sha1(path: Path) -> Optional[str]:
    try:
        if not path.exists() or not path.is_file():
            return None
        digest = hashlib.sha1()
        with path.open("rb") as fh:
            while True:
                chunk = fh.read(1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)
        return digest.hexdigest()
    except Exception:
        return None


@contextmanager
def exclusive_file_lock(lock_path: Path) -> Iterator[None]:
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+", encoding="utf-8") as handle:
        if fcntl is not None:
            try:
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            except Exception:
                pass
        try:
            yield
        finally:
            if fcntl is not None:
                try:
                    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
                except Exception:
                    pass


def git_info(repo_root: Path) -> Dict[str, Any]:
    info: Dict[str, Any] = {"repo": None, "branch": None, "commit": None, "dirty": None}

    def _run(cmd: list[str]) -> Optional[str]:
        try:
            res = subprocess.run(
                cmd,
                cwd=str(repo_root),
                capture_output=True,
                text=True,
                check=False,
            )
            if res.returncode == 0:
                out = res.stdout.strip()
                return out or None
        except Exception:
            return None
        return None

    try:
        info["repo"] = _run(["git", "config", "--get", "remote.origin.url"])
        info["branch"] = _run(["git", "rev-parse", "--abbrev-ref", "HEAD"])
        info["commit"] = _run(["git", "rev-parse", "HEAD"])
        dirty_output = _run(["git", "status", "--porcelain"])
        info["dirty"] = bool(dirty_output) if dirty_output is not None else None
    except Exception:
        logging.warning(
            "[run-manifest] Failed to gather git metadata root=%s",
            repo_root,
            exc_info=True,
        )
    return info


def resources_snapshot() -> Dict[str, Any]:
    host = None
    try:
        host = os.uname().nodename
    except Exception:
        try:
            host = platform.node()
        except Exception:
            host = None

    try:
        n_cores = os.cpu_count()
    except Exception:
        n_cores = None

    ram_gb = None
    try:
        page_size = os.sysconf("SC_PAGE_SIZE")
        phys_pages = os.sysconf("SC_PHYS_PAGES")
        ram_bytes = float(page_size) * float(phys_pages)
        ram_gb = round(ram_bytes / (1024**3), 2)
    except Exception:
        ram_gb = None

    gpu = None
    try:
        cuda_visible = os.environ.get("CUDA_VISIBLE_DEVICES")
        if cuda_visible:
            gpu = cuda_visible
    except Exception:
        gpu = None

    slurm = {
        "job_id": os.environ.get("SLURM_JOB_ID"),
        "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
        "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        "array_task_count": os.environ.get("SLURM_ARRAY_TASK_COUNT"),
        "cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
        "cpus_on_node": os.environ.get("SLURM_CPUS_ON_NODE"),
        "job_name": os.environ.get("SLURM_JOB_NAME"),
    }
    slurm = {key: value for key, value in slurm.items() if value}

    return {
        "host": host,
        "n_cores": n_cores,
        "ram_gb": ram_gb,
        "gpu": gpu,
        "slurm": slurm,
    }


def normalize_for_hash(obj: Any) -> Any:
    if isinstance(obj, Mapping):
        return {
            str(k): normalize_for_hash(v)
            for k, v in sorted(obj.items(), key=lambda kv: str(kv[0]))
        }
    if isinstance(obj, (list, tuple, set)):
        return [normalize_for_hash(v) for v in obj]
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return repr(obj)


def compute_config_hash(cfg: Mapping[str, Any]) -> str:
    try:
        normalized = normalize_for_hash(dict(cfg))
        payload = json.dumps(normalized, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()
    except Exception:
        try:
            return hashlib.sha256(repr(cfg).encode("utf-8")).hexdigest()
        except Exception:
            return "UNKNOWN"


def extract_error_from_fail_log(fail_log_path: Optional[Path]) -> str:
    try:
        if not fail_log_path or not fail_log_path.exists():
            return "<unknown error>"
        with fail_log_path.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                stripped = line.lstrip()
                if stripped.startswith("exception     ="):
                    parts = stripped.split("=", 1)
                    if len(parts) == 2:
                        return parts[1].strip()
        return "<unknown error>"
    except Exception:
        return "<unknown error>"


def coerce_vec3(value: Any) -> Optional[list[float]]:
    if value is None:
        return None
    try:
        seq = list(value)
    except Exception:
        return None
    if len(seq) < 3:
        return None
    try:
        return [float(seq[0]), float(seq[1]), float(seq[2])]
    except Exception:
        return None
