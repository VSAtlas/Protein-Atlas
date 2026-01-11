"""
Best-effort run manifest writer for Atlas runs.

All helpers must remain non-fatal: any error (missing PyYAML, IO, git)
is logged at WARNING level and silently skipped so scientific behavior
is unchanged.
"""

from __future__ import annotations

import datetime
import json
import logging
import os
import platform
import subprocess
import tempfile
import threading
import time
from dataclasses import asdict, dataclass
from pathlib import Path
import hashlib
from typing import Any, Dict, Mapping, MutableMapping, Optional, Sequence, Tuple

try:  # PyYAML is optional; we degrade gracefully if unavailable
    import yaml
except Exception:  # pragma: no cover - import guard
    yaml = None  # type: ignore[assignment]


STAGE_KEYS = ("prep", "pocket_detection", "docking", "postprocessing")
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


def _utc_now_iso() -> str:
    """UTC timestamp without microseconds for stable manifest entries."""

    return datetime.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"


def _require_yaml() -> Optional[Any]:
    """Return the yaml module or log once and skip writes if missing."""

    global _missing_yaml_logged
    if yaml is None:
        if not _missing_yaml_logged:
            logging.warning(
                "[run-manifest] PyYAML not available; manifest writes disabled"
            )
            _missing_yaml_logged = True
        return None
    return yaml


def _default_stage_entry() -> Dict[str, Any]:
    return {
        "status": "pending",
        "error": None,
        "timing": {"started_at": None, "finished_at": None, "wall_time_sec": None},
        # Arbitrary stage-specific metadata (e.g., pocket detection info)
        "details": {},
    }


def _default_protein_entry() -> Dict[str, Any]:
    return {
        "pdb_id": None,
        "library": None,
        "variant": None,
        "ph": None,
        "status": "pending",
        "error": None,
        "timing": {"started_at": None, "finished_at": None, "wall_time_sec": None},
        "stages": {k: _default_stage_entry() for k in STAGE_KEYS},
    }


def _normalize_pdb_id_token(token: Any) -> Optional[str]:
    """
    Best-effort normalization of arbitrary tokens into 4-char uppercase PDB IDs.
    Mirrors the CLI helper in main.py without importing it to avoid cycles.
    """
    if token is None:
        return None
    try:
        t = str(token).strip()
    except Exception:
        return None
    if not t:
        return None
    while t.startswith("-"):
        t = t[1:]
    t = os.path.basename(t)
    if t.lower().endswith(".pdb"):
        t = t[:-4]
    t = t.replace("_cleaned", "")
    t = t.upper()
    if len(t) >= 4:
        cand = t[:4]
        return cand if cand.isalnum() else None
    return None


def _normalize_string_token(value: Any) -> Optional[str]:
    if value is None:
        return None
    try:
        s = str(value).strip()
    except Exception:
        return None
    if len(s) >= 2 and s[0] == s[-1] and s[0] in {'"', "'"}:
        s = s[1:-1].strip()
    return s or None


def _coerce_bool_token(value: Any) -> Optional[bool]:
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
    truthy = {"1", "true", "yes", "on", "y", "t"}
    falsy = {"0", "false", "no", "off", "n", "f"}
    if token in truthy:
        return True
    if token in falsy:
        return False
    return None


def _refresh_summary(manifest: MutableMapping[str, Any]) -> None:
    summary = manifest.get("summary")
    if not isinstance(summary, MutableMapping):
        summary = {}
        manifest["summary"] = summary

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

    try:
        # Phase A: scan entries and classify base vs pH-tagged
        ph_pairs: set[tuple[str, str]] = set()
        base_entries: dict[tuple[str, str], str] = {}
        for key, entry in proteins.items():
            raw_key = str(key)
            parts = raw_key.split("|", 2)
            if len(parts) != 3:
                continue
            pdb_part, variant_part, ph_token = parts
            if ph_token and ph_token != "base":
                ph_pairs.add((pdb_part, variant_part))
            elif ph_token == "base":
                base_entries[(pdb_part, variant_part)] = raw_key

        # Phase B: promote pocket-detection details from base to pH entries
        for pdb_part, variant_part in ph_pairs:
            base_key = base_entries.get((pdb_part, variant_part))
            if not base_key:
                continue

            base_entry = proteins.get(base_key)
            if not isinstance(base_entry, MutableMapping):
                continue

            base_stages = base_entry.get("stages") or {}
            if not isinstance(base_stages, MutableMapping):
                base_stages = {}

            base_pocket = base_stages.get("pocket_detection") or {}
            if not isinstance(base_pocket, MutableMapping):
                base_pocket = {}
            base_details = base_pocket.get("details") or {}
            if not isinstance(base_details, MutableMapping) or not base_details:
                continue

            promoted_ph_tags: list[str] = []
            for key2, entry2 in proteins.items():
                raw_key2 = str(key2)
                parts2 = raw_key2.split("|", 2)
                if len(parts2) != 3:
                    continue
                pdb2, var2, ph2 = parts2
                if (pdb2, var2) != (pdb_part, variant_part):
                    continue
                if ph2 == "base":
                    continue

                if not isinstance(entry2, MutableMapping):
                    continue

                stages2 = entry2.setdefault("stages", {})
                if not isinstance(stages2, MutableMapping):
                    stages2 = {}
                    entry2["stages"] = stages2

                pocket2 = stages2.get("pocket_detection")
                if not isinstance(pocket2, MutableMapping):
                    pocket2 = _default_stage_entry()
                    stages2["pocket_detection"] = pocket2

                details2 = pocket2.get("details")
                if not isinstance(details2, MutableMapping):
                    details2 = {}
                pocket2["details"] = details2

                for dk, dv in base_details.items():
                    details2.setdefault(dk, dv)

                if base_pocket.get("status") == "completed":
                    pocket2.setdefault("status", "completed")
                    pocket2.setdefault("error", None)

                promoted_ph_tags.append(ph2)

            if promoted_ph_tags:
                try:
                    logging.debug(
                        "[run-manifest.pocket_detection.promote] pdb=%s variant=%s ph_tags=%s",
                        pdb_part,
                        variant_part,
                        sorted(set(promoted_ph_tags)),
                    )
                except Exception:
                    logging.warning(
                        "[run-manifest] Failed to log pocket_detection promotion pdb=%s variant=%s",
                        pdb_part,
                        variant_part,
                    )

        # Phase C: drop base entries when pH-tagged entries exist
        base_keys_to_drop: list[str] = []
        for pdb_part, variant_part in ph_pairs:
            base_key = base_entries.get((pdb_part, variant_part))
            if base_key:
                base_keys_to_drop.append(base_key)

        for k in base_keys_to_drop:
            proteins.pop(k, None)

        total = len(proteins)
        completed = 0
        failed = 0
        for entry in proteins.values():
            status = entry.get("status") if isinstance(entry, Mapping) else None
            if status == "completed":
                completed += 1
            elif status == "failed":
                failed += 1

        if "total_proteins_scheduled" not in summary:
            summary["total_proteins_scheduled"] = total
        if "total_protein_list" not in summary:
            fallback_ids: list[str] = []
            for entry in proteins.values():
                if not isinstance(entry, Mapping):
                    continue
                nid = entry.get("pdb_id")
                canon = _normalize_pdb_id_token(nid)
                if canon:
                    fallback_ids.append(canon)
            summary["total_protein_list"] = sorted(set(fallback_ids))
        summary["total_proteins_completed"] = completed
        summary["total_proteins_failed"] = failed
    except Exception:
        logging.warning(
            "[run-manifest] Failed to refresh summary (promotion + dedupe)",
            exc_info=True,
        )


def _load_manifest(manifest_path: Path) -> Optional[Dict[str, Any]]:
    yaml_mod = _require_yaml()
    if yaml_mod is None:
        return None

    if not manifest_path.exists():
        return {}

    try:
        with manifest_path.open("r", encoding="utf-8") as fh:
            data = yaml_mod.safe_load(fh) or {}
        if not isinstance(data, dict):
            return {}
        return data
    except Exception:
        logging.warning(
            "[run-manifest] Failed to load manifest path=%s",
            manifest_path,
            exc_info=True,
        )
        return {}


def _write_manifest(manifest_path: Path, manifest: Dict[str, Any]) -> None:
    yaml_mod = _require_yaml()
    if yaml_mod is None:
        return

    tmp_path: Optional[Path] = None
    try:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)

        # Unique temp per process/thread to avoid replace collisions under parallel writes
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=manifest_path.name + ".",
            suffix=f".{os.getpid()}.{threading.get_ident()}.tmp",
            dir=manifest_path.parent,
            delete=False,
        ) as fh:
            tmp_path = Path(fh.name)
            yaml_mod.safe_dump(manifest, fh, default_flow_style=False, sort_keys=False)
            fh.flush()
            os.fsync(fh.fileno())

        tmp_path.replace(manifest_path)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to write manifest path=%s",
            manifest_path,
            exc_info=True,
        )
    finally:
        try:
            if tmp_path and tmp_path.exists():
                tmp_path.unlink()
        except Exception:
            pass


def _protein_key(
    pdb_id: str, variant_label: Optional[str], ph_tag: Optional[str]
) -> str:
    variant_norm = (variant_label or "legacy").strip().upper() or "LEGACY"
    ph_norm = (ph_tag or "base").strip()
    ph_token = ph_norm if ph_norm else "base"
    return f"{str(pdb_id).upper()}|{variant_norm}|{ph_token}"


def get_manifest_paths(cfg: Mapping[str, Any], run_id: str) -> Tuple[Path, Path]:
    """
    Return (manifest_dir, manifest_path) for the current run and ensure the
    directory exists. Never raises.
    """

    try:
        repo_root = Path(cfg.get("OVERALL_DIR", "."))
        manifest_dir = repo_root / "manifests" / str(run_id)
        manifest_dir.mkdir(parents=True, exist_ok=True)
        manifest_path = manifest_dir / "run_manifest.yaml"
        return manifest_dir, manifest_path
    except Exception:
        logging.warning(
            "[run-manifest] Failed to prepare manifest paths run_id=%s",
            run_id,
            exc_info=True,
        )
        fallback_dir = Path(cfg.get("OVERALL_DIR", ".")) / "manifests"
        return fallback_dir, fallback_dir / "run_manifest.yaml"


def _git_info(repo_root: Path) -> Dict[str, Any]:
    info = {"repo": None, "branch": None, "commit": None, "dirty": None}
    try:

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


def update_manifest_for_druggability_and_engine_plan(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: Optional[str],
    *,
    ph_tag: Optional[str],
    tier: Optional[str],
    use_gnina: bool,
    use_ledock: bool,
    use_dock6: bool,
) -> None:
    """
    Record fpocket druggability tier (A/B/C) and planned engines for (pdb, variant, pH).

    Best-effort only: errors are logged and ignored.
    """
    ph_label = ph_tag if ph_tag is not None else "base"
    try:
        if not run_id:
            logging.debug(
                "[run-manifest.druggability-plan.skip] no run_id pdb=%s variant=%s ph=%s",
                pdb_id,
                variant_label,
                ph_label,
            )
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            logging.warning(
                "[run-manifest.druggability-plan.skip] manifest missing run_id=%s path=%s pdb=%s variant=%s ph=%s",
                run_id,
                manifest_path,
                pdb_id,
                variant_label,
                ph_label,
            )
            return

        entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
        stages = entry.setdefault("stages", {})

        pocket_stage = stages.get("pocket_detection")
        if not isinstance(pocket_stage, MutableMapping):
            pocket_stage = _default_stage_entry()
            stages["pocket_detection"] = pocket_stage

        pocket_details = pocket_stage.get("details")
        if not isinstance(pocket_details, MutableMapping):
            pocket_details = {}
        pocket_stage["details"] = pocket_details

        tier_norm = None
        if tier is not None:
            try:
                t = str(tier).strip().upper()
            except Exception:
                t = ""
            if t in {"A", "B", "C"}:
                tier_norm = t

        if tier_norm is not None:
            pocket_details["druggability_tier"] = tier_norm

        docking_stage = stages.get("docking")
        if not isinstance(docking_stage, MutableMapping):
            docking_stage = _default_stage_entry()
            stages["docking"] = docking_stage

        docking_details = docking_stage.get("details")
        if not isinstance(docking_details, MutableMapping):
            docking_details = {}
        docking_stage["details"] = docking_details

        engine_plan = docking_details.get("engine_plan")
        if not isinstance(engine_plan, MutableMapping):
            engine_plan = {}
        docking_details["engine_plan"] = engine_plan

        if tier_norm is not None:
            engine_plan["tier"] = tier_norm

        engines = {
            "vina": True,
            "gnina": bool(use_gnina),
            "ledock": bool(use_ledock),
            "dock6": bool(use_dock6),
        }
        engine_plan["engines"] = engines

        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)

        logging.debug(
            "[run-manifest.druggability-plan.ok] run_id=%s pdb=%s variant=%s ph=%s tier=%s engines=%r",
            run_id,
            pdb_id,
            variant_label,
            ph_label,
            tier_norm,
            engines,
        )
    except Exception:
        logging.warning(
            "[run-manifest.druggability-plan.error] run_id=%s pdb=%s variant=%s ph=%s",
            run_id,
            pdb_id,
            variant_label,
            ph_label,
            exc_info=True,
        )


def _resources_snapshot() -> Dict[str, Any]:
    host = None
    try:
        host = os.uname().nodename
    except Exception:
        try:
            host = platform.node()
        except Exception:
            host = None

    n_cores = None
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

    return {"host": host, "n_cores": n_cores, "ram_gb": ram_gb, "gpu": gpu}


def init_run_manifest(
    cfg: Mapping[str, Any], run_id: str, argv: list[str], log_path: str
) -> None:
    try:
        manifest_dir, manifest_path = get_manifest_paths(cfg, run_id)
        yaml_mod = _require_yaml()
        if yaml_mod is None:
            return

        run_dir = Path(cfg.get("CONFIG_RUN_DIR", ""))
        overall_dir = Path(cfg.get("OVERALL_DIR", "."))
        selection_mode = None
        test_mode = str(cfg.get("TEST_MODE_ENABLE", "off"))

        if test_mode != "off":
            selection_mode = "TEST_LIBRARY_MAP"
        elif cfg.get("_EFFECTIVE_SPECIFIED_PROTEINS"):
            selection_mode = "SPECIFIED_PROTEINS"
        else:
            selection_mode = "AUTO"

        pdb_list = None
        if selection_mode == "SPECIFIED_PROTEINS":
            pdb_list = list(cfg.get("_EFFECTIVE_SPECIFIED_PROTEINS", []) or [])

        manifest = {
            "run_id": str(run_id),
            "status": "running",
            "command": {
                "argv": " ".join(argv),
                "pdb_selection_mode": selection_mode,
                "TEST_MODE_ENABLE": test_mode,
                "PH_ENSEMBLE_ENABLE": bool(cfg.get("PH_ENSEMBLE", False)),
            },
            "git": _git_info(overall_dir),
            "paths": {
                "run_dir": str(run_dir) if run_dir else None,
                "log_file": str(log_path),
                "docked_dir": str(cfg.get("DOCKED_DIR", overall_dir / "docked")),
                "input_pdb_dir": str(cfg.get("INPUT_DIR", overall_dir / "input_pdbs")),
                "processed_pdb_dir": str(
                    cfg.get("OUTPUT_DIR", overall_dir / "processed_pdbs")
                ),
                "prepped_ligands_dir": str(
                    cfg.get("PREPPED_LIGANDS_DIR", overall_dir / "prepped_ligands")
                ),
                "config_file": "run_config.yaml",
            },
            "timing": {
                "created_at": _utc_now_iso(),
                "started_at": _utc_now_iso(),
                "finished_at": None,
                "wall_time_sec": None,
            },
            "resources": _resources_snapshot(),
            "summary": {
                "total_proteins_scheduled": 0,
                "total_proteins_completed": 0,
                "total_proteins_failed": 0,
            },
            "proteins": {},
        }

        if pdb_list:
            manifest["command"]["pdb_list"] = pdb_list

        manifest["command"]["config_hash"] = compute_config_hash(cfg)

        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to initialize run manifest run_id=%s",
            run_id,
            exc_info=True,
        )


def _ensure_protein(
    manifest: MutableMapping[str, Any],
    pdb_id: str,
    variant_label: Optional[str],
    ph_tag: Optional[str],
) -> Dict[str, Any]:
    proteins = manifest.setdefault("proteins", {})
    if not isinstance(proteins, MutableMapping):
        manifest["proteins"] = {}
        proteins = manifest["proteins"]

    key = _protein_key(pdb_id, variant_label, ph_tag)
    entry = proteins.get(key)
    if not isinstance(entry, MutableMapping):
        entry = _default_protein_entry()
        proteins[key] = entry
    else:
        # Ensure stages/timing keys exist
        entry.setdefault("pdb_id", str(pdb_id))
        entry.setdefault("library", None)
        entry.setdefault("variant", None)
        entry.setdefault("ph", ph_tag)
        entry.setdefault("status", "pending")
        entry.setdefault("error", None)
        entry.setdefault(
            "timing", {"started_at": None, "finished_at": None, "wall_time_sec": None}
        )
        stages = entry.setdefault("stages", {})
        for stage_key in STAGE_KEYS:
            stage_entry = stages.get(stage_key)
            if not isinstance(stage_entry, MutableMapping):
                stage_entry = _default_stage_entry()
                stages[stage_key] = stage_entry
            stage_entry.setdefault("status", "pending")
            stage_entry.setdefault("error", None)
            stage_entry.setdefault(
                "timing",
                {"started_at": None, "finished_at": None, "wall_time_sec": None},
            )
            stage_entry.setdefault("details", {})

    entry["pdb_id"] = str(pdb_id).upper()
    entry["variant"] = (variant_label or "legacy").strip().upper() or "LEGACY"
    entry["ph"] = ph_tag

    return entry  # type: ignore[return-value]


def _normalize_for_hash(obj: Any) -> Any:
    if isinstance(obj, Mapping):
        return {
            str(k): _normalize_for_hash(v)
            for k, v in sorted(obj.items(), key=lambda kv: str(kv[0]))
        }
    if isinstance(obj, (list, tuple, set)):
        return [_normalize_for_hash(v) for v in obj]
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return repr(obj)


def compute_config_hash(cfg: Mapping[str, Any]) -> str:
    try:
        normalized = _normalize_for_hash(dict(cfg))
        payload = json.dumps(
            normalized,
            sort_keys=True,
            separators=(",", ":"),
        )
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()
    except Exception:
        try:
            return hashlib.sha256(repr(cfg).encode("utf-8")).hexdigest()
        except Exception:
            return "UNKNOWN"


def update_manifest_for_scheduled_proteins(
    cfg: Mapping[str, Any],
    run_id: str,
    scheduled_pdb_ids: Sequence[Any],
) -> None:
    """
    Record the canonical list of PDB IDs queued for this run.
    """
    try:
        if not run_id:
            logging.debug("[run-manifest.scheduled.skip] reason=missing_run_id")
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            logging.warning(
                "[run-manifest.scheduled.skip] manifest missing run_id=%s path=%s",
                run_id,
                manifest_path,
            )
            return

        summary = manifest.get("summary")
        if not isinstance(summary, MutableMapping):
            summary = {}
            manifest["summary"] = summary

        norm_ids: set[str] = set()
        for raw in scheduled_pdb_ids:
            nid = _normalize_pdb_id_token(raw)
            if nid:
                norm_ids.add(nid)

        sorted_ids = sorted(norm_ids)
        summary["total_proteins_scheduled"] = len(sorted_ids)
        summary["total_protein_list"] = sorted_ids

        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest.scheduled.error] run_id=%s",
            run_id,
            exc_info=True,
        )


def update_manifest_for_run_config(
    cfg: Mapping[str, Any],
    run_id: str,
) -> None:
    """
    Record resolved apo/holo mode and water-policy knobs in the manifest.
    """
    try:
        if not run_id:
            logging.debug("[run-manifest.run-config.skip] reason=missing_run_id")
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            logging.warning(
                "[run-manifest.run-config.skip] manifest missing run_id=%s path=%s",
                run_id,
                manifest_path,
            )
            return

        cmd = manifest.get("command")
        if not isinstance(cmd, MutableMapping):
            cmd = {}
            manifest["command"] = cmd

        apo_mode = _normalize_string_token(cfg.get("_RESOLVED_APO_HOLO_MODE"))
        if apo_mode is None:
            apo_mode = _normalize_string_token(cfg.get("APO_HOLO_MODE"))
        if apo_mode is not None:
            cmd["APO_HOLO_MODE"] = apo_mode

        if "REMOVE_WATERS" in cfg:
            rw_raw = cfg.get("REMOVE_WATERS")
            rw_val = _coerce_bool_token(rw_raw)
            cmd["REMOVE_WATERS"] = rw_val if rw_val is not None else rw_raw

        if "KEEP_WATERS_WITHIN_A" in cfg:
            val = cfg.get("KEEP_WATERS_WITHIN_A")
            try:
                cmd["KEEP_WATERS_WITHIN_A"] = float(val) if val is not None else None
            except Exception:
                cmd["KEEP_WATERS_WITHIN_A"] = val

        if "WATER_KEEP_POLICY" in cfg:
            val = _normalize_string_token(cfg.get("WATER_KEEP_POLICY"))
            cmd["WATER_KEEP_POLICY"] = (
                val if val is not None else cfg.get("WATER_KEEP_POLICY")
            )

        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest.run-config.error] run_id=%s",
            run_id,
            exc_info=True,
        )


def update_manifest_for_config_hash(
    cfg: Mapping[str, Any],
    run_id: str,
) -> None:
    """
    Record the resolved config hash under command.config_hash.
    """
    try:
        if not run_id:
            logging.debug("[run-manifest.config-hash.skip] missing run_id")
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            logging.warning(
                "[run-manifest.config-hash.skip] manifest missing run_id=%s path=%s",
                run_id,
                manifest_path,
            )
            return

        cmd = manifest.get("command")
        if not isinstance(cmd, MutableMapping):
            cmd = {}
            manifest["command"] = cmd

        cmd["config_hash"] = compute_config_hash(cfg)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest.config-hash.error] run_id=%s",
            run_id,
            exc_info=True,
        )


def load_run_manifest(cfg: Mapping[str, Any], run_id: str) -> Optional[Dict[str, Any]]:
    """
    Best-effort loader for run_manifest.yaml for a given run_id.

    Returns the manifest dict, or None if it can't be loaded.
    """
    try:
        if not run_id:
            logging.debug("[run-manifest.load.skip] reason=missing_run_id")
            return None

        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            logging.warning(
                "[run-manifest.load.missing] run_id=%s path=%s",
                run_id,
                manifest_path,
            )
        return manifest
    except Exception:
        logging.warning(
            "[run-manifest.load.error] run_id=%s",
            run_id,
            exc_info=True,
        )
        return None


def update_manifest_for_protein_start(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: str,
    library: Optional[str],
    *,
    ph_tag: Optional[str] = None,
) -> None:
    try:
        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            return

        entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
        entry["variant"] = (variant_label or "legacy").upper()
        if library is None:
            library = cfg.get("LIBRARY_SUBDIR_DEFAULT")
        entry["library"] = library
        entry["status"] = "running"
        entry.setdefault("ph", ph_tag)
        timing = entry.get("timing", {})
        if isinstance(timing, MutableMapping):
            if not timing.get("started_at"):
                timing["started_at"] = _utc_now_iso()
        else:
            entry["timing"] = {
                "started_at": _utc_now_iso(),
                "finished_at": None,
                "wall_time_sec": None,
            }

        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to record start pdb=%s variant=%s ph=%s run_id=%s",
            pdb_id,
            variant_label,
            ph_tag if ph_tag is not None else "base",
            run_id,
            exc_info=True,
        )


def update_manifest_for_protein_success(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: str,
    elapsed_sec: float,
    *,
    ph_tag: Optional[str] = None,
) -> None:
    try:
        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            return

        entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
        entry["status"] = "completed"
        entry["error"] = None

        timing = entry.get("timing", {})
        if isinstance(timing, MutableMapping):
            timing.setdefault("started_at", _utc_now_iso())
            timing["finished_at"] = _utc_now_iso()
            timing["wall_time_sec"] = round(float(elapsed_sec), 3)
        else:
            entry["timing"] = {
                "started_at": _utc_now_iso(),
                "finished_at": _utc_now_iso(),
                "wall_time_sec": round(float(elapsed_sec), 3),
            }

        stages = entry.get("stages", {})
        if isinstance(stages, MutableMapping):
            for stage_key in STAGE_KEYS:
                stage_entry = stages.get(stage_key)
                if not isinstance(stage_entry, MutableMapping):
                    stages[stage_key] = _default_stage_entry()
                    stage_entry = stages[stage_key]
                stage_entry["status"] = "completed"
                stage_entry.setdefault("error", None)
        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to record success pdb=%s variant=%s ph=%s run_id=%s",
            pdb_id,
            variant_label,
            ph_tag if ph_tag is not None else "base",
            run_id,
            exc_info=True,
        )


def _extract_error_from_fail_log(fail_log_path: Optional[Path]) -> str:
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


def _coerce_vec3(value: Any) -> Optional[list[float]]:
    """
    Best-effort: normalize a 3-vector (center/box) into a JSON/YAML-friendly
    [x, y, z] list of floats. Returns None on failure.
    """
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


def update_manifest_for_protein_failure(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: str,
    fail_log_path: Optional[Path],
    *,
    ph_tag: Optional[str] = None,
) -> None:
    try:
        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            return

        entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
        entry["variant"] = (variant_label or "legacy").upper()
        entry["status"] = "failed"
        entry["error"] = _extract_error_from_fail_log(fail_log_path)
        entry["fail_log"] = str(fail_log_path) if fail_log_path is not None else None

        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to record failure pdb=%s variant=%s ph=%s run_id=%s",
            pdb_id,
            variant_label,
            ph_tag if ph_tag is not None else "base",
            run_id,
            exc_info=True,
        )


def emit_pocket_detection_event(
    cfg: Mapping[str, Any], event: PocketDetectionEvent
) -> None:
    """
    Append a pocket-detection event to a JSONL log for later replay by the
    main process. Best-effort only; failures are logged as warnings.
    """
    ph_label = event.ph_tag if event.ph_tag is not None else "base"
    try:
        if not event.run_id:
            logging.debug(
                "[run-manifest.pocket_detection.event.skip] reason=missing_run_id pdb=%s variant=%s ph=%s",
                event.pdb_id,
                event.variant_label,
                ph_label,
            )
            return

        _, manifest_path = get_manifest_paths(cfg, event.run_id)
        events_path = manifest_path.with_name("run_manifest_events.jsonl")
        payload = {"type": "pocket_detection", **asdict(event)}
        events_path.parent.mkdir(parents=True, exist_ok=True)
        with events_path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(payload, sort_keys=True) + "\n")

        logging.debug(
            "[run-manifest.pocket_detection.event] run_id=%s pdb=%s variant=%s ph=%s method=%s",
            event.run_id,
            event.pdb_id,
            event.variant_label,
            ph_label,
            event.method,
        )
    except Exception:
        logging.warning(
            "[run-manifest.pocket_detection.event.error] run_id=%s pdb=%s variant=%s ph=%s",
            event.run_id,
            event.pdb_id,
            event.variant_label,
            ph_label,
            exc_info=True,
        )


def apply_pocket_detection_events(cfg: Mapping[str, Any], run_id: str) -> None:
    """
    Replay pocket-detection events for a run and write them to the manifest in
    a single-threaded main-process context. Best-effort; errors are logged and
    ignored.
    """
    if not run_id:
        return

    try:
        _, manifest_path = get_manifest_paths(cfg, run_id)
        events_path = manifest_path.with_name("run_manifest_events.jsonl")
        if not events_path.exists() or events_path.stat().st_size == 0:
            return

        try:
            with events_path.open("r", encoding="utf-8") as fh:
                lines = fh.readlines()
        except Exception:
            logging.warning(
                "[run-manifest.pocket_detection.events.read.error] run_id=%s path=%s",
                run_id,
                events_path,
                exc_info=True,
            )
            return

        events: list[dict[str, Any]] = []
        for raw in lines:
            line = raw.strip()
            if not line:
                continue
            try:
                payload = json.loads(line)
            except Exception:
                continue
            if payload.get("type") != "pocket_detection":
                continue
            if str(payload.get("run_id", "")) != str(run_id):
                continue
            events.append(payload)

        if not events:
            try:
                events_path.unlink()
            except Exception:
                pass
            return

        for payload in events:
            try:
                update_manifest_for_pocket_detection(
                    cfg,
                    run_id=str(payload.get("run_id") or run_id),
                    pdb_id=payload.get("pdb_id"),
                    variant_label=payload.get("variant_label"),
                    ph_tag=payload.get("ph_tag"),
                    method=payload.get("method"),
                    center=payload.get("center"),
                    box_size=payload.get("box_size"),
                )
            except Exception:
                logging.warning(
                    "[run-manifest.pocket_detection.events.apply.error] run_id=%s pdb=%s variant=%s ph=%s",
                    run_id,
                    payload.get("pdb_id"),
                    payload.get("variant_label"),
                    payload.get("ph_tag"),
                    exc_info=True,
                )

        try:
            events_path.unlink()
        except Exception:
            try:
                events_path.write_text("", encoding="utf-8")
            except Exception:
                pass
    except Exception:
        logging.warning(
            "[run-manifest.pocket_detection.events.error] run_id=%s",
            run_id,
            exc_info=True,
        )


def update_manifest_for_pocket_detection(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: Optional[str],
    *,
    ph_tag: Optional[str],
    method: Optional[str],
    center: Optional[Any],
    box_size: Optional[Any],
) -> None:
    """
    Record pocket-detection method + geometry under the pocket_detection stage
    for a given (pdb, variant, pH) context.

    Best-effort only: any error logs a WARNING and is otherwise ignored.
    """
    ph_label = ph_tag if ph_tag is not None else "base"
    try:
        if not run_id:
            logging.debug(
                "[run-manifest.pocket_detection.skip] no run_id pdb=%s variant=%s ph=%s",
                pdb_id,
                variant_label,
                ph_label,
            )
            return

        logging.debug(
            "[run-manifest.pocket_detection.request] run_id=%s pdb=%s variant=%s ph=%s method=%s center=%r box=%r",
            run_id,
            pdb_id,
            variant_label,
            ph_label,
            method,
            center,
            box_size,
        )

        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            logging.warning(
                "[run-manifest.pocket_detection.skip] manifest missing run_id=%s path=%s pdb=%s variant=%s ph=%s",
                run_id,
                manifest_path,
                pdb_id,
                variant_label,
                ph_label,
            )
            return

        entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)

        stages = entry.setdefault("stages", {})
        pocket_stage = stages.get("pocket_detection")
        if not isinstance(pocket_stage, MutableMapping):
            pocket_stage = _default_stage_entry()
            stages["pocket_detection"] = pocket_stage

        details = pocket_stage.get("details")
        if not isinstance(details, MutableMapping):
            details = {}
        pocket_stage["details"] = details

        existing_method = details.get("method")
        method_raw = (method or "").strip()
        if not method_raw and existing_method is not None:
            try:
                method_raw = str(existing_method).strip()
            except Exception:
                method_raw = ""
        method_str = method_raw or None
        center_vec = _coerce_vec3(center)
        box_vec = _coerce_vec3(box_size)

        if method_str is not None:
            details["method"] = method_str
        if center_vec is not None:
            details["center"] = [round(float(x), 3) for x in center_vec]
        if box_vec is not None:
            details["box_size"] = [round(float(x), 1) for x in box_vec]

        if method_str and center_vec and box_vec:
            pocket_stage.setdefault("status", "completed")
            pocket_stage.setdefault("error", None)
            logging.info(
                "[run-manifest.pocket_detection.ok] run_id=%s pdb=%s variant=%s ph=%s method=%s center=%s box=%s",
                run_id,
                pdb_id,
                variant_label,
                ph_label,
                method_str,
                details.get("center"),
                details.get("box_size"),
            )
        else:
            logging.info(
                "[run-manifest.pocket_detection.partial] run_id=%s pdb=%s variant=%s ph=%s method=%s center=%r box=%r",
                run_id,
                pdb_id,
                variant_label,
                ph_label,
                method_str,
                center,
                box_size,
            )

        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)

    except Exception:
        logging.warning(
            "[run-manifest.pocket_detection.error] run_id=%s pdb=%s variant=%s ph=%s",
            run_id,
            pdb_id,
            variant_label,
            ph_label,
            exc_info=True,
        )


def _ensure_stage_timing(
    stage_entry: MutableMapping[str, Any],
) -> MutableMapping[str, Any]:
    timing = stage_entry.get("timing")
    if not isinstance(timing, MutableMapping):
        timing = {"started_at": None, "finished_at": None, "wall_time_sec": None}
    else:
        timing.setdefault("started_at", None)
        timing.setdefault("finished_at", None)
        timing.setdefault("wall_time_sec", None)
    stage_entry["timing"] = timing
    return timing


def update_manifest_for_docking_overall(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: Optional[str],
    *,
    ph_tag: Optional[str] = None,
    event: str = "start",
    elapsed_sec: Optional[float] = None,
    error: Optional[str] = None,
) -> None:
    """
    Record overall docking timing + status for a single (pdb, variant, pH).

    event == "start" : mark running + capture started_at if unset.
    event == "end"   : mark completed + capture finished_at/wall_time_sec.
    event == "fail"  : mark failed + capture finished_at/wall_time_sec/error.
    """
    ph_label = ph_tag if ph_tag is not None else "base"
    try:
        if not run_id:
            logging.debug(
                "[run-manifest.docking-overall.skip] no run_id pdb=%s variant=%s ph=%s event=%s",
                pdb_id,
                variant_label,
                ph_label,
                event,
            )
            return

        now = _utc_now_iso()
        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            logging.warning(
                "[run-manifest.docking-overall.skip] manifest missing run_id=%s path=%s pdb=%s variant=%s ph=%s",
                run_id,
                manifest_path,
                pdb_id,
                variant_label,
                ph_label,
            )
            return

        entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
        stages = entry.setdefault("stages", {})
        docking_stage = stages.get("docking")
        if not isinstance(docking_stage, MutableMapping):
            docking_stage = _default_stage_entry()
            stages["docking"] = docking_stage

        timing = _ensure_stage_timing(docking_stage)

        if event == "start":
            docking_stage["status"] = "running"
            docking_stage["error"] = None
            if not timing.get("started_at"):
                timing["started_at"] = now
        elif event in ("end", "fail"):
            if not timing.get("started_at"):
                timing["started_at"] = now
            timing["finished_at"] = now
            if elapsed_sec is not None:
                timing["wall_time_sec"] = round(float(elapsed_sec), 3)
            docking_stage["status"] = "failed" if event == "fail" else "completed"
            if error is not None:
                docking_stage["error"] = str(error)
        else:
            logging.debug(
                "[run-manifest.docking-overall.skip] unknown_event=%s run_id=%s pdb=%s variant=%s ph=%s",
                event,
                run_id,
                pdb_id,
                variant_label,
                ph_label,
            )
            return

        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest.docking-overall.error] run_id=%s pdb=%s variant=%s ph=%s event=%s",
            run_id,
            pdb_id,
            variant_label,
            ph_label,
            event,
            exc_info=True,
        )


def update_manifest_for_docking_stage(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: Optional[str],
    stage_name: str,
    status: str,
    *,
    ph_tag: Optional[str] = None,
    elapsed_sec: Optional[float] = None,
    error: Optional[str] = None,
) -> None:
    """
    Record status/timing for a single docking stage (e.g., "stage1", "stage2").
    """
    ph_label = ph_tag if ph_tag is not None else "base"
    stage_token = (stage_name or "stage").strip() or "stage"

    # --- normalize stage token into subrun + engine + base stage name ---
    lib_prefix = None
    core = stage_token
    if core.startswith("dud_"):
        lib_prefix = "dud"
        core = core[len("dud_") :] or core
    elif core.startswith("hmdb_"):
        lib_prefix = "hmdb"
        core = core[len("hmdb_") :] or core

    engine = None
    if core.startswith("vina_"):
        engine = "vina"
        core = core[len("vina_") :] or core
    elif core.startswith("gnina_"):
        engine = "gnina"
        core = core[len("gnina_") :] or core
    elif core.startswith("ledock_"):
        engine = "ledock"
        core = core[len("ledock_") :] or core

    # Handle engine-prefixed strings that still carry a library prefix (e.g., gnina_dud_stage1)
    if lib_prefix is None:
        if core.startswith("dud_"):
            lib_prefix = "dud"
            core = core[len("dud_") :] or core
        elif core.startswith("hmdb_"):
            lib_prefix = "hmdb"
            core = core[len("hmdb_") :] or core

    base_stage_name = core or stage_token
    engine = engine or "vina"
    subrun_label = lib_prefix if lib_prefix else "primary"

    if lib_prefix:
        raw_name = f"{lib_prefix}_{engine}_{base_stage_name}"
    else:
        raw_name = f"{engine}_{base_stage_name}"
    try:
        if not run_id:
            logging.debug(
                "[run-manifest.docking-stage.skip] no run_id pdb=%s variant=%s ph=%s stage=%s",
                pdb_id,
                variant_label,
                ph_label,
                stage_token,
            )
            return

        now = _utc_now_iso()
        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            logging.warning(
                "[run-manifest.docking-stage.skip] manifest missing run_id=%s path=%s pdb=%s variant=%s ph=%s stage=%s",
                run_id,
                manifest_path,
                pdb_id,
                variant_label,
                ph_label,
                stage_token,
            )
            return

        entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
        stages = entry.setdefault("stages", {})
        docking_stage = stages.get("docking")
        if not isinstance(docking_stage, MutableMapping):
            docking_stage = _default_stage_entry()
            stages["docking"] = docking_stage

        details = docking_stage.get("details")
        if not isinstance(details, MutableMapping):
            details = {}
        docking_stage["details"] = details

        per_stage = details.get("per_stage")
        if not isinstance(per_stage, MutableMapping):
            per_stage = {}
        details["per_stage"] = per_stage

        stage_entry = per_stage.get(raw_name)
        if not isinstance(stage_entry, MutableMapping):
            stage_entry = {
                "status": "pending",
                "error": None,
                "timing": {
                    "started_at": None,
                    "finished_at": None,
                    "wall_time_sec": None,
                },
                "subrun": subrun_label,
                "stage_base_name": base_stage_name,
            }
        per_stage[raw_name] = stage_entry
        stage_entry["subrun"] = subrun_label
        stage_entry["stage_base_name"] = base_stage_name

        timing = stage_entry.get("timing")
        if not isinstance(timing, MutableMapping):
            timing = {"started_at": None, "finished_at": None, "wall_time_sec": None}
        else:
            timing.setdefault("started_at", None)
            timing.setdefault("finished_at", None)
            timing.setdefault("wall_time_sec", None)
        stage_entry["timing"] = timing

        by_subrun = details.get("by_subrun")
        if not isinstance(by_subrun, MutableMapping):
            by_subrun = {}
        details["by_subrun"] = by_subrun

        subrun_map = by_subrun.get(subrun_label)
        if not isinstance(subrun_map, MutableMapping):
            subrun_map = {}
        by_subrun[subrun_label] = subrun_map
        subrun_map[raw_name] = stage_entry

        if status == "running":
            if not timing.get("started_at"):
                timing["started_at"] = now
            stage_entry["status"] = "running"
            stage_entry["error"] = None
        elif status in ("completed", "failed"):
            if not timing.get("started_at"):
                timing["started_at"] = now
            timing["finished_at"] = now
            if elapsed_sec is not None:
                timing["wall_time_sec"] = round(float(elapsed_sec), 3)
            stage_entry["status"] = status
            if error is not None:
                stage_entry["error"] = str(error)
        else:
            logging.debug(
                "[run-manifest.docking-stage.skip] unknown_status=%s run_id=%s pdb=%s variant=%s ph=%s stage=%s",
                status,
                run_id,
                pdb_id,
                variant_label,
                ph_label,
                stage_token,
            )
            return

        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest.docking-stage.error] run_id=%s pdb=%s variant=%s ph=%s stage=%s status=%s",
            run_id,
            pdb_id,
            variant_label,
            ph_label,
            stage_token,
            status,
            exc_info=True,
        )


def finalize_run_manifest(
    cfg: Mapping[str, Any],
    run_id: str,
    start_time: float,
    failed_entries: list[tuple[str, str, str, str, str]],
) -> None:
    try:
        _, manifest_path = get_manifest_paths(cfg, run_id)
        manifest = _load_manifest(manifest_path)
        if manifest is None:
            return

        manifest_status = "completed" if len(failed_entries) == 0 else "failed"
        manifest["status"] = manifest_status

        timing = manifest.setdefault("timing", {})
        if isinstance(timing, MutableMapping):
            timing.setdefault("created_at", _utc_now_iso())
            timing.setdefault("started_at", _utc_now_iso())
            timing["finished_at"] = _utc_now_iso()
            timing["wall_time_sec"] = round(time.time() - float(start_time), 3)
        else:
            manifest["timing"] = {
                "created_at": _utc_now_iso(),
                "started_at": _utc_now_iso(),
                "finished_at": _utc_now_iso(),
                "wall_time_sec": round(time.time() - float(start_time), 3),
            }

        _refresh_summary(manifest)
        _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to finalize manifest run_id=%s",
            run_id,
            exc_info=True,
        )
