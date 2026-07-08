#!/usr/bin/env python3
# ruff: noqa: E402
from __future__ import annotations

import argparse
import csv
import datetime
import json
import logging
import os
import sys
import time

import yaml  # type: ignore[import-untyped]
from pathlib import Path
from typing import Any

# Ensure repo-root and src imports resolve when launched as `python tools/...`.
REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
for _path in (REPO_ROOT, SRC_ROOT):
    _path_str = str(_path)
    if _path_str not in sys.path:
        sys.path.insert(0, _path_str)

import sitecustomize  # noqa: F401

from config.runtime_config import load_inputs, validate_config
from cli.postrun_hooks_runtime import (
    _log_rescore_verification,
    _maybe_run_dud_eval,
    _maybe_run_master_schema_export,
    _maybe_run_report_generation,
    _maybe_run_throughput_integrity,
    _write_run_efficiency_report,
)
from cli.run_manifest_runtime import (
    finalize_run_manifest,
    load_run_manifest,
    reconcile_distributed_manifest_state,
)
from cli.qol.manifest_repair import repair_run_manifest
from config.output_paths import output_root, runtime_root
from cli.run_context import ConfigDict, _apply_resume_config_from_snapshot


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Finalize a distributed run after array workers complete. "
            "Merges distributed manifest shards, finalizes run_manifest.yaml, "
            "and runs run-level post hooks."
        )
    )
    parser.add_argument("--run-id", required=True, help="Atlas run id to finalize.")
    parser.add_argument(
        "--no-hooks",
        action="store_true",
        help="Finalize manifest only; skip run-level post hooks.",
    )
    parser.add_argument(
        "--reconcile-only",
        action="store_true",
        help=(
            "Reconcile distributed per-protein state into run_manifest.yaml and exit "
            "(skip manifest status finalization and post hooks)."
        ),
    )
    return parser.parse_args()


def _load_effective_cfg(run_id: str) -> ConfigDict:
    os.environ["ATLAS_RUN_ID"] = str(run_id)
    base_cfg = ConfigDict(load_inputs())
    snap_cfg = _apply_resume_config_from_snapshot(base_cfg, run_id)
    cfg = ConfigDict(snap_cfg)
    cfg["RUN_ID"] = str(run_id)
    validate_config(cfg)
    return cfg


def _marker_dir(cfg: dict[str, Any], run_id: str) -> Path:
    return _manifest_run_dir(cfg, run_id) / "distributed" / "markers"


def _distributed_proteins_dir(cfg: dict[str, Any], run_id: str) -> Path:
    return _manifest_run_dir(cfg, run_id) / "distributed" / "proteins"


def _manifest_run_dir(cfg: dict[str, Any], run_id: str) -> Path:
    return runtime_root(cfg, "MANIFESTS_DIR", "manifests") / str(run_id)


def _to_int(value: Any) -> int | None:
    try:
        if value is None:
            return None
        return int(str(value).strip())
    except Exception:
        return None


def _collect_marker_stats(cfg: dict[str, Any], run_id: str) -> dict[str, Any]:
    marker_root = _marker_dir(cfg, run_id)
    payloads: list[dict[str, Any]] = []
    present_task_ids: set[int] = set()
    task_count_hints: set[int] = set()
    marker_files = 0

    if marker_root.exists():
        for marker_path in sorted(marker_root.glob("task_*.json")):
            if not marker_path.is_file():
                continue
            marker_files += 1
            try:
                payload = json.loads(marker_path.read_text(encoding="utf-8")) or {}
            except Exception:
                logging.warning(
                    "[dist-finalize.marker.read] action=skip path=%s",
                    marker_path,
                    exc_info=True,
                )
                continue
            if not isinstance(payload, dict):
                continue
            payloads.append(payload)
            tid = _to_int(payload.get("task_id"))
            if tid is not None and tid >= 0:
                present_task_ids.add(tid)
            hinted = _to_int(payload.get("task_count"))
            if hinted is not None and hinted > 0:
                task_count_hints.add(hinted)

    expected_task_count: int | None = None
    if task_count_hints:
        expected_task_count = max(task_count_hints)

    missing_task_ids: list[int] = []
    if expected_task_count and present_task_ids:
        base = min(present_task_ids)
        expected_ids = set(range(base, base + int(expected_task_count)))
        missing_task_ids = sorted(expected_ids - present_task_ids)

    return {
        "marker_files": marker_files,
        "marker_payloads": len(payloads),
        "present_task_ids": sorted(present_task_ids),
        "expected_task_count": expected_task_count,
        "missing_task_ids": missing_task_ids,
    }


def _collect_failed_entries_from_markers(
    cfg: dict[str, Any], run_id: str
) -> list[tuple[str, str, str, str, str]]:
    out: list[tuple[str, str, str, str, str]] = []
    marker_root = _marker_dir(cfg, run_id)
    if not marker_root.exists():
        return out

    for marker_path in sorted(marker_root.glob("task_*.json")):
        if not marker_path.is_file():
            continue
        try:
            payload = json.loads(marker_path.read_text(encoding="utf-8")) or {}
        except Exception:
            logging.warning(
                "[dist-finalize.marker.read] action=skip path=%s", marker_path, exc_info=True
            )
            continue
        entries = payload.get("failed_entries") or []
        if not isinstance(entries, list):
            continue
        for raw in entries:
            if not isinstance(raw, (list, tuple)) or len(raw) < 5:
                continue
            out.append(
                (
                    str(raw[0]),
                    str(raw[1]),
                    str(raw[2]),
                    str(raw[3]),
                    str(raw[4]),
                )
            )
    return out


def _completed_manifest_entry_keys(
    cfg: dict[str, Any], run_id: str
) -> set[tuple[str, str, str]]:
    manifest = load_run_manifest(cfg, run_id) or {}
    proteins = manifest.get("proteins") if isinstance(manifest, dict) else None
    if not isinstance(proteins, dict):
        return set()

    completed: set[tuple[str, str, str]] = set()
    for raw_key, entry in proteins.items():
        if not isinstance(entry, dict):
            continue
        status = str(entry.get("status") or "").strip().lower()
        if status != "completed":
            continue
        key_bits = str(raw_key).split("|")
        pdb_id = str(entry.get("pdb_id") or (key_bits[0] if key_bits else "")).strip()
        variant = str(
            entry.get("variant") or (key_bits[1] if len(key_bits) > 1 else "LEGACY")
        ).strip()
        ph = str(entry.get("ph") or (key_bits[2] if len(key_bits) > 2 else "base")).strip()
        if not pdb_id:
            continue
        completed.add(
            (
                pdb_id.upper(),
                _normalize_variant_token(variant),
                _normalize_ph_token(ph),
            )
        )
    return completed


def _drop_resolved_failure_entries(
    cfg: dict[str, Any],
    run_id: str,
    entries: list[tuple[str, str, str, str, str]],
) -> tuple[list[tuple[str, str, str, str, str]], int]:
    if not entries:
        return entries, 0
    completed = _completed_manifest_entry_keys(cfg, run_id)
    if not completed:
        return entries, 0

    unresolved: list[tuple[str, str, str, str, str]] = []
    resolved = 0
    for entry in entries:
        key = (
            str(entry[0]).strip().upper(),
            _normalize_variant_token(entry[1]),
            _normalize_ph_token(entry[2]),
        )
        if key in completed:
            resolved += 1
            continue
        unresolved.append(entry)
    return unresolved, resolved


def _collect_scorch_queue_marker_failures(
    cfg: dict[str, Any],
    run_id: str,
) -> list[tuple[str, str, int]]:
    out: list[tuple[str, str, int]] = []
    phase_root = _manifest_run_dir(cfg, run_id) / "distributed" / "phase_markers"
    if not phase_root.exists():
        return out

    fields = ("scorch_queue_failed", "scorch_reconcile_failed")
    for phase in ("scorch_queue_closed", "scorch_drain"):
        marker_dir = phase_root / phase
        if not marker_dir.exists():
            continue
        for marker_path in sorted(marker_dir.glob("task_*.json")):
            if not marker_path.is_file():
                continue
            try:
                payload = json.loads(marker_path.read_text(encoding="utf-8")) or {}
            except Exception:
                logging.warning(
                    "[dist-finalize.scorch-markers] action=skip path=%s",
                    marker_path,
                    exc_info=True,
                )
                continue
            if not isinstance(payload, dict):
                continue
            for field in fields:
                value = _to_int(payload.get(field))
                if value is not None and value > 0:
                    out.append((phase, marker_path.name, int(value)))
    return out


def _collect_incomplete_entries_from_manifest(
    cfg: dict[str, Any], run_id: str
) -> list[tuple[str, str, str, str, str]]:
    manifest = load_run_manifest(cfg, run_id) or {}
    proteins = manifest.get("proteins") if isinstance(manifest, dict) else None
    if not isinstance(proteins, dict):
        return []

    out: list[tuple[str, str, str, str, str]] = []
    terminal = {"completed", "failed", "skipped"}
    chunk_coverage = _load_distributed_chunk_coverage(cfg, run_id)
    for raw_key, entry in sorted(proteins.items()):
        if not isinstance(entry, dict):
            continue
        status = str(entry.get("status") or "").strip().lower()
        if status in {"failed", "skipped"}:
            continue

        key_bits = str(raw_key).split("|")
        pdb_id = str(entry.get("pdb_id") or (key_bits[0] if key_bits else "")).strip()
        variant = str(
            entry.get("variant") or (key_bits[1] if len(key_bits) > 1 else "LEGACY")
        ).strip()
        ph = str(entry.get("ph") or (key_bits[2] if len(key_bits) > 2 else "base")).strip()
        if not pdb_id:
            continue
        coverage_key = (
            str(pdb_id).strip().upper(),
            _normalize_variant_token(variant),
            _normalize_ph_token(ph),
        )
        complete, coverage_reason = _consensus_output_complete(
            cfg,
            run_id,
            pdb_id,
            variant,
            ph,
            chunk_coverage=chunk_coverage,
        )
        chunk_plans_present = bool(chunk_coverage)
        has_coverage = coverage_key in chunk_coverage
        if complete and status in terminal and (not chunk_plans_present or has_coverage):
            continue
        if complete and not has_coverage and chunk_plans_present:
            coverage_reason = "missing_chunk_plan_coverage"
        elif complete and status not in terminal and not has_coverage:
            coverage_reason = "coverage_unresolved_nonterminal"
        if status in terminal:
            reason = coverage_reason or "incomplete_consensus_coverage"
        else:
            reason = f"nonterminal_status={status or 'missing'}"
            if coverage_reason:
                reason = f"{reason};{coverage_reason}"
        out.append((pdb_id, variant or "LEGACY", ph or "base", "DistributedIncomplete", reason))
    return out


def _normalize_variant_token(raw: Any) -> str:
    token = str(raw or "").strip().upper()
    if token in {"", "BASE", "NONE", "NULL"}:
        return "LEGACY"
    return token


def _normalize_ph_token(raw: Any) -> str:
    token = str(raw or "").strip()
    if token.lower() in {"", "base", "none", "null"}:
        return "base"
    return token


def _distributed_manifest_dirs(cfg: dict[str, Any], run_id: str) -> list[Path]:
    candidates = [
        runtime_root(cfg, "MANIFESTS_DIR", "manifests") / str(run_id) / "distributed"
    ]
    overall_dir = Path(str(cfg.get("OVERALL_DIR", ".") or "."))
    candidates.append(output_root(overall_dir, "manifests") / str(run_id) / "distributed")
    out: list[Path] = []
    seen: set[str] = set()
    for path in candidates:
        token = str(path)
        if token in seen:
            continue
        seen.add(token)
        out.append(path)
    return out


def _load_distributed_chunk_coverage(
    cfg: dict[str, Any],
    run_id: str,
) -> dict[tuple[str, str, str], dict[str, set[str]]]:
    coverage: dict[tuple[str, str, str], dict[str, set[str]]] = {}
    for dist_dir in _distributed_manifest_dirs(cfg, run_id):
        chunk_results: dict[str, str] = {}
        for result_path in sorted((dist_dir / "chunk_results").glob("*.json")):
            if not result_path.is_file():
                continue
            try:
                payload = json.loads(result_path.read_text(encoding="utf-8")) or {}
            except Exception:
                continue
            if not isinstance(payload, dict):
                continue
            chunk_id = str(payload.get("chunk_id") or result_path.stem).strip()
            status = str(payload.get("status") or "").strip().lower()
            if chunk_id:
                chunk_results[chunk_id] = status

        for plan_path in sorted(dist_dir.glob("combo_chunks_*.json")):
            if not plan_path.is_file():
                continue
            try:
                payload = json.loads(plan_path.read_text(encoding="utf-8")) or {}
            except Exception:
                continue
            chunks = payload.get("chunks") if isinstance(payload, dict) else payload
            if not isinstance(chunks, list):
                continue
            for chunk in chunks:
                if not isinstance(chunk, dict):
                    continue
                pdb_id = str(chunk.get("pdb_id") or "").strip().upper()
                if not pdb_id:
                    pdb_file = str(chunk.get("pdb_file") or "").strip()
                    if pdb_file:
                        pdb_id = Path(pdb_file).stem.strip().upper()
                if not pdb_id:
                    continue
                key = (
                    pdb_id,
                    _normalize_variant_token(chunk.get("variant_label")),
                    _normalize_ph_token(chunk.get("ph_tag")),
                )
                slot = coverage.setdefault(
                    key,
                    {
                        "chunks": set(),
                        "completed_chunks": set(),
                        "failed_chunks": set(),
                        "expected_ligands": set(),
                    },
                )
                chunk_id = str(chunk.get("chunk_id") or "").strip()
                if chunk_id:
                    slot["chunks"].add(chunk_id)
                    chunk_status = chunk_results.get(chunk_id)
                    if chunk_status == "completed":
                        slot["completed_chunks"].add(chunk_id)
                    elif chunk_status in {"failed", "terminal_failed"}:
                        slot["failed_chunks"].add(chunk_id)
                ligands_raw = chunk.get("ligand_bases") or []
                if isinstance(ligands_raw, list):
                    slot["expected_ligands"].update(
                        str(item).strip() for item in ligands_raw if str(item).strip()
                    )
    return coverage


def _consensus_output_candidates(
    cfg: dict[str, Any],
    run_id: str,
    pdb_id: str,
    variant: str,
    ph: str,
) -> list[Path]:
    docked_root = runtime_root(cfg, "DOCKED_DIR", "outputs/docked") / str(run_id)
    pdb = str(pdb_id).strip()
    variant_token = str(variant or "").strip()
    ph_token = str(ph or "").strip()
    ph_is_base = ph_token.lower() in {"", "base", "none", "null"}
    variant_is_flat = variant_token.lower() in {"", "legacy", "base", "none", "null"}

    candidates = [docked_root / pdb / "consensus_docking_scores.csv"]
    if variant_token and not variant_is_flat:
        candidates.append(docked_root / pdb / variant_token / "consensus_docking_scores.csv")
        if ph_token and not ph_is_base:
            candidates.append(
                docked_root
                / pdb
                / variant_token
                / ph_token
                / "consensus_docking_scores.csv"
            )
    if ph_token and not ph_is_base:
        candidates.append(docked_root / pdb / ph_token / "consensus_docking_scores.csv")
    return candidates


def _count_csv_rows(path: Path) -> int:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            return sum(1 for _ in csv.DictReader(handle))
    except Exception:
        return 0


def _consensus_output_complete(
    cfg: dict[str, Any],
    run_id: str,
    pdb_id: str,
    variant: str,
    ph: str,
    *,
    chunk_coverage: dict[tuple[str, str, str], dict[str, set[str]]] | None = None,
) -> tuple[bool, str]:
    found_path: Path | None = None
    for path in _consensus_output_candidates(cfg, run_id, pdb_id, variant, ph):
        try:
            if path.is_file() and path.stat().st_size > 0:
                found_path = path
                break
        except Exception:
            continue
    if found_path is None:
        return False, "missing_consensus_csv"

    key = (
        str(pdb_id).strip().upper(),
        _normalize_variant_token(variant),
        _normalize_ph_token(ph),
    )
    coverage = (chunk_coverage or {}).get(key)
    if not coverage:
        return True, ""

    chunk_ids = coverage.get("chunks") or set()
    completed = coverage.get("completed_chunks") or set()
    failed = coverage.get("failed_chunks") or set()
    if chunk_ids:
        if failed:
            return False, f"chunk_plan_failed={len(failed)}/{len(chunk_ids)}"
        if len(completed) < len(chunk_ids):
            return False, f"chunk_plan_incomplete={len(completed)}/{len(chunk_ids)}"

    expected_ligands = coverage.get("expected_ligands") or set()
    if expected_ligands:
        row_count = _count_csv_rows(found_path)
        if row_count < len(expected_ligands):
            return False, f"consensus_rows={row_count}/{len(expected_ligands)}"
    return True, ""


def _consensus_output_exists(
    cfg: dict[str, Any],
    run_id: str,
    pdb_id: str,
    variant: str,
    ph: str,
) -> bool:
    for path in _consensus_output_candidates(cfg, run_id, pdb_id, variant, ph):
        try:
            if path.is_file() and path.stat().st_size > 0:
                return True
        except Exception:
            continue
    return False


def _resolve_manifest_started_epoch(cfg: dict[str, Any], run_id: str) -> float:
    manifest = load_run_manifest(cfg, run_id) or {}
    timing = manifest.get("timing") if isinstance(manifest, dict) else None
    started_at = timing.get("started_at") if isinstance(timing, dict) else None
    if not started_at:
        return time.time()
    try:
        started_dt = datetime.datetime.fromisoformat(
            str(started_at).replace("Z", "+00:00")
        )
        return float(started_dt.timestamp())
    except Exception:
        return time.time()


def _scheduled_pdb_ids(cfg: dict[str, Any], run_id: str) -> list[str]:
    manifest = load_run_manifest(cfg, run_id) or {}
    summary = manifest.get("summary") if isinstance(manifest, dict) else None
    listed = summary.get("total_protein_list") if isinstance(summary, dict) else None
    if not isinstance(listed, list):
        return []
    return [str(x) for x in listed if str(x).strip()]


def _requires_strict_throughput(cfg: dict[str, Any], run_id: str) -> bool:
    manifest = load_run_manifest(cfg, run_id) or {}
    command = manifest.get("command") if isinstance(manifest, dict) else None
    argv = str(command.get("argv") or "") if isinstance(command, dict) else ""
    argv_l = argv.lower()
    run_id_l = str(run_id or "").lower()
    if "-bench-small" in argv_l or "-bench2" in argv_l or "-bench" in argv_l:
        return True
    return "bench" in run_id_l


def _truthy_cfg(value: Any) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y", "on"}


def _scorch_repair_summary_incomplete(require_scorch: bool, summary: Any) -> bool:
    if not require_scorch:
        return False
    return bool(
        int(getattr(summary, "targets_running", 0) or 0) > 0
        or int(getattr(summary, "targets_failed", 0) or 0) > 0
        or int(getattr(summary, "targets_completed", 0) or 0)
        < int(getattr(summary, "targets", 0) or 0)
    )


def _repair_summary_complete(require_scorch: bool, summary: Any) -> bool:
    targets = int(getattr(summary, "targets", 0) or 0)
    planned_chunks = int(getattr(summary, "planned_chunks", 0) or 0)
    if targets <= 0 or planned_chunks <= 0:
        return False
    if int(getattr(summary, "targets_completed", 0) or 0) < targets:
        return False
    if int(getattr(summary, "targets_running", 0) or 0) > 0:
        return False
    if int(getattr(summary, "targets_failed", 0) or 0) > 0:
        return False
    if int(getattr(summary, "completed_chunks", 0) or 0) < planned_chunks:
        return False
    if int(getattr(summary, "missing_chunks", 0) or 0) > 0:
        return False
    if int(getattr(summary, "failed_chunks", 0) or 0) > 0:
        return False
    if require_scorch and int(getattr(summary, "scorch_reranked_files", 0) or 0) < targets:
        return False
    return True


def _all_dirs_root_from_cfg(cfg: dict[str, Any]) -> Path:
    manifests_root = runtime_root(cfg, "MANIFESTS_DIR", "manifests")
    if manifests_root.name == "manifests" and manifests_root.parent.name == "outputs":
        return manifests_root.parent.parent
    return runtime_root(cfg, "OVERALL_DIR", ".")


def _chunk_library_overrides(
    cfg: dict[str, Any],
    run_id: str,
) -> dict[tuple[str, str, str], str]:
    dist_dir = _manifest_run_dir(cfg, run_id) / "distributed"
    if not dist_dir.exists():
        return {}

    by_key: dict[tuple[str, str, str], set[str]] = {}
    for plan_path in sorted(dist_dir.glob("combo_chunks_*.json")):
        if not plan_path.is_file():
            continue
        try:
            payload = json.loads(plan_path.read_text(encoding="utf-8")) or {}
        except Exception:
            continue
        chunks = payload.get("chunks") if isinstance(payload, dict) else payload
        if not isinstance(chunks, list):
            continue
        for chunk in chunks:
            if not isinstance(chunk, dict):
                continue
            library_name = str(chunk.get("library_name") or "").strip()
            if not library_name:
                continue
            pdb_id = str(chunk.get("pdb_id") or "").strip().upper()
            if not pdb_id:
                pdb_file = str(chunk.get("pdb_file") or "").strip()
                if pdb_file:
                    pdb_id = Path(pdb_file).stem.strip().upper()
            variant = str(chunk.get("variant_label") or "").strip().upper()
            ph = str(chunk.get("ph_tag") or "").strip()
            if not (pdb_id and variant and ph):
                continue
            by_key.setdefault((pdb_id, variant, ph), set()).add(library_name)

    out: dict[tuple[str, str, str], str] = {}
    for key, choices in by_key.items():
        if len(choices) == 1:
            out[key] = next(iter(choices))
    return out


def _backfill_manifest_libraries_from_chunks(
    cfg: dict[str, Any],
    run_id: str,
) -> int:
    overrides = _chunk_library_overrides(cfg, run_id)
    if not overrides:
        return 0

    manifest_path = _manifest_run_dir(cfg, run_id) / "run_manifest.yaml"
    if not manifest_path.exists():
        return 0

    try:
        manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
    except Exception:
        return 0
    if not isinstance(manifest, dict):
        return 0

    proteins = manifest.get("proteins")
    if not isinstance(proteins, dict):
        return 0

    changed = 0
    for raw_key, entry in proteins.items():
        if not isinstance(entry, dict):
            continue

        pdb_id = str(entry.get("pdb_id") or "").strip().upper()
        variant = str(entry.get("variant") or "").strip().upper()
        ph = str(entry.get("ph") or "").strip()
        if not (pdb_id and variant and ph):
            bits = str(raw_key).split("|")
            if len(bits) >= 3:
                pdb_id = pdb_id or bits[0].strip().upper()
                variant = variant or bits[1].strip().upper()
                ph = ph or bits[2].strip()
        if not (pdb_id and variant and ph):
            continue

        library_name = overrides.get((pdb_id, variant, ph))
        if not library_name:
            continue
        current_library = str(entry.get("library") or "").strip()
        if current_library == library_name:
            continue

        entry["library"] = library_name
        changed += 1

    if changed > 0:
        manifest_path.write_text(
            yaml.safe_dump(manifest, sort_keys=False),
            encoding="utf-8",
        )
    return changed


def main() -> int:
    args = _parse_args()
    run_id = str(args.run_id).strip()
    if not run_id:
        print("ERROR: --run-id is required", file=sys.stderr)
        return 2

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s - %(message)s",
    )

    cfg = _load_effective_cfg(run_id)
    marker_stats = _collect_marker_stats(cfg, run_id)
    protein_state_count = len(
        [
            p
            for p in _distributed_proteins_dir(cfg, run_id).glob("*.json")
            if p.is_file()
        ]
    )
    logging.info(
        "[dist-finalize.inputs] run_id=%s marker_files=%d marker_payloads=%d state_files=%d expected_task_count=%s present_tasks=%s",
        run_id,
        int(marker_stats.get("marker_files", 0)),
        int(marker_stats.get("marker_payloads", 0)),
        protein_state_count,
        marker_stats.get("expected_task_count"),
        marker_stats.get("present_task_ids"),
    )
    missing_ids = marker_stats.get("missing_task_ids") or []
    if missing_ids:
        logging.warning(
            "[dist-finalize.inputs] run_id=%s marker_gap missing_task_ids=%s",
            run_id,
            missing_ids,
        )

    reconcile_stats = reconcile_distributed_manifest_state(cfg, run_id, force=True)
    logging.info(
        "[dist-finalize.reconcile] run_id=%s action=%s state_files=%s merged_entries=%s stale=%s",
        run_id,
        reconcile_stats.get("action"),
        reconcile_stats.get("state_files"),
        reconcile_stats.get("merged_entries"),
        reconcile_stats.get("manifest_stale"),
    )
    library_backfilled = _backfill_manifest_libraries_from_chunks(cfg, run_id)
    if library_backfilled > 0:
        logging.info(
            "[dist-finalize.library-backfill] run_id=%s entries=%d",
            run_id,
            int(library_backfilled),
        )
    if args.reconcile_only:
        logging.info(
            "[dist-finalize.reconcile] run_id=%s action=done (--reconcile-only)",
            run_id,
        )
        return 0

    require_scorch = _truthy_cfg(cfg.get("USE_SCORCH"))
    scorch_repair_incomplete = False
    repair_complete = False
    try:
        repair_summary = repair_run_manifest(
            run_id,
            root=_all_dirs_root_from_cfg(cfg),
            dry_run=False,
            require_scorch=require_scorch,
            clear_incomplete_scorch_outputs=False,
        )
        logging.info(
            "[dist-finalize.manifest-repair] run_id=%s require_scorch=%s targets=%d completed=%d running=%d failed=%d chunks=%d/%d scorch=%d",
            run_id,
            str(require_scorch).lower(),
            int(repair_summary.targets),
            int(repair_summary.targets_completed),
            int(repair_summary.targets_running),
            int(repair_summary.targets_failed),
            int(repair_summary.completed_chunks),
            int(repair_summary.planned_chunks),
            int(repair_summary.scorch_reranked_files),
        )
        manifest_after = load_run_manifest(cfg, run_id) or {}
        manifest_status = (
            manifest_after.get("status", "unknown")
            if isinstance(manifest_after, dict)
            else "unknown"
        )
        logging.info(
            "[dist-finalize.manifest] run_id=%s status=%s source=repair",
            run_id,
            manifest_status,
        )
        scorch_repair_incomplete = _scorch_repair_summary_incomplete(
            require_scorch,
            repair_summary,
        )
        repair_complete = _repair_summary_complete(require_scorch, repair_summary)
    except Exception:
        logging.warning(
            "[dist-finalize.manifest-repair] run_id=%s action=skip reason=unexpected_exception",
            run_id,
            exc_info=True,
        )
    if require_scorch and scorch_repair_incomplete:
        logging.error(
            "[dist-finalize.manifest-repair] run_id=%s require_scorch=true status=incomplete",
            run_id,
        )
        return 3

    failed_entries = _collect_failed_entries_from_markers(cfg, run_id)
    failed_entries, resolved_failures = _drop_resolved_failure_entries(
        cfg,
        run_id,
        failed_entries,
    )
    if resolved_failures:
        logging.info(
            "[dist-finalize.failures] run_id=%s resolved_by_manifest_repair=%d",
            run_id,
            int(resolved_failures),
        )
    if repair_complete:
        incomplete_entries = []
        logging.info(
            "[dist-finalize.incomplete] run_id=%s action=skip reason=repair_summary_complete",
            run_id,
        )
    else:
        incomplete_entries = _collect_incomplete_entries_from_manifest(cfg, run_id)
    if incomplete_entries:
        seen = set(failed_entries)
        for entry in incomplete_entries:
            if entry not in seen:
                failed_entries.append(entry)
                seen.add(entry)
        logging.warning(
            "[dist-finalize.incomplete] run_id=%s nonterminal_entries=%d examples=%s",
            run_id,
            len(incomplete_entries),
            ",".join("|".join(x[:3]) for x in incomplete_entries[:8]),
        )
    logging.info(
        "[dist-finalize.failures] run_id=%s failed_entries=%d",
        run_id,
        len(failed_entries),
    )

    started_epoch = _resolve_manifest_started_epoch(cfg, run_id)
    finalize_run_manifest(
        cfg,
        run_id,
        start_time=started_epoch,
        failed_entries=failed_entries,
    )
    manifest_after = load_run_manifest(cfg, run_id) or {}
    manifest_status = (
        manifest_after.get("status", "unknown")
        if isinstance(manifest_after, dict)
        else "unknown"
    )
    logging.info("[dist-finalize.manifest] run_id=%s status=%s", run_id, manifest_status)

    scorch_marker_failures = _collect_scorch_queue_marker_failures(cfg, run_id)
    if scorch_marker_failures:
        if repair_complete:
            logging.info(
                "[dist-finalize.scorch-markers] run_id=%s action=ignore reason=repair_summary_complete stale_failures=%d examples=%s",
                run_id,
                len(scorch_marker_failures),
                ",".join(
                    f"{phase}/{name}:{count}"
                    for phase, name, count in scorch_marker_failures[:8]
                ),
            )
        else:
            logging.error(
                "[dist-finalize.scorch-markers] run_id=%s status=failed failures=%d examples=%s",
                run_id,
                len(scorch_marker_failures),
                ",".join(
                    f"{phase}/{name}:{count}"
                    for phase, name, count in scorch_marker_failures[:8]
                ),
            )
            return 4

    if args.no_hooks:
        logging.info("[dist-finalize.hooks] run_id=%s action=skip (--no-hooks)", run_id)
        return 0

    pdb_ids = _scheduled_pdb_ids(cfg, run_id)
    try:
        _maybe_run_dud_eval(cfg, run_id, pdb_ids)
    except Exception:
        logging.warning(
            "[dist-finalize.dud-eval] run_id=%s action=skip reason=unexpected_exception",
            run_id,
            exc_info=True,
        )
    try:
        _maybe_run_master_schema_export(cfg, run_id)
    except Exception:
        logging.warning(
            "[dist-finalize.master-export] run_id=%s action=skip reason=unexpected_exception",
            run_id,
            exc_info=True,
        )
    try:
        _maybe_run_report_generation(cfg, run_id)
    except Exception:
        logging.warning(
            "[dist-finalize.report] run_id=%s action=skip reason=unexpected_exception",
            run_id,
            exc_info=True,
        )
    try:
        _write_run_efficiency_report(cfg, run_id)
    except Exception:
        logging.warning(
            "[dist-finalize.run-efficiency] run_id=%s action=skip reason=unexpected_exception",
            run_id,
            exc_info=True,
        )
    strict_throughput = _requires_strict_throughput(cfg, run_id)
    throughput_ok = True
    try:
        throughput_ok = _maybe_run_throughput_integrity(
            cfg, run_id, strict=strict_throughput
        )
    except Exception:
        throughput_ok = False
        logging.warning(
            "[dist-finalize.throughput-integrity] run_id=%s action=skip reason=unexpected_exception",
            run_id,
            exc_info=True,
        )
    if strict_throughput and not throughput_ok:
        logging.error(
            "[dist-finalize.throughput-integrity] run_id=%s strict=true status=failed",
            run_id,
        )
        return 2
    try:
        _log_rescore_verification(run_id, cfg)
    except Exception:
        logging.warning(
            "[dist-finalize.post-check] run_id=%s action=skip reason=unexpected_exception",
            run_id,
            exc_info=True,
        )
    logging.info("[dist-finalize.hooks] run_id=%s status=done", run_id)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
