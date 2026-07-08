# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
)

import yaml  # type: ignore[import-untyped]
from config.output_paths import output_root, run_output_dir

_prod_canonical_ligand_base: Optional[Callable[[str], str]] = None
try:
    from src.post_docking.rescoring.rescore_reranker import (
        canonical_ligand_base as _prod_canonical_ligand_base_imported,
    )
    _prod_canonical_ligand_base = _prod_canonical_ligand_base_imported
except Exception:  # pragma: no cover - fallback only
    pass

COMPONENT = "[throughput-integrity]"

_TRUE_VALUES = {"1", "true", "yes", "on", "y", "t"}
_TS_PREFIX_RE = re.compile(r"^(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3})")
_FAILED_MARKER_SKEW_SEC = 60.0
_ALLOW_EMPTY_EXPECTED_RUN_IDS: Set[str] = set()
_COMPLETED_STATUSES = {"completed", "success", "succeeded"}


def _configure_logging(verbose: bool) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")
    return logging.getLogger("throughput-integrity")


def _as_bool(value: Any) -> bool:
    token = str(value or "").strip().lower()
    return token in _TRUE_VALUES


def _chunk_plan_checks_enabled() -> bool:
    return not _as_bool(os.environ.get("ATLAS_DISABLE_COMBO_CHUNKS", "0"))


def _canonical_ligand_base(name: str) -> str:
    s = Path(str(name or "")).name.strip()
    s = re.sub(r"\.failed\.txt$", "", s, flags=re.IGNORECASE)
    if _prod_canonical_ligand_base is not None:
        try:
            canon = str(_prod_canonical_ligand_base(s)).strip()
            if canon:
                return canon
        except Exception:
            pass
    s = re.sub(r"\.(pdbqt|mol2|sdf|pdb|txt|csv)$", "", s, flags=re.IGNORECASE)
    s = s.replace(".sanitized", "")
    s = re.sub(r"(_gnina_stage\d+)$", "", s)
    s = re.sub(r"(_dock6_stage\d+)$", "", s)
    s = re.sub(r"(_ledock_stage\d+)$", "", s)
    s = re.sub(r"(_stage\d+)$", "", s)
    s = re.sub(r"(__gnina_stage\d+)$", "", s)
    s = re.sub(r"(__dock6_stage\d+)$", "", s)
    s = re.sub(r"(__ledock_stage\d+)$", "", s)
    s = s.replace("__", "_")
    s = re.sub(r"_+$", "", s)
    return s


def _parse_combo_key(key: str) -> Optional[Tuple[str, str, str]]:
    parts = str(key).split("|")
    if len(parts) < 3:
        return None
    pdb_id = parts[0].strip().upper()
    variant = parts[1].strip().upper()
    ph = parts[2].strip()
    if not pdb_id or not variant or not ph:
        return None
    return pdb_id, variant, ph


def _read_manifest(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    try:
        payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _env_run_dir(env_key: str, fallback: Path, run_id: str) -> Path:
    explicit = os.environ.get(env_key)
    if explicit and str(explicit).strip():
        return Path(str(explicit)).expanduser() / str(run_id)
    return fallback


def _env_root(env_key: str, fallback: Path) -> Path:
    explicit = os.environ.get(env_key)
    if explicit and str(explicit).strip():
        return Path(str(explicit)).expanduser()
    return fallback


def _entry_completed(entry: Dict[str, Any]) -> bool:
    status = str(entry.get("status") or "").strip().lower()
    return status in _COMPLETED_STATUSES


def _read_library_manifest(library_dir: Path) -> Tuple[Set[str], Set[str]]:
    """
    Returns:
      expected: all ligand bases declared in manifest
      prepped:  declared ligand bases whose files exist
    """
    expected: Set[str] = set()
    prepped: Set[str] = set()
    manifest_path = library_dir / "_manifest.json"
    if manifest_path.exists():
        try:
            payload = json.loads(manifest_path.read_text(encoding="utf-8")) or {}
        except Exception:
            payload = {}
        entries = payload.get("entries")
        if isinstance(entries, dict):
            for key, rel in entries.items():
                base = _canonical_ligand_base(str(key))
                if not base:
                    continue
                expected.add(base)
                rel_path = library_dir / str(rel)
                if rel_path.exists():
                    prepped.add(base)

    if expected:
        return expected, prepped

    # Fallback for libraries without manifest metadata.
    for p in library_dir.rglob("*.pdbqt"):
        if not p.is_file():
            continue
        base = _canonical_ligand_base(p.name)
        if base:
            expected.add(base)
            prepped.add(base)
    return expected, prepped


def _load_stage1_configs(combo_config_dir: Path) -> Set[str]:
    out: Set[str] = set()
    stage1_dir = combo_config_dir / "stage1"
    if not stage1_dir.exists():
        return out
    for cfg in stage1_dir.glob("*_stage1.txt"):
        base = _canonical_ligand_base(cfg.stem)
        if base:
            out.add(base)
    return out


def _resolve_docking_summary_csv(combo_docked_dir: Path, library: str) -> Path:
    suffix = "_docking_score_summary.csv"
    default_path = combo_docked_dir / "docking_score_summary.csv"
    prefixed = sorted(
        p
        for p in combo_docked_dir.glob(f"*{suffix}")
        if p.name != "docking_score_summary.csv"
    )

    if not prefixed:
        return default_path

    lib = str(library or "").strip()
    if lib:
        wanted = {lib.lower(), _canonical_ligand_base(lib).lower()}
        for candidate in prefixed:
            prefix = candidate.name[: -len(suffix)].strip().lower()
            if prefix in wanted:
                return candidate

    if default_path.exists():
        return default_path
    return prefixed[0]


def _load_distributed_chunk_plans(
    manifest_path: Path,
) -> tuple[
    Dict[Tuple[str, str, str, str], Set[str]],
    Dict[Tuple[str, str, str], Set[str]],
    Dict[Tuple[str, str, str, str], Set[str]],
    Dict[Tuple[str, str, str], Set[str]],
]:
    by_combo4_expected: Dict[Tuple[str, str, str, str], Set[str]] = {}
    by_combo3_expected: Dict[Tuple[str, str, str], Set[str]] = {}
    by_combo4_chunks: Dict[Tuple[str, str, str, str], Set[str]] = {}
    by_combo3_chunks: Dict[Tuple[str, str, str], Set[str]] = {}

    dist_dir = manifest_path.parent / "distributed"
    if not dist_dir.exists():
        return by_combo4_expected, by_combo3_expected, by_combo4_chunks, by_combo3_chunks

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
            variant = str(chunk.get("variant_label") or "").strip().upper()
            ph = str(chunk.get("ph_tag") or "").strip() or "base"
            library = str(chunk.get("library_name") or "").strip()
            if not (pdb_id and variant and ph):
                continue
            ligands_raw = chunk.get("ligand_bases") or []
            if not isinstance(ligands_raw, list):
                ligands_raw = []
            ligands = {
                _canonical_ligand_base(str(x))
                for x in ligands_raw
                if _canonical_ligand_base(str(x))
            }
            chunk_id = str(chunk.get("chunk_id") or "").strip()

            key4 = (pdb_id, variant, ph, library)
            key3 = (pdb_id, variant, ph)
            if ligands:
                by_combo4_expected.setdefault(key4, set()).update(ligands)
                by_combo3_expected.setdefault(key3, set()).update(ligands)
            if chunk_id:
                by_combo4_chunks.setdefault(key4, set()).add(chunk_id)
                by_combo3_chunks.setdefault(key3, set()).add(chunk_id)

    return by_combo4_expected, by_combo3_expected, by_combo4_chunks, by_combo3_chunks


def _load_distributed_chunk_results(manifest_path: Path) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    result_dir = manifest_path.parent / "distributed" / "chunk_results"
    if not result_dir.exists():
        return out
    for path in sorted(result_dir.glob("*.json")):
        if not path.is_file():
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8")) or {}
        except Exception:
            continue
        if not isinstance(payload, dict):
            continue
        chunk_id = str(payload.get("chunk_id") or "").strip()
        if not chunk_id:
            chunk_id = path.stem
        if chunk_id:
            out[chunk_id] = payload
    return out


def _coverage_snapshot_path(
    manifest_path: Path, *, pdb_id: str, variant: str, ph: str
) -> Path:
    filename = f"{str(pdb_id).upper()}__{str(variant).upper()}__{str(ph)}.json"
    return manifest_path.parent / "coverage" / filename


def _load_combo_coverage_snapshot(
    manifest_path: Path, *, pdb_id: str, variant: str, ph: str
) -> Dict[str, Any]:
    path = _coverage_snapshot_path(
        manifest_path, pdb_id=pdb_id, variant=variant, ph=ph
    )
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8")) or {}
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _snapshot_keyset(payload: Mapping[str, Any], key: str) -> Set[str]:
    raw = payload.get(key)
    if not isinstance(raw, list):
        return set()
    out: Set[str] = set()
    for item in raw:
        val = _canonical_ligand_base(str(item))
        if val:
            out.add(val)
    return out


def _resolve_chunk_expected_for_combo(
    *,
    pdb_id: str,
    variant: str,
    ph: str,
    library: str,
    by_combo4_expected: Dict[Tuple[str, str, str, str], Set[str]],
    by_combo3_expected: Dict[Tuple[str, str, str], Set[str]],
    by_combo4_chunks: Dict[Tuple[str, str, str, str], Set[str]],
    by_combo3_chunks: Dict[Tuple[str, str, str], Set[str]],
) -> tuple[Set[str], Set[str], str]:
    key4 = (str(pdb_id).upper(), str(variant).upper(), str(ph), str(library or ""))
    if key4 in by_combo4_expected or key4 in by_combo4_chunks:
        return (
            set(by_combo4_expected.get(key4, set())),
            set(by_combo4_chunks.get(key4, set())),
            "distributed_chunk_plan_exact",
        )

    key3 = (str(pdb_id).upper(), str(variant).upper(), str(ph))
    if key3 in by_combo3_expected or key3 in by_combo3_chunks:
        return (
            set(by_combo3_expected.get(key3, set())),
            set(by_combo3_chunks.get(key3, set())),
            "distributed_chunk_plan_fallback",
        )

    return set(), set(), ""


def _chunk_plan_counts(
    chunk_ids: Set[str],
    chunk_results: Dict[str, Dict[str, Any]],
) -> tuple[int, int, int]:
    completed = 0
    failed = 0
    missing = 0
    for chunk_id in chunk_ids:
        payload = chunk_results.get(str(chunk_id))
        if not isinstance(payload, dict):
            missing += 1
            continue
        status = str(payload.get("status") or "").strip().lower()
        if status == "completed":
            completed += 1
        elif status in {"failed", "terminal_failed"}:
            failed += 1
        else:
            missing += 1
    return completed, failed, missing


def _load_docking_sets(
    docking_summary_csv: Path,
) -> Tuple[Set[str], Set[str]]:
    attempted: Set[str] = set()
    scored: Set[str] = set()
    if not docking_summary_csv.exists():
        return attempted, scored

    with docking_summary_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        ligand_col = next(
            (
                c
                for c in ("Ligand", "ligand", "ligand_file", "ligand_id", "Ligand_ID")
                if c in fieldnames
            ),
            None,
        )
        stage_cols = [c for c in fieldnames if c.lower().startswith("stage")]
        for row in reader:
            ligand_raw = str(row.get(ligand_col or "", "")).strip()
            if not ligand_raw:
                continue
            base = _canonical_ligand_base(ligand_raw)
            if not base:
                continue
            attempted.add(base)
            if not stage_cols:
                scored.add(base)
                continue
            if any(str(row.get(col, "")).strip() for col in stage_cols):
                scored.add(base)
    return attempted, scored


def _marker_stage(path: Path) -> str:
    for token in reversed(path.parts):
        part = token.lower()
        if "stage1" in part:
            return "stage1"
        if "stage2" in part:
            return "stage2"
        if "stage3" in part:
            return "stage3"
    if "post_docked" in str(path).lower():
        return "post"
    return "unknown"


def _marker_within_run_window(
    marker: Path,
    start_ts: Optional[float],
    end_ts: Optional[float],
) -> bool:
    if start_ts is None or end_ts is None:
        return True
    try:
        ts = float(marker.stat().st_mtime)
    except Exception:
        return False
    return (start_ts - _FAILED_MARKER_SKEW_SEC) <= ts <= (
        end_ts + _FAILED_MARKER_SKEW_SEC
    )


def _load_failed_markers(
    combo_docked_dir: Path,
    combo_post_dir: Path,
    *,
    run_start_ts: Optional[float],
    run_end_ts: Optional[float],
) -> Dict[str, Set[str]]:
    out: Dict[str, Set[str]] = {
        "stage1": set(),
        "stage2": set(),
        "stage3": set(),
        "post": set(),
        "unknown": set(),
    }
    roots: List[Tuple[Path, Optional[str]]] = [
        (combo_docked_dir, None),
        (combo_post_dir, "post"),
    ]
    for root, default_stage in roots:
        if not root.exists():
            continue
        for marker in root.rglob("*.failed.txt"):
            if not _marker_within_run_window(marker, run_start_ts, run_end_ts):
                continue
            base = _canonical_ligand_base(marker.name)
            if not base:
                continue
            stage = default_stage or _marker_stage(marker)
            out.setdefault(stage, set()).add(base)
    return out


def _load_post_sets(
    master_rows_csv: Path, pdb_id: str, variant: str, ph: str
) -> Tuple[Set[str], Set[str], Set[str]]:
    """
    Returns (present_non_control, scored_non_control, control_bases)
    """
    present: Set[str] = set()
    scored: Set[str] = set()
    controls: Set[str] = set()
    if not master_rows_csv.exists():
        return present, scored, controls

    with master_rows_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if str(row.get("pdb_id", "")).strip().upper() != pdb_id:
                continue
            if str(row.get("variant", "")).strip().upper() != variant:
                continue
            if str(row.get("ph_label", "")).strip() != ph:
                continue

            base = str(row.get("ligand_base", "")).strip()
            if not base:
                base = _canonical_ligand_base(str(row.get("ligand_file", "")))
            if not base:
                continue

            if _as_bool(row.get("is_control")):
                controls.add(base)
                continue

            present.add(base)
            if str(row.get("final_score", "")).strip():
                scored.add(base)
    return present, scored, controls


def _sample(items: Iterable[str], n: int = 12) -> List[str]:
    return sorted(set(items))[:n]


def _pipeline_window_from_log(
    path: Path,
) -> Tuple[Optional[dt.datetime], Optional[dt.datetime], Optional[float]]:
    if not path.exists():
        return None, None, None
    fmt = "%Y-%m-%d %H:%M:%S,%f"
    first: Optional[dt.datetime] = None
    last: Optional[dt.datetime] = None
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            match = _TS_PREFIX_RE.match(line)
            if not match:
                continue
            try:
                ts = dt.datetime.strptime(match.group(1), fmt)
            except Exception:
                continue
            if first is None:
                first = ts
            last = ts
    if first is None or last is None:
        return first, last, None
    return first, last, max(0.0, float((last - first).total_seconds()))


def _prep_required_paths(
    processed_root: Path, pdb_id: str, variant: str, ph: str
) -> List[Path]:
    if variant.upper() in {"LEGACY", "BASE", ""}:
        receptor_root = processed_root / pdb_id / "receptor"
    else:
        receptor_root = processed_root / pdb_id / variant / "receptor"
    ph_root = receptor_root / "ph_ensemble"
    return [
        ph_root / f"{pdb_id}_{ph}.withH.pdb",
        ph_root / f"{pdb_id}_{ph}.pdbqt",
        receptor_root / f"{pdb_id}_cleaned.pdb",
    ]


@dataclass(frozen=True)
class ComboResult:
    pdb_id: str
    variant: str
    ph: str
    library: str
    status: str
    expected_source: str
    expected_resolved: bool
    expected_resolution_reason: str
    expected_n: int
    prep_n: int
    queued_stage1_n: int
    chunk_expected_mode: str
    chunk_plan_total_n: int
    chunk_plan_completed_n: int
    chunk_plan_failed_n: int
    chunk_plan_missing_result_n: int
    prep_artifacts_required_n: int
    prep_artifacts_found_n: int
    prep_artifacts_missing: List[str]
    docking_attempted_n: int
    docking_scored_n: int
    post_present_n: int
    post_scored_n: int
    controls_post_n: int
    known_failed_stage1_n: int
    known_failed_stage2_n: int
    known_failed_stage3_n: int
    known_failed_post_n: int
    known_failed_n: int
    missing_prep_unexplained_n: int
    missing_docking_unexplained_n: int
    missing_post_unexplained_n: int
    pass_integrity: bool
    failure_reasons: List[str]
    missing_examples: Dict[str, List[str]]


def compute_throughput_integrity(
    *,
    repo_root: Path,
    run_id: str,
    logger: logging.Logger,
) -> Dict[str, Any]:
    manifest_dir = _env_run_dir(
        "MANIFESTS_DIR",
        run_output_dir(repo_root, "manifests", run_id),
        run_id,
    )
    manifest_path = manifest_dir / "run_manifest.yaml"
    manifest = _read_manifest(manifest_path)
    proteins = manifest.get("proteins") if isinstance(manifest, dict) else None
    proteins = proteins if isinstance(proteins, dict) else {}

    paths = manifest.get("paths") if isinstance(manifest, dict) else None
    paths = paths if isinstance(paths, dict) else {}
    run_dir = Path(
        str(
            paths.get("run_dir")
            or run_output_dir(repo_root, "configs", run_id)
        )
    )
    prepped_root = Path(
        str(paths.get("prepped_ligands_dir") or (repo_root / "prepped_ligands"))
    )
    docked_root = Path(str(paths.get("docked_dir") or output_root(repo_root, "docked")))
    processed_root = Path(
        str(
            paths.get("processed_pdb_dir")
            or run_output_dir(repo_root, "processed_pdbs", run_id)
        )
    )
    post_root = Path(
        str(
            paths.get("post_docked_dir")
            or _env_root("POST_DOCKED_DIR", output_root(repo_root, "post_docked"))
        )
    )

    data_dir = _env_run_dir("DATA_DIR", run_output_dir(repo_root, "data", run_id), run_id)
    master_rows_csv = data_dir / "master_rows.csv"
    pipeline_log = docked_root / run_id / "logs" / "pipeline.log"
    run_start, run_end, wall_seconds = _pipeline_window_from_log(pipeline_log)
    run_start_ts = run_start.timestamp() if run_start is not None else None
    run_end_ts = run_end.timestamp() if run_end is not None else None

    if _chunk_plan_checks_enabled():
        (
            chunk_plan_by_combo4_expected,
            chunk_plan_by_combo3_expected,
            chunk_plan_by_combo4_ids,
            chunk_plan_by_combo3_ids,
        ) = _load_distributed_chunk_plans(manifest_path)
        chunk_results = _load_distributed_chunk_results(manifest_path)
    else:
        logger.info(
            "%s action=chunk_plan_checks status=disabled reason=ATLAS_DISABLE_COMBO_CHUNKS",
            COMPONENT,
        )
        chunk_plan_by_combo4_expected = {}
        chunk_plan_by_combo3_expected = {}
        chunk_plan_by_combo4_ids = {}
        chunk_plan_by_combo3_ids = {}
        chunk_results = {}

    combos: List[ComboResult] = []

    for key, entry_any in sorted(proteins.items()):
        combo = _parse_combo_key(str(key))
        if combo is None:
            continue
        pdb_id, variant, ph = combo
        entry = entry_any if isinstance(entry_any, dict) else {}
        entry_status = str(entry.get("status") or "").strip().lower()
        combo_completed = _entry_completed(entry)
        library = str(entry.get("library") or "").strip()

        combo_config_dir = run_dir / pdb_id / variant / ph
        combo_docked_dir = docked_root / run_id / pdb_id / variant / ph
        combo_post_dir = post_root / run_id / pdb_id / variant / ph
        docking_summary_csv = _resolve_docking_summary_csv(combo_docked_dir, library)

        expected_library: Set[str] = set()
        prepped_library: Set[str] = set()
        expected_source = "stage1_configs"
        if library:
            expected_library, prepped_library = _read_library_manifest(
                prepped_root / library
            )
        stage1_queued = _load_stage1_configs(combo_config_dir)
        coverage_snapshot = _load_combo_coverage_snapshot(
            manifest_path,
            pdb_id=pdb_id,
            variant=variant,
            ph=ph,
        )
        snapshot_expected = _snapshot_keyset(coverage_snapshot, "expected_ligands")
        snapshot_docking = _snapshot_keyset(
            coverage_snapshot, "docking_scored_ligands"
        )
        snapshot_post = _snapshot_keyset(coverage_snapshot, "post_scored_ligands")

        chunk_expected, chunk_ids, chunk_reason = _resolve_chunk_expected_for_combo(
            pdb_id=pdb_id,
            variant=variant,
            ph=ph,
            library=library,
            by_combo4_expected=chunk_plan_by_combo4_expected,
            by_combo3_expected=chunk_plan_by_combo3_expected,
            by_combo4_chunks=chunk_plan_by_combo4_ids,
            by_combo3_chunks=chunk_plan_by_combo3_ids,
        )

        if chunk_expected:
            expected = set(chunk_expected)
            prep_set = (prepped_library & expected) if prepped_library else set(expected)
            expected_source = "distributed_chunk_plan"
            expected_resolution_reason = chunk_reason or "distributed_chunk_plan"
            chunk_expected_mode = "planned_full"
        elif expected_library:
            expected = expected_library
            prep_set = prepped_library
            expected_source = "library_manifest"
            expected_resolution_reason = "library_manifest"
            chunk_expected_mode = "none"
        elif stage1_queued:
            expected = stage1_queued
            prep_set = set(stage1_queued)
            expected_source = "stage1_configs"
            expected_resolution_reason = "stage1_configs"
            chunk_expected_mode = "none"
        elif snapshot_expected:
            expected = set(snapshot_expected)
            prep_set = set(snapshot_expected)
            expected_source = "coverage_snapshot"
            expected_resolution_reason = "coverage_snapshot"
            chunk_expected_mode = "none"
        else:
            expected = set()
            prep_set = set()
            expected_source = "observed_fallback"
            expected_resolution_reason = (
                "library_unresolved_and_stage1_missing" if library else "library_and_stage1_missing"
            )
            chunk_expected_mode = "none"
        expected_resolved = len(expected) > 0

        chunk_plan_total_n = len(chunk_ids)
        (
            chunk_plan_completed_n,
            chunk_plan_failed_n,
            chunk_plan_missing_result_n,
        ) = _chunk_plan_counts(chunk_ids, chunk_results)

        live_docking_attempted_all, live_docking_scored_all = _load_docking_sets(
            docking_summary_csv
        )
        if coverage_snapshot:
            docking_attempted_all = set(snapshot_docking)
            docking_scored_all = set(snapshot_docking)
            # Coverage snapshots can be written before late chunk aggregation updates
            # the canonical docking summary. Prefer the summary CSV when it has
            # strictly better evidence, mirroring the post/master_rows refresh below.
            if live_docking_scored_all and (
                not docking_scored_all
                or len(live_docking_scored_all) > len(docking_scored_all)
                or not live_docking_scored_all.issubset(docking_scored_all)
            ):
                docking_attempted_all = set(live_docking_attempted_all)
                docking_scored_all = set(live_docking_scored_all)
        else:
            docking_attempted_all = set(live_docking_attempted_all)
            docking_scored_all = set(live_docking_scored_all)
        if expected:
            docking_attempted = docking_attempted_all & expected
            docking_scored = docking_scored_all & expected
        else:
            docking_attempted = docking_attempted_all
            docking_scored = docking_scored_all

        master_post_present, master_post_scored, post_controls = _load_post_sets(
            master_rows_csv, pdb_id, variant, ph
        )
        if not expected and (docking_scored_all or master_post_scored or post_controls):
            expected = set(docking_scored_all) | set(master_post_scored) | set(
                post_controls
            )
            prep_set = set(expected)
            expected_source = "observed_outputs"
            expected_resolution_reason = "observed_docking_or_post_outputs"
            expected_resolved = len(expected) > 0
        if coverage_snapshot:
            post_present_all = set(snapshot_post)
            post_scored_all = set(snapshot_post)
            # Coverage snapshots are best-effort and may be generated before the final
            # canonical post CSV is written. Prefer master_rows when it has strictly
            # better evidence for scored post ligands.
            if master_post_scored and (
                not post_scored_all
                or len(master_post_scored) > len(post_scored_all)
                or not master_post_scored.issubset(post_scored_all)
            ):
                post_present_all = set(master_post_present)
                post_scored_all = set(master_post_scored)
        else:
            post_present_all = set(master_post_present)
            post_scored_all = set(master_post_scored)
        if expected:
            post_present = post_present_all & expected
            post_scored = post_scored_all & expected
        else:
            post_present = post_present_all
            post_scored = post_scored_all

        failed_markers = _load_failed_markers(
            combo_docked_dir,
            combo_post_dir,
            run_start_ts=run_start_ts,
            run_end_ts=run_end_ts,
        )

        missing_prep = expected - prep_set
        missing_docking = expected - docking_scored
        post_accounted = post_scored | post_controls
        missing_post = expected - post_accounted

        known_failed_stage1 = missing_docking & failed_markers.get("stage1", set())
        known_failed_stage2 = missing_docking & failed_markers.get("stage2", set())
        known_failed_stage3 = missing_docking & failed_markers.get("stage3", set())
        known_failed_docking = (
            known_failed_stage1 | known_failed_stage2 | known_failed_stage3
        )
        known_failed_post = missing_post & failed_markers.get("post", set())
        known_failed = known_failed_docking | known_failed_post
        missing_prep_unexplained = missing_prep
        missing_docking_unexplained = missing_docking - known_failed_docking
        missing_post_unexplained = missing_post - known_failed_post

        prep_required_paths = _prep_required_paths(processed_root, pdb_id, variant, ph)
        prep_artifacts_missing = [str(p) for p in prep_required_paths if not p.exists()]
        prep_artifacts_found_n = len(prep_required_paths) - len(prep_artifacts_missing)

        failure_reasons: List[str] = []
        if not combo_completed:
            failure_reasons.append("combo_not_completed")
        if (
            combo_completed
            and not expected_resolved
            and run_id not in _ALLOW_EMPTY_EXPECTED_RUN_IDS
        ):
            failure_reasons.append("expected_set_unresolved")
        if combo_completed and prep_artifacts_missing:
            failure_reasons.append("prep_artifacts_missing")
        if missing_prep_unexplained:
            failure_reasons.append("missing_prep_unexplained")
        if missing_docking_unexplained:
            failure_reasons.append("missing_docking_unexplained")
        if missing_post_unexplained:
            failure_reasons.append("missing_post_unexplained")
        if combo_completed and chunk_plan_total_n > 0:
            if chunk_plan_completed_n < chunk_plan_total_n:
                failure_reasons.append("chunk_plan_incomplete")
            if chunk_plan_failed_n > 0:
                failure_reasons.append("chunk_plan_failed")
            if chunk_plan_missing_result_n > 0:
                failure_reasons.append("chunk_results_missing")

        pass_integrity = len(failure_reasons) == 0

        combos.append(
            ComboResult(
                pdb_id=pdb_id,
                variant=variant,
                ph=ph,
                library=library,
                status=entry_status,
                expected_source=expected_source,
                expected_resolved=expected_resolved,
                expected_resolution_reason=expected_resolution_reason,
                expected_n=len(expected),
                prep_n=len(prep_set & expected if expected else prep_set),
                queued_stage1_n=len(stage1_queued),
                chunk_expected_mode=chunk_expected_mode,
                chunk_plan_total_n=chunk_plan_total_n,
                chunk_plan_completed_n=chunk_plan_completed_n,
                chunk_plan_failed_n=chunk_plan_failed_n,
                chunk_plan_missing_result_n=chunk_plan_missing_result_n,
                prep_artifacts_required_n=len(prep_required_paths),
                prep_artifacts_found_n=prep_artifacts_found_n,
                prep_artifacts_missing=prep_artifacts_missing,
                docking_attempted_n=len(docking_attempted),
                docking_scored_n=len(docking_scored),
                post_present_n=len(post_present),
                post_scored_n=len(post_scored),
                controls_post_n=len(post_controls),
                known_failed_stage1_n=len(known_failed_stage1),
                known_failed_stage2_n=len(known_failed_stage2),
                known_failed_stage3_n=len(known_failed_stage3),
                known_failed_post_n=len(known_failed_post),
                known_failed_n=len(known_failed),
                missing_prep_unexplained_n=len(missing_prep_unexplained),
                missing_docking_unexplained_n=len(missing_docking_unexplained),
                missing_post_unexplained_n=len(missing_post_unexplained),
                pass_integrity=pass_integrity,
                failure_reasons=failure_reasons,
                missing_examples={
                    "prep": _sample(missing_prep_unexplained),
                    "docking": _sample(missing_docking_unexplained),
                    "post": _sample(missing_post_unexplained),
                    "known_failed": _sample(known_failed),
                },
            )
        )

    expected_total = sum(c.expected_n for c in combos)
    prep_total = sum(c.prep_n for c in combos)
    docking_total = sum(c.docking_scored_n for c in combos)
    post_total = sum(c.post_scored_n for c in combos)
    controls_total = sum(c.controls_post_n for c in combos)
    planned_chunks_total = sum(c.chunk_plan_total_n for c in combos)
    planned_chunks_completed = sum(c.chunk_plan_completed_n for c in combos)
    planned_chunks_failed = sum(c.chunk_plan_failed_n for c in combos)
    planned_chunks_missing = sum(c.chunk_plan_missing_result_n for c in combos)

    unexplained_total = sum(
        c.missing_prep_unexplained_n
        + c.missing_docking_unexplained_n
        + c.missing_post_unexplained_n
        for c in combos
    )
    unresolved_expected_combos = sum(1 for c in combos if not c.expected_resolved)
    missing_prep_artifact_combos = sum(1 for c in combos if c.prep_artifacts_missing)
    failed_combos = sum(1 for c in combos if not c.pass_integrity)

    integrity_pass = failed_combos == 0
    ligands_per_second_expected = (
        (expected_total / wall_seconds) if wall_seconds and wall_seconds > 0 else None
    )
    ligands_per_second_post = (
        (post_total / wall_seconds) if wall_seconds and wall_seconds > 0 else None
    )

    summary = {
        "run_id": run_id,
        "generated_at_utc": dt.datetime.utcnow().replace(microsecond=0).isoformat()
        + "Z",
        "integrity_pass": bool(integrity_pass),
        "wall_seconds": wall_seconds,
        "totals": {
            "combo_count": len(combos),
            "expected_ligands": expected_total,
            "prepped_ligands": prep_total,
            "docking_scored_ligands": docking_total,
            "post_scored_ligands": post_total,
            "post_controls": controls_total,
            "unexplained_missing_total": unexplained_total,
            "failed_combos": failed_combos,
            "unresolved_expected_combos": unresolved_expected_combos,
            "missing_prep_artifact_combos": missing_prep_artifact_combos,
            "planned_chunks_total": planned_chunks_total,
            "planned_chunks_completed": planned_chunks_completed,
            "planned_chunks_failed": planned_chunks_failed,
            "planned_chunks_missing": planned_chunks_missing,
            "chunk_coverage_pass": (
                (planned_chunks_total == planned_chunks_completed)
                if planned_chunks_total > 0
                else True
            ),
            "ligands_per_second_expected": ligands_per_second_expected,
            "ligands_per_second_post": ligands_per_second_post,
        },
        "failure_reasons": [
            {
                "pdb_id": c.pdb_id,
                "variant": c.variant,
                "ph": c.ph,
                "reasons": c.failure_reasons,
            }
            for c in combos
            if c.failure_reasons
        ],
        "combos": [
            {
                "pdb_id": c.pdb_id,
                "variant": c.variant,
                "ph": c.ph,
                "library": c.library,
                "status": c.status,
                "expected_source": c.expected_source,
                "expected_resolved": c.expected_resolved,
                "expected_resolution_reason": c.expected_resolution_reason,
                "expected_n": c.expected_n,
                "prep_n": c.prep_n,
                "queued_stage1_n": c.queued_stage1_n,
                "chunk_expected_mode": c.chunk_expected_mode,
                "chunk_plan_total_n": c.chunk_plan_total_n,
                "chunk_plan_completed_n": c.chunk_plan_completed_n,
                "chunk_plan_failed_n": c.chunk_plan_failed_n,
                "chunk_plan_missing_result_n": c.chunk_plan_missing_result_n,
                "prep_artifacts_required_n": c.prep_artifacts_required_n,
                "prep_artifacts_found_n": c.prep_artifacts_found_n,
                "prep_artifacts_missing": c.prep_artifacts_missing,
                "docking_attempted_n": c.docking_attempted_n,
                "docking_scored_n": c.docking_scored_n,
                "post_present_n": c.post_present_n,
                "post_scored_n": c.post_scored_n,
                "controls_post_n": c.controls_post_n,
                "known_failed_stage1_n": c.known_failed_stage1_n,
                "known_failed_stage2_n": c.known_failed_stage2_n,
                "known_failed_stage3_n": c.known_failed_stage3_n,
                "known_failed_post_n": c.known_failed_post_n,
                "known_failed_n": c.known_failed_n,
                "missing_prep_unexplained_n": c.missing_prep_unexplained_n,
                "missing_docking_unexplained_n": c.missing_docking_unexplained_n,
                "missing_post_unexplained_n": c.missing_post_unexplained_n,
                "pass_integrity": c.pass_integrity,
                "failure_reasons": c.failure_reasons,
                "missing_examples": c.missing_examples,
            }
            for c in combos
        ],
    }

    logger.info(
        "%s action=compute run_id=%s combos=%d expected=%d post=%d integrity_pass=%s",
        COMPONENT,
        run_id,
        len(combos),
        expected_total,
        post_total,
        str(integrity_pass).lower(),
    )
    return summary


def _write_summary_csv(out_csv: Path, summary: Dict[str, Any]) -> None:
    rows = summary.get("combos") or []
    fieldnames = [
        "run_id",
        "pdb_id",
        "variant",
        "ph",
        "library",
        "status",
        "expected_source",
        "expected_resolved",
        "expected_resolution_reason",
        "expected_n",
        "prep_n",
        "queued_stage1_n",
        "chunk_expected_mode",
        "chunk_plan_total_n",
        "chunk_plan_completed_n",
        "chunk_plan_failed_n",
        "chunk_plan_missing_result_n",
        "prep_artifacts_required_n",
        "prep_artifacts_found_n",
        "docking_attempted_n",
        "docking_scored_n",
        "post_present_n",
        "post_scored_n",
        "controls_post_n",
        "known_failed_stage1_n",
        "known_failed_stage2_n",
        "known_failed_stage3_n",
        "known_failed_post_n",
        "known_failed_n",
        "missing_prep_unexplained_n",
        "missing_docking_unexplained_n",
        "missing_post_unexplained_n",
        "pass_integrity",
        "failure_reasons",
    ]
    with out_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            payload = {"run_id": summary.get("run_id")}
            for k in fieldnames:
                if k == "run_id":
                    continue
                if k == "failure_reasons":
                    payload[k] = "|".join(str(x) for x in (row.get(k) or []))
                else:
                    payload[k] = row.get(k)
            writer.writerow(payload)


def run(
    *,
    repo_root: Path,
    run_id: str,
    strict: bool,
    logger: logging.Logger,
) -> int:
    summary = compute_throughput_integrity(
        repo_root=repo_root,
        run_id=run_id,
        logger=logger,
    )
    data_dir = _env_run_dir("DATA_DIR", output_root(repo_root, "data") / run_id, run_id)
    data_dir.mkdir(parents=True, exist_ok=True)
    out_json = data_dir / "throughput_integrity.json"
    out_csv = data_dir / "throughput_integrity.csv"
    out_json.write_text(json.dumps(summary, indent=2, sort_keys=False), encoding="utf-8")
    _write_summary_csv(out_csv, summary)

    integrity_pass = bool(summary.get("integrity_pass", False))
    logger.info(
        "%s action=write run_id=%s json=%s csv=%s strict=%s integrity_pass=%s",
        COMPONENT,
        run_id,
        out_json,
        out_csv,
        str(strict).lower(),
        str(integrity_pass).lower(),
    )
    if strict and not integrity_pass:
        return 2
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate ligand-throughput integrity across prep, docking, and post-docking."
    )
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[2]))
    parser.add_argument("--strict", dest="strict", action="store_true", default=True)
    parser.add_argument("--no-strict", dest="strict", action="store_false")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(list(argv) if argv is not None else None)

    logger = _configure_logging(args.verbose)
    repo_root = Path(args.repo_root).resolve()
    return run(repo_root=repo_root, run_id=str(args.run_id), strict=bool(args.strict), logger=logger)


if __name__ == "__main__":
    raise SystemExit(main())
