from __future__ import annotations

import argparse
import csv
import datetime as _dt
import json
import logging
import os
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, MutableMapping, Sequence

from cli.run_manifest_io import protein_key
from cli.run_manifest_support import (
    default_protein_entry,
    default_stage_entry,
    refresh_summary,
    require_yaml,
)
from post_docking.rescoring.scorch_coverage import (
    combo_coverage_summary_path,
    evaluate_scorch_coverage,
    normalize_combo,
    remove_done_sentinel,
    resolve_top_fraction,
    write_coverage_summary,
)
from post_docking.rescoring.rescore_reranker import (
    find_consensus_csv,
    rerank_consensus_with_scorch,
)
from post_docking.rescoring.scorch_selection import load_control_bases
from post_docking.rescoring.scorch_shards import (
    _output_valid as _scorch_shard_output_valid,
    planned_completed_output_bases,
    planned_completed_output_csvs,
    read_all_shard_records,
)
from post_docking.rescoring.scorch_postprocess import aggregate_combo
from post_docking.rescoring.scorch_types import StageSpec


@dataclass(frozen=True)
class _Scope:
    pdb_id: str
    variant: str
    ph: str


@dataclass(frozen=True)
class ManifestRepairSummary:
    run_id: str
    all_dirs: Path
    manifest_path: Path
    dry_run: bool
    require_scorch: bool
    planned_chunks: int
    completed_chunks: int
    failed_chunks: int
    missing_chunks: int
    planned_ligand_assignments: int
    completed_ligand_assignments: int
    targets: int
    targets_completed: int
    targets_running: int
    targets_failed: int
    total_proteins_scheduled: int
    total_protein_list: list[str]
    total_proteins_completed: int
    total_proteins_failed: int
    consensus_files: int
    scorch_reranked_files: int
    repair_scorch_selection: bool
    consensus_files_checked: int
    consensus_files_sorted: int
    consensus_files_sort_failed: int
    scorch_output_files_cleared: int
    scorch_done_sentinels_cleared: int
    scorch_coverage_summaries_cleared: int
    scorch_input_dirs_cleared: int
    scorch_shard_summaries_cleared: int
    backup_path: Path | None


def _normalize_variant(raw: Any) -> str:
    token = str(raw or "").strip().upper()
    if token in {"", "BASE", "NONE", "NULL"}:
        return "LEGACY"
    return token


def _normalize_ph(raw: Any) -> str:
    token = str(raw or "").strip()
    if token.lower() in {"", "base", "none", "null"}:
        return "base"
    return token


def _read_json(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _read_yaml(path: Path) -> dict[str, Any]:
    yaml_mod = require_yaml()
    if yaml_mod is None:
        raise RuntimeError("PyYAML is required to repair run manifests.")
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        payload = yaml_mod.safe_load(handle) or {}
    if not isinstance(payload, dict):
        return {}
    return payload


def _write_yaml(path: Path, payload: Mapping[str, Any]) -> None:
    yaml_mod = require_yaml()
    if yaml_mod is None:
        raise RuntimeError("PyYAML is required to repair run manifests.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        yaml_mod.safe_dump(dict(payload), handle, default_flow_style=False, sort_keys=False)


def _resolve_all_dirs(raw: str | None) -> Path:
    token = str(raw or "").strip()
    if token and token.lower() != "none":
        return Path(token).expanduser().resolve()
    env_root = str(os.environ.get("ATLAS_ALL_DIRS", "")).strip()
    if env_root:
        return Path(env_root).expanduser().resolve()
    return Path(__file__).resolve().parents[3]


def _manifest_dir(all_dirs: Path, run_id: str) -> Path:
    return all_dirs / "outputs" / "manifests" / str(run_id)


def _docked_dir(all_dirs: Path, run_id: str) -> Path:
    return all_dirs / "outputs" / "docked" / str(run_id)


def _post_docked_dir(all_dirs: Path, run_id: str) -> Path:
    return all_dirs / "outputs" / "post_docked" / str(run_id)


def _processed_root(all_dirs: Path) -> Path:
    return all_dirs / "outputs" / "processed_pdbs"


def _iter_plan_chunks(dist_dir: Path) -> list[dict[str, Any]]:
    chunks: list[dict[str, Any]] = []
    for plan_path in sorted(dist_dir.glob("combo_chunks_*.json")):
        payload = _read_json(plan_path)
        raw_chunks = payload.get("chunks") if isinstance(payload, dict) else payload
        if not isinstance(raw_chunks, list):
            continue
        for idx, raw in enumerate(raw_chunks):
            if isinstance(raw, dict):
                chunk = dict(raw)
                chunk["_plan_path"] = str(plan_path)
                chunk["_plan_index"] = idx
                chunks.append(chunk)
    return chunks


def _chunk_scope(chunk: Mapping[str, Any]) -> _Scope | None:
    pdb_id = str(chunk.get("pdb_id") or "").strip().upper()
    if not pdb_id:
        pdb_file = str(chunk.get("pdb_file") or "").strip()
        if pdb_file:
            pdb_id = Path(pdb_file).stem.strip().upper()
    if not pdb_id:
        return None
    return _Scope(
        pdb_id=pdb_id,
        variant=_normalize_variant(chunk.get("variant_label")),
        ph=_normalize_ph(chunk.get("ph_tag")),
    )


def _chunk_ligand_count(chunk: Mapping[str, Any]) -> int:
    try:
        return int(chunk.get("ligand_count") or 0)
    except Exception:
        pass
    ligands = chunk.get("ligand_bases") or []
    return len(ligands) if isinstance(ligands, list) else 0


def _completed_ligand_count(
    result: Mapping[str, Any] | None,
    chunk: Mapping[str, Any],
) -> int:
    if isinstance(result, Mapping):
        for key in ("ligand_count", "completed_ligands", "n_ligands"):
            try:
                value = int(result.get(key) or 0)
            except Exception:
                value = 0
            if value > 0:
                return value
    return _chunk_ligand_count(chunk)


def _load_chunk_results(dist_dir: Path) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for result_path in sorted((dist_dir / "chunk_results").glob("*.json")):
        try:
            payload = _read_json(result_path)
        except Exception:
            continue
        if not isinstance(payload, dict):
            continue
        chunk_id = str(payload.get("chunk_id") or result_path.stem).split(".", 1)[0]
        if chunk_id:
            out[chunk_id] = payload
    return out


def _chunk_id(chunk: Mapping[str, Any]) -> str:
    token = str(chunk.get("chunk_id") or "").strip()
    if token:
        return token
    plan_name = Path(str(chunk.get("_plan_path") or "combo_chunks.json")).name
    return f"missing_id:{plan_name}:{chunk.get('_plan_index', 0)}"


def _utc_now_iso() -> str:
    return _dt.datetime.utcnow().replace(microsecond=0).isoformat() + "Z"


def _backup_path(manifest_path: Path) -> Path:
    stamp = _dt.datetime.utcnow().strftime("%Y%m%d_%H%M%S")
    candidate = manifest_path.with_suffix(manifest_path.suffix + f".bak.{stamp}")
    suffix = 1
    while candidate.exists():
        candidate = manifest_path.with_suffix(
            manifest_path.suffix + f".bak.{stamp}.{suffix}"
        )
        suffix += 1
    return candidate


def _count_named_files(root: Path, filename: str) -> int:
    if not root.exists():
        return 0
    return sum(1 for path in root.rglob(filename) if path.is_file())


def _scorch_required_auto(post_root: Path) -> bool:
    if not post_root.exists():
        return False
    return any(post_root.glob("*/scorch_scores_all.csv")) or any(
        post_root.glob("*/consensus_reranked_scorch.csv")
    )


def _has_consensus(docked_root: Path, scope: _Scope) -> bool:
    candidates = [docked_root / scope.pdb_id / "consensus_docking_scores.csv"]
    if scope.variant != "LEGACY":
        candidates.append(docked_root / scope.pdb_id / scope.variant / "consensus_docking_scores.csv")
        if scope.ph != "base":
            candidates.append(
                docked_root / scope.pdb_id / scope.variant / scope.ph / "consensus_docking_scores.csv"
            )
    if scope.ph != "base":
        candidates.append(docked_root / scope.pdb_id / scope.ph / "consensus_docking_scores.csv")
    return any(path.is_file() and path.stat().st_size > 0 for path in candidates)


def _scorch_csv_candidates(post_root: Path, scope: _Scope) -> list[Path]:
    candidates = [post_root / scope.pdb_id / "consensus_reranked_scorch.csv"]
    candidates.append(
        post_root / scope.pdb_id / scope.variant / scope.ph / "consensus_reranked_scorch.csv"
    )
    if scope.variant != "LEGACY":
        candidates.append(post_root / scope.pdb_id / scope.variant / "consensus_reranked_scorch.csv")
        if scope.ph != "base":
            candidates.append(
                post_root / scope.pdb_id / scope.variant / scope.ph / "consensus_reranked_scorch.csv"
            )
    if scope.ph != "base":
        candidates.append(post_root / scope.pdb_id / scope.ph / "consensus_reranked_scorch.csv")
    return candidates


def _find_scorch_csv(post_root: Path, scope: _Scope) -> Path | None:
    for path in _scorch_csv_candidates(post_root, scope):
        if path.is_file() and path.stat().st_size > 0:
            return path
    return None


def _coverage_control_bases(processed_root: Path | None, scope: _Scope) -> set[str]:
    if processed_root is None:
        return set()
    try:
        return set(
            load_control_bases(
                Path(processed_root),
                scope.pdb_id,
                logging.getLogger("manifest_repair.scorch_coverage"),
                component="[scorch-coverage]",
            )
        )
    except Exception:
        return set()


def _scope_combo(scope: _Scope) -> tuple[str, str, str]:
    variant = "" if scope.variant in {"", "LEGACY", "BASE"} else scope.variant
    ph = "" if scope.ph in {"", "base", "BASE"} else scope.ph
    return normalize_combo(scope.pdb_id, variant, ph)


def _manifest_cfg_for_shards(manifest_run_root: Path | None) -> tuple[dict[str, str], str] | None:
    if manifest_run_root is None:
        return None
    run_root = Path(manifest_run_root)
    if not run_root.name:
        return None
    return {"MANIFESTS_DIR": str(run_root.parent)}, run_root.name


def _reaggregate_scorch_from_shards(
    post_root: Path,
    scope: _Scope,
    *,
    docked_root: Path | None,
    manifest_run_root: Path | None,
    dry_run: bool,
) -> str:
    if dry_run or docked_root is None:
        return ""
    manifest_ctx = _manifest_cfg_for_shards(manifest_run_root)
    if manifest_ctx is None:
        return ""
    cfg, run_id = manifest_ctx
    combo = _scope_combo(scope)
    records = read_all_shard_records(cfg, run_id)
    specs_by_mode: dict[str, dict[tuple[str, str, str], StageSpec]] = {}
    decoy_prefix = "decoys"
    for record in records:
        raw_combo = record.get("combo") if isinstance(record.get("combo"), list) else []
        record_combo = normalize_combo(
            raw_combo[0] if len(raw_combo) > 0 else "",
            raw_combo[1] if len(raw_combo) > 1 else "",
            raw_combo[2] if len(raw_combo) > 2 else "",
        )
        if record_combo != combo:
            continue
        mode = str(record.get("run_mode") or "").strip().lower()
        if mode not in {"fda", "dud"}:
            continue
        prefix = str(record.get("decoy_prefix") or "").strip()
        if prefix:
            decoy_prefix = prefix
        spec_raw = record.get("spec") if isinstance(record.get("spec"), Mapping) else {}
        source = str(spec_raw.get("source") or "").strip()
        stage_dir = str(spec_raw.get("stage_dir") or "").strip()
        output_name = str(spec_raw.get("output_name") or "").strip()
        if not (source and stage_dir and output_name):
            continue
        specs_by_mode.setdefault(mode, {})[
            (source, stage_dir, output_name)
        ] = StageSpec(source, stage_dir, output_name)
    if not specs_by_mode:
        return ""

    logger = logging.getLogger("manifest_repair.scorch_reaggregate")
    aggregate_paths: list[Path] = []
    for mode, spec_map in sorted(specs_by_mode.items()):
        specs = list(spec_map.values())
        planned_csvs = planned_completed_output_csvs(
            cfg,
            run_id,
            combo,
            specs,
            mode,
            decoy_prefix,
        )
        if planned_csvs is None:
            continue
        if not planned_csvs:
            logger.info(
                "[manifest-repair.scorch] action=reaggregate status=skip reason=incomplete_plan_bound_parts run_id=%s pdb_id=%s mode=%s",
                run_id,
                scope.pdb_id,
                mode,
            )
            return ""
        planned_bases = planned_completed_output_bases(
            cfg,
            run_id,
            combo,
            specs,
            mode,
            decoy_prefix,
        )
        output_name = (
            "scorch_scores_all.csv"
            if mode == "fda"
            else f"{decoy_prefix}_scorch_scores_all.csv"
        )
        agg_path = aggregate_combo(
            post_root,
            specs,
            combo,
            logger,
            component="[manifest-repair.scorch]",
            run_mode=mode,
            output_name=output_name,
            decoy_prefix=decoy_prefix,
            input_csvs=planned_csvs,
            allowed_ligand_bases=planned_bases,
        )
        if agg_path is None:
            return ""
        aggregate_paths.append(agg_path)

    if not aggregate_paths:
        return ""
    dock_combo_dir = Path(docked_root) / combo[0] / combo[1] / combo[2]
    consensus_csv = find_consensus_csv(dock_combo_dir)
    if consensus_csv is None:
        logger.info(
            "[manifest-repair.scorch] action=rerank status=skip reason=missing_consensus run_id=%s pdb_id=%s",
            run_id,
            scope.pdb_id,
        )
        return ""
    out_csv = aggregate_paths[0].parent / "consensus_reranked_scorch.csv"
    ok = rerank_consensus_with_scorch(
        consensus_csv,
        aggregate_paths,
        out_csv,
        logger,
        overwrite=True,
        decoy_prefix=decoy_prefix,
    )
    if not ok:
        return ""
    logger.info(
        "[manifest-repair.scorch] action=reaggregate status=ok run_id=%s pdb_id=%s modes=%s",
        run_id,
        scope.pdb_id,
        ",".join(sorted(specs_by_mode)),
    )
    return decoy_prefix


def _has_scorch(
    post_root: Path,
    scope: _Scope,
    *,
    docked_root: Path | None = None,
    manifest_run_root: Path | None = None,
    processed_root: Path | None = None,
    require_coverage: bool = False,
    dry_run: bool = False,
    trust_existing_coverage: bool = False,
) -> bool:
    if not require_coverage:
        return _find_scorch_csv(post_root, scope) is not None
    if docked_root is None:
        return False
    if trust_existing_coverage and _existing_scorch_coverage_complete(post_root, scope):
        return True
    decoy_prefix = _reaggregate_scorch_from_shards(
        post_root,
        scope,
        docked_root=docked_root,
        manifest_run_root=manifest_run_root,
        dry_run=dry_run,
    )
    summary = evaluate_scorch_coverage(
        docked_run_root=docked_root,
        post_run_root=post_root,
        pdb_id=scope.pdb_id,
        variant=scope.variant,
        ph=scope.ph,
        top_fraction=resolve_top_fraction({}),
        decoy_prefix=decoy_prefix or "dud",
        control_bases=_coverage_control_bases(processed_root, scope),
        manifest_run_root=manifest_run_root,
    )
    if not dry_run:
        write_coverage_summary(post_root, summary)
    if not summary.all_complete and not dry_run:
        remove_done_sentinel(post_root, summary)
    return bool(summary.all_complete)


def _canonical_scorch_done_dir(post_root: Path, scope: _Scope) -> Path:
    combo_dir = post_root / scope.pdb_id
    if scope.variant != "LEGACY":
        combo_dir = combo_dir / scope.variant
    if scope.ph != "base":
        combo_dir = combo_dir / scope.ph
    return combo_dir / "scorch"


def _stream_coverage_complete(raw: Any) -> bool:
    if not isinstance(raw, Mapping):
        return False
    if raw.get("complete") is not True:
        return False
    try:
        invalid = int(raw.get("invalid_rows") or 0)
        expected = int(raw.get("expected_rows") or 0)
        raw_rows = int(raw.get("raw_rows") or 0)
        final_rows = int(raw.get("final_rescored_rows") or 0)
    except Exception:
        return False
    return invalid == 0 and raw_rows >= expected and final_rows >= expected


def _existing_scorch_coverage_complete(post_root: Path, scope: _Scope) -> bool:
    if _find_scorch_csv(post_root, scope) is None:
        return False
    summary_path = combo_coverage_summary_path(
        post_root,
        normalize_combo(scope.pdb_id, scope.variant, scope.ph),
    )
    payload = _read_json_maybe(summary_path)
    if not payload or payload.get("all_complete") is not True:
        return False
    streams = payload.get("streams") if isinstance(payload.get("streams"), Mapping) else payload
    return _stream_coverage_complete(streams.get("fda")) and _stream_coverage_complete(
        streams.get("dud")
    )


def _ensure_scorch_done_sentinel(
    post_root: Path,
    scope: _Scope,
    *,
    docked_root: Path | None = None,
    manifest_run_root: Path | None = None,
    processed_root: Path | None = None,
    require_coverage: bool = False,
    dry_run: bool = False,
    trust_existing_coverage: bool = False,
) -> bool:
    if require_coverage and not _has_scorch(
        post_root,
        scope,
        docked_root=docked_root,
        manifest_run_root=manifest_run_root,
        processed_root=processed_root,
        require_coverage=True,
        dry_run=dry_run,
        trust_existing_coverage=trust_existing_coverage,
    ):
        return False
    scorch_csv = _find_scorch_csv(post_root, scope)
    if scorch_csv is None:
        return False
    done_dirs = {
        scorch_csv.parent / "scorch",
        _canonical_scorch_done_dir(post_root, scope),
    }
    wrote = False
    for done_dir in done_dirs:
        sentinel = done_dir / "_DONE"
        if sentinel.exists():
            continue
        done_dir.mkdir(parents=True, exist_ok=True)
        sentinel.write_text("ok\n", encoding="utf-8")
        wrote = True
    return wrote


def _read_json_maybe(path: Path) -> dict[str, Any] | None:
    try:
        payload = _read_json(path)
    except Exception:
        return None
    return payload if isinstance(payload, dict) else None


def _scope_from_scorch_payload(payload: Mapping[str, Any]) -> _Scope | None:
    pdb_id = str(payload.get("pdb_id") or "").strip().upper()
    if not pdb_id:
        return None
    variant_raw = str(payload.get("variant") or "LEGACY").strip()
    variant = "LEGACY" if variant_raw in {"", "*"} else _normalize_variant(variant_raw)
    return _Scope(
        pdb_id=pdb_id,
        variant=variant,
        ph=_normalize_ph(payload.get("ph")),
    )


def _scope_from_scorch_plan(payload: Mapping[str, Any]) -> _Scope | None:
    combo_raw = payload.get("combo")
    if not isinstance(combo_raw, (list, tuple)):
        return None
    if not combo_raw:
        return None
    pdb_id = str(combo_raw[0] if len(combo_raw) > 0 else "").strip().upper()
    if not pdb_id:
        return None
    return _Scope(
        pdb_id=pdb_id,
        variant=_normalize_variant(combo_raw[1] if len(combo_raw) > 1 else ""),
        ph=_normalize_ph(combo_raw[2] if len(combo_raw) > 2 else ""),
    )


def _shard_ids_from_plan(payload: Mapping[str, Any]) -> set[str]:
    shards = payload.get("shards")
    if not isinstance(shards, (list, tuple)):
        return set()
    ids: set[str] = set()
    for record in shards:
        if not isinstance(record, Mapping):
            continue
        shard_id = str(record.get("shard_id") or "").strip()
        if shard_id:
            ids.add(shard_id)
    return ids


def _scorch_result_has_valid_output(payload: Mapping[str, Any]) -> bool:
    status = str(payload.get("status") or "").strip().lower()
    if status != "completed":
        return False
    if bool(payload.get("no_output", False)):
        return True
    output_raw = str(payload.get("output_csv") or "").strip()
    if not output_raw:
        return False
    output_path = Path(output_raw)
    try:
        return output_path.is_file() and output_path.stat().st_size > 0
    except Exception:
        return False


def _scorch_record_expected_rows(record: Mapping[str, Any]) -> int:
    explicit = record.get("expected_output_rows")
    if explicit is not None:
        try:
            return max(0, int(explicit or 0))
        except Exception:
            pass
    allowed_raw = record.get("allowed_bases")
    control_raw = record.get("control_bases")
    if isinstance(allowed_raw, list):
        allowed = {str(item) for item in allowed_raw if str(item).strip()}
        controls = (
            {str(item) for item in control_raw if str(item).strip()}
            if isinstance(control_raw, list)
            else set()
        )
        return max(0, len(allowed - controls))
    try:
        return max(0, int(record.get("allowed_count") or 0))
    except Exception:
        return 0


def _scorch_ligand_key(value: Any) -> str:
    return Path(str(value or "").strip()).stem.strip()


def _scorch_csv_rows_and_ligands(path: Path) -> tuple[int, set[str]]:
    rows = 0
    ligands: set[str] = set()
    try:
        with path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                rows += 1
                for key in (
                    "Ligand_ID",
                    "ligand",
                    "ligand_id",
                    "Ligand",
                    "ligand_base",
                    "base",
                ):
                    token = _scorch_ligand_key(row.get(key))
                    if token:
                        ligands.add(token)
                        break
    except Exception:
        return 0, set()
    return rows, ligands


def _scorch_record_output_valid(record: Mapping[str, Any]) -> tuple[bool, int]:
    return _scorch_shard_output_valid(record)


def _write_repaired_scorch_result(
    result_path: Path,
    *,
    run_id: str,
    record: Mapping[str, Any],
    row_count: int,
    dry_run: bool,
) -> bool:
    if dry_run:
        return True
    shard_id = str(record.get("shard_id") or "").strip()
    if not shard_id:
        return False
    result_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "atlas.scorch_shard_result.v1",
        "run_id": str(run_id),
        "shard_id": shard_id,
        "status": "completed",
        "completed_at": float(time.time()),
        "combo": list(record.get("combo", []) or []),
        "source": str((record.get("spec") or {}).get("source", ""))
        if isinstance(record.get("spec"), Mapping)
        else "",
        "run_mode": str(record.get("run_mode") or ""),
        "chunk_tag": record.get("chunk_tag"),
        "output_csv": str(record.get("output_csv") or ""),
        "row_count": int(row_count),
        "repaired_from_output": True,
    }
    with result_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, sort_keys=True)
    return True


def _repair_scorch_shard_result_ledgers(
    dist_dir: Path,
    *,
    run_id: str,
    dry_run: bool,
) -> dict[str, int]:
    counts = {
        "scorch_shard_results_repaired": 0,
        "scorch_invalid_output_files_cleared": 0,
    }
    plans_dir = dist_dir / "scorch_shard_plans"
    results_dir = dist_dir / "scorch_shard_results"
    if not plans_dir.exists():
        return counts

    for plan_path in sorted(plans_dir.glob("*.json")):
        payload = _read_json_maybe(plan_path)
        if not payload:
            continue
        shards = payload.get("shards")
        if not isinstance(shards, list):
            continue
        for record in shards:
            if not isinstance(record, Mapping):
                continue
            shard_id = str(record.get("shard_id") or "").strip()
            if not shard_id:
                continue
            result_path = results_dir / f"{shard_id}.json"
            result = _read_json_maybe(result_path) if result_path.exists() else None
            valid, row_count = _scorch_record_output_valid(record)
            status = str((result or {}).get("status") or "").strip().lower()
            if result is None:
                if valid and _write_repaired_scorch_result(
                    result_path,
                    run_id=run_id,
                    record=record,
                    row_count=row_count,
                    dry_run=dry_run,
                ):
                    counts["scorch_shard_results_repaired"] += 1
                continue
            if status == "completed" and not valid:
                output_raw = str(record.get("output_csv") or "").strip()
                if output_raw and _unlink_if_present(Path(output_raw), dry_run=dry_run):
                    counts["scorch_invalid_output_files_cleared"] += 1
                _unlink_if_present(result_path, dry_run=dry_run)
    return counts


def _unlink_if_present(path: Path, *, dry_run: bool) -> bool:
    if not path.exists():
        return False
    if dry_run:
        return True
    try:
        path.unlink()
        return True
    except Exception:
        return False


def _dedup_paths(paths: Sequence[Path]) -> list[Path]:
    out: list[Path] = []
    seen: set[str] = set()
    for path in paths:
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        out.append(path)
    return out


def _scope_candidate_dirs(root: Path, scope: _Scope) -> list[Path]:
    pdb_root = root / scope.pdb_id
    candidates = [pdb_root]
    if scope.variant != "LEGACY":
        variant_dir = pdb_root / scope.variant
        candidates.append(variant_dir)
        if scope.ph != "base":
            candidates.append(variant_dir / scope.ph)
    else:
        candidates.append(pdb_root / "LEGACY")
        candidates.append(pdb_root / "LEGACY" / "base")
        if scope.ph != "base":
            candidates.append(pdb_root / "LEGACY" / scope.ph)
    if scope.ph != "base":
        candidates.append(pdb_root / scope.ph)
    return _dedup_paths(candidates)


def _target_roots(root: Path, scopes: Sequence[_Scope] | None) -> list[Path]:
    if not scopes:
        return [root]
    return _dedup_paths([root / scope.pdb_id for scope in scopes])


def _consensus_paths_for_scope(docked_root: Path, scope: _Scope) -> list[Path]:
    paths: list[Path] = []
    for candidate_dir in _scope_candidate_dirs(docked_root, scope):
        if not candidate_dir.is_dir():
            continue
        paths.extend(
            path
            for path in sorted(candidate_dir.glob("*consensus_docking_scores.csv"))
            if path.is_file()
        )
    return _dedup_paths(paths)


def _score_value(row: Mapping[str, Any]) -> float:
    try:
        return float(str(row.get("consensus_score") or "").strip())
    except Exception:
        return float("-inf")


def _consensus_needs_sort(rows: Sequence[Mapping[str, Any]]) -> bool:
    previous = float("inf")
    for row in rows:
        current = _score_value(row)
        if current > previous:
            return True
        previous = current
    return False


def _sort_consensus_outputs(
    docked_root: Path,
    *,
    scopes: Sequence[_Scope] | None = None,
    dry_run: bool,
) -> dict[str, int]:
    counts = {
        "consensus_files_checked": 0,
        "consensus_files_sorted": 0,
        "consensus_files_sort_failed": 0,
    }
    if not docked_root.exists():
        return counts
    if scopes:
        paths = _dedup_paths(
            [
                path
                for scope in sorted(scopes, key=lambda item: (item.pdb_id, item.variant, item.ph))
                for path in _consensus_paths_for_scope(docked_root, scope)
            ]
        )
    else:
        paths = sorted(docked_root.rglob("*consensus_docking_scores.csv"))
    for csv_path in paths:
        if not csv_path.is_file():
            continue
        counts["consensus_files_checked"] += 1
        try:
            with csv_path.open("r", newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                fieldnames = list(reader.fieldnames or [])
                if "consensus_score" not in fieldnames:
                    continue
                rows = list(reader)
        except Exception:
            counts["consensus_files_sort_failed"] += 1
            continue
        if not _consensus_needs_sort(rows):
            continue
        counts["consensus_files_sorted"] += 1
        if dry_run:
            continue
        tmp_path = csv_path.with_name(f".{csv_path.name}.repair.tmp")
        try:
            with tmp_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=fieldnames,
                    extrasaction="ignore",
                )
                writer.writeheader()
                writer.writerows(
                    sorted(rows, key=_score_value, reverse=True)
                )
            tmp_path.replace(csv_path)
        except Exception:
            counts["consensus_files_sort_failed"] += 1
            try:
                tmp_path.unlink()
            except Exception:
                pass
    return counts


_SCORCH_POST_FILE_PATTERNS = (
    "*scorch_scores*.csv",
    "*reranked_scorch*.csv",
)


def _clear_scorch_post_outputs(
    post_root: Path,
    *,
    scopes: Sequence[_Scope] | None = None,
    dry_run: bool,
) -> dict[str, int]:
    counts = {
        "scorch_output_files_cleared": 0,
        "scorch_done_sentinels_cleared": 0,
        "scorch_coverage_summaries_cleared": 0,
        "scorch_input_dirs_cleared": 0,
    }
    if not post_root.exists():
        return counts
    if scopes is not None and len(scopes) == 0:
        return counts

    for root in _target_roots(post_root, scopes):
        if not root.exists():
            continue
        for current, dirnames, filenames in os.walk(root, topdown=True):
            current_path = Path(current)
            if ".scorch_inputs" in dirnames:
                input_dir = current_path / ".scorch_inputs"
                counts["scorch_input_dirs_cleared"] += 1
                if not dry_run:
                    try:
                        shutil.rmtree(input_dir)
                    except Exception:
                        pass
                dirnames.remove(".scorch_inputs")
            for filename in filenames:
                path = current_path / filename
                if (
                    filename.endswith(".csv")
                    and (
                        "scorch_scores" in filename
                        or "reranked_scorch" in filename
                    )
                ):
                    if _unlink_if_present(path, dry_run=dry_run):
                        counts["scorch_output_files_cleared"] += 1
                    continue
                if filename == "scorch_coverage_summary.json":
                    if _unlink_if_present(path, dry_run=dry_run):
                        counts["scorch_coverage_summaries_cleared"] += 1
                    continue
                if filename == "_DONE" and current_path.name == "scorch":
                    if _unlink_if_present(path, dry_run=dry_run):
                        counts["scorch_done_sentinels_cleared"] += 1

    return counts


def _count_named_files_scoped(
    root: Path,
    filename: str,
    *,
    scopes: Sequence[_Scope] | None = None,
    prune_dirs: Sequence[str] = (),
) -> int:
    if not root.exists():
        return 0
    count = 0
    pruned = set(prune_dirs)
    for target_root in _target_roots(root, scopes):
        if not target_root.exists():
            continue
        for _current, dirnames, filenames in os.walk(target_root, topdown=True):
            if pruned:
                dirnames[:] = [name for name in dirnames if name not in pruned]
            count += sum(1 for name in filenames if name == filename)
    return count


def _count_files_under(path: Path) -> int:
    if not path.exists():
        return 0
    if path.is_file():
        return 1
    return sum(1 for child in path.rglob("*") if child.is_file() or child.is_symlink())


def _clear_dir_contents(path: Path, *, dry_run: bool) -> int:
    count = _count_files_under(path)
    if dry_run or count <= 0 or not path.exists():
        return count
    try:
        shutil.rmtree(path)
        path.mkdir(parents=True, exist_ok=True)
    except Exception:
        return 0
    return count


def _clear_scorch_selection_state(
    dist_dir: Path,
    *,
    dry_run: bool,
) -> dict[str, int]:
    counts = {
        "scorch_scope_states_cleared": 0,
        "scorch_scope_claims_cleared": 0,
        "scorch_shard_plans_cleared": 0,
        "scorch_shard_results_cleared": 0,
        "scorch_shard_claims_cleared": 0,
        "scorch_shard_summaries_cleared": 0,
    }
    mapping = {
        "scorch_scope_state": "scorch_scope_states_cleared",
        "scorch_scope_claims": "scorch_scope_claims_cleared",
        "scorch_shard_plans": "scorch_shard_plans_cleared",
        "scorch_shard_results": "scorch_shard_results_cleared",
        "scorch_shard_claims": "scorch_shard_claims_cleared",
    }
    for dirname, key in mapping.items():
        counts[key] += _clear_dir_contents(dist_dir / dirname, dry_run=dry_run)
    summary_path = dist_dir / "scorch_shard_summary.json"
    if _unlink_if_present(summary_path, dry_run=dry_run):
        counts["scorch_shard_summaries_cleared"] += 1
    return counts


def _cleanup_scorch_distributed_state(
    *,
    dist_dir: Path,
    docked_root: Path,
    post_root: Path,
    manifest_run_root: Path,
    processed_root: Path,
    require_scorch: bool,
    dry_run: bool,
) -> dict[str, int]:
    counts = {
        "scorch_scope_states_cleared": 0,
        "scorch_scope_claims_cleared": 0,
        "scorch_shard_plans_cleared": 0,
        "scorch_shard_results_cleared": 0,
        "scorch_shard_claims_cleared": 0,
    }

    scope_state_dir = dist_dir / "scorch_scope_state"
    scope_claim_dir = dist_dir / "scorch_scope_claims"
    if scope_state_dir.exists():
        for state_path in sorted(scope_state_dir.glob("*.json")):
            payload = _read_json_maybe(state_path)
            if not payload:
                continue
            status = str(payload.get("status") or "").strip().lower()
            if status == "completed":
                continue
            if status in {"running", ""}:
                if _unlink_if_present(state_path, dry_run=dry_run):
                    counts["scorch_scope_states_cleared"] += 1

    if scope_claim_dir.exists():
        for claim_path in sorted(scope_claim_dir.glob("*.claim.json")):
            if _unlink_if_present(claim_path, dry_run=dry_run):
                counts["scorch_scope_claims_cleared"] += 1

    shard_results_dir = dist_dir / "scorch_shard_results"
    shard_claims_dir = dist_dir / "scorch_shard_claims"
    completed_valid_shards: set[str] = set()
    if shard_results_dir.exists():
        for result_path in sorted(shard_results_dir.glob("*.json")):
            payload = _read_json_maybe(result_path)
            if not payload:
                continue
            shard_id = str(payload.get("shard_id") or result_path.stem)
            if _scorch_result_has_valid_output(payload):
                completed_valid_shards.add(shard_id)
                continue
            status = str(payload.get("status") or "").strip().lower()
            if status in {"running", ""}:
                if _unlink_if_present(result_path, dry_run=dry_run):
                    counts["scorch_shard_results_cleared"] += 1

    if shard_claims_dir.exists():
        for claim_path in sorted(shard_claims_dir.glob("*.claim.json")):
            shard_id = claim_path.name[: -len(".claim.json")]
            if shard_id not in completed_valid_shards:
                if _unlink_if_present(claim_path, dry_run=dry_run):
                    counts["scorch_shard_claims_cleared"] += 1

    return counts


def _summary_from_report(report: Mapping[str, Any]) -> ManifestRepairSummary:
    backup_raw = report.get("backup_path")
    return ManifestRepairSummary(
        run_id=str(report.get("run_id") or ""),
        all_dirs=Path(str(report.get("all_dirs") or ".")),
        manifest_path=Path(str(report.get("manifest_path") or ".")),
        dry_run=bool(report.get("dry_run")),
        require_scorch=bool(report.get("require_scorch")),
        planned_chunks=int(report.get("planned_chunks") or 0),
        completed_chunks=int(report.get("completed_chunks") or 0),
        failed_chunks=int(report.get("failed_chunks") or 0),
        missing_chunks=int(report.get("missing_chunks") or 0),
        planned_ligand_assignments=int(report.get("planned_ligand_assignments") or 0),
        completed_ligand_assignments=int(report.get("completed_ligand_assignments") or 0),
        targets=int(report.get("targets") or 0),
        targets_completed=int(report.get("targets_completed") or 0),
        targets_running=int(report.get("targets_running") or 0),
        targets_failed=int(report.get("targets_failed") or 0),
        total_proteins_scheduled=int(report.get("total_proteins_scheduled") or 0),
        total_protein_list=[str(item) for item in (report.get("total_protein_list") or [])],
        total_proteins_completed=int(report.get("total_proteins_completed") or 0),
        total_proteins_failed=int(report.get("total_proteins_failed") or 0),
        consensus_files=int(report.get("consensus_files") or 0),
        scorch_reranked_files=int(report.get("scorch_reranked_files") or 0),
        repair_scorch_selection=bool(report.get("repair_scorch_selection")),
        consensus_files_checked=int(report.get("consensus_files_checked") or 0),
        consensus_files_sorted=int(report.get("consensus_files_sorted") or 0),
        consensus_files_sort_failed=int(report.get("consensus_files_sort_failed") or 0),
        scorch_output_files_cleared=int(report.get("scorch_output_files_cleared") or 0),
        scorch_done_sentinels_cleared=int(report.get("scorch_done_sentinels_cleared") or 0),
        scorch_coverage_summaries_cleared=int(report.get("scorch_coverage_summaries_cleared") or 0),
        scorch_input_dirs_cleared=int(report.get("scorch_input_dirs_cleared") or 0),
        scorch_shard_summaries_cleared=int(report.get("scorch_shard_summaries_cleared") or 0),
        backup_path=Path(str(backup_raw)) if backup_raw else None,
    )


def _fresh_stage(status: str, details: Mapping[str, Any]) -> dict[str, Any]:
    stage = default_stage_entry()
    stage["status"] = status
    stage["details"] = dict(details)
    return stage


def repair_manifest(
    *,
    run_id: str,
    all_dirs: Path,
    dry_run: bool = False,
    require_scorch: bool | None = None,
    repair_scorch_selection: bool = False,
    clear_incomplete_scorch_outputs: bool = False,
    trust_existing_scorch_coverage: bool = False,
) -> dict[str, Any]:
    manifest_dir = _manifest_dir(all_dirs, run_id)
    manifest_path = manifest_dir / "run_manifest.yaml"
    dist_dir = manifest_dir / "distributed"
    docked_root = _docked_dir(all_dirs, run_id)
    post_root = _post_docked_dir(all_dirs, run_id)
    processed_root = _processed_root(all_dirs)
    if not dist_dir.exists():
        raise RuntimeError(f"Distributed state directory not found: {dist_dir}")

    chunks = _iter_plan_chunks(dist_dir)
    if not chunks:
        raise RuntimeError(f"No distributed chunk plans found under: {dist_dir}")
    results = _load_chunk_results(dist_dir)
    if require_scorch is None:
        require_scorch = _scorch_required_auto(post_root)
    if repair_scorch_selection:
        require_scorch = True

    existing = _read_yaml(manifest_path)
    old_proteins = existing.get("proteins") if isinstance(existing.get("proteins"), dict) else {}
    by_scope: dict[_Scope, list[dict[str, Any]]] = {}
    for chunk in chunks:
        scope = _chunk_scope(chunk)
        if scope is not None:
            by_scope.setdefault(scope, []).append(chunk)
    scopes = sorted(by_scope, key=lambda item: (item.pdb_id, item.variant, item.ph))

    consensus_sort_counts = (
        _sort_consensus_outputs(
            docked_root,
            scopes=scopes,
            dry_run=bool(dry_run),
        )
        if repair_scorch_selection
        else {
            "consensus_files_checked": 0,
            "consensus_files_sorted": 0,
            "consensus_files_sort_failed": 0,
        }
    )
    scorch_selection_post_cleanup = (
        _clear_scorch_post_outputs(
            post_root,
            scopes=scopes,
            dry_run=bool(dry_run),
        )
        if repair_scorch_selection
        else {
            "scorch_output_files_cleared": 0,
            "scorch_done_sentinels_cleared": 0,
            "scorch_coverage_summaries_cleared": 0,
            "scorch_input_dirs_cleared": 0,
        }
    )
    scorch_selection_state_cleanup = (
        _clear_scorch_selection_state(dist_dir, dry_run=bool(dry_run))
        if repair_scorch_selection
        else {
            "scorch_scope_states_cleared": 0,
            "scorch_scope_claims_cleared": 0,
            "scorch_shard_plans_cleared": 0,
            "scorch_shard_results_cleared": 0,
            "scorch_shard_claims_cleared": 0,
            "scorch_shard_summaries_cleared": 0,
        }
    )
    scorch_ledger_repair = (
        _repair_scorch_shard_result_ledgers(
            dist_dir,
            run_id=run_id,
            dry_run=bool(dry_run),
        )
        if bool(require_scorch) and not repair_scorch_selection
        else {
            "scorch_shard_results_repaired": 0,
            "scorch_invalid_output_files_cleared": 0,
        }
    )

    proteins: dict[str, Any] = {}
    repaired_at = _utc_now_iso()
    status_counts = {"completed": 0, "running": 0, "failed": 0}
    missing_chunks_total = 0
    completed_chunks_total = 0
    failed_chunks_total = 0
    planned_ligands_total = 0
    completed_ligands_total = 0
    consensus_count = 0
    scorch_count = 0
    incomplete_scorch_scopes: list[_Scope] = []

    for scope in scopes:
        scope_chunks = by_scope[scope]
        chunk_ids = [_chunk_id(ch) for ch in scope_chunks]
        completed_ids = [
            cid
            for cid in chunk_ids
            if str((results.get(cid) or {}).get("status") or "").strip().lower()
            == "completed"
        ]
        failed_ids = [
            cid
            for cid in chunk_ids
            if str((results.get(cid) or {}).get("status") or "").strip().lower()
            in {"failed", "terminal_failed"}
        ]
        missing_ids = [cid for cid in chunk_ids if cid not in results]
        planned_ligands = sum(_chunk_ligand_count(ch) for ch in scope_chunks)
        chunks_by_id = {_chunk_id(ch): ch for ch in scope_chunks}
        completed_ligands = sum(
            _completed_ligand_count(results.get(cid), chunks_by_id.get(cid, {}))
            for cid in completed_ids
        )
        has_consensus = _has_consensus(docked_root, scope)
        has_scorch = (
            False
            if repair_scorch_selection
            else _has_scorch(
                post_root,
                scope,
                docked_root=docked_root,
                manifest_run_root=manifest_dir,
                processed_root=processed_root,
                require_coverage=bool(require_scorch),
                dry_run=dry_run,
                trust_existing_coverage=bool(trust_existing_scorch_coverage),
            )
        )
        if has_scorch:
            if not dry_run:
                _ensure_scorch_done_sentinel(
                    post_root,
                    scope,
                    docked_root=docked_root,
                    manifest_run_root=manifest_dir,
                    processed_root=processed_root,
                    require_coverage=bool(require_scorch),
                    dry_run=dry_run,
                    trust_existing_coverage=bool(trust_existing_scorch_coverage),
                )
        elif require_scorch:
            incomplete_scorch_scopes.append(scope)
        consensus_count += int(has_consensus)
        scorch_count += int(has_scorch)

        if failed_ids:
            status = "failed"
            reason = "chunk_results_failed"
        elif missing_ids:
            status = "running"
            reason = "chunk_results_missing"
        elif require_scorch and not has_scorch:
            status = "running"
            reason = "missing_or_incomplete_scorch_coverage"
        else:
            status = "completed"
            reason = "complete"

        status_counts[status] += 1
        missing_chunks_total += len(missing_ids)
        completed_chunks_total += len(completed_ids)
        failed_chunks_total += len(failed_ids)
        planned_ligands_total += planned_ligands
        completed_ligands_total += completed_ligands

        key = protein_key(scope.pdb_id, scope.variant, scope.ph)
        entry = dict(old_proteins.get(key) or default_protein_entry())
        entry["pdb_id"] = scope.pdb_id
        entry["variant"] = scope.variant
        entry["ph"] = None if scope.ph == "base" else scope.ph
        libraries = sorted({str(ch.get("library_name") or "").strip() for ch in scope_chunks if str(ch.get("library_name") or "").strip()})
        run_modes = sorted({str(ch.get("run_mode") or "").strip() for ch in scope_chunks if str(ch.get("run_mode") or "").strip()})
        entry["library"] = "+".join(libraries) if libraries else entry.get("library")
        entry["status"] = status
        entry["error"] = None if status != "failed" else reason
        stages = entry.get("stages")
        if not isinstance(stages, MutableMapping):
            stages = {}
        details = {
            "repair_reason": reason,
            "planned_chunks": len(chunk_ids),
            "completed_chunks": len(completed_ids),
            "failed_chunks": len(failed_ids),
            "missing_result_chunks": len(missing_ids),
            "planned_ligands": planned_ligands,
            "completed_ligands": completed_ligands,
            "libraries": libraries,
            "run_modes": run_modes,
            "has_consensus_csv": bool(has_consensus),
            "require_scorch": bool(require_scorch),
            "has_scorch_reranked_csv": bool(has_scorch),
            "repaired_at": repaired_at,
        }
        stages["distributed_chunk_coverage"] = _fresh_stage(status, details)
        docking_status = "failed" if failed_ids else ("completed" if not missing_ids else "running")
        stages["docking"] = _fresh_stage(docking_status, details)
        if require_scorch:
            stages["postprocessing"] = _fresh_stage(
                "completed" if has_scorch else "running",
                {
                    "has_scorch_reranked_csv": bool(has_scorch),
                    "repair_reason": "complete"
                    if has_scorch
                    else "missing_or_incomplete_scorch_coverage",
                    "repaired_at": repaired_at,
                },
            )
        entry["stages"] = dict(stages)
        proteins[key] = entry

    repaired = dict(existing)
    repaired["run_id"] = run_id
    command = repaired.get("command") if isinstance(repaired.get("command"), dict) else {}
    command.setdefault("RUN_ID", run_id)
    command.setdefault("ATLAS_ALL_DIRS", str(all_dirs))
    command.setdefault("DISTRIBUTED_MODE", "slurm_array")
    repaired["command"] = command
    paths = repaired.get("paths") if isinstance(repaired.get("paths"), dict) else {}
    paths.setdefault("all_dirs", str(all_dirs))
    paths.setdefault("manifest_dir", str(manifest_dir))
    paths.setdefault("docked_dir", str(docked_root))
    paths.setdefault("post_docked_dir", str(post_root))
    repaired["paths"] = paths
    repaired["proteins"] = proteins
    repaired["status"] = (
        "failed"
        if status_counts["failed"]
        else ("completed" if status_counts["running"] == 0 else "running")
    )
    repaired["summary"] = {}
    refresh_summary(repaired)
    summary = repaired.get("summary") if isinstance(repaired.get("summary"), dict) else {}
    target_count = int(summary.get("total_proteins_scheduled") or len(proteins))
    targets_completed = int(summary.get("total_proteins_completed") or 0)
    targets_failed = int(summary.get("total_proteins_failed") or 0)
    timing = repaired.get("timing") if isinstance(repaired.get("timing"), dict) else {}
    timing.setdefault("started_at", repaired_at)
    timing["manifest_repaired_at"] = repaired_at
    repaired["timing"] = timing
    if (
        bool(require_scorch)
        and not repair_scorch_selection
        and bool(clear_incomplete_scorch_outputs)
    ):
        stale_post_cleanup = _clear_scorch_post_outputs(
            post_root,
            scopes=incomplete_scorch_scopes,
            dry_run=bool(dry_run),
        )
        for key, value in stale_post_cleanup.items():
            scorch_selection_post_cleanup[key] = int(
                scorch_selection_post_cleanup.get(key) or 0
            ) + int(value or 0)
    if repair_scorch_selection:
        scorch_cleanup = {
            "scorch_scope_states_cleared": 0,
            "scorch_scope_claims_cleared": 0,
            "scorch_shard_plans_cleared": 0,
            "scorch_shard_results_cleared": 0,
            "scorch_shard_claims_cleared": 0,
        }
    elif bool(require_scorch):
        scorch_cleanup = _cleanup_scorch_distributed_state(
            dist_dir=dist_dir,
            docked_root=docked_root,
            post_root=post_root,
            manifest_run_root=manifest_dir,
            processed_root=processed_root,
            require_scorch=bool(require_scorch),
            dry_run=bool(dry_run),
        )
    else:
        scorch_cleanup = {
            "scorch_scope_states_cleared": 0,
            "scorch_scope_claims_cleared": 0,
            "scorch_shard_plans_cleared": 0,
            "scorch_shard_results_cleared": 0,
            "scorch_shard_claims_cleared": 0,
        }
    for key, value in scorch_selection_state_cleanup.items():
        scorch_cleanup[key] = int(scorch_cleanup.get(key) or 0) + int(value or 0)
    repair_report = {
        "run_id": run_id,
        "all_dirs": str(all_dirs),
        "manifest_path": str(manifest_path),
        "dry_run": bool(dry_run),
        "require_scorch": bool(require_scorch),
        "repair_scorch_selection": bool(repair_scorch_selection),
        "planned_chunks": len(chunks),
        "completed_chunks": completed_chunks_total,
        "failed_chunks": failed_chunks_total,
        "missing_chunks": missing_chunks_total,
        "planned_ligand_assignments": planned_ligands_total,
        "completed_ligand_assignments": completed_ligands_total,
        "targets": target_count,
        "targets_completed": targets_completed,
        "targets_running": max(0, target_count - targets_completed - targets_failed),
        "targets_failed": targets_failed,
        "total_proteins_scheduled": target_count,
        "total_protein_list": list(summary.get("total_protein_list") or []),
        "total_proteins_completed": targets_completed,
        "total_proteins_failed": targets_failed,
        "consensus_files": _count_named_files_scoped(
            docked_root,
            "consensus_docking_scores.csv",
            scopes=scopes,
        ),
        "scorch_reranked_files": (
            0
            if repair_scorch_selection and not dry_run
            else _count_named_files_scoped(
                post_root,
                "consensus_reranked_scorch.csv",
                scopes=scopes,
                prune_dirs=(".scorch_inputs",),
            )
        ),
        **consensus_sort_counts,
        **scorch_selection_post_cleanup,
        **scorch_cleanup,
        **scorch_ledger_repair,
        "backup_path": None,
    }
    repaired["repair"] = repair_report | {"repaired_at": repaired_at}

    if not dry_run:
        if manifest_path.exists():
            backup_path = _backup_path(manifest_path)
            shutil.copy2(manifest_path, backup_path)
            repair_report["backup_path"] = str(backup_path)
            repaired["repair"]["backup_path"] = str(backup_path)
        _write_yaml(manifest_path, repaired)

    return repair_report


def repair_run_manifest(
    run_id: str,
    *,
    root: str | Path,
    dry_run: bool = False,
    require_scorch: bool | None = None,
    repair_scorch_selection: bool = False,
    clear_incomplete_scorch_outputs: bool = False,
    trust_existing_scorch_coverage: bool = False,
) -> ManifestRepairSummary:
    report = repair_manifest(
        run_id=str(run_id),
        all_dirs=Path(root).expanduser().resolve(),
        dry_run=bool(dry_run),
        require_scorch=require_scorch,
        repair_scorch_selection=bool(repair_scorch_selection),
        clear_incomplete_scorch_outputs=bool(clear_incomplete_scorch_outputs),
        trust_existing_scorch_coverage=bool(trust_existing_scorch_coverage),
    )
    return _summary_from_report(report)


def _print_repair_report(report: Mapping[str, Any]) -> None:
    fields = [
        "run_id",
        "manifest_path",
        "dry_run",
        "require_scorch",
        "repair_scorch_selection",
        "planned_chunks",
        "completed_chunks",
        "missing_chunks",
        "failed_chunks",
        "planned_ligand_assignments",
        "completed_ligand_assignments",
        "targets",
        "targets_completed",
        "targets_running",
        "targets_failed",
        "total_proteins_scheduled",
        "total_proteins_completed",
        "total_proteins_failed",
        "consensus_files",
        "scorch_reranked_files",
        "consensus_files_checked",
        "consensus_files_sorted",
        "consensus_files_sort_failed",
        "scorch_output_files_cleared",
        "scorch_done_sentinels_cleared",
        "scorch_coverage_summaries_cleared",
        "scorch_input_dirs_cleared",
        "scorch_scope_states_cleared",
        "scorch_scope_claims_cleared",
        "scorch_shard_plans_cleared",
        "scorch_shard_results_cleared",
        "scorch_shard_claims_cleared",
        "scorch_shard_summaries_cleared",
        "scorch_shard_results_repaired",
        "scorch_invalid_output_files_cleared",
        "backup_path",
    ]
    for key in fields:
        print(f"{key}={report.get(key)}")
    protein_list = ",".join(str(item) for item in (report.get("total_protein_list") or []))
    print(f"total_protein_list={protein_list}")
    print(
        "chunks planned={planned_chunks} completed={completed_chunks} missing={missing_chunks} failed={failed_chunks}".format(
            **report
        )
    )
    print(
        "outputs consensus_files={consensus_files} scorch_reranked_files={scorch_reranked_files}".format(
            **report
        )
    )


def add_repair_manifest_parser(subparsers: argparse._SubParsersAction[Any]) -> None:
    parser = subparsers.add_parser(
        "repair-manifest",
        help="Repair run_manifest.yaml from distributed chunk plan/results.",
    )
    parser.add_argument("run_id")
    parser.add_argument(
        "--all-dirs",
        default="",
        help="Relocated run root passed to tools/run_relocated_mode.py --all-dirs.",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--repair-scorch-selection",
        action="store_true",
        help=(
            "Sort consensus CSVs by consensus_score and clear stale SCORCH "
            "outputs/state so resume recomputes score-based SCORCH selection. "
            "Implies --require-scorch."
        ),
    )
    parser.add_argument(
        "--preserve-incomplete-scorch-outputs",
        action="store_true",
        help=(
            "Repair manifest state without deleting partial SCORCH CSVs for "
            "targets that fail strict SCORCH coverage."
        ),
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--require-scorch",
        action="store_true",
        help="Only mark a target completed when reranked SCORCH output exists.",
    )
    group.add_argument(
        "--no-require-scorch",
        action="store_true",
        help="Ignore SCORCH outputs when computing target completion.",
    )


def cmd_repair_manifest(args: argparse.Namespace) -> int:
    try:
        require_scorch: bool | None = None
        if bool(getattr(args, "require_scorch", False)):
            require_scorch = True
        elif bool(getattr(args, "no_require_scorch", False)):
            require_scorch = False
        report = repair_manifest(
            run_id=str(args.run_id),
            all_dirs=_resolve_all_dirs(str(args.all_dirs)),
            dry_run=bool(getattr(args, "dry_run", False)),
            require_scorch=require_scorch,
            repair_scorch_selection=bool(
                getattr(args, "repair_scorch_selection", False)
            ),
            clear_incomplete_scorch_outputs=bool(require_scorch)
            and not bool(getattr(args, "preserve_incomplete_scorch_outputs", False)),
        )
    except Exception as exc:
        print(f"atlas slurm repair-manifest: {exc}", flush=True)
        return 1
    _print_repair_report(report)
    return 0


def _cmd_repair_manifest(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(
        prog="atlas slurm repair-manifest",
        description="Repair run_manifest.yaml from distributed chunk plan/results.",
    )
    parser.add_argument("run_id")
    parser.add_argument("--all-dirs", "--root", dest="all_dirs", default="")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--repair-scorch-selection", action="store_true")
    parser.add_argument(
        "--preserve-incomplete-scorch-outputs",
        action="store_true",
        help=(
            "Repair manifest state without deleting partial SCORCH CSVs for "
            "incomplete targets. Useful before resumable SCORCH recovery."
        ),
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--require-scorch", action="store_true")
    group.add_argument("--no-require-scorch", action="store_true")
    args = parser.parse_args(list(argv))
    rc = cmd_repair_manifest(args)
    return rc


__all__ = [
    "add_repair_manifest_parser",
    "_cmd_repair_manifest",
    "cmd_repair_manifest",
    "ManifestRepairSummary",
    "repair_manifest",
    "repair_run_manifest",
]
