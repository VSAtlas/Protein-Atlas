from __future__ import annotations

import logging
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from post_docking.rescoring import scorch_housekeeping as _scorch_housekeeping_mod
from post_docking.rescoring import scorch_postprocess as _scorch_postprocess_mod
from post_docking.rescoring import scorch_prefix_and_stages as _scorch_prefix_mod
from post_docking.rescoring import scorch_provisional_cache as _scorch_cache_mod
from post_docking.rescoring import scorch_selection as _scorch_selection_mod
from post_docking.rescoring import scorch_stage_io as _scorch_stage_io_mod
from post_docking.rescoring import scorch_stage_runner as _scorch_stage_runner_mod
from post_docking.rescoring import rescoring_scorch_support as _scorch_support_mod
from post_docking.rescoring.scorch_events import emit_task_event
from post_docking.rescoring import scorch_runtime_state as _scorch_runtime_mod
from post_docking.rescoring.scorch_runtime_state import (
    COMPONENT,
    SCORCH_DONE_DIR,
    SCORCH_DONE_SENTINEL,
)
from post_docking.rescoring.scorch_types import (
    PoseSelectionResult,
    SelectionResult,
    StageSpec,
)

_materialize_inputs = _scorch_stage_io_mod.materialize_inputs
_annotate_csv = _scorch_stage_io_mod.annotate_csv
_decoy_prefix_value = _scorch_support_mod._decoy_prefix_value
_decoy_engine_stage_dirs_all = _scorch_support_mod._decoy_engine_stage_dirs_all


def _set_decoy_prefix(value: str, logger: Optional[logging.Logger] = None) -> str:
    return _scorch_support_mod._set_decoy_prefix(value, logger)


def _decoy_prefixes_from_test_mode(cfg, override):
    return _scorch_support_mod._decoy_prefixes_from_test_mode(cfg, override)


def _parse_test_mode_tokens(cfg):
    return _scorch_support_mod._parse_test_mode_tokens(cfg)


def _resolve_decoy_prefix_override(args):
    return _scorch_support_mod._resolve_decoy_prefix_override(args)


def _decoy_stage_dirs(order: Sequence[int]) -> Tuple[str, ...]:
    return _scorch_support_mod._decoy_stage_dirs(order)


def _decoy_engine_stage_dirs(engine: str, order: Sequence[int]) -> Tuple[str, ...]:
    return _scorch_support_mod._decoy_engine_stage_dirs(engine, order)


def _decoy_engine_stage_dirs_legacy(engine: str, order: Sequence[int]) -> Tuple[str, ...]:
    return _scorch_support_mod._decoy_engine_stage_dirs_legacy(engine, order)


def _decoy_post_stage_dirs() -> Tuple[str, ...]:
    return _scorch_support_mod._decoy_post_stage_dirs()


def _vina_stage_dirs() -> Tuple[str, ...]:
    return _scorch_support_mod._vina_stage_dirs()


def _gnina_stage_dirs() -> Tuple[str, ...]:
    return _scorch_support_mod._gnina_stage_dirs()


def _dock6_stage_dirs() -> Tuple[str, ...]:
    return _scorch_support_mod._dock6_stage_dirs()


def _ledock_stage_dirs() -> Tuple[str, ...]:
    return _scorch_support_mod._ledock_stage_dirs()


def _post_stage_dirs() -> Tuple[str, ...]:
    return _scorch_support_mod._post_stage_dirs()


def stage_dir_candidates(
    source: str, mode: str, root: Optional[Path] = None
) -> Tuple[str, ...]:
    return _scorch_support_mod.stage_dir_candidates(source, mode, root)


def _discover_mode_dirs(combo, run_root, post_root, specs, logger):
    return _scorch_support_mod._discover_mode_dirs(
        combo,
        run_root,
        post_root,
        specs,
        logger,
    )


def _collect_combo_from_rel(parts: Sequence[str]) -> Optional[Tuple[str, str, str]]:
    return _scorch_prefix_mod._collect_combo_from_rel(parts)


def discover_combos(
    run_root: Path,
    post_root: Path,
    *,
    pdb_id: Optional[str] = None,
    variant: Optional[str] = None,
    ph: Optional[str] = None,
) -> Set[Tuple[str, str, str]]:
    return _scorch_prefix_mod.discover_combos(
        run_root,
        post_root,
        decoy_prefix=_decoy_prefix_value(),
        pdb_id=pdb_id,
        variant=variant,
        ph=ph,
    )


def _normalized_filter_token(value: Optional[str]) -> Optional[str]:
    return _scorch_prefix_mod.normalized_filter_token(value)


def _combo_matches_filters(
    combo: Tuple[str, str, str],
    *,
    pdb_id: Optional[str] = None,
    variant: Optional[str] = None,
    ph: Optional[str] = None,
) -> bool:
    return _scorch_prefix_mod.combo_matches_filters(
        combo,
        pdb_id=pdb_id,
        variant=variant,
        ph=ph,
    )


def _filter_combos(
    combos: Set[Tuple[str, str, str]],
    *,
    pdb_id: Optional[str] = None,
    variant: Optional[str] = None,
    ph: Optional[str] = None,
) -> Set[Tuple[str, str, str]]:
    return _scorch_prefix_mod.filter_combos(
        combos,
        pdb_id=pdb_id,
        variant=variant,
        ph=ph,
    )


def _done_sentinel_path(post_run_root: Path, combo: Tuple[str, str, str]) -> Path:
    return _scorch_prefix_mod.done_sentinel_path(
        post_run_root,
        combo,
        done_dir=SCORCH_DONE_DIR,
        done_sentinel=SCORCH_DONE_SENTINEL,
    )


def _filter_done_combos(
    combos: Set[Tuple[str, str, str]],
    post_run_root: Path,
    overwrite: bool,
    logger: logging.Logger,
) -> Set[Tuple[str, str, str]]:
    pending, _stats = _filter_done_combos_with_stats(
        combos,
        post_run_root,
        overwrite,
        logger,
    )
    return pending


def _filter_done_combos_with_stats(
    combos: Set[Tuple[str, str, str]],
    post_run_root: Path,
    overwrite: bool,
    logger: logging.Logger,
) -> Tuple[Set[Tuple[str, str, str]], Dict[str, int]]:
    return _scorch_prefix_mod.filter_done_combos_with_stats(
        combos,
        post_run_root,
        overwrite,
        logger,
        component=COMPONENT,
        done_dir=SCORCH_DONE_DIR,
        done_sentinel=SCORCH_DONE_SENTINEL,
    )


def _log_no_combos_after_filters(
    logger: logging.Logger,
    *,
    run_id: str,
    decoy_prefix: str,
    pdb_id_filter: str,
    variant_filter: str,
    ph_filter: str,
    idempotent_only: bool,
    combos_input: int = 0,
    combos_skipped: int = 0,
    filtered_scope: bool = False,
) -> None:
    _scorch_prefix_mod.log_no_combos_after_filters(
        logger,
        component=COMPONENT,
        run_id=run_id,
        decoy_prefix=decoy_prefix,
        pdb_id_filter=pdb_id_filter,
        variant_filter=variant_filter,
        ph_filter=ph_filter,
        idempotent_only=idempotent_only,
        combos_input=combos_input,
        combos_skipped=combos_skipped,
        filtered_scope=filtered_scope,
    )


def _mark_combo_done(
    post_run_root: Path, combo: Tuple[str, str, str], logger: logging.Logger
) -> None:
    _scorch_prefix_mod.mark_combo_done(
        post_run_root,
        combo,
        logger,
        component=COMPONENT,
        done_dir=SCORCH_DONE_DIR,
        done_sentinel=SCORCH_DONE_SENTINEL,
    )


def _collect_stage_pdbqts(ph_root: Path, stage_dir: str) -> List[Path]:
    return _scorch_selection_mod.collect_stage_pdbqts(ph_root, stage_dir)


def _pose_base_from_path(p: Path) -> str:
    return _scorch_selection_mod.pose_base_from_path(
        p,
        decoy_prefix=_decoy_prefix_value(),
    )


def _control_base_from_path(p: Path) -> str:
    return _scorch_selection_mod.control_base_from_path(p)


def _load_control_bases(
    processed_root: Path, pdb_id: str, logger: logging.Logger
) -> Set[str]:
    return _scorch_selection_mod.load_control_bases(
        processed_root,
        pdb_id,
        logger,
        component=COMPONENT,
    )


def _stage_priority(stage_dir: str, pdbqt: Optional[Path] = None) -> int:
    return _scorch_selection_mod.stage_priority(stage_dir, pdbqt)


def _collect_best_pose_per_base(
    ph_root: Path,
    stage_dirs: Sequence[str],
    allowed_bases: Optional[Set[str]],
    logger: logging.Logger,
    preferred_stage_by_base: Optional[Dict[str, str]] = None,
) -> PoseSelectionResult:
    return _scorch_selection_mod.collect_best_pose_per_base(
        ph_root,
        stage_dirs,
        allowed_bases,
        logger,
        component=COMPONENT,
        decoy_prefix=_decoy_prefix_value(),
        preferred_stage_by_base=preferred_stage_by_base,
    )


def _load_consensus_top_bases(
    consensus_csv: Path,
    frac: float,
    logger: logging.Logger,
    control_bases: Optional[Set[str]] = None,
) -> Tuple[Set[str], Set[str], int, int, int, int, int]:
    return _scorch_selection_mod.load_consensus_top_bases(
        consensus_csv,
        frac,
        logger,
        component=COMPONENT,
        control_bases=control_bases,
        decoy_prefix=_decoy_prefix_value(),
    )


def _select_top_bases_from_score_csv(
    score_csv: Path,
    frac: float,
    control_bases: Optional[Set[str]],
    *,
    higher_is_better: bool,
    logger: logging.Logger,
    score_cols: Optional[Sequence[str]] = None,
) -> SelectionResult:
    return _scorch_selection_mod.select_top_bases_from_score_csv(
        score_csv,
        frac,
        control_bases,
        higher_is_better=higher_is_better,
        logger=logger,
        component=COMPONENT,
        score_cols=score_cols,
        decoy_prefix=_decoy_prefix_value(),
    )


def _score_csv_for_spec(
    run_root: Path,
    combo: Tuple[str, str, str],
    spec: StageSpec,
    run_mode: str,
) -> Tuple[Optional[Path], List[str], bool]:
    return _scorch_selection_mod.score_csv_for_spec(
        run_root,
        combo,
        spec,
        run_mode,
        decoy_prefix=_decoy_prefix_value(),
    )


def _log_score_csv_coverage(
    source: str,
    combo: Tuple[str, str, str],
    score_csv: Path,
    allowed_bases: Set[str],
    logger: logging.Logger,
) -> None:
    _scorch_selection_mod.log_score_csv_coverage(
        source,
        combo,
        score_csv,
        allowed_bases,
        logger,
        component=COMPONENT,
        decoy_prefix=_decoy_prefix_value(),
    )


def _posebusters_missing(post_root: Path, combos: Set[Tuple[str, str, str]]) -> bool:
    return _scorch_housekeeping_mod.posebusters_missing(post_root, combos)


def _prep_missing(
    run_root: Path, post_root: Path, combos: Set[Tuple[str, str, str]]
) -> bool:
    return _scorch_housekeeping_mod.prep_missing(
        run_root,
        post_root,
        combos,
        decoy_prefix=_decoy_prefix_value(),
        decoy_engine_stage_dirs_all=_decoy_engine_stage_dirs_all,
    )


def _run_pose_bust(
    run_id: str,
    repo_root: Path,
    docked_root: Path,
    post_docked_root: Path,
    overwrite: bool,
    max_workers: int,
    logger: logging.Logger,
) -> bool:
    return _scorch_housekeeping_mod.run_pose_bust(
        run_id,
        repo_root,
        docked_root,
        post_docked_root,
        overwrite,
        max_workers,
        logger,
        component=COMPONENT,
        sys_executable=str(sys.executable),
        run_subprocess=subprocess.run,
    )


def _launch_pose_bust_async(
    run_id: str,
    repo_root: Path,
    docked_root: Path,
    post_docked_root: Path,
    overwrite: bool,
    max_workers: int,
    logger: logging.Logger,
) -> int:
    return _scorch_housekeeping_mod.launch_pose_bust_async(
        run_id,
        repo_root,
        docked_root,
        post_docked_root,
        overwrite,
        max_workers,
        logger,
        component=COMPONENT,
        sys_executable=str(sys.executable),
        popen_subprocess=subprocess.Popen,
    )


def _run_prep_for_scorch(
    run_id: str,
    repo_root: Path,
    docked_root: Path,
    post_docked_root: Path,
    overwrite: bool,
    logger: logging.Logger,
    decoy_prefix: Optional[str] = None,
    pdb_id: Optional[str] = None,
    variant: Optional[str] = None,
    ph: Optional[str] = None,
    max_workers: int = 1,
) -> bool:
    return _scorch_housekeeping_mod.run_prep_for_scorch(
        run_id,
        repo_root,
        docked_root,
        post_docked_root,
        overwrite,
        logger,
        component=COMPONENT,
        sys_executable=str(sys.executable),
        run_subprocess=subprocess.run,
        decoy_prefix=decoy_prefix,
        pdb_id=pdb_id,
        variant=variant,
        ph=ph,
        max_workers=max_workers,
    )


def _scorch_command(
    receptor: Path,
    ligands: Path,
    threads: int,
    cfg: Mapping[str, object] | None = None,
) -> List[str]:
    cfg_env_prefix = None
    cfg_env_name = None
    if cfg is not None:
        raw_prefix = cfg.get("SCORCH_ENV_PREFIX")
        raw_env = cfg.get("SCORCH_ENV")
        if raw_prefix and str(raw_prefix).strip():
            cfg_env_prefix = Path(str(raw_prefix).strip()).expanduser().resolve()
        if raw_env and str(raw_env).strip():
            cfg_env_name = str(raw_env).strip()
    env_prefix = (
        cfg_env_prefix
        if cfg_env_prefix is not None
        else _scorch_runtime_mod.SCORCH_ENV_PREFIX
    )
    env_name = cfg_env_name if cfg_env_name is not None else _scorch_runtime_mod.SCORCH_ENV
    scorch_script = _scorch_runtime_mod.SCORCH_SCRIPT
    if scorch_script is None:
        raw_script = (cfg or {}).get("SCORCH") or (cfg or {}).get("SCORCH_SCRIPT")
        if raw_script and str(raw_script).strip():
            scorch_script = Path(str(raw_script)).expanduser().resolve()
    if scorch_script is None:
        raise RuntimeError("SCORCH script is not initialized; run _preflight first")
    if _scorch_runtime_mod.SCORCH_USE_MICROMAMBA:
        cmd = ["micromamba", "run"]
        if env_prefix is not None:
            cmd.extend(["-p", str(env_prefix)])
        else:
            cmd.extend(["-n", env_name])
        cmd.extend(["python", str(scorch_script)])
    else:
        cmd = [
            str(_scorch_runtime_mod.SCORCH_PYTHON or Path(sys.executable).resolve()),
            str(scorch_script),
        ]
    cmd.extend(
        [
            "--receptor",
            str(receptor),
            "--ligand",
            str(ligands),
            "--out",
            "{out}",
            "--threads",
            str(threads),
            "--verbose",
        ]
    )
    return cmd


def _score_stage(
    cfg: Dict[str, object],
    spec: StageSpec,
    combo: Tuple[str, str, str],
    run_root: Path,
    post_root: Path,
    receptor: Path,
    threads: int,
    overwrite: bool,
    logger: logging.Logger,
    allowed_bases: Optional[Set[str]] = None,
    control_bases: Optional[Set[str]] = None,
    run_mode: str = "fda",
    stage_dirs_override: Optional[Sequence[str]] = None,
    score_csv: Optional[Path] = None,
    selected_stage_by_base: Optional[Dict[str, str]] = None,
    selected_score_by_base: Optional[Dict[str, float]] = None,
    chunk_tag: Optional[str] = None,
    hedge_group_id: Optional[str] = None,
    hedge_role: str = "primary",
    rechunk_generation: int = 0,
) -> Tuple[bool, Optional[Path]]:
    return _scorch_stage_runner_mod.score_stage(
        cfg,
        spec,
        combo,
        run_root,
        post_root,
        receptor,
        threads,
        overwrite,
        logger,
        decoy_prefix=_decoy_prefix_value(),
        component=COMPONENT,
        emit_task_event=emit_task_event,
        materialize_inputs=_materialize_inputs,
        scorch_command=_scorch_command,
        annotate_csv=_annotate_csv,
        scorch_root=_scorch_runtime_mod.SCORCH_ROOT,
        scorch_python=_scorch_runtime_mod.SCORCH_PYTHON,
        scorch_script=_scorch_runtime_mod.SCORCH_SCRIPT,
        allowed_bases=allowed_bases,
        control_bases=control_bases,
        run_mode=run_mode,
        stage_dirs_override=stage_dirs_override,
        score_csv=score_csv,
        selected_stage_by_base=selected_stage_by_base,
        selected_score_by_base=selected_score_by_base,
        chunk_tag=chunk_tag,
        hedge_group_id=hedge_group_id,
        hedge_role=hedge_role,
        rechunk_generation=rechunk_generation,
        artifact_post_root=(
            _scorch_cache_mod.artifact_post_root(cfg, run_root.name) / run_root.name
            if _scorch_cache_mod.provisional_run_enabled(cfg)
            else None
        ),
        cache_read=(
            _scorch_cache_mod.final_reuse_enabled(cfg)
            and not _scorch_cache_mod.cache_write_enabled(cfg)
            and not overwrite
        ),
        cache_write=_scorch_cache_mod.cache_write_enabled(cfg),
        cache_phase=(
            "provisional"
            if _scorch_cache_mod.provisional_run_enabled(cfg)
            else "final"
        ),
        run_id=run_root.name,
    )


def _aggregate_combo(
    post_root: Path,
    specs: List[StageSpec],
    combo: Tuple[str, str, str],
    logger: logging.Logger,
    *,
    run_mode: str = "fda",
    output_name: str = "scorch_scores_all.csv",
    input_csvs: Optional[Sequence[Path | str]] = None,
    allowed_ligand_bases: Optional[Sequence[str] | set[str]] = None,
) -> Optional[Path]:
    return _scorch_postprocess_mod.aggregate_combo(
        post_root,
        specs,
        combo,
        logger,
        component=COMPONENT,
        run_mode=run_mode,
        output_name=output_name,
        decoy_prefix=_decoy_prefix_value(),
        input_csvs=input_csvs,
        allowed_ligand_bases=allowed_ligand_bases,
    )


def _compute_best_composites(rows: List[Dict[str, str]]) -> Dict[str, float]:
    return _scorch_postprocess_mod.compute_best_composites(
        rows,
        pose_base_from_path=_pose_base_from_path,
    )


def _annotate_scorch_file(
    csv_path: Path,
    best_map: Dict[str, float],
    mu_decoy: Optional[float],
    sigma_decoy: Optional[float],
    n_decoys: int,
) -> None:
    _scorch_postprocess_mod.annotate_scorch_file(
        csv_path,
        best_map,
        mu_decoy,
        sigma_decoy,
        n_decoys,
        pose_base_from_path=_pose_base_from_path,
    )


def annotate_scorch_z_scores(
    fda_csv: Path, dud_csv: Path, logger: logging.Logger
) -> None:
    _scorch_postprocess_mod.annotate_scorch_z_scores(
        fda_csv,
        dud_csv,
        logger,
        pose_base_from_path=_pose_base_from_path,
    )


def annotate_scorch_t_scores(
    fda_csv: Path, dud_csv: Path, logger: logging.Logger
) -> None:
    _scorch_postprocess_mod.annotate_scorch_t_scores(
        fda_csv,
        dud_csv,
        logger,
        pose_base_from_path=_pose_base_from_path,
    )
