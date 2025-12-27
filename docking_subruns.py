from __future__ import annotations

import logging
import math
import os
import re
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from apo_holo_mode import _record_apo_holo_usage
from checkpoints import (
    checkpoint_invalidate_from,
    checkpoint_mark_done,
    checkpoint_should_skip,
)
from docking_centering import CenterSelector
from docking_utils import (
    _fingerprint_stage,
    _write_audit_json,
    early_recenter_decision,
    final_pose_validation_and_screenshots,
    run_completion_audit,
    norm,
)
from fallback_recenter import (
    GlobalCenterGuard,
    RecenterParams,
    fallback_recentering_if_empty,
)
from input_and_export_functions import (
    annotate_fda_long_csv_with_t_scores_vs_decoys,
    record_score,
    _to_bool,
)
from path_router import Paths, docked_dir, receptor_file, load_ph_tags
from ph_ensemble_docking import (
    enumerate_ligands_for_ph_context,
    init_ph_tags_and_manifest,
    prewarm_ph_ligand_microstates,
)
import prep_dock6
from run_manifest import (
    update_manifest_for_docking_stage,
    update_manifest_for_protein_failure,
    update_manifest_for_protein_start,
    update_manifest_for_protein_success,
    update_manifest_for_druggability_and_engine_plan,
)
from druggability_orchestrator import EnginePolicy, decide_engine_policy
from record_data import compute_ligand_efficiency, record_le
from docking_gnina import (
    run_gnina_for_stage,
    write_gnina_scores_csv,
    annotate_gnina_fda_long_csv_with_t_scores_vs_decoys,
)
from docking_ledock import (
    run_ledock_for_stage,
    should_run_ledock_for_target,
    write_ledock_scores_csv,
    annotate_ledock_fda_long_csv_with_t_scores_vs_decoys,
)
from docking_dock6 import (
    run_dock6_for_stage,
    should_run_dock6_for_target,
    write_dock6_scores_csv,
)
from docking_consensus_score import compute_consensus_for_variant_ph
from docking_vina import emit_vina_config, write_scores_csv
from docking_ligands import (
    compute_stage_membership_from_scores,
    _count_heavy_atoms_from_pdbqt,
    _lib_roots_for_pdb,
    _resolve_test_mode,
    prepare_and_filter_ligands,
    select_ligands_for_next,
)
from run_vina import run_docking_task
from single_ligand_index import (
    _ensure_single_ligand_index,
    _load_fda_name_map,
    _resolve_single_ligand,
)


def _apply_force_carry_and_doping(
    cfg: Dict[str, Any],
    docking_mode: str,
    stage_index: int,
    stages_for_run: List[Dict[str, Any]],
    scores: Dict[str, float],
    logger: logging.Logger,
    *,
    stage1_original: List[str],
    forced_extracted_for_stage3: Optional[set] = None,
    invalids: Optional[Dict[str, Tuple[Optional[float], str]]] = None,
    next_stage_is_stage3: bool = False,
    control_norms: Optional[set[str]] = None,
) -> List[str]:
    """
    Mirror the existing Vina selection + doping flow:
      - percentile selection via select_ligands_for_next
      - optional force-carry of extracted controls into stage3
        (filtered back out here; controls are injected at Stage3 run time)
      - optional rescue of self-RMSD near-miss ligands
    """
    use_stage1_base = (docking_mode == "polypharmacology" and stage_index == 1)

    rescue: List[str] = []
    if invalids and stage_index < len(stages_for_run) - 1:
        for lig, (sc, reason) in invalids.items():
            if sc is not None and "self_rmsd_" in str(reason).lower() and sc <= float(
                cfg.get("RESCUE_SELF_RMSD_SCORE_MAX", -8.0)
            ):
                rescue.append((sc, lig))
        rescue = [lig for _, lig in sorted(rescue)[: int(cfg.get("RESCUE_SELF_RMSD_TOP_N", 10))]]

    selected_raw = select_ligands_for_next(
        docking_mode,
        stage_index,
        stages_for_run,
        scores,
        logger,
        base_pool_n=(len(stage1_original) if use_stage1_base else None),
        force_include=(forced_extracted_for_stage3 if next_stage_is_stage3 else None),
    )

    def _filter_controls(seq: List[str]) -> List[str]:
        if not control_norms:
            return list(seq)
        filt = []
        for lig in seq:
            if norm(lig) in control_norms:
                continue
            filt.append(lig)
        return filt

    selected = _filter_controls(selected_raw)
    removed_selected = len(selected_raw) - len(selected)
    if removed_selected and next_stage_is_stage3:
        logger.info(
            "[Force-carry] Dropped %d control ligand(s) from selection; controls are merged at Stage3 docking time.",
            removed_selected,
        )

    if rescue:
        sel_set = set(selected)
        rescue_unique = [r for r in _filter_controls(rescue) if r not in sel_set]
        ligands = rescue_unique + selected
    else:
        ligands = selected
    return ligands


def _use_ledock(cfg: Dict[str, Any]) -> bool:
    try:
        return should_run_ledock_for_target(cfg)
    except Exception:
        return False


def _use_dock6(cfg: Dict[str, Any]) -> bool:
    try:
        return should_run_dock6_for_target(cfg)
    except Exception:
        return False


def _is_stage3(stage_name: str) -> bool:
    s = str(stage_name).lower()
    return re.search(r"(?:^|_)stage3(?:$|_)", s) is not None


def _split_controls_and_noncontrols(
    ligands: List[str],
    prepped_root: Optional[Path],
    *,
    allowed_subdirs: Optional[set[str]] = None,
    control_stems: Optional[set[str]] = None,
) -> tuple[list[str], list[str]]:
    controls: list[str] = []
    noncontrols: list[str] = []
    allow = allowed_subdirs or {"controls", "reference"}
    ctrl_stems_lower = {s.lower() for s in (control_stems or set())}

    def _is_control_path(p: Path) -> bool:
        if not prepped_root:
            return False
        try:
            rel = p.resolve().relative_to(prepped_root.resolve())
        except Exception:
            return False
        if not rel.parts:
            return False
        return rel.parts[0] in allow

    def _is_control_stem(p: Path) -> bool:
        stem = p.stem.split("_stage")[0].lower()
        return stem in ctrl_stems_lower

    for lig in ligands:
        p = Path(lig)
        if _is_control_path(p) or _is_control_stem(p):
            controls.append(str(p))
        else:
            noncontrols.append(str(p))

    return controls, noncontrols


def _interleave_controls(noncontrols: List[str], controls: List[str]) -> List[str]:
    merged: list[str] = []
    i = 0
    max_len = max(len(noncontrols), len(controls))
    while i < max_len:
        if i < len(noncontrols):
            merged.append(noncontrols[i])
        if i < len(controls):
            merged.append(controls[i])
        i += 1
    return merged


@dataclass
class SubrunSpec:
    run_mode: Optional[str]
    csv_prefix: str
    stage_name_prefix: str


@dataclass
class ProteinDockingContext:
    cfg: Dict[str, Any]
    paths: Paths
    logger: logging.Logger
    pdb_id: str

    variant_env: str
    variant_token: Optional[str]
    variant_label: str
    legacy_mode: bool

    cleaned_pdb: Optional[str]
    receptor_pdbqt: Optional[str]
    center: Optional[Tuple[float, float, float]]
    box_size: Optional[Tuple[float, float, float]]

    stages: List[Dict[str, Any]]
    recenter_params: RecenterParams

    control_stems: List[str]
    control_lookup: Dict[str, Path]


def subruns_for_test_mode(test_mode: str) -> List[SubrunSpec]:
    """
    Return a list of subrun descriptors with keys:
      - run_mode: None | "dud" | "fda" | "hmdb"
      - csv_prefix: str
      - stage_name_prefix: str
    """
    if test_mode in (None, "", "off", "dud"):
        return [
            SubrunSpec(
                run_mode=None,
                csv_prefix="",
                stage_name_prefix="",
            )
        ]
    if test_mode == "fda":
        return [
            SubrunSpec(
                run_mode="fda",
                csv_prefix="",
                stage_name_prefix="",
            )
        ]
    if test_mode == "fda+dud":
        return [
            SubrunSpec(run_mode="dud", csv_prefix="dud_", stage_name_prefix="dud_"),
            SubrunSpec(run_mode="fda", csv_prefix="", stage_name_prefix=""),
        ]
    if test_mode == "hmdb":
        return [
            SubrunSpec(run_mode="hmdb", csv_prefix="hmdb_", stage_name_prefix="hmdb_"),
        ]
    if test_mode == "hmdb+dud":
        return [
            SubrunSpec(run_mode="hmdb", csv_prefix="hmdb_", stage_name_prefix="hmdb_"),
            SubrunSpec(run_mode="dud", csv_prefix="dud_", stage_name_prefix="dud_"),
        ]
    if test_mode == "hmdb+fda":
        return [
            SubrunSpec(run_mode="hmdb", csv_prefix="hmdb_", stage_name_prefix="hmdb_"),
            SubrunSpec(run_mode="fda", csv_prefix="", stage_name_prefix=""),
        ]
    if test_mode == "fda+dud+hmdb":
        return [
            SubrunSpec(run_mode="dud", csv_prefix="dud_", stage_name_prefix="dud_"),
            SubrunSpec(run_mode="hmdb", csv_prefix="hmdb_", stage_name_prefix="hmdb_"),
            SubrunSpec(run_mode="fda", csv_prefix="", stage_name_prefix=""),
        ]
    return [
        SubrunSpec(
            run_mode=None,
            csv_prefix="",
            stage_name_prefix="",
        )
    ]


def run_ligand_pipeline_subrun(ctx: ProteinDockingContext, subrun: SubrunSpec) -> None:
    """
    Run the existing ligands + multi-stage docking pipeline once,
    but parameterized by:
      - run_mode: None | "dud" | "fda" | "hmdb"
      - csv_prefix: "" or "dud_" or "hmdb_"
      - stage_name_prefix: "" or "dud_" or "hmdb_"
    Variant and pH behavior must remain unchanged: only the final stage
    component gets the prefix.
    """
    from docking import RetryManager, run_one_stage  # late import to avoid circular dependency
    from prep_for_ledock import ensure_mol2_for_ledock

    cfg = ctx.cfg
    paths = ctx.paths
    logger = ctx.logger
    pdb_id = ctx.pdb_id
    variant_env = ctx.variant_env
    variant_token = ctx.variant_token
    variant_label = ctx.variant_label
    legacy_mode = ctx.legacy_mode
    cleaned_pdb = ctx.cleaned_pdb
    center = ctx.center
    box_size = ctx.box_size
    stages = ctx.stages
    params = ctx.recenter_params
    control_stems = ctx.control_stems
    control_lookup = ctx.control_lookup

    run_mode = subrun.run_mode
    csv_prefix = subrun.csv_prefix
    stage_name_prefix = subrun.stage_name_prefix

    cfg.setdefault("_EFFECTIVE_SINGLE_LIGAND", "")
    single_ligand_hit: Optional[Path] = None
    cfg.pop("_SINGLE_RESOLVED_PATH", None)
    if cfg["_EFFECTIVE_SINGLE_LIGAND"]:
        _ensure_single_ligand_index(cfg, paths, logger)
        cfg.setdefault("paths", {})
        cfg["paths"]["prepped_ligands_dir"] = str(paths.prepped_ligands_dir)

        hit = _resolve_single_ligand(cfg["_EFFECTIVE_SINGLE_LIGAND"], pdb_id, cfg, logger)
        if hit:
            cfg["_SINGLE_RESOLVED_PATH"] = str(hit)
            single_ligand_hit = hit
        else:
            selector_token = cfg["_EFFECTIVE_SINGLE_LIGAND"]
            suggestions: list[str] = []
            try:
                import difflib

                fda_map = _load_fda_name_map(cfg, logger)
                suggestions = difflib.get_close_matches(
                    selector_token,
                    list(fda_map.keys()),
                    n=5,
                    cutoff=0.7,
                )
            except Exception:
                suggestions = []
            if suggestions:
                logger.error("[single.miss.suggest] did_you_mean=%s", ", ".join(suggestions))

            allow_flag = os.environ.get("ALLOW_FDA_FALLBACK")
            if allow_flag is None:
                allow_flag = cfg.get("ALLOW_FDA_FALLBACK", False)
            if not _to_bool(allow_flag):
                logger.error(
                    "[single.block] selector '%s' not found in fda_library via FDA_MAPPING_CSV; aborting instead of fallback.",
                    selector_token,
                )
                raise SystemExit(2)
            logger.warning(
                "[single.block] selector '%s' not found; ALLOW_FDA_FALLBACK enabled, continuing with fallback flow.",
                selector_token,
            )
            cfg["_EFFECTIVE_SINGLE_LIGAND"] = ""

    if cfg.get("_EFFECTIVE_SINGLE_LIGAND"):
        if not single_ligand_hit:
            return
        ligands = [str(single_ligand_hit)]
        ha = _count_heavy_atoms_from_pdbqt(single_ligand_hit)
        heavy_atom_counts = {str(single_ligand_hit): ha}
        pains_flags = {}
        logger.info(f"[single] Active ? docking only: {single_ligand_hit.name} (heavy={ha})")

        if cfg.get("PH_LIGAND_MODE", "").lower() == "context_window" and cfg.get("PH_ENSEMBLE_IN_PREP"):
            try:
                from prep_ligands import enumerate_ligands_for_docking

                ph_values = [6.0, 8.0]
                if "_PH_CONTEXT_VALUES" in cfg:
                    ph_values = cfg["_PH_CONTEXT_VALUES"]
                ligand_window = sorted(
                    {round(p, 1) for ph in ph_values for p in (float(ph) - 1.0, float(ph), float(ph) + 1.0)}
                )
                logger.info(f"[single.ph_ligand] Using ligand window {ligand_window}")

                if run_mode in {"dud", "hmdb"}:
                    ph_override_single = run_mode
                elif run_mode == "fda":
                    ph_override_single = "off"
                else:
                    ph_override_single = None
                _, noncontrol_roots_single = _lib_roots_for_pdb(
                    cfg,
                    paths.pdb_id.upper(),
                    paths,
                    logger,
                    test_mode_override=ph_override_single,
                )
                ph_root_path = noncontrol_roots_single[0] if noncontrol_roots_single else None
                ph_root_cfg = str(ph_root_path) if ph_root_path else ""

                logger.info(
                    "[single.ph_ligand.bridge] ph_root_cfg=%s ph_root_path=%s exists=%s",
                    ph_root_cfg,
                    str(ph_root_path) if ph_root_path is not None else "",
                    ph_root_path.exists() if ph_root_path is not None else False,
                )
                if ph_root_path is not None and ph_root_path.exists():
                    enumerate_ligands_for_docking(
                        requested_ph_values=ligand_window,
                        root_dir=ph_root_path,
                        microstate_dedup=True,
                        force=False,
                    )
                else:
                    logger.info(
                        "[single.ph_ligand.bridge.skip] no valid ph_ligand_root; "
                        "skipping microstate priming for single-ligand mode"
                    )
            except Exception as e:
                logger.warning(f"[single.ph_ligand.skip] Could not run PH-ligand window for single mode: {e}")
    else:
        ligands, heavy_atom_counts, pains_flags = prepare_and_filter_ligands(
            cfg,
            paths,
            logger,
            run_mode=run_mode,
        )

    def _norm_dedupe(seq):
        seen = set()
        out = []
        for p in seq:
            pn = norm(p)
            if pn not in seen:
                seen.add(pn)
                out.append(pn)
        return out

    control_stems_lower_init = {s.lower() for s in control_stems}
    control_pool_raw, noncontrol_pool_raw = _split_controls_and_noncontrols(
        ligands,
        paths.prepped_ligands_dir,
        allowed_subdirs={"controls", "reference"},
        control_stems=control_stems_lower_init,
    )

    if not cfg.get("_EFFECTIVE_SINGLE_LIGAND"):
        ctrl_stems_lower = {s.lower() for s in control_stems}

        prepped_control_pdbqts = []
        scan_roots = [paths.prepped_ligands_dir]
        if cfg.get("OUTPUT_LIGANDS_DIR"):
            try:
                out_root = Path(cfg["OUTPUT_LIGANDS_DIR"])
                if out_root.exists():
                    scan_roots.append(out_root)
            except Exception:
                pass

        for root in scan_roots:
            if root and root.exists():
                for p in root.glob("*.pdbqt"):
                    stem0 = p.stem.split("_stage")[0].lower()
                    if stem0 in ctrl_stems_lower:
                        prepped_control_pdbqts.append(p)

        lig_set = {norm(x) for x in control_pool_raw + noncontrol_pool_raw}
        missing_controls = [p for p in prepped_control_pdbqts if norm(p) not in lig_set]

        if missing_controls:
            logger.info(f"[Controls] Adding {len(missing_controls)} prepared control(s) to control pool for Stage3.")
            control_pool_raw.extend(str(p) for p in missing_controls)
            for p in missing_controls:
                try:
                    heavy_atom_counts.setdefault(str(p), _count_heavy_atoms_from_pdbqt(p))
                except Exception:
                    heavy_atom_counts.setdefault(str(p), 0)
                pains_flags.setdefault(str(p), False)

    control_pool = _norm_dedupe(control_pool_raw)
    control_norms = {norm(p) for p in control_pool}
    noncontrol_pool = [p for p in _norm_dedupe(noncontrol_pool_raw) if norm(p) not in control_norms]

    ligands = noncontrol_pool[:]

    base_controls = control_pool[:]
    base_ligands = ligands[:]
    base_heavy_atoms = dict(heavy_atom_counts)
    base_pains_flags = dict(pains_flags)
    base_center = tuple(center)
    base_box = tuple(box_size)

    ph_log = logging.getLogger("ph_ensemble")
    ph_enabled = bool(cfg.get("PH_ENSEMBLE"))
    plan_only = os.environ.get("A2_PLAN_ONLY") == "1"

    manifest_run_id = cfg.get("RUN_ID")
    library_for_manifest = None
    try:
        mode_for_manifest = _resolve_test_mode(cfg)
        if mode_for_manifest != "off":
            lib_map = cfg.get("_TEST_LIBRARY_CANONICAL", {}) or {}
            library_for_manifest = lib_map.get(paths.pdb_id.upper())
        if not library_for_manifest:
            library_for_manifest = cfg.get("LIBRARY_SUBDIR_DEFAULT")
    except Exception:
        library_for_manifest = cfg.get("LIBRARY_SUBDIR_DEFAULT")

    ph_tags = init_ph_tags_and_manifest(cfg, paths.pdb_id, variant_token, legacy_mode)
    ph_tags = [str(tag).strip() for tag in ph_tags if tag is not None and str(tag).strip()]
    if not ph_tags:
        fallback_ph = load_ph_tags(paths.pdb_id, variant=variant_token) or []
        ph_tags = [str(tag).strip() for tag in fallback_ph if tag is not None and str(tag).strip()]
    if not ph_tags:
        ph_tags = ["base"]  # ensure downstream logging/manifest updates occur
    if ph_enabled and not ph_tags:
        logger.info(
            "[subrun.ph] run_mode=%s ph_enabled=True but no ph_tags; skipping PH run",
            run_mode or "None",
        )
        return

    if run_mode in {"dud", "hmdb"}:
        ph_test_mode_override = run_mode
    elif run_mode == "fda":
        ph_test_mode_override = "off"
    else:
        ph_test_mode_override = None
    _ctrl_roots_ph, noncontrol_roots_ph = _lib_roots_for_pdb(
        cfg,
        paths.pdb_id.upper(),
        paths,
        logger,
        test_mode_override=ph_test_mode_override,
    )
    ph_ligand_root = noncontrol_roots_ph[0] if noncontrol_roots_ph else None

    logger.info(
        "[subrun.ph] run_mode=%s ph_tags=%s ph_ligand_root=%s override=%s",
        run_mode or "None",
        ",".join(ph_tags) if ph_tags else "(none)",
        str(ph_ligand_root) if ph_ligand_root else "(none)",
        ph_test_mode_override or "(none)",
    )
    prewarm_ph_ligand_microstates(cfg, ph_tags, ph_ligand_root)

    stages_for_run = [
        {**stage, "name": f"{stage_name_prefix}{stage['name']}"}
        for stage in stages
    ]

    for ph_label in ph_tags:
        ph_start_ts = time.time()
        try:
            update_manifest_for_protein_start(
                cfg,
                manifest_run_id or "",
                paths.pdb_id,
                variant_label,
                library_for_manifest,
                ph_tag=ph_label,
            )
        except Exception:
            ph_log.warning(
                "[run-manifest.skip] pdb=%s variant=%s ph=%s reason=start",
                paths.pdb_id,
                variant_label,
                ph_label if ph_label else "base",
                exc_info=True,
            )

        rec_path = receptor_file(paths.pdb_id, variant=variant_token, ph_tag=ph_label, legacy=legacy_mode)
        out_root = docked_dir(paths.pdb_id, variant=variant_token, ph_tag=ph_label, legacy=legacy_mode)
        ph_print = ph_label or "(none)"
        rec_exists = rec_path.exists()
        logger.info(
            "[router] pdb=%s variant=%s ph=%s receptor_file=%s docked_dir=%s exists=%s",
            paths.pdb_id,
            variant_label,
            ph_print,
            str(rec_path),
            str(out_root),
            rec_exists,
        )

        if plan_only:
            print(
                f"pdb={paths.pdb_id} variant={variant_label} ph={ph_print} "
                f"receptor_file={rec_path} docked_dir={out_root} exists={rec_exists}"
            )
            continue

        _record_apo_holo_usage(cfg, paths.pdb_id, variant_token, ph_label, rec_path)

        if not rec_exists:
            ph_log.warning(
                "[ph_ensemble.dock.stage] pdb_id=%s variant=%s ph=%s receptor_missing=%s",
                paths.pdb_id,
                variant_label,
                ph_print,
                str(rec_path),
            )
            try:
                update_manifest_for_protein_failure(
                    cfg,
                    manifest_run_id or "",
                    paths.pdb_id,
                    variant_label,
                    rec_path,
                    ph_tag=ph_label,
                )
            except Exception:
                ph_log.warning(
                    "[run-manifest.skip] pdb=%s variant=%s ph=%s reason=receptor-missing",
                    paths.pdb_id,
                    variant_label,
                    ph_label if ph_label else "base",
                    exc_info=True,
                )
            continue

        try:

            ligands = base_ligands[:]
            controls_for_run = base_controls[:]
            heavy_atom_counts = dict(base_heavy_atoms)
            pains_flags = dict(base_pains_flags)
            center = tuple(base_center)
            box_size = tuple(base_box)
            receptor_pdbqt = str(rec_path)

            enumerated = enumerate_ligands_for_ph_context(
                cfg=cfg,
                pdb_id=paths.pdb_id,
                ph_label=ph_label,
                ph_ligand_root=ph_ligand_root,
            )

            if enumerated:
                ligands = [str(p) for p in enumerated]
                heavy_atom_counts = {
                    str(p): _count_heavy_atoms_from_pdbqt(p) for p in enumerated
                }
                pains_flags = {
                    k: base_pains_flags.get(
                        k,
                        base_pains_flags.get(Path(k).stem, False),
                    )
                    for k in ligands
                }
                for c in controls_for_run:
                    if c not in heavy_atom_counts and c in base_heavy_atoms:
                        heavy_atom_counts[c] = base_heavy_atoms[c]
                    pains_flags.setdefault(
                        c,
                        base_pains_flags.get(c, base_pains_flags.get(Path(c).stem, False)),
                    )
            else:
                ligands = base_ligands[:]
                heavy_atom_counts = dict(base_heavy_atoms)
                pains_flags = dict(base_pains_flags)

            ctrl_stems_lower = control_stems_lower_init
            ctrl_blacklist = {t.strip().upper() for t in str(cfg.get("CONTROL_BLACKLIST", "")).split(",") if t.strip()}
            min_ha = int(cfg.get("CONTROL_MIN_HEAVY_ATOMS", 10))

            def _eligible_control_path(p: str) -> bool:
                stem = Path(p).stem.split("_stage")[0].lower()
                if stem.upper() in ctrl_blacklist:
                    return False
                if stem not in ctrl_stems_lower:
                    return False
                ha = heavy_atom_counts.get(p)
                return (ha is None) or (ha >= min_ha)

            controls_for_run = [c for c in controls_for_run if _eligible_control_path(c)]
            control_norms_for_run = {norm(c) for c in controls_for_run}

            if not ligands and controls_for_run:
                logger.info(
                    "[Controls] No non-control ligands selected; controls will run in Stage3 only."
                )
            if not ligands and not controls_for_run:
                logger.warning("No valid ligands after filtering; skipping protein.")
                continue

            selector = CenterSelector(cfg, logger, control_stems, heavy_atom_counts, center)
            guard = GlobalCenterGuard(
                max_global_switches=int(cfg.get("MAX_GLOBAL_CENTER_SWITCHES", 2))
            )

            stage1_original = [l for l in ligands if norm(l) not in control_norms_for_run]

            # --- NO_LIBRARY_DOCKING: only plan ligands, do not dock -----------
            if bool(cfg.get("NO_LIBRARY_DOCKING", False)):
                # How many ligands would be docked for this PH / mode?
                preview_n = int(cfg.get("NO_DOCKING_PREVIEW_N", 10))
                preview_names = ", ".join(
                    Path(l).name for l in stage1_original[:preview_n]
                )

                logger.info(
                    "[no-docking-planned] pdb=%s variant=%s ph=%s mode=%s "
                    "n_stage1=%d preview=[%s]",
                    paths.pdb_id,
                    variant_label,
                    ph_label if ph_label else "base",
                    run_mode or "(unspecified)",
                    len(stage1_original),
                    preview_names,
                )

                # Also drop a simple text file with the full Stage1 ligand list
                try:
                    out_dir = docked_dir(
                        paths.pdb_id,
                        variant=variant_token,
                        ph_tag=ph_label,
                        legacy=legacy_mode,
                    )
                    out_dir.mkdir(parents=True, exist_ok=True)
                    suffix = run_mode or "run"
                    out_txt = out_dir / f"planned_ligands_{suffix}.txt"

                    with out_txt.open("w") as fh:
                        for lig in stage1_original:
                            fh.write(f"{lig}\n")

                    logger.info(
                        "[no-docking-planned] Wrote planned ligands to %s",
                        out_txt,
                    )
                except Exception as e:
                    logger.warning(
                        "[no-docking-planned] Failed to write planned ligand list: %s",
                        e,
                    )

                # Skip all docking stages for this PH context.
                continue

            forced_extracted_for_stage3 = set(controls_for_run)
            logger.info(
                "[Force-carry] Stage3 control pool size=%d noncontrols_stage1=%d",
                len(forced_extracted_for_stage3),
                len(stage1_original),
            )

            score_history: Dict[str, Dict[str, Dict]] = defaultdict(dict)
            score_history_gnina: Dict[str, Dict[str, Dict]] = defaultdict(dict)
            gnina_metrics_by_stage_best: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(dict)
            gnina_metrics_by_stage_full: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(dict)
            ledock_metrics_by_stage: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(dict)
            dock6_metrics_by_stage: Dict[str, Dict[Path, Dict[str, Any]]] = defaultdict(dict)
            gnina_records: List[Dict[str, Any]] = []
            validated_ligands_last: List[str] = []
            recenter_attempts = 0
            docking_mode = cfg.get("DOCKING_MODE", "discovery").lower()

            retry_mgr = RetryManager()
            gnina_jobs: List[Dict[str, Any]] = []
            ledock_jobs: List[Dict[str, Any]] = []
            dock6_jobs: List[Dict[str, Any]] = []
            ledock_enabled = _use_ledock(cfg)
            dock6_enabled = _use_dock6(cfg)

            i = 0
            while i < len(stages_for_run):
                guard.reset_stage()
                stage = stages_for_run[i]
                stage_is_stage3 = _is_stage3(stage["name"])
                stage_controls = sorted(controls_for_run, key=lambda p: Path(p).name) if stage_is_stage3 else []
                stage_noncontrols = [
                    l
                    for l in ligands
                    if norm(l) not in control_norms_for_run
                    and Path(l).stem.split("_stage")[0].lower() not in ctrl_stems_lower
                ]
                stage_ligands = (
                    _interleave_controls(stage_noncontrols, stage_controls)
                    if stage_is_stage3
                    else stage_noncontrols
                )

                if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                    fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
                    if checkpoint_should_skip(
                        cfg,
                        paths.pdb_id,
                        stage["name"],
                        fp,
                        ph_label=ph_label,
                        variant=variant_env or None,
                        engine="vina",
                    ):
                        logger.info(f"[Checkpoint] Skipping {stage['name']} (fingerprint matched).")
                        i += 1
                        continue

                if not stage_ligands:
                    logger.warning(f"No ligands to dock at {stage['name']}; skipping this stage.")
                    i += 1
                    continue

                stage_ligands_for_audit = [norm(l) for l in stage_ligands]
                stage_dir = paths.docked_stage_dir(variant_env or None, stage["name"], ph_label)

                logger.info(
                    "[ligand.stage] stage=%s ph=%s controls=%d noncontrols=%d docking_ligands=%d",
                    stage["name"],
                    ph_label if ph_label else "base",
                    len(stage_controls),
                    len(stage_noncontrols),
                    len(stage_ligands),
                )
                if ph_label:
                    ph_log.info(
                        "[ph_ensemble.dock.stage] pdb_id=%s variant=%s ph=%s stage=%s receptor=%s out=%s",
                        paths.pdb_id,
                        variant_label,
                        ph_label,
                        stage["name"],
                        receptor_pdbqt,
                        str(stage_dir),
                    )

                scores, validated, distances, raw_docked, invalids = run_one_stage(
                    cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage,
                    stage_ligands, logger, retry_mgr, control_lookup, ph_label=ph_label
                )

                validated_ligands_last = validated

                ctrl_blacklist_set = ctrl_blacklist

                def _is_control(lig: str) -> bool:
                    lig_norm = norm(lig)
                    stem = Path(lig).stem.split("_stage")[0].lower()
                    if stem.upper() in ctrl_blacklist_set:
                        return False
                    if lig_norm not in control_norms_for_run and stem not in ctrl_stems_lower:
                        return False
                    ha = heavy_atom_counts.get(lig)
                    if ha is not None and ha < min_ha:
                        return False
                    return True

                control_anchor_hit = any(_is_control(lig) for lig in validated)

                lock_score_max = float(cfg.get("CONTROL_LOCK_SCORE_MAX", float("inf")))
                lock_min_hits = int(cfg.get("CONTROL_LOCK_MIN_HITS", 1))
                lock_center_max = float(cfg.get("CONTROL_LOCK_CENTER_MAX_DIST", 4.0))

                qualified_controls = []
                for lig in validated:
                    if not _is_control(lig):
                        continue
                    sc = scores.get(lig)
                    if sc is None or not np.isfinite(sc):
                        continue
                    if sc > lock_score_max:
                        continue
                    pose_path = raw_docked.get(lig)
                    if not pose_path:
                        continue
                    c = CenterSelector._pdbqt_centroid(pose_path)
                    if c is None:
                        continue
                    if np.linalg.norm(c - np.array(center, float)) <= lock_center_max:
                        qualified_controls.append(lig)

                if len(qualified_controls) >= lock_min_hits and not guard.locked:
                    guard.lock()
                    logger.info(
                        "[CONTROL-LOCK] Control(s) validated with strong confidence "
                        f"(n={len(qualified_controls)}, score<={lock_score_max}, dist<={lock_center_max} Ang); "
                        "center is now anchored; future center switches are disabled."
                    )

                try:
                    processed = {norm(x) for x in stage_ligands}
                    valid_set = {norm(x) for x in scores.keys()}
                    invalid_set = {norm(x) for x in invalids.keys()}
                    both = valid_set & invalid_set
                    missing = processed - (valid_set | invalid_set)
                    if both or missing:
                        logger.error(
                            f"Invariant violation at {stage['name']}: both={len(both)}, missing={len(missing)}"
                        )
                        if both:
                            logger.error(
                                "Ligands marked both valid & invalid: "
                                + ", ".join(os.path.basename(x) for x in list(both)[:10])
                            )
                        if missing:
                            logger.error(
                                "Ligands missing from results: "
                                + ", ".join(os.path.basename(x) for x in list(missing)[:10])
                            )
                            for lig_m in missing:
                                invalids[lig_m] = (None, "not_processed")
                except Exception as _e:
                    logger.warning(f"Invariant check failed: {_e}")

                for lig, sc in scores.items():
                    record_score(score_history, stage['name'], lig, sc, True)
                    record_le(score_history, stage['name'], lig, sc, heavy_atom_counts)
                for lig, (sc, reason) in invalids.items():
                    record_score(score_history, stage['name'], lig, sc, False, reason=reason)
                    record_le(score_history, stage['name'], lig, sc, heavy_atom_counts)

                promoted_this_stage = False
                try:
                    decision = selector.consider_switch(
                        stage['name'], scores, validated, raw_docked, receptor_pdbqt, center, guard
                    )
                    if decision.promoted and decision.new_center is not None:
                        old = center
                        center = decision.new_center
                        promoted_this_stage = True
                        guard.mark_switch()
                        if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                            checkpoint_invalidate_from(
                                cfg,
                                paths.pdb_id,
                                stages_for_run,
                                start_index=i,
                                ph_label=ph_label,
                                variant=variant_env or None,
                            )
                        logger.info(
                            f"[CENTER] Switched from {old} -> {center} ({decision.reason}, "
                            f"SwitchScore={decision.switchscore:.2f}) [global switch]"
                        )
                except Exception as e:
                    logger.warning(f"CenterSelector failed gracefully: {e}")

                if ph_label:
                    ph_log.info(
                        "[ph_ensemble.dock.scores] pdb_id=%s variant=%s ph=%s stage=%s valid=%d invalid=%d",
                        paths.pdb_id,
                        variant_label,
                        ph_label,
                        stage["name"],
                        len(scores),
                        len(invalids),
                    )

                if not promoted_this_stage:
                    restart, center, box_size, redo_ligands, recenter_attempts = early_recenter_decision(
                        i, scores, distances, box_size, center, stage1_original, recenter_attempts, params,
                        cfg, paths.pdb_id, receptor_pdbqt, logger, raw_docked, guard, control_anchor_hit
                    )
                    if restart:
                        ligands = redo_ligands
                        if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                            checkpoint_invalidate_from(
                                cfg,
                                paths.pdb_id,
                                stages_for_run,
                                start_index=0,
                                ph_label=ph_label,
                                variant=variant_env or None,
                            )
                        gnina_jobs.clear()
                        i = 0
                        continue

                try:
                    if cfg.get("ADAPTIVE_SHRINK_ENABLE", True) and validated:
                        med = (
                            float(np.median([d for d in distances if isinstance(d, (int, float))]))
                            if distances
                            else None
                        )
                        if (med is not None) and (med < float(cfg.get("ADAPTIVE_SHRINK_MEDIAN_MAX", 4.0))):
                            dec = float(cfg.get("ADAPTIVE_SHRINK_DEC", 4.0))
                            min_box = float(cfg.get("ADAPTIVE_SHRINK_MIN_BOX", 14.0))
                            new_box = tuple(max(min_box, s - dec) for s in box_size)
                            if new_box != box_size:
                                logger.info(
                                    f"Adaptive shrink: median dist {med:.2f} A -> box {box_size} -> {new_box}"
                                )
                                box_size = new_box
                except Exception as _e:
                    logger.warning(f"Adaptive shrink skipped: {_e}")

                if i < len(stages_for_run) - 1:
                    if not scores:
                        restart, center, box_size, redo_ligands = fallback_recentering_if_empty(
                            cfg, paths.pdb_id, stage['name'], scores, raw_docked,
                            receptor_pdbqt, center, box_size, stage1_original, logger, guard, control_anchor_hit
                        )
                        if restart:
                            ligands = redo_ligands
                            if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                                checkpoint_invalidate_from(
                                    cfg,
                                    paths.pdb_id,
                                    stages_for_run,
                                    start_index=0,
                                    ph_label=ph_label,
                                    variant=variant_env or None,
                                )
                            gnina_jobs.clear()
                            i = 0
                            continue
                        else:
                            break

                    next_stage_is_stage3 = _is_stage3(stages_for_run[i + 1]["name"])
                    scores_for_selection = {
                        lig: sc for lig, sc in scores.items() if norm(lig) not in control_norms_for_run
                    }
                    invalids_for_selection = {
                        lig: val for lig, val in invalids.items() if norm(lig) not in control_norms_for_run
                    }
                    ligands = _apply_force_carry_and_doping(
                        cfg,
                        docking_mode,
                        i,
                        stages_for_run,
                        scores_for_selection,
                        logger,
                        stage1_original=stage1_original,
                        forced_extracted_for_stage3=forced_extracted_for_stage3,
                        invalids=invalids_for_selection,
                        next_stage_is_stage3=next_stage_is_stage3,
                        control_norms=control_norms_for_run,
                    )
                    if not ligands:
                        logger.warning(
                            f"No ligands selected for {stages_for_run[i + 1]['name']}; stopping."
                        )
                        break

                # Engine follow-ups controlled by fpocket druggability tiers (or legacy flags).
                if validated:
                    policy = decide_engine_policy(
                        cfg=cfg,
                        pdb_id=paths.pdb_id,
                        variant=variant_env or None,
                        ph_label=ph_label,
                        center=center,
                        ledock_enabled=ledock_enabled,
                        dock6_enabled=dock6_enabled,
                        logger=logger,
                    )

                    run_id_token = str(cfg.get("RUN_ID") or "")
                    if run_id_token:
                        try:
                            update_manifest_for_druggability_and_engine_plan(
                                cfg=cfg,
                                run_id=run_id_token,
                                pdb_id=paths.pdb_id,
                                variant_label=variant_env or None,
                                ph_tag=ph_label,
                                tier=policy.tier,
                                use_gnina=policy.use_gnina,
                                use_ledock=policy.use_ledock,
                                use_dock6=policy.use_dock6,
                            )
                        except Exception:
                            logger.warning(
                                "[run-manifest.druggability-plan.skip] run_id=%s pdb=%s variant=%s ph=%s",
                                run_id_token,
                                paths.pdb_id,
                                variant_env,
                                ph_label,
                                exc_info=True,
                            )

                    if policy.use_gnina:
                        gnina_jobs.append(
                            {
                                "stage_name": stage["name"],
                                "stage_info": dict(stage),
                                "ligands": list(validated),
                                "center": tuple(center) if center is not None else None,
                                "box_size": tuple(box_size) if box_size is not None else None,
                            }
                        )
                    else:
                        logger.info(
                            "[gnina.skip] pdb=%s ph=%s tier=%s reason=%s",
                            paths.pdb_id,
                            ph_label if ph_label else "base",
                            policy.tier,
                            policy.reason,
                        )

                    if policy.use_ledock:
                        ledock_jobs.append(
                            {
                                "stage_name": stage["name"],
                                "stage_info": dict(stage),
                                "ligands": list(validated),
                                "center": tuple(center) if center is not None else None,
                                "box_size": tuple(box_size) if box_size is not None else None,
                            }
                        )
                    if policy.use_dock6:
                        dock6_jobs.append(
                            {
                                "stage_name": stage["name"],
                                "stage_info": dict(stage),
                                "ligands": list(validated),
                                "center": tuple(center) if center is not None else None,
                                "box_size": tuple(box_size) if box_size is not None else None,
                                "variant": variant_env or None,
                                "ph_label": ph_label,
                            }
                        )

                # Stage-level completion audit ensures every ligand has an output
                # artifact (pose or failure marker) before marking the checkpoint.
                threads_per_vina = int(cfg.get("THREADS_PER_VINA", 1))

                def _expected_vina_path(lig: str) -> Path:
                    return stage_dir / f"{Path(lig).stem}_{stage['name']}.pdbqt"

                def _rerun_vina_missing(lig: str) -> tuple[bool, Optional[str], Optional[Path]]:
                    stage_for_cfg = dict(stage)
                    stage_for_cfg["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))
                    if cfg.get("FAST_MODE"):
                        stage_for_cfg["exhaustiveness"] = 1
                    conf_path, out_path = emit_vina_config(
                        cfg,
                        paths.pdb_id,
                        receptor_pdbqt,
                        center,
                        box_size,
                        lig,
                        stage["name"],
                        stage_for_cfg,
                        threads_per_vina,
                        logger,
                        variant=variant_env or None,
                        ph_token=ph_label,
                        legacy=legacy_mode,
                    )
                    try:
                        Path(conf_path).resolve().relative_to(Path(cfg["CONFIG_RUN_DIR"]).resolve())
                    except Exception:
                        return False, "rerun_config_outside_run_dir", None
                    try:
                        run_docking_task(
                            cfg["VINA_EXE"],
                            conf_path,
                            lig,
                            out_path,
                            write_failure_marker_flag=True,
                        )
                    except Exception as exc:
                        return False, f"rerun_error:{exc}", None
                    if Path(out_path).exists() and Path(out_path).stat().st_size > 0:
                        return True, "rerun_ok", None
                    return False, "rerun_no_output", None

                completion_report_vina = run_completion_audit(
                    engine="vina",
                    pdb_id=paths.pdb_id,
                    stage_name=stage["name"],
                    ligands=stage_ligands_for_audit,
                    expected_output_path=_expected_vina_path,
                    rerun_one=_rerun_vina_missing,
                    stage_dir=stage_dir,
                    cfg=cfg,
                    logger=logger,
                    retries=1,
                    ph_label=ph_label,
                    variant=variant_env or None,
                )
                missing_after_vina = completion_report_vina.get("missing_ligands_after") or []
                if missing_after_vina:
                    for lig_miss in missing_after_vina:
                        if lig_miss not in score_history.get(stage["name"], {}):
                            record_score(
                                score_history,
                                stage["name"],
                                lig_miss,
                                None,
                                False,
                                reason="completion_missing",
                            )
                            record_le(score_history, stage["name"], lig_miss, None, heavy_atom_counts)

                if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                    try:
                        fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
                        if completion_report_vina.get("success", False):
                            checkpoint_mark_done(
                                cfg,
                                paths.pdb_id,
                                stage["name"],
                                fp,
                                ph_label=ph_label,
                                variant=variant_env or None,
                            )
                        else:
                            logger.warning(
                                "[Checkpoint] defer mark stage=%s pdb=%s reason=completion_missing",
                                stage["name"],
                                paths.pdb_id,
                            )
                    except Exception:
                        pass

                i += 1

            if gnina_jobs:
                logger.info(
                    "[gnina.scheduler] pdb=%s variant=%s ph=%s policy=after_all_vina_stages n_jobs=%d",
                    paths.pdb_id,
                    variant_label,
                    ph_label if ph_label else "base",
                    len(gnina_jobs),
                )
                for job in gnina_jobs:
                    gnina_stage_name = f"gnina_{job['stage_name']}"
                    gnina_fp = None
                    if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                        try:
                            gnina_fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, job["stage_info"])
                            if checkpoint_should_skip(
                                cfg,
                                paths.pdb_id,
                                gnina_stage_name,
                                gnina_fp,
                                ph_label=ph_label,
                                variant=variant_env or None,
                                engine="gnina",
                            ):
                                logger.info(f"[Checkpoint] Skipping {gnina_stage_name} (fingerprint matched).")
                                if manifest_run_id:
                                    try:
                                        update_manifest_for_docking_stage(
                                            cfg,
                                            manifest_run_id,
                                            paths.pdb_id,
                                            variant_label,
                                            gnina_stage_name,
                                            status="completed",
                                            ph_tag=ph_label,
                                            elapsed_sec=0.0,
                                        )
                                    except Exception:
                                        logger.warning(
                                            "[run-manifest.docking-stage] failed to record GNINA checkpoint skip pdb=%s variant=%s ph=%s stage=%s",
                                            paths.pdb_id,
                                            variant_label,
                                            ph_label if ph_label is not None else "base",
                                            gnina_stage_name,
                                            exc_info=True,
                                        )
                                continue
                        except Exception:
                            gnina_fp = None
                    gnina_start_ts = time.time()
                    if manifest_run_id:
                        try:
                            update_manifest_for_docking_stage(
                                cfg,
                                manifest_run_id,
                                paths.pdb_id,
                                variant_label,
                                gnina_stage_name,
                                status="running",
                                ph_tag=ph_label,
                            )
                        except Exception:
                            logger.warning(
                                "[run-manifest.docking-stage] failed to record GNINA start pdb=%s variant=%s ph=%s stage=%s",
                                paths.pdb_id,
                                variant_label,
                                ph_label if ph_label is not None else "base",
                                gnina_stage_name,
                                exc_info=True,
                            )

                    try:
                        scores_gnina, gnina_metrics, gnina_completion = run_gnina_for_stage(
                            cfg=cfg,
                            paths=paths,
                            pdb_id=paths.pdb_id,
                            variant=variant_env or None,
                            ph_label=ph_label,
                            stage_name=job["stage_name"],
                            stage_info=job["stage_info"],
                            ligands=job["ligands"],
                            center=job["center"],
                            box_size=job["box_size"],
                            logger=logger,
                            receptor_pdbqt=receptor_pdbqt,
                            control_lookup=control_lookup,
                                )
                        gnina_elapsed = time.time() - gnina_start_ts
                        if manifest_run_id:
                            try:
                                update_manifest_for_docking_stage(
                                    cfg,
                                    manifest_run_id,
                                    paths.pdb_id,
                                    variant_label,
                                    gnina_stage_name,
                                    status="completed",
                                    ph_tag=ph_label,
                                    elapsed_sec=gnina_elapsed,
                                )
                            except Exception:
                                logger.warning(
                                    "[run-manifest.docking-stage] failed to record GNINA end pdb=%s variant=%s ph=%s stage=%s",
                                    paths.pdb_id,
                                    variant_label,
                                    ph_label if ph_label is not None else "base",
                                    gnina_stage_name,
                                    exc_info=True,
                                )
                        if bool(cfg.get("CHECKPOINT_ENABLE", True)) and gnina_fp is not None:
                            try:
                                if gnina_completion.get("success", False):
                                    checkpoint_mark_done(
                                        cfg,
                                        paths.pdb_id,
                                        gnina_stage_name,
                                        gnina_fp,
                                        ph_label=ph_label,
                                        variant=variant_env or None,
                                    )
                                else:
                                    logger.warning(
                                        "[Checkpoint] defer mark gnina_stage=%s pdb=%s reason=completion_missing",
                                        gnina_stage_name,
                                        paths.pdb_id,
                                    )
                            except Exception:
                                pass
                    except Exception as e:
                        if manifest_run_id:
                            try:
                                update_manifest_for_docking_stage(
                                    cfg,
                                    manifest_run_id,
                                    paths.pdb_id,
                                    variant_label,
                                    gnina_stage_name,
                                    status="failed",
                                    ph_tag=ph_label,
                                    error=str(e),
                                    elapsed_sec=time.time() - gnina_start_ts,
                                )
                            except Exception:
                                logger.warning(
                                    "[run-manifest.docking-stage] failed to record GNINA failure pdb=%s variant=%s ph=%s stage=%s",
                                    paths.pdb_id,
                                    variant_label,
                                    ph_label if ph_label is not None else "base",
                                    gnina_stage_name,
                                    exc_info=True,
                                )
                        raise
                    gnina_stage_key = f"gnina_{job['stage_name']}"
                    for lig, sc in scores_gnina.items():
                        metrics_entry = gnina_metrics.get(
                            lig,
                            {
                                "minimized_affinity_kcal": None,
                                "cnn_score": None,
                                "cnn_affinity_pK": None,
                                "gnina_primary_score": sc,
                                "valid": False,
                                "reason": "gnina_failed" if sc is None else "",
                                "self_rmsd": None,
                            },
                        )
                        primary_score = metrics_entry.get("gnina_primary_score", sc)
                        minimized_affinity = metrics_entry.get("minimized_affinity_kcal")
                        cnn_score = metrics_entry.get("cnn_score")
                        cnn_affinity = metrics_entry.get("cnn_affinity_pK")
                        valid_flag = bool(metrics_entry.get("valid", False))
                        reason_str = metrics_entry.get("reason", "") or ("gnina_failed" if not valid_flag else "")
                        self_rmsd_val = metrics_entry.get("self_rmsd", None)

                        ha_val = heavy_atom_counts.get(lig)
                        le_val = compute_ligand_efficiency(primary_score, ha_val)
                        pains_val = pains_flags.get(lig, pains_flags.get(Path(lig).stem, False))

                        gnina_metrics_by_stage_full[gnina_stage_key][lig] = {
                            "minimized_affinity_kcal": minimized_affinity,
                            "cnn_score": cnn_score,
                            "cnn_affinity_pK": cnn_affinity,
                            "gnina_primary_score": primary_score,
                            "valid": bool(valid_flag),
                            "reason": reason_str,
                            "heavy_atoms": int(ha_val) if isinstance(ha_val, (int, float)) else None,
                            "le": le_val,
                            "self_rmsd": self_rmsd_val,
                            "pains_flag": bool(pains_val) if pains_val is not None else False,
                        }

                        gnina_records.append(
                            {
                                "stage_name": gnina_stage_key,
                                "ligand": lig,
                                "primary_score": primary_score,
                                "minimized_affinity_kcal": minimized_affinity,
                                "cnn_score": cnn_score,
                                "cnn_affinity_pK": cnn_affinity,
                                "valid": bool(valid_flag),
                                "reason": reason_str or "",
                                "heavy_atoms": int(ha_val) if isinstance(ha_val, (int, float)) else None,
                                "le": le_val,
                                "self_rmsd": self_rmsd_val,
                                "pains_flag": bool(pains_val) if pains_val is not None else False,
                            }
                        )
            else:
                logger.info(
                    "[gnina.scheduler] pdb=%s variant=%s ph=%s action=skip reason=no_jobs",
                    paths.pdb_id,
                    variant_label,
                    ph_label if ph_label else "base",
                )

            if gnina_records:
                gnina_primary_scores: Dict[str, float] = {}
                for rec in gnina_records:
                    lig = rec["ligand"]
                    sc = rec.get("primary_score")
                    if not isinstance(sc, (int, float)) or not math.isfinite(sc):
                        continue
                    prev = gnina_primary_scores.get(lig)
                    if prev is None or sc > prev:
                        gnina_primary_scores[lig] = sc

                gnina_stage_membership: Dict[int, List[str]] = {}
                gnina_stage_index_for_ligand: Dict[str, int] = {}
                if gnina_primary_scores:
                    try:
                        gnina_stage_membership = compute_stage_membership_from_scores(
                            cfg=cfg,
                            docking_mode=docking_mode,
                            scores=gnina_primary_scores,
                            higher_is_better=True,
                            n_stages=len(stages_for_run),
                        )
                    except Exception as e:
                        logger.warning("[gnina.staging] failed to compute stage membership: %s", e)
                        gnina_stage_membership = {}

                for stage_idx, lig_list in (gnina_stage_membership or {}).items():
                    for lig in lig_list:
                        gnina_stage_index_for_ligand[lig] = stage_idx

                best_rec_for_lig: Dict[str, Dict[str, Any]] = {}
                for rec in gnina_records:
                    lig = rec["ligand"]
                    sc = rec.get("primary_score")
                    prev = best_rec_for_lig.get(lig)
                    if prev is None:
                        best_rec_for_lig[lig] = rec
                        continue
                    prev_sc = prev.get("primary_score")
                    if (
                        isinstance(sc, (int, float))
                        and math.isfinite(sc)
                        and (not isinstance(prev_sc, (int, float)) or not math.isfinite(prev_sc) or sc > prev_sc)
                    ):
                        best_rec_for_lig[lig] = rec

                for lig, rec in best_rec_for_lig.items():
                    stage_idx = gnina_stage_index_for_ligand.get(lig)
                    if stage_idx is not None:
                        gnina_stage_name = f"gnina_stage{stage_idx}"
                    else:
                        gnina_stage_name = rec.get("stage_name") or "gnina_stage1"

                    record_score(
                        score_history_gnina,
                        gnina_stage_name,
                        lig,
                        rec["primary_score"],
                        rec["valid"],
                        reason=rec["reason"] or None,
                    )
                    record_le(score_history_gnina, gnina_stage_name, lig, rec["primary_score"], heavy_atom_counts)

                    gnina_metrics_by_stage_best[gnina_stage_name][lig] = {
                        "minimized_affinity_kcal": rec["minimized_affinity_kcal"],
                        "cnn_score": rec["cnn_score"],
                        "cnn_affinity_pK": rec["cnn_affinity_pK"],
                        "gnina_primary_score": rec["primary_score"],
                        "valid": rec["valid"],
                        "reason": rec["reason"],
                        "heavy_atoms": rec["heavy_atoms"],
                        "le": rec["le"],
                        "self_rmsd": rec["self_rmsd"],
                        "pains_flag": rec["pains_flag"],
                    }

            if gnina_metrics_by_stage_full:
                write_gnina_scores_csv(
                    cfg,
                    paths.pdb_id,
                    gnina_metrics_by_stage_full,
                    ph_label=ph_label,
                    variant=variant_env or None,
                    csv_prefix=csv_prefix,
                )

            if dock6_jobs and dock6_enabled:
                try:
                    logger.info(
                        "[DOCK6_PREP_ONCE] ensuring DOCK6 site for pdb=%s variant=%s ph=%s",
                        paths.pdb_id,
                        variant_label,
                        ph_label if ph_label else "base",
                    )
                    prep_dock6.ensure_dock6_site(
                        cfg=cfg,
                        pdb_id=paths.pdb_id,
                        variant=variant_env or None,
                        ph_label=ph_label,
                        logger=logger,
                    )
                except Exception as exc:
                    logger.warning(
                        "[dock6.surface.skip] pdb=%s variant=%s ph=%s reason=%s",
                        paths.pdb_id,
                        variant_label,
                        ph_label if ph_label else "base",
                        exc,
                    )

            if dock6_jobs:
                logger.info(
                    "[dock6.scheduler] pdb=%s variant=%s ph=%s jobs=%d",
                    paths.pdb_id,
                    variant_label,
                    ph_label if ph_label else "base",
                    len(dock6_jobs),
                )
                for job in dock6_jobs:
                    stage_name = job["stage_name"]
                    stage_info = job["stage_info"]
                    ligands = [Path(l) for l in job.get("ligands", [])]
                    dock6_stage_name = f"dock6_{stage_name}"
                    dock6_start_ts = time.time()
                    if manifest_run_id:
                        try:
                            update_manifest_for_docking_stage(
                                cfg,
                                manifest_run_id,
                                paths.pdb_id,
                                variant_label,
                                dock6_stage_name,
                                status="running",
                                elapsed_sec=None,
                                ph_tag=ph_label,
                            )
                        except Exception:
                            logger.warning(
                                "[run-manifest.docking-stage] failed to record DOCK6 start pdb=%s variant=%s ph=%s stage=%s",
                                paths.pdb_id,
                                variant_label,
                                ph_label if ph_label is not None else "base",
                                dock6_stage_name,
                                exc_info=True,
                            )

                    try:
                        dock6_scores, dock6_metrics = run_dock6_for_stage(
                            cfg=cfg,
                            paths=paths,
                            pdb_id=paths.pdb_id,
                            variant=variant_env or None,
                            ph_label=ph_label,
                            stage_name=stage_name,
                            stage_info=stage_info,
                            ligands=ligands,
                            center=job.get("center"),
                            box_size=job.get("box_size"),
                            logger=logger,
                        )
                        dock6_metrics_by_stage[stage_name] = dock6_metrics
                        valid_count = sum(1 for rec in dock6_metrics.values() if rec.get("valid"))
                        invalid_count = max(len(dock6_metrics) - valid_count, 0)
                        dock6_elapsed = time.time() - dock6_start_ts
                        if manifest_run_id:
                            try:
                                update_manifest_for_docking_stage(
                                    cfg,
                                    manifest_run_id,
                                    paths.pdb_id,
                                    variant_label,
                                    dock6_stage_name,
                                    status="completed",
                                    elapsed_sec=dock6_elapsed,
                                    ph_tag=ph_label,
                                )
                            except Exception:
                                logger.warning(
                                    "[run-manifest.docking-stage] failed to record DOCK6 completion pdb=%s variant=%s ph=%s stage=%s",
                                    paths.pdb_id,
                                    variant_label,
                                    ph_label if ph_label is not None else "base",
                                    dock6_stage_name,
                                    exc_info=True,
                                )
                        logger.info(
                            "[dock6.done] pdb=%s stage=%s variant=%s ph=%s valid=%d invalid=%d elapsed_sec=%.2f",
                            paths.pdb_id,
                            stage_name,
                            variant_label,
                            ph_label if ph_label else "base",
                            valid_count,
                            invalid_count,
                            dock6_elapsed,
                        )
                    except Exception as e:
                        logger.warning(
                            "[dock6.error] pdb=%s stage=%s variant=%s ph=%s reason=%s",
                            paths.pdb_id,
                            stage_name,
                            variant_label,
                            ph_label if ph_label else "base",
                            e,
                            exc_info=True,
                        )
                        dock6_metrics_by_stage[stage_name] = {}
                        if manifest_run_id:
                            try:
                                update_manifest_for_docking_stage(
                                    cfg,
                                    manifest_run_id,
                                    paths.pdb_id,
                                    variant_label,
                                    dock6_stage_name,
                                    status="failed",
                                    ph_tag=ph_label,
                                    error=str(e),
                                    elapsed_sec=time.time() - dock6_start_ts,
                                )
                            except Exception:
                                logger.warning(
                                    "[run-manifest.docking-stage] failed to record DOCK6 failure pdb=%s variant=%s ph=%s stage=%s",
                                    paths.pdb_id,
                                    variant_label,
                                    ph_label if ph_label is not None else "base",
                                    dock6_stage_name,
                                    exc_info=True,
                                )

            elif dock6_enabled:
                logger.info(
                    "[dock6.scheduler] pdb=%s variant=%s ph=%s action=skip reason=no_jobs",
                    paths.pdb_id,
                    variant_label,
                    ph_label if ph_label else "base",
                )

            if ledock_jobs:
                logger.info(
                    "[ledock.scheduler] pdb=%s variant=%s ph=%s policy=after_gnina n_jobs=%d",
                    paths.pdb_id,
                    variant_label,
                    ph_label if ph_label else "base",
                    len(ledock_jobs),
                )
                for job in ledock_jobs:
                    stage_name = job["stage_name"]
                    ledock_stage_name = f"ledock_{stage_name}"
                    ledock_fp = None

                    if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                        try:
                            ledock_fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, job["stage_info"])
                            if checkpoint_should_skip(
                                cfg,
                                paths.pdb_id,
                                ledock_stage_name,
                                ledock_fp,
                                ph_label=ph_label,
                                variant=variant_env or None,
                                engine="ledock",
                            ):
                                logger.info(f"[Checkpoint] Skipping {ledock_stage_name} (fingerprint matched).")
                                if manifest_run_id:
                                    try:
                                        update_manifest_for_docking_stage(
                                            cfg,
                                            manifest_run_id,
                                            paths.pdb_id,
                                            variant_label,
                                            ledock_stage_name,
                                            status="completed",
                                            ph_tag=ph_label,
                                            elapsed_sec=0.0,
                                        )
                                    except Exception:
                                        logger.warning(
                                            "[run-manifest.docking-stage] failed to record LeDock checkpoint skip pdb=%s variant=%s ph=%s stage=%s",
                                            paths.pdb_id,
                                            variant_label,
                                            ph_label if ph_label is not None else "base",
                                            ledock_stage_name,
                                            exc_info=True,
                                        )
                                continue
                        except Exception:
                            ledock_fp = None

                    ledock_start_ts = time.time()
                    if manifest_run_id:
                        try:
                            update_manifest_for_docking_stage(
                                cfg,
                                manifest_run_id,
                                paths.pdb_id,
                                variant_label,
                                ledock_stage_name,
                                status="running",
                                elapsed_sec=None,
                                ph_tag=ph_label,
                            )
                        except Exception:
                            logger.warning(
                                "[run-manifest.docking-stage] failed to record LeDock start pdb=%s variant=%s ph=%s stage=%s",
                                paths.pdb_id,
                                variant_label,
                                ph_label if ph_label is not None else "base",
                                ledock_stage_name,
                                exc_info=True,
                            )
                    try:
                        scores_ledock, ledock_metrics, ledock_completion = run_ledock_for_stage(
                            cfg=cfg,
                            paths=paths,
                            pdb_id=paths.pdb_id,
                            variant=variant_env or None,
                            ph_label=ph_label,
                            stage_name=stage_name,
                            stage_info=job["stage_info"],
                            ligands=job["ligands"],
                            center=job["center"],
                            box_size=job["box_size"],
                            logger=logger,
                        )
                        ledock_metrics_by_stage[stage_name] = ledock_metrics
                        valid_count = sum(1 for rec in ledock_metrics.values() if rec.get("valid"))
                        invalid_count = max(len(ledock_metrics) - valid_count, 0)
                        ledock_elapsed = time.time() - ledock_start_ts
                        if manifest_run_id:
                            try:
                                update_manifest_for_docking_stage(
                                    cfg,
                                    manifest_run_id,
                                    paths.pdb_id,
                                    variant_label,
                                    ledock_stage_name,
                                    status="completed",
                                    elapsed_sec=ledock_elapsed,
                                    ph_tag=ph_label,
                                )
                            except Exception:
                                logger.warning(
                                    "[run-manifest.docking-stage] failed to record LeDock end pdb=%s variant=%s ph=%s stage=%s",
                                    paths.pdb_id,
                                    variant_label,
                                    ph_label if ph_label is not None else "base",
                                    ledock_stage_name,
                                    exc_info=True,
                                )
                        if bool(cfg.get("CHECKPOINT_ENABLE", True)) and ledock_fp is not None:
                            try:
                                if ledock_completion.get("success", False):
                                    checkpoint_mark_done(
                                        cfg,
                                        paths.pdb_id,
                                        ledock_stage_name,
                                        ledock_fp,
                                        ph_label=ph_label,
                                        variant=variant_env or None,
                                    )
                                else:
                                    logger.warning(
                                        "[Checkpoint] defer mark ledock_stage=%s pdb=%s reason=completion_missing",
                                        ledock_stage_name,
                                        paths.pdb_id,
                                    )
                            except Exception:
                                pass
                        logger.info(
                            "[ledock.done] pdb=%s stage=%s variant=%s ph=%s valid=%d invalid=%d elapsed_sec=%.2f",
                            paths.pdb_id,
                            stage_name,
                            variant_label,
                            ph_label if ph_label else "base",
                            valid_count,
                            invalid_count,
                            time.time() - ledock_start_ts,
                        )
                    except Exception as e:
                        logger.warning(
                            "[ledock.error] pdb=%s stage=%s variant=%s ph=%s reason=%s",
                            paths.pdb_id,
                            stage_name,
                            variant_label,
                            ph_label if ph_label else "base",
                            e,
                            exc_info=True,
                        )
                        ledock_metrics_by_stage[stage_name] = {}
                        if manifest_run_id:
                            try:
                                update_manifest_for_docking_stage(
                                    cfg,
                                    manifest_run_id,
                                    paths.pdb_id,
                                    variant_label,
                                    ledock_stage_name,
                                    status="failed",
                                    ph_tag=ph_label,
                                    error=str(e),
                                    elapsed_sec=time.time() - ledock_start_ts,
                                )
                            except Exception:
                                logger.warning(
                                    "[run-manifest.docking-stage] failed to record LeDock failure pdb=%s variant=%s ph=%s stage=%s",
                                    paths.pdb_id,
                                    variant_label,
                                    ph_label if ph_label else "base",
                                    ledock_stage_name,
                                    exc_info=True,
                                )
            else:
                logger.info(
                    "[ledock.scheduler] pdb=%s variant=%s ph=%s action=skip reason=no_jobs",
                    paths.pdb_id,
                    variant_label,
                    ph_label if ph_label else "base",
                )

            if ledock_metrics_by_stage:
                write_ledock_scores_csv(
                    cfg,
                    paths.pdb_id,
                    ledock_metrics_by_stage,
                    ph_label=ph_label,
                    variant=variant_env or None,
                    csv_prefix=csv_prefix,
                )

            if dock6_metrics_by_stage:
                write_dock6_scores_csv(
                    cfg,
                    paths.pdb_id,
                    dock6_metrics_by_stage,
                    ph_label=ph_label,
                    variant=variant_env or None,
                    csv_prefix=csv_prefix,
                )

            final_pose_validation_and_screenshots(
                cfg, paths.pdb_id, stages_for_run, receptor_pdbqt, center, validated_ligands_last,
                score_history, cleaned_pdb, docking_mode, logger, ph_label
            )

            csv_path = write_scores_csv(
                cfg,
                paths.pdb_id,
                score_history,
                ph_label=ph_label,
                variant=variant_env or None,
                csv_prefix=csv_prefix,
            )
            logger.info(
                "[Scores] ph_label=%s summary=%s",
                ph_label if ph_label else "base",
                csv_path,
            )

            # Optionally annotate FDA long CSV with T-scores vs DUD decoys.
            try:
                test_mode_now = _resolve_test_mode(cfg)
            except Exception:
                test_mode_now = None

            if run_mode == "fda" and test_mode_now in ("fda+dud", "fda+dud+hmdb"):
                try:
                    annotate_fda_long_csv_with_t_scores_vs_decoys(
                        cfg,
                        paths.pdb_id,
                        ph_label=ph_label,
                        logger=logger,
                    )
                except Exception as e:
                    logger.warning(
                        "[t-score.warn] pdb_id=%s ph=%s reason=%s",
                        paths.pdb_id,
                        ph_label if ph_label else "base",
                        e,
                    )
                try:
                    annotate_gnina_fda_long_csv_with_t_scores_vs_decoys(
                        cfg,
                        paths.pdb_id,
                        ph_label=ph_label,
                        logger=logger,
                    )
                except Exception as e:
                    logger.warning(
                        "[gnina.t-score.warn] pdb_id=%s ph=%s reason=%s",
                        paths.pdb_id,
                        ph_label if ph_label else "base",
                        e,
                    )
                try:
                    annotate_ledock_fda_long_csv_with_t_scores_vs_decoys(
                        cfg,
                        paths.pdb_id,
                        ph_label=ph_label,
                        logger=logger,
                    )
                except Exception as e:
                    logger.warning(
                        "[ledock.t-score.warn] pdb_id=%s ph=%s reason=%s",
                        paths.pdb_id,
                        ph_label if ph_label else "base",
                        e,
                    )

            try:
                compute_consensus_for_variant_ph(
                    cfg=cfg,
                    paths=paths,
                    ph_label=ph_label,
                    variant_env=variant_env,
                    variant_label=variant_label,
                    csv_prefix=csv_prefix,
                    logger=logger,
                )
            except Exception:
                logger.exception(
                    "[consensus.error] Failed to compute consensus scores for pdb_id=%s variant=%s ph=%s",
                    paths.pdb_id,
                    variant_label,
                    ph_label,
                )

            if _use_ledock(cfg):
                try:
                    logger.info(
                        "[ledock.mol2] starting_mol2_prep pdb=%s variant=%s ph=%s",
                        paths.pdb_id,
                        variant_label,
                        ph_label if ph_label else "base",
                    )
                    ligands_for_mol2 = _norm_dedupe(stage1_original + list(controls_for_run))
                    ensure_mol2_for_ledock(cfg, ligands_for_mol2, logger)
                    logger.info(
                        "[ledock.mol2] completed_mol2_prep pdb=%s variant=%s ph=%s",
                        paths.pdb_id,
                        variant_label,
                        ph_label if ph_label else "base",
                    )
                except Exception as e:
                    logger.warning(
                        "[ledock.mol2.warn] pdb=%s variant=%s ph=%s reason=%s",
                        paths.pdb_id,
                        variant_label,
                        ph_label if ph_label else "base",
                        e,
                        exc_info=True,
                    )
            else:
                logger.info(
                    "[ledock.mol2] skip_mol2_prep pdb=%s variant=%s ph=%s reason=use_ledock_disabled",
                    paths.pdb_id,
                    variant_label,
                    ph_label if ph_label else "base",
                )

            difficulty_info = None
            should_eval_difficulty = False
            if test_mode_now:
                should_eval_difficulty = "dud" in str(test_mode_now).lower()
            if should_eval_difficulty and run_mode not in (None, "dud"):
                should_eval_difficulty = False

            if should_eval_difficulty:
                logger.info(
                    "[difficulty] pdb=%s action=skip reason=druggability_orchestrator_enabled",
                    paths.pdb_id,
                )

            try:
                summary = {
                    "pdb_id": paths.pdb_id,
                    "center": tuple(map(float, center)) if center else None,
                    "box_size": tuple(map(float, box_size)) if box_size else None,
                    "n_ligands_stage1": len(stage1_original),
                    "n_valid_last_stage": len(validated_ligands_last),
                    "switch_history": getattr(selector, "switch_history", []),
                    "global_switches": guard.global_switches,
                    "stages": [s["name"] for s in stages_for_run],
                    "ph_label": ph_label,
                }
                if difficulty_info:
                    summary["difficulty"] = difficulty_info.difficulty
                    summary["roc_auc"] = difficulty_info.roc_auc
                    summary["difficulty_N"] = difficulty_info.N
                    summary["difficulty_n_actives"] = difficulty_info.n_actives
                _write_audit_json(cfg, paths.pdb_id, summary, ph_label=ph_label, variant=variant_env or None)
            except Exception as _e:
                logger.warning(f"Audit JSON write failed: {_e}")

            try:
                update_manifest_for_protein_success(
                    cfg,
                    manifest_run_id or "",
                    paths.pdb_id,
                    variant_label,
                    time.time() - ph_start_ts,
                    ph_tag=ph_label,
                )
            except Exception:
                ph_log.warning(
                    "[run-manifest.skip] pdb=%s variant=%s ph=%s reason=success",
                    paths.pdb_id,
                    variant_label,
                    ph_label if ph_label else "base",
                    exc_info=True,
                )
        except Exception:
            try:
                update_manifest_for_protein_failure(
                    cfg,
                    manifest_run_id or "",
                    paths.pdb_id,
                    variant_label,
                    None,
                    ph_tag=ph_label,
                )
            except Exception:
                ph_log.warning(
                    "[run-manifest.skip] pdb=%s variant=%s ph=%s reason=fail",
                    paths.pdb_id,
                    variant_label,
                    ph_label if ph_label else "base",
                    exc_info=True,
                )
            raise
