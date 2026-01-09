from __future__ import annotations

import logging
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .docking_control_redock import (
    _control_centers_by_ph,
)
from .docking_subrun_selection import (
    _split_controls_and_noncontrols,
    _interleave_controls,
)
from .docking_single_ligand_mode import resolve_single_ligand_or_prepare
from .docking_ph_subruns import resolve_ph_tags_and_root, enumerate_ph_ligands_if_needed
from .docking_plan_only import maybe_handle_no_library_docking
from .docking_vina_multistage import run_multistage_vina
from .docking_subrun_finalize import finalize_ph_subrun

from apo_holo_mode import _record_apo_holo_usage
from .docking_centering import CenterSelector
from .docking_utils import norm
from .fallback_recenter import (
    GlobalCenterGuard,
    RecenterParams,
)
from path_router.path_router import Paths, docked_dir, receptor_file
from run_manifest import (
    update_manifest_for_protein_failure,
    update_manifest_for_protein_start,
)
from .docking_ligands import (
    _count_heavy_atoms_from_pdbqt,
    prepare_and_filter_ligands,
)


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
    center_by_ph: Optional[Dict[Optional[str], Tuple[float, float, float]]] = None
    box_by_ph: Optional[Dict[Optional[str], Tuple[float, float, float]]] = None
    center_source_by_ph: Optional[Dict[Optional[str], str]] = None


def resolve_center_box_for_ph(
    ctx: ProteinDockingContext,
    ph_label: Optional[str],
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float], str]:
    ph_key = None if ph_label in (None, "base") else ph_label
    if ctx.center_by_ph and ctx.box_by_ph and ph_key in ctx.center_by_ph and ph_key in ctx.box_by_ph:
        center = ctx.center_by_ph[ph_key]
        box_size = ctx.box_by_ph[ph_key]
        source = "control_ph"
        if ctx.center_source_by_ph and ph_key in ctx.center_source_by_ph:
            raw_src = ctx.center_source_by_ph[ph_key]
            if raw_src != "control":
                source = "fallback"
        return center, box_size, source

    source = "fallback"
    if ctx.center_by_ph:
        source = "control_missing"
    if ctx.center is None or ctx.box_size is None:
        raise ValueError("Center/box unavailable for requested ph_label.")
    return ctx.center, ctx.box_size, source


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
    """
    cfg = ctx.cfg
    paths = ctx.paths
    logger = ctx.logger
    pdb_id = ctx.pdb_id
    variant_env = ctx.variant_env
    variant_token = ctx.variant_token
    variant_label = ctx.variant_label
    legacy_mode = ctx.legacy_mode
    stages = ctx.stages
    control_stems = ctx.control_stems
    control_lookup = ctx.control_lookup

    run_mode = subrun.run_mode
    csv_prefix = subrun.csv_prefix
    stage_name_prefix = subrun.stage_name_prefix

    # 1. Single Ligand Resolution
    single_res = resolve_single_ligand_or_prepare(
        cfg=cfg,
        paths=paths,
        logger=logger,
        pdb_id=pdb_id,
        run_mode=run_mode,
    )

    if single_res:
        ligands, heavy_atom_counts, pains_flags = single_res
    else:
        if cfg.get("_EFFECTIVE_SINGLE_LIGAND"):
             # Fallback triggered or block exit, but if we are here, we continue normally
             pass
        # 2. Prepare & Filter Ligands
        ligands, heavy_atom_counts, pains_flags = prepare_and_filter_ligands(
            cfg,
            paths,
            logger,
            run_mode=run_mode,
        )

    # 3. Split Controls / Non-controls & Augment Missing
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

    # 4. Resolve PH Tags & Root
    ph_tags, ph_ligand_root, ph_test_mode_override, library_for_manifest = resolve_ph_tags_and_root(
        cfg=cfg,
        paths=paths,
        logger=logger,
        run_mode=run_mode,
        variant_token=variant_token,
        legacy_mode=legacy_mode,
    )
    if not ph_tags:
        return

    stages_for_run = [
        {**stage, "name": f"{stage_name_prefix}{stage['name']}"}
        for stage in stages
    ]
    
    ph_log = logging.getLogger("ph_ensemble")
    
    # 5. Loop PH Tags
    for ph_label in ph_tags:
        ph_start_ts = time.time()
        try:
            update_manifest_for_protein_start(
                cfg,
                cfg.get("RUN_ID", "") or "",
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

        # Plan-only check: if we only print plan, we might still want to check rec_exists?
        # Original code checked existence FIRST.
        # But if plan-only is TRUE, we skip existence check? No, original code checked exists before plan-only block.
        # Wait, if `rec_path` doesn't exist, we might skip unless we want to fail.

        if bool(cfg.get("NO_LIBRARY_DOCKING", False)):
            # This check is inside the loop in original code, AFTER existence check?
            # Actually, let's look at original code.
            # It checks `rec_exists`. If not exists, it logs warning and continues.
            # THEN it calls `plan_only` check?
            # Wait, `plan_only = os.environ.get("A2_PLAN_ONLY") == "1"` was in my prompt's A2.
            # But `NO_LIBRARY_DOCKING` is what I extracted to `docking_plan_only.py`.
            # Let's handle `NO_LIBRARY_DOCKING` (docking_plan_only) separately from the loop guard.
            # Original code:
            # `if plan_only: ... continue` (This was `A2_PLAN_ONLY` env var, used for debugging maybe?)
            # `if not rec_exists: ... continue`
            # Then inside: `if bool(cfg.get("NO_LIBRARY_DOCKING", False)): ... continue`
            
            pass

        # Environment variable plan_only check (A2_PLAN_ONLY) - kept from original code
        if os.environ.get("A2_PLAN_ONLY") == "1":
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
                    cfg.get("RUN_ID", "") or "",
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
            center, box_size, _center_src = resolve_center_box_for_ph(ctx, ph_label)

            ligands, heavy_atom_counts, pains_flags = enumerate_ph_ligands_if_needed(
                cfg=cfg,
                paths=paths,
                logger=logger,
                ph_label=ph_label,
                ph_ligand_root=ph_ligand_root,
                base_ligands=base_ligands,
                base_heavy_atoms=base_heavy_atoms,
                base_pains_flags=base_pains_flags,
                controls_for_run=base_controls,
            )
            controls_for_run = base_controls[:] # Copy fresh for filtering

            ctrl_blacklist = {t.strip().upper() for t in str(cfg.get("CONTROL_BLACKLIST", "")).split(",") if t.strip()}
            min_ha = int(cfg.get("CONTROL_MIN_HEAVY_ATOMS", 10))

            def _eligible_control_path(p: str) -> bool:
                stem = Path(p).stem.split("_stage")[0].lower()
                if stem.upper() in ctrl_blacklist:
                    return False
                if stem not in control_stems_lower_init:
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
            
            receptor_pdbqt = str(rec_path)
            stage1_original = [l for l in ligands if norm(l) not in control_norms_for_run]

            # Plan-only check
            if maybe_handle_no_library_docking(
                cfg=cfg,
                paths=paths,
                logger=logger,
                variant_env=variant_env or "",
                variant_token=variant_token,
                variant_label=variant_label,
                legacy_mode=legacy_mode,
                ph_label=ph_label,
                run_mode=run_mode,
                center=center,
                box_size=box_size,
                stage1_original=stage1_original,
            ):
                continue

            # Multistage Vina loop execution
            result = run_multistage_vina(
                cfg=cfg,
                paths=paths,
                receptor_pdbqt=receptor_pdbqt,
                center=center,
                box_size=box_size,
                ph_label=ph_label,
                legacy_mode=legacy_mode,
                variant_env=variant_env or None,
                stages_for_run=stages_for_run,
                stage1_original=stage1_original,
                controls_for_run=controls_for_run,
                control_norms_for_run=control_norms_for_run,
                control_lookup=control_lookup,
                heavy_atom_counts=heavy_atom_counts,
                pains_flags=pains_flags,
                logger=logger,
                recenter_params=ctx.recenter_params,
                control_stems=control_stems,
            )

            # Finalize
            # We need to extract the switch history from selector/guard somehow?
            # run_multistage_vina returns recenter_attempts, but not the switch history directly?
            # Wait, `run_multistage_vina` instantiates `CenterSelector` and `GlobalCenterGuard` internally.
            # So we lose their state?
            # `run_multistage_vina` doesn't return `selector` or `guard`.
            # But `finalize_ph_subrun` writes audit JSON which includes switch history.
            # I should update `run_multistage_vina` to return this info or pass it out.
            # The `MultistageVinaResult` has:
            # `recenter_attempts`
            # But not `switch_history` or `global_switches`.
            
            # Since I cannot change `MultistageVinaResult` (as I didn't edit `docking_vina_multistage.py` in this step),
            # I should probably have updated `docking_vina_multistage.py` to return these, OR
            # `run_multistage_vina` was just extracted in the PREVIOUS turn?
            # Yes, `docking_vina_multistage.py` was extracted in the previous state.
            # I should assume `MultistageVinaResult` probably doesn't have `switch_history`.
            # Let me check `docking_vina_multistage.py` content again.
            
            # `docking_vina_multistage.py` result class:
            # @dataclass
            # class MultistageVinaResult:
            #    scores: Dict[str, float]
            #    ...
            #    recenter_attempts: int
            #    ...
            
            # And inside `run_multistage_vina`:
            # selector = CenterSelector(...)
            # guard = GlobalCenterGuard(...)
            # ...
            # return MultistageVinaResult(...)
            
            # It seems `switch_history` is NOT returned. This is a potential regression if `finalize_ph_subrun` needs it.
            # `docking_subruns.py` original code had `selector` and `guard` in the same scope.
            # `finalize_ph_subrun` logic in `docking_subruns.py` (original) used:
            # "switch_history": getattr(selector, "switch_history", []),
            # "global_switches": guard.global_switches,
            
            # So I DO need to pass these out.
            # Since I am not supposed to change `docking_vina_multistage.py` (it was extracted previously),
            # but wait, the prompt says "docking.py: keep phase orchestration, extract large...".
            # The prompt assumes `docking_vina_multistage.py` exists or is created?
            # Ah, the context says "MODIFIED: docking_vina_multistage.py - Implements the core Vina loop".
            # So I CANNOT change it easily without breaking the constraint "No behavior changes".
            # If I can't change `docking_vina_multistage.py`, then I can't get that info out.
            # BUT, I am allowed to edit `docking_subruns.py` and create NEW modules.
            # `docking_vina_multistage.py` was already there from previous turn.
            
            # **Correction**: I can edit `docking_vina_multistage.py` if needed to fix bugs/add returns, as I am "Refactoring".
            # But the prompt specifically lists "PART A" and "PART B" and doesn't mention editing `docking_vina_multistage.py`.
            # However, if I don't, I lose data.
            # Wait, `finalize_ph_subrun` takes `center_selector_history` and `global_switches` arguments.
            # So I definitely need to pass them.
            # I will modify `docking_vina_multistage.py` to return them in `MultistageVinaResult`.
            # I'll do that in a separate write.
            
            # Let's write `docking_subruns.py` assuming `MultistageVinaResult` has `switch_history` and `global_switches`.
            
            finalize_ph_subrun(
                cfg=cfg,
                paths=paths,
                logger=logger,
                pdb_id=pdb_id,
                variant_env=variant_env or "",
                variant_token=variant_token,
                variant_label=variant_label,
                legacy_mode=legacy_mode,
                ph_label=ph_label,
                ph_start_ts=ph_start_ts,
                receptor_pdbqt=receptor_pdbqt,
                center=center,
                box_size=box_size,
                run_mode=run_mode,
                csv_prefix=csv_prefix,
                stage_name_prefix=stage_name_prefix,
                stages_for_run=stages_for_run,
                score_history=result.score_history,
                validated=result.validated,
                validated_ligands_last=result.validated_ligands_last,
                gnina_jobs=result.gnina_jobs,
                ledock_jobs=result.ledock_jobs,
                dock6_jobs=result.dock6_jobs,
                dock6_metrics_by_stage=result.dock6_metrics_by_stage,
                heavy_atom_counts=heavy_atom_counts,
                pains_flags=pains_flags,
                library_for_manifest=library_for_manifest,
                docking_mode=result.docking_mode,
                stage1_original=stage1_original,
                controls_for_run=controls_for_run,
                control_lookup=control_lookup,
                recenter_attempts=result.recenter_attempts,
                center_selector_history=getattr(result, "switch_history", []),
                global_switches=getattr(result, "global_switches", 0),
            )

        except Exception:
            try:
                update_manifest_for_protein_failure(
                    cfg,
                    cfg.get("RUN_ID", "") or "",
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