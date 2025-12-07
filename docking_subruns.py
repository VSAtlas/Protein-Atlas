from __future__ import annotations

import logging
import os
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
from path_router import Paths, docked_dir, receptor_file
from ph_ensemble_docking import (
    enumerate_ligands_for_ph_context,
    init_ph_tags_and_manifest,
    prewarm_ph_ligand_microstates,
)
from run_manifest import (
    update_manifest_for_protein_failure,
    update_manifest_for_protein_start,
    update_manifest_for_protein_success,
)
from record_data import record_le, write_scores_csv
from docking_ligands import (
    _count_heavy_atoms_from_pdbqt,
    _lib_roots_for_pdb,
    _resolve_test_mode,
    prepare_and_filter_ligands,
    select_ligands_for_next,
)
from single_ligand_index import (
    _ensure_single_ligand_index,
    _load_fda_name_map,
    _resolve_single_ligand,
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

        lig_set = {norm(x) for x in ligands}
        missing_controls = [p for p in prepped_control_pdbqts if norm(p) not in lig_set]

        if missing_controls:
            logger.info(f"[Controls] Adding {len(missing_controls)} prepared control(s) to Stage1.")
            ligands = [str(p) for p in missing_controls] + ligands
            for p in missing_controls:
                try:
                    heavy_atom_counts.setdefault(str(p), _count_heavy_atoms_from_pdbqt(p))
                except Exception:
                    heavy_atom_counts.setdefault(str(p), 0)

    def _norm_dedupe(seq):
        seen = set()
        out = []
        for p in seq:
            pn = norm(p)
            if pn not in seen:
                seen.add(pn)
                out.append(pn)
        return out

    ligands = _norm_dedupe(ligands)

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
            else:
                ligands = base_ligands[:]
                heavy_atom_counts = dict(base_heavy_atoms)
                pains_flags = dict(base_pains_flags)

            ctrl_stems_lower = {s.lower() for s in control_stems}
            ctrl_blacklist = {t.strip().upper() for t in str(cfg.get("CONTROL_BLACKLIST", "")).split(",") if t.strip()}
            min_ha = int(cfg.get("CONTROL_MIN_HEAVY_ATOMS", 10))

            def _is_control_path(p: str) -> bool:
                stem = Path(p).stem.split("_stage")[0]
                if stem.upper() in ctrl_blacklist:
                    return False
                if stem.lower() not in ctrl_stems_lower:
                    return False
                ha = heavy_atom_counts.get(p)
                return (ha is None) or (ha >= min_ha)

            ctrls = [p for p in ligands if _is_control_path(p)]
            non_ctrls = [p for p in ligands if not _is_control_path(p)]
            if ctrls:
                ligands = ctrls + non_ctrls
                logger.info(
                    f"[Controls] Front-loading {len(ctrls)} controls. "
                    f"First wave: {[Path(x).name for x in ligands[:int(cfg.get('MAX_PARALLEL_JOBS', 1))]]}"
                )

            present_ctrls = [
                Path(l).stem.split("_stage")[0].lower()
                for l in ligands
                if Path(l).stem.split("_stage")[0].lower() in ctrl_stems_lower
            ]

            if not present_ctrls:
                logger.warning(
                    "[Controls] No control ligands present in Stage1 ligand list -- "
                    "self-RMSD/locking will not be possible. (Check prep errors above.)"
                )
            if not ligands:
                logger.warning("No valid ligands after filtering; skipping protein.")
                continue

            selector = CenterSelector(cfg, logger, control_stems, heavy_atom_counts, center)
            guard = GlobalCenterGuard(
                max_global_switches=int(cfg.get("MAX_GLOBAL_CENTER_SWITCHES", 2))
            )

            stage1_original = ligands[:]

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

            control_stems_lower = {s.lower() for s in control_stems}
            forced_extracted_for_stage3 = {
                lig for lig in stage1_original
                if Path(lig).stem.split("_stage")[0].lower() in control_stems_lower
            }
            logger.info(
                f"[Force-carry] Extracted ligands earmarked for Stage3: {len(forced_extracted_for_stage3)}"
            )

            score_history: Dict[str, Dict[str, Dict]] = defaultdict(dict)
            validated_ligands_last: List[str] = []
            recenter_attempts = 0
            docking_mode = cfg.get("DOCKING_MODE", "discovery").lower()

            retry_mgr = RetryManager()

            i = 0
            while i < len(stages_for_run):
                guard.reset_stage()
                stage = stages_for_run[i]

                if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                    fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
                    if checkpoint_should_skip(
                        cfg,
                        paths.pdb_id,
                        stage["name"],
                        fp,
                        ph_label=ph_label,
                        variant=variant_env or None,
                    ):
                        logger.info(f"[Checkpoint] Skipping {stage['name']} (fingerprint matched).")
                        i += 1
                        continue

                if not ligands:
                    logger.warning(f"No ligands to dock at {stage['name']}; stopping for this protein.")
                    break

                logger.info(f"Starting {stage['name']} with {len(ligands)} ligands...")
                if ph_label:
                    stage_dir = paths.docked_stage_dir(variant_env or None, stage["name"], ph_label)
                    ph_log.info(
                        "[ph_ensemble.dock.stage] pdb_id=%s variant=%s ph=%s stage=%s receptor=%s out=%s",
                        paths.pdb_id,
                        variant_label,
                        ph_label,
                        stage["name"],
                        receptor_pdbqt,
                        str(stage_dir),
                    )

                if i == 0 and ctrls and non_ctrls:
                    logger.info(
                        f"Stage1 two-wave: {len(ctrls)} controls first, then {len(non_ctrls)} others."
                    )

                    s1, v1, d1, rd1, inv1 = run_one_stage(
                        cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage,
                        ctrls, logger, retry_mgr, control_lookup, ph_label=ph_label
                    )

                    try:
                        dec = selector.consider_switch(stage['name'], s1, v1, rd1, receptor_pdbqt, center, guard)
                        if dec.promoted and dec.new_center is not None:
                            old = center
                            center = dec.new_center
                            guard.mark_switch()
                            logger.info(
                                f"[CENTER] Switched before library run: {old} -> {center} ({dec.reason}) [global switch]"
                            )
                    except Exception as e:
                        logger.warning(f"CenterSelector (controls-only) failed gracefully: {e}")

                    lock_score_max = float(cfg.get("CONTROL_LOCK_SCORE_MAX", -6.0))
                    lock_min_hits = int(cfg.get("CONTROL_LOCK_MIN_HITS", 1))
                    lock_center_max = float(cfg.get("CONTROL_LOCK_CENTER_MAX_DIST", 4.0))
                    qualified_controls = []

                    for lig in v1:
                        stem = Path(lig).stem.split("_stage")[0].lower()
                        ha = heavy_atom_counts.get(lig)
                        if stem in control_stems_lower and (ha is None or ha >= min_ha):
                            sc = s1.get(lig)
                            if sc is not None and np.isfinite(sc) and sc <= lock_score_max:
                                pose_path = rd1.get(lig)
                                c = CenterSelector._pdbqt_centroid(pose_path) if pose_path else None
                                if c is not None and np.linalg.norm(c - np.array(center, float)) <= lock_center_max:
                                    qualified_controls.append(lig)

                    if len(qualified_controls) >= lock_min_hits and not guard.locked:
                        guard.lock()
                        logger.info(
                            "[CONTROL-LOCK] Early lock from controls-only wave "
                            f"(n={len(qualified_controls)}, score<={lock_score_max}, dist<={lock_center_max} A); "
                            "future center switches disabled."
                        )

                    s2, v2, d2, rd2, inv2 = run_one_stage(
                        cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage,
                        non_ctrls, logger, retry_mgr, control_lookup, ph_label=ph_label
                    )

                    scores, validated, distances = ({**s1, **s2}, v1 + v2, d1 + d2)
                    raw_docked = {**rd1, **rd2}
                    invalids = {**inv1, **inv2}
                else:
                    scores, validated, distances, raw_docked, invalids = run_one_stage(
                        cfg, paths.pdb_id, receptor_pdbqt, center, box_size, stage,
                        ligands, logger, retry_mgr, control_lookup, ph_label=ph_label
                    )

                validated_ligands_last = validated

                def _is_control(lig: str) -> bool:
                    stem = Path(lig).stem.split("_stage")[0].lower()
                    if stem.upper() in {s.strip().upper() for s in cfg.get("CONTROL_BLACKLIST", "").split(",") if s.strip()}:
                        return False
                    ha = heavy_atom_counts.get(lig)
                    if ha is not None and ha < int(cfg.get("CONTROL_MIN_HEAVY_ATOMS", 10)):
                        return False
                    return stem in {s.lower() for s in control_stems}

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
                    processed = {norm(x) for x in ligands}
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
                            i = 0
                            continue
                        else:
                            break

                    use_stage1_base = (docking_mode == "polypharmacology" and i == 1)

                    rescue = []
                    if i < len(stages_for_run) - 1:
                        for lig, (sc, reason) in invalids.items():
                            if sc is not None and "self_rmsd_" in str(reason).lower() and sc <= float(
                                    cfg.get("RESCUE_SELF_RMSD_SCORE_MAX", -8.0)):
                                rescue.append((sc, lig))
                        rescue = [lig for _, lig in sorted(rescue)[:int(cfg.get("RESCUE_SELF_RMSD_TOP_N", 10))]]

                    selected = select_ligands_for_next(
                        docking_mode,
                        i,
                        stages_for_run,
                        scores,
                        logger,
                        base_pool_n=(len(stage1_original) if use_stage1_base else None),
                        force_include=(forced_extracted_for_stage3 if use_stage1_base else None)
                    )

                    if rescue:
                        sel_set = set(selected)
                        rescue_unique = [r for r in rescue if r not in sel_set]
                        ligands = rescue_unique + selected
                    else:
                        ligands = selected
                    if not ligands:
                        logger.warning(f"No ligands selected for {stages_for_run[i + 1]['name']}; stopping.")
                        break

                if bool(cfg.get("CHECKPOINT_ENABLE", True)):
                    try:
                        fp = _fingerprint_stage(cfg, receptor_pdbqt, center, box_size, stage)
                        checkpoint_mark_done(
                            cfg,
                            paths.pdb_id,
                            stage["name"],
                            fp,
                            ph_label=ph_label,
                            variant=variant_env or None,
                        )
                    except Exception:
                        pass

                i += 1

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
