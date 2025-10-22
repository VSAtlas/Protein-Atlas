
def general_fallback_recenter_if_needed(cfg, pdb_id, stage_name, scores, raw_docked_ligands,
                                        receptor_pdbqt, center, box_size, logger, ligands_stage1_original):
    """
    If a stage yields no valid ligands, try fallback recentering and redo stage1.
    Returns (should_restart_stage1, new_center, new_box_size, ligands_to_redock)
    """
    if scores:
        return False, center, box_size, None

    logger.warning(f"No valid ligands in {stage_name}. Attempting fallback recentering...")
    fb_pose, new_center, best_score, chosen_ligand = attempt_fallback_recenter(
        fallback_ligands=raw_docked_ligands,
        receptor_pdbqt=receptor_pdbqt,
        docking_dir=os.path.join(cfg["DOCKED_DIR"], pdb_id),
        stage_name=stage_name,
        pocket_center=center,
        logger=logger,
        exclude_basenames=set(),
    )
    if not fb_pose or new_center is None:
        logger.warning("Fallback recovery failed. Ending docking for this protein.")
        return False, center, box_size, None

    max_box = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
    center = new_center
    box_size = tuple(min(max_box, s) for s in box_size)
    logger.info("Re-running stage1 with new center after no-valid fallback.")
    return True, center, box_size, ligands_stage1_original[:]


def early_recenter(stage_index, all_distances, scores, recenter_knobs, box_size, center,
                         ligands_stage1_original, tried_centers, recenter_attempts, max_recenter_attempts,
                         raw_docked_ligands, cfg, pdb_id, receptor_pdbqt, logger):
    """
    Stage-1 heuristic: decide whether to expand/recenter box early when results are systematically off.
    Returns (should_redo_stage1, new_center, new_box_size, recenter_attempts, tried_centers, ligands_to_redock)
    """
    if stage_index != 0:
        return False, center, box_size, recenter_attempts, tried_centers, None

    evaluated = len(all_distances)
    valid_count = len(scores)
    if evaluated < recenter_knobs["EARLY_RECENTER_MIN_EVAL"]:
        logger.info(f"Early recenter skipped: evaluated={evaluated} < {recenter_knobs['EARLY_RECENTER_MIN_EVAL']}.")
        return False, center, box_size, recenter_attempts, tried_centers, None

    far = sum(1 for d in all_distances if isinstance(d, (int, float)) and d > recenter_knobs["EARLY_RECENTER_FAR_A"])
    far_ratio = far / evaluated if evaluated else 0.0
    med_dist = float(np.median(all_distances)) if all_distances else 0.0

    # One-time box expansion for borderline cases
    if (recenter_knobs["ALLOW_BOX_EXPAND"] and 0.5 <= far_ratio < recenter_knobs["EARLY_RECENTER_RATIO"]
            and 8.0 <= med_dist < recenter_knobs["EARLY_RECENTER_MEDIAN_A"] and valid_count == 0):
        max_box = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
        new_box = tuple(min(max_box, s + 4.0) for s in box_size)
        if new_box != box_size:
            logger.info(f"Borderline far_ratio={far_ratio:.2f}, median={med_dist:.1f} Å → expand box to {new_box} and redo stage1.")
            return True, center, new_box, recenter_attempts, tried_centers, ligands_stage1_original[:]

    # Recenter when clearly off and no valid poses
    if (far_ratio >= recenter_knobs["EARLY_RECENTER_RATIO"]) and (med_dist >= recenter_knobs["EARLY_RECENTER_MEDIAN_A"]) and (valid_count == 0):
        logger.warning(f"Early recenter trigger: far_ratio={far_ratio:.2f}, median={med_dist:.1f} Å, valid=0 → recentering.")
        remaining = max(0, max_recenter_attempts - recenter_attempts)
        while remaining > 0:
            fb_pose, new_center, best_score, chosen_ligand = attempt_fallback_recenter(
                fallback_ligands=raw_docked_ligands,
                receptor_pdbqt=receptor_pdbqt,
                docking_dir=os.path.join(cfg['DOCKED_DIR'], pdb_id),
                stage_name="stage1",
                pocket_center=center,
                logger=logger,
                exclude_basenames=set()
            )
            if new_center is None:
                logger.warning("Fallback couldn’t produce a new center from remaining candidates.")
                break

            center_key = tuple(round(c, 1) for c in new_center)
            if center_key in tried_centers:
                logger.warning("Proposed center equals a previously tried center; trying next candidate...")
                remaining -= 1
                continue

            tried_centers.add(center_key)
            recenter_attempts += 1
            center = new_center
            max_box = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
            box_size = tuple(min(max_box, s) for s in box_size)
            logger.info("Re-running stage1 with new center and tightened box.")
            return True, center, box_size, recenter_attempts, tried_centers, ligands_stage1_original[:]

        logger.warning("Early recenter did not yield a new center; proceeding without recenter.")
    return False, center, box_size, recenter_attempts, tried_centers, None