"""Protein cleaning pipeline orchestration helpers."""

from __future__ import annotations

import logging
import json
import os
from pathlib import Path
import re
import shutil
from typing import Any, Optional, Union

from protein_prep.automate_protein_prep import (
    _cfg_env_or_default,
    _helium_postwrite_counter,
    _set_clean_provenance,
    assert_no_metal_in_peptidic,
    build_missing_loops,
    config,
    count_waters_within,
    file_contains_hydrogens,
    filter_invalid_chains,
    quick_element_histogram,
    run_phenix_pdbtools,
)
from path_router import make_paths  # type: ignore[attr-defined]
from protein_prep.pdb_fixer_runtime import (
    assert_no_helium_in_hydrogen_names,
    fix_element_columns_in_file,
    fix_pdb_elements,
    scan_helium_counts,
)
from protein_prep.aliases_policy import RULES, _flatten_semicolons
from protein_prep.altloc_filter import filter_altlocs
from protein_prep.network.apply import (
    apply_binding_site_network_recommendations,
)
from protein_prep.chain_prune import (
    _chains_to_keep,
    _prune_chains_conservative,
    prune_to_chains,
    select_chains_to_keep,
)
from protein_prep.hydrogen_cleanup import clean_hydrogens
from protein_prep.geometry.constructive import (
    apply_targeted_binding_site_geometry_policy,
)
from protein_prep.ion_audit import (
    _IonAuditManager,
    _collect_monoatomic_records,
    _diff_detail_records,
    _emit_ion_breadcrumb,
    _format_ion_hist,
    _log_ions_probe,
    _pop_ion_audit_manager,
    _push_ion_audit_manager,
    _reset_ion_probe,
)
from protein_prep.het_state_selection import select_het_states
from protein_prep.layout import canon_paths
from protein_prep.ligand_extract import extract_ligands_from_filtered
from protein_prep.metal.duplicate_sites import resolve_duplicate_metal_sites
from protein_prep.metal.retention import rescue_retained_hets
from protein_prep.openmm_repair import get_target_ph_for_prep, repair_with_pdbfixer
from protein_prep.prep_utils import _cfg_bool, _resolve_variant_token
from protein_prep.protonation import (
    REDUCE_EXE,
    _protonate_with_pdb2pqr_if_available,
    assign_protonation_states,
)
from protein_prep.strip_nsr import strip_nonstandard_residues
from protein_prep.water.policy import audit_water_policy, write_supported_water_receptor
from protein_prep.water.records import compute_control_centroids, filter_waters_near_points


_ION_ELEMENTS = {
    "LI",
    "NA",
    "K",
    "RB",
    "CS",
    "MG",
    "CA",
    "SR",
    "BA",
    "ZN",
    "MN",
    "FE",
    "CO",
    "NI",
    "CU",
    "AL",
    "CD",
    "HG",
    "PB",
    "AG",
    "AU",
    "PT",
    "MO",
    "RU",
    "RH",
    "PD",
    "CL",
    "BR",
    "I",
    "F",
    "AT",
    "LA",
    "CE",
    "PR",
    "ND",
    "PM",
    "SM",
    "EU",
    "GD",
    "TB",
    "DY",
    "HO",
    "ER",
    "TM",
    "YB",
    "LU",
    "U",
    "TH",
}


def _line_element_token(line: str) -> str:
    element = line[76:78].strip().upper()
    if element:
        return element
    guess_chars = []
    for char in line[12:16].strip():
        if char.isalpha():
            guess_chars.append(char)
        else:
            break
    guess = "".join(guess_chars).upper()
    if len(guess) >= 2 and guess[:2] in _ION_ELEMENTS:
        return guess[:2]
    return guess[:1]


def _count_ion_like_records(pdb_path: Union[str, Path]) -> int:
    path = Path(pdb_path)
    if not path.exists():
        return 0
    count = 0
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if (
                line.startswith(("ATOM  ", "HETATM"))
                and _line_element_token(line) in _ION_ELEMENTS
            ):
                count += 1
    return count


def _log_ion_count(
    log: logging.Logger,
    ordinal: int,
    stage: str,
    phase: str,
    pdb_path: Union[str, Path],
) -> int:
    count = _count_ion_like_records(pdb_path)
    log.info("(%d) %s %s: metals=%d file=%s", ordinal, stage, phase, count, pdb_path)
    return count


def _largest_ligand_pdb(ligand_dir: Path) -> Path | None:
    candidates: list[tuple[int, Path]] = []
    if not ligand_dir.exists():
        return None
    for path in ligand_dir.glob("*.pdb"):
        heavy_atoms = 0
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                element = (line[76:78].strip() or line[12:16].strip()[:1]).upper()
                if element != "H":
                    heavy_atoms += 1
        if heavy_atoms > 0:
            candidates.append((heavy_atoms, path))
    if not candidates:
        return None
    return max(candidates, key=lambda item: item[0])[1]


def _selected_water_policy_enabled() -> bool:
    policy = (
        os.environ.get("ATLAS_WATER_BENCH_POLICY")
        or config.get("WATER_POLICY", "site_only")
        or "site_only"
    ).strip().lower()
    return policy in {"site_only", "selected", "selected_water", "supported"}


def _summary_int(summary: dict[str, object], key: str) -> int:
    try:
        return int(str(summary.get(key, 0) or 0))
    except (TypeError, ValueError):
        return 0


def _apply_binding_site_network_stage(
    *,
    receptor_pdb: Path,
    filtered_pdb: Path,
    ligand_dir: Path,
    work_dir: Path,
    pdb_id: str,
    target_ph: float,
    log: logging.Logger,
) -> Path:
    if not _cfg_bool("APPLY_BINDING_SITE_NETWORK", True):
        log.info("[binding_site_network] stage=skip reason=config_disabled")
        return receptor_pdb

    ligand_pdb = _largest_ligand_pdb(ligand_dir)
    selection_sidecar = work_dir / f"{pdb_id}_het_state_selection.json"
    try:
        selection_summary, _ = select_het_states(
            filtered_pdb,
            work_dir,
            target_ph,
            sidecar_path=selection_sidecar,
        )
        log.info(
            "[binding_site_network] selection_status=%s sidecar=%s",
            selection_summary.get("het_state_selection_status", ""),
            selection_sidecar,
        )
    except Exception as exc:
        log.warning("[binding_site_network] selection_failed err=%s", exc)
        return receptor_pdb

    protein_applied = work_dir / f"{pdb_id}_binding_site_network_applied.pdb"
    protein_summary = apply_binding_site_network_recommendations(
        receptor_pdb,
        protein_applied,
        selection_sidecar=selection_sidecar,
        ligand_pdb=ligand_pdb,
        apply_protein_modes=True,
        apply_water_hydrogens=False,
        audit_path=work_dir / f"{pdb_id}_binding_site_network_application.json",
    )
    if _summary_int(protein_summary, "network_application_applied_count") > 0:
        shutil.copyfile(protein_applied, receptor_pdb)
        log.info(
            "[binding_site_network] protein_modes_applied=%s audit=%s",
            protein_summary.get("network_application_applied_count", 0),
            protein_summary.get("network_application_audit_path", ""),
        )
    else:
        log.info(
            "[binding_site_network] protein_modes_applied=0 skipped_low_confidence=%s",
            protein_summary.get("network_application_low_confidence_skipped_count", 0),
        )

    if ligand_pdb is None or not _selected_water_policy_enabled():
        return receptor_pdb

    try:
        _, water_rows = audit_water_policy(
            filtered_pdb,
            receptor_pdb,
            ligand_pdb,
            sidecar_path=work_dir / f"{pdb_id}_water_policy_audit.json",
        )
        selected_water_pdb = work_dir / f"{pdb_id}_selected_water_receptor.pdb"
        selection = write_supported_water_receptor(
            receptor_pdb,
            selected_water_pdb,
            water_rows,
            source_water_pdb=filtered_pdb,
        )
        (work_dir / f"{pdb_id}_selected_water_policy.json").write_text(
            json.dumps(
                {"selection": selection, "rows": water_rows},
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )
        if selected_water_pdb.exists():
            shutil.copyfile(selected_water_pdb, receptor_pdb)
        water_applied = work_dir / f"{pdb_id}_selected_water_network_applied.pdb"
        water_summary = apply_binding_site_network_recommendations(
            receptor_pdb,
            water_applied,
            selection_sidecar=selection_sidecar,
            ligand_pdb=ligand_pdb,
            water_policy_rows=water_rows,
            apply_protein_modes=False,
            apply_water_hydrogens=True,
            audit_path=work_dir / f"{pdb_id}_water_network_application.json",
        )
        if _summary_int(water_summary, "network_application_applied_count") > 0:
            shutil.copyfile(water_applied, receptor_pdb)
        log.info(
            "[binding_site_network] selected_waters=%s water_h_added=%s",
            selection.get("water_policy_selected_water_count", 0),
            water_summary.get("network_application_water_hydrogen_added_count", 0),
        )
    except Exception as exc:
        log.warning("[binding_site_network] water_application_failed err=%s", exc)
    return receptor_pdb


def _apply_binding_site_geometry_policy_stage(
    *,
    receptor_pdb: Path,
    reference_pdb: Path,
    ligand_dir: Path,
    work_dir: Path,
    pdb_id: str,
    target_ph: float,
    log: logging.Logger,
) -> Path:
    if not _cfg_bool("APPLY_BINDING_SITE_GEOMETRY_REPAIR", True):
        log.info("[binding_site_geometry] stage=skip reason=config_disabled")
        return receptor_pdb
    ligand_pdb = _largest_ligand_pdb(ligand_dir)
    if ligand_pdb is None:
        log.info("[binding_site_geometry] stage=skip reason=no_extracted_ligand")
        return receptor_pdb
    repaired = work_dir / f"{pdb_id}_binding_site_geometry_policy.pdb"
    try:
        summary = apply_targeted_binding_site_geometry_policy(
            receptor_pdb,
            repaired,
            reference_pdb=reference_pdb,
            ligand_pdb=ligand_pdb,
            work_dir=work_dir,
            target_ph=target_ph,
            audit_prefix=pdb_id,
            logger=log,
        )
    except Exception as exc:
        log.warning("[binding_site_geometry] policy_failed err=%s", exc)
        return receptor_pdb
    status = str(summary.get("binding_site_geometry_policy_status", ""))
    if repaired.exists() and status not in {"", "unchanged"}:
        shutil.copyfile(repaired, receptor_pdb)
    log.info(
        "[binding_site_geometry] status=%s residual_gate=%s pre_clashes=%s final_clashes=%s site_sidechain_atoms=%s audit=%s",
        status,
        summary.get("binding_site_geometry_residual_gate_status", ""),
        summary.get("binding_site_geometry_pre_clash_count", ""),
        summary.get("binding_site_geometry_final_clash_count", ""),
        summary.get("binding_site_geometry_site_sidechain_atom_count", ""),
        summary.get("binding_site_geometry_policy_audit_path", ""),
    )
    return receptor_pdb


def _apply_duplicate_metal_site_policy_stage(
    *,
    receptor_pdb: Path,
    work_dir: Path,
    pdb_id: str,
    log: logging.Logger,
) -> Path:
    if not _cfg_bool("APPLY_DUPLICATE_METAL_SITE_POLICY", True):
        log.info("[duplicate_metal_site] stage=skip reason=config_disabled")
        return receptor_pdb
    repaired = work_dir / f"{pdb_id}_duplicate_metal_sites_resolved.pdb"
    sidecar = work_dir / f"{pdb_id}_duplicate_metal_sites.json"
    summary = resolve_duplicate_metal_sites(
        receptor_pdb,
        repaired,
        sidecar_path=sidecar,
    )
    status = str(summary.get("duplicate_metal_site_status", ""))
    if bool(summary.get("duplicate_metal_site_fail_for_review")):
        log.error(
            "[duplicate_metal_site] status=%s action=fail_for_review clusters=%s audit=%s",
            status,
            summary.get("duplicate_metal_site_review_cluster_count", 0),
            sidecar,
        )
        raise RuntimeError(
            "duplicate_metal_site_review_required:"
            f"{summary.get('duplicate_metal_site_review_cluster_count', 0)}"
        )
    dropped_count = int(str(summary.get("duplicate_metal_site_dropped_atom_count", 0) or 0))
    if repaired.exists() and dropped_count:
        shutil.copyfile(repaired, receptor_pdb)
    log.info(
        "[duplicate_metal_site] status=%s clusters=%s dropped=%s audit=%s",
        status,
        summary.get("duplicate_metal_site_cluster_count", 0),
        summary.get("duplicate_metal_site_dropped_atom_count", 0),
        sidecar,
    )
    return receptor_pdb


def _run_clean_pdb_original(
    pdb_file: Union[str, Path],
    output_root: Union[str, Path],
    logger: Optional[logging.Logger] = None,
) -> Optional[str]:
    """Run the full cleaning pipeline and return path to final cleaned PDB (receptor)."""
    prepare_before_count: Optional[int] = None

    output_root = str(output_root)
    Path(output_root).mkdir(parents=True, exist_ok=True)

    raw_stem = os.path.splitext(os.path.basename(str(pdb_file)))[0]
    pdb_id = re.sub(r"(_nolig(_cleaned)?|_cleaned)$", "", raw_stem, flags=re.I).upper()
    _set_clean_provenance("automate_protein_prep.clean_pdb")
    _reset_ion_probe(pdb_id)
    log = logger or logging.getLogger(pdb_id)
    log.info("[prep.id] clean_pdb stem=%s -> base_id=%s", raw_stem, pdb_id)
    variant_token = _resolve_variant_token(config)
    variant_label = variant_token or "legacy"
    paths = canon_paths(pdb_id, output_root, variant=variant_token)
    log.info(
        "[prep.paths] protein_root=%s receptor=%s nolig=%s work=%s",
        paths["protein_root"],
        paths["receptor"],
        paths["nolig"],
        paths["work"],
    )
    receptor_target = paths["receptor"] / f"{pdb_id}_cleaned.pdb"
    log.info(
        "[receptor.clean.location] path=%s variant=%s variant_scoped=%s",
        receptor_target,
        variant_label,
        bool(variant_token),
    )

    log.info(
        "[proteinprep] entering clean_pdb pdb_file=%s output_root=%s",
        pdb_file,
        output_root,
    )
    input_path_for_metals = Path(pdb_file)
    if os.path.exists(str(input_path_for_metals)):
        _log_ion_count(log, 1, "elemfix", "BEFORE", input_path_for_metals)

    for d in ["protein_root", "raw", "work", "ligands_raw", "nolig", "receptor"]:
        paths[d].mkdir(parents=True, exist_ok=True)

    ion_audit: Any
    try:
        ion_audit = _IonAuditManager(pdb_id, paths["work"], variant_token)
    except NameError:

        class _NoOpIonAuditManager:
            def probe(self, *_args, **_kwargs):
                return None

        ion_audit = _NoOpIonAuditManager()
    _push_ion_audit_manager(ion_audit)

    try:
        ion_audit.probe("input", pdb_file)
        _emit_ion_breadcrumb("input", pdb_file)

        # (1) Working copy → raw/
        working_pdb = paths["raw"] / f"{pdb_id}_working.pdb"
        shutil.copyfile(str(pdb_file), working_pdb)
        _helium_postwrite_counter("copy_working", working_pdb)
        # (2) AltLoc filtering → raw/filtered.pdb
        filtered_pdb = paths["raw"] / f"{pdb_id}_filtered.pdb"
        filter_altlocs(working_pdb, filtered_pdb)

        # (2a) EARLY text-level element fix (YAML-driven), before any heavy tools
        try:
            _helium_postwrite_counter("elemfix_filtered", filtered_pdb)
            log.info("Early text-level element fix applied to %s", filtered_pdb)
        except Exception as e:
            log.warning(
                "Early text-level element fix skipped for %s: %s", filtered_pdb, e
            )

        # (3) Extract ligands now (controls live here), with YAML element repair per-file
        _ = extract_ligands_from_filtered(filtered_pdb, paths["ligands_raw"])
        _log_ions_probe(pdb_id, "extract_ligands", filtered_pdb)

        if os.environ.get("EARLY_CHAIN_PRUNE", "1").lower() not in {"0", "false", "no"}:
            try:
                with open(filtered_pdb, "r", encoding="utf-8", errors="ignore") as fh:
                    lines = fh.readlines()
                keep = _chains_to_keep(
                    lines, pocket_center=None
                )  # center can be wired later
                orig = {ln[21] for ln in lines if ln.startswith(("ATOM", "HETATM"))}
                if keep and keep != orig:
                    pruned = _prune_chains_conservative(lines, keep)
                    with open(filtered_pdb, "w", encoding="utf-8") as out:
                        out.writelines(pruned)
                    log.info(
                        "[chains] early-pruned chains keep=%s drop=%s",
                        "".join(sorted(keep)),
                        "".join(sorted(orig - keep)),
                    )
            except Exception as e:
                log.warning("[chains] early prune skipped: %s", e)
        # (3a) Optional early chain-prune (conservative, pocket-aware)
        source_for_strip = filtered_pdb
        if _cfg_bool("CHAIN_PRUNE", False):
            try:
                kept = select_chains_to_keep(filtered_pdb, paths["ligands_raw"], config)
                if kept:
                    pruned_filtered = paths["raw"] / f"{pdb_id}_filtered_pruned.pdb"
                    prune_to_chains(filtered_pdb, kept, pruned_filtered)
                    log.info(
                        "[chain_prune] using pruned source for step (4): %s",
                        pruned_filtered,
                    )
                    source_for_strip = pruned_filtered
                else:
                    log.info(
                        "[chain_prune] not applied (kept empty or conservative_abort); using unpruned file"
                    )
            except Exception as e:
                log.warning("[chain_prune] skipped due to exception: %s", e)

        # (4) Strip nonstandard from protein (policy aware) → work/stripped.pdb
        stripped_pdb = paths["work"] / f"{pdb_id}_stripped.pdb"
        ion_audit.probe("strip_nsr_before", source_for_strip)
        _emit_ion_breadcrumb("strip_nsr_before", source_for_strip)
        _log_ions_probe(pdb_id, "strip_nonstandard", source_for_strip, phase="before")
        if source_for_strip and os.path.exists(str(source_for_strip)):
            _log_ion_count(log, 5, "strip_nonstandard", "BEFORE", source_for_strip)
        removed_count, _out = strip_nonstandard_residues(
            source_for_strip,
            stripped_pdb,
            variant=variant_token,
        )
        _log_ions_probe(pdb_id, "strip_nonstandard", stripped_pdb, phase="after")
        ion_audit.probe("strip_nsr_after", stripped_pdb)
        _emit_ion_breadcrumb("strip_nsr_after", stripped_pdb)
        log.info("Removed %d nonstandard residue lines.", removed_count)
        if os.path.exists(str(stripped_pdb)):
            _log_ion_count(log, 5, "strip_nonstandard", "AFTER ", stripped_pdb)

        # (5) Element fix → open structural repair → element fix again (PDB only)
        elemfix_pdb = paths["work"] / f"{pdb_id}_elemfix.pdb"

        fix_pdb_elements(stripped_pdb, elemfix_pdb)
        ion_audit.probe("elemfix", elemfix_pdb)
        _log_ions_probe(pdb_id, "elemfix", elemfix_pdb)
        quick_element_histogram(elemfix_pdb)
        _helium_postwrite_counter("elemfix_before_modeller", elemfix_pdb)
        if os.path.exists(str(elemfix_pdb)):
            _log_ion_count(log, 1, "elemfix", "AFTER ", elemfix_pdb)
            _log_ion_count(log, 2, "openmm_repair", "BEFORE", elemfix_pdb)

        target_ph = get_target_ph_for_prep(elemfix_pdb)
        openmm_repaired_pdb = paths["work"] / f"{pdb_id}_openmm_repaired.pdb"
        openmm_ok = repair_with_pdbfixer(
            elemfix_pdb,
            openmm_repaired_pdb,
            target_ph=target_ph,
            cfg=config,
            logger=log,
            add_hydrogens=False,
        )
        use_modeller = _cfg_bool("USE_MODELLER", False) and not openmm_ok
        if use_modeller:
            loop_fixed_pdb = build_missing_loops(elemfix_pdb, paths["work"])
            repair_stage = "modeller"
        elif openmm_ok:
            loop_fixed_pdb = str(openmm_repaired_pdb)
            repair_stage = "openmm_repair"
        else:
            loop_fixed_pdb = str(elemfix_pdb)
            repair_stage = "copy"

        fix_pdb_elements(loop_fixed_pdb, loop_fixed_pdb)
        ion_audit.probe(repair_stage, loop_fixed_pdb)
        _emit_ion_breadcrumb(repair_stage, loop_fixed_pdb)
        _log_ions_probe(pdb_id, repair_stage, loop_fixed_pdb)
        quick_element_histogram(loop_fixed_pdb)
        _helium_postwrite_counter(f"elemfix_after_{repair_stage}", loop_fixed_pdb)
        if os.path.exists(str(loop_fixed_pdb)):
            _log_ion_count(log, 2, repair_stage, "AFTER ", loop_fixed_pdb)

        modeller_ok = repair_stage == "modeller" and os.path.isfile(loop_fixed_pdb)

        receptor_pdb = paths["receptor"] / f"{pdb_id}_cleaned.pdb"
        # (6) Optional external Phenix polish (non-fatal if missing)
        phenix_ok = False
        _use_phenix = str(
            config.get("use_phenix", config.get("USE_PHENIX", "false"))
        ).strip().lower() in ("1", "true", "yes")
        if _use_phenix:  # Water policy & radius
            _remove_waters = str(
                _cfg_env_or_default("REMOVE_WATERS", "true")
            ).strip().lower() in ("1", "true", "yes")
            _policy = (
                (_cfg_env_or_default("WATER_KEEP_POLICY", "none") or "none")
                .strip()
                .lower()
            )
            _keep_R = float(_cfg_env_or_default("KEEP_WATERS_WITHIN_A", "6.0") or 6.0)

            # Default: feed Phenix the loop-fixed input
            _phenix_in: Union[str, Path] = loop_fixed_pdb

            if _remove_waters and _policy != "none":
                # Policy active: derive reference points
                ref_pts: list[tuple[float, float, float]] = []
                # Prefer chosen center when available; here we’re early, so fall back to control centroids
                try:
                    ref_pts = compute_control_centroids(paths["ligands_raw"])
                except Exception:
                    ref_pts = []
                if ref_pts:
                    _phenix_in = paths["work"] / f"{pdb_id}_prefiltered_waters.pdb"
                    water_kept_count = filter_waters_near_points(
                        loop_fixed_pdb, _phenix_in, ref_pts, _keep_R
                    )
                    log.info(
                        "[waters] policy=%s kept=%d within %.1f Å of %d centers",
                        _policy,
                        water_kept_count,
                        _keep_R,
                        len(ref_pts),
                    )
                else:
                    log.info(
                        "[waters] policy=%s but no reference points found; skipping prefilter",
                        _policy,
                    )

            # Blanket removal only when policy is 'none'
            _phenix_remove = bool(_remove_waters and _policy == "none")
            if os.path.exists(str(_phenix_in)):
                _log_ion_count(log, 6, "phenix_clean", "BEFORE", _phenix_in)
            phenix_ok = run_phenix_pdbtools(
                input_pdb=_phenix_in,
                output_pdb=receptor_pdb,
                remove_waters=_phenix_remove,
            )

            # (6b) Dry-run sanity: count HOH within 8 Å of control-centroid center in final receptor
            try:
                ref_pts = compute_control_centroids(paths["ligands_raw"])
                center0 = None
                if ref_pts:
                    # quick average as an approximate center for the dry-run note
                    cx = sum(p[0] for p in ref_pts) / len(ref_pts)
                    cy = sum(p[1] for p in ref_pts) / len(ref_pts)
                    cz = sum(p[2] for p in ref_pts) / len(ref_pts)
                    center0 = (cx, cy, cz)
                if center0:
                    kept8 = count_waters_within(receptor_pdb, center0, 8.0)
                    log.info(
                        "[waters] dry-run kept_within_8A=%d center=(%.2f,%.2f,%.2f) file=%s",
                        kept8,
                        center0[0],
                        center0[1],
                        center0[2],
                        receptor_pdb,
                    )
            except Exception as _e:
                log.debug("[waters] dry-run check skipped: %s", _e)

        if not phenix_ok:
            shutil.copyfile(loop_fixed_pdb, receptor_pdb)
        _helium_postwrite_counter("phenix_or_copy_receptor", receptor_pdb)
        if os.path.exists(str(receptor_pdb)):
            _log_ion_count(log, 6, "phenix_clean", "AFTER ", receptor_pdb)

        reduce_deferred = True

        log.info(
            "[proteinprep] steps: Reduce=%s Phenix=%s MODELLER=%s",
            "deferred" if reduce_deferred else "applied",
            str(phenix_ok),
            str(modeller_ok),
        )
        if modeller_ok:
            sz = os.path.getsize(loop_fixed_pdb)
            log.info("[proteinprep] modeller_out=%s size=%d", loop_fixed_pdb, sz)
        elif openmm_ok:
            sz = os.path.getsize(loop_fixed_pdb)
            log.info("[proteinprep] openmm_repair_out=%s size=%d", loop_fixed_pdb, sz)
        else:
            log.info("[proteinprep] modeller_out=none (kept %s)", elemfix_pdb)

        if phenix_ok:
            sz = os.path.getsize(receptor_pdb)
            log.info("[proteinprep] phenix_applied_to=%s size=%d", receptor_pdb, sz)

        # (7) Hydrogen cleanup & chain validation
        debulked_pdb = paths["work"] / f"{pdb_id}_debulked.pdb"
        shutil.copyfile(receptor_pdb, debulked_pdb)
        clean_hydrogens(debulked_pdb, use_conect_if_reliable=True, conect_min_cov=0.6)

        chain_validated_pdb = paths["work"] / f"{pdb_id}_validated.pdb"
        if os.path.exists(str(debulked_pdb)):
            _log_ion_count(log, 3, "validate", "BEFORE", debulked_pdb)
        filter_invalid_chains(debulked_pdb, chain_validated_pdb)
        _helium_postwrite_counter("chain_validate", chain_validated_pdb)
        ion_audit.probe("altloc_validate", chain_validated_pdb)

        if os.path.exists(str(chain_validated_pdb)):
            _log_ion_count(log, 3, "validate", "AFTER ", chain_validated_pdb)

        try:
            _txt_before = Path(chain_validated_pdb).read_text(
                encoding="utf-8", errors="ignore"
            )
        except Exception:
            _txt_before = ""

        try:
            # Rewrite PDB element columns (77–78) using the unified rules (also for ATOM when rewrite_atoms=True)
            fix_element_columns_in_file(
                chain_validated_pdb, chain_validated_pdb, rewrite_atoms=True
            )

            # Secondary invariant: if any H-named atom still carries He, fix and summarize
            _txt_after = Path(chain_validated_pdb).read_text(
                encoding="utf-8", errors="ignore"
            )
            _fixed_text, _nname = assert_no_helium_in_hydrogen_names(_txt_after)
            if _nname > 0:
                Path(chain_validated_pdb).write_text(_fixed_text, encoding="utf-8")

            # Grep-friendly one-liner with He→H delta
            _before = scan_helium_counts(_txt_before)
            _after_text = (
                _fixed_text
                if _nname > 0
                else _txt_after
            )
            _after = scan_helium_counts(_after_text)
            _delta = max(0, _before - _after)
            log.info(
                f"[elem-fix] file={Path(chain_validated_pdb).name} stage=preflight He->H={_delta}"
            )
        except Exception as _e:
            log.warning(
                f"[elements] receptor preflight failed for {Path(chain_validated_pdb).name}: {_e}"
            )

        quick_element_histogram(chain_validated_pdb)

        # (8) Protonation (Reduce when safe; else Open Babel fallback)
        def _present_resnames(pdb_path: Path) -> set[str]:
            res = set()
            with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
                for ln in fh:
                    if ln.startswith(("ATOM  ", "HETATM")):
                        res.add(ln[17:20].strip().upper())
            return res

        present_resnames = _present_resnames(chain_validated_pdb)
        NUC_LIKE = set(_flatten_semicolons(RULES.get("nucleotide_like_resnames", [])))
        use_reduce = not any(r in NUC_LIKE for r in present_resnames)
        # ---- Prefer PDB2PQR (with PROPKA) at pipeline pH; fall back to Reduce ----
        # This uses the helper already defined earlier in this file.
        pdb_for_reduce, used_pdb2pqr, _pk_log = _protonate_with_pdb2pqr_if_available(
            str(
                chain_validated_pdb
            ),  # protonate the validated, ligand-free coordinates
            Path(
                paths["work"]
            ),  # write PROPKA/PDB2PQR artifacts into the work directory
            logging,
            variant=variant_token,
        )
        if used_pdb2pqr and pdb_for_reduce and Path(pdb_for_reduce).exists():
            _log_ions_probe(pdb_id, "pdb2pqr", pdb_for_reduce)
            quick_element_histogram(pdb_for_reduce)

        # If PDB2PQR succeeded, pdb_for_reduce now has hydrogens and titration states.
        # assign_protonation_states() will detect H presence and run Reduce WITHOUT -BUILD,
        # i.e., do flips/cleanup only. If PDB2PQR failed, pdb_for_reduce == chain_validated_pdb
        # and Reduce will run with -BUILD as needed.

        if not use_reduce:
            log.info(
                "[protonation] Skipping Reduce due to detected nucleotides; using OpenBabel path."
            )

        reduced_pdb = paths["work"] / f"{pdb_id}_reduced.pdb"
        reduce_input_path = pdb_for_reduce if pdb_for_reduce else chain_validated_pdb
        if reduce_input_path and os.path.exists(str(reduce_input_path)):
            _log_ion_count(log, 4, "reduce", "BEFORE", reduce_input_path)
        assign_protonation_states(
            pdb_for_reduce,
            reduced_pdb,
            reduce_exe=REDUCE_EXE
            if use_reduce
            else None,  # skip Reduce for nucleotide cofactors
        )
        rescue_retained_hets(
            source_pdb=chain_validated_pdb,
            target_pdb=reduced_pdb,
            stage="post_reduce",
            cfg=config,
            logger=log,
        )

        _helium_postwrite_counter("reduce_or_fallback", reduced_pdb)
        _log_ions_probe(pdb_id, "reduce", reduced_pdb)
        ion_audit.probe("reduce", reduced_pdb)
        print(
            f"[proteinprep] Reduce/alt_protonation wrote={Path(reduced_pdb).is_file()} -> {reduced_pdb}"
        )
        if os.path.exists(str(reduced_pdb)):
            _log_ion_count(log, 4, "reduce", "AFTER ", reduced_pdb)

        # (9) Final element fix and sanity on the protonated file
        fix_pdb_elements(reduced_pdb)
        _helium_postwrite_counter("elemfix_after_reduce", reduced_pdb)

        quick_element_histogram(reduced_pdb)

        assert_no_metal_in_peptidic(reduced_pdb)

        # (10) Move to receptor. Element normalization is already done on the
        # reduced PDB and is repeated once at the Meeko boundary.
        # [ions] receptor_write audit
        variant_env = (
            variant_token or (os.environ.get("APO_HOLO_VARIANT") or "").strip().upper()
        )
        variant_label = variant_env if variant_env else "legacy"
        before_counts, before_detail = _collect_monoatomic_records(reduced_pdb)
        log.info(
            "[ions.policy] stage=receptor_write variant=%s source=%s target=%s reason=copy_reduced_to_cleaned",
            variant_label,
            reduced_pdb,
            receptor_pdb,
        )
        log.info(
            "[ions.counts.before] stage=receptor_write file=%s metals=%s",
            reduced_pdb,
            _format_ion_hist(before_counts),
        )

        shutil.copyfile(reduced_pdb, receptor_pdb)
        _helium_postwrite_counter("promote_receptor_copy", receptor_pdb)
        _log_ions_probe(pdb_id, "receptor_write", receptor_pdb)

        after_counts, after_detail = _collect_monoatomic_records(receptor_pdb)
        log.info(
            "[ions.counts.after] stage=receptor_write file=%s metals=%s",
            receptor_pdb,
            _format_ion_hist(after_counts),
        )
        diff_list = _diff_detail_records(before_detail, after_detail)
        kept_total = sum(after_counts.values())
        stripped_total = max(0, sum(before_counts.values()) - kept_total)
        log.info(
            "[ions.receptor.copy] action=write_cleaned kept=%d stripped=%d changed=%d",
            kept_total,
            stripped_total,
            len(diff_list),
        )
        if diff_list:
            log.info("[ions.diff.reduced→cleaned] lost=%s", ",".join(diff_list))

        _helium_postwrite_counter("final_receptor", receptor_pdb, fix_elements=False)
        ion_audit.probe("final_cleaned", receptor_pdb)

        receptor_pdb = _apply_binding_site_network_stage(
            receptor_pdb=receptor_pdb,
            filtered_pdb=filtered_pdb,
            ligand_dir=paths["ligands_raw"],
            work_dir=paths["work"],
            pdb_id=pdb_id,
            target_ph=target_ph,
            log=log,
        )
        _helium_postwrite_counter("binding_site_network_stage", receptor_pdb)
        ion_audit.probe("binding_site_network_stage", receptor_pdb)

        receptor_pdb = _apply_duplicate_metal_site_policy_stage(
            receptor_pdb=receptor_pdb,
            work_dir=paths["work"],
            pdb_id=pdb_id,
            log=log,
        )
        _helium_postwrite_counter("duplicate_metal_site_policy", receptor_pdb)
        ion_audit.probe("duplicate_metal_site_policy", receptor_pdb)

        receptor_pdb = _apply_binding_site_geometry_policy_stage(
            receptor_pdb=receptor_pdb,
            reference_pdb=filtered_pdb,
            ligand_dir=paths["ligands_raw"],
            work_dir=paths["work"],
            pdb_id=pdb_id,
            target_ph=target_ph,
            log=log,
        )
        _helium_postwrite_counter("binding_site_geometry_policy", receptor_pdb)
        ion_audit.probe("binding_site_geometry_policy", receptor_pdb)

        quick_element_histogram(receptor_pdb)
        assert file_contains_hydrogens(
            receptor_pdb
        ), f"[FATAL] Cleaned file lost hydrogens: {receptor_pdb}"

        log.info("Cleaned receptor: %s", receptor_pdb)
        if os.path.exists(str(receptor_pdb)):
            prepare_before_count = _log_ion_count(
                log, 7, "prepare_receptor", "BEFORE", receptor_pdb
            )
            log.info(
                "(7) prepare_receptor AFTER : metals=%d file=%s (pdbqt not parsed)",
                prepare_before_count,
                receptor_pdb,
            )
        try:
            path_cfg = dict(config)
            path_cfg["OUTPUT_DIR"] = str(output_root)
            run_id = (
                os.environ.get("ATLAS_RUN_ID") or str(path_cfg.get("RUN_ID", ""))
            ).strip()
            if run_id:
                path_cfg["RUN_ID"] = run_id
            router_paths = make_paths(path_cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
            receptor_pdbqt_path = router_paths.receptor_pdbqt(
                variant_token, ph_token=None
            )
            receptor_dir_path = router_paths.receptor_dir(variant_token)
            ph_dir = receptor_dir_path / "ph_ensemble"
            variant_log = (variant_token or "NONE").upper()
            log.info(
                "[receptor.path.final] variant=%s receptor_pdbqt=%s",
                variant_log,
                receptor_pdbqt_path,
            )
            ph_enabled = 0
            try:
                if ph_dir.exists():
                    next(ph_dir.iterdir())
                    ph_enabled = 1
            except StopIteration:
                ph_enabled = 0
            except Exception:
                ph_enabled = 1 if ph_dir.exists() else 0
            log.info(
                "[receptor.path.ensemble] enabled=%d dir=%s",
                ph_enabled,
                ph_dir,
            )
        except Exception as exc:
            log.warning(
                "[receptor.path.final] variant=%s action=skip reason=%s",
                (variant_token or "NONE").upper(),
                exc,
            )
        summary_count = prepare_before_count
        if summary_count is None and os.path.exists(str(receptor_pdb)):
            summary_count = _count_ion_like_records(receptor_pdb)
        if summary_count is not None:
            summary_variant = variant_token or "legacy"
            log.info(
                f"[ions.clean.counts] pdb={pdb_id} variant={summary_variant} file={str(receptor_pdb)} present_pdb=metals:{summary_count}"
            )
        print(
            f"[proteinprep] cleaned receptor exists={Path(receptor_pdb).is_file()} -> {receptor_pdb}"
        )
        return str(receptor_pdb)
    finally:
        _pop_ion_audit_manager(ion_audit)

def clean_pdb(
    pdb_file: Union[str, Path],
    output_root: Union[str, Path],
    logger: Optional[logging.Logger] = None,
) -> Optional[str]:
    return _run_clean_pdb_original(
        pdb_file=pdb_file,
        output_root=output_root,
        logger=logger,
    )
