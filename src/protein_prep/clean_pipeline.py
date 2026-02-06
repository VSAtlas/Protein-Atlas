"""Protein cleaning pipeline orchestration helpers."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Optional, Union


def _hydrate_legacy_globals() -> None:
    import automate_protein_prep as legacy

    g = globals()
    for name, value in legacy.__dict__.items():
        g.setdefault(name, value)


def _run_clean_pdb_original(
    pdb_file: Union[str, Path],
    output_root: Union[str, Path],
    logger: Optional[logging.Logger] = None,
) -> Optional[str]:
    """Run the full cleaning pipeline and return path to final cleaned PDB (receptor)."""
    ion_elements = {
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
        elemfix_before_count = 0
        with open(
            input_path_for_metals, "r", encoding="utf-8", errors="ignore"
        ) as _handle:
            for _line in _handle:
                if not _line.startswith(("ATOM  ", "HETATM")):
                    continue
                _element = _line[76:78].strip().upper()
                if not _element:
                    _name_field = _line[12:16].strip()
                    _guess = []
                    for _ch in _name_field:
                        if _ch.isalpha():
                            _guess.append(_ch)
                        else:
                            break
                    _guess_text = "".join(_guess).upper()
                    if len(_guess_text) >= 2 and _guess_text[:2] in ion_elements:
                        _element = _guess_text[:2]
                    elif _guess_text[:1] in ion_elements:
                        _element = _guess_text[:1]
                if _element in ion_elements:
                    elemfix_before_count += 1
        # (1) elemfix BEFORE: metals=<elemfix_before_count>
        log.info(
            f"(1) elemfix BEFORE: metals={elemfix_before_count} file={str(input_path_for_metals)}"
        )

    for d in ["protein_root", "raw", "work", "ligands_raw", "nolig", "receptor"]:
        paths[d].mkdir(parents=True, exist_ok=True)

    try:
        ion_audit = _IonAuditManager(pdb_id, paths["work"], variant_token)
    except NameError:

        class _NoOpIonAuditManager:
            enabled = False

            def probe(self, *args, **kwargs):
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
            fix_element_columns_in_file(filtered_pdb, filtered_pdb, rewrite_atoms=True)
            _helium_postwrite_counter("elemfix_filtered", filtered_pdb)
            log.info("Early text-level element fix applied to %s", filtered_pdb)
        except Exception as e:
            log.warning(
                "Early text-level element fix skipped for %s: %s", filtered_pdb, e
            )

        # (3) Extract ligands now (controls live here), with YAML element repair per-file
        _ = extract_ligands_from_filtered(filtered_pdb, paths["ligands_raw"])
        _log_ions_probe(pdb_id, "extract_ligands", filtered_pdb)

        chain_pruned = False
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
                    chain_pruned = True
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
                    chain_pruned = True
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
            strip_nonstandard_before_count = 0
            with open(
                source_for_strip, "r", encoding="utf-8", errors="ignore"
            ) as _handle:
                for _line in _handle:
                    if not _line.startswith(("ATOM  ", "HETATM")):
                        continue
                    _element = _line[76:78].strip().upper()
                    if not _element:
                        _name_field = _line[12:16].strip()
                        _guess = []
                        for _ch in _name_field:
                            if _ch.isalpha():
                                _guess.append(_ch)
                            else:
                                break
                        _guess_text = "".join(_guess).upper()
                        if len(_guess_text) >= 2 and _guess_text[:2] in ion_elements:
                            _element = _guess_text[:2]
                        elif _guess_text[:1] in ion_elements:
                            _element = _guess_text[:1]
                    if _element in ion_elements:
                        strip_nonstandard_before_count += 1
            # (5) strip_nonstandard BEFORE: metals=<strip_nonstandard_before_count>
            log.info(
                f"(5) strip_nonstandard BEFORE: metals={strip_nonstandard_before_count} file={str(source_for_strip)}"
            )
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
            strip_nonstandard_after_count = 0
            with open(stripped_pdb, "r", encoding="utf-8", errors="ignore") as _handle:
                for _line in _handle:
                    if not _line.startswith(("ATOM  ", "HETATM")):
                        continue
                    _element = _line[76:78].strip().upper()
                    if not _element:
                        _name_field = _line[12:16].strip()
                        _guess = []
                        for _ch in _name_field:
                            if _ch.isalpha():
                                _guess.append(_ch)
                            else:
                                break
                        _guess_text = "".join(_guess).upper()
                        if len(_guess_text) >= 2 and _guess_text[:2] in ion_elements:
                            _element = _guess_text[:2]
                        elif _guess_text[:1] in ion_elements:
                            _element = _guess_text[:1]
                    if _element in ion_elements:
                        strip_nonstandard_after_count += 1
            # (5) strip_nonstandard AFTER:  metals=<strip_nonstandard_after_count>
            log.info(
                f"(5) strip_nonstandard AFTER:  metals={strip_nonstandard_after_count} file={str(stripped_pdb)}"
            )

        # (5) Element fix → MODELLER → element fix again (PDB only)
        elemfix_pdb = paths["work"] / f"{pdb_id}_elemfix.pdb"

        fix_pdb_elements(stripped_pdb, elemfix_pdb)
        ion_audit.probe("elemfix", elemfix_pdb)
        _log_ions_probe(pdb_id, "elemfix", elemfix_pdb)
        quick_element_histogram(elemfix_pdb)
        _helium_postwrite_counter("elemfix_before_modeller", elemfix_pdb)
        if os.path.exists(str(elemfix_pdb)):
            elemfix_after_count = 0
            with open(elemfix_pdb, "r", encoding="utf-8", errors="ignore") as _handle:
                for _line in _handle:
                    if not _line.startswith(("ATOM  ", "HETATM")):
                        continue
                    _element = _line[76:78].strip().upper()
                    if not _element:
                        _name_field = _line[12:16].strip()
                        _guess = []
                        for _ch in _name_field:
                            if _ch.isalpha():
                                _guess.append(_ch)
                            else:
                                break
                        _guess_text = "".join(_guess).upper()
                        if len(_guess_text) >= 2 and _guess_text[:2] in ion_elements:
                            _element = _guess_text[:2]
                        elif _guess_text[:1] in ion_elements:
                            _element = _guess_text[:1]
                    if _element in ion_elements:
                        elemfix_after_count += 1
            # (1) elemfix AFTER:  metals=<elemfix_after_count>
            log.info(
                f"(1) elemfix AFTER:  metals={elemfix_after_count} file={str(elemfix_pdb)}"
            )
            modeller_before_count = elemfix_after_count
            # (2) modeller_fill BEFORE: metals=<modeller_before_count>
            log.info(
                f"(2) modeller_fill BEFORE: metals={modeller_before_count} file={str(elemfix_pdb)}"
            )
        loop_fixed_pdb = build_missing_loops(elemfix_pdb, paths["work"])
        fix_pdb_elements(loop_fixed_pdb, loop_fixed_pdb)
        ion_audit.probe("modeller", loop_fixed_pdb)
        _emit_ion_breadcrumb("modeller", loop_fixed_pdb)
        _log_ions_probe(pdb_id, "modeller", loop_fixed_pdb)
        quick_element_histogram(loop_fixed_pdb)
        _helium_postwrite_counter("elemfix_after_modeller", loop_fixed_pdb)
        if os.path.exists(str(loop_fixed_pdb)):
            modeller_after_count = 0
            with open(
                loop_fixed_pdb, "r", encoding="utf-8", errors="ignore"
            ) as _handle:
                for _line in _handle:
                    if not _line.startswith(("ATOM  ", "HETATM")):
                        continue
                    _element = _line[76:78].strip().upper()
                    if not _element:
                        _name_field = _line[12:16].strip()
                        _guess = []
                        for _ch in _name_field:
                            if _ch.isalpha():
                                _guess.append(_ch)
                            else:
                                break
                        _guess_text = "".join(_guess).upper()
                        if len(_guess_text) >= 2 and _guess_text[:2] in ion_elements:
                            _element = _guess_text[:2]
                        elif _guess_text[:1] in ion_elements:
                            _element = _guess_text[:1]
                    if _element in ion_elements:
                        modeller_after_count += 1
            # (2) modeller_fill AFTER:  metals=<modeller_after_count>
            log.info(
                f"(2) modeller_fill AFTER:  metals={modeller_after_count} file={str(loop_fixed_pdb)}"
            )

        # MODELLER (detect whether a new file was actually produced)
        modeller_ok = os.path.basename(
            loop_fixed_pdb
        ) == "modeller_filled.pdb" and os.path.isfile(loop_fixed_pdb)

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
            _phenix_in = loop_fixed_pdb

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
                    kept = filter_waters_near_points(
                        loop_fixed_pdb, _phenix_in, ref_pts, _keep_R
                    )
                    log.info(
                        "[waters] policy=%s kept=%d within %.1f Å of %d centers",
                        _policy,
                        kept,
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
                phenix_before_count = 0
                with open(
                    _phenix_in, "r", encoding="utf-8", errors="ignore"
                ) as _handle:
                    for _line in _handle:
                        if not _line.startswith(("ATOM  ", "HETATM")):
                            continue
                        _element = _line[76:78].strip().upper()
                        if not _element:
                            _name_field = _line[12:16].strip()
                            _guess = []
                            for _ch in _name_field:
                                if _ch.isalpha():
                                    _guess.append(_ch)
                                else:
                                    break
                            _guess_text = "".join(_guess).upper()
                            if (
                                len(_guess_text) >= 2
                                and _guess_text[:2] in ion_elements
                            ):
                                _element = _guess_text[:2]
                            elif _guess_text[:1] in ion_elements:
                                _element = _guess_text[:1]
                        if _element in ion_elements:
                            phenix_before_count += 1
                # (6) phenix_clean BEFORE: metals=<phenix_before_count>
                log.info(
                    f"(6) phenix_clean BEFORE: metals={phenix_before_count} file={str(_phenix_in)}"
                )
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
            phenix_after_count = 0
            with open(receptor_pdb, "r", encoding="utf-8", errors="ignore") as _handle:
                for _line in _handle:
                    if not _line.startswith(("ATOM  ", "HETATM")):
                        continue
                    _element = _line[76:78].strip().upper()
                    if not _element:
                        _name_field = _line[12:16].strip()
                        _guess = []
                        for _ch in _name_field:
                            if _ch.isalpha():
                                _guess.append(_ch)
                            else:
                                break
                        _guess_text = "".join(_guess).upper()
                        if len(_guess_text) >= 2 and _guess_text[:2] in ion_elements:
                            _element = _guess_text[:2]
                        elif _guess_text[:1] in ion_elements:
                            _element = _guess_text[:1]
                    if _element in ion_elements:
                        phenix_after_count += 1
            # (6) phenix_clean AFTER:  metals=<phenix_after_count>
            log.info(
                f"(6) phenix_clean AFTER:  metals={phenix_after_count} file={str(receptor_pdb)}"
            )

        # choose the file to pass downstream
        pdb_for_reduce = loop_fixed_pdb if modeller_ok else elemfix_pdb

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
            validate_before_count = 0
            with open(debulked_pdb, "r", encoding="utf-8", errors="ignore") as _handle:
                for _line in _handle:
                    if not _line.startswith(("ATOM  ", "HETATM")):
                        continue
                    _element = _line[76:78].strip().upper()
                    if not _element:
                        _name_field = _line[12:16].strip()
                        _guess = []
                        for _ch in _name_field:
                            if _ch.isalpha():
                                _guess.append(_ch)
                            else:
                                break
                        _guess_text = "".join(_guess).upper()
                        if len(_guess_text) >= 2 and _guess_text[:2] in ion_elements:
                            _element = _guess_text[:2]
                        elif _guess_text[:1] in ion_elements:
                            _element = _guess_text[:1]
                    if _element in ion_elements:
                        validate_before_count += 1
            # (3) validate BEFORE: metals=<validate_before_count>
            log.info(
                f"(3) validate BEFORE: metals={validate_before_count} file={str(debulked_pdb)}"
            )
        filter_invalid_chains(debulked_pdb, chain_validated_pdb)
        _helium_postwrite_counter("chain_validate", chain_validated_pdb)
        ion_audit.probe("altloc_validate", chain_validated_pdb)

        if os.path.exists(str(chain_validated_pdb)):
            validate_after_count = 0
            with open(
                chain_validated_pdb, "r", encoding="utf-8", errors="ignore"
            ) as _handle:
                for _line in _handle:
                    if not _line.startswith(("ATOM  ", "HETATM")):
                        continue
                    _element = _line[76:78].strip().upper()
                    if not _element:
                        _name_field = _line[12:16].strip()
                        _guess = []
                        for _ch in _name_field:
                            if _ch.isalpha():
                                _guess.append(_ch)
                            else:
                                break
                        _guess_text = "".join(_guess).upper()
                        if len(_guess_text) >= 2 and _guess_text[:2] in ion_elements:
                            _element = _guess_text[:2]
                        elif _guess_text[:1] in ion_elements:
                            _element = _guess_text[:1]
                    if _element in ion_elements:
                        validate_after_count += 1
            # (3) validate AFTER:  metals=<validate_after_count>
            log.info(
                f"(3) validate AFTER:  metals={validate_after_count} file={str(chain_validated_pdb)}"
            )

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
            _after = scan_helium_counts(
                Path(chain_validated_pdb).read_text(encoding="utf-8", errors="ignore")
            )
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
        pdb_for_reduce, used_pdb2pqr, pk_log = _protonate_with_pdb2pqr_if_available(
            str(
                chain_validated_pdb
            ),  # protonate the validated, ligand-free coordinates
            str(
                paths["work"]
            ),  # write PROPKA/PDB2PQR artifacts into the work directory
            logging,
            variant=variant_token,
        )
        if used_pdb2pqr and pdb_for_reduce and Path(pdb_for_reduce).exists():
            _log_ions_probe(pdb_id, "pdb2pqr", pdb_for_reduce)
            quick_element_histogram(pdb_for_reduce)

        validated_source = (
            Path(pdb_for_reduce) if pdb_for_reduce else Path(chain_validated_pdb)
        )

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
            reduce_before_count = 0
            with open(
                reduce_input_path, "r", encoding="utf-8", errors="ignore"
            ) as _handle:
                for _line in _handle:
                    if not _line.startswith(("ATOM  ", "HETATM")):
                        continue
                    _element = _line[76:78].strip().upper()
                    if not _element:
                        _name_field = _line[12:16].strip()
                        _guess = []
                        for _ch in _name_field:
                            if _ch.isalpha():
                                _guess.append(_ch)
                            else:
                                break
                        _guess_text = "".join(_guess).upper()
                        if len(_guess_text) >= 2 and _guess_text[:2] in ion_elements:
                            _element = _guess_text[:2]
                        elif _guess_text[:1] in ion_elements:
                            _element = _guess_text[:1]
                    if _element in ion_elements:
                        reduce_before_count += 1
            # (4) reduce BEFORE: metals=<reduce_before_count>
            log.info(
                f"(4) reduce BEFORE: metals={reduce_before_count} file={str(reduce_input_path)}"
            )
        assign_protonation_states(
            pdb_for_reduce,
            reduced_pdb,
            reduce_exe=REDUCE_EXE
            if use_reduce
            else None,  # skip Reduce for nucleotide cofactors
        )

        _helium_postwrite_counter("reduce_or_fallback", reduced_pdb)
        _log_ions_probe(pdb_id, "reduce", reduced_pdb)
        ion_audit.probe("reduce", reduced_pdb)
        print(
            f"[proteinprep] Reduce/alt_protonation wrote={Path(reduced_pdb).is_file()} -> {reduced_pdb}"
        )
        if os.path.exists(str(reduced_pdb)):
            reduce_after_count = 0
            with open(reduced_pdb, "r", encoding="utf-8", errors="ignore") as _handle:
                for _line in _handle:
                    if not _line.startswith(("ATOM  ", "HETATM")):
                        continue
                    _element = _line[76:78].strip().upper()
                    if not _element:
                        _name_field = _line[12:16].strip()
                        _guess = []
                        for _ch in _name_field:
                            if _ch.isalpha():
                                _guess.append(_ch)
                            else:
                                break
                        _guess_text = "".join(_guess).upper()
                        if len(_guess_text) >= 2 and _guess_text[:2] in ion_elements:
                            _element = _guess_text[:2]
                        elif _guess_text[:1] in ion_elements:
                            _element = _guess_text[:1]
                    if _element in ion_elements:
                        reduce_after_count += 1
            # (4) reduce AFTER:  metals=<reduce_after_count>
            log.info(
                f"(4) reduce AFTER:  metals={reduce_after_count} file={str(reduced_pdb)}"
            )

        # (9) Final element fix and sanity on the protonated file
        fix_pdb_elements(reduced_pdb)
        _helium_postwrite_counter("elemfix_after_reduce", reduced_pdb)

        quick_element_histogram(reduced_pdb)

        assert_no_metal_in_peptidic(reduced_pdb)

        # (10) Move to receptor and re-fix (post-step edits)
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

        fix_pdb_elements(receptor_pdb)
        _helium_postwrite_counter("elemfix_final_receptor", receptor_pdb)
        ion_audit.probe("final_cleaned", receptor_pdb)

        quick_element_histogram(receptor_pdb)
        assert file_contains_hydrogens(
            receptor_pdb
        ), f"[FATAL] Cleaned file lost hydrogens: {receptor_pdb}"

        log.info("Cleaned receptor: %s", receptor_pdb)
        if os.path.exists(str(receptor_pdb)):
            prepare_before_count = 0
            with open(receptor_pdb, "r", encoding="utf-8", errors="ignore") as _handle:
                for _line in _handle:
                    if not _line.startswith(("ATOM  ", "HETATM")):
                        continue
                    _element = _line[76:78].strip().upper()
                    if not _element:
                        _name_field = _line[12:16].strip()
                        _guess = []
                        for _ch in _name_field:
                            if _ch.isalpha():
                                _guess.append(_ch)
                            else:
                                break
                        _guess_text = "".join(_guess).upper()
                        if len(_guess_text) >= 2 and _guess_text[:2] in ion_elements:
                            _element = _guess_text[:2]
                        elif _guess_text[:1] in ion_elements:
                            _element = _guess_text[:1]
                    if _element in ion_elements:
                        prepare_before_count += 1
            # (7) prepare_receptor4 BEFORE: metals=<prepare_before_count>
            log.info(
                f"(7) prepare_receptor4 BEFORE: metals={prepare_before_count} file={str(receptor_pdb)}"
            )
            # (7) prepare_receptor4 AFTER:  metals=<prepare_before_count>
            log.info(
                f"(7) prepare_receptor4 AFTER:  metals={prepare_before_count} file={str(receptor_pdb)} (pdbqt not parsed)"
            )
        try:
            router_paths = make_paths(config, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
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
            summary_count = 0
            with open(receptor_pdb, "r", encoding="utf-8", errors="ignore") as _handle:
                for _line in _handle:
                    if not _line.startswith(("ATOM  ", "HETATM")):
                        continue
                    _element = _line[76:78].strip().upper()
                    if not _element:
                        _name_field = _line[12:16].strip()
                        _guess = []
                        for _ch in _name_field:
                            if _ch.isalpha():
                                _guess.append(_ch)
                            else:
                                break
                        _guess_text = "".join(_guess).upper()
                        if len(_guess_text) >= 2 and _guess_text[:2] in ion_elements:
                            _element = _guess_text[:2]
                        elif _guess_text[:1] in ion_elements:
                            _element = _guess_text[:1]
                    if _element in ion_elements:
                        summary_count += 1
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



def init_workspace_and_paths(ctx: dict[str, Any]) -> dict[str, Any]:
    return ctx


def step_altloc_filter(ctx: dict[str, Any]) -> dict[str, Any]:
    return ctx


def step_strip_nonstandard(ctx: dict[str, Any]) -> dict[str, Any]:
    return ctx


def step_modeller_loopfill(ctx: dict[str, Any]) -> dict[str, Any]:
    return ctx


def step_phenix_clean(ctx: dict[str, Any]) -> dict[str, Any]:
    return ctx


def step_reduce_protonation(ctx: dict[str, Any]) -> dict[str, Any]:
    return ctx


def step_prepare_receptor(ctx: dict[str, Any]) -> dict[str, Any]:
    _hydrate_legacy_globals()
    ctx["result"] = _run_clean_pdb_original(
        pdb_file=ctx["pdb_file"],
        output_root=ctx["output_root"],
        logger=ctx.get("logger"),
    )
    return ctx


def finalize_outputs_and_logs(ctx: dict[str, Any]) -> Optional[str]:
    return ctx.get("result")


def clean_pdb(
    pdb_file: Union[str, Path],
    output_root: Union[str, Path],
    logger: Optional[logging.Logger] = None,
) -> Optional[str]:
    ctx: dict[str, Any] = {
        "pdb_file": pdb_file,
        "output_root": output_root,
        "logger": logger,
    }
    for step in (
        init_workspace_and_paths,
        step_altloc_filter,
        step_strip_nonstandard,
        step_modeller_loopfill,
        step_phenix_clean,
        step_reduce_protonation,
        step_prepare_receptor,
    ):
        ctx = step(ctx)
    return finalize_outputs_and_logs(ctx)
