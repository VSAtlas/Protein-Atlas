from __future__ import annotations

import json
import logging
import os
import shutil
from distutils.util import strtobool
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Tuple

from path_router import Paths, docked_dir

from docking_controls import (
    receptor_sanity_check,
    _ph_values_from_context,
    _ph_ligand_mode,
    _resolve_ph_scope,
    select_center_via_control_redock,
    detect_pocket,
    _summarize_ions_file,
    build_control_lookup,
    extract_ligands_to_nolig,
)

from docking_ligands import _lib_roots_for_pdb
from prep_ligands.prep_ligands_microstates import enumerate_ligands_for_docking

# NOTE: `automate_protein_prep` is imported inside prepare_receptor itself.
# There is no need to import it at the top level.


def norm(p: str | Path) -> str:
    """
    Normalize path to a clean, forward-slash string for logs & keys.

    This is intentionally duplicated from docking.py to avoid a circular
    import (docking imports prepare_receptor from this module).
    """
    return os.path.abspath(str(p)).replace("\\", "/")


def prepare_receptor(
    cfg: Dict,
    paths: Paths,
    logger: logging.Logger,
    center: Optional[Tuple[float, float, float]] = None,
    box_size: Optional[Tuple[float, float, float]] = None,
) -> Tuple[Optional[str], Optional[str]]:
    import automate_protein_prep
    from distutils.util import strtobool

    prepare_receptor.last_provenance = "unknown"
    force_reprocess = bool(strtobool(str(cfg.get("FORCE_REPROCESS", False))))
    log = logging.getLogger("ph_ensemble")
    var = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    variant = var if var else None
    ph_token = None
    cleaned_pdb_path = paths.receptor_cleaned_pdb(variant)
    receptor_pdbqt_path = paths.receptor_pdbqt(variant, ph_token)
    logger.info(
        f"FORCE_REPROCESS={force_reprocess} | "
        f"cleaned_exists={cleaned_pdb_path.exists()} "
        f"receptor_exists={receptor_pdbqt_path.exists()}"
    )

    def _build_ph_ensemble(cleaned_path: str) -> Optional[str]:
        if not bool(cfg.get("PH_ENSEMBLE")):
            return None
        cleaned_path = str(cleaned_path)
        log.info("[ph_ensemble.anchor] cleaned_receptor_pdb=%s", cleaned_path)
        log.info("[ph_ensemble.begin] pdb_id=%s path=%s", paths.pdb_id, cleaned_path)

        def _collect_dock_targets(manifest_path: str) -> Optional[list[tuple[str, str]]]:
            manifest_file = Path(manifest_path)
            try:
                payload = json.loads(manifest_file.read_text())
            except Exception as exc:
                log.error("[ph_ensemble.manifest.read.error] path=%s err=%s", manifest_path, exc)
                return None

            members = payload.get("members") or []
            # Diagnostics: enumerate keys and canonical counts.
            try:
                key_universe = sorted({k for m in members for k in (m.keys() if isinstance(m, dict) else [])})
            except Exception:
                key_universe = []
            log.info(
                "[ph_ensemble.manifest.stats] members=%d canonical=%d keys=%s",
                len(members),
                sum(1 for m in members if isinstance(m, dict) and bool(m.get("canonical", False))),
                ",".join(key_universe),
            )

            canonical = [m for m in members if isinstance(m, dict) and bool(m.get("canonical", False))]
            if canonical:
                members = canonical

            prefix = f"{paths.pdb_id}_"
            targets: list[tuple[str, str]] = []
            for entry in members:
                if not isinstance(entry, dict):
                    continue
                receptor_path = (
                    entry.get("pdbqt")
                    or entry.get("receptor_pdbqt")
                    or entry.get("output_pdbqt")
                    or entry.get("path")
                    or entry.get("receptor")
                )
                if not receptor_path:
                    log.warning("[ph_ensemble.manifest.entry.missing_pdbqt] keys=%s", list(entry.keys()))
                    continue

                ph_label = entry.get("label") or entry.get("ph_label")
                if not ph_label:
                    stem = Path(receptor_path).stem
                    ph_label = stem[len(prefix):] if stem.startswith(prefix) else stem

                targets.append((str(ph_label), str(receptor_path)))

            if not targets:
                log.error("[ph_ensemble.manifest.no_targets] path=%s members=%d", manifest_path, len(members))
            return targets

        def _bridge_manifest_targets(manifest_path: str) -> None:
            targets = _collect_dock_targets(manifest_path)
            if targets is None:
                return
            cfg.setdefault("_PH_ENSEMBLE_CANONICAL", {})[paths.pdb_id] = targets
            # --- mapping debug (anchor: [ph_ensemble.map]) ---
            log.info(
                "[ph_ensemble.map] pdb_id=%s canonical=%s",
                paths.pdb_id,
                ";".join(f"{lbl}:{Path(p).name}" for lbl, p in targets) if targets else ""
            )
            variant_label = variant or "legacy"
            log.info(
                "[ph_ensemble.dock.begin] pdb_id=%s variant=%s n=%d labels=%s",
                paths.pdb_id,
                variant_label,
                len(targets),
                ",".join(lbl for lbl, _ in targets) or "",
            )

            if not targets:
                log.info(
                    "[ph_ensemble.dock.done] pdb_id=%s variant=%s n=0 targets=",
                    paths.pdb_id,
                    variant_label,
                )
                return
            legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
            for ph_label, _ in targets:
                ph_root = docked_dir(paths.pdb_id, variant=variant, ph_tag=ph_label, legacy=legacy_mode)
                log.info(
                    "[ph_ensemble.dock.root] pdb_id=%s variant=%s ph=%s dock_root=%s",
                    paths.pdb_id,
                    variant_label,
                    ph_label,
                    str(ph_root),
                )
            summary = ";".join(f"{ph}:{rec}" for ph, rec in targets)
            log.info(
                "[ph_ensemble.dock.done] pdb_id=%s variant=%s n=%d targets=%s",
                paths.pdb_id,
                variant_label,
                len(targets),
                summary,
            )

        try:
            from path_router.context_ph import select_ph_from_pdb
            ctx_result = select_ph_from_pdb(cleaned_path)
        except Exception as exc:
            log.warning("[ph_ensemble.ctx.error] %s", exc)
            ctx_result = {"target_pH": 7.0, "ensemble": None}
        target_pH = float(ctx_result.get("target_pH", 7.0) or 7.0)
        ensemble_from_context = ctx_result.get("ensemble")
        log.info("[ph_ensemble.ctx] target_pH=%.2f raw_ensemble=%s", target_pH, repr(ensemble_from_context))
        raw_values = list(ensemble_from_context or [target_pH])
        ph_values: list[float] = []
        for value in raw_values:
            try:
                ph = float(value)
            except Exception:
                continue
            ph = round(ph, 1)
            if ph < 3.0:
                ph = 3.0
            if ph > 10.5:
                ph = 10.5
            ph_values.append(ph)
        if not ph_values:
            fallback = round(target_pH, 1)
            if fallback < 3.0:
                fallback = 3.0
            if fallback > 10.5:
                fallback = 10.5
            ph_values = [fallback]
        ph_values = sorted({round(p, 1) for p in ph_values})
        log.info("[ph_ensemble.list] canonical=%s", ",".join(f"{p:.1f}" for p in ph_values))
        radius_nominal = getattr(cfg, "PH_RADIUS", 10.0)
        try:
            radius_nominal = float(radius_nominal)
        except Exception:
            radius_nominal = 10.0
        scope_cfg = getattr(cfg, "PH_SCOPE", "")
        scope, center, eff_radius = _resolve_ph_scope(scope_cfg, radius_nominal, paths, cleaned_path, log)
        log.info("[ph_ensemble.pick.scope] scope=%s", scope)
        if scope == "pocket":
            log.info("[ph_ensemble.pick.center] center=(%.3f,%.3f,%.3f) radius=%.1f", center[0], center[1], center[2], radius_nominal)
        else:
            log.info("[ph_ensemble.pick.center] center=GLOBAL radius=ALL")
        log.info("[ph_ensemble.call] building pH ensemble for %s", paths.pdb_id)
        prev_cfg = getattr(automate_protein_prep, "config", None)
        try:
            automate_protein_prep.config = cfg
        except Exception:
            prev_cfg = None
        try:
            import ph_ensemble
            manifest_path = ph_ensemble.build_ph_ensemble(
                pdb_id=paths.pdb_id,
                cleaned_receptor_pdb=cleaned_path,
                out_dir=str(Path(cfg["OUTPUT_DIR"])),
                center=center,
                radius=eff_radius,
                ph_values=ph_values,
                variant=variant,
                legacy=bool(cfg.get("_ROUTER_LEGACY", False)),
            )

            log.info("[ph_ensemble.manifest] path=%s", manifest_path)
            log.info("[ph.manifest.json] written=%s", manifest_path)

            # --- DEBUG: measure map size before/after bridge ---
            try:
                _pre = len((cfg.get("_PH_ENSEMBLE_CANONICAL") or {}).get(paths.pdb_id, []))
            except Exception:
                _pre = -1
            log.info("[ph_ensemble.debug] before-bridge map_len[%s]=%d", paths.pdb_id, _pre)

            if manifest_path:
                _bridge_manifest_targets(str(manifest_path))

            try:
                _post = len((cfg.get("_PH_ENSEMBLE_CANONICAL") or {}).get(paths.pdb_id, []))
            except Exception:
                _post = -1
            log.info("[ph_ensemble.debug] after-bridge map_len[%s]=%d", paths.pdb_id, _post)
            # ---------------------------------------------------

            log.info("[ph_ensemble.done] ok=True")
            return manifest_path




        except Exception as exc:
            log.error("[ph_ensemble.error] %s", exc)
            log.info("[ph_ensemble.done] ok=False")
            return None
        finally:
            if prev_cfg is not None:
                automate_protein_prep.config = prev_cfg

    if cleaned_pdb_path.exists() and receptor_pdbqt_path.exists() and not force_reprocess:
        logger.info("Reusing existing cleaned PDB and receptor PDBQT.")
        try:
            if bool(cfg.get("RECEPTOR_SANITY_CHECK", True)) and not receptor_sanity_check(str(receptor_pdbqt_path)):
                logger.warning("Receptor sanity check failed (cached receptor).")
                prepare_receptor.last_provenance = "cache_reuse_failed"
                return None, None
        except Exception as _e:
            logger.warning(f"Receptor sanity check skipped due to error: {_e}")
        cleaned_norm = norm(cleaned_pdb_path)
        receptor_norm = norm(receptor_pdbqt_path)
        if bool(cfg.get("PH_ENSEMBLE_IN_PREP", False)):
            _build_ph_ensemble(cleaned_norm)
            mode = _ph_ligand_mode(cfg)
            if mode != "off":
                try:
                    context_ph_values = _ph_values_from_context(cleaned_pdb_path)
                except Exception as e:
                    logger.warning("[ph_ligand] failed to load context pH values; skipping ligand enumeration: %s", e)
                    context_ph_values = []

                if context_ph_values:
                    window: set[float] = set()
                    for p in context_ph_values:
                        if p is None:
                            continue
                        try:
                            p_val = float(p)
                        except Exception:
                            continue
                        for delta in (-1.0, 0.0, +1.0):
                            v = p_val + delta
                            if v < 3.0 or v > 10.5:
                                continue
                            window.add(round(v, 1))

                    ligand_ph_values = sorted(window)
                    if ligand_ph_values:
                        logger.info(
                            "[ph_ligand] mode=%s context_pH=%s ligand_pH_window=%s",
                            mode,
                            ",".join(f"{p:.1f}" for p in sorted(context_ph_values)),
                            ",".join(f"{p:.1f}" for p in ligand_ph_values),
                        )
                        _ctrl_roots, noncontrol_roots = _lib_roots_for_pdb(
                            cfg, paths.pdb_id.upper(), paths, logger
                        )
                        ph_ligand_root = noncontrol_roots[0] if noncontrol_roots else None
                        logger.info(
                            "[ph_ligand.context.bridge] pdb=%s root_dir=%s requested_ph=%s",
                            paths.pdb_id,
                            ph_ligand_root,
                            ligand_ph_values,
                        )
                        logger.debug(
                            "[ph_ligand.debug.call_enum] pdb_id=%s library_root=%s n_ph=%d ph_values=%s",
                            paths.pdb_id,
                            ph_ligand_root,
                            len(ligand_ph_values),
                            ",".join(f"{p:.2f}" for p in ligand_ph_values),
                        )
                        if ph_ligand_root is None or not ph_ligand_root.exists():
                            logger.info(
                                "[ph_ligand.context.bridge.skip] no valid ligand root; skipping microstate priming",
                            )
                        else:
                            try:
                                enumerate_ligands_for_docking(
                                    requested_ph_values=ligand_ph_values,
                                    root_dir=ph_ligand_root,
                                    microstate_dedup=True,
                                    force=False,
                                    cfg=cfg,
                                    pdb_id=paths.pdb_id,
                                )
                            except Exception as e:
                                logger.warning("[ph_ligand] ligand enumeration failed (non-fatal): %s", e)
                    else:
                        logger.info("[ph_ligand] context pH values present but window is empty after clamping; skipping ligand enumeration")
                else:
                    logger.info("[ph_ligand] no context pH values available; ligand enumeration skipped")
        prepare_receptor.last_provenance = "cache_reuse"
        return cleaned_norm, receptor_norm
    
    # --- PH_ENSEMBLE gating of legacy protonation ---
    if bool(cfg.get("PH_ENSEMBLE", False)):
        os.environ["A2_SKIP_PDB2PQR"] = "1"
        logger.info("[ph_ensemble] enabling ensemble mode: A2_SKIP_PDB2PQR=1 for cleaning stage")
    else:
        os.environ.pop("A2_SKIP_PDB2PQR", None)
        logger.info("[ph_ensemble] disabled; legacy cleaning path unchanged")
    # ------------------------------------------------
    # Fresh prep path: clean PDB then create PDBQT into the variant-aware target
    try:
        cleaned_pdb = automate_protein_prep.clean_pdb(
            pdb_file=str(paths.input_pdb_path),
            output_root=str(Path(cfg["OUTPUT_DIR"])),  # processed_pdbs root; module lays out subdirs
            logger=logger,
        )
    except Exception as e:
        logger.warning(f"Protein cleaning failed: {e}")
        prepare_receptor.last_provenance = "clean_failed"
        return None, None

    if not cleaned_pdb or not Path(cleaned_pdb).exists():
        logger.warning("Protein cleaning did not produce a cleaned PDB.")
        prepare_receptor.last_provenance = "clean_failed"
        return None, None

    # Relocate cleaned PDB into the variant receptor dir if needed
    try:
        if Path(cleaned_pdb).resolve() != cleaned_pdb_path.resolve():
            cleaned_pdb_path.parent.mkdir(parents=True, exist_ok=True)
            from shutil import copy2
            copy2(str(cleaned_pdb), str(cleaned_pdb_path))
            cleaned_pdb = str(cleaned_pdb_path)
        else:
            cleaned_pdb = str(cleaned_pdb_path)
    except Exception as e:
        logger.warning(f"Could not relocate cleaned PDB: {e}")

    # Build protein pH ensemble if requested
    ph_ensemble_in_prep = bool(cfg.get("PH_ENSEMBLE_IN_PREP", False))
    if ph_ensemble_in_prep:
        _build_ph_ensemble(cleaned_pdb)

    # Optional ligand microstate priming based on context pH
    mode = _ph_ligand_mode(cfg)
    if ph_ensemble_in_prep and mode != "off":
        try:
            # Use the existing context-pH helper to recover the canonical pH values
            # used to build the ensemble, instead of parsing any filenames.
            context_ph_values = _ph_values_from_context(cleaned_pdb_path)
        except Exception as e:
            logger.warning("[ph_ligand] failed to load context pH values; skipping ligand enumeration: %s", e)
            context_ph_values = []

        if context_ph_values:
            # Build a ±1 pH window union across all context pHs, clamped to [3.0, 10.5]
            window: set[float] = set()
            for p in context_ph_values:
                if p is None:
                    continue
                try:
                    p_val = float(p)
                except Exception:
                    continue
                for delta in (-1.0, 0.0, +1.0):
                    v = p_val + delta
                    if v < 3.0 or v > 10.5:
                        continue
                    # round to 1 decimal place to match microstate registry convention
                    window.add(round(v, 1))

            ligand_ph_values = sorted(window)
            if ligand_ph_values:
                logger.info(
                    "[ph_ligand] mode=%s context_pH=%s ligand_pH_window=%s",
                    mode,
                    ",".join(f"{p:.1f}" for p in sorted(context_ph_values)),
                    ",".join(f"{p:.1f}" for p in ligand_ph_values),
                )
                _ctrl_roots, noncontrol_roots = _lib_roots_for_pdb(cfg, paths.pdb_id.upper(), paths, logger)
                ph_ligand_root = noncontrol_roots[0] if noncontrol_roots else None

                if ligand_ph_values and ph_ligand_root is not None and ph_ligand_root.exists():
                    logger.info(
                        "[ph_ligand.context.bridge] pdb=%s root_dir=%s requested_ph=%s",
                        paths.pdb_id,
                        ph_ligand_root,
                        ligand_ph_values,
                    )
                    logger.debug(
                        "[ph_ligand.debug.call_enum] pdb_id=%s library_root=%s n_ph=%d ph_values=%s",
                        paths.pdb_id,
                        ph_ligand_root,
                        len(ligand_ph_values),
                        ",".join(f"{p:.2f}" for p in ligand_ph_values),
                    )
                    try:
                        enumerate_ligands_for_docking(
                            requested_ph_values=ligand_ph_values,
                            root_dir=ph_ligand_root,
                            microstate_dedup=True,
                            force=False,
                        )
                    except Exception as e:
                        logger.warning("[ph_ligand] ligand enumeration failed (non-fatal): %s", e)
                else:
                    logger.info(
                        "[ph_ligand.context.bridge.skip] no valid ligand root or empty window; skipping microstate priming"
                    )
            else:
                logger.info("[ph_ligand] context pH values present but window is empty after clamping; skipping ligand enumeration")
        else:
            logger.info("[ph_ligand] no context pH values available; ligand enumeration skipped")

    try:
        provenance = getattr(automate_protein_prep, "get_clean_provenance", lambda: "clean_pdb")()
    except Exception:
        provenance = "clean_pdb"
    prepare_receptor.last_provenance = provenance

    # Generate receptor PDBQT directly at the variant-aware path
    try:
        ok = automate_protein_prep.run_prepare_receptor(
            input_pdb=cleaned_pdb,
            output_pdbqt=str(receptor_pdbqt_path),
            cfg=cfg
        )
        receptor_pdbqt = str(receptor_pdbqt_path) if ok else None
    except Exception as e:
        logger.warning(f"Receptor PDBQT prep failed: {e}")
        receptor_pdbqt = None

    if not receptor_pdbqt or not Path(receptor_pdbqt).exists():
        logger.warning("Receptor PDBQT was not created.")
        prepare_receptor.last_provenance = "clean_failed"
        return None, None

    try:
        if Path(receptor_pdbqt).resolve() != receptor_pdbqt_path.resolve():
            from shutil import copy2
            receptor_pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
            copy2(receptor_pdbqt, receptor_pdbqt_path)
            receptor_pdbqt = str(receptor_pdbqt_path)
    except Exception as e:
        logger.warning(f"Could not relocate receptor PDBQT: {e}")

    try:
        if bool(cfg.get("RECEPTOR_SANITY_CHECK", True)):
            ok = receptor_sanity_check(receptor_pdbqt)
            if not ok:
                logger.warning("Receptor sanity check failed (too few atoms or zero coords).")
                prepare_receptor.last_provenance = "clean_failed"
                return None, None
    except Exception as _e:
        logger.warning(f"Receptor sanity check skipped due to error: {_e}")


    return norm(cleaned_pdb), norm(receptor_pdbqt)
