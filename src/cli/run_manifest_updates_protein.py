"""Run-manifest protein-level update entrypoints."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Mapping, Optional

from cli.run_manifest_distributed import distributed_mutate_protein_entry as _distributed_mutate_protein_entry
from cli.run_manifest_io import (
    get_manifest_paths,
    load_manifest as _load_manifest,
    manifest_lock_path as _manifest_lock_path,
    write_manifest as _write_manifest,
)
from cli.run_manifest_mutations import (
    ensure_protein as _ensure_protein,
    mutate_entry_druggability_and_engine_plan as _mutate_entry_druggability_and_engine_plan,
    mutate_entry_protein_failure as _mutate_entry_protein_failure,
    mutate_entry_protein_start as _mutate_entry_protein_start,
    mutate_entry_protein_success as _mutate_entry_protein_success,
)
from cli.run_manifest_support import exclusive_file_lock as _exclusive_file_lock, refresh_summary as _refresh_summary

def update_manifest_for_druggability_and_engine_plan(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: Optional[str],
    *,
    ph_tag: Optional[str],
    tier: Optional[str],
    use_gnina: bool,
    use_ledock: bool,
    use_dock6: bool,
) -> None:
    """
    Record protein difficulty tier (A/B/C) and planned engines for (pdb, variant, pH).

    Best-effort only: errors are logged and ignored.
    """
    ph_label = ph_tag if ph_tag is not None else "base"
    try:
        if not run_id:
            logging.debug(
                "[run-manifest.druggability-plan.skip] no run_id pdb=%s variant=%s ph=%s",
                pdb_id,
                variant_label,
                ph_label,
            )
            return

        if _distributed_mutate_protein_entry(
            cfg,
            run_id,
            pdb_id,
            variant_label,
            ph_tag,
            lambda entry: _mutate_entry_druggability_and_engine_plan(
                entry,
                tier=tier,
                use_gnina=use_gnina,
                use_ledock=use_ledock,
                use_dock6=use_dock6,
            ),
        ):
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        with _exclusive_file_lock(_manifest_lock_path(manifest_path)):
            manifest = _load_manifest(manifest_path)
            if manifest is None:
                logging.warning(
                    "[run-manifest.druggability-plan.skip] manifest missing run_id=%s path=%s pdb=%s variant=%s ph=%s",
                    run_id,
                    manifest_path,
                    pdb_id,
                    variant_label,
                    ph_label,
                )
                return

            entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
            tier_norm, engines = _mutate_entry_druggability_and_engine_plan(
                entry,
                tier=tier,
                use_gnina=use_gnina,
                use_ledock=use_ledock,
                use_dock6=use_dock6,
            )

            _refresh_summary(manifest)
            _write_manifest(manifest_path, manifest)

        logging.debug(
            "[run-manifest.druggability-plan.ok] run_id=%s pdb=%s variant=%s ph=%s tier=%s engines=%r",
            run_id,
            pdb_id,
            variant_label,
            ph_label,
            tier_norm,
            engines,
        )
    except Exception:
        logging.warning(
            "[run-manifest.druggability-plan.error] run_id=%s pdb=%s variant=%s ph=%s",
            run_id,
            pdb_id,
            variant_label,
            ph_label,
            exc_info=True,
        )

def update_manifest_for_protein_start(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: str,
    library: Optional[str],
    *,
    ph_tag: Optional[str] = None,
) -> None:
    try:
        if _distributed_mutate_protein_entry(
            cfg,
            run_id,
            pdb_id,
            variant_label,
            ph_tag,
            lambda entry: _mutate_entry_protein_start(entry, cfg, variant_label, library, ph_tag),
        ):
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        lock_path = _manifest_lock_path(manifest_path)
        with _exclusive_file_lock(lock_path):
            manifest = _load_manifest(manifest_path)
            if manifest is None:
                return

            entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
            _mutate_entry_protein_start(entry, cfg, variant_label, library, ph_tag)

            _refresh_summary(manifest)
            _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to record start pdb=%s variant=%s ph=%s run_id=%s",
            pdb_id,
            variant_label,
            ph_tag if ph_tag is not None else "base",
            run_id,
            exc_info=True,
        )


def update_manifest_for_protein_success(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: str,
    elapsed_sec: float,
    *,
    ph_tag: Optional[str] = None,
) -> None:
    try:
        if _distributed_mutate_protein_entry(
            cfg,
            run_id,
            pdb_id,
            variant_label,
            ph_tag,
            lambda entry: _mutate_entry_protein_success(entry, elapsed_sec),
        ):
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        lock_path = _manifest_lock_path(manifest_path)
        with _exclusive_file_lock(lock_path):
            manifest = _load_manifest(manifest_path)
            if manifest is None:
                return

            entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
            _mutate_entry_protein_success(entry, elapsed_sec)
            _refresh_summary(manifest)
            _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to record success pdb=%s variant=%s ph=%s run_id=%s",
            pdb_id,
            variant_label,
            ph_tag if ph_tag is not None else "base",
            run_id,
            exc_info=True,
        )

def update_manifest_for_protein_failure(
    cfg: Mapping[str, Any],
    run_id: str,
    pdb_id: str,
    variant_label: str,
    fail_log_path: Optional[Path],
    *,
    ph_tag: Optional[str] = None,
) -> None:
    try:
        if _distributed_mutate_protein_entry(
            cfg,
            run_id,
            pdb_id,
            variant_label,
            ph_tag,
            lambda entry: _mutate_entry_protein_failure(entry, variant_label, fail_log_path),
        ):
            return

        _, manifest_path = get_manifest_paths(cfg, run_id)
        lock_path = _manifest_lock_path(manifest_path)
        with _exclusive_file_lock(lock_path):
            manifest = _load_manifest(manifest_path)
            if manifest is None:
                return

            entry = _ensure_protein(manifest, pdb_id, variant_label, ph_tag)
            _mutate_entry_protein_failure(entry, variant_label, fail_log_path)

            _refresh_summary(manifest)
            _write_manifest(manifest_path, manifest)
    except Exception:
        logging.warning(
            "[run-manifest] Failed to record failure pdb=%s variant=%s ph=%s run_id=%s",
            pdb_id,
            variant_label,
            ph_tag if ph_tag is not None else "base",
            run_id,
            exc_info=True,
        )
