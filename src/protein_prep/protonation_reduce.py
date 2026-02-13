"""Backward-compatible import shim for the protonation module."""

from __future__ import annotations

import subprocess

from protein_prep.protonation import (
    DEFAULT_HET,
    OPENBABEL_PATH,
    PHENIX_DIR,
    REDUCE_EXE,
    _REPO_ROOT,
    _het_dict_path,
    _maybe_get_target_ph,
    _pick_reduce_exe,
    _protonate_with_pdb2pqr_if_available,
    assign_protonation_states,
    conect_coverage,
    hydrogenation_status,
    run_openbabel_add_h,
)

__all__ = [
    "PHENIX_DIR",
    "OPENBABEL_PATH",
    "REDUCE_EXE",
    "DEFAULT_HET",
    "_REPO_ROOT",
    "_pick_reduce_exe",
    "_het_dict_path",
    "hydrogenation_status",
    "conect_coverage",
    "_maybe_get_target_ph",
    "_protonate_with_pdb2pqr_if_available",
    "assign_protonation_states",
    "run_openbabel_add_h",
    "subprocess",
]
