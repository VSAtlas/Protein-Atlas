# -*- coding: utf-8 -*-
from __future__ import annotations

import os
from typing import Dict, Sequence

from input_and_export_functions import _to_bool


def _normalize_flag_name(tok: str) -> str:
    return str(tok).lstrip("-").lower()


def _parse_force_deepcoy_flag(argv: Sequence[str]) -> bool:
    try:
        return any(_normalize_flag_name(tok) == "force-deepcoy" for tok in argv if str(tok).startswith("-"))
    except Exception:
        return False


def apply_deepcoy_defaults(cfg: Dict[str, str]) -> None:
    cfg.setdefault("DEEPCOY_DUDS_DIR", "/home/michael/atlas/code/protein_automation/DeepCoy_duds")
    cfg.setdefault("DEEPCOY_PYTHON", "/home/michael/atlas/anaconda3/envs/DeepCoy-env/bin/python")
    cfg.setdefault("DEEPCOY_INPUT_PDB_DIR", "/home/michael/atlas/code/protein_automation/input_pdbs")
    cfg.setdefault(
        "DEEPCOY_OUT_ROOT",
        os.path.join(
            cfg.get("LIGAND_EXTRACTED_DIR", "/home/michael/atlas/code/protein_automation/extracted_ligands"),
            "deepcoy",
        ),
    )
    cfg.setdefault("DEEPCOY_WORK_ROOT", os.path.join(cfg["DEEPCOY_DUDS_DIR"], "deepcoy_work"))
    cfg.setdefault("DEEPCOY_ENABLE_AUTOGEN_SDF", "on")
    cfg.setdefault("DEEPCOY_ENABLE_AUTOGEN_PDBQT", "on")
    cfg.setdefault("DEEPCOY_FALLBACK_SMILES", "")
    cfg.setdefault("DEEPCOY_RESTRICT_DATA", 0)
    cfg.setdefault("DEEPCOY_PREPPED_SUBDIR", "deepcoy")
    cfg.setdefault("DEEPCOY_LIGPREP_FORCE", "off")
    cfg.setdefault("DEEPCOY_USE_AS_DUD_LIBRARY", "on")
    cfg.setdefault("DEEPCOY_FORCE", "off")
    env_deepcoy_toggle = os.environ.get("DEEPCOY_ENABLE_AUTOGEN_SDF")
    if env_deepcoy_toggle is not None:
        cfg["DEEPCOY_ENABLE_AUTOGEN_SDF"] = env_deepcoy_toggle


def apply_deepcoy_cli_overrides(cfg: Dict[str, str], argv: Sequence[str]) -> None:
    if _parse_force_deepcoy_flag(argv):
        cfg["DEEPCOY_FORCE"] = "on"
