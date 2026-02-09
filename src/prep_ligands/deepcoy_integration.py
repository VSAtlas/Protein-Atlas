# -*- coding: utf-8 -*-
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, Sequence


def _normalize_flag_name(tok: str) -> str:
    return str(tok).lstrip("-").lower()


def _parse_force_deepcoy_flag(argv: Sequence[str]) -> bool:
    try:
        return any(
            _normalize_flag_name(tok) == "force-deepcoy"
            for tok in argv
            if str(tok).startswith("-")
        )
    except Exception:
        return False


def apply_deepcoy_defaults(cfg: Dict[str, object]) -> None:
    overall_dir = Path(
        str(
            cfg.get("OVERALL_DIR")
            or os.environ.get("PROTEIN_AUTOMATION_DIR")
            or Path(__file__).resolve().parents[2]
        )
    )
    deepcoy_duds_dir = os.environ.get("DEEPCOY_DUDS_DIR") or str(
        overall_dir / "DeepCoy_duds"
    )
    input_pdb_dir = (
        os.environ.get("DEEPCOY_INPUT_PDB_DIR")
        or cfg.get("INPUT_DIR")
        or str(overall_dir / "input_pdbs")
    )
    extracted_root = str(
        cfg.get("EXTRACTED_LIGANDS_DIR") or (overall_dir / "extracted_ligands")
    )
    deepcoy_python = os.environ.get("DEEPCOY_PYTHON") or sys.executable

    cfg.setdefault("USE_DEEPCOY", "on")
    cfg.setdefault("DEEPCOY_DUDS_DIR", deepcoy_duds_dir)
    cfg.setdefault("DEEPCOY_PYTHON", deepcoy_python)
    cfg.setdefault("DEEPCOY_INPUT_PDB_DIR", input_pdb_dir)
    cfg.setdefault(
        "DEEPCOY_OUT_ROOT",
        os.path.join(extracted_root, "deepcoy"),
    )
    cfg.setdefault(
        "DEEPCOY_WORK_ROOT",
        os.path.join(str(cfg["DEEPCOY_DUDS_DIR"]), "deepcoy_work"),
    )
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


def apply_deepcoy_cli_overrides(cfg: Dict[str, object], argv: Sequence[str]) -> None:
    if _parse_force_deepcoy_flag(argv):
        cfg["DEEPCOY_FORCE"] = "on"
