from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Mapping


METHOD_KEYS = (
    "MMGBSA_PROTOCOL",
    "MMGBSA_MD_ENABLED",
    "MMGBSA_MD_SOLVENT_MODEL",
    "MMGBSA_MD_NREPLICATES",
    "MMGBSA_MD_HEAT_PS",
    "MMGBSA_MD_EQUIL_PS",
    "MMGBSA_MD_PROD_PS",
    "MMGBSA_MD_FRAME_STRIDE_PS",
    "MMGBSA_ANALYSIS_START_PS",
    "MMGBSA_ANALYSIS_END_PS",
    "MMGBSA_ANALYSIS_INTERVAL",
    "MMGBSA_GB_IGB",
    "MMGBSA_GB_SALTCON",
    "MMGBSA_FRAME_AGG",
    "MMGBSA_MD_REP_AGG",
    "MMGBSA_KEEP_WATERS",
    "MMGBSA_KEEP_METALS",
    "MMGBSA_CHEMISTRY_REVIEW_REQUIRED",
    "MMGBSA_CHEMISTRY_REVIEW_STATUS",
    "MMGBSA_CHEMISTRY_REVIEW_NOTES",
    "MMGBSA_CHEMISTRY_SUBGATES_REQUIRED",
    "MMGBSA_PROTONATION_REVIEW_STATUS",
    "MMGBSA_TAUTOMER_REVIEW_STATUS",
    "MMGBSA_STEREOCHEMISTRY_REVIEW_STATUS",
    "MMGBSA_NET_CHARGE_REVIEW_STATUS",
    "MMGBSA_PARAMETER_REVIEW_STATUS",
    "MMGBSA_LIGAND_AT",
    "MMGBSA_LIGAND_CHARGE_METHOD",
    "MMGBSA_LIGAND_PRIMARY_CHARGE_METHOD",
    "MMGBSA_LIGAND_FALLBACK_CHARGE_METHOD",
    "MMGBSA_LIGAND_CHEMISTRY_STRICT",
    "MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE",
    "MMGBSA_INPUT_CHEMISTRY_SOURCE",
    "MMGBSA_INPUT_CHEMISTRY_NOTES",
)


def write_methods_json(
    path: Path,
    cfg: Mapping[str, Any],
    *,
    run_id: str,
    pdb_id: str,
    variant: str,
    ph_label: str,
) -> None:
    payload = {
        "run_id": run_id,
        "pdb_id": pdb_id,
        "variant": variant,
        "ph_label": ph_label,
        "config": {key: cfg.get(key) for key in METHOD_KEYS if key in cfg},
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".part")
    with tmp_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp_path, path)
