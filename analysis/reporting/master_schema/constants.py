# -*- coding: utf-8 -*-
import csv
import re

# Try to import canonical_ligand_base from rescore_reranker
try:
    from post_docking.rescoring.rescore_reranker import canonical_ligand_base
except ImportError:
    # Fallback if import fails (should not happen if rescore_reranker is in path)
    def canonical_ligand_base(lig: str) -> str:
        s = str(lig or "").strip()
        s = s.replace(".sanitized", "")
        s = re.sub(r"(_dud_gnina_stage\d+)$", "", s)
        s = re.sub(r"(_gnina_dud_stage\d+)$", "", s)
        s = re.sub(r"(_dock6_dud_stage\d+)$", "", s)
        s = re.sub(r"(_dud_dock6_stage\d+)$", "", s)
        s = re.sub(r"(_dud_stage\d+)$", "", s)
        s = re.sub(r"(__dock6_dud_stage\d+)$", "", s)
        s = re.sub(r"(__ledock_stage\d+)$", "", s)
        s = re.sub(r"(__dock6_stage\d+)$", "", s)
        s = re.sub(r"(_gnina_stage\d+)$", "", s)
        s = re.sub(r"(_stage\d+)$", "", s)
        s = re.sub(r"\.(mol2|pdbqt)$", "", s, flags=re.IGNORECASE)
        s = s.replace("__", "_")
        s = re.sub(r"_+$", "", s)
        return s


COMPONENT = "[master-export]"
DECOY_PREFIX_KEY = "DECOY_PREFIX"
DUD_PREFIX_KEY = "DUD_PREFIX"
DECOY_PREFIX_DEFAULT = "dud"
_SCHEMA_CSV_ERRORS: tuple[type[BaseException], ...] = (
    OSError,
    UnicodeDecodeError,
    csv.Error,
)
_SCHEMA_WRITE_ERRORS: tuple[type[BaseException], ...] = _SCHEMA_CSV_ERRORS + (
    ValueError,
    KeyError,
)
_DECOY_RE = re.compile(r"\bdecoys?_", re.IGNORECASE)
MIN_DECOYS_FOR_FDR = 200
MIN_UNIQUE_DECOY_SCORES = 10
