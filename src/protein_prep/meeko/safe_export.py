"""Conservative Meeko receptor-export rescue inputs.

These helpers keep the prepared PDB intact, but create a Meeko-only receptor
input when retained HET records, unsupported metals, or malformed local residue
geometry prevent PDBQT export.
"""

from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path
from typing import Iterable

from protein_prep.geometry.repair import repair_terminal_oxt_geometry_in_pdb

WATER_NAMES = {"HOH", "WAT", "H2O", "DOD"}
PROTEIN_RECORD_RESNAMES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}
PROTEIN_MEEKO_ALIASES = {
    "ASH": "ASP",
    "CYM": "CYS",
    "CYX": "CYS",
    "GLH": "GLU",
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "LYN": "LYS",
}

ResidueDropKey = tuple[str, str]


def parse_unmatched_residue_keys(error_text: str) -> set[ResidueDropKey]:
    """Parse Meeko residue_key='A:123' template failures."""

    keys: set[ResidueDropKey] = set()
    for match in re.finditer(r"residue_key='([^:']+):([^']+)'", error_text):
        chain = match.group(1).strip() or "-"
        resseq = match.group(2).strip()
        if resseq:
            keys.add((chain, resseq))
    return keys


def write_meeko_safe_receptor_copy(
    input_pdb: Path,
    output_pdb: Path,
    *,
    drop_residue_keys: Iterable[ResidueDropKey] = (),
) -> dict[str, object]:
    """Write a Meeko-only ATOM receptor input and audit every dropped record."""

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    drop_keys = set(drop_residue_keys)
    dropped: list[dict[str, str]] = []
    kept: list[str] = []
    with input_pdb.open("r", encoding="utf-8", errors="ignore") as src:
        for line in src:
            normalized, reason = _meeko_safe_line(line, drop_keys)
            if reason:
                dropped.append(_drop_row(line, reason))
                continue
            if normalized is not None:
                kept.append(normalized)
    output_pdb.write_text("".join(kept).rstrip() + "\nEND\n", encoding="utf-8")
    repair_terminal_oxt_geometry_in_pdb(
        output_pdb,
        audit_path=output_pdb.with_suffix(".oxt_repair_audit.json"),
    )
    summary = _summary(input_pdb, output_pdb, dropped, drop_keys)
    output_pdb.with_suffix(".drop_audit.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return summary


def meeko_safe_stage_label(summary: dict[str, object]) -> str:
    dropped_atoms = _int_summary_value(summary.get("dropped_atom_count"))
    dropped_residues = _int_summary_value(summary.get("dropped_residue_count"))
    label = f"protein_only_sanitized_retry:dropped_atoms={dropped_atoms}"
    if dropped_residues:
        label += f",dropped_residues={dropped_residues}"
    return label


def _int_summary_value(value: object) -> int:
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.strip().isdigit():
        return int(value)
    return 0


def _meeko_safe_line(
    line: str,
    drop_residue_keys: set[ResidueDropKey],
) -> tuple[str | None, str]:
    if line.startswith("HETATM"):
        resname = line[17:20].strip().upper()
        return None, "meeko_safe_drop_water" if resname in WATER_NAMES else "meeko_safe_drop_hetatm"
    if line.startswith(("CONECT", "USER  MOD")):
        return None, "meeko_safe_drop_connectivity_record"
    if not line.startswith("ATOM  "):
        return (line if line.startswith("TER") else None), ""
    if _line_element(line) in {"H", "D"}:
        return None, "meeko_safe_drop_hydrogen"
    key = (line[21:22].strip() or "-", line[22:26].strip())
    if key in drop_residue_keys:
        return None, "meeko_safe_drop_template_failed_residue"
    resname = _normalize_resname_for_meeko(line[17:21].strip())
    if resname not in PROTEIN_RECORD_RESNAMES:
        return None, "meeko_safe_drop_nonstandard_atom_residue"
    return _rewrite_resname(line, resname), ""


def _normalize_resname_for_meeko(resname: str) -> str:
    token = resname.strip().upper()
    if len(token) == 4 and token[0] in {"N", "C"}:
        token = token[1:]
    return PROTEIN_MEEKO_ALIASES.get(token, token)


def _rewrite_resname(line: str, resname: str) -> str:
    return f"{line[:17]}{resname.rjust(3)}{line[20:]}"


def _line_element(line: str) -> str:
    return (line[76:78].strip() or line[12:16].strip()[:1]).upper()


def _drop_row(line: str, reason: str) -> dict[str, str]:
    return {
        "record": line[:6].strip(),
        "chain": line[21:22].strip() or "-",
        "resseq": line[22:26].strip(),
        "icode": line[26:27].strip(),
        "resname": line[17:21].strip(),
        "atom_name": line[12:16].strip(),
        "element": _line_element(line) if line.startswith(("ATOM", "HETATM")) else "",
        "reason": reason,
    }


def _summary(
    input_pdb: Path,
    output_pdb: Path,
    dropped: list[dict[str, str]],
    drop_keys: set[ResidueDropKey],
) -> dict[str, object]:
    reasons = Counter(row["reason"] for row in dropped)
    return {
        "input_pdb": str(input_pdb),
        "output_pdb": str(output_pdb),
        "dropped_atom_count": len(dropped),
        "dropped_residue_count": len(drop_keys),
        "dropped_residue_keys": [f"{chain}:{resseq}" for chain, resseq in sorted(drop_keys)],
        "drop_reasons": dict(sorted(reasons.items())),
        "dropped_atoms": dropped,
    }


__all__ = [
    "meeko_safe_stage_label",
    "parse_unmatched_residue_keys",
    "write_meeko_safe_receptor_copy",
]
