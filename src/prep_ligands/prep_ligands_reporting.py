"""Ligand-prep reporting and audit helpers shared across prep paths."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Union, cast

from rdkit import Chem

from protein_prep import pdb_fixer_runtime as _pdb_fixer


MALFORMED_LOG = Path("malformed_ligands.txt")
MALFORMED_DIR: Path | None = None


def _log_malformed(p: Path, reason: str, log_dir: Path | None = None) -> None:
    """Append a one-line malformed reason into the prep output root."""
    try:
        line = f"{p.name}\t{reason}\n"
        base = (
            Path(log_dir)
            if log_dir
            else (Path(MALFORMED_DIR) if MALFORMED_DIR else Path.cwd())
        )
        dest = base / MALFORMED_LOG.name
        dest.parent.mkdir(parents=True, exist_ok=True)
        with open(dest, "a", encoding="utf-8") as out:
            out.write(line)
    except Exception:
        pass


def _append_prep_status(
    status_log_path: Path,
    ligand_name: str,
    status: str,
    reason: str = "",
    relpath: str = "",
    *,
    stage: str = "",
    failure_code: str = "",
    failure_detail: str = "",
    fixes_count: int = 0,
    rules_ver: str = "",
    writer_final: str = "",
    rescue_used: str = "",
    torsion_root_rule: str = "",
    polarH: str = "",
    ad4_types_ok: str = "",
    charges_ok: str = "",
    torsdof: str = "",
) -> None:
    """Append a backward-compatible prep TSV row with provenance columns."""
    header = (
        "ligand\tstatus\treason\tpdbqt_rel\t"
        "stage\tfailure_code\tfailure_detail\tfixes_count\trules_version\t"
        "writer_final\trescue_used\ttorsion_root_rule\tpolarH\tad4_types_ok\tcharges_ok\ttorsdof\n"
    )
    if not status_log_path.exists():
        status_log_path.write_text(header, encoding="utf-8")
    with status_log_path.open("a", encoding="utf-8") as fh:
        fh.write(
            "\t".join(
                [
                    ligand_name,
                    status,
                    reason,
                    relpath,
                    stage,
                    failure_code,
                    failure_detail,
                    str(fixes_count),
                    rules_ver,
                    writer_final,
                    rescue_used,
                    torsion_root_rule,
                    polarH,
                    ad4_types_ok,
                    charges_ok,
                    str(torsdof or ""),
                ]
            )
            + "\n"
        )


def _log_elem_fix_summary(
    file_path: Path,
    stage: str,
    before_text: Optional[str] = None,
    after_text: Optional[str] = None,
) -> None:
    """Emit a compact He->H summary line for a rewritten file."""
    try:
        if before_text is None:
            before_text = Path(file_path).read_text(encoding="utf-8", errors="ignore")
        n_before = cast(int, _pdb_fixer.scan_helium_counts(before_text.splitlines()))
        if after_text is None:
            after_text = Path(file_path).read_text(encoding="utf-8", errors="ignore")
        n_after = cast(int, _pdb_fixer.scan_helium_counts(after_text.splitlines()))
        logging.info(
            "[elem-fix] file=%s stage=%s He->H=%d",
            file_path.name,
            stage,
            max(0, n_before - n_after),
        )
    except Exception as exc:
        logging.warning("[elem-fix] summary failed for %s: %s", file_path, exc)


def _audit_protonation_metrics(
    tag: str, mol_or_path: Union["Chem.Mol", str, Path], context: str = "pre"
) -> Dict[str, Any]:
    """Collect lightweight protonation/charge/aromaticity metrics for debug logs."""
    notes: List[str] = []
    h_count = 0
    formal_charge = 0
    has_partial = False
    aromatic_atoms = 0
    ad4_ok: Optional[bool] = None

    mol: Optional[Chem.Mol] = None
    try:
        if hasattr(mol_or_path, "GetNumAtoms"):
            mol = cast(Chem.Mol, mol_or_path)
        else:
            path_str = str(mol_or_path)
            ext = os.path.splitext(path_str)[1].lower()
            if ext == ".sdf":
                suppl = Chem.SDMolSupplier(path_str, sanitize=False, removeHs=False)
                mol = next((m for m in suppl if m), None)
            elif ext == ".mol2":
                mol = Chem.MolFromMol2File(path_str, sanitize=False, removeHs=False)
            elif ext == ".pdb":
                mol = Chem.MolFromPDBFile(path_str, sanitize=False, removeHs=False)
        if mol is not None:
            try:
                Chem.SanitizeMol(mol, catchErrors=True)
            except Exception:
                notes.append("sanitize_warn")
            try:
                formal_charge = Chem.GetFormalCharge(mol)
            except Exception:
                notes.append("formal_charge_calc_warn")
            try:
                h_count = sum(
                    atom.GetNumExplicitHs() + atom.GetTotalNumHs()
                    for atom in mol.GetAtoms()
                )
            except Exception:
                notes.append("h_count_warn")
            try:
                aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
            except Exception:
                notes.append("aroma_warn")
    except Exception:
        notes.append("rdkit_missing")

    if mol is None and isinstance(mol_or_path, (str, Path)):
        try:
            has_q = False
            ad4_seen = False
            h_count_fallback = 0
            with open(mol_or_path, "r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if not line.startswith(("ATOM", "HETATM")):
                        continue
                    atom_name = line[12:16].strip().upper()
                    elem = line[76:78].strip().upper()
                    if atom_name.startswith("H") or elem == "H":
                        h_count_fallback += 1
                    tokens = line.split()
                    if len(tokens) >= 9:
                        try:
                            float(tokens[8])
                            has_q = True
                        except Exception:
                            pass
                    if "A" in line[-6:]:
                        ad4_seen = True
            h_count = h_count or h_count_fallback
            has_partial = has_q
            ad4_ok = ad4_seen
        except Exception:
            notes.append("pdbqt_parse_warn")

    return {
        "ok": True,
        "h_count": int(h_count),
        "formal_charge": int(formal_charge),
        "has_partial_charges": bool(has_partial),
        "aromatic_atoms": int(aromatic_atoms),
        "ad4_types_ok": ad4_ok,
        "notes": notes,
        "tag": tag,
        "context": context,
    }


__all__ = [
    "MALFORMED_DIR",
    "MALFORMED_LOG",
    "_append_prep_status",
    "_audit_protonation_metrics",
    "_log_elem_fix_summary",
    "_log_malformed",
]
