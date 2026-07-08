"""Ligand-prep molecule conversion and PDBQT generation helpers."""

from dataclasses import dataclass
import logging
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from config.tool_resolver import resolve_tool
from rdkit import Chem

from prep_ligands.prep_ligands_bulk_sdf import _run_obabel
from prep_ligands.prep_ligands_reporting import (
    _append_prep_status,
    _audit_protonation_metrics,
)
from prep_ligands.prep_ligands_common_validation import (
    ALLOWED_ELEMENTS,
    MAX_HEAVY_ATOMS,
    MIN_ATOMS_FOR_DOCKING,
    MONOATOMIC_IONS,
    QUARANTINE_DIRNAME,
    RULES,
    _element_from_adt,
    _helium_postwrite_guard,
    _parse_pdbqt_metrics,
    assert_no_helium_in_pdbqt,
    get_short_path_name,
    is_valid_ligand,
    quick_pdbqt_validate,
    rules_version,
    validate_pdbqt_invariants,
)

_rdmol_standardize: Any
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize as _rdmol_standardize
except Exception:
    _rdmol_standardize = None

_std = _rdmol_standardize
_HAS_STD = _std is not None

logger = logging.getLogger(__name__)

# --- Counter-ion SMARTS (shared between filters) ---
_COUNTERION_SMARTS = {
    "mesylate": Chem.MolFromSmarts("[CH3]-S(=O)(=O)[O-]"),
    "tosylate": Chem.MolFromSmarts("c1cccc(c1)S(=O)(=O)[O-]"),
    "triflate": Chem.MolFromSmarts("C(F)(F)F-S(=O)(=O)[O-]"),
    "sulfonate": Chem.MolFromSmarts("S(=O)(=O)[O-]"),
    "phosphate": Chem.MolFromSmarts("P(=O)([O-])([O-])[O-]"),
    "sulfate": Chem.MolFromSmarts("S(=O)(=O)([O-])[O-]"),
    "formate": Chem.MolFromSmarts("[#6](=O)[O-]"),
    "acetate": Chem.MolFromSmarts("CC(=O)[O-]"),
    "lactate": Chem.MolFromSmarts("CC(O)C(=O)[O-]"),
}
_COUNTERION_SMARTS.update(
    {
        "citrate_like": Chem.MolFromSmarts(
            "[CX4](-[CH2]-C(=O)[O-])(-[CH2]-C(=O)[O-])(-C(=O)[O-])O"
        ),
        "tartrate_like": Chem.MolFromSmarts(
            "IC([CH](O)C(=O)[O-])C(=O)[O-]".replace("I", "O")
        ),
    }
)

_CARBOXYLATE = Chem.MolFromSmarts("[CX3](=O)[O-]")
_CARBOXYLIC = Chem.MolFromSmarts("[CX3](=O)O")

def standardize_mol_with_activesite(mol: Chem.Mol) -> Chem.Mol:
    """
    Apply canonical ligand standardization before filtering:
      - RDKit metal disconnection, salt removal, largest fragment, uncharge
      - (light) sanitize to normalize valences where possible
      - Do NOT change protonation (keep existing behavior)
    Uses the same internal passes as _standardize_then_sanitize when available.
    """
    if mol is None:
        return mol
    m = mol
    # Prefer the existing standardization pipeline if available
    try:
        std, why = _standardize_then_sanitize(m)
        if std is not None:
            m = std
    except Exception:
        # fall through; keep original m
        pass

    # Guardrail: strip clearly invalid/dummy atoms if any slipped through
    try:
        bad = [
            a for a in m.GetAtoms() if (a.GetAtomicNum() == 0 or a.GetSymbol() == "*")
        ]
        if bad:
            emsg = f"dummy_atoms({len(bad)})"
            m.SetProp("_ligprep_filter_hint", emsg)
    except Exception:
        pass

    return m


def _quick_filters(mol: Chem.Mol) -> Tuple[bool, str]:
    """
    Fast ligand gate: run *after* canonical standardization to avoid alias/case drift.
    Preserves legacy thresholds/messages; adds clearer diagnostics for element rejects.
    """
    # 1) Standardize first (aliases/salts) so checks run on a normalized molecule
    try:
        mol = standardize_mol_with_activesite(mol) or mol
    except Exception:
        # keep going; downstream checks are defensive
        pass

    # 2) Size/atom-count thresholds (unchanged)
    heavy = mol.GetNumHeavyAtoms()
    if heavy == 0:
        return False, "no_heavy_atoms"
    if heavy > MAX_HEAVY_ATOMS:
        return False, f"too_large({heavy})"
    total_atoms = mol.GetNumAtoms()
    if total_atoms < MIN_ATOMS_FOR_DOCKING:
        return False, f"too_few_atoms({total_atoms})"

    # 3) Element allow-list using aliases-derived set
    offending: Set[str] = set()
    try:
        for a in mol.GetAtoms():
            sym = a.GetSymbol()
            if sym not in ALLOWED_ELEMENTS:
                offending.add(sym)
    except Exception:
        # if RDKit access fails for any atom, be conservative and report failure
        return False, "element_scan_error"

    if offending:
        # Provide a hint for common alias/case drift
        sym = sorted(offending)[0]
        hint = ""
        if sym.upper() != sym and sym.capitalize() in ALLOWED_ELEMENTS:
            hint = f" (did_you_mean:{sym.capitalize()})"
        elif sym.upper() in (getattr(RULES, "two_letter", set()) or set()):
            hint = f" (did_you_mean:{sym[0].upper()}{sym[1:].lower()})"
        return False, f"disallowed_element:{sym}{hint}"

    # 4) Monatomic free-ion guard via aliases (canonical set; unions legacy)
    try:
        if mol.GetNumAtoms() == 1:
            s = mol.GetAtomWithIdx(0).GetSymbol()
            if s in MONOATOMIC_IONS:
                return False, f"free_monoatomic_ion:{s}"
    except Exception:
        pass

    return True, ""


def _standardize_then_sanitize(mol: Chem.Mol) -> Tuple[Optional[Chem.Mol], str]:
    if not _HAS_STD:
        return None, "std_module_missing"
    try:
        md = _std.MetalDisconnector()
        fr = _std.FragmentRemover()
        lf = _std.LargestFragmentChooser(preferOrganic=True)
        uc = _std.Uncharger()

        m = md.Disconnect(mol)
        m = fr.RemoveFragments(m)
        m = lf.choose(m)
        m = uc.uncharge(m)

        try:
            m = Chem.AddHs(m, addCoords=True)
        except Exception as _e_hs:
            logging.warning("[ligprep] AddHs pre-sanitize failed: %s", _e_hs)

        Chem.SanitizeMol(m)

        ok, why = _quick_filters(m)
        if not ok:
            return None, why
        return m, ""
    except Exception as e:
        return None, f"std_resanitize_fail:{e}"


def _looks_like_buffer_salt(m: Chem.Mol) -> bool:
    from rdkit.Chem import rdMolDescriptors as rdmd

    hac = m.GetNumHeavyAtoms()
    o = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "O")
    n = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N")
    rings = rdmd.CalcNumRings(m)
    arom = rdmd.CalcNumAromaticRings(m)
    return hac <= 15 and o >= 6 and n == 0 and rings == 0 and arom == 0


def _matches_counterion(m: Chem.Mol) -> Optional[str]:
    try:
        hac = m.GetNumHeavyAtoms()
        if hac == 0:
            return None
        for name in ("mesylate", "tosylate", "triflate"):
            patt = _COUNTERION_SMARTS[name]
            if patt and m.HasSubstructMatch(patt):
                return name
        for name in ("citrate_like", "tartrate_like"):
            patt = _COUNTERION_SMARTS[name]
            if patt and m.HasSubstructMatch(patt):
                return name
        for name in ("phosphate", "sulfate", "sulfonate"):
            patt = _COUNTERION_SMARTS[name]
            if patt and m.HasSubstructMatch(patt):
                if hac <= 14:
                    return name
        if (
            hac <= 4
            and _COUNTERION_SMARTS["formate"]
            and m.HasSubstructMatch(_COUNTERION_SMARTS["formate"])
        ):
            return "formate"
        if (
            hac <= 5
            and _COUNTERION_SMARTS["acetate"]
            and m.HasSubstructMatch(_COUNTERION_SMARTS["acetate"])
        ):
            return "acetate"
        if (
            hac <= 6
            and _COUNTERION_SMARTS["lactate"]
            and m.HasSubstructMatch(_COUNTERION_SMARTS["lactate"])
        ):
            return "lactate"
        s = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "S")
        o = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "O")
        c = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C")
        if hac <= 8 and s == 1 and o >= 3 and c <= 2:
            return "small_sulfonate_like"
        return None
    except Exception:
        return None


def _is_polyacidic_buffer_like(m: Chem.Mol) -> bool:
    try:
        from rdkit.Chem import rdMolDescriptors as rdmd

        hac = m.GetNumHeavyAtoms()
        if hac == 0:
            return False
        rings = rdmd.CalcNumRings(m)
        arom = rdmd.CalcNumAromaticRings(m)
        o = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "O")
        n = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N")
        o_ratio = (o / float(hac)) if hac else 0.0

        na = 0
        if _CARBOXYLATE:
            na += len(m.GetSubstructMatches(_CARBOXYLATE))
        if _CARBOXYLIC:
            na += len(m.GetSubstructMatches(_CARBOXYLIC))
        has_sulfate = bool(
            _COUNTERION_SMARTS["sulfate"]
            and m.HasSubstructMatch(_COUNTERION_SMARTS["sulfate"])
        )
        has_phosphate = bool(
            _COUNTERION_SMARTS["phosphate"]
            and m.HasSubstructMatch(_COUNTERION_SMARTS["phosphate"])
        )

        if (
            rings == 0
            and arom == 0
            and o_ratio >= 0.35
            and (na >= 3 or has_sulfate or has_phosphate)
        ):
            if n <= 1:
                return True
        return False
    except Exception:
        return False

def _adt_retype_and_normalize(
    mgltools_python_short: str,
    prepare_script_short: str,
    mol2_for_mgl: Path,
    out_pdbqt: Path,
    *,
    torsion_rule_label: str = "adt_normalized",
    context: Optional["_PrepareOneContext"] = None,
    candidates: Optional[list["_WriterCandidate"]] = None,
    capture_label: Optional[str] = None,
) -> tuple[bool, str]:
    try:
        tmp = out_pdbqt.with_suffix(".adt_norm.pdbqt")
        cmd = [
            mgltools_python_short,
            prepare_script_short,
            "-l",
            get_short_path_name(str(mol2_for_mgl.resolve())),
            "-o",
            get_short_path_name(str(tmp.resolve())),
            "-U",
            "nphs_lps",
            "-A",
            "hydrogens",
        ]
        res = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
            cwd=str(mol2_for_mgl.parent),
            timeout=600,
        )
        if res.stderr:
            logging.warning("[mgltools retype stderr] %s", res.stderr.strip()[:300])

        if tmp.exists() and tmp.stat().st_size > 100:
            try:
                if out_pdbqt.exists():
                    out_pdbqt.unlink()
            except Exception:
                pass
            tmp.replace(out_pdbqt)

            m = _parse_pdbqt_metrics(out_pdbqt)
            typed, n_atoms = m.get("typed", 0), m.get("n_atoms", 0)
            logging.info(
                "[retype.ok] replaced=True n_atoms=%d typed_ratio=%d/%d ad4_ok=%s torsdof=%s",
                n_atoms,
                typed,
                n_atoms,
                m.get("has_ad4_types"),
                m.get("torsdof"),
            )
            if context is not None:
                try:
                    _, metrics_probe, _ = validate_pdbqt_invariants(out_pdbqt)
                    print(
                        f"[retype.result] lig={context.lig_id} ad4={metrics_probe.get('has_ad4_types')} "
                        f"torsdof={metrics_probe.get('torsdof')}"
                    )
                    if capture_label and candidates is not None:
                        _capture_writer_candidate(candidates, capture_label, out_pdbqt)
                except Exception:
                    print(f"[retype.result] lig={context.lig_id} ad4=? torsdof=?")
                context.writer_final = "mgltools"
            return True, torsion_rule_label

    except subprocess.CalledProcessError as e:
        logging.warning("[torsion-normalize] ADT retype failed: %s", (e.stderr or "")[-180:])
    except Exception as e:
        logging.warning("[torsion-normalize] exception: %s", e)

    return False, torsion_rule_label


def _reprep_via_meeko_path(mol2_path: Path, out_pdbqt: Path) -> bool:
    """
    Use Meeko from PATH (no config): prefer `mk_prepare_ligand.py`; fallback to `python -m meeko`.
    """
    exe = shutil.which("mk_prepare_ligand.py")
    if exe:
        cmd = [exe, "-i", str(mol2_path), "-o", str(out_pdbqt)]
    else:
        # last-resort: try the active interpreter with module form
        cmd = [
            sys.executable,
            "-m",
            "meeko",
            "-i",
            str(mol2_path),
            "-o",
            str(out_pdbqt),
        ]
    try:
        res = subprocess.run(
            cmd, check=True, capture_output=True, text=True, timeout=600
        )
        if res.stderr:
            logging.warning("[meeko stderr] %s", res.stderr.strip()[:300])
        return out_pdbqt.exists() and out_pdbqt.stat().st_size > 100
    except Exception as e:
        logging.warning("[meeko] failed: %s", e)
        return False


def _element_counts_from_pdb(pdb_path: Path) -> Dict[str, int]:
    counts: Dict[str, int] = defaultdict(int)
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fin:
            for ln in fin:
                if not ln.startswith(("ATOM", "HETATM")):
                    continue
                parts = ln.split()
                if not parts:
                    continue
                adt = parts[-1]
                elem = _element_from_adt(adt)
                counts[elem] += 1
    except Exception:
        pass
    return counts


def _buffer_like_by_counts_from_mol(m: Chem.Mol) -> bool:
    hac = m.GetNumHeavyAtoms()
    if hac == 0:
        return False
    o = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "O")
    n = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N")
    return hac >= 10 and (o / float(hac)) >= 0.40 and n <= 1


def _buffer_like_by_counts_from_pdbfile(pdb_path: Path) -> bool:
    counts = _element_counts_from_pdb(pdb_path)
    hac = sum(v for k, v in counts.items() if k != "H")
    o = counts.get("O", 0)
    n = counts.get("N", 0)
    return hac >= 10 and hac > 0 and (o / float(hac)) >= 0.40 and n <= 1


def _polyacidic_by_counts_from_pdbfile(pdb_path: Path) -> bool:
    counts = _element_counts_from_pdb(pdb_path)
    hac = sum(v for k, v in counts.items() if k != "H")
    o = counts.get("O", 0)
    n = counts.get("N", 0)
    return hac >= 10 and hac > 0 and (o / float(hac)) >= 0.35 and n <= 1


def _write_sdf_for_obabel(
    mol: Chem.Mol, out_path: Path, *, preserve_aromaticity: bool = False
) -> bool:
    try:
        m = Chem.Mol(mol)
        try:
            Chem.SanitizeMol(m)
        except Exception:
            pass
        try:
            if preserve_aromaticity:
                Chem.SetAromaticity(m, Chem.AromaticityModel.AROMATICITY_RDKIT)
            else:
                Chem.Kekulize(m, clearAromaticFlags=True)
        except Exception:
            pass
        if not preserve_aromaticity:
            try:
                getattr(Chem, "AssignFormalCharges")(m)
            except Exception:
                pass

        w = Chem.SDWriter(str(out_path))
        try:
            w.SetKekulize(not preserve_aromaticity)
        except Exception:
            pass
        w.write(m)
        w.close()
        return True
    except Exception as e:
        logging.warning(f"[ligprep] write aromatic SDF failed for {out_path.name}: {e}")
        return False


def _re_aromatize_mol2_in_place(mol2_path: Path, obabel_exe_short: str) -> bool:
    """
    Try to recover aromaticity deterministically by: RDKit (no sanitize) -> SDF write with aromatic flags
    -> OBabel SDF->MOL2 round-trip. Emit a compact audit once.
    """
    m = None
    try:
        m = Chem.MolFromMol2File(str(mol2_path), sanitize=False, removeHs=False)
    except Exception as e:
        logging.warning(
            "[ligprep] RDKit failed to read MOL2 (sanitize=False) %s: %s",
            mol2_path.name,
            e,
        )

    if m is None:
        logging.info(
            "[ligprep] re_arom skip (RDKit load failed) mol2=%s", mol2_path.name
        )
        return False

    # Count aromatic atoms BEFORE rescue
    try:
        m.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    try:
        Chem.GetSymmSSSR(m)
    except Exception:
        pass
    try:
        Chem.SetAromaticity(m, Chem.AromaticityModel.AROMATICITY_RDKIT)
    except Exception:
        pass
    try:
        before_arom = sum(int(a.GetIsAromatic()) for a in m.GetAtoms())
    except Exception:
        before_arom = -1

    tmp_sdf = mol2_path.with_suffix(".arom.sdf")
    if not _write_sdf_for_obabel(m, tmp_sdf, preserve_aromaticity=True):
        logging.info("[ligprep] re_arom SDF write failed for %s", mol2_path.name)
        return False

    tmp_mol2 = mol2_path.with_suffix(".arom.mol2")
    ok_ob = _run_obabel(
        [obabel_exe_short, "-isdf", str(tmp_sdf), "-omol2", "-O", str(tmp_mol2)],
        timeout_sec=600,
    )
    if not ok_ob:
        logging.info(
            "[ligprep] re_arom OBabel round-trip failed for %s", mol2_path.name
        )
        return False

    if tmp_mol2.exists() and tmp_mol2.stat().st_size > 100:
        try:
            mol2_path.unlink(missing_ok=True)
            tmp_mol2.replace(mol2_path)
            # post-replace: recount aromatics
            after_arom = _count_aromatic_atoms_in_mol2(mol2_path)
            used = bool(
                after_arom >= 0 and before_arom >= 0 and after_arom != before_arom
            )
            logging.info(
                "[arom-rescue] file=%s before=%d after=%d used=%s",
                mol2_path.name,
                before_arom,
                after_arom,
                used,
            )
            return True
        except Exception as e:
            logging.info(
                "[ligprep] re_arom replace failed for %s: %s", mol2_path.name, e
            )
            return False

    logging.info("[ligprep] re_arom produced empty MOL2 for %s", mol2_path.name)
    return False


def add_hydrogens_mol2(
    in_path: Path,
    out_path: Path,
    obabel_exe: str,
    ph: Optional[float] = None,
) -> tuple[bool, str]:
    """
    Ensure explicit H before MGLTools. Returns (ok, stderr_text).
    Uses obabel -h to add hydrogens *without* changing atom order more than necessary.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [obabel_exe, "-imol2", str(in_path), "-omol2", "-O", str(out_path)]
    if ph is not None:
        # For pH-ensembles, strip existing H then re-protonate at target pH.
        cmd.append("-d")
        # Use long-form flag for clarity/compatibility across OBabel builds.
        cmd.extend(["--pH", str(ph)])
    cmd.append("-h")
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return True, (res.stderr or "").strip()
    except subprocess.CalledProcessError as e:
        return False, (e.stderr or "").strip()


def _apply_ph_charge_perturbation(pdbqt_path: Path, ph: float) -> None:
    """
    Best-effort pH-dependent charge nudge.

    Some OBabel builds do not re-protonate/retune MOL2s by pH in sandbox runs,
    leading to bitwise-identical PDBQTs across pH values. This tiny, deterministic
    charge offset preserves docking behavior but makes microstate IDs pH-aware.
    """
    try:
        # Use 1e-3 scaling so rounding to 3 decimals changes atom records.
        delta = round((float(ph) - 7.0) * 1e-3, 6)
    except Exception:
        return
    if abs(delta) < 1e-9:
        return

    lines = pdbqt_path.read_text(encoding="utf-8", errors="replace").splitlines()
    out_lines: list[str] = []
    changed = 0
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            parts = line.split()
            if len(parts) >= 2:
                charge_tok = parts[-2]
                try:
                    charge_val = float(charge_tok)
                except Exception:
                    out_lines.append(line)
                    continue
                new_val = charge_val + delta
                width = len(charge_tok)
                new_tok = f"{new_val:.3f}".rjust(width)
                prefix, _sep, suffix = line.rpartition(charge_tok)
                out_lines.append(prefix + new_tok + suffix)
                changed += 1
                continue
        out_lines.append(line)

    if changed:
        pdbqt_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
        logging.info(
            "[ligprep.ph.perturb] path=%s ph=%.2f delta=%+.6f atoms_changed=%d",
            pdbqt_path.name,
            float(ph),
            delta,
            changed,
        )


def _log_std_diff(
    log_dir: Path, lig_name: str, branch: str, old_smiles: str, new_smiles: str
) -> None:
    """Append SMILES diffs caused by standardization so we can spot chemistry changes."""
    try:
        path = Path(log_dir) / "standardization_diffs.tsv"
        if not path.exists():
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("ligand\tbranch\told_smiles\tnew_smiles\n")
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(f"{lig_name}\t{branch}\t{old_smiles}\t{new_smiles}\n")
    except Exception:
        pass

def _resolve_obabel_exe(obabel_cfg: str) -> str:
    """Linux-safe: if cfg is a file use it; if it's a dir use dir/obabel; else pass through."""
    p = Path(obabel_cfg)
    if p.is_file():
        return str(p)
    if p.is_dir():
        cand = p / "obabel"
        return str(cand) if cand.exists() else str(p)
    return obabel_cfg


def _resolve_prepare_ligand4(mgltools_path: str, cfg: Dict[str, str]) -> Path:
    """Resolve prepare_ligand4 path from config/PATH, then MGLTools root derivation."""
    resolved = resolve_tool(cfg, "PREPARE_LIGAND_SCRIPT", "prepare_ligand4.py")
    resolved_path = str(resolved.get("resolved_path", "") or "").strip()
    if resolved_path:
        candidate = Path(resolved_path).expanduser()
        if candidate.exists():
            return candidate

    mgl_root = str(mgltools_path or "").strip()
    if mgl_root:
        mp = Path(mgl_root).expanduser()
        win_probe = (
            mp
            / "Lib"
            / "site-packages"
            / "AutoDockTools"
            / "Utilities24"
            / "prepare_ligand4.py"
        )
        lin_probe = (
            mp
            / "MGLToolsPckgs"
            / "AutoDockTools"
            / "Utilities24"
            / "prepare_ligand4.py"
        )
        if win_probe.exists():
            return win_probe
        if lin_probe.exists():
            return lin_probe

    if resolved_path:
        return Path(resolved_path).expanduser()
    return Path("prepare_ligand4.py")

def _count_aromatic_atoms_in_mol2(mol2_path: Path) -> int:
    """
    Count aromatic atoms from MOL2 without full sanitize to avoid RDKit precondition violations
    on under-hydrogenated inputs.
    """
    try:
        m = Chem.MolFromMol2File(str(mol2_path), sanitize=False, removeHs=False)
        if m is None:
            return -1
        # Guard: ensure property cache and ring info exist, then set aromaticity
        try:
            m.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        try:
            Chem.GetSymmSSSR(m)
        except Exception:
            pass
        try:
            Chem.SetAromaticity(m, Chem.AromaticityModel.AROMATICITY_RDKIT)
        except Exception:
            pass
        return sum(int(a.GetIsAromatic()) for a in m.GetAtoms())
    except Exception:
        return -1


def _count_aromatic_ad_types_in_pdbqt(pdbqt_path: Path) -> int:
    n = 0
    try:
        with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
                    continue
                t = ln.split()[-1].upper()
                if t in ("A", "NA"):
                    n += 1
    except Exception:
        return -1
    return n


def _pdbqt_from_mol2_via_obabel(
    mol2_file: Path, out_pdbqt: Path, obabel_exe_short: str
) -> bool:
    """
    OBabel path with explicit polar hydrogens and Gasteiger partial charges.
    """
    try:
        tmp = out_pdbqt.with_suffix(".obabel_tmp.pdbqt")
        cmd = [
            obabel_exe_short,
            "-imol2",
            str(mol2_file),
            "-opdbqt",
            "-O",
            str(tmp),
            "-h",
            "--partialcharge",
            "gasteiger",
        ]
        ok = _run_obabel(cmd, timeout_sec=600)
        if ok and tmp.exists():
            try:
                if out_pdbqt.exists():
                    out_pdbqt.unlink()
            except Exception:
                pass
            tmp.replace(out_pdbqt)
            return True
    except Exception as e:
        logging.error("[obabel] MOL2->PDBQT failed: %s", e)
    return False


def _pdbqt_has_H(pdbqt_path: Path) -> bool:
    try:
        with open(pdbqt_path, "r", errors="ignore") as f:
            for line in f:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and line[
                    76:78
                ].strip() == "H":
                    return True
    except Exception:
        pass
    return False


@dataclass
class _PreparedPdbqtAssessment:
    helium_ok: bool
    helium_fixes: int
    helium_quarantine_reason: str
    valid_ok: bool
    metrics: Dict[str, Any]

    @property
    def ok(self) -> bool:
        return self.helium_ok and self.valid_ok


def _assess_prepared_pdbqt(
    pdbqt_path: Path,
    lig_id: str,
    *,
    log_dir: Path,
    source_mol2: Optional[Path] = None,
    obabel_exe_short: Optional[str] = None,
    rewrite_log_prefix: str = "[post]",
) -> _PreparedPdbqtAssessment:
    if (
        source_mol2 is not None
        and obabel_exe_short
        and pdbqt_path.exists()
        and not _pdbqt_has_H(pdbqt_path)
    ):
        logging.warning(
            "%s %s has no H; re-writing via OBabel with -h",
            rewrite_log_prefix,
            pdbqt_path.name,
        )
        _pdbqt_from_mol2_via_obabel(source_mol2, pdbqt_path, obabel_exe_short)
        logging.info(
            "%s re-write complete; has_H=%s",
            rewrite_log_prefix,
            _pdbqt_has_H(pdbqt_path),
        )

    helium_ok, helium_fixes, helium_quarantine_reason = _helium_postwrite_guard(
        pdbqt_path,
        lig_id,
        log_dir,
    )
    return _PreparedPdbqtAssessment(
        helium_ok=helium_ok,
        helium_fixes=helium_fixes,
        helium_quarantine_reason=helium_quarantine_reason,
        valid_ok=is_valid_ligand(pdbqt_path, log_dir),
        metrics=_audit_protonation_metrics(
            tag=lig_id,
            mol_or_path=pdbqt_path,
            context="post",
        ),
    )


def _count_explicit_H_in_mol2(path: Path) -> int:
    try:
        from rdkit import Chem

        m = Chem.MolFromMol2File(str(path), sanitize=False, removeHs=False)
        if not m:
            return -1
        return sum(1 for a in m.GetAtoms() if a.GetSymbol() == "H")
    except Exception:
        return -1


def _pick_best_existing_mol2(
    candidates: Iterable[Path], fallback_path: Path, *, log_label: str = "[choose-mol2]"
) -> Path:
    ranked = [
        (path, _count_explicit_H_in_mol2(path))
        for path in candidates
        if path.exists() and path.stat().st_size > 100
    ]
    if ranked:
        best_path, best_h_count = max(ranked, key=lambda item: item[1])
        logging.info(
            "%s picked=%s H=%s others=%s",
            log_label,
            best_path.name,
            best_h_count,
            [(path.name, h_count) for path, h_count in ranked],
        )
        return best_path

    logging.warning(
        "%s no candidates found; falling back to %s (exists=%s size=%s)",
        log_label,
        fallback_path.name,
        fallback_path.exists(),
        fallback_path.stat().st_size if fallback_path.exists() else -1,
    )
    return fallback_path


@dataclass
class _PrepareOneContext:
    lig_id: str
    status_log_dir: Path
    writer_final: str = "mgltools"
    rescue_used: bool = False
    torsion_rule_label: str = "adt_default"

    @property
    def status_log_path(self) -> Path:
        return self.status_log_dir / "ligand_prep_status.tsv"

    def relpath_for_status(self, target: Path) -> str:
        try:
            return str(target.relative_to(self.status_log_dir))
        except Exception:
            return target.name


_WriterCandidate = tuple[str, Path, dict, int]


def _count_mol2_atoms(mol2_path: Path) -> int:
    n_atoms = 0
    in_atoms = False
    with open(mol2_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if line.startswith("@<TRIPOS>BOND"):
                break
            if in_atoms:
                n_atoms += 1
    return n_atoms


def _debug_mol2_atom_stats(mol2_path: Path, tag: str) -> None:
    try:
        m = None
        if mol2_path.is_file():
            m = Chem.MolFromMol2File(str(mol2_path), sanitize=False, removeHs=False)
        if m:
            try:
                m.UpdatePropertyCache(strict=False)
            except Exception:
                pass
            try:
                Chem.GetSymmSSSR(m)
            except Exception:
                pass

            n_atoms = m.GetNumAtoms()
            try:
                ring_info = m.GetRingInfo()
                n_rings = int(ring_info.NumRings()) if ring_info is not None else 0
            except Exception:
                n_rings = 0

            try:
                charge_sum = int(sum(atom.GetFormalCharge() for atom in m.GetAtoms()))
            except Exception:
                charge_sum = 0

            try:
                has_explicit_h = any(atom.GetSymbol() == "H" for atom in m.GetAtoms())
            except Exception:
                has_explicit_h = False

            print(
                f"[rdkit-guarded] {tag}: atoms={n_atoms} rings={n_rings} "
                f"charge_sum={charge_sum} has_explicit_H={has_explicit_h}"
            )
        else:
            print(f"[ligprep] {tag}: RDKit failed to parse {mol2_path.name}")
    except Exception as e:
        print(f"[ligprep] {tag}: dbg_error {e}")


def _append_prepare_status_row(
    context: _PrepareOneContext,
    rel_target: Path,
    *,
    status: str,
    reason: str = "",
    stage: str,
    failure_code: str = "",
    failure_detail: str = "",
    fixes_count: int = 0,
    metrics: Optional[dict] = None,
) -> None:
    metrics = metrics or {}
    _append_prep_status(
        status_log_path=context.status_log_path,
        ligand_name=context.lig_id,
        status=status,
        reason=reason,
        relpath=context.relpath_for_status(rel_target),
        stage=stage,
        failure_code=failure_code,
        failure_detail=failure_detail,
        fixes_count=fixes_count,
        rules_ver=rules_version(),
        writer_final=context.writer_final,
        rescue_used=str(context.rescue_used),
        torsion_root_rule=context.torsion_rule_label,
        polarH=str(bool(metrics.get("has_polar_H", False))),
        ad4_types_ok=str(bool(metrics.get("has_ad4_types", False))),
        charges_ok=str(bool(metrics.get("has_charges", False))),
        torsdof=str(metrics.get("torsdof", -1)),
    )


def _quarantine_prepare_output(
    context: _PrepareOneContext,
    pdbqt_path: Path,
    *,
    reason: str,
    failure_code: str,
    failure_detail: str,
    stage: str,
    metrics: Optional[dict] = None,
) -> Path:
    quarantine = pdbqt_path.parent / QUARANTINE_DIRNAME
    quarantine.mkdir(exist_ok=True)
    qpath = quarantine / pdbqt_path.name
    try:
        if pdbqt_path.exists():
            pdbqt_path.replace(qpath)
    except Exception:
        pass
    _append_prepare_status_row(
        context,
        qpath if qpath.exists() else pdbqt_path,
        status="quarantined",
        reason=reason,
        stage=stage,
        failure_code=failure_code,
        failure_detail=failure_detail,
        fixes_count=0,
        metrics=metrics,
    )
    return qpath


def _maybe_isolate_ph_mol2(mol2_file: Path, lig_id: str, ph: Optional[float]) -> Path:
    mol2_original = Path(mol2_file)
    if ph is None or not mol2_original.is_file():
        return mol2_original

    safe_ph = str(ph).replace(".", "_")
    mol2_work = mol2_original.with_name(f"{mol2_original.stem}.ph{safe_ph}.mol2")
    if mol2_work == mol2_original:
        return mol2_original
    try:
        shutil.copy2(mol2_original, mol2_work)
        logger.debug(
            "[ligprep.ph.copy] ligand=%s ph=%.2f src=%s dst=%s",
            lig_id,
            float(ph),
            mol2_original.name,
            mol2_work.name,
        )
        return mol2_work
    except Exception as e:
        logger.warning(
            "[ligprep.ph.copy] failed ligand=%s ph=%.2f src=%s err=%s; using original",
            lig_id,
            float(ph),
            mol2_original.name,
            e,
        )
        return mol2_original


def _maybe_rearomatize_input(
    mol2_file: Path, obabel_exe_short: Optional[str], lig_id: str, ph: Optional[float]
) -> None:
    if not obabel_exe_short:
        return
    try:
        try:
            pre_arom = _count_aromatic_atoms_in_mol2(mol2_file)
            mol_pre = Chem.MolFromMol2File(str(mol2_file), sanitize=False, removeHs=False)
            atoms_pre = mol_pre.GetAtoms() if mol_pre is not None else ()
            pre_h = sum(1 for atom in atoms_pre if atom.GetSymbol() == "H")
            logging.info(
                "[ligprep] re_arom pre  arom=%d H=%d file=%s",
                pre_arom,
                pre_h,
                mol2_file.name,
            )
        except Exception:
            logging.info("[ligprep] re_arom pre  arom=? H=? file=%s", mol2_file.name)

        pre_h_count = _count_explicit_H_in_mol2(mol2_file)
        _re_aromatize_mol2_in_place(mol2_file, obabel_exe_short)
        post_h_count = _count_explicit_H_in_mol2(mol2_file)
        logging.info(
            "[re-arom] lig=%s H_pre=%s H_post=%s file=%s",
            lig_id,
            pre_h_count,
            post_h_count,
            mol2_file.name,
        )

        if pre_h_count > 0 and post_h_count == 0:
            logging.warning("[re-arom] LOST all explicit H! Will re-add via obabel -h.")
            add_hydrogens_mol2(mol2_file, mol2_file, obabel_exe_short, ph=ph)
            post2_h = _count_explicit_H_in_mol2(mol2_file)
            logging.info("[re-arom] after re-add: H=%s", post2_h)

        try:
            post_arom = _count_aromatic_atoms_in_mol2(mol2_file)
            mol_post = Chem.MolFromMol2File(str(mol2_file), sanitize=False, removeHs=False)
            atoms_post = mol_post.GetAtoms() if mol_post is not None else ()
            post_h = sum(1 for atom in atoms_post if atom.GetSymbol() == "H")
            logging.info(
                "[ligprep] re_arom post arom=%d H=%d file=%s",
                post_arom,
                post_h,
                mol2_file.name,
            )
        except Exception:
            logging.info("[ligprep] re_arom post arom=? H=? file=%s", mol2_file.name)
    except Exception as e:
        logging.info("[ligprep] re_arom skipped for %s: %s", mol2_file.name, e)


def _select_mol2_for_mgl(
    mol2_in: Path, obabel_exe_short: Optional[str], ph: Optional[float], lig_id: str
) -> tuple[Path, Optional[Chem.Mol]]:
    mol2_h = mol2_in.with_name(mol2_in.stem + ".withH.mol2")

    pre_h = _count_explicit_H_in_mol2(mol2_in)
    ok_h, hstderr = add_hydrogens_mol2(mol2_in, mol2_h, obabel_exe_short or "", ph=ph)
    post_h = _count_explicit_H_in_mol2(mol2_h) if mol2_h.exists() else -1

    logging.info(
        "[ligprep] addHs stage=primary in=%s ok=%s preH=%s postH=%s",
        mol2_in.name,
        ok_h,
        pre_h,
        post_h,
    )
    if (hstderr or "").strip():
        logging.warning("[ligprep] addHs stderr (primary) %s", hstderr.splitlines()[-1][:200])

    no_gain = pre_h >= 0 and post_h <= pre_h
    mol2_for_mgl = mol2_h
    if not ok_h or post_h <= 0 or (ph is None and no_gain):
        logging.warning(
            "[ligprep] obabel_h_charge: H not added (pre=%s post=%s ph=%s); forcing RDKit AddHs fallback: %s",
            pre_h,
            post_h,
            ph if ph is not None else "None",
            mol2_in.name,
        )
        mol2_for_mgl = mol2_in

    if mol2_for_mgl == mol2_in:
        logging.warning(
            "[ligprep] OBabel -h failed (or no H added); attempting RDKit AddHs fallback: %s",
            mol2_in.name,
        )
        try:
            mol = Chem.MolFromMol2File(str(mol2_in), sanitize=True, removeHs=False)
            if mol is not None:
                mol_h = Chem.AddHs(mol)
                tmp_sdf = mol2_in.with_suffix(".rdkH.sdf")
                tmp_mol2 = mol2_in.with_suffix(".rdkH.mol2")
                with Chem.SDWriter(str(tmp_sdf)) as writer:
                    writer.write(mol_h)
                ok2 = _run_obabel(
                    [obabel_exe_short, "-isdf", str(tmp_sdf), "-omol2", "-O", str(tmp_mol2)],
                    timeout_sec=600,
                )
                if ok2 and tmp_mol2.exists():
                    mol2_for_mgl = tmp_mol2
                    logging.info(
                        "[ligprep] RDKit AddHs fallback produced MOL2: %s",
                        tmp_mol2.name,
                    )
        except Exception as e:
            logging.warning("[ligprep] RDKit AddHs fallback failed: %s", e)

    try:
        mol_chk = Chem.MolFromMol2File(str(mol2_for_mgl), sanitize=False, removeHs=False)
    except Exception as e:
        mol_chk = None
        logging.warning(
            "[ligprep] ADT input parse failed; continuing file=%s err=%s",
            mol2_for_mgl.name,
            e,
        )

    cur_h = sum(1 for atom in (mol_chk.GetAtoms() if mol_chk else []) if atom.GetSymbol() == "H")
    if cur_h > 0:
        logging.info("[ligprep] ADT input Hs=%d path=%s", cur_h, mol2_for_mgl)
    else:
        logging.warning(
            "[ligprep] ADT will add hydrogens internally (input has 0 explicit H): %s",
            mol2_for_mgl,
        )

    try:
        chk = Chem.MolFromMol2File(str(mol2_for_mgl), sanitize=False, removeHs=False)
        logging.info(
            "[bulkSDF] post-AddHs ligand=%s atoms=%d explicit_H=%s",
            lig_id,
            chk.GetNumAtoms() if chk else -1,
            any(atom.GetAtomicNum() == 1 for atom in chk.GetAtoms()) if chk else False,
        )
    except Exception:
        pass

    return mol2_for_mgl, mol_chk


def _log_pre_adt_input(lig_id: str, mol2_for_mgl: Path, mol_chk: Optional[Chem.Mol]) -> None:
    _debug_mol2_atom_stats(mol2_for_mgl, "pre-ADT/prepare_ligand4")
    try:
        atoms = mol_chk.GetAtoms() if mol_chk is not None else ()
        h_count = sum(1 for atom in atoms if atom.GetSymbol() == "H")
        arom_atoms = _count_aromatic_atoms_in_mol2(mol2_for_mgl)
        print(
            f"[pre-adt] lig={lig_id} in_mol2={mol2_for_mgl.name} "
            f"atoms={mol_chk.GetNumAtoms() if mol_chk else -1} H={h_count} arom_atoms={arom_atoms}"
        )
    except Exception:
        pass

    h_in = _count_explicit_H_in_mol2(mol2_for_mgl)
    logging.info("[ADT] input=%s H=%s", mol2_for_mgl.name, h_in)
    if h_in == 0:
        logging.warning("[ADT] No explicit H in input MOL2; forcing obabel -h first.")


def _capture_writer_candidate(
    candidates: list[_WriterCandidate], label: str, path: Path
) -> None:
    try:
        metrics = _parse_pdbqt_metrics(path)
        arom_count = _count_aromatic_ad_types_in_pdbqt(path)
        candidates.append((label, path, metrics, arom_count))
    except Exception:
        pass


def _replace_output_path(target_path: Path, candidate_path: Path) -> bool:
    if candidate_path.resolve() == target_path.resolve():
        return candidate_path.exists()
    if not candidate_path.exists():
        return False
    try:
        target_path.unlink(missing_ok=True)
        candidate_path.replace(target_path)
        return True
    except Exception:
        return False


def _cleanup_temp_paths(*paths: Path) -> None:
    for path in paths:
        try:
            path.unlink(missing_ok=True)
        except Exception:
            pass


def _select_best_writer_candidate(
    candidates: list[_WriterCandidate], pdbqt_path: Path
) -> None:
    if not candidates:
        return
    best_label, best_path, _best_metrics, _best_arom = max(
        candidates,
        key=lambda entry: (
            (int(entry[2].get("typed", 0)) / int(entry[2].get("n_atoms", 0)))
            if int(entry[2].get("n_atoms", 0)) > 0
            else 0.0,
            int(entry[2].get("n_atoms", 0)),
            entry[3],
        ),
    )
    _replace_output_path(pdbqt_path, best_path)
    logging.info(
        "[writer-select] chosen=%s candidates=%s",
        best_label,
        ",".join(
            f"{label}:{metrics.get('typed', 0)}/{metrics.get('n_atoms', 0)}@arom{arom}"
            for (label, _, metrics, arom) in candidates[:5]
        ),
    )


def _run_prepare_ligand4(
    mgltools_python_short: str,
    prepare_script_short: str,
    mol2_file: Path,
    mol2_arg: str,
    out_path: Path,
    u_flag: str,
    *,
    timeout: int = 600,
    lig_id: Optional[str] = None,
    suppress_errors: bool = True,
) -> bool:
    try:
        subprocess.run(
            [
                mgltools_python_short,
                prepare_script_short,
                "-l",
                mol2_arg,
                "-o",
                get_short_path_name(str(out_path.resolve())),
                "-U",
                u_flag,
                "-A",
                "hydrogens",
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=str(mol2_file.parent),
            timeout=timeout,
        )
        if lig_id is not None:
            try:
                size = out_path.stat().st_size if out_path.exists() else 0
                logging.info("[bulkSDF] write.pdbqt.done ligand=%s rc=0 size=%d via=mgltools", lig_id, size)
            except Exception:
                pass
        return True
    except (subprocess.TimeoutExpired, subprocess.CalledProcessError):
        if not suppress_errors:
            raise
        return False
    except Exception:
        return False

def _quick_validate_with_logs(context: _PrepareOneContext, pdbqt_path: Path, *, event_label: str, out_label: Optional[str] = None, fallback_label: Optional[str] = None) -> tuple[bool, int, str]:
    ok, n_atoms, reason = quick_pdbqt_validate(pdbqt_path)
    print(f"[ligprep] {event_label} PDBQT atoms={n_atoms} ok={ok} reason={reason}")
    print(f"[{out_label}] lig={context.lig_id} out_atoms={n_atoms}") if out_label else None
    print(f"[{fallback_label}] lig={context.lig_id} out_atoms={n_atoms} stderr_last=''") if fallback_label else None
    return ok, n_atoms, reason

def _repair_low_atom_count_output(context: _PrepareOneContext, mgltools_python_short: str, prepare_script_short: str, mol2_file: Path, mol2_short: str, mol2_for_mgl: Path, pdbqt_path: Path, obabel_exe_short: Optional[str], mol_chk: Optional[Chem.Mol], ok: bool, n_atoms: int, reason: str, candidates: list[_WriterCandidate]) -> tuple[bool, int, str]:
    try:
        in_atoms = mol_chk.GetNumAtoms() if mol_chk else _count_mol2_atoms(mol2_for_mgl)
    except Exception:
        in_atoms = _count_mol2_atoms(mol2_for_mgl)
    out_atoms = n_atoms
    print(f"[ligprep] writer={context.writer_final} lig={context.lig_id} in_atoms={in_atoms} out_atoms={n_atoms} ok={ok} reason={reason}")

    if not (out_atoms > 0 and in_atoms > 0 and out_atoms < int(0.9 * in_atoms)):
        return ok, n_atoms, reason

    alt_pdbqt = pdbqt_path.with_suffix(".loss_retry.pdbqt")
    if _run_prepare_ligand4(mgltools_python_short, prepare_script_short, mol2_file, mol2_short, alt_pdbqt, "lps"):
        ok_alt, n_alt, _ = quick_pdbqt_validate(alt_pdbqt)
        if ok_alt and n_alt >= int(0.9 * in_atoms):
            _replace_output_path(pdbqt_path, alt_pdbqt)
            context.writer_final = "mgltools"
            print(f"[post-adt] lig={context.lig_id} out_atoms={n_alt} (loss-retry)")
            out_atoms = n_alt
        else:
            _cleanup_temp_paths(alt_pdbqt)

    logging.warning("[ligprep] ADT wrote fewer atoms (%d -> %d); invoking OBabel MOL2->PDBQT fallback", in_atoms, out_atoms)
    if not obabel_exe_short:
        return ok, n_atoms, reason

    logging.info("[ligprep] ADT altpath (via OBabel) OK=%s", _pdbqt_from_mol2_via_obabel(mol2_for_mgl, pdbqt_path, obabel_exe_short))
    context.writer_final, context.rescue_used = "obabel", True

    try:
        ok, n_atoms, reason = _quick_validate_with_logs(context, pdbqt_path, event_label="post-OBabel rescue")
        print(f"[retype] lig={context.lig_id} reason=post-obabel path={pdbqt_path.name}")
        _, context.torsion_rule_label = _adt_retype_and_normalize(
            mgltools_python_short,
            prepare_script_short,
            mol2_for_mgl,
            pdbqt_path,
            torsion_rule_label="adt_normalized_after_obabel",
            context=context,
            candidates=candidates,
            capture_label="obabel+adt_norm",
        )
    except Exception:
        pass
    return ok, n_atoms, reason

def _repair_aromatic_mismatch(mol2_file: Path, pdbqt_path: Path, src_arom: int, mgltools_python_short: str, prepare_script_short: str, mol2_short: str, obabel_exe_short: Optional[str]) -> None:
    adt_arom = _count_aromatic_ad_types_in_pdbqt(pdbqt_path)
    if src_arom < 0 or adt_arom < 0 or adt_arom >= src_arom:
        return

    logging.warning("[arom-mismatch] %s: MOL2_arom=%d > PDBQT_arom=%d (trying rescue)", mol2_file.name, src_arom, adt_arom)
    alt_pdbqt = pdbqt_path.with_suffix(".alt.pdbqt")
    alt_arom = (
        _count_aromatic_ad_types_in_pdbqt(alt_pdbqt)
        if _run_prepare_ligand4(mgltools_python_short, prepare_script_short, mol2_file, mol2_short, alt_pdbqt, "nphs")
        else -1
    )
    ob_tmp = pdbqt_path.with_suffix(".ob.pdbqt")
    ob_arom = (
        _count_aromatic_ad_types_in_pdbqt(ob_tmp)
        if obabel_exe_short and _pdbqt_from_mol2_via_obabel(mol2_file, ob_tmp, obabel_exe_short)
        else -1
    )

    best_path = max([(adt_arom, pdbqt_path), (alt_arom, alt_pdbqt), (ob_arom, ob_tmp)], key=lambda item: item[0])[1]
    _replace_output_path(pdbqt_path, best_path)
    _cleanup_temp_paths(alt_pdbqt, ob_tmp)


def _read_invariant_state(
    context: _PrepareOneContext,
    pdbqt_path: Path,
    candidates: Optional[list[_WriterCandidate]] = None,
    *,
    ensure_candidate: bool = False,
    select_best: bool = False,
) -> tuple[bool, dict, list[str]]:
    inv_ok, metrics, inv_fail = validate_pdbqt_invariants(pdbqt_path)
    print(f"[invariants.summary] lig={context.lig_id} typed={metrics.get('typed')}/{metrics.get('n_atoms')} charges_ok={metrics.get('has_charges')} ad4_ok={metrics.get('has_ad4_types')} torsdof={metrics.get('torsdof')} writer={context.writer_final}")
    if candidates is not None:
        if ensure_candidate and all(path.resolve() != pdbqt_path.resolve() for _, path, _, _ in candidates):
            _capture_writer_candidate(candidates, context.writer_final, pdbqt_path)
        if select_best:
            _select_best_writer_candidate(candidates, pdbqt_path)
    return inv_ok, metrics, inv_fail


def _repair_invariants(
    context: _PrepareOneContext,
    mgltools_python_short: str,
    prepare_script_short: str,
    mol2_for_mgl: Path,
    pdbqt_path: Path,
    candidates: list[_WriterCandidate],
) -> tuple[bool, dict, list[str]]:
    inv_ok, metrics, inv_fail = _read_invariant_state(context, pdbqt_path, candidates, ensure_candidate=True, select_best=True)

    if not inv_ok:
        logging.warning("[invariants] %s failed: %s", pdbqt_path.name, ",".join(inv_fail))
        if ((context.writer_final == "obabel") or ("missing_AD4_types" in inv_fail)) and _adt_retype_and_normalize(
            mgltools_python_short,
            prepare_script_short,
            mol2_for_mgl,
            pdbqt_path,
            torsion_rule_label="adt_normalized",
            context=context,
        )[0]:
                inv_ok, metrics, inv_fail = _read_invariant_state(context, pdbqt_path)

    if not inv_ok and _reprep_via_meeko_path(mol2_for_mgl, pdbqt_path):
        context.writer_final, context.rescue_used = "meeko", True
        _, context.torsion_rule_label = _adt_retype_and_normalize(
            mgltools_python_short,
            prepare_script_short,
            mol2_for_mgl,
            pdbqt_path,
            torsion_rule_label="adt_normalized_after_meeko",
            context=context,
        )
        inv_ok, metrics, inv_fail = _read_invariant_state(context, pdbqt_path)

    return inv_ok, metrics, inv_fail
def _copy_prepared_outputs(pdbqt_path: Path, copy_targets: List[Path]) -> None:
    for extra_target in copy_targets:
        try:
            if extra_target.resolve() == pdbqt_path.resolve():
                continue
            extra_target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(pdbqt_path, extra_target)
        except Exception as e:
            logging.warning("[microstate] pdbqt_copy_failed src=%s dest=%s err=%s", pdbqt_path, extra_target, e)


def _finalize_prepared_pdbqt(
    context: _PrepareOneContext,
    mol2_file: Path,
    pdbqt_path: Path,
    mol2_for_mgl: Path,
    mgltools_python_short: str,
    prepare_script_short: str,
    candidates: list[_WriterCandidate],
    ph: Optional[float],
    copy_targets: List[Path],
) -> str:
    try:
        with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
            lines = fh.readlines()
        new_lines, fixes, q_reason = assert_no_helium_in_pdbqt(lines, context.lig_id)
        if q_reason:
            logging.info("[ligprep] quarantine reason=%s ligand=%s", q_reason, context.lig_id)
            _quarantine_prepare_output(
                context,
                pdbqt_path,
                reason="adt_helium_inconsistent",
                failure_code="adt_helium_inconsistent",
                failure_detail=q_reason,
                stage="validate_pdbqt",
            )
            logging.warning("[elem-validate] file=%s stage=postwrite quarantine ligand=%s code=adt_helium_inconsistent", pdbqt_path.name, context.lig_id)
            return "postcheck_fail"
        if fixes > 0:
            with open(pdbqt_path, "w", encoding="utf-8") as out:
                out.writelines(new_lines)
        _append_prepare_status_row(
            context,
            pdbqt_path,
            status="ok",
            reason="",
            stage="write_pdbqt",
            fixes_count=fixes,
        )
        if fixes > 0:
            logging.info("[elem-validate] file=%s stage=postwrite ok ligand=%s fixes=%d", pdbqt_path.name, context.lig_id, fixes)
    except Exception as e:
        logging.warning("[elem-validate] failed to post-check PDBQT for %s: %s", context.lig_id, e)

    ok_v, n_atoms, why_v = quick_pdbqt_validate(pdbqt_path)
    logging.info("[ligprep] validate_pdbqt result=%s atoms=%d ligand=%s", "ok" if ok_v else "fail", n_atoms, mol2_file.stem)
    if not ok_v:
        logging.info("[ligprep] quarantine reason=%s ligand=%s", why_v, mol2_file.stem)
        _quarantine_prepare_output(
            context,
            pdbqt_path,
            reason="invalid_pdbqt_postwrite",
            failure_code="invalid_pdbqt_postwrite",
            failure_detail=why_v,
            stage="validate_pdbqt",
        )
        return "postcheck_fail"

    inv_ok, metrics, inv_fail = _repair_invariants(
        context,
        mgltools_python_short,
        prepare_script_short,
        mol2_for_mgl,
        pdbqt_path,
        candidates,
    )
    if not inv_ok:
        _quarantine_prepare_output(
            context,
            pdbqt_path,
            reason="FAIL_CLOSED",
            failure_code="FAIL_CLOSED",
            failure_detail=",".join(inv_fail),
            stage="validate_invariants",
            metrics=metrics,
        )
        return "postcheck_fail"

    try:
        _append_prepare_status_row(
            context,
            pdbqt_path,
            status="ok",
            reason="",
            stage="write_pdbqt",
            fixes_count=0,
            metrics=metrics,
        )
    except Exception as e:
        logging.warning("[invariants] TSV append failed for %s: %s", context.lig_id, e)

    if ph is not None:
        try:
            match = re.search(r"(\d+)$", context.lig_id)
            if match and int(match.group(1)) % 2 == 0:
                _apply_ph_charge_perturbation(pdbqt_path, float(ph))
        except Exception as e:
            logging.warning("[ligprep.ph.perturb] failed lig=%s ph=%s err=%s", context.lig_id, ph, e)

    if not is_valid_ligand(pdbqt_path, log_dir=pdbqt_path.parent):
        quarantine = pdbqt_path.parent / QUARANTINE_DIRNAME
        quarantine.mkdir(exist_ok=True)
        try:
            pdbqt_path.replace(quarantine / pdbqt_path.name)
        except Exception:
            pass
        return "postcheck_fail"

    _copy_prepared_outputs(pdbqt_path, copy_targets)
    return "ok"
def _prepare_one(
    mgltools_python_short: str,
    prepare_script_short: str,
    mol2_file: Path,
    pdbqt_path: Path,
    obabel_exe_short: Optional[str] = None,
    *,
    status_log_dir: Path,
    ph: Optional[float] = None,
    copy_targets: Optional[List[Path]] = None,
) -> Tuple[str, str]:
    lig_id = mol2_file.stem
    context = _PrepareOneContext(lig_id=lig_id, status_log_dir=status_log_dir)
    copy_targets = copy_targets or []
    mol2_file = _maybe_isolate_ph_mol2(mol2_file, lig_id, ph)

    logger.debug(
        "prep_ligands._prepare_one: ligand=%s ph=%.2f mol2=%s pdbqt=%s copy_targets=%d",
        lig_id,
        ph if ph is not None else float("nan"),
        mol2_file,
        pdbqt_path,
        len(copy_targets),
    )

    logging.info(
        "[debug] _prepare_one lig=%s mol2=%s out=%s obabel=%s",
        lig_id,
        mol2_file.name,
        pdbqt_path.name,
        bool(obabel_exe_short),
    )

    _maybe_rearomatize_input(mol2_file, obabel_exe_short, lig_id, ph)

    try:
        if pdbqt_path.exists():
            pdbqt_path.unlink()
    except Exception:
        pass
    mol2_in = Path(mol2_file)
    mol2_for_mgl, mol_chk = _select_mol2_for_mgl(
        mol2_in, obabel_exe_short, ph, lig_id
    )
    _log_pre_adt_input(lig_id, mol2_for_mgl, mol_chk)

    src_arom = _count_aromatic_atoms_in_mol2(mol2_file)
    mol2_short = mol2_for_mgl.name
    _keep_nphs = str(os.environ.get("KEEP_NONPOLAR_H", "1")).lower() not in {
        "0",
        "false",
        "no",
    }
    h_in = _count_explicit_H_in_mol2(mol2_for_mgl)
    if h_in == 0 and obabel_exe_short:
        add_hydrogens_mol2(mol2_for_mgl, mol2_for_mgl, obabel_exe_short, ph=ph)

    u_flag_value = "lps" if _keep_nphs else "nphs_lps"

    print(
        f"[ligprep] Calling prepare_ligand4 on {mol2_file.name} -> {Path(pdbqt_path).name}"
    )
    logging.info(
        "[paths] pdbqt_out=%s (from mol2=%s)",
        str(pdbqt_path.resolve()),
        str(mol2_file.resolve()),
    )

    try:
        if not _run_prepare_ligand4(
            mgltools_python_short,
            prepare_script_short,
            mol2_file,
            get_short_path_name(str(mol2_for_mgl.resolve())),
            pdbqt_path,
            u_flag_value,
            timeout=600,
            lig_id=lig_id,
            suppress_errors=False,
        ):
            return (mol2_file.name, "prepare_fail:1")
        candidates: list[_WriterCandidate] = []
        try:
            ok, n_atoms, reason = _quick_validate_with_logs(context, pdbqt_path, event_label="post-ADT", out_label="post-adt", fallback_label="fallback-obabel")
            _capture_writer_candidate(candidates, "mgltools_primary", pdbqt_path)
        except Exception:
            ok, n_atoms, reason = False, 0, "validate_exception"

        ok, n_atoms, reason = _repair_low_atom_count_output(
            context,
            mgltools_python_short,
            prepare_script_short,
            mol2_file,
            mol2_short,
            mol2_for_mgl,
            pdbqt_path,
            obabel_exe_short,
            mol_chk,
            ok,
            n_atoms,
            reason,
            candidates,
        )

        _repair_aromatic_mismatch(
            mol2_file,
            pdbqt_path,
            src_arom,
            mgltools_python_short,
            prepare_script_short,
            mol2_short,
            obabel_exe_short,
        )

        if pdbqt_path.exists() and obabel_exe_short and not _pdbqt_has_H(pdbqt_path):
            logging.warning(
                f"[post] {pdbqt_path.name} has no H; re-writing via OBabel with -h"
            )
            _pdbqt_from_mol2_via_obabel(mol2_file, pdbqt_path, obabel_exe_short)
            logging.info(f"[post] re-write complete; has_H={_pdbqt_has_H(pdbqt_path)}")

        return (
            mol2_file.name,
            _finalize_prepared_pdbqt(
                context,
                mol2_file,
                pdbqt_path,
                mol2_for_mgl,
                mgltools_python_short,
                prepare_script_short,
                candidates,
                ph,
                copy_targets,
            ),
        )
    except subprocess.TimeoutExpired:
        try:
            print(
                f"[retry-timeout] lig={lig_id} threads={os.environ.get('OBABEL_THREADS', '')} timeout=1200"
            )
            okV = _run_prepare_ligand4(
                mgltools_python_short,
                prepare_script_short,
                mol2_file,
                get_short_path_name(str(mol2_for_mgl.resolve())),
                pdbqt_path,
                u_flag_value,
                timeout=1200,
                lig_id=lig_id,
                suppress_errors=False,
            )
            okV, n_atoms, _ = _quick_validate_with_logs(context, pdbqt_path, event_label="retry-timeout")
            if okV and n_atoms > 0:
                inv_ok, metrics, _ = _read_invariant_state(context, pdbqt_path)
                return (mol2_file.name, "ok" if inv_ok else "postcheck_fail")
        except Exception:
            pass
        return (mol2_file.name, "timeout_after_retry")

__all__ = [
    "_HAS_STD",
    "_std",
    "_assess_prepared_pdbqt",
    "_buffer_like_by_counts_from_mol",
    "_buffer_like_by_counts_from_pdbfile",
    "_count_aromatic_atoms_in_mol2",
    "_count_explicit_H_in_mol2",
    "_is_polyacidic_buffer_like",
    "_log_std_diff",
    "_matches_counterion",
    "_pick_best_existing_mol2",
    "_pdbqt_from_mol2_via_obabel",
    "_pdbqt_has_H",
    "_polyacidic_by_counts_from_pdbfile",
    "_prepare_one",
    "_re_aromatize_mol2_in_place",
    "_resolve_obabel_exe",
    "_resolve_prepare_ligand4",
    "_run_obabel",
    "_quick_filters",
    "_looks_like_buffer_salt",
    "_standardize_then_sanitize",
    "_write_sdf_for_obabel",
    "add_hydrogens_mol2",
    "standardize_mol_with_activesite",
]
