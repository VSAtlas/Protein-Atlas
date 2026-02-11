import os
from pathlib import Path
import subprocess
import logging
import concurrent.futures
import importlib
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Optional, Tuple, Dict, Set, Any
import re
from datetime import datetime
import shutil
from collections import Counter

from config.tool_resolver import resolve_tool
from input_and_export_functions import load_config, validate_config
from path_router import make_paths
from prep_ligands.prep_ligands_crystal import prep_ligands_from_pdb

# --- RDKit / Standardization imports ---
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import Crippen
from rdkit.Chem import rdmolops as rdMolOps  # new: shared rdMolOps import

from prep_ligands.prep_ligands_bulk_sdf import (
    LIGPREP_OBABEL_TIMEOUT_SEC,
    _sdf_to_mol2,
    convert_sdf_to_mol2_split_parallel,
    rename_mol2_with_prefix,
)
from prep_ligands.prep_ligands_microstates import (
    _collect_only_from_env_and_cli,
    enumerate_ligands_for_docking,
    compute_microstate_id,
    load_microstate_registry,
    save_microstate_registry,
)
import prep_ligands.prep_ligands_common as ligprep_common
from prep_ligands.prep_ligands_common import (
    EXCLUDE_CRYSTAL_ADDITIVES,
    MIN_ATOMS_FOR_DOCKING,
    QUARANTINE_DIRNAME,
    STANDARD_AMINO_ACIDS,
    _HAS_STD,
    _append_prep_status,
    _audit_protonation_metrics,
    _buffer_like_by_counts_from_mol,
    _buffer_like_by_counts_from_pdbfile,
    _count_aromatic_atoms_in_mol2,
    _count_explicit_H_in_mol2,
    _is_polyacidic_buffer_like,
    _log_elem_fix_summary,
    _log_malformed,
    _log_std_diff,
    _matches_counterion,
    _pdbqt_from_mol2_via_obabel,
    _pdbqt_has_H,
    _polyacidic_by_counts_from_pdbfile,
    _prepare_one,
    _re_aromatize_mol2_in_place,
    collapse_sanitized_path,
    _resolve_obabel_exe,
    _resolve_prepare_ligand4,
    _run_obabel,
    _std,
    _looks_like_buffer_salt,
    _write_aromatic_sdf,
    get_short_path_name,
    is_valid_ligand,
    read_config,
    standardize_mol_with_activesite,
    _quick_filters,
    _standardize_then_sanitize,
)

logger = logging.getLogger(__name__)

try:
    from rdkit.Chem.MolStandardize import rdMolStandardize as _std  # unified handle

    _HAS_STD = True
except Exception:
    _std = None
    _HAS_STD = False


LIGPREP_PH = 7.4
KEEP_NONPOLAR_H = 1

# --- Salvage / logging config ---
RUN_TAG = datetime.now().strftime("%Y%m%d_%H%M%S")

# --- Resume mode toggle ---
RESUME_SKIP = False


MIN_PARENT_HEAVY = 4

# ----------- tuning switches -----------
USE_RDKIT_FOR_3D = True
OBABEL_THREADS = 50
CHUNK_SIZE = 200


# --------------------------------------

# --- salt remover (fallback path also uses this) ---
_REM = SaltRemover()
# =========================
# Alias/rules centralization
# =========================

# Effective, alias-derived sets/maps populated at import:
RULES = None
ALLOWED_ELEMENTS: Set[str] = set()
MONOATOMIC_IONS: Set[str] = set()
# AD4_TYPES already exists later in the file; we will reconcile it here.


# =========================
# Utility / logging helpers
# =========================

# --- PDBQT metrics + invariants ---------------------------------------------


def _elem_hist(m):
    try:
        return dict(Counter(a.GetSymbol() for a in m.GetAtoms()))
    except Exception:
        return {}


def _formal_charge(m):
    try:
        return sum(a.GetFormalCharge() for a in m.GetAtoms())
    except Exception:
        return 0


def _has_explicit_Hs(m):
    try:
        return any(a.GetAtomicNum() == 1 for a in m.GetAtoms())
    except Exception:
        return False


def _dbg_atom_stats(m: Chem.Mol, tag: str) -> None:
    """
    Lightweight, fully-guarded debug helper for RDKit molecules.
    Prints element histogram, ring count, charge sum, and explicit-H flag.
    Never raises.
    """
    try:
        if m is None:
            print(f"[ligprep][{tag}] DEBUG_STATS: mol=None")
            return
        atoms = list(m.GetAtoms())
        elems = Counter(a.GetSymbol() for a in atoms)
        # Ring count (guarded)
        try:
            ri = m.GetRingInfo()
            ring_count = ri.NumRings() if ri is not None else 0
        except Exception:
            ring_count = -1
        # Formal charge sum
        try:
            chg_sum = sum(int(a.GetFormalCharge()) for a in atoms)
        except Exception:
            chg_sum = 0
        # Explicit H presence
        try:
            has_exp_h = any(a.GetNumExplicitHs() > 0 for a in atoms)
        except Exception:
            has_exp_h = False
        print(
            f"[ligprep][{tag}] atoms={len(atoms)} elems={dict(elems)} "
            f"rings={ring_count} charge_sum={chg_sum} explicitH={int(has_exp_h)}"
        )
    except Exception as e:
        print(f"[ligprep][{tag}] DEBUG_STATS_ERROR: {e}")


def collapse_sanitized_once(p: Path) -> Path:
    """Collapse repeated '.sanitized' tokens in the basename (single pass)."""
    return collapse_sanitized_path(p)


def _ph_label(ph: float) -> str:
    """Canonical label for ligand pH subdirectories, e.g. 7.0 -> 'pH7_0'."""
    return f"pH{ph:.1f}".replace(".", "_")


# =========================
# Parent chooser
# =========================


def _to_parent_mol(m: Chem.Mol) -> Optional[Chem.Mol]:
    try:
        if _HAS_STD and _std is not None:
            p = _std.Cleanup(m)
            chooser = _std.LargestFragmentChooser(preferOrganic=True)
            p = chooser.choose(p)
            p = _std.ChargeParent(p)
            if (
                _looks_like_buffer_salt(p)
                or _matches_counterion(p)
                or _is_polyacidic_buffer_like(p)
                or _buffer_like_by_counts_from_mol(p)
            ):
                frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True)
                p = max(frags, key=lambda x: x.GetNumHeavyAtoms())
                p = _std.ChargeParent(p)
        else:
            p = _REM.StripMol(m, dontRemoveEverything=True)
            frags = Chem.GetMolFrags(p, asMols=True, sanitizeFrags=True)
            p = max(frags, key=lambda x: x.GetNumHeavyAtoms())

        Chem.SanitizeMol(p)
        try:
            Chem.SetAromaticity(p, Chem.AromaticityModel.AROMATICITY_RDKIT)
        except Exception:
            pass

        if not any(a.GetSymbol() == "C" for a in p.GetAtoms()):
            return None
        if _is_polyacidic_buffer_like(p) or _buffer_like_by_counts_from_mol(p):
            return None
        if p.GetNumHeavyAtoms() < MIN_PARENT_HEAVY:
            return None
        if _looks_like_buffer_salt(p) or _matches_counterion(p):
            return None
        return p
    except Exception:
        return None




# Optional pKa/logD dependency: try to import, else fall back to cLogP
try:
    _pka = importlib.import_module("pkasolver").pkasolver
    _HAS_PKASOLVER = True
except Exception:
    _pka = None
    _HAS_PKASOLVER = False

# SCAM substructure alerts
SCAM_SMARTS: Dict[str, str] = {
    # electrophiles / reactive
    "epoxide": "[OX2r3]",
    "aziridine": "[NX3r3]",
    "alkyl_halide": "[CX4;H0,H1,H2][Cl,Br,I,F]",
    "michael_acceptor": "[C,c]=[C,c]-[C,S](=O)[O,N,S] | [C,c]=[C,c]-C(=O)[O,N,S]",
    "acrylamide": "C=CC(=O)N",
    "isothiocyanate": "N=C=S",
    "sulfonyl_fluoride": "S(=O)(=O)F",
    # redox / interference
    "p_quinone": "O=C1C=CC(=O)C=C1",
    "o_quinone": "O=c1ccc(=O)[cH][cH]1",
    "phenothiazine_like": "n2c1ccccn1Sc3ccccc23",
    # chelators / aggregators
    "catechol": "c1cc(O)c(O)cc1",
    "hydroxamate": "C(=O)N[OH]",
    "8_hydroxyquinoline": "Oc1cccc2ncccc12",
    "tannin_polyphenol": "c(O)c(O)c(O)",
    # nucleophiles / potentially reactive
    "hydrazine": "NN",
    "hydroxylamine": "N[OH]",
    "thiol": "[SH]",
    "dithiol": "SCCS",
    # rhodanine / known hitters
    "rhodanine": "O=C1NC(=S)SC1",
    "barbiturate_like": "O=C1NC(=O)NC(=O)1",
}
SCAM_QUERIES = {name: Chem.MolFromSmarts(s) for name, s in SCAM_SMARTS.items()}


def scam_flags(mol: Chem.Mol) -> List[str]:
    flags = []
    for name, patt in SCAM_QUERIES.items():
        if patt is not None and mol.HasSubstructMatch(patt):
            flags.append(f"hard: {name}")
    return flags


def predict_logD(mol: Chem.Mol, ph: float = 7.4) -> float:
    """
    Prefer pkasolver logD if available, otherwise gracefully fall back to cLogP.
    """
    try:
        if _HAS_PKASOLVER and _pka is not None:
            smiles = Chem.MolToSmiles(mol)
            return float(_pka.calculate_logd(smiles, ph=ph))
    except Exception as e:
        print(f"[WARN] logD prediction failed for {Chem.MolToSmiles(mol)}: {e}")
    return float(Crippen.MolLogP(mol))


def annotate_ligand_with_scam(lig_path: str, ligand_record: Dict) -> Dict:
    """
    Annotates the ligand record with SCAM filter results (cLogP, logD7.4, SCAM_Flags, SMILES).
    """
    try:
        # load molecule from file
        if lig_path.endswith(".sdf") or lig_path.endswith(".mol"):
            mol = Chem.MolFromMolFile(lig_path, sanitize=True)
        elif lig_path.endswith(".mol2"):
            mol = Chem.MolFromMol2File(lig_path, sanitize=True)
        else:
            raise ValueError(f"Unsupported ligand format: {lig_path}")

        if mol is None:
            ligand_record["SCAM_Flags"] = "InvalidMol"
            return ligand_record

        # get ligand descriptors
        smiles = Chem.MolToSmiles(mol)
        clogp = Crippen.MolLogP(mol)
        logd = predict_logD(mol, ph=7.4)

        # check for typical SCAM substructures
        flags = scam_flags(mol)

        # soft flags based on property cutoffs
        if clogp is not None and clogp > 3.5:
            flags.append(f"soft:high_logP({clogp:.2f})")

        if logd is not None:
            if logd < 0:
                flags.append(f"soft:low_logD({logd:.2f})")
            elif logd > 3.5:
                flags.append(f"soft:high_logD({logd:.2f})")
            if logd > 5:
                flags.append(f"hard:very_high_logD({logd:.2f})")

        flag_str = ";".join(flags) if flags else "None"

        ligand_record.update(
            {"SMILES": smiles, "cLogP": clogp, "logD7.4": logd, "SCAM_Flags": flag_str}
        )

    except Exception as e:
        ligand_record["SCAM_Flags"] = f"Error:{e}"

    return ligand_record


# =========================
# OBabel-friendly SDF writer
# =========================


def _write_obabel_friendly_sdf(mol: Chem.Mol, out_path: Path) -> bool:
    try:
        m = Chem.Mol(mol)
        Chem.SanitizeMol(m)
        try:
            # Helps fix some O valence issues due to charge misassignment
            Chem.Kekulize(m, clearAromaticFlags=True)
        except Exception:
            pass
        try:
            rdMolOps.AssignFormalCharges(m)
        except Exception:
            pass
        w = Chem.SDWriter(str(out_path))
        w.write(m)
        w.close()
        return True
    except Exception:
        try:
            m2 = Chem.Mol(mol)
            Chem.SanitizeMol(m2)
            Chem.SetAromaticity(m2, Chem.AromaticityModel.AROMATICITY_RDKIT)
            w = Chem.SDWriter(str(out_path))
            try:
                w.SetKekulize(False)
            except Exception:
                pass
            w.write(m2)
            w.close()
            return True
        except Exception:
            return False


# =========================
# OBabel helpers
# =========================


def _reserialize_mol_via_obabel(mol, obabel_exe_short: str, target_mol2: Path):
    tmp_sdf = target_mol2.with_suffix(".std.sdf")
    if not _write_obabel_friendly_sdf(mol, tmp_sdf):
        _log_malformed(
            target_mol2,
            "rdkit_sdf_write_fail:kekulize_or_aromaticity",
            log_dir=ligprep_common.MALFORMED_DIR,
        )
        return None

    fresh_mol2 = target_mol2.with_suffix(".std.mol2")
    cmd = [
        obabel_exe_short,
        "-isdf",
        str(tmp_sdf),
        "--gen3d",
        "-omol2",
        "-O",
        str(fresh_mol2),
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        _log_malformed(
            target_mol2, f"obabel_reserialize_fail:{getattr(e, 'stderr', '')[:200]}"
        )
        try:
            tmp_sdf.unlink(missing_ok=True)
        except Exception:
            pass
        return None

    try:
        tmp_sdf.unlink(missing_ok=True)
    except Exception:
        pass

    if not fresh_mol2.exists() or fresh_mol2.stat().st_size < 100:
        _log_malformed(fresh_mol2, "fresh_mol2_empty_or_small")
        return None
    return fresh_mol2


# =========================
# Open Babel runners
# =========================


def rdkit_embed_sdf_to_mol2(
    sdf_in: Path,
    mol2_out_dir: Path,
    obabel_exe: str,
    max_workers: int = 8,
    only_set: Optional[Set[str]] = None,
) -> List[Path]:
    # Ensure both imports exist in this scope
    from rdkit import Chem
    from rdkit.Chem import AllChem

    sdf_stem = sdf_in.stem
    name_prefix = sdf_stem

    allowed_indices: Optional[Set[int]] = None
    if only_set:
        try:
            allowed_indices = {
                int(tok.split("_", 1)[1])
                for tok in only_set
                if tok.startswith("rdk_") and tok.split("_", 1)[1].isdigit()
            }
        except Exception:
            allowed_indices = None

    suppl = Chem.SDMolSupplier(str(sdf_in), removeHs=False, sanitize=False)

    mols: List[Tuple[int, Chem.Mol]] = []
    raw_count = 0
    for i, m in enumerate(suppl):
        raw_count += 1
        if m is None:
            continue
        if allowed_indices is not None and i not in allowed_indices:
            continue
        mols.append((i, m))

    if allowed_indices is not None:
        print(
            f"[test-mode.rdkit] RDKit supplier read {raw_count} records from {sdf_in.name}; "
            f"kept {len(mols)} based on ONLY filter (requested={len(allowed_indices)})"
        )
    print(f"RDKit: loaded {len(mols)} molecules from {sdf_in.name}")

    parent_names: List[str] = []
    debug_limit = 5
    if mols:
        max_idx = max(i for i, _ in mols)
        parent_names = ["" for _ in range(max_idx + 1)]
        for idx, mol in mols:
            raw_name = ""
            parent_name = f"{name_prefix}_{idx + 1:05d}"
            parent_names[idx] = parent_name
            if idx < debug_limit:
                print(
                    f"[rdkit-debug] idx={idx} raw_name={raw_name!r} "
                    f"fallback_base={name_prefix!r} stored_parent={parent_name!r}"
                )

    sdf_tmp_dir = mol2_out_dir / "_rdkit_embedded_sdf"
    sdf_tmp_dir.mkdir(parents=True, exist_ok=True)

    def _embed_one(i_m):
        i, m = i_m


        logging.info(
            "[bulkSDF] pre-standardize ligand=rdk_%07d atoms=%d charge=%+d elements=%s has_conf=%s",
            i,
            (m.GetNumAtoms() if m else -1),
            _formal_charge(m),
            _elem_hist(m),
            bool(m and m.GetNumConformers() > 0),
        )

        try:
            Chem.SanitizeMol(m)
            logging.info(
                "[bulkSDF] post-sanitize ligand=rdk_%07d ok=True charge=%+d",
                i,
                _formal_charge(m),
            )
        except Exception as e:
            print(f"[ligprep] SanitizeMol hard fail: {e} (trying staged flags)")
            try:
                rdMolOps.SanitizeMol(
                    m,
                    sanitizeOps=rdMolOps.SanitizeFlags.SANITIZE_FINDRADICALS
                    | rdMolOps.SanitizeFlags.SANITIZE_KEKULIZE
                    | rdMolOps.SanitizeFlags.SANITIZE_SETAROMATICITY
                    | rdMolOps.SanitizeFlags.SANITIZE_ADJUSTHS,
                )
                print("[ligprep] staged-sanitize succeeded")
                logging.info(
                    "[bulkSDF] post-sanitize ligand=rdk_%07d ok=True(staged) charge=%+d",
                    i,
                    _formal_charge(m),
                )
            except Exception as e2:
                _dbg_atom_stats(m, "sanitize-fail-snapshot")
                logging.info(
                    "[bulkSDF] post-sanitize ligand=rdk_%07d ok=False err=%s",
                    i,
                    str(e2)[:180],
                )
                raise  # keep failing fast, but with context

        orig_heavy = m.GetNumHeavyAtoms()
        p = _to_parent_mol(m)
        if p is None:
            _log_malformed(
                Path(f"rdk_{i:07d}"), "no_parent_or_too_small_after_desalting"
            )
            return None
        if p is not m and p.GetNumHeavyAtoms() != orig_heavy:
            logging.info(
                f"[parent-pick] rdk_{i:07d}: {orig_heavy}?{p.GetNumHeavyAtoms()} heavy atoms"
            )

        m = p

        logging.info(
            "[bulkSDF] post-standardize ligand=rdk_%07d atoms=%d charge=%+d elements=%s",
            i,
            m.GetNumAtoms(),
            _formal_charge(m),
            _elem_hist(m),
        )

        try:
            params = AllChem.ETKDGv3()
            params.randomSeed = 0xC0FFEE

            ok = AllChem.EmbedMolecule(m, params)
            if ok != 0:
                return None
            try:
                AllChem.UFFOptimizeMolecule(m, maxIters=200)
            except Exception:
                pass
            out = sdf_tmp_dir / f"rdk_{i:07d}.sdf"
            if _write_obabel_friendly_sdf(m, out):
                logging.info("[paths] sdf_out=%s", str(out.resolve()))
                return out

            else:
                _log_malformed(out, "write_obabel_friendly_sdf_failed")
                return None
        except Exception:
            return None

    paths: List[Path] = []
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        for res in ex.map(_embed_one, mols):
            if res:
                paths.append(res)

    print(f"RDKit: embedded {len(paths)} molecules; converting to MOL2")
    out_files: List[Path] = []

    def _stem_to_rdk_index(stem: str) -> Optional[int]:
        if stem.startswith("rdk_"):
            try:
                return int(stem.split("_", 1)[1])
            except Exception:
                return None
        return None

    for pth in paths:
        out = mol2_out_dir / (pth.stem + ".mol2")
        print(
            f"[ligprep] pre-MOL2-write: rdkit_embed sdf={pth.name} -> mol2={out.name}"
        )
        logging.info(
            "[bulkSDF] write.mol2 ligand=%s via=obabel out=%s",
            pth.stem,
            str(out.resolve()),
        )
        parent_name = ""
        idx = _stem_to_rdk_index(pth.stem)
        if idx is not None and 0 <= idx < len(parent_names):
            parent_name = parent_names[idx]
        name_path = out.with_suffix(".name")
        try:
            to_write = (parent_name or "").rstrip("\r\n")
            name_path.write_text(to_write + "\n", encoding="utf-8")
        except Exception as e:
            logging.warning("[bulkSDF] name_write_failed stem=%s err=%s", pth.stem, e)
        if idx is not None and idx < debug_limit:
            exists = name_path.exists()
            try:
                contents = name_path.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                contents = "<unreadable>"
            print(
                f"[rdkit-debug] mol2={out.name} name_path={name_path.name} "
                f"exists={exists} contents={contents!r}"
            )
        if idx is not None and idx < debug_limit:
            print(
                f"[rdkit-debug] writing .name for {out.name}: parent_name={parent_name!r}"
            )

        if out.exists() and out.stat().st_size > 100:
            out_files.append(out)
            try:
                size = out.stat().st_size
                logging.info(
                    "[bulkSDF] write.mol2.done ligand=%s rc=0(size_only) size=%d",
                    pth.stem,
                    size,
                )
            except Exception:
                pass
            continue

        cmd = [obabel_exe, "-isdf", str(pth), "-omol2", "-O", str(out)]
        if _run_obabel(cmd, timeout_sec=120):
            out_files.append(out)
            try:
                size = out.stat().st_size if out.exists() else 0
                logging.info(
                    "[bulkSDF] write.mol2.done ligand=%s rc=%s size=%d",
                    pth.stem,
                    "0" if size > 0 else "unknown",
                    size,
                )
            except Exception:
                pass

            try:
                n_mol2_atoms = 0
                with open(out, "r", errors="ignore") as fh:
                    in_atoms = False
                    for ln in fh:
                        s = ln.strip()
                        if s.startswith("@<TRIPOS>ATOM"):
                            in_atoms = True
                            continue
                        if s.startswith("@<TRIPOS>") and in_atoms:
                            break
                        if in_atoms and s and s[0].isdigit():
                            n_mol2_atoms += 1
                print(f"[ligprep] MOL2_written≈{n_mol2_atoms} file={out.name}")
            except Exception as e:
                print(f"[ligprep] MOL2_count_error {out.name}: {e}")

    print(f"RDKit: wrote {len(out_files)} MOL2 files")
    return out_files


# =========================
# Optional: pre-clean entire SDF to parents only
# =========================


def _write_parent_only_sdf(src_sdf: Path, dst_sdf: Path) -> int:
    suppl = Chem.SDMolSupplier(str(src_sdf), removeHs=False, sanitize=False)
    w = Chem.SDWriter(str(dst_sdf))
    try:
        w.SetKekulize(False)
    except Exception:
        pass
    kept = 0
    for m in suppl:
        if m is None:
            continue
        try:
            Chem.SanitizeMol(m)
        except Exception:
            continue
        p = _to_parent_mol(m)
        if p is None:
            continue
        w.write(p)
        kept += 1
    w.close()
    return kept


# =========================
# MGLTools: prepare_ligand4 & validation
# =========================


# === AROMATICITY AUDIT & RESCUE ===


def load_mol2_lenient(path, logger):
    mol = Chem.MolFromMol2File(str(path), sanitize=False, removeHs=False)
    if mol is None:
        _log_malformed(Path(path), "rdkit_read_fail")
        return None
    try:
        # --- DEBUG: pre-sanitize fingerprint ---

        Chem.SanitizeMol(mol)
        _dbg_atom_stats(mol, "post-sanitize")
        ok, why = _quick_filters(mol)
        if not ok:
            _log_malformed(Path(path), why)
            return None
        return mol
    except Exception as e:
        std_mol, why = _standardize_then_sanitize(mol)
        if std_mol is None:
            _log_malformed(Path(path), f"sanitize_fail:{e}|{why}")
            if logger:
                logger.warning(f"RDKit failed to sanitize: {path} ({e})")
                # On sanitize/valence error: dump quick forensics
                try:
                    smi = (
                        Chem.MolToSmiles(mol, isomericSmiles=True)
                        if mol is not None
                        else "None"
                    )
                except Exception:
                    smi = "MolToSmiles_failed"
                try:
                    from collections import Counter

                    elem_hist = (
                        dict(Counter(a.GetSymbol() for a in mol.GetAtoms()))
                        if mol is not None
                        else {}
                    )
                    fcharge = sum(
                        int(a.GetFormalCharge())
                        for a in (mol.GetAtoms() if mol else [])
                    )
                except Exception:
                    elem_hist, fcharge = {}, 0
                logging.warning(
                    "[ligprep] sanitize_failed elem=%s formal_charge=%d smiles=%s",
                    elem_hist,
                    fcharge,
                    smi,
                )
            return None
        return std_mol


# =========================
# Main SDF ? MOL2 ? PDBQT pipeline
# =========================


def _valid_pdbqt(path: Path, log_dir: Path) -> bool:
    """
    Lightweight sanity check for PDBQT files during docking enumeration.

    We assume heavy validation has already happened at prep time.
    Here we only enforce existence and a minimum size threshold.
    """
    try:
        return path.exists() and path.stat().st_size > 100
    except Exception:
        return False


def _env_bool(name: str, default: Any) -> bool:
    v = os.environ.get(name)
    return (
        (str(v).strip().lower() in {"1", "true", "yes", "on"})
        if v is not None
        else bool(default)
    )


def _env_int(name: str, default: Any) -> int:
    v = os.environ.get(name)
    return int(v) if v not in (None, "") else int(default)


def _env_float(name: str, default: Any) -> float:
    v = os.environ.get(name)
    return float(v) if v not in (None, "") else float(default)


def _sanitize_ligand_name_for_filename(name: str) -> str:
    name = name.strip()
    name = name.replace(" ", "_")
    name = re.sub(r"[^A-Za-z0-9_.+-]+", "_", name)
    return name[:80] or "ligand"


def _relative_to_output(path: Path, output_dir: Path) -> str:
    try:
        return str(path.relative_to(output_dir))
    except Exception:
        return path.name


def _bulk_init_config_and_paths(
    *,
    force: bool,
    microstate_dedup: bool,
    ph_values: Optional[List[float]],
    in_sdf: Optional[Path],
    in_sdf_dir: Optional[Path],
    in_pdb_dir: Optional[Path],
    mol2_dir: Optional[Path],
    out_pdbqt_dir: Optional[Path],
    status_log: Optional[Path],
    root_dir: Optional[Path],
) -> Dict[str, Any]:
    """Phase 1: initialize config, effective paths, and environment-driven overrides."""
    ctx: Dict[str, Any] = {}
    force_flag = force
    if not force_flag:
        env_force = os.environ.get("LIGPREP_FORCE", "").strip().lower()
        force_flag = env_force in {"1", "true", "yes", "y"}
    ctx["force"] = force_flag

    print("Starting ligand preparation")

    if microstate_dedup:
        logger.info(
            "prep_ligands: microstate_dedup=True (stage 4: microstate registry + canonical PDBQT shadow)"
        )

    cfg = load_config()
    validate_config(cfg)
    ctx["cfg"] = cfg
    pdb_token = (
        os.environ.get("PDB_ID")
        or cfg.get("PDB_ID")
        or cfg.get("TARGET_PDB")
        or cfg.get("PDB")
        or cfg.get("INPUT_PDB")
        or cfg.get("PDB_FILE")
        or ""
    )
    pdb_id = Path(str(pdb_token)).stem.upper() if pdb_token else "LIGPREP"
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ctx["paths"] = paths

    ctx["in_sdf_env"] = (os.environ.get("LIGPREP_IN_SDF", "") or "").strip() or None
    ctx["in_sdf_dir_env"] = (
        os.environ.get("LIGPREP_IN_SDF_DIR", "") or ""
    ).strip() or None
    ctx["in_pdb_dir_env"] = (
        os.environ.get("LIGPREP_IN_PDB_DIR", "") or ""
    ).strip() or None
    ctx["mol2_dir_env"] = (os.environ.get("LIGPREP_MOL2_DIR", "") or "").strip() or None
    ctx["out_dir_env"] = (os.environ.get("LIGPREP_OUT_DIR", "") or "").strip() or None
    ctx["status_log_env"] = (
        os.environ.get("LIGPREP_STATUS_LOG", "") or ""
    ).strip() or None
    rename_prefix = (os.environ.get("LIGPREP_RENAME_PREFIX", "") or "").strip()
    try:
        rename_pad = int(os.environ.get("LIGPREP_RENAME_PAD", "5"))
    except Exception:
        rename_pad = 5
    try:
        rename_start = int(os.environ.get("LIGPREP_RENAME_START", "1"))
    except Exception:
        rename_start = 1
    rename_force = (
        os.environ.get("LIGPREP_RENAME_FORCE", "") or ""
    ).strip().lower() in {"1", "true", "yes", "y"}
    ctx["rename_prefix"] = rename_prefix
    ctx["rename_pad"] = max(1, rename_pad)
    ctx["rename_start"] = max(1, rename_start)
    ctx["rename_force"] = rename_force
    ctx["rename_active"] = bool(rename_prefix or rename_force)

    ligands_raw_dir = paths.ligand_output_dir
    prepped_lig_dir = paths.prepped_ligands_dir
    ligands_mol2_root = paths.ligands_mol2_dir

    library_env = (os.environ.get("LIGPREP_LIBRARY", "") or "").strip()
    library_base = ""

    if root_dir is not None:
        root_dir = Path(root_dir).resolve()
        output_ligands_dir = root_dir
        prepped_ligands_dir = output_ligands_dir
        library_hint = prepped_ligands_dir.name.lower()
    else:
        output_ligands_dir = (
            Path(ctx["out_dir_env"]).resolve()
            if ctx["out_dir_env"]
            else prepped_lig_dir
        )
        prepped_ligands_dir = output_ligands_dir
        library_hint = prepped_ligands_dir.name.lower()

    if library_env:
        library_base = library_env.lower()
    elif ctx["in_sdf_env"]:
        try:
            in_sdf_path = Path(ctx["in_sdf_env"]).resolve()
            detected = ""
            parts = list(in_sdf_path.parts)
            for idx, part in enumerate(parts):
                if part == "extracted_ligands" and idx + 1 < len(parts):
                    detected = parts[idx + 1]
                    break
            library_base = (detected or in_sdf_path.stem).lower()
        except Exception:
            library_base = ""
    else:
        library_base = library_hint
    library_base = (library_base or library_hint or "ligprep").lower()
    library_name = library_base
    ctx["library_name"] = library_name

    if root_dir is not None and not ctx["in_sdf_dir_env"]:
        project_root = None
        for parent in root_dir.parents:
            if parent.name == "prepped_ligands":
                project_root = parent.parent
                break
        if project_root is None:
            project_root = root_dir.parent.parent

        candidate_extracted = project_root / "extracted_ligands" / library_base
        if candidate_extracted.is_dir():
            ligand_extracted_dir = candidate_extracted
            logger.info(
                "prep_ligands: using extracted_ligands source path=%s",
                candidate_extracted,
            )
        else:
            ligand_extracted_dir = ligands_raw_dir
    else:
        ligand_extracted_dir = (
            Path(ctx["in_sdf_dir_env"]).resolve()
            if ctx["in_sdf_dir_env"]
            else ligands_raw_dir
        )
    ctx["ligand_extracted_dir"] = ligand_extracted_dir

    if ctx["mol2_dir_env"]:
        ligands_mol2_dir = Path(ctx["mol2_dir_env"]).resolve()
    else:
        ligands_mol2_dir = Path(ligands_mol2_root) / library_base
    ctx["ligands_mol2_dir"] = ligands_mol2_dir

    if not library_env and root_dir is None:
        os.environ["LIGPREP_LIBRARY"] = library_base

    ctx["is_fda_library"] = library_hint == "fda"

    ligprep_common.MALFORMED_DIR = prepped_ligands_dir

    output_ligands_dir.mkdir(parents=True, exist_ok=True)
    ligands_mol2_dir.mkdir(parents=True, exist_ok=True)

    ctx["library_out_dir"] = prepped_ligands_dir
    if microstate_dedup:
        microstates_dir = prepped_ligands_dir / "microstates"
        microstates_dir.mkdir(parents=True, exist_ok=True)
        microstate_registry, microstate_index = load_microstate_registry(
            prepped_ligands_dir, library_name
        )
        microstate_registry_dirty = False

        alias_index: Dict[Tuple[str, str, float], dict] = {}
        for entry in microstate_registry.get("microstates", []) or []:
            for alias in entry.get("aliases", []) or []:
                try:
                    key = (
                        str(alias.get("ligand_stem", "")),
                        str(alias.get("ph_label", "")),
                        float(alias.get("ph_value", 0.0)),
                    )
                except Exception:
                    continue
                alias_index[key] = entry
        ctx["microstates_dir"] = microstates_dir
        ctx["microstate_registry"] = microstate_registry
        ctx["microstate_index"] = microstate_index
        ctx["alias_index"] = alias_index
        ctx["microstate_registry_dirty"] = microstate_registry_dirty
    else:
        ctx["microstates_dir"] = None
        ctx["microstate_registry"] = None
        ctx["microstate_index"] = None
        ctx["alias_index"] = None
        ctx["microstate_registry_dirty"] = False

    if (
        microstate_dedup
        and ph_values is not None
        and ctx.get("microstate_registry") is not None
    ):
        logger.info(
            "prep_ligands: microstate_dedup init library=%s out_dir=%s existing_microstates=%d",
            library_name,
            prepped_ligands_dir,
            len(ctx["microstate_registry"].get("microstates", []) or []),
        )

    def _maybe_save_microstate_registry(force: bool = False) -> None:
        if (
            not microstate_dedup
            or ctx.get("microstate_registry") is None
            or ctx.get("library_out_dir") is None
        ):
            return
        microstate_registry_local = ctx["microstate_registry"]
        microstate_list = microstate_registry_local.get("microstates") or []
        alias_count = sum(len(entry.get("aliases") or []) for entry in microstate_list)
        if microstate_dedup:
            if microstate_list:
                logger.info(
                    "prep_ligands: microstate_dedup summary library=%s out_dir=%s microstates=%d aliases=%d",
                    library_name,
                    ctx["library_out_dir"],
                    len(microstate_list),
                    alias_count,
                )
            else:
                logger.warning(
                    "prep_ligands: microstate_dedup summary library=%s out_dir=%s microstates=%d aliases=%d (EMPTY REGISTRY)",
                    library_name,
                    ctx["library_out_dir"],
                    len(microstate_list),
                    alias_count,
                )
        if ctx.get("microstate_registry_dirty") or force:
            save_microstate_registry(ctx["library_out_dir"], microstate_registry_local)
            ctx["microstate_registry_dirty"] = False

    ctx["maybe_save_microstate_registry"] = _maybe_save_microstate_registry

    print(
        "[paths.effective]"
        f" in_sdf={ctx['in_sdf_env'] or 'None'}"
        f" in_sdf_dir={ligand_extracted_dir}"
        f" in_pdb_dir={ctx['in_pdb_dir_env'] or 'None'}"
        f" mol2_dir={ligands_mol2_dir}"
        f" out_pdbqt_dir={output_ligands_dir}"
        f" status_log={(ctx['status_log_env'] or cfg.get('LIGAND_STATUS_LOG_BASENAME', 'ligand_prep_status.tsv'))}"
    )

    ctx["USE_RDKIT_FOR_3D"] = _env_bool(
        "USE_RDKIT_FOR_3D", cfg.get("USE_RDKIT_FOR_3D", True)
    )
    ctx["OBABEL_THREADS"] = _env_int("OBABEL_THREADS", cfg.get("OBABEL_THREADS", 50))
    ctx["OBABEL_TIMEOUT_S"] = _env_int(
        "OBABEL_TIMEOUT_S", cfg.get("OBABEL_TIMEOUT_S", 900)
    )
    ctx["CHUNK_SIZE"] = _env_int(
        "LIGPREP_CHUNK_SIZE", cfg.get("LIGPREP_CHUNK_SIZE", 200)
    )
    ctx["LIGPREP_PH"] = _env_float("LIGPREP_PH", cfg.get("LIGPREP_PH", 7.4))
    ctx["KEEP_NONPOLAR_H"] = _env_int("KEEP_NONPOLAR_H", cfg.get("KEEP_NONPOLAR_H", 1))
    ctx["MAX_HEAVY_ATOMS"] = _env_int(
        "MAX_HEAVY_ATOMS", cfg.get("MAX_HEAVY_ATOMS", 1200)
    )
    ctx["MIN_ATOMS_FOR_DOCKING"] = _env_int(
        "MIN_ATOMS_FOR_DOCKING", cfg.get("MIN_ATOMS_FOR_DOCKING", 5)
    )
    ctx["MIN_PARENT_HEAVY"] = _env_int(
        "MIN_PARENT_HEAVY", cfg.get("MIN_PARENT_HEAVY", 8)
    )
    ctx["MIN_TORS_DOF"] = _env_int("MIN_TORS_DOF", cfg.get("MIN_TORS_DOF", 0))

    if ph_values is not None and len(ph_values) > 0:
        eff_ph_values = [float(x) for x in ph_values]
    else:
        eff_ph_values = [float(ctx["LIGPREP_PH"])]
    ctx["eff_ph_values"] = eff_ph_values

    multiple_ph = ph_values is not None and len(eff_ph_values) > 1
    use_ph_subdirs = ph_values is not None
    ph_source = "python" if ph_values is not None else "config"
    ph_summary = ",".join(f"{ph:.1f}" for ph in eff_ph_values)
    print(
        f"[ligprep] ph_schedule source={ph_source} values={ph_summary} multiple={multiple_ph}"
    )
    logger.info("prep_ligands: pH schedule eff_ph_values=%s", eff_ph_values)

    output_ligands_dir.mkdir(parents=True, exist_ok=True)
    if "prepped_ligands" in str(ligands_mol2_dir):
        suggested = Path(
            str(ligands_mol2_dir).replace("prepped_ligands", "ligands_mol2")
        ).resolve()
        logging.warning(
            "[compat] LIGANDS_MOL2_DIR points at prepped_ligands; redirecting to %s",
            suggested,
        )
        ligands_mol2_dir = suggested
        ctx["ligands_mol2_dir"] = ligands_mol2_dir

    ligands_mol2_dir.mkdir(parents=True, exist_ok=True)

    print(
        f"[ligprep] sdf_dir={ligands_mol2_dir / '_rdkit_embedded_sdf'} "
        f"mol2_dir={ligands_mol2_dir} "
        f"pdbqt_dir={output_ligands_dir}"
    )

    mgltools_python = str(
        resolve_tool(cfg, "MGLTOOLS_PYTHON", "pythonsh").get("resolved_path", "") or ""
    ).strip()
    mgltools_path = str(cfg.get("MGLTOOLS_PATH", "") or "").strip()
    obabel_exe = str(
        resolve_tool(cfg, "OPENBABEL_PATH", "obabel").get("resolved_path", "") or ""
    ).strip()

    for label, p in [
        ("MGLTOOLS_PYTHON", mgltools_python),
        ("OPENBABEL_PATH", obabel_exe),
        ("EXTRACTED_LIGANDS_DIR", ligand_extracted_dir),
        ("LIGANDS_MOL2_DIR", ligands_mol2_dir),
        ("PREPPED_LIGANDS_DIR", output_ligands_dir),
    ]:
        if not str(p).strip():
            raise RuntimeError(f"Config value missing/empty: {label}")

    if not os.environ.get("BABEL_DATADIR"):
        obabel_dir = Path(obabel_exe).resolve().parent
        data_dir = obabel_dir / "data"
        if data_dir.exists():
            os.environ["BABEL_DATADIR"] = str(data_dir)

    obabel_exe_short = get_short_path_name(obabel_exe)
    mgltools_python_short = get_short_path_name(mgltools_python)
    prepare_script = _resolve_prepare_ligand4(mgltools_path, cfg)
    if not prepare_script.exists():
        raise FileNotFoundError(
            f"prepare_ligand4.py not found at {prepare_script} (set PREPARE_LIGAND_SCRIPT in config.txt)"
        )
    prepare_script_short = get_short_path_name(str(prepare_script.resolve()))

    status_log_path = (
        Path(ctx["status_log_env"]).resolve()
        if (ctx["status_log_env"] and Path(ctx["status_log_env"]).suffix)
        else (
            output_ligands_dir
            / (
                ctx["status_log_env"]
                or cfg.get("LIGAND_STATUS_LOG_BASENAME", "ligand_prep_status.tsv")
            )
        )
    )

    ctx.update(
        {
            "mgltools_python_short": mgltools_python_short,
            "prepare_script_short": prepare_script_short,
            "obabel_exe_short": obabel_exe_short,
            "status_log": status_log_path,
            "output_ligands_dir": output_ligands_dir,
            "prepped_ligands_dir": prepped_ligands_dir,
            "library_hint": library_hint,
            "microstate_dedup": microstate_dedup,
            "use_ph_subdirs": use_ph_subdirs,
            "ph_values": ph_values,
            "in_sdf": in_sdf,
            "in_sdf_dir": in_sdf_dir,
            "in_pdb_dir": in_pdb_dir,
            "mol2_dir": mol2_dir,
            "out_pdbqt_dir": out_pdbqt_dir,
        }
    )
    return ctx


def _bulk_select_unit_sdfs(
    ctx: Dict[str, Any], only_set: Optional[Set[str]]
) -> Tuple[List[Path], Optional[Set[str]], bool]:
    """Phase 2: inspect per-ligand SDFs, handle true test-mode branch, and return unit sdf metadata."""
    ligands_mol2_dir = ctx["ligands_mol2_dir"]
    output_ligands_dir = ctx["output_ligands_dir"]
    rdkit_unit_sdf_dir = ligands_mol2_dir / "_rdkit_embedded_sdf"
    unit_sdfs = (
        sorted(rdkit_unit_sdf_dir.glob("*.sdf")) if rdkit_unit_sdf_dir.is_dir() else []
    )

    in_sdf_env = ctx["in_sdf_env"]
    in_pdb_dir_env = ctx["in_pdb_dir_env"]
    has_only = bool(only_set)
    microstate_dedup = ctx["microstate_dedup"]
    ph_values = ctx.get("ph_values")

    if only_set:
        logging.info(
            "[ligprep] ONLY-set active for library=%s count=%d sample=%s",
            ctx["library_name"],
            len(only_set),
            ",".join(sorted(list(only_set))[:10]),
        )
    else:
        logging.info(
            "[ligprep] ONLY-set empty; full library will be processed for library=%s",
            ctx["library_name"],
        )

    test_mode_allowed = (in_sdf_env is None) and (in_pdb_dir_env is None)
    use_unit_sdf_test_mode = (
        unit_sdfs
        and test_mode_allowed
        and has_only
        and not (microstate_dedup and ph_values is not None)
    )

    if use_unit_sdf_test_mode:
        print(f"[test-mode] per-ligand SDFs detected dir={rdkit_unit_sdf_dir}")
        print(
            "[test-mode] using legacy unit-SDF pipeline (microstate_dedup=False or no ph_values); "
            "bulk SDF path is disabled."
        )
        selected_sdfs: List[Path]
        if only_set:
            requested_ids = sorted(only_set)
            selected_sdfs = [
                rdkit_unit_sdf_dir / f"{rid}.sdf"
                for rid in requested_ids
                if (rdkit_unit_sdf_dir / f"{rid}.sdf").exists()
            ]
            missing = [
                rid
                for rid in requested_ids
                if not (rdkit_unit_sdf_dir / f"{rid}.sdf").exists()
            ]
            print(
                f"[test-mode] enabled only_count={len(requested_ids)} "
                f"sample={','.join(requested_ids[:10])}"
            )
            print(
                f"[test-mode] selected_sdf_count={len(selected_sdfs)} missing={len(missing)}"
                + (f" first_missing={','.join(missing[:10])}" if missing else "")
            )
            if len(selected_sdfs) == 0:
                print(
                    "[test-mode] No requested per-ligand SDFs found. "
                    "Nothing to do; exiting cleanly."
                )
                ctx["maybe_save_microstate_registry"](force=True)
                return unit_sdfs, only_set, True
        else:
            selected_sdfs = unit_sdfs

        clean_env = os.environ.get("LIGPREP_CLEAN", "").strip().lower() in {
            "1",
            "true",
            "yes",
            "y",
        }
        if clean_env and selected_sdfs:
            cleaned = []
            for sdf_in in selected_sdfs:
                out_pdbqt = output_ligands_dir / (sdf_in.stem + ".pdbqt")
                if out_pdbqt.exists():
                    try:
                        out_pdbqt.unlink()
                        cleaned.append(out_pdbqt.name)
                    except Exception:
                        pass
            if cleaned:
                logging.info(
                    "[clean] removed %d existing PDBQTs: %s",
                    len(cleaned),
                    ",".join(cleaned[:20]),
                )

        rename_active = ctx.get("rename_active", False)
        rename_prefix = ctx.get("rename_prefix", "")
        rename_pad = int(ctx.get("rename_pad", 5))
        rename_start = int(ctx.get("rename_start", 1))

        mol2_files: List[Path] = []
        for idx, sdf_in in enumerate(selected_sdfs, start=rename_start):
            sdf_name = sdf_in.name
            mol2_stem = sdf_in.stem
            if rename_active:
                mol2_stem = f"{rename_prefix}{idx:0{rename_pad}d}"
            mol2_out = ligands_mol2_dir / f"{mol2_stem}.mol2"
            print(
                f"[ligprep] pre-MOL2-write: per-ligand sdf={sdf_name} -> mol2={mol2_out.name}"
            )
            ok, stderr_text = _sdf_to_mol2(
                sdf_in,
                mol2_out,
                ctx["obabel_exe_short"],
                timeout_sec=max(60, LIGPREP_OBABEL_TIMEOUT_SEC),
            )
            if not ok:
                print(
                    f"[ligprep] WARNING: Open Babel SDF->MOL2 failed for {sdf_name}: {stderr_text[:200]}"
                )
                continue
            mol2_files.append(mol2_out)

        if only_set:
            sample = ",".join(sorted({Path(p).stem for p in mol2_files})[:10])
            print(f"[test-mode] scheduling_only={len(mol2_files)} sample={sample}")

        print(
            f"Preparing {len(mol2_files)} MOL2 files with MGLTools (parallel) "
            f"into {output_ligands_dir}"
        )
        max_workers = max(1, int(os.environ.get("CPU", "8")))
        futures = []
        ok_count = 0
        fail_count = 0
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            for idx, mol2_file in enumerate(mol2_files, start=rename_start):
                pdbqt_stem = mol2_file.stem
                if rename_active:
                    pdbqt_stem = f"{rename_prefix}{idx:0{rename_pad}d}"
                pdbqt_path = output_ligands_dir / f"{pdbqt_stem}.pdbqt"
                if pdbqt_path.exists():
                    try:
                        pdbqt_path.unlink()
                        logging.info(
                            "[test-mode] overwriting existing PDBQT: %s",
                            pdbqt_path.name,
                        )
                    except Exception as e:
                        logging.warning(
                            "[test-mode] unable to remove existing PDBQT %s: %s",
                            pdbqt_path.name,
                            e,
                        )
                futures.append(
                    ex.submit(
                        _prepare_one,
                        ctx["mgltools_python_short"],
                        ctx["prepare_script_short"],
                        mol2_file,
                        pdbqt_path,
                        obabel_exe_short=ctx["obabel_exe_short"],
                        status_log_dir=output_ligands_dir,
                        ph=ctx["eff_ph_values"][0],
                    )
                )

            for fut in as_completed(futures):
                try:
                    _, status = fut.result()
                    if status == "ok":
                        ok_count += 1
                    else:
                        fail_count += 1
                except Exception:
                    fail_count += 1
        print(
            f"[test-mode] finished: ok={ok_count} fail={fail_count} "
            f"scheduled={len(futures)}"
        )
        ctx["maybe_save_microstate_registry"](force=True)
        return unit_sdfs, only_set, True

    if (
        unit_sdfs
        and test_mode_allowed
        and has_only
        and microstate_dedup
        and ph_values is not None
    ):
        print(
            "[test-mode] per-ligand SDFs detected but microstate_dedup=True with ph_values; "
            "using bulk RDKit SDF pipeline with ONLY filter instead of legacy unit-SDF path."
        )
    elif unit_sdfs and test_mode_allowed and not has_only:
        print(
            "[test-mode] per-ligand SDFs present but ONLY not set; using bulk ligprep path instead"
        )
    elif unit_sdfs and not test_mode_allowed:
        print(
            "[test-mode] per-ligand SDFs present at "
            f"{rdkit_unit_sdf_dir} but skipping because explicit input was provided: "
            f"in_sdf={in_sdf_env!r}, in_pdb_dir={in_pdb_dir_env!r}"
        )
    return unit_sdfs, only_set, False


def _bulk_generate_mol2_files(
    ctx: Dict[str, Any], only_set: Optional[Set[str]]
) -> List[Path]:
    """Phase 3: convert bulk SDFs (or explicit overrides) into MOL2 files."""
    ligand_extracted_dir = ctx["ligand_extracted_dir"]
    in_sdf_env = ctx["in_sdf_env"]
    ligands_mol2_dir = ctx["ligands_mol2_dir"]
    obabel_exe_short = ctx["obabel_exe_short"]
    output_ligands_dir = ctx["output_ligands_dir"]
    _maybe_save = ctx["maybe_save_microstate_registry"]

    sdf_files = (
        [Path(in_sdf_env).resolve()]
        if in_sdf_env
        else list(ligand_extracted_dir.glob("*.sdf"))
    )
    initial_count = len(sdf_files)
    sdf_files = [p for p in sdf_files if "docked" not in p.name.lower()]
    skipped = initial_count - len(sdf_files)
    if skipped > 0:
        logging.info(
            "[ligprep.skip] skipped %d SDF(s) with 'docked' in name under %s",
            skipped,
            ligand_extracted_dir,
        )

    print(f"Found {len(sdf_files)} SDF file(s)")

    if not sdf_files:
        _maybe_save(force=True)
        return []

    all_mol2_files: List[Path] = []
    use_rdkit = ctx["USE_RDKIT_FOR_3D"]
    obabel_threads = ctx["OBABEL_THREADS"]
    obabel_timeout = ctx["OBABEL_TIMEOUT_S"]
    chunk_size = ctx["CHUNK_SIZE"]

    for sdf_file in sdf_files:
        print(f"=== Processing SDF: {sdf_file.name} ===")
        norm_sdf = collapse_sanitized_once(Path(sdf_file))
        if norm_sdf != Path(sdf_file):
            try:
                norm_sdf.write_bytes(Path(sdf_file).read_bytes())
                logging.info(
                    "[tidy] normalized double-sanitized SDF -> %s", norm_sdf.name
                )
                sdf_file = norm_sdf
            except Exception as e:
                logging.warning("[tidy] unable to normalize SDF %s: %s", sdf_file, e)
        sdf_abs = Path(sdf_file).resolve()
        sdf_abs = sdf_file.resolve()
        mol2_files: List[Path]
        max_workers = max(1, int(os.environ.get("CPU", "8")))
        if use_rdkit:
            print(
                "Using RDKit ETKDG for 3D with parent-picking; OBabel only for format conversion "
            )
            mol2_files = rdkit_embed_sdf_to_mol2(
                sdf_abs,
                ligands_mol2_dir,
                obabel_exe=obabel_exe_short,
                max_workers=max_workers,
                only_set=only_set,
            )
        else:
            print("Using OBabel --gen3d; pre-cleaning SDF to parent-only ")
            cleaned_sdf = ligands_mol2_dir / (sdf_abs.stem + "_parents.sdf")
            n_kept = _write_parent_only_sdf(sdf_abs, cleaned_sdf)
            print(f"Parent-only SDF kept {n_kept} records")
            if n_kept == 0:
                print("No parent molecules survived desalting; skipping.")
                continue
            mol2_files = convert_sdf_to_mol2_split_parallel(
                cleaned_sdf,
                ligands_mol2_dir,
                obabel_exe_short,
                threads=obabel_threads,
                timeout_sec=obabel_timeout,
                chunk_size=chunk_size,
            )

        if not mol2_files:
            for _root in (ligands_mol2_dir, output_ligands_dir):
                probe = sorted(_root.glob("*.mol2"))
                if probe:
                    logging.warning(
                        "[compat] Found pre-existing MOL2s under %s (DEPRECATED layout); continuing.",
                        _root,
                    )
                    mol2_files = probe
                    break

        if not mol2_files:
            print("No MOL2 files produced; skipping this SDF.")
            continue

        if only_set:
            stems = {Path(p).stem for p in mol2_files}
            requested = set(only_set)
            mol2_files = [p for p in mol2_files if Path(p).stem in requested]
            sample = sorted(list(requested))[:10]
            remain = len(mol2_files)
            logging.info(
                "[test-mode] enabled count=%d remain=%d sample=%s",
                len(requested),
                remain,
                ",".join(sample),
            )
            print(
                f"[test-mode] enabled count={len(requested)} remain={remain} sample={','.join(sample)}"
            )
            missing = sorted(list(requested - stems))
            if missing:
                logging.warning(
                    "[test-mode] requested ligands not found among MOL2s: %s",
                    ",".join(missing[:20]) + ("..." if len(missing) > 20 else ""),
                )
            if not mol2_files:
                print(
                    "[test-mode] No requested ligands were found. Nothing to do; exiting cleanly."
                )
                _maybe_save(force=True)
                return []

        all_mol2_files.extend(mol2_files)
    return all_mol2_files


def _bulk_prepare_pdbqts(ctx: Dict[str, Any], mol2_files: List[Path]) -> List[Path]:
    """Phase 4: normalize MOL2s, run MGLTools/ADT prep, and update microstate metadata."""
    if not mol2_files:
        return []
    output_ligands_dir = ctx["output_ligands_dir"]
    status_log = ctx["status_log"]
    mgltools_python_short = ctx["mgltools_python_short"]
    prepare_script_short = ctx["prepare_script_short"]
    obabel_exe_short = ctx["obabel_exe_short"]
    is_fda_library = ctx["is_fda_library"]
    microstate_dedup = ctx["microstate_dedup"]
    ph_values = ctx.get("ph_values")
    eff_ph_values = ctx["eff_ph_values"]
    use_ph_subdirs = ctx["use_ph_subdirs"]
    library_out_dir = ctx["library_out_dir"]
    microstate_registry = ctx.get("microstate_registry")
    microstate_index = ctx.get("microstate_index")
    alias_index = ctx.get("alias_index")
    microstates_dir = ctx.get("microstates_dir")
    library_name = ctx["library_name"]
    force = ctx["force"]
    def relative_output(path: Path) -> str:
        return _relative_to_output(path, output_ligands_dir)

    rename_active = ctx.get("rename_active", False)
    rename_prefix = ctx.get("rename_prefix", "")
    rename_pad = int(ctx.get("rename_pad", 5))
    rename_start = int(ctx.get("rename_start", 1))
    rename_force = bool(ctx.get("rename_force", False))

    print(
        f"Preparing {len(mol2_files)} MOL2 files with MGLTools (parallel) into {output_ligands_dir}"
    )
    max_workers = max(1, int(os.environ.get("CPU", "8")))
    futures = []
    future_relpaths: Dict[Any, str] = {}
    future_meta: Dict[concurrent.futures.Future, dict] = {}
    resume_skips = 0
    resume_examples: List[str] = []
    ligprep_debug_seen = 0
    ordered_mol2_files = (
        sorted(mol2_files, key=lambda p: p.name) if rename_active else mol2_files
    )
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        for seq, mol2_file in enumerate(ordered_mol2_files, start=rename_start):
            norm_mol2 = collapse_sanitized_once(Path(mol2_file))
            if norm_mol2 != mol2_file:
                try:
                    norm_mol2.write_bytes(Path(mol2_file).read_bytes())
                    mol2_file = norm_mol2
                    logging.info(
                        "[tidy] normalized double-sanitized -> %s", mol2_file.name
                    )
                except Exception as e:
                    logging.warning("[tidy] unable to normalize %s: %s", mol2_file, e)

            mol2_input = Path(mol2_file)
            rdk_stem = mol2_input.stem
            lig_stem = rdk_stem
            name_path = mol2_input.with_suffix(".name")
            name_exists = name_path.exists()
            parent_name_raw: str = ""
            try:
                if name_exists:
                    parent_name_raw = name_path.read_text(
                        encoding="utf-8", errors="ignore"
                    ).strip()
            except Exception as e:
                logging.warning(
                    "[ligprep] unable to read .name for %s: %s", mol2_input.name, e
                )

            if rename_active:
                lig_stem = f"{rename_prefix}{seq:0{rename_pad}d}"
                if rename_force or mol2_input.stem != lig_stem:
                    try:
                        mol2_input = rename_mol2_with_prefix(mol2_input, lig_stem)
                    except Exception as e:
                        logging.warning(
                            "[ligprep] rename failed for %s: %s", mol2_input, e
                        )
            elif not is_fda_library and parent_name_raw:
                lig_stem = _sanitize_ligand_name_for_filename(parent_name_raw)

            if ligprep_debug_seen < 5:
                print(
                    f"[ligprep-debug] mol2={mol2_input.name} rename_active={rename_active} "
                    f"name_path_exists={name_exists} parent_name_raw={parent_name_raw!r} "
                    f"lig_stem={lig_stem}"
                )
                ligprep_debug_seen += 1

            for ph_value in eff_ph_values:
                ph_label = _ph_label(ph_value) if use_ph_subdirs else None
                ph_out_dir = (
                    output_ligands_dir / ph_label if ph_label else output_ligands_dir
                )
                ph_out_dir.mkdir(parents=True, exist_ok=True)

                if use_ph_subdirs and ph_label:
                    base, _, rest = lig_stem.partition("_")
                    if rest:
                        lig_stem_ph = f"{base}{ph_label}_{rest}"
                    else:
                        lig_stem_ph = f"{lig_stem}{ph_label}"
                    pdbqt_name = f"{lig_stem_ph}.pdbqt"
                else:
                    pdbqt_name = f"{lig_stem}.pdbqt"

                pdbqt_path = ph_out_dir / pdbqt_name

                try:
                    pdbqt_rel = relative_output(pdbqt_path)
                except Exception:
                    pdbqt_rel = pdbqt_path.name

                resume_skip = (
                    not force
                    and pdbqt_path.exists()
                    and pdbqt_path.stat().st_size > 100
                    and is_valid_ligand(pdbqt_path, log_dir=output_ligands_dir)
                )
                if resume_skip:
                    logging.info(
                        "[resume] Valid PDBQT exists, skipping: %s ph=%.2f",
                        pdbqt_path,
                        ph_value,
                    )
                    resume_skips += 1
                    if len(resume_examples) < 10:
                        resume_examples.append(pdbqt_rel)
                    continue

                if force and pdbqt_path.exists():
                    logging.info("[force] Overwriting existing PDBQT: %s", pdbqt_rel)

                logger.debug(
                    "prep_ligands: scheduling ligand prep stem=%s ph=%.2f",
                    lig_stem,
                    ph_value,
                )
                fut = ex.submit(
                    _prepare_one,
                    mgltools_python_short,
                    prepare_script_short,
                    mol2_input,
                    pdbqt_path,
                    obabel_exe_short=obabel_exe_short,
                    status_log_dir=output_ligands_dir,
                    ph=ph_value,
                )
                futures.append(fut)
                future_relpaths[fut] = pdbqt_rel
                future_meta[fut] = {
                    "lig_stem": lig_stem,
                    "ph_label": ph_label or "",
                    "ph_value": float(ph_value),
                    "pdbqt_path": pdbqt_path,
                }

    total = len(futures)
    result_paths: List[Path] = []
    for i, fut in enumerate(as_completed(futures), 1):
        name, status = fut.result()
        lig_stem_fallback = Path(name).stem
        relpath = future_relpaths.get(fut, f"{lig_stem_fallback}.pdbqt")
        result_paths.append(Path(output_ligands_dir) / relpath)
        if status == "ok":
            try:
                _append_prep_status(status_log, lig_stem_fallback, "OK", "", relpath)
            except Exception:
                pass
        else:
            try:
                _append_prep_status(
                    status_log, lig_stem_fallback, "FAIL", status, relpath
                )
            except Exception:
                pass

        meta = future_meta.get(fut) or {}
        lig_stem = meta.get("lig_stem", lig_stem_fallback)
        ph_label = meta.get("ph_label") or ""
        try:
            ph_value = float(meta.get("ph_value", 0.0))
        except Exception:
            ph_value = 0.0
        pdbqt_meta = meta.get("pdbqt_path") or (output_ligands_dir / relpath)
        pdbqt_path = pdbqt_meta if isinstance(pdbqt_meta, Path) else Path(pdbqt_meta)

        if (
            status == "ok"
            and microstate_dedup
            and ph_values is not None
            and microstates_dir is not None
            and microstate_registry is not None
            and microstate_index is not None
            and pdbqt_path.is_file()
        ):
            microstate_id = compute_microstate_id(pdbqt_path)
            if not microstate_id:
                logger.debug(
                    "[microstate] skip_dedup_empty_id ligand=%s ph_label=%s ph_value=%.2f path=%s",
                    lig_stem,
                    ph_label,
                    ph_value,
                    pdbqt_path,
                )
            else:
                microstate_entry = microstate_index.get(microstate_id)

                if microstate_entry is None:
                    canonical_name = f"{lig_stem}__ms_{microstate_id}.pdbqt"
                    canonical_rel_path = f"microstates/{canonical_name}"
                    canonical_path = microstates_dir / canonical_name

                    try:
                        canonical_path.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(pdbqt_path, canonical_path)
                    except Exception as e:
                        logger.warning(
                            "[microstate] copy_to_canonical_failed ligand=%s src=%s dst=%s err=%s",
                            lig_stem,
                            pdbqt_path,
                            canonical_path,
                            e,
                        )

                    microstate_entry = {
                        "microstate_id": microstate_id,
                        "pdbqt_path": canonical_rel_path,
                        "aliases": [],
                    }
                    microstate_registry["microstates"].append(microstate_entry)
                    microstate_index[microstate_id] = microstate_entry
                    ctx["microstate_registry_dirty"] = True

                    logger.debug(
                        "[microstate] new_microstate library=%s ligand=%s microstate_id=%s canonical=%s",
                        library_name,
                        lig_stem,
                        microstate_id,
                        canonical_name,
                    )
                else:
                    canonical_rel_path = microstate_entry.get("pdbqt_path")
                    if canonical_rel_path:
                        canonical_path = library_out_dir / canonical_rel_path
                        if not canonical_path.exists():
                            try:
                                canonical_path.parent.mkdir(parents=True, exist_ok=True)
                                shutil.copy2(pdbqt_path, canonical_path)
                            except Exception as e:
                                logger.warning(
                                    "[microstate] canonical_missing_copy_failed ligand=%s src=%s dst=%s err=%s",
                                    lig_stem,
                                    pdbqt_path,
                                    canonical_path,
                                    e,
                                )

                    logger.debug(
                        "[microstate] reuse_microstate library=%s ligand=%s microstate_id=%s canonical=%s",
                        library_name,
                        lig_stem,
                        microstate_id,
                        canonical_rel_path,
                    )

                alias = {
                    "ligand_stem": lig_stem,
                    "ph_label": ph_label,
                    "ph_value": ph_value,
                }
                if alias not in microstate_entry.get("aliases", []):
                    microstate_entry.setdefault("aliases", []).append(alias)
                    ctx["microstate_registry_dirty"] = True
                    logger.debug(
                        "microstate_dedup: add_alias stem=%s ph_label=%s ph_value=%.2f microstate_id=%s",
                        lig_stem,
                        ph_label,
                        ph_value,
                        microstate_entry["microstate_id"],
                    )

                if alias_index is not None:
                    alias_key = (lig_stem, ph_label, ph_value)
                    alias_index[alias_key] = microstate_entry

        if i % 100 == 0 or status != "ok":
            print(f"[{i}/{total}] {name}: {status}")
    try:
        print(
            f"[ligprep] scheduled={len(futures)} resume_skips={resume_skips} force={force} "
            f"examples_skipped={resume_examples[:5]}"
        )
    except Exception:
        pass
    return result_paths


def prep_ligands_with_mgltools(
    *,
    force: bool = False,
    only: Optional[Set[str]] = None,
    ph_values: Optional[List[float]] = None,
    microstate_dedup: bool = False,
    in_sdf: Optional[Path] = None,
    in_sdf_dir: Optional[Path] = None,
    in_pdb_dir: Optional[Path] = None,
    mol2_dir: Optional[Path] = None,
    out_pdbqt_dir: Optional[Path] = None,
    status_log: Optional[Path] = None,
    root_dir: Optional[Path] = None,
) -> None:
    ctx = _bulk_init_config_and_paths(
        force=force,
        microstate_dedup=microstate_dedup,
        ph_values=ph_values,
        in_sdf=in_sdf,
        in_sdf_dir=in_sdf_dir,
        in_pdb_dir=in_pdb_dir,
        mol2_dir=mol2_dir,
        out_pdbqt_dir=out_pdbqt_dir,
        status_log=status_log,
        root_dir=root_dir,
    )

    unit_sdfs, only_set, test_mode_done = _bulk_select_unit_sdfs(ctx, only)
    if test_mode_done:
        return
    if ctx["in_pdb_dir_env"] and not ctx["in_sdf_env"]:
        result = prep_ligands_from_pdb(
            Path(ctx["in_pdb_dir_env"]).resolve(),
            ctx["ligands_mol2_dir"],
            ctx["output_ligands_dir"],
        )
        ctx["maybe_save_microstate_registry"](force=True)
        return result
    mol2_files = _bulk_generate_mol2_files(ctx, only_set)
    if not mol2_files:
        ctx["maybe_save_microstate_registry"](force=True)
        return
    _bulk_prepare_pdbqts(ctx, mol2_files)
    ctx["maybe_save_microstate_registry"](force=True)

    def _as_path(p) -> Path:
        return p if isinstance(p, Path) else Path(p)

    def _cfg_env_or_default(key: str, default: Optional[str] = None) -> Optional[str]:
        """Resolve config from env first, then shared read_config(), then default."""
        v = os.environ.get(key)
        if v:
            return v
        cfg = read_config()
        val = cfg.get(key)
        if val not in (None, ""):
            return str(val)
        return default

    def _canon_base(output_root: Path, pdb_id: str) -> Path:
        """Canonical per-protein base dir: processed_pdbs/<PDB>"""
        output_root = _as_path(output_root)
        return (output_root / pdb_id.upper()).resolve()

    def _merge_dir(src: Path, dst: Path) -> None:
        """Merge src directory into dst (mkdirs as needed); removes src if emptied."""
        src, dst = src.resolve(), dst.resolve()
        if not src.exists():
            return
        dst.mkdir(parents=True, exist_ok=True)
        for root, dirs, files in os.walk(src):
            r = Path(root)
            rel = r.relative_to(src)
            (dst / rel).mkdir(parents=True, exist_ok=True)
            for d in dirs:
                (dst / rel / d).mkdir(parents=True, exist_ok=True)
            for f in files:
                s = r / f
                t = dst / rel / f
                if t.exists():
                    # prefer keeping existing canonical artifacts; only overwrite if target is missing
                    try:
                        # If same file, skip; else overwrite (safe in our case)
                        if s.stat().st_size == t.stat().st_size:
                            continue
                    except Exception:
                        pass
                shutil.move(str(s), str(t))
        # try to remove empty src tree
        try:
            shutil.rmtree(src)
        except Exception:
            pass

        def fold_legacy_layout(pdb_id: str, output_root) -> None:
            """
            Override with extended handling to pull uppercase legacy dirs into the canonical tree:
              <PDB>_NOLIG            -> processed_pdbs/<PDB>/nolig
              <PDB>_CLEANED_LIGANDS  -> processed_pdbs/<PDB>/ligands_raw
              <PDB>_nolig(.pdb)      -> processed_pdbs/<PDB>/nolig/<PDB>_nolig_phenix_clean.pdb
            """
            try:
                root = _as_path(output_root).resolve()
                pdb_idU = pdb_id.upper()
                base = _canon_base(root, pdb_idU)
                (base / "nolig").mkdir(parents=True, exist_ok=True)
                (base / "ligands_raw").mkdir(parents=True, exist_ok=True)

                # legacy directories (both lower and UPPER)
                legacy_dirs = [
                    (root / f"{pdb_id}_nolig", base / "nolig"),
                    (root / f"{pdb_id.lower()}_nolig", base / "nolig"),
                    (root / f"{pdb_idU}_NOLIG", base / "nolig"),
                    (root / f"{pdb_id}_cleaned_ligands", base / "ligands_raw"),
                    (root / f"{pdb_id.lower()}_cleaned_ligands", base / "ligands_raw"),
                    (root / f"{pdb_idU}_CLEANED_LIGANDS", base / "ligands_raw"),
                ]
                for src, dst in legacy_dirs:
                    if src.exists():
                        logging.info("Migrating legacy directory %s -> %s", src, dst)
                        _merge_dir(src, dst)

                # legacy files
                candidates_files = [
                    (
                        root / f"{pdb_id}_nolig.pdb",
                        base / "nolig" / f"{pdb_idU}_nolig_phenix_clean.pdb",
                    ),
                    (
                        root / f"{pdb_id.lower()}_nolig.pdb",
                        base / "nolig" / f"{pdb_idU}_nolig_phenix_clean.pdb",
                    ),
                ]
                for src, dst in candidates_files:
                    if src.exists() and not dst.exists():
                        logging.info("Moving legacy file %s -> %s", src, dst)
                        dst.parent.mkdir(parents=True, exist_ok=True)
                        shutil.move(str(src), str(dst))
            except Exception as e:
                logging.warning(
                    "fold_legacy_layout (extended) failed for %s: %s", pdb_id, e
                )

        def _expose_ligand_intermediates_for_debug(
            pdb_id: str, ligands_raw: Path
        ) -> None:
            """
            If PREPPED_LIGANDS_DIR is configured, create/update:
              <PREPPED_LIGANDS_DIR>/<PDB>/intermediates -> <processed_pdbs>/<PDB>/ligands_raw
            so intermediates (sanitized PDB, MOL2) are visible next to final PDBQTs.
            """
            try:
                prepped_root = _cfg_env_or_default("PREPPED_LIGANDS_DIR", "")
                if not prepped_root:
                    return
                dst_dir = _as_path(prepped_root) / pdb_id.upper()
                dst_dir.mkdir(parents=True, exist_ok=True)
                link = dst_dir / "intermediates"
                if link.is_symlink() or link.exists():
                    try:
                        if link.is_dir() and not link.is_symlink():
                            shutil.rmtree(link)
                        else:
                            link.unlink()
                    except Exception:
                        pass
                link.symlink_to(ligands_raw.resolve(), target_is_directory=True)
                logging.info("Debug symlink: %s -> %s", link, ligands_raw)
            except Exception as e:
                logging.warning("Could not create debug symlink for %s: %s", pdb_id, e)

        def _element_fix_all_in_dir(ligands_raw: Path) -> None:
            """
            Run your element-column fixer on every PDB in ligands_raw if the function exists.
            This reduces 'Unknown atom name F29 -> C' noise and improves template matches downstream.
            """
            try:
                fixer = globals().get("fix_element_columns_in_file", None)
                if fixer is None:
                    return
                for p in ligands_raw.glob("*.pdb"):
                    try:
                        fixer(p, p)
                    except Exception as e:
                        logging.warning("Element-fix skipped for %s: %s", p.name, e)
            except Exception as e:
                logging.warning("Bulk element-fix failed in %s: %s", ligands_raw, e)

        # Wrap/extend clean_pdb to: (1) fold legacy for this PDB, (2) expose intermediates, (3) element-fix newly extracted ligands.
        if "clean_pdb" in globals():
            _orig_clean_pdb = globals()["clean_pdb"]  # type: ignore[misc]
            def clean_pdb(pdb_file, output_root, logger=None):
                pdb_file = _as_path(pdb_file)
                pdb_id = pdb_file.stem.upper()
                # 1) Make sure any legacy dirs for this PDB are migrated before we proceed
                fold_legacy_layout(pdb_id, output_root)
                # 2) Run the original pipeline
                try:
                    out = _orig_clean_pdb(pdb_file, output_root, logger=logger)
                except TypeError:
                    out = _orig_clean_pdb(pdb_file, output_root)
                fold_legacy_layout(pdb_id, output_root)

                # 3) Expose intermediates via a symlink next to prepped_ligands/<PDB>
                try:
                    # Use the same canonical path helper your code already uses, if present
                    if "canon_paths" in globals():
                        paths = globals()["canon_paths"](pdb_id, output_root)  # type: ignore[index]
                        lig_raw = paths.get("ligands_raw", None)
                        if lig_raw:
                            _element_fix_all_in_dir(lig_raw)
                            _expose_ligand_intermediates_for_debug(pdb_id, lig_raw)
                except Exception as e:
                    logging.warning(
                        "post-clean_pdb expose failed for %s: %s", pdb_id, e
                    )
                return out


__all__ = [
    "prep_ligands_with_mgltools",
    "enumerate_ligands_for_docking",
    "_collect_only_from_env_and_cli",
    "_valid_pdbqt",
    "_log_malformed",
    "is_valid_ligand",
    "_prepare_one",
    "_pdbqt_from_mol2_via_obabel",
    "_audit_protonation_metrics",
    "read_config",
    "compute_microstate_id",
    "load_microstate_registry",
    "save_microstate_registry",
    "get_short_path_name",
    "_run_obabel",
    "_count_explicit_H_in_mol2",
    "_re_aromatize_mol2_in_place",
    "_pdbqt_has_H",
    "_log_elem_fix_summary",
    "standardize_mol_with_activesite",
    "_log_std_diff",
    "_write_aromatic_sdf",
    "_count_aromatic_atoms_in_mol2",
    "_buffer_like_by_counts_from_mol",
    "_matches_counterion",
    "_looks_like_buffer_salt",
    "_is_polyacidic_buffer_like",
    "_buffer_like_by_counts_from_pdbfile",
    "_polyacidic_by_counts_from_pdbfile",
    "MIN_ATOMS_FOR_DOCKING",
    "QUARANTINE_DIRNAME",
    "EXCLUDE_CRYSTAL_ADDITIVES",
    "STANDARD_AMINO_ACIDS",
    "_resolve_obabel_exe",
    "_resolve_prepare_ligand4",
]
