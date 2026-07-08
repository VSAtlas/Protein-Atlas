import logging
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Optional, Tuple, Dict, Set
from collections import Counter


# --- RDKit / Standardization imports ---
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import Crippen
from rdkit.Chem import rdmolops as rdMolOps  # new: shared rdMolOps import

from prep_ligands.prep_ligands_bulk_sdf import (
    LIGPREP_OBABEL_TIMEOUT_SEC,
    _sdf_to_mol2,
    convert_sdf_to_mol2_split_parallel,
)
from prep_ligands.prep_ligands_bulk_meeko import prep_ligands_with_meeko
from prep_ligands.prep_ligands_bulk_mol2 import bulk_generate_mol2_files
from prep_ligands.prep_ligands_bulk_prepare import bulk_prepare_pdbqts
from prep_ligands.prep_ligands_microstates import (
    _collect_only_from_env_and_cli,
    compute_microstate_id,
    enumerate_ligands_for_docking,
    load_microstate_registry,
    save_microstate_registry,
)
from prep_ligands.prep_ligands_reporting import (
    _audit_protonation_metrics,
    _log_elem_fix_summary,
    _log_malformed,
)
from prep_ligands.prep_ligands_common import (
    EXCLUDE_CRYSTAL_ADDITIVES,
    MIN_ATOMS_FOR_DOCKING,
    QUARANTINE_DIRNAME,
    STANDARD_AMINO_ACIDS,
    _HAS_STD,
    _buffer_like_by_counts_from_mol,
    _buffer_like_by_counts_from_pdbfile,
    _count_aromatic_atoms_in_mol2,
    _count_explicit_H_in_mol2,
    _is_polyacidic_buffer_like,
    _log_std_diff,
    _matches_counterion,
    _pdbqt_from_mol2_via_obabel,
    _pdbqt_has_H,
    _polyacidic_by_counts_from_pdbfile,
    _prepare_one,
    _re_aromatize_mol2_in_place,
    _resolve_obabel_exe,
    _resolve_prepare_ligand4,
    _run_obabel,
    _std,
    _looks_like_buffer_salt,
    _write_sdf_for_obabel,
    get_short_path_name,
    is_valid_ligand,
    read_config,
    standardize_mol_with_activesite,
    _quick_filters,
    _standardize_then_sanitize,
)
from prep_ligands.prep_ligands_runtime import (
    BulkContext,
)

logger = logging.getLogger(__name__)

# --- salt remover (fallback path also uses this) ---
_REM = SaltRemover()
MIN_PARENT_HEAVY = 4


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


def _ph_label(ph: float) -> str:
    """Canonical label for ligand pH subdirectories, e.g. 7.0 -> 'pH7_0'."""
    return f"pH{ph:.1f}".replace(".", "_")


# =========================
# Parent chooser
# =========================


def _to_parent_mol(m: Chem.Mol) -> Optional[Chem.Mol]:
    try:
        if _HAS_STD and _std is not None:
            try:
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
            except Exception as exc:
                logging.debug(
                    "[ligprep] rdMolStandardize parent chooser failed; falling back to salt remover (%s)",
                    exc,
                )
                p = None
        else:
            p = None

        if p is None:
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
    from pkasolver import pkasolver as _pka

    _HAS_PKASOLVER = True
except ImportError:
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
    return float(getattr(Crippen, "MolLogP")(mol))


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
        clogp = float(getattr(Crippen, "MolLogP")(mol))
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
            if _write_sdf_for_obabel(m, out):
                logging.info("[paths] sdf_out=%s", str(out.resolve()))
                return out

            else:
                _log_malformed(out, "write_sdf_for_obabel_failed")
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
        rdk_idx = _stem_to_rdk_index(pth.stem)
        if rdk_idx is not None and 0 <= rdk_idx < len(parent_names):
            parent_name = parent_names[rdk_idx]
        name_path = out.with_suffix(".name")
        try:
            to_write = (parent_name or "").rstrip("\r\n")
            name_path.write_text(to_write + "\n", encoding="utf-8")
        except Exception as e:
            logging.warning("[bulkSDF] name_write_failed stem=%s err=%s", pth.stem, e)
        if rdk_idx is not None and rdk_idx < debug_limit:
            exists = name_path.exists()
            try:
                contents = name_path.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                contents = "<unreadable>"
            print(
                f"[rdkit-debug] mol2={out.name} name_path={name_path.name} "
                f"exists={exists} contents={contents!r}"
            )
        if rdk_idx is not None and rdk_idx < debug_limit:
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


def _bulk_select_unit_sdfs(
    ctx: BulkContext, only_set: Optional[Set[str]]
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
    ctx: BulkContext, only_set: Optional[Set[str]]
) -> List[Path]:
    """Phase 3: convert bulk SDFs (or explicit overrides) into MOL2 files."""
    return bulk_generate_mol2_files(
        ctx,
        only_set,
        rdkit_embed_fn=rdkit_embed_sdf_to_mol2,
        write_parent_sdf_fn=_write_parent_only_sdf,
        convert_parallel_fn=convert_sdf_to_mol2_split_parallel,
    )


def _bulk_prepare_pdbqts(ctx: BulkContext, mol2_files: List[Path]) -> List[Path]:
    """Phase 4: normalize MOL2s, run MGLTools/ADT prep, and update microstate metadata."""
    return bulk_prepare_pdbqts(ctx, mol2_files)


prep_ligands_with_mgltools = prep_ligands_with_meeko


__all__ = [
    "prep_ligands_with_meeko",
    "prep_ligands_with_mgltools",
    "enumerate_ligands_for_docking",
    "_collect_only_from_env_and_cli",
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
    "_write_sdf_for_obabel",
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
