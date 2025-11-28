# dedup_sdfs.py
from __future__ import annotations
import csv, sys, re
from pathlib import Path
from typing import Iterable, List, Dict, Tuple, Any
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
# --- robust imports for standardization ---
Standardize = None
try:
    from rdkit.Chem import rdMolStandardize as Standardize  # type: ignore
except Exception:
    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize as Standardize  # type: ignore
    except Exception:
        Standardize = None
from rdkit.Chem.SaltRemover import SaltRemover
from input_and_export_functions import load_config, validate_config
_cfg = load_config("config.txt") or {}
try:
    validate_config(_cfg)
except Exception:
    pass

_ID_PREFIX = str(_cfg.get("DEDUP_ID_PREFIX", "fda_"))
_ID_WIDTH  = int(_cfg.get("DEDUP_ID_WIDTH", 4))
_FDA_RE = re.compile(rf"^{re.escape(_ID_PREFIX)}\d{{{_ID_WIDTH}}}$", re.IGNORECASE)




def _iter_mols(paths: Iterable[Path]):
    """Yield RDKit molecules from one or more SDF files with provenance props set."""
    for p in paths:
        suppl = Chem.SDMolSupplier(str(p), removeHs=False, sanitize=True)
        for i, m in enumerate(suppl):
            if m is None:
                continue
            m.SetProp("_source_file", str(p))
            m.SetProp("_source_idx", str(i))
            if not m.HasProp("_Name") or not m.GetProp("_Name").strip():
                m.SetProp("_Name", f"unnamed_{i}")
            yield m


def _parent_neutralize(m: Chem.Mol) -> Chem.Mol:
    """Return desalted, cleaned, neutralized, canonical-tautomer parent (when possible)."""
    if Standardize:
        lf = Standardize.LargestFragmentChooser()
        m = lf.choose(m)
        m = Standardize.Cleanup(m)
        try:
            m = Standardize.Uncharger().uncharge(m)
        except Exception:
            pass
        try:
            m = Standardize.TautomerEnumerator().Canonicalize(m)
        except Exception:
            pass
        return m
    else:
        # Minimal fallback: remove common salts and keep largest fragment
        remover = SaltRemover(defnData=SaltRemover.SaltRemoverParams().defnData)
        m = remover.StripMol(m, dontRemoveEverything=True)
        frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True)
        m = max(frags, key=lambda x: x.GetNumHeavyAtoms()) if frags else m
        return m


def _make_keys(m_parent: Chem.Mol) -> Tuple[str, str]:
    """Return (full InChIKey or SMILES, parent14)."""
    try:
        ik = rdMolDescriptors.InchiKey(m_parent)
    except Exception:
        ik = Chem.MolToSmiles(m_parent, isomericSmiles=True)
    parent14 = ik.split("-")[0] if "-" in ik else ik[:14]
    return ik, parent14


def _choose_rep(existing: Chem.Mol, candidate: Chem.Mol) -> Chem.Mol:
    """Prefer higher heavy-atom count; then 3D conformer."""
    def score(m: Chem.Mol):
        ha = m.GetNumHeavyAtoms()
        has3d = False
        try:
            has3d = m.GetConformer().Is3D()
        except Exception:
            pass
        return (ha, 1 if has3d else 0)

    return candidate if score(candidate) > score(existing) else existing


def _pick_rep_name(members: List[str], fallback_name: str) -> str:
    """If any member already has a minted ID, keep that exact name as representative."""
    for nm in members:
        if _FDA_RE.match(nm):
            return nm
    return fallback_name



def dedup_sdfs(
    input_files: List[str | Path],
    out_file: str | Path,
    map_file: str | Path,
    dups_file: str | Path,
    id_seed: int = 1,
) -> Dict[str, Any]:
    """
    De-duplicate molecules across one or more SDFs, choosing a representative per exact structure,
    preserving existing 'fda_####' names when present, otherwise minting new IDs.

    Returns a small summary dict.
    """
    # Normalize paths
    inputs = [Path(p) for p in input_files]
    out_file = Path(out_file)
    map_file = Path(map_file)
    dups_file = Path(dups_file)

    reps: Dict[str, Chem.Mol] = {}          # key_full -> representative RDKit Mol
    groups: Dict[str, List[Dict[str, str]]] = {}   # key_full -> list of member metadata
    names_by_key: Dict[str, List[str]] = {}  # key_full -> original names

    # Build groups and choose representatives
    for m in _iter_mols(inputs):
        orig_name = m.GetProp("_Name")
        src = m.GetProp("_source_file")
        idx = m.GetProp("_source_idx")

        parent = _parent_neutralize(m)
        key_full, key_parent14 = _make_keys(parent)

        groups.setdefault(key_full, []).append(
            {"orig_name": orig_name, "source_file": src, "source_idx": idx, "parent14": key_parent14}
        )
        names_by_key.setdefault(key_full, []).append(orig_name)

        if key_full in reps:
            reps[key_full] = _choose_rep(reps[key_full], m)
        else:
            reps[key_full] = m

    # Assign representative names: keep existing fda_#### else mint
    next_id = int(id_seed)
    key_to_repname: Dict[str, str] = {}
    used_ids = set()
    # First pass: reserve any pre-existing minted IDs
    for key, rep in reps.items():
        rep_name = _pick_rep_name(names_by_key[key], rep.GetProp("_Name"))
        if _FDA_RE.match(rep_name):
            rep_name = rep_name.lower()
            key_to_repname[key] = rep_name
            used_ids.add(rep_name)

    def mint() -> str:
        nonlocal next_id
        while True:
            cand = f"{_ID_PREFIX}{next_id:0{_ID_WIDTH}d}"
            next_id += 1
            if cand not in used_ids:
                used_ids.add(cand)
                return cand

    for key in reps:
        if key not in key_to_repname:
            key_to_repname[key] = mint()

    # Write merged SDF with representative names & keys
    out_file.parent.mkdir(parents=True, exist_ok=True)
    w = Chem.SDWriter(str(out_file))
    for key, mol in reps.items():
        repname = key_to_repname[key]
        mol.SetProp("_Name", repname)
        mol.SetProp("DEDUPE_KEY", key)
        p14 = key.split("-")[0] if "-" in key else key[:14]
        mol.SetProp("PARENT14", p14)
        mol.SetProp("ORIGINAL_NAMES", ";".join(sorted(set(names_by_key[key]))))
        w.write(mol)
    w.close()

    # Write mapping (all originals -> chosen rep)
    map_file.parent.mkdir(parents=True, exist_ok=True)
    with map_file.open("w", newline="") as f:
        wcsv = csv.writer(f)
        wcsv.writerow(
            ["dedupe_key_full", "parent14", "rep_fda_id", "member_name", "member_source_file", "member_source_idx"]
        )
        for key, members in groups.items():
            rep_fda = key_to_repname[key]
            for m in members:
                wcsv.writerow([key, m["parent14"], rep_fda, m["orig_name"], m["source_file"], m["source_idx"]])

    # Parent-level variant report
    parent_groups: Dict[str, set] = {}
    for key in groups:
        p14 = key.split("-")[0] if "-" in key else key[:14]
        parent_groups.setdefault(p14, set()).add(key)

    with dups_file.open("w", newline=True) as f:
        wcsv = csv.writer(f)
        wcsv.writerow(["parent14", "num_exact_variants", "example_full_keys"])
        for p14, fullset in sorted(parent_groups.items(), key=lambda x: -len(x[1])):
            wcsv.writerow([p14, len(fullset), ";".join(list(fullset)[:5])])

    if not Standardize:
        print(
            "WARNING: rdMolStandardize not found. Used fallback (SaltRemover + largest fragment). "
            "Install a newer RDKit for full Cleanup/Uncharger/Tautomer canonicalization.",
            file=sys.stderr,
        )

    return {
        "n_unique_exact": len(reps),
        "n_parent_groups": len(parent_groups),
        "last_id_used": max(int(x.split("_")[1]) for x in used_ids if _FDA_RE.match(x)) if used_ids else None,
        "out_file": str(out_file),
        "map_file": str(map_file),
        "dups_file": str(dups_file),
    }


# Optional CLI wrapper so this file also works as a script if you ever need it.
def _main_cli():
    import argparse
    ap = argparse.ArgumentParser(description="De-duplicate SDFs and mint/keep IDs.")
    ap.add_argument("--in",  dest="inputs",   nargs="+", required=True)
    ap.add_argument("--out", dest="out_file", required=True)
    ap.add_argument("--map", dest="map_file", required=True)
    ap.add_argument("--dups",dest="dups_file",required=True)
    # [ANCHOR CLI default from cfg; user-provided still wins]
    ap.add_argument("--id-seed", type=int, default=int(_cfg.get("DEDUP_ID_SEED", 1)))
    args = ap.parse_args()
    summary = dedup_sdfs(
        input_files=args.inputs,
        out_file=args.out_file,
        map_file=args.map_file,
        dups_file=args.dups_file,
        id_seed=args.id_seed,
    )
    print(f"Done. Summary: {summary}")

if __name__ == "__main__":
    _main_cli()
