"""Ligand extraction and debug-intermediate helpers."""

from __future__ import annotations

import json
import logging
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import List, Union

from protein_prep.aliases_policy import (
    ALIASES,
    _RETAIN_VARIANT,
    _RETAIN_VARIANT_CANONICAL,
    _WATER_NAMES,
    _normalize_resname,
)

try:
    from rdkit import Chem  # type: ignore

    _HAS_RDKIT = True
except Exception:
    _HAS_RDKIT = False


def _hydrate_legacy_globals() -> None:
    import automate_protein_prep as legacy

    g = globals()
    for name, value in legacy.__dict__.items():
        g.setdefault(name, value)


def collapse_sanitized_once(p: Union[str, Path]) -> Path:
    """
    Return a Path with a sanitized basename (single pass; no filesystem changes).
    Keeps directory the same, replaces non [A-Za-z0-9._-] with underscores.
    """
    _hydrate_legacy_globals()
    p = Path(p)
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", p.name)
    return p.with_name(safe)


def _write_pristine_reference(pdb_lig_path: Path) -> None:
    """Write pristine ligand copy next to ligands_raw (SDF + light atom map JSON)."""
    _hydrate_legacy_globals()
    ref_dir = pdb_lig_path.parent.parent / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)
    ref_sdf = ref_dir / (pdb_lig_path.stem + ".sdf")
    ref_json = ref_dir / (pdb_lig_path.stem + ".map.json")

    if _HAS_RDKIT:
        m = Chem.MolFromPDBFile(str(pdb_lig_path), sanitize=False, removeHs=False)
        if m is not None:
            w = Chem.SDWriter(str(ref_sdf))
            w.write(m)
            w.close()
            amap = []
            with open(pdb_lig_path, "r", encoding="utf-8", errors="ignore") as f:
                for ln in f:
                    if ln.startswith(("ATOM", "HETATM")):
                        amap.append(
                            {
                                "name": ln[12:16].strip(),
                                "element": (
                                    ln[76:78].strip() or ln[12:16].strip()[:1].upper()
                                ),
                            }
                        )
            json.dump({"atoms": amap}, open(ref_json, "w"), indent=2)
            return

    # Fallback minimal SDF stub
    with open(ref_sdf, "w", encoding="utf-8") as out:
        out.write(f"{pdb_lig_path.stem}\n  -Pristine-\n\n")
        out.write("$$$$\n")


def extract_ligands_from_filtered(
    filtered_pdb: Union[str, Path], out_dir: Union[str, Path]
) -> List[Path]:
    """Extract non-water HETATM residues into individual PDBs (RES_CHAINRESI.pdb),
    and run text-level element repair on each (YAML rules)."""
    _hydrate_legacy_globals()
    outd = Path(out_dir)
    outd.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    cur_key = None
    bucket: List[str] = []

    def flush():
        nonlocal bucket, cur_key, written
        if not bucket or cur_key is None:
            return
        resname, chain, resseq, icode = cur_key
        name = f"{resname}_{chain}{int(resseq)}"
        outp = outd / f"{name}.pdb"
        with open(outp, "w") as w:
            w.write(f"REMARK Extracted {name}\n")
            for ln in bucket:
                w.write(ln)
            w.write("TER\nEND\n")

        # normalize once (no-op if not needed)
        norm = collapse_sanitized_once(outp)
        if norm.name != outp.name:
            try:
                outp.replace(norm)
                outp = norm
            except Exception:
                pass

        # YAML-backed element repair (PDB only)
        try:
            fix_element_columns_in_file(outp, outp)
        except Exception as e:
            logging.warning("Element-fix skipped for %s: %s", outp.name, e)
        written.append(outp)
        bucket = []
        cur_key = None
        try:
            _write_pristine_reference(outp)
        except Exception as e:
            logging.warning(
                "Could not write pristine reference for %s: %s", outp.name, e
            )

    lines: List[str] = []
    try:
        with open(filtered_pdb, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
    except Exception as exc:
        logging.warning("[extract.check] file=%s err=%s", filtered_pdb, exc)
        lines = []

    candidate_resnames = {
        ln[17:20].strip().upper()
        for ln in lines
        if ln.startswith("HETATM") and ln[17:20].strip()
    }
    logging.info(
        "[extract.check] file=%s candidate_resnames=%s",
        filtered_pdb,
        ",".join(sorted(candidate_resnames)) if candidate_resnames else "none",
    )
    retain_allow = {tok.upper() for tok in getattr(ALIASES, "retain_resnames", [])}
    retain_allow_canonical = {
        _normalize_resname(tok)
        for tok in getattr(ALIASES, "retain_resnames", [])
        if _normalize_resname(tok)
    }
    overlap = sorted(
        res
        for res in candidate_resnames
        if (res in retain_allow) or (_normalize_resname(res) in retain_allow_canonical)
    )
    if overlap:
        logging.warning(
            "[extract.violation] file=%s will_extract_retain_list=%s",
            filtered_pdb,
            ",".join(overlap),
        )

    for ln in lines:
        if not ln.startswith("HETATM"):
            continue
        resname = ln[17:20].strip().upper()
        if resname in _WATER_NAMES:
            continue  # never extract waters
        canonical_res = _normalize_resname(resname)
        if (resname in _RETAIN_VARIANT) or (canonical_res in _RETAIN_VARIANT_CANONICAL):
            continue  # don't extract cofactors/metals/ions you keep with protein
        chain = ln[21]
        resseq = ln[22:26].strip() or "0"
        icode = ln[26]
        key = (resname, chain, resseq, icode)
        if cur_key is None:
            cur_key = key
        if key != cur_key:
            flush()
            cur_key = key
        bucket.append(ln)
    flush()
    logging.info("Extracted %d ligand residues to %s", len(written), out_dir)
    return written


def element_fix_all_in_dir(
    dir_path: Union[str, Path], rewrite_atoms: bool = False
) -> int:
    """
    Run the text-level element column fixer on every *.pdb under dir_path.
    Returns the number of files rewritten.
    """
    _hydrate_legacy_globals()
    d = Path(dir_path)
    if not d.exists():
        return 0
    n = 0
    for p in d.rglob("*.pdb"):
        try:
            fix_element_columns_in_file(p, dst_path=p, rewrite_atoms=rewrite_atoms)
            n += 1
        except Exception as e:
            logging.warning("element_fix_all_in_dir skipped %s: %s", p, e)
    logging.info("Element column sweep fixed %d PDB files under %s", n, d)
    return n


def expose_ligand_intermediates_for_debug(
    src_dir: Union[str, Path], link_dir: Union[str, Path]
) -> None:
    """
    Create/refresh a symlink 'link_dir' -> 'src_dir' for easy browsing.
    On Windows, tries directory junction fallback if symlink fails.
    Quietly skips if the source doesn't exist yet.
    """
    _hydrate_legacy_globals()
    src = Path(src_dir).resolve()
    dst = Path(link_dir)

    # Avoid noise if source not present yet during early prep
    if not src.exists():
        logging.debug("intermediates: skip symlink, source missing: %s", src)
        return

    # Ensure link's parent exists
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logging.warning("intermediates: could not ensure parent dir for %s: %s", dst, e)
        return

    try:
        # Remove existing link/dir
        if dst.is_symlink() or dst.exists():
            try:
                if dst.is_symlink():
                    dst.unlink()
                else:
                    shutil.rmtree(dst)
            except Exception:
                pass

        # Create link (with Windows junction fallback)
        if os.name == "nt":
            try:
                os.symlink(src, dst, target_is_directory=True)
            except OSError:
                cmd = ["cmd", "/c", "mklink", "/J", str(dst), str(src)]
                subprocess.run(cmd, check=True)
        else:
            os.symlink(src, dst, target_is_directory=True)

        logging.info("Exposed ligand intermediates: %s -> %s", dst, src)
    except Exception as e:
        logging.warning(
            "intermediates: could not create symlink %s -> %s: %s", dst, src, e
        )
