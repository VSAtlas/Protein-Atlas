from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

from path_router.path_router import make_paths  # type: ignore[import-not-found]
from prep_ligands.library_index import LibraryIndex  # type: ignore[import-not-found]

import analysis.dud_eval_types as dud_eval_types
from analysis.dud_eval_log import dbg
from analysis.dud_eval_types import TargetSpec

_LIB_INDEX_CACHE: Dict[Tuple[str, str], LibraryIndex] = {}


def _get_library_index(root: Path, manifest_filename: str) -> LibraryIndex:
    """
    Return a cached LibraryIndex for a given prepped-ligands root + manifest filename.
    This avoids re-reading the same manifest for every target.
    """
    key = (str(root.resolve()), manifest_filename)
    idx = _LIB_INDEX_CACHE.get(key)
    if idx is None:
        idx = LibraryIndex(manifest_filename=manifest_filename)
        idx.load([root])
        _LIB_INDEX_CACHE[key] = idx
    return idx


def _is_probable_run_id_dirname(name: str) -> bool:
    """
    Heuristic to detect run-id style directory names (e.g., 20251205_225826, 2025-12-18T...).
    Avoid matching classic 4-char PDB IDs.
    """
    if not name or len(name) <= 6:
        return False
    upper = name.upper()
    if len(upper) == 4 and upper.isalnum():
        return False
    if "_" in name:
        return True
    if re.fullmatch(r"\d{8,}", name):
        return True
    if "T" in name and re.match(r"\d{4}-?\d{2}-?\d{2}", name):
        return True
    return False


def _resolve_scan_roots(docked_root: Path, run_id: Optional[str]) -> List[Path]:
    """
    Determine which roots to scan for docking outputs.

    Priority:
      - If run_id is provided and docked_root/run_id exists, scan that first,
        then fall back to docked_root for legacy layouts.
      - If run_id provided but missing, fall back to docked_root.
      - If no run_id: scan run-like subdirs (mtime desc) if present, then docked_root.
      - If docked_root is not a dir, just return [docked_root].
    """
    if not docked_root.is_dir():
        return [docked_root]

    if run_id:
        candidate = docked_root / run_id
        if candidate.is_dir():
            roots: List[Path] = [candidate]
            if docked_root != candidate:
                roots.append(docked_root)
            return roots
        return [docked_root]

    run_subdirs = [
        p
        for p in docked_root.iterdir()
        if p.is_dir() and _is_probable_run_id_dirname(p.name)
    ]
    if run_subdirs:
        run_subdirs.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        return run_subdirs + [docked_root]

    return [docked_root]


def make_target_key(pdb_id: str, variant: Optional[str], ph_tag: Optional[str]) -> str:
    """
    Build a path-safe identifier for a target.
    Backward compatibility: if both variant and pH are blank, return the pdb_id unchanged.
    """
    variant_norm = (variant or "").strip()
    ph_norm = (ph_tag or "").strip()
    if not variant_norm and not ph_norm:
        return pdb_id
    variant_norm = variant_norm.upper().replace(" ", "_")
    ph_norm = ph_norm.replace(" ", "_")
    for sep in (os.sep, os.altsep):
        if sep:
            variant_norm = variant_norm.replace(sep, "_")
            ph_norm = ph_norm.replace(sep, "_")
    return f"{pdb_id}__{variant_norm}__{ph_norm}"


def _dedup(seq: Iterable[str]) -> List[str]:
    """Preserve order while removing duplicates."""
    seen: Set[str] = set()
    out: List[str] = []
    for item in seq:
        if item and item not in seen:
            seen.add(item)
            out.append(item)
    return out


def _normalize_pdb_ids(tokens: Optional[Iterable[str]]) -> List[str]:
    """
    Normalize tokens to uppercase 4-character PDB IDs.
    Strips extensions and non-alphanumerics; drops too-short tokens.
    """
    if tokens is None:
        return []
    normalized: List[str] = []
    for raw in tokens:
        if raw is None:
            continue
        text = str(raw).strip()
        if not text:
            continue
        text = re.sub(r"\.pdb(?:\.gz)?$", "", text, flags=re.IGNORECASE)
        alnum = re.sub(r"[^A-Za-z0-9]", "", text).upper()
        if len(alnum) < 4:
            dbg("WARN", "pdb-filter", f"token={text} reason=too_short")
            continue
        normalized.append(alnum[:4])
    return _dedup(normalized)


def _is_under(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except Exception:
        return False


def _resolve_reranked_scorch_path(
    spec: TargetSpec,
    docked_root: Path,
    post_docked_root: Path,
    run_id: Optional[str],
    run_dir: Optional[Path],
) -> Optional[Path]:
    """
    Best-effort resolution of consensus_reranked_scorch.csv for a target.
    Prefers new layout (post_docked/<run_id>/...) but falls back to legacy.
    """
    candidates: List[Path] = []
    seen: Set[str] = set()

    if spec.csv_path:
        try:
            rel = spec.csv_path.parent.relative_to(docked_root)
            cand = post_docked_root / rel / dud_eval_types.RERANKED_SCORCH_BASENAME
            if str(cand) not in seen:
                candidates.append(cand)
                seen.add(str(cand))
        except Exception:
            pass

    if run_dir and spec.csv_path and _is_under(spec.csv_path, run_dir):
        try:
            rel_new = spec.csv_path.parent.relative_to(run_dir)
            base = post_docked_root / (run_id or "")
            cand = base / rel_new / dud_eval_types.RERANKED_SCORCH_BASENAME
            if str(cand) not in seen:
                candidates.append(cand)
                seen.add(str(cand))
        except Exception:
            pass

    if run_id:
        base = post_docked_root / run_id
        combos = [
            base
            / spec.pdb_id
            / (spec.variant or "")
            / (spec.ph_tag or "")
            / dud_eval_types.RERANKED_SCORCH_BASENAME,
            base
            / spec.pdb_id
            / (spec.variant or "")
            / dud_eval_types.RERANKED_SCORCH_BASENAME,
            base / spec.pdb_id / dud_eval_types.RERANKED_SCORCH_BASENAME,
        ]
        for cand in combos:
            if str(cand) not in seen:
                candidates.append(cand)
                seen.add(str(cand))

    legacy_base = post_docked_root
    legacy_candidates = [
        legacy_base
        / spec.pdb_id
        / (spec.variant or "")
        / (spec.ph_tag or "")
        / dud_eval_types.RERANKED_SCORCH_BASENAME,
        legacy_base
        / spec.pdb_id
        / (spec.variant or "")
        / dud_eval_types.RERANKED_SCORCH_BASENAME,
        legacy_base / spec.pdb_id / dud_eval_types.RERANKED_SCORCH_BASENAME,
    ]
    for cand in legacy_candidates:
        if str(cand) not in seen:
            candidates.append(cand)
            seen.add(str(cand))

    for cand in candidates:
        if cand.exists() and cand.is_file():
            return cand
    return None


def extract_compnd_molecules(pdb_lines: List[str]) -> List[str]:
    """Pull COMPND.MOLECULE entries from a PDB header."""
    molecules: List[str] = []
    pending: Optional[str] = None
    last_key: Optional[str] = None

    def _flush() -> None:
        nonlocal pending
        if pending:
            molecules.append(pending.strip().rstrip(","))
            pending = None

    for line in pdb_lines:
        if not line.startswith("COMPND"):
            continue
        payload = line[10:].strip()
        if not payload:
            continue
        segments = [seg.strip() for seg in payload.split(";")]
        for seg in segments:
            if not seg:
                continue
            if ":" in seg:
                key, val = seg.split(":", 1)
                key = key.strip().upper()
                val = val.strip()
                if key != last_key and last_key == "MOLECULE":
                    _flush()
                last_key = key
                if key == "MOLECULE":
                    pending = f"{pending} {val}".strip() if pending else val
                elif pending and key != "MOLECULE":
                    _flush()
            else:
                if last_key == "MOLECULE":
                    pending = f"{pending} {seg}".strip() if pending else seg
    if last_key != "MOLECULE":
        _flush()
    _flush()
    return _dedup(molecules)


def extract_uniprot_from_dbref(pdb_lines: List[str]) -> Tuple[List[str], List[str]]:
    """Return (entry_names, accessions) discovered from DBREF UNP rows."""
    entry_names: List[str] = []
    accessions: List[str] = []
    for line in pdb_lines:
        if not line.startswith("DBREF"):
            continue
        tokens = line.split()
        if len(tokens) < 7:
            continue
        db = tokens[5].upper()
        if db not in {"UNP", "UNIPROT"}:
            continue
        accession = tokens[6].strip()
        entry = tokens[7].strip() if len(tokens) >= 8 else ""
        if accession:
            accessions.append(accession)
        if entry:
            entry_names.append(entry)
    return _dedup(entry_names), _dedup(accessions)


def choose_target_name(
    compnd_mols: List[str],
    uniprot_entries: List[str],
    uniprot_accessions: List[str],
    prefer: str = "auto",
) -> str:
    """Pick the best target label respecting preference order."""
    prefer_key = (prefer or "auto").lower()
    compnd_choice = ", ".join(compnd_mols) if compnd_mols else ""
    uniprot_entry = uniprot_entries[0] if uniprot_entries else ""
    uniprot_acc = uniprot_accessions[0] if uniprot_accessions else ""

    if prefer_key == "compnd":
        return compnd_choice or uniprot_entry or uniprot_acc
    if prefer_key == "uniprot":
        return uniprot_entry or uniprot_acc or compnd_choice
    # auto
    return compnd_choice or uniprot_entry or uniprot_acc


def _read_pdb_header_lines(pdb_path: Path) -> List[str]:
    lines: List[str] = []
    try:
        with pdb_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith(("ATOM", "HETATM")):
                    break
                lines.append(line.rstrip("\n"))
    except Exception:
        return []
    return lines


def _infer_library_name_via_manifest(
    target_id: str,
    ligand_basenames: Set[str],
    roots: List[Path],
    cfg: Dict,
) -> str:
    """
    Fast path: infer library name using a pre-built manifest under prepped_ligands.

    Strategy:
      - For each existing root that has a manifest file (default: _manifest.json),
        load a LibraryIndex for that root.
      - For a bounded sample of ligand basenames, look up each basename in the
        manifest via LibraryIndex.lookup_filename.
      - Each hit yields a relative path like "abl1/ABL1_active_10.pdbqt";
        treat the first path component ("abl1") as the library name and count votes.
      - If there is a unique best library with >0 votes, return it.
      - Otherwise, return "" so the caller can fall back to the legacy scan.
    """
    if not ligand_basenames:
        return ""

    manifest_filename = cfg.get("LIBRARY_MANIFEST_FILENAME", "_manifest.json")
    if not manifest_filename:
        manifest_filename = "_manifest.json"
    manifest_filename = str(manifest_filename)

    sample_names = sorted(ligand_basenames)
    max_sample = 200
    if len(sample_names) > max_sample:
        sample_names = sample_names[:max_sample]

    overall_counts: Dict[str, int] = {}

    for root in roots:
        if not root or not root.exists():
            continue
        manifest_path = root / manifest_filename
        if not manifest_path.exists():
            dbg("DEBUG", "library", f"pdb={target_id} manifest_missing={manifest_path}")
            continue

        try:
            idx = _get_library_index(root, manifest_filename)
        except Exception as exc:
            dbg(
                "WARN",
                "library",
                f"pdb={target_id} manifest_load_failed root={root} err={exc}",
            )
            continue

        for name in sample_names:
            base = os.path.basename(name)
            if not base:
                continue
            hit = idx.lookup_filename(base, roots=[root])
            if not hit:
                continue
            try:
                rel = hit.relative_to(root)
            except Exception:
                rel = hit
            parts = rel.as_posix().split("/")
            if not parts:
                continue
            lib_name = parts[0].strip()
            if not lib_name:
                continue
            overall_counts[lib_name] = overall_counts.get(lib_name, 0) + 1

    if not overall_counts:
        return ""

    best_lib, best_count = max(overall_counts.items(), key=lambda kv: (kv[1], kv[0]))
    n_best = sum(1 for c in overall_counts.values() if c == best_count)
    if best_count <= 0 or n_best != 1:
        dbg(
            "WARN",
            "library",
            f"pdb={target_id} manifest_votes ambiguous best_count={best_count} "
            f"candidates={len(overall_counts)}",
        )
        return ""

    dbg(
        "INFO",
        "library",
        f"pdb={target_id} picked='{best_lib}' method='manifest_entries' "
        f"matches={best_count}",
    )
    return best_lib


def infer_library_name(
    target_id: str,
    ligand_basenames: Set[str],
    *,
    prepped_root_override: Optional[Path],
    cfg: Dict,
) -> str:
    """
    Infer library folder name when libraries live at:
        <prepped_root>/<library>/*files*
    No per-PDB subdirectory is expected.

    Preference:
      1) Global prepped-ligands manifest via LibraryIndex (fast path).
      2) manifest.json at library root with "library_name"/"name"/"library"/"label"
      3) filename-overlap between docked basenames and files in <library> (non-recursive; optional shallow)
    """

    # --- gather candidate roots from CLI override + config
    candidate_roots: List[Path] = []
    if prepped_root_override:
        candidate_roots.append(prepped_root_override)
    if cfg:
        try:
            paths = make_paths(cfg, base_id=target_id, pdb_file=f"{target_id}.pdb")
            candidate_roots.append(paths.prepped_root)
        except Exception:
            pass

    roots: List[Path] = []
    seen: Set[Path] = set()
    for r in candidate_roots:
        if r and r not in seen:
            seen.add(r)
            roots.append(r)

    # DEBUG: show roots + docked basename sample
    if roots:
        dbg(
            "DEBUG",
            "library",
            f"pdb={target_id} search_roots={';'.join(str(r) for r in roots)}",
        )
    else:
        dbg(
            "WARN",
            "library",
            f"pdb={target_id} no candidate roots for library inference",
        )

    if ligand_basenames:
        dbg(
            "DEBUG",
            "library",
            f"pdb={target_id} docked_basenames_n={len(ligand_basenames)} "
            f"sample={list(sorted(ligand_basenames))[:3]}",
        )
    else:
        dbg("DEBUG", "library", f"pdb={target_id} docked_basenames_n=0")

    # 1) Fast path: try prepped-ligands manifest via LibraryIndex (if present).
    manifest_choice = _infer_library_name_via_manifest(
        target_id,
        ligand_basenames,
        roots,
        cfg,
    )
    if manifest_choice:
        return manifest_choice

    def _manifest_label(manifest_path: Path) -> Optional[str]:
        try:
            with manifest_path.open("r", encoding="utf-8") as fh:
                data = json.load(fh)
        except Exception:
            return None
        for key in ("library_name", "name", "library", "label"):
            val = data.get(key)
            if isinstance(val, str) and val.strip():
                return val.strip()
        return None

    def _collect_basenames_top(folder: Path) -> Set[str]:
        """Non-recursive: only files directly under <library>."""
        allowed = {".pdbqt", ".sdf", ".mol2"}
        out: Set[str] = set()
        try:
            for e in folder.iterdir():
                if e.is_file() and e.suffix.lower() in allowed:
                    out.add(e.name)
        except FileNotFoundError:
            pass
        return out

    # (Optional) shallow peek—uncomment if you want actives/decoys support without full recursion
    def _collect_basenames_shallow(folder: Path) -> Set[str]:
        allowed = {".pdbqt", ".sdf", ".mol2"}
        out: Set[str] = set()
        # top level
        out |= _collect_basenames_top(folder)
        # shallow subdirs commonly used
        for sub in ("actives", "decoys"):
            subdir = folder / sub
            if subdir.is_dir():
                try:
                    for e in subdir.iterdir():
                        if e.is_file() and e.suffix.lower() in allowed:
                            out.add(e.name)
                except FileNotFoundError:
                    pass
        return out

    best_name: Optional[str] = None
    best_score = -1
    have_tie = False

    for root in roots:
        dbg("DEBUG", "library", f"pdb={target_id} scanning_root={root}")
        if not root.exists():
            dbg("WARN", "library", f"pdb={target_id} root_missing={root}")
            continue

        # iterate each <library> folder at the root
        for library_dir in sorted(root.iterdir()):
            if not library_dir.is_dir():
                continue

            dbg(
                "DEBUG",
                "library",
                f"pdb={target_id} candidate_library={library_dir.name}",
            )

            # 1) manifest preference
            manifest = library_dir / "manifest.json"
            if manifest.exists():
                label = _manifest_label(manifest)
                dbg(
                    "DEBUG",
                    "library",
                    f"pdb={target_id} manifest_found={manifest} label={label or '<none>'}",
                )
                if label:
                    dbg(
                        "INFO",
                        "library",
                        f"pdb={target_id} picked='{label}' method='manifest' root={library_dir}",
                    )
                    return label

            # 2) filename-overlap (non-recursive; flip to _collect_basenames_shallow if needed)
            basenames = _collect_basenames_top(library_dir)
            dbg(
                "DEBUG",
                "library",
                f"pdb={target_id} library={library_dir.name} files_seen={len(basenames)} "
                f"sample={list(sorted(basenames))[:3] if basenames else []}",
            )

            if not basenames:
                # Uncomment next two lines to consider shallow subdirs if top-level empty
                # basenames = _collect_basenames_shallow(library_dir)
                # dbg("DEBUG", "library", f"pdb={target_id} library={library_dir.name} shallow_files_seen={len(basenames)}")
                if not basenames:
                    continue

            score = len(ligand_basenames & basenames) if ligand_basenames else 0
            dbg(
                "DEBUG",
                "library",
                f"pdb={target_id} library={library_dir.name} overlap_score={score}",
            )
            if score > best_score:
                best_score = score
                best_name = library_dir.name
                have_tie = False
            elif score == best_score and score >= 0:
                have_tie = True

    if best_name and best_score > 0 and not have_tie:
        dbg(
            "INFO",
            "library",
            f"pdb={target_id} picked='{best_name}' method='filename_overlap' score={best_score}",
        )
        return best_name
    if best_name and have_tie:
        dbg(
            "WARN",
            "library",
            f"pdb={target_id} filename_overlap ties score={best_score}",
        )
    dbg("WARN", "library", f"pdb={target_id} no library match found in roots")
    return ""


def derive_target_name(
    target_id: str, *, prefer: str, pdb_root_override: Optional[Path], cfg: Dict
) -> str:
    pdb_path = Path("input_pdbs") / f"{target_id}.pdb"
    exists = pdb_path.exists()
    print(
        f"[dbg.target.source] pdb={target_id} path=input_pdbs/{target_id}.pdb exists={str(exists)}"
    )
    if not exists:
        return ""

    header_lines = _read_pdb_header_lines(pdb_path)
    compnd = extract_compnd_molecules(header_lines)
    entries, accessions = extract_uniprot_from_dbref(header_lines)
    compnd_label = compnd[0] if compnd else ""
    uniprot_label = accessions[0] if accessions else (entries[0] if entries else "")
    print(
        f"[dbg.target.extract] pdb={target_id} compnd='{compnd_label}' uniprot='{uniprot_label}' prefer={prefer}"
    )
    choice = choose_target_name(compnd, entries, accessions, prefer=prefer)
    source = ""
    compnd_choice = ", ".join(compnd) if compnd else ""
    if choice:
        if compnd_choice and choice == compnd_choice:
            source = "PDB:COMPND"
        elif entries and choice == entries[0]:
            source = "PDB:DBREF_ENTRY"
        elif accessions and choice == accessions[0]:
            source = "PDB:DBREF_ACCESSION"
    source_label = source if source else "none"
    dbg("DEBUG", "target", f"pdb={target_id} name='{choice}' source='{source_label}'")
    return choice
