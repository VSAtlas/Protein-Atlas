"""Microstate registry helpers and ligand enumeration logic."""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import shutil
from pathlib import Path
from typing import Any, Collection, Dict, List, Optional, Set, Union

from path_router import make_paths, load_ph_tags
from rdkit import Chem

from docking.ph_ensemble_docking import _parse_ph_values_from_label

logger = logging.getLogger(__name__)


def compute_microstate_id_from_pdbqt(pdbqt_path: Union[str, Path]) -> str:
    """
    Compute a microstate identifier directly from a prepared PDBQT file.

    We derive the ID from the ATOM/HETATM records only, so it depends on:
      - Atom identity / order
      - Atom types and partial charges
      - Bonding / connectivity as encoded in the atom records

    and ignores comments, REMARKs, and ROOT/ENDROOT blocks.

    Any change in protonation, tautomer, formal charge, or atom typing should
    change the ID. Bitwise-identical PDBQTs (e.g., identical microstate at
    different pH) will share an ID.
    """
    path = Path(pdbqt_path)
    try:
        with path.open("rt", encoding="utf-8", errors="replace") as fh:
            atom_lines: list[str] = []
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    atom_lines.append(line.rstrip("\r\n"))
    except Exception as e:
        logger.warning(
            "[microstate] compute_microstate_id_from_pdbqt_failed path=%s err=%s",
            path,
            e,
        )
        return ""

    if not atom_lines:
        # Fallback: hash entire file if no ATOM/HETATM lines are present
        try:
            payload = path.read_bytes()
        except Exception as e:
            logger.warning(
                "[microstate] compute_microstate_id_from_pdbqt_readbytes_failed path=%s err=%s",
                path,
                e,
            )
            return ""
    else:
        payload = ("\n".join(atom_lines)).encode("utf-8")

    h = hashlib.sha1(payload).hexdigest()
    microstate_id = h[:16]
    logger.debug(
        "[microstate] pdbqt_id path=%s id=%s n_atom_lines=%d",
        path.name,
        microstate_id,
        len(atom_lines),
    )
    return microstate_id


def compute_microstate_id(pdbqt_path: Union[str, Path]) -> str:
    """
    Backwards-compatible wrapper used by prep_ligands_bulk.

    Note: The signature has changed: we now accept a PDBQT path instead of an
    RDKit Mol and derive the microstate ID from the PDBQT content.
    """
    return compute_microstate_id_from_pdbqt(pdbqt_path)


_ONLY_TOKEN_RE = re.compile(r"[0-9]{1,7}")


def _normalize_only_token(tok: str) -> Optional[str]:
    """
    Accepts things like 'rdk_0004931', 'rdk_4931', '0004931', '4931'.
    Returns canonical 'rdk_0004931' or None if it can't be parsed.
    """
    if not tok:
        return None
    t = tok.strip().lower().replace(",", " ")
    if not t:
        return None
    m = _ONLY_TOKEN_RE.search(t)
    if not m:
        return None
    n = m.group(0)
    try:
        i = int(n)
    except Exception:
        return None
    if i < 0 or i > 9_999_999:
        return None
    return f"rdk_{i:07d}"


def _collect_only_from_env_and_cli(cli_only: Optional[List[str]] = None) -> Set[str]:
    """
    Merge LIGPREP_ONLY (env), LIGPREP_ONLY_FILE (env path), and --only (CLI list).
    Normalize all tokens to 'rdk_0000000'. Empty/invalid tokens are ignored.
    """
    out: Set[str] = set()

    env_only = os.environ.get("LIGPREP_ONLY", "")
    if env_only:
        for raw in re.split(r"[,\s]+", env_only.strip()):
            norm = _normalize_only_token(raw)
            if norm:
                out.add(norm)

    env_file = os.environ.get("LIGPREP_ONLY_FILE", "")
    if env_file:
        try:
            with open(env_file, "r", encoding="utf-8") as fh:
                for line in fh:
                    for raw in re.split(r"[,\s]+", line.strip()):
                        norm = _normalize_only_token(raw)
                        if norm:
                            out.add(norm)
        except Exception as e:
            logging.warning(
                "[test-mode] could not read LIGPREP_ONLY_FILE=%s: %s", env_file, e
            )

    if cli_only:
        for raw in cli_only:
            for tok in re.split(r"[,\s]+", raw.strip()):
                norm = _normalize_only_token(tok)
                if norm:
                    out.add(norm)

    return out


def load_microstate_registry(
    library_out_dir: Path, library_name: str
) -> tuple[dict, dict]:
    """
    Load or initialize the microstate registry for a given library.
    Returns a dict with keys: version, library, microstates (list) and an index mapping microstate_id -> entry.
    """
    registry_path = library_out_dir / "microstates.json"
    default_registry = {
        "version": 1,
        "library": library_name,
        "microstates": [],
    }

    registry = default_registry.copy()
    if registry_path.exists():
        try:
            registry = json.loads(registry_path.read_text(encoding="utf-8"))
        except Exception as e:
            logging.warning(
                "[microstate] registry_load_failed path=%s err=%s", registry_path, e
            )
            registry = default_registry.copy()
        else:
            if registry.get("library") and registry["library"] != library_name:
                logging.warning(
                    "[microstate] registry_library_mismatch path=%s expected=%s found=%s",
                    registry_path,
                    library_name,
                    registry.get("library"),
                )

    microstate_index: Dict[str, dict] = {}
    for entry in registry.get("microstates", []) or []:
        microstate_id = entry.get("microstate_id")
        if microstate_id:
            microstate_index[microstate_id] = entry

    return registry, microstate_index


def _rehydrate_microstate_registry_from_ph_pdbqts(
    library_out_dir: Path,
    library_name: str,
) -> dict:
    """
    Rebuild a microstate registry from existing pH-specific PDBQT files when the
    registry is empty (e.g., after a previous partial write or resume skip).
    """
    registry, index = load_microstate_registry(library_out_dir, library_name)
    if registry.get("microstates"):
        return registry

    ph_dirs = [
        d
        for d in library_out_dir.iterdir()
        if d.is_dir() and d.name.lower().startswith("ph")
    ]
    if not ph_dirs:
        return registry

    rebuilt = {
        "version": registry.get("version", 1),
        "library": registry.get("library", library_name),
        "microstates": [],
    }
    ms_index: Dict[str, dict] = {}

    for ph_dir in sorted(ph_dirs):
        ph_label = ph_dir.name
        ph_vals = _parse_ph_values_from_label(ph_label)
        ph_value = ph_vals[0] if ph_vals else None

        for pdbqt_path in ph_dir.glob("*.pdbqt"):
            microstate_id = compute_microstate_id_from_pdbqt(pdbqt_path)
            if not microstate_id:
                continue

            base = pdbqt_path.stem.replace(ph_label, "")
            while "__" in base:
                base = base.replace("__", "_")
            base = base.strip("_")
            ligand_stem = base if base else pdbqt_path.stem

            entry = ms_index.get(microstate_id)
            if entry is None:
                canonical_name = f"{ligand_stem}__ms_{microstate_id}.pdbqt"
                canonical_rel = f"microstates/{canonical_name}"
                canonical_path = library_out_dir / canonical_rel
                try:
                    canonical_path.parent.mkdir(parents=True, exist_ok=True)
                    if not canonical_path.exists():
                        shutil.copy2(pdbqt_path, canonical_path)
                except Exception as e:
                    logger.warning(
                        "[microstate.rebuild] copy_failed ligand=%s src=%s dst=%s err=%s",
                        ligand_stem,
                        pdbqt_path,
                        canonical_path,
                        e,
                    )

                entry = {
                    "microstate_id": microstate_id,
                    "pdbqt_path": canonical_rel,
                    "aliases": [],
                }
                rebuilt["microstates"].append(entry)
                ms_index[microstate_id] = entry

            alias = {"ligand_stem": ligand_stem, "ph_label": ph_label}
            if ph_value is not None:
                alias["ph_value"] = ph_value
            aliases = entry.setdefault("aliases", [])
            if alias not in aliases:
                aliases.append(alias)

    save_microstate_registry(library_out_dir, rebuilt)
    return rebuilt


def save_microstate_registry(library_out_dir: Path, registry: dict) -> None:
    """
    Save the microstate registry to microstates.json atomically.

    Merges with the existing registry under a simple file lock to avoid
    clobbering entries when multiple workers write in parallel.
    """
    import errno
    import time

    final_path = library_out_dir / "microstates.json"
    tmp_path = library_out_dir / "microstates.json.tmp"
    lock_path = library_out_dir / "microstates.json.lock"

    fd = None
    while True:
        try:
            fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            break
        except OSError as e:
            if e.errno == errno.EEXIST:
                time.sleep(0.05)
                continue
            raise

    try:
        library_name = registry.get("library") or "ligprep"
        try:
            existing_registry, _ = load_microstate_registry(
                library_out_dir, library_name
            )
        except Exception:
            existing_registry = {
                "version": registry.get("version", 1),
                "library": library_name,
                "microstates": [],
            }

        existing_microstates = list(existing_registry.get("microstates") or [])
        index: Dict[str, dict] = {}
        for entry in existing_microstates:
            mid = entry.get("microstate_id")
            if mid:
                index[mid] = entry

        new_microstates = registry.get("microstates") or []
        for entry in new_microstates:
            mid = entry.get("microstate_id")
            if not mid:
                continue
            existing_entry = index.get(mid)
            if existing_entry is None:
                existing_microstates.append(entry)
                index[mid] = entry
                continue

            # Merge aliases for existing microstates so parallel writers don't drop metadata
            existing_aliases = existing_entry.get("aliases") or []
            new_aliases = entry.get("aliases") or []
            for alias in new_aliases:
                if alias and alias not in existing_aliases:
                    existing_aliases.append(alias)
            existing_entry["aliases"] = existing_aliases

        merged = existing_registry
        merged["microstates"] = existing_microstates
        if "version" not in merged and "version" in registry:
            merged["version"] = registry["version"]
        if "library" not in merged and registry.get("library"):
            merged["library"] = registry["library"]

        payload = json.dumps(merged, indent=2, sort_keys=True)
        tmp_path.write_text(payload, encoding="utf-8")
        os.replace(tmp_path, final_path)
    finally:
        try:
            if fd is not None:
                os.close(fd)
        finally:
            try:
                os.remove(lock_path)
            except FileNotFoundError:
                pass


def _collect_ligand_stems_with_microstates(
    registry: Dict[str, Any],
    requested_ph_values: Optional[Set[float]] = None,
) -> Set[str]:
    """
    Return the set of canonical ligand stems present in microstates.json.
    """

    stems: Set[str] = set()
    microstates = registry.get("microstates") or []
    if not microstates:
        return stems

    for entry in microstates:
        aliases = entry.get("aliases") or []
        if not aliases:
            continue

        if requested_ph_values:
            for a in aliases:
                try:
                    ph_val = float(a.get("ph_value", -999.0))
                except Exception:
                    continue
                if ph_val in requested_ph_values:
                    stem = a.get("ligand_stem")
                    if stem:
                        stems.add(stem)
                    break
        else:
            for a in aliases:
                stem = a.get("ligand_stem")
                if stem:
                    stems.add(stem)
                    break

    return stems


def _collect_expected_ligand_stems_from_library(
    library_out_dir: Path,
    library_name: str,
) -> Set[str]:
    """
    Infer canonical ligand stems expected for a library by parsing extracted SDFs.
    """

    expected: Set[str] = set()

    project_root = None
    try:
        for parent in library_out_dir.parents:
            if parent.name == "prepped_ligands":
                project_root = parent.parent
                break
        if project_root is None:
            project_root = library_out_dir.parent.parent
    except Exception:
        project_root = library_out_dir.parent.parent

    extracted_root = project_root / "extracted_ligands" / library_name
    if not extracted_root.is_dir():
        return expected

    sdf_paths: List[Path] = []
    for p in sorted(extracted_root.rglob("*.sdf")):
        sdf_paths.append(p)
    for p in sorted(extracted_root.rglob("*.sdf.gz")):
        sdf_paths.append(p)

    for sdf_path in sdf_paths:
        try:
            stem_base = sdf_path.stem
            if stem_base.endswith(".sdf"):
                stem_base = stem_base[:-4]
            supplier = None
            if sdf_path.suffix == ".gz":
                import gzip

                with gzip.open(sdf_path, "rb") as fh:
                    _ = fh.read()
                continue
            else:
                supplier = Chem.SDMolSupplier(
                    str(sdf_path), removeHs=False, sanitize=False
                )

            if supplier is None:
                continue

            for idx, mol in enumerate(supplier):
                if mol is None or mol.GetNumAtoms() == 0:
                    continue
                stem = f"{stem_base}_{idx + 1:05d}"
                expected.add(stem)
        except Exception:
            continue

    return expected


def enumerate_ligands_for_docking(
    requested_ph_values: Optional[Collection[float]] = None,
    *,
    cfg: Optional[Dict] = None,
    pdb_id: Optional[str] = None,
    root_dir: Optional[str | Path] = None,
    microstate_dedup: bool = True,
    force: bool = False,
) -> List[Path]:
    """Enumerate canonical ligand PDBQTs for docking, deduplicated at the microstate level."""

    from prep_ligands.prep_ligands_bulk import (
        prep_ligands_with_mgltools,
    )  # lazy import to avoid circular

    _ = cfg
    pdb_label = (pdb_id or "LIGPREP").upper()

    root_dir_path = Path(root_dir).resolve() if root_dir is not None else None
    out_dir_env = (os.environ.get("LIGPREP_OUT_DIR", "") or "").strip() or None

    library_hint: Optional[str] = None
    if root_dir_path is not None:
        library_hint = root_dir_path.name.lower()
    elif out_dir_env:
        library_hint = Path(out_dir_env).name.lower()

    library_env = (os.environ.get("LIGPREP_LIBRARY", "") or "").strip()
    library_cfg = (cfg.get("LIGPREP_LIBRARY", "") or "").strip() if cfg else ""

    library_base: Optional[str] = library_hint
    if not library_base:
        if library_env:
            library_base = library_env.lower()
        elif library_cfg:
            library_base = library_cfg.lower()

    library_name = (library_base or "ligprep").lower()

    prepped_root_env = (os.environ.get("PREPPED_LIGANDS_DIR", "") or "").strip()
    if prepped_root_env:
        prepped_root = Path(prepped_root_env).resolve()
    else:
        try:
            paths = (
                make_paths(cfg, base_id=pdb_label, pdb_file=f"{pdb_label}.pdb")
                if cfg is not None
                else None
            )
            prepped_root = (
                paths.prepped_ligands_dir.parent
                if paths is not None
                else Path("prepped_ligands").resolve()
            )
        except Exception:
            prepped_root = Path("prepped_ligands").resolve()

    if root_dir_path is not None:
        library_out_dir = root_dir_path
        library_source = "root_dir"
    else:
        library_out_dir = (
            Path(out_dir_env).resolve() if out_dir_env else prepped_root / library_name
        )
        library_source = "config/env"

    registry_path = library_out_dir / "microstates.json"

    requested_set: Optional[Set[float]] = None
    if requested_ph_values:
        requested_set = {float(ph) for ph in requested_ph_values}

    # Auto-derive pH window from pH ensemble labels when none is provided explicitly
    if requested_set is None and cfg is not None and pdb_id:
        try:
            mode = str(cfg.get("PH_LIGAND_MODE", "off")).lower()
            ph_ensemble_enabled = bool(cfg.get("PH_ENSEMBLE"))
            if ph_ensemble_enabled and mode not in ("off", ""):
                ph_tags = load_ph_tags(pdb_id, variant=None) or []

                context_phs: list[float] = []
                for tag in ph_tags:
                    for ph in _parse_ph_values_from_label(tag):
                        context_phs.extend((ph - 1.0, ph, ph + 1.0))

                ligand_window = sorted(
                    {round(ph, 1) for ph in context_phs if 0.0 < ph < 15.0}
                )

                if ligand_window:
                    requested_set = set(ligand_window)
                    logger.info(
                        "enumerate_ligands_for_docking: auto pH window from ensemble pdb=%s mode=%s ph_values=%s",
                        (pdb_id or "").upper(),
                        mode,
                        ligand_window,
                    )
        except Exception as e:
            logger.warning(
                "enumerate_ligands_for_docking: failed to derive pH window from ensemble for %s: %s",
                (pdb_id or "").upper(),
                e,
            )

    logger.info(
        "enumerate_ligands_for_docking: root_dir=%s ph_values=%s microstate_dedup=%s force=%s",
        str(root_dir_path) if root_dir_path is not None else "",
        sorted(requested_set) if requested_set is not None else [],
        microstate_dedup,
        force,
    )

    logger.info(
        "enumerate_ligands_for_docking: pdb=%s library_hint=%s library_env=%s library_cfg=%s resolved_library=%s "
        "prepped_root=%s library_out_dir=%s registry=%s source=%s",
        pdb_label,
        library_hint,
        library_env,
        library_cfg,
        library_name,
        str(prepped_root),
        library_out_dir,
        registry_path,
        library_source,
    )

    only_set = _collect_only_from_env_and_cli(None)
    only_for_microstate: Optional[Set[str]] = only_set or None

    def _run_microstate_prep_for_phs(ph_values: Collection[float]) -> None:
        if not ph_values:
            return
        logger.info(
            "enumerate_ligands_for_docking: running microstate prep for missing pH values: %s",
            ",".join(f"{ph:.2f}" for ph in sorted(set(ph_values))),
        )
        prep_ligands_with_mgltools(
            force=False,
            only=only_for_microstate,
            ph_values=sorted(set(float(ph) for ph in ph_values)),
            microstate_dedup=True,
            root_dir=library_out_dir,
        )

    if requested_set is not None:
        if not registry_path.exists():
            _run_microstate_prep_for_phs(requested_set)
        registry, _microstate_index = load_microstate_registry(
            library_out_dir, library_name
        )
        if not (registry.get("microstates") or []):
            try:
                registry = _rehydrate_microstate_registry_from_ph_pdbqts(
                    library_out_dir, library_name
                )
                _, _microstate_index = load_microstate_registry(
                    library_out_dir, library_name
                )
            except Exception as e:
                logger.warning(
                    "[microstate.rebuild.skip] library=%s reason=%s", library_name, e
                )

        existing_ph_values: Set[float] = set()
        for entry in registry.get("microstates", []) or []:
            for alias in entry.get("aliases", []) or []:
                ph_value = alias.get("ph_value")
                if ph_value is None:
                    continue
                try:
                    existing_ph_values.add(float(ph_value))
                except Exception:
                    continue

        logger.info(
            "enumerate_ligands_for_docking: library=%s requested_ph=%s existing_ph=%s",
            library_name,
            ",".join(f"{p:.2f}" for p in sorted(requested_set)),
            ",".join(f"{p:.2f}" for p in sorted(existing_ph_values)),
        )

        missing = requested_set - existing_ph_values
        if missing:
            logger.info(
                "enumerate_ligands_for_docking: library=%s missing_ph=%s (will trigger microstate prep)",
                library_name,
                ",".join(f"{p:.2f}" for p in sorted(missing)),
            )
        else:
            logger.info(
                "enumerate_ligands_for_docking: library=%s no missing_ph; microstates already cover requested pH values",
                library_name,
            )
        if missing:
            logger.info(
                "enumerate_ligands_for_docking: requested_ph=%s missing_ph=%s",
                sorted(requested_set),
                sorted(missing),
            )
            _run_microstate_prep_for_phs(missing)
            registry, _microstate_index = load_microstate_registry(
                library_out_dir, library_name
            )
    else:
        registry, _microstate_index = load_microstate_registry(
            library_out_dir, library_name
        )

    microstate_entries = registry.get("microstates") or []

    if microstate_dedup and microstate_entries and requested_set is not None:
        expected_ligands = _collect_expected_ligand_stems_from_library(
            library_out_dir=library_out_dir,
            library_name=library_name,
        )

        if expected_ligands:
            covered_ligands = _collect_ligand_stems_with_microstates(
                registry,
                requested_set,
            )

            missing_ligands = expected_ligands - covered_ligands
            coverage_fraction = len(covered_ligands) / max(len(expected_ligands), 1)

            logger.info(
                "[microstate.coverage] library=%s requested_ph=%s expected=%d covered=%d missing=%d coverage=%.3f",
                library_name,
                sorted(requested_set),
                len(expected_ligands),
                len(covered_ligands),
                len(missing_ligands),
                coverage_fraction,
            )

    if (not registry_path.exists()) or not (registry.get("microstates") or []):
        if not registry_path.exists():
            logger.warning(
                "enumerate_ligands_for_docking: no microstate registry found at %s",
                registry_path,
            )
        else:
            logger.warning(
                "enumerate_ligands_for_docking: microstate registry %s has no entries; attempting manifest fallback",
                registry_path,
            )

        manifest_path = library_out_dir / "_manifest.json"
        if manifest_path.exists():
            try:
                manifest = json.loads(manifest_path.read_text())
            except Exception as e:
                logger.warning(
                    "enumerate_ligands_for_docking: failed to read manifest %s: %s; returning empty list",
                    manifest_path,
                    e,
                )
                return []

            entries = manifest.get("entries") or {}
            pdbqt_paths: List[Path] = []
            for key, rel_name in entries.items():
                p = library_out_dir / rel_name
                try:
                    if p.exists() and p.stat().st_size > 100:
                        pdbqt_paths.append(p)
                    else:
                        logger.debug(
                            "enumerate_ligands_for_docking: skipping manifest entry %s -> %s (missing or too small)",
                            key,
                            p,
                        )
                except OSError:
                    logger.debug(
                        "enumerate_ligands_for_docking: skipping manifest entry %s -> %s (stat failed)",
                        key,
                        p,
                    )

            if pdbqt_paths:
                logger.info(
                    "enumerate_ligands_for_docking: manifest fallback returning %d PDBQT(s) (ignoring requested_ph_values=%s)",
                    len(pdbqt_paths),
                    "any" if requested_set is None else sorted(requested_set),
                )
                return sorted(pdbqt_paths, key=lambda p: p.name)

        logger.warning(
            "enumerate_ligands_for_docking: no usable microstates or manifest entries; returning empty list",
        )
        return []

    result: Set[Path] = set()
    total_microstates = len(registry.get("microstates", []) or [])

    sample_limit = 20
    sample_count = 0

    for entry in registry.get("microstates", []) or []:
        rel_path = entry.get("pdbqt_path")
        if not rel_path:
            continue

        aliases = entry.get("aliases", []) or []
        canonical_path = library_out_dir / rel_path

        include = requested_set is None
        selected_alias: Optional[dict] = None
        if requested_set is not None:
            for alias in aliases:
                ph_value = alias.get("ph_value")
                if ph_value is None:
                    continue
                try:
                    pv = float(ph_value)
                except Exception:
                    continue
                if pv in requested_set:
                    include = True
                    selected_alias = alias
                    break
        else:
            if aliases:
                selected_alias = aliases[0]

        if not include:
            continue

        if canonical_path.exists() and canonical_path.stat().st_size > 100:
            if sample_count < sample_limit:
                logger.info(
                    "enumerate_ligands_for_docking.selection: stem=%s ph_label=%s ph_value=%s microstate_id=%s pdbqt=%s",
                    (selected_alias or {}).get("ligand_stem"),
                    (selected_alias or {}).get("ph_label"),
                    (selected_alias or {}).get("ph_value"),
                    entry.get("microstate_id"),
                    canonical_path,
                )
                sample_count += 1
            result.add(canonical_path)

    if requested_set is None:
        logger.info(
            "enumerate_ligands_for_docking: mode=all_ph microstates=%d included=%d",
            total_microstates,
            len(result),
        )
    else:
        logger.info(
            "enumerate_ligands_for_docking: mode=filtered ph_values=%s microstates=%d included=%d",
            sorted(requested_set),
            total_microstates,
            len(result),
        )

    return sorted(result, key=lambda p: p.name)


__all__ = [
    "load_microstate_registry",
    "save_microstate_registry",
    "_collect_ligand_stems_with_microstates",
    "_collect_expected_ligand_stems_from_library",
    "_collect_only_from_env_and_cli",
    "enumerate_ligands_for_docking",
]
