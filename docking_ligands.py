from __future__ import annotations

import json
import logging
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import math
from rdkit import Chem
from rdkit.Chem import FilterCatalog, rdFMCS, rdMolAlign
from rdkit.Chem.MolStandardize import rdMolStandardize

from library_index import LibraryIndex
from path_router import Paths
from prep_ligands import enumerate_ligands_for_docking, prep_ligands_from_pdb


def norm(p: str | Path) -> str:
    """Normalize path to a clean, forward-slash string for logs & keys."""
    return os.path.abspath(str(p)).replace("\\", "/")


def _selection_schedule(docking_mode: str, n_stages: Optional[int]) -> List[float]:
    """
    Selection percentages per stage.
    Mirrors the existing Vina staging defaults.
    """
    default_disc = [1.0, 0.1, 0.01, 0.001, 0.001]
    default_poly = [1.0, 0.05, 0.005]
    sched = {
        "discovery": default_disc,
        "polypharmacology": default_poly,
    }.get(docking_mode, [1.0] * n_stages if n_stages is not None else [1.0])
    # Trim/extend to requested number of stages to avoid index errors.
    if n_stages is None:
        return list(sched)
    if len(sched) < n_stages:
        sched = sched + [sched[-1]] * (n_stages - len(sched))
    elif len(sched) > n_stages:
        sched = sched[:n_stages]
    return sched


def _is_pytest_context(cfg: dict) -> bool:
    """
    Detect pytest-driven runs so we can swap HMDB roots for lightweight fixtures.
    """
    try:
        selection_mode = str(cfg.get("PDB_SELECTION_MODE", "")).strip().lower()
    except Exception:
        selection_mode = ""
    env_selection = str(os.environ.get("PDB_SELECTION_MODE", "")).strip().lower()
    return (
        "pytest" in os.environ.get("PYTEST_CURRENT_TEST", "")
        or selection_mode == "test_library_map"
        or env_selection == "test_library_map"
    )


def _iter_pdbqt_dirfirst(root: Path, allowed_subdirs: Optional[set[str]] = None):
    """
    Yield .pdbqt files with a directory-first strategy:
      - files directly under `root`
      - then files under first-level subdirs (filter via `allowed_subdirs` if provided)
    Falls back to recursive glob if any listing fails to keep things robust.
    """
    try:
        if not root or not root.exists():
            return
        for p in root.glob("*.pdbqt"):
            yield p
        for d in root.iterdir():
            if not d.is_dir():
                continue
            if allowed_subdirs is not None and d.name not in allowed_subdirs:
                continue
            for p in d.glob("*.pdbqt"):
                yield p
    except Exception:
        for p in root.rglob("*.pdbqt"):
            yield p


def select_ligands_for_next(
        docking_mode: str,
        i: int,
        stages: List[Dict],
        scores: Dict[str, float],
        logger: logging.Logger,
        base_pool_n: Optional[int] = None,          #  if provided, select % of this
        force_include: Optional[set] = None         #  always add these
) -> List[str]:
    if not scores:
        # Still allow force-carry if provided and next stage exists
        return sorted(force_include) if force_include else []

    schedule = _selection_schedule(docking_mode, len(stages))
    pct = schedule[i + 1] if i + 1 < len(schedule) else 0.01

    # Use provided base if given (e.g., Stage1 pool size) -- otherwise fall back to valid-count
    pool_n = base_pool_n if (base_pool_n is not None) else len(scores)

    # Select K by the base pool, but cap at the number of valid scores available
    k_target = max(1, int(pool_n * pct))
    k = max(1, min(k_target, len(scores)))

    # take best k from valid scores
    next_list = [l for l, _ in sorted(scores.items(), key=lambda kv: kv[1])[:k]]

    # Force-carry: add any requested ligands (e.g., extracted controls) to the next stage
    if force_include:
        # maintain stable order: extend with any forced ligands not already selected
        in_set = set(next_list)
        forced_add = [l for l in sorted(force_include) if l not in in_set]
        next_list.extend(forced_add)
        if forced_add:
            logger.info(f"[Force-carry] Added {len(forced_add)} extracted ligands to next stage.")

    logger.info(
        f"Selected {k} by score (+{len(force_include or [])} forced) "
        f"= {len(next_list)} total ({pct * 100:.5f}% of base={pool_n})."
    )
    return next_list


def compute_stage_membership_from_scores(
    cfg: Dict[str, Any],
    docking_mode: str,
    scores: Dict[str, float],
    *,
    higher_is_better: bool,
    n_stages: Optional[int] = None,
) -> Dict[int, List[str]]:
    """
    Compute stage membership using the same percentile logic as Vina staging.

    Returns {stage_index (1-based): [ligands]}.
    """
    del cfg  # cfg reserved for future config-driven selection; staging matches current Vina defaults.
    if not scores:
        return {}

    # Filter to numeric scores only
    valid_scores = {l: s for l, s in scores.items() if isinstance(s, (int, float)) and math.isfinite(s)}
    if not valid_scores:
        return {}

    order = [l for l, _ in sorted(valid_scores.items(), key=lambda kv: kv[1], reverse=higher_is_better)]
    total = len(order)
    schedule = _selection_schedule(docking_mode, n_stages)

    # Assign most selective stages first to avoid duplicates, then return in ascending stage order.
    remaining = list(order)
    stage_membership: Dict[int, List[str]] = {}
    for stage_idx in range(len(schedule), 0, -1):
        pct = schedule[stage_idx - 1]
        target = int(math.ceil(total * pct))
        target = max(1, target)
        if stage_idx == 1:
            assign_count = len(remaining)
        else:
            assign_count = min(target, len(remaining))
        stage_membership[stage_idx] = remaining[:assign_count]
        remaining = remaining[assign_count:]

    return {k: stage_membership[k] for k in sorted(stage_membership)}


def _count_heavy_atoms_from_pdbqt(pdbqt_path: Path) -> int:
    heavy = 0
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            element = line[76:78].strip() if len(line) >= 78 else ""
            if element:
                if element.upper() != "H":
                    heavy += 1
            else:
                atom_name = line[12:16].strip()
                if not atom_name.upper().startswith("H"):
                    heavy += 1
    return heavy


def _load_mol_any(pdbqt_path: Path, obabel_exe: str | None) -> Chem.Mol | None:
    base = pdbqt_path.with_suffix("")
    # prefer SDF, then MOL2, then PDB
    sdf = base.with_suffix(".sdf"); mol2 = base.with_suffix(".mol2"); pdb = base.with_suffix(".pdb")
    if sdf.exists():
        supp = Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=True)
        for m in supp:
            if m: return m
    for fp, reader in [(mol2, Chem.MolFromMol2File), (pdb, Chem.MolFromPDBFile)]:
        if fp.exists():
            m = reader(str(fp), sanitize=True, removeHs=False)
            if m: return m
    # fallback: PDBQT -> SDF via obabel (Linux-friendly)
    if obabel_exe:
        with tempfile.TemporaryDirectory() as td:
            out = Path(td, "tmp.sdf")
            try:
                subprocess.check_call([obabel_exe, "-ipdbqt", str(pdbqt_path), "-osdf", "-O", str(out), "--retype", "--addh"])
                supp = Chem.SDMolSupplier(str(out), removeHs=False, sanitize=True)
                for m in supp:
                    if m: return m
            except Exception:
                return None
    return None


def _standardize(m: Chem.Mol) -> Chem.Mol:
    parent = rdMolStandardize.ChargeParent(m)   # neutralize/parent
    rdMolStandardize.Normalize(parent)          # FG normalization
    Chem.SanitizeMol(parent)
    return parent


# Build catalog with PAINS A/B/C
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
pains_catalog = FilterCatalog.FilterCatalog(params)


def _coerce_test_map(m) -> Dict[str, str]:
    import json as _json, ast as _ast

    if isinstance(m, dict):
        return {str(k).upper(): str(v) for k, v in m.items()}
    s = str(m).strip()
    if not s:
        return {}
    parsed = None
    try:
        parsed = _json.loads(s)
    except Exception:
        try:
            parsed = _ast.literal_eval(s)
        except Exception:
            parsed = {}
    return {str(k).upper(): str(v) for k, v in (parsed if isinstance(parsed, dict) else {}).items()}


def _resolve_test_mode(cfg) -> str:
    """
    Normalize TEST_MODE_ENABLE to one of:
      "off", "dud", "fda", "fda+dud", "hmdb", "hmdb+dud", "hmdb+fda", "fda+dud+hmdb".

    Accepts:
      - Booleans / bool-like strings for backwards compatibility:
          True  / "true" / "yes" / "on" / "1"  -> "fda+dud"
          False / "false" / "no"  / "off" / "0" / "" / None -> "off"
      - Explicit string modes (case-insensitive, tolerant of "_" / "+" / "-"):
          "off"                 -> "off"
          "dud"                 -> "dud"
          "fda"                 -> "fda"
          "fda+dud"/"both"      -> "fda+dud"
          "hmdb"                -> "hmdb"
          "hmdb+dud"            -> "hmdb+dud"
          "hmdb+fda"            -> "hmdb+fda"
          "fda+dud+hmdb"        -> "fda+dud+hmdb"

    Any unrecognized string logs a warning and falls back to "off".
    """
    #env has priority over config
    raw = os.environ.get("TEST_MODE_ENABLE", cfg.get("TEST_MODE_ENABLE", "off"))

    if isinstance(raw, bool):
        return "fda+dud" if raw else "off"

    s = str(raw).strip()
    if not s:
        return "off"

    s_lower = s.lower()
    if s_lower in ("", "0", "false", "no", "off", "none", "null"):
        return "off"
    if s_lower in ("true", "yes", "on", "1"):
        return "fda+dud"

    direct_aliases = {
        "fda": "fda",
        "both": "fda+dud",
        "fda_dud": "fda+dud",
        "fda+dud": "fda+dud",
        "dud+fda": "fda+dud",
    }
    if s_lower in direct_aliases:
        return direct_aliases[s_lower]

    normalized = (
        s_lower.replace("_", "+")
        .replace("-", "+")
        .replace(" ", "")
        .replace("only", "")
    ).strip("+")
    if normalized in direct_aliases:
        return direct_aliases[normalized]

    tokens = [tok for tok in normalized.split("+") if tok and tok != "and"]
    token_set = set(tokens)
    if token_set == {"dud"}:
        return "dud"
    if token_set == {"fda"}:
        return "fda"
    if token_set == {"hmdb"}:
        return "hmdb"
    if token_set == {"fda", "dud"}:
        return "fda+dud"
    if token_set == {"hmdb", "dud"}:
        return "hmdb+dud"
    if token_set == {"hmdb", "fda"}:
        return "hmdb+fda"
    if token_set == {"fda", "dud", "hmdb"}:
        return "fda+dud+hmdb"
    if normalized in ("fdaadud", "fdaanddud"):
        return "fda+dud"

    print(f"[test-mode] WARNING: Unknown TEST_MODE_ENABLE={raw!r}; treating as 'off'.")
    return "off"

def _dedup_index_roots(seq: list[Path]) -> list[Path]:
    deduped: list[Path] = []
    seen: set[str] = set()
    for candidate in seq:
        if not candidate:
            continue
        path_obj = Path(candidate)
        try:
            key = str(path_obj.resolve())
        except Exception:
            key = str(path_obj)
        if key not in seen:
            seen.add(key)
            deduped.append(path_obj)
    return deduped


def _lib_roots_for_pdb(
    cfg: Dict,
    pdb_id: str,
    paths: Paths,
    logger: logging.Logger,
    *,
    test_mode_override: Optional[str] = None,
) -> tuple[list[Path], list[Path]]:
    extra_dirs = str(cfg.get("LIBRARY_EXTRA_DIRS", "")).strip()
    extra_paths: list[Path] = []
    if extra_dirs:
        for d in extra_dirs.split(";"):
            d = d.strip()
            if not d:
                continue
            p = Path(d)
            if p.exists():
                extra_paths.append(p)

    subdir_default = str(cfg.get("LIBRARY_SUBDIR_DEFAULT", "fda_library"))
    hmdb_subdir_raw = str(cfg.get("HMDB_LIBRARY_SUBDIR", "hmdb"))
    hmdb_test_subdir = str(cfg.get("HMDB_TEST_LIBRARY_SUBDIR", "hmdb_test_library_10"))
    test_mode = test_mode_override or _resolve_test_mode(cfg)
    pytest_mode = _is_pytest_context(cfg)
    hmdb_use_test = pytest_mode and ("hmdb" in str(test_mode))
    hmdb_subdir = hmdb_test_subdir if hmdb_use_test else hmdb_subdir_raw

    maybe_map = cfg.get("TEST_LIBRARY_MAP", {})
    test_map: Dict[str, str] = {}
    test_map = _coerce_test_map(maybe_map)
    logger.info(f"[lib-roots.map] raw_type={type(maybe_map).__name__} keys={len(test_map)}")
    mapped_value = (test_map or {}).get(pdb_id)

    base_root = Path(
        cfg.get("OUTPUT_LIGANDS_DIR")
        or cfg.get("PREPPED_LIGANDS_ROOT")
        or "prepped_ligands"
    )

    if hmdb_use_test and hmdb_subdir != hmdb_subdir_raw:
        logger.info(
            "[ligands.hmdb-test] test_mode=%s pytest=%s subdir=%s raw=%s",
            test_mode,
            pytest_mode,
            hmdb_subdir,
            hmdb_subdir_raw,
        )

    hmdb_root = base_root / hmdb_subdir
    dud_root = (base_root / mapped_value) if mapped_value else None

    need_dud_map = test_mode in ("dud", "fda+dud", "hmdb+dud", "fda+dud+hmdb")
    if need_dud_map and not mapped_value:
        logger.warning(
            "[test-mode] PDB %s missing from TEST_LIBRARY_MAP; using default library=%s",
            pdb_id,
            subdir_default,
        )
    dud_roots: list[Path] = []
    if dud_root:
        dud_roots = [dud_root]
    elif need_dud_map:
        dud_roots = [base_root / subdir_default]

    roots_for_mode: list[Path]
    if test_mode == "off":
        roots_for_mode = [base_root / subdir_default]
    elif test_mode == "dud":
        roots_for_mode = dud_roots or [base_root / subdir_default]
    elif test_mode == "fda+dud":
        roots_for_mode = dud_roots + [base_root / subdir_default] if dud_roots else [base_root / subdir_default]
    elif test_mode == "hmdb":
        roots_for_mode = [hmdb_root]
    elif test_mode == "hmdb+dud":
        roots_for_mode = [hmdb_root] + (dud_roots or [base_root / subdir_default])
    elif test_mode == "hmdb+fda":
        roots_for_mode = [hmdb_root, base_root / subdir_default]
    elif test_mode == "fda+dud+hmdb":
        if dud_roots:
            roots_for_mode = dud_roots + [hmdb_root, base_root / subdir_default]
        else:
            roots_for_mode = [hmdb_root, base_root / subdir_default]
    else:
        roots_for_mode = [base_root / subdir_default]

    deduped_roots = _dedup_index_roots(roots_for_mode)
    allowed_noncontrol_roots: list[Path] = []
    for r in deduped_roots:
        if r.exists():
            allowed_noncontrol_roots.append(r)
        else:
            logger.warning("[ligands.test-roots] missing=%s", r)

    extra_paths_cfg: list[str] = cfg.get("EXTRA_LIGAND_ROOTS", []) or []
    for d in extra_paths_cfg:
        p = Path(d)
        if p.exists():
            allowed_noncontrol_roots.append(p)
        else:
            logger.warning("[ligands.extra-roots] missing=%s", p)

    allowed_noncontrol_roots.extend(extra_paths)
    allowed_noncontrol_roots = _dedup_index_roots(allowed_noncontrol_roots)

    logger.info(
        "[ligands.allowed-roots] test_mode=%s pdb=%s roots=%d",
        test_mode,
        pdb_id,
        len(allowed_noncontrol_roots),
    )
    cfg["_ALLOWED_NONCONTROL_ROOTS"] = [str(p) for p in allowed_noncontrol_roots]

    primary_root_str = str(allowed_noncontrol_roots[0]) if allowed_noncontrol_roots else None
    logger.info(
        "[ph_ligand.roots] test_mode=%s pdb=%s ph_root=%s noncontrol_roots=%s",
        test_mode,
        pdb_id,
        primary_root_str,
        cfg.get("_ALLOWED_NONCONTROL_ROOTS"),
    )
    cfg["_TEST_MODE_EFFECTIVE"] = test_mode
    return [], allowed_noncontrol_roots


def prepare_and_filter_ligands(
    cfg: Dict,
    paths: Paths,
    logger: logging.Logger,
    *,
    run_mode: Optional[str] = None,
) -> Tuple[List[str], Dict[str, int], Dict[str, bool]]:
    """
    Gathers candidate ligands, keeps existing validation/PAINS logic, and
    filters the *non-control* pool to allowed library roots:

      - TEST_MODE_ENABLE="off"      -> OUTPUT_LIGANDS_DIR/<LIBRARY_SUBDIR_DEFAULT>
      - TEST_MODE_ENABLE="dud"      -> OUTPUT_LIGANDS_DIR/<mapped_subdir>
      - TEST_MODE_ENABLE="fda+dud"  -> OUTPUT_LIGANDS_DIR/<mapped_subdir> + OUTPUT_LIGANDS_DIR/<LIBRARY_SUBDIR_DEFAULT>
      - TEST_MODE_ENABLE="hmdb"     -> OUTPUT_LIGANDS_DIR/<HMDB_LIBRARY_SUBDIR>
      - TEST_MODE_ENABLE="hmdb+dud" -> OUTPUT_LIGANDS_DIR/<HMDB_LIBRARY_SUBDIR> + OUTPUT_LIGANDS_DIR/<mapped_subdir>
      - TEST_MODE_ENABLE="hmdb+fda" -> OUTPUT_LIGANDS_DIR/<HMDB_LIBRARY_SUBDIR> + OUTPUT_LIGANDS_DIR/<LIBRARY_SUBDIR_DEFAULT>
      - TEST_MODE_ENABLE="fda+dud+hmdb" -> OUTPUT_LIGANDS_DIR/<mapped_subdir> + OUTPUT_LIGANDS_DIR/<HMDB_LIBRARY_SUBDIR> + OUTPUT_LIGANDS_DIR/<LIBRARY_SUBDIR_DEFAULT>

    run_mode:
      - None: preserve TEST_MODE_ENABLE semantics (off / dud / fda+dud / hmdb / combinations)
      - "dud": force DUD-only by using the mapped subdir only
      - "fda": force FDA-only by using the default library only
      - "hmdb": force HMDB-only by using HMDB_LIBRARY_SUBDIR only

    Controls are *never* filtered out here.
    LIBRARY_EXTRA_DIRS remain included (unchanged).
    """
    overall_mode = _resolve_test_mode(cfg)
    if run_mode == "dud":
        effective_mode = "dud"
    elif run_mode == "fda":
        effective_mode = "off"
    elif run_mode == "hmdb":
        effective_mode = "hmdb"
    else:
        effective_mode = overall_mode

    # Keep existing prep step for extracted controls (harmless if nothing to do)
    prep_ligands_from_pdb(
        ligand_output_dir=paths.ligand_output_dir,
        ligands_mol2_dir=paths.ligands_mol2_dir,
        prepped_ligands_dir=paths.prepped_ligands_dir,
    )

    _control_roots, allowed_noncontrol_roots = _lib_roots_for_pdb(
        cfg,
        paths.pdb_id.upper(),
        paths,
        logger,
        test_mode_override=effective_mode,
    )
    per_index_roots: list[Path] = []
    if paths.prepped_ligands_dir:
        per_index_roots.append(paths.prepped_ligands_dir)

    index_library_roots: list[Path] = _dedup_index_roots(allowed_noncontrol_roots)
    per_index_roots = _dedup_index_roots(per_index_roots)
    index_roots = _dedup_index_roots(per_index_roots + index_library_roots)

    manifest_filename = str(cfg.get("LIBRARY_MANIFEST_FILENAME", "_manifest.json"))
    lib_index = cfg.get("_LIB_INDEX")
    if not isinstance(lib_index, LibraryIndex):
        lib_index = LibraryIndex(manifest_filename=manifest_filename, logger=logger)
        cfg["_LIB_INDEX"] = lib_index
    if index_roots:
        lib_index.load(index_roots)
    cfg["_LIB_INDEX_PER_ROOTS"] = [str(p) for p in per_index_roots]
    cfg["_LIB_INDEX_LIBRARY_ROOTS"] = [str(p) for p in index_library_roots]

    def _under(p: Path, root: Path) -> bool:
        try:
            p.resolve().relative_to(root.resolve())
            return True
        except Exception:
            return False

    def _enumerate_noncontrol_candidates_via_index(
        cfg: Dict,
        allowed_roots: list[Path],
        logger: logging.Logger,
    ) -> list[Path]:
        resolved_roots = _dedup_index_roots([Path(r) for r in allowed_roots if r])
        resolved_existing = [r for r in resolved_roots if r.exists()]
        logger.info(
            "[lib-roots] non-control roots = %s",
            [str(p) for p in resolved_existing],
        )

        candidates: list[Path] = []
        microstate_roots = [r for r in resolved_existing if (r / "microstates.json").exists()]
        if microstate_roots:
            try:
                all_ms: list[Path] = []
                for root in microstate_roots:
                    ms_paths = enumerate_ligands_for_docking(
                        requested_ph_values=None,
                        cfg=cfg,
                        pdb_id=paths.pdb_id.upper(),
                        root_dir=root,
                        microstate_dedup=True,
                        force=False,
                    )
                    all_ms.extend(ms_paths)
                ms_filtered = [p for p in all_ms if any(_under(p, root) for root in resolved_existing)]
                candidates.extend(ms_filtered)
                logger.info(
                    "[lib-index.microstate] roots=%d ligands=%d",
                    len(microstate_roots),
                    len(ms_filtered),
                )
            except Exception as exc:
                logger.warning("[lib-index.microstate] error=%s", exc)

        manifest_candidates: list[Path] = []
        if isinstance(cfg.get("_LIB_INDEX"), LibraryIndex):
            index_obj: LibraryIndex = cfg["_LIB_INDEX"]
            for root in resolved_existing:
                manifest = getattr(index_obj, "_cache", {}).get(Path(root))
                if not manifest:
                    continue
                for rel in manifest.entries.values():
                    manifest_candidates.append(Path(root) / rel)
        if manifest_candidates:
            candidates.extend(manifest_candidates)
            logger.info(
                "[lib-index.manifest] roots=%d ligands=%d",
                len(resolved_existing),
                len(manifest_candidates),
            )

        if not candidates:
            logger.warning(
                "[lib-index] no usable index detected; falling back to filesystem scan under %d roots (this is slow)",
                len(resolved_existing),
            )
            seen: set[str] = set()
            for root in resolved_existing:
                for p in _iter_pdbqt_dirfirst(root):
                    pn = norm(p)
                    if pn not in seen:
                        seen.add(pn)
                        candidates.append(p)

        deduped: list[Path] = []
        seen_keys: set[str] = set()
        for p in candidates:
            key = norm(p)
            if key not in seen_keys:
                seen_keys.add(key)
                deduped.append(Path(p))
        return deduped

    scan_roots: list[Path] = []
    prepped_lig_root = paths.prepped_ligands_dir
    if prepped_lig_root and prepped_lig_root.exists():
        scan_roots.append(prepped_lig_root)

    scan_roots.extend(allowed_noncontrol_roots)
    scan_roots = _dedup_index_roots(scan_roots)

    single_selector = cfg.get("_EFFECTIVE_SINGLE_LIGAND")
    crawl_allowed = not bool(single_selector)
    logger.info(
        "[ligands.scan.guard] roots=%d crawl_allowed=%s selector=%s",
        len(scan_roots),
        crawl_allowed,
        single_selector,
    )

    seen: set[str] = set()
    controls: list[Path] = []
    per_protein_allow = {"controls", "reference"}
    if prepped_lig_root and prepped_lig_root.exists():
        logger.info(
            "[ligands.scan.root] root=%s allowed_subdirs=%s crawl=%s",
            prepped_lig_root,
            per_protein_allow,
            crawl_allowed,
        )
        for p in _iter_pdbqt_dirfirst(prepped_lig_root, allowed_subdirs=per_protein_allow):
            pn = norm(p)
            if pn not in seen:
                seen.add(pn)
                controls.append(p)

    noncontrol_candidates = _enumerate_noncontrol_candidates_via_index(cfg, allowed_noncontrol_roots, logger)

    all_pdbqt_paths: list[Path] = []
    all_pdbqt_paths.extend(controls)
    all_pdbqt_paths.extend(noncontrol_candidates)

    logger.info("[ligands.scan] roots=%d found=%d", len(scan_roots), len(all_pdbqt_paths))
    cfg["ALL_LIGAND_PATHS"] = [str(p) for p in all_pdbqt_paths]

    if not all_pdbqt_paths:
        logger.warning("No .pdbqt ligands were found under the configured roots.")
        return [], {}, {}

    global_root = Path(cfg["OUTPUT_LIGANDS_DIR"]) if cfg.get("OUTPUT_LIGANDS_DIR") else None

    valid_pdbqt: Dict[str, Path] = {}
    for p in all_pdbqt_paths:
        try:
            if p.exists() and p.stat().st_size > 100:
                valid_pdbqt[norm(p)] = p
            else:
                logger.debug(f"Excluded malformed ligand (missing or too small): {p}")
        except Exception:
            logger.debug(f"Excluded malformed ligand (exception): {p}")

    cfg["ALL_LIGAND_PATHS_VALID"] = [str(p) for p in valid_pdbqt.values()]
    logger.info("[ligands.valid] count=%d", len(valid_pdbqt))

    controls_valid: list[Path] = []
    noncontrols: list[Path] = []
    for p in valid_pdbqt.values():
        if prepped_lig_root and _under(p, prepped_lig_root):
            controls_valid.append(p)
        else:
            noncontrols.append(p)

    filtered_noncontrols: list[Path] = []
    for p in noncontrols:
        keep = False
        for root in allowed_noncontrol_roots:
            if root.exists() and _under(p, root):
                keep = True
                break
        if keep:
            filtered_noncontrols.append(p)

    # When pH ligand mode is disabled, drop pH-annotated microstate ligands
    # from the non-control pool. We detect these by 'pH' in the filename stem.
    ph_mode_raw = str(cfg.get("PH_LIGAND_MODE", "off")).strip().lower()
    ph_ligand_mode_on = ph_mode_raw not in ("", "off", "none", "false", "0")
    if not ph_ligand_mode_on:
        before = len(filtered_noncontrols)
        filtered_noncontrols = [
            p for p in filtered_noncontrols
            if "pH" not in p.stem
            and not any(part.lower() == "microstates" for part in p.parts)
        ]
        removed = before - len(filtered_noncontrols)
        if removed:
            logger.info(
                "[ligands.ph-filter] PH_LIGAND_MODE=False -> removed %d ligands with 'pH' in stem or under microstates/ from non-control pool",
                removed,
            )

    # --- Optional: build library manifests from scan results ---
    if cfg.get("LIBRARY_MANIFEST_BUILD_ON_SCAN"):
        try:
            manifest_filename = str(cfg.get("LIBRARY_MANIFEST_FILENAME", "_manifest.json"))
            by_root: Dict[Path, list[Path]] = {}

            for root in allowed_noncontrol_roots:
                root = Path(root)
                if not root.exists():
                    continue
                for lig in filtered_noncontrols:
                    lig_path = Path(lig)
                    if not lig_path.exists():
                        continue
                    if not _under(lig_path, root):
                        continue
                    by_root.setdefault(root, []).append(lig_path)

            for root, ligs in by_root.items():
                manifest_path = root / manifest_filename
                if manifest_path.exists():
                    logger.info(
                        "[lib-manifest.scan.skip] root=%s reason=exists path=%s",
                        str(root),
                        str(manifest_path),
                    )
                    continue

                if not ligs:
                    continue

                logger.info(
                    "[lib-manifest.scan.build] root=%s ligands=%d manifest=%s",
                    str(root),
                    len(ligs),
                    str(manifest_path),
                )

                entries: list[str] = []
                for lig in ligs:
                    try:
                        rel = Path(lig).resolve().relative_to(root.resolve())
                        entries.append(rel.as_posix())
                    except Exception:
                        entries.append(os.path.relpath(str(lig), str(root)))

                tmp_path = manifest_path.with_suffix(manifest_path.suffix + ".tmp")

                try:
                    lib_index = LibraryIndex(
                        manifest_filename=manifest_filename,
                        logger=logger,
                    )
                    if hasattr(lib_index, "write_manifest_for_root"):
                        lib_index.write_manifest_for_root(root, entries, tmp_path)
                    else:
                        import json

                        tmp_path.parent.mkdir(parents=True, exist_ok=True)
                        data = {"root": str(root), "entries": entries}
                        tmp_path.write_text(json.dumps(data, indent=2, sort_keys=True))

                    tmp_path.replace(manifest_path)
                except Exception:
                    logger.exception(
                        "[lib-manifest.scan.error] root=%s manifest=%s",
                        str(root),
                        str(manifest_path),
                    )
        except Exception:
            logger.exception("[lib-manifest.scan.error] unexpected failure during build_on_scan")
    # ------------------------------------------------------------

    # Merge back: controls (unaltered) + filtered non-controls
    if cfg.get("_EFFECTIVE_SINGLE_LIGAND") and cfg.get("_SINGLE_RESOLVED_PATH"):
        resolved_path = Path(cfg["_SINGLE_RESOLVED_PATH"])
        final_paths = controls_valid + [resolved_path]
        filtered_noncontrols = [resolved_path]
        logger.info("[single.fuel] resolved=%s controls=%d (blocking non-control pool)", cfg["_SINGLE_RESOLVED_PATH"], len(controls_valid))
    else:
        final_paths = controls_valid + filtered_noncontrols

    # --- PAINS flags (keep as before; default to {}) ---
    pains_flags: Dict[str, bool] = {}
    try:
        from rdkit import Chem
        from rdkit.Chem import FilterCatalog, rdMolStandardize
        # Build catalog once at module-level if you prefer; safe inline here too
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
        pains_catalog = FilterCatalog.FilterCatalog(params)

        def _has_pains(pdbqt_path: Path) -> bool:
            try:
                # Try to locate a mol2 (or similar) neighbor if your original logic requires it.
                # Fallback: False (non-blocking).
                return False
            except Exception:
                return False

        for p in final_paths:
            pains_flags[p.stem] = _has_pains(p)
    except Exception:
        pains_flags = {}

    # Heavy atom counts (reuse your existing helper)
    heavy_atom_counts: Dict[str, int] = {}
    for p in final_paths:
        try:
            heavy_atom_counts[p.stem] = _count_heavy_atoms_from_pdbqt(p)
        except Exception:
            heavy_atom_counts[p.stem] = 0

    # Final return (stringify paths)
    ligands = [str(p) for p in final_paths]
    non_control_count = max(0, len(final_paths) - len(controls))
    logger.info(f"Selected ligands -> controls={len(controls)} + non-controls={non_control_count} = total={len(ligands)}")
    # GUARD: enforce .pdbqt-only pool
    bad = [p for p in ligands if not str(p).lower().endswith(".pdbqt")]
    if bad:
        raise ValueError(f"Ligand is not a .pdbqt file: {bad[0]}")
    return ligands, heavy_atom_counts, pains_flags


def _read_any_lig(path: str):
    """
    Load ligand from SDF/MOL2/PDB with consistent settings.
    Returns an RDKit Mol or None.
    """
    mol = None
    loader = "unknown"
    sanitize = True
    sanitize_failed = False
    load_err = None
    ext = os.path.splitext(path)[1].lower()
    _rlog = logging.getLogger("rmsd")

    try:
        if ext in (".sdf", ".sd"):
            loader = "SDMolSupplier"
            sanitize = True
            suppl = Chem.SDMolSupplier(path, removeHs=False, sanitize=True)
            mol = next((m for m in suppl if m is not None), None)
        elif ext in (".mol2",):
            loader = "MolFromMol2File"
            sanitize = True
            mol = Chem.MolFromMol2File(path, sanitize=True, removeHs=False)
        elif ext in (".pdbqt",):
            loader = "MolFromPDBFile(pdbqt)"
            sanitize = False
            mol = Chem.MolFromPDBFile(path, sanitize=False, removeHs=False, proximityBonding=True)
            if mol is not None:
                try:
                    Chem.SanitizeMol(mol)
                except Exception as e:
                    sanitize_failed = True
                    if _rlog:
                        _rlog.warning(f"[read_any] sanitize failed for path='{path}' err={e!r}")
        elif ext in (".pdb",):
            loader = "MolFromPDBFile"
            sanitize = True
            # If you use proximityBonding or flavor flags elsewhere, keep them consistent here.
            try:
                mol = Chem.MolFromPDBFile(path, sanitize=True, removeHs=False, proximityBonding=True)
            except Exception as e:
                load_err = e
                mol = None
            if mol is None:
                try:
                    with open(path, "rt", errors="ignore") as fh:
                        head = fh.read(1024)
                except Exception:
                    head = ""
                looks_like_pdbqt = ("REMARK VINA" in head) or ("TORSDOF" in head) or ("ROOT" in head)
                if looks_like_pdbqt:
                    if _rlog:
                        _rlog.warning(
                            f"[read_any] PDBQT fallback for path='{path}' (sanitize=True failed; retry sanitize=False)"
                        )
                    loader = "MolFromPDBFile(pdbqt-fallback)"
                    sanitize = False
                    try:
                        mol = Chem.MolFromPDBFile(path, sanitize=False, removeHs=False, proximityBonding=True)
                    except Exception as e:
                        load_err = e
                        mol = None
                    if mol is not None:
                        try:
                            Chem.SanitizeMol(mol)
                        except Exception as e:
                            sanitize_failed = True
                            if _rlog:
                                _rlog.warning(f"[read_any] sanitize failed for path='{path}' err={e!r}")
        else:
            loader = "auto"
            sanitize = True
            mol = Chem.MolFromMolFile(path, sanitize=True, removeHs=False)  # last-ditch; or return None
    except Exception as e:
        load_err = e if load_err is None else load_err
        mol = None
    if (mol is None) and (load_err is None):
        load_err = "load_returned_None"

    # ──  single debug line about what we actually loaded ───────────────────
    try:
        _rlog = logging.getLogger("rmsd")
        if _rlog and mol is not None:
            from rdkit.Chem import rdMolDescriptors
            # formula = e.g., "C20H25N3O"
            formula = rdMolDescriptors.CalcMolFormula(mol)
            # InChIKey may be unavailable if RDKit was built without InChI; guard it.
            try:
                from rdkit.Chem import inchi
                inchikey = inchi.MolToInchiKey(mol)
            except Exception:
                inchikey = "NA"
            extra = " sanitize_failed=True" if sanitize_failed else ""
            _rlog.info(f"[read_any] loader={loader} sanitize={sanitize} path='{path}'{extra} "
                       f"atoms={mol.GetNumAtoms()} heavy={mol.GetNumHeavyAtoms()} formula={formula} inchikey={inchikey}")
        elif _rlog:
            if (mol is None) and (load_err is not None):
                _rlog.info(f"[read_any] loader={loader} sanitize={sanitize} path='{path}' mol=None err={load_err!r}")
            else:
                _rlog.info(f"[read_any] loader={loader} sanitize={sanitize} path='{path}' mol=None")
    except Exception:
        pass
    # ───────────────────────────────────────────────────────────────────────────

    return mol


def _is_readable_ref(pth: Path) -> bool:
    try:
        m = _read_any_lig(str(pth))
        return (m is not None) and (m.GetNumHeavyAtoms() > 0)
    except Exception:
        return False


def compute_rmsd(ref_path: str, docked_path: str) -> float:
    """Heavy-atom RMSD using best mapping; supports PDB/SDF/MOL2 refs and adds a minimal MCS fallback."""
    _rlog = logging.getLogger("rmsd")
    ref = _read_any_lig(ref_path)
    dock = _read_any_lig(docked_path)
    if not ref or not dock:
        return float("inf")

    # 1) Fast path: RDKit best alignment
    try:
        rmsd = float(rdMolAlign.GetBestRMS(ref, dock))
        if _rlog:
            _rlog.info(f"[rmsd.best] ref='{ref_path}' dock='{docked_path}' rmsd={rmsd:.3f}A")
        return rmsd
    except Exception as e:
        if _rlog:
            _rlog.warning(
                f"[rmsd.best] failed ref='{ref_path}' dock='{docked_path}' err={e!r}; falling back to MCS"
            )

    # 2) Tiny, robust fallback via MCS
    try:
        # Ensure ring info for potentially unsanitized PDBQT-derived mols
        Chem.FastFindRings(ref)
        Chem.FastFindRings(dock)

        mcs = rdFMCS.FindMCS([ref, dock],
                             ringMatchesRingOnly=True,
                             completeRingsOnly=True,
                             matchValences=True)
        patt = Chem.MolFromSmarts(mcs.smartsString)
        if patt is None:
            if _rlog:
                _rlog.warning(f"[rmsd.mcs] no SMARTS pattern for ref='{ref_path}' dock='{docked_path}'")
            return float("inf")
        ref_match = ref.GetSubstructMatch(patt)
        dock_match = dock.GetSubstructMatch(patt)
        if not ref_match or not dock_match or (len(ref_match) != len(dock_match)):
            if _rlog:
                _rlog.warning(f"[rmsd.mcs] match-fail ref='{ref_path}' dock='{docked_path}'")
            return float("inf")
        amap = list(zip(dock_match, ref_match))  # (probe->ref)
        rmsd = float(rdMolAlign.AlignMol(dock, ref, atomMap=amap))
        if _rlog:
            _rlog.info(f"[rmsd.mcs] ref='{ref_path}' dock='{docked_path}' rmsd={rmsd:.3f}A atoms={len(amap)}")
        return rmsd
    except Exception as e:
        if _rlog:
            _rlog.warning(f"[rmsd.mcs] failed ref='{ref_path}' dock='{docked_path}' err={e!r}")
        return float("inf")


def validate_ligand(
        ligand_name: str,
        docked_path: str,
        crystal_path: str = None,
        rmsd_thresh: float = 2.0,
        self_rmsd: float = None,
        logger=None
) -> bool:
    """
    Validate ligand docking.
      * If crystal structure available ? use redocking RMSD (hard gate).
      * Otherwise (non-controls) ? self-RMSD is *log-only* (never reject).
    """
    if crystal_path and Path(crystal_path).exists():
        redock_rmsd = compute_rmsd(crystal_path, docked_path)
        if logger:
            sr = f"{self_rmsd:.2f}" if isinstance(self_rmsd, (int, float)) else "n/a"
            logger.info(f"[validate] {ligand_name}: redock_RMSD={redock_rmsd:.2f} A, self_RMSD={sr}")
        if redock_rmsd <= rmsd_thresh:
            return True
        else:
            if logger:
                logger.warning(
                    f"[validate] {ligand_name}: redocking failed (RMSD {redock_rmsd:.2f} A > {rmsd_thresh:.2f})"
                )
            return False

    # Non-controls: log self-RMSD but do not gate on it
    try:
        sr_val = float(self_rmsd) if self_rmsd is not None else None
    except Exception:
        sr_val = None
    if logger:
        sr_txt = f"{sr_val:.2f}" if isinstance(sr_val, (int, float)) else "n/a"
        logger.info(f"[validate] {ligand_name}: self_RMSD={sr_txt} A (LOG-ONLY)")
    return True


__all__ = [
    "select_ligands_for_next",
    "prepare_and_filter_ligands",
    "compute_rmsd",
    "validate_ligand",
    "_coerce_test_map",
    "_resolve_test_mode",
    "_lib_roots_for_pdb",
    "_dedup_index_roots",
    "_count_heavy_atoms_from_pdbqt",
    "_read_any_lig",
    "_is_readable_ref",
]
