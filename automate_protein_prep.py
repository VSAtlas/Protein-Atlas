# -*- coding: utf-8 -*-
"""
Automated protein preparation pipeline (Windows/WSL-friendly)

This module prepares a receptor structure from a raw PDB by performing:
  1) Alternate-conformation filtering (altLoc) with deterministic selection
  2) Removal of nonstandard residues (policy-aware keep/drop of cofactors/metals/waters)
  3) Element-column fixes to ensure PDB compliance (via activesite/YAML)
  4) Optional loop/residue completion via MODELLER
  5) Sanity checks for incomplete residues
  6) Phenix cleaning passes (auto-detect Linux Phenix, Windows Phenix, or skip—WSL-safe)
  7) Hydrogen sanity cleanup (CONECT- and geometry-based)
  8) Chain validation (ensure at least one CA-containing chain)
  9) Protonation via Reduce with automatic fallbacks (temp retry/OpenBabel)
 10) Final polish via phenix.pdbtools when available; otherwise degrade gracefully
 11) Receptor PDBQT preparation (Meeko preferred; ADT fallback)

Design notes
------------
• Single source of truth for element handling: **activesite** helpers (YAML-driven). No local element inference.
• WSL-aware external calls.
• Logging instead of prints; short, numbered steps for easier debugging.
• Never touch PDBQT with any PDB element-fixing logic.

Prerequisites (recommended)
---------------------------
• Phenix (Linux: phenix.pdbtools; Windows: phenix.pdbtools.bat, phenix.python.bat)
• Reduce (reduce.exe or reduce) — often bundled with Phenix
• MODELLER (licensed) — optional but recommended for loop filling
• Open Babel (obabel) — fallback for hydrogen addition
• Biopython — for PDB parsing
• MGLTools (for prepare_receptor4.py) and/or Meeko (mk_prepare_receptor)
"""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import (
    Collection,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
    Union,
)
import os
import shutil
import logging

# >>> PATHS IMPORT START
from path_router import make_paths

# >>> PATHS IMPORT END
# NOTE:
# These imports intentionally expose a broad legacy surface for helper modules
# that hydrate globals from automate_protein_prep.
from protein_prep.altloc_filter import filter_altlocs  # noqa: F401
from protein_prep.aliases_policy import (  # noqa: F401
    ALIASES,
    RULES,
    _ALIASES_BIND_LOGGED,
    _COFACTOR_CANONICAL,
    _COFACTOR_DROP_LOGGED,
    _COFACTOR_NAMES,
    _COFACTOR_RAW_ALL,
    _ELEM_CANON,
    _ELEMENT_TOKENS_RAW,
    _ION_AUDIT_ALIAS_MAP,
    _ION_AUDIT_DISABLED_VALUES,
    _ION_AUDIT_ENABLED_VALUES,
    _ION_AUDIT_METALS,
    _ION_AUDIT_SIMPLE_IONS,
    _ION_AUDIT_WATERS,
    _ION_BREADCRUMB_METAL_ORDER,
    _ION_BREADCRUMB_SIMPLE_ORDER,
    _ION_KEEP_LOGGED,
    _IONS_CFG_CACHE,
    _IONS_CFG_LOGGED,
    _POLICY_MODE,
    _RETAIN_VARIANT,
    _RETAIN_VARIANT_CANONICAL,
    _WATER_NAMES,
    _flatten_semicolons,
    _normalize_resname,
    _to_upper_set,
    load_aliases as _aliases_policy_load_aliases,
)
from protein_prep.element_guard import (  # noqa: F401
    _meeko_preflight_or_fail,
    _post_write_element_guard,
)
from protein_prep.os_utils import (  # noqa: F401
    _as_path,
    _powershell,
    _run_and_log,
    _win_path,
)
from protein_prep.prep_utils import (  # noqa: F401
    _cfg,
    _cfg_bool,
    _cfg_chain_keep_list,
    _cfg_float,
    _cfg_int,
    _load_retain_allowlist,
    _persist_subproc,
    _resolve_variant_token,
    _short_path_for_log,
    set_runtime_config,
)
from protein_prep.protonation import (  # noqa: F401
    REDUCE_EXE,
    _protonate_with_pdb2pqr_if_available,
    assign_protonation_states,
    conect_coverage,
    hydrogenation_status,
)
from protein_prep.receptor_prep import (  # noqa: F401
    _MEEKO_DROP_IONS,
    _classify_and_rename_histidines,
    _clean_receptor_pdbqt,
    _distance_from_center,
    _drop_free_ions_for_meeko,
    _first_last_lines,
    _ion_pairs_from_pdb,
    _ion_pairs_from_pdbqt,
    _log_ion_diff,
    _log_pdb_pdbqt_counts_diff,
    _normalize_ion_policy,
    _ok_receptor_file,
    _rewrite_his_default,
    _strip_reduce_user_lines,
    run_prepare_receptor,
)
from protein_prep.tool_runners import (  # noqa: F401
    build_adt_prepare_receptor_cmd,
    build_meeko_base_cmd,
    build_meeko_legacy_cmd,
    build_meeko_modern_cmd,
    run_subprocess_capture,
)

from installation import load_config
import pdb_fixer as _activesite_mod  # noqa: F401
from pdb_fixer import (  # noqa: F401
    assert_no_helium_in_hydrogen_names,
    fix_pdb_elements,
    fix_element_columns_in_file,
    get_atom_rules,
    load_canonical_cofactors,
    load_canonical_metals,
    load_canonical_waters,
    scan_helium_counts,
)
_HE_POSTWRITE_VERBOSE = os.environ.get("HELIUM_POSTWRITE_VERBOSE", "1") != "0"


def load_aliases():
    return _aliases_policy_load_aliases()


_ION_PIPE_AUDIT: dict = {}
_ION_PIPE_WARNED = False


def _format_histogram(counter: Counter[str]) -> str:
    from protein_prep.ion_audit import _format_histogram as _format_histogram_impl

    return _format_histogram_impl(counter)



def _format_token_list(tokens: Iterable[str], *, limit: int = 8) -> str:
    from protein_prep.ion_audit import _format_token_list as _format_token_list_impl

    return _format_token_list_impl(tokens, limit=limit)



def _summarize_ions_file(file_path: Union[str, Path]) -> dict[str, object]:
    from protein_prep.ion_audit import _summarize_ions_file as _summarize_ions_file_impl

    return _summarize_ions_file_impl(file_path)



def _format_ion_pairs(pairs: Iterable[tuple[str, str]]) -> str:
    seq = sorted(pairs)
    if not seq:
        return "none"
    return ",".join(f"{res}:{loc}" for res, loc in seq)


def _ion_audit_enabled() -> bool:
    from protein_prep.ion_audit import _ion_audit_enabled as _ion_audit_enabled_impl

    return _ion_audit_enabled_impl()



def _canon_ion_resname(resname: str) -> str:
    from protein_prep.ion_audit import _canon_ion_resname as _canon_ion_resname_impl

    return _canon_ion_resname_impl(resname)



def _format_breadcrumb_counts(counts: dict[str, int], order: Sequence[str]) -> str:
    from protein_prep.ion_audit import (
        _format_breadcrumb_counts as _format_breadcrumb_counts_impl,
    )

    return _format_breadcrumb_counts_impl(counts, order)



def _breadcrumbs_enabled() -> bool:
    from protein_prep.ion_audit import _breadcrumbs_enabled as _breadcrumbs_enabled_impl

    return _breadcrumbs_enabled_impl()



def _emit_ion_breadcrumb(stage: str, file_path: Union[str, Path]) -> None:
    from protein_prep.ion_audit import _emit_ion_breadcrumb as _emit_ion_breadcrumb_impl

    _emit_ion_breadcrumb_impl(stage, file_path)



def _serialize_counts(counts: dict[str, int]) -> str:
    from protein_prep.ion_audit import _serialize_counts as _serialize_counts_impl

    return _serialize_counts_impl(counts)



def _legacy_ion_global_missing() -> None:
    from protein_prep.ion_audit import (
        _legacy_ion_global_missing as _legacy_ion_global_missing_impl,
    )

    _legacy_ion_global_missing_impl()



def _gather_ion_counts(path: Path) -> tuple[dict[str, int], dict[str, int], int, int]:
    from protein_prep.ion_audit import _gather_ion_counts as _gather_ion_counts_impl

    return _gather_ion_counts_impl(path)



def audit_ions(
    pdb_path: Union[str, Path],
    pdb_id: str,
    stage: str,
    variant: Optional[str],
    logger: logging.Logger | None = None,
) -> Optional[dict[str, object]]:
    from protein_prep.ion_audit import audit_ions as _audit_ions_impl

    return _audit_ions_impl(
        pdb_path=pdb_path,
        pdb_id=pdb_id,
        stage=stage,
        variant=variant,
        logger=logger,
    )



def diff_ions(
    prev: dict[str, object],
    curr: dict[str, object],
    logger: logging.Logger | None = None,
) -> dict[str, object]:
    from protein_prep.ion_audit import diff_ions as _diff_ions_impl

    return _diff_ions_impl(prev=prev, curr=curr, logger=logger)



class _IonAuditManager:
    def __new__(cls, *args, **kwargs):
        from protein_prep.ion_audit import _IonAuditManager as _IonAuditManagerImpl

        return _IonAuditManagerImpl(*args, **kwargs)



_ION_AUDIT_STACK: list[_IonAuditManager] = []


def _push_ion_audit_manager(manager: _IonAuditManager) -> None:
    from protein_prep.ion_audit import _push_ion_audit_manager as _push_impl

    _push_impl(manager)



def _pop_ion_audit_manager(manager: _IonAuditManager) -> None:
    from protein_prep.ion_audit import _pop_ion_audit_manager as _pop_impl

    _pop_impl(manager)



def _current_ion_audit_manager() -> Optional[_IonAuditManager]:
    from protein_prep.ion_audit import _current_ion_audit_manager as _current_impl

    return _current_impl()



def emit_ion_audit_probe(
    stage: str, file_path: Union[str, Path], *, variant: Optional[str] = None
) -> None:
    from protein_prep.ion_audit import emit_ion_audit_probe as _emit_probe_impl

    _emit_probe_impl(stage=stage, file_path=file_path, variant=variant)



_LAST_CLEAN_PROVENANCE = "unknown"


def get_clean_provenance() -> str:
    return _LAST_CLEAN_PROVENANCE


def _set_clean_provenance(label: str) -> None:
    global _LAST_CLEAN_PROVENANCE
    _LAST_CLEAN_PROVENANCE = label or "unknown"


def _reset_ion_probe(pdb_id: str) -> None:
    from protein_prep.ion_audit import _reset_ion_probe as _reset_probe_impl

    _reset_probe_impl(pdb_id)



def _ion_candidate_tokens(rules_obj=ALIASES) -> set[str]:
    tokens: set[str] = set()
    canonical = getattr(rules_obj, "elem_tokens_canonical", None) or set()
    for tok in canonical:
        if tok is None:
            continue
        text = str(tok).strip()
        if text:
            tokens.add(text.upper())
    alias_sets = getattr(rules_obj, "alias_sets", None)
    if alias_sets and getattr(alias_sets, "elem_tokens_canonical", None):
        for tok in getattr(alias_sets, "elem_tokens_canonical", set()):
            if tok is None:
                continue
            text = str(tok).strip()
            if text:
                tokens.add(text.upper())
    return tokens


def _scan_metal_map(path: Union[str, Path], rules=ALIASES) -> dict[str, int]:
    file_path = Path(path)
    if not file_path.exists():
        return {}
    tokens = _ion_candidate_tokens(rules)
    counts: Counter[str] = Counter()
    try:
        with file_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.startswith(("HETATM", "ATOM  ")):
                    continue
                resname = line[17:20].strip().upper()
                if not resname:
                    continue
                if resname in tokens:
                    counts[resname] += 1
    except Exception as exc:
        logging.warning("[ions.probe] stage=scan_fail file=%s err=%s", file_path, exc)
        return dict(counts)
    return dict(counts)


def _format_probe_counts(counts: dict[str, int], *, limit: int | None = None) -> str:
    if not counts:
        return "none"
    items = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    if limit is not None:
        items = items[:limit]
    return ",".join(f"{res}:{cnt}" for res, cnt in items)


def _log_ions_probe(
    pdb_id: str,
    stage: str,
    file_path: Union[str, Path],
    *,
    phase: str | None = None,
) -> dict[str, int]:
    from protein_prep.ion_audit import _log_ions_probe as _log_probe_impl

    return _log_probe_impl(
        pdb_id=pdb_id,
        stage=stage,
        file_path=file_path,
        phase=phase,
    )



def _format_diff_map(data: dict[str, int]) -> str:
    if not data:
        return "none"
    filtered = [(res, count) for res, count in data.items() if count > 0]
    if not filtered:
        return "none"
    items = sorted(filtered, key=lambda kv: (-kv[1], kv[0]))
    return ",".join(f"{res}:{cnt}" for res, cnt in items)


def get_ion_probe_map(pdb_id: str) -> dict[str, dict[str, int]]:
    from protein_prep.ion_audit import get_ion_probe_map as _get_probe_map_impl

    return _get_probe_map_impl(pdb_id)



def _bucket_counts(counter: Counter) -> dict[str, int]:
    keys = ["ZN", "MG", "NA", "K", "CA", "MN", "FE", "CL"]
    out = {k: 0 for k in keys}
    other = 0
    for resn, count in counter.items():
        token = resn.upper()
        if token in out:
            out[token] += count
        else:
            other += count
    out["OTHER"] = other
    return out


def _format_counts(counter: Counter) -> str:
    bucketed = _bucket_counts(counter)
    keys = ["ZN", "MG", "NA", "K", "CA", "MN", "FE", "CL", "OTHER"]
    parts = [f"{k}={bucketed.get(k, 0)}" for k in keys]
    return " ".join(parts)


# [ions] monoatomic audit helpers
def _collect_monoatomic_records(pdb_path: Union[str, Path]) -> tuple[Counter, Counter]:
    from protein_prep.ion_audit import _collect_monoatomic_records as _collect_impl

    return _collect_impl(pdb_path)



def _format_ion_hist(counter: Counter) -> str:
    from protein_prep.ion_audit import _format_ion_hist as _format_hist_impl

    return _format_hist_impl(counter)



def _diff_detail_records(before: Counter, after: Counter) -> list[str]:
    from protein_prep.ion_audit import _diff_detail_records as _diff_detail_impl

    return _diff_detail_impl(before, after)



def _format_ion_pairs(pairs: Iterable[tuple[str, str]]) -> str:
    sorted_pairs = sorted(pairs)
    if not sorted_pairs:
        return "none"
    return ",".join(f"{res}:{loc}" for res, loc in sorted_pairs)


def _parse_atoms_from_pdb_like_lines(
    lines: list[str], canonical_metals: Collection[str] | None = None
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    from protein_prep.metal_site_audit import (
        _parse_atoms_from_pdb_like_lines as _parse_atoms_impl,
    )

    return _parse_atoms_impl(lines, canonical_metals)



def _find_metal_donors(
    metal: Mapping[str, object],
    parsed_atoms: Sequence[Mapping[str, object]],
    canonical_waters: Collection[str] | None = None,
) -> list[dict[str, object]]:
    from protein_prep.metal_site_audit import _find_metal_donors as _find_donors_impl

    return _find_donors_impl(metal, parsed_atoms, canonical_waters)



def run_metal_site_audit(
    *,
    pdb_id: str,
    router_paths,
    input_pdb_path: str,
    receptor_pdb_path: str | None,
    receptor_pdbqt_path: str | None,
    center: Optional[tuple[float, float, float]] = None,
    variant_label: str | None = None,
    ph_label: str | None = None,
) -> None:
    from protein_prep.metal_site_audit import run_metal_site_audit as _metal_audit_impl

    _metal_audit_impl(
        pdb_id=pdb_id,
        router_paths=router_paths,
        input_pdb_path=input_pdb_path,
        receptor_pdb_path=receptor_pdb_path,
        receptor_pdbqt_path=receptor_pdbqt_path,
        center=center,
        variant_label=variant_label,
        ph_label=ph_label,
    )



def _maybe_strip_ions(
    pdb_path: Union[str, Path],
    cfg: Optional[dict] = None,
    *,
    variant: Optional[str] = None,
    pocket_center: Optional[tuple[float, float, float]] = None,
    extra_keep: Optional[Iterable[str]] = None,
) -> int:
    from protein_prep.strip_nsr import _maybe_strip_ions as _maybe_strip_ions_impl

    return _maybe_strip_ions_impl(
        pdb_path=pdb_path,
        cfg=cfg,
        variant=variant,
        pocket_center=pocket_center,
        extra_keep=extra_keep,
    )



def _is_element_token(sym):
    from protein_prep.strip_nsr import _is_element_token as _is_element_token_impl

    return _is_element_token_impl(sym)



# Treat as "ion-like" if it normalizes to a canonical element token retained by policy
def _is_retained_ion(resname: str) -> bool:
    from protein_prep.strip_nsr import _is_retained_ion as _is_retained_ion_impl

    return _is_retained_ion_impl(resname)


# =============================
# Configuration helpers
# =============================
# Single unified config load for this module:
# Priority: environment overrides > config.txt keys > legacy key aliases > defaults.
try:
    from input_and_export_functions import load_config, validate_config

    _CFG = load_config("config.txt")
    validate_config(_CFG)
except Exception:
    # Fallback for older setups
    from installation import load_config as _legacy_load_config

    _CFG = _legacy_load_config()

# Keep both names so existing call-sites continue to work
config = _CFG
_CONFIG = _CFG
set_runtime_config(_CFG)


PHENIX_DIR = _cfg("PHENIX_DIR", "", "phenix_dir")
PHENIX_LIB_PATH = _cfg("PHENIX_LIB_PATH", "", "phenix_lib_path")
PHENIX_CLEAN_SCRIPT = _cfg("PHENIX_CLEAN_SCRIPT", "", "phenix_clean_script")
INPUT_DIR = _cfg("INPUT_DIR", ".")

MGLTOOLS_PYTHON = _cfg("MGLTOOLS_PYTHON", "")
PREPARE_RECEPTOR_SCRIPT = _cfg("PREPARE_RECEPTOR_SCRIPT", "")
USE_MEEKO = str(_cfg("USE_MEEKO", "")).lower() in ("1", "true", "yes")

# Windows .bat wrappers (if using Windows Phenix)
PDBTOOLS_BAT = str(Path(PHENIX_DIR) / "phenix.pdbtools.bat") if PHENIX_DIR else ""
MOLPROBITY_BAT = str(Path(PHENIX_DIR) / "phenix.molprobity.bat") if PHENIX_DIR else ""
PHENIX_PYTHON_BAT = str(Path(PHENIX_DIR) / "phenix.python.bat") if PHENIX_DIR else ""


# --- HOLO-only restore helper ---
def _holo_restore_from_input_if_needed(
    pdb_id: str,
    cleaned_pdb: Union[str, Path],
    output_pdbqt: Union[str, Path],
    config: dict,
    center: Optional[tuple[float, float, float]] = None,
    box_size: Optional[tuple[float, float, float]] = None,
) -> tuple[int, int, bool]:
    from protein_prep.holo_restore import _holo_restore_from_input_if_needed as _holo_restore_impl

    return _holo_restore_impl(
        pdb_id=pdb_id,
        cleaned_pdb=cleaned_pdb,
        output_pdbqt=output_pdbqt,
        config=config,
        center=center,
        box_size=box_size,
    )


def _helium_postwrite_counter(step_name: str, pdb_path: str | Path) -> None:
    """
    Run the element fixer on the freshly written receptor PDB and emit:
      [helium] stage=post_write step=<name> file=<path> He->H=<n> residual_He=<n>
    Fail fast if residual_He > 0.
    """
    p = Path(pdb_path)
    try:
        before = p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        before = ""

    # Count BEFORE
    he_before = scan_helium_counts(before.splitlines()) if before else 0

    # Fix elements in-place on the exact file we pass forward
    fix_element_columns_in_file(p, p, rewrite_atoms=True)

    try:
        after = p.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        after = ""

    # Count AFTER and compute delta
    he_after = scan_helium_counts(after.splitlines()) if after else 0
    delta = max(0, he_before - he_after)

    if _HE_POSTWRITE_VERBOSE:
        logging.info(
            "[helium] stage=post_write step=%s file=%s He->H=%d residual_He=%d",
            step_name,
            str(p),
            delta,
            he_after,
        )

    if he_after > 0:
        try:
            bad = [
                ln
                for ln in after.splitlines()
                if (" He" in ln) or (len(ln) >= 78 and ln[76:78].strip() == "HE")
            ][:3]
            logging.warning("[helium] residual examples: %r", bad)
        except Exception:
            pass
        # Fail fast so we can see which step leaked helium
        raise RuntimeError("helium_residual_post_write")


try:

    _HAS_RDKIT = True
except Exception:
    _HAS_RDKIT = False

def collapse_sanitized_once(p: Union[str, Path]) -> Path:
    from protein_prep.ligand_extract import collapse_sanitized_once as _collapse_impl

    return _collapse_impl(p)



def compute_control_centroids(
    ligands_dir: Union[str, Path],
) -> list[tuple[float, float, float]]:
    from protein_prep.waters import compute_control_centroids as _compute_centroids_impl

    return _compute_centroids_impl(ligands_dir)



def filter_waters_near_points(
    src_pdb: Union[str, Path],
    dst_pdb: Union[str, Path],
    points: list[tuple[float, float, float]],
    radius_A: float,
) -> int:
    from protein_prep.waters import (
        filter_waters_near_points as _filter_waters_near_points_impl,
    )

    return _filter_waters_near_points_impl(
        src_pdb=src_pdb,
        dst_pdb=dst_pdb,
        points=points,
        radius_A=radius_A,
    )



def count_waters_within(
    pdb_path: Union[str, Path], point_xyz: tuple[float, float, float], radius_A: float
) -> int:
    """Count distinct HOH residues within radius_A of point_xyz."""
    R2 = float(radius_A) * float(radius_A)
    seen: set[tuple[str, str]] = set()
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if not ln.startswith(("ATOM  ", "HETATM")) or ln[17:20] != "HOH":
                    continue
                try:
                    x = float(ln[30:38])
                    y = float(ln[38:46])
                    z = float(ln[46:54])
                except Exception:
                    continue
                dx = x - point_xyz[0]
                dy = y - point_xyz[1]
                dz = z - point_xyz[2]
                if (dx * dx + dy * dy + dz * dz) <= R2:
                    seen.add((ln[21], ln[22:27]))
    except Exception:
        return 0
    return len(seen)


def _cfg_env_or_default(key: str, default: Optional[str] = None) -> Optional[str]:
    """Prefer env var, then normalized load_config(), else default."""
    v = os.environ.get(key)
    if v:
        return v
    try:
        root = Path(__file__).resolve().parent
        cfg = load_config(config_path=str(root / "config.txt"), base_dir=root)
        val = cfg.get(key)
        if val not in (None, ""):
            return str(val).strip()
    except Exception:
        pass
    return default


def _canon_base(output_root: Path, pdb_id: str) -> Path:
    from protein_prep.layout import _canon_base as _impl

    return _impl(output_root, pdb_id)


def _merge_dir(src: Path, dst: Path) -> None:
    from protein_prep.layout import _merge_dir as _impl

    _impl(src, dst)


def fold_legacy_layout(pdb_id: str, output_root) -> None:
    from protein_prep.layout import fold_legacy_layout as _impl

    _impl(pdb_id, output_root)


# =============================
# Ligand pristine reference (optional; for RMSD downstream)
# =============================


def _write_pristine_reference(pdb_lig_path: Path) -> None:
    from protein_prep.ligand_extract import (
        _write_pristine_reference as _write_pristine_reference_impl,
    )

    _write_pristine_reference_impl(pdb_lig_path)



# =============================
# Canonical per-protein directory layout
# =============================


def canon_paths(
    pdb_id: str,
    output_root: Union[str, Path],
    *,
    variant: Optional[str] = None,
) -> Dict[str, Path]:
    from protein_prep.layout import canon_paths as _impl

    return _impl(
        pdb_id=pdb_id,
        output_root=output_root,
        variant=variant,
    )


# =============================
# Ligand extraction (single source lives here)
# =============================


def extract_ligands_from_filtered(
    filtered_pdb: Union[str, Path], out_dir: Union[str, Path]
) -> List[Path]:
    from protein_prep.ligand_extract import (
        extract_ligands_from_filtered as _extract_ligands_impl,
    )

    return _extract_ligands_impl(filtered_pdb, out_dir)



def element_fix_all_in_dir(
    dir_path: Union[str, Path], rewrite_atoms: bool = False
) -> int:
    from protein_prep.ligand_extract import element_fix_all_in_dir as _fix_all_impl

    return _fix_all_impl(dir_path, rewrite_atoms=rewrite_atoms)



def expose_ligand_intermediates_for_debug(
    src_dir: Union[str, Path], link_dir: Union[str, Path]
) -> None:
    from protein_prep.ligand_extract import (
        expose_ligand_intermediates_for_debug as _expose_impl,
    )

    _expose_impl(src_dir, link_dir)



# =============================
# Element & Hydrogen utilities (PDB only)
# =============================


def file_contains_hydrogens(pdb_path: Union[str, Path]) -> bool:
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                el = line[76:78].strip().upper()
                if el in {"H", "D", "T"}:
                    return True
    except Exception as e:
        logging.warning("Could not read %s to check for H atoms: %s", pdb_path, e)
    return False


def strip_monoatomic_ions_inplace(
    pdb_path: Union[str, Path],
    keep_resnames: Optional[Set[str]] = None,
    *,
    cfg: Optional[dict] = None,
    variant: Optional[str] = None,
    pocket_center: Optional[tuple[float, float, float]] = None,
    force_policy: Optional[str] = None,
) -> int:
    """Legacy wrapper that now delegates to the policy-aware ion stripping."""

    base_cfg = cfg if cfg is not None else config
    if isinstance(base_cfg, dict):
        cfg_obj: dict = dict(base_cfg)
    else:
        cfg_obj = {}

    if force_policy:
        cfg_obj = cfg_obj or {}
        cfg_obj["ION_STRIP_POLICY"] = force_policy

    extra_keep = set(keep_resnames or []) or None

    removed = _maybe_strip_ions(
        pdb_path,
        cfg=cfg_obj,
        variant=variant,
        pocket_center=pocket_center,
        extra_keep=extra_keep,
    )

    if removed:
        logging.info("Stripped %d monoatomic ions from %s", removed, pdb_path)

    return removed or 0


# Compatibility alias for older callers
def _strip_monoatomic_ions_inplace(*args, **kwargs):
    return strip_monoatomic_ions_inplace(*args, **kwargs)


def detect_catalytic_metals(pdb_path: Union[str, Path]) -> Set[str]:
    """
    Return the set of retained ion-like resnames present in the PDB (e.g., ZN, MG).
    Uses YAML retain list + element tokens. For diagnostics only.
    """
    RETAIN = set(_RETAIN_VARIANT)

    ionic = set()
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if not ln.startswith("HETATM"):
                continue
            resname = ln[17:20].strip().upper()
            if (resname in RETAIN) and _is_element_token(resname):
                ionic.add(resname)

    if ionic:
        logging.info("Catalytic/retained ions present: %s", sorted(ionic))
    return ionic


def compare_ion_presence_between_pdb_and_pdbqt(
    clean_pdb: Union[str, Path], receptor_pdbqt: Union[str, Path]
) -> None:
    """
    Log a warning if an ion present in the cleaned PDB is missing in receptor PDBQT.
    Uses YAML retain list + element tokens as "ions" definition.
    """

    def ions_in_pdb(p: Union[str, Path]) -> Set[Tuple[str, str]]:
        ions = set()
        with open(p, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if not ln.startswith("HETATM"):
                    continue
                resname = ln[17:20].strip().upper()
                chain = ln[21]
                resi = ln[22:26].strip()
                if _is_element_token(resname):
                    ions.add((resname, f"{chain}:{resi}"))
        return ions

    def ions_in_pdbqt(p: Union[str, Path]) -> Set[Tuple[str, str]]:
        ions = set()
        if not Path(p).exists():
            return ions
        with open(p, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if not ln.startswith(("ATOM  ", "HETATM")):
                    continue
                resname = ln[17:20].strip().upper()
                chain = ln[21]
                resi = ln[22:26].strip()
                if _is_element_token(resname):
                    ions.add((resname, f"{chain}:{resi}"))
        return ions

    pdb_ions = ions_in_pdb(clean_pdb)
    pdbqt_ions = ions_in_pdbqt(receptor_pdbqt)
    missing = pdb_ions - pdbqt_ions
    if missing:
        logging.warning("[iondiff.warn] retained_missing=%s", sorted(missing))


def log_possible_metal_mislabels(pdb_path: Union[str, Path]) -> None:
    """
    Heuristic: flag HET residues whose resname is not an element token, but
    whose element column shows an element token and the residue has very few atoms.
    """
    suspects = []
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        by_res = {}
        for ln in f:
            if ln.startswith(("ATOM  ", "HETATM")):
                key = (ln[21], ln[22:26], ln[26], ln[17:20].strip().upper())
                by_res.setdefault(key, []).append(ln)

    for (chain, resi, icode, resname), atms in by_res.items():
        if len(atms) > 4:
            continue
        if _is_element_token(resname):
            continue
        elem_tokens = {ln[76:78].strip().upper() for ln in atms if len(ln) >= 78}
        if any(_is_element_token(e) for e in elem_tokens):
            suspects.append((resname, f"{chain}:{resi}", sorted(elem_tokens)))

    if suspects:
        logging.warning(
            "Possible metal mislabels (check resname vs element cols): %s", suspects
        )


# =============================
# Structural cleaning utilities
# =============================


def build_missing_loops(
    input_pdb: Union[str, Path], output_dir: Union[str, Path]
) -> str:
    """Fill missing loops/residues using MODELLER; return output PDB path."""
    try:
        from modeller import environ, log
        from modeller.scripts import complete_pdb
    except Exception as e:
        logging.warning("MODELLER not available; skipping loop completion: %s", e)
        return str(input_pdb)

    output_pdb = os.path.join(str(output_dir), "modeller_filled.pdb")
    log.none()
    logging.info("Running MODELLER to complete missing parts of %s", input_pdb)

    try:
        env = environ()
        env.io.hetatm = True
        env.io.water = True
        env.libs.topology.read(file="$(LIB)/top_heav.lib")
        env.libs.parameters.read(file="$(LIB)/par.lib")
        # MODELLER loop building
        mdl = complete_pdb(env, str(input_pdb))
        mdl.write(file=output_pdb)
        _post_write_element_guard("MODELLER", output_pdb)

        if os.path.exists(output_pdb):
            logging.info("MODELLER filled PDB saved to %s", output_pdb)
            return output_pdb
        logging.warning("MODELLER did not produce expected output: %s", output_pdb)
        return str(input_pdb)
    except Exception as e:
        logging.warning("MODELLER failed on %s: %s", input_pdb, e)
        return str(input_pdb)


def filter_invalid_chains(
    pdb_path: Union[str, Path], output_path: Union[str, Path]
) -> None:
    """Remove entire chains lacking backbone atoms (CA/N/C/O)."""
    chains: Dict[str, List[str]] = defaultdict(list)
    valid_chains: Set[str] = set()

    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                chain_id = line[21]
                atom_name = line[12:16].strip()
                chains[chain_id].append(line)
                if atom_name in {"CA", "N", "C", "O"}:
                    valid_chains.add(chain_id)
            else:
                chains["HEADER"].append(line)

    with open(output_path, "w", encoding="utf-8") as f:
        for chain_id in chains:
            if chain_id == "HEADER" or chain_id in valid_chains:
                f.writelines(chains[chain_id])
            else:
                logging.warning("Skipping invalid chain '%s' (no CA atoms)", chain_id)


def _has_backbone_atoms(lines):
    from protein_prep.chain_prune import _has_backbone_atoms as _impl

    return _impl(lines)



def _group_by_chain(lines):
    from protein_prep.chain_prune import _group_by_chain as _impl

    return _impl(lines)



def _prune_chains_conservative(lines, keep_chains):
    from protein_prep.chain_prune import _prune_chains_conservative as _impl

    return _impl(lines, keep_chains)



def _chains_to_keep(lines, pocket_center=None, r=12.0):
    from protein_prep.chain_prune import _chains_to_keep as _impl

    return _impl(lines, pocket_center=pocket_center, r=r)



# --- Hydrogen cleanup (geometry + CONECT) ---


def _parse_xyz(line: str):
    try:
        return (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except Exception:
        return None


def remove_unbonded_atoms(pdb_path: Union[str, Path]) -> None:
    from protein_prep.hydrogen_cleanup import remove_unbonded_atoms as _impl

    _impl(pdb_path)


def remove_implausible_hydrogens_by_distance(pdb_path: Union[str, Path]) -> None:
    from protein_prep.hydrogen_cleanup import (
        remove_implausible_hydrogens_by_distance as _impl,
    )

    _impl(pdb_path)


def clean_hydrogens(
    pdb_path: Union[str, Path],
    use_conect_if_reliable: bool = True,
    conect_min_cov: float = 0.6,
) -> None:
    from protein_prep.hydrogen_cleanup import clean_hydrogens as _impl

    _impl(
        pdb_path,
        use_conect_if_reliable=use_conect_if_reliable,
        conect_min_cov=conect_min_cov,
    )


# =============================
# Cofactors / waters policy
# =============================


def _is_metal(resname: str) -> bool:
    """Treat as 'metal/ion' iff the name normalizes to a canonical element token."""
    canonical = _normalize_resname(resname)
    return bool(canonical) and (canonical in _ELEM_CANON)


def _cofactor_policy_keep(resname: str) -> bool:
    from protein_prep.strip_nsr import _cofactor_policy_keep as _cofactor_keep_impl

    return _cofactor_keep_impl(resname)



def strip_nonstandard_residues(
    input_pdb: Union[str, Path],
    output_pdb: Union[str, Path],
    *,
    variant: Optional[str] = None,
) -> Tuple[int, str]:
    from protein_prep.strip_nsr import (
        strip_nonstandard_residues as _strip_nonstandard_residues_impl,
    )

    return _strip_nonstandard_residues_impl(
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        variant=variant,
    )



def quick_element_histogram(pdb_path: Union[str, Path]) -> None:
    cnt = Counter()
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM", "HETATM")):
                el = ln[76:78].strip().upper()
                cnt[el or ""] += 1
    logging.info(
        "[Elem histogram %s] %s",
        os.path.basename(str(pdb_path)),
        dict(sorted(cnt.items())),
    )


def assert_no_metal_in_peptidic(pdb_path: Union[str, Path]) -> None:
    # Build a peptide-like name set from YAML
    pep_name_lines = RULES["element_sets"].get("peptide_like_names", [])
    peptidey = set()
    for line in pep_name_lines:
        for tok in str(line).split(";"):
            t = tok.strip().upper()
            if t:
                peptidey.add(t)

    # PTM whitelist: YAML override if present; else default PTR/SEP/TPO
    ptm_yaml = set(
        _flatten_semicolons(RULES.get("element_sets", {}).get("ptm_resnames", []))
    )
    ptm_resnames = ptm_yaml or {"PTR", "SEP", "TPO"}

    # Element tokens considered "ionic" from YAML context (retain + element list)
    rules_snapshot = get_atom_rules()
    ionic_tokens = {
        str(tok).strip().upper()
        for tok in getattr(rules_snapshot, "elem_tokens_canonical", set())
        if str(tok).strip()
    }
    if not ionic_tokens:
        ionic_tokens = set(_ELEM_CANON)

    res_atoms = defaultdict(list)
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM  ", "HETATM")):
                key = (ln[21], ln[22:26], ln[26], ln[17:20].strip().upper())
                res_atoms[key].append(ln)

    offenders = []
    for key, lines in res_atoms.items():
        chain, resi, icode, resname = key
        # Skip known PTMs entirely
        if resname in ptm_resnames:
            continue
        names = [ln[12:16].strip().upper() for ln in lines]
        if not names:
            continue
        hits = sum(
            (n in peptidey) or (n[:2] in {"OE", "NE", "OD", "ND", "SD"}) for n in names
        )
        pep_like = hits >= max(4, 0.6 * len(names))
        if not pep_like:
            continue
        # if any atom's element is an ionic token (from YAML), flag
        if any((ln[76:78].strip().upper() in ionic_tokens) for ln in lines):
            offenders.append(key)

    # Compact info line to make the warning greppable/noisy only when meaningful
    logging.info(
        "[peptide-ion check] offenders=%s peptidey_n=%d ionic_tokens_n=%d",
        offenders,
        len(peptidey),
        len(ionic_tokens),
    )

    if offenders:
        logging.warning(
            "Peptide-like residues contain ionic elements (check labeling): %s",
            offenders,
        )


# DELETE CHAINS HELPERS
def detect_pocket_center_from_ligands(
    filtered_pdb: Union[str, Path], ligands_dir: Union[str, Path]
) -> Optional[Tuple[float, float, float]]:
    from protein_prep.chain_prune import (
        detect_pocket_center_from_ligands as _detect_pocket_center_impl,
    )

    return _detect_pocket_center_impl(filtered_pdb, ligands_dir)



def score_chain_contacts(
    pdb_path: Union[str, Path], keep_chains: set[str]
) -> Dict[str, int]:
    from protein_prep.chain_prune import score_chain_contacts as _score_contacts_impl

    return _score_contacts_impl(pdb_path, keep_chains)



def select_chains_to_keep(
    filtered_pdb: Union[str, Path], ligands_dir: Union[str, Path], cfg=None
) -> set[str]:
    from protein_prep.chain_prune import select_chains_to_keep as _select_chains_impl

    return _select_chains_impl(filtered_pdb, ligands_dir, cfg=cfg)



def prune_to_chains(
    input_pdb: Union[str, Path], kept_chains: set[str], output_pdb: Union[str, Path]
) -> None:
    from protein_prep.chain_prune import prune_to_chains as _prune_to_chains_impl

    _prune_to_chains_impl(input_pdb, kept_chains, output_pdb)



# ============================
# End-to-end Cleaning Pipeline
# =============================
def clean_pdb(
    pdb_file: Union[str, Path],
    output_root: Union[str, Path],
    logger: Optional[logging.Logger] = None,
) -> Optional[str]:
    from protein_prep.clean_pipeline import clean_pdb as _clean_pipeline_pdb

    return _clean_pipeline_pdb(
        pdb_file=pdb_file,
        output_root=output_root,
        logger=logger,
    )


def _phenix_detect() -> tuple[str, list[str] | None, dict]:
    from protein_prep.phenix_tools import _phenix_detect as _impl

    return _impl()


def run_phenix_pdbtools(
    input_pdb: Union[str, Path],
    output_pdb: Union[str, Path],
    remove_waters: bool = True,
) -> bool:
    from protein_prep.phenix_tools import run_phenix_pdbtools as _impl

    return _impl(
        input_pdb=input_pdb,
        output_pdb=output_pdb,
        remove_waters=remove_waters,
    )


def run_windows_phenix_clean_script(
    loop_fixed_pdb: Union[str, Path], nolig_dir: Union[str, Path]
) -> int:
    from protein_prep.phenix_tools import run_windows_phenix_clean_script as _impl

    return _impl(loop_fixed_pdb, nolig_dir)


def run_molprobity_validate(
    pdb_path: Union[str, Path], work_dir: Optional[Union[str, Path]] = None
) -> int:
    from protein_prep.phenix_tools import run_molprobity_validate as _impl

    return _impl(pdb_path=pdb_path, work_dir=work_dir)


# =============================
# Module Entrypoint
# =============================


def main(
    pdb_filename: str,
    output_dir: Union[str, Path] = r"./processed_pdbs",
    center: Optional[Tuple[float, float, float]] = None,
    box_size: Optional[Tuple[float, float, float]] = None,
):
    """High-level wrapper: clean a PDB and prepare the receptor PDBQT."""
    try:
        # Resolve input path robustly:
        # 1) If the provided path already points to a file (absolute or relative), use it as-is.
        # 2) Otherwise, fall back to INPUT_DIR/pdb_filename.
        cand = Path(pdb_filename)
        if cand.is_file():
            pdb_path = str(cand.resolve())
        else:
            pdb_path = str(Path(_cfg("INPUT_DIR", ".")).joinpath(pdb_filename))

        logging.info(
            "[prep] argv pdb_filename=%s resolved=%s output_dir=%s cwd=%s",
            pdb_filename,
            pdb_path,
            str(output_dir),
            os.getcwd(),
        )

        if not os.path.isfile(pdb_path):
            logging.error("ERROR: File does not exist: %s", pdb_path)
            raise FileNotFoundError(pdb_path)

        # Run cleaning → returns the final cleaned receptor PDB path
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        cleaned_pdb = clean_pdb(pdb_path, output_dir)
        if not cleaned_pdb:
            logging.error("ERROR: Cleaning failed for %s", pdb_filename)
            raise RuntimeError(f"cleaning_failed: {pdb_filename}")

        # Prepare receptor PDBQT next to the cleaned tree
        raw_stem = Path(pdb_filename).stem
        pdb_id = re.sub(
            r"(_nolig(_cleaned)?|_cleaned)$", "", raw_stem, flags=re.I
        ).upper()
        logging.info("[prep.id] main stem=%s -> base_id=%s", raw_stem, pdb_id)
        legacy_paths = canon_paths(pdb_id, output_dir)
        logging.info(
            "[prep.paths] protein_root=%s receptor=%s nolig=%s work=%s",
            legacy_paths["protein_root"],
            legacy_paths["receptor"],
            legacy_paths["nolig"],
            legacy_paths["work"],
        )

        paths = make_paths(config, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
        cleaned_pdb_out = str(paths.receptor_cleaned_pdb(None))
        receptor_pdbqt_out = str(paths.receptor_pdbqt(None, ph_token=None))

        cleaned_target = Path(cleaned_pdb_out)
        cleaned_target.parent.mkdir(parents=True, exist_ok=True)
        if Path(cleaned_pdb).resolve() != cleaned_target.resolve():
            shutil.copyfile(cleaned_pdb, str(cleaned_target))
        cleaned_pdb = str(cleaned_target)
        receptor_target = Path(receptor_pdbqt_out)
        receptor_target.parent.mkdir(parents=True, exist_ok=True)
        output_pdbqt = str(receptor_target)
        fix_element_columns_in_file(cleaned_pdb, cleaned_pdb, rewrite_atoms=True)

        if not run_prepare_receptor(cleaned_pdb, output_pdbqt, config):
            logging.error("ERROR: Failed to prepare receptor PDBQT for %s", pdb_id)
            try:
                rp = Path(output_pdbqt)
                exists = rp.exists()
                size = rp.stat().st_size if exists else 0
            except Exception:
                exists = False
                size = 0
            logging.info(
                "[receptor-summary]\n"
                "cleaned_pdb=%s\n"
                "receptor_pdbqt=%s exists=%s size=%d\n"
                "meeko_attempts=(see work/*.cmd.txt | *.stderr.txt)",
                cleaned_pdb,
                output_pdbqt,
                exists,
                size,
            )
            raise RuntimeError(f"receptor_pdbqt_failed: {pdb_id}")

        # Run base receptor PDBQT build
        if not run_prepare_receptor(cleaned_pdb, output_pdbqt, config):
            logging.error("ERROR: Failed to prepare receptor PDBQT for %s", pdb_id)
            try:
                rp = Path(output_pdbqt)
                exists = rp.exists()
                size = rp.stat().st_size if exists else 0
            except Exception:
                exists = False
                size = 0
            logging.info(
                "[receptor-summary]\n"
                "cleaned_pdb=%s\n"
                "receptor_pdbqt=%s exists=%s size=%d\n"
                "meeko_attempts=(see work/*.cmd.txt | *.stderr.txt)",
                cleaned_pdb,
                output_pdbqt,
                exists,
                size,
            )
            raise RuntimeError(f"receptor_pdbqt_failed: {pdb_id}")

        # --- HOLO-only stash & restore (moved to process_one_protein) ---
        # NOTE: _holo_restore_from_input_if_needed is now invoked from main.process_one_protein
        # after pocket detection, once the docking box center/size are known.
        # This block is intentionally left as a no-op to avoid double-restore.

        # Success path summary
        try:
            rp = Path(output_pdbqt)
            exists = rp.exists()
            size = rp.stat().st_size if exists else 0
        except Exception:
            exists = False
            size = 0
        logging.info(
            "[receptor-summary]\n"
            "cleaned_pdb=%s\n"
            "receptor_pdbqt=%s exists=%s size=%d\n"
            "meeko_attempts=(see work/*.cmd.txt | *.stderr.txt)",
            cleaned_pdb,
            output_pdbqt,
            exists,
            size,
        )

        try:
            _log_pdb_pdbqt_counts_diff(cleaned_pdb, output_pdbqt, tool="final")
        except Exception as exc:
            logging.warning("[iondiff.pdb_pdbqt] action=skip reason=%s", exc)

        logging.info(
            "[prep.return] cleaned=%s receptor_pdbqt=%s", cleaned_pdb, output_pdbqt
        )
        logging.info("Prepared receptor PDBQT: %s", output_pdbqt)
        return cleaned_pdb, output_pdbqt

    except Exception as e:
        logging.exception("[FATAL] automate_protein_prep.main() failed: %s", e)
        raise
