from __future__ import annotations

import hashlib
import json
import logging
import os
import subprocess
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Optional, Tuple, Mapping, TypedDict, cast


from docking.docking_vina import emit_vina_config
from path_router.path_router import Paths, docked_dir, receptor_file
from docking.active_site_detection import detect_active_site
from docking.run_vina import run_docking_task
from docking.global_scheduler import (
    acquire_global_cores,
    get_global_docking_sem,
)

from docking.ligand_metrics import _is_readable_ref, compute_rmsd
from docking.pose_validation import compute_redock_rmsd
from protein_prep.pdb_records import centroid as _coord_centroid
from protein_prep.pdb_records import line_xyz


class _ControlRedockResult(TypedDict):
    order: int
    base: str
    center: tuple[float, float, float]
    ligand_name: str
    best_pdb: str | None
    best_e: float | None
    score: float | None
    out_path: str


def _canonical_ctrl_base_from_stem(stem: str) -> str:
    """
    Normalize control ligand stems so variants like `LIG_A301.sanitized.protoB` map to `LIG_A301`.
    Drops any _stage suffix and trims trailing .sanitized chains.
    """
    stem = stem.split("_stage")[0]
    return stem.split(".sanitized")[0]


def _cfg_bool_value(cfg: Mapping[str, Any], key: str, default: bool) -> bool:
    raw = cfg.get(key, default)
    if isinstance(raw, bool):
        return raw
    return str(raw).strip().lower() in {"1", "true", "yes", "y", "on"}


def _control_redock_stage_info(cfg: Mapping[str, Any]) -> dict[str, Any]:
    info: dict[str, Any] = {
        "exhaustiveness": int(cfg.get("CTRL_REDOCK_EXHAUSTIVENESS", 64)),
        "num_modes": int(cfg.get("CTRL_REDOCK_NMODES", 9)),
        "energy_range": float(cfg.get("CTRL_REDOCK_ENERGY_RANGE", 6)),
        "verbosity": int(cfg.get("VINA_VERBOSITY", 0)),
        "seed": int(cfg.get("CTRL_REDOCK_SEED", 1337)),
    }
    if cfg.get("FAST_MODE") and _cfg_bool_value(
        cfg,
        "CTRL_REDOCK_RESPECT_FAST_MODE",
        False,
    ):
        info["exhaustiveness"] = 1
        info["num_modes"] = 1
    return info


# [ions] audit classification tokens
_ION_AUDIT_METALS = {
    "ZN",
    "MG",
    "MN",
    "FE",
    "CU",
    "CO",
    "NI",
    "CA",
    "HG",
}
_ION_AUDIT_SALTS = {
    "NA",
    "K",
    "CL",
    "BR",
    "I",
    "LI",
    "RB",
    "CS",
}


# --------- control ligand lookup (prefer SDF > MOL2 > PDB; search ligands_raw + reference) ---------
def build_control_lookup(paths: Paths) -> dict:
    """
    Map base extracted-ligand stem -> crystal file path.
    Prefer a readable .sdf > .mol2 > .pdb, searching ligands_raw and optional reference folder.
    """
    prefs = [".sdf", ".mol2", ".pdb"]
    by_base: dict[str, dict[str, Path]] = {}

    search_dirs = [
        paths.ligand_output_dir,
        paths.ligand_output_dir.parent / "reference",
    ]
    for root in search_dirs:
        if not root.exists():
            continue
        for p in root.rglob("*"):
            ext = p.suffix.lower()
            if ext not in prefs:
                continue
            base = _canonical_ctrl_base_from_stem(p.stem)
            by_base.setdefault(base, {})
            by_base[base][ext] = p

    chosen: dict[str, Path] = {}
    for base, candidates in by_base.items():
        picked = None
        for ext in prefs:
            candidate_path = candidates.get(ext)
            if not candidate_path:
                continue
            # Prefer RDKit-readable refs; for PDB keep a non-empty file as last resort
            # so control redock can still run when SDF/MOL2 parsing is unavailable.
            if _is_readable_ref(candidate_path):
                picked = candidate_path
                break
            if ext == ".pdb":
                try:
                    if candidate_path.exists() and candidate_path.stat().st_size > 0:
                        picked = candidate_path
                        break
                except Exception:
                    continue
        if picked:
            chosen[base] = picked

    return chosen


def receptor_sanity_check(receptor_pdbqt: str, min_atoms: int = 10) -> bool:
    try:
        atoms = 0
        any_nonzero = False
        with open(receptor_pdbqt, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                xyz = line_xyz(ln)
                if xyz is not None:
                    atoms += 1
                    if (abs(xyz[0]) + abs(xyz[1]) + abs(xyz[2])) > 0.0:
                        any_nonzero = True
        return (atoms >= min_atoms) and any_nonzero
    except Exception:
        return False


# ======================
# Phase 1 5: Prep steps
# ======================
def extract_ligands_to_nolig(paths: Paths, logger: logging.Logger) -> Tuple[int, set]:
    """
    Run crystallographic ligand extraction into processed_pdbs/<PDB>/ligands_raw/,
    create/update nolig/<PDB>_nolig.pdb, and return (n_controls, control_stems).
    Compatible with multiple legacy signatures of activesite.extract_and_remove_ligands().
    """
    from docking.pockets.ligand_pocket_runtime import extract_and_remove_ligands  # authoritative extractor

    # Ensure intermediate dirs exist and clear the malformed log for a fresh run
    paths.ligand_output_dir.mkdir(parents=True, exist_ok=True)
    paths.ligands_mol2_dir.mkdir(parents=True, exist_ok=True)
    malformed_log = paths.ligands_mol2_dir / "malformed_ligands.txt"
    if malformed_log.exists():
        malformed_log.unlink()

    src_pdb = paths.input_pdb_path
    nolig_dst = paths.nolig_pdb_path
    ligands_dir = paths.ligand_output_dir

    logger.debug(
        "[extract.debug] in=%s nolig=%s ldir=%s", src_pdb, nolig_dst, ligands_dir
    )

    ligands_dict, _ = extract_and_remove_ligands(
        str(src_pdb), str(nolig_dst), str(ligands_dir)
    )
    logger.info(f"Extracted {len(ligands_dict)} ligands -> {ligands_dir}")
    try:
        counts = {".pdb": 0, ".mol2": 0, ".sdf": 0}
        samples: list[str] = []
        for ext in (".pdb", ".mol2", ".sdf"):
            for p in paths.ligand_output_dir.glob(f"*{ext}"):
                counts[ext] += 1
                if len(samples) < 6:
                    samples.append(p.name)
        logger.info("[extract.audit] counts=%s samples=%s", counts, samples)
    except Exception as e:
        logger.debug("[extract.audit] listing failed: %s", e)

    # Build control stems from what was actually written
    control_stems: set[str] = set()
    for ext in (".mol2", ".pdb", ".sdf"):
        for p in paths.ligand_output_dir.rglob(f"*{ext}"):
            control_stems.add(Path(p).stem)

    logger.info(f"Extracted {len(control_stems)} ligands -> {paths.ligand_output_dir}")
    return len(control_stems), control_stems


def _ph_control_centroid(paths: Paths) -> Optional[Tuple[float, float, float]]:
    try:
        primary = paths.ligand_output_dir
    except Exception:
        primary = None
    roots = []
    if primary is not None:
        roots.append(primary)
        legacy = primary.parent.parent / f"{paths.pdb_id}_NOLIG" / "ligands_raw"
        roots.append(legacy)
    else:
        roots.append(paths.root_pdb_dir / "ligands_raw")
        roots.append(
            paths.root_pdb_dir.parent / f"{paths.pdb_id}_NOLIG" / "ligands_raw"
        )
    centroids = []
    seen = set()
    for root in roots:
        if not root:
            continue
        root = Path(root)
        key = str(root)
        if key in seen or not root.exists():
            continue
        seen.add(key)
        for pdb_path in sorted(root.glob("*.pdb")):
            coords: list[tuple[float, float, float]] = []
            try:
                with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
                    for line in fh:
                        xyz = line_xyz(line)
                        if xyz is not None:
                            coords.append(xyz)
            except Exception:
                continue
            center = _coord_centroid(coords)
            if center is not None:
                centroids.append(center)
    if centroids:
        n = float(len(centroids))
        return (
            sum(x for x, _, _ in centroids) / n,
            sum(y for _, y, _ in centroids) / n,
            sum(z for _, _, z in centroids) / n,
        )
    return None


# pH helpers  ----------------------
def _resolve_ph_scope(
    scope_cfg: str,
    radius_nominal: float,
    paths: Paths,
    cleaned_pdb: str,
    log: logging.Logger,
) -> Tuple[str, Tuple[float, float, float], float]:
    scope = (scope_cfg or "").strip().lower()
    try:
        radius = float(radius_nominal)
    except Exception:
        radius = 10.0
    if radius <= 0:
        radius = 10.0
    if scope == "pocket":
        center = _ph_control_centroid(paths)
        if center is None:
            log.warning(
                "[ph_ensemble.scope] pocket requested but no control center found; fallback=global"
            )
            return "global", (0.0, 0.0, 0.0), 1_000_000.0
        cx, cy, cz = (float(center[0]), float(center[1]), float(center[2]))
        return "pocket", (cx, cy, cz), radius
    return "global", (0.0, 0.0, 0.0), 1_000_000.0


def _ph_values_from_context(pdb_path: str) -> list[float]:
    """
    Ask context_ph for the full pH list (ensemble if present, else [target]),
    then round to 0.1 and clamp to [3.0, 10.5].
    """
    vals = []
    try:
        from path_router.context_ph import select_ph_values_for_protonation

        raw = select_ph_values_for_protonation(pdb_path)  # returns ensemble or [target]
        logging.info(f"[ph.ctx.list] taken_from_context={raw}")

        for x in raw or []:
            # round & clamp
            v = max(3.0, min(10.5, round(float(x), 1)))
            vals.append(v)
        # dedupe + sort for stability
        vals = sorted({v for v in vals})
    except Exception as e:
        logging.warning(f"[ph.context] failed to resolve; falling back to [7.0]: {e}")
        vals = [7.0]
    logging.info(f"[ph.list] n={len(vals)} values={vals}")
    return vals


def _ph_ligand_mode(cfg: Mapping[str, Any]) -> str:
    """
    Interpret PH_LIGAND_MODE from config.

    Recognized values (case-insensitive):
    - off / none / false / 0 / "" -> "off"
    - anything else -> "context_window"

    This keeps the behavior opt-in while allowing future modes later.
    """
    raw = str(cfg.get("PH_LIGAND_MODE", "off")).strip().lower()
    if raw in ("", "off", "none", "false", "0"):
        return "off"
    return "context_window"


# ----------------------


# ----------------------


def _ensure_ctrl_vina_manifest(
    cfg, paths, stage_name, variant_token, ph_label, legacy_mode, logger
):
    run_id = cfg["RUN_ID"]
    stage_dir = paths.configs_stage_dir(run_id, variant_token, stage_name, ph_label)
    manifest_path = stage_dir / "vina.json"
    stage_dir.mkdir(parents=True, exist_ok=True)
    precreated = 0
    if not manifest_path.exists():
        payload = {
            "run_id": run_id,
            "pdb_id": paths.pdb_id,
            "stage": stage_name,
            "variant": variant_token,
            "ph": ph_label,
            "legacy": bool(legacy_mode),
            "entries": [],
        }
        tmp = manifest_path.with_suffix(".part")
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, manifest_path)
        precreated = 1
    logger.info("[ctrl.vina.json] precreated=%d path=%s", precreated, manifest_path)
    return manifest_path


def _finalize_ctrl_vina_manifest(
    cfg,
    paths,
    stage_name,
    variant_token,
    ph_label,
    legacy_mode,
    manifest_path,
    jobs,
):
    if manifest_path is None:
        return
    run_id = cfg["RUN_ID"]
    stage_dir = manifest_path.parent
    receptor_path = receptor_file(
        paths.pdb_id,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    out_root = docked_dir(
        paths.pdb_id,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    out_root.mkdir(parents=True, exist_ok=True)
    out_dir = out_root / stage_name
    out_dir.mkdir(parents=True, exist_ok=True)

    entries = []
    for job in jobs:
        lig_base = Path(job["ligand_path"]).stem
        entry = {
            "ligand": lig_base,
            "config": str(stage_dir / f"{lig_base}_{stage_name}.txt"),
            "out": str(out_dir / f"{lig_base}_{stage_name}.pdbqt"),
            "receptor": str(receptor_path),
        }
        entries.append(entry)

    entries.sort(key=lambda e: e["ligand"])

    payload = {
        "run_id": run_id,
        "pdb_id": paths.pdb_id,
        "stage": stage_name,
        "variant": variant_token,
        "ph": ph_label,
        "legacy": bool(legacy_mode),
        "entries": entries,
    }

    tmp = manifest_path.with_suffix(".part")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, manifest_path)


def _ctrl_quick_file_sig(pth: str) -> str:
    try:
        p = Path(pth)
        if not p.exists():
            return "exists=False size=-1 sha=0000000000"
        sz = p.stat().st_size
        h = hashlib.sha1()
        with open(p, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if ln.startswith(("ATOM", "HETATM")):
                    h.update(ln[12:54].encode("utf-8", "ignore"))
        return f"exists=True size={sz} sha={h.hexdigest()[:10]}"
    except Exception:
        return "sig=unavailable"


def _ctrl_best_model_to_pdb(pdbqt_file: str, obabel: str):
    from pathlib import Path as _Path
    import tempfile
    import shutil as _sh

    best_e = None
    best_chunk = None
    pdbqt_path = _Path(pdbqt_file)
    try:
        if (not pdbqt_path.exists()) or (pdbqt_path.stat().st_size <= 0):
            return None, None
    except Exception:
        return None, None
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
        chunk = []
        in_model = False
        for ln in fh:
            u = ln.strip().upper()
            if u.startswith("MODEL"):
                chunk = [ln]
                in_model = True
            elif u.startswith("ENDMDL"):
                chunk.append(ln)
                in_model = False
                for cl in chunk:
                    if "REMARK VINA RESULT" in cl.upper():
                        try:
                            e = float(cl.strip().split()[3])
                            if (best_e is None) or (e < best_e):
                                best_e = e
                                best_chunk = chunk[:]
                        except Exception:
                            pass
            else:
                if in_model:
                    chunk.append(ln)
    if best_chunk is None:
        try:
            with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as fh:
                for ln in fh:
                    if "REMARK VINA RESULT" in ln.upper():
                        best_e = float(ln.strip().split()[3])
                        break
            best_chunk = None
        except Exception:
            return None, None

    td = _Path(tempfile.mkdtemp(prefix="ctrl_redock_"))
    best_pdbqt = td / "best.pdbqt"
    if best_chunk:
        with open(best_pdbqt, "w", encoding="utf-8") as out:
            out.writelines(best_chunk)
    else:
        _sh.copy2(pdbqt_path, best_pdbqt)
    out_pdb = td / "best.pdb"
    try:
        subprocess.run(
            [obabel, "-ipdbqt", str(best_pdbqt), "-opdb", "-O", str(out_pdb)],
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception:
        return None, best_e
    return (out_pdb if out_pdb.exists() else None), best_e


def _ctrl_redock_job(payload: dict) -> dict:
    cfg = payload["cfg"]
    stage_info = dict(payload["stage_info"])
    conf_path, out_path = emit_vina_config(
        cfg,
        payload["pdb_id"],
        payload["receptor_pdbqt"],
        tuple(payload["center"]),
        tuple(payload["box_size"]),
        payload["ligand_path"],
        payload["stage_name"],
        stage_info,
        int(payload["threads_per_job"]),
        logger=None,
        variant=payload.get("variant"),
        ph_token=payload.get("ph"),
        legacy=bool(payload.get("legacy")),
    )
    job_logger = logging.getLogger("control-redock")
    job_logger.info(
        "[control-redock.job] base=%s order=%d lig=%s center=%s box=%s conf=%s out=%s",
        payload.get("base"),
        int(payload.get("order", -1)),
        payload.get("ligand_name"),
        payload.get("center"),
        payload.get("box_size"),
        conf_path,
        out_path,
    )
    try:
        _, score = run_docking_task(
            payload["vina_exe"],
            str(conf_path),
            payload["ligand_path"],
            str(out_path),
        )
    except Exception:
        score = None
    best_pdb, best_e = _ctrl_best_model_to_pdb(str(out_path), payload["obabel"])
    return {
        "order": payload["order"],
        "base": payload["base"],
        "center": tuple(payload["center"]),
        "ligand_name": payload["ligand_name"],
        "best_pdb": str(best_pdb) if best_pdb else None,
        "best_e": best_e,
        "score": score,
        "out_path": str(out_path),
    }


def select_center_via_control_redock(
    cfg,
    paths,
    receptor_pdbqt,
    logger,
    *,
    variant: Optional[str] = None,
    ph_token: Optional[str] = None,
    legacy: bool = False,
):
    """
    Returns (center_tuple, (24.0,24.0,24.0)) or (None, None).
    - If multiple controls and max pairwise centroid distance <= CONTROL_CENTER_CLOSE_MAX_A -> average (consensus).
    - Else (far apart), redock each prepped control and pick lowest RMSD vs its crystal.
    - If no controls or all redocks fail, returns (None, None) to signal P2Rank fallback.
    """
    import numpy as _np
    import math as _math
    from pathlib import Path as _Path
    from docking.run_vina import run_docking_task as _run_dock

    def _find_control_pdbs(d: _Path) -> list[_Path]:
        return sorted([p for p in d.glob("*.pdb") if p.is_file()])

    def _centroid_from_pdb(p: _Path):
        coords: list[tuple[float, float, float]] = []
        try:
            with open(p, "r", encoding="utf-8", errors="ignore") as f:
                for ln in f:
                    xyz = line_xyz(ln)
                    if xyz is not None:
                        coords.append(xyz)
        except Exception:
            return None
        return _coord_centroid(coords)

    # discover crystal controls in preferred locations (current + legacy sibling)
    ctrl_pdbs = _find_control_pdbs(paths.ligand_output_dir)
    logger.info(
        f"[control-redock] controls_found={len(ctrl_pdbs)} dir={paths.ligand_output_dir}"
    )
    if not ctrl_pdbs:
        legacy_dir = Path(
            paths.ligand_output_dir.parent.parent
            / f"{paths.pdb_id}_NOLIG"
            / "ligands_raw"
        )
        if legacy_dir.exists():
            ctrl_pdbs = _find_control_pdbs(legacy_dir)

    if not ctrl_pdbs:
        logger.warning(
            "[control-redock] No extracted control PDBs present; skipping redock."
        )
        return None, None  # let caller go to P2Rank directly if no control
    # Map base control name -> extracted crystal PDB (for Kabsch RMSD)
    ctrl_pdb_map: dict[str, _Path] = {}
    for p in ctrl_pdbs:
        try:
            base = _canonical_ctrl_base_from_stem(p.stem)
            # First hit wins; avoids ambiguity if multiple variants exist.
            if base not in ctrl_pdb_map:
                ctrl_pdb_map[base] = p
        except Exception:
            continue

    def _ref_for_base(base: str) -> Optional[_Path]:
        return control_lookup.get(base) or ctrl_pdb_map.get(base)

    policy = str(cfg.get("CONTROL_CENTER_POLICY", "best_redock")).lower().strip()
    variant_env = (
        variant
        if variant is not None
        else (os.environ.get("APO_HOLO_VARIANT", "") or "")
    )
    variant_token = str(variant_env).strip().upper() or None
    ph_label = str(ph_token).strip() if ph_token is not None else None
    legacy_mode = bool(legacy)
    thr = float(cfg.get("CONTROL_CENTER_CLOSE_MAX_A", 8.0))

    # compute centroids + pairwise spread
    centroids = {}
    for p in ctrl_pdbs:
        c = _centroid_from_pdb(p)
        if c:
            base = _canonical_ctrl_base_from_stem(p.stem)
            centroids[base] = c

    bases = list(centroids.keys())
    coords = [centroids[b] for b in bases]

    def _dist(a, b):
        return float(
            ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) ** 0.5
        )

    max_delta = 0.0
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            d = _dist(coords[i], coords[j])
            if d > max_delta:
                max_delta = d

    if not coords:
        raise RuntimeError(
            "[control-centers] No control centroids available; cannot select center."
        )
    logger.info(
        f"[control-centers] n={len(coords)} max?={max_delta:.2f}A policy={policy} thr={thr:.2f}A"
    )
    for b, c in zip(bases, coords):
        logger.debug(f"[control-centers] {b}: ({c[0]:.3f},{c[1]:.3f},{c[2]:.3f})")

    # single-control or simple policies
    if len(coords) == 1 and policy != "best_redock":
        center = coords[0]
        logger.info(
            f"[Control-center] chosen={bases[0]} center=({center[0]:.3f},{center[1]:.3f},{center[2]:.3f}) box=(24,24,24)"
        )
        return center, (24.0, 24.0, 24.0)

    if policy == "first":
        center = coords[0]
        logger.info(
            f"[Control-center] chosen={bases[0]} center=({center[0]:.3f},{center[1]:.3f},{center[2]:.3f}) box=(24,24,24)"
        )
        return center, (24.0, 24.0, 24.0)
    # if best_redock, skip consensus short-circuit:
    if policy == "best_redock":
        pass  # fall through to redock block below

    elif max_delta <= thr:
        # consensus average when controls are close
        c = (
            float(_np.mean([x for x, _, _ in coords])),
            float(_np.mean([y for _, y, _ in coords])),
            float(_np.mean([z for _, _, z in coords])),
        )
        logger.info(
            f"[Control-center] chosen=consensus center=({c[0]:.3f},{c[1]:.3f},{c[2]:.3f}) box=(24,24,24)"
        )
        return c, (24.0, 24.0, 24.0)

    elif policy == "average_when_close":
        # far apart ? fall back to first per spec
        center = coords[0]
        logger.info(
            f"[Control-center] chosen={bases[0]} center=({center[0]:.3f},{center[1]:.3f},{center[2]:.3f}) box=(24,24,24)"
        )
        return center, (24.0, 24.0, 24.0)

    # best_redock path (controls far apart)
    control_lookup = build_control_lookup(paths)  # base -> reference path
    # collect prepped control pdbqts (per-protein dir and global output dir)
    prepped_dirs = [
        paths.prepped_ligands_dir,
        _Path(str(cfg.get("PREPPED_LIGANDS_DIR", ""))) / paths.pdb_id,
    ]
    cand_pdbqts = []
    logger.info(
        f"[control-redock] search_prepped_dirs={[str(d) for d in prepped_dirs]}"
    )

    logger.info(f"[control-redock.debug] centroids_keys={list(centroids.keys())}")
    logger.info(
        f"[control-redock.debug] control_lookup_keys={list(control_lookup.keys())}"
    )

    seen = set()
    for root in prepped_dirs:
        if not root:
            continue
        root_path = _Path(root)
        if not root_path.exists():
            continue
        pdbqt_files = sorted(root_path.glob("*.pdbqt"))
        logger.info(
            "[control-redock] prepped_pdbqt_count dir=%s count=%d",
            root_path,
            len(pdbqt_files),
        )
        # Find all PDBQTs, then group by base to pick the best variant
        candidates_by_base: dict[str, list[_Path]] = {}
        for p in pdbqt_files:
            base = _canonical_ctrl_base_from_stem(p.stem)
            if base in centroids and _ref_for_base(base) is not None:
                candidates_by_base.setdefault(base, []).append(p)

        for base, paths_list in candidates_by_base.items():
            if base in seen:
                continue
            # Prefer the one with exactly one ".sanitized" in the stem
            best_p = paths_list[0]
            for p in paths_list:
                if p.stem.count(".sanitized") == 1:
                    best_p = p
                    break
            cand_pdbqts.append(best_p)
            seen.add(base)
    max_controls = max(0, int(cfg.get("CTRL_REDOCK_MAX_CONTROLS", 12) or 0))
    if max_controls and len(cand_pdbqts) > max_controls:
        logger.info(
            "[control-redock] limiting candidates old=%d new=%d policy=deterministic",
            len(cand_pdbqts),
            max_controls,
        )
        cand_pdbqts = cand_pdbqts[:max_controls]
    logger.info(f"[control-redock] candidates={len(cand_pdbqts)}")

    if not cand_pdbqts:
        logger.warning(
            "[control-redock] No prepped control PDBQTs found; redock impossible (will fall back)."
        )
        return None, None

    def _build_ctrl_stage_info() -> dict[str, Any]:
        return _control_redock_stage_info(cfg)

    if cfg.get("FAST_MODE") and not _cfg_bool_value(
        cfg,
        "CTRL_REDOCK_RESPECT_FAST_MODE",
        False,
    ):
        stage_preview = _build_ctrl_stage_info()
        logger.info(
            "[control-redock] fast_mode_ignored_for_controls=true exhaustiveness=%s num_modes=%s",
            stage_preview.get("exhaustiveness"),
            stage_preview.get("num_modes"),
        )

    threads_per_vina = int(
        cfg.get("THREADS_PER_VINA_CTRL", 8)
    )  # control redock uses its own threads default=8
    vina_exe = str(cfg.get("VINA_EXE") or cfg.get("VINA_PATH") or "vina")
    obabel = str(cfg.get("OPENBABEL_PATH") or "obabel")

    best: tuple[float, float | None, str, tuple[float, float, float]] | None = None

    cpu_total = int(cfg.get("CPU", os.cpu_count() or 1) or 1)
    cpu_total = max(1, cpu_total)
    workers = max(1, min(cpu_total, len(cand_pdbqts)))
    job_threads = max(1, min(threads_per_vina, max(1, cpu_total // workers)))
    total_threads = workers * job_threads
    logger.info(
        "[ctrl.parallel] CPU=%d workers=%d threads_per_job=%d total_threads=%d",
        cpu_total,
        workers,
        job_threads,
        total_threads,
    )
    if job_threads < threads_per_vina:
        logger.info(
            "[ctrl.parallel] adjust threads_per_vina old=%d new=%d",
            threads_per_vina,
            job_threads,
        )

    stage_name = "ctrl_redock"
    parallel_enabled = workers > 1 and len(cand_pdbqts) > 1
    if get_global_docking_sem(cfg) is not None and parallel_enabled:
        # Cross-process pools cannot coordinate with the in-process semaphore token pool.
        parallel_enabled = False
        logger.info(
            "[ctrl.parallel] global_scheduler_active=true forcing serial_ctrl_redock=true"
        )

    if not parallel_enabled:
        for order, lig_pdbqt in enumerate(cand_pdbqts):
            stem = lig_pdbqt.stem.split("_stage")[0]
            base = _canonical_ctrl_base_from_stem(stem)
            center = centroids.get(base)
            if not center:
                ref = _ref_for_base(base)
                if ref and ref.suffix.lower() == ".pdb":
                    center = _centroid_from_pdb(ref)
            if not center:
                continue

            stage_info = _build_ctrl_stage_info()
            logger.info(
                "[emit.debug] pdb=%s stage=%s variant=%s ph=%s",
                paths.pdb_id,
                stage_name,
                variant_token or "None",
                ph_label or "None",
            )
            conf_path, out_path = emit_vina_config(
                cfg,
                paths.pdb_id,
                receptor_pdbqt,
                center,
                (24.0, 24.0, 24.0),
                str(lig_pdbqt),
                stage_name,
                stage_info,
                job_threads,
                logger=None,
                variant=variant_token,
                ph_token=ph_label,
                legacy=legacy_mode,
            )
            logger.info(
                "[control-redock.job] base=%s order=%d lig=%s center=%s box=%s conf=%s out=%s",
                base,
                order,
                lig_pdbqt.name,
                tuple(center),
                (24.0, 24.0, 24.0),
                conf_path,
                out_path,
            )
            try:
                with acquire_global_cores(
                    cfg,
                    cores=max(1, int(job_threads)),
                    min_cores=1,
                    priority=3,
                    task_id=f"{paths.pdb_id}:ctrl:{base}:{order}",
                    task_type="dock_stage_ligand",
                ) as granted:
                    _, score = _run_dock(
                        vina_exe,
                        conf_path,
                        lig_pdbqt.name,
                        out_path,
                        cpu_override=max(1, int(granted)),
                    )
            except Exception:
                score = None

            best_pdb, best_e = _ctrl_best_model_to_pdb(str(out_path), obabel)
            ref_path = _ref_for_base(base)
            rmsd = float("inf")

            if best_pdb and ref_path:
                logger.info(
                    f"[rmsd.debug] ref={ref_path} | {_ctrl_quick_file_sig(str(ref_path))}"
                )
                logger.info(
                    f"[rmsd.debug] dock={best_pdb} | {_ctrl_quick_file_sig(str(best_pdb))}"
                )
                same_file = Path(ref_path).resolve() == Path(best_pdb).resolve()
                if same_file:
                    logger.warning(
                        "[rmsd.debug] ref and dock paths resolve to the same file! RMSD=0.0 is expected."
                    )
                try:
                    rmsd = compute_rmsd(str(ref_path), str(best_pdb))
                except Exception as e:
                    rmsd = float("inf")
                    logger.exception(f"[rmsd.debug] compute_rmsd failed: {e}")

            e_print = (
                best_e
                if (best_e is not None)
                else (score if score is not None else float("nan"))
            )
            logger.info(
                f"[control-redock] lig={lig_pdbqt.name} rmsd={rmsd:.2f}A score={e_print if e_print is not None else float('nan')} kcal/mol"
            )

            if _math.isfinite(rmsd):
                if (
                    (best is None)
                    or (rmsd < best[0])
                    or (
                        rmsd == best[0]
                        and (e_print is not None)
                        and (best[1] is None or e_print < best[1])
                    )
                ):
                    best = (
                        rmsd,
                        e_print if e_print is not None else None,
                        base,
                        center,
                    )
    else:
        cfg_payload = dict(cfg)
        jobs = []
        for order, lig_pdbqt in enumerate(cand_pdbqts):
            stem = lig_pdbqt.stem.split("_stage")[0]
            base = _canonical_ctrl_base_from_stem(stem)
            center = centroids.get(base)
            if not center:
                ref = _ref_for_base(base)
                if ref and ref.suffix.lower() == ".pdb":
                    center = _centroid_from_pdb(ref)
            if not center:
                continue
            stage_info = _build_ctrl_stage_info()
            logger.info(
                "[emit.debug] pdb=%s stage=%s variant=%s ph=%s",
                paths.pdb_id,
                stage_name,
                variant_token or "None",
                ph_label or "None",
            )
            jobs.append(
                {
                    "order": order,
                    "cfg": cfg_payload,
                    "pdb_id": paths.pdb_id,
                    "receptor_pdbqt": str(receptor_pdbqt),
                    "center": tuple(center),
                    "box_size": (24.0, 24.0, 24.0),
                    "ligand_path": str(lig_pdbqt),
                    "ligand_name": lig_pdbqt.name,
                    "stage_name": stage_name,
                    "stage_info": stage_info,
                    "threads_per_job": job_threads,
                    "variant": variant_token,
                    "ph": ph_label,
                    "legacy": legacy_mode,
                    "vina_exe": vina_exe,
                    "obabel": obabel,
                    "base": base,
                }
            )

        if not jobs:
            return None, None

        results: list[_ControlRedockResult] = []
        if jobs:
            try:
                with ProcessPoolExecutor(max_workers=workers) as pool:
                    future_map = {
                        pool.submit(_ctrl_redock_job, job): job for job in jobs
                    }
                    for fut in as_completed(future_map):
                        res = fut.result()
                        results.append(cast(_ControlRedockResult, res))
            except Exception as e:
                logger.warning(
                    "[control-redock] ProcessPoolExecutor unavailable for pdb=%s; "
                    "falling back to serial redock jobs (err=%s)",
                    paths.pdb_id,
                    e,
                )
                results = []
                for job in jobs:
                    try:
                        results.append(cast(_ControlRedockResult, _ctrl_redock_job(job)))
                    except Exception as serial_exc:
                        logger.warning(
                            "[control-redock] serial job failed for pdb=%s lig=%s err=%s",
                            paths.pdb_id,
                            str(job.get("ligand_name", "?")),
                            serial_exc,
                        )
        else:
            logger.warning(
                "[control-redock] No jobs constructed for ctrl_redock; skipping redock."
            )
            return None, None

        if not results:
            logger.warning(
                "[control-redock.summary] pdb=%s no results from ctrl_redock jobs (results empty); will fall back.",
                paths.pdb_id,
            )
            return None, None

        results.sort(key=lambda r: int(r["order"]))
        for result_row in results:
            base = result_row["base"]
            center = tuple(float(v) for v in result_row["center"])

            # Legacy control lookup (SDF/MOL2/PDB); kept only for logging / debugging.
            ref_path = _ref_for_base(base)

            # crystal PDB taken from ligands_raw (or legacy NOLIG) via ctrl_pdb_map
            crystal_pdb = ctrl_pdb_map.get(base)

            # Docked control redock PDBQT (full multi-model file that Vina wrote)
            dock_pdbqt = result_row.get("out_path")

            # Also keep best_pdb around for debugging if needed
            best_pdb = result_row["best_pdb"]

            rmsd = float("inf")

            # Prefer Kabsch RMSD using crystal PDB vs docked PDBQT
            if crystal_pdb and dock_pdbqt:
                logger.info(
                    f"[rmsd.debug] crystal_pdb={crystal_pdb} | {_ctrl_quick_file_sig(str(crystal_pdb))}"
                )
                # best_pdb is a PDB-converted best pose; still useful to log for sanity
                if best_pdb:
                    logger.info(
                        f"[rmsd.debug] best_pdb={best_pdb} | {_ctrl_quick_file_sig(str(best_pdb))}"
                    )
                # log ref_path for comparison, even though we don't use it for RMSD now
                if ref_path:
                    logger.info(f"[rmsd.debug] ref_lookup={ref_path}")

                kabsch_rmsd = compute_redock_rmsd(str(crystal_pdb), str(dock_pdbqt))
                if kabsch_rmsd is None:
                    logger.warning(
                        "[rmsd.debug] compute_redock_rmsd returned None for base=%s; treating RMSD=inf",
                        base,
                    )
                else:
                    rmsd = float(kabsch_rmsd)
            else:
                logger.warning(
                    "[rmsd.debug] missing crystal_pdb or dock_pdbqt for base=%s; crystal=%s dock_pdbqt=%s",
                    base,
                    crystal_pdb,
                    dock_pdbqt,
                )

            best_e = result_row["best_e"]
            score = result_row["score"]
            e_print = (
                best_e
                if (best_e is not None)
                else (score if score is not None else float("nan"))
            )
            logger.info(
                f"[control-redock] lig={result_row['ligand_name']} rmsd={rmsd:.2f}A "
                f"score={e_print if e_print is not None else float('nan')} kcal/mol"
            )

            if _math.isfinite(rmsd):
                if (
                    (best is None)
                    or (rmsd < best[0])
                    or (
                        rmsd == best[0]
                        and (e_print is not None)
                        and (best[1] is None or e_print < best[1])
                    )
                ):
                    best = (
                        rmsd,
                        e_print if e_print is not None else None,
                        base,
                        center,
                    )

    if best is None:
        # Fail-soft: when control redock produced no usable poses, keep docking
        # by falling back to a deterministic control centroid.
        fallback_base = bases[0] if bases else None
        fallback_center = centroids.get(fallback_base) if fallback_base else None
        if fallback_center is None:
            return None, None
        logger.warning(
            "[control-redock] no finite redock RMSD; fallback centroid base=%s center=(%.3f,%.3f,%.3f)",
            fallback_base,
            float(fallback_center[0]),
            float(fallback_center[1]),
            float(fallback_center[2]),
        )
        return fallback_center, (24.0, 24.0, 24.0)

    chosen_center = best[3]

    logger.info(
        f"[Control-center] chosen={best[2]} center=({chosen_center[0]:.3f},{chosen_center[1]:.3f},{chosen_center[2]:.3f})"
    )
    # BOX SIZE SPECIFIED HERE, NEED TO EDIT THIS TO CALCULATE BOX SIZE, LARGE BOX  SIZES DECREASE VINA  ACCURACY
    return chosen_center, (24.0, 24.0, 24.0)


def detect_pocket(
    cleaned_pdb: str, ligand_dir: Path, logger: logging.Logger
) -> Tuple[
    Optional[Tuple[float, float, float]],
    Optional[Tuple[float, float, float]],
    str,  # source ("control" | "ligand_top" | "p2rank" | "activesite" | "none")
]:
    """
    Prefer control ligands for docking center/box. If none, fall back to the active-site module.
    """

    def _his_counts_within(
        pdb_path: str, center_xyz: tuple[float, float, float], r: float = 6.0
    ) -> tuple[int, int, int]:
        HID = HIE = HIP = 0
        try:
            with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
                seen = set()
                cx, cy, cz = center_xyz
                for ln in fh:
                    if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
                        continue
                    res = ln[17:20].strip().upper()  # residue name
                    if res not in {"HID", "HIE", "HIP"}:
                        continue
                    xyz = line_xyz(ln)
                    if xyz is None:
                        continue
                    x, y, z = xyz
                    if (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= r * r:
                        # key by (chain, resseq, resname) so we count each residue once
                        key = (ln[21].strip(), ln[22:26].strip(), res)
                        if key in seen:
                            continue
                        seen.add(key)
                        if res == "HID":
                            HID += 1
                        elif res == "HIE":
                            HIE += 1
                        elif res == "HIP":
                            HIP += 1
        except Exception:
            pass
        return HID, HIE, HIP

    def find_control_pdbs(d: Path) -> list[Path]:
        return sorted([p for p in d.glob("*.pdb") if p.is_file()])

    # 1) Controls check in canonical ligands_raw
    ctrl_files = find_control_pdbs(ligand_dir)

    # --- AUDIT: summarize control centroids & policy ---
    def _centroid_of_pdb(p: Path) -> tuple[float, float, float] | None:
        coords: list[tuple[float, float, float]] = []
        try:
            with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                for ln in fh:
                    xyz = line_xyz(ln)
                    if xyz is not None:
                        coords.append(xyz)
        except Exception:
            return None
        return _coord_centroid(coords)

    centers = [c for c in (_centroid_of_pdb(p) for p in ctrl_files) if c]
    dmax = 0.0
    if len(centers) >= 2:
        for i in range(len(centers)):
            for j in range(i + 1, len(centers)):
                dx = centers[i][0] - centers[j][0]
                dy = centers[i][1] - centers[j][1]
                dz = centers[i][2] - centers[j][2]
                d = float((dx * dx + dy * dy + dz * dz) ** 0.5)
                if d > dmax:
                    dmax = d
    policy = "first" if ctrl_files else "activesite"
    logger.info(
        "[control-centers] n=%d max?=%.2f A policy=%s", len(ctrl_files), dmax, policy
    )

    # Back-compat (read-only): if none found, check legacy sibling <PDB>_NOLIG/ligands_raw
    if not ctrl_files:
        # ligand_dir = .../<PDB>/ligands_raw
        pdb_root = ligand_dir.parent  # .../<PDB>
        legacy = pdb_root.parent / f"{pdb_root.name}_NOLIG" / "ligands_raw"
        if legacy.exists():
            ctrl_files = find_control_pdbs(legacy)
            if ctrl_files:
                logger.info(
                    f"[Control-center] Found controls in legacy sibling: {legacy}"
                )

    if ctrl_files:
        p = ctrl_files[0]
        coords: list[tuple[float, float, float]] = []
        with open(p, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                xyz = line_xyz(ln)
                if xyz is not None:
                    coords.append(xyz)
        ctrl_center = _coord_centroid(coords)
        if ctrl_center is not None:
            box_size = (24.0, 24.0, 24.0)
            hid, hie, hip = _his_counts_within(cleaned_pdb, ctrl_center, r=6.0)
            logger.info(
                "[reduce] his={'HID':%d,'HIE':%d,'HIP':%d} flips_near_box=%d",
                hid,
                hie,
                hip,
                0,
            )
            return ctrl_center, box_size, "control"

    # 2) Fallback to active-site module (ligand-based or P2Rank)
    active_center: tuple[float, float, float] | None = None
    active_box_size: tuple[float, float, float] | None = None
    src = None

    try:
        active_center, active_box_size, src = detect_active_site(cleaned_pdb)
    except Exception as exc:
        logger.error("Active-site detection raised exception: %s", exc)
        active_center = None
        active_box_size = None
        src = None

    if active_center and active_box_size:
        final_box_size: tuple[float, float, float] = (
            min(28.0, float(active_box_size[0])),
            min(28.0, float(active_box_size[1])),
            min(28.0, float(active_box_size[2])),
        )
        source = str(src or "activesite")
        logger.info(
            "[active-site] Using center %s with box %s source=%s",
            active_center,
            final_box_size,
            source,
        )
        hid, hie, hip = _his_counts_within(cleaned_pdb, active_center, r=6.0)
        logger.info(
            "[reduce] his={'HID':%d,'HIE':%d,'HIP':%d} flips_near_box=%d",
            hid,
            hie,
            hip,
            0,
        )
        return active_center, final_box_size, source
    else:
        logger.error(
            "Active-site detection failed (no controls, active-site returned None)."
        )
        return None, None, "none"


def _summarize_ions_file(file_path: Path | str) -> dict[str, object]:
    path = Path(file_path)
    if not path.exists():
        return {
            "hist": "missing",
            "counts": {},
            "metals_present": False,
            "salts_present": False,
            "error": "missing",
        }
    counts: Counter[str] = Counter()
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for ln in fh:
                if not ln.startswith("HETATM"):
                    continue
                res = ln[17:20].strip().upper()
                elem = (ln[76:78].strip() or res).upper()
                token = elem if elem.isalpha() and 1 <= len(elem) <= 2 else res
                if token and token.isalpha() and len(token) <= 3:
                    counts[token] += 1
    except Exception as exc:  # pragma: no cover - diagnostics
        return {
            "hist": "error",
            "counts": {},
            "metals_present": False,
            "salts_present": False,
            "error": str(exc),
        }

    hist = (
        ",".join(f"{tok}:{counts[tok]}" for tok in sorted(counts)) if counts else "none"
    )
    metals_present = any(
        token in _ION_AUDIT_METALS and counts[token] > 0 for token in counts
    )
    salts_present = any(
        token in _ION_AUDIT_SALTS and counts[token] > 0 for token in counts
    )
    return {
        "hist": hist,
        "counts": dict(counts),
        "metals_present": metals_present,
        "salts_present": salts_present,
        "error": None,
    }
