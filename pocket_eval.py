from __future__ import annotations

import csv
import hashlib
import json
import logging
import math
import os
import shutil
import random
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
from sklearn.model_selection import StratifiedKFold

from analysis import dud_eval
from calibrator.chemdbl_calibrator import (
    DEFAULT_ACTIVITY_TYPES,
    DEFAULT_CHEMBL_MAX_PHASE,
    DEFAULT_LABEL_THRESHOLDS,
    DEFAULT_OUT_ROOT,
    DEFAULT_RETRIES,
    DEFAULT_TIMEOUT,
    run_calibrator_for_pdb,
)
from calibrator.uniprot_resolver import resolve_chain_uniprot_segments, select_primary_chain
from DeepCoy_duds.external_sources import fetch_chembl_labeled_smiles
from DeepCoy_duds.generate_dud_library import convert_smi_to_sdf
from docking.run_vina import run_docking_task
from ligand_pocket import compute_box_from_ligand_coords
from path_router.path_router import make_paths
from prep_ligands.prep_ligands_bulk import prep_ligands_with_mgltools
try:
    from rdkit import Chem  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    Chem = None

REPO_ROOT = Path(__file__).resolve().parent

_LOG = logging.getLogger(__name__)


def _sanitize_ligand_name_for_filename(name: str) -> str:
    name = name.strip()
    name = name.replace(" ", "_")
    name = re.sub(r"[^A-Za-z0-9_.+-]+", "_", name)
    return name[:80] or "ligand"


def pocket_id_from_entry(entry: Dict[str, Any]) -> str:
    pocket_id = entry.get("pocket_id")
    if pocket_id:
        return str(pocket_id)
    resname = (entry.get("ligand_resname") or entry.get("het_name") or "UNK").strip()
    chain = (entry.get("ligand_chain") or entry.get("protein_chain") or "X").strip()
    resnum = str(entry.get("ligand_resnum") or entry.get("resnum") or "0").strip()
    return f"{resname}_{chain}_{resnum}"


def _center_and_box_from_pocket(entry: Dict[str, Any]) -> Tuple[Optional[tuple], Optional[tuple]]:
    coords = entry.get("coords") or []
    if coords:
        parsed = []
        for xyz in coords:
            if not isinstance(xyz, (list, tuple)) or len(xyz) < 3:
                continue
            try:
                parsed.append((float(xyz[0]), float(xyz[1]), float(xyz[2])))
            except Exception:
                continue
        if parsed:
            return compute_box_from_ligand_coords(parsed)

    center = entry.get("center")
    bounds = entry.get("bounds") or {}
    min_b = bounds.get("min") if isinstance(bounds, dict) else None
    max_b = bounds.get("max") if isinstance(bounds, dict) else None
    if center and min_b and max_b:
        try:
            cx, cy, cz = (float(center[0]), float(center[1]), float(center[2]))
            sx = float(max_b[0]) - float(min_b[0]) + 5.0
            sy = float(max_b[1]) - float(min_b[1]) + 5.0
            sz = float(max_b[2]) - float(min_b[2]) + 5.0
            return (cx, cy, cz), (sx, sy, sz)
        except Exception:
            return None, None
    return None, None


def _load_calibrator_cache(path: Path, logger: logging.Logger) -> Optional[List[Dict[str, Any]]]:
    if not path.exists():
        return None
    try:
        if path.stat().st_size == 0:
            return None
    except Exception:
        return None

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        logger.warning("[pocket-eval] calibrator_cache unreadable=%s err=%s", path, exc)
        return None

    rows: Optional[List[Dict[str, Any]]] = None
    if isinstance(payload, list):
        rows = payload
    elif isinstance(payload, dict):
        rows = payload.get("rows") or payload.get("calibrators")

    if not isinstance(rows, list):
        return None

    cleaned: List[Dict[str, Any]] = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        ligand_id = row.get("ligand_id")
        smiles = row.get("smiles")
        label = row.get("label")
        if not ligand_id or not smiles or not label:
            continue
        cleaned.append(
            {"ligand_id": str(ligand_id), "smiles": str(smiles), "label": str(label)}
        )

    return cleaned if cleaned else None


def _write_calibrator_cache(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _load_smiles_lines(path: Path) -> List[str]:
    if not path.exists():
        return []
    lines: List[str] = []
    try:
        with path.open("r", encoding="utf-8") as handle:
            for raw in handle:
                stripped = raw.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                parts = stripped.split()
                if not parts:
                    continue
                lines.append(parts[0])
    except Exception:
        return []
    return lines


def _load_extracted_calibrator_labels(
    extracted_dir: Path,
) -> Optional[Dict[str, List[str]]]:
    if not extracted_dir.exists():
        return None
    labels: Dict[str, List[str]] = {}
    found_any = False
    for label in ("strong", "weak", "non"):
        smi_path = extracted_dir / f"{label}_binders.smi"
        if smi_path.exists():
            found_any = True
            labels[label] = _load_smiles_lines(smi_path)
    if not found_any:
        return None
    for label in ("strong", "weak", "non"):
        labels.setdefault(label, [])
    if not any(labels.values()):
        return None
    return labels


def _smiles_hash(smiles: str) -> str:
    return hashlib.sha1(smiles.encode("utf-8")).hexdigest()[:10]


def _build_rows_from_labels(labels: Dict[str, List[str]]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for label in ("strong", "weak", "non"):
        entries = labels.get(label) or []
        for smi in sorted(set(entries)):
            lid = f"chembl:{label}_{_smiles_hash(smi)}"
            rows.append({"ligand_id": lid, "smiles": smi, "label": label})
    return rows


def _filter_supported_calibrators(
    rows: List[Dict[str, Any]], logger: logging.Logger
) -> List[Dict[str, Any]]:
    if Chem is None:
        return rows

    allowed = {"C", "H", "N", "O", "S", "P", "F", "CL", "BR", "I"}
    filtered: List[Dict[str, Any]] = []
    removed = 0
    for row in rows:
        smi = row.get("smiles")
        if not smi:
            removed += 1
            continue
        try:
            mol = Chem.MolFromSmiles(str(smi))
        except Exception:
            mol = None
        if mol is None:
            removed += 1
            continue
        symbols = {atom.GetSymbol().upper() for atom in mol.GetAtoms()}
        if symbols and symbols.issubset(allowed):
            filtered.append(row)
        else:
            removed += 1

    if removed:
        logger.info(
            "[pocket-eval] calibrator_filter unsupported_atoms removed=%d kept=%d",
            removed,
            len(filtered),
        )
    return filtered


_VINA_ALLOWED_TYPES = {
    "C",
    "A",
    "N",
    "NA",
    "O",
    "OA",
    "S",
    "SA",
    "H",
    "HD",
    "F",
    "CL",
    "BR",
    "I",
    "P",
    "SI",
    "SE",
    "ZN",
    "MG",
    "CA",
    "MN",
    "FE",
    "K",
    "CU",
    "CO",
    "NI",
    "AL",
    "AG",
    "AU",
    "PT",
    "LI",
    "BA",
    "SR",
    "CS",
    "RB",
}


def _pdbqt_has_only_vina_types(
    pdbqt_path: Path, logger: logging.Logger
) -> tuple[bool, Optional[str]]:
    try:
        with pdbqt_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for raw in handle:
                if not raw.startswith(("ATOM", "HETATM")):
                    continue
                parts = raw.split()
                if len(parts) < 3:
                    continue
                adt = parts[-1].strip().upper()
                while adt and not adt[-1].isalnum():
                    adt = adt[:-1]
                if adt in {"NA+", "NA_", "NA"}:
                    adt = "NA"
                if adt and adt not in _VINA_ALLOWED_TYPES:
                    logger.info(
                        "[pocket-eval] calibrator_pdbqt_invalid type=%s path=%s",
                        adt,
                        pdbqt_path,
                    )
                    return False, adt
    except Exception:
        return False, "read_error"
    return True, None


def _sample_rows(
    rows: List[Dict[str, Any]], n_keep: int, rng: random.Random
) -> List[Dict[str, Any]]:
    if n_keep <= 0 or len(rows) <= n_keep:
        return rows
    indices = sorted(rng.sample(range(len(rows)), n_keep))
    return [rows[i] for i in indices]


def _limit_calibrators(
    strong_rows: List[Dict[str, Any]],
    non_rows: List[Dict[str, Any]],
    max_total: Optional[int],
    seed: int,
    logger: logging.Logger,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    if not max_total or max_total <= 0:
        return strong_rows, non_rows
    total = len(strong_rows) + len(non_rows)
    if total <= max_total:
        return strong_rows, non_rows

    rng = random.Random(seed)
    strong_rows = sorted(strong_rows, key=lambda r: r["ligand_id"])
    non_rows = sorted(non_rows, key=lambda r: r["ligand_id"])

    target_strong = min(len(strong_rows), max_total // 2)
    target_non = min(len(non_rows), max_total - target_strong)

    if target_strong == 0 and strong_rows:
        target_strong = 1
    if target_non == 0 and non_rows:
        target_non = 1

    while target_strong + target_non > max_total:
        if target_strong >= target_non and target_strong > 0:
            target_strong -= 1
        elif target_non > 0:
            target_non -= 1
        else:
            break

    while target_strong + target_non < max_total:
        if len(strong_rows) - target_strong > len(non_rows) - target_non:
            if target_strong < len(strong_rows):
                target_strong += 1
            else:
                break
        else:
            if target_non < len(non_rows):
                target_non += 1
            else:
                break

    logger.info(
        "[pocket-eval] calibrator_sample max_total=%d keep_strong=%d keep_non=%d",
        max_total,
        target_strong,
        target_non,
    )
    return (
        _sample_rows(strong_rows, target_strong, rng),
        _sample_rows(non_rows, target_non, rng),
    )


def get_calibration_set_for_pdb(
    pdb_id: str, cfg: Dict[str, Any], logger: Optional[logging.Logger] = None
) -> List[Dict[str, Any]]:
    log = logger or _LOG
    pdb_norm = pdb_id.upper()
    paths = make_paths(cfg, base_id=pdb_norm, pdb_file=f"{pdb_norm}.pdb")
    dock_root = paths.docked_pdb_root() / "pocket_eval"
    dock_root.mkdir(parents=True, exist_ok=True)
    cache_path = dock_root / "calibrator_cache.json"
    cache_dir = Path(cfg.get("CALIBRATOR_CACHE_DIR", REPO_ROOT / "calibrator" / ".cache"))
    force_calibrator = bool(cfg.get("FORCE_CALIBRATOR", False))

    cached_rows = None
    if force_calibrator:
        log.info(
            "[pocket-eval] calibrator_cache=force_invalidate path=%s", cache_path
        )
        try:
            cache_path.unlink(missing_ok=True)
        except Exception:
            try:
                cache_path.unlink()
            except Exception:
                pass
    else:
        cached_rows = _load_calibrator_cache(cache_path, log)
        if cached_rows:
            log.info(
                "[pocket-eval] calibrator_cache=hit rows=%d path=%s",
                len(cached_rows),
                cache_path,
            )
            return cached_rows

    extracted_root = Path(cfg.get("LIGAND_EXTRACTED_DIR", DEFAULT_OUT_ROOT))
    extracted_dir = extracted_root / f"{pdb_norm}_calibrator"
    labels = None
    if force_calibrator:
        log.info(
            "[pocket-eval] calibrator_extracted=force_invalidate dir=%s", extracted_dir
        )
        shutil.rmtree(extracted_dir, ignore_errors=True)
    else:
        labels = _load_extracted_calibrator_labels(extracted_dir)
    if labels is None:
        try:
            run_calibrator_for_pdb(
                pdb_norm,
                out_root=extracted_root,
                cache_dir=cache_dir,
                timeout=DEFAULT_TIMEOUT,
                retries=DEFAULT_RETRIES,
                activity_types=list(DEFAULT_ACTIVITY_TYPES),
                chembl_max_phase=DEFAULT_CHEMBL_MAX_PHASE,
                label_thresholds=DEFAULT_LABEL_THRESHOLDS,
                run_tag=cfg.get("RUN_ID") or os.environ.get("ATLAS_RUN_ID"),
                force_refresh=force_calibrator,
            )
        except Exception as exc:
            log.warning(
                "[pocket-eval] calibrator_run_failed pdb=%s err=%s", pdb_norm, exc
            )
        labels = _load_extracted_calibrator_labels(extracted_dir)

    if labels is not None:
        rows = _build_rows_from_labels(labels)
        payload = {
            "pdb_id": pdb_norm,
            "rows": rows,
            "labels": labels,
            "source": "extracted_ligands",
        }
        _write_calibrator_cache(cache_path, payload)
        log.info(
            "[pocket-eval] calibrator_cache=extracted rows=%d path=%s",
            len(rows),
            cache_path,
        )
        return rows

    chain_map = resolve_chain_uniprot_segments(pdb_norm, write_file=False)
    _, segment = select_primary_chain(chain_map)
    uniprot_id = segment.get("uniprot") if isinstance(segment, dict) else None
    if not uniprot_id:
        log.warning("[pocket-eval] calibrator_uniprot_missing pdb=%s", pdb_norm)
        return []

    unp_start = segment.get("unp_start") if isinstance(segment, dict) else None
    unp_end = segment.get("unp_end") if isinstance(segment, dict) else None
    labels, chembl_meta = fetch_chembl_labeled_smiles(
        uniprot_id,
        pdb_norm,
        cache_dir,
        DEFAULT_TIMEOUT,
        DEFAULT_RETRIES,
        list(DEFAULT_ACTIVITY_TYPES),
        DEFAULT_CHEMBL_MAX_PHASE,
        DEFAULT_LABEL_THRESHOLDS,
        unp_start=unp_start,
        unp_end=unp_end,
        force_refresh=force_calibrator,
    )
    rows = _build_rows_from_labels(labels)

    payload = {
        "pdb_id": pdb_norm,
        "uniprot": uniprot_id,
        "unp_start": unp_start,
        "unp_end": unp_end,
        "rows": rows,
        "labels": labels,
        "meta": chembl_meta,
    }
    _write_calibrator_cache(cache_path, payload)
    log.info(
        "[pocket-eval] calibrator_cache=miss rows=%d path=%s", len(rows), cache_path
    )
    return rows


def prepare_calibrator_ligands(
    cal_rows: List[Dict[str, Any]],
    out_dir: Path,
    cfg: Dict[str, Any],
    logger: Optional[logging.Logger] = None,
    inputs_dir_override: Optional[Path] = None,
) -> List[Dict[str, Any]]:
    log = logger or _LOG
    out_dir.mkdir(parents=True, exist_ok=True)
    inputs_dir = Path(inputs_dir_override) if inputs_dir_override else out_dir / "inputs"
    prepped_dir = out_dir
    inputs_dir.mkdir(parents=True, exist_ok=True)
    prepped_dir.mkdir(parents=True, exist_ok=True)
    smi_path = inputs_dir / "calibrators.smi"
    sdf_path = inputs_dir / "calibrators.sdf"

    force = bool(cfg.get("FORCE_REPROCESS", False))
    sorted_rows = sorted(cal_rows, key=lambda r: r["ligand_id"])

    pdbqt_files = sorted(prepped_dir.glob("*.pdbqt"))
    cache_hit = bool(
        not force and pdbqt_files and len(pdbqt_files) >= len(sorted_rows)
    )
    if cache_hit:
        log.info(
            "[pocket-eval] ligand_prep=cache_hit n=%d dir=%s",
            len(pdbqt_files),
            prepped_dir,
        )
        if (not smi_path.exists()) or (not sdf_path.exists()):
            with smi_path.open("w", encoding="utf-8") as handle:
                for row in sorted_rows:
                    handle.write(f"{row['smiles']} {row['ligand_id']}\n")
            convert_smi_to_sdf(smi_path, sdf_path, require_output=True)
    else:
        with smi_path.open("w", encoding="utf-8") as handle:
            for row in sorted_rows:
                handle.write(f"{row['smiles']} {row['ligand_id']}\n")
        convert_smi_to_sdf(smi_path, sdf_path, require_output=True)
        prev_in_sdf = os.environ.get("LIGPREP_IN_SDF")
        os.environ["LIGPREP_IN_SDF"] = str(sdf_path)
        try:
            prep_ligands_with_mgltools(
                force=force,
                microstate_dedup=False,
                ph_values=None,
                root_dir=prepped_dir,
            )
        finally:
            if prev_in_sdf is None:
                os.environ.pop("LIGPREP_IN_SDF", None)
            else:
                os.environ["LIGPREP_IN_SDF"] = prev_in_sdf
        pdbqt_files = sorted(prepped_dir.glob("*.pdbqt"))

    prepared: List[Dict[str, Any]] = []
    missing = 0
    expected_by_id = {
        row["ligand_id"]: prepped_dir
        / f"{_sanitize_ligand_name_for_filename(row['ligand_id'])}.pdbqt"
        for row in sorted_rows
    }
    if expected_by_id and all(p.exists() for p in expected_by_id.values()):
        mapping = expected_by_id
    elif pdbqt_files and len(pdbqt_files) >= len(sorted_rows):
        mapping = {
            row["ligand_id"]: pdbqt_files[idx]
            for idx, row in enumerate(sorted_rows)
        }
    else:
        mapping = expected_by_id

    for row in sorted_rows:
        path = mapping.get(row["ligand_id"])
        if not path or not path.exists():
            missing += 1
            continue
        prepared.append(
            {
                "ligand_id": row["ligand_id"],
                "label": row["label"],
                "pdbqt_path": str(path),
            }
        )

    invalid = 0
    if prepared:
        validated: List[Dict[str, Any]] = []
        for row in prepared:
            pdbqt_path = Path(row["pdbqt_path"])
            ok, _bad = _pdbqt_has_only_vina_types(pdbqt_path, log)
            if ok:
                validated.append(row)
            else:
                invalid += 1
        prepared = validated

    log.info(
        "[pocket-eval] ligand_prep complete prepared=%d missing=%d invalid=%d dir=%s",
        len(prepared),
        missing,
        invalid,
        prepped_dir,
    )
    return prepared


def _write_vina_config(
    config_path: Path,
    receptor_pdbqt: Path,
    ligand_path: Path,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    *,
    cpu: int,
    exhaustiveness: int,
    num_modes: int,
    energy_range: int = 4,
    verbosity: int = 0,
    seed: Optional[int] = None,
    out_path: Optional[Path] = None,
) -> None:
    lines = [
        f"receptor = {receptor_pdbqt}",
        f"ligand   = {ligand_path}",
        f"center_x = {center[0]:.3f}",
        f"center_y = {center[1]:.3f}",
        f"center_z = {center[2]:.3f}",
        f"size_x   = {box_size[0]:.3f}",
        f"size_y   = {box_size[1]:.3f}",
        f"size_z   = {box_size[2]:.3f}",
        f"cpu      = {int(cpu)}",
        f"exhaustiveness = {int(exhaustiveness)}",
        f"energy_range   = {int(energy_range)}",
        f"num_modes      = {int(num_modes)}",
        f"verbosity      = {int(verbosity)}",
    ]
    if out_path is not None:
        lines.append(f"out = {out_path}")
    if seed is not None:
        lines.append(f"seed = {int(seed)}")

    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _read_scores_csv(path: Path) -> Optional[List[Dict[str, Any]]]:
    if not path.exists():
        return None
    try:
        if path.stat().st_size == 0:
            return None
    except Exception:
        return None

    rows = []
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                ligand_id = row.get("ligand_id")
                label = row.get("label")
                score_raw = row.get("score")
                if not ligand_id or not label or score_raw in (None, ""):
                    continue
                try:
                    score = float(score_raw)
                except Exception:
                    continue
                rows.append({"ligand_id": ligand_id, "label": label, "score": score})
    except Exception:
        return None
    return rows if rows else None


def _write_scores_csv(path: Path, rows: Iterable[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["ligand_id", "score", "label"])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "ligand_id": row["ligand_id"],
                    "score": row["score"],
                    "label": row["label"],
                }
            )


def _compute_fold_metrics(
    rows: List[Dict[str, Any]], folds: int, seed: int
) -> Dict[str, Any]:
    y_true = np.array([1 if r["label"] == "strong" else 0 for r in rows], dtype=int)
    scores = np.array([float(r["score"]) for r in rows], dtype=float)
    n_strong = int(y_true.sum())
    n_non = int(len(y_true) - n_strong)

    metrics: Dict[str, Any] = {
        "n_strong": n_strong,
        "n_non": n_non,
        "n_total_used": int(len(y_true)),
    }

    min_class = min(n_strong, n_non)
    if min_class <= 0:
        metrics["error"] = "missing_class"
        return metrics

    if min_class < 2 or folds < 2:
        try:
            metrics["auc_mean"] = float(dud_eval.roc_auc_score(y_true, -scores))
        except Exception:
            metrics["auc_mean"] = float("nan")
        metrics["auc_std"] = 0.0
        ef = dud_eval.ef_at_fractions(y_true, scores, fractions=(0.01,))
        metrics["ef1_mean"] = float(ef.get("EF@1%", float("nan")))
        metrics["ef1_std"] = 0.0
        return metrics

    effective_folds = min(int(folds), min_class)
    skf = StratifiedKFold(n_splits=effective_folds, shuffle=True, random_state=seed)
    auc_vals: List[float] = []
    ef_vals: List[float] = []

    for _, test_idx in skf.split(scores, y_true):
        y_fold = y_true[test_idx]
        score_low = scores[test_idx]
        try:
            auc_val = float(dud_eval.roc_auc_score(y_fold, -score_low))
        except Exception:
            auc_val = float("nan")
        auc_vals.append(auc_val)

        ef = dud_eval.ef_at_fractions(y_fold, score_low, fractions=(0.01,))
        ef_val = float(ef.get("EF@1%", float("nan")))
        ef_vals.append(ef_val)

    metrics["auc_mean"] = float(np.nanmean(auc_vals)) if auc_vals else float("nan")
    metrics["auc_std"] = float(np.nanstd(auc_vals)) if auc_vals else float("nan")
    metrics["ef1_mean"] = float(np.nanmean(ef_vals)) if ef_vals else float("nan")
    metrics["ef1_std"] = float(np.nanstd(ef_vals)) if ef_vals else float("nan")
    return metrics


def _metric_or_default(value: Any, default: float) -> float:
    if value is None:
        return default
    try:
        val = float(value)
    except Exception:
        return default
    if math.isnan(val):
        return default
    return val


def _select_best_pocket(pocket_rows: List[Dict[str, Any]]) -> Tuple[Optional[Dict[str, Any]], str]:
    if not pocket_rows:
        return None, "no_pockets"

    def _sort_key(entry: Dict[str, Any]) -> Tuple[float, float, float, str]:
        metrics = entry.get("metrics") or {}
        auc = _metric_or_default(metrics.get("auc_mean"), float("-inf"))
        ef = _metric_or_default(metrics.get("ef1_mean"), float("-inf"))
        std = _metric_or_default(metrics.get("auc_std"), float("inf"))
        pocket_id = str(entry.get("pocket_id") or "")
        return (-auc, -ef, std, pocket_id)

    ranked = sorted(pocket_rows, key=_sort_key)
    if not ranked:
        return None, "no_pockets"
    best = ranked[0]
    best_auc = _metric_or_default(
        (best.get("metrics") or {}).get("auc_mean"), float("-inf")
    )
    if best_auc == float("-inf"):
        return None, "no_valid_metrics"
    return best, "max_auc_mean"


def _write_performance_error(
    path: Path,
    pdb_id: str,
    reason: str,
    *,
    seed: int,
    folds: int,
    pockets: Optional[List[Dict[str, Any]]] = None,
) -> None:
    payload = {
        "pdb_id": pdb_id,
        "mode": "pocket_eval",
        "labels_used": {"positive": "strong", "negative": "non", "ignored": ["weak"]},
        "split": {"seed": seed, "folds": folds},
        "pockets": pockets or [],
        "selected_pocket_id": None,
        "selected_reason": None,
        "error": reason,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _select_calibrator_rows_for_eval_and_prep(
    cal_rows: List[Dict[str, Any]],
    max_total: Optional[int],
    seed: int,
    logger: logging.Logger,
) -> Tuple[List[Dict[str, Any]], int]:
    strong_rows = [r for r in cal_rows if r.get("label") == "strong"]
    non_rows = [r for r in cal_rows if r.get("label") == "non"]
    weak_rows = [r for r in cal_rows if r.get("label") == "weak"]

    if not strong_rows or not non_rows:
        return [], len(weak_rows)

    strong_rows_used, non_rows_used = _limit_calibrators(
        strong_rows, non_rows, max_total, seed, logger
    )
    rows_used = strong_rows_used + non_rows_used
    return rows_used, len(weak_rows)


def select_pocket_with_eval(
    *,
    pdb_id: str,
    cfg: Dict[str, Any],
    logger: Optional[logging.Logger] = None,
    pockets_json_path: Path,
    receptor_pdbqt_path: Path,
) -> Optional[Dict[str, Any]]:
    log = logger or _LOG
    pdb_norm = pdb_id.upper()
    paths = make_paths(cfg, base_id=pdb_norm, pdb_file=f"{pdb_norm}.pdb")
    dock_root = paths.docked_pdb_root() / "pocket_eval"
    dock_root.mkdir(parents=True, exist_ok=True)
    perf_path = dock_root / "pocket_performance.json"
    seed = int(cfg.get("POCKET_EVAL_SEED", 0))
    folds = int(cfg.get("POCKET_EVAL_FOLDS", 5))
    max_cal = cfg.get("POCKET_EVAL_MAX_CALIBRATORS")

    log.info(
        "[pocket-eval] start pdb=%s pockets_json=%s", pdb_norm, pockets_json_path
    )

    cal_rows = get_calibration_set_for_pdb(pdb_norm, cfg, log)
    if not cal_rows:
        log.warning("[pocket-eval] no_calibrators pdb=%s", pdb_norm)
        _write_performance_error(perf_path, pdb_norm, "calibrators_empty", seed=seed, folds=folds)
        return None
    cal_rows = _filter_supported_calibrators(cal_rows, log)
    if not cal_rows:
        log.warning("[pocket-eval] calibrators_filtered_empty pdb=%s", pdb_norm)
        _write_performance_error(perf_path, pdb_norm, "calibrators_empty", seed=seed, folds=folds)
        return None

    rows_used, weak_ignored = _select_calibrator_rows_for_eval_and_prep(
        cal_rows, max_cal, seed, log
    )
    strong_used = sum(1 for r in rows_used if r.get("label") == "strong")
    non_used = sum(1 for r in rows_used if r.get("label") == "non")
    if not strong_used or not non_used:
        reason = "calibrators_missing_class"
        log.warning(
            "[pocket-eval] %s strong=%d non=%d weak=%d",
            reason,
            strong_used,
            non_used,
            weak_ignored,
        )
        _write_performance_error(perf_path, pdb_norm, reason, seed=seed, folds=folds)
        return None

    cal_rows_used = rows_used
    prep_rows = rows_used
    log.info(
        "[pocket-eval] calibrators used strong=%d non=%d weak_ignored=%d",
        strong_used,
        non_used,
        weak_ignored,
    )

    prep_root = paths.prepped_ligands_dir.parent / f"{pdb_norm}_calibrator"
    prep_dir = prep_root
    extracted_root = Path(cfg.get("LIGAND_EXTRACTED_DIR", DEFAULT_OUT_ROOT))
    extracted_dir = extracted_root / f"{pdb_norm}_calibrator"
    extracted_dir.mkdir(parents=True, exist_ok=True)
    prepared = prepare_calibrator_ligands(
        prep_rows, prep_dir, cfg, log, inputs_dir_override=extracted_dir
    )
    used_ids = {r["ligand_id"] for r in cal_rows_used}
    prepared_scoring = [r for r in prepared if r.get("ligand_id") in used_ids]
    strong_prepped = [r for r in prepared_scoring if r.get("label") == "strong"]
    non_prepped = [r for r in prepared_scoring if r.get("label") == "non"]
    if not strong_prepped or not non_prepped:
        reason = "calibrators_prep_missing_class"
        log.warning(
            "[pocket-eval] %s strong=%d non=%d",
            reason,
            len(strong_prepped),
            len(non_prepped),
        )
        _write_performance_error(perf_path, pdb_norm, reason, seed=seed, folds=folds)
        return None

    if not pockets_json_path.exists():
        reason = "pockets_json_missing"
        log.warning("[pocket-eval] %s path=%s", reason, pockets_json_path)
        _write_performance_error(perf_path, pdb_norm, reason, seed=seed, folds=folds)
        return None

    try:
        pockets_payload = json.loads(pockets_json_path.read_text(encoding="utf-8"))
    except Exception as exc:
        reason = "pockets_json_unreadable"
        log.warning("[pocket-eval] %s err=%s", reason, exc)
        _write_performance_error(perf_path, pdb_norm, reason, seed=seed, folds=folds)
        return None

    pockets = pockets_payload.get("pockets") if isinstance(pockets_payload, dict) else None
    if not pockets:
        reason = "pockets_empty"
        log.warning("[pocket-eval] %s", reason)
        _write_performance_error(perf_path, pdb_norm, reason, seed=seed, folds=folds)
        return None

    vina_exe = cfg.get("VINA_EXE") or cfg.get("VINA_PATH") or "vina"
    threads = max(1, int(cfg.get("THREADS_PER_VINA", 1)))
    exhaustiveness = 1
    num_modes = 1
    verbosity = int(cfg.get("VINA_VERBOSITY", 0))

    pocket_rows: List[Dict[str, Any]] = []
    for pocket in pockets:
        if not isinstance(pocket, dict):
            continue
        pocket_id = pocket_id_from_entry(pocket)
        center, box_size = _center_and_box_from_pocket(pocket)
        if not center or not box_size:
            log.warning(
                "[pocket-eval] pocket_skip missing_box pocket_id=%s", pocket_id
            )
            continue

        dock_dir = dock_root / pocket_id
        scores_path = dock_dir / "calibration_scores.csv"
        scores_rows = None
        if not cfg.get("FORCE_REPROCESS"):
            cached = _read_scores_csv(scores_path)
            if cached:
                expected_ids = {r["ligand_id"] for r in prepared_scoring}
                cached_ids = {r["ligand_id"] for r in cached}
                if expected_ids.issubset(cached_ids):
                    scores_rows = cached
                    log.info(
                        "[pocket-eval] scores_cache_hit pocket_id=%s n=%d",
                        pocket_id,
                        len(scores_rows),
                    )

        if scores_rows is None:
            scores_rows = []
            for entry in prepared_scoring:
                ligand_id = entry["ligand_id"]
                ligand_path = Path(entry["pdbqt_path"])
                safe_id = _sanitize_ligand_name_for_filename(ligand_id)
                config_path = dock_dir / "configs" / f"{safe_id}.txt"
                out_path = dock_dir / f"{safe_id}.pdbqt"

                _write_vina_config(
                    config_path,
                    receptor_pdbqt_path,
                    ligand_path,
                    center,
                    box_size,
                    cpu=threads,
                    exhaustiveness=exhaustiveness,
                    num_modes=num_modes,
                    verbosity=verbosity,
                    seed=seed,
                    out_path=out_path,
                )
                _, score = run_docking_task(
                    str(vina_exe), str(config_path), ligand_id, str(out_path)
                )
                if score is None:
                    continue
                scores_rows.append(
                    {"ligand_id": ligand_id, "label": entry["label"], "score": score}
                )

            _write_scores_csv(scores_path, scores_rows)

        metrics = _compute_fold_metrics(scores_rows, folds, seed)
        pocket_rows.append(
            {
                "pocket_id": pocket_id,
                "ligand_resname": pocket.get("ligand_resname"),
                "chain": pocket.get("ligand_chain") or pocket.get("protein_chain"),
                "resnum": pocket.get("ligand_resnum"),
                "metrics": metrics,
                "artifacts": {
                    "dock_dir": str(dock_dir),
                    "scores_path": str(scores_path),
                },
                "center": center,
                "box_size": box_size,
            }
        )
        log.info(
            "[pocket-eval] pocket=%s auc_mean=%s ef1_mean=%s",
            pocket_id,
            metrics.get("auc_mean"),
            metrics.get("ef1_mean"),
        )

    selected, reason = _select_best_pocket(pocket_rows)
    if not selected:
        _write_performance_error(perf_path, pdb_norm, "no_valid_pockets", seed=seed, folds=folds, pockets=pocket_rows)
        return None

    performance_payload = {
        "pdb_id": pdb_norm,
        "mode": "pocket_eval",
        "labels_used": {"positive": "strong", "negative": "non", "ignored": ["weak"]},
        "split": {"seed": seed, "folds": folds},
        "pockets": [
            {
                "pocket_id": p["pocket_id"],
                "ligand_resname": p.get("ligand_resname"),
                "chain": p.get("chain"),
                "resnum": p.get("resnum"),
                "metrics": p.get("metrics"),
                "artifacts": p.get("artifacts"),
            }
            for p in pocket_rows
        ],
        "selected_pocket_id": selected.get("pocket_id"),
        "selected_reason": reason,
    }
    perf_path.parent.mkdir(parents=True, exist_ok=True)
    perf_path.write_text(json.dumps(performance_payload, indent=2), encoding="utf-8")

    log.info(
        "[pocket-eval] selected pocket_id=%s reason=%s perf=%s",
        selected.get("pocket_id"),
        reason,
        perf_path,
    )
    return {
        "pocket_id": selected.get("pocket_id"),
        "center": selected.get("center"),
        "box_size": selected.get("box_size"),
        "selection_reason": reason,
        "performance_path": str(perf_path),
    }
