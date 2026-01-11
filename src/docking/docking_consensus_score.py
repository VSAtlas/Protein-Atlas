# Consensus equations (rank-based)
# N = number of candidate ligands for this (variant, pH)
# E = set of engines/signals (e.g. vina, gnina_energy, gnina_cnn, ledock, dock6)
# r_e(l) = rank of ligand l under engine e, in {1, 2, ..., N} with 1 = best, N = worst
# Rank percentile in [0, 1]:
# p_e(l) = 1 - (r_e(l) - 1) / (N - 1)
# = (N - r_e(l)) / (N - 1)
# Best ligand under engine e: r_e(l) = 1 -> p_e(l) = 1
# Worst or missing ligand: r_e(l) = N -> p_e(l) = 0
# Consensus score (higher is better):
# S(l) = sum_e ( w_e * p_e(l) ) where weights w_e >= 0 and the sum over participating engines is 1.
# (Implementation groups engines into energy, ledock, and cnn families.)
# Energy-family combination (per ligand):
# alpha_vina = 0.4
# alpha_dock6 = 0.6
# alpha_gnina_energy = 0.4
# p_energy(l) = (alpha_vina * p_vina(l)
#              + alpha_dock6 * p_dock6(l)
#              + alpha_gnina_energy * p_gnina_energy(l)) / denom(l)
# denom(l) = sum of alpha_x for engines x in {vina, dock6, gnina_energy} that have a valid score for ligand l.
# If no energy engines have a score for l, set p_energy(l) = 0.
# Full runs (CNN available):
# w_energy = 0.25
# w_ledock = 0.30
# w_cnn = 0.45
# S(l) = w_energy * p_energy(l)
# + w_ledock * p_ledock(l)
# + w_cnn * p_cnn(l)
# Fast runs (CNN missing):
# w_energy = 0.45
# w_ledock = 0.55
# w_cnn = 0.0
# In all cases, only engines that actually have data for this run participate:
# - If an engine has no usable scores for this (variant, pH), its base weight is set to 0
#   and the remaining weights are renormalized to sum to 1.
# - For a specific ligand l, if an engine e has no valid score, p_e(l) = 0 (worst rank).

from __future__ import annotations

import csv
import logging
import math
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence


def _normalize_ligand_id(lig: str) -> str:
    """Normalize ligand identifiers across engines by stripping one known extension."""
    name = Path(lig).name
    for ext in (".pdbqt", ".mol2", ".dok", ".sdf", ".pdb"):
        if name.lower().endswith(ext):
            return name[: -len(ext)]
    return name


def _to_bool_flag(raw: Any) -> bool:
    if raw is None:
        return False
    if isinstance(raw, bool):
        return raw
    token = str(raw).strip().lower()
    return token in {"1", "true", "t", "yes", "y"}


FDA_PREFIX_DEFAULTS = ["fda_"]
DECOY_PREFIX_DEFAULTS = ["dud_", "deepcoy_", "decoy_", "decoys"]


def _library_prefixes(
    cfg: Dict[str, Any], key: str, defaults: Sequence[str]
) -> List[str]:
    raw = cfg.get(key)
    if raw is None:
        raw = cfg.get(key.lower())
    if raw is None:
        return list(defaults)
    if isinstance(raw, str):
        parts = [p.strip() for p in re.split(r"[;,]", raw) if p.strip()]
        return parts or list(defaults)
    try:
        return [str(p).strip() for p in raw if str(p).strip()]
    except Exception:
        return list(defaults)


def _infer_library(
    ligand_display: str, fda_prefixes: Sequence[str], decoy_prefixes: Sequence[str]
) -> str:
    name = Path(ligand_display).name.lower()
    for prefix in fda_prefixes:
        if name.startswith(prefix.lower()):
            return "FDA"
    for prefix in decoy_prefixes:
        if name.startswith(prefix.lower()):
            return "DECOY"
    return "UNKNOWN"


def _collect_ligands_from_csv(
    csv_path: Path,
    ligand_ids: List[str],
    display_name_by_id: Dict[str, str],
    seen: set[str],
    logger: logging.Logger,
    *,
    label: str,
) -> None:
    if not csv_path.exists():
        return
    overlaps = 0
    added = 0
    try:
        with csv_path.open("r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                lig = row.get("ligand")
                if not lig:
                    continue
                lig_id = _normalize_ligand_id(lig)
                if not lig_id:
                    continue
                if lig_id in seen:
                    overlaps += 1
                    continue
                seen.add(lig_id)
                ligand_ids.append(lig_id)
                display_name_by_id[lig_id] = lig
                added += 1
    except Exception as exc:
        logger.warning(
            "[consensus.skip] reason=read_error path=%s label=%s error=%s",
            str(csv_path),
            label,
            exc,
        )
        return
    if overlaps:
        logger.warning(
            "[consensus.decoy_overlap] path=%s label=%s overlaps=%d action=prefer_existing",
            str(csv_path),
            label,
            overlaps,
        )
    if added:
        logger.info(
            "[consensus.extend] path=%s label=%s added=%d", str(csv_path), label, added
        )


def _merge_scores(
    base_scores: Dict[str, float],
    incoming: Dict[str, float],
    *,
    higher_is_better: bool,
    logger: logging.Logger,
    label: str,
) -> Dict[str, float]:
    if not incoming:
        return base_scores
    overlaps = set(base_scores).intersection(incoming)
    if overlaps:
        logger.warning(
            "[consensus.decoy_overlap] label=%s overlaps=%d action=keep_best",
            label,
            len(overlaps),
        )
    for lig, val in incoming.items():
        existing = base_scores.get(lig)
        if existing is None:
            base_scores[lig] = val
            continue
        if higher_is_better:
            if val > existing:
                base_scores[lig] = val
        else:
            if val < existing:
                base_scores[lig] = val
    return base_scores


def _detect_cnn_key(csv_path: Path) -> Optional[str]:
    if not csv_path.exists():
        return None
    try:
        with csv_path.open("r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            fields = reader.fieldnames or []
            if "gnina_cnn_affinity_pK" in fields:
                return "gnina_cnn_affinity_pK"
            if "gnina_cnn_score" in fields:
                return "gnina_cnn_score"
    except Exception:
        return None
    return None


def _load_best_scores_per_ligand(
    csv_path: Path,
    value_key: Optional[str],
    valid_key: Optional[str],
    *,
    higher_is_better: bool,
    logger: logging.Logger,
) -> Dict[str, float]:
    scores: Dict[str, float] = {}
    if not csv_path.exists() or value_key is None:
        logger.info(
            "[consensus.skip] reason=missing_csv_or_key path=%s key=%s",
            str(csv_path),
            value_key,
        )
        return scores

    with csv_path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ligand_raw = row.get("ligand")
            lig_id = _normalize_ligand_id(ligand_raw) if ligand_raw else ""
            if not lig_id:
                continue

            if valid_key is not None and not _to_bool_flag(row.get(valid_key)):
                continue

            raw = row.get(value_key)
            if raw is None or str(raw).strip() == "":
                continue
            try:
                val = float(raw)
            except Exception:
                continue
            if not math.isfinite(val):
                continue

            if lig_id not in scores:
                scores[lig_id] = val
            else:
                if higher_is_better:
                    if val > scores[lig_id]:
                        scores[lig_id] = val
                else:
                    if val < scores[lig_id]:
                        scores[lig_id] = val

    logger.info(
        "[consensus.load] path=%s key=%s n_ligands=%d",
        str(csv_path),
        value_key,
        len(scores),
    )
    return scores


def _percentiles_from_scores(
    ligands: List[str],
    scores: Dict[str, float],
    *,
    higher_is_better: bool,
) -> Dict[str, float]:
    n = len(ligands)
    if n <= 1:
        return {lig: 1.0 for lig in ligands}

    ranked_items = []
    for lig, sc in scores.items():
        if isinstance(sc, (int, float)) and math.isfinite(sc):
            ranked_items.append((sc, lig))

    if not ranked_items:
        return {lig: 0.0 for lig in ligands}

    ranked_items.sort(key=lambda x: x[0], reverse=higher_is_better)
    ranks: Dict[str, int] = {}
    for idx, (_, lig) in enumerate(ranked_items, start=1):
        ranks[lig] = idx

    denom = max(n - 1, 1)
    percentiles: Dict[str, float] = {}
    for lig in ligands:
        r_e = ranks.get(lig, n)
        percentiles[lig] = 1.0 - (r_e - 1) / denom
    return percentiles


def _pick_best_engine(
    lig_id: str,
    p_vina: Dict[str, float],
    p_ledock: Dict[str, float],
    p_dock6: Dict[str, float],
    p_gnina_energy: Dict[str, float],
    p_cnn: Dict[str, float],
    vina_scores: Dict[str, float],
    ledock_scores: Dict[str, float],
    dock6_scores: Dict[str, float],
    gnina_energy_scores: Dict[str, float],
    cnn_scores: Dict[str, float],
) -> (str, str):
    priority = {
        "gnina_cnn": 5,
        "ledock": 4,
        "dock6": 3,
        "vina": 2,
        "gnina_energy": 1,
    }

    candidates = []
    if lig_id in cnn_scores:
        candidates.append(
            (p_cnn.get(lig_id, 0.0), priority["gnina_cnn"], "gnina", "gnina_cnn")
        )
    if lig_id in ledock_scores:
        candidates.append(
            (p_ledock.get(lig_id, 0.0), priority["ledock"], "ledock", "ledock")
        )
    if lig_id in dock6_scores:
        candidates.append(
            (p_dock6.get(lig_id, 0.0), priority["dock6"], "dock6", "dock6")
        )
    if lig_id in vina_scores:
        candidates.append((p_vina.get(lig_id, 0.0), priority["vina"], "vina", "vina"))
    if lig_id in gnina_energy_scores:
        candidates.append(
            (
                p_gnina_energy.get(lig_id, 0.0),
                priority["gnina_energy"],
                "gnina",
                "gnina_energy",
            )
        )

    if not candidates:
        return "", ""

    best_p, _, best_engine, best_signal = max(candidates, key=lambda x: (x[0], x[1]))
    return (
        best_engine if best_p is not None else "",
        best_signal if best_p is not None else "",
    )


def compute_consensus_for_variant_ph(
    cfg: Dict[str, Any],
    paths: Any,
    ph_label: Optional[str],
    *,
    variant_env: Optional[str],
    variant_label: str,
    csv_prefix: str,
    logger: Optional[logging.Logger] = None,
) -> Optional[Path]:
    log = logger or logging.getLogger(__name__)
    variant_root = Path(paths.docked_variant_root(variant_env, ph_label))

    vina_csv = variant_root / f"{csv_prefix}docking_score_long.csv"
    gnina_csv = variant_root / f"{csv_prefix}gnina_docking_score_long.csv"
    ledock_csv = variant_root / f"{csv_prefix}ledock_docking_score_long.csv"
    dock6_csv = variant_root / f"{csv_prefix}dock6_docking_score_long.csv"
    decoy_vina_csv = variant_root / f"{csv_prefix}dud_docking_score_long.csv"
    decoy_gnina_csv = variant_root / f"{csv_prefix}dud_gnina_docking_score_long.csv"
    decoy_ledock_csv = variant_root / f"{csv_prefix}dud_ledock_docking_score_long.csv"
    decoy_dock6_csv = variant_root / f"{csv_prefix}dud_dock6_docking_score_long.csv"

    if not vina_csv.exists():
        log.info("[consensus.skip] reason=missing_vina_csv path=%s", str(vina_csv))
        return None

    ligand_ids: List[str] = []
    display_name_by_id: Dict[str, str] = {}
    seen = set()
    with vina_csv.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            lig = row.get("ligand")
            if not lig:
                continue
            lig_id = _normalize_ligand_id(lig)
            if lig_id and lig_id not in seen:
                seen.add(lig_id)
                ligand_ids.append(lig_id)
                display_name_by_id[lig_id] = lig
    include_decoys_raw = cfg.get("CONSENSUS_INCLUDE_DECOYS")
    if include_decoys_raw is None:
        include_decoys_raw = cfg.get("consensus_include_decoys")
    include_decoys = (
        True if include_decoys_raw is None else _to_bool_flag(include_decoys_raw)
    )
    fda_prefixes = _library_prefixes(cfg, "FDA_LIGAND_PREFIXES", FDA_PREFIX_DEFAULTS)
    decoy_prefixes = _library_prefixes(
        cfg, "DECOY_LIGAND_PREFIXES", DECOY_PREFIX_DEFAULTS
    )

    if include_decoys:
        _collect_ligands_from_csv(
            decoy_vina_csv,
            ligand_ids,
            display_name_by_id,
            seen,
            log,
            label="vina_decoys",
        )
        _collect_ligands_from_csv(
            decoy_gnina_csv,
            ligand_ids,
            display_name_by_id,
            seen,
            log,
            label="gnina_decoys",
        )
        _collect_ligands_from_csv(
            decoy_ledock_csv,
            ligand_ids,
            display_name_by_id,
            seen,
            log,
            label="ledock_decoys",
        )
        _collect_ligands_from_csv(
            decoy_dock6_csv,
            ligand_ids,
            display_name_by_id,
            seen,
            log,
            label="dock6_decoys",
        )
    if not ligand_ids:
        log.warning("[consensus.skip] reason=no_ligands path=%s", str(vina_csv))
        return None

    run_id = cfg.get("RUN_ID", "")

    vina_scores = _load_best_scores_per_ligand(
        vina_csv,
        value_key="score",
        valid_key="valid",
        higher_is_better=False,
        logger=log,
    )
    gnina_energy_scores = _load_best_scores_per_ligand(
        gnina_csv,
        value_key="gnina_minimized_affinity_kcal",
        valid_key="valid",
        higher_is_better=False,
        logger=log,
    )

    cnn_key: Optional[str] = _detect_cnn_key(gnina_csv)
    cnn_scores = _load_best_scores_per_ligand(
        gnina_csv,
        value_key=cnn_key,
        valid_key="valid",
        higher_is_better=True,
        logger=log,
    )

    ledock_scores = _load_best_scores_per_ligand(
        ledock_csv,
        value_key="ledock_best_score_kcal",
        valid_key=None,
        higher_is_better=False,
        logger=log,
    )
    dock6_scores = _load_best_scores_per_ligand(
        dock6_csv,
        value_key="dock6_best_score_kcal",
        valid_key=None,
        higher_is_better=False,
        logger=log,
    )
    if not dock6_scores and dock6_csv.exists():
        dock6_scores = _load_best_scores_per_ligand(
            dock6_csv,
            value_key="dock6_grid_score",
            valid_key=None,
            higher_is_better=False,
            logger=log,
        )

    if include_decoys:
        decoy_vina_scores = _load_best_scores_per_ligand(
            decoy_vina_csv,
            value_key="score",
            valid_key="valid",
            higher_is_better=False,
            logger=log,
        )
        vina_scores = _merge_scores(
            vina_scores,
            decoy_vina_scores,
            higher_is_better=False,
            logger=log,
            label="vina_decoys",
        )

        decoy_gnina_energy_scores = _load_best_scores_per_ligand(
            decoy_gnina_csv,
            value_key="gnina_minimized_affinity_kcal",
            valid_key="valid",
            higher_is_better=False,
            logger=log,
        )
        gnina_energy_scores = _merge_scores(
            gnina_energy_scores,
            decoy_gnina_energy_scores,
            higher_is_better=False,
            logger=log,
            label="gnina_energy_decoys",
        )
        decoy_cnn_key = _detect_cnn_key(decoy_gnina_csv)
        if decoy_cnn_key:
            decoy_cnn_scores = _load_best_scores_per_ligand(
                decoy_gnina_csv,
                value_key=decoy_cnn_key,
                valid_key="valid",
                higher_is_better=True,
                logger=log,
            )
            cnn_scores = _merge_scores(
                cnn_scores,
                decoy_cnn_scores,
                higher_is_better=True,
                logger=log,
                label="gnina_cnn_decoys",
            )

        decoy_ledock_scores = _load_best_scores_per_ligand(
            decoy_ledock_csv,
            value_key="ledock_best_score_kcal",
            valid_key=None,
            higher_is_better=False,
            logger=log,
        )
        ledock_scores = _merge_scores(
            ledock_scores,
            decoy_ledock_scores,
            higher_is_better=False,
            logger=log,
            label="ledock_decoys",
        )

        decoy_dock6_scores = _load_best_scores_per_ligand(
            decoy_dock6_csv,
            value_key="dock6_best_score_kcal",
            valid_key=None,
            higher_is_better=False,
            logger=log,
        )
        if not decoy_dock6_scores and decoy_dock6_csv.exists():
            decoy_dock6_scores = _load_best_scores_per_ligand(
                decoy_dock6_csv,
                value_key="dock6_grid_score",
                valid_key=None,
                higher_is_better=False,
                logger=log,
            )
        dock6_scores = _merge_scores(
            dock6_scores,
            decoy_dock6_scores,
            higher_is_better=False,
            logger=log,
            label="dock6_decoys",
        )

    p_vina = _percentiles_from_scores(ligand_ids, vina_scores, higher_is_better=False)
    p_gnina_energy = _percentiles_from_scores(
        ligand_ids, gnina_energy_scores, higher_is_better=False
    )
    p_ledock = _percentiles_from_scores(
        ligand_ids, ledock_scores, higher_is_better=False
    )
    p_cnn = _percentiles_from_scores(ligand_ids, cnn_scores, higher_is_better=True)
    p_dock6 = _percentiles_from_scores(ligand_ids, dock6_scores, higher_is_better=False)

    has_cnn = any(v > 0.0 for v in p_cnn.values())
    has_ledock = any(v > 0.0 for v in p_ledock.values())
    has_dock6 = any(v > 0.0 for v in p_dock6.values())

    alpha_vina = 0.4
    alpha_dock6 = 0.6
    alpha_gnina_energy = 0.4
    p_energy: Dict[str, float] = {}
    for lig in ligand_ids:
        numerator = 0.0
        denom = 0.0
        if lig in vina_scores:
            numerator += alpha_vina * p_vina.get(lig, 0.0)
            denom += alpha_vina
        if lig in dock6_scores:
            numerator += alpha_dock6 * p_dock6.get(lig, 0.0)
            denom += alpha_dock6
        if lig in gnina_energy_scores:
            numerator += alpha_gnina_energy * p_gnina_energy.get(lig, 0.0)
            denom += alpha_gnina_energy
        p_energy[lig] = numerator / denom if denom > 0 else 0.0

    if has_cnn:
        w_energy_base = 0.25
        w_cnn_base = 0.45
        w_ledock_base = 0.30 if has_ledock else 0.0
    else:
        w_energy_base = 0.45
        w_cnn_base = 0.0
        w_ledock_base = 0.55 if has_ledock else 0.0

    total = w_energy_base + w_cnn_base + w_ledock_base
    if total <= 0:
        w_energy = 1.0
        w_cnn = 0.0
        w_ledock = 0.0
    else:
        w_energy = w_energy_base / total
        w_cnn = w_cnn_base / total
        w_ledock = w_ledock_base / total

    log.info(
        "[consensus.weights] has_cnn=%s has_ledock=%s w_energy=%.3f w_ledock=%.3f w_cnn=%.3f",
        has_cnn,
        has_ledock,
        w_energy,
        w_ledock,
        w_cnn,
    )

    out_rows = []
    for lig in ligand_ids:
        consensus_score = (
            w_energy * p_energy.get(lig, 0.0)
            + w_ledock * p_ledock.get(lig, 0.0)
            + w_cnn * p_cnn.get(lig, 0.0)
        )
        n_engines = sum(
            1
            for store in (
                vina_scores,
                gnina_energy_scores,
                ledock_scores,
                dock6_scores,
                cnn_scores,
            )
            if lig in store
        )
        best_engine, best_signal = _pick_best_engine(
            lig,
            p_vina,
            p_ledock,
            p_dock6,
            p_gnina_energy,
            p_cnn,
            vina_scores,
            ledock_scores,
            dock6_scores,
            gnina_energy_scores,
            cnn_scores,
        )
        ligand_display = display_name_by_id.get(lig, lig)
        library = _infer_library(ligand_display, fda_prefixes, decoy_prefixes)
        out_rows.append(
            {
                "run_id": run_id,
                "pdb_id": paths.pdb_id,
                "variant": variant_label,
                "ph_label": ph_label or "",
                "ligand": ligand_display,
                "library": library,
                "consensus_score": consensus_score,
                "p_energy": p_energy.get(lig, 0.0),
                "p_ledock": p_ledock.get(lig, 0.0),
                "p_cnn": p_cnn.get(lig, 0.0),
                "p_vina": p_vina.get(lig, 0.0),
                "p_gnina_energy": p_gnina_energy.get(lig, 0.0),
                "p_dock6": p_dock6.get(lig, 0.0),
                "n_engines_with_data": n_engines,
                "best_engine": best_engine,
                "best_signal": best_signal,
                "t_vs_decoys_consensus": "",
            }
        )

    decoy_scores = [
        float(r["consensus_score"])
        for r in out_rows
        if r.get("library") == "DECOY"
        and isinstance(r.get("consensus_score"), (int, float))
        and math.isfinite(r.get("consensus_score"))
    ]
    if decoy_scores:
        mu_decoy = sum(decoy_scores) / len(decoy_scores)
        variance = sum((v - mu_decoy) ** 2 for v in decoy_scores) / len(decoy_scores)
        sigma_decoy = math.sqrt(variance)
        if sigma_decoy > 0 and math.isfinite(mu_decoy) and math.isfinite(sigma_decoy):
            for row in out_rows:
                val = row.get("consensus_score")
                if not isinstance(val, (int, float)) or not math.isfinite(val):
                    row["t_vs_decoys_consensus"] = ""
                    continue
                t_val = (val - mu_decoy) / sigma_decoy
                row["t_vs_decoys_consensus"] = f"{t_val:.6g}"
            log.info(
                "[t-score.consensus] pdb=%s variant=%s ph=%s n_decoys=%d mu=%.6g sigma=%.6g",
                paths.pdb_id,
                variant_label,
                ph_label or "",
                len(decoy_scores),
                mu_decoy,
                sigma_decoy,
            )
        else:
            log.warning(
                "[t-score.skip] pdb=%s variant=%s ph=%s reason=sigma_zero_or_nonfinite n_decoys=%d",
                paths.pdb_id,
                variant_label,
                ph_label or "",
                len(decoy_scores),
            )
    else:
        log.info(
            "[t-score.skip] pdb=%s variant=%s ph=%s reason=no_decoy_scores n_decoys=0",
            paths.pdb_id,
            variant_label,
            ph_label or "",
        )

    out_path = variant_root / "consensus_docking_scores.csv"
    if out_rows:
        fieldnames = list(out_rows[0].keys())
        with out_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(out_rows)
        log.info(
            "[consensus.write] pdb=%s variant=%s ph=%s n_rows=%d out=%s",
            paths.pdb_id,
            variant_label,
            ph_label or "",
            len(out_rows),
            str(out_path),
        )
    return out_path
