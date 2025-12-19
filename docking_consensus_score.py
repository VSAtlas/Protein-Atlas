# Consensus equations (rank-based)
# N = number of candidate ligands for this (variant, pH)
# E = set of engines/signals (e.g. vina, gnina_energy, gnina_cnn, ledock)
# r_e(l) = rank of ligand l under engine e, in {1, 2, ..., N} with 1 = best, N = worst
# Rank percentile in [0, 1]:
# p_e(l) = 1 - (r_e(l) - 1) / (N - 1)
# = (N - r_e(l)) / (N - 1)
# Best ligand under engine e: r_e(l) = 1 -> p_e(l) = 1
# Worst or missing ligand: r_e(l) = N -> p_e(l) = 0
# Consensus score (higher is better):
# S(l) = sum_e ( w_e * p_e(l) )
# with weights w_e >= 0 and sum_e w_e = 1
# Energy-family combination:
# p_energy(l) = 0.6 * p_vina(l) + 0.4 * p_gnina_energy(l)
# If gnina_energy is not available: p_energy(l) = p_vina(l)
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
from pathlib import Path
from typing import Any, Dict, List, Optional


def _to_bool_flag(raw: Any) -> bool:
    if raw is None:
        return False
    if isinstance(raw, bool):
        return raw
    token = str(raw).strip().lower()
    return token in {"1", "true", "t", "yes", "y"}


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
        logger.info("[consensus.skip] reason=missing_csv_or_key path=%s key=%s", str(csv_path), value_key)
        return scores

    with csv_path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ligand = row.get("ligand")
            if not ligand:
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

            if ligand not in scores:
                scores[ligand] = val
            else:
                if higher_is_better:
                    if val > scores[ligand]:
                        scores[ligand] = val
                else:
                    if val < scores[ligand]:
                        scores[ligand] = val

    logger.info("[consensus.load] path=%s key=%s n_ligands=%d", str(csv_path), value_key, len(scores))
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

    if not vina_csv.exists():
        log.info("[consensus.skip] reason=missing_vina_csv path=%s", str(vina_csv))
        return None

    ligands: List[str] = []
    seen = set()
    with vina_csv.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            lig = row.get("ligand")
            if lig and lig not in seen:
                seen.add(lig)
                ligands.append(lig)
    if not ligands:
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

    cnn_key: Optional[str] = None
    if gnina_csv.exists():
        with gnina_csv.open("r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            fields = reader.fieldnames or []
            if "gnina_cnn_affinity_pK" in fields:
                cnn_key = "gnina_cnn_affinity_pK"
            elif "gnina_cnn_score" in fields:
                cnn_key = "gnina_cnn_score"
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

    p_vina = _percentiles_from_scores(ligands, vina_scores, higher_is_better=False)
    p_gnina_energy = _percentiles_from_scores(ligands, gnina_energy_scores, higher_is_better=False)
    p_ledock = _percentiles_from_scores(ligands, ledock_scores, higher_is_better=False)
    p_cnn = _percentiles_from_scores(ligands, cnn_scores, higher_is_better=True)

    has_cnn = any(v > 0.0 for v in p_cnn.values())
    has_ledock = any(v > 0.0 for v in p_ledock.values())

    if any(v > 0.0 for v in p_gnina_energy.values()):
        p_energy = {lig: 0.6 * p_vina.get(lig, 0.0) + 0.4 * p_gnina_energy.get(lig, 0.0) for lig in ligands}
    else:
        p_energy = {lig: p_vina.get(lig, 0.0) for lig in ligands}

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
    for lig in ligands:
        consensus_score = (
            w_energy * p_energy.get(lig, 0.0)
            + w_ledock * p_ledock.get(lig, 0.0)
            + w_cnn * p_cnn.get(lig, 0.0)
        )
        n_engines = sum(
            1
            for store in (vina_scores, gnina_energy_scores, ledock_scores, cnn_scores)
            if lig in store
        )
        out_rows.append(
            {
                "run_id": run_id,
                "pdb_id": paths.pdb_id,
                "variant": variant_label,
                "ph_label": ph_label or "",
                "ligand": lig,
                "consensus_score": consensus_score,
                "p_energy": p_energy.get(lig, 0.0),
                "p_ledock": p_ledock.get(lig, 0.0),
                "p_cnn": p_cnn.get(lig, 0.0),
                "p_vina": p_vina.get(lig, 0.0),
                "p_gnina_energy": p_gnina_energy.get(lig, 0.0),
                "n_engines_with_data": n_engines,
            }
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
