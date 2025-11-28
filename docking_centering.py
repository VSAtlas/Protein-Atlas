from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from fallback_recenter import GlobalCenterGuard


@dataclass
class CenterDecision:
    new_center: Optional[Tuple[float, float, float]]
    reason: str = ""
    switchscore: float = 0.0
    promoted: bool = False


class CenterSelector:
    """
    Decides when to keep the current center (control-anchored) vs. switch to a newly discovered pocket.
    Uses clustering of valid pose centroids + SwitchScore, now with stronger control anchoring and guard.
    """

    def __init__(self, cfg: Dict, logger: logging.Logger,
                 control_stems: set,
                 heavy_atom_counts: Dict[str, int],
                 initial_center: Tuple[float, float, float]):
        self.cfg = cfg
        self.log = logger
        self.control_stems = {s.lower() for s in control_stems}
        self.heavy = heavy_atom_counts
        self.eps = float(cfg.get("CLUSTER_EPS_ANG", 3.5))
        self.mode = cfg.get("CENTER_MODE", "control-first").lower()
        self.ctrl_blacklist = {s.strip().upper() for s in cfg.get("CONTROL_BLACKLIST", "").split(",") if s.strip()}
        self.ctrl_min_heavy = int(cfg.get("CONTROL_MIN_HEAVY_ATOMS", 10))
        self.allow_switch_from_control = bool(cfg.get("ALLOW_SWITCH_FROM_CONTROL", True))
        self.require_control_failure = bool(cfg.get("REQUIRE_CONTROL_FAILURE_TO_SWITCH", False))
        self.threshold = float(cfg.get("SWITCH_SCORE_THRESHOLD", 0.7))
        self.away_from_control_boost = float(cfg.get("SWITCH_AWAY_FROM_CONTROL_MIN_BOOST", 2.5))
        self.hysteresis = float(cfg.get("SWITCH_SCORE_HYSTERESIS", 0.5))
        self.switch_history = []  # keep last few SwitchScores
        self.current_center = np.array(initial_center, float)

        self.curr_stats = {"median_score": None, "median_le": None, "valid_rate": 0.0}
        self.last_decision_had_control_anchor = False

    @staticmethod
    def _strip_stage(name: str) -> str:
        stem = Path(name).stem
        return stem.split("_stage")[0].lower()

    @staticmethod
    def _pdbqt_centroid(pdbqt_path: str) -> Optional[np.ndarray]:
        if not pdbqt_path or not os.path.exists(pdbqt_path):
            return None
        xs, ys, zs = [], [], []
        try:
            with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if line.startswith(("ATOM", "HETATM")):
                        try:
                            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                        except Exception:
                            parts = line.split()
                            if len(parts) >= 8:
                                x = float(parts[5]); y = float(parts[6]); z = float(parts[7])
                            else:
                                continue
                        xs.append(x); ys.append(y); zs.append(z)
            if xs:
                return np.array([np.mean(xs), np.mean(ys), np.mean(zs)], dtype=float)
        except Exception:
            return None
        return None

    def _is_blacklisted_control_name(self, stem_upper: str) -> bool:
        return any(stem_upper.startswith(bad) for bad in self.ctrl_blacklist)

    def _is_control(self, lig_path: str) -> bool:
        stem = self._strip_stage(os.path.basename(lig_path))
        if stem.upper() and self._is_blacklisted_control_name(stem.upper()):
            return False
        ha = self.heavy.get(lig_path)
        if ha is not None and ha < self.ctrl_min_heavy:
            return False
        return (stem in self.control_stems)

    def _cluster(self, points: List[np.ndarray]) -> List[List[int]]:
        clusters: List[List[int]] = []
        for i, p in enumerate(points):
            placed = False
            for cl in clusters:
                c = np.mean([points[j] for j in cl], axis=0)
                if np.linalg.norm(p - c) <= self.eps:
                    cl.append(i); placed = True; break
            if not placed:
                clusters.append([i])
        return clusters

    def _median(self, arr):
        return float(np.median(arr)) if arr else None

    def _cluster_stats(self,
                       member_idxs: List[int],
                       ligs: List[str],
                       scores_map: Dict[str, float],
                       centroids: List[np.ndarray],
                       total_docked: int) -> dict:
        members = [ligs[i] for i in member_idxs]
        scores = [scores_map[m] for m in members if m in scores_map]
        les = []
        for m in members:
            s = scores_map.get(m)
            if s is None:
                continue
            ha = self.heavy.get(m)
            if ha and ha > 0:
                les.append((-s) / ha)  # s is negative kcal/mol
        center = np.mean([centroids[i] for i in member_idxs], axis=0)
        valid_rate = len(members) / max(1, total_docked)
        ctrl_hits = sum(1 for m in members if self._is_control(m))
        return {
            "center": tuple(center.tolist()),
            "n": len(members),
            "valid_rate": float(valid_rate),
            "median_score": self._median(scores),
            "median_le": self._median(les),
            "control_hits": ctrl_hits,
        }

    def _switch_score(self, valid_rate, score_boost_kcal, le_gain, pocketability=0.5) -> float:
        consensus_boost = max(0.0, min(1.0, score_boost_kcal / 3.0))
        le_norm = max(0.0, min(1.0, le_gain / 0.05))
        return 0.40*valid_rate + 0.20*consensus_boost + 0.20*pocketability + 0.20*le_norm

    def consider_switch(self,
                        stage_name: str,
                        scores: Dict[str, float],
                        validated_ligands: List[str],
                        raw_docked: Dict[str, str],
                        receptor_pdbqt: str,
                        current_center: Tuple[float, float, float],
                        guard: GlobalCenterGuard) -> CenterDecision:
        # If the guard forbids a switch this stage or we're capped out, exit early
        if not guard.can_switch():
            return CenterDecision(None, "guard_disallowed")

        if not validated_ligands:
            self.last_decision_had_control_anchor = False
            return CenterDecision(None, "no_validated_poses")

        ligs = sorted(set(validated_ligands))
        centroids = []
        ligs_kept = []
        for lig in ligs:
            pose = raw_docked.get(lig)
            c = self._pdbqt_centroid(pose) if pose else None
            if c is not None:
                centroids.append(c); ligs_kept.append(lig)

        if len(centroids) < 3:
            self.last_decision_had_control_anchor = False
            return CenterDecision(None, "too_few_centroids")

        clusters = self._cluster(centroids)
        total_docked = max(1, len(raw_docked))
        stats = [self._cluster_stats(cl, ligs_kept, scores, centroids, total_docked) for cl in clusters]

        curr_idx = None
        curr_center = np.array(current_center, float)
        for idx, st in enumerate(stats):
            if np.linalg.norm(np.array(st["center"]) - curr_center) <= self.eps:
                curr_idx = idx; break

        curr_median = stats[curr_idx]["median_score"] if curr_idx is not None else None
        curr_le = stats[curr_idx]["median_le"] if curr_idx is not None else None
        curr_valid = stats[curr_idx]["valid_rate"] if curr_idx is not None else 0.0

        curr_ctrl_hits = stats[curr_idx]["control_hits"] if curr_idx is not None else 0
        self.curr_stats = {"median_score": curr_median, "median_le": curr_le, "valid_rate": curr_valid}

        # Strong control anchoring (compute for THIS stage before using it):
        self.last_decision_had_control_anchor = (
                curr_ctrl_hits > 0 and curr_valid >= float(self.cfg.get("CONTROL_ANCHOR_MIN_VALID_RATE", 0.10))
        )

        # If policy requires control failure to switch, stop here while control anchors
        if self.require_control_failure and self.last_decision_had_control_anchor:
            return CenterDecision(None, "control_anchor_lock")

        # Choose best alternative cluster
        best_cand = None
        best_score = -1.0
        for idx, st in enumerate(stats):
            if idx == curr_idx:
                continue
            score_boost = 0.0
            le_gain = 0.0
            if curr_median is not None and st["median_score"] is not None:
                score_boost = (curr_median - st["median_score"])
            if curr_le is not None and st["median_le"] is not None:
                le_gain = st["median_le"] - curr_le

            if st["valid_rate"] < float(self.cfg.get("SWITCH_VALID_RATE_MIN", 0.40)):
                continue
            # Require bigger improvement if leaving a control-anchored cluster
            required_boost = float(self.cfg.get("SWITCH_SCORE_IMPROVE_MIN", 1.5))
            if self.last_decision_had_control_anchor:
                required_boost = max(required_boost, self.away_from_control_boost)

            if score_boost < required_boost:
                continue
            if le_gain < float(self.cfg.get("SWITCH_LE_GAIN_MIN", 0.02)):
                continue

            sw = self._switch_score(st["valid_rate"], score_boost, le_gain, pocketability=0.5)
            if sw > best_score:
                best_score = sw
                best_cand = st

        if not best_cand:
            return CenterDecision(None, "no_candidate_passed_gates")

        self.switch_history.append(best_score)
        promoted = best_score >= self.threshold
        if not promoted:
            return CenterDecision(None, f"below_threshold:{best_score:.2f}", best_score, False)

        # mark guard switch outside (caller), but flag that we want to promote
        return CenterDecision(new_center=tuple(best_cand["center"]),
                              reason=f"promote_new_center score={best_score:.2f}",
                              switchscore=best_score,
                              promoted=True)


__all__ = ["CenterDecision", "CenterSelector"]
