from __future__ import annotations

import logging
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Dict, Iterator, List, Optional, Tuple

from path_router.path_router import make_paths
from docking.pose_validation import validate_pose_pdbqt, attempt_fallback_recenter


class BudgetGuard:
    """Simple per-ligand wall-clock guard for retries/validation."""

    def __init__(
        self,
        max_seconds: float | None = None,
        *,
        dock_max_seconds: float | None = None,
        idle_max_seconds: float | None = None,
        log: logging.Logger | None = None,
        **_ignored: object,
    ) -> None:
        effective_max = max_seconds if max_seconds is not None else dock_max_seconds
        if effective_max is None:
            effective_max = 300.0

        self.max_seconds = float(effective_max)
        self.idle_max_seconds = (
            float(idle_max_seconds) if idle_max_seconds is not None else None
        )
        self.log = log
        self._deadline: float = 0.0
        self._lig_cfg: Optional[str] = None
        self._out_pdbqt: Optional[str] = None
        self.reset()

    def reset(self) -> None:
        """Reset the wall-clock deadline from now."""

        self._deadline = time.time() + self.max_seconds

    def bind(self, lig_cfg: str | None, out_pdbqt: str | None) -> None:
        """
        Optionally remember which config/output this guard is associated with.
        This is for debugging/telemetry only; no logic elsewhere depends on it.
        """

        self._lig_cfg = lig_cfg
        self._out_pdbqt = out_pdbqt

    def expired(self) -> bool:
        """Return True if the wall-clock budget has been exceeded."""

        return time.time() >= self._deadline


def _iter_pdbqt_models(pdbqt_path: str) -> Iterator[str]:
    """
    Yield individual MODEL..ENDMDL blocks from a (possibly multi-model) PDBQT.
    If no MODEL/ENDMDL markers exist, yield the whole file once.
    """
    buf: List[str] = []
    saw_model = False
    with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith("MODEL"):
                if buf:
                    yield "".join(buf)
                    buf = []
                saw_model = True
                buf.append(ln)
            elif ln.startswith("ENDMDL"):
                buf.append(ln)
                yield "".join(buf)
                buf = []
                saw_model = True
            else:
                if saw_model:
                    buf.append(ln)

    # If we never saw a MODEL block, treat the whole file as one model
    if not saw_model:
        try:
            with open(pdbqt_path, "r", encoding="utf-8", errors="ignore") as f2:
                yield f2.read()
        except Exception:
            yield ""


def validate_first_valid_pose(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    pocket_center: tuple[float, float, float],
    surface_coords,
    max_models: int = 3,
    clash_threshold: float = 2.0,
    clash_tol: int = 3,
    dist_surf: float = 6.0,
    dist_centroid: float = 4.5,
):
    """
    Validate poses in order and return as soon as one passes.
    Falls back to the last invalid result if none pass.
    """
    tmp_dir = Path(ligand_pdbqt).parent
    best_invalid = None
    count = 0

    for idx, model_text in enumerate(_iter_pdbqt_models(ligand_pdbqt)):
        if max_models and count >= int(max_models):
            break
        count += 1

        tmp = Path(ligand_pdbqt)
        try:
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                errors="ignore",
                dir=str(tmp_dir),
                prefix=f"{Path(ligand_pdbqt).stem}.m{idx}.",
                suffix=".tmp.pdbqt",
                delete=False,
            ) as handle:
                handle.write(model_text)
                tmp = Path(handle.name)
        except Exception:
            tmp = Path(ligand_pdbqt)

        try:
            res = validate_pose_pdbqt(
                protein_pdbqt=receptor_pdbqt,
                ligand_pdbqt=str(tmp),
                pocket_center=pocket_center,
                clash_threshold=clash_threshold,
                CLASH_TOLERANCE=clash_tol,
                DIST_THRESHOLD_SURFACE=dist_surf,
                DIST_THRESHOLD_CENTROID=dist_centroid,
                surface_atom_coords=surface_coords,
            )
        finally:
            if tmp.name.endswith(".tmp.pdbqt"):
                try:
                    tmp.unlink(missing_ok=True)
                except Exception:
                    pass

        if res.get("valid", False):
            return res  # early exit on first valid

        best_invalid = res  # keep the last invalid for diagnostics

    return best_invalid or {"valid": False, "reason": "no_poses"}


@dataclass
class RecenterParams:
    """Thresholds for early/fallback recenter heuristics."""

    EARLY_RECENTER_RATIO: float = 0.70
    EARLY_RECENTER_MIN_EVAL: int = 10
    EARLY_RECENTER_FAR_A: float = 15.0
    EARLY_RECENTER_MEDIAN_A: float = 10.0
    ALLOW_BOX_EXPAND: bool = True
    MAX_RECENTER_ATTEMPTS: int = 1  # tightened: fewer early recenter tries


def load_recenter_params(cfg: Mapping[str, Any]) -> RecenterParams:
    """Read early/fallback recentering knobs from runtime config."""
    return RecenterParams(
        EARLY_RECENTER_RATIO=float(cfg.get("EARLY_RECENTER_RATIO", 0.70)),
        EARLY_RECENTER_MIN_EVAL=int(cfg.get("EARLY_RECENTER_MIN_EVAL", 10)),
        EARLY_RECENTER_FAR_A=float(cfg.get("EARLY_RECENTER_FAR_A", 15.0)),
        EARLY_RECENTER_MEDIAN_A=float(cfg.get("EARLY_RECENTER_MEDIAN_A", 10.0)),
        ALLOW_BOX_EXPAND=bool(cfg.get("ALLOW_BOX_EXPAND", True)),
        MAX_RECENTER_ATTEMPTS=int(cfg.get("MAX_RECENTER_ATTEMPTS", 1)),
    )


@dataclass
class GlobalCenterGuard:
    """
    Gatekeeper for ANY global center change (early recenter, empty-stage fallback,
    CenterSelector promotions). Supports a hard 'lock' after a validated control ligand.
    """

    max_global_switches: int = 2
    global_switches: int = 0
    switched_this_stage: bool = False
    locked: bool = False

    def reset_stage(self) -> None:
        """Reset per-stage switch flag (call at the start of each stage)."""
        self.switched_this_stage = False

    def can_switch(self) -> bool:
        """
        True if a center change is allowed right now.
        Respects: hard lock, one-per-stage, and global cap.
        """
        return (
            (not self.locked)
            and (not self.switched_this_stage)
            and (self.global_switches < self.max_global_switches)
        )

    def mark_switch(self) -> None:
        """Record that a center change just happened this stage."""
        self.global_switches += 1
        self.switched_this_stage = True

    def lock(self) -> None:
        """Hard-lock: disallow any further center changes for the remainder of the run."""
        self.locked = True


def fallback_recentering_if_empty(
    cfg: Dict,
    pdb_id: str,
    stage_name: str,
    scores: Dict[str, float],
    raw_docked_ligands: Dict[str, str],
    receptor_pdbqt: str,
    center: Tuple[float, float, float],
    box_size: Tuple[float, float, float],
    stage1_original: List[str],
    logger: logging.Logger,
    guard: GlobalCenterGuard,
    control_anchor_hit: bool,
) -> Tuple[bool, Tuple[float, float, float], Tuple[float, float, float], List[str]]:
    """
    If a stage yields no valid ligands, attempt a fallback recenter and restart stage1.
    De-duped: will not fire if early recenter or CenterSelector already switched this stage,
    or if a control anchor validated in this stage, or if global cap reached.
    """
    if scores:
        return False, center, box_size, []
    if control_anchor_hit:
        logger.info(
            "Empty-stage fallback skipped: control-anchored validation present earlier."
        )
        return False, center, box_size, []
    if not guard.can_switch():
        logger.info(
            "Empty-stage fallback skipped: global switch guard disallows further switches."
        )
        return False, center, box_size, []

    logger.warning(
        f"No valid ligands in {stage_name}. Attempting fallback recentering..."
    )
    # >>> DOCKED PATHS PATCH START
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    # >>> DOCKED PATHS PATCH END
    fb_pose, new_center, _best_score, _chosen_lig = attempt_fallback_recenter(
        fallback_ligands=raw_docked_ligands,
        receptor_pdbqt=receptor_pdbqt,
        docking_dir=str(paths.docked_pdb_root()),
        stage_name=stage_name,
        pocket_center=center,
        logger=logger,
        exclude_basenames=set(),
    )
    if new_center is None:
        logger.warning("Fallback recovery failed.")
        return False, center, box_size, []

    new_center_tuple = (
        float(new_center[0]),
        float(new_center[1]),
        float(new_center[2]),
    )
    box_cap = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
    new_box = (
        min(box_cap, float(box_size[0])),
        min(box_cap, float(box_size[1])),
        min(box_cap, float(box_size[2])),
    )
    guard.mark_switch()  # counts as a global switch
    logger.info(
        "Re-running stage1 with new center after no-valid fallback. [global switch]"
    )
    return True, new_center_tuple, new_box, stage1_original[:]
