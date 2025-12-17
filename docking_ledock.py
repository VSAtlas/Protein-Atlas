from __future__ import annotations

import json
import logging
import os
import re
import subprocess
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple

from input_and_export_functions import _to_bool
from path_router import ph_ensemble_dir
from prep_ligands_mol2_for_ledock import map_pdbqt_to_mol2_path

LEDOCK_STAGE_PARAMS = {
    "stage1": {"rmsd": 1.5, "n_poses": 10},
    "stage2": {"rmsd": 1.0, "n_poses": 20},
    "stage3": {"rmsd": 0.5, "n_poses": 40},
}

_CLUSTER_RE = re.compile(
    r"REMARK\s+Cluster\s+\d+\s+of\s+Poses:\s*(\d+)\s+Score:\s*([-0-9.]+)",
    re.IGNORECASE,
)


def should_run_ledock_for_target(cfg: Dict[str, Any]) -> bool:
    """
    Return True if LeDock should run for this target, based on USE_LEDOCK/use_ledock flags
    and PH_ENSEMBLE enablement.
    """
    raw_flag = None
    if isinstance(cfg, dict):
        file_cfg = cfg.get("_FILE_CFG") if isinstance(cfg, dict) else None
        raw_flag = cfg.get("USE_LEDOCK", cfg.get("use_ledock", None))
        if raw_flag is None and isinstance(file_cfg, dict):
            raw_flag = file_cfg.get("USE_LEDOCK", file_cfg.get("use_ledock", None))

    if raw_flag is None:
        return False
    if not _to_bool(raw_flag):
        logging.getLogger(__name__).info(
            "[ledock.disabled] USE_LEDOCK=%r -> LeDock docking is globally disabled",
            raw_flag,
        )
        return False

    if not bool(cfg.get("PH_ENSEMBLE")):
        return False
    canonical = cfg.get("_PH_ENSEMBLE_CANONICAL")
    if isinstance(canonical, dict) and not canonical:
        return False
    return True


def _normalize_ph_label(ph_label: Optional[str]) -> Optional[str]:
    if ph_label is None:
        return None
    token = str(ph_label).strip()
    return token or None


def _variant_for_ph(variant: Optional[str], legacy_mode: bool) -> Optional[str]:
    v = (str(variant).strip().upper() or None) if variant is not None else None
    if v:
        return v
    if legacy_mode:
        return None
    # Default to HOLO to match pH-ensemble receptor layout when variant is unset.
    return "HOLO"


def _label_from_manifest_entry(entry: Dict[str, Any], pdb_id: str) -> Optional[str]:
    if not isinstance(entry, dict):
        return None
    for key in ("label", "ph_label"):
        if entry.get(key):
            return str(entry[key])

    tag = entry.get("tag")
    if tag:
        tag_str = str(tag)
        prefix = f"{pdb_id}_"
        return tag_str[len(prefix):] if tag_str.startswith(prefix) else tag_str

    for key in ("withH", "pdbqt", "receptor_pdbqt", "output_pdbqt", "path", "receptor"):
        raw = entry.get(key)
        if not raw:
            continue
        stem = Path(str(raw)).stem
        if stem.endswith(".withH"):
            stem = stem[: -len(".withH")]
        prefix = f"{pdb_id}_"
        return stem[len(prefix):] if stem.startswith(prefix) else stem
    return None


def _resolve_ledock_receptor(
    cfg: Dict[str, Any],
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    logger: logging.Logger,
) -> Optional[Path]:
    label = _normalize_ph_label(ph_label)
    if not label:
        return None

    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    variant_token = _variant_for_ph(variant, legacy_mode)
    ensemble_dir = ph_ensemble_dir(pdb_id, variant=variant_token, legacy=legacy_mode)
    manifest_path = ensemble_dir / "ensemble.json"

    if manifest_path.exists():
        try:
            payload = json.loads(manifest_path.read_text(encoding="utf-8", errors="ignore"))
        except Exception as exc:
            logger.warning(
                "[ledock.manifest.error] path=%s reason=%s",
                str(manifest_path),
                exc,
            )
        else:
            members = payload.get("members") if isinstance(payload, dict) else None
            if isinstance(members, list):
                for entry in members:
                    entry_label = _label_from_manifest_entry(entry, str(pdb_id).upper())
                    if not entry_label:
                        continue
                    if entry_label.strip().lower() != label.lower():
                        continue
                    withh = entry.get("withH") or entry.get("withH_pdb") or entry.get("withH_path")
                    if withh:
                        path = Path(str(withh))
                        if not path.is_absolute():
                            path = (manifest_path.parent / path)
                        return path

    prefix = f"{str(pdb_id).upper()}_"
    label_token = label
    if not label_token.startswith(prefix):
        label_token = f"{prefix}{label_token}"
    # Fallback to naming convention when manifest is missing or incomplete.
    return ensemble_dir / f"{label_token}.withH.pdb"


def _write_text_atomic(path: Path, text: str) -> None:
    tmp = path.with_suffix(path.suffix + ".part")
    with open(tmp, "w", encoding="utf-8") as f:
        f.write(text)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def emit_ledock_config(
    receptor_pdb: Path,
    ligands_list_path: Path,
    center: Tuple[float, float, float],
    box_size: Tuple[float, float, float],
    rmsd: float,
    n_poses: int,
    logger: logging.Logger,
) -> Path:
    """Write a dock.in-style config file and return its path."""
    xmin = center[0] - box_size[0] / 2.0
    xmax = center[0] + box_size[0] / 2.0
    ymin = center[1] - box_size[1] / 2.0
    ymax = center[1] + box_size[1] / 2.0
    zmin = center[2] - box_size[2] / 2.0
    zmax = center[2] + box_size[2] / 2.0

    config_path = ligands_list_path.with_name("dock.in")
    lines = [
        "//",
        "Receptor",
        str(receptor_pdb),
        "",
        "RMSD",
        f"{float(rmsd):.3f}",
        "",
        "Binding pocket",
        f"{xmin:.3f} {xmax:.3f}",
        f"{ymin:.3f} {ymax:.3f}",
        f"{zmin:.3f} {zmax:.3f}",
        "",
        "Number of binding poses",
        str(int(n_poses)),
        "",
        "Ligands list",
        str(ligands_list_path),
        "//",
    ]
    _write_text_atomic(config_path, "\n".join(lines) + "\n")

    logger.info(
        "[ledock.cfg] receptor=%s ligands=%s rmsd=%.3f n_poses=%d dock_in=%s",
        str(receptor_pdb),
        str(ligands_list_path),
        float(rmsd),
        int(n_poses),
        str(config_path),
    )
    return config_path


def _parse_ledock_dok(dok_path: Path) -> Dict[str, Any]:
    if not dok_path.exists():
        return {
            "best_score_kcal": None,
            "n_poses": 0,
            "cluster_count": 0,
            "valid": False,
            "reason": "ledock_missing_dok",
        }

    try:
        lines = dok_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return {
            "best_score_kcal": None,
            "n_poses": 0,
            "cluster_count": 0,
            "valid": False,
            "reason": "ledock_parse_error",
        }

    has_atom = False
    cluster_scores: list[float] = []
    cluster_count = 0
    n_poses = 0

    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            has_atom = True
        m = _CLUSTER_RE.search(line)
        if not m:
            continue
        cluster_count += 1
        try:
            poses = int(m.group(1))
            if poses > n_poses:
                n_poses = poses
        except Exception:
            pass
        try:
            score = float(m.group(2))
            cluster_scores.append(score)
        except Exception:
            pass

    if not has_atom or cluster_count == 0 or not cluster_scores:
        return {
            "best_score_kcal": None,
            "n_poses": 0,
            "cluster_count": 0,
            "valid": False,
            "reason": "ledock_parse_error",
        }

    best_score = min(cluster_scores)
    return {
        "best_score_kcal": float(best_score),
        "n_poses": int(n_poses),
        "cluster_count": int(cluster_count),
        "valid": True,
        "reason": None,
    }


def _resolve_stage_params(stage_name: str) -> Tuple[Optional[str], Optional[Dict[str, Any]]]:
    if stage_name in LEDOCK_STAGE_PARAMS:
        return stage_name, dict(LEDOCK_STAGE_PARAMS[stage_name])
    for key in LEDOCK_STAGE_PARAMS:
        if str(stage_name).endswith(key):
            return key, dict(LEDOCK_STAGE_PARAMS[key])
    return None, None


def run_ledock_for_stage(
    cfg: Dict[str, Any],
    paths: Any,
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    stage_name: str,
    stage_info: Dict[str, Any],
    ligands: Iterable[Path],
    center: Optional[Tuple[float, float, float]],
    box_size: Optional[Tuple[float, float, float]],
    logger: logging.Logger,
    receptor_pdb: Optional[Path] = None,
) -> Tuple[Dict[Path, float], Dict[Path, Dict[str, Any]]]:
    stage_info = stage_info or {}
    scores: Dict[Path, float] = {}
    ledock_metrics: Dict[Path, Dict[str, Any]] = {}

    stage_key, stage_params = _resolve_stage_params(stage_name)
    if not stage_params:
        logger.warning("[ledock.skip] reason=unknown_stage stage=%s", stage_name)
        return scores, ledock_metrics

    if center is None or box_size is None:
        logger.info("[ledock.skip] reason=missing_center_or_box")
        return scores, ledock_metrics

    if not bool(cfg.get("PH_ENSEMBLE")):
        logger.info("[ledock.skip] reason=ph_ensemble_disabled_or_missing_receptor")
        return scores, ledock_metrics

    if not should_run_ledock_for_target(cfg):
        logger.info("[ledock.skip] reason=use_ledock_disabled")
        return scores, ledock_metrics

    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    variant_token = (str(variant).strip().upper() or None) if variant is not None else None
    variant_for_ph = _variant_for_ph(variant_token, legacy_mode)
    ph_token = _normalize_ph_label(ph_label)

    receptor_path = Path(receptor_pdb) if receptor_pdb else None
    if receptor_path is None:
        receptor_path = _resolve_ledock_receptor(cfg, pdb_id, variant_token, ph_label, logger)

    if receptor_path:
        logger.info(
            "[ledock.receptor] pdb=%s variant=%s ph=%s receptor_pdb=%s",
            pdb_id,
            variant_for_ph or "HOLO",
            ph_token or "ph_ensemble",
            str(receptor_path),
        )
    if receptor_path is None or not receptor_path.exists():
        logger.info("[ledock.skip] reason=ph_ensemble_disabled_or_missing_receptor")
        return scores, ledock_metrics

    ligand_paths = [Path(lig) for lig in ligands]
    if not ligand_paths:
        logger.info("[ledock.skip] reason=no_ligands")
        return scores, ledock_metrics

    mol2_map: Dict[Path, Path] = {}
    for lig_path in ligand_paths:
        mol2_path = map_pdbqt_to_mol2_path(lig_path, logger=logger)
        if mol2_path is None or not mol2_path.exists():
            ledock_metrics[lig_path] = {
                "best_score_kcal": None,
                "n_poses": 0,
                "cluster_count": 0,
                "valid": False,
                "reason": "ledock_missing_mol2",
            }
            continue
        mol2_map[lig_path] = mol2_path

    if not mol2_map:
        logger.info("[ledock.skip] reason=no_mol2_ligands")
        return scores, ledock_metrics

    stage_params["rmsd"] = float(stage_info.get("rmsd", stage_params["rmsd"]))
    stage_params["n_poses"] = int(stage_info.get("n_poses", stage_params["n_poses"]))

    dock_root = paths.docked_variant_root(variant_for_ph, ph_token)
    ledock_root = dock_root / "ledock"
    stage_root = ledock_root / stage_name
    stage_root.mkdir(parents=True, exist_ok=True)

    mol2_paths = [mol2_map[lig] for lig in ligand_paths if lig in mol2_map]
    ligands_list_path = stage_root / f"ligands_{stage_key}.list"
    _write_text_atomic(
        ligands_list_path,
        "\n".join(str(p) for p in mol2_paths) + "\n",
    )

    config_path = emit_ledock_config(
        receptor_pdb=receptor_path,
        ligands_list_path=ligands_list_path,
        center=center,
        box_size=box_size,
        rmsd=stage_params["rmsd"],
        n_poses=stage_params["n_poses"],
        logger=logger,
    )

    logger.info(
        "[ledock.stage] pdb=%s stage=%s variant=%s ph=%s n_lig=%d",
        pdb_id,
        stage_name,
        variant_for_ph or "HOLO",
        ph_token or "ph_ensemble",
        len(mol2_paths),
    )

    cmd = ["ledock", str(config_path)]
    logger.info("[ledock.cmd] cmd=%s cwd=%s", " ".join(cmd), str(stage_root))
    try:
        subprocess.run(cmd, check=True, cwd=stage_root)
    except Exception as exc:
        logger.warning(
            "[ledock.error] pdb=%s stage=%s reason=%s",
            pdb_id,
            stage_name,
            exc,
        )
        for lig in mol2_map:
            if lig in ledock_metrics:
                continue
            ledock_metrics[lig] = {
                "best_score_kcal": None,
                "n_poses": 0,
                "cluster_count": 0,
                "valid": False,
                "reason": "ledock_run_failed",
            }
        return scores, ledock_metrics

    for lig_path, mol2_path in mol2_map.items():
        dok_path = stage_root / f"{mol2_path.stem}.dok"
        metrics = _parse_ledock_dok(dok_path)
        ledock_metrics[lig_path] = metrics
        if metrics.get("valid"):
            score = metrics.get("best_score_kcal")
            if isinstance(score, (int, float)):
                scores[lig_path] = float(score)
            logger.info(
                "[ledock.score] ligand=%s score=%s rmsd=%.3f",
                mol2_path.name,
                score if score is not None else "None",
                float(stage_params["rmsd"]),
            )
        else:
            logger.warning(
                "[ledock.parse] ligand=%s reason=%s dok=%s",
                mol2_path.name,
                metrics.get("reason"),
                str(dok_path),
            )

    return scores, ledock_metrics
