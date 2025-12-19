from __future__ import annotations

import json
import logging
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Tuple

import pandas as pd

from dud_eval import compute_decoy_stats_from_long_csv, guess_ligfile_col, guess_score_col
from input_and_export_functions import _to_bool, write_score_summary_to_csv
from path_router import make_paths, ph_ensemble_dir
from prep_for_ledock import ensure_ledock_receptor, map_pdbqt_to_mol2_path, symlink_mol2_for_stage

LEDOCK_STAGE_PARAMS = {
    "stage1": {"rmsd": 1.5, "n_poses": 10},
    "stage2": {"rmsd": 1.0, "n_poses": 20},
    "stage3": {"rmsd": 0.5, "n_poses": 40},
}

_CLUSTER_RE = re.compile(
    r"REMARK\s+Cluster\s+\d+\s+of\s+Poses:\s*(\d+)\s+Score:\s*([-0-9.]+)",
    re.IGNORECASE,
)


def annotate_ledock_fda_long_csv_with_t_scores_vs_decoys(
    cfg: Dict[str, Any],
    pdb_id: str,
    ph_label: Optional[str] = None,
    logger=None,
) -> Optional[str]:
    """
    Annotate LeDock FDA long CSV with T-scores vs decoys, mirroring Vina/GNINA helpers.
    """
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    var = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
    ph_token = (ph_label or "").strip() or None
    variant_root = Path(paths.docked_variant_root(var, ph_token))

    dud_csv = variant_root / "dud_ledock_docking_score_long.csv"
    fda_csv = variant_root / "ledock_docking_score_long.csv"

    if not dud_csv.exists() or not fda_csv.exists():
        if logger:
            logger.info(
                "[ledock.t-score.skip] pdb_id=%s ph=%s reason=missing_csv dud=%s fda=%s",
                pdb_id,
                ph_label or "base",
                str(dud_csv),
                str(fda_csv),
            )
        return None

    mu, sigma, n_decoys = compute_decoy_stats_from_long_csv(dud_csv)
    if (not n_decoys) or (not math.isfinite(mu)) or (not math.isfinite(sigma)) or sigma == 0.0:
        if logger:
            logger.info(
                "[ledock.t-score.skip] pdb_id=%s ph=%s reason=degenerate_stats n=%s mu=%s sigma=%s",
                pdb_id,
                ph_label or "base",
                n_decoys,
                mu,
                sigma,
            )
        return None

    df = pd.read_csv(fda_csv)
    if df.empty:
        return None

    lig_col = guess_ligfile_col(df, None)
    score_col = guess_score_col(df, None)
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")

    best = df.groupby(lig_col, as_index=False).agg(best_score=(score_col, "min"))
    best["ledock_t_vs_decoys"] = (mu - best["best_score"]) / sigma
    t_map = dict(zip(best[lig_col], best["ledock_t_vs_decoys"]))
    df["ledock_t_vs_decoys"] = df[lig_col].map(t_map)

    df.to_csv(fda_csv, index=False)
    if logger:
        logger.info(
            "[ledock.t-score.ok] pdb_id=%s ph=%s n_decoys=%s mean=%.3f std=%.3f out=%s",
            pdb_id,
            ph_label or "base",
            n_decoys,
            mu,
            sigma,
            str(fda_csv),
        )
    return str(fda_csv)


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
    variant_token = (str(variant).strip().upper() or None) if variant is not None else None
    variant_for_ph = _variant_for_ph(variant_token, legacy_mode)
    ensemble_dir = ph_ensemble_dir(pdb_id, variant=variant_for_ph, legacy=legacy_mode)
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
            logger.info(
                "[ledock.manifest.miss] pdb=%s variant=%s ph=%s path=%s",
                pdb_id,
                variant_for_ph or "HOLO",
                label,
                str(manifest_path),
            )
    else:
        logger.info(
            "[ledock.manifest.miss] pdb=%s variant=%s ph=%s path=%s",
            pdb_id,
            variant_for_ph or "HOLO",
            label,
            str(manifest_path),
        )

    prefix = f"{str(pdb_id).upper()}_"
    label_token = label
    if not label_token.startswith(prefix):
        label_token = f"{prefix}{label_token}"
    fallback = ensemble_dir / f"{label_token}.withH.pdb"
    if not fallback.exists():
        logger.info(
            "[ledock.receptor.missing] pdb=%s variant=%s ph=%s path=%s",
            pdb_id,
            variant_for_ph or "HOLO",
            label,
            str(fallback),
        )
        return None
    return fallback


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


def resolve_ledock_stage_params(
    cfg: Mapping[str, Any],
    stage_name: str,
) -> Tuple[Optional[str], Optional[Dict[str, Any]]]:
    """
    Wrapper around _resolve_stage_params that also applies FAST_MODE overrides.
    """
    stage_key, params = _resolve_stage_params(stage_name)
    if not params:
        return stage_key, params
    if cfg.get("FAST_MODE"):
        fast_rmsd = float(cfg.get("LEDOCK_FAST_RMSD", 1.5))
        fast_n_poses = int(cfg.get("LEDOCK_FAST_N_POSES", 1))
        params["rmsd"] = fast_rmsd
        params["n_poses"] = fast_n_poses
    return stage_key, params


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

    stage_key, stage_params = resolve_ledock_stage_params(cfg, stage_name)
    if not stage_params:
        logger.warning("[ledock.skip] reason=unknown_stage stage=%s", stage_name)
        return scores, ledock_metrics

    if center is None or box_size is None:
        logger.info("[ledock.skip] reason=missing_center_or_box")
        return scores, ledock_metrics

    if not should_run_ledock_for_target(cfg):
        return scores, ledock_metrics

    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    variant_token = (str(variant).strip().upper() or None) if variant is not None else None
    variant_for_ph = _variant_for_ph(variant_token, legacy_mode)
    ph_token = _normalize_ph_label(ph_label)

    receptor_path = Path(receptor_pdb) if receptor_pdb else None
    if receptor_path is None:
        receptor_path = ensure_ledock_receptor(cfg, pdb_id, variant_token, ph_label, logger)

    if receptor_path:
        logger.info(
            "[ledock.receptor] pdb=%s variant=%s ph=%s receptor_pdb=%s",
            pdb_id,
            variant_for_ph or "HOLO",
            ph_token or "ph_ensemble",
            str(receptor_path),
        )
    if receptor_path is None or not receptor_path.exists():
        logger.info("[ledock.skip] reason=missing_receptor")
        return scores, ledock_metrics

    ligand_paths = [Path(lig) for lig in ligands]
    if not ligand_paths:
        logger.info("[ledock.skip] reason=no_ligands")
        return scores, ledock_metrics

    ligand_pairs: list[Tuple[Path, Path]] = []
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
        ligand_pairs.append((lig_path, mol2_path))

    if not ligand_pairs:
        logger.info("[ledock.skip] reason=no_mol2_ligands")
        return scores, ledock_metrics

    stage_params["rmsd"] = float(stage_info.get("rmsd", stage_params["rmsd"]))
    stage_params["n_poses"] = int(stage_info.get("n_poses", stage_params["n_poses"]))

    if cfg.get("FAST_MODE"):
        fast_rmsd = float(cfg.get("LEDOCK_FAST_RMSD", 1.5))
        fast_n_poses = int(cfg.get("LEDOCK_FAST_N_POSES", 1))
        stage_params["rmsd"] = fast_rmsd
        stage_params["n_poses"] = fast_n_poses
        logger.info(
            "[ledock.fast] enabled=True pdb=%s stage=%s rmsd=%.3f n_poses=%d",
            pdb_id,
            stage_name,
            fast_rmsd,
            fast_n_poses,
        )

    ensemble_dir = ph_ensemble_dir(pdb_id, variant=variant_for_ph, legacy=legacy_mode)
    ledock_root = ensemble_dir / "ledock"
    stage_root = ledock_root / stage_name
    stage_root.mkdir(parents=True, exist_ok=True)
    for stale in stage_root.glob("*.dok"):
        try:
            stale.unlink()
        except FileNotFoundError:
            pass
    for stale in stage_root.glob("*.mol2"):
        try:
            stale.unlink()
        except FileNotFoundError:
            pass

    alias_map = symlink_mol2_for_stage(ligand_pairs, stage_root, logger)
    if not alias_map:
        logger.info("[ledock.skip] reason=no_symlinked_ligands")
        return scores, ledock_metrics

    ligands_list_path = stage_root / f"ligands_{stage_key}.list"
    ligand_lines = [alias.name for alias in alias_map.values()]
    _write_text_atomic(
        ligands_list_path,
        "\n".join(ligand_lines) + "\n",
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
        len(alias_map),
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
        for lig_path, _mol2_path in ligand_pairs:
            if lig_path in ledock_metrics:
                continue
            ledock_metrics[lig_path] = {
                "best_score_kcal": None,
                "n_poses": 0,
                "cluster_count": 0,
                "valid": False,
                "reason": "ledock_run_failed",
            }
        return scores, ledock_metrics

    for lig_path, alias_path in alias_map.items():
        dok_path = stage_root / f"{alias_path.stem}.dok"
        metrics = _parse_ledock_dok(dok_path)
        ledock_metrics[lig_path] = metrics
        if metrics.get("valid"):
            score = metrics.get("best_score_kcal")
            if isinstance(score, (int, float)):
                scores[lig_path] = float(score)
            logger.info(
                "[ledock.score] stage=%s ligand=%s score=%s",
                stage_name,
                alias_path.name,
                score if score is not None else "None",
            )
        else:
            logger.warning(
                "[ledock.parse] stage=%s ligand=%s reason=%s dok=%s",
                stage_name,
                alias_path.name,
                metrics.get("reason"),
                str(dok_path),
            )

    dock_root = paths.docked_variant_root(variant_for_ph, ph_token)
    if stage_key == "stage1":
        dok_dest = dock_root / "ledock_stage1"
    elif stage_key == "stage2":
        dok_dest = dock_root / "ledock_stage2"
    elif stage_key == "stage3":
        dok_dest = dock_root / "ledock_stage3"
    else:
        dok_dest = dock_root / f"ledock_{stage_name}"
    dok_dest.mkdir(parents=True, exist_ok=True)
    for stale in dok_dest.glob("*.dok"):
        try:
            stale.unlink()
        except FileNotFoundError:
            pass

    moved_any = False
    for alias_path in alias_map.values():
        src = stage_root / f"{alias_path.stem}.dok"
        if not src.exists():
            continue
        dest = dok_dest / src.name
        try:
            if dest.exists() or dest.is_symlink():
                dest.unlink()
        except FileNotFoundError:
            pass
        shutil.move(str(src), dest)
        moved_any = True

    if moved_any:
        logger.info(
            "[ledock.move_dok] stage=%s source=%s dest=%s",
            stage_name,
            str(stage_root),
            str(dok_dest),
        )

    return scores, ledock_metrics


def write_ledock_scores_csv(
    cfg: Dict,
    pdb_id: str,
    ledock_metrics_by_stage: Dict[str, Dict[str, Dict[str, Any]]],
    *,
    ph_label: Optional[str] = None,
    variant: Optional[str] = None,
    csv_prefix: str = "",
) -> str:
    """
    Emit LeDock score CSVs:

    - Wide summary: ledock_docking_score_summary.csv
        One row per ligand, one column per LeDock stage.
        Cell values are best LeDock energies (kcal/mol), formatted to 2 decimals.

    - Long format: <prefix>ledock_docking_score_long.csv
        One row per ligand per LeDock stage, with:
        run_id, optional variant, stage, ligand,
        ledock_best_score_kcal, ledock_cluster_count, ledock_n_poses.

    NOTE: For now, we ignore valid/invalid flags; no 'valid' or 'reason' columns.
    """
    import csv as _csv
    import logging
    import os
    from pathlib import Path

    logger = logging.getLogger(__name__)

    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_token = (ph_label or "").strip() or None

    variant_env = (variant or os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    variant_token = variant_env or None

    dock_dir = Path(paths.docked_variant_root(variant_token, ph_token))
    dock_dir.mkdir(parents=True, exist_ok=True)

    run_id_value = str(cfg.get("RUN_ID") or "")
    variant_value = variant_env
    include_variant = bool(variant_value)

    summary_name = "ledock_docking_score_summary.csv"
    long_name = f"{csv_prefix}ledock_docking_score_long.csv"

    csv_out_wide = str(dock_dir / summary_name)
    csv_out_long = str(dock_dir / long_name)

    flat: Dict[str, Dict[str, str]] = {}

    for stage_name, stage_map in ledock_metrics_by_stage.items():
        stage_dict: Dict[str, str] = {}
        for lig, rec in stage_map.items():
            lig_key = os.path.basename(str(lig))
            best = rec.get("best_score_kcal")
            if isinstance(best, (int, float)):
                stage_dict[lig_key] = f"{best:.2f}"
            else:
                stage_dict[lig_key] = ""
        flat[stage_name] = stage_dict

    write_score_summary_to_csv(
        flat,
        output_path=csv_out_wide,
        run_id=run_id_value,
        variant=variant_value if include_variant else None,
    )

    with open(csv_out_long, "w", newline="", encoding="utf-8") as f:
        writer = _csv.writer(f)

        header = ["run_id"]
        if include_variant:
            header.append("variant")
        header.extend(
            [
                "stage",
                "ligand",
                "ledock_best_score_kcal",
                "ledock_cluster_count",
                "ledock_n_poses",
            ]
        )
        writer.writerow(header)

        def _fmt(val: Any, places: int = 2) -> str:
            return f"{val:.{places}f}" if isinstance(val, (int, float)) else ""

        for stage_name, stage_map in ledock_metrics_by_stage.items():
            for lig, rec in stage_map.items():
                lig_key = os.path.basename(str(lig))
                best = rec.get("best_score_kcal")
                cluster_count = rec.get("cluster_count", None)
                n_poses = rec.get("n_poses", None)

                row = [run_id_value]
                if include_variant:
                    row.append(variant_value)
                row.extend(
                    [
                        stage_name,
                        lig_key,
                        _fmt(best),
                        cluster_count if cluster_count is not None else "",
                        n_poses if n_poses is not None else "",
                    ]
                )
                writer.writerow(row)

    logger.info(
        "[ledock.csv] pdb=%s ph=%s variant=%s summary=%s long=%s",
        pdb_id,
        ph_token,
        variant_token,
        csv_out_wide,
        csv_out_long,
    )

    return csv_out_wide
