from __future__ import annotations

import json
import logging
import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

import pandas as pd

from chemdb.target_difficulty import TargetDifficulty
from config.normalize import _to_bool
from docking.score_io import write_score_summary_to_csv
from analysis.dud_eval import guess_ligfile_col, guess_score_col, parse_name_and_label
from path_router.path_router import (
    make_paths,
    config_dir as router_config_dir,
    config_file as router_config_file,
    docked_dir as router_docked_dir,
    receptor_file as router_receptor_file,
)

HARD_BUCKETS = {"hard", "degenerate"}


def _default_retry_recipes() -> Dict[str, List[Dict[str, Any]]]:
    return {
        "too_far_from_pocket": [
            {"recenter": True, "box_pad_delta": +1.0, "num_modes": 4},
            {
                "recenter": True,
                "box_pad_delta": +2.0,
                "exhaustiveness": 6,
                "num_modes": 4,
            },
        ],
        "no_valid_pose": [
            {
                "exhaustiveness": 6,
                "num_modes": 5,
                "seed_jitter": True,
                "energy_range": 6,
            },
            {
                "recenter": True,
                "box_pad_delta": +1.0,
                "exhaustiveness": 6,
                "num_modes": 5,
                "seed_jitter": True,
                "energy_range": 6,
            },
        ],
        "timeout": [
            {"exhaustiveness": 3, "num_modes": 3, "seed_jitter": True},
            {"exhaustiveness": 2, "num_modes": 2},
        ],
        "malformed": [],
    }


def _compute_gnina_decoy_stats_from_long_csv(
    csv_path: Path | str,
    lig_col_override: Optional[str] = None,
    score_col_override: Optional[str] = "gnina_primary_score",
    *,
    assume_unlabeled_decoys: bool = False,
) -> tuple[float, float, int]:
    csv_path = Path(csv_path)
    if not csv_path.exists():
        return float("nan"), float("nan"), 0

    df = pd.read_csv(csv_path)
    if df.empty:
        return float("nan"), float("nan"), 0

    lig_col = guess_ligfile_col(df, lig_col_override)
    score_col = (
        score_col_override
        if score_col_override in df.columns
        else guess_score_col(df, score_col_override)
    )
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")

    parsed = df[lig_col].astype(str).apply(parse_name_and_label)
    df["lig_id"] = parsed.apply(lambda t: t[0])
    df["is_active"] = parsed.apply(lambda t: t[1])
    if assume_unlabeled_decoys and int(df["is_active"].notna().sum()) == 0:
        df["is_active"] = 0

    best = df.groupby("lig_id", as_index=False).agg(
        best_score=(score_col, "max"),
        is_active=("is_active", "max"),
    )
    decoys = best[best["is_active"] == 0]
    n_decoys = int(len(decoys))
    if n_decoys == 0:
        return float("nan"), float("nan"), 0

    mu = float(decoys["best_score"].mean())
    sigma = float(decoys["best_score"].std(ddof=1))
    return mu, sigma, n_decoys


def annotate_gnina_fda_long_csv_with_z_scores_vs_decoys(
    cfg: Dict[str, Any],
    pdb_id: str,
    ph_label: Optional[str] = None,
    *,
    csv_prefix: str = "",
    decoy_csv_prefix: str = "dud_",
    logger=None,
) -> Optional[str]:
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    var = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
    ph_token = (ph_label or "").strip() or None
    variant_root = Path(paths.docked_variant_root(var, ph_token))

    dud_csv = variant_root / f"{decoy_csv_prefix}gnina_docking_score_long.csv"
    fda_csv = variant_root / f"{csv_prefix}gnina_docking_score_long.csv"

    if not dud_csv.exists() or not fda_csv.exists():
        if logger:
            logger.info(
                "[gnina.z-score.skip] pdb_id=%s ph=%s reason=missing_csv dud=%s fda=%s",
                pdb_id,
                ph_label or "base",
                str(dud_csv),
                str(fda_csv),
            )
        return None

    mu, sigma, n_decoys = _compute_gnina_decoy_stats_from_long_csv(
        dud_csv,
        score_col_override="gnina_primary_score",
        assume_unlabeled_decoys=True,
    )
    if (
        not n_decoys
        or not math.isfinite(mu)
        or not math.isfinite(sigma)
        or sigma == 0.0
    ):
        if logger:
            logger.info(
                "[gnina.z-score.skip] pdb_id=%s ph=%s reason=degenerate_stats n=%s mu=%s sigma=%s",
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
    score_col = (
        "gnina_primary_score"
        if "gnina_primary_score" in df.columns
        else guess_score_col(df, None)
    )
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")

    best = df.groupby(lig_col, as_index=False).agg(best_score=(score_col, "max"))
    best["gnina_z_vs_decoys"] = (best["best_score"] - mu) / sigma
    z_map = dict(zip(best[lig_col], best["gnina_z_vs_decoys"]))
    df["gnina_z_vs_decoys"] = df[lig_col].map(z_map)
    df["gnina_t_vs_decoys"] = df["gnina_z_vs_decoys"]
    df.to_csv(fda_csv, index=False)

    if logger:
        logger.info(
            "[gnina.z-score.ok] pdb_id=%s ph=%s n_decoys=%s mean=%.3f std=%.3f out=%s",
            pdb_id,
            ph_label or "base",
            n_decoys,
            mu,
            sigma,
            str(fda_csv),
        )
    return str(fda_csv)


def annotate_gnina_fda_long_csv_with_t_scores_vs_decoys(
    cfg: Dict[str, Any],
    pdb_id: str,
    ph_label: Optional[str] = None,
    *,
    csv_prefix: str = "",
    decoy_csv_prefix: str = "dud_",
    logger=None,
) -> Optional[str]:
    return annotate_gnina_fda_long_csv_with_z_scores_vs_decoys(
        cfg,
        pdb_id,
        ph_label=ph_label,
        csv_prefix=csv_prefix,
        decoy_csv_prefix=decoy_csv_prefix,
        logger=logger,
    )


def emit_gnina_config(
    cfg: Dict[str, Any],
    pdb_id: str,
    receptor_pdbqt: str,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    ligand_path: str,
    stage_name: str,
    stage_info: Dict[str, Any],
    cpu_per_job: int,
    logger: Optional[Any] = None,
    *,
    variant: Optional[str] = None,
    ph_token: Optional[str] = None,
    legacy: bool = False,
):
    make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")

    variant_token = (
        (str(variant).strip().upper() or None) if variant is not None else None
    )
    ph_label = (str(ph_token).strip() or None) if ph_token is not None else None
    legacy_mode = bool(legacy)

    lig_base = Path(ligand_path).stem
    run_id = cfg["RUN_ID"]

    cfg_dir = router_config_dir(
        run_id,
        pdb_id,
        stage_name,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    cfg_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = router_config_file(
        run_id,
        pdb_id,
        stage_name,
        variant=variant_token,
        ph_tag=ph_label,
        name="gnina.json",
        legacy=legacy_mode,
    )

    cfg_path = cfg_dir / f"{lig_base}_{stage_name}.txt"

    stage_root = router_docked_dir(
        pdb_id,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    out_dir = stage_root / stage_name
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{lig_base}_{stage_name}.pdbqt"

    expected_receptor = router_receptor_file(
        pdb_id,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    receptor_for_config = str(expected_receptor)
    receptor_exists = expected_receptor.exists()

    variant_display = variant_token or "None"
    ph_display = ph_label or "None"
    breadcrumb = (
        "[cfg.emit] run=%s pdb=%s stage=%s variant=%s ph=%s\n"
        "           cfg_dir=%s receptor=%s out_root=%s"
    )
    breadcrumb_args = (
        run_id,
        pdb_id,
        stage_name,
        variant_display,
        ph_display,
        str(cfg_dir),
        receptor_for_config,
        str(stage_root),
    )
    if logger:
        logger.info(breadcrumb, *breadcrumb_args)
    else:
        print(breadcrumb % breadcrumb_args)

    if not receptor_exists:
        msg = (
            f"[router.error] missing receptor for pdb={pdb_id} variant={variant_display} "
            f"ph={ph_display} -> {expected_receptor}"
        )
        if logger:
            logger.error(msg)
        else:
            print(msg)

    lines = [
        f"receptor = {receptor_for_config}",
        f"ligand   = {ligand_path}",
        f"center_x = {center[0]:.3f}",
        f"center_y = {center[1]:.3f}",
        f"center_z = {center[2]:.3f}",
        f"size_x   = {box_size[0]:.3f}",
        f"size_y   = {box_size[1]:.3f}",
        f"size_z   = {box_size[2]:.3f}",
        f"cpu      = {int(cpu_per_job)}",
        f"exhaustiveness = {int(stage_info.get('exhaustiveness', 8))}",
        f"energy_range   = {int(stage_info.get('energy_range', 4))}",
        f"num_modes      = {int(stage_info.get('num_modes', 4))}",
        f"verbosity      = {int(stage_info.get('verbosity', 0))}",
        f"out = {out_path}",
    ]
    if "seed" in stage_info:
        lines.append(f"seed = {int(stage_info['seed'])}")

    payload = ("\n".join(lines)).encode("utf-8")
    overwrite = cfg_path.exists()
    tmp = cfg_path.with_suffix(".part")
    with open(tmp, "wb") as f:
        f.write(payload)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, cfg_path)

    manifest_data: Dict[str, Any]
    entries_map: Dict[str, Dict[str, Any]] = {}
    if manifest_path.exists():
        try:
            manifest_data = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception:
            manifest_data = {}
    else:
        manifest_data = {}

    entries = manifest_data.get("entries") if isinstance(manifest_data, dict) else None
    if isinstance(entries, list):
        for item in entries:
            if isinstance(item, dict):
                lig = str(item.get("ligand", ""))
                if lig:
                    entries_map[lig] = item

    entries_map[lig_base] = {
        "ligand": lig_base,
        "config": str(cfg_path),
        "out": str(out_path),
        "receptor": receptor_for_config,
    }
    manifest_data = {
        "run_id": run_id,
        "pdb_id": pdb_id,
        "stage": stage_name,
        "variant": variant_token,
        "ph": ph_label,
        "legacy": legacy_mode,
        "entries": [entries_map[k] for k in sorted(entries_map.keys())],
    }

    manifest_tmp = manifest_path.with_suffix(".part")
    try:
        manifest_tmp.parent.mkdir(parents=True, exist_ok=True)
        with open(manifest_tmp, "w", encoding="utf-8") as fh:
            json.dump(manifest_data, fh, indent=2, sort_keys=True)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(manifest_tmp, manifest_path)
    except FileNotFoundError:
        if logger:
            logger.warning(
                "[cfg.emit] manifest_tmp missing; skipping manifest update tmp=%s dest=%s stage=%s",
                str(manifest_tmp),
                str(manifest_path),
                stage_name,
            )
    except Exception:
        if logger:
            logger.exception(
                "[cfg.emit] manifest update failed; continuing without manifest stage=%s path=%s",
                stage_name,
                str(manifest_path),
            )

    emit_msg = (
        "[cfg.emit] run=%s pdb=%s variant=%s ph=%s stage=%s ligand=%s "
        "cfg_dir=%s docked_root=%s path=%s overwrite=%s bytes=%d"
    )
    emit_args = (
        run_id,
        pdb_id,
        variant_display,
        ph_display,
        stage_name,
        lig_base,
        str(cfg_dir),
        str(stage_root),
        str(cfg_path),
        str(overwrite).lower(),
        len(payload),
    )
    if logger:
        logger.info(emit_msg, *emit_args)
    else:
        print(emit_msg % emit_args)

    return str(cfg_path), str(out_path)


def _map_reason_to_category(reason: str) -> str:
    if not reason:
        return "no_valid_pose"
    r = str(reason).lower()
    if "timeout" in r:
        return "timeout"
    if "too far" in r or "distance" in r or "centroid" in r:
        return "too_far_from_pocket"
    if "malformed" in r or "parse" in r or "format" in r:
        return "malformed"
    if "no pose" in r or "no_valid" in r or "all_poses_invalid" in r:
        return "no_valid_pose"
    return "no_valid_pose"


@dataclass
class RetryManager:
    max_retries: int = 2
    recipes: Dict[str, List[Dict[str, Any]]] = field(
        default_factory=_default_retry_recipes
    )

    def apply(
        self, base_params: Dict[str, Any], err_type: str, attempt: int
    ) -> Optional[Dict[str, Any]]:
        if err_type not in self.recipes or attempt >= len(self.recipes[err_type]):
            return None
        params = base_params.copy()
        for key, value in self.recipes[err_type][attempt].items():
            if key.endswith("_delta"):
                target_key = key.replace("_delta", "")
                params[target_key] = params.get(target_key, 0.0) + value
            else:
                params[key] = value
        return params


def should_run_gnina_for_target(
    td: Optional[TargetDifficulty], cfg: Mapping[str, Any]
) -> bool:
    raw_flag = None
    if cfg is not None:
        file_cfg = cfg.get("_FILE_CFG") if isinstance(cfg, dict) else None
        raw_flag = cfg.get("USE_GNINA", cfg.get("use_gnina", None))
        if raw_flag is None and isinstance(file_cfg, dict):
            raw_flag = file_cfg.get("USE_GNINA", file_cfg.get("use_gnina", None))
    if raw_flag is not None and not _to_bool(raw_flag):
        logging.getLogger(__name__).info(
            "[gnina.disabled] USE_GNINA=%r -> GNINA follow-up docking is globally disabled",
            raw_flag,
        )
        return False

    if not cfg.get("ENABLE_GNINA", True):
        return False
    if not (cfg.get("GNINA_EXE") or "").strip():
        return False
    if td is None:
        return False
    diff = (getattr(td, "difficulty", "") or "").strip().lower()
    return diff in HARD_BUCKETS


def _resolve_gnina_stage_params(
    cfg: Mapping[str, Any],
    stage: Mapping[str, Any],
) -> Dict[str, Any]:
    stage_for_gnina = dict(stage)
    stage_for_gnina["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))
    if cfg.get("FAST_MODE"):
        stage_for_gnina["exhaustiveness"] = 1
        stage_for_gnina["num_modes"] = 1
    return stage_for_gnina


def write_gnina_scores_csv(
    cfg: Dict[str, Any],
    pdb_id: str,
    gnina_metrics_by_stage: Dict[str, Dict[str, Dict[str, Any]]],
    *,
    ph_label: Optional[str] = None,
    variant: Optional[str] = None,
    csv_prefix: str = "",
) -> str:
    import csv as _csv

    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_token = (ph_label or "").strip() or None
    variant_env = (
        (variant or os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    )
    variant_token = variant_env or None
    dock_dir = paths.docked_variant_root(variant_token, ph_token)
    dock_dir.mkdir(parents=True, exist_ok=True)

    run_id_value = str(cfg.get("RUN_ID") or "")
    variant_value = variant_env
    include_variant = bool(variant_value)

    csv_out_wide = str(dock_dir / f"{csv_prefix}gnina_docking_score_summary.csv")
    flat: Dict[str, Dict[str, Any]] = {}
    for stage_name, stage_map in gnina_metrics_by_stage.items():
        flat[stage_name] = {}
        for lig, rec in stage_map.items():
            lig_key = os.path.basename(lig)
            primary = rec.get("gnina_primary_score")
            if rec.get("valid"):
                flat[stage_name][lig_key] = (
                    f"{primary:.2f}" if isinstance(primary, (int, float)) else ""
                )
            elif isinstance(primary, (int, float)):
                flat[stage_name][lig_key] = f"{primary:.2f} (invalid)"
            else:
                flat[stage_name][lig_key] = "(invalid)"

    write_score_summary_to_csv(
        flat,
        output_path=csv_out_wide,
        run_id=run_id_value,
        variant=variant_value if include_variant else None,
    )

    csv_out_long = str(dock_dir / f"{csv_prefix}gnina_docking_score_long.csv")
    with open(csv_out_long, "w", newline="", encoding="utf-8") as f:
        writer = _csv.writer(f)
        header = ["run_id"]
        if include_variant:
            header.append("variant")
        header.extend(
            [
                "stage",
                "ligand",
                "valid",
                "reason",
                "heavy_atoms",
                "le",
                "self_rmsd",
                "pains_flag",
                "gnina_minimized_affinity_kcal",
                "gnina_cnn_score",
                "gnina_cnn_affinity_pK",
                "gnina_primary_score",
            ]
        )
        writer.writerow(header)

        for stage_name, stage_map in gnina_metrics_by_stage.items():
            for lig, rec in stage_map.items():
                lig_key = os.path.basename(lig)
                valid = bool(rec.get("valid", False))
                reason = rec.get("reason", "")
                ha = rec.get("heavy_atoms", None)
                le = rec.get("le", None)
                sr = rec.get("self_rmsd", None)
                pains_hit = rec.get("pains_flag", False)
                minimized_affinity = rec.get("minimized_affinity_kcal")
                cnn_score = rec.get("cnn_score")
                cnn_affinity = rec.get("cnn_affinity_pK")
                primary = rec.get("gnina_primary_score")

                def _fmt(val: Any, places: int = 2) -> str:
                    return f"{val:.{places}f}" if isinstance(val, (int, float)) else ""

                row: List[str] = [run_id_value]
                if include_variant:
                    row.append(variant_value)
                row.extend(
                    [
                        stage_name,
                        lig_key,
                        str(int(valid)),
                        str(reason) if reason is not None else "",
                        str(int(ha)) if isinstance(ha, (int, float)) else "",
                        f"{le:.4f}" if isinstance(le, (int, float)) else "",
                        f"{sr:.2f}"
                        if isinstance(sr, (int, float)) and math.isfinite(sr)
                        else "",
                        str(int(bool(pains_hit))),
                        _fmt(minimized_affinity),
                        _fmt(cnn_score, places=3),
                        _fmt(cnn_affinity, places=3),
                        _fmt(primary, places=3),
                    ]
                )
                writer.writerow(row)

    return csv_out_wide


__all__ = [
    "RetryManager",
    "_default_retry_recipes",
    "_map_reason_to_category",
    "_resolve_gnina_stage_params",
    "annotate_gnina_fda_long_csv_with_t_scores_vs_decoys",
    "annotate_gnina_fda_long_csv_with_z_scores_vs_decoys",
    "emit_gnina_config",
    "should_run_gnina_for_target",
    "write_gnina_scores_csv",
]
