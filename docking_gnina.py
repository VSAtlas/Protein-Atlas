from __future__ import annotations

import logging
import json
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

import pandas as pd
from tqdm import tqdm

from input_and_export_functions import _to_bool, extract_gnina_scores, write_score_summary_to_csv
from chemdb.target_difficulty import TargetDifficulty
from docking_ligands import compute_rmsd, validate_ligand
from pose_validation import (
    compute_self_rmsd,
    extract_surface_atoms,
    filter_and_rewrite_poses_by_rmsd,
    _pose_centroid_from_pdbqt,
)
from fallback_recenter import BudgetGuard, validate_first_valid_pose
from dataclasses import dataclass, field
from docking_utils import run_completion_audit
from path_router import (
    make_paths,
    config_dir as router_config_dir,
    config_file as router_config_file,
    docked_dir as router_docked_dir,
    receptor_file as router_receptor_file,
)
from dud_eval import guess_ligfile_col, guess_score_col, parse_name_and_label

HARD_BUCKETS = {"hard", "degenerate"}


def _compute_gnina_decoy_stats_from_long_csv(
    csv_path: Path | str,
    lig_col_override: Optional[str] = None,
    score_col_override: Optional[str] = "gnina_primary_score",
) -> tuple[float, float, int]:
    """
    Compute (mean, std, n) of best decoy scores from a GNINA long CSV.

    GNINA scores are higher-is-better, so we take the max per ligand.
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        return float("nan"), float("nan"), 0

    df = pd.read_csv(csv_path)
    if df.empty:
        return float("nan"), float("nan"), 0

    lig_col = guess_ligfile_col(df, lig_col_override)
    score_col = score_col_override if score_col_override in df.columns else guess_score_col(df, score_col_override)
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")

    parsed = df[lig_col].astype(str).apply(parse_name_and_label)
    df["lig_id"] = parsed.apply(lambda t: t[0])
    df["is_active"] = parsed.apply(lambda t: t[1])

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


def annotate_gnina_fda_long_csv_with_t_scores_vs_decoys(
    cfg: Dict[str, Any],
    pdb_id: str,
    ph_label: Optional[str] = None,
    logger=None,
) -> Optional[str]:
    """
    Mirror the Vina T-score annotation, but using GNINA long CSVs:
      - dud_gnina_docking_score_long.csv provides decoy stats (higher-is-better).
      - gnina_docking_score_long.csv gets annotated with gnina_t_vs_decoys.
    """
    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    var = (os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper() or None
    ph_token = (ph_label or "").strip() or None
    variant_root = Path(paths.docked_variant_root(var, ph_token))

    dud_csv = variant_root / "dud_gnina_docking_score_long.csv"
    fda_csv = variant_root / "gnina_docking_score_long.csv"

    if not dud_csv.exists() or not fda_csv.exists():
        if logger:
            logger.info(
                "[gnina.t-score.skip] pdb_id=%s ph=%s reason=missing_csv dud=%s fda=%s",
                pdb_id,
                ph_label or "base",
                str(dud_csv),
                str(fda_csv),
            )
        return None

    mu, sigma, n_decoys = _compute_gnina_decoy_stats_from_long_csv(
        dud_csv, score_col_override="gnina_primary_score"
    )
    if (
        not n_decoys
        or not math.isfinite(mu)
        or not math.isfinite(sigma)
        or sigma == 0.0
    ):
        if logger:
            logger.info(
                "[gnina.t-score.skip] pdb_id=%s ph=%s reason=degenerate_stats n=%s mu=%s sigma=%s",
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
    score_col = "gnina_primary_score" if "gnina_primary_score" in df.columns else guess_score_col(df, None)
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")

    best = df.groupby(lig_col, as_index=False).agg(best_score=(score_col, "max"))
    best["gnina_t_vs_decoys"] = (best["best_score"] - mu) / sigma
    t_map = dict(zip(best[lig_col], best["gnina_t_vs_decoys"]))
    df["gnina_t_vs_decoys"] = df[lig_col].map(t_map)

    df.to_csv(fda_csv, index=False)

    if logger:
        logger.info(
            "[gnina.t-score.ok] pdb_id=%s ph=%s n_decoys=%s mean=%.3f std=%.3f out=%s",
            pdb_id,
            ph_label or "base",
            n_decoys,
            mu,
            sigma,
            str(fda_csv),
        )

    return str(fda_csv)


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
    """
    GNINA config writer mirroring emit_vina_config.

    Key differences:
    - Writes directory manifest as gnina.json (not vina.json).
    - Keeps the same router path layout so downstream indexing remains stable.
    """
    make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")

    variant_token = (str(variant).strip().upper() or None) if variant is not None else None
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
    receptor_exists = expected_receptor.exists()
    receptor_for_config = str(expected_receptor)

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
    if logger:
        cx, cy, cz = center
        sx, sy, sz = box_size
        lig_name = os.path.basename(str(ligand_path))
        logger.info(
            "[gnina.cfg] lig=%s center=(%.3f,%.3f,%.3f) size=(%.1f,%.1f,%.1f)",
            lig_name,
            cx,
            cy,
            cz,
            sx,
            sy,
            sz,
        )

    payload = ("\n".join(lines)).encode("utf-8")
    overwrite = cfg_path.exists()

    tmp = cfg_path.with_suffix(".part")
    with open(tmp, "wb") as f:
        f.write(payload)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, cfg_path)

    manifest_data: Dict[str, Any]
    entries_map: Dict[str, Dict[str, Any]]
    if manifest_path.exists():
        try:
            manifest_data = json.loads(manifest_path.read_text(encoding="utf-8"))
        except Exception:
            manifest_data = {}
    else:
        manifest_data = {}

    entries = manifest_data.get("entries") if isinstance(manifest_data, dict) else None
    entries_map = {}
    if isinstance(entries, list):
        for item in entries:
            if isinstance(item, dict):
                lig = str(item.get("ligand", ""))
                if lig:
                    entries_map[lig] = item

    entry = {
        "ligand": lig_base,
        "config": str(cfg_path),
        "out": str(out_path),
        "receptor": receptor_for_config,
    }
    entries_map[lig_base] = entry

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
    recipes: Dict[str, List[Dict[str, Any]]] = field(default_factory=lambda: {
        "too_far_from_pocket": [
            {"recenter": True, "box_pad_delta": +1.0, "num_modes": 4},
            {"recenter": True, "box_pad_delta": +2.0, "exhaustiveness": 6, "num_modes": 4},
        ],
        "no_valid_pose": [
            {"exhaustiveness": 6, "num_modes": 5, "seed_jitter": True, "energy_range": 6},
            {"recenter": True, "box_pad_delta": +1.0, "exhaustiveness": 6, "num_modes": 5, "seed_jitter": True, "energy_range": 6},
        ],
        "timeout": [
            {"exhaustiveness": 3, "num_modes": 3, "seed_jitter": True},
            {"exhaustiveness": 2, "num_modes": 2},
        ],
        "malformed": [],
    })

    def apply(self, base_params: Dict[str, Any], err_type: str, attempt: int) -> Optional[Dict[str, Any]]:
        if err_type not in self.recipes or attempt >= len(self.recipes[err_type]):
            return None
        p = base_params.copy()
        for k, v in self.recipes[err_type][attempt].items():
            if k.endswith("_delta"):
                key = k.replace("_delta", "")
                p[key] = p.get(key, 0.0) + v
            else:
                p[k] = v
        return p


def should_run_gnina_for_target(td: Optional[TargetDifficulty], cfg: Mapping[str, Any]) -> bool:
    """
    Decide whether to run GNINA based on difficulty and config.

    Honors USE_GNINA/use_gnina (including _FILE_CFG), ENABLE_GNINA, GNINA_EXE,
    and only runs for hard/degenerate targets.
    """
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
    """
    Return GNINA-specific stage parameters derived from the base Vina stage.
    """
    stage_for_gnina = dict(stage)
    stage_for_gnina["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))
    if cfg.get("FAST_MODE"):
        stage_for_gnina["exhaustiveness"] = 1
        stage_for_gnina["num_modes"] = 1
    return stage_for_gnina


def write_gnina_scores_csv(
    cfg: Dict,
    pdb_id: str,
    gnina_metrics_by_stage: Dict[str, Dict[str, Dict[str, Any]]],
    *,
    ph_label: Optional[str] = None,
    variant: Optional[str] = None,
    csv_prefix: str = "",
) -> str:
    import csv as _csv
    import math as _math

    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_token = (ph_label or "").strip() or None
    variant_env = (variant or os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    variant_token = variant_env or None
    dock_dir = paths.docked_variant_root(variant_token, ph_token)
    dock_dir.mkdir(parents=True, exist_ok=True)

    run_id_value = str(cfg.get("RUN_ID") or "")
    variant_value = variant_env
    include_variant = bool(variant_value)

    summary_name = "gnina_docking_score_summary.csv"
    long_name = f"{csv_prefix}gnina_docking_score_long.csv"
    csv_out_wide = str(dock_dir / summary_name)

    flat: Dict[str, Dict[str, Any]] = {}
    for stage_name, stage_map in gnina_metrics_by_stage.items():
        flat[stage_name] = {}
        for lig, rec in stage_map.items():
            lig_key = os.path.basename(lig)
            primary = rec.get("gnina_primary_score")
            if rec.get("valid"):
                flat[stage_name][lig_key] = f"{primary:.2f}" if isinstance(primary, (int, float)) else ""
            else:
                if isinstance(primary, (int, float)):
                    flat[stage_name][lig_key] = f"{primary:.2f} (invalid)"
                else:
                    flat[stage_name][lig_key] = "(invalid)"

    write_score_summary_to_csv(
        flat,
        output_path=csv_out_wide,
        run_id=run_id_value,
        variant=variant_value if include_variant else None,
    )

    csv_out_long = str(dock_dir / long_name)
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

                ha_str = str(int(ha)) if isinstance(ha, (int, float)) else ""
                le_str = f"{le:.4f}" if isinstance(le, (int, float)) else ""
                sr_str = f"{sr:.2f}" if isinstance(sr, (int, float)) and _math.isfinite(sr) else ""
                reason_str = str(reason) if reason is not None else ""

                def _fmt(val: Any, places: int = 2) -> str:
                    return f"{val:.{places}f}" if isinstance(val, (int, float)) else ""

                row = [run_id_value]
                if include_variant:
                    row.append(variant_value)
                row.extend(
                    [
                        stage_name,
                        lig_key,
                        int(valid),
                        reason_str,
                        ha_str,
                        le_str,
                        sr_str,
                        int(bool(pains_hit)),
                        _fmt(minimized_affinity),
                        _fmt(cnn_score, places=3),
                        _fmt(cnn_affinity, places=3),
                        _fmt(primary, places=3),
                    ]
                )
                writer.writerow(row)

    return csv_out_wide


def _run_gnina_for_ligand(
    cfg: Dict[str, Any],
    pdb_id: str,
    receptor_pdbqt: str,
    lig: str,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    stage: Dict[str, Any],
    variant_token: Optional[str],
    ph_label: Optional[str],
    legacy_mode: bool,
    logger: logging.Logger,
    threads_per_vina: int,
) -> Optional[str]:
    """
    Minimal GNINA follow-up docking for a single ligand.

    - Writes config under gnina_<stage_name> and runs GNINA with --no_gpu --cpu 1.
    - Returns the GNINA out_path on success, or None on failure.
    """
    gnina_exe = (cfg.get("GNINA_EXE") or "").strip()
    if not gnina_exe:
        logger.debug("[gnina.skip] GNINA_EXE not set; skipping GNINA follow-up.")
        return None

    # For fast test targets that use the small test_library_10, allow GNINA to use
    # more CPU threads to keep CI/local tests quick. This must be reflected both
    # in the emitted config and the CLI invocation.
    test_pdbs = {"TEST", "T3MP", "TEMP"}
    cpu_threads = 10 if pdb_id.strip().upper() in test_pdbs else 1

    # Build GNINA-specific stage params without mutating the shared stage dict.
    stage_for_gnina = _resolve_gnina_stage_params(cfg, stage)

    base_stage_name = stage.get("name") or "stage"
    gnina_stage_name = f"gnina_{base_stage_name}"

    logger.info(
        "[gnina.emit] pdb=%s stage=%s variant=%s ph=%s lig=%s",
        pdb_id,
        gnina_stage_name,
        variant_token or "None",
        ph_label or "None",
        os.path.basename(lig),
    )

    conf_path, out_path = emit_gnina_config(
        cfg,
        pdb_id,
        receptor_pdbqt,
        center,
        box_size,
        lig,
        gnina_stage_name,
        stage_for_gnina,
        cpu_threads,  # GNINA threads (test targets may use more)
        logger,
        variant=variant_token,
        ph_token=ph_label,
        legacy=legacy_mode,
    )

    # GNINA accepts a mostly Vina-compatible config, but rejects some Vina-only
    # options like energy_range/verbosity. Strip unsupported keys in-place to
    # keep router paths identical while ensuring GNINA runs.
    try:
        conf_p = Path(conf_path)
        raw_lines = conf_p.read_text(encoding="utf-8").splitlines()
        filtered_lines: List[str] = []
        removed_keys: List[str] = []
        for line in raw_lines:
            key = line.split("=", 1)[0].strip().lower()
            if key in {"energy_range", "verbosity"}:
                removed_keys.append(key)
                continue
            filtered_lines.append(line)
        if removed_keys:
            tmp = conf_p.with_suffix(".part")
            tmp.write_text("\n".join(filtered_lines) + "\n", encoding="utf-8")
            os.replace(tmp, conf_p)
            logger.debug(
                "[gnina.cfg.strip] removed=%s path=%s",
                ",".join(sorted(set(removed_keys))),
                conf_path,
            )
    except Exception as e:
        logger.warning("[gnina.cfg.strip] failed path=%s err=%s", conf_path, e)

    try:
        Path(conf_path).resolve().relative_to(Path(cfg["CONFIG_RUN_DIR"]).resolve())
    except Exception:
        raise RuntimeError(
            f"Refusing to launch GNINA with config outside current RUN_DIR: {conf_path}"
        )

    # Build GNINA CLI. In FAST_MODE, explicitly disable CNN scoring to keep runs fast
    # and easier to test, and sort poses by energy instead of CNN score.
    cmd = [gnina_exe, "--no_gpu", "--cpu", str(cpu_threads), "--config", conf_path]
    if cfg.get("FAST_MODE"):
        cmd.extend(["--cnn_scoring", "none", "--pose_sort_order", "energy"])

    logger.info(
        "[gnina.call] exe=%s config=%s cmd=%s",
        gnina_exe,
        conf_path,
        " ".join(cmd),
    )

    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except Exception as e:
        logger.warning(
            "[gnina.error] pdb=%s lig=%s stage=%s err=%s",
            pdb_id,
            os.path.basename(lig),
            gnina_stage_name,
            e,
        )
        return None

    return out_path


def run_gnina_for_stage(
    *,
    cfg: Dict[str, Any],
    paths: Any,
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    stage_name: str,
    stage_info: Mapping[str, Any],
    ligands: List[str],
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    logger: logging.Logger,
    receptor_pdbqt: Optional[str] = None,
    control_lookup: Optional[Dict[str, Path]] = None,
) -> tuple[Dict[str, Optional[float]], Dict[str, Dict[str, Optional[float]]], Dict[str, Any]]:
    """
    Run GNINA for the given stage and ligand list.

    - Writes GNINA config(s) into configs/.../gnina_<stage>/gnina.json
    - Writes docked poses into docked/.../gnina_<stage>/<ligand>.pdbqt
    - Uses GNINA_EXE from cfg
    - Uses --no_gpu and --cpu 1
    - Returns (primary_scores, metrics) where:
        primary_scores: {ligand_path: cnn_affinity or None}
        metrics: {
            ligand_path: {
                minimized_affinity_kcal, cnn_score, cnn_affinity_pK, gnina_primary_score
            }
        }
    """
    threads_per_vina = int(cfg.get("THREADS_PER_VINA", 1))
    cpu = int(cfg.get("CPU", os.cpu_count() or 1))
    max_jobs = int(cfg.get("MAX_PARALLEL_JOBS", cpu))
    max_workers = min(max_jobs, len(ligands)) if ligands else 1
    if max_workers < 1:
        max_workers = 1

    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    variant_token = (str(variant).strip().upper() or None) if variant is not None else None

    base_stage_name = stage_name
    gnina_stage_name = f"gnina_{base_stage_name}"
    stage_dir = paths.docked_stage_dir(variant_token, gnina_stage_name, ph_label)
    logger.info(
        "[gnina.stage] pdb=%s stage=%s variant=%s ph=%s n_lig=%d",
        pdb_id,
        gnina_stage_name,
        variant_token or "None",
        ph_label or "None",
        len(ligands),
    )

    receptor_pdbqt_path = receptor_pdbqt or str(paths.receptor_pdbqt(variant_token, ph_label))
    try:
        surface_coords = extract_surface_atoms(pdbqt_path=receptor_pdbqt_path, center=center)
    except Exception:
        surface_coords = None

    def _best_pose_pdb_from_pdbqt(pdbqt_path: str, obabel_path: Optional[str] = None) -> Optional[str]:
        """Convert first model of PDBQT -> PDB (no hydrogens) using OpenBabel."""
        try:
            from shutil import which

            obabel = obabel_path or which("obabel")
            if not obabel:
                return None

            out_pdb = Path(pdbqt_path).with_suffix(".best.pdb")
            cmd = [obabel, "-ipdbqt", pdbqt_path, "-opdb", "-O", str(out_pdb), "-d", "--first"]
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            return str(out_pdb)
        except Exception:
            return None

    scores: Dict[str, Optional[float]] = {}
    gnina_metrics: Dict[str, Dict[str, Optional[float]]] = {}
    gnina_futures: Dict[Any, str] = {}

    default_budget_seconds = float(cfg.get("MAX_RETRY_SECONDS_PER_LIGAND", cfg.get("BENCH_MAX_SECONDS", 300.0)))
    guards_for_ligand: Dict[str, BudgetGuard] = {}
    retry_mgr = RetryManager()

    def _dock_once(
        lig: str,
        stage_params: Dict[str, Any],
        center_use: tuple[float, float, float],
        box_use: tuple[float, float, float],
        stage_label: str,
    ) -> tuple[Optional[float], Dict[str, Optional[float]], Dict[str, Any], Optional[str]]:
        """
        Run GNINA once, parse scores, and perform pose validation/self-RMSD/control gate.
        Returns (primary_score, metrics, validation_result, out_path).
        """
        out_path = _run_gnina_for_ligand(
            cfg,
            pdb_id,
            str(paths.receptor_pdbqt(variant_token, ph_label)),
            lig,
            center_use,
            box_use,
            stage_params,
            variant_token,
            ph_label,
            legacy_mode,
            logger,
            threads_per_vina,
        )

        metrics: Dict[str, Optional[float]] = {
            "minimized_affinity_kcal": None,
            "cnn_score": None,
            "cnn_affinity_pK": None,
            "gnina_primary_score": None,
            "valid": False,
            "reason": "",
            "self_rmsd": None,
        }
        validation_result: Dict[str, Any] = {"valid": False, "reason": "gnina_failed"}

        if not out_path:
            metrics["reason"] = "gnina_failed"
            return None, metrics, validation_result, None

        try:
            parsed = extract_gnina_scores(out_path)
        except Exception:
            parsed = {
                "minimized_affinity_kcal": None,
                "cnn_score": None,
                "cnn_affinity_pK": None,
            }

        minimized_affinity = parsed.get("minimized_affinity_kcal")
        cnn_score = parsed.get("cnn_score")
        cnn_affinity = parsed.get("cnn_affinity_pK")
        # Primary score is CNN affinity when available; otherwise fall back to the
        # minimized empirical affinity (e.g. when FAST_MODE disables CNN scoring).
        primary_score = cnn_affinity if cnn_affinity is not None else minimized_affinity

        metrics.update(
            {
                "minimized_affinity_kcal": minimized_affinity,
                "cnn_score": cnn_score,
                "cnn_affinity_pK": cnn_affinity,
                "gnina_primary_score": primary_score,
            }
        )
        logger.info(
            "[gnina.out] pdb=%s stage=%s lig=%s out=%s cnn_affinity=%s cnn_score=%s minimized_affinity=%s",
            pdb_id,
            stage_label,
            os.path.basename(lig),
            out_path,
            cnn_affinity if cnn_affinity is not None else "None",
            cnn_score if cnn_score is not None else "None",
            minimized_affinity if minimized_affinity is not None else "None",
        )

        try:
            kept, removed = filter_and_rewrite_poses_by_rmsd(
                out_path,
                rmsd_tol=float(cfg.get("RMSD_FILTER_ANG", 2.0)),
                max_models=int(cfg.get("RMSD_MAX_MODELS", 3)),
            )
            if removed > 0:
                logger.info(f"{os.path.basename(out_path)}: RMSD filter kept {kept}, removed {removed}")
        except Exception as e:
            logger.warning(f"RMSD filtering failed for {os.path.basename(out_path)}: {e}")

        try:
            validation_result = validate_first_valid_pose(
                receptor_pdbqt=receptor_pdbqt_path,
                ligand_pdbqt=out_path,
                pocket_center=center_use,
                surface_coords=surface_coords,
                max_models=int(cfg.get("EARLY_EXIT_MAX_MODELS", 3)),
                clash_threshold=2.0,
                clash_tol=3,
                dist_surf=6.0,
                dist_centroid=4.5,
            )
        except Exception:
            validation_result = {"valid": False, "reason": "pose_invalid"}

        valid_pose = bool(validation_result.get("valid", False))
        reason_str = validation_result.get("reason", "") or ""

        self_rmsd_val = None
        if bool(cfg.get("LOG_SELF_RMSD", True)):
            try:
                self_rmsd_val = compute_self_rmsd(out_path)
                logger.info(f"[gnina.self-rmsd] lig={os.path.basename(lig)} rmsd={self_rmsd_val}")
            except Exception as _e:
                logger.warning(f"[gnina.self-rmsd] failed for {os.path.basename(lig)}: {_e}")

        ctrl_ref = None
        if control_lookup:
            base = Path(lig).stem.split("_stage")[0]
            ctrl_ref = control_lookup.get(base)
        if ctrl_ref:
            best_pdb = _best_pose_pdb_from_pdbqt(out_path, obabel_path=cfg.get("OPENBABEL_PATH"))
            if not best_pdb:
                valid_pose = False
                reason_str = "no_best_pose_for_rmsd"
            else:
                ok = validate_ligand(
                    ligand_name=os.path.basename(lig),
                    docked_path=best_pdb,
                    crystal_path=str(ctrl_ref),
                    rmsd_thresh=float(cfg.get("CONTROL_RMSD_MAX_ANG", 2.0)),
                    self_rmsd=self_rmsd_val,
                    logger=logger,
                )
                if not ok:
                    valid_pose = False
                    reason_str = "rmsd_fail"

        metrics["valid"] = bool(valid_pose)
        metrics["reason"] = reason_str
        metrics["self_rmsd"] = self_rmsd_val
        return primary_score, metrics, validation_result, out_path

    def _near_miss_retry(
        lig: str,
        base_stage: Dict[str, Any],
        result: Dict[str, Any],
        guard: BudgetGuard,
        last_out: Optional[str],
    ) -> Optional[Dict[str, Optional[float]]]:
        if not bool(cfg.get("RETRY_NEAR_MISS", True)):
            return None

        reason = result.get("reason", "") or ""
        near_miss = ("clash" in reason) or (result.get("distance_to_surface") or 0.0) < 6.5 or (
            result.get("distance_to_centroid") or 0.0
        ) < 4.0
        if not near_miss:
            return None

        stage_retry = dict(base_stage)
        stage_retry["name"] = f"{base_stage.get('name', base_stage_name)}_retry"
        stage_retry["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))
        if cfg.get("FAST_MODE"):
            stage_retry["exhaustiveness"] = 1

        try:
            base_seed = int(stage_retry.get("seed", 0)) if "seed" in stage_retry else 0
        except Exception:
            base_seed = 0
        stage_retry["seed"] = base_seed + 137
        stage_retry["num_modes"] = int(cfg.get("NEAR_MISS_NUM_MODES", stage_retry.get("num_modes", 4)))
        stage_retry["energy_range"] = float(cfg.get("NEAR_MISS_ENERGY_RANGE", stage_retry.get("energy_range", 4.0)))

        center_nm = center
        box_nm = box_size
        if bool(cfg.get("NEAR_MISS_RECENTER", True)) and last_out:
            try:
                cent = _pose_centroid_from_pdbqt(str(last_out))
            except Exception as _e:
                logger.warning(f"[near-miss] failed to compute centroid for {os.path.basename(lig)}: {_e}")
                cent = None
            if cent and isinstance(cent, (list, tuple)) and len(cent) == 3:
                try:
                    center_nm = tuple(float(x) for x in cent)
                    max_box = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
                    box_nm = tuple(min(max_box, s + 1.0) for s in box_size)
                except Exception:
                    center_nm = center
                    box_nm = box_size

        if guard.expired():
            return {
                "minimized_affinity_kcal": None,
                "cnn_score": None,
                "cnn_affinity_pK": None,
                "gnina_primary_score": None,
                "valid": False,
                "reason": "budget_exceeded",
                "self_rmsd": None,
            }

        primary_retry, metrics_retry, _, _ = _dock_once(
            lig,
            stage_retry,
            center_nm,
            box_nm,
            f"{stage_retry['name']}",
        )
        if metrics_retry.get("valid", False):
            logger.info(
                "%s | %s score: %s (gnina near-miss rescued)",
                os.path.basename(lig),
                stage_retry["name"],
                f"{primary_retry:.2f}" if isinstance(primary_retry, (int, float)) else "None",
            )
        return metrics_retry

    def _structured_retries(
        lig: str,
        base_stage: Dict[str, Any],
        result: Dict[str, Any],
        guard: BudgetGuard,
        last_out: Optional[str],
    ) -> Optional[Dict[str, Optional[float]]]:
        err_cat = _map_reason_to_category(result.get("reason", ""))
        attempt = 0
        metrics_candidate: Optional[Dict[str, Optional[float]]] = None
        last_out_path = last_out

        while attempt < retry_mgr.max_retries:
            if guard.expired():
                metrics_candidate = {
                    "minimized_affinity_kcal": None,
                    "cnn_score": None,
                    "cnn_affinity_pK": None,
                    "gnina_primary_score": None,
                    "valid": False,
                    "reason": "budget_exceeded",
                    "self_rmsd": None,
                }
                break

            recipe = retry_mgr.apply(base_stage, err_cat, attempt)
            if not recipe:
                break

            stage_retry2 = dict(base_stage)
            stage_retry2["name"] = f"{base_stage.get('name', base_stage_name)}_r{attempt + 1}"
            stage_retry2["verbosity"] = int(cfg.get("VINA_VERBOSITY", 0))

            if "exhaustiveness" in recipe:
                stage_retry2["exhaustiveness"] = recipe["exhaustiveness"]
            if "num_modes" in recipe:
                stage_retry2["num_modes"] = recipe["num_modes"]
            if recipe.get("energy_range") is not None:
                stage_retry2["energy_range"] = recipe["energy_range"]
            if cfg.get("FAST_MODE"):
                stage_retry2["exhaustiveness"] = 1
            if recipe.get("seed_jitter", False):
                try:
                    base_seed = int(stage_retry2.get("seed", 0)) if "seed" in stage_retry2 else 0
                except Exception:
                    base_seed = 0
                stage_retry2["seed"] = base_seed + (attempt + 1) * 137

            retry_center = center
            retry_box = box_size
            try:
                if recipe.get("recenter", False) and last_out_path:
                    cent = _pose_centroid_from_pdbqt(str(last_out_path))
                    if cent and isinstance(cent, (list, tuple)) and len(cent) == 3:
                        retry_center = tuple(float(x) for x in cent)
                if "box_pad_delta" in recipe and isinstance(recipe["box_pad_delta"], (int, float)):
                    dx = float(recipe["box_pad_delta"])
                    box_cap = float(cfg.get("BOX_SIZE_MAX_A", 28.0))
                    retry_box = tuple(min(box_cap, s + dx) for s in box_size)
            except Exception as _e:
                logger.warning(f"Retry recenter/box tweak failed: {_e}")

            primary_retry, metrics_retry, result_retry, out_retry = _dock_once(
                lig,
                stage_retry2,
                retry_center,
                retry_box,
                f"{stage_retry2['name']}",
            )
            last_out_path = out_retry or last_out_path
            metrics_candidate = metrics_retry
            if metrics_retry.get("valid", False):
                logger.info(
                    "%s | %s score: %s (gnina structured retry rescued)",
                    os.path.basename(lig),
                    stage_retry2["name"],
                    f"{primary_retry:.2f}" if isinstance(primary_retry, (int, float)) else "None",
                )
                break
            result = result_retry
            attempt += 1

        return metrics_candidate

    submit_ligands: List[str] = []
    for lig in ligands:
        guard = guards_for_ligand.get(lig) or BudgetGuard(default_budget_seconds, log=logger)
        guards_for_ligand[lig] = guard
        if guard.expired():
            scores[lig] = None
            gnina_metrics[lig] = {
                "minimized_affinity_kcal": None,
                "cnn_score": None,
                "cnn_affinity_pK": None,
                "gnina_primary_score": None,
                "valid": False,
                "reason": "budget_exceeded",
                "self_rmsd": None,
            }
            logger.info("[gnina.budget] pdb=%s lig=%s stage=%s reason=budget_exceeded", pdb_id, lig, gnina_stage_name)
            continue
        submit_ligands.append(lig)

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        for lig in submit_ligands:
            gnina_futures[
                pool.submit(
                    _dock_once,
                    lig,
                    dict(stage_info),
                    center,
                    box_size,
                    gnina_stage_name,
                )
            ] = lig

        with tqdm(
            total=len(gnina_futures),
            desc=f"GNINA ({base_stage_name})",
            unit="ligand",
            position=2,
            dynamic_ncols=True,
            mininterval=0.2,
            leave=True,
            file=sys.stdout,
        ) as pbar:
            for fut in as_completed(gnina_futures):
                lig = gnina_futures[fut]
                lig_name = os.path.basename(lig)
                try:
                    primary, metrics, result, out_gnina = fut.result()
                except Exception as e:
                    logger.warning(
                        "[gnina.error] pdb=%s lig=%s stage=%s err=%s",
                        pdb_id,
                        lig_name,
                        gnina_stage_name,
                        e,
                    )
                    primary, metrics, result, out_gnina = None, {
                        "minimized_affinity_kcal": None,
                        "cnn_score": None,
                        "cnn_affinity_pK": None,
                        "gnina_primary_score": None,
                        "valid": False,
                        "reason": "gnina_failed",
                        "self_rmsd": None,
                    }, {"valid": False, "reason": "gnina_failed"}, None

                guard = guards_for_ligand.get(lig) or BudgetGuard(default_budget_seconds, log=logger)
                guards_for_ligand[lig] = guard

                final_metrics = metrics
                if not metrics.get("valid", False):
                    if guard.expired():
                        final_metrics = dict(metrics)
                        final_metrics["reason"] = "budget_exceeded"
                    else:
                        near_miss_metrics = _near_miss_retry(lig, dict(stage_info), result, guard, out_gnina)
                        if near_miss_metrics and near_miss_metrics.get("valid", False):
                            final_metrics = near_miss_metrics
                        else:
                            structured_metrics = _structured_retries(lig, dict(stage_info), result, guard, out_gnina)
                            if structured_metrics and structured_metrics.get("valid", False):
                                final_metrics = structured_metrics
                            elif structured_metrics:
                                final_metrics = structured_metrics

                scores[lig] = final_metrics.get("gnina_primary_score")
                gnina_metrics[lig] = final_metrics
                pbar.update(1)

    # Completion audit: rerun missing GNINA outputs once before checkpointing.
    stage_params_for_gnina = _resolve_gnina_stage_params(cfg, stage_info)

    def _expected_gnina_path(lig: str) -> Path:
        return stage_dir / f"{Path(lig).stem}_{gnina_stage_name}.pdbqt"

    def _rerun_gnina_missing(lig: str) -> tuple[bool, Optional[str], Optional[Path]]:
        try:
            primary_r, metrics_r, _result_r, out_path = _dock_once(
                lig,
                dict(stage_params_for_gnina),
                center,
                box_size,
                gnina_stage_name,
            )
            gnina_metrics[lig] = metrics_r
            scores[lig] = metrics_r.get("gnina_primary_score")
        except Exception as exc:
            return False, f"rerun_error:{exc}", None
        if out_path and Path(out_path).exists() and Path(out_path).stat().st_size > 0:
            return True, "rerun_ok", None
        return False, "rerun_no_output", None

    completion_report_gnina = run_completion_audit(
        engine="gnina",
        pdb_id=pdb_id,
        stage_name=gnina_stage_name,
        ligands=ligands,
        expected_output_path=_expected_gnina_path,
        rerun_one=_rerun_gnina_missing,
        stage_dir=stage_dir,
        cfg=cfg,
        logger=logger,
        retries=1,
        ph_label=ph_label,
        variant=variant_token,
    )
    missing_after = completion_report_gnina.get("missing_ligands_after") or []
    if missing_after:
        for lig in missing_after:
            if lig not in gnina_metrics:
                gnina_metrics[lig] = {
                    "minimized_affinity_kcal": None,
                    "cnn_score": None,
                    "cnn_affinity_pK": None,
                    "gnina_primary_score": None,
                    "valid": False,
                    "reason": "completion_missing",
                    "self_rmsd": None,
                }
            scores.setdefault(lig, gnina_metrics[lig].get("gnina_primary_score"))

    return scores, gnina_metrics, completion_report_gnina
