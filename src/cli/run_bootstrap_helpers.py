from __future__ import annotations

import errno
import logging
import os
import shutil
import time
from pathlib import Path
from typing import Any, Dict

from config.output_paths import output_root
from config.runtime_config import load_inputs
from cli.run_context import ConfigDict
from docking.apo_holo_runtime import resolve_apo_holo_mode
from docking.docking_vina import emit_vina_config
from path_router.path_router import make_paths, receptor_file


def init_config_run_dir(cfg, run_id=None, reset=None, logger=None):
    root = Path(
        cfg.get("CONFIGS_DIR", output_root(Path(cfg["OVERALL_DIR"]), "configs"))
    )
    rid = run_id or cfg.get("RUN_ID")
    if not rid:
        rid = time.strftime("%Y%m%d_%H%M%S")
    cfg["RUN_ID"] = rid
    run_dir = root / rid
    cfg["CONFIG_RUN_DIR"] = str(run_dir)

    reset = bool(cfg.get("RESET_CONFIGS", True)) if reset is None else bool(reset)
    removed = 0
    if reset and run_dir.exists():
        for p in run_dir.rglob("*"):
            removed += 1
        try:
            shutil.rmtree(run_dir, ignore_errors=False)
        except FileNotFoundError:
            # Concurrent workers may race on reset; treat already-removed as success.
            pass
        except OSError as exc:
            if exc.errno == errno.ENOENT:
                # File vanished during traversal; safe to continue.
                pass
            else:
                raise
    run_dir.mkdir(parents=True, exist_ok=True)

    msg = f"[cfg.reset] run_dir={run_dir} removed={removed}"
    (logger.info(msg) if logger else print(msg))


# -------------------------
# P2Rank helpers (optional)
# -------------------------
def p2rank_out_dir(cfg: Dict[str, Any]) -> Path:
    """
    Central place to define where P2Rank writes results.
    Default: OVERALL_DIR/p2rank_out, but user may override P2RANK_OUTPUT_DIR.
    """
    out = Path(cfg.get("P2RANK_OUTPUT_DIR") or Path(cfg["OVERALL_DIR"]) / "p2rank_out")
    out.mkdir(parents=True, exist_ok=True)
    return out


# -------------------------
# Docking stages
# -------------------------
def define_docking_stages(mode="discovery"):
    if mode == "discovery":
        return [
            {"name": "stage1", "num_modes": 1, "energy_range": 2, "exhaustiveness": 2},
            {"name": "stage2", "num_modes": 3, "energy_range": 3, "exhaustiveness": 4},
            {"name": "stage3", "num_modes": 5, "energy_range": 4, "exhaustiveness": 6},
            {"name": "stage4", "num_modes": 9, "energy_range": 6, "exhaustiveness": 8},
            {
                "name": "stage5",
                "num_modes": 20,
                "energy_range": 9,
                "exhaustiveness": 20,
            },
        ]
    elif mode == "polypharmacology":
        return [
            {"name": "stage1", "num_modes": 3, "energy_range": 2, "exhaustiveness": 4},
            {
                "name": "stage2",
                "num_modes": 10,
                "energy_range": 6,
                "exhaustiveness": 12,
            },
            {
                "name": "stage3",
                "num_modes": 20,
                "energy_range": 9,
                "exhaustiveness": 24,
            },
        ]
    else:
        raise ValueError(f"Unknown docking mode: {mode}")


def smoke_emit_config_demo(repo_root: Path | None = None) -> None:
    """Emit a small config to exercise router paths in isolation."""
    smoke_log = logging.getLogger("smoke")
    old_variant = os.environ.get("APO_HOLO_VARIANT")
    old_run_env = os.environ.get("ATLAS_RUN_ID")
    try:
        base_cfg = ConfigDict(load_inputs())
    except Exception as exc:
        smoke_log.warning("[smoke.emit.skip] reason=%s", exc)
        return

    try:
        cfg = ConfigDict(base_cfg.copy())
        root = (
            Path(repo_root).resolve()
            if repo_root is not None
            else Path(__file__).resolve().parents[2]
        )
        smoke_root = root / "analysis" / "_smoke"
        overrides = {
            "OVERALL_DIR": smoke_root,
            "INPUT_DIR": smoke_root / "input_pdbs",
            "OUTPUT_DIR": smoke_root / "processed_pdbs",
            "DOCKED_DIR": smoke_root / "docked",
            "PREPPED_LIGANDS_DIR": smoke_root / "prepped_ligands",
            "LIGANDS_MOL2_DIR": smoke_root / "ligands_mol2",
            "CONFIGS_DIR": smoke_root / "configs",
        }
        for key, path_value in overrides.items():
            cfg[key] = str(path_value)
            Path(path_value).mkdir(parents=True, exist_ok=True)

        cfg["RUN_ID"] = "smoke_demo"
        cfg["RESET_CONFIGS"] = False
        init_config_run_dir(cfg, run_id=cfg["RUN_ID"], reset=False, logger=smoke_log)

        mode, variants = resolve_apo_holo_mode(cfg)
        smoke_log.info("[smoke.emit] mode=%s variants=%s", mode, variants)

        pdb_id = "3CS9"
        paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
        lig_dir = paths.prepped_ligands_dir
        lig_dir.mkdir(parents=True, exist_ok=True)
        lig_path = lig_dir / "smoke_ligand.pdbqt"
        if not lig_path.exists():
            lig_path.write_text("SMOKE", encoding="utf-8")

        ph_label = "pH6_7"
        os.environ["APO_HOLO_VARIANT"] = "APO"

        receptor_path = receptor_file(
            paths.pdb_id, variant="APO", ph_tag=ph_label, legacy=False
        )
        receptor_path.parent.mkdir(parents=True, exist_ok=True)
        if not receptor_path.exists():
            receptor_path.write_text("RECEPTOR", encoding="utf-8")

        stage_info = {
            "name": "smoke_stage",
            "exhaustiveness": 8,
            "num_modes": 9,
            "verbosity": 0,
        }
        conf_path, out_path = emit_vina_config(
            cfg,
            paths.pdb_id,
            str(receptor_path),
            (0.0, 0.0, 0.0),
            (20.0, 20.0, 20.0),
            str(lig_path),
            stage_info["name"],
            stage_info,
            1,
            smoke_log,
            variant="APO",
            ph_token=ph_label,
            legacy=False,
        )
        smoke_log.info("[smoke.emit.done] config=%s out=%s", conf_path, out_path)
    except Exception as exc:
        smoke_log.warning("[smoke.emit.skip] reason=%s", exc)
    finally:
        if old_run_env is None:
            os.environ.pop("ATLAS_RUN_ID", None)
        else:
            os.environ["ATLAS_RUN_ID"] = old_run_env
        if old_variant is None:
            os.environ.pop("APO_HOLO_VARIANT", None)
        else:
            os.environ["APO_HOLO_VARIANT"] = old_variant


# -------------------------
# Score I/O
# -------------------------
