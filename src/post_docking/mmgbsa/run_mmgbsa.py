from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path
from typing import Dict, Optional, Sequence

from installation import load_config


_COMPONENT = "mmgbsa.run"


def _get_logger() -> logging.Logger:
    logger = logging.getLogger(_COMPONENT)
    if logger.handlers:
        return logger

    logger.setLevel(logging.DEBUG)
    stream = logging.StreamHandler()
    stream.setLevel(logging.INFO)
    stream.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(stream)
    logger.propagate = False
    return logger


def _log(
    logger: logging.Logger, level: str, kvs: Dict[str, object], msg: str | None = None
) -> None:
    level_map = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
    }
    parts = [f"[{level}]", f"[{_COMPONENT}]"]
    if kvs:
        parts.append(" ".join(f"{k}={v}" for k, v in kvs.items()))
    if msg:
        parts.append(f"| {msg}")
    line = " ".join(parts)
    logger.log(level_map[level], line)


def _tmp_path(final_path: Path) -> Path:
    return final_path.with_name(f"{final_path.name}.tmp.{uuid.uuid4().hex}")


def _write_text_atomic(path: Path, text: str) -> None:
    tmp_path = _tmp_path(path)
    with tmp_path.open("w", encoding="utf-8") as handle:
        handle.write(text)
    if not tmp_path.exists() or tmp_path.stat().st_size == 0:
        if tmp_path.exists():
            tmp_path.unlink()
        raise RuntimeError(f"output not created: {tmp_path}")
    # Atomic replace avoids partial files on failures.
    os.replace(tmp_path, path)


def _to_bool(val: object, default: bool = False) -> bool:
    if isinstance(val, bool):
        return val
    if val is None:
        return default
    text = str(val).strip().lower()
    if text in {"1", "true", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "no", "n", "off"}:
        return False
    return default


def _to_int(val: object, default: int) -> int:
    try:
        return int(val)
    except Exception:
        return default


def _to_float(val: object, default: float) -> float:
    try:
        return float(val)
    except Exception:
        return default


def _cfg_get(cfg: object | None, key: str, default: object) -> object:
    if cfg is None:
        return default
    if isinstance(cfg, dict):
        return cfg.get(key, default)
    try:
        return cfg.get(key, default)
    except Exception:
        return getattr(cfg, key, default)


def _cfg_first(cfg: object | None, keys: Sequence[str], default: object) -> object:
    for key in keys:
        val = _cfg_get(cfg, key, None)
        if val is None:
            continue
        if str(val).strip() == "":
            continue
        return val
    return default


def _cfg_has_value(cfg: object | None, key: str) -> bool:
    if cfg is None:
        return False
    if isinstance(cfg, dict):
        if key not in cfg:
            return False
        return str(cfg.get(key)).strip() != ""
    try:
        val = cfg.get(key, None)
    except Exception:
        val = getattr(cfg, key, None)
    if val is None:
        return False
    return str(val).strip() != ""


def _mmgbsa_effective_md_enabled(cfg: object | None) -> bool:
    if _cfg_has_value(cfg, "MMGBSA_MD_ENABLED"):
        return _to_bool(_cfg_get(cfg, "MMGBSA_MD_ENABLED", False), default=False)
    traj_mode = str(_cfg_get(cfg, "MMGBSA_TRAJ_MODE", "") or "").strip().upper()
    if traj_mode == "IMPLICIT_MD":
        return True
    return _to_bool(_cfg_get(cfg, "MMGBSA_MD_RUN", False), default=False)


def parse_mmpbsa_delta_total(csv_path: Path) -> Optional[float]:
    """
    Lightweight parser for MMPBSA FINAL_RESULTS_MMPBSA.csv.
    Looks for the DELTA Energy Terms block and returns the first DELTA TOTAL value.
    """
    if not csv_path.exists():
        return None
    try:
        lines = csv_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return None

    in_delta = False
    header: Optional[list[str]] = None
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("DELTA Energy Terms"):
            in_delta = True
            header = None
            continue
        if not in_delta:
            continue
        if header is None:
            header = [t.strip() for t in stripped.split(",")]
            continue
        values = [t.strip() for t in stripped.split(",")]
        if not header:
            break
        try:
            idx = header.index("DELTA TOTAL")
        except ValueError:
            break
        if idx >= len(values):
            continue
        try:
            return float(values[idx])
        except Exception:
            return None
    return None


def write_aggregated_mmpbsa_results(
    work_dir: Path,
    mean_score: float,
    std_score: float,
    cfg: object,
    force: bool = False,
) -> Dict[str, str]:
    """
    Write aggregate FINAL_RESULTS files in a minimal MMPBSA-compatible format so downstream
    readers (frame aggregators) can parse DELTA TOTAL.
    """
    out_dat_name = str(
        _cfg_get(cfg, "MMGBSA_MMPBSA_OUT_DAT", "FINAL_RESULTS_MMPBSA.dat")
        or "FINAL_RESULTS_MMPBSA.dat"
    )
    out_csv_name = str(
        _cfg_get(cfg, "MMGBSA_MMPBSA_OUT_CSV", "FINAL_RESULTS_MMPBSA.csv")
        or "FINAL_RESULTS_MMPBSA.csv"
    )

    work_dir.mkdir(parents=True, exist_ok=True)
    out_dat_path = work_dir / out_dat_name
    out_csv_path = work_dir / out_csv_name

    if not force and out_dat_path.exists() and out_csv_path.exists():
        return {"out_dat": str(out_dat_path), "out_csv": str(out_csv_path)}

    dat_lines = [
        "FINAL RESULTS",
        f"DELTA TOTAL       {mean_score:.6f}   (stddev {std_score:.6f})",
        "",
    ]
    _write_text_atomic(out_dat_path, "\n".join(dat_lines))

    csv_lines = [
        "DELTA Energy Terms",
        "VDWAALS,EELEC,EGB,ESURF,EPB,ECAVITY,DELTA TOTAL,STDDEV",
        f"0,0,0,0,0,0,{mean_score:.6f},{std_score:.6f}",
    ]
    _write_text_atomic(out_csv_path, "\n".join(csv_lines) + "\n")

    return {"out_dat": str(out_dat_path), "out_csv": str(out_csv_path)}


def _looks_like_pkgs_cache(prefix: Path) -> bool:
    return "/micromamba/pkgs" in prefix.as_posix()


def _prefix_has_tool(prefix: Path, tool: str) -> bool:
    return (prefix / "bin" / tool).is_file()


def _resolve_micromamba() -> Optional[Path]:
    env_mm = os.environ.get("MICROMAMBA_EXE")
    if env_mm:
        env_path = Path(env_mm).expanduser()
        if env_path.is_file():
            return env_path
    found = shutil.which("micromamba")
    if found:
        return Path(found)
    home_candidate = Path.home() / "micromamba" / "bin" / "micromamba"
    if home_candidate.is_file():
        return home_candidate
    return None


def _common_ambertools_prefixes(cfg: object | None) -> list[Path]:
    candidates: list[Path] = []
    env_candidates = [
        os.environ.get("MMGBSA_AMBERTOOLS_PREFIX"),
        os.environ.get("AMBERTOOLS_PREFIX"),
        os.environ.get("CONDA_PREFIX"),
    ]
    for raw in env_candidates:
        if raw:
            candidates.append(Path(raw).expanduser())

    overall = _cfg_get(cfg, "OVERALL_DIR", None)
    if overall:
        base = Path(str(overall)).expanduser().resolve().parent
        candidates.append(base / "tools" / "envs" / "ambertools")

    mamba_root = os.environ.get("MAMBA_ROOT_PREFIX")
    if mamba_root:
        candidates.append(Path(mamba_root).expanduser() / "envs" / "AmberTools25")

    home = Path.home()
    candidates.extend(
        [
            home / "tools" / "envs" / "ambertools",
            home / "micromamba" / "envs" / "AmberTools25",
            home / "miniconda3" / "envs" / "AmberTools25",
            home / "anaconda3" / "envs" / "AmberTools25",
        ]
    )

    seen: set[str] = set()
    out: list[Path] = []
    for path in candidates:
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        out.append(path)
    return out


def _resolve_ambertools_prefix(cfg: object | None) -> Optional[str]:
    env_prefix = os.environ.get("AMBERTOOLS_PREFIX")
    if env_prefix:
        return env_prefix
    cfg_prefix = _cfg_get(cfg, "MMGBSA_AMBERTOOLS_PREFIX", None)
    if cfg_prefix:
        return str(cfg_prefix)
    fallback = _cfg_get(cfg, "AMBERTOOLS_PREFIX", None)
    if fallback:
        return str(fallback)
    return None


def build_mmpbsa_input(
    startframe: int,
    endframe: int,
    interval: int,
    verbose: int,
    igb: int,
    saltcon: float,
) -> str:
    saltcon_text = f"{float(saltcon):.3f}"
    lines = [
        "&general",
        f" startframe={int(startframe)}, endframe={int(endframe)}, interval={int(interval)},",
        f" verbose={int(verbose)},",
        "/",
        "&gb",
        f" igb={int(igb)}, saltcon={saltcon_text},",
        "/",
        "",
    ]
    return "\n".join(lines)


def write_mmpbsa_input(text: str, out_path: str, force: bool = False) -> str:
    path = Path(out_path)
    if path.exists() and path.stat().st_size > 0 and not force:
        return str(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    _write_text_atomic(path, text)
    return str(path)


def resolve_mmpbsa_runner(cfg: object) -> Dict[str, object]:
    logger = _get_logger()
    prefix_value = _resolve_ambertools_prefix(cfg)
    if prefix_value:
        prefix = Path(prefix_value).expanduser()
        if _looks_like_pkgs_cache(prefix):
            _log(
                logger,
                "WARNING",
                {"reason": "pkgs_cache_prefix", "fallback": "PATH_or_common_prefix"},
                "ignoring_AMBERTOOLS_PREFIX",
            )
        elif _prefix_has_tool(prefix, "MMPBSA.py"):
            micromamba = _resolve_micromamba()
            if not micromamba:
                raise FileNotFoundError(
                    "micromamba not found; cannot run MMPBSA.py from prefix"
                )
            runner = [str(micromamba), "run", "-p", str(prefix)]
            return {
                "runner": runner,
                "exe": "MMPBSA.py",
                "source": f"prefix:{prefix}",
            }
        else:
            _log(
                logger,
                "WARNING",
                {"reason": "missing_mmpbsa", "fallback": "PATH_or_common_prefix"},
                "ignoring_AMBERTOOLS_PREFIX",
            )

    mmpbsa_path = shutil.which("MMPBSA.py")
    if mmpbsa_path:
        return {
            "runner": [],
            "exe": mmpbsa_path,
            "source": "PATH",
        }

    common_prefixes = _common_ambertools_prefixes(cfg)
    for prefix in common_prefixes:
        if _prefix_has_tool(prefix, "MMPBSA.py") and not _looks_like_pkgs_cache(prefix):
            micromamba = _resolve_micromamba()
            if not micromamba:
                raise FileNotFoundError(
                    "micromamba not found; cannot run MMPBSA.py from prefix"
                )
            runner = [str(micromamba), "run", "-p", str(prefix)]
            return {
                "runner": runner,
                "exe": "MMPBSA.py",
                "source": f"prefix:{prefix}",
            }

    raise FileNotFoundError("could not locate MMPBSA.py; set AMBERTOOLS_PREFIX or PATH")


def _ensure_nonempty(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    if path.stat().st_size == 0:
        raise FileNotFoundError(f"{label} is empty: {path}")


def run_mmgbsa(
    complex_prmtop: str,
    receptor_prmtop: str,
    ligand_prmtop: str,
    trajectory_path: str,
    work_dir: str,
    cfg: object,
    force: bool = False,
    run: bool = True,
) -> dict:
    logger = _get_logger()
    work_path = Path(work_dir)

    enabled = _to_bool(_cfg_get(cfg, "MMGBSA_MMPBSA_ENABLED", True), default=True)
    if not enabled:
        _log(logger, "INFO", {"enabled": False}, "mmpbsa_disabled")
        return {
            "enabled": False,
            "run": False,
            "work_dir": str(work_path),
            "input_path": "",
            "log_path": "",
            "out_dat": "",
            "out_csv": "",
            "trajectory": "",
            "cmd_preview": "",
            "skipped": True,
            "returncode": None,
        }

    work_path.mkdir(parents=True, exist_ok=True)

    startframe = _to_int(
        _cfg_first(cfg, ["MMGBSA_GENERAL_STARTFRAME", "MMGBSA_MMPBSA_STARTFRAME"], 1), 1
    )
    endframe = _to_int(
        _cfg_first(cfg, ["MMGBSA_GENERAL_ENDFRAME", "MMGBSA_MMPBSA_ENDFRAME"], 1), 1
    )
    interval = _to_int(
        _cfg_first(cfg, ["MMGBSA_GENERAL_INTERVAL", "MMGBSA_MMPBSA_INTERVAL"], 1), 1
    )
    verbose = _to_int(
        _cfg_first(cfg, ["MMGBSA_GENERAL_VERBOSE", "MMGBSA_MMPBSA_VERBOSE"], 2), 2
    )
    igb = _to_int(_cfg_get(cfg, "MMGBSA_GB_IGB", 5), 5)
    saltcon = _to_float(_cfg_get(cfg, "MMGBSA_GB_SALTCON", 0.150), 0.150)

    if _cfg_has_value(cfg, "MMGBSA_MMPBSA_USE_TRAJ_FRAMES"):
        use_traj_frames = _to_bool(
            _cfg_get(cfg, "MMGBSA_MMPBSA_USE_TRAJ_FRAMES", True), default=True
        )
    else:
        use_traj_frames = _mmgbsa_effective_md_enabled(cfg)
    if use_traj_frames:
        startframe = 1
        endframe = 999999
        interval = 1

    input_name = str(
        _cfg_get(cfg, "MMGBSA_MMPBSA_INPUT_NAME", "mmpbsa.in") or "mmpbsa.in"
    )
    log_name = str(
        _cfg_get(cfg, "MMGBSA_MMPBSA_LOG_NAME", "mmpbsa.log") or "mmpbsa.log"
    )
    out_dat_name = str(
        _cfg_get(cfg, "MMGBSA_MMPBSA_OUT_DAT", "FINAL_RESULTS_MMPBSA.dat")
        or "FINAL_RESULTS_MMPBSA.dat"
    )
    out_csv_name = str(
        _cfg_get(cfg, "MMGBSA_MMPBSA_OUT_CSV", "FINAL_RESULTS_MMPBSA.csv")
        or "FINAL_RESULTS_MMPBSA.csv"
    )
    default_traj = str(_cfg_get(cfg, "MMGBSA_DEFAULT_TRAJ_NAME", "mdcrd") or "mdcrd")

    traj_value = trajectory_path or str(work_path / default_traj)
    traj_path = Path(traj_value)

    input_path = work_path / input_name
    log_path = work_path / log_name
    out_dat_path = work_path / out_dat_name
    out_csv_path = work_path / out_csv_name

    mmpbsa_text = build_mmpbsa_input(
        startframe=startframe,
        endframe=endframe,
        interval=interval,
        verbose=verbose,
        igb=igb,
        saltcon=saltcon,
    )
    write_mmpbsa_input(mmpbsa_text, str(input_path), force=force)

    out_dat_ok = out_dat_path.exists() and out_dat_path.stat().st_size > 0
    out_csv_ok = out_csv_path.exists() and out_csv_path.stat().st_size > 0
    outputs_ok = out_dat_ok and out_csv_ok

    run_flag = _to_bool(_cfg_get(cfg, "MMGBSA_MMPBSA_RUN", True), default=True)
    if not run:
        run_flag = False

    if outputs_ok and not force:
        run_flag = False

    cmd_preview = (
        "MMPBSA.py -O -i {input} -cp {cp} -rp {rp} -lp {lp} -y {traj} -o {out_dat} -eo {out_csv}"
    ).format(
        input=input_name,
        cp=Path(complex_prmtop).name,
        rp=Path(receptor_prmtop).name,
        lp=Path(ligand_prmtop).name,
        traj=traj_path.name,
        out_dat=out_dat_name,
        out_csv=out_csv_name,
    )

    _log(
        logger,
        "INFO",
        {
            "work_dir": work_path,
            "input": input_path,
            "out_dat": out_dat_path,
            "out_csv": out_csv_path,
            "traj": traj_path,
            "frames": f"{startframe}-{endframe}:{interval}",
            "use_traj_frames": use_traj_frames,
            "run": run_flag,
        },
        "mmpbsa_plan",
    )
    _log(logger, "INFO", {"cmd": cmd_preview}, "cmd")

    if not run_flag:
        return {
            "enabled": True,
            "run": False,
            "work_dir": str(work_path),
            "input_path": str(input_path),
            "log_path": str(log_path),
            "out_dat": str(out_dat_path),
            "out_csv": str(out_csv_path),
            "trajectory": str(traj_path),
            "cmd_preview": cmd_preview,
            "skipped": True,
            "returncode": None,
        }

    _ensure_nonempty(Path(complex_prmtop), "complex_prmtop")
    _ensure_nonempty(Path(receptor_prmtop), "receptor_prmtop")
    _ensure_nonempty(Path(ligand_prmtop), "ligand_prmtop")
    _ensure_nonempty(traj_path, "trajectory")

    runner_info = resolve_mmpbsa_runner(cfg)
    runner = list(runner_info["runner"])
    exe = str(runner_info["exe"])

    cp_rel = os.path.relpath(Path(complex_prmtop), work_path)
    rp_rel = os.path.relpath(Path(receptor_prmtop), work_path)
    lp_rel = os.path.relpath(Path(ligand_prmtop), work_path)
    traj_rel = os.path.relpath(traj_path, work_path)

    cmd = runner + [
        exe,
        "-O",
        "-i",
        input_name,
        "-cp",
        cp_rel,
        "-rp",
        rp_rel,
        "-lp",
        lp_rel,
        "-y",
        traj_rel,
        "-o",
        out_dat_name,
        "-eo",
        out_csv_name,
    ]

    _log(
        logger,
        "INFO",
        {"cmd": " ".join(cmd), "source": runner_info["source"]},
        "run",
    )

    with log_path.open("w", encoding="utf-8") as handle:
        proc = subprocess.run(
            cmd, cwd=str(work_path), stdout=handle, stderr=subprocess.STDOUT
        )

    result = {
        "enabled": True,
        "run": True,
        "work_dir": str(work_path),
        "input_path": str(input_path),
        "log_path": str(log_path),
        "out_dat": str(out_dat_path),
        "out_csv": str(out_csv_path),
        "trajectory": str(traj_path),
        "cmd_preview": " ".join(cmd),
        "skipped": False,
        "returncode": proc.returncode,
    }

    if proc.returncode != 0:
        _log(
            logger,
            "ERROR",
            {"work_dir": work_path, "rc": proc.returncode, "input": input_path},
            "mmpbsa_failed",
        )
        raise RuntimeError(f"MMPBSA.py failed (code={proc.returncode}); see {log_path}")

    _log(
        logger,
        "INFO",
        {
            "work_dir": work_path,
            "ok": "true",
            "rc": proc.returncode,
            "input": input_path,
            "out_dat": out_dat_path,
            "out_csv": out_csv_path,
            "traj": traj_path,
        },
        "done",
    )

    return result


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run MMGBSA using AmberTools MMPBSA.py."
    )
    parser.add_argument(
        "--work-dir",
        required=True,
        help="Working directory for mmpbsa input/log/output.",
    )
    parser.add_argument("--complex-prmtop", required=True, help="Complex prmtop path.")
    parser.add_argument(
        "--receptor-prmtop", required=True, help="Receptor prmtop path."
    )
    parser.add_argument("--ligand-prmtop", required=True, help="Ligand prmtop path.")
    parser.add_argument(
        "--traj",
        default=None,
        help="Trajectory path (defaults to work_dir/MMGBSA_DEFAULT_TRAJ_NAME).",
    )
    parser.add_argument(
        "--force", action="store_true", help="Overwrite outputs if present."
    )
    parser.add_argument(
        "--no-run", action="store_true", help="Write mmpbsa input only."
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_arg_parser()
    args = parser.parse_args(argv)
    logger = _get_logger()

    try:
        cfg = load_config()
        default_traj = str(
            _cfg_get(cfg, "MMGBSA_DEFAULT_TRAJ_NAME", "mdcrd") or "mdcrd"
        )
        traj_path = args.traj or str(Path(args.work_dir) / default_traj)
        run_mmgbsa(
            complex_prmtop=args.complex_prmtop,
            receptor_prmtop=args.receptor_prmtop,
            ligand_prmtop=args.ligand_prmtop,
            trajectory_path=traj_path,
            work_dir=args.work_dir,
            cfg=cfg,
            force=args.force,
            run=not args.no_run,
        )
    except Exception as exc:
        _log(
            logger,
            "ERROR",
            {"reason": type(exc).__name__, "detail": str(exc)},
            "failed",
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
