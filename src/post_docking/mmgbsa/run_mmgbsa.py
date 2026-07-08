from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional, Sequence, cast

from config.value_access import cfg_first, cfg_get, cfg_has_value, to_bool, to_float, to_int
from config.tool_resolver import select_ambertools_tool
from config.runtime_config import load_config
from post_docking.mmgbsa._atomic_io import write_text_atomic
from post_docking.mmgbsa.mmgbsa_mpi_env import (
    apply_blas_single_thread_env_defaults,
    resolve_mmpbsa_mpi_ranks,
    select_mpi_launcher,
)
from post_docking.mmgbsa.operation_exceptions import (
    FLOAT_COERCE_ERRORS,
    MMGBSA_CLI_ERRORS,
    TEXT_FILE_READ_ERRORS,
)


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


def _write_text_atomic(path: Path, text: str) -> None:
    write_text_atomic(path, text, require_nonempty=True)


def _to_bool(val: object, default: bool = False) -> bool:
    return bool(to_bool(val, default=default))


def _to_int(val: object, default: int) -> int:
    coerced = to_int(val, default=default)
    return int(default if coerced is None else coerced)


def _to_float(val: object, default: float) -> float:
    coerced = to_float(val, default=default)
    return float(default if coerced is None else coerced)


def _cfg_get(cfg: object | None, key: str, default: object) -> object:
    return cfg_get(cfg, key, default)


def _cfg_first(cfg: object | None, keys: Sequence[str], default: object) -> object:
    return cfg_first(cfg, keys, default)


def _cfg_has_value(cfg: object | None, key: str) -> bool:
    return cfg_has_value(cfg, key)


def _mmgbsa_effective_md_enabled(cfg: object | None) -> bool:
    if _cfg_has_value(cfg, "MMGBSA_MD_ENABLED"):
        return _to_bool(_cfg_get(cfg, "MMGBSA_MD_ENABLED", False), default=False)
    traj_mode = str(_cfg_get(cfg, "MMGBSA_TRAJ_MODE", "") or "").strip().upper()
    if traj_mode == "IMPLICIT_MD":
        return True
    return _to_bool(_cfg_get(cfg, "MMGBSA_MD_RUN", False), default=False)


def _analysis_frame_window(
    cfg: object | None, startframe: int, endframe: int, interval: int
) -> tuple[int, int, int, dict[str, object]]:
    frame_stride_ps = _to_float(_cfg_get(cfg, "MMGBSA_MD_FRAME_STRIDE_PS", 0.0), 0.0)
    start_ps = _to_float(_cfg_get(cfg, "MMGBSA_ANALYSIS_START_PS", 0.0), 0.0)
    end_ps = _to_float(_cfg_get(cfg, "MMGBSA_ANALYSIS_END_PS", 0.0), 0.0)
    explicit_start = _to_int(_cfg_get(cfg, "MMGBSA_ANALYSIS_STARTFRAME", 0), 0)
    explicit_end = _to_int(_cfg_get(cfg, "MMGBSA_ANALYSIS_ENDFRAME", 0), 0)
    explicit_interval = _to_int(
        _cfg_get(cfg, "MMGBSA_ANALYSIS_INTERVAL", interval), interval
    )

    if frame_stride_ps > 0 and start_ps > 0:
        startframe = max(1, int(start_ps // frame_stride_ps) + 1)
    if frame_stride_ps > 0 and end_ps > 0:
        endframe = max(startframe, int(end_ps // frame_stride_ps))
    if explicit_start > 0:
        startframe = explicit_start
    if explicit_end > 0:
        endframe = explicit_end
    interval = max(1, explicit_interval)
    meta: dict[str, object] = {
        "analysis_start_ps": start_ps,
        "analysis_end_ps": end_ps,
        "frame_stride_ps": frame_stride_ps,
        "analysis_startframe": startframe,
        "analysis_endframe": endframe,
        "analysis_interval": interval,
    }
    return startframe, endframe, interval, meta


def _read_mmpbsa_csv_lines(csv_path: Path) -> list[str]:
    try:
        return csv_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except TEXT_FILE_READ_ERRORS as exc:
        logging.getLogger(_COMPONENT).debug(
            "[mmgbsa.parse_delta] action=csv_read_failed path=%s err=%s",
            csv_path,
            exc,
            exc_info=True,
        )
    return []


def _parse_delta_total_from_lines(lines: Sequence[str]) -> Optional[float]:
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
            header = [token.strip() for token in stripped.split(",")]
            continue
        return _parse_delta_total_row(header, stripped)
    return None


def _parse_delta_total_row(header: Sequence[str], row: str) -> Optional[float]:
    try:
        idx = header.index("DELTA TOTAL")
    except ValueError:
        return None
    values = [token.strip() for token in row.split(",")]
    if idx >= len(values):
        return None
    try:
        return float(values[idx])
    except FLOAT_COERCE_ERRORS:
        return None


def parse_mmpbsa_delta_total(csv_path: Path) -> Optional[float]:
    """
    Lightweight parser for MMPBSA FINAL_RESULTS_MMPBSA.csv.
    Looks for the DELTA Energy Terms block and returns the first DELTA TOTAL value.
    """
    if not csv_path.exists():
        return None
    return _parse_delta_total_from_lines(_read_mmpbsa_csv_lines(csv_path))


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


def resolve_mmpbsa_runner(cfg: object, *, prefer_mpi: bool = True) -> Dict[str, object]:
    mpi_enabled = prefer_mpi and _to_bool(
        _cfg_get(cfg, "MMGBSA_MMPBSA_MPI_ENABLED", True), True
    )
    if mpi_enabled:
        try:
            selection = select_ambertools_tool(
                cfg,
                "MMPBSA.py.MPI",
                prefix_exe_style="absolute",
            )
            launcher, launcher_source = select_mpi_launcher(
                selection.prefix,
                probe_common_prefixes=True,
                cfg_for_common_prefixes=None,
            )
            if launcher:
                ranks = resolve_mmpbsa_mpi_ranks(cfg)
                runner = list(selection.runner) + [launcher, "-np", str(ranks)]
                return {
                    "runner": runner,
                    "exe": selection.exe,
                    "source": f"{selection.source};mpi:{launcher_source};ranks:{ranks}",
                    "mpi": True,
                    "mpi_ranks": ranks,
                    "amberhome": str(selection.prefix) if selection.prefix else "",
                }
        except FileNotFoundError:
            pass
    selection = select_ambertools_tool(cfg, "MMPBSA.py")
    return {
        "runner": selection.runner,
        "exe": selection.exe,
        "source": selection.source,
        "mpi": False,
        "mpi_ranks": 1,
        "amberhome": str(selection.prefix) if selection.prefix else "",
    }


def _mmpbsa_mpi_import_failed(log_path: Path) -> bool:
    try:
        text = log_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return False
    markers = (
        "Could not import mpi4py package",
        "ModuleNotFoundError: No module named 'mpi4py'",
    )
    return any(marker in text for marker in markers)


def _ensure_nonempty(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    if path.stat().st_size == 0:
        raise FileNotFoundError(f"{label} is empty: {path}")


def _disabled_mmpbsa_result(work_path: Path) -> dict:
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


def _resolve_mmpbsa_frames(cfg: object) -> tuple[int, int, int, bool, dict[str, object]]:
    startframe = _to_int(
        _cfg_first(cfg, ["MMGBSA_GENERAL_STARTFRAME", "MMGBSA_MMPBSA_STARTFRAME"], 1), 1
    )
    endframe = _to_int(
        _cfg_first(cfg, ["MMGBSA_GENERAL_ENDFRAME", "MMGBSA_MMPBSA_ENDFRAME"], 1), 1
    )
    interval = _to_int(
        _cfg_first(cfg, ["MMGBSA_GENERAL_INTERVAL", "MMGBSA_MMPBSA_INTERVAL"], 1), 1
    )
    if _cfg_has_value(cfg, "MMGBSA_MMPBSA_USE_TRAJ_FRAMES"):
        use_traj_frames = _to_bool(
            _cfg_get(cfg, "MMGBSA_MMPBSA_USE_TRAJ_FRAMES", True), default=True
        )
    else:
        use_traj_frames = _mmgbsa_effective_md_enabled(cfg)
    if use_traj_frames:
        startframe, endframe, interval, analysis_meta = _analysis_frame_window(
            cfg, 1, 999999, 1
        )
        return startframe, endframe, interval, use_traj_frames, analysis_meta
    analysis_meta = {
        "analysis_startframe": startframe,
        "analysis_endframe": endframe,
        "analysis_interval": interval,
    }
    return startframe, endframe, interval, use_traj_frames, analysis_meta


def _mmpbsa_output_names(cfg: object) -> dict[str, str]:
    return {
        "input": str(_cfg_get(cfg, "MMGBSA_MMPBSA_INPUT_NAME", "mmpbsa.in") or "mmpbsa.in"),
        "log": str(_cfg_get(cfg, "MMGBSA_MMPBSA_LOG_NAME", "mmpbsa.log") or "mmpbsa.log"),
        "out_dat": str(
            _cfg_get(cfg, "MMGBSA_MMPBSA_OUT_DAT", "FINAL_RESULTS_MMPBSA.dat")
            or "FINAL_RESULTS_MMPBSA.dat"
        ),
        "out_csv": str(
            _cfg_get(cfg, "MMGBSA_MMPBSA_OUT_CSV", "FINAL_RESULTS_MMPBSA.csv")
            or "FINAL_RESULTS_MMPBSA.csv"
        ),
        "default_traj": str(_cfg_get(cfg, "MMGBSA_DEFAULT_TRAJ_NAME", "mdcrd") or "mdcrd"),
    }


def _mmpbsa_command(
    runner_info: Dict[str, object],
    *,
    input_name: str,
    cp_rel: str,
    rp_rel: str,
    lp_rel: str,
    traj_rel: str,
    out_dat_name: str,
    out_csv_name: str,
) -> list[str]:
    runner = list(cast(Sequence[str], runner_info["runner"]))
    return runner + [
        str(runner_info["exe"]),
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


def _mmpbsa_env(runner_info: Dict[str, object]) -> dict[str, str]:
    env = os.environ.copy()
    apply_blas_single_thread_env_defaults(env)
    amberhome = str(runner_info.get("amberhome", "") or "").strip()
    if amberhome:
        env.setdefault("AMBERHOME", amberhome)
    return env


def _run_command_to_log(
    cmd: Sequence[str],
    *,
    work_path: Path,
    log_path: Path,
    env: dict[str, str],
    append: bool = False,
) -> subprocess.CompletedProcess[str]:
    mode = "a" if append else "w"
    with log_path.open(mode, encoding="utf-8") as handle:
        if append:
            handle.write("\n\n# Retrying with serial MMPBSA.py after MPI import failure.\n")
        return subprocess.run(
            list(cmd),
            cwd=str(work_path),
            stdout=handle,
            stderr=subprocess.STDOUT,
            env=env,
            text=True,
        )


def _run_mmpbsa_with_serial_fallback(
    *,
    cfg: object,
    runner_info: Dict[str, object],
    cmd: Sequence[str],
    work_path: Path,
    log_path: Path,
    logger: logging.Logger,
    command_parts: dict[str, str],
) -> tuple[subprocess.CompletedProcess[str], Dict[str, object], list[str]]:
    env = _mmpbsa_env(runner_info)
    proc = _run_command_to_log(cmd, work_path=work_path, log_path=log_path, env=env)
    if not (proc.returncode != 0 and runner_info.get("mpi") and _mmpbsa_mpi_import_failed(log_path)):
        return proc, runner_info, list(cmd)

    fallback_info = resolve_mmpbsa_runner(cfg, prefer_mpi=False)
    fallback_cmd = _mmpbsa_command(fallback_info, **command_parts)
    _log(
        logger,
        "WARNING",
        {
            "reason": "mpi4py_missing",
            "fallback": "serial",
            "cmd": " ".join(fallback_cmd),
        },
        "mmpbsa_mpi_fallback",
    )
    fallback_env = _mmpbsa_env(fallback_info)
    proc = _run_command_to_log(
        fallback_cmd,
        work_path=work_path,
        log_path=log_path,
        env=fallback_env,
        append=True,
    )
    return proc, fallback_info, fallback_cmd


def _should_run_mmpbsa(cfg: object, *, requested: bool, outputs_ok: bool, force: bool) -> bool:
    run_flag = _to_bool(_cfg_get(cfg, "MMGBSA_MMPBSA_RUN", True), default=True)
    return bool(run_flag and requested and (force or not outputs_ok))


def _skip_mmpbsa_result(
    *,
    work_path: Path,
    input_path: Path,
    log_path: Path,
    out_dat_path: Path,
    out_csv_path: Path,
    traj_path: Path,
    cmd_preview: str,
    analysis_meta: dict[str, object],
) -> dict:
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
        "analysis": analysis_meta,
    }


def _completed_mmpbsa_result(
    *,
    work_path: Path,
    input_path: Path,
    log_path: Path,
    out_dat_path: Path,
    out_csv_path: Path,
    traj_path: Path,
    cmd: Sequence[str],
    proc: subprocess.CompletedProcess[str],
    analysis_meta: dict[str, object],
    runner_info: Dict[str, object],
) -> dict:
    return {
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
        "analysis": analysis_meta,
        "mpi": runner_info.get("mpi", False),
        "mpi_ranks": runner_info.get("mpi_ranks", 1),
    }


def _log_mmpbsa_plan(
    logger: logging.Logger,
    *,
    work_path: Path,
    input_path: Path,
    out_dat_path: Path,
    out_csv_path: Path,
    traj_path: Path,
    startframe: int,
    endframe: int,
    interval: int,
    use_traj_frames: bool,
    run_flag: bool,
    analysis_meta: dict[str, object],
    cmd_preview: str,
) -> None:
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
            "analysis_start_ps": analysis_meta.get("analysis_start_ps", ""),
            "analysis_end_ps": analysis_meta.get("analysis_end_ps", ""),
        },
        "mmpbsa_plan",
    )
    _log(logger, "INFO", {"cmd": cmd_preview}, "cmd")


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
        return _disabled_mmpbsa_result(work_path)

    work_path.mkdir(parents=True, exist_ok=True)
    verbose = _to_int(
        _cfg_first(cfg, ["MMGBSA_GENERAL_VERBOSE", "MMGBSA_MMPBSA_VERBOSE"], 2), 2
    )
    igb = _to_int(_cfg_get(cfg, "MMGBSA_GB_IGB", 5), 5)
    saltcon = _to_float(_cfg_get(cfg, "MMGBSA_GB_SALTCON", 0.150), 0.150)
    startframe, endframe, interval, use_traj_frames, analysis_meta = _resolve_mmpbsa_frames(cfg)
    names = _mmpbsa_output_names(cfg)
    input_name = names["input"]
    out_dat_name = names["out_dat"]
    out_csv_name = names["out_csv"]

    traj_value = trajectory_path or str(work_path / names["default_traj"])
    traj_path = Path(traj_value)

    input_path = work_path / input_name
    log_path = work_path / names["log"]
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
    run_flag = _should_run_mmpbsa(
        cfg, requested=run, outputs_ok=out_dat_ok and out_csv_ok, force=force
    )

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

    _log_mmpbsa_plan(
        logger,
        work_path=work_path,
        input_path=input_path,
        out_dat_path=out_dat_path,
        out_csv_path=out_csv_path,
        traj_path=traj_path,
        startframe=startframe,
        endframe=endframe,
        interval=interval,
        use_traj_frames=use_traj_frames,
        run_flag=run_flag,
        analysis_meta=analysis_meta,
        cmd_preview=cmd_preview,
    )

    if not run_flag:
        return _skip_mmpbsa_result(
            work_path=work_path,
            input_path=input_path,
            log_path=log_path,
            out_dat_path=out_dat_path,
            out_csv_path=out_csv_path,
            traj_path=traj_path,
            cmd_preview=cmd_preview,
            analysis_meta=analysis_meta,
        )

    _ensure_nonempty(Path(complex_prmtop), "complex_prmtop")
    _ensure_nonempty(Path(receptor_prmtop), "receptor_prmtop")
    _ensure_nonempty(Path(ligand_prmtop), "ligand_prmtop")
    _ensure_nonempty(traj_path, "trajectory")

    runner_info = resolve_mmpbsa_runner(cfg, prefer_mpi=use_traj_frames)

    cp_rel = os.path.relpath(Path(complex_prmtop), work_path)
    rp_rel = os.path.relpath(Path(receptor_prmtop), work_path)
    lp_rel = os.path.relpath(Path(ligand_prmtop), work_path)
    traj_rel = os.path.relpath(traj_path, work_path)

    command_parts = {
        "input_name": input_name,
        "cp_rel": cp_rel,
        "rp_rel": rp_rel,
        "lp_rel": lp_rel,
        "traj_rel": traj_rel,
        "out_dat_name": out_dat_name,
        "out_csv_name": out_csv_name,
    }
    cmd = _mmpbsa_command(runner_info, **command_parts)

    _log(
        logger,
        "INFO",
        {
            "cmd": " ".join(cmd),
            "source": runner_info["source"],
            "mpi": runner_info.get("mpi", False),
            "mpi_ranks": runner_info.get("mpi_ranks", 1),
        },
        "run",
    )

    proc, runner_info, cmd = _run_mmpbsa_with_serial_fallback(
        cfg=cfg,
        runner_info=runner_info,
        cmd=cmd,
        work_path=work_path,
        log_path=log_path,
        logger=logger,
        command_parts=command_parts,
    )

    result = _completed_mmpbsa_result(
        work_path=work_path,
        input_path=input_path,
        log_path=log_path,
        out_dat_path=out_dat_path,
        out_csv_path=out_csv_path,
        traj_path=traj_path,
        cmd=cmd,
        proc=proc,
        analysis_meta=analysis_meta,
        runner_info=runner_info,
    )

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
    except MMGBSA_CLI_ERRORS as exc:
        _log(
            logger,
            "ERROR",
            {
                "reason": type(exc).__name__,
                "detail": str(exc),
                "action": "cli_run_failed",
            },
            "failed",
        )
        logging.getLogger(_COMPONENT).debug(
            "[mmgbsa.cli] action=run_failed", exc_info=True
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
