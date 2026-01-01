from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from installation import load_config


_COMPONENT = "mmgbsa.cpptraj"


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


def _log(logger: logging.Logger, level: str, kvs: Dict[str, object], msg: str | None = None) -> None:
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


def _looks_like_pkgs_cache(prefix: Path) -> bool:
    return "/micromamba/pkgs" in prefix.as_posix()


def _prefix_has_tool(prefix: Path, tool: str) -> bool:
    return (prefix / "bin" / tool).is_file()


def _resolve_micromamba() -> Optional[Path]:
    preferred = Path("/home/michael/atlas/micromamba/bin/micromamba")
    if preferred.is_file():
        return preferred
    found = shutil.which("micromamba")
    return Path(found) if found else None


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


def _select_cpptraj_runner(logger: logging.Logger, cfg: object | None) -> Tuple[Sequence[str], str, str]:
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
        elif _prefix_has_tool(prefix, "cpptraj"):
            micromamba = _resolve_micromamba()
            if not micromamba:
                raise FileNotFoundError("micromamba not found; cannot run cpptraj from prefix")
            runner = [str(micromamba), "run", "-p", str(prefix)]
            return runner, "cpptraj", f"prefix:{prefix}"
        else:
            _log(
                logger,
                "WARNING",
                {"reason": "missing_cpptraj", "fallback": "PATH_or_common_prefix"},
                "ignoring_AMBERTOOLS_PREFIX",
            )

    cpptraj_path = shutil.which("cpptraj")
    if cpptraj_path:
        return [], cpptraj_path, "PATH"

    common_prefixes = [
        Path("/home/atlas/micromamba/envs/AmberTools25"),
        Path("/home/michael/atlas/tools/envs/ambertools"),
        Path("/home/michael/atlas/micromamba/envs/AmberTools25"),
    ]
    for prefix in common_prefixes:
        if _prefix_has_tool(prefix, "cpptraj") and not _looks_like_pkgs_cache(prefix):
            micromamba = _resolve_micromamba()
            if not micromamba:
                raise FileNotFoundError("micromamba not found; cannot run cpptraj from prefix")
            runner = [str(micromamba), "run", "-p", str(prefix)]
            return runner, "cpptraj", f"prefix:{prefix}"

    raise FileNotFoundError("could not locate cpptraj; set AMBERTOOLS_PREFIX or PATH")


def _select_sander_runner(logger: logging.Logger, cfg: object | None) -> Tuple[Sequence[str], str, str]:
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
        elif _prefix_has_tool(prefix, "sander"):
            micromamba = _resolve_micromamba()
            if not micromamba:
                raise FileNotFoundError("micromamba not found; cannot run sander from prefix")
            runner = [str(micromamba), "run", "-p", str(prefix)]
            return runner, "sander", f"prefix:{prefix}"
        else:
            _log(
                logger,
                "WARNING",
                {"reason": "missing_sander", "fallback": "PATH_or_common_prefix"},
                "ignoring_AMBERTOOLS_PREFIX",
            )

    sander_path = shutil.which("sander")
    if sander_path:
        return [], sander_path, "PATH"

    common_prefixes = [
        Path("/home/atlas/micromamba/envs/AmberTools25"),
        Path("/home/michael/atlas/tools/envs/ambertools"),
        Path("/home/michael/atlas/micromamba/envs/AmberTools25"),
    ]
    for prefix in common_prefixes:
        if _prefix_has_tool(prefix, "sander") and not _looks_like_pkgs_cache(prefix):
            micromamba = _resolve_micromamba()
            if not micromamba:
                raise FileNotFoundError("micromamba not found; cannot run sander from prefix")
            runner = [str(micromamba), "run", "-p", str(prefix)]
            return runner, "sander", f"prefix:{prefix}"

    raise FileNotFoundError("could not locate sander; set AMBERTOOLS_PREFIX or PATH")


def _format_restraint_block(cfg: object | None) -> List[str]:
    restrain = _to_bool(_cfg_get(cfg, "MMGBSA_MD_RESTRAIN_PROTEIN_HEAVY", True), default=True)
    if not restrain:
        return []
    wt = _to_float(_cfg_get(cfg, "MMGBSA_MD_RESTRAINT_WT", 5.0), 5.0)
    mask = str(_cfg_get(cfg, "MMGBSA_MD_RESTRAINT_MASK", ":1-999999 & !@H=") or ":1-999999 & !@H=")
    return [
        " ntr=1,",
        f" restraint_wt={wt:.3f},",
        f" restraintmask='{mask}',",
    ]


def _compute_steps(ps: float, dt_ps: float) -> int:
    if dt_ps <= 0:
        dt_ps = 0.002
    steps = int(round(ps / dt_ps))
    return max(1, steps)


def build_sander_inputs(cfg: object | None) -> Dict[str, str]:
    dt_ps = _to_float(_cfg_get(cfg, "MMGBSA_MD_DT_PS", 0.002), 0.002)
    heat_ps = _to_float(_cfg_get(cfg, "MMGBSA_MD_HEAT_PS", 10.0), 10.0)
    equil_ps = _to_float(_cfg_get(cfg, "MMGBSA_MD_EQUIL_PS", 20.0), 20.0)
    prod_ps = _to_float(_cfg_get(cfg, "MMGBSA_MD_PROD_PS", 100.0), 100.0)
    temp0 = _to_float(_cfg_get(cfg, "MMGBSA_MD_TEMP0", 300.0), 300.0)
    ntt = _to_int(_cfg_get(cfg, "MMGBSA_MD_NTT", 3), 3)
    gamma_ln = _to_float(_cfg_get(cfg, "MMGBSA_MD_GAMMA_LN", 2.0), 2.0)
    ntc = _to_int(_cfg_get(cfg, "MMGBSA_MD_NTC", 2), 2)
    ntf = _to_int(_cfg_get(cfg, "MMGBSA_MD_NTF", 2), 2)
    igb = _to_int(_cfg_get(cfg, "MMGBSA_MD_IGB", 5), 5)
    saltcon = _to_float(_cfg_get(cfg, "MMGBSA_MD_SALTCON", 0.150), 0.150)
    frame_stride_ps = _to_float(_cfg_get(cfg, "MMGBSA_MD_FRAME_STRIDE_PS", 2.0), 2.0)
    traj_format = str(_cfg_get(cfg, "MMGBSA_MD_TRAJ_FORMAT", "nc") or "nc").strip().lower()

    heat_steps = _compute_steps(heat_ps, dt_ps)
    equil_steps = _compute_steps(equil_ps, dt_ps)
    prod_steps = _compute_steps(prod_ps, dt_ps)
    ntwx = max(1, int(round(frame_stride_ps / dt_ps))) if frame_stride_ps > 0 else 1

    restraint_lines = _format_restraint_block(cfg)
    saltcon_text = f"{saltcon:.3f}"

    min_lines = [
        "min",
        "&cntrl",
        " imin=1,",
        " maxcyc=500,",
        " ncyc=250,",
        " ntb=0,",
        " cut=999.0,",
        f" igb={igb},",
        f" saltcon={saltcon_text},",
    ]
    min_lines.extend(restraint_lines)
    min_lines.append("/")

    heat_lines = [
        "heat",
        "&cntrl",
        " imin=0,",
        " irest=0,",
        " ntx=1,",
        f" nstlim={heat_steps},",
        f" dt={dt_ps:.6f},",
        " tempi=0.0,",
        f" temp0={temp0:.2f},",
        " ntb=0,",
        f" ntt={ntt},",
        f" gamma_ln={gamma_ln:.3f},",
        f" ntc={ntc},",
        f" ntf={ntf},",
        " ntpr=100,",
        f" igb={igb},",
        f" saltcon={saltcon_text},",
    ]
    heat_lines.extend(restraint_lines)
    heat_lines.append(" nmropt=1,")
    heat_lines.append("/")
    heat_lines.append(
        f"&wt type='TEMP0', istep1=0, istep2={heat_steps}, value1=0.0, value2={temp0:.2f} /"
    )
    heat_lines.append("&wt type='END' /")

    equil_lines = [
        "equil",
        "&cntrl",
        " imin=0,",
        " irest=1,",
        " ntx=5,",
        f" nstlim={equil_steps},",
        f" dt={dt_ps:.6f},",
        f" temp0={temp0:.2f},",
        " ntb=0,",
        f" ntt={ntt},",
        f" gamma_ln={gamma_ln:.3f},",
        f" ntc={ntc},",
        f" ntf={ntf},",
        " ntpr=100,",
        f" igb={igb},",
        f" saltcon={saltcon_text},",
    ]
    equil_lines.extend(restraint_lines)
    equil_lines.append("/")

    prod_lines = [
        "prod",
        "&cntrl",
        " imin=0,",
        " irest=1,",
        " ntx=5,",
        f" nstlim={prod_steps},",
        f" dt={dt_ps:.6f},",
        f" temp0={temp0:.2f},",
        " ntb=0,",
        f" ntt={ntt},",
        f" gamma_ln={gamma_ln:.3f},",
        f" ntc={ntc},",
        f" ntf={ntf},",
        " ntpr=100,",
        f" ntwx={ntwx},",
        f" ntwr={ntwx},",
        f" igb={igb},",
        f" saltcon={saltcon_text},",
    ]
    if traj_format == "nc":
        prod_lines.append(" ioutfm=1,")
    prod_lines.extend(restraint_lines)
    prod_lines.append("/")

    return {
        "min": "\n".join(min_lines) + "\n",
        "heat": "\n".join(heat_lines) + "\n",
        "equil": "\n".join(equil_lines) + "\n",
        "prod": "\n".join(prod_lines) + "\n",
    }


def build_cpptraj_input(
    prmtop_path: str,
    trajin_path: str,
    trajin_format: str,
    startframe: int,
    endframe: int,
    interval: int,
    trajout_path: str,
    trajout_format: str,
) -> str:
    trajin_line = f"trajin {trajin_path} {startframe} {endframe}"
    if interval != 1:
        trajin_line = f"{trajin_line} {interval}"
    if trajin_format and trajin_format.lower() not in {"inpcrd", ""}:
        trajin_line = f"{trajin_line} {trajin_format}"

    trajout_line = f"trajout {trajout_path}"
    if trajout_format and trajout_format.lower() not in {"mdcrd", ""}:
        trajout_line = f"{trajout_line} {trajout_format}"

    lines = [
        f"parm {prmtop_path}",
        trajin_line,
        trajout_line,
        "run",
        "quit",
        "",
    ]
    return "\n".join(lines)


def write_cpptraj_file(text: str, out_path: str, force: bool = False) -> str:
    path = Path(out_path)
    if path.exists() and path.stat().st_size > 0 and not force:
        return str(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    _write_text_atomic(path, text)
    return str(path)


def _infer_trajout_path(cpptraj_in: Path, work_dir: Path) -> Optional[Path]:
    try:
        lines = cpptraj_in.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return None
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.lower().startswith("trajout "):
            tokens = stripped.split()
            if len(tokens) >= 2:
                path = Path(tokens[1])
                if not path.is_absolute():
                    return work_dir / path
                return path
    return None


def _run_cpptraj_impl(
    cpptraj_in: str,
    work_dir: str,
    runner: Optional[Callable[[Sequence[str], str, str], object]] = None,
) -> dict:
    logger = _get_logger()
    cpptraj_path = Path(cpptraj_in)
    work_path = Path(work_dir)

    if not cpptraj_path.exists():
        raise FileNotFoundError(f"cpptraj input not found: {cpptraj_path}")

    work_path.mkdir(parents=True, exist_ok=True)

    cfg = load_config()
    cmd_prefix, cpptraj_bin, source = _select_cpptraj_runner(logger, cfg)
    cpptraj_arg = cpptraj_path.name if cpptraj_path.parent == work_path else str(cpptraj_path)
    cmd = list(cmd_prefix) + [cpptraj_bin, "-i", cpptraj_arg]
    log_path = work_path / "cpptraj.log"

    _log(
        logger,
        "INFO",
        {"cpptraj": " ".join(cmd), "source": source, "work_dir": work_path},
        "run",
    )

    result: Dict[str, object] = {
        "args": " ".join(cmd),
        "returncode": 0,
        "log_path": str(log_path),
        "ok": True,
    }

    if runner is not None:
        runner_result = runner(cmd, str(work_path), str(log_path))
        if isinstance(runner_result, dict):
            result.update(runner_result)
        elif hasattr(runner_result, "returncode"):
            result["returncode"] = int(getattr(runner_result, "returncode"))
            result["ok"] = result["returncode"] == 0
        elif isinstance(runner_result, int):
            result["returncode"] = runner_result
            result["ok"] = runner_result == 0
    else:
        with log_path.open("w", encoding="utf-8") as handle:
            proc = subprocess.run(cmd, cwd=str(work_path), stdout=handle, stderr=subprocess.STDOUT)
        result["returncode"] = proc.returncode
        result["ok"] = proc.returncode == 0

    trajout_path = _infer_trajout_path(cpptraj_path, work_path)
    result["trajout_path"] = str(trajout_path) if trajout_path else ""

    if not result["ok"]:
        _log(logger, "ERROR", {"returncode": result["returncode"], "log": log_path}, "cpptraj_failed")
        raise RuntimeError(f"cpptraj failed (code={result['returncode']}); see {log_path}")

    _log(logger, "INFO", {"ok": "true", "log": log_path}, "cpptraj_done")
    return result


def run_cpptraj(
    cpptraj_in: str,
    work_dir: str,
    runner: Optional[Callable[[Sequence[str], str, str], object]] = None,
) -> dict:
    return _run_cpptraj_impl(cpptraj_in=cpptraj_in, work_dir=work_dir, runner=runner)


def _inject_seed(text: str, seed: int) -> str:
    lines = text.splitlines()
    for idx, line in enumerate(lines):
        marker = line.strip().lower()
        if marker in {"/", "&end"}:
            lines.insert(idx, f" ig={int(seed)},")
            break
    return "\n".join(lines) + "\n"


def run_implicit_md(
    complex_prmtop: str,
    complex_inpcrd: str,
    out_dir: str,
    cfg: object,
    replicate_index: int,
    seed: int,
    force: bool = False,
    run: bool = True,
) -> dict:
    logger = _get_logger()
    prmtop_path = Path(complex_prmtop)
    inpcrd_path = Path(complex_inpcrd)
    base_dir = Path(out_dir)

    engine = str(_cfg_get(cfg, "MMGBSA_MD_ENGINE", "sander") or "sander").strip().lower()
    if engine != "sander":
        raise ValueError(f"MMGBSA_MD_ENGINE must be sander (got {engine})")

    if not prmtop_path.exists():
        raise FileNotFoundError(f"prmtop not found: {prmtop_path}")
    if not inpcrd_path.exists():
        raise FileNotFoundError(f"inpcrd not found: {inpcrd_path}")

    rep_dir = base_dir / f"rep{replicate_index}"
    md_dir = rep_dir / "md"
    md_dir.mkdir(parents=True, exist_ok=True)

    traj_name = str(_cfg_get(cfg, "MMGBSA_MD_TRAJ_NAME", "prod.nc") or "prod.nc")
    traj_format = str(_cfg_get(cfg, "MMGBSA_MD_TRAJ_FORMAT", "nc") or "nc").strip().lower()
    traj_path = md_dir / traj_name

    dt_ps = _to_float(_cfg_get(cfg, "MMGBSA_MD_DT_PS", 0.002), 0.002)
    prod_ps = _to_float(_cfg_get(cfg, "MMGBSA_MD_PROD_PS", 100.0), 100.0)
    frame_stride_ps = _to_float(_cfg_get(cfg, "MMGBSA_MD_FRAME_STRIDE_PS", 2.0), 2.0)
    prod_steps = _compute_steps(prod_ps, dt_ps)
    ntwx = max(1, int(round(frame_stride_ps / dt_ps))) if frame_stride_ps > 0 else 1
    n_frames_est = max(1, int(prod_steps // ntwx))

    inputs = build_sander_inputs(cfg)
    inputs["heat"] = _inject_seed(inputs["heat"], seed)

    for name, text in inputs.items():
        _write_text_atomic(md_dir / f"{name}.in", text)

    traj_exists = traj_path.exists() and traj_path.stat().st_size > 0
    if traj_exists and not force:
        _log(
            logger,
            "INFO",
            {"traj": traj_path, "replicate": replicate_index, "skipped": True},
            "md_exists",
        )
        return {
            "ok": True,
            "run": False,
            "skipped": True,
            "traj_path": str(traj_path),
            "traj_format": traj_format,
            "n_frames_est": n_frames_est,
            "replicate": replicate_index,
            "seed": seed,
            "md_dir": str(md_dir),
            "replicate_dir": str(rep_dir),
        }

    if not run:
        _log(
            logger,
            "INFO",
            {"traj": traj_path, "replicate": replicate_index, "run": False},
            "md_planned",
        )
        return {
            "ok": False,
            "run": False,
            "skipped": True,
            "traj_path": str(traj_path),
            "traj_format": traj_format,
            "n_frames_est": n_frames_est,
            "replicate": replicate_index,
            "seed": seed,
            "md_dir": str(md_dir),
            "replicate_dir": str(rep_dir),
        }

    cmd_prefix, sander_bin, source = _select_sander_runner(logger, cfg)
    restrain = _to_bool(_cfg_get(cfg, "MMGBSA_MD_RESTRAIN_PROTEIN_HEAVY", True), default=True)
    prmtop_rel = os.path.relpath(prmtop_path, md_dir)
    inpcrd_rel = os.path.relpath(inpcrd_path, md_dir)

    steps = [
        ("min", inpcrd_rel, "min.rst7"),
        ("heat", "min.rst7", "heat.rst7"),
        ("equil", "heat.rst7", "equil.rst7"),
        ("prod", "equil.rst7", "prod.rst7"),
    ]

    for step_name, coord_in, coord_out in steps:
        out_file = f"{step_name}.out"
        cmd = list(cmd_prefix) + [
            sander_bin,
            "-O",
            "-i",
            f"{step_name}.in",
            "-p",
            prmtop_rel,
            "-c",
            coord_in,
            "-o",
            out_file,
            "-r",
            coord_out,
        ]
        if step_name == "prod":
            cmd.extend(["-x", traj_name])
        if restrain:
            cmd.extend(["-ref", inpcrd_rel])

        _log(
            logger,
            "INFO",
            {"cmd": " ".join(cmd), "source": source, "step": step_name, "replicate": replicate_index},
            "md_run",
        )

        proc = subprocess.run(cmd, cwd=str(md_dir), stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        if proc.returncode != 0:
            _log(
                logger,
                "ERROR",
                {"step": step_name, "rc": proc.returncode, "log": md_dir / out_file},
                "md_failed",
            )
            return {
                "ok": False,
                "run": True,
                "skipped": False,
                "traj_path": str(traj_path),
                "traj_format": traj_format,
                "n_frames_est": n_frames_est,
                "replicate": replicate_index,
                "seed": seed,
                "md_dir": str(md_dir),
                "replicate_dir": str(rep_dir),
            }

    if not traj_path.exists() or traj_path.stat().st_size == 0:
        _log(
            logger,
            "ERROR",
            {"traj": traj_path, "replicate": replicate_index},
            "md_missing_trajectory",
        )
        return {
            "ok": False,
            "run": True,
            "skipped": False,
            "traj_path": str(traj_path),
            "traj_format": traj_format,
            "n_frames_est": n_frames_est,
            "replicate": replicate_index,
            "seed": seed,
            "md_dir": str(md_dir),
            "replicate_dir": str(rep_dir),
        }

    _log(
        logger,
        "INFO",
        {"traj": traj_path, "replicate": replicate_index, "ok": "true"},
        "md_done",
    )
    return {
        "ok": True,
        "run": True,
        "skipped": False,
        "traj_path": str(traj_path),
        "traj_format": traj_format,
        "n_frames_est": n_frames_est,
        "replicate": replicate_index,
        "seed": seed,
        "md_dir": str(md_dir),
        "replicate_dir": str(rep_dir),
    }


def make_mmgbsa_trajectory(
    complex_prmtop: str,
    complex_inpcrd: str,
    out_dir: str,
    cfg: object,
    force: bool = False,
    run_cpptraj: bool = True,
) -> dict:
    logger = _get_logger()
    prmtop_path = Path(complex_prmtop)
    inpcrd_path = Path(complex_inpcrd)
    out_path = Path(out_dir)

    if not prmtop_path.exists():
        raise FileNotFoundError(f"prmtop not found: {prmtop_path}")
    if not inpcrd_path.exists():
        raise FileNotFoundError(f"inpcrd not found: {inpcrd_path}")

    enabled = _to_bool(_cfg_get(cfg, "MMGBSA_CPPTRAJ_ENABLED", True), default=True)
    if not enabled:
        _log(logger, "INFO", {"enabled": False}, "cpptraj_disabled")
        return {
            "enabled": False,
            "cpptraj_in": "",
            "log_path": "",
            "trajout_path": "",
            "cpptraj": None,
        }

    out_path.mkdir(parents=True, exist_ok=True)

    trajout_name = str(_cfg_get(cfg, "MMGBSA_TRAJOUT_NAME", "mdcrd") or "mdcrd")
    trajout_format = str(_cfg_get(cfg, "MMGBSA_TRAJOUT_FORMAT", "mdcrd") or "mdcrd")
    trajin_source = str(_cfg_get(cfg, "MMGBSA_TRAJIN_SOURCE", "INPCRD") or "INPCRD").strip().upper()
    trajin_format = str(_cfg_get(cfg, "MMGBSA_TRAJIN_FORMAT", "inpcrd") or "inpcrd")

    startframe = _to_int(_cfg_get(cfg, "MMGBSA_TRAJ_STARTFRAME", 1), 1)
    endframe = _to_int(_cfg_get(cfg, "MMGBSA_TRAJ_ENDFRAME", 1), 1)
    interval = _to_int(_cfg_get(cfg, "MMGBSA_TRAJ_INTERVAL", 1), 1)

    if trajin_source == "INPCRD":
        trajin_path = inpcrd_path
    elif trajin_source == "EXTERNAL":
        trajin_value = str(_cfg_get(cfg, "MMGBSA_TRAJIN_PATH", "") or "")
        if not trajin_value:
            raise ValueError("MMGBSA_TRAJIN_PATH is required when MMGBSA_TRAJIN_SOURCE=EXTERNAL")
        trajin_path = Path(trajin_value)
        if not trajin_path.exists():
            raise FileNotFoundError(f"external trajin not found: {trajin_path}")
    else:
        raise ValueError(f"MMGBSA_TRAJIN_SOURCE must be INPCRD or EXTERNAL (got {trajin_source})")

    cpptraj_in_path = out_path / "cpptraj_mmgbsa.in"
    trajout_path = out_path / trajout_name

    prmtop_rel = os.path.relpath(prmtop_path, out_path)
    trajin_rel = os.path.relpath(trajin_path, out_path)
    trajout_rel = os.path.relpath(trajout_path, out_path)

    cpptraj_text = build_cpptraj_input(
        prmtop_path=prmtop_rel,
        trajin_path=trajin_rel,
        trajin_format=trajin_format,
        startframe=startframe,
        endframe=endframe,
        interval=interval,
        trajout_path=trajout_rel,
        trajout_format=trajout_format,
    )

    cpptraj_in_written = write_cpptraj_file(cpptraj_text, str(cpptraj_in_path), force=force)

    run_flag = _to_bool(_cfg_get(cfg, "MMGBSA_CPPTRAJ_RUN", True), default=True)
    if not run_cpptraj:
        run_flag = False

    trajout_exists = trajout_path.exists() and trajout_path.stat().st_size > 0
    if trajout_exists and not force:
        run_flag = False

    _log(
        logger,
        "INFO",
        {
            "prmtop": prmtop_path,
            "trajin": trajin_path,
            "trajout": trajout_path,
            "run": run_flag,
            "force": force,
        },
        "cpptraj_plan",
    )

    cpptraj_result = None
    if run_flag:
        cpptraj_result = _run_cpptraj_impl(cpptraj_in_written, str(out_path))

    _log(logger, "INFO", {"ok": "true", "trajout": trajout_path}, "done")

    return {
        "enabled": True,
        "cpptraj_in": str(cpptraj_in_path),
        "log_path": str(out_path / "cpptraj.log"),
        "trajout_path": str(trajout_path),
        "cpptraj": cpptraj_result,
    }


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate an MMGBSA trajectory via cpptraj.")
    parser.add_argument("--prmtop", required=True, help="Complex prmtop path.")
    parser.add_argument("--inpcrd", required=True, help="Complex inpcrd path.")
    parser.add_argument("--out-dir", required=True, help="Output directory for cpptraj input and trajectory.")
    parser.add_argument("--force", action="store_true", help="Overwrite outputs if present.")
    parser.add_argument("--no-run", action="store_true", help="Write cpptraj input only.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_arg_parser()
    args = parser.parse_args(argv)
    logger = _get_logger()

    try:
        cfg = load_config()
        make_mmgbsa_trajectory(
            complex_prmtop=args.prmtop,
            complex_inpcrd=args.inpcrd,
            out_dir=args.out_dir,
            cfg=cfg,
            force=args.force,
            run_cpptraj=not args.no_run,
        )
    except Exception as exc:
        _log(logger, "ERROR", {"reason": type(exc).__name__, "detail": str(exc)}, "failed")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
