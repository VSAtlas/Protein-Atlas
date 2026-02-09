from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import uuid
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple


_COMPONENT = "mmgbsa.prep"


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


def _ensure_file_handler(logger: logging.Logger, log_file: Path) -> None:
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            if (
                Path(getattr(handler, "baseFilename", "")).resolve()
                == log_file.resolve()
            ):
                return
    log_file.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_file, mode="a", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(fh)


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


def _log_custom(
    logger: logging.Logger,
    component: str,
    level: str,
    kvs: Dict[str, object],
    msg: str | None = None,
) -> None:
    level_map = {
        "DEBUG": logging.DEBUG,
        "INFO": logging.INFO,
        "WARNING": logging.WARNING,
        "ERROR": logging.ERROR,
    }
    parts = [f"[{level}]", f"[{component}]"]
    if kvs:
        parts.append(" ".join(f"{k}={v}" for k, v in kvs.items()))
    if msg:
        parts.append(f"| {msg}")
    line = " ".join(parts)
    logger.log(level_map[level], line)


def _cfg_get(cfg: object | None, key: str, default: object) -> object:
    if cfg is None:
        return default
    if isinstance(cfg, dict):
        return cfg.get(key, default)
    try:
        return cfg.get(key, default)
    except Exception:
        return getattr(cfg, key, default)


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


def _parse_int_list(value: object) -> List[int]:
    if value is None:
        return []
    text = str(value).strip()
    if not text:
        return []
    out: List[int] = []
    for chunk in text.replace(",", " ").split():
        token = chunk.strip()
        if not token:
            continue
        try:
            out.append(int(token))
        except Exception:
            continue
    return out


def _resolve_post_docked_paths(sdf_path: Path) -> Tuple[str, str, str, str, str, Path]:
    try:
        resolved = sdf_path.resolve(strict=False)
    except Exception:
        resolved = sdf_path

    parts = resolved.parts
    try:
        idx = parts.index("post_docked")
    except ValueError:
        raise ValueError(
            "sdf_path must live under post_docked/<runid>/<pdb>/<variant>/<pH>/<stage_dir>/file.sdf"
        )

    if len(parts) <= idx + 5:
        raise ValueError(
            "sdf_path does not include runid/pdb/variant/pH/stage_dir components"
        )

    runid = parts[idx + 1]
    pdb = parts[idx + 2]
    variant = parts[idx + 3]
    ph_label = parts[idx + 4]
    expected_ph_dir = Path(*parts[: idx + 5])

    stage_dir = resolved.parent.name
    if resolved.parent.parent != expected_ph_dir:
        raise ValueError(
            "sdf_path must be under post_docked/<runid>/<pdb>/<variant>/<pH>/<stage_dir>/file.sdf"
        )

    mmgbsa_dir = expected_ph_dir / "mmgbsa"
    return runid, pdb, variant, ph_label, stage_dir, mmgbsa_dir


def _read_rdkit_info(
    sdf_path: Path,
    logger: logging.Logger,
) -> Tuple[Optional[int], Optional[int]]:
    try:
        from rdkit import Chem
    except Exception as exc:
        _log(
            logger,
            "WARNING",
            {
                "reason": "rdkit_unavailable",
                "detail": str(exc),
                "fallback": "skip_rdkit",
            },
            "rdkit_unavailable",
        )
        return None, None

    try:
        supplier = Chem.SDMolSupplier(str(sdf_path), sanitize=True, removeHs=False)
        mol = next((m for m in supplier if m is not None), None)
    except Exception:
        mol = None

    if mol is None:
        try:
            supplier = Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)
            mol = next((m for m in supplier if m is not None), None)
        except Exception:
            mol = None

    if mol is None:
        _log(
            logger,
            "WARNING",
            {"reason": "rdkit_parse_failed", "fallback": "skip_rdkit"},
            "rdkit_parse_failed",
        )
        return None, None

    try:
        formal_charge = int(Chem.GetFormalCharge(mol))
    except Exception:
        formal_charge = None
    try:
        radical_electrons = int(
            sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
        )
    except Exception:
        radical_electrons = None

    return formal_charge, radical_electrons


def _compute_net_charge_rdkit(sdf_path: Path, logger: logging.Logger) -> Optional[int]:
    formal_charge, _ = _read_rdkit_info(sdf_path, logger)
    return formal_charge


def _looks_like_pkgs_cache(prefix: Path) -> bool:
    prefix_str = prefix.as_posix()
    return "/micromamba/pkgs" in prefix_str


def _prefix_has_tools(prefix: Path) -> bool:
    return (prefix / "bin" / "antechamber").is_file() and (
        prefix / "bin" / "parmchk2"
    ).is_file()


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


def _common_ambertools_prefixes() -> list[Path]:
    candidates: list[Path] = []
    env_candidates = [
        os.environ.get("MMGBSA_AMBERTOOLS_PREFIX"),
        os.environ.get("AMBERTOOLS_PREFIX"),
        os.environ.get("CONDA_PREFIX"),
    ]
    for raw in env_candidates:
        if raw:
            candidates.append(Path(raw).expanduser())

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


def _select_ambertools(
    amber_prefix: Optional[str],
    logger: logging.Logger,
) -> Tuple[List[str], str, str, str]:
    if amber_prefix:
        prefix = Path(amber_prefix).expanduser()
        if _looks_like_pkgs_cache(prefix):
            raise ValueError(
                f"refusing amber_prefix under pkgs cache: {prefix}; use an environment prefix instead"
            )
        if not _prefix_has_tools(prefix):
            raise FileNotFoundError(
                f"amber_prefix missing AmberTools binaries: {prefix}"
            )
        micromamba = _resolve_micromamba()
        if not micromamba:
            raise FileNotFoundError(
                "micromamba not found; cannot run AmberTools from prefix"
            )
        runner = [str(micromamba), "run", "-p", str(prefix)]
        return runner, "antechamber", "parmchk2", f"prefix:{prefix}"

    env_prefix = os.environ.get("AMBERTOOLS_PREFIX")
    if env_prefix:
        prefix = Path(env_prefix).expanduser()
        if _looks_like_pkgs_cache(prefix):
            _log(
                logger,
                "WARNING",
                {"reason": "pkgs_cache_prefix", "fallback": "PATH_or_common_prefix"},
                "ignoring_AMBERTOOLS_PREFIX",
            )
        elif _prefix_has_tools(prefix):
            micromamba = _resolve_micromamba()
            if not micromamba:
                raise FileNotFoundError(
                    "micromamba not found; cannot run AmberTools from prefix"
                )
            runner = [str(micromamba), "run", "-p", str(prefix)]
            return runner, "antechamber", "parmchk2", f"prefix:{prefix}"
        else:
            _log(
                logger,
                "WARNING",
                {"reason": "missing_tools", "fallback": "PATH_or_common_prefix"},
                "ignoring_AMBERTOOLS_PREFIX",
            )

    antechamber_path = shutil.which("antechamber")
    parmchk2_path = shutil.which("parmchk2")
    if antechamber_path and parmchk2_path:
        return [], antechamber_path, parmchk2_path, "PATH"

    common_prefixes = _common_ambertools_prefixes()
    for prefix in common_prefixes:
        if _prefix_has_tools(prefix) and not _looks_like_pkgs_cache(prefix):
            micromamba = _resolve_micromamba()
            if not micromamba:
                raise FileNotFoundError(
                    "micromamba not found; cannot run AmberTools from prefix"
                )
            runner = [str(micromamba), "run", "-p", str(prefix)]
            return runner, "antechamber", "parmchk2", f"prefix:{prefix}"

    raise FileNotFoundError(
        "could not locate AmberTools antechamber/parmchk2; set AMBERTOOLS_PREFIX or PATH"
    )


def _tmp_path(final_path: Path) -> Path:
    return final_path.with_name(f"{final_path.name}.tmp.{uuid.uuid4().hex}")


def _write_json_atomic(path: Path, payload: Dict[str, object]) -> None:
    tmp_path = _tmp_path(path)
    with tmp_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    if not tmp_path.exists() or tmp_path.stat().st_size == 0:
        if tmp_path.exists():
            tmp_path.unlink()
        raise RuntimeError(f"output not created: {tmp_path}")
    os.replace(tmp_path, path)


def _tail_log(path: Path, max_lines: int = 50) -> str:
    try:
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return ""
    tail = lines[-max_lines:]
    return "\n".join(tail)


def _resolve_nominal_charge(
    sdf_path: Path,
    logger: logging.Logger,
    cfg: object | None,
    override: Optional[int],
) -> int:
    if override is not None:
        return int(override)

    cfg_value = _cfg_get(cfg, "MMGBSA_LIGAND_NOMINAL_NET_CHARGE", None)
    if cfg_value is None or str(cfg_value).strip() == "":
        cfg_value = _cfg_get(cfg, "MMGBSA_LIGAND_NET_CHARGE", None)

    if cfg_value is not None and str(cfg_value).strip() != "":
        text = str(cfg_value).strip().lower()
        if text in {"auto", "rdkit"}:
            computed = _compute_net_charge_rdkit(sdf_path, logger)
            return computed if computed is not None else 0
        return _to_int(cfg_value, 0)

    computed = _compute_net_charge_rdkit(sdf_path, logger)
    return computed if computed is not None else 0


def _resolve_net_charge(
    sdf_path: Path,
    logger: logging.Logger,
    cfg: object | None,
    net_charge: Optional[int],
) -> int:
    return _resolve_nominal_charge(sdf_path, logger, cfg, net_charge)


def _run_command_to_log(
    cmd: List[str],
    log_path: Path,
    logger: logging.Logger,
    cwd: Optional[Path] = None,
) -> int:
    _log(logger, "DEBUG", {"cmd": " ".join(cmd), "log": log_path})
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as handle:
        proc = subprocess.run(
            cmd, cwd=str(cwd) if cwd else None, stdout=handle, stderr=subprocess.STDOUT
        )
    return int(proc.returncode)


def _log_attempt_failure(
    logger: logging.Logger,
    cmd: List[str],
    rc: int,
    log_path: Path,
    context: str,
) -> None:
    _log(logger, "ERROR", {"cmd": " ".join(cmd), "rc": rc, "log": log_path}, context)
    tail = _tail_log(log_path)
    if tail:
        tail_text = tail.replace("\n", "\\n")
        _log(logger, "ERROR", {"tail": tail_text}, f"{context}_tail")


def _build_bcc_sweep(cfg: object | None) -> List[int]:
    sweep = _parse_int_list(_cfg_get(cfg, "MMGBSA_LIGAND_BCC_CHARGE_SWEEP", "-1,1"))
    include_pm2 = _to_bool(
        _cfg_get(cfg, "MMGBSA_LIGAND_BCC_SWEEP_INCLUDE_PLUSMINUS2", False),
        default=False,
    )
    if include_pm2:
        for val in (-2, 2):
            if val not in sweep:
                sweep.append(val)
    return sweep


def parameterize_ligand_with_fallback(
    sdf_path: str,
    out_mol2: str,
    out_frcmod: str,
    cfg: object,
    stage_dir: str,
    ligand_stem: str,
    force: bool = False,
) -> dict:
    logger = _get_logger()
    sdf = Path(sdf_path)
    mol2_path = Path(out_mol2)
    frcmod_path = Path(out_frcmod)

    if not sdf.exists():
        raise FileNotFoundError(f"sdf_path not found: {sdf}")

    mmgbsa_dir = mol2_path.parent.parent.parent
    log_dir = mmgbsa_dir / "logs" / stage_dir / ligand_stem
    log_dir.mkdir(parents=True, exist_ok=True)

    metadata_path = mol2_path.with_suffix(".mmgbsa_prep.json")
    work_dir = mmgbsa_dir / "ambertools_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    rdkit_validate = _to_bool(
        _cfg_get(cfg, "MMGBSA_RDKit_VALIDATE", True), default=True
    )
    rdkit_formal_charge: Optional[int] = None
    rdkit_radical_electrons: Optional[int] = None
    low_confidence = False
    if rdkit_validate:
        rdkit_formal_charge, rdkit_radical_electrons = _read_rdkit_info(sdf, logger)
        if rdkit_radical_electrons is not None:
            threshold = _to_int(
                _cfg_get(cfg, "MMGBSA_RDKit_RADICAL_LOWCONF_THRESHOLD", 1), 1
            )
            if threshold <= 0:
                threshold = 1
            if rdkit_radical_electrons >= threshold:
                low_confidence = True
                _log_custom(
                    logger,
                    "mmgbsa.ligand.radical",
                    "INFO",
                    {
                        "ligand": ligand_stem,
                        "radicals": rdkit_radical_electrons,
                        "formal_charge": rdkit_formal_charge,
                    },
                )

    nominal_charge = _resolve_nominal_charge(sdf, logger, cfg, None)

    ligand_at = (
        str(_cfg_get(cfg, "MMGBSA_LIGAND_AT", "gaff2") or "gaff2").strip() or "gaff2"
    )
    primary_method = (
        str(
            _cfg_get(
                cfg,
                "MMGBSA_LIGAND_PRIMARY_CHARGE_METHOD",
                _cfg_get(cfg, "MMGBSA_LIGAND_CHARGE_METHOD", "bcc"),
            )
            or "bcc"
        )
        .strip()
        .lower()
    )
    if primary_method not in {"bcc", "gas"}:
        primary_method = "bcc"
    fallback_method = (
        str(_cfg_get(cfg, "MMGBSA_LIGAND_FALLBACK_CHARGE_METHOD", "gas") or "gas")
        .strip()
        .lower()
    )
    if fallback_method not in {"bcc", "gas"}:
        fallback_method = "gas"

    sqm_level_val = _to_int(_cfg_get(cfg, "MMGBSA_LIGAND_SQM_LEVEL", 2), 2)

    outputs_ok = (
        mol2_path.exists()
        and mol2_path.stat().st_size > 0
        and frcmod_path.exists()
        and frcmod_path.stat().st_size > 0
    )
    if outputs_ok and not force:
        metadata = {
            "ok": True,
            "input_sdf": str(sdf),
            "output_mol2": str(mol2_path),
            "output_frcmod": str(frcmod_path),
            "charge_method": primary_method,
            "net_charge_used": nominal_charge,
            "bcc_attempts": [],
            "gas_fallback_used": False,
            "rdkit_formal_charge": rdkit_formal_charge,
            "rdkit_radical_electrons": rdkit_radical_electrons,
            "low_confidence": low_confidence,
            "notes": "skipped_existing",
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "skipped_existing": True,
        }
        _write_json_atomic(metadata_path, metadata)
        return {
            "ok": True,
            "skipped_existing": True,
            "output_mol2": str(mol2_path),
            "output_frcmod": str(frcmod_path),
            "metadata_path": str(metadata_path),
            "charge_method": primary_method,
            "net_charge_used": nominal_charge,
            "low_confidence": low_confidence,
            "commands_run": [],
        }

    amber_prefix = _cfg_get(cfg, "MMGBSA_AMBERTOOLS_PREFIX", None)
    if amber_prefix is None or str(amber_prefix).strip() == "":
        amber_prefix = _cfg_get(cfg, "AMBERTOOLS_PREFIX", None)
    amber_prefix = str(amber_prefix).strip() if amber_prefix else None

    runner, antechamber_bin, parmchk2_bin, tool_source = _select_ambertools(
        amber_prefix, logger
    )
    antechamber_label = (
        " ".join(runner + [antechamber_bin]) if runner else antechamber_bin
    )
    parmchk2_label = " ".join(runner + [parmchk2_bin]) if runner else parmchk2_bin
    _log(
        logger,
        "INFO",
        {
            "antechamber": antechamber_label,
            "parmchk2": parmchk2_label,
            "source": tool_source,
        },
    )

    commands_run: List[str] = []
    bcc_attempts: List[dict] = []
    gas_fallback_used = False
    charge_method_used: Optional[str] = None
    net_charge_used: Optional[int] = None

    def _attempt_antechamber(method: str, charge: int) -> bool:
        attempt_log = log_dir / f"antechamber_attempt_{method}_nc{charge}.log"
        mol2_tmp = _tmp_path(mol2_path)
        cmd = runner + [
            antechamber_bin,
            "-i",
            str(sdf),
            "-fi",
            "sdf",
            "-o",
            str(mol2_tmp),
            "-fo",
            "mol2",
            "-at",
            ligand_at,
            "-c",
            method,
            "-nc",
            str(charge),
            "-s",
            str(sqm_level_val),
        ]
        rc = _run_command_to_log(cmd, attempt_log, logger, cwd=work_dir)
        ok = rc == 0 and mol2_tmp.exists() and mol2_tmp.stat().st_size > 0
        if ok:
            os.replace(mol2_tmp, mol2_path)
        else:
            if mol2_tmp.exists():
                mol2_tmp.unlink()
            _log_attempt_failure(logger, cmd, rc, attempt_log, "antechamber_failed")
        commands_run.append(" ".join(cmd))
        if method == "bcc":
            bcc_attempts.append(
                {"nc": charge, "ok": ok, "rc": rc, "log": str(attempt_log)}
            )
        return ok

    succeeded = False

    if primary_method == "bcc":
        if _attempt_antechamber("bcc", nominal_charge):
            succeeded = True
            charge_method_used = "bcc"
            net_charge_used = nominal_charge
        else:
            sweep = _build_bcc_sweep(cfg)
            for charge in sweep:
                if charge == nominal_charge:
                    continue
                if _attempt_antechamber("bcc", charge):
                    succeeded = True
                    charge_method_used = "bcc"
                    net_charge_used = charge
                    break

    if not succeeded and fallback_method == "gas":
        gas_fallback_used = True
        if _attempt_antechamber("gas", nominal_charge):
            succeeded = True
            charge_method_used = "gas"
            net_charge_used = nominal_charge
        else:
            charge_method_used = "gas"
            net_charge_used = nominal_charge

    parmchk2_ok = False
    parmchk2_log = log_dir / "parmchk2.log"
    if succeeded:
        frcmod_tmp = _tmp_path(frcmod_path)
        parmchk2_cmd = runner + [
            parmchk2_bin,
            "-i",
            str(mol2_path),
            "-f",
            "mol2",
            "-o",
            str(frcmod_tmp),
        ]
        rc = _run_command_to_log(parmchk2_cmd, parmchk2_log, logger, cwd=work_dir)
        parmchk2_ok = rc == 0 and frcmod_tmp.exists() and frcmod_tmp.stat().st_size > 0
        if parmchk2_ok:
            os.replace(frcmod_tmp, frcmod_path)
        else:
            if frcmod_tmp.exists():
                frcmod_tmp.unlink()
            _log_attempt_failure(
                logger, parmchk2_cmd, rc, parmchk2_log, "parmchk2_failed"
            )
        commands_run.append(" ".join(parmchk2_cmd))

    ok = succeeded and parmchk2_ok
    notes = ""
    if not ok:
        notes = "param_failed" if succeeded else "antechamber_failed"

    metadata = {
        "ok": ok,
        "input_sdf": str(sdf),
        "output_mol2": str(mol2_path),
        "output_frcmod": str(frcmod_path),
        "charge_method": charge_method_used or primary_method,
        "net_charge_used": int(net_charge_used)
        if net_charge_used is not None
        else nominal_charge,
        "bcc_attempts": bcc_attempts,
        "gas_fallback_used": gas_fallback_used,
        "rdkit_formal_charge": rdkit_formal_charge,
        "rdkit_radical_electrons": rdkit_radical_electrons,
        "low_confidence": low_confidence,
        "notes": notes,
        "timestamp": datetime.utcnow().isoformat() + "Z",
    }
    _write_json_atomic(metadata_path, metadata)

    return {
        "ok": ok,
        "output_mol2": str(mol2_path),
        "output_frcmod": str(frcmod_path),
        "metadata_path": str(metadata_path),
        "charge_method": charge_method_used or primary_method,
        "net_charge_used": int(net_charge_used)
        if net_charge_used is not None
        else nominal_charge,
        "low_confidence": low_confidence,
        "commands_run": commands_run,
        "bcc_attempts": bcc_attempts,
        "gas_fallback_used": gas_fallback_used,
    }


def prep_mmgbsa_from_sdf(
    sdf_path: str,
    net_charge: Optional[int] = None,
    amber_prefix: Optional[str] = None,
    dry_run: bool = False,
    force: bool = False,
    ligand_at: Optional[str] = None,
    charge_method: Optional[str] = None,
    sqm_level: Optional[int] = None,
    cfg: object | None = None,
) -> dict:
    logger = _get_logger()
    path = Path(sdf_path)

    if path.suffix.lower() != ".sdf":
        raise ValueError(f"sdf_path must end with .sdf, got: {path}")
    if not path.exists():
        raise FileNotFoundError(f"sdf_path not found: {path}")

    _log(logger, "INFO", {"input": path})

    runid, pdb, variant, ph_label, stage_dir, mmgbsa_dir = _resolve_post_docked_paths(
        path
    )
    if not dry_run:
        _ensure_file_handler(logger, mmgbsa_dir / "mmgbsa_prep.log")

    _log(logger, "INFO", {"stage_dir": stage_dir, "mmgbsa_dir": mmgbsa_dir})
    work_dir = mmgbsa_dir / "ambertools_work"
    _log(logger, "INFO", {"work_dir": work_dir})
    _log(
        logger,
        "DEBUG",
        {"runid": runid, "pdb": pdb, "variant": variant, "ph": ph_label},
    )

    ligand_base = path.stem
    mol2_path = mmgbsa_dir / "mol2" / stage_dir / f"{ligand_base}.mol2"
    frcmod_path = mmgbsa_dir / "frcmod" / stage_dir / f"{ligand_base}.frcmod"

    _log(logger, "INFO", {"mol2": mol2_path, "frcmod": frcmod_path})

    cfg_local: dict = dict(cfg) if isinstance(cfg, dict) else {}
    if amber_prefix is not None and str(amber_prefix).strip():
        cfg_local["MMGBSA_AMBERTOOLS_PREFIX"] = str(amber_prefix).strip()
    if net_charge is not None:
        cfg_local["MMGBSA_LIGAND_NOMINAL_NET_CHARGE"] = int(net_charge)
        cfg_local["MMGBSA_LIGAND_NET_CHARGE"] = int(net_charge)
    if ligand_at:
        cfg_local["MMGBSA_LIGAND_AT"] = str(ligand_at)
    if charge_method:
        cfg_local["MMGBSA_LIGAND_PRIMARY_CHARGE_METHOD"] = str(charge_method)
    if sqm_level is not None:
        cfg_local["MMGBSA_LIGAND_SQM_LEVEL"] = int(sqm_level)

    force_val = force or _to_bool(
        _cfg_get(cfg_local, "MMGBSA_LIGAND_FORCE", False), default=False
    )
    if _to_bool(_cfg_get(cfg_local, "MMGBSA_FORCE", False), default=False):
        force_val = True

    if dry_run:
        _log(logger, "INFO", {"dry_run": True}, "skip_parameterization")
        return {
            "sdf_path": str(path),
            "stage_dir": stage_dir,
            "mmgbsa_dir": str(mmgbsa_dir),
            "mol2_path": str(mol2_path),
            "frcmod_path": str(frcmod_path),
            "net_charge_used": None,
            "commands_run": [],
            "dry_run": True,
        }

    mol2_path.parent.mkdir(parents=True, exist_ok=True)
    frcmod_path.parent.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    result = parameterize_ligand_with_fallback(
        sdf_path=str(path),
        out_mol2=str(mol2_path),
        out_frcmod=str(frcmod_path),
        cfg=cfg_local,
        stage_dir=stage_dir,
        ligand_stem=ligand_base,
        force=force_val,
    )

    if not result.get("ok"):
        raise RuntimeError(result.get("notes") or "mmgbsa_prep_failed")

    return {
        "sdf_path": str(path),
        "stage_dir": stage_dir,
        "mmgbsa_dir": str(mmgbsa_dir),
        "mol2_path": str(mol2_path),
        "frcmod_path": str(frcmod_path),
        "net_charge_used": result.get("net_charge_used"),
        "charge_method": result.get("charge_method"),
        "low_confidence": result.get("low_confidence"),
        "gas_fallback_used": result.get("gas_fallback_used"),
        "metadata_path": result.get("metadata_path"),
        "commands_run": result.get("commands_run", []),
    }


def prep_mmgbsa_from_sdfs(
    sdf_paths: List[str],
    cfg: object | None = None,
    max_ligands: Optional[int] = None,
    net_charge: Optional[int] = None,
    amber_prefix: Optional[str] = None,
    dry_run: bool = False,
    force: Optional[bool] = None,
) -> List[dict]:
    logger = _get_logger()
    paths = [str(p) for p in (sdf_paths or []) if str(p).strip()]

    cfg_max = _to_int(_cfg_get(cfg, "MMGBSA_MAX_LIGANDS", 0), 0)
    if max_ligands is None and cfg_max > 0:
        max_ligands = cfg_max

    if max_ligands is not None and max_ligands > 0:
        paths = paths[:max_ligands]

    cfg_force = _to_bool(_cfg_get(cfg, "MMGBSA_FORCE", False), default=False)
    cfg_force = cfg_force or _to_bool(
        _cfg_get(cfg, "MMGBSA_LIGAND_FORCE", False), default=False
    )
    force_val = cfg_force if force is None else force

    _log(
        logger,
        "INFO",
        {
            "total": len(sdf_paths or []),
            "selected": len(paths),
            "max_ligands": max_ligands if max_ligands is not None else "all",
            "force": force_val,
        },
        "batch_start",
    )

    results: List[dict] = []
    for sdf_path in paths:
        try:
            result = prep_mmgbsa_from_sdf(
                sdf_path=sdf_path,
                net_charge=net_charge,
                amber_prefix=amber_prefix,
                dry_run=dry_run,
                force=force_val,
                cfg=cfg,
            )
            results.append(result)
        except Exception as exc:
            _log(
                logger,
                "ERROR",
                {"sdf": sdf_path, "reason": type(exc).__name__, "detail": str(exc)},
                "prep_failed",
            )
            results.append(
                {
                    "sdf_path": sdf_path,
                    "error": str(exc),
                    "ok": False,
                }
            )

    _log(logger, "INFO", {"ok": "true", "processed": len(results)}, "batch_done")
    return results


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare MMGBSA mol2/frcmod files from a docked SDF."
    )
    parser.add_argument(
        "--sdf",
        required=True,
        help="Input SDF under post_docked/<runid>/<pdb>/<variant>/<pH>/<stage_dir>/.",
    )
    parser.add_argument(
        "--net-charge",
        type=int,
        default=None,
        help="Net charge override for antechamber.",
    )
    parser.add_argument(
        "--amber-prefix",
        default=None,
        help="AmberTools environment prefix (contains bin/antechamber).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running AmberTools.",
    )
    parser.add_argument(
        "--force", action="store_true", help="Regenerate mol2/frcmod even if present."
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_arg_parser()
    args = parser.parse_args(argv)

    logger = _get_logger()
    try:
        prep_mmgbsa_from_sdf(
            sdf_path=args.sdf,
            net_charge=args.net_charge,
            amber_prefix=args.amber_prefix,
            dry_run=args.dry_run,
            force=args.force,
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
