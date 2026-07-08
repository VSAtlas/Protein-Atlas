from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from config.value_access import cfg_get, to_bool, to_int
from config.tool_resolver import (
    common_ambertools_prefixes as _common_ambertools_prefixes_shared,
    looks_like_pkgs_cache as _looks_like_pkgs_cache_shared,
    prefix_has_tools as _prefix_has_tools_shared,
    resolve_micromamba as _resolve_micromamba_shared,
    select_ambertools_tools,
)
from post_docking.mmgbsa._atomic_io import tmp_path, write_json_atomic
from post_docking.mmgbsa.mmgbsa_ligand_state import (
    analyze_ligand_state,
    suspicious_blocks_publication,
)


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
    return cfg_get(cfg, key, default)


def _to_bool(val: object, default: bool = False) -> bool:
    return bool(to_bool(val, default=default))


def _to_int(val: object, default: int) -> int:
    return int(to_int(val, default=default))


def _dict_int(payload: Dict[str, object], key: str) -> int:
    return _to_int(payload.get(key, 0), 0)


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

    tail_len = len(parts) - idx - 1
    if tail_len < 5:
        raise ValueError(
            "sdf_path does not include runid/pdb/pH/stage_dir components"
        )

    runid = parts[idx + 1]
    pdb = parts[idx + 2]
    if tail_len == 5:
        variant = "LEGACY"
        ph_label = parts[idx + 3]
        expected_ph_dir = Path(*parts[: idx + 4])
    else:
        variant = parts[idx + 3]
        ph_label = parts[idx + 4]
        expected_ph_dir = Path(*parts[: idx + 5])

    stage_dir = resolved.parent.name
    if resolved.parent.parent != expected_ph_dir:
        raise ValueError(
            "sdf_path must be under post_docked/<runid>/<pdb>/<pH>/<stage_dir>/file.sdf "
            "or post_docked/<runid>/<pdb>/<variant>/<pH>/<stage_dir>/file.sdf"
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


def _rdkit_sanitize_error(sdf_path: Path, logger: logging.Logger) -> str:
    try:
        from rdkit import Chem
    except Exception:
        return ""
    try:
        supplier = Chem.SDMolSupplier(str(sdf_path), sanitize=True, removeHs=False)
        mol = next((m for m in supplier if m is not None), None)
        if mol is not None:
            return ""
    except Exception as exc:
        return str(exc)
    try:
        supplier = Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)
        mol = next((m for m in supplier if m is not None), None)
    except Exception as exc:
        return str(exc)
    if mol is None:
        return ""
    try:
        Chem.SanitizeMol(mol)
        return ""
    except Exception as exc:
        _log(
            logger,
            "ERROR",
            {"sdf": sdf_path, "detail": str(exc)},
            "rdkit_sanitize_failed",
        )
        return str(exc)


def _compute_net_charge_rdkit(sdf_path: Path, logger: logging.Logger) -> Optional[int]:
    formal_charge, _ = _read_rdkit_info(sdf_path, logger)
    return formal_charge


def _looks_like_pkgs_cache(prefix: Path) -> bool:
    return _looks_like_pkgs_cache_shared(prefix)


def _prefix_has_tools(prefix: Path) -> bool:
    return _prefix_has_tools_shared(prefix, ("antechamber", "parmchk2"))


def _resolve_micromamba() -> Optional[Path]:
    return _resolve_micromamba_shared()


def _common_ambertools_prefixes() -> list[Path]:
    return _common_ambertools_prefixes_shared(None)


def _select_ambertools(
    amber_prefix: Optional[str],
    logger: logging.Logger,
) -> Tuple[List[str], str, str, str]:
    del logger
    runner, executables, source = select_ambertools_tools(
        None,
        ("antechamber", "parmchk2"),
        explicit_prefix=amber_prefix,
        strict_explicit_prefix=bool(amber_prefix),
    )
    return runner, executables["antechamber"], executables["parmchk2"], source


def _tmp_path(final_path: Path) -> Path:
    return tmp_path(final_path)


def _write_json_atomic(path: Path, payload: Dict[str, object]) -> None:
    write_json_atomic(path, payload, require_nonempty=True)


def _tail_log(path: Path, max_lines: int = 50) -> str:
    try:
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return ""
    tail = lines[-max_lines:]
    return "\n".join(tail)


def _review_status(cfg: object | None, key: str) -> str:
    return str(_cfg_get(cfg, key, "unreviewed") or "unreviewed").strip().lower()


def _review_ok(status: str) -> bool:
    return status in {"reviewed", "approved", "ok"}


def _chemistry_review_payload(cfg: object | None) -> Dict[str, object]:
    keys = (
        "MMGBSA_CHEMISTRY_REVIEW_STATUS",
        "MMGBSA_PROTONATION_REVIEW_STATUS",
        "MMGBSA_TAUTOMER_REVIEW_STATUS",
        "MMGBSA_STEREOCHEMISTRY_REVIEW_STATUS",
        "MMGBSA_NET_CHARGE_REVIEW_STATUS",
        "MMGBSA_PARAMETER_REVIEW_STATUS",
    )
    statuses = {key: _review_status(cfg, key) for key in keys}
    required = _to_bool(
        _cfg_get(cfg, "MMGBSA_CHEMISTRY_SUBGATES_REQUIRED", False), default=False
    )
    missing = [key for key, status in statuses.items() if not _review_ok(status)]
    return {
        "statuses": statuses,
        "subgates_required": required,
        "subgates_ok": not missing,
        "missing_review_keys": missing,
        "notes": str(_cfg_get(cfg, "MMGBSA_CHEMISTRY_REVIEW_NOTES", "") or ""),
    }


def _input_chemistry_payload(cfg: object | None) -> Dict[str, object]:
    authoritative = _to_bool(
        _cfg_get(cfg, "MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE", False), default=False
    )
    return {
        "authoritative": authoritative,
        "source": str(
            _cfg_get(cfg, "MMGBSA_INPUT_CHEMISTRY_SOURCE", "unreviewed_sdf")
            or "unreviewed_sdf"
        ),
        "notes": str(_cfg_get(cfg, "MMGBSA_INPUT_CHEMISTRY_NOTES", "") or ""),
    }


def _read_existing_metadata(path: Path) -> Dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _mol2_atom_type_warnings(mol2_path: Path) -> Dict[str, object]:
    try:
        lines = mol2_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return {"unknown_atom_type_count": 0, "unknown_atom_type_examples": []}

    in_atoms = False
    unknown: List[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if stripped.startswith("@<TRIPOS>") and in_atoms:
            break
        if not in_atoms or not stripped:
            continue
        parts = stripped.split()
        if len(parts) >= 6 and parts[5].upper() in {"DU", "DUM", "XX"}:
            unknown.append(parts[5])
    return {
        "unknown_atom_type_count": len(unknown),
        "unknown_atom_type_examples": unknown[:10],
    }


def _frcmod_warning_summary(frcmod_path: Path) -> Dict[str, object]:
    try:
        lines = frcmod_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return {
            "frcmod_warning_count": 0,
            "frcmod_du_count": 0,
            "frcmod_warning_examples": [],
        }

    warnings: List[str] = []
    du_count = 0
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        upper = line.upper()
        if upper.strip() == "REMARK LINE GOES HERE":
            continue
        if "ATTN" in upper or "NEED" in upper or "REMARK" in upper:
            warnings.append(stripped)
        if "DU" in upper:
            du_count += 1
    return {
        "frcmod_warning_count": len(warnings),
        "frcmod_du_count": du_count,
        "frcmod_warning_examples": warnings[:20],
    }


def _parameter_warning_summary(
    mol2_path: Path, frcmod_path: Path
) -> Tuple[Dict[str, object], Dict[str, object]]:
    mol2_warnings = (
        _mol2_atom_type_warnings(mol2_path)
        if mol2_path.exists()
        else {"unknown_atom_type_count": 0, "unknown_atom_type_examples": []}
    )
    frcmod_warnings = (
        _frcmod_warning_summary(frcmod_path)
        if frcmod_path.exists()
        else {
            "frcmod_warning_count": 0,
            "frcmod_du_count": 0,
            "frcmod_warning_examples": [],
        }
    )
    return mol2_warnings, frcmod_warnings


def _publication_ready_gate(
    *,
    parameterized: bool,
    low_confidence: bool,
    input_chemistry: Dict[str, object],
    charge_method: str,
    gas_fallback_used: bool,
    chemistry_review: Dict[str, object],
    suspicious_chemistry: Dict[str, object],
    mol2_warnings: Dict[str, object],
    frcmod_warnings: Dict[str, object],
) -> bool:
    charge_method_normalized = charge_method.lower()
    return (
        parameterized
        and not low_confidence
        and bool(input_chemistry["authoritative"])
        and charge_method_normalized not in {"gas", "mul"}
        and not gas_fallback_used
        and not suspicious_blocks_publication(suspicious_chemistry)
        and _dict_int(mol2_warnings, "unknown_atom_type_count") == 0
        and _dict_int(frcmod_warnings, "frcmod_warning_count") == 0
        and _dict_int(frcmod_warnings, "frcmod_du_count") == 0
        and (
            not chemistry_review["subgates_required"]
            or bool(chemistry_review["subgates_ok"])
        )
    )


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
    timeout: Optional[int] = None,
) -> int:
    _log(logger, "DEBUG", {"cmd": " ".join(cmd), "log": log_path})
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as handle:
        try:
            proc = subprocess.run(
                cmd,
                cwd=str(cwd) if cwd else None,
                stdout=handle,
                stderr=subprocess.STDOUT,
                timeout=timeout,
            )
            return int(proc.returncode)
        except TypeError as exc:
            if "timeout" not in str(exc):
                raise
            proc = subprocess.run(
                cmd,
                cwd=str(cwd) if cwd else None,
                stdout=handle,
                stderr=subprocess.STDOUT,
            )
            return int(proc.returncode)
        except subprocess.TimeoutExpired as exc:
            handle.write(
                f"\nTIMEOUT: command exceeded {timeout} seconds and was stopped.\n"
            )
            if exc.stdout:
                handle.write(str(exc.stdout))
            if exc.stderr:
                handle.write(str(exc.stderr))
            _log(logger, "ERROR", {"cmd": " ".join(cmd), "timeout": timeout}, "command_timeout")
            return 124


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
        _cfg_get(cfg, "MMGBSA_LIGAND_BCC_SWEEP_INCLUDE_PLUSMINUS2", True),
        default=True,
    )
    if include_pm2:
        for val in (-2, 2):
            if val not in sweep:
                sweep.append(val)
    return sweep


def _resp_mol2_candidates(sdf: Path) -> List[Path]:
    stem = sdf.stem
    parent = sdf.parent
    names = [
        f"{stem}.resp.mol2",
        f"{stem}.qm_resp.mol2",
        f"{stem}.qm.mol2",
        f"{stem}.precharged.mol2",
        f"{stem}.mol2",
    ]
    candidates = [parent / name for name in names]
    candidates.extend(
        parent / subdir / name
        for subdir in ("resp", "qm", "mol2")
        for name in names
    )
    return candidates


def _find_resp_mol2(sdf: Path) -> Optional[Path]:
    for candidate in _resp_mol2_candidates(sdf):
        try:
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
        except Exception:
            continue
    return None


def _resp_like_charge_methods() -> List[str]:
    return ["resp", "mul"]


def _mol2_total_charge(mol2_path: Path) -> Optional[int]:
    try:
        lines = mol2_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return None
    in_atoms = False
    total = 0.0
    seen_charge = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if stripped.startswith("@<TRIPOS>") and in_atoms:
            break
        if not in_atoms or not stripped:
            continue
        parts = stripped.split()
        if len(parts) < 9:
            continue
        try:
            total += float(parts[-1])
            seen_charge = True
        except Exception:
            continue
    return int(round(total)) if seen_charge else None


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
    chemistry_review = _chemistry_review_payload(cfg)
    input_chemistry = _input_chemistry_payload(cfg)
    strict_chemistry = _to_bool(
        _cfg_get(cfg, "MMGBSA_LIGAND_CHEMISTRY_STRICT", False), default=False
    )
    ligand_state = analyze_ligand_state(sdf, cfg if isinstance(cfg, dict) else {})
    ligand_state_provenance = dict(ligand_state.get("ligand_state_provenance", {}))
    suspicious_chemistry = dict(ligand_state.get("suspicious_chemistry", {}))
    if suspicious_blocks_publication(suspicious_chemistry):
        low_confidence = True
    if rdkit_validate:
        sanitize_error = _rdkit_sanitize_error(sdf, logger)
        if sanitize_error:
            return {
                "ok": False,
                "error": "rdkit_sanitize_failed",
                "detail": sanitize_error,
                "output_mol2": str(mol2_path),
                "output_frcmod": str(frcmod_path),
                "metadata_path": str(metadata_path),
                "charge_method": "",
                "actual_charge_method": "",
                "net_charge_used": None,
                "publication_ready": False,
                "low_confidence": True,
                "ligand_state_provenance": ligand_state_provenance,
                "suspicious_chemistry": suspicious_chemistry,
                "commands_run": [],
                "notes": f"rdkit_sanitize_failed:{sanitize_error}",
            }
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
    requested_primary_method = (
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
    resp_methods = {"resp", "qm", "qm_resp", "precharged"}
    primary_method = requested_primary_method
    if primary_method in resp_methods:
        primary_method = "resp"
    elif primary_method not in {"bcc", "gas"}:
        primary_method = "bcc"
    fallback_method = (
        str(_cfg_get(cfg, "MMGBSA_LIGAND_FALLBACK_CHARGE_METHOD", "gas") or "gas")
        .strip()
        .lower()
    )
    if fallback_method not in {"bcc", "gas"}:
        fallback_method = "gas"

    sqm_level_val = _to_int(_cfg_get(cfg, "MMGBSA_LIGAND_SQM_LEVEL", 2), 2)
    antechamber_timeout = _to_int(
        _cfg_get(cfg, "MMGBSA_LIGAND_ANTECHAMBER_TIMEOUT_SEC", 300),
        300,
    )
    if antechamber_timeout <= 0:
        antechamber_timeout = 300

    outputs_ok = (
        mol2_path.exists()
        and mol2_path.stat().st_size > 0
        and frcmod_path.exists()
        and frcmod_path.stat().st_size > 0
    )
    if outputs_ok and not force:
        existing_metadata = _read_existing_metadata(metadata_path)
        existing_charge_method = str(
            existing_metadata.get("charge_method", primary_method) or primary_method
        ).lower()
        existing_gas_fallback = _to_bool(
            existing_metadata.get("gas_fallback_used", False), default=False
        )
        mol2_warnings, frcmod_warnings = _parameter_warning_summary(
            mol2_path, frcmod_path
        )
        publication_ready = _publication_ready_gate(
            parameterized=True,
            low_confidence=low_confidence,
            input_chemistry=input_chemistry,
            charge_method=existing_charge_method,
            gas_fallback_used=existing_gas_fallback,
            chemistry_review=chemistry_review,
            suspicious_chemistry=suspicious_chemistry,
            mol2_warnings=mol2_warnings,
            frcmod_warnings=frcmod_warnings,
        )
        existing_requested_charge_workflow = str(
            existing_metadata.get("requested_charge_workflow", requested_primary_method)
            or requested_primary_method
        ).lower()
        if (
            existing_requested_charge_workflow in resp_methods
            and existing_charge_method != "resp"
        ):
            publication_ready = False
        strict_block = strict_chemistry and not publication_ready
        ok = not strict_block
        notes = (
            "chemistry_not_publication_ready"
            if strict_block
            else "skipped_existing"
        )
        metadata = {
            "ok": ok,
            "publication_ready": publication_ready,
            "input_sdf": str(sdf),
            "output_mol2": str(mol2_path),
            "output_frcmod": str(frcmod_path),
            "input_format": "sdf",
            "input_chemistry_authoritative": bool(input_chemistry["authoritative"]),
            "input_chemistry_source": input_chemistry["source"],
            "input_chemistry_notes": input_chemistry["notes"],
            "force_field": ligand_at,
            "requested_charge_workflow": existing_metadata.get(
                "requested_charge_workflow", requested_primary_method
            ),
            "charge_method": existing_charge_method,
            "actual_charge_method": existing_metadata.get(
                "actual_charge_method", existing_charge_method
            ),
            "internal_charge_fitting": _to_bool(
                existing_metadata.get("internal_charge_fitting", False),
                default=False,
            ),
            "charge_command_logs": existing_metadata.get("charge_command_logs", []),
            "charge_limitations": existing_metadata.get("charge_limitations", []),
            "net_charge_used": nominal_charge,
            "bcc_attempts": existing_metadata.get("bcc_attempts", []),
            "gas_fallback_used": existing_gas_fallback,
            "rdkit_formal_charge": rdkit_formal_charge,
            "rdkit_radical_electrons": rdkit_radical_electrons,
            "low_confidence": low_confidence,
            "ligand_state_provenance": ligand_state_provenance,
            "suspicious_chemistry": suspicious_chemistry,
            "chemistry_review": chemistry_review,
            "mol2_warnings": mol2_warnings,
            "frcmod_warnings": frcmod_warnings,
            "strict_chemistry": strict_chemistry,
            "notes": notes,
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "skipped_existing": True,
        }
        _write_json_atomic(metadata_path, metadata)
        return {
            "ok": ok,
            "skipped_existing": True,
            "output_mol2": str(mol2_path),
            "output_frcmod": str(frcmod_path),
            "metadata_path": str(metadata_path),
            "charge_method": existing_charge_method,
            "actual_charge_method": existing_metadata.get(
                "actual_charge_method", existing_charge_method
            ),
            "net_charge_used": nominal_charge,
            "publication_ready": publication_ready,
            "low_confidence": low_confidence,
            "ligand_state_provenance": ligand_state_provenance,
            "suspicious_chemistry": suspicious_chemistry,
            "commands_run": [],
            "notes": notes,
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
    charge_command_logs: List[dict] = []
    charge_limitations: List[str] = []
    bcc_attempts: List[dict] = []
    resp_attempts: List[dict] = []
    gas_fallback_used = False
    charge_method_used: Optional[str] = None
    net_charge_used: Optional[int] = None
    resp_mol2_used: Optional[Path] = None
    internal_charge_fitting = False

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
        rc = _run_command_to_log(
            cmd,
            attempt_log,
            logger,
            cwd=work_dir,
            timeout=antechamber_timeout,
        )
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
        elif method in _resp_like_charge_methods():
            resp_attempts.append(
                {
                    "method": method,
                    "nc": charge,
                    "ok": ok,
                    "rc": rc,
                    "log": str(attempt_log),
                }
            )
        charge_command_logs.append(
            {
                "method": method,
                "nc": charge,
                "ok": ok,
                "rc": rc,
                "log": str(attempt_log),
            }
        )
        return ok

    succeeded = False

    if primary_method == "resp":
        resp_mol2 = _find_resp_mol2(sdf)
        if resp_mol2 is not None:
            shutil.copy2(resp_mol2, mol2_path)
            commands_run.append(f"copy_precharged_mol2 {resp_mol2} {mol2_path}")
            succeeded = True
            charge_method_used = "resp"
            net_charge_used = _mol2_total_charge(resp_mol2) or nominal_charge
            resp_mol2_used = resp_mol2
        else:
            _log(
                logger, "INFO", {"ligand": ligand_stem, "sdf": sdf}, "resp_mol2_not_found"
            )
            if _attempt_antechamber("resp", nominal_charge):
                succeeded = True
                charge_method_used = "resp"
                net_charge_used = nominal_charge
                internal_charge_fitting = True
            elif _attempt_antechamber("mul", nominal_charge):
                succeeded = True
                charge_method_used = "mul"
                net_charge_used = nominal_charge
                internal_charge_fitting = True
                charge_limitations.append(
                    "AmberTools -c mul was used after -c resp failed; this is "
                    "QM-derived Mulliken charge assignment, not true RESP, and is "
                    "not publication-ready without independent charge validation."
                )
                _log(
                    logger,
                    "WARNING",
                    {"ligand": ligand_stem, "method": "mul"},
                    "resp_internal_fallback_not_publication_ready",
                )
            else:
                charge_limitations.append(
                    "No precharged RESP/QM MOL2 was found and internal AmberTools "
                    "RESP-like charge generation failed."
                )
                _log(
                    logger,
                    "ERROR",
                    {"ligand": ligand_stem, "sdf": sdf},
                    "resp_internal_charge_failed",
                )

    if primary_method == "bcc" or (
        primary_method == "resp" and not succeeded and fallback_method == "bcc"
    ):
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
        charge_limitations.append(
            "Gas charges were used after higher-quality charge assignment failed; "
            "this is screening-only and is not publication-ready."
        )
        if _attempt_antechamber("gas", nominal_charge):
            succeeded = True
            charge_method_used = "gas"
            net_charge_used = nominal_charge
        else:
            charge_method_used = "gas"
            net_charge_used = nominal_charge

    if primary_method == "resp" and charge_method_used not in {None, "resp", "mul"}:
        charge_limitations.append(
            f"Requested RESP-like workflow fell back to {charge_method_used}; this "
            "is not a RESP/QM charge fit and is not publication-ready for RESP use."
        )

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

    mol2_warnings, frcmod_warnings = _parameter_warning_summary(
        mol2_path, frcmod_path
    )
    publication_ready = _publication_ready_gate(
        parameterized=succeeded and parmchk2_ok,
        low_confidence=low_confidence,
        input_chemistry=input_chemistry,
        charge_method=charge_method_used or primary_method,
        gas_fallback_used=gas_fallback_used,
        chemistry_review=chemistry_review,
        suspicious_chemistry=suspicious_chemistry,
        mol2_warnings=mol2_warnings,
        frcmod_warnings=frcmod_warnings,
    )
    if primary_method == "resp" and charge_method_used != "resp":
        publication_ready = False
    strict_block = strict_chemistry and not publication_ready
    ok = succeeded and parmchk2_ok and not strict_block
    notes = ""
    if not ok:
        if strict_block:
            notes = "chemistry_not_publication_ready"
        else:
            notes = "param_failed" if succeeded else "antechamber_failed"

    metadata = {
        "ok": ok,
        "publication_ready": publication_ready,
        "input_sdf": str(sdf),
        "output_mol2": str(mol2_path),
        "output_frcmod": str(frcmod_path),
        "input_format": "sdf",
        "input_chemistry_authoritative": bool(input_chemistry["authoritative"]),
        "input_chemistry_source": input_chemistry["source"],
        "input_chemistry_notes": input_chemistry["notes"],
        "force_field": ligand_at,
        "requested_charge_workflow": requested_primary_method,
        "charge_method": charge_method_used or primary_method,
        "actual_charge_method": charge_method_used or primary_method,
        "internal_charge_fitting": internal_charge_fitting,
        "charge_command_logs": charge_command_logs,
        "charge_limitations": charge_limitations,
        "resp_mol2_used": str(resp_mol2_used) if resp_mol2_used else "",
        "net_charge_used": int(net_charge_used)
        if net_charge_used is not None
        else nominal_charge,
        "bcc_attempts": bcc_attempts,
        "resp_attempts": resp_attempts,
        "gas_fallback_used": gas_fallback_used,
        "rdkit_formal_charge": rdkit_formal_charge,
        "rdkit_radical_electrons": rdkit_radical_electrons,
        "low_confidence": low_confidence,
        "ligand_state_provenance": ligand_state_provenance,
        "suspicious_chemistry": suspicious_chemistry,
        "chemistry_review": chemistry_review,
        "mol2_warnings": mol2_warnings,
        "frcmod_warnings": frcmod_warnings,
        "strict_chemistry": strict_chemistry,
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
        "actual_charge_method": charge_method_used or primary_method,
        "requested_charge_workflow": requested_primary_method,
        "internal_charge_fitting": internal_charge_fitting,
        "charge_command_logs": charge_command_logs,
        "charge_limitations": charge_limitations,
        "resp_mol2_used": str(resp_mol2_used) if resp_mol2_used else "",
        "net_charge_used": int(net_charge_used)
        if net_charge_used is not None
        else nominal_charge,
        "publication_ready": publication_ready,
        "low_confidence": low_confidence,
        "ligand_state_provenance": ligand_state_provenance,
        "suspicious_chemistry": suspicious_chemistry,
        "commands_run": commands_run,
        "bcc_attempts": bcc_attempts,
        "resp_attempts": resp_attempts,
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
    cfg_local["MMGBSA_INPUT_PH_LABEL"] = ph_label

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
        "publication_ready": result.get("publication_ready"),
        "ligand_state_provenance": result.get("ligand_state_provenance", {}),
        "suspicious_chemistry": result.get("suspicious_chemistry", {}),
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
