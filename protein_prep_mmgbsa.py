from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from installation import load_config
from activesite import (
    load_canonical_cofactors,
    load_canonical_metals,
    load_canonical_waters,
    get_atom_rules,
)
from prep_for_mmgbsa import prep_mmgbsa_from_sdf


_COMPONENT_RECEPTOR = "mmgbsa.receptor"
_COMPONENT_LEAP = "mmgbsa.leap"
_COMPONENT_LEAP_PREP = "mmgbsa.leap.prep"
_NTERM_H_ALLOWLIST = {"H", "H1", "H2", "H3", "HT1", "HT2", "HT3", "HN"}


def _get_logger() -> logging.Logger:
    logger = logging.getLogger(_COMPONENT_RECEPTOR)
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
            if Path(getattr(handler, "baseFilename", "")).resolve() == log_file.resolve():
                return
    log_file.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_file, mode="a", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(fh)


def _log_component(
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


def _log_receptor(logger: logging.Logger, level: str, kvs: Dict[str, object], msg: str | None = None) -> None:
    _log_component(logger, _COMPONENT_RECEPTOR, level, kvs, msg)


def _log_leap(logger: logging.Logger, level: str, kvs: Dict[str, object], msg: str | None = None) -> None:
    _log_component(logger, _COMPONENT_LEAP, level, kvs, msg)


def _log_leap_prep(logger: logging.Logger, level: str, kvs: Dict[str, object], msg: str | None = None) -> None:
    _log_component(logger, _COMPONENT_LEAP_PREP, level, kvs, msg)


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


def _to_float(val: object, default: float) -> float:
    try:
        return float(val)
    except Exception:
        return default


def _resolve_post_docked_context(pdb_path: Path) -> Tuple[str, str, str, str, Path]:
    try:
        resolved = pdb_path.resolve(strict=False)
    except Exception:
        resolved = pdb_path

    parts = resolved.parts
    try:
        idx = parts.index("post_docked")
    except ValueError:
        raise ValueError(
            "pdb_path must live under post_docked/<runid>/<pdb>/<variant>/<pH>/..."
        )

    if len(parts) <= idx + 4:
        raise ValueError("pdb_path does not include runid/pdb/variant/pH components")

    runid = parts[idx + 1]
    pdb = parts[idx + 2]
    variant = parts[idx + 3]
    ph_label = parts[idx + 4]
    base_dir = Path(*parts[: idx + 5])
    return runid, pdb, variant, ph_label, base_dir


def _parse_center(text: str) -> Tuple[float, float, float]:
    tokens = [t.strip() for t in text.split(",") if t.strip()]
    if len(tokens) != 3:
        raise ValueError("center must be provided as 'x,y,z'")
    return float(tokens[0]), float(tokens[1]), float(tokens[2])


def _parse_token_list(value: object) -> List[str]:
    if value is None:
        return []
    text = str(value).strip()
    if not text:
        return []
    tokens: List[str] = []
    for chunk in text.replace(",", " ").split():
        token = chunk.strip().upper()
        if token:
            tokens.append(token)
    return tokens


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
    os.replace(tmp_path, path)


def _residue_key(line: str) -> str:
    chain = (line[21:22] or "-").strip() or "-"
    resseq = (line[22:26] or "0").strip() or "0"
    icode = (line[26:27] or "-").strip() or "-"
    return f"{chain}:{resseq}:{icode}"


def _parse_coords(line: str) -> Optional[Tuple[float, float, float]]:
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except Exception:
        return None
    return x, y, z


def _distance_sq(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return dx * dx + dy * dy + dz * dz


def _sanitize_receptor_lines(
    kept_lines: List[str],
    strip_nterm_h: bool,
    insert_ter_on_chainbreak: bool,
    chainbreak_cn_max_a: float,
) -> Tuple[List[str], Dict[str, object]]:
    if not kept_lines or not (strip_nterm_h or insert_ter_on_chainbreak):
        return kept_lines, {
            "inserted_TER_count": 0,
            "nterm_h_stripped_count": 0,
            "nterm_h_stripped_residues": [],
        }

    residue_first_idx: Dict[str, int] = {}
    residue_chain_order: Dict[str, List[str]] = {}
    residue_coords: Dict[str, Dict[str, Tuple[float, float, float]]] = {}

    for idx, line in enumerate(kept_lines):
        if not line.startswith("ATOM"):
            continue
        key = _residue_key(line)
        chain = (line[21:22] or "-").strip() or "-"
        if key not in residue_first_idx:
            residue_first_idx[key] = idx
            residue_chain_order.setdefault(chain, []).append(key)

        atom_name = line[12:16].strip().upper()
        if atom_name in {"C", "N"}:
            coords = _parse_coords(line)
            if coords is not None:
                residue_coords.setdefault(key, {})[atom_name] = coords

    break_before: set[str] = set()
    if insert_ter_on_chainbreak:
        threshold_sq = float(chainbreak_cn_max_a) ** 2
        for chain, keys in residue_chain_order.items():
            for idx in range(1, len(keys)):
                prev_key = keys[idx - 1]
                next_key = keys[idx]
                prev_coords = residue_coords.get(prev_key, {})
                next_coords = residue_coords.get(next_key, {})
                if "C" in prev_coords and "N" in next_coords:
                    if _distance_sq(prev_coords["C"], next_coords["N"]) > threshold_sq:
                        break_before.add(next_key)

    segment_start_keys: set[str] = set()
    if strip_nterm_h:
        for chain, keys in residue_chain_order.items():
            if keys:
                segment_start_keys.add(keys[0])
        pending_segment_start = False
        for line in kept_lines:
            if line.startswith("TER"):
                pending_segment_start = True
                continue
            if pending_segment_start and line.startswith("ATOM"):
                segment_start_keys.add(_residue_key(line))
                pending_segment_start = False
        if break_before:
            segment_start_keys |= break_before

    sanitized_lines: List[str] = []
    inserted_ter_count = 0
    nterm_h_stripped_count = 0
    stripped_residues: set[str] = set()
    seen_residues: set[str] = set()

    for line in kept_lines:
        if line.startswith("ATOM"):
            key = _residue_key(line)
            if key not in seen_residues:
                if insert_ter_on_chainbreak and key in break_before:
                    if not sanitized_lines or not sanitized_lines[-1].startswith("TER"):
                        sanitized_lines.append("TER\n")
                        inserted_ter_count += 1
                seen_residues.add(key)

            if strip_nterm_h and key in segment_start_keys:
                atom_name = line[12:16].strip().upper()
                if atom_name in _NTERM_H_ALLOWLIST:
                    element = line[76:78].strip().upper()
                    if (element and element == "H") or (not element and atom_name.startswith("H")):
                        nterm_h_stripped_count += 1
                        stripped_residues.add(key)
                        continue

        sanitized_lines.append(line)

    return sanitized_lines, {
        "inserted_TER_count": inserted_ter_count,
        "nterm_h_stripped_count": nterm_h_stripped_count,
        "nterm_h_stripped_residues": sorted(stripped_residues),
    }


def _strip_receptor_h_for_leap(receptor_path: Path, strip_all_h: bool) -> Tuple[Path, Dict[str, object]]:
    if not strip_all_h:
        hoh_residues: set[str] = set()
        with receptor_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                resn_raw = line[17:20].strip().upper()
                if resn_raw == "HOH":
                    hoh_residues.add(_residue_key(line))
        return receptor_path, {
            "strip_all_h": False,
            "removed_h": 0,
            "hoh_residue_count": len(hoh_residues),
        }

    out_path = receptor_path.with_name(f"{receptor_path.stem}.noH{receptor_path.suffix}")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    kept_lines: List[str] = []
    hoh_residues: set[str] = set()
    removed_h = 0

    with receptor_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resn_raw = line[17:20].strip().upper()
                if resn_raw == "HOH":
                    hoh_residues.add(_residue_key(line))

                atom_name = line[12:16].strip().upper()
                element = line[76:78].strip().upper()
                if element == "H" or atom_name.startswith("H"):
                    removed_h += 1
                    continue

            kept_lines.append(line)

    _write_text_atomic(out_path, "".join(kept_lines))

    return out_path, {
        "strip_all_h": True,
        "removed_h": removed_h,
        "hoh_residue_count": len(hoh_residues),
    }


def _leap_water_lines(water_model: str, map_hoh_to_wat: bool, hoh_residue_count: int) -> List[str]:
    lines: List[str] = []
    if water_model:
        lines.append(f"source leaprc.water.{water_model}")
    if map_hoh_to_wat and hoh_residue_count > 0:
        lines.append('addPdbResMap { { "WAT" "HOH" } }')
    return lines


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


def _resolve_ambertools_prefix(cfg: Optional[object]) -> Optional[str]:
    env_prefix = os.environ.get("AMBERTOOLS_PREFIX")
    if env_prefix:
        return env_prefix
    if cfg is None:
        return None
    try:
        cfg_prefix = cfg.get("MMGBSA_AMBERTOOLS_PREFIX") or cfg.get("AMBERTOOLS_PREFIX")
    except Exception:
        cfg_prefix = getattr(cfg, "MMGBSA_AMBERTOOLS_PREFIX", None) or getattr(cfg, "AMBERTOOLS_PREFIX", None)
    if cfg_prefix:
        return str(cfg_prefix)
    return None


def _select_tleap_runner(logger: logging.Logger) -> Tuple[List[str], str, str]:
    cfg = load_config()
    prefix_value = _resolve_ambertools_prefix(cfg)
    if prefix_value:
        prefix = Path(prefix_value).expanduser()
        if _looks_like_pkgs_cache(prefix):
            _log_leap(
                logger,
                "WARNING",
                {"reason": "pkgs_cache_prefix", "fallback": "PATH_or_common_prefix"},
                "ignoring_AMBERTOOLS_PREFIX",
            )
        elif _prefix_has_tool(prefix, "tleap"):
            micromamba = _resolve_micromamba()
            if not micromamba:
                raise FileNotFoundError("micromamba not found; cannot run tleap from prefix")
            runner = [str(micromamba), "run", "-p", str(prefix)]
            return runner, "tleap", f"prefix:{prefix}"
        else:
            _log_leap(
                logger,
                "WARNING",
                {"reason": "missing_tleap", "fallback": "PATH_or_common_prefix"},
                "ignoring_AMBERTOOLS_PREFIX",
            )

    tleap_path = shutil.which("tleap")
    if tleap_path:
        return [], tleap_path, "PATH"

    common_prefixes = [
        Path("/home/atlas/micromamba/envs/AmberTools25"),
        Path("/home/michael/atlas/tools/envs/ambertools"),
        Path("/home/michael/atlas/micromamba/envs/AmberTools25"),
    ]
    for prefix in common_prefixes:
        if _prefix_has_tool(prefix, "tleap") and not _looks_like_pkgs_cache(prefix):
            micromamba = _resolve_micromamba()
            if not micromamba:
                raise FileNotFoundError("micromamba not found; cannot run tleap from prefix")
            runner = [str(micromamba), "run", "-p", str(prefix)]
            return runner, "tleap", f"prefix:{prefix}"

    raise FileNotFoundError("could not locate tleap; set AMBERTOOLS_PREFIX or PATH")


def _run_tleap_impl(leap_file_path: str, work_dir: str, force: bool = False) -> dict:
    logger = _get_logger()
    leap_path = Path(leap_file_path)
    work_path = Path(work_dir)

    if not leap_path.exists():
        raise FileNotFoundError(f"tleap input not found: {leap_path}")

    work_path.mkdir(parents=True, exist_ok=True)

    runner, tleap_bin, source = _select_tleap_runner(logger)
    leap_arg = leap_path.name if leap_path.parent == work_path else str(leap_path)
    cmd = runner + [tleap_bin, "-f", leap_arg]
    log_path = work_path / "tleap.log"

    _log_leap(
        logger,
        "INFO",
        {"tleap": " ".join(cmd), "source": source, "work_dir": work_path, "force": force},
        "run",
    )

    with log_path.open("w", encoding="utf-8") as handle:
        proc = subprocess.run(cmd, cwd=str(work_path), stdout=handle, stderr=subprocess.STDOUT)

    result = {
        "args": " ".join(cmd),
        "returncode": proc.returncode,
        "log_path": str(log_path),
        "ok": proc.returncode == 0,
    }

    if proc.returncode != 0:
        _log_leap(logger, "ERROR", {"returncode": proc.returncode, "log": log_path}, "tleap_failed")
        raise RuntimeError(f"tleap failed (code={proc.returncode}); see {log_path}")

    _log_leap(logger, "INFO", {"ok": "true", "log": log_path}, "tleap_done")
    return result


def run_tleap(leap_file_path: str, work_dir: str, force: bool = False) -> dict:
    return _run_tleap_impl(leap_file_path=leap_file_path, work_dir=work_dir, force=force)


def write_leap_receptor_only(
    receptor_pdb_path: str,
    out_dir: str,
    force: bool = False,
    water_model: str = "tip3p",
    map_hoh_to_wat: bool = False,
    hoh_residue_count: int = 0,
) -> dict:
    logger = _get_logger()
    receptor_path = Path(receptor_pdb_path)
    out_path = Path(out_dir)

    if not receptor_path.exists():
        raise FileNotFoundError(f"receptor PDB not found: {receptor_path}")

    out_path.mkdir(parents=True, exist_ok=True)

    leap_file = out_path / "build_receptor.leap"
    prmtop_path = out_path / "receptor.prmtop"
    inpcrd_path = out_path / "receptor.inpcrd"

    receptor_rel = os.path.relpath(receptor_path, out_path)
    leap_lines = ["source leaprc.protein.ff14SB"]
    leap_lines.extend(_leap_water_lines(water_model=water_model, map_hoh_to_wat=map_hoh_to_wat, hoh_residue_count=hoh_residue_count))
    leap_lines.extend(
        [
            f"REC = loadpdb {receptor_rel}",
            "saveamberparm REC receptor.prmtop receptor.inpcrd",
            "quit",
            "",
        ]
    )
    leap_text = "\n".join(leap_lines)
    _write_text_atomic(leap_file, leap_text)

    outputs_ok = (
        prmtop_path.exists()
        and prmtop_path.stat().st_size > 0
        and inpcrd_path.exists()
        and inpcrd_path.stat().st_size > 0
    )
    skip_tleap = outputs_ok and not force

    return {
        "leap_file": str(leap_file),
        "prmtop_path": str(prmtop_path),
        "inpcrd_path": str(inpcrd_path),
        "skip_tleap": skip_tleap,
    }


def write_leap_for_ligand(
    receptor_pdb_path: str,
    ligand_mol2: str,
    ligand_frcmod: str,
    out_dir: str,
    force: bool = False,
    water_model: str = "tip3p",
    map_hoh_to_wat: bool = False,
    hoh_residue_count: int = 0,
) -> dict:
    logger = _get_logger()
    receptor_path = Path(receptor_pdb_path)
    mol2_path = Path(ligand_mol2)
    frcmod_path = Path(ligand_frcmod)
    out_path = Path(out_dir)

    if not receptor_path.exists():
        raise FileNotFoundError(f"receptor PDB not found: {receptor_path}")
    if not mol2_path.exists():
        raise FileNotFoundError(f"ligand mol2 not found: {mol2_path}")
    if not frcmod_path.exists():
        raise FileNotFoundError(f"ligand frcmod not found: {frcmod_path}")

    out_path.mkdir(parents=True, exist_ok=True)

    leap_file = out_path / "build.leap"
    ligand_prmtop = out_path / "ligand.prmtop"
    ligand_inpcrd = out_path / "ligand.inpcrd"
    receptor_prmtop = out_path / "receptor.prmtop"
    receptor_inpcrd = out_path / "receptor.inpcrd"
    complex_prmtop = out_path / "complex.prmtop"
    complex_inpcrd = out_path / "complex.inpcrd"

    receptor_rel = os.path.relpath(receptor_path, out_path)
    mol2_rel = os.path.relpath(mol2_path, out_path)
    frcmod_rel = os.path.relpath(frcmod_path, out_path)

    leap_lines = [
        "source leaprc.protein.ff14SB",
        "source leaprc.gaff2",
    ]
    leap_lines.extend(_leap_water_lines(water_model=water_model, map_hoh_to_wat=map_hoh_to_wat, hoh_residue_count=hoh_residue_count))
    leap_lines.extend(
        [
            f"loadamberparams {frcmod_rel}",
            f"LIG = loadmol2 {mol2_rel}",
            f"REC = loadpdb {receptor_rel}",
            "COM = combine { REC LIG }",
            "saveamberparm REC receptor.prmtop receptor.inpcrd",
            "saveamberparm LIG ligand.prmtop ligand.inpcrd",
            "saveamberparm COM complex.prmtop complex.inpcrd",
            "quit",
            "",
        ]
    )
    leap_text = "\n".join(leap_lines)
    _write_text_atomic(leap_file, leap_text)

    outputs_ok = complex_prmtop.exists() and complex_prmtop.stat().st_size > 0
    skip_tleap = outputs_ok and not force

    return {
        "leap_file": str(leap_file),
        "receptor_prmtop": str(receptor_prmtop),
        "receptor_inpcrd": str(receptor_inpcrd),
        "ligand_prmtop": str(ligand_prmtop),
        "ligand_inpcrd": str(ligand_inpcrd),
        "complex_prmtop": str(complex_prmtop),
        "complex_inpcrd": str(complex_inpcrd),
        "skip_tleap": skip_tleap,
    }


def write_leap_for_ligands(
    receptor_pdb_path: str,
    ligand_mol2_frcmod_pairs: List[dict],
    out_dir_base: str,
    force: bool = False,
    water_model: str = "tip3p",
    map_hoh_to_wat: bool = False,
    hoh_residue_count: int = 0,
) -> List[dict]:
    results: List[dict] = []
    base_dir = Path(out_dir_base)
    base_dir.mkdir(parents=True, exist_ok=True)

    for entry in ligand_mol2_frcmod_pairs:
        stage_dir = str(entry.get("stage_dir", "")).strip() or "stage"
        ligand_base = str(entry.get("ligand_base", "")).strip()
        mol2_path = str(entry.get("mol2_path", "")).strip()
        frcmod_path = str(entry.get("frcmod_path", "")).strip()
        if not ligand_base or not mol2_path or not frcmod_path:
            raise ValueError(f"invalid ligand entry for leap generation: {entry}")

        ligand_out_dir = base_dir / stage_dir / ligand_base
        topo = write_leap_for_ligand(
            receptor_pdb_path=receptor_pdb_path,
            ligand_mol2=mol2_path,
            ligand_frcmod=frcmod_path,
            out_dir=str(ligand_out_dir),
            force=force,
            water_model=water_model,
            map_hoh_to_wat=map_hoh_to_wat,
            hoh_residue_count=hoh_residue_count,
        )
        topo.update(
            {
                "stage_dir": stage_dir,
                "ligand_base": ligand_base,
                "output_dir": str(ligand_out_dir),
            }
        )
        results.append(topo)

    return results


def prep_mmgbsa_receptor(
    pdb_path: str,
    runid: str,
    center: Tuple[float, float, float],
    radius: float | None,
    force: bool = False,
) -> dict:
    logger = _get_logger()
    path = Path(pdb_path)

    if not path.exists():
        raise FileNotFoundError(f"pdb_path not found: {path}")

    _log_receptor(logger, "INFO", {"input": path})

    inferred_runid, pdb_id, variant, ph_label, base_dir = _resolve_post_docked_context(path)
    if runid and runid != inferred_runid:
        raise ValueError(f"runid mismatch: arg={runid} path={inferred_runid}")
    runid = runid or inferred_runid

    out_dir = base_dir / "mmgbsa_receptor"
    out_path = out_dir / f"{pdb_id}_mmgbsa_receptor.pdb"
    out_dir.mkdir(parents=True, exist_ok=True)
    _ensure_file_handler(logger, out_dir / "mmgbsa_receptor_prep.log")

    _log_receptor(logger, "INFO", {"runid": runid, "pdb": pdb_id, "variant": variant, "ph": ph_label})
    _log_receptor(logger, "INFO", {"output": out_path})

    cfg = load_config()
    file_cfg = cfg.get("_FILE_CFG", {}) if isinstance(cfg, dict) else {}

    receptor_force = _to_bool(cfg.get("MMGBSA_RECEPTOR_FORCE", False), default=False)
    if receptor_force:
        force = True

    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        _log_receptor(logger, "INFO", {"skipped": True}, "output_present")
        _log_receptor(logger, "INFO", {"ok": "true"}, "done")
        return {
            "pdb_path": str(path),
            "output_path": str(out_path),
            "skipped": True,
            "metal_lines_removed": 0,
            "water_residues_total": 0,
            "water_residues_kept": 0,
            "water_residues_removed": 0,
            "water_lines_removed": 0,
            "water_lines_kept": 0,
        }

    keep_metals = _to_bool(cfg.get("MMGBSA_KEEP_METALS", False), default=False)
    if "MMGBSA_KEEP_METALS" not in file_cfg and "MMGBSA_STRIP_METALS" in file_cfg:
        keep_metals = not _to_bool(cfg.get("MMGBSA_STRIP_METALS", True), default=True)
    strip_metals = not keep_metals

    keep_waters = _to_bool(cfg.get("MMGBSA_KEEP_WATERS", True), default=True)
    water_policy_raw = str(cfg.get("MMGBSA_WATER_POLICY", "ACTIVE_SITE") or "ACTIVE_SITE")
    water_policy = water_policy_raw.strip().upper().replace("-", "_").replace(" ", "_")
    if "MMGBSA_KEEP_WATERS" in file_cfg or "MMGBSA_WATER_POLICY" not in file_cfg:
        water_policy = "ACTIVE_SITE" if keep_waters else "NONE"

    water_radius_cfg = cfg.get("MMGBSA_WATER_KEEP_RADIUS", cfg.get("MMGBSA_WATER_KEEP_RADIUS_A", 6.0))
    water_radius = _to_float(water_radius_cfg, 6.0)
    if radius is not None:
        water_radius = float(radius)

    if water_policy not in {"NONE", "ACTIVE_SITE", "ALL"}:
        raise ValueError(f"MMGBSA_WATER_POLICY must be NONE, ACTIVE_SITE, or ALL (got {water_policy_raw})")

    if water_policy == "ACTIVE_SITE" and center is None:
        raise ValueError("center is required when MMGBSA_WATER_POLICY=ACTIVE_SITE")

    _log_receptor(
        logger,
        "INFO",
        {
            "strip_metals": strip_metals,
            "keep_metals": keep_metals,
            "keep_waters": keep_waters,
            "water_policy": water_policy,
            "water_radius": water_radius,
        },
        "policy",
    )

    water_override = _parse_token_list(cfg.get("MMGBSA_WATER_RETAIN_TOKENS", ""))
    metal_override = _parse_token_list(cfg.get("MMGBSA_METAL_RETAIN_TOKENS", ""))

    use_water_aliases = _to_bool(cfg.get("MMGBSA_WATER_USE_ALIASES", True), default=True)
    use_metal_aliases = _to_bool(cfg.get("MMGBSA_METAL_USE_ALIASES", True), default=True)

    rules = get_atom_rules()
    metals_set = {rules.normalize_resname(t) or t for t in metal_override} if metal_override else load_canonical_metals(None)
    waters_set = {rules.normalize_resname(t) or t for t in water_override} if water_override else load_canonical_waters(None)
    cofactors_set = load_canonical_cofactors(None)

    if not use_water_aliases or not use_metal_aliases:
        _log_receptor(
            logger,
            "WARNING",
            {
                "reason": "aliases_required",
                "water_aliases": use_water_aliases,
                "metal_aliases": use_metal_aliases,
            },
            "using_alias_sets",
        )

    _log_receptor(
        logger,
        "INFO",
        {
            "metals": len(metals_set),
            "waters": len(waters_set),
            "cofactors": len(cofactors_set),
            "water_override": bool(water_override),
            "metal_override": bool(metal_override),
        },
        "aliases_loaded",
    )

    lines: List[str] = []
    water_info: Dict[str, Dict[str, Optional[Tuple[float, float, float]]]] = {}

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            lines.append(line)
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            resn_raw = line[17:20].strip().upper()
            if not resn_raw:
                continue
            resn_norm = rules.normalize_resname(resn_raw) or resn_raw
            if resn_norm not in waters_set:
                continue

            key = _residue_key(line)
            info = water_info.setdefault(key, {"oxygen": None})
            atom_name = line[12:16].strip().upper()
            if atom_name.startswith("O"):
                coords = _parse_coords(line)
                if coords is not None:
                    info["oxygen"] = coords

    water_keys = set(water_info.keys())
    water_keep_keys: set[str] = set()
    water_remove_keys: set[str] = set()
    missing_oxygen_keys: set[str] = set()

    if water_policy == "ALL":
        water_keep_keys = set(water_keys)
    elif water_policy == "NONE":
        water_remove_keys = set(water_keys)
    else:
        center_vec = center
        radius_sq = float(water_radius) ** 2
        for key, info in water_info.items():
            coords = info.get("oxygen")
            if coords is None:
                missing_oxygen_keys.add(key)
                water_remove_keys.add(key)
                continue
            if _distance_sq(coords, center_vec) <= radius_sq:
                water_keep_keys.add(key)
            else:
                water_remove_keys.add(key)

    if missing_oxygen_keys:
        _log_receptor(
            logger,
            "WARNING",
            {"reason": "water_missing_oxygen", "count": len(missing_oxygen_keys), "action": "removed"},
        )

    kept_lines: List[str] = []
    metal_lines_removed = 0
    water_lines_removed = 0
    water_lines_kept = 0

    for line in lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            kept_lines.append(line)
            continue

        resn_raw = line[17:20].strip().upper()
        if not resn_raw:
            kept_lines.append(line)
            continue
        resn_norm = rules.normalize_resname(resn_raw) or resn_raw

        if strip_metals:
            element_raw = line[76:78].strip().upper()
            element_norm = rules.normalize_resname(element_raw) if element_raw else element_raw
            if resn_norm in metals_set or (element_norm and element_norm in metals_set):
                metal_lines_removed += 1
                continue

        if water_policy != "ALL" and resn_norm in waters_set:
            key = _residue_key(line)
            if key in water_remove_keys:
                water_lines_removed += 1
                continue
            water_lines_kept += 1

        kept_lines.append(line)

    strip_nterm_h = _to_bool(cfg.get("MMGBSA_TLEAP_STRIP_NTERM_H", True), default=True)
    insert_ter_on_break = _to_bool(cfg.get("MMGBSA_TLEAP_INSERT_TER_ON_CHAINBREAK", True), default=True)
    chainbreak_cn_max = _to_float(cfg.get("MMGBSA_TLEAP_CHAINBREAK_CN_MAX_A", 2.2), 2.2)

    if strip_nterm_h or insert_ter_on_break:
        kept_lines, sanitize_info = _sanitize_receptor_lines(
            kept_lines=kept_lines,
            strip_nterm_h=strip_nterm_h,
            insert_ter_on_chainbreak=insert_ter_on_break,
            chainbreak_cn_max_a=chainbreak_cn_max,
        )
        log_kvs = {
            "inserted_TER_count": sanitize_info.get("inserted_TER_count", 0),
            "nterm_h_stripped_count": sanitize_info.get("nterm_h_stripped_count", 0),
        }
        stripped_res = sanitize_info.get("nterm_h_stripped_residues") or []
        if stripped_res:
            shown = stripped_res[:10]
            extra = len(stripped_res) - len(shown)
            if extra > 0:
                log_kvs["nterm_h_stripped_residues"] = ",".join(shown) + f",+{extra} more"
            else:
                log_kvs["nterm_h_stripped_residues"] = ",".join(shown)
        _log_receptor(logger, "INFO", log_kvs, "sanitize")

    tmp_path = _tmp_path(out_path)
    with tmp_path.open("w", encoding="utf-8") as handle:
        handle.writelines(kept_lines)
    if not tmp_path.exists() or tmp_path.stat().st_size == 0:
        if tmp_path.exists():
            tmp_path.unlink()
        raise RuntimeError(f"output not created: {tmp_path}")
    os.replace(tmp_path, out_path)

    _log_receptor(
        logger,
        "INFO",
        {
            "metal_lines_removed": metal_lines_removed,
            "water_residues_total": len(water_keys),
            "water_residues_kept": len(water_keep_keys),
            "water_residues_removed": len(water_remove_keys),
            "water_lines_removed": water_lines_removed,
            "water_lines_kept": water_lines_kept,
        },
        "counts",
    )
    _log_receptor(logger, "INFO", {"ok": "true"}, "done")

    return {
        "pdb_path": str(path),
        "output_path": str(out_path),
        "metal_lines_removed": metal_lines_removed,
        "water_residues_total": len(water_keys),
        "water_residues_kept": len(water_keep_keys),
        "water_residues_removed": len(water_remove_keys),
        "water_lines_removed": water_lines_removed,
        "water_lines_kept": water_lines_kept,
    }


def prep_mmgbsa_receptor_and_topologies(
    pdb_path: str,
    runid: str,
    center: Tuple[float, float, float],
    radius: float | None,
    ligand_sdf_paths: List[str],
    force: bool = False,
    run_tleap: bool = True,
) -> dict:
    logger = _get_logger()
    cfg = load_config()

    topo_enabled = _to_bool(cfg.get("MMGBSA_TLEAP_ENABLED", cfg.get("MMGBSA_TOPOLOGY_PREP_ENABLED", True)), default=True)
    if not topo_enabled:
        _log_leap(logger, "INFO", {"enabled": False}, "topology_prep_disabled")
        receptor_result = prep_mmgbsa_receptor(pdb_path, runid, center, radius, force=force)
        return {
            "receptor": receptor_result,
            "topology_enabled": False,
            "receptor_topology": None,
            "ligand_topologies": [],
        }

    strip_all_h_for_leap = _to_bool(cfg.get("MMGBSA_TLEAP_STRIP_ALL_H", True), default=True)
    water_model = str(cfg.get("MMGBSA_TLEAP_WATER_MODEL", "tip3p") or "tip3p").strip().lower()
    map_hoh_to_wat = _to_bool(cfg.get("MMGBSA_TLEAP_MAP_HOH_TO_WAT", True), default=True)
    if water_model not in {"tip3p"}:
        raise ValueError(f"MMGBSA_TLEAP_WATER_MODEL supports tip3p only (got {water_model})")

    inferred_runid, pdb_id, variant, ph_label, base_dir = _resolve_post_docked_context(Path(pdb_path))
    if runid and runid != inferred_runid:
        raise ValueError(f"runid mismatch: arg={runid} path={inferred_runid}")

    topology_dirname = str(cfg.get("MMGBSA_TOPOLOGY_DIRNAME", "mmgbsa_topologies") or "mmgbsa_topologies")
    top_root = base_dir / topology_dirname

    if _to_bool(cfg.get("MMGBSA_TLEAP_FORCE", False), default=False):
        force = True

    receptor_result = prep_mmgbsa_receptor(pdb_path, runid, center, radius, force=force)
    receptor_pdb = receptor_result["output_path"]
    receptor_for_leap_path, strip_info = _strip_receptor_h_for_leap(Path(receptor_pdb), strip_all_h_for_leap)
    hoh_residue_count = strip_info.get("hoh_residue_count", 0)
    _log_leap_prep(
        logger,
        "INFO",
        {
            "strip_all_h": strip_all_h_for_leap,
            "input": receptor_pdb,
            "output": str(receptor_for_leap_path),
            "removed_H": strip_info.get("removed_h", 0),
            "hoh_residues": hoh_residue_count,
        },
    )
    _log_leap_prep(
        logger,
        "INFO",
        {"map_hoh_to_wat": map_hoh_to_wat, "water_model": water_model, "hoh_residues": hoh_residue_count},
        "water_mapping",
    )

    if not _to_bool(cfg.get("MMGBSA_TLEAP_RUN", True), default=True):
        run_tleap = False

    _log_leap(
        logger,
        "INFO",
        {
            "runid": inferred_runid,
            "pdb": pdb_id,
            "variant": variant,
            "ph": ph_label,
            "output_dir": top_root,
            "tleap_run": run_tleap,
        },
        "topology_prep",
    )

    receptor_top_dir = top_root / "receptor"
    receptor_topo = write_leap_receptor_only(
        str(receptor_for_leap_path),
        str(receptor_top_dir),
        force=force,
        water_model=water_model,
        map_hoh_to_wat=map_hoh_to_wat,
        hoh_residue_count=hoh_residue_count,
    )
    receptor_leap_log = receptor_top_dir / "tleap.log"
    _log_leap(
        logger,
        "INFO",
        {
            "output_dir": receptor_top_dir,
            "receptor_pdb_used": receptor_for_leap_path,
            "build_leap_path": receptor_topo["leap_file"],
            "tleap_log_path": receptor_leap_log,
            "tleap_run": run_tleap,
            "skip": receptor_topo["skip_tleap"],
        },
        "receptor_topology",
    )

    receptor_tleap_result = None
    if run_tleap and not receptor_topo["skip_tleap"]:
        receptor_tleap_result = _run_tleap_impl(
            leap_file_path=receptor_topo["leap_file"],
            work_dir=str(receptor_top_dir),
            force=force,
        )

    ligand_results: List[dict] = []

    for ligand_sdf in ligand_sdf_paths:
        sdf_path = Path(ligand_sdf)
        if not sdf_path.exists():
            raise FileNotFoundError(f"ligand sdf not found: {sdf_path}")

        lig_runid, lig_pdb, lig_variant, lig_ph, lig_base_dir = _resolve_post_docked_context(sdf_path)
        if lig_base_dir != base_dir:
            raise ValueError(f"ligand path does not match receptor base dir: {sdf_path}")
        if lig_runid != inferred_runid or lig_pdb != pdb_id or lig_variant != variant or lig_ph != ph_label:
            raise ValueError(f"ligand path context mismatch for {sdf_path}")

        stage_dir = sdf_path.parent.name
        ligand_base = sdf_path.stem

        prep_result = prep_mmgbsa_from_sdf(
            sdf_path=str(sdf_path),
            net_charge=None,
            amber_prefix=None,
            dry_run=False,
            force=force,
        )

        mol2_path = prep_result["mol2_path"]
        frcmod_path = prep_result["frcmod_path"]

        ligand_out_dir = top_root / stage_dir / ligand_base
        ligand_topo = write_leap_for_ligand(
            receptor_pdb_path=str(receptor_for_leap_path),
            ligand_mol2=mol2_path,
            ligand_frcmod=frcmod_path,
            out_dir=str(ligand_out_dir),
            force=force,
            water_model=water_model,
            map_hoh_to_wat=map_hoh_to_wat,
            hoh_residue_count=hoh_residue_count,
        )
        ligand_leap_log = ligand_out_dir / "tleap.log"

        _log_leap(
            logger,
            "INFO",
            {
                "stage_dir": stage_dir,
                "ligand": ligand_base,
                "output_dir": ligand_out_dir,
                "receptor_pdb_used": receptor_for_leap_path,
                "build_leap_path": ligand_topo["leap_file"],
                "tleap_log_path": ligand_leap_log,
                "tleap_run": run_tleap,
                "skip": ligand_topo["skip_tleap"],
            },
            "ligand_topology",
        )

        ligand_tleap_result = None
        if run_tleap and not ligand_topo["skip_tleap"]:
            ligand_tleap_result = _run_tleap_impl(
                leap_file_path=ligand_topo["leap_file"],
                work_dir=str(ligand_out_dir),
                force=force,
            )

        ligand_results.append(
            {
                "sdf_path": str(sdf_path),
                "stage_dir": stage_dir,
                "ligand_base": ligand_base,
                "prep_result": prep_result,
                "topology": ligand_topo,
                "tleap": ligand_tleap_result,
            }
        )

    _log_leap(logger, "INFO", {"ok": "true"}, "topology_done")

    return {
        "receptor": receptor_result,
        "topology_enabled": True,
        "topology_root": str(top_root),
        "receptor_topology": receptor_topo,
        "receptor_tleap": receptor_tleap_result,
        "ligand_topologies": ligand_results,
    }


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Prepare an MMGBSA receptor PDB from a post_docked receptor.")
    parser.add_argument("--pdb", required=True, help="Input receptor PDB under post_docked/<runid>/<pdb>/<variant>/<pH>/." )
    parser.add_argument("--runid", default="", help="Run ID for validation if path contains post_docked.")
    parser.add_argument("--center", required=True, help="Active-site center as 'x,y,z'.")
    parser.add_argument("--radius", type=float, default=None, help="Active-site radius (A). Overrides config.")
    parser.add_argument("--ligand-sdf", action="append", default=[], help="Ligand SDF under post_docked/<runid>/<pdb>/<variant>/<pH>/<stage_dir>/. Can be repeated.")
    parser.add_argument("--no-topology-prep", action="store_true", help="Skip MMGBSA topology prep.")
    parser.add_argument("--no-tleap", action="store_true", help="Write tleap scripts but do not run tleap.")
    parser.add_argument("--force", action="store_true", help="Overwrite outputs if present.")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_arg_parser()
    args = parser.parse_args(argv)
    logger = _get_logger()

    try:
        center = _parse_center(args.center)
        cfg = load_config()

        topo_enabled = _to_bool(cfg.get("MMGBSA_TOPOLOGY_PREP_ENABLED", True), default=True)
        if args.no_topology_prep:
            topo_enabled = False

        run_tleap_flag = _to_bool(cfg.get("MMGBSA_TLEAP_RUN", True), default=True)
        if args.no_tleap:
            run_tleap_flag = False

        if topo_enabled:
            prep_mmgbsa_receptor_and_topologies(
                pdb_path=args.pdb,
                runid=args.runid,
                center=center,
                radius=args.radius,
                ligand_sdf_paths=list(args.ligand_sdf or []),
                force=args.force,
                run_tleap=run_tleap_flag,
            )
        else:
            _log_leap(logger, "INFO", {"enabled": False}, "topology_prep_disabled")
            prep_mmgbsa_receptor(
                pdb_path=args.pdb,
                runid=args.runid,
                center=center,
                radius=args.radius,
                force=args.force,
            )

    except Exception as exc:
        _log_receptor(logger, "ERROR", {"reason": type(exc).__name__, "detail": str(exc)}, "failed")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
