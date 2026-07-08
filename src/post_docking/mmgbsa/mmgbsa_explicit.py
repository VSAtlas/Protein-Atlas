from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from config.tool_resolver import select_ambertools_tool
from post_docking.mmgbsa._atomic_io import write_json_atomic, write_text_atomic
from post_docking.mmgbsa.mmgbsa_mpi_env import (
    apply_blas_single_thread_env_defaults,
    resolve_md_mpi_ranks,
    select_mpi_launcher,
)
from post_docking.mmgbsa.mmgbsa_openmm import run_openmm_step
from protein_prep.pdb_records import line_xyz as _pdb_xyz

def _cfg_get(cfg: object | None, key: str, default: object = None) -> object:
    if cfg is None:
        return default
    if isinstance(cfg, Mapping):
        return cfg.get(key, default)
    getter = getattr(cfg, "get", None)
    if callable(getter):
        try:
            return getter(key, default)
        except Exception:
            pass
    return getattr(cfg, key, default)


def _to_bool(value: object, default: bool = False) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "no", "n", "off"}:
        return False
    return default


def _to_float(value: object, default: float) -> float:
    try:
        if isinstance(value, (int, float, str, bytes, bytearray)):
            return float(value)
    except (TypeError, ValueError, OverflowError):
        pass
    return default


def _to_int(value: object, default: int) -> int:
    try:
        if isinstance(value, (int, float, str, bytes, bytearray)):
            return int(float(value))
    except (TypeError, ValueError, OverflowError):
        pass
    return default


def _write_text_atomic(path: Path, text: str) -> None:
    write_text_atomic(path, text, require_nonempty=True)


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    write_json_atomic(path, payload, require_nonempty=True)


def _select_amber_tool(cfg: object | None, tool: str) -> tuple[list[str], str, str]:
    selection = select_ambertools_tool(cfg, tool)
    return selection.runner, selection.exe, selection.source


def _select_sander_command(cfg: object | None) -> tuple[list[str], str]:
    mpi_enabled = _to_bool(_cfg_get(cfg, "MMGBSA_MD_MPI_ENABLED", True), True)
    if mpi_enabled:
        try:
            selection = select_ambertools_tool(cfg, "sander.MPI", prefix_exe_style="absolute")
            launcher, launcher_source = select_mpi_launcher(
                selection.prefix, probe_common_prefixes=False
            )
            if launcher:
                ranks = resolve_md_mpi_ranks(cfg)
                runner = list(selection.runner)
                runner.extend([launcher, "-np", str(ranks), selection.exe])
                source = f"{selection.source};mpi:{launcher_source};ranks:{ranks}"
                return runner, source
        except FileNotFoundError:
            pass
    runner, sander_bin, source = _select_amber_tool(cfg, "sander")
    return list(runner) + [sander_bin], source


def _md_engine(cfg: object | None) -> str:
    engine = str(_cfg_get(cfg, "MMGBSA_MD_ENGINE", "sander") or "sander").strip().lower()
    if engine in {"amber", "sander", "sander.mpi"}:
        return "sander"
    if engine in {"openmm", "omm", "hip", "opencl"}:
        return "openmm"
    raise ValueError(f"Unsupported MMGBSA_MD_ENGINE={engine!r}; use sander or openmm")


def _publication_protocol(cfg: object | None) -> bool:
    token = str(_cfg_get(cfg, "MMGBSA_PROTOCOL", "") or "").strip().lower()
    return token in {"production", "publication", "publish"}


def _openmm_start_stage(cfg: object | None) -> str:
    raw = str(_cfg_get(cfg, "MMGBSA_OPENMM_START_STAGE", "prod") or "prod").strip().lower()
    aliases = {
        "all": "min1",
        "min": "min1",
        "minimize": "min1",
        "minimization": "min1",
        "production": "prod",
    }
    stage = aliases.get(raw, raw)
    valid = {"min1", "min2", "min_rescue1", "min_rescue2", "heat", "density", "equil", "prod"}
    if stage not in valid:
        raise ValueError(
            f"Unsupported MMGBSA_OPENMM_START_STAGE={raw!r}; "
            "use min1, min2, heat, density, equil, or prod"
        )
    return stage


def _engine_for_stage(engine: str, stage: str, openmm_start_stage: str) -> str:
    if engine != "openmm":
        return engine
    def rank(name: str) -> int:
        if name.startswith("heat"):
            return 4
        return {
        "min1": 0,
        "min2": 1,
        "min_rescue1": 2,
        "min_rescue2": 3,
        "density": 5,
        "equil": 6,
        "prod": 7,
        }.get(name, 7)

    return "openmm" if rank(stage) >= rank(openmm_start_stage) else "sander"


def _water_box_command(model: str, buffer_a: float) -> str:
    model_token = model.strip().lower()
    if model_token != "tip3p":
        raise ValueError(f"explicit solvent currently supports tip3p only, got {model}")
    return f"solvateBox COM TIP3PBOX {buffer_a:.3f}"


def _mol2_xyz(line: str, in_atoms: bool) -> tuple[float, float, float] | None:
    if not in_atoms:
        return None
    parts = line.strip().split()
    if len(parts) < 5:
        return None
    try:
        return float(parts[2]), float(parts[3]), float(parts[4])
    except Exception:
        return None


def _path_xyz(path: Path) -> list[tuple[float, float, float]]:
    try:
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return []
    coords: list[tuple[float, float, float]] = []
    in_mol2_atoms = False
    for line in lines:
        pdb_coord = _pdb_xyz(line)
        if pdb_coord is not None:
            coords.append(pdb_coord)
            continue
        stripped = line.strip()
        if stripped.startswith("@<TRIPOS>ATOM"):
            in_mol2_atoms = True
            continue
        if stripped.startswith("@<TRIPOS>") and in_mol2_atoms:
            in_mol2_atoms = False
            continue
        mol2_coord = _mol2_xyz(stripped, in_mol2_atoms)
        if mol2_coord is not None:
            coords.append(mol2_coord)
    return coords


def _coord_bounds(paths: Sequence[Path]) -> tuple[float, float, float]:
    coords: list[tuple[float, float, float]] = []
    for path in paths:
        coords.extend(_path_xyz(path))
    if not coords:
        return 0.0, 0.0, 0.0
    xs, ys, zs = zip(*coords)
    return max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)


def _estimate_salt_pairs(
    paths: Sequence[Path], salt_molar: float, water_buffer_a: float
) -> int:
    if salt_molar <= 0:
        return 0
    dx, dy, dz = _coord_bounds(paths)
    if min(dx, dy, dz) <= 0:
        return 0
    volume_a3 = (dx + 2.0 * water_buffer_a) * (dy + 2.0 * water_buffer_a) * (
        dz + 2.0 * water_buffer_a
    )
    return max(0, int(round(float(salt_molar) * 0.000602214076 * volume_a3)))


def _explicit_plan_dict(
    *,
    out_path: Path,
    leap_path: Path,
    dry_complex_prmtop: str | Path,
    dry_receptor_prmtop: str | Path,
    dry_ligand_prmtop: str | Path,
) -> dict[str, str]:
    return {
        "work_dir": str(out_path),
        "solvated_leap": str(leap_path),
        "solvated_prmtop": str(out_path / "complex_solvated.prmtop"),
        "solvated_inpcrd": str(out_path / "complex_solvated.inpcrd"),
        "dry_complex_prmtop": str(dry_complex_prmtop),
        "dry_receptor_prmtop": str(dry_receptor_prmtop),
        "dry_ligand_prmtop": str(dry_ligand_prmtop),
        "stripped_traj": str(out_path / "prod_stripped.nc"),
        "production_traj": str(out_path / "prod.nc"),
        "methods_json": str(out_path / "explicit_solvent_methods.json"),
    }


def _solvated_tleap_lines(
    *,
    out_dir: str | Path,
    receptor_path: Path,
    mol2_path: Path,
    frcmod_path: Path,
    water_model: str,
    water_buffer_a: float,
    neutralize: bool,
    salt_molar: float,
) -> tuple[list[str], int]:
    out_path = Path(out_dir)
    receptor_rel = os.path.relpath(receptor_path, out_path)
    mol2_rel = os.path.relpath(mol2_path, out_path)
    frcmod_rel = os.path.relpath(frcmod_path, out_path)
    lines = [
        "source leaprc.protein.ff14SB",
        "source leaprc.gaff2",
        "source leaprc.water.tip3p",
        f"loadamberparams {frcmod_rel}",
        f"LIG = loadmol2 {mol2_rel}",
        f"REC = loadpdb {receptor_rel}",
        "COM = combine { REC LIG }",
        _water_box_command(water_model, water_buffer_a),
    ]
    if neutralize:
        lines.extend(["addIonsRand COM Na+ 0", "addIonsRand COM Cl- 0"])
    salt_pairs = _estimate_salt_pairs(
        (receptor_path, mol2_path), salt_molar, water_buffer_a
    )
    if salt_pairs > 0:
        lines.append(f"addIonsRand COM Na+ {salt_pairs}")
        lines.append(f"addIonsRand COM Cl- {salt_pairs}")
    lines.extend(
        [
            "saveamberparm COM complex_solvated.prmtop complex_solvated.inpcrd",
            "quit",
            "",
        ]
    )
    return lines, salt_pairs


def write_solvated_tleap(
    *,
    receptor_pdb: str | Path,
    ligand_mol2: str | Path,
    ligand_frcmod: str | Path,
    dry_complex_prmtop: str | Path,
    dry_receptor_prmtop: str | Path,
    dry_ligand_prmtop: str | Path,
    out_dir: str | Path,
    water_model: str = "tip3p",
    water_buffer_a: float = 10.0,
    neutralize: bool = True,
    salt_molar: float = 0.150,
    force: bool = False,
) -> dict[str, Any]:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    receptor_path = Path(receptor_pdb)
    mol2_path = Path(ligand_mol2)
    frcmod_path = Path(ligand_frcmod)
    for required in (receptor_path, mol2_path, frcmod_path):
        if not required.exists():
            raise FileNotFoundError(f"explicit-solvent input not found: {required}")

    leap_path = out_path / "build_explicit.leap"
    lines, salt_pairs = _solvated_tleap_lines(
        out_dir=out_path,
        receptor_path=receptor_path,
        mol2_path=mol2_path,
        frcmod_path=frcmod_path,
        water_model=water_model,
        water_buffer_a=water_buffer_a,
        neutralize=neutralize,
        salt_molar=salt_molar,
    )
    if force or not leap_path.exists():
        _write_text_atomic(leap_path, "\n".join(lines))

    plan = _explicit_plan_dict(
        out_path=out_path,
        leap_path=leap_path,
        dry_complex_prmtop=dry_complex_prmtop,
        dry_receptor_prmtop=dry_receptor_prmtop,
        dry_ligand_prmtop=dry_ligand_prmtop,
    )
    return {
        **plan,
        "water_model": water_model,
        "water_buffer_a": water_buffer_a,
        "neutralize": neutralize,
        "salt_molar": salt_molar,
        "estimated_salt_pairs": salt_pairs,
        "skip_tleap": Path(plan["solvated_prmtop"]).exists()
        and Path(plan["solvated_prmtop"]).stat().st_size > 0
        and Path(plan["solvated_inpcrd"]).exists()
        and Path(plan["solvated_inpcrd"]).stat().st_size > 0
        and not force,
    }


def _sander_input(
    *,
    title: str,
    imin: int,
    nstlim: int,
    dt: float,
    temp0: float,
    ntb: int,
    ntp: int,
    cut: float,
    ntwx: int,
    irest: int = 0,
    ntx: int = 1,
    restraint_wt: float = 0.0,
    restraintmask: str = "@CA,C,N,O",
    seed: Optional[int] = None,
) -> str:
    lines = [
        title,
        " &cntrl",
        f"  imin={imin}, irest={irest}, ntx={ntx},",
        f"  ntb={ntb}, ntp={ntp}, cut={cut:.3f},",
        f"  ntpr=500, ntwx={ntwx}, ntwr=5000,",
    ]
    if imin == 1:
        lines.append(f"  maxcyc={nstlim}, ncyc={max(1, nstlim // 2)}, ntmin=1,")
    else:
        lines.extend(
            [
                f"  nstlim={nstlim}, dt={dt:.6f},",
                f"  tempi=10.0, temp0={temp0:.2f},",
                "  ntt=3, gamma_ln=2.0,",
                "  ntc=2, ntf=2,",
            ]
        )
    if seed is not None:
        lines.append(f"  ig={int(seed)},")
    if restraint_wt > 0:
        lines.extend(
            [
                "  ntr=1,",
                f"  restraint_wt={restraint_wt:.3f},",
                f"  restraintmask='{restraintmask}',",
            ]
        )
    lines.extend([" /", ""])
    return "\n".join(lines)


def _mdout_float(token: str) -> float | None:
    text = token.strip()
    if not text or "*" in text:
        return None
    try:
        return float(text)
    except (TypeError, ValueError, OverflowError):
        return None


def _parse_minimization_out(path: Path) -> dict[str, Any]:
    try:
        lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception as exc:
        return {"ok": False, "path": str(path), "reason": f"unreadable:{exc}"}

    nstep_rows: list[dict[str, Any]] = []
    components: dict[str, float | None] = {}
    nonfinite_marker_count = 0
    row_pattern = re.compile(
        r"^\s*(\d+)\s+([*+\-0-9.Ee]+)\s+([*+\-0-9.Ee]+)\s+"
        r"([*+\-0-9.Ee]+)\s+(\S+)\s+(\d+)"
    )
    component_pattern = re.compile(r"\b([A-Z0-9\- ]+?)\s*=\s*([*+\-0-9.Ee]+)")
    for line in lines:
        if "*" in line:
            nonfinite_marker_count += 1
        row_match = row_pattern.match(line)
        if row_match:
            nstep_rows.append(
                {
                    "step": int(row_match.group(1)),
                    "energy": _mdout_float(row_match.group(2)),
                    "rms": _mdout_float(row_match.group(3)),
                    "gmax": _mdout_float(row_match.group(4)),
                    "gmax_atom": row_match.group(5),
                    "gmax_atom_number": int(row_match.group(6)),
                }
            )
            continue
        for key, raw in component_pattern.findall(line):
            components[key.strip().replace(" ", "_").lower()] = _mdout_float(raw)
    return {
        "ok": bool(nstep_rows),
        "path": str(path),
        "last_row": nstep_rows[-1] if nstep_rows else {},
        "components": components,
        "nonfinite_marker_count": nonfinite_marker_count,
    }


def _preheat_qc_thresholds(cfg: object | None) -> dict[str, float]:
    return {
        "max_gmax": _to_float(_cfg_get(cfg, "MMGBSA_PREHEAT_QC_MAX_GMAX", 5000.0), 5000.0),
        "max_rms": _to_float(_cfg_get(cfg, "MMGBSA_PREHEAT_QC_MAX_RMS", 100.0), 100.0),
        "max_abs_energy": _to_float(
            _cfg_get(cfg, "MMGBSA_PREHEAT_QC_MAX_ABS_ENERGY", 1.0e9),
            1.0e9,
        ),
        "max_abs_vdw": _to_float(_cfg_get(cfg, "MMGBSA_PREHEAT_QC_MAX_ABS_VDW", 1.0e8), 1.0e8),
    }


def _preheat_qc_problems(
    *,
    gmax: object,
    rms: object,
    energy: object,
    vdw: object,
    thresholds: Mapping[str, float],
) -> list[str]:
    checks = [
        ("high_gmax", gmax, thresholds["max_gmax"], False),
        ("high_rms", rms, thresholds["max_rms"], False),
        ("high_abs_energy", energy, thresholds["max_abs_energy"], True),
        ("high_abs_vdw", vdw, thresholds["max_abs_vdw"], True),
    ]
    problems: list[str] = []
    for label, value, threshold, use_abs in checks:
        if value is None:
            problems.append(label)
            continue
        if not isinstance(value, (int, float, str)):
            problems.append(label)
            continue
        parsed = abs(float(value)) if use_abs else float(value)
        if parsed > threshold:
            problems.append(label)
    return problems


def _explicit_preheat_qc(
    *,
    work_dir: Path,
    cfg: object | None,
    step_name: str = "min2",
) -> dict[str, Any]:
    parsed = _parse_minimization_out(work_dir / f"{step_name}.out")
    if not parsed.get("ok"):
        return {**parsed, "ok": False, "qc_pass": False, "problems": ["missing_minimization_rows"]}

    last_row = parsed.get("last_row", {})
    components = parsed.get("components", {})
    gmax = last_row.get("gmax")
    rms = last_row.get("rms")
    energy = last_row.get("energy")
    vdw = components.get("vdwaals")
    thresholds = _preheat_qc_thresholds(cfg)
    problems = _preheat_qc_problems(
        gmax=gmax,
        rms=rms,
        energy=energy,
        vdw=vdw,
        thresholds=thresholds,
    )
    marker_count = int(parsed.get("nonfinite_marker_count", 0) or 0)
    return {
        **parsed,
        "qc_pass": not problems,
        "problems": problems,
        "thresholds": thresholds,
        "early_nonfinite_marker_count": marker_count,
    }


def _rst7_coordinates(path: Path) -> list[tuple[float, float, float]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    if len(lines) < 3:
        return []
    try:
        natom = int(lines[1].split()[0])
    except (IndexError, ValueError):
        return []
    values: list[float] = []
    for line in lines[2:]:
        for token in line.split():
            try:
                values.append(float(token))
            except ValueError:
                pass
        if len(values) >= natom * 3:
            break
    coords: list[tuple[float, float, float]] = []
    for offset in range(0, min(len(values), natom * 3), 3):
        if offset + 2 < len(values):
            coords.append((values[offset], values[offset + 1], values[offset + 2]))
    return coords


def _close_contact_qc(work_dir: Path, *, restart_name: str = "min2.rst7") -> dict[str, Any]:
    coords = _rst7_coordinates(work_dir / restart_name)
    cutoff_a = 0.65
    cell = cutoff_a
    grid: dict[tuple[int, int, int], list[tuple[int, tuple[float, float, float]]]] = {}
    severe_contacts = 0
    min_distance = float("inf")
    for idx, xyz in enumerate(coords):
        key = tuple(int(value // cell) for value in xyz)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for other_idx, other_xyz in grid.get((key[0] + dx, key[1] + dy, key[2] + dz), []):
                        if abs(idx - other_idx) <= 3:
                            continue
                        dist2 = sum((xyz[axis] - other_xyz[axis]) ** 2 for axis in range(3))
                        if dist2 < cutoff_a * cutoff_a:
                            severe_contacts += 1
                            min_distance = min(min_distance, dist2**0.5)
        grid.setdefault(key, []).append((idx, xyz))
    ok = bool(coords) and severe_contacts == 0
    return {
        "ok": ok,
        "restart": str(work_dir / restart_name),
        "atom_count": len(coords),
        "severe_contact_cutoff_a": cutoff_a,
        "severe_contact_count": severe_contacts,
        "min_distance_a": None if min_distance == float("inf") else min_distance,
    }


def _add_close_contact_qc(preheat_qc: dict[str, Any], work_dir: Path) -> dict[str, Any]:
    contact_qc = _close_contact_qc(work_dir)
    problems = list(preheat_qc.get("problems", []))
    if not contact_qc.get("ok"):
        problems.append("severe_close_contacts")
    return {
        **preheat_qc,
        "qc_pass": bool(preheat_qc.get("qc_pass")) and bool(contact_qc.get("ok")),
        "problems": problems,
        "close_contact_qc": contact_qc,
    }


def write_explicit_mdin_files(
    *,
    out_dir: str | Path,
    minimization_steps: int = 5000,
    heat_ps: float = 100.0,
    equil_ps: float = 1000.0,
    production_ps: float = 10000.0,
    frame_stride_ps: float = 20.0,
    dt_ps: float = 0.002,
    temperature_k: float = 300.0,
    seed: int = 1,
    force: bool = False,
) -> dict[str, str]:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    ntwx = max(1, int(round(frame_stride_ps / dt_ps)))
    min_steps = max(1, int(minimization_steps))
    heat_steps = max(1, int(round(max(0.0, float(heat_ps)) / dt_ps)))
    density_steps = max(1, int(round(max(0.0, float(equil_ps)) * 0.5 / dt_ps)))
    equil_steps = max(1, int(round(max(0.0, float(equil_ps)) * 0.5 / dt_ps)))
    prod_steps = max(1, int(round(production_ps / dt_ps)))
    inputs = {
        "min1": _sander_input(
            title="Explicit solvent restrained minimization",
            imin=1,
            nstlim=min_steps,
            dt=dt_ps,
            temp0=temperature_k,
            ntb=1,
            ntp=0,
            cut=10.0,
            ntwx=ntwx,
            restraint_wt=10.0,
        ),
        "min2": _sander_input(
            title="Explicit solvent unrestrained minimization",
            imin=1,
            nstlim=min_steps,
            dt=dt_ps,
            temp0=temperature_k,
            ntb=1,
            ntp=0,
            cut=10.0,
            ntwx=ntwx,
        ),
        "heat": _sander_input(
            title="Explicit solvent heating",
            imin=0,
            nstlim=heat_steps,
            dt=dt_ps,
            temp0=temperature_k,
            ntb=1,
            ntp=0,
            cut=10.0,
            ntwx=ntwx,
            irest=0,
            ntx=1,
            restraint_wt=5.0,
            seed=seed,
        ),
        "density": _sander_input(
            title="Explicit solvent density equilibration",
            imin=0,
            nstlim=density_steps,
            dt=dt_ps,
            temp0=temperature_k,
            ntb=2,
            ntp=1,
            cut=10.0,
            ntwx=ntwx,
            irest=1,
            ntx=5,
            restraint_wt=2.0,
            seed=seed,
        ),
        "equil": _sander_input(
            title="Explicit solvent unrestrained equilibration",
            imin=0,
            nstlim=equil_steps,
            dt=dt_ps,
            temp0=temperature_k,
            ntb=2,
            ntp=1,
            cut=10.0,
            ntwx=ntwx,
            irest=1,
            ntx=5,
            seed=seed,
        ),
        "prod": _sander_input(
            title="Explicit solvent production",
            imin=0,
            nstlim=prod_steps,
            dt=dt_ps,
            temp0=temperature_k,
            ntb=2,
            ntp=1,
            cut=10.0,
            ntwx=ntwx,
            irest=1,
            ntx=5,
            seed=seed,
        ),
    }
    written: dict[str, str] = {}
    for name, text in inputs.items():
        path = out_path / f"{name}.in"
        _write_text_atomic(path, text)
        written[name] = str(path)
    return written


def _write_rescue_mdin_files(
    *,
    out_dir: Path,
    rescue_steps: int,
    dt_ps: float,
    temperature_k: float,
    ntwx: int,
) -> dict[str, str]:
    steps = max(1, int(rescue_steps))
    inputs = {
        "min_rescue1": _sander_input(
            title="Explicit solvent high-energy rescue restrained minimization",
            imin=1,
            nstlim=steps,
            dt=dt_ps,
            temp0=temperature_k,
            ntb=1,
            ntp=0,
            cut=10.0,
            ntwx=ntwx,
            restraint_wt=25.0,
        ),
        "min_rescue2": _sander_input(
            title="Explicit solvent high-energy rescue unrestrained minimization",
            imin=1,
            nstlim=steps,
            dt=dt_ps,
            temp0=temperature_k,
            ntb=1,
            ntp=0,
            cut=10.0,
            ntwx=ntwx,
        ),
    }
    written: dict[str, str] = {}
    for name, text in inputs.items():
        path = out_dir / f"{name}.in"
        _write_text_atomic(path, text)
        written[name] = str(path)
    return written


def run_tleap(plan: Mapping[str, Any], cfg: object | None, force: bool = False) -> dict[str, Any]:
    prmtop = Path(str(plan["solvated_prmtop"]))
    inpcrd = Path(str(plan["solvated_inpcrd"]))
    if prmtop.exists() and prmtop.stat().st_size > 0 and inpcrd.exists() and inpcrd.stat().st_size > 0 and not force:
        return {"ok": True, "skipped": True, "returncode": 0}
    runner, tleap_bin, source = _select_amber_tool(cfg, "tleap")
    work_dir = Path(str(plan["work_dir"]))
    leap_path = Path(str(plan["solvated_leap"]))
    cmd = list(runner) + [tleap_bin, "-f", os.path.relpath(leap_path, work_dir)]
    log_path = work_dir / "tleap_explicit.log"
    with log_path.open("w", encoding="utf-8") as handle:
        proc = subprocess.run(cmd, cwd=str(work_dir), stdout=handle, stderr=subprocess.STDOUT)
    ok = proc.returncode == 0 and prmtop.exists() and prmtop.stat().st_size > 0
    return {"ok": ok, "skipped": False, "returncode": proc.returncode, "cmd": cmd, "log": str(log_path), "source": source}


def _run_sander_step(
    *,
    work_dir: Path,
    cfg: object | None,
    prmtop: Path,
    mdin: str,
    coord_in: str,
    restart_out: str,
    traj_out: str = "",
    ref_coord: str = "",
) -> dict[str, Any]:
    sander_cmd, source = _select_sander_command(cfg)
    cmd = list(sander_cmd) + [
        "-O",
        "-i",
        f"{mdin}.in",
        "-p",
        os.path.relpath(prmtop, work_dir),
        "-c",
        coord_in,
        "-o",
        f"{mdin}.out",
        "-r",
        restart_out,
    ]
    if traj_out:
        cmd.extend(["-x", traj_out])
    if ref_coord:
        cmd.extend(["-ref", ref_coord])
    cmd_path = work_dir / f"{mdin}.cmd.txt"
    cmd_path.write_text(" ".join(cmd), encoding="utf-8")
    log_path = work_dir / f"{mdin}.launch.log"
    env = os.environ.copy()
    apply_blas_single_thread_env_defaults(env)
    with log_path.open("w", encoding="utf-8") as handle:
        proc = subprocess.run(
            cmd,
            cwd=str(work_dir),
            stdout=handle,
            stderr=subprocess.STDOUT,
            env=env,
        )
    return {"step": mdin, "ok": proc.returncode == 0, "returncode": proc.returncode, "cmd": str(cmd_path), "log": str(log_path), "source": source}


def _run_md_step(
    *,
    engine: str,
    work_dir: Path,
    cfg: object | None,
    prmtop: Path,
    mdin: str,
    coord_in: str,
    restart_out: str,
    traj_out: str,
    ref_coord: str,
    nsteps: int,
    dt_ps: float,
    temperature_k: float,
    ntwx: int,
    pressure: bool,
    restraint_wt: float,
    restraint_scope: str,
    seed: int,
    minimize: bool,
) -> dict[str, Any]:
    if engine == "openmm":
        return run_openmm_step(
            work_dir=work_dir,
            cfg=cfg,
            prmtop=prmtop,
            coord_in=coord_in,
            restart_out=restart_out,
            step_name=mdin,
            nsteps=nsteps,
            dt_ps=dt_ps,
            temperature_k=temperature_k,
            ntwx=ntwx,
            traj_out=traj_out,
            pressure=pressure,
            restraint_wt=restraint_wt,
            restraint_scope=restraint_scope,
            seed=seed,
            minimize=minimize,
        )
    return _run_sander_step(
        work_dir=work_dir,
        cfg=cfg,
        prmtop=prmtop,
        mdin=mdin,
        coord_in=coord_in,
        restart_out=restart_out,
        traj_out=traj_out,
        ref_coord=ref_coord,
    )


def _stage_resume_ready(
    work_dir: Path,
    *,
    mdin: str,
    restart_out: str,
    traj_out: str,
    force: bool,
) -> bool:
    if force:
        return False
    restart_path = work_dir / restart_out
    if not restart_path.exists() or restart_path.stat().st_size == 0:
        return False
    if traj_out:
        traj_path = work_dir / traj_out
        if not traj_path.exists() or traj_path.stat().st_size == 0:
            return False
    if mdin in {"min1", "min2", "min_rescue1", "min_rescue2"}:
        out_path = work_dir / f"{mdin}.out"
        return out_path.exists() and out_path.stat().st_size > 0
    return True


def _normalize_step_spec(
    spec: Sequence[Any],
    *,
    default_temperature_k: float,
    default_dt_ps: float,
) -> tuple[str, str, str, str, str, int, bool, float, bool, str, float, float]:
    mdin, coord_in, restart_out, traj_out, ref_coord, nsteps, pressure, restraint_wt, minimize = spec[:9]
    restraint_scope = str(spec[9]) if len(spec) > 9 else "backbone"
    temperature_k = float(spec[10]) if len(spec) > 10 else default_temperature_k
    dt_ps = float(spec[11]) if len(spec) > 11 else default_dt_ps
    return (
        str(mdin),
        str(coord_in),
        str(restart_out),
        str(traj_out),
        str(ref_coord),
        int(nsteps),
        bool(pressure),
        float(restraint_wt),
        bool(minimize),
        restraint_scope,
        temperature_k,
        dt_ps,
    )


def _stage_resume_result(
    work_dir: Path,
    *,
    mdin: str,
    restart_out: str,
    traj_out: str,
    engine: str,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "ok": True,
        "skipped": True,
        "resume": True,
        "step": mdin,
        "engine": engine,
        "restart": str(work_dir / restart_out),
    }
    if traj_out:
        payload["traj"] = str(work_dir / traj_out)
    return payload


def _run_explicit_step_sequence(
    *,
    steps: Sequence[tuple[str, str, str, str, str, int, bool, float, bool]],
    engine: str,
    openmm_start_stage: str,
    work_dir: Path,
    cfg: object | None,
    prmtop: Path,
    dt_ps: float,
    temperature_k: float,
    ntwx: int,
    seed: int,
    force: bool,
    allow_min2_failure: bool = False,
) -> tuple[list[dict[str, Any]], str]:
    results: list[dict[str, Any]] = []
    for spec in steps:
        (
            mdin,
            coord_in,
            restart_out,
            traj_out,
            ref_coord,
            nsteps,
            pressure,
            restraint_wt,
            minimize,
            restraint_scope,
            step_temperature_k,
            step_dt_ps,
        ) = _normalize_step_spec(
            spec,
            default_temperature_k=temperature_k,
            default_dt_ps=dt_ps,
        )
        step_engine = _engine_for_stage(engine, mdin, openmm_start_stage)
        if _stage_resume_ready(
            work_dir,
            mdin=mdin,
            restart_out=restart_out,
            traj_out=traj_out,
            force=force,
        ):
            result = _stage_resume_result(
                work_dir,
                mdin=mdin,
                restart_out=restart_out,
                traj_out=traj_out,
                engine=step_engine,
            )
        else:
            result = _run_md_step(
                engine=step_engine,
                work_dir=work_dir,
                cfg=cfg,
                prmtop=prmtop,
                mdin=mdin,
                coord_in=coord_in,
                restart_out=restart_out,
                traj_out=traj_out,
                ref_coord=ref_coord,
                nsteps=nsteps,
                dt_ps=step_dt_ps,
                temperature_k=step_temperature_k,
                ntwx=ntwx,
                pressure=pressure,
                restraint_wt=restraint_wt,
                restraint_scope=restraint_scope,
                seed=seed,
                minimize=minimize,
            )
        results.append(result)
        if not result["ok"] and not (allow_min2_failure and mdin == "min2"):
            return results, mdin
    return results, ""


def _run_explicit_rescue_steps(
    *,
    engine: str,
    openmm_start_stage: str,
    work_dir: Path,
    cfg: object | None,
    prmtop: Path,
    rescue_steps: int,
    dt_ps: float,
    temperature_k: float,
    ntwx: int,
    seed: int,
    force: bool,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rescue_runs: list[dict[str, Any]] = []
    rescue_result: dict[str, Any] = {"enabled": True, "attempted": True}
    rescue_specs = [
        ("min_rescue1", "min1.rst7", "min_rescue1.rst7", "", "min1.rst7"),
        ("min_rescue2", "min_rescue1.rst7", "min2.rst7", "", ""),
    ]
    step_specs = [
        (mdin, coord_in, restart_out, traj_out, ref_coord, rescue_steps, False, 10.0 if mdin == "min_rescue1" else 0.0, True)
        for mdin, coord_in, restart_out, traj_out, ref_coord in rescue_specs
    ]
    rescue_runs, _failed = _run_explicit_step_sequence(
        steps=step_specs,
        engine=engine,
        openmm_start_stage=openmm_start_stage,
        work_dir=work_dir,
        cfg=cfg,
        prmtop=prmtop,
        dt_ps=dt_ps,
        temperature_k=temperature_k,
        ntwx=ntwx,
        seed=seed,
        force=force,
    )
    rescue_result["steps"] = rescue_runs
    return rescue_runs, rescue_result


def _maybe_rescue_preheat(
    *,
    results: list[dict[str, Any]],
    preheat_qc: dict[str, Any],
    rescue_enabled: bool,
    work_dir: Path,
    cfg: object | None,
    prmtop: Path,
    engine: str,
    openmm_start_stage: str,
    rescue_steps: int,
    dt_ps: float,
    frame_stride_ps: float,
    temperature_k: float,
    ntwx: int,
    seed: int,
    force: bool,
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    rescue_result: dict[str, Any] = {"enabled": rescue_enabled, "attempted": False}
    if results[-1].get("ok") and preheat_qc.get("qc_pass"):
        return results, preheat_qc, rescue_result
    if not rescue_enabled:
        return results, preheat_qc, rescue_result
    _write_rescue_mdin_files(
        out_dir=work_dir,
        rescue_steps=rescue_steps,
        dt_ps=dt_ps,
        temperature_k=_to_float(_cfg_get(cfg, "MMGBSA_MD_TEMP0", 300.0), 300.0),
        ntwx=max(1, int(round(frame_stride_ps / dt_ps))),
    )
    rescue_runs, rescue_result = _run_explicit_rescue_steps(
        engine=engine,
        openmm_start_stage=openmm_start_stage,
        work_dir=work_dir,
        cfg=cfg,
        prmtop=prmtop,
        rescue_steps=rescue_steps,
        dt_ps=dt_ps,
        temperature_k=temperature_k,
        ntwx=ntwx,
        seed=seed,
        force=force,
    )
    results.extend(rescue_runs)
    preheat_qc = _add_close_contact_qc(
        _explicit_preheat_qc(work_dir=work_dir, cfg=cfg, step_name="min2"),
        work_dir,
    )
    return results, preheat_qc, rescue_result


def _heat_rescue_settings(
    cfg: object | None,
    *,
    dt_ps: float,
    heat_ps: float,
    temperature_k: float,
) -> dict[str, float]:
    rescue_dt_ps = _to_float(
        _cfg_get(cfg, "MMGBSA_MD_HEAT_RESCUE_DT_PS", min(dt_ps * 0.5, 0.0005)),
        min(dt_ps * 0.5, 0.0005),
    )
    rescue_heat_ps = _to_float(
        _cfg_get(cfg, "MMGBSA_MD_HEAT_RESCUE_PS", max(20.0, heat_ps * 5.0)),
        max(20.0, heat_ps * 5.0),
    )
    restraint_wt = _to_float(
        _cfg_get(cfg, "MMGBSA_MD_HEAT_RESCUE_RESTRAINT_WT", 10.0),
        10.0,
    )
    return {
        "dt_ps": max(0.0001, min(float(dt_ps), float(rescue_dt_ps))),
        "heat_ps": max(float(heat_ps), float(rescue_heat_ps)),
        "temperature_k": float(temperature_k),
        "restraint_wt": max(0.0, float(restraint_wt)),
    }


def _maybe_rescue_heat_failure(
    *,
    results: list[dict[str, Any]],
    failed_step: str,
    rescue_result: dict[str, Any],
    preheat_qc: dict[str, Any],
    rescue_enabled: bool,
    work_dir: Path,
    cfg: object | None,
    prmtop: Path,
    engine: str,
    openmm_start_stage: str,
    rescue_steps: int,
    dt_ps: float,
    heat_ps: float,
    frame_stride_ps: float,
    temperature_k: float,
    ntwx: int,
    density_steps: int,
    equil_steps: int,
    prod_steps: int,
    traj_path: Path,
    seed: int,
    force: bool,
) -> tuple[list[dict[str, Any]], str, dict[str, Any], dict[str, Any]]:
    if not failed_step.startswith("heat") or not rescue_enabled:
        return results, failed_step, preheat_qc, rescue_result

    rescue_result["enabled"] = True
    rescue_result["heat_attempted"] = True
    rescue_result["heat_reason"] = "heat_failed"
    _write_rescue_mdin_files(
        out_dir=work_dir,
        rescue_steps=rescue_steps,
        dt_ps=dt_ps,
        temperature_k=temperature_k,
        ntwx=ntwx,
    )
    rescue_runs, heat_rescue_result = _run_explicit_rescue_steps(
        engine=engine,
        openmm_start_stage=openmm_start_stage,
        work_dir=work_dir,
        cfg=cfg,
        prmtop=prmtop,
        rescue_steps=rescue_steps,
        dt_ps=dt_ps,
        temperature_k=temperature_k,
        ntwx=ntwx,
        seed=seed,
        force=True,
    )
    results.extend(rescue_runs)
    rescue_result["heat_relaxation"] = heat_rescue_result
    if any(not step.get("ok") for step in rescue_runs):
        rescue_result["heat_ok"] = False
        return results, "heat_rescue_relaxation", preheat_qc, rescue_result

    preheat_qc = _add_close_contact_qc(
        _explicit_preheat_qc(work_dir=work_dir, cfg=cfg, step_name="min2"),
        work_dir,
    )
    if not preheat_qc.get("qc_pass"):
        rescue_result["heat_ok"] = False
        return results, "heat_rescue_qc", preheat_qc, rescue_result

    rescue = _heat_rescue_settings(
        cfg,
        dt_ps=dt_ps,
        heat_ps=heat_ps,
        temperature_k=temperature_k,
    )
    rescue_heat_steps = max(1, int(round(rescue["heat_ps"] / rescue["dt_ps"])))
    rescue_ntwx = max(1, int(round(frame_stride_ps / rescue["dt_ps"])))
    heat_retry_steps = [
        (
            "heat",
            "min2.rst7",
            "heat.rst7",
            "",
            "min2.rst7",
            rescue_heat_steps,
            False,
            rescue["restraint_wt"],
            False,
        )
    ]
    heat_retry, heat_failed = _run_explicit_step_sequence(
        steps=heat_retry_steps,
        engine=engine,
        openmm_start_stage=openmm_start_stage,
        work_dir=work_dir,
        cfg=cfg,
        prmtop=prmtop,
        dt_ps=rescue["dt_ps"],
        temperature_k=rescue["temperature_k"],
        ntwx=rescue_ntwx,
        seed=seed,
        force=True,
    )
    results.extend(heat_retry)
    rescue_result["heat_retry"] = {
        "dt_ps": rescue["dt_ps"],
        "heat_ps": rescue["heat_ps"],
        "restraint_wt": rescue["restraint_wt"],
        "steps": heat_retry,
    }
    if heat_failed:
        rescue_result["heat_ok"] = False
        return results, heat_failed, preheat_qc, rescue_result

    remaining_steps = [
        ("density", "heat.rst7", "density.rst7", "", "heat.rst7", density_steps, True, 2.0, False),
        ("equil", "density.rst7", "equil.rst7", "", "", equil_steps, True, 0.0, False),
        ("prod", "equil.rst7", "prod.rst7", traj_path.name, "", prod_steps, True, 0.0, False),
    ]
    remaining_results, remaining_failed = _run_explicit_step_sequence(
        steps=remaining_steps,
        engine=engine,
        openmm_start_stage=openmm_start_stage,
        work_dir=work_dir,
        cfg=cfg,
        prmtop=prmtop,
        dt_ps=dt_ps,
        temperature_k=temperature_k,
        ntwx=ntwx,
        seed=seed,
        force=force,
    )
    results.extend(remaining_results)
    rescue_result["heat_ok"] = not bool(remaining_failed)
    return results, remaining_failed, preheat_qc, rescue_result


def _explicit_md_preflight(
    *, prmtop: Path, inpcrd: Path, traj_path: Path, force: bool, run: bool
) -> dict[str, Any] | None:
    if traj_path.exists() and traj_path.stat().st_size > 0 and not force:
        return {"ok": True, "skipped": True, "traj_path": str(traj_path)}
    if not run:
        return {"ok": False, "skipped": True, "planned": True, "traj_path": str(traj_path)}
    missing_topology = not prmtop.exists() or not inpcrd.exists()
    if missing_topology:
        raise FileNotFoundError("explicit solvated topology is missing; run tleap first")
    return None


def run_explicit_md(plan: Mapping[str, Any], cfg: object | None, *, seed: int, force: bool = False, run: bool = True) -> dict[str, Any]:
    work_dir = Path(str(plan["work_dir"]))
    prmtop = Path(str(plan["solvated_prmtop"]))
    inpcrd = Path(str(plan["solvated_inpcrd"]))
    traj_path = Path(str(plan["production_traj"]))
    dt_ps = _to_float(_cfg_get(cfg, "MMGBSA_MD_DT_PS", 0.002), 0.002)
    min_steps = int(_to_float(_cfg_get(cfg, "MMGBSA_MD_MIN_STEPS", 10000), 10000))
    rescue_enabled = _to_bool(_cfg_get(cfg, "MMGBSA_MD_RESCUE_ENABLED", True), True)
    rescue_steps = _to_int(
        _cfg_get(cfg, "MMGBSA_MD_RESCUE_STEPS", max(50000, min_steps * 5)),
        max(50000, min_steps * 5),
    )
    heat_ps = _to_float(_cfg_get(cfg, "MMGBSA_MD_HEAT_PS", 100.0), 100.0)
    equil_ps = _to_float(_cfg_get(cfg, "MMGBSA_MD_EQUIL_PS", 1000.0), 1000.0)
    prod_ps = _to_float(_cfg_get(cfg, "MMGBSA_MD_PROD_PS", 10000.0), 10000.0)
    frame_stride_ps = _to_float(_cfg_get(cfg, "MMGBSA_MD_FRAME_STRIDE_PS", 20.0), 20.0)
    temperature_k = _to_float(_cfg_get(cfg, "MMGBSA_MD_TEMP0", 300.0), 300.0)
    engine = _md_engine(cfg)
    openmm_start_stage = _openmm_start_stage(cfg)
    publication = _publication_protocol(cfg)
    ntwx = max(1, int(round(frame_stride_ps / dt_ps)))
    heat_steps = max(1, int(round(max(0.0, heat_ps) / dt_ps)))
    density_steps = max(1, int(round(max(0.0, equil_ps) * 0.5 / dt_ps)))
    equil_steps = max(1, int(round(max(0.0, equil_ps) * 0.5 / dt_ps)))
    prod_steps = max(1, int(round(prod_ps / dt_ps)))
    write_explicit_mdin_files(
        out_dir=work_dir,
        minimization_steps=min_steps,
        heat_ps=heat_ps,
        equil_ps=equil_ps,
        production_ps=prod_ps,
        frame_stride_ps=frame_stride_ps,
        dt_ps=dt_ps,
        seed=seed,
        temperature_k=temperature_k,
        force=force,
    )
    preflight = _explicit_md_preflight(
        prmtop=prmtop, inpcrd=inpcrd, traj_path=traj_path, force=force, run=run
    )
    if preflight is not None:
        return preflight

    if publication and engine == "openmm":
        heat_chunk = max(1, heat_steps // 3)
        preheat_steps = [
            (
                "min1",
                os.path.relpath(inpcrd, work_dir),
                "min1.rst7",
                "",
                os.path.relpath(inpcrd, work_dir),
                min_steps,
                False,
                25.0,
                True,
                "solute_heavy",
            ),
            ("min2", "min1.rst7", "min2.rst7", "", "min1.rst7", min_steps, False, 10.0, True, "solute_heavy"),
        ]
        production_steps = [
            ("heat_100", "min2.rst7", "heat_100.rst7", "", "min2.rst7", heat_chunk, False, 10.0, False, "solute_heavy", 100.0),
            ("heat_200", "heat_100.rst7", "heat_200.rst7", "", "heat_100.rst7", heat_chunk, False, 7.5, False, "solute_heavy", 200.0),
            (
                "heat",
                "heat_200.rst7",
                "heat.rst7",
                "",
                "heat_200.rst7",
                max(1, heat_steps - 2 * heat_chunk),
                False,
                5.0,
                False,
                "solute_heavy",
                temperature_k,
            ),
            ("density", "heat.rst7", "density.rst7", "", "heat.rst7", density_steps, True, 2.0, False, "solute_heavy"),
            ("equil", "density.rst7", "equil.rst7", "", "density.rst7", equil_steps, True, 0.5, False, "solute_heavy"),
            ("prod", "equil.rst7", "prod.rst7", traj_path.name, "", prod_steps, True, 0.0, False, "backbone"),
        ]
    else:
        preheat_steps = [
            (
                "min1",
                os.path.relpath(inpcrd, work_dir),
                "min1.rst7",
                "",
                os.path.relpath(inpcrd, work_dir),
                min_steps,
                False,
                10.0,
                True,
            ),
            ("min2", "min1.rst7", "min2.rst7", "", "", min_steps, False, 0.0, True),
        ]
        production_steps = [
            ("heat", "min2.rst7", "heat.rst7", "", "min2.rst7", heat_steps, False, 5.0, False),
            ("density", "heat.rst7", "density.rst7", "", "heat.rst7", density_steps, True, 2.0, False),
            ("equil", "density.rst7", "equil.rst7", "", "", equil_steps, True, 0.0, False),
            ("prod", "equil.rst7", "prod.rst7", traj_path.name, "", prod_steps, True, 0.0, False),
        ]
    results, failed_step = _run_explicit_step_sequence(
        steps=preheat_steps,
        engine=engine,
        openmm_start_stage=openmm_start_stage,
        work_dir=work_dir,
        cfg=cfg,
        prmtop=prmtop,
        dt_ps=dt_ps,
        temperature_k=temperature_k,
        ntwx=ntwx,
        seed=seed,
        force=force,
        allow_min2_failure=True,
    )
    if failed_step:
        return {"ok": False, "skipped": False, "step_failed": failed_step, "steps": results, "traj_path": str(traj_path)}
    preheat_qc = _add_close_contact_qc(
        _explicit_preheat_qc(work_dir=work_dir, cfg=cfg, step_name="min2"),
        work_dir,
    )
    results, preheat_qc, rescue_result = _maybe_rescue_preheat(
        results=results,
        preheat_qc=preheat_qc,
        rescue_enabled=rescue_enabled,
        work_dir=work_dir,
        cfg=cfg,
        prmtop=prmtop,
        engine=engine,
        openmm_start_stage=openmm_start_stage,
        rescue_steps=rescue_steps,
        dt_ps=dt_ps,
        frame_stride_ps=frame_stride_ps,
        temperature_k=temperature_k,
        ntwx=ntwx,
        seed=seed,
        force=force,
    )
    if not preheat_qc.get("qc_pass"):
        return {
            "ok": False,
            "skipped": False,
            "step_failed": "preheat_qc",
            "steps": results,
            "preheat_qc": preheat_qc,
            "rescue": rescue_result,
            "traj_path": str(traj_path),
        }

    production_results, failed_step = _run_explicit_step_sequence(
        steps=production_steps,
        engine=engine,
        openmm_start_stage=openmm_start_stage,
        work_dir=work_dir,
        cfg=cfg,
        prmtop=prmtop,
        dt_ps=dt_ps,
        temperature_k=temperature_k,
        ntwx=ntwx,
        seed=seed,
        force=force,
    )
    results.extend(production_results)
    if failed_step:
        results, failed_step, preheat_qc, rescue_result = _maybe_rescue_heat_failure(
            results=results,
            failed_step=failed_step,
            rescue_result=rescue_result,
            preheat_qc=preheat_qc,
            rescue_enabled=rescue_enabled,
            work_dir=work_dir,
            cfg=cfg,
            prmtop=prmtop,
            engine=engine,
            openmm_start_stage=openmm_start_stage,
            rescue_steps=rescue_steps,
            dt_ps=dt_ps,
            heat_ps=heat_ps,
            frame_stride_ps=frame_stride_ps,
            temperature_k=temperature_k,
            ntwx=ntwx,
            density_steps=density_steps,
            equil_steps=equil_steps,
            prod_steps=prod_steps,
            traj_path=traj_path,
            seed=seed,
            force=force,
        )
    if failed_step:
        return {
            "ok": False,
            "skipped": False,
            "step_failed": failed_step,
            "steps": results,
            "preheat_qc": preheat_qc,
            "rescue": rescue_result,
            "traj_path": str(traj_path),
        }
    return {
        "ok": traj_path.exists() and traj_path.stat().st_size > 0,
        "skipped": False,
        "steps": results,
        "preheat_qc": preheat_qc,
        "rescue": rescue_result,
        "traj_path": str(traj_path),
        "engine": engine,
        "openmm_start_stage": openmm_start_stage,
    }


def write_strip_cpptraj(plan: Mapping[str, Any], *, force: bool = False) -> str:
    work_dir = Path(str(plan["work_dir"]))
    input_path = work_dir / "strip_explicit_for_mmpbsa.in"
    text = "\n".join(
        [
            f"parm {plan['solvated_prmtop']}",
            f"trajin {plan['production_traj']}",
            "autoimage",
            "strip :WAT,HOH,Na+,Cl-,K+",
            f"trajout {plan['stripped_traj']} netcdf",
            "run",
            "",
        ]
    )
    if force or not input_path.exists():
        _write_text_atomic(input_path, text)
    return str(input_path)


def run_cpptraj_strip(plan: Mapping[str, Any], cfg: object | None, *, force: bool = False, run: bool = True) -> dict[str, Any]:
    stripped = Path(str(plan["stripped_traj"]))
    input_path = Path(write_strip_cpptraj(plan, force=force))
    if stripped.exists() and stripped.stat().st_size > 0 and not force:
        return {"ok": True, "skipped": True, "stripped_traj": str(stripped), "input": str(input_path)}
    if not run:
        return {"ok": False, "skipped": True, "planned": True, "stripped_traj": str(stripped), "input": str(input_path)}
    runner, cpptraj_bin, source = _select_amber_tool(cfg, "cpptraj")
    work_dir = Path(str(plan["work_dir"]))
    cmd = list(runner) + [cpptraj_bin, "-i", os.path.relpath(input_path, work_dir)]
    log_path = work_dir / "strip_explicit_for_mmpbsa.log"
    with log_path.open("w", encoding="utf-8") as handle:
        proc = subprocess.run(cmd, cwd=str(work_dir), stdout=handle, stderr=subprocess.STDOUT)
    ok = proc.returncode == 0 and stripped.exists() and stripped.stat().st_size > 0
    return {"ok": ok, "skipped": False, "returncode": proc.returncode, "cmd": cmd, "log": str(log_path), "source": source, "stripped_traj": str(stripped), "input": str(input_path)}


def _series_svg(values: Sequence[tuple[float, float]], title: str) -> str:
    width = 760
    height = 280
    pad = 36
    if not values:
        return "<svg xmlns='http://www.w3.org/2000/svg' width='760' height='280'></svg>\n"
    xs, ys = zip(*values)
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    if x_min == x_max:
        x_max = x_min + 1.0
    if y_min == y_max:
        y_min -= 1.0
        y_max += 1.0
    points = []
    for x_val, y_val in values:
        x = pad + (x_val - x_min) * (width - 2 * pad) / (x_max - x_min)
        y = height - pad - (y_val - y_min) * (height - 2 * pad) / (y_max - y_min)
        points.append(f"{x:.2f},{y:.2f}")
    return (
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
        f"viewBox='0 0 {width} {height}'>\n"
        "<rect width='100%' height='100%' fill='white'/>\n"
        f"<text x='{pad}' y='22' font-family='sans-serif' font-size='14'>{title}</text>\n"
        f"<line x1='{pad}' y1='{height-pad}' x2='{width-pad}' y2='{height-pad}' stroke='#333'/>\n"
        f"<line x1='{pad}' y1='{pad}' x2='{pad}' y2='{height-pad}' stroke='#333'/>\n"
        f"<polyline fill='none' stroke='#1f77b4' stroke-width='2' points='{' '.join(points)}'/>\n"
        f"<text x='{pad}' y='{height-8}' font-family='sans-serif' font-size='11'>frame/residue</text>\n"
        "</svg>\n"
    )


def _read_numeric_xy(path: Path) -> list[tuple[float, float]]:
    values: list[tuple[float, float]] = []
    if not path.exists():
        return values
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "@")):
            continue
        parts = stripped.split()
        if len(parts) < 2:
            continue
        try:
            values.append((float(parts[0]), float(parts[-1])))
        except (TypeError, ValueError, OverflowError):
            continue
    return values


def _write_qc_plot(dat_path: Path, title: str) -> str:
    svg_path = dat_path.with_suffix(".svg")
    svg_path.write_text(_series_svg(_read_numeric_xy(dat_path), title), encoding="utf-8")
    return str(svg_path)


def write_explicit_qc_cpptraj(plan: Mapping[str, Any], cfg: object | None, *, force: bool = False) -> str:
    work_dir = Path(str(plan["work_dir"]))
    qc_dir = work_dir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    input_path = qc_dir / "cpptraj_md_qc.in"
    ligand_mask = str(_cfg_get(cfg, "MMGBSA_QC_LIGAND_MASK", ":MOL") or ":MOL")
    complex_rmsd = qc_dir / "complex_rmsd.dat"
    ligand_rmsd = qc_dir / "ligand_rmsd.dat"
    protein_rmsf = qc_dir / "protein_rmsf_byres.dat"
    text = "\n".join(
        [
            f"parm {plan['dry_complex_prmtop']}",
            f"trajin {plan['stripped_traj']}",
            "autoimage",
            f"rms first out {complex_rmsd} @CA,C,N,O",
            f"rms first out {ligand_rmsd} {ligand_mask}",
            f"atomicfluct out {protein_rmsf} byres :1-999999&!{ligand_mask}",
            "run",
            "",
        ]
    )
    if force or not input_path.exists():
        _write_text_atomic(input_path, text)
    return str(input_path)


def run_explicit_qc_cpptraj(plan: Mapping[str, Any], cfg: object | None, *, force: bool = False) -> dict[str, Any]:
    stripped = Path(str(plan["stripped_traj"]))
    input_path = Path(write_explicit_qc_cpptraj(plan, cfg, force=force))
    if not stripped.exists() or stripped.stat().st_size == 0:
        return {"ok": False, "skipped": True, "reason": "missing_stripped_traj"}
    runner, cpptraj_bin, source = _select_amber_tool(cfg, "cpptraj")
    work_dir = Path(str(plan["work_dir"]))
    log_path = work_dir / "qc" / "cpptraj_md_qc.log"
    cmd = list(runner) + [cpptraj_bin, "-i", os.path.relpath(input_path, work_dir)]
    with log_path.open("w", encoding="utf-8") as handle:
        proc = subprocess.run(cmd, cwd=str(work_dir), stdout=handle, stderr=subprocess.STDOUT)
    qc_dir = work_dir / "qc"
    complex_rmsd = qc_dir / "complex_rmsd.dat"
    ligand_rmsd = qc_dir / "ligand_rmsd.dat"
    protein_rmsf = qc_dir / "protein_rmsf_byres.dat"
    ok = proc.returncode == 0
    return {
        "ok": ok,
        "skipped": False,
        "returncode": proc.returncode,
        "cmd": cmd,
        "log": str(log_path),
        "source": source,
        "input": str(input_path),
        "complex_rmsd": str(complex_rmsd),
        "ligand_rmsd": str(ligand_rmsd),
        "protein_rmsf": str(protein_rmsf),
        "complex_rmsd_svg": _write_qc_plot(complex_rmsd, "complex RMSD") if ok else "",
        "ligand_rmsd_svg": _write_qc_plot(ligand_rmsd, "ligand RMSD") if ok else "",
        "protein_rmsf_svg": _write_qc_plot(protein_rmsf, "protein RMSF by residue") if ok else "",
    }


def write_explicit_methods(plan: Mapping[str, Any], cfg: object | None, extra: Mapping[str, Any] | None = None) -> str:
    payload = {
        "schema_version": 1,
        "workflow": "explicit_solvent_mmgbsa",
        "md_engine": _cfg_get(cfg, "MMGBSA_MD_ENGINE", "sander"),
        "openmm_platform": _cfg_get(cfg, "MMGBSA_OPENMM_PLATFORM", "auto"),
        "openmm_device_index": _cfg_get(cfg, "MMGBSA_OPENMM_DEVICE_INDEX", ""),
        "openmm_start_stage": _cfg_get(cfg, "MMGBSA_OPENMM_START_STAGE", "prod"),
        "water_model": _cfg_get(cfg, "MMGBSA_TLEAP_WATER_MODEL", "tip3p"),
        "salt_molar": _cfg_get(cfg, "MMGBSA_EXPLICIT_SALT_MOLAR", 0.150),
        "water_buffer_a": _cfg_get(cfg, "MMGBSA_EXPLICIT_WATER_BUFFER_A", 10.0),
        "md_prod_ps": _cfg_get(cfg, "MMGBSA_MD_PROD_PS", 10000.0),
        "md_heat_ps": _cfg_get(cfg, "MMGBSA_MD_HEAT_PS", 100.0),
        "md_equil_ps": _cfg_get(cfg, "MMGBSA_MD_EQUIL_PS", 1000.0),
        "md_min_steps": _cfg_get(cfg, "MMGBSA_MD_MIN_STEPS", 10000),
        "md_rescue_enabled": _cfg_get(cfg, "MMGBSA_MD_RESCUE_ENABLED", True),
        "md_rescue_steps": _cfg_get(cfg, "MMGBSA_MD_RESCUE_STEPS", ""),
        "md_frame_stride_ps": _cfg_get(cfg, "MMGBSA_MD_FRAME_STRIDE_PS", 20.0),
        "preheat_qc": {
            "max_gmax": _cfg_get(cfg, "MMGBSA_PREHEAT_QC_MAX_GMAX", 5000.0),
            "max_rms": _cfg_get(cfg, "MMGBSA_PREHEAT_QC_MAX_RMS", 100.0),
            "max_abs_energy": _cfg_get(cfg, "MMGBSA_PREHEAT_QC_MAX_ABS_ENERGY", 1.0e9),
            "max_abs_vdw": _cfg_get(cfg, "MMGBSA_PREHEAT_QC_MAX_ABS_VDW", 1.0e8),
        },
        "strip_mask": ":WAT,HOH,Na+,Cl-,K+",
        "plan": dict(plan),
    }
    if extra:
        payload["extra"] = dict(extra)
    path = Path(str(plan["methods_json"]))
    _write_json_atomic(path, payload)
    return str(path)


def _heat_like_failure(md_result: Mapping[str, Any]) -> bool:
    failed = str(md_result.get("step_failed", "") or "")
    return failed.startswith("heat")


def build_and_run_explicit_workflow(
    *,
    receptor_pdb: str | Path,
    ligand_mol2: str | Path,
    ligand_frcmod: str | Path,
    dry_complex_prmtop: str | Path,
    dry_receptor_prmtop: str | Path,
    dry_ligand_prmtop: str | Path,
    out_dir: str | Path,
    cfg: object | None,
    seed: int,
    force: bool = False,
    run: bool = True,
) -> dict[str, Any]:
    plan = write_solvated_tleap(
        receptor_pdb=receptor_pdb,
        ligand_mol2=ligand_mol2,
        ligand_frcmod=ligand_frcmod,
        dry_complex_prmtop=dry_complex_prmtop,
        dry_receptor_prmtop=dry_receptor_prmtop,
        dry_ligand_prmtop=dry_ligand_prmtop,
        out_dir=out_dir,
        water_model=str(_cfg_get(cfg, "MMGBSA_TLEAP_WATER_MODEL", "tip3p") or "tip3p"),
        water_buffer_a=_to_float(_cfg_get(cfg, "MMGBSA_EXPLICIT_WATER_BUFFER_A", 10.0), 10.0),
        neutralize=_to_bool(_cfg_get(cfg, "MMGBSA_EXPLICIT_NEUTRALIZE", True), True),
        salt_molar=_to_float(_cfg_get(cfg, "MMGBSA_EXPLICIT_SALT_MOLAR", 0.150), 0.150),
        force=force,
    )
    tleap_result = run_tleap(plan, cfg, force=force) if run else {"ok": False, "planned": True}
    md_result = {"ok": False, "skipped": True}
    strip_result = {"ok": False, "skipped": True}
    qc_result = {"ok": False, "skipped": True}
    if run and tleap_result.get("ok"):
        md_result = run_explicit_md(plan, cfg, seed=seed, force=force, run=True)
        if (
            not md_result.get("ok")
            and _publication_protocol(cfg)
            and _heat_like_failure(md_result)
            and _to_bool(_cfg_get(cfg, "MMGBSA_MD_RESCUE_ENABLED", True), True)
        ):
            water_buffer_a = float(plan.get("water_buffer_a", _cfg_get(cfg, "MMGBSA_EXPLICIT_WATER_BUFFER_A", 10.0)) or 10.0)
            rebuild_plan = write_solvated_tleap(
                receptor_pdb=receptor_pdb,
                ligand_mol2=ligand_mol2,
                ligand_frcmod=ligand_frcmod,
                dry_complex_prmtop=dry_complex_prmtop,
                dry_receptor_prmtop=dry_receptor_prmtop,
                dry_ligand_prmtop=dry_ligand_prmtop,
                out_dir=out_dir,
                water_model=str(_cfg_get(cfg, "MMGBSA_TLEAP_WATER_MODEL", "tip3p") or "tip3p"),
                water_buffer_a=water_buffer_a + 2.0,
                neutralize=_to_bool(_cfg_get(cfg, "MMGBSA_EXPLICIT_NEUTRALIZE", True), True),
                salt_molar=_to_float(_cfg_get(cfg, "MMGBSA_EXPLICIT_SALT_MOLAR", 0.150), 0.150),
                force=True,
            )
            rebuild_tleap = run_tleap(rebuild_plan, cfg, force=True)
            retry_md = {"ok": False, "skipped": True}
            if rebuild_tleap.get("ok"):
                retry_md = run_explicit_md(rebuild_plan, cfg, seed=seed + 7919, force=True, run=True)
            md_result = {
                **retry_md,
                "initial_failure": md_result,
                "solvent_rebuild": {
                    "attempted": True,
                    "water_buffer_a": water_buffer_a + 2.0,
                    "tleap": rebuild_tleap,
                },
                "structurally_suspect": not bool(retry_md.get("ok")),
            }
            tleap_result = rebuild_tleap
            plan = rebuild_plan
        if md_result.get("ok"):
            strip_result = run_cpptraj_strip(plan, cfg, force=force, run=True)
            if strip_result.get("ok"):
                qc_result = run_explicit_qc_cpptraj(plan, cfg, force=force)
    methods_json = write_explicit_methods(
        plan,
        cfg,
        {
            "tleap": tleap_result,
            "md": md_result,
            "strip": strip_result,
            "qc": qc_result,
        },
    )
    return {
        "ok": bool(strip_result.get("ok")),
        "plan": plan,
        "tleap": tleap_result,
        "md": md_result,
        "strip": strip_result,
        "qc": qc_result,
        "traj_path": strip_result.get("stripped_traj", plan["stripped_traj"]),
        "methods_json": methods_json,
    }


__all__ = [
    "build_and_run_explicit_workflow",
    "run_cpptraj_strip",
    "run_explicit_qc_cpptraj",
    "run_explicit_md",
    "write_explicit_qc_cpptraj",
    "write_explicit_mdin_files",
    "write_explicit_methods",
    "write_solvated_tleap",
    "write_strip_cpptraj",
]
