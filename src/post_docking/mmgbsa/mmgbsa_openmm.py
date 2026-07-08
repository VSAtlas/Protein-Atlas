from __future__ import annotations

import math
import os
import re
import subprocess
from pathlib import Path
from typing import Any, Mapping


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


def _to_float(value: object, default: float) -> float:
    try:
        if isinstance(value, (int, float, str, bytes, bytearray)):
            return float(value)
    except (TypeError, ValueError):
        pass
    return default


def _to_int(value: object, default: int) -> int:
    try:
        if isinstance(value, (int, float, str, bytes, bytearray)):
            return int(float(value))
    except (TypeError, ValueError):
        pass
    return default


def _available_platforms() -> list[str]:
    import openmm as mm

    return [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]


def _resolve_platform_name(cfg: object | None) -> str:
    requested = str(_cfg_get(cfg, "MMGBSA_OPENMM_PLATFORM", "auto") or "auto").strip()
    platforms = _available_platforms()
    if requested and requested.lower() != "auto":
        for name in platforms:
            if name.lower() == requested.lower():
                return name
        raise ValueError(f"OpenMM platform {requested!r} is not available; found {platforms}")
    for preferred in ("HIP", "OpenCL", "CUDA", "CPU"):
        if preferred in platforms:
            return preferred
    return platforms[0]


def _discover_rocm_gpu_indices() -> str:
    tokens = _visible_device_tokens()
    if tokens:
        return ",".join(tokens)
    for env_name in ("HIP_VISIBLE_DEVICES", "ROCR_VISIBLE_DEVICES"):
        value = os.environ.get(env_name, "").strip()
        if value and value.lower() != "all":
            return value
    rocm_smi = Path("/opt/rocm/bin/rocm-smi")
    if not rocm_smi.exists():
        return "0"
    try:
        proc = subprocess.run(
            [str(rocm_smi), "--showproductname"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return "0"
    indices = sorted({int(match.group(1)) for match in re.finditer(r"GPU\[(\d+)\]", proc.stdout)})
    return ",".join(str(idx) for idx in indices) if indices else "0"


def _split_device_tokens(raw: str) -> list[str]:
    return [part.strip() for part in raw.replace(";", ",").split(",") if part.strip() and part.strip().lower() != "all"]


def _visible_device_tokens() -> list[str]:
    for env_name in ("GPU_DEVICE_ORDINAL", "ROCR_VISIBLE_DEVICES", "HIP_VISIBLE_DEVICES", "CUDA_VISIBLE_DEVICES"):
        value = os.environ.get(env_name, "").strip()
        tokens = _split_device_tokens(value)
        if tokens:
            return tokens
    return []


def _resolve_device_index(cfg: object | None, platform_name: str) -> str:
    raw = str(_cfg_get(cfg, "MMGBSA_OPENMM_DEVICE_INDEX", "") or "").strip()
    if raw.lower() == "all":
        return _discover_rocm_gpu_indices()
    if raw:
        return raw
    if platform_name in {"HIP", "CUDA", "OpenCL"}:
        visible = _visible_device_tokens()
        if visible:
            # Visibility masks such as ROCR_VISIBLE_DEVICES/HIP_VISIBLE_DEVICES
            # expose a dense local device namespace to HIP/OpenMM. If a worker
            # masks itself to physical GPU "2", OpenMM must usually select
            # logical DeviceIndex "0" inside that process.
            return "0"
        return str(_cfg_get(cfg, "MMGBSA_OPENMM_DEFAULT_DEVICE_INDEX", "0") or "0")
    return ""


def _platform_properties(cfg: object | None, platform_name: str) -> dict[str, str]:
    props: dict[str, str] = {}
    precision = str(_cfg_get(cfg, "MMGBSA_OPENMM_PRECISION", "mixed") or "mixed").strip()
    if precision and platform_name in {"HIP", "CUDA", "OpenCL"}:
        props["Precision"] = precision
    device_index = _resolve_device_index(cfg, platform_name)
    if device_index and platform_name in {"HIP", "CUDA", "OpenCL"}:
        props["DeviceIndex"] = device_index
    if platform_name == "CPU":
        threads = str(_cfg_get(cfg, "MMGBSA_OPENMM_CPU_THREADS", "1") or "1").strip()
        if threads:
            props["Threads"] = threads
    return props


class AmberNetCDFReporter:
    def __init__(self, file: str | Path, report_interval: int, atom_count: int) -> None:
        self._path = Path(file)
        self._interval = int(report_interval)
        self._atom_count = int(atom_count)
        self._frame = 0
        self._nc: Any | None = None

    def describeNextReport(self, simulation: Any) -> dict[str, Any]:  # noqa: N802
        steps = self._interval - simulation.currentStep % self._interval
        return {
            "steps": steps,
            "periodic": True,
            "include": ["positions"],
        }

    def report(self, simulation: Any, state: Any) -> None:
        import numpy as np
        from openmm import unit
        from scipy.io import netcdf_file

        self.describeNextReport
        nc = self._nc
        if nc is None:
            self._path.parent.mkdir(parents=True, exist_ok=True)
            nc = netcdf_file(str(self._path), "w", version=2)
            self._nc = nc
            setattr(nc, "Conventions", "AMBER")
            setattr(nc, "ConventionVersion", "1.0")
            setattr(nc, "program", "OpenMM")
            setattr(nc, "programVersion", "8")
            nc.createDimension("frame", None)
            nc.createDimension("atom", self._atom_count)
            nc.createDimension("spatial", 3)
            nc.createDimension("cell_spatial", 3)
            nc.createDimension("cell_angular", 3)
            nc.createDimension("label", 5)
            nc.createVariable("spatial", "c", ("spatial",))
            nc.variables["spatial"][:] = np.array([b"x", b"y", b"z"], dtype="S1")
            nc.createVariable("cell_spatial", "c", ("cell_spatial",))
            nc.variables["cell_spatial"][:] = np.array([b"a", b"b", b"c"], dtype="S1")
            nc.createVariable("cell_angular", "c", ("cell_angular", "label"))
            nc.variables["cell_angular"][:] = np.array(
                [[b"a", b"l", b"p", b"h", b"a"], [b"b", b"e", b"t", b"a", b" "], [b"g", b"a", b"m", b"m", b"a"]],
                dtype="S1",
            )
            nc.createVariable("coordinates", "f", ("frame", "atom", "spatial"))
            setattr(nc.variables["coordinates"], "units", "angstrom")
            nc.createVariable("time", "f", ("frame",))
            setattr(nc.variables["time"], "units", "picosecond")
            nc.createVariable("cell_lengths", "f", ("frame", "cell_spatial"))
            setattr(nc.variables["cell_lengths"], "units", "angstrom")
            nc.createVariable("cell_angles", "f", ("frame", "cell_angular"))
            setattr(nc.variables["cell_angles"], "units", "degree")

        coords = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        nc.variables["coordinates"][self._frame, :, :] = coords
        nc.variables["time"][self._frame] = state.getTime().value_in_unit(unit.picosecond)
        vectors = state.getPeriodicBoxVectors()
        lengths, angles = _box_lengths_angles(vectors)
        nc.variables["cell_lengths"][self._frame, :] = lengths
        nc.variables["cell_angles"][self._frame, :] = angles
        self._frame += 1
        nc.flush()

    def __del__(self) -> None:
        nc = getattr(self, "_nc", None)
        if nc is not None:
            nc.close()


def _box_lengths_angles(vectors: Any) -> tuple[list[float], list[float]]:
    from openmm import Vec3, unit

    vals = [vec.value_in_unit(unit.angstrom) for vec in vectors]
    a, b, c = vals

    def norm(v: Vec3) -> float:
        return math.sqrt(float(v.x * v.x + v.y * v.y + v.z * v.z))

    def angle(u: Vec3, v: Vec3) -> float:
        dot = float(u.x * v.x + u.y * v.y + u.z * v.z)
        denom = max(norm(u) * norm(v), 1.0e-12)
        return math.degrees(math.acos(max(-1.0, min(1.0, dot / denom))))

    return [norm(a), norm(b), norm(c)], [angle(b, c), angle(a, c), angle(a, b)]


def _write_amber_rst7(
    path: Path,
    positions: Any,
    box_vectors: Any,
    *,
    title: str,
    time_ps: float,
) -> None:
    from openmm import unit

    coords = positions.value_in_unit(unit.angstrom)
    lengths, angles = _box_lengths_angles(box_vectors)
    values = [coord for pos in coords for coord in (float(pos.x), float(pos.y), float(pos.z))]
    lines = [title, f"{len(coords):6d}{float(time_ps):15.7f}"]
    for i in range(0, len(values), 6):
        lines.append("".join(f"{value:12.7f}" for value in values[i : i + 6]))
    lines.append("".join(f"{value:12.7f}" for value in [*lengths, *angles]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _is_solvent_or_ion(residue_name: str) -> bool:
    return residue_name.upper() in {"WAT", "HOH", "Na+", "Cl-", "K+", "NA", "CL"}


def _is_hydrogen(atom: Any) -> bool:
    element = getattr(atom, "element", None)
    symbol = str(getattr(element, "symbol", "") or "").upper()
    return symbol == "H" or str(getattr(atom, "name", "") or "").upper().startswith("H")


def _restrain_atom(atom: Any, scope: str) -> bool:
    token = str(scope or "backbone").strip().lower()
    residue_name = str(getattr(getattr(atom, "residue", None), "name", "") or "")
    if token == "solute_heavy":
        return not _is_solvent_or_ion(residue_name) and not _is_hydrogen(atom)
    if token == "heavy":
        return not _is_hydrogen(atom)
    return str(getattr(atom, "name", "") or "") in {"CA", "C", "N", "O"}


def _add_positional_restraints(
    system: Any, positions: Any, topology: Any, weight: float, scope: str = "backbone"
) -> None:
    if weight <= 0:
        return
    import openmm as mm
    from openmm import unit

    force = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter("k", float(weight) * 418.4 * unit.kilojoule_per_mole / unit.nanometer**2)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    pos_nm = positions.value_in_unit(unit.nanometer)
    for atom in topology.atoms():
        if _restrain_atom(atom, scope):
            ref = pos_nm[int(atom.index)]
            force.addParticle(int(atom.index), [float(ref.x), float(ref.y), float(ref.z)])
    if force.getNumParticles() > 0:
        system.addForce(force)


def _write_minimization_out(path: Path, *, step: int, energy_kcal: float, gmax: float, rms: float) -> None:
    text = "\n".join(
        [
            "OpenMM minimization summary",
            "",
            "   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER",
            f"{step:7d} {energy_kcal:14.4E} {rms:12.4E} {gmax:12.4E}     O         1",
            "",
            " VDWAALS =        0.0000  EEL     =        0.0000  RESTRAINT  =        0.0000",
            f" EAMBER  = {energy_kcal:14.4f}",
            "",
        ]
    )
    path.write_text(text, encoding="utf-8")


def _force_stats(state: Any) -> tuple[float, float]:
    from openmm import unit

    forces = state.getForces(asNumpy=True).value_in_unit(
        unit.kilocalorie_per_mole / unit.angstrom
    )
    mags = [math.sqrt(float(f[0] * f[0] + f[1] * f[1] + f[2] * f[2])) for f in forces]
    if not mags:
        return 0.0, 0.0
    return max(mags), math.sqrt(sum(v * v for v in mags) / len(mags))


def _build_simulation(
    *,
    cfg: object | None,
    prmtop: Path,
    coord_path: Path,
    dt_ps: float,
    temperature_k: float,
    pressure: bool,
    restraint_wt: float,
    restraint_scope: str,
    seed: int,
) -> tuple[Any, Any, str, dict[str, str]]:
    import openmm as mm
    from openmm import unit
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile, HBonds, PME, Simulation

    platform_name = _resolve_platform_name(cfg)
    platform = mm.Platform.getPlatformByName(platform_name)
    properties = _platform_properties(cfg, platform_name)
    prmtop_obj = AmberPrmtopFile(str(prmtop))
    inpcrd_obj = AmberInpcrdFile(str(coord_path))
    system = prmtop_obj.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=10.0 * unit.angstrom,
        constraints=HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
    )
    _add_positional_restraints(
        system,
        inpcrd_obj.positions,
        prmtop_obj.topology,
        restraint_wt,
        scope=restraint_scope,
    )
    if pressure:
        system.addForce(mm.MonteCarloBarostat(1.0 * unit.atmosphere, temperature_k * unit.kelvin, 25))
    integrator = mm.LangevinMiddleIntegrator(
        temperature_k * unit.kelvin,
        2.0 / unit.picosecond,
        dt_ps * unit.picoseconds,
    )
    integrator.setRandomNumberSeed(int(seed))
    simulation = Simulation(prmtop_obj.topology, system, integrator, platform, properties)
    simulation.context.setPositions(inpcrd_obj.positions)
    if inpcrd_obj.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd_obj.boxVectors)
    simulation.context.setVelocitiesToTemperature(temperature_k * unit.kelvin, int(seed))
    return simulation, prmtop_obj, platform_name, properties


def _write_command(
    path: Path,
    *,
    platform_name: str,
    properties: Mapping[str, str],
    step_name: str,
    nsteps: int,
) -> None:
    path.write_text(
        " ".join(
            [
                "openmm",
                f"--platform={platform_name}",
                f"--properties={dict(properties)}",
                f"--step={step_name}",
                f"--nsteps={nsteps}",
            ]
        ),
        encoding="utf-8",
    )


def _run_minimization(simulation: Any, out_path: Path, nsteps: int) -> None:
    from openmm import unit

    simulation.minimizeEnergy(maxIterations=max(1, int(nsteps)))
    state = simulation.context.getState(getEnergy=True, getForces=True, getPositions=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    gmax, rms = _force_stats(state)
    _write_minimization_out(
        out_path,
        step=max(1, int(nsteps)),
        energy_kcal=float(energy),
        gmax=float(gmax),
        rms=float(rms),
    )


def _run_dynamics(
    *,
    simulation: Any,
    out_handle: Any,
    work_dir: Path,
    traj_out: str,
    ntwx: int,
    nsteps: int,
    atom_count: int,
    temperature_k: float,
    heat_ramp: bool = False,
) -> None:
    from openmm import unit
    from openmm.app import StateDataReporter

    simulation.reporters.append(
        StateDataReporter(
            out_handle,
            max(1, min(500, int(nsteps))),
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            speed=True,
            separator=" ",
        )
    )
    if traj_out:
        simulation.reporters.append(
            AmberNetCDFReporter(work_dir / traj_out, max(1, int(ntwx)), atom_count)
        )
    total_steps = max(1, int(nsteps))
    if not heat_ramp or total_steps < 10:
        simulation.step(total_steps)
        return
    integrator = simulation.integrator
    chunks = min(20, total_steps)
    done = 0
    for chunk_index in range(chunks):
        chunk_steps = total_steps // chunks
        if chunk_index < total_steps % chunks:
            chunk_steps += 1
        done += chunk_steps
        fraction = done / total_steps
        target_k = max(25.0, float(temperature_k) * fraction)
        integrator.setTemperature(target_k * unit.kelvin)
        simulation.step(chunk_steps)


def _write_final_restart(simulation: Any, path: Path, step_name: str) -> None:
    from openmm import unit

    final_state = simulation.context.getState(getPositions=True)
    _write_amber_rst7(
        path,
        final_state.getPositions(),
        final_state.getPeriodicBoxVectors(),
        title=f"OpenMM {step_name} restart",
        time_ps=final_state.getTime().value_in_unit(unit.picosecond),
    )


def _result_payload(
    *,
    step_name: str,
    ok: bool,
    cmd_path: Path,
    log_path: Path,
    platform_name: str,
    properties: Mapping[str, str],
    error: str = "",
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "step": step_name,
        "ok": ok,
        "returncode": 0 if ok else 1,
        "cmd": str(cmd_path),
        "log": str(log_path),
        "source": f"openmm:{platform_name}",
        "platform": platform_name,
        "platform_properties": dict(properties),
    }
    if error:
        payload["error"] = error
    return payload


def run_openmm_step(
    *,
    work_dir: Path,
    cfg: object | None,
    prmtop: Path,
    coord_in: str,
    restart_out: str,
    step_name: str,
    nsteps: int,
    dt_ps: float,
    temperature_k: float,
    ntwx: int,
    traj_out: str = "",
    pressure: bool = False,
    restraint_wt: float = 0.0,
    restraint_scope: str = "backbone",
    seed: int = 1,
    minimize: bool = False,
) -> dict[str, Any]:
    simulation, prmtop_obj, platform_name, properties = _build_simulation(
        cfg=cfg,
        prmtop=prmtop,
        coord_path=work_dir / coord_in,
        dt_ps=dt_ps,
        temperature_k=temperature_k,
        pressure=pressure,
        restraint_wt=restraint_wt,
        restraint_scope=restraint_scope,
        seed=seed,
    )

    cmd_path = work_dir / f"{step_name}.cmd.txt"
    _write_command(
        cmd_path,
        platform_name=platform_name,
        properties=properties,
        step_name=step_name,
        nsteps=nsteps,
    )
    out_path = work_dir / f"{step_name}.out"
    log_path = work_dir / f"{step_name}.launch.log"
    try:
        with out_path.open("w", encoding="utf-8") as out_handle:
            if minimize:
                _run_minimization(simulation, out_path, nsteps)
            else:
                if step_name.startswith("heat") and coord_in == "min2.rst7":
                    from openmm import unit

                    simulation.context.setVelocitiesToTemperature(25.0 * unit.kelvin, int(seed))
                _run_dynamics(
                    simulation=simulation,
                    out_handle=out_handle,
                    work_dir=work_dir,
                    traj_out=traj_out,
                    ntwx=ntwx,
                    nsteps=nsteps,
                    atom_count=prmtop_obj.topology.getNumAtoms(),
                    temperature_k=temperature_k,
                    heat_ramp=step_name.startswith("heat"),
                )
            _write_final_restart(simulation, work_dir / restart_out, step_name)
        log_path.write_text(
            f"platform={platform_name}\nproperties={properties}\nstep={step_name}\n",
            encoding="utf-8",
        )
        ok = (work_dir / restart_out).exists() and (work_dir / restart_out).stat().st_size > 0
        if traj_out:
            ok = ok and (work_dir / traj_out).exists() and (work_dir / traj_out).stat().st_size > 0
        return _result_payload(
            step_name=step_name,
            ok=ok,
            cmd_path=cmd_path,
            log_path=log_path,
            platform_name=platform_name,
            properties=properties,
        )
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}"
        log_path.write_text(f"OpenMM {step_name} failed: {error}\n", encoding="utf-8")
        return _result_payload(
            step_name=step_name,
            ok=False,
            cmd_path=cmd_path,
            log_path=log_path,
            platform_name=platform_name,
            properties=properties,
            error=error,
        )
