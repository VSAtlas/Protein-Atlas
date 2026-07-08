"""AutoDock4Zn planning and optional execution helpers.

AutoDock4Zn is a zinc-specific AD4 map workflow. These helpers only route
eligible Zn sites to that workflow; they do not invent generic pseudo-atoms for
other metals.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Mapping, Sequence


AUTODOCK4ZN_TOOLS: tuple[str, ...] = (
    "zinc_pseudo.py",
    "prepare_gpf4zn.py",
    "AD4Zn.dat",
    "pythonsh",
    "autogrid4",
    "vina",
)


def build_autodock4zn_plan(
    *,
    receptor_pdbqt: Path | None,
    metal_rows: Sequence[Mapping[str, object]],
    work_dir: Path,
    ligand_pdbqt: Path | None = None,
    center: tuple[float, float, float] | None = None,
    box_size: float | None = None,
    exhaustiveness: int = 32,
) -> dict[str, object]:
    """Return a runnable AutoDock4Zn plan for eligible zinc sites."""

    candidates = _ad4zn_candidates(metal_rows)
    non_zinc_excluded = _non_zinc_parameterized_count(metal_rows)
    tools = _tool_paths(AUTODOCK4ZN_TOOLS)
    missing = [name for name, path in tools.items() if not path]
    status = _plan_status(
        candidate_count=len(candidates),
        missing_tools=missing,
        receptor_ready=_path_ready(receptor_pdbqt),
        ligand_ready=ligand_pdbqt is None or ligand_pdbqt.exists(),
    )
    work_dir.mkdir(parents=True, exist_ok=True)
    outputs = _ad4zn_outputs(work_dir)
    commands = _autodock4zn_commands(
        receptor_pdbqt=receptor_pdbqt,
        receptor_tz=outputs["receptor_tz_pdbqt"],
        gpf=outputs["gpf"],
        glg=outputs["glg"],
        maps_prefix=outputs["maps_prefix"],
        ligand_pdbqt=ligand_pdbqt,
        center=center,
        box_size=box_size,
        exhaustiveness=exhaustiveness,
    )
    return {
        "status": status,
        "candidate_count": len(candidates),
        "non_zinc_excluded_count": non_zinc_excluded,
        "tool_paths": tools,
        "missing_tools": missing,
        "receptor_pdbqt": str(receptor_pdbqt) if receptor_pdbqt else "",
        "ligand_pdbqt": str(ligand_pdbqt) if ligand_pdbqt else "",
        "work_dir": str(work_dir),
        "outputs": {name: str(path) for name, path in outputs.items()},
        "commands": commands,
        "metals": [_metal_plan_row(row) for row in candidates],
        "note": (
            "AutoDock4Zn applies only to eligible Zn sites and requires AD4Zn maps. "
            "It is not a generic metal formal-charge or PDBQT parameter solution."
        ),
    }


def write_autodock4zn_plan(path: Path, plan: Mapping[str, object]) -> None:
    """Write an AutoDock4Zn plan sidecar."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(plan, indent=2, sort_keys=True), encoding="utf-8")


def summarize_autodock4zn_plan(plan: Mapping[str, object]) -> dict[str, object]:
    """Flatten an AutoDock4Zn plan for benchmark CSV output."""

    missing = plan.get("missing_tools", [])
    return {
        "autodock4zn_status": str(plan.get("status", "")),
        "autodock4zn_candidate_count": _safe_int(plan.get("candidate_count", 0)),
        "autodock4zn_non_zinc_excluded_count": _safe_int(
            plan.get("non_zinc_excluded_count", 0)
        ),
        "autodock4zn_missing_tools": ",".join(missing)
        if isinstance(missing, list)
        else str(missing or ""),
    }


def run_autodock4zn_plan(plan: Mapping[str, object]) -> dict[str, object]:
    """Execute the receptor/map-generation commands from a ready plan."""

    if plan.get("status") != "ready":
        return {
            "status": "skipped",
            "reason": f"plan_not_ready:{plan.get('status')}",
            "steps": [],
        }
    steps: list[dict[str, object]] = []
    _stage_ad4zn_dat(plan, steps)
    for command in _command_steps(plan):
        result = subprocess.run(
            command,
            cwd=str(plan.get("work_dir") or "."),
            text=True,
            capture_output=True,
            check=False,
        )
        steps.append(
            {
                "command": command,
                "returncode": result.returncode,
                "stdout_tail": result.stdout[-4000:],
                "stderr_tail": result.stderr[-4000:],
            }
        )
        if result.returncode != 0:
            return {"status": "failed", "steps": steps}
    return {"status": "ok", "steps": steps}


def _stage_ad4zn_dat(
    plan: Mapping[str, object],
    steps: list[dict[str, object]],
) -> None:
    tool_paths = plan.get("tool_paths", {})
    if not isinstance(tool_paths, Mapping):
        return
    source = str(tool_paths.get("AD4Zn.dat", ""))
    if not source:
        return
    work_dir = Path(str(plan.get("work_dir") or "."))
    target = work_dir / "AD4Zn.dat"
    if Path(source).resolve() != target.resolve():
        shutil.copyfile(source, target)
    steps.append(
        {
            "command": ["stage_ad4zn_dat", source, str(target)],
            "returncode": 0,
            "stdout_tail": "",
            "stderr_tail": "",
        }
    )


def _tool_paths(names: Sequence[str]) -> dict[str, str]:
    return {name: _resolve_autodock4zn_tool(name) for name in names}


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _staged_autodock4zn_dir() -> Path:
    return Path(
        str(os.environ.get("AUTODOCK4ZN_DIR") or _repo_root() / "tools" / "autodock4zn")
    ).expanduser()


def _resolve_autodock4zn_tool(name: str) -> str:
    if name == "pythonsh":
        return shutil.which("pythonsh") or _mgltools_pythonsh()
    if name == "AD4Zn.dat":
        configured = os.environ.get("AD4ZN_DAT")
        if configured and Path(configured).expanduser().is_file():
            return str(Path(configured).expanduser())
        staged = _staged_autodock4zn_dir() / "AD4Zn.dat"
        return str(staged) if staged.is_file() else ""
    staged = _staged_autodock4zn_dir() / name
    if staged.is_file():
        return str(staged)
    return shutil.which(name) or ""


def _mgltools_pythonsh() -> str:
    candidates = (
        Path("/stor/system/opt/mgltools_x86_64Linux2_1.5.7/bin/pythonsh"),
        Path.home() / "mgltools_x86_64Linux2_1.5.7" / "bin" / "pythonsh",
    )
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    return ""


def _ad4zn_candidates(
    metal_rows: Sequence[Mapping[str, object]],
) -> list[Mapping[str, object]]:
    return [
        row for row in metal_rows if row.get("metal_treatment") == "autodock4zn_candidate"
    ]


def _non_zinc_parameterized_count(metal_rows: Sequence[Mapping[str, object]]) -> int:
    return sum(
        1
        for row in metal_rows
        if row.get("requires_parameterization")
        and row.get("metal_treatment") != "autodock4zn_candidate"
    )


def _path_ready(path: Path | None) -> bool:
    return path is not None and path.exists()


def _ad4zn_outputs(work_dir: Path) -> dict[str, Path]:
    return {
        "receptor_tz_pdbqt": work_dir / "receptor_ad4zn_tz.pdbqt",
        "gpf": work_dir / "receptor_ad4zn_tz.gpf",
        "glg": work_dir / "receptor_ad4zn_tz.glg",
        "maps_prefix": work_dir / "receptor_ad4zn_tz",
    }


def _plan_status(
    *,
    candidate_count: int,
    missing_tools: Sequence[str],
    receptor_ready: bool,
    ligand_ready: bool,
) -> str:
    if candidate_count == 0:
        return "not_applicable"
    if not receptor_ready:
        return "missing_receptor"
    if not ligand_ready:
        return "missing_ligand"
    if missing_tools:
        return "missing_tools"
    return "ready"


def _autodock4zn_commands(
    *,
    receptor_pdbqt: Path | None,
    receptor_tz: Path,
    gpf: Path,
    glg: Path,
    maps_prefix: Path,
    ligand_pdbqt: Path | None,
    center: tuple[float, float, float] | None,
    box_size: float | None,
    exhaustiveness: int,
) -> list[dict[str, object]]:
    commands: list[dict[str, object]] = []
    tools = _tool_paths(AUTODOCK4ZN_TOOLS)
    if receptor_pdbqt is not None:
        commands.append(
            {
                "name": "add_tz_pseudoatoms",
                "argv": [
                    sys.executable,
                    tools.get("zinc_pseudo.py", "zinc_pseudo.py"),
                    "-r",
                    str(receptor_pdbqt),
                    "-o",
                    str(receptor_tz),
                ],
            }
        )
    gpf_args = [
        tools.get("pythonsh", "pythonsh"),
        tools.get("prepare_gpf4zn.py", "prepare_gpf4zn.py"),
    ]
    if ligand_pdbqt is not None:
        gpf_args.extend(["-l", str(ligand_pdbqt)])
    gpf_args.extend(["-r", str(receptor_tz), "-o", str(gpf)])
    if box_size and center:
        npts = max(1, round(float(box_size) / 0.375))
        gpf_args.extend(["-p", f"npts={npts},{npts},{npts}"])
        gpf_args.extend(["-p", f"gridcenter={center[0]:.3f},{center[1]:.3f},{center[2]:.3f}"])
    gpf_args.extend(["-p", f"parameter_file={tools.get('AD4Zn.dat', 'AD4Zn.dat')}"])
    commands.append({"name": "write_ad4zn_gpf", "argv": gpf_args})
    commands.append(
        {
            "name": "run_autogrid4",
            "argv": [tools.get("autogrid4", "autogrid4"), "-p", str(gpf), "-l", str(glg)],
        }
    )
    if ligand_pdbqt is not None:
        commands.append(
            {
                "name": "dock_with_vina_ad4",
                "argv": [
                    tools.get("vina", "vina"),
                    "--ligand",
                    str(ligand_pdbqt),
                    "--maps",
                    str(maps_prefix),
                    "--scoring",
                    "ad4",
                    "--exhaustiveness",
                    str(int(exhaustiveness)),
                    "--out",
                    str(maps_prefix.with_name("ligand_ad4zn_out.pdbqt")),
                ],
            }
        )
    return commands


def _command_steps(plan: Mapping[str, object]) -> list[list[str]]:
    raw_commands = plan.get("commands", [])
    if not isinstance(raw_commands, list):
        return []
    commands: list[list[str]] = []
    for raw in raw_commands:
        if not isinstance(raw, Mapping):
            continue
        argv = raw.get("argv", [])
        if isinstance(argv, list) and all(isinstance(part, str) for part in argv):
            commands.append(list(argv))
    return commands


def _metal_plan_row(row: Mapping[str, object]) -> dict[str, object]:
    return {
        "id": str(row.get("id", "")),
        "element": str(row.get("element", "")),
        "coordination_number": row.get("coordination_number", ""),
        "coordination_geometry": str(row.get("coordination_geometry", "")),
        "formal_charge_primary": str(row.get("formal_charge_primary", "")),
    }


def _safe_int(value: object) -> int:
    try:
        return int(float(str(value or 0)))
    except Exception:
        return 0
