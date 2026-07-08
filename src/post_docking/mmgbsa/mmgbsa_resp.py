from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any, Mapping

from config.tool_resolver import select_ambertools_tool
from post_docking.mmgbsa._atomic_io import write_json_atomic, write_text_atomic

GAUSSIAN_ROUTE = "# HF/6-31G* Pop=MK IOp(6/33=2,6/41=10,6/42=17) SCF=Tight"


def _write_text_atomic(path: Path, text: str) -> None:
    write_text_atomic(path, text)


def _write_json_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    write_json_atomic(path, payload)


def _select_amber_tool(cfg: object | None, tool: str) -> tuple[list[str], str, str]:
    selection = select_ambertools_tool(cfg, tool)
    return selection.runner, selection.exe, selection.source


def _suffix_file_type(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".sdf":
        return "sdf"
    if suffix == ".mol2":
        return "mol2"
    if suffix in {".pdb", ".ent"}:
        return "pdb"
    raise ValueError(f"unsupported ligand input for RESP planning: {path}")


def _qm_output_type(qm_output: Path, qm_engine: str) -> str:
    engine = qm_engine.strip().lower()
    if engine in {"gaussian", "g16", "g09"}:
        return "gout"
    if engine in {"gamess", "gamess-us"}:
        return "gamout"
    suffix = qm_output.suffix.lower()
    if suffix in {".log", ".out"}:
        return "gout"
    raise ValueError(f"unsupported QM engine/output type: {qm_engine}")


def _run_command(
    cmd: list[str],
    *,
    cwd: Path,
    log_path: Path,
) -> dict[str, Any]:
    with log_path.open("w", encoding="utf-8") as handle:
        proc = subprocess.run(
            cmd,
            cwd=str(cwd),
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return {
        "cmd": cmd,
        "cwd": str(cwd),
        "log": str(log_path),
        "returncode": proc.returncode,
        "ok": proc.returncode == 0,
    }


def _shell_join(parts: list[str]) -> str:
    return " ".join(_quote_shell(part) for part in parts)


def _quote_shell(value: str) -> str:
    if not value:
        return "''"
    safe = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_./:-=+")
    if all(ch in safe for ch in value):
        return value
    return "'" + value.replace("'", "'\"'\"'") + "'"


def _make_qm_input_command(
    *,
    antechamber: str,
    ligand_path: Path,
    qm_input: Path,
    input_type: str,
    net_charge: int,
    residue_name: str,
    atom_type: str,
) -> list[str]:
    return [
        antechamber,
        "-i",
        str(ligand_path),
        "-fi",
        input_type,
        "-o",
        str(qm_input),
        "-fo",
        "gcrt",
        "-c",
        "gas",
        "-nc",
        str(net_charge),
        "-rn",
        residue_name,
        "-at",
        atom_type,
    ]


def _finalize_resp_command(
    *,
    antechamber: str,
    qm_output: Path,
    resp_mol2: Path,
    qm_output_type: str,
    net_charge: int,
    residue_name: str,
    atom_type: str,
) -> list[str]:
    return [
        antechamber,
        "-i",
        str(qm_output),
        "-fi",
        qm_output_type,
        "-o",
        str(resp_mol2),
        "-fo",
        "mol2",
        "-c",
        "resp",
        "-nc",
        str(net_charge),
        "-rn",
        residue_name,
        "-at",
        atom_type,
    ]


def _parmchk2_command(
    *,
    parmchk2: str,
    resp_mol2: Path,
    frcmod: Path,
    atom_type: str,
) -> list[str]:
    return [
        parmchk2,
        "-i",
        str(resp_mol2),
        "-f",
        "mol2",
        "-o",
        str(frcmod),
        "-s",
        atom_type,
    ]


def _resp_output_paths(qm_path: Path, out_path: Path) -> tuple[Path, Path]:
    stem = qm_path.name
    for suffix in (".log", ".out", ".gout"):
        if stem.lower().endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return out_path / f"{stem}.resp.mol2", out_path / f"{stem}.resp.frcmod"


def _outputs_ok(resp_mol2: Path, frcmod: Path) -> bool:
    return (
        resp_mol2.exists()
        and resp_mol2.stat().st_size > 0
        and frcmod.exists()
        and frcmod.stat().st_size > 0
    )


def _amber_tool_bins(cfg: object | None) -> tuple[list[str], str, str, str, str]:
    runner, antechamber_bin, source = _select_amber_tool(cfg, "antechamber")
    _, parmchk2_bin, parmchk_source = _select_amber_tool(cfg, "parmchk2")
    antechamber = antechamber_bin if not runner else str(Path(antechamber_bin))
    parmchk2 = parmchk2_bin if not runner else str(Path(parmchk2_bin))
    return runner, antechamber, parmchk2, source, parmchk_source


def _finalize_command_pair(
    *,
    runner: list[str],
    antechamber: str,
    parmchk2: str,
    qm_path: Path,
    resp_mol2: Path,
    frcmod: Path,
    qm_engine: str,
    net_charge: int,
    residue_name: str,
    atom_type: str,
) -> tuple[list[str], list[str]]:
    finalize_cmd = list(runner) + _finalize_resp_command(
        antechamber=antechamber,
        qm_output=qm_path,
        resp_mol2=resp_mol2,
        qm_output_type=_qm_output_type(qm_path, qm_engine),
        net_charge=net_charge,
        residue_name=residue_name,
        atom_type=atom_type,
    )
    parmchk_cmd = list(runner) + _parmchk2_command(
        parmchk2=parmchk2,
        resp_mol2=resp_mol2,
        frcmod=frcmod,
        atom_type=atom_type,
    )
    return finalize_cmd, parmchk_cmd


def _run_finalize_commands(
    *,
    finalize_cmd: list[str],
    parmchk_cmd: list[str],
    out_path: Path,
    resp_mol2: Path,
    frcmod: Path,
    force: bool,
) -> list[dict[str, Any]]:
    executed: list[dict[str, Any]] = []
    if force or not resp_mol2.exists():
        executed.append(
            _run_command(
                finalize_cmd,
                cwd=out_path,
                log_path=out_path / "finalize_resp.log",
            )
        )
    if resp_mol2.exists() and resp_mol2.stat().st_size > 0 and (
        force or not frcmod.exists()
    ):
        executed.append(
            _run_command(
                parmchk_cmd,
                cwd=out_path,
                log_path=out_path / "parmchk2_resp.log",
            )
        )
    return executed


def _write_resp_scripts(
    *,
    out_dir: Path,
    make_qm_cmd: list[str],
    finalize_cmd: list[str],
    parmchk_cmd: list[str],
    qm_output: Path,
) -> dict[str, str]:
    make_script = out_dir / "01_make_qm_input.sh"
    finalize_script = out_dir / "02_finalize_resp.sh"
    _write_text_atomic(
        make_script,
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                _shell_join(make_qm_cmd),
                "",
                "# Review the generated QM input, run Gaussian/GAMESS externally,",
                f"# then place the completed output at: {qm_output}",
                "",
            ]
        ),
    )
    _write_text_atomic(
        finalize_script,
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                f"test -s {_quote_shell(str(qm_output))}",
                _shell_join(finalize_cmd),
                _shell_join(parmchk_cmd),
                "",
            ]
        ),
    )
    make_script.chmod(0o755)
    finalize_script.chmod(0o755)
    return {"make_qm_input": str(make_script), "finalize_resp": str(finalize_script)}


def plan_resp_workflow(
    *,
    ligand_path: str | Path,
    out_dir: str | Path,
    net_charge: int,
    cfg: object | None = None,
    residue_name: str = "LIG",
    atom_type: str = "gaff2",
    qm_engine: str = "gaussian",
    qm_route: str | None = GAUSSIAN_ROUTE,
    run_antechamber: bool = False,
    force: bool = False,
) -> dict[str, Any]:
    """Create a reviewed QM/RESP handoff plan from authoritative ligand input."""
    ligand = Path(ligand_path)
    if not ligand.exists():
        raise FileNotFoundError(f"ligand input not found: {ligand}")
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    input_type = _suffix_file_type(ligand)
    runner, antechamber_bin, source = _select_amber_tool(cfg, "antechamber")
    _, parmchk2_bin, parmchk_source = _select_amber_tool(cfg, "parmchk2")
    ligand_stem = ligand.stem
    qm_input = out_path / f"{ligand_stem}.resp.gcrt"
    qm_output = out_path / f"{ligand_stem}.resp.log"
    resp_mol2 = out_path / f"{ligand_stem}.resp.mol2"
    frcmod = out_path / f"{ligand_stem}.resp.frcmod"
    antechamber = antechamber_bin if not runner else str(Path(antechamber_bin))
    parmchk2 = parmchk2_bin if not runner else str(Path(parmchk2_bin))
    make_cmd = list(runner) + _make_qm_input_command(
        antechamber=antechamber,
        ligand_path=ligand,
        qm_input=qm_input,
        input_type=input_type,
        net_charge=net_charge,
        residue_name=residue_name,
        atom_type=atom_type,
    )
    finalize_cmd = list(runner) + _finalize_resp_command(
        antechamber=antechamber,
        qm_output=qm_output,
        resp_mol2=resp_mol2,
        qm_output_type=_qm_output_type(qm_output, qm_engine),
        net_charge=net_charge,
        residue_name=residue_name,
        atom_type=atom_type,
    )
    parmchk_cmd = list(runner) + _parmchk2_command(
        parmchk2=parmchk2,
        resp_mol2=resp_mol2,
        frcmod=frcmod,
        atom_type=atom_type,
    )
    scripts = _write_resp_scripts(
        out_dir=out_path,
        make_qm_cmd=make_cmd,
        finalize_cmd=finalize_cmd,
        parmchk_cmd=parmchk_cmd,
        qm_output=qm_output,
    )
    commands: list[dict[str, Any]] = []
    if run_antechamber and (force or not qm_input.exists()):
        commands.append(
            _run_command(
                make_cmd,
                cwd=out_path,
                log_path=out_path / "make_qm_input.log",
            )
        )
    payload = {
        "schema": "atlas_mmgbsa_resp_workflow_v1",
        "status": "planned",
        "ligand_input": str(ligand),
        "ligand_input_type": input_type,
        "net_charge": int(net_charge),
        "residue_name": residue_name,
        "atom_type": atom_type,
        "qm_engine": qm_engine,
        "qm_route": qm_route or GAUSSIAN_ROUTE,
        "qm_input": str(qm_input),
        "qm_output_expected": str(qm_output),
        "resp_mol2": str(resp_mol2),
        "frcmod": str(frcmod),
        "scripts": scripts,
        "commands": {
            "make_qm_input": make_cmd,
            "finalize_resp": finalize_cmd,
            "parmchk2": parmchk_cmd,
        },
        "executed": commands,
        "antechamber_source": source,
        "parmchk2_source": parmchk_source,
        "publication_notes": [
            "Review protonation, tautomer, stereochemistry, net charge, and QM settings before running the QM job.",
            "Atlas finalizes RESP from completed Gaussian/GAMESS output; it does not replace chemistry review.",
        ],
    }
    _write_json_atomic(out_path / "resp_workflow.json", payload)
    return payload


def finalize_resp_workflow(
    *,
    qm_output: str | Path,
    out_dir: str | Path,
    net_charge: int,
    cfg: object | None = None,
    residue_name: str = "LIG",
    atom_type: str = "gaff2",
    qm_engine: str = "gaussian",
    force: bool = False,
) -> dict[str, Any]:
    """Create RESP MOL2/frcmod files from a completed QM output file."""
    qm_path = Path(qm_output)
    if not qm_path.exists():
        raise FileNotFoundError(f"QM output not found: {qm_path}")
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    runner, antechamber, parmchk2, source, parmchk_source = _amber_tool_bins(cfg)
    resp_mol2, frcmod = _resp_output_paths(qm_path, out_path)
    finalize_cmd, parmchk_cmd = _finalize_command_pair(
        runner=runner,
        antechamber=antechamber,
        parmchk2=parmchk2,
        qm_path=qm_path,
        resp_mol2=resp_mol2,
        frcmod=frcmod,
        qm_engine=qm_engine,
        net_charge=net_charge,
        residue_name=residue_name,
        atom_type=atom_type,
    )
    executed = _run_finalize_commands(
        finalize_cmd=finalize_cmd,
        parmchk_cmd=parmchk_cmd,
        out_path=out_path,
        resp_mol2=resp_mol2,
        frcmod=frcmod,
        force=force,
    )
    payload = {
        "schema": "atlas_mmgbsa_resp_finalize_v1",
        "ok": _outputs_ok(resp_mol2, frcmod),
        "qm_output": str(qm_path),
        "qm_output_type": _qm_output_type(qm_path, qm_engine),
        "net_charge": int(net_charge),
        "residue_name": residue_name,
        "atom_type": atom_type,
        "resp_mol2": str(resp_mol2),
        "frcmod": str(frcmod),
        "commands": {
            "finalize_resp": finalize_cmd,
            "parmchk2": parmchk_cmd,
        },
        "executed": executed,
        "antechamber_source": source,
        "parmchk2_source": parmchk_source,
    }
    _write_json_atomic(out_path / "resp_finalize.json", payload)
    return payload
