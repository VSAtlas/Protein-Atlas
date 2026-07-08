from __future__ import annotations

import argparse
import importlib
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, TextIO

from config.normalize import _to_bool
from config.tool_resolver import resolve_tool
from cli.doctor_display import (
    combo_to_json,
    display_path,
    print_doctor_report,
    row_to_json,
)


@dataclass(frozen=True)
class _DoctorRow:
    name: str
    required: bool
    ok: bool
    detail: str
    fix: str = "none"


@dataclass(frozen=True)
class _PdbqtCheck:
    path: Path
    kind: str
    ok: bool
    severity: str
    atom_count: int
    model_count: int
    issues: tuple[str, ...]
    fix: str


@dataclass(frozen=True)
class _ComboDoctorReport:
    pdb_id: str
    input_pdb: Path | None
    receptor_pdbqt: Path | None
    ligand_root: Path | None
    ligand_files_checked: int
    ligand_files_total: int
    checks: tuple[_PdbqtCheck, ...]
    missing: tuple[str, ...]


_REQUIRED_IMPORTS: tuple[tuple[str, str, bool], ...] = (
    ("yaml", "yaml", True),
    ("pandas", "pandas", True),
    ("rdkit", "rdkit", True),
    ("meeko", "meeko", True),
    ("numpy", "numpy", True),
    ("pyarrow", "pyarrow", True),
    ("Bio", "Bio", True),
    ("openbabel", "openbabel", False),
)

_TOOL_VERIFY_TARGETS: tuple[tuple[str, str, bool], ...] = (
    ("VINA_EXE", "vina", True),
    ("OPENBABEL_PATH", "obabel", False),
    ("P2RANK_PATH", "prank", True),
    ("SCORCH", "scorch.py", False),
)


_PDBQT_ALLOWED_RECORDS = {
    "ATOM",
    "HETATM",
    "TER",
    "END",
    "REMARK",
    "ROOT",
    "ENDROOT",
    "BRANCH",
    "ENDBRANCH",
    "TORSDOF",
    "MODEL",
    "ENDMDL",
}


def cli_wants_doctor(argv: list[str]) -> bool:
    if "--doctor" in argv:
        return True
    return len(argv) > 1 and str(argv[1]).strip().lower() == "doctor"


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="atlas doctor",
        description="Check Atlas setup and optionally validate a target/ligand docking input combo.",
    )
    parser.add_argument(
        "target",
        nargs="?",
        help="Optional PDB ID shorthand; equivalent to --pdb.",
    )
    parser.add_argument("--pdb", "--target", dest="pdb_id")
    parser.add_argument(
        "--receptor-pdbqt",
        help="Prepared receptor PDBQT to validate instead of resolving from --pdb.",
    )
    parser.add_argument(
        "--input-pdb",
        help="Raw input receptor PDB to require/check instead of resolving INPUT_DIR/<PDB>.pdb.",
    )
    parser.add_argument(
        "--ligands",
        "--library",
        dest="ligand_library",
        help="Prepared ligand library token or subdirectory, e.g. fda, chembl, hmdb.",
    )
    parser.add_argument(
        "--ligand-dir",
        help="Directory containing prepared ligand PDBQT files.",
    )
    parser.add_argument(
        "--ligand-file",
        action="append",
        default=[],
        help="Specific ligand PDBQT file to validate. Repeat for multiple files.",
    )
    parser.add_argument(
        "--limit-ligands",
        type=int,
        default=50,
        help="Maximum ligand PDBQT files to inspect from a directory/library.",
    )
    parser.add_argument(
        "--config",
        help="Optional Vina config file; referenced receptor/ligand paths are checked.",
    )
    parser.add_argument("--json", action="store_true", dest="as_json")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Treat warnings as a failing doctor result.",
    )
    return parser


def _parse_doctor_args(argv: list[str]) -> argparse.Namespace:
    args = list(argv[1:])
    if args and args[0].strip().lower() == "doctor":
        args = args[1:]
    elif "--doctor" in args:
        args.remove("--doctor")
    parser = _build_arg_parser()
    return parser.parse_args(args)


def _detect_env_name() -> str:
    for key in ("CONDA_DEFAULT_ENV", "MAMBA_DEFAULT_ENV", "VIRTUAL_ENV"):
        raw = str(os.environ.get(key, "")).strip()
        if not raw:
            continue
        if key == "VIRTUAL_ENV":
            return Path(raw).name or raw
        return raw
    conda_prefix = str(os.environ.get("CONDA_PREFIX", "")).strip()
    if conda_prefix:
        return Path(conda_prefix).name or conda_prefix
    return "unknown"


def _tool_specs_for_mode(cfg: Mapping[str, Any]) -> list[tuple[str, str, bool]]:
    specs: list[tuple[str, str, bool]] = []
    scorch_required = _to_bool(cfg.get("USE_SCORCH", False))
    for key, fallback_cmd, required in _TOOL_VERIFY_TARGETS:
        required_now = required or (key == "SCORCH" and scorch_required)
        specs.append((key, fallback_cmd, required_now))
    return specs


def _check_python_version() -> _DoctorRow:
    version = sys.version_info
    ok = (version.major, version.minor) >= (3, 10)
    detail = f"{version.major}.{version.minor}.{version.micro}"
    fix = "none" if ok else "use the docking-env Python 3.10 environment"
    return _DoctorRow(
        name="python>=3.10",
        required=True,
        ok=ok,
        detail=detail,
        fix=fix,
    )


def _check_import(module_name: str, label: str, required: bool) -> _DoctorRow:
    try:
        importlib.import_module(module_name)
    except Exception as exc:
        return _DoctorRow(
            name=f"import:{label}",
            required=required,
            ok=False,
            detail=str(exc),
            fix=_import_fix(module_name, required),
        )
    return _DoctorRow(
        name=f"import:{label}",
        required=required,
        ok=True,
        detail="ok",
    )


def _check_tool_resolution(cfg: Mapping[str, Any]) -> list[_DoctorRow]:
    rows: list[_DoctorRow] = []
    for key, fallback_cmd, required in _tool_specs_for_mode(cfg):
        resolved = resolve_tool(cfg, key, fallback_cmd)
        resolved_path = str(resolved.get("resolved_path", "")).strip()
        source = str(resolved.get("source", "missing") or "missing")
        ok = bool(resolved_path)
        detail = f"{display_path(resolved_path)} ({source})" if resolved_path else "unresolved"
        fix = "none" if ok else _tool_fix(key, required)
        rows.append(
            _DoctorRow(
                name=f"tool:{key}",
                required=required,
                ok=ok,
                detail=detail,
                fix=fix,
            )
        )
    return rows


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _resolve_path(raw: str | os.PathLike[str], *, base: Path | None = None) -> Path:
    path = Path(raw).expanduser()
    if path.is_absolute():
        return path
    return (base or _repo_root()) / path


def _cfg_path(cfg: Mapping[str, Any], key: str, default: str) -> Path:
    raw = str(cfg.get(key) or default)
    return _resolve_path(raw)


def _resolve_pdb_id(args: argparse.Namespace) -> str:
    raw = args.pdb_id or args.target or ""
    return str(raw).strip().upper()


def _resolve_input_pdb(
    cfg: Mapping[str, Any],
    args: argparse.Namespace,
    pdb_id: str,
) -> Path | None:
    if args.input_pdb:
        return _resolve_path(args.input_pdb)
    if not pdb_id:
        return None
    return _cfg_path(cfg, "INPUT_DIR", "input_pdbs") / f"{pdb_id}.pdb"


def _resolve_receptor_pdbqt(
    cfg: Mapping[str, Any],
    args: argparse.Namespace,
    pdb_id: str,
) -> Path | None:
    if args.receptor_pdbqt:
        return _resolve_path(args.receptor_pdbqt)
    if not pdb_id:
        return None
    processed_root = _cfg_path(cfg, "OUTPUT_DIR", "outputs/processed_pdbs")
    candidates = [
        processed_root / pdb_id / "receptor" / f"{pdb_id}.pdbqt",
        processed_root / pdb_id / "HOLO" / "receptor" / f"{pdb_id}.pdbqt",
        processed_root / pdb_id / "APO" / "receptor" / f"{pdb_id}.pdbqt",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def _resolve_ligand_root(
    cfg: Mapping[str, Any],
    args: argparse.Namespace,
) -> Path | None:
    if args.ligand_dir:
        return _resolve_path(args.ligand_dir)
    if not args.ligand_library:
        return None
    try:
        from prep_ligands.ligand_library_manager import paths_for_source, resolve_source

        return paths_for_source(dict(cfg), resolve_source(args.ligand_library)).library_dir
    except Exception:
        prepped_root = _cfg_path(cfg, "PREPPED_LIGANDS_DIR", "prepped_ligands")
        return prepped_root / str(args.ligand_library).strip()


def _collect_ligand_files(args: argparse.Namespace, ligand_root: Path | None) -> list[Path]:
    files = [_resolve_path(path) for path in args.ligand_file]
    if files:
        return files
    if ligand_root is None or not ligand_root.exists():
        return []
    limit = max(1, int(args.limit_ligands or 50))
    try:
        return sorted(ligand_root.rglob("*.pdbqt"))[:limit]
    except OSError:
        return []


def _count_ligand_files(ligand_root: Path | None, explicit_files: list[str]) -> int:
    if explicit_files:
        return len(explicit_files)
    if ligand_root is None or not ligand_root.exists():
        return 0
    try:
        return sum(1 for _path in ligand_root.rglob("*.pdbqt"))
    except OSError:
        return 0


def _parse_vina_config(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        values[key.strip().lower()] = value.strip().strip('"').strip("'")
    return values


def _config_path_issue(config_path: Path, key: str, value: str) -> str | None:
    if not value:
        return None
    candidate = _resolve_path(value, base=config_path.parent)
    if candidate.exists():
        return None
    return f"config {key} path does not exist: {display_path(candidate)}"


def _same_path(left: Path, right: Path) -> bool:
    try:
        return left.resolve() == right.resolve()
    except OSError:
        return left.absolute() == right.absolute()


def _check_config_paths(
    args: argparse.Namespace,
    *,
    receptor_pdbqt: Path | None,
    ligand_root: Path | None,
    ligand_files: list[Path],
) -> tuple[str, ...]:
    if not args.config:
        return ()
    config_path = _resolve_path(args.config)
    if not config_path.exists():
        return (f"config file does not exist: {display_path(config_path)}",)
    values = _parse_vina_config(config_path)
    return (
        *_config_missing_path_issues(config_path, values),
        *_config_receptor_mismatch(config_path, values, receptor_pdbqt),
        *_config_ligand_mismatch(config_path, values, ligand_root, ligand_files),
    )


def _config_missing_path_issues(
    config_path: Path,
    values: Mapping[str, str],
) -> tuple[str, ...]:
    issues = []
    for key in ("receptor", "ligand", "flex"):
        issue = _config_path_issue(config_path, key, values.get(key, ""))
        if issue:
            issues.append(issue)
    return tuple(issues)


def _config_receptor_mismatch(
    config_path: Path,
    values: Mapping[str, str],
    receptor_pdbqt: Path | None,
) -> tuple[str, ...]:
    receptor_value = values.get("receptor", "")
    if not receptor_value or receptor_pdbqt is None:
        return ()
    config_receptor = _resolve_path(receptor_value, base=config_path.parent)
    if not config_receptor.exists() or _same_path(config_receptor, receptor_pdbqt):
        return ()
    return (
        "config receptor mismatch: "
        f"{display_path(config_receptor)} does not match resolved receptor {display_path(receptor_pdbqt)}",
    )


def _config_ligand_mismatch(
    config_path: Path,
    values: Mapping[str, str],
    ligand_root: Path | None,
    ligand_files: list[Path],
) -> tuple[str, ...]:
    ligand_value = values.get("ligand", "")
    if not ligand_value:
        return ()
    config_ligand = _resolve_path(ligand_value, base=config_path.parent)
    if not config_ligand.exists():
        return ()
    if ligand_files:
        return _config_selected_ligand_mismatch(config_ligand, ligand_files)
    if ligand_root is not None:
        return _config_ligand_root_mismatch(config_ligand, ligand_root)
    return ()


def _config_selected_ligand_mismatch(
    config_ligand: Path,
    ligand_files: list[Path],
) -> tuple[str, ...]:
    selected = any(_same_path(config_ligand, ligand) for ligand in ligand_files)
    if selected:
        return ()
    return (
        "config ligand mismatch: "
        f"{display_path(config_ligand)} is not one of the selected ligand files",
    )


def _config_ligand_root_mismatch(
    config_ligand: Path,
    ligand_root: Path,
) -> tuple[str, ...]:
    try:
        config_ligand.relative_to(ligand_root.resolve())
    except (OSError, ValueError):
        return (
            "config ligand mismatch: "
            f"{display_path(config_ligand)} is outside resolved ligand root {display_path(ligand_root)}",
        )
    return ()


def _check_pdbqt_file(path: Path, *, kind: str) -> _PdbqtCheck:
    if not path.exists():
        return _pdbqt_error_check(path, kind, f"{kind} PDBQT does not exist")
    if not path.is_file():
        return _pdbqt_error_check(
            path,
            kind,
            f"{kind} PDBQT is not a file",
            fix="provide a concrete .pdbqt file path",
        )
    lines, read_error = _read_pdbqt_lines(path)
    if read_error:
        return _pdbqt_error_check(
            path,
            kind,
            read_error,
            fix="check file permissions or regenerate the PDBQT",
        )
    atom_count, model_count, bad_records = _scan_pdbqt_lines(lines)
    issues = _pdbqt_issues(atom_count, model_count, bad_records)
    severity = _pdbqt_severity(issues)
    return _PdbqtCheck(
        path=path,
        kind=kind,
        ok=severity == "pass",
        severity=severity,
        atom_count=atom_count,
        model_count=model_count,
        issues=tuple(issues),
        fix=_pdbqt_fix(kind, issues),
    )


def _pdbqt_error_check(
    path: Path,
    kind: str,
    issue: str,
    *,
    fix: str | None = None,
) -> _PdbqtCheck:
    return _PdbqtCheck(
        path=path,
        kind=kind,
        ok=False,
        severity="error",
        atom_count=0,
        model_count=0,
        issues=(issue,),
        fix=fix or _missing_pdbqt_fix(kind),
    )


def _read_pdbqt_lines(path: Path) -> tuple[list[str], str]:
    try:
        return path.read_text(encoding="utf-8", errors="replace").splitlines(), ""
    except OSError as exc:
        return [], f"cannot read file: {exc}"


def _scan_pdbqt_lines(lines: list[str]) -> tuple[int, int, list[str]]:
    atom_count = 0
    model_count = 0
    bad_records: list[str] = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        record = stripped.split(maxsplit=1)[0].upper()
        if record in {"ATOM", "HETATM"}:
            atom_count += 1
        elif record == "MODEL":
            model_count += 1
        if record not in _PDBQT_ALLOWED_RECORDS and record not in bad_records:
            bad_records.append(record)
    return atom_count, model_count, bad_records


def _pdbqt_issues(
    atom_count: int,
    model_count: int,
    bad_records: list[str],
) -> list[str]:
    issues: list[str] = []
    if bad_records:
        issues.append("unknown PDBQT records: " + ", ".join(bad_records[:8]))
    if model_count > 1:
        issues.append("unexpected multi-MODEL PDBQT input")
    if atom_count == 0:
        issues.append("no atoms found")
    elif atom_count <= 1:
        issues.append("one-atom ligand/receptor is likely a failed conversion")
    return issues


def _pdbqt_severity(issues: list[str]) -> str:
    if not issues:
        return "pass"
    if any(
        token in issue
        for issue in issues
        for token in ("unknown PDBQT records", "multi-MODEL", "no atoms", "one-atom")
    ):
        return "error"
    return "warn"


def _missing_pdbqt_fix(kind: str) -> str:
    if kind == "receptor":
        return "prepare the target first, e.g. `atlas --pdb <PDB> --no-docking`"
    return "run `atlas ligands install <library>` or provide --ligand-file"


def _pdbqt_fix(kind: str, issues: list[str]) -> str:
    joined = " ".join(issues).lower()
    if not issues:
        return "none"
    if "multi-model" in joined:
        return "split poses with vina_split or provide one prepared ligand PDBQT per file"
    if "unknown pdbqt records" in joined:
        return "regenerate with Meeko; avoid raw PDB/Open Babel PDBQT shortcuts"
    if "one-atom" in joined or "no atoms" in joined:
        return _missing_pdbqt_fix(kind)
    return "inspect or regenerate the PDBQT"


def _build_combo_report(
    cfg: Mapping[str, Any],
    args: argparse.Namespace,
) -> _ComboDoctorReport:
    pdb_id = _resolve_pdb_id(args)
    input_pdb = _resolve_input_pdb(cfg, args, pdb_id)
    receptor_pdbqt = _resolve_receptor_pdbqt(cfg, args, pdb_id)
    ligand_root = _resolve_ligand_root(cfg, args)
    ligand_files = _collect_ligand_files(args, ligand_root)
    ligand_total = _count_ligand_files(ligand_root, list(args.ligand_file))
    missing: list[str] = []
    checks: list[_PdbqtCheck] = []

    if input_pdb is not None and not input_pdb.exists():
        missing.append(f"input PDB not found: {display_path(input_pdb)}")
    if receptor_pdbqt is not None:
        checks.append(_check_pdbqt_file(receptor_pdbqt, kind="receptor"))
    if args.ligand_library or args.ligand_dir:
        if ligand_root is None or not ligand_root.exists():
            missing.append(f"ligand library directory not found: {display_path(ligand_root)}")
        elif ligand_total == 0:
            missing.append(f"no ligand PDBQT files found under: {display_path(ligand_root)}")
    for ligand_file in ligand_files:
        checks.append(_check_pdbqt_file(ligand_file, kind="ligand"))
    missing.extend(
        _check_config_paths(
            args,
            receptor_pdbqt=receptor_pdbqt,
            ligand_root=ligand_root,
            ligand_files=ligand_files,
        )
    )

    return _ComboDoctorReport(
        pdb_id=pdb_id,
        input_pdb=input_pdb,
        receptor_pdbqt=receptor_pdbqt,
        ligand_root=ligand_root,
        ligand_files_checked=len(ligand_files),
        ligand_files_total=ligand_total,
        checks=tuple(checks),
        missing=tuple(missing),
    )


def _doctor_args_include_combo(args: argparse.Namespace) -> bool:
    return bool(
        args.target
        or args.pdb_id
        or args.receptor_pdbqt
        or args.input_pdb
        or args.ligand_library
        or args.ligand_dir
        or args.ligand_file
        or args.config
    )


def _import_fix(module_name: str, required: bool) -> str:
    conda_names = {
        "yaml": "pyyaml",
        "Bio": "biopython",
        "rdkit": "rdkit",
        "openbabel": "openbabel",
        "meeko": "meeko",
    }
    package = conda_names.get(module_name, module_name)
    prefix = "rerun one-command installer" if required else "optional"
    return f"{prefix}; install conda-forge::{package} if needed"


def _tool_fix(key: str, required: bool) -> str:
    fixes = {
        "VINA_EXE": "rerun one-command installer, or set VINA_EXE=/path/to/vina",
        "OPENBABEL_PATH": (
            "optional; rerun installer, or set OPENBABEL_PATH=/path/to/obabel"
        ),
        "P2RANK_PATH": "BYOL/optional: register P2Rank; set P2RANK_PATH",
        "SCORCH": "optional: register SCORCH if rescoring is enabled",
    }
    fallback = "set the matching config.txt path"
    if not required and key not in fixes:
        fallback = "optional; configure only if you need this feature"
    return fixes.get(key, fallback)


def _combo_has_failures(report: _ComboDoctorReport, *, strict: bool) -> bool:
    if report.missing:
        return True
    for check in report.checks:
        if check.severity == "error":
            return True
        if strict and check.severity == "warn":
            return True
    return False


def _build_doctor_rows(cfg: Mapping[str, Any]) -> list[_DoctorRow]:
    rows: list[_DoctorRow] = [_check_python_version()]
    for module_name, label, required in _REQUIRED_IMPORTS:
        rows.append(_check_import(module_name, label, required))
    rows.extend(_check_tool_resolution(cfg))
    return rows


def _required_failures(rows: list[_DoctorRow]) -> list[_DoctorRow]:
    return [row for row in rows if row.required and not row.ok]


def _doctor_payload(
    rows: list[_DoctorRow],
    combo_report: _ComboDoctorReport | None,
    *,
    strict: bool,
) -> dict[str, Any]:
    combo_failed = combo_report is not None and _combo_has_failures(
        combo_report,
        strict=strict,
    )
    return {
        "python_executable": sys.executable,
        "environment": _detect_env_name(),
        "checks": [row_to_json(row) for row in rows],
        "combo": combo_to_json(combo_report) if combo_report is not None else None,
        "ok": not _required_failures(rows) and not combo_failed,
    }


def _print_doctor_report(
    rows: list[_DoctorRow],
    combo_report: _ComboDoctorReport | None,
    stream: TextIO,
) -> None:
    print_doctor_report(
        rows,
        combo_report,
        env_name=_detect_env_name(),
        stream=stream,
    )


def _doctor_exit_code(
    rows: list[_DoctorRow],
    combo_report: _ComboDoctorReport | None,
    *,
    strict: bool,
    stream: TextIO | None = None,
) -> int:
    required_failures = _required_failures(rows)
    if required_failures:
        if stream is not None:
            print(
                f"doctor result: FAIL ({len(required_failures)} required checks)",
                file=stream,
            )
        return 2
    if combo_report is not None and _combo_has_failures(combo_report, strict=strict):
        if stream is not None:
            print("doctor result: FAIL (target/ligand preflight)", file=stream)
        return 2
    if stream is not None:
        print("doctor result: PASS", file=stream)
    return 0


def run_doctor(cfg: Mapping[str, Any], argv: list[str], stream: TextIO | None = None) -> int:
    out = stream or sys.stdout
    args = _parse_doctor_args(argv)
    rows = _build_doctor_rows(cfg)
    combo_report = _build_combo_report(cfg, args) if _doctor_args_include_combo(args) else None
    if args.as_json:
        payload = _doctor_payload(rows, combo_report, strict=bool(args.strict))
        print(json.dumps(payload, indent=2, sort_keys=True), file=out)
        return 0 if payload["ok"] else 2

    _print_doctor_report(rows, combo_report, out)
    return _doctor_exit_code(rows, combo_report, strict=bool(args.strict), stream=out)
