from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any, TextIO


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def path_display_bases() -> tuple[tuple[Path, str], ...]:
    return (
        (repo_root(), "."),
        (Path(sys.prefix), "$PYTHON_PREFIX"),
        (Path(os.environ.get("CONDA_PREFIX", "") or ""), "$CONDA_PREFIX"),
        (Path(os.environ.get("VIRTUAL_ENV", "") or ""), "$VIRTUAL_ENV"),
        (Path.home(), "~"),
    )


def relative_display_for_base(path: Path, base: Path, label: str) -> str | None:
    if not str(base):
        return None
    try:
        rel = path.relative_to(base.expanduser())
    except ValueError:
        return None
    rel_text = rel.as_posix()
    if label == ".":
        return rel_text or "."
    return f"{label}/{rel_text}" if rel_text else label


def abbreviate_absolute_path(path: Path) -> str:
    if path.parts[:2] in {("/", "bin"), ("/", "usr")}:
        return path.as_posix()
    tail = "/".join(path.parts[-3:])
    return f".../{tail}" if tail else path.name


def display_path(value: str | os.PathLike[str] | None) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    path = Path(text).expanduser()
    if not path.is_absolute():
        return text
    for base, label in path_display_bases():
        shown = relative_display_for_base(path, base, label)
        if shown is not None:
            return shown
    return abbreviate_absolute_path(path)


def print_table(rows: list[Any], stream: TextIO) -> None:
    name_width = max([len("check"), *(len(row.name) for row in rows)])
    required_width = len("required")
    status_width = len("status")
    detail_width = max([len("detail"), *(min(len(row.detail), 48) for row in rows)])
    print(
        f"{'check':<{name_width}}  {'required':<{required_width}}  "
        f"{'status':<{status_width}}  {'detail':<{detail_width}}  fix",
        file=stream,
    )
    print(
        f"{'-' * name_width}  {'-' * required_width}  "
        f"{'-' * status_width}  {'-' * detail_width}  {'-' * 16}",
        file=stream,
    )
    for row in rows:
        required = "yes" if row.required else "no"
        if row.ok:
            status = "PASS"
        elif row.required:
            status = "FAIL"
        else:
            status = "WARN"
        detail = row.detail if len(row.detail) <= 48 else f"{row.detail[:45]}..."
        print(
            f"{row.name:<{name_width}}  {required:<{required_width}}  "
            f"{status:<{status_width}}  {detail:<{detail_width}}  {row.fix}",
            file=stream,
        )


def print_combo_report(report: Any, stream: TextIO) -> None:
    print("", file=stream)
    print("target/ligand preflight", file=stream)
    print_combo_summary(report, stream)
    for item in report.missing:
        print(f"ERROR missing: {item}", file=stream)
    if not report.checks:
        print_no_combo_checks_hint(stream)
        return

    for check in report.checks:
        print_pdbqt_check(check, stream)


def print_combo_summary(report: Any, stream: TextIO) -> None:
    if report.pdb_id:
        print(f"pdb_id: {report.pdb_id}", file=stream)
    if report.input_pdb is not None:
        print(f"input_pdb: {display_path(report.input_pdb)}", file=stream)
    if report.receptor_pdbqt is not None:
        print(f"receptor_pdbqt: {display_path(report.receptor_pdbqt)}", file=stream)
    if report.ligand_root is not None:
        print(f"ligand_root: {display_path(report.ligand_root)}", file=stream)
        print(
            f"ligands_checked: {report.ligand_files_checked}/{report.ligand_files_total}",
            file=stream,
        )


def print_no_combo_checks_hint(stream: TextIO) -> None:
    print("No PDBQT files were selected for validation.", file=stream)
    print("hint: use --pdb/--receptor-pdbqt and --ligands/--ligand-file", file=stream)


def print_pdbqt_check(check: Any, stream: TextIO) -> None:
    status = "PASS" if check.severity == "pass" else check.severity.upper()
    print(
        f"{status} {check.kind}: {display_path(check.path)} atoms={check.atom_count} models={check.model_count}",
        file=stream,
    )
    for issue in check.issues:
        print(f"  issue: {issue}", file=stream)
    if check.fix != "none":
        print(f"  fix: {check.fix}", file=stream)


def row_to_json(row: Any) -> dict[str, Any]:
    return {
        "name": row.name,
        "required": row.required,
        "ok": row.ok,
        "detail": row.detail,
        "fix": row.fix,
    }


def check_to_json(check: Any) -> dict[str, Any]:
    return {
        "path": str(check.path),
        "kind": check.kind,
        "ok": check.ok,
        "severity": check.severity,
        "atom_count": check.atom_count,
        "model_count": check.model_count,
        "issues": list(check.issues),
        "fix": check.fix,
    }


def combo_to_json(report: Any) -> dict[str, Any]:
    return {
        "pdb_id": report.pdb_id,
        "input_pdb": str(report.input_pdb) if report.input_pdb is not None else "",
        "receptor_pdbqt": str(report.receptor_pdbqt)
        if report.receptor_pdbqt is not None
        else "",
        "ligand_root": str(report.ligand_root) if report.ligand_root is not None else "",
        "ligand_files_checked": report.ligand_files_checked,
        "ligand_files_total": report.ligand_files_total,
        "missing": list(report.missing),
        "checks": [check_to_json(check) for check in report.checks],
    }


def print_doctor_report(
    rows: list[Any],
    combo_report: Any | None,
    *,
    env_name: str,
    stream: TextIO,
) -> None:
    print("atlas doctor", file=stream)
    print(f"python_executable: {display_path(sys.executable)}", file=stream)
    print(f"environment: {env_name}", file=stream)
    print_table(rows, stream=stream)
    if combo_report is not None:
        print_combo_report(combo_report, stream)
