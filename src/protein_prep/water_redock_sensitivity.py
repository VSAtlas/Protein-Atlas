"""Water-policy redocking sensitivity checks for the prep benchmark."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Callable, Sequence

RedockRunner = Callable[..., tuple[float | None, str, str]]
WaterVariant = tuple[str, Path]


def audit_water_redock_sensitivity(
    *,
    ligand_pdb: Path,
    ligand_resname: str,
    variants: Sequence[WaterVariant],
    out_dir: Path,
    sidecar_path: Path,
    redock_runner: RedockRunner,
    box_size: float,
    exhaustiveness: int,
    redock_rmsd_max: float,
) -> dict[str, object]:
    """Run seeded redocking against water-policy receptor variants."""

    rows = [
        _run_variant(
            label,
            receptor_pdbqt,
            ligand_pdb=ligand_pdb,
            ligand_resname=ligand_resname,
            out_dir=out_dir,
            redock_runner=redock_runner,
            box_size=box_size,
            exhaustiveness=exhaustiveness,
        )
        for label, receptor_pdbqt in variants
    ]
    summary = _summarize(rows, redock_rmsd_max=redock_rmsd_max)
    sidecar_path.parent.mkdir(parents=True, exist_ok=True)
    sidecar_path.write_text(
        json.dumps(
            {
                "ligand_pdb": str(ligand_pdb),
                "summary": summary,
                "variants": rows,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    summary["water_redock_sensitivity_path"] = str(sidecar_path)
    return summary


def empty_water_redock_sensitivity_summary(status: str = "disabled") -> dict[str, object]:
    """Return stable row fields when sensitivity redocking is not run."""

    return {
        "water_redock_sensitivity_status": status,
        "water_redock_variant_count": 0,
        "water_redock_ok_variant_count": 0,
        "water_redock_best_variant": "",
        "water_redock_best_rmsd_a": "",
        "water_redock_baseline_rmsd_a": "",
        "water_redock_selected_rmsd_a": "",
        "water_redock_remove_all_rmsd_a": "",
        "water_redock_rmsd_range_a": "",
        "water_redock_selected_vs_remove_delta_a": "",
        "water_redock_selected_supported": False,
        "water_redock_sensitivity_path": "",
    }


def _run_variant(
    label: str,
    receptor_pdbqt: Path,
    *,
    ligand_pdb: Path,
    ligand_resname: str,
    out_dir: Path,
    redock_runner: RedockRunner,
    box_size: float,
    exhaustiveness: int,
) -> dict[str, object]:
    if not receptor_pdbqt.exists():
        return {
            "variant": label,
            "receptor_pdbqt": str(receptor_pdbqt),
            "status": "missing_receptor_pdbqt",
            "rmsd_a": "",
            "ligand_prep_status": "not_run",
        }
    variant_dir = out_dir / _safe_variant_dir(label)
    rmsd, status, ligand_prep_status = redock_runner(
        ligand_pdb,
        receptor_pdbqt,
        variant_dir,
        ligand_resname=ligand_resname,
        box_size=box_size,
        exhaustiveness=exhaustiveness,
    )
    return {
        "variant": label,
        "receptor_pdbqt": str(receptor_pdbqt),
        "status": status,
        "rmsd_a": _rounded_or_blank(rmsd),
        "ligand_prep_status": ligand_prep_status,
    }


def _summarize(
    rows: Sequence[dict[str, object]],
    *,
    redock_rmsd_max: float,
) -> dict[str, object]:
    ok_rows = _ok_redock_rows(rows)
    finite_values = _finite_rmsd_values(ok_rows)
    best = _best_row(ok_rows)
    selected = _variant_rmsd(ok_rows, "selected")
    baseline = _variant_rmsd(ok_rows, "baseline")
    remove_all = _variant_rmsd(ok_rows, "remove_all")
    selected_supported = _selected_supported(
        selected,
        baseline,
        remove_all,
        redock_rmsd_max=redock_rmsd_max,
    )
    status = _sensitivity_status(rows, ok_rows, selected_supported)
    return {
        "water_redock_sensitivity_status": status,
        "water_redock_variant_count": len(rows),
        "water_redock_ok_variant_count": len(ok_rows),
        "water_redock_best_variant": _best_variant(best),
        "water_redock_best_rmsd_a": _rounded_or_blank(_row_rmsd(best)),
        "water_redock_baseline_rmsd_a": _rounded_or_blank(baseline),
        "water_redock_selected_rmsd_a": _rounded_or_blank(selected),
        "water_redock_remove_all_rmsd_a": _rounded_or_blank(remove_all),
        "water_redock_rmsd_range_a": _rounded_or_blank(_rmsd_range(finite_values)),
        "water_redock_selected_vs_remove_delta_a": _rounded_or_blank(
            _delta(selected, remove_all)
        ),
        "water_redock_selected_absolute_pass": (
            selected is not None and selected <= redock_rmsd_max
        ),
        "water_redock_selected_not_worse": _selected_not_worse(
            selected,
            baseline,
            remove_all,
        ),
        "water_redock_selected_supported": selected_supported,
    }


def _ok_redock_rows(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    return [
        row
        for row in rows
        if _row_rmsd(row) is not None and row.get("status") == "ok"
    ]


def _finite_rmsd_values(rows: Sequence[dict[str, object]]) -> list[float]:
    return [value for row in rows if (value := _row_rmsd(row)) is not None]


def _best_row(rows: Sequence[dict[str, object]]) -> dict[str, object] | None:
    if not rows:
        return None
    return min(rows, key=lambda row: _row_rmsd(row) or math.inf)


def _best_variant(row: dict[str, object] | None) -> str:
    return str(row.get("variant", "")) if row else ""


def _selected_supported(
    selected: float | None,
    baseline: float | None,
    remove_all: float | None,
    *,
    redock_rmsd_max: float,
) -> bool:
    del redock_rmsd_max
    return _selected_not_worse(selected, baseline, remove_all)


def _selected_not_worse(
    selected: float | None,
    baseline: float | None,
    remove_all: float | None,
    *,
    tolerance_a: float = 0.25,
) -> bool:
    if selected is None:
        return False
    references = [value for value in (baseline, remove_all) if value is not None]
    return not references or selected <= min(references) + tolerance_a


def _sensitivity_status(
    rows: Sequence[dict[str, object]],
    ok_rows: Sequence[dict[str, object]],
    selected_supported: bool,
) -> str:
    if not rows:
        return "no_variants"
    if len(ok_rows) < 2:
        return "insufficient_successful_variants"
    if selected_supported:
        return "selected_supported_by_redocking"
    return "sensitivity_review"


def _variant_rmsd(rows: Sequence[dict[str, object]], variant: str) -> float | None:
    for row in rows:
        if row.get("variant") == variant:
            return _row_rmsd(row)
    return None


def _row_rmsd(row: dict[str, object] | None) -> float | None:
    if row is None:
        return None
    try:
        value = float(str(row.get("rmsd_a", "")))
    except Exception:
        return None
    return value if math.isfinite(value) else None


def _rmsd_range(values: Sequence[float]) -> float | None:
    if len(values) < 2:
        return None
    return max(values) - min(values)


def _delta(a: float | None, b: float | None) -> float | None:
    if a is None or b is None:
        return None
    return a - b


def _rounded_or_blank(value: float | None) -> float | str:
    if value is None or not math.isfinite(value):
        return ""
    return round(float(value), 3)


def _safe_variant_dir(label: str) -> str:
    safe = "".join(char if char.isalnum() or char in {"_", "-"} else "_" for char in label)
    return safe or "variant"
