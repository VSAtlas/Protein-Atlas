#!/usr/bin/env python3
"""Summarize large protein-prep stress runs from structured sidecars."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

from path_router.path_router import run_scoped_root


def _read_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _as_int(value: Any) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def _resolve_run_root(processed_root: Path, run_id: str) -> Path:
    scoped = run_scoped_root(processed_root.expanduser(), run_id)
    if scoped.exists():
        return scoped
    return processed_root.expanduser()


def _geometry_summary(work_dir: Path, pdb_id: str) -> dict[str, Any]:
    data = _read_json(work_dir / f"{pdb_id}_geometry_clash_audit_final.json")
    summary = data.get("summary", {})
    return summary if isinstance(summary, dict) else {}


def _residual_geometry_summary(work_dir: Path, pdb_id: str) -> dict[str, Any]:
    data = _read_json(work_dir / f"{pdb_id}_geometry_residual_escalation.json")
    return data if isinstance(data, dict) else {}


def _metal_summary(work_dir: Path) -> dict[str, int]:
    audits = sorted(work_dir.glob("*.retained_hets_audit.json"), key=lambda p: p.stat().st_mtime)
    if not audits:
        return {
            "metal_before": 0,
            "metal_after": 0,
            "metal_missing_after": 0,
            "cofactor_before": 0,
            "cofactor_after": 0,
        }
    data = _read_json(audits[-1])
    before = data.get("before", {})
    after = data.get("after", {})
    if not isinstance(before, dict):
        before = {}
    if not isinstance(after, dict):
        after = {}
    return {
        "metal_before": _as_int(before.get("metals")),
        "metal_after": _as_int(after.get("metals")),
        "metal_missing_after": _as_int(data.get("missing_after_count")),
        "cofactor_before": _as_int(before.get("cofactors")),
        "cofactor_after": _as_int(after.get("cofactors")),
    }


def _drop_summary(work_dir: Path) -> dict[str, int]:
    atoms = 0
    residues = 0
    for path in work_dir.glob("*.drop_audit.json"):
        data = _read_json(path)
        atoms += _as_int(data.get("dropped_atom_count"))
        residues += _as_int(data.get("dropped_residue_count"))
    return {"dropped_atoms": atoms, "dropped_residues": residues}


def _duplicate_metal_summary(work_dir: Path, pdb_id: str) -> dict[str, Any]:
    data = _read_json(work_dir / f"{pdb_id}_duplicate_metal_sites.json")
    return data if isinstance(data, dict) else {}


def _p2rank_summary(
    target_dir: Path,
    processed_root: Path,
    pdb_id: str,
) -> dict[str, Any]:
    candidates = [target_dir]
    direct_target = processed_root / pdb_id
    if direct_target not in candidates:
        candidates.append(direct_target)
    prediction_path: Path | None = None
    pockets_path: Path | None = None
    for candidate in candidates:
        if prediction_path is None:
            path = (
                candidate
                / "receptor"
                / "_p2rank"
                / f"{pdb_id}_cleaned.pdb_predictions.csv"
            )
            if path.exists():
                prediction_path = path
        if pockets_path is None:
            path = candidate / "pockets" / "pockets.json"
            if path.exists():
                pockets_path = path
    row_count = 0
    top_score = ""
    top_probability = ""
    top_center = ""
    if prediction_path is not None:
        try:
            with prediction_path.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                rows = [
                    {_clean_csv_key(k): v for k, v in row.items()}
                    for row in reader
                ]
            row_count = len(rows)
            if rows:
                top = rows[0]
                top_score = str(top.get("score", "")).strip()
                top_probability = str(top.get("probability", "")).strip()
                center = [
                    str(top.get("center_x", "")).strip(),
                    str(top.get("center_y", "")).strip(),
                    str(top.get("center_z", "")).strip(),
                ]
                top_center = ",".join(value for value in center if value)
        except Exception:
            row_count = 0
    pockets_count = 0
    if pockets_path is not None:
        data = _read_json(pockets_path)
        pockets = data.get("pockets", []) if isinstance(data, dict) else []
        pockets_count = len(pockets) if isinstance(pockets, list) else 0
    return {
        "p2rank_predictions": prediction_path is not None,
        "p2rank_pocket_count": row_count,
        "p2rank_pockets_json": pockets_path is not None,
        "p2rank_pockets_json_count": pockets_count,
        "p2rank_top_score": top_score,
        "p2rank_top_probability": top_probability,
        "p2rank_top_center": top_center,
        "p2rank_predictions_path": str(prediction_path) if prediction_path else "",
        "p2rank_pockets_path": str(pockets_path) if pockets_path else "",
    }


def _clean_csv_key(key: str | None) -> str:
    return "" if key is None else key.strip()


def _ligand_count(target_dir: Path) -> int:
    candidates = [target_dir / "ligands_raw"]
    if target_dir.parent != target_dir:
        candidates.append(target_dir.parent / "ligands_raw")
    for ligand_dir in candidates:
        if ligand_dir.exists():
            return sum(1 for path in ligand_dir.glob("*.pdb") if path.is_file())
    return 0


def _center_sources(configs_root: Path, pdb_id: str, variant: str = "") -> Counter[str]:
    counter: Counter[str] = Counter()
    pdb_cfg = configs_root / pdb_id
    if variant:
        variant_cfg = pdb_cfg / variant
        if variant_cfg.exists():
            pdb_cfg = variant_cfg
    if not pdb_cfg.exists():
        return counter
    for path in pdb_cfg.glob("**/*ctrl_redock*.txt"):
        counter["control_config"] += 1
    return counter


def _target_row(
    target_dir: Path,
    configs_root: Path,
    processed_root: Path,
    pdb_id: str,
    variant: str,
) -> dict[str, Any]:
    receptor_dir = target_dir / "receptor"
    work_dir = target_dir / "work"
    cleaned = receptor_dir / f"{pdb_id}_cleaned.pdb"
    pdbqt = receptor_dir / f"{pdb_id}.pdbqt"
    meeko = _read_json(receptor_dir / f"{pdb_id}.meeko_input_stage.json")
    if not meeko:
        meeko = _read_json(work_dir / "meeko_input_stage.json")
    geometry = _geometry_summary(work_dir, pdb_id)
    residual_geometry = _residual_geometry_summary(work_dir, pdb_id)
    metal = _metal_summary(work_dir)
    drops = _drop_summary(work_dir)
    duplicate_metals = _duplicate_metal_summary(work_dir, pdb_id)
    p2rank = _p2rank_summary(target_dir, processed_root, pdb_id)
    water = _read_json(work_dir / f"{pdb_id}_selected_water_policy.json")
    het = _read_json(work_dir / f"{pdb_id}_het_state_selection.json")
    center_counts = _center_sources(configs_root, pdb_id, variant=variant)
    ligand_count = _ligand_count(target_dir)
    meeko_ok = bool(meeko.get("ok")) if meeko else pdbqt.exists()
    geometry_status = str(geometry.get("geometry_status", "missing" if not geometry else "unknown"))
    residual_gate = str(residual_geometry.get("geometry_residual_gate_status", ""))
    duplicate_metal_review = bool(
        duplicate_metals.get("duplicate_metal_site_fail_for_review", False)
    )
    metal_ok = metal["metal_missing_after"] == 0
    geometry_ok = geometry_status == "ok" or residual_gate == "offsite_geometry_review"
    has_control_ligand = ligand_count > 0
    p2rank_site_ok = _as_int(p2rank.get("p2rank_pocket_count")) > 0
    site_ok = has_control_ligand or p2rank_site_ok
    passed = (
        cleaned.exists()
        and pdbqt.exists()
        and meeko_ok
        and not duplicate_metal_review
        and metal_ok
        and geometry_ok
        and site_ok
    )
    root_cause = "passed"
    if not cleaned.exists():
        root_cause = "cleaned_receptor_missing"
    elif duplicate_metal_review:
        root_cause = "duplicate_metal_site_review"
    elif not pdbqt.exists():
        root_cause = "receptor_pdbqt_missing"
    elif not meeko_ok:
        root_cause = "meeko_export_failed"
    elif not metal_ok:
        root_cause = "metal_or_cofactor_loss"
    elif geometry_status == "missing":
        root_cause = "geometry_audit_missing"
    elif not geometry_ok:
        if residual_gate == "binding_site_geometry_review":
            root_cause = "geometry_binding_site_clashes"
        elif residual_gate == "severe_geometry_review":
            root_cause = "geometry_severe_clashes"
        else:
            root_cause = f"geometry_{geometry_status}"
    elif not has_control_ligand and p2rank.get("p2rank_predictions") and not p2rank_site_ok:
        root_cause = "p2rank_no_predicted_pockets"
    elif not site_ok:
        root_cause = "control_ligand_missing"
    return {
        "pdb_id": pdb_id,
        "variant": variant,
        "passed": passed,
        "root_cause": root_cause,
        "cleaned_receptor": cleaned.exists(),
        "receptor_pdbqt": pdbqt.exists(),
        "meeko_ok": meeko_ok,
        "meeko_stage": meeko.get("stage", ""),
        "meeko_review_required": bool(meeko.get("review_required", False)),
        "meeko_publication_review_required": bool(
            meeko.get("publication_review_required", False)
        ),
        "control_ligand_count": ligand_count,
        "active_site_identification_route": "control_ligand"
        if has_control_ligand
        else "p2rank_predicted_pocket"
        if p2rank_site_ok
        else "missing",
        "p2rank_predictions": bool(p2rank.get("p2rank_predictions")),
        "p2rank_pocket_count": _as_int(p2rank.get("p2rank_pocket_count")),
        "p2rank_pockets_json": bool(p2rank.get("p2rank_pockets_json")),
        "p2rank_pockets_json_count": _as_int(
            p2rank.get("p2rank_pockets_json_count")
        ),
        "p2rank_top_score": p2rank.get("p2rank_top_score", ""),
        "p2rank_top_probability": p2rank.get("p2rank_top_probability", ""),
        "p2rank_top_center": p2rank.get("p2rank_top_center", ""),
        "duplicate_metal_site_status": duplicate_metals.get(
            "duplicate_metal_site_status", ""
        ),
        "duplicate_metal_site_review_required": bool(
            duplicate_metals.get("duplicate_metal_site_review_required", False)
        ),
        "duplicate_metal_site_cluster_count": _as_int(
            duplicate_metals.get("duplicate_metal_site_cluster_count")
        ),
        "duplicate_metal_site_dropped_atom_count": _as_int(
            duplicate_metals.get("duplicate_metal_site_dropped_atom_count")
        ),
        "geometry_status": geometry_status,
        "geometry_residual_gate_status": residual_gate,
        "geometry_residual_clash_count": _as_int(
            residual_geometry.get("geometry_residual_clash_count")
        ),
        "geometry_residual_near_ligand_clash_count": _as_int(
            residual_geometry.get("geometry_residual_near_ligand_clash_count")
        ),
        "geometry_residual_severe_clash_count": _as_int(
            residual_geometry.get("geometry_residual_severe_clash_count")
        ),
        "geometry_residual_offsite_clash_count": _as_int(
            residual_geometry.get("geometry_residual_offsite_clash_count")
        ),
        "geometry_residual_repair_status": residual_geometry.get(
            "geometry_residual_repair_status", ""
        ),
        "geometry_residual_repair_accepted": bool(
            residual_geometry.get("geometry_residual_repair_accepted", False)
        ),
        "geometry_clash_count": _as_int(geometry.get("geometry_clash_count")),
        "geometry_ignored_close_contact_count": _as_int(
            geometry.get("geometry_ignored_close_contact_count")
        ),
        "metal_before": metal["metal_before"],
        "metal_after": metal["metal_after"],
        "metal_missing_after": metal["metal_missing_after"],
        "cofactor_before": metal["cofactor_before"],
        "cofactor_after": metal["cofactor_after"],
        "dropped_atoms": drops["dropped_atoms"],
        "dropped_residues": drops["dropped_residues"],
        "water_policy_status": water.get("status", ""),
        "het_state_status": het.get("status", ""),
        "control_redock_config_count": center_counts["control_config"],
    }


def _iter_target_output_dirs(run_root: Path) -> list[tuple[Path, str, str]]:
    rows: list[tuple[Path, str, str]] = []
    variant_names = {"APO", "HOLO", "LEGACY"}
    for pdb_dir in sorted(path for path in run_root.iterdir() if path.is_dir()):
        variant_rows: list[tuple[Path, str, str]] = []
        for variant_dir in sorted(path for path in pdb_dir.iterdir() if path.is_dir()):
            variant = variant_dir.name.upper()
            if variant not in variant_names:
                continue
            receptor_dir = variant_dir / "receptor"
            work_dir = variant_dir / "work"
            if receptor_dir.exists() or work_dir.exists():
                variant_rows.append((variant_dir, pdb_dir.name, variant))
        if variant_rows:
            rows.extend(variant_rows)
            continue
        direct_receptor = pdb_dir / "receptor"
        direct_work = pdb_dir / "work"
        if direct_receptor.exists() or direct_work.exists():
            rows.append((pdb_dir, pdb_dir.name, ""))
    return rows


def build_summary(processed_root: Path, configs_root: Path, run_id: str) -> dict[str, Any]:
    expanded_processed_root = processed_root.expanduser()
    run_root = _resolve_run_root(expanded_processed_root, run_id)
    rows = [
        _target_row(
            path,
            configs_root / run_id,
            expanded_processed_root,
            pdb_id=pdb_id,
            variant=variant,
        )
        for path, pdb_id, variant in _iter_target_output_dirs(run_root)
    ]
    root_causes = Counter(str(row["root_cause"]) for row in rows)
    meeko_stages = Counter(str(row["meeko_stage"]) for row in rows)
    geometry_statuses = Counter(str(row["geometry_status"]) for row in rows)
    return {
        "run_id": run_id,
        "run_root": str(run_root),
        "target_count": len(rows),
        "passed_count": sum(1 for row in rows if row["passed"]),
        "failed_count": sum(1 for row in rows if not row["passed"]),
        "root_causes": dict(root_causes.most_common()),
        "meeko_stages": dict(meeko_stages.most_common()),
        "geometry_statuses": dict(geometry_statuses.most_common()),
        "rows": rows,
    }


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0]) if rows else ["pdb_id", "passed", "root_cause"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _print_rows(label: str, rows: list[dict[str, Any]], *, limit: int) -> None:
    if not rows:
        print(f"{label}=none")
        return
    shown = rows[:limit] if limit > 0 else rows
    names = [
        f"{row['pdb_id']}:{row.get('root_cause', '')}"
        for row in shown
    ]
    suffix = "" if len(shown) == len(rows) else f" (+{len(rows) - len(shown)} more)"
    print(f"{label}={','.join(names)}{suffix}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--processed-root", required=True, type=Path)
    parser.add_argument("--configs-root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument(
        "--print-clashes",
        action="store_true",
        help="Print PDB IDs with residual geometry clash root causes.",
    )
    parser.add_argument(
        "--print-failures",
        action="store_true",
        help="Print non-passing PDB IDs with root-cause labels.",
    )
    parser.add_argument(
        "--print-limit",
        type=int,
        default=80,
        help="Maximum IDs to print for --print-clashes/--print-failures; 0 prints all.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    summary = build_summary(args.processed_root, args.configs_root, args.run_id)
    out_dir = args.out_dir.expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / f"{args.run_id}_prep_stress_summary.json"
    csv_path = out_dir / f"{args.run_id}_prep_stress_targets.csv"
    rows = list(summary.pop("rows"))
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    _write_csv(csv_path, rows)
    print(
        "targets={target_count} passed={passed_count} failed={failed_count}".format(
            **summary
        )
    )
    print(f"summary_json={json_path}")
    print(f"target_csv={csv_path}")
    if args.print_clashes:
        clash_rows = [
            row for row in rows if str(row.get("root_cause", "")).startswith("geometry_")
        ]
        _print_rows("clash_targets", clash_rows, limit=args.print_limit)
    if args.print_failures:
        failure_rows = [row for row in rows if not row.get("passed")]
        _print_rows("failure_targets", failure_rows, limit=args.print_limit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
