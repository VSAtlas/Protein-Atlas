from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from post_docking.mmgbsa._atomic_io import write_csv_rows_atomic, write_json_atomic
from post_docking.mmgbsa.mmgbsa_wang2016_benchmark import (
    WANG2016_TABLE1,
    WANG2016_TABLE8,
    _protocol_gap_report,
    write_representative_complexes,
    write_static_targets,
)


def _pose_source(stage_dir: object, ligand_id: object) -> str:
    text = f"{stage_dir} {ligand_id}".lower()
    if "crystal" in text or "cocrystal" in text or "co_crystal" in text:
        return "crystal_pose"
    return "redocked_pose"


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return payload if isinstance(payload, dict) else {}


def _family_for_pdb(pdb_id: object, representative_csv: Path) -> str:
    if not representative_csv.exists():
        return ""
    target = str(pdb_id or "").strip().upper()
    for line in representative_csv.read_text(encoding="utf-8").splitlines()[1:]:
        parts = [part.strip() for part in line.split(",")]
        if len(parts) >= 2 and parts[0].upper() == target:
            return parts[1]
    return ""


def _md_status_for_replicate(replicate: Mapping[str, Any]) -> dict[str, Any]:
    status_path = Path(str(replicate.get("status_json", "")))
    if not status_path.exists():
        return {}
    return _read_json(status_path)


def _qc_paths(md_status: Mapping[str, Any]) -> dict[str, Any]:
    result = md_status.get("result", {})
    if not isinstance(result, Mapping):
        return {}
    qc = result.get("qc", {})
    if isinstance(qc, Mapping):
        return dict(qc)
    explicit = result.get("md", {})
    if isinstance(explicit, Mapping):
        nested = explicit.get("qc", {})
        if isinstance(nested, Mapping):
            return dict(nested)
    return {}


def _summary_rows(
    *,
    staged_manifest: Mapping[str, Any],
    representative_csv: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    result_rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    for aggregate in staged_manifest.get("aggregate_jobs", []):
        if not isinstance(aggregate, Mapping):
            continue
        summary_path = Path(str(aggregate.get("summary_json", "")))
        summary = _read_json(summary_path) if summary_path.exists() else {}
        ligand = aggregate.get("ligand", {})
        if not isinstance(ligand, Mapping):
            ligand = {}
        best = summary.get("best_replicate", {})
        if not isinstance(best, Mapping):
            best = {}
        pdb_id = ligand.get("pdb_id", "")
        family = _family_for_pdb(pdb_id, representative_csv)
        row = {
            "pdb_id": pdb_id,
            "family": family,
            "ligand_id": ligand.get("ligand_id", ""),
            "ligand_stem": ligand.get("ligand_stem", ""),
            "stage_dir": ligand.get("stage_dir", ""),
            "pose_source": _pose_source(ligand.get("stage_dir", ""), ligand.get("ligand_id", "")),
            "validation_tier": aggregate.get("validation_tier", ""),
            "ok": summary.get("ok", False),
            "delta_total": summary.get("agg_delta", ""),
            "frames_used": best.get("frames_used", ""),
            "frame_mean": best.get("frame_mean", ""),
            "frame_sd": best.get("frame_sd", ""),
            "frame_sem": best.get("frame_sem", ""),
            "frame_ci95": best.get("frame_ci95", ""),
            "frame_qc_pass": best.get("qc_pass", ""),
            "frame_qc_json": best.get("qc_json", ""),
            "results_csv": best.get("results_csv", ""),
            "summary_json": str(summary_path),
        }
        result_rows.append(row)
        for replicate in aggregate.get("replicates", []):
            if not isinstance(replicate, Mapping):
                continue
            md_status = _md_status_for_replicate(replicate)
            qc = _qc_paths(md_status)
            qc_rows.append(
                {
                    "pdb_id": pdb_id,
                    "family": family,
                    "ligand_id": ligand.get("ligand_id", ""),
                    "replicate": replicate.get("replicate", ""),
                    "md_status_json": replicate.get("status_json", ""),
                    "md_ok": md_status.get("ok", ""),
                    "complex_rmsd": qc.get("complex_rmsd", ""),
                    "complex_rmsd_svg": qc.get("complex_rmsd_svg", ""),
                    "ligand_rmsd": qc.get("ligand_rmsd", ""),
                    "ligand_rmsd_svg": qc.get("ligand_rmsd_svg", ""),
                    "protein_rmsf": qc.get("protein_rmsf", ""),
                    "protein_rmsf_svg": qc.get("protein_rmsf_svg", ""),
                    "cpptraj_qc_ok": qc.get("ok", ""),
                    "cpptraj_qc_log": qc.get("log", ""),
                }
            )
    return result_rows, qc_rows


def _target_rows(results: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    families = {str(row.get("family", "")) for row in results if row.get("family")}
    rows: list[dict[str, Any]] = []
    for family, setting, rmsd, pearson_r, slope in WANG2016_TABLE1:
        if family in families:
            rows.append(
                {
                    "family": family,
                    "benchmark_kind": "relative_binding_affinity",
                    "setting": setting,
                    "wang_rmsd_kcal_mol": rmsd,
                    "wang_pearson_r": pearson_r,
                    "wang_slope": slope,
                    "atlas_lite_note": "family-level target; not a per-complex energy reference",
                }
            )
    for family, setting, rmsd in WANG2016_TABLE8:
        if family in families:
            rows.append(
                {
                    "family": family,
                    "benchmark_kind": "absolute_binding_affinity",
                    "setting": setting,
                    "wang_rmsd_kcal_mol": rmsd,
                    "wang_pearson_r": "",
                    "wang_slope": "",
                    "atlas_lite_note": "family-level target; not a per-complex energy reference",
                }
            )
    return rows


def _float_or_none(value: object) -> float | None:
    try:
        text = str(value).strip()
        return float(text) if text else None
    except (TypeError, ValueError):
        return None


def _completed_pose_rows(
    results: Sequence[Mapping[str, Any]],
) -> dict[str, dict[str, Mapping[str, Any]]]:
    by_pdb: dict[str, dict[str, Mapping[str, Any]]] = {}
    for row in results:
        if str(row.get("ok", "")).lower() != "true":
            continue
        pdb_id = str(row.get("pdb_id", "")).upper()
        pose_source = str(row.get("pose_source", ""))
        if not pdb_id or pose_source not in {"crystal_pose", "redocked_pose"}:
            continue
        by_pdb.setdefault(pdb_id, {})[pose_source] = row
    return by_pdb


def _pose_delta_row(
    pdb_id: str,
    crystal: Mapping[str, Any],
    redocked: Mapping[str, Any],
) -> dict[str, Any]:
    crystal_delta = _float_or_none(crystal.get("delta_total"))
    redocked_delta = _float_or_none(redocked.get("delta_total"))
    delta_shift: float | str = ""
    if crystal_delta is not None and redocked_delta is not None:
        delta_shift = redocked_delta - crystal_delta
    return {
        "pdb_id": pdb_id,
        "family": crystal.get("family") or redocked.get("family", ""),
        "crystal_ligand_id": crystal.get("ligand_id", ""),
        "redocked_ligand_id": redocked.get("ligand_id", ""),
        "crystal_delta_total": crystal_delta if crystal_delta is not None else "",
        "redocked_delta_total": redocked_delta if redocked_delta is not None else "",
        "redock_minus_crystal_delta_total": delta_shift,
        "crystal_frames_used": crystal.get("frames_used", ""),
        "redocked_frames_used": redocked.get("frames_used", ""),
        "crystal_results_csv": crystal.get("results_csv", ""),
        "redocked_results_csv": redocked.get("results_csv", ""),
    }


def _pose_comparison_rows(results: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for pdb_id, rows in sorted(_completed_pose_rows(results).items()):
        crystal = rows.get("crystal_pose")
        redocked = rows.get("redocked_pose")
        if not crystal or not redocked:
            continue
        out.append(_pose_delta_row(pdb_id, crystal, redocked))
    return out


def _crystal_pose_row(row: Mapping[str, str], run_id: str) -> dict[str, str] | None:
    ligand_id = _crystal_ligand_id(row)
    source_sdf, source_mol2 = _crystal_source_paths(row)
    if not _crystal_row_is_usable(row, ligand_id, source_sdf, source_mol2):
        return None
    out_row = dict(row)
    out_row["run_id"] = run_id or str(row.get("run_id", "") or "").strip()
    out_row["ligand_id"] = f"{ligand_id}_crystal"
    out_row["stage_dir"] = "crystal"
    out_row["source_sdf"] = source_sdf
    out_row["source_mol2"] = source_mol2
    out_row["score_note"] = "crystal_pose_pdbbind"
    return out_row


def _crystal_ligand_id(row: Mapping[str, str]) -> str:
    return str(row.get("ligand_id", "") or row.get("ligand_base", "") or "").strip()


def _crystal_source_paths(row: Mapping[str, str]) -> tuple[str, str]:
    return (
        str(row.get("pdbbind_sdf", "") or "").strip(),
        str(row.get("pdbbind_mol2", "") or "").strip(),
    )


def _crystal_row_is_usable(
    row: Mapping[str, str],
    ligand_id: str,
    source_sdf: str,
    source_mol2: str,
) -> bool:
    pdb_id = str(row.get("pdb_id", "") or "").strip()
    return bool(pdb_id and ligand_id and (source_sdf or source_mol2))


def _crystal_report_fields(fieldnames: Sequence[str]) -> list[str]:
    output_fields = list(fieldnames)
    for field in (
        "source_mol2",
        "score_note",
        "run_id",
        "ligand_id",
        "stage_dir",
        "source_sdf",
    ):
        if field not in output_fields:
            output_fields.append(field)
    return output_fields


def write_wang_crystal_pose_report(
    *, report: str | Path, out_csv: str | Path, run_id: str = ""
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    with Path(report).open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            out_row = _crystal_pose_row(row, run_id)
            if out_row:
                rows.append(out_row)
    output_fields = _crystal_report_fields(fieldnames)
    out_path = Path(out_csv)
    ordered_rows = [{field: row.get(field, "") for field in output_fields} for row in rows]
    write_csv_rows_atomic(out_path, ordered_rows)
    return {"ok": bool(rows), "out_csv": str(out_path), "rows": len(rows)}


def build_wang_lite_report(*, staged_dir: str | Path, out_dir: str | Path) -> dict[str, Any]:
    return build_wang_lite_combined_report(staged_dirs=[staged_dir], out_dir=out_dir)


def build_wang_lite_combined_report(
    *, staged_dirs: Sequence[str | Path], out_dir: str | Path
) -> dict[str, Any]:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    target_paths = write_static_targets(out_path)
    representative_csv = Path(write_representative_complexes(out_path))
    result_rows: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    manifests: list[str] = []
    staged_paths = [Path(path) for path in staged_dirs]
    for staged_path in staged_paths:
        manifest_path = staged_path / "staged_manifest.json"
        manifests.append(str(manifest_path))
        manifest = _read_json(manifest_path)
        rows, qc = _summary_rows(
            staged_manifest=manifest,
            representative_csv=representative_csv,
        )
        result_rows.extend(rows)
        qc_rows.extend(qc)
    results_csv = out_path / "wang_lite_results.csv"
    qc_csv = out_path / "wang_lite_qc_index.csv"
    target_csv = out_path / "wang_lite_family_targets.csv"
    pose_comparison_csv = out_path / "wang_lite_crystal_vs_redocked.csv"
    write_csv_rows_atomic(results_csv, result_rows)
    write_csv_rows_atomic(qc_csv, qc_rows)
    ok_rows = [row for row in result_rows if str(row.get("ok")).lower() == "true"]
    write_csv_rows_atomic(target_csv, _target_rows(ok_rows))
    write_csv_rows_atomic(pose_comparison_csv, _pose_comparison_rows(ok_rows))
    payload = {
        "ok": bool(ok_rows),
        "schema": "atlas_mmgbsa_wang_lite_gate_v1",
        "staged_dirs": [str(path) for path in staged_paths],
        "staged_manifests": manifests,
        "staged_dir": str(staged_paths[0]) if len(staged_paths) == 1 else "",
        "staged_manifest": manifests[0] if len(manifests) == 1 else "",
        "results_csv": str(results_csv),
        "qc_index_csv": str(qc_csv),
        "family_targets_csv": str(target_csv),
        "pose_comparison_csv": str(pose_comparison_csv),
        "representative_complexes_csv": str(representative_csv),
        "wang_static_targets": target_paths,
        "n_results": len(result_rows),
        "n_ok": len(ok_rows),
        "protocol_gap_report": _protocol_gap_report(),
        "remaining_publication_gaps": [
            "Wang-lite is one replicate and shorter production MD than Wang 2016.",
            "Family target rows are scaffolded; the Wang supplement does not provide per-complex full MMPBSA rows.",
            "Chemistry/water/protonation review remains required before publication comparison.",
        ],
    }
    report_json = out_path / "wang_lite_gate_report.json"
    write_json_atomic(report_json, payload)
    return {**payload, "report_json": str(report_json)}


__all__ = [
    "build_wang_lite_combined_report",
    "build_wang_lite_report",
    "write_wang_crystal_pose_report",
]
