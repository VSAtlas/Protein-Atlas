#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import yaml  # type: ignore[import-untyped]

from config.output_paths import run_output_dir

_REFERENCE_DELTA_LIMIT = 0.05
_BENCH2_TARGET_KEYS: tuple[str, ...] = (
    "BNJS|HOLO|pH7_0",
    "BNNQ|HOLO|pH5_6",
    "BOJG|HOLO|pH7_2",
    "BOJG|HOLO|pH7_7",
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _safe_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    token = str(value).strip().replace(",", "")
    if not token:
        return None
    try:
        return float(token)
    except Exception:
        return None


def _safe_int(value: Any) -> Optional[int]:
    if value is None:
        return None
    token = str(value).strip().replace(",", "")
    if not token:
        return None
    try:
        return int(float(token))
    except Exception:
        return None


def _read_json(path: Path) -> Mapping[str, Any]:
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8")) or {}
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _read_yaml(path: Path) -> Mapping[str, Any]:
    if not path.exists():
        return {}
    try:
        payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def _parse_iso_utc(value: Any) -> Optional[dt.datetime]:
    if value is None:
        return None
    token = str(value).strip()
    if not token:
        return None
    token = token.replace("Z", "+00:00")
    try:
        return dt.datetime.fromisoformat(token)
    except Exception:
        return None


def _fmt_iso_utc(value: Optional[dt.datetime]) -> Optional[str]:
    if value is None:
        return None
    if value.tzinfo is None:
        value = value.replace(tzinfo=dt.timezone.utc)
    return value.astimezone(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def _count_planned_chunks(distributed_dir: Path) -> int:
    total = 0
    for plan_path in sorted(distributed_dir.glob("combo_chunks_*.json")):
        payload = _read_json(plan_path)
        chunk_list = payload.get("chunks")
        if isinstance(chunk_list, list):
            count = 0
            for item in chunk_list:
                if isinstance(item, Mapping) and isinstance(item.get("chunk_id"), str):
                    count += 1
                elif isinstance(item, str):
                    count += 1
            total += count
            continue
        meta_count: Optional[int] = _safe_int(payload.get("chunk_count"))
        if meta_count is None:
            meta_count = _safe_int(payload.get("chunks_total"))
        if meta_count is not None and meta_count > 0:
            total += int(meta_count)
            continue
        try:
            text = plan_path.read_text(encoding="utf-8")
            total += int(text.count('"chunk_id"'))
        except Exception:
            pass
    return int(total)


def _count_csv_rows(path: Path) -> int:
    if not path.exists():
        return 0
    rows = 0
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle)
        try:
            next(reader)
        except StopIteration:
            return 0
        for _ in reader:
            rows += 1
    return rows


def _choose_ligand_column(
    fieldnames: Sequence[str],
    preferred: Sequence[str],
) -> Optional[str]:
    fields = list(fieldnames)
    lower_map = {f.lower(): f for f in fields}
    for cand in preferred:
        hit = lower_map.get(cand.lower())
        if hit is not None:
            return hit
    return fields[0] if fields else None


def _duplicate_stats_for_pattern(
    root: Path,
    pattern: str,
    preferred_cols: Sequence[str],
) -> Dict[str, int]:
    files = sorted(root.rglob(pattern))
    total_rows = 0
    duplicate_rows = 0
    files_with_duplicates = 0
    for path in files:
        with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames:
                continue
            ligand_col = _choose_ligand_column(reader.fieldnames, preferred_cols)
            if ligand_col is None:
                continue
            seen: set[str] = set()
            local_dup = 0
            for row in reader:
                token = str(row.get(ligand_col, "")).strip()
                if not token:
                    continue
                total_rows += 1
                if token in seen:
                    local_dup += 1
                else:
                    seen.add(token)
            duplicate_rows += local_dup
            if local_dup > 0:
                files_with_duplicates += 1
    return {
        "files": int(len(files)),
        "rows": int(total_rows),
        "duplicate_rows": int(duplicate_rows),
        "files_with_duplicates": int(files_with_duplicates),
    }


def _extract_combo_from_path(path: Path, run_root: Path) -> Optional[str]:
    try:
        rel = path.resolve().relative_to(run_root.resolve())
    except Exception:
        return None
    parts = list(rel.parts)
    if len(parts) < 4:
        return None
    pdb_id, variant, ph = parts[0], parts[1], parts[2]
    if not pdb_id or not variant or not ph:
        return None
    return f"{pdb_id.upper()}|{variant.upper()}|{ph}"


def _stage_counts_from_docking_summary(path: Path) -> Dict[str, int]:
    counts = {"rows": 0, "stage1": 0, "stage2": 0, "stage3": 0}
    if not path.exists():
        return counts
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            counts["rows"] += 1
            for stage in ("stage1", "stage2", "stage3"):
                token = str(row.get(stage, "")).strip()
                if token and token.lower() not in {"nan", "none"}:
                    counts[stage] += 1
    return counts


def _collect_stage_screen_counts(
    docked_run_root: Path,
) -> Dict[str, Dict[str, int]]:
    completion_expected = {"stage1": 0, "stage2": 0, "stage3": 0}
    summary_scored = {"rows": 0, "stage1": 0, "stage2": 0, "stage3": 0}
    combo_completion: Dict[str, Dict[str, int]] = {}
    summary_seen: set[str] = set()

    for completion_path in sorted(docked_run_root.rglob("completion_vina.json")):
        stage = completion_path.parent.name.lower().strip()
        if stage not in {"stage1", "stage2", "stage3"}:
            continue
        payload = _read_json(completion_path)
        expected_n = int(_safe_int(payload.get("expected_count")) or 0)
        completion_expected[stage] += expected_n

        combo_key = _extract_combo_from_path(completion_path, docked_run_root)
        if combo_key:
            slot = combo_completion.setdefault(
                combo_key,
                {"stage1": 0, "stage2": 0, "stage3": 0},
            )
            slot[stage] += expected_n
            if combo_key not in summary_seen:
                summary_seen.add(combo_key)
                parts = combo_key.split("|", 2)
                if len(parts) == 3:
                    summary_path = (
                        docked_run_root
                        / parts[0]
                        / parts[1]
                        / parts[2]
                        / "docking_score_summary.csv"
                    )
                    s_counts = _stage_counts_from_docking_summary(summary_path)
                    for key in ("rows", "stage1", "stage2", "stage3"):
                        summary_scored[key] += int(s_counts.get(key, 0))

    return {
        "completion_expected": completion_expected,
        "docking_summary_scored": summary_scored,
        "combo_count": {"value": int(len(combo_completion))},
    }


def _collect_stage_order(manifest: Mapping[str, Any]) -> Dict[str, Any]:
    proteins = manifest.get("proteins")
    if not isinstance(proteins, Mapping):
        return {"rows": {}, "ok_count": 0, "fail_count": 0, "missing_count": 0}

    rows: Dict[str, Dict[str, Any]] = {}
    ok_count = 0
    fail_count = 0
    missing_count = 0

    for key, payload in proteins.items():
        if not isinstance(payload, Mapping):
            continue
        combo_key = str(key)
        pdb_id = str(payload.get("pdb_id") or "").strip().upper()
        variant = str(payload.get("variant") or "").strip().upper()
        ph = str(payload.get("ph") or "").strip()
        if pdb_id and variant and ph:
            combo_key = f"{pdb_id}|{variant}|{ph}"

        stages = payload.get("stages")
        docking = stages.get("docking") if isinstance(stages, Mapping) else None
        docking_details = docking.get("details") if isinstance(docking, Mapping) else None
        per_stage = (
            docking_details.get("per_stage") if isinstance(docking_details, Mapping) else None
        )

        def _stage_finished(stage_name: str) -> Optional[dt.datetime]:
            if not isinstance(per_stage, Mapping):
                return None
            block = per_stage.get(stage_name)
            if not isinstance(block, Mapping):
                return None
            timing = block.get("timing")
            if not isinstance(timing, Mapping):
                return None
            return _parse_iso_utc(timing.get("finished_at"))

        s1 = _stage_finished("vina_stage1")
        s2 = _stage_finished("vina_stage2")
        s3 = _stage_finished("vina_stage3")

        postprocessing = stages.get("postprocessing") if isinstance(stages, Mapping) else None
        post_timing = postprocessing.get("timing") if isinstance(postprocessing, Mapping) else None
        scorch_start = (
            _parse_iso_utc(post_timing.get("started_at"))
            if isinstance(post_timing, Mapping)
            else None
        )

        missing_fields = [
            name
            for name, value in (
                ("stage1_finished", s1),
                ("stage2_finished", s2),
                ("stage3_finished", s3),
                ("scorch_started", scorch_start),
            )
            if value is None
        ]
        order_ok = (
            not missing_fields
            and s1 is not None
            and s2 is not None
            and s3 is not None
            and scorch_start is not None
            and s1 <= s2 <= s3 <= scorch_start
        )
        if missing_fields:
            missing_count += 1
        elif order_ok:
            ok_count += 1
        else:
            fail_count += 1

        rows[combo_key] = {
            "stage1_finished": _fmt_iso_utc(s1),
            "stage2_finished": _fmt_iso_utc(s2),
            "stage3_finished": _fmt_iso_utc(s3),
            "scorch_started": _fmt_iso_utc(scorch_start),
            "missing_fields": missing_fields,
            "order_ok": bool(order_ok),
        }

    return {
        "rows": rows,
        "ok_count": int(ok_count),
        "fail_count": int(fail_count),
        "missing_count": int(missing_count),
    }


def _load_table_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        return []
    header_idx = 0
    for i, line in enumerate(lines):
        token = line.strip().lower()
        if token.startswith("run_id") or token.startswith("variant"):
            header_idx = i
            break
    data_lines = lines[header_idx:]
    delimiter = "\t" if "\t" in data_lines[0] else ","
    reader = csv.DictReader(data_lines, delimiter=delimiter)
    return [dict(row) for row in reader]


def _load_post_docked_summary_table(path: Path) -> Dict[str, Dict[str, Optional[float]]]:
    if not path.exists():
        return {}
    rows = _load_table_rows(path)
    out: Dict[str, Dict[str, Optional[float]]] = {}
    for row in rows:
        pdb_id = str(row.get("pdb_id", "")).strip().upper()
        variant = str(row.get("variant", "")).strip().upper()
        ph = str(row.get("pH", row.get("ph", ""))).strip()
        if not pdb_id or not variant or not ph:
            continue
        key = f"{pdb_id}|{variant}|{ph}"
        out[key] = {
            "ROC_AUC": _safe_float(row.get("ROC_AUC")),
            "PR_AUC": _safe_float(row.get("PR_AUC")),
            "N": _safe_float(row.get("N")),
            "n_actives": _safe_float(row.get("n_actives")),
        }
    return out


def _find_post_docked_summary_file(root: Path, run_id: str) -> Optional[Path]:
    base = root / "analysis" / "dud_eval" / run_id / "post_docked"
    if not base.exists():
        return None
    candidates = sorted(base.glob("consensus_reranked_scorch_summary*.tsv"))
    return candidates[0] if candidates else None


def _load_reference_table(path: Optional[Path]) -> Dict[str, Dict[str, Optional[float]]]:
    if path is None or not path.exists():
        return {}
    rows = _load_table_rows(path)
    out: Dict[str, Dict[str, Optional[float]]] = {}
    for row in rows:
        pdb_id = str(row.get("pdb_id", "")).strip().upper()
        variant = str(row.get("variant", "")).strip().upper()
        ph = str(row.get("pH", row.get("ph", ""))).strip()
        if not pdb_id or not variant or not ph:
            continue
        key = f"{pdb_id}|{variant}|{ph}"
        out[key] = {
            "ROC_AUC": _safe_float(row.get("ROC_AUC")),
            "PR_AUC": _safe_float(row.get("PR_AUC")),
        }
    return out


def _collect_target_dud_metrics(root: Path, run_id: str) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    summary_path = _find_post_docked_summary_file(root, run_id)
    if summary_path is not None:
        summary_table = _load_post_docked_summary_table(summary_path)
        for key in _BENCH2_TARGET_KEYS:
            row = summary_table.get(key)
            if row is None:
                continue
            out[key] = {
                "source": str(summary_path),
                **row,
            }

    base = root / "analysis" / "dud_eval" / run_id
    for key in _BENCH2_TARGET_KEYS:
        if key in out:
            continue
        pdb, variant, ph = key.split("|", 2)
        for sub in ("", "consensus", "post_docked", "post_docked_scorch"):
            metrics = (
                base / f"{pdb}__{variant}__{ph}" / "metrics.tsv"
                if not sub
                else base / sub / f"{pdb}__{variant}__{ph}" / "metrics.tsv"
            )
            rows = _load_table_rows(metrics)
            if not rows:
                continue
            row_metrics = rows[0]
            out[key] = {
                "source": str(metrics),
                "ROC_AUC": _safe_float(row_metrics.get("ROC_AUC")),
                "PR_AUC": _safe_float(row_metrics.get("PR_AUC")),
                "N": _safe_float(row_metrics.get("N")),
                "n_actives": _safe_float(row_metrics.get("n_actives")),
            }
            break
    return out


def _run_strict_throughput_integrity(root: Path, run_id: str) -> Dict[str, Any]:
    cmd = [
        sys.executable,
        "-m",
        "analysis.cli.throughput_integrity",
        "--run-id",
        str(run_id),
        "--repo-root",
        str(root),
        "--strict",
    ]
    proc = subprocess.run(
        cmd,
        cwd=str(root),
        check=False,
        capture_output=True,
        text=True,
    )
    payload = _read_json(run_output_dir(root, "data", run_id) / "throughput_integrity.json")
    return {
        "returncode": int(proc.returncode),
        "integrity_pass": bool(payload.get("integrity_pass"))
        if payload.get("integrity_pass") is not None
        else None,
    }


def _quality_delta_rows(
    lhs: Mapping[str, Mapping[str, Optional[float]]],
    rhs: Mapping[str, Mapping[str, Optional[float]]],
    rhs_label: str,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for key in sorted(lhs.keys()):
        if key not in rhs:
            continue
        left = lhs[key]
        right = rhs[key]
        l_roc = _safe_float(left.get("ROC_AUC"))
        r_roc = _safe_float(right.get("ROC_AUC"))
        l_pr = _safe_float(left.get("PR_AUC"))
        r_pr = _safe_float(right.get("PR_AUC"))
        if l_roc is None or r_roc is None or l_pr is None or r_pr is None:
            continue
        delta_roc = l_roc - r_roc
        delta_pr = l_pr - r_pr
        pdb_id, variant, ph = key.split("|", 2)
        rows.append(
            {
                "key": key,
                "pdb_id": pdb_id,
                "variant": variant,
                "ph": ph,
                "compare_to": rhs_label,
                "roc_auc": l_roc,
                "pr_auc": l_pr,
                "delta_roc_auc": delta_roc,
                "delta_pr_auc": delta_pr,
                "pass_roc": abs(delta_roc) <= _REFERENCE_DELTA_LIMIT,
                "pass_pr": abs(delta_pr) <= _REFERENCE_DELTA_LIMIT,
                "pass_both": abs(delta_roc) <= _REFERENCE_DELTA_LIMIT
                and abs(delta_pr) <= _REFERENCE_DELTA_LIMIT,
            }
        )
    return rows


def _write_delim(
    rows: Sequence[Mapping[str, Any]],
    path: Path,
    *,
    delimiter: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(row))


def _write_csv(rows: Sequence[Mapping[str, Any]], path: Path) -> None:
    _write_delim(rows, path, delimiter=",")


def _write_tsv(rows: Sequence[Mapping[str, Any]], path: Path) -> None:
    _write_delim(rows, path, delimiter="\t")


def _bool_text(value: bool) -> str:
    return "selected" if value else "not-selected"


def _write_optimization_backlog(
    out_path: Path,
    *,
    snapshot: Mapping[str, Any],
    run_efficiency: Mapping[str, Any],
) -> None:
    chunk_planned = int(snapshot.get("chunk_planned") or 0)
    chunk_completed = int(snapshot.get("chunk_completed") or 0)
    stage_counts_raw = snapshot.get("stage_screen_counts")
    stage_counts = stage_counts_raw if isinstance(stage_counts_raw, Mapping) else {}
    completion_raw = (
        stage_counts.get("completion_expected")
        if isinstance(stage_counts.get("completion_expected"), Mapping)
        else {}
    )
    expected_stage1 = int(
        _safe_int(
            completion_raw.get("stage1") if isinstance(completion_raw, Mapping) else 0
        )
        or 0
    )
    stage_order = snapshot.get("per_combo_stage_order")
    stage_fail = 0
    if isinstance(stage_order, Mapping):
        stage_fail = int(stage_order.get("fail_count") or 0)

    scorch_share = _safe_float(run_efficiency.get("scorch_share_of_combo_wall"))
    avg_ligands_per_chunk = float(expected_stage1) / max(1.0, float(chunk_planned or 1))
    chunk_density_high = bool(chunk_planned >= 100 and avg_ligands_per_chunk <= 48.0)

    items = [
        {
            "item": "Chunk execution tail imbalance",
            "expected_gain": "Medium to high on tails with idle workers",
            "risk": "Low",
            "complexity": "Low",
            "selected": bool(chunk_density_high),
            "reason": (
                f"planned_chunks={chunk_planned} "
                f"avg_ligands_per_chunk={avg_ligands_per_chunk:.1f}"
            ),
        },
        {
            "item": "SCORCH tail occupancy and housekeeping overhead",
            "expected_gain": "High when SCORCH dominates wall time",
            "risk": "Low to medium",
            "complexity": "Low",
            "selected": bool((scorch_share or 0.0) >= 0.35),
            "reason": f"scorch_share={float(scorch_share or 0.0):.3f}",
        },
        {
            "item": "Control-plane overhead from distributed state scans",
            "expected_gain": "Medium for large chunk counts",
            "risk": "Low",
            "complexity": "Medium",
            "selected": bool(chunk_planned >= 100),
            "reason": (
                f"chunks={chunk_planned} completed={chunk_completed} "
                f"avg_ligands_per_chunk={avg_ligands_per_chunk:.1f}"
            ),
        },
        {
            "item": "High-cardinality logging hotspots",
            "expected_gain": "Low to medium",
            "risk": "Low",
            "complexity": "Low",
            "selected": bool(stage_fail == 0),
            "reason": f"stage_order_failures={stage_fail}",
        },
    ]

    lines = [
        "# Optimization Backlog",
        "",
        "| Item | Expected Gain | Regression Risk | Complexity | Decision | Basis |",
        "|---|---|---|---|---|---|",
    ]
    for item in items:
        lines.append(
            "| {item} | {expected_gain} | {risk} | {complexity} | {decision} | {reason} |".format(
                item=item["item"],
                expected_gain=item["expected_gain"],
                risk=item["risk"],
                complexity=item["complexity"],
                decision=_bool_text(bool(item["selected"])),
                reason=item["reason"],
            )
        )
    lines.extend(
        [
            "",
            "Selection rule: selected items are high-ROI, low behavior-risk items that do not change ligand selection/scoring semantics.",
        ]
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@dataclass
class RunSnapshot:
    run_id: str
    status: str
    started_at: Optional[str]
    finished_at: Optional[str]
    elapsed_sec: Optional[float]
    chunk_planned: int
    chunk_completed: int
    chunk_claims: int
    combos_scheduled: Optional[int]
    combos_completed: Optional[int]
    combos_failed: Optional[int]
    proteins_scheduled: Optional[int]
    proteins_completed: Optional[int]
    proteins_failed: Optional[int]
    throughput_integrity_pass: Optional[bool]
    throughput_integrity_strict_pass: Optional[bool]
    throughput_integrity_strict_returncode: Optional[int]
    rescored_rows_all: int
    stage_screen_counts: Dict[str, Dict[str, int]]
    per_combo_stage_order: Dict[str, Any]
    consensus_duplicates: Dict[str, int]
    scorch_duplicates: Dict[str, int]
    post_docked_scorch_summary: Dict[str, Dict[str, Optional[float]]]
    dud_post_docked_target_metrics: Dict[str, Dict[str, Any]]


def _collect_snapshot(root: Path, run_id: str, *, run_strict: bool) -> RunSnapshot:
    manifests_dir = run_output_dir(root, "manifests", run_id)
    manifest_path = manifests_dir / "run_manifest.yaml"
    manifest = _read_yaml(manifest_path)
    summary_raw = manifest.get("summary")
    summary: Mapping[str, Any] = summary_raw if isinstance(summary_raw, Mapping) else {}
    status = str(manifest.get("status") or "")

    timing = manifest.get("timing") if isinstance(manifest.get("timing"), Mapping) else {}
    started = _parse_iso_utc(timing.get("started_at")) if isinstance(timing, Mapping) else None
    finished = _parse_iso_utc(timing.get("finished_at")) if isinstance(timing, Mapping) else None
    elapsed_sec = (
        float((finished - started).total_seconds())
        if started is not None and finished is not None
        else _safe_float(timing.get("wall_time_sec"))
        if isinstance(timing, Mapping)
        else None
    )

    distributed = manifests_dir / "distributed"
    chunk_planned = _count_planned_chunks(distributed)
    chunk_completed = len(list((distributed / "chunk_results").glob("*.json")))
    chunk_claims = len(list((distributed / "chunk_claims").glob("*.json")))

    ti_path = run_output_dir(root, "data", run_id) / "throughput_integrity.json"
    ti_payload = _read_json(ti_path)
    integrity_pass_raw = ti_payload.get("integrity_pass")
    integrity_pass = bool(integrity_pass_raw) if integrity_pass_raw is not None else None

    strict_info: Dict[str, Any] = {}
    if run_strict:
        strict_info = _run_strict_throughput_integrity(root, run_id)

    post_root = run_output_dir(root, "post_docked", run_id)
    rescored_rows_all = 0
    for score_csv in sorted(post_root.rglob("scorch_scores_all.csv")):
        rescored_rows_all += _count_csv_rows(score_csv)

    docked_run_root = run_output_dir(root, "docked", run_id)
    stage_counts = _collect_stage_screen_counts(docked_run_root)
    stage_order = _collect_stage_order(manifest)

    consensus_dup_stats = _duplicate_stats_for_pattern(
        post_root,
        "consensus_reranked_scorch.csv",
        ("ligand", "Ligand", "Ligand_ID", "lig_id", "ligand_file"),
    )
    scorch_dup_stats = _duplicate_stats_for_pattern(
        post_root,
        "scorch_scores_all.csv",
        ("Ligand_ID", "ligand", "Ligand", "lig_id", "ligand_file"),
    )

    summary_file = _find_post_docked_summary_file(root, run_id)
    post_summary = (
        _load_post_docked_summary_table(summary_file) if summary_file is not None else {}
    )

    target_metrics = _collect_target_dud_metrics(root, run_id)

    return RunSnapshot(
        run_id=run_id,
        status=status,
        started_at=_fmt_iso_utc(started),
        finished_at=_fmt_iso_utc(finished),
        elapsed_sec=elapsed_sec,
        chunk_planned=int(chunk_planned),
        chunk_completed=int(chunk_completed),
        chunk_claims=int(chunk_claims),
        combos_scheduled=_safe_int(summary.get("total_combos_scheduled")),
        combos_completed=_safe_int(summary.get("total_combos_completed")),
        combos_failed=_safe_int(summary.get("total_combos_failed")),
        proteins_scheduled=_safe_int(summary.get("total_proteins_scheduled")),
        proteins_completed=_safe_int(summary.get("total_proteins_completed")),
        proteins_failed=_safe_int(summary.get("total_proteins_failed")),
        throughput_integrity_pass=integrity_pass,
        throughput_integrity_strict_pass=(
            bool(strict_info.get("integrity_pass"))
            if strict_info.get("integrity_pass") is not None
            else None
        ),
        throughput_integrity_strict_returncode=(
            _safe_int(strict_info.get("returncode"))
        ),
        rescored_rows_all=int(rescored_rows_all),
        stage_screen_counts=stage_counts,
        per_combo_stage_order=stage_order,
        consensus_duplicates=consensus_dup_stats,
        scorch_duplicates=scorch_dup_stats,
        post_docked_scorch_summary=post_summary,
        dud_post_docked_target_metrics=target_metrics,
    )


def _collect_acceptance(
    *,
    comparison: Mapping[str, Any],
    quality_vs_reference: Sequence[Mapping[str, Any]],
) -> Dict[str, Any]:
    speedup_pct = _safe_float(comparison.get("speedup_pct"))
    perf_ok = speedup_pct is None or speedup_pct >= -2.0

    quality_rows = list(quality_vs_reference)
    quality_ok = True
    for row in quality_rows:
        if not bool(row.get("pass_both")):
            quality_ok = False
            break

    stage_parity = bool(comparison.get("stage_screen_count_parity"))
    rescored_parity = bool(comparison.get("rescored_row_parity"))
    no_dups = (
        int(comparison.get("consensus_duplicate_rows") or 0) == 0
        and int(comparison.get("scorch_duplicate_rows") or 0) == 0
    )

    return {
        "performance_ok": perf_ok,
        "parity_ok": stage_parity and rescored_parity,
        "duplicates_ok": no_dups,
        "quality_vs_reference_ok": quality_ok,
        "overall_ok": perf_ok and stage_parity and rescored_parity and no_dups and quality_ok,
    }


def _capture_run(
    root: Path,
    run_id: str,
    reference: Optional[Path],
    *,
    run_strict: bool,
) -> Dict[str, Any]:
    snap = _collect_snapshot(root, run_id, run_strict=run_strict)
    out_dir = root / "analysis" / "benchmarks" / run_id
    out_dir.mkdir(parents=True, exist_ok=True)

    payload: Dict[str, Any] = {"snapshot": asdict(snap)}
    ref_table = _load_reference_table(reference)
    if ref_table:
        rows = _quality_delta_rows(snap.post_docked_scorch_summary, ref_table, "reference")
        payload["quality_vs_reference"] = rows
        _write_csv(rows, out_dir / "quality_compare_reference.csv")

    run_eff = _read_json(run_output_dir(root, "manifests", run_id) / "run_efficiency.json")
    payload["run_efficiency"] = dict(run_eff)

    _write_optimization_backlog(
        out_dir / "optimization_backlog.md",
        snapshot=asdict(snap),
        run_efficiency=run_eff,
    )

    out_path = out_dir / "baseline_summary.json"
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "path": str(out_path),
        "rows_quality_vs_reference": len(payload.get("quality_vs_reference", [])),
        "backlog_path": str(out_dir / "optimization_backlog.md"),
    }


def _compare_runs(
    root: Path,
    baseline_run_id: str,
    optimized_run_id: str,
    reference: Optional[Path],
    *,
    run_strict: bool,
) -> Dict[str, Any]:
    baseline = _collect_snapshot(root, baseline_run_id, run_strict=run_strict)
    optimized = _collect_snapshot(root, optimized_run_id, run_strict=run_strict)

    elapsed_delta = None
    speedup_pct = None
    if baseline.elapsed_sec and optimized.elapsed_sec:
        elapsed_delta = float(optimized.elapsed_sec - baseline.elapsed_sec)
        if baseline.elapsed_sec > 0:
            speedup_pct = float(
                (baseline.elapsed_sec - optimized.elapsed_sec) * 100.0 / baseline.elapsed_sec
            )

    stage_parity = baseline.stage_screen_counts == optimized.stage_screen_counts
    rescored_parity = baseline.rescored_rows_all == optimized.rescored_rows_all

    quality_vs_baseline = _quality_delta_rows(
        optimized.post_docked_scorch_summary,
        baseline.post_docked_scorch_summary,
        "baseline",
    )
    ref_table = _load_reference_table(reference)
    quality_vs_reference = (
        _quality_delta_rows(optimized.post_docked_scorch_summary, ref_table, "reference")
        if ref_table
        else []
    )

    comparison_block: Dict[str, Any] = {
        "elapsed_delta_sec": elapsed_delta,
        "speedup_pct": speedup_pct,
        "stage_screen_count_parity": stage_parity,
        "rescored_row_parity": rescored_parity,
        "consensus_duplicate_rows": optimized.consensus_duplicates.get("duplicate_rows", 0),
        "consensus_duplicate_files": optimized.consensus_duplicates.get("files_with_duplicates", 0),
        "scorch_duplicate_rows": optimized.scorch_duplicates.get("duplicate_rows", 0),
        "scorch_duplicate_files": optimized.scorch_duplicates.get("files_with_duplicates", 0),
    }

    payload: Dict[str, Any] = {
        "baseline": asdict(baseline),
        "optimized": asdict(optimized),
        "comparison": comparison_block,
        "quality_vs_baseline": quality_vs_baseline,
        "quality_vs_reference": quality_vs_reference,
        "acceptance": _collect_acceptance(
            comparison=comparison_block,
            quality_vs_reference=quality_vs_reference,
        ),
    }

    out_dir = root / "analysis" / "benchmarks" / optimized_run_id
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "comparison.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    _write_tsv(quality_vs_baseline, out_dir / "quality_compare_baseline.tsv")
    if quality_vs_reference:
        _write_csv(quality_vs_reference, out_dir / "quality_compare_reference.csv")
    return {
        "path": str(out_dir / "comparison.json"),
        "quality_vs_baseline_rows": len(quality_vs_baseline),
        "quality_vs_reference_rows": len(quality_vs_reference),
    }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Capture and compare bench2 benchmark outputs and quality metrics.",
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    cap = sub.add_parser("capture", help="Capture baseline summary for one run.")
    cap.add_argument("--run-id", required=True)
    cap.add_argument("--reference", help="Optional CSV/TSV reference file for quality deltas.")
    cap.add_argument(
        "--no-strict-integrity",
        action="store_true",
        help="Skip strict throughput-integrity execution during snapshot capture.",
    )

    cmp = sub.add_parser("compare", help="Compare optimized run against baseline.")
    cmp.add_argument("--baseline-run-id", required=True)
    cmp.add_argument("--optimized-run-id", required=True)
    cmp.add_argument("--reference", help="Optional CSV/TSV reference file for quality deltas.")
    cmp.add_argument(
        "--no-strict-integrity",
        action="store_true",
        help="Skip strict throughput-integrity execution during snapshot capture.",
    )
    return parser


def main() -> int:
    args = _build_parser().parse_args()
    root = _repo_root()
    reference = Path(args.reference).resolve() if getattr(args, "reference", None) else None
    run_strict = not bool(getattr(args, "no_strict_integrity", False))

    if args.cmd == "capture":
        result = _capture_run(
            root,
            str(args.run_id),
            reference,
            run_strict=run_strict,
        )
    elif args.cmd == "compare":
        result = _compare_runs(
            root,
            str(args.baseline_run_id),
            str(args.optimized_run_id),
            reference,
            run_strict=run_strict,
        )
    else:
        raise RuntimeError(f"Unsupported command: {args.cmd}")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
