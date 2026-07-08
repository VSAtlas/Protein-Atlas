from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, List, Optional

from config.output_paths import run_output_dir


_COMPLETED_STATUSES = {"completed", "success", "succeeded"}


def _load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(str(path))
    return json.loads(path.read_text(encoding="utf-8"))


def _is_completed(combo: dict[str, Any]) -> bool:
    status = str(combo.get("status") or "").strip().lower()
    return status in _COMPLETED_STATUSES


def _hard_integrity_check(summary: dict[str, Any]) -> tuple[bool, list[str]]:
    reasons: list[str] = []
    if not bool(summary.get("integrity_pass", False)):
        reasons.append("integrity_pass=false")
    combos = summary.get("combos") or []
    for combo in combos:
        if not isinstance(combo, dict):
            continue
        if not _is_completed(combo):
            continue
        key = "|".join(
            (
                str(combo.get("pdb_id") or ""),
                str(combo.get("variant") or ""),
                str(combo.get("ph") or ""),
            )
        )
        expected = int(combo.get("expected_n") or 0)
        docking = int(combo.get("docking_scored_n") or 0)
        post = int(combo.get("post_scored_n") or 0)
        expected_resolved = bool(combo.get("expected_resolved", False))
        known_failed = int(combo.get("known_failed_n") or 0)
        missing_docking = int(combo.get("missing_docking_unexplained_n") or 0)
        missing_post = int(combo.get("missing_post_unexplained_n") or 0)
        missing_prep = int(combo.get("missing_prep_unexplained_n") or 0)
        if not expected_resolved or expected <= 0:
            reasons.append(f"{key}:expected_unresolved")
            continue
        if known_failed > 0:
            reasons.append(f"{key}:known_failed={known_failed}")
        if missing_prep > 0:
            reasons.append(f"{key}:missing_prep_unexplained={missing_prep}")
        if missing_docking > 0:
            reasons.append(f"{key}:missing_docking_unexplained={missing_docking}")
        if missing_post > 0:
            reasons.append(f"{key}:missing_post_unexplained={missing_post}")
        if docking != expected:
            reasons.append(f"{key}:docking_scored={docking} expected={expected}")
        if post != expected:
            reasons.append(f"{key}:post_scored={post} expected={expected}")
    return len(reasons) == 0, reasons


def _load_util_summary(path: Optional[Path]) -> Optional[dict[str, Any]]:
    if path is None:
        return None
    return _load_json(path)


def _sum_frag_dispatch(summary: dict[str, Any]) -> float:
    attribution = summary.get("attribution") or {}
    frag = float((attribution.get("fragmentation") or {}).get("idle_fraction") or 0.0)
    disp = float((attribution.get("dispatch_gap") or {}).get("idle_fraction") or 0.0)
    return frag + disp


def _build_report(
    *,
    baseline: dict[str, Any],
    candidate: dict[str, Any],
    baseline_util: Optional[dict[str, Any]],
    candidate_util: Optional[dict[str, Any]],
    target_wall_pct: float,
    target_idle_core_pct: float,
    target_idle_frac_abs: float,
    target_frag_disp_pct: float,
) -> dict[str, Any]:
    base_wall = float(baseline.get("wall_seconds") or 0.0)
    cand_wall = float(candidate.get("wall_seconds") or 0.0)
    wall_improve_pct = (
        ((base_wall - cand_wall) / base_wall) * 100.0 if base_wall > 0.0 else 0.0
    )

    hard_ok, hard_reasons = _hard_integrity_check(candidate)

    util_metrics: dict[str, Any] = {"enabled": False}
    util_ok = True
    if baseline_util is not None and candidate_util is not None:
        util_metrics["enabled"] = True
        base_idle_core = float(baseline_util.get("IdleCoreSeconds") or 0.0)
        cand_idle_core = float(candidate_util.get("IdleCoreSeconds") or 0.0)
        base_idle_frac = float(baseline_util.get("IdleFrac") or 0.0)
        cand_idle_frac = float(candidate_util.get("IdleFrac") or 0.0)
        base_frag_disp = _sum_frag_dispatch(baseline_util)
        cand_frag_disp = _sum_frag_dispatch(candidate_util)

        idle_core_improve_pct = (
            ((base_idle_core - cand_idle_core) / base_idle_core) * 100.0
            if base_idle_core > 0.0
            else 0.0
        )
        idle_frac_abs_drop = base_idle_frac - cand_idle_frac
        frag_disp_improve_pct = (
            ((base_frag_disp - cand_frag_disp) / base_frag_disp) * 100.0
            if base_frag_disp > 0.0
            else 0.0
        )
        util_metrics.update(
            {
                "base_idle_core_seconds": base_idle_core,
                "candidate_idle_core_seconds": cand_idle_core,
                "idle_core_improve_pct": idle_core_improve_pct,
                "base_idle_frac": base_idle_frac,
                "candidate_idle_frac": cand_idle_frac,
                "idle_frac_abs_drop": idle_frac_abs_drop,
                "base_frag_plus_dispatch": base_frag_disp,
                "candidate_frag_plus_dispatch": cand_frag_disp,
                "frag_plus_dispatch_improve_pct": frag_disp_improve_pct,
            }
        )
        util_ok = (
            idle_core_improve_pct >= target_idle_core_pct
            and idle_frac_abs_drop >= target_idle_frac_abs
            and frag_disp_improve_pct >= target_frag_disp_pct
        )

    wall_ok = wall_improve_pct >= target_wall_pct
    accepted = bool(hard_ok and wall_ok and util_ok)

    return {
        "accepted": accepted,
        "checks": {
            "hard_integrity_ok": hard_ok,
            "wall_ok": wall_ok,
            "util_ok": util_ok,
        },
        "targets": {
            "wall_improvement_pct": target_wall_pct,
            "idle_core_improvement_pct": target_idle_core_pct,
            "idle_frac_abs_drop": target_idle_frac_abs,
            "frag_plus_dispatch_improvement_pct": target_frag_disp_pct,
        },
        "results": {
            "wall_improvement_pct": wall_improve_pct,
            "baseline_wall_seconds": base_wall,
            "candidate_wall_seconds": cand_wall,
            "util": util_metrics,
        },
        "hard_integrity_failures": hard_reasons,
    }


def run(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Compare throughput runs with hard integrity gates and KPI targets."
    )
    ap.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[2]))
    ap.add_argument("--baseline-run-id", required=True)
    ap.add_argument("--candidate-run-id", required=True)
    ap.add_argument(
        "--util-baseline-summary",
        default="",
        help="Optional util bench baseline summary.json path",
    )
    ap.add_argument(
        "--util-candidate-summary",
        default="",
        help="Optional util bench candidate summary.json path",
    )
    ap.add_argument("--target-wall-improvement-pct", type=float, default=10.0)
    ap.add_argument("--target-idle-core-improvement-pct", type=float, default=15.0)
    ap.add_argument("--target-idle-frac-abs", type=float, default=0.03)
    ap.add_argument("--target-frag-dispatch-improvement-pct", type=float, default=10.0)
    ap.add_argument("--out-json", default="")
    args = ap.parse_args(argv)

    repo_root = Path(args.repo_root).resolve()
    baseline_dir = run_output_dir(repo_root, "data", str(args.baseline_run_id))
    candidate_dir = run_output_dir(repo_root, "data", str(args.candidate_run_id))
    baseline_path = baseline_dir / "throughput_integrity.json"
    candidate_path = candidate_dir / "throughput_integrity.json"
    baseline = _load_json(baseline_path)
    candidate = _load_json(candidate_path)

    util_baseline = _load_util_summary(
        Path(args.util_baseline_summary).resolve()
        if str(args.util_baseline_summary).strip()
        else None
    )
    util_candidate = _load_util_summary(
        Path(args.util_candidate_summary).resolve()
        if str(args.util_candidate_summary).strip()
        else None
    )

    report = _build_report(
        baseline=baseline,
        candidate=candidate,
        baseline_util=util_baseline,
        candidate_util=util_candidate,
        target_wall_pct=float(args.target_wall_improvement_pct),
        target_idle_core_pct=float(args.target_idle_core_improvement_pct),
        target_idle_frac_abs=float(args.target_idle_frac_abs),
        target_frag_disp_pct=float(args.target_frag_dispatch_improvement_pct),
    )

    out_json = str(args.out_json).strip()
    if out_json:
        out_path = Path(out_json).resolve()
    else:
        out_path = (
            repo_root
            / "data"
            / str(args.candidate_run_id)
            / "throughput_acceptance.json"
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=False))
    return 0 if bool(report.get("accepted", False)) else 2


def main(argv: Optional[List[str]] = None) -> int:
    return run(argv)


if __name__ == "__main__":
    raise SystemExit(main())
