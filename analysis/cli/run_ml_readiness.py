from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from analysis.cli.run_ml_training_pass import EXPERT_DEFAULTS, main as training_main
from analysis.ml.leaderboard import write_ml_leaderboard
from analysis.ml.preflight import build_ml_doctor_report, write_doctor_report
from analysis.ml.split_locking import lock_run_splits


def _locked_split_overrides(split_manifest: dict[str, Any]) -> dict[str, Path]:
    overrides: dict[str, Path] = {}
    for row in split_manifest.get("results", []):
        if row.get("status") != "locked":
            continue
        expert = str(row.get("expert") or "")
        default = EXPERT_DEFAULTS.get(expert)
        if not default or row.get("split_mode") != default.get("split"):
            continue
        metadata = row.get("split_metadata")
        manifest = row.get("split_manifest")
        if metadata:
            overrides[expert] = Path(str(metadata)).parent
        elif manifest:
            overrides[expert] = Path(str(manifest)).parent
    return overrides


def _sample_training_status(sample_out: Path, exit_code: int) -> dict[str, Any]:
    status: dict[str, Any] = {"exit_code": int(exit_code or 0), "out_dir": str(sample_out)}
    manifest_path = sample_out / "ml_training_pass_manifest.json"
    if not manifest_path.exists():
        status["status"] = "failed"
        status["error"] = "sample training did not write ml_training_pass_manifest.json"
        return status
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception as exc:
        status["status"] = "failed"
        status["error"] = f"could not read sample training manifest: {exc}"
        return status
    experts = list(manifest.get("experts") or [])
    trained = [row for row in experts if row.get("status") == "trained"]
    failed = [row for row in experts if row.get("status") == "failed"]
    status["trained_experts"] = [row.get("expert") for row in trained]
    status["failed_experts"] = [row.get("expert") for row in failed]
    if exit_code or failed or not trained:
        status["status"] = "failed"
        if not trained:
            status["error"] = "sample training produced no trained experts"
        elif failed:
            status["error"] = "one or more sample training experts failed"
        return status
    status["status"] = "ok"
    return status


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Prepare an imported Atlas run for fast ML iteration.")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--ml-root", required=True, type=Path)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--claim-mode", choices=["exploratory", "publication"], default="exploratory")
    parser.add_argument("--sample-rows", type=int, default=5000)
    parser.add_argument("--sample-seed", type=int, default=42)
    parser.add_argument("--skip-splits", action="store_true")
    parser.add_argument("--skip-sample-train", action="store_true")
    args = parser.parse_args(argv)
    ml_root = args.ml_root.resolve()
    readiness_dir = ml_root / "readiness"
    readiness_dir.mkdir(parents=True, exist_ok=True)
    steps: list[dict[str, Any]] = []
    doctor = build_ml_doctor_report(
        run_id=args.run_id,
        run_dir=args.run_dir,
        claim_mode=args.claim_mode,
        strict=args.claim_mode == "publication",
    )
    doctor_path = write_doctor_report(doctor, readiness_dir / "ml_doctor_report.json")
    steps.append({"name": "doctor", "status": doctor.get("status"), "path": str(doctor_path)})
    locked_split_overrides: dict[str, Path] = {}
    if not args.skip_splits:
        try:
            split_manifest = lock_run_splits(run_id=args.run_id, run_dir=args.run_dir, out_dir=ml_root / "splits", seed=args.sample_seed)
            locked_split_overrides = _locked_split_overrides(split_manifest)
            steps.append({"name": "split_lock", "status": "written", "path": str(ml_root / "splits" / "split_lock_manifest.json"), "locked": len([r for r in split_manifest.get("results", []) if r.get("status") == "locked"])})
            if locked_split_overrides:
                steps.append(
                    {
                        "name": "locked_split_training_overrides",
                        "status": "wired",
                        "experts": {key: str(value) for key, value in locked_split_overrides.items()},
                    }
                )
        except Exception as exc:
            steps.append({"name": "split_lock", "status": "failed", "error": str(exc)})
    if not args.skip_sample_train and doctor.get("status") != "blocked":
        sample_out = ml_root / "sample_training_pass"
        train_args = [
            "--repo-root", str(args.repo_root),
            "--run-id", args.run_id,
            "--run-dir", str(args.run_dir),
            "--out-dir", str(sample_out),
            "--claim-mode", args.claim_mode,
            "--sample-rows", str(args.sample_rows),
            "--sample-seed", str(args.sample_seed),
            "--bootstraps", "0",
        ]
        if locked_split_overrides:
            train_args.extend(
                [
                    "--expert-split-manifest",
                    *[f"{expert}={path}" for expert, path in sorted(locked_split_overrides.items())],
                ]
            )
        try:
            code = training_main(train_args)
            sample_status = _sample_training_status(sample_out, int(code or 0))
            steps.append({"name": "sample_train", **sample_status})
        except Exception as exc:
            steps.append({"name": "sample_train", "status": "failed", "error": str(exc), "out_dir": str(sample_out)})
    elif doctor.get("status") == "blocked":
        steps.append({"name": "sample_train", "status": "skipped", "reason": "doctor_blocked"})
    leaderboard = write_ml_leaderboard(ml_root, ml_root / "leaderboard")
    steps.append({"name": "leaderboard", "status": leaderboard.get("status"), "path": leaderboard.get("csv"), "rows": leaderboard.get("rows")})
    sample_failed = any(step.get("name") == "sample_train" and step.get("status") == "failed" for step in steps)
    status = (
        "blocked"
        if doctor.get("status") == "blocked" or sample_failed
        else "ready_with_warnings"
        if doctor.get("status") == "ready_with_warnings"
        else "ready"
    )
    manifest = {"run_id": args.run_id, "run_dir": str(args.run_dir), "ml_root": str(ml_root), "status": status, "steps": steps}
    (readiness_dir / "ml_readiness_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str), encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return 1 if status == "blocked" else 0


if __name__ == "__main__":
    raise SystemExit(main())
