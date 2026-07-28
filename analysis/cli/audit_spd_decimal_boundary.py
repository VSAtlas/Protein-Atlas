from __future__ import annotations

import argparse
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path
import platform
import shlex
import subprocess
import sys
from typing import Any

import pandas as pd

from analysis.external.spd import (
    SPD_LABEL_POLICY_VERSION,
    _spd_label_status,
)
from analysis.ml.spd_four_expert_tables import BINDING_LABEL_POLICY_VERSION
from analysis.spd_activity_policy import (
    SPD_ACTIVITY_PARSER_VERSION,
    label_spd_binding_observation,
    label_spd_exposure_observation,
    parse_spd_activity_interval,
)


EXPECTED_GIT_SHA = "3554269ce088b709344f6a22dd9a9cd8958c6105"
EXPECTED_INPUT_HASHES = {
    "phase_b_raw": "e3753933e65698383a8f9560c849c38c4ee1129cb66b3292fc379b18c58e7a74",
    "phase_b_pairs": "3d9a871506dff5356079fb296fe4bc1579f8ccb567def3407adcd20845f4001d",
    "phase_b_manifest": "2f8ea7a729101e92eb6caf37bb05236db6b2b3473a14b3769c88a566394edb7c",
    "mapping_manifest": "cb5c1096a57542dedf6fe9ef775839e5da7fbb8a3e5dc6941ef85c579f3489d6",
    "spd_workbook": "6b7681c4cc51670740ca8a20a6317927e47fe3672ae768b5ac4ab2512f06c640",
}
EXPECTED_COUNTS = {
    "raw_binding": {"positive": 3_008, "negative": 112_800, "unknown": 5_289},
    "pair_binding": {"positive": 2_484, "negative": 88_345, "unknown": 4_683},
    "raw_exposure": {"positive": 1_680, "negative": 50_095, "unknown": 69_322},
    "pair_exposure": {"positive": 1_397, "negative": 37_684, "unknown": 56_431},
}
COMPARISON_COLUMNS = [
    "comparison_level",
    "source_row_id",
    "drug_identity",
    "target_identity",
    "raw_relation",
    "raw_activity_value",
    "activity_unit",
    "free_cmax_um",
    "phase_b_binding_label",
    "prompt2c_binding_label",
    "binding_prompt2c_changed",
    "phase_b_exposure_label",
    "prompt2c_exposure_label",
    "exposure_prompt2c_changed",
    "phase_b_exposure_status",
    "prompt2c_exposure_status",
    "exposure_status_prompt2c_changed",
    "legacy_binding_label",
    "legacy_exposure_label",
    "legacy_to_prompt2c_binding_changed",
    "legacy_to_prompt2c_exposure_changed",
]
BOUNDARY_COLUMNS = [
    "case_id",
    "activity_value",
    "activity_relation",
    "activity_unit",
    "free_cmax_um",
    "expected_binding_label",
    "actual_binding_label",
    "expected_exposure_label",
    "actual_exposure_label",
    "expected_margin_interval",
    "actual_margin_lower",
    "actual_margin_upper",
    "passed",
]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git(repo_root: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=repo_root, text=True
    ).strip()


def _label_name(value: Any) -> str:
    if value is None:
        return "unknown"
    try:
        if bool(pd.isna(value)):
            return "unknown"
    except (TypeError, ValueError):
        pass
    text = str(value).strip().lower()
    if text in {"1", "1.0", "positive", "true"}:
        return "positive"
    if text in {"0", "0.0", "negative", "false"}:
        return "negative"
    return "unknown"


def _status_label(status: str) -> str:
    if status == "labeled_relevant":
        return "positive"
    if status.startswith("labeled_"):
        return "negative"
    return "unknown"


def _is_pass_status(value: Any) -> bool:
    return str(value or "").lower() in {"pass", "passed"}


def _decimal_or_nan(value: Any, factor: str = "1") -> Decimal | float:
    try:
        if bool(pd.isna(value)):
            return float("nan")
    except (TypeError, ValueError):
        pass
    try:
        return Decimal(str(value)) * Decimal(factor)
    except (InvalidOperation, TypeError, ValueError):
        return float("nan")


def _canonical_labels(
    activity: Any,
    relation: Any,
    unit: str,
    free_cmax_um: Any,
) -> tuple[str, str]:
    interval = parse_spd_activity_interval(activity, relation, unit)
    binding = label_spd_binding_observation(interval)
    exposure = label_spd_exposure_observation(interval, free_cmax_um)
    return _label_name(binding.numeric_label), _label_name(exposure.numeric_label)


def _boundary_cases() -> pd.DataFrame:
    specifications = [
        ("exposure_07_open", "0.7", ">", "uM", "0.07", "unknown", "negative", "(10,+inf)"),
        ("exposure_07_closed", "0.7", ">=", "uM", "0.07", "unknown", "unknown", "[10,+inf)"),
        ("exposure_07_exact", "0.7", "=", "uM", "0.07", "positive", "positive", "[10,10]"),
        ("exposure_03_open", "0.3", ">", "uM", "0.03", "unknown", "negative", "(10,+inf)"),
        ("exposure_03_closed", "0.3", ">=", "uM", "0.03", "unknown", "unknown", "[10,+inf)"),
        ("exposure_03_exact", "0.3", "=", "uM", "0.03", "positive", "positive", "[10,10]"),
        ("exposure_nm_open", "700", ">", "nM", "0.07", "unknown", "negative", "(10,+inf)"),
        ("exposure_nm_closed", "700", ">=", "nM", "0.07", "unknown", "unknown", "[10,+inf)"),
        ("exposure_nm_exact", "700", "=", "nM", "0.07", "positive", "positive", "[10,10]"),
        ("binding_one_um", "1", "=", "uM", None, "positive", "unknown", "not_tested"),
        ("binding_one_um_nm", "1000", "=", "nM", None, "positive", "unknown", "not_tested"),
        ("binding_ten_um", "10", "=", "uM", None, "negative", "unknown", "not_tested"),
        ("binding_ten_um_nm", "10000", "=", "nM", None, "negative", "unknown", "not_tested"),
    ]
    rows: list[dict[str, Any]] = []
    for (
        case_id,
        activity,
        relation,
        unit,
        free_cmax,
        expected_binding,
        expected_exposure,
        expected_interval,
    ) in specifications:
        interval = parse_spd_activity_interval(activity, relation, unit)
        binding = label_spd_binding_observation(interval)
        if free_cmax is None:
            exposure_label = "unknown"
            margin_lower = None
            margin_upper = None
        else:
            exposure = label_spd_exposure_observation(interval, free_cmax)
            exposure_label = _label_name(exposure.numeric_label)
            margin_lower = exposure.exposure_margin_lower
            margin_upper = exposure.exposure_margin_upper
        actual_binding = _label_name(binding.numeric_label)
        rows.append(
            {
                "case_id": case_id,
                "activity_value": activity,
                "activity_relation": relation,
                "activity_unit": unit,
                "free_cmax_um": free_cmax,
                "expected_binding_label": expected_binding,
                "actual_binding_label": actual_binding,
                "expected_exposure_label": expected_exposure,
                "actual_exposure_label": exposure_label,
                "expected_margin_interval": expected_interval,
                "actual_margin_lower": margin_lower,
                "actual_margin_upper": margin_upper,
                "passed": (
                    actual_binding == expected_binding
                    and exposure_label == expected_exposure
                ),
            }
        )
    return pd.DataFrame(rows, columns=BOUNDARY_COLUMNS)


def _raw_comparison(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row in frame.itertuples(index=False):
        binding, exposure = _canonical_labels(
            row.raw_activity_value,
            row.raw_relation,
            "uM",
            row.free_cmax_um,
        )
        status = _spd_label_status(
            {
                "ac50_nM": _decimal_or_nan(row.raw_activity_value, "1000"),
                "free_cmax_nM": _decimal_or_nan(row.free_cmax_um, "1000"),
                "activity_relation": row.raw_relation,
            },
            10.0,
            100.0,
        )
        production_exposure = _status_label(status)
        if production_exposure != exposure:
            raise RuntimeError(
                f"raw canonical/production mismatch at source row {row.source_row_id}"
            )
        phase_b_binding = _label_name(row.production_binding_label)
        phase_b_exposure = _label_name(row.production_exposure_label)
        legacy_binding = _label_name(row.legacy_binding_label)
        legacy_exposure = _label_name(row.legacy_exposure_label)
        rows.append(
            {
                "comparison_level": "raw_source_observation",
                "source_row_id": row.source_row_id,
                "drug_identity": row.drug_identity,
                "target_identity": row.target_identity,
                "raw_relation": row.raw_relation,
                "raw_activity_value": row.raw_activity_value,
                "activity_unit": "uM",
                "free_cmax_um": row.free_cmax_um,
                "phase_b_binding_label": phase_b_binding,
                "prompt2c_binding_label": binding,
                "binding_prompt2c_changed": phase_b_binding != binding,
                "phase_b_exposure_label": phase_b_exposure,
                "prompt2c_exposure_label": exposure,
                "exposure_prompt2c_changed": phase_b_exposure != exposure,
                "phase_b_exposure_status": row.production_exposure_status,
                "prompt2c_exposure_status": status,
                "exposure_status_prompt2c_changed": (
                    str(row.production_exposure_status) != status
                ),
                "legacy_binding_label": legacy_binding,
                "legacy_exposure_label": legacy_exposure,
                "legacy_to_prompt2c_binding_changed": legacy_binding != binding,
                "legacy_to_prompt2c_exposure_changed": legacy_exposure != exposure,
            }
        )
    return pd.DataFrame(rows, columns=COMPARISON_COLUMNS)


def _pair_comparison(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row in frame.itertuples(index=False):
        free_cmax_nm = _decimal_or_nan(row.free_cmax_nM)
        free_cmax_um = (
            free_cmax_nm / Decimal("1000")
            if isinstance(free_cmax_nm, Decimal)
            else float("nan")
        )
        binding, exposure = _canonical_labels(
            row.ac50_nM,
            row.spd_activity_relation,
            "nM",
            free_cmax_um,
        )
        status = _spd_label_status(
            {
                "ac50_nM": row.ac50_nM,
                "free_cmax_nM": row.free_cmax_nM,
                "activity_relation": row.spd_activity_relation,
            },
            10.0,
            100.0,
        )
        production_exposure = _status_label(status)
        if production_exposure != exposure:
            raise RuntimeError(
                f"pair canonical/production mismatch at {row.drug_id}/{row.target_id}"
            )
        phase_b_binding = _label_name(row.production_binding_label)
        phase_b_exposure = _label_name(row.production_exposure_label)
        legacy_binding = _label_name(row.legacy_binding_label)
        legacy_exposure = _label_name(row.legacy_exposure_label)
        rows.append(
            {
                "comparison_level": "selected_drug_target_pair",
                "source_row_id": row.source_row_id,
                "drug_identity": row.drug_id,
                "target_identity": row.target_id,
                "raw_relation": row.spd_activity_relation,
                "raw_activity_value": row.ac50_nM,
                "activity_unit": "nM",
                "free_cmax_um": free_cmax_um,
                "phase_b_binding_label": phase_b_binding,
                "prompt2c_binding_label": binding,
                "binding_prompt2c_changed": phase_b_binding != binding,
                "phase_b_exposure_label": phase_b_exposure,
                "prompt2c_exposure_label": exposure,
                "exposure_prompt2c_changed": phase_b_exposure != exposure,
                "phase_b_exposure_status": row.spd_label_status,
                "prompt2c_exposure_status": status,
                "exposure_status_prompt2c_changed": str(row.spd_label_status) != status,
                "legacy_binding_label": legacy_binding,
                "legacy_exposure_label": legacy_exposure,
                "legacy_to_prompt2c_binding_changed": legacy_binding != binding,
                "legacy_to_prompt2c_exposure_changed": legacy_exposure != exposure,
            }
        )
    return pd.DataFrame(rows, columns=COMPARISON_COLUMNS)


def _counts(series: pd.Series) -> dict[str, int]:
    values = series.map(_label_name).value_counts()
    return {
        "positive": int(values.get("positive", 0)),
        "negative": int(values.get("negative", 0)),
        "unknown": int(values.get("unknown", 0)),
    }


def run_audit(args: argparse.Namespace) -> dict[str, Any]:
    started = datetime.now(timezone.utc)
    args.out_dir.mkdir(parents=True, exist_ok=False)
    inputs = {
        "phase_b_raw": args.phase_b_raw,
        "phase_b_pairs": args.phase_b_pairs,
        "phase_b_manifest": args.phase_b_manifest,
        "mapping_manifest": args.mapping_manifest,
        "spd_workbook": args.spd_workbook,
    }
    input_hashes = {name: _sha256(path) for name, path in inputs.items()}
    failures = [
        f"{name}_hash={digest}"
        for name, digest in input_hashes.items()
        if digest != EXPECTED_INPUT_HASHES[name]
    ]
    phase_b_manifest = json.loads(args.phase_b_manifest.read_text())
    mapping_manifest = json.loads(args.mapping_manifest.read_text())
    raw = pd.read_csv(args.phase_b_raw, low_memory=False)
    pairs = pd.read_csv(args.phase_b_pairs, low_memory=False)
    raw_comparison = _raw_comparison(raw)
    pair_comparison = _pair_comparison(pairs)
    comparison = pd.concat([raw_comparison, pair_comparison], ignore_index=True)
    cases = _boundary_cases()

    counts = {
        "raw_binding": _counts(raw_comparison["prompt2c_binding_label"]),
        "pair_binding": _counts(pair_comparison["prompt2c_binding_label"]),
        "raw_exposure": _counts(raw_comparison["prompt2c_exposure_label"]),
        "pair_exposure": _counts(pair_comparison["prompt2c_exposure_label"]),
    }
    changed = {
        "raw_binding": int(raw_comparison["binding_prompt2c_changed"].sum()),
        "raw_exposure": int(raw_comparison["exposure_prompt2c_changed"].sum()),
        "raw_exposure_status": int(
            raw_comparison["exposure_status_prompt2c_changed"].sum()
        ),
        "pair_binding": int(pair_comparison["binding_prompt2c_changed"].sum()),
        "pair_exposure": int(pair_comparison["exposure_prompt2c_changed"].sum()),
        "pair_exposure_status": int(
            pair_comparison["exposure_status_prompt2c_changed"].sum()
        ),
    }
    transitions = {
        "raw_exposure_unknown_to_negative": int(
            raw_comparison["legacy_to_prompt2c_exposure_changed"].sum()
        ),
        "pair_exposure_unknown_to_negative": int(
            pair_comparison["legacy_to_prompt2c_exposure_changed"].sum()
        ),
    }
    if len(raw_comparison) != 121_097 or len(pair_comparison) != 95_512:
        failures.append(
            f"row_counts=raw:{len(raw_comparison)},pairs:{len(pair_comparison)}"
        )
    if counts != EXPECTED_COUNTS:
        failures.append(f"label_counts={counts}")
    if any(changed.values()):
        failures.append(f"prompt2c_changes={changed}")
    if transitions != {
        "raw_exposure_unknown_to_negative": 77,
        "pair_exposure_unknown_to_negative": 46,
    }:
        failures.append(f"phase_b_transitions={transitions}")
    if not bool(cases["passed"].all()):
        failures.append("boundary_cases_failed")
    if not _is_pass_status(phase_b_manifest.get("status")):
        failures.append("phase_b_manifest_not_passed")
    if not _is_pass_status(mapping_manifest.get("status")):
        failures.append("mapping_manifest_not_passed")
    mapping_checks = mapping_manifest.get("checks", {})
    for check in (
        "strict_model_ready_excludes_6huj_8hcq",
        "strict_model_ready_retains_eligible_9i52",
        "legacy_no_row_multiplication",
        "strict_no_row_multiplication",
    ):
        if mapping_checks.get(check) is not True:
            failures.append(f"mapping_check_failed:{check}")
    if SPD_ACTIVITY_PARSER_VERSION != "spd_activity_interval_v3":
        failures.append(f"parser_version={SPD_ACTIVITY_PARSER_VERSION}")
    if BINDING_LABEL_POLICY_VERSION != "spd_binding_interval_aware_v2":
        failures.append(f"binding_policy={BINDING_LABEL_POLICY_VERSION}")
    if SPD_LABEL_POLICY_VERSION != "spd_interval_censor_aware_v3":
        failures.append(f"exposure_policy={SPD_LABEL_POLICY_VERSION}")

    cases_path = args.out_dir / "decimal_boundary_cases.csv"
    comparison_path = args.out_dir / "decimal_current_data_comparison.csv"
    summary_path = args.out_dir / "decimal_boundary_summary.md"
    manifest_path = args.out_dir / "decimal_boundary_manifest.json"
    cases.to_csv(cases_path, index=False)
    comparison.to_csv(comparison_path, index=False)
    summary_path.write_text(
        "\n".join(
            [
                "# SPD exact-decimal boundary audit",
                "",
                f"- Gate: **{'PASS' if not failures else 'FAIL'}**",
                f"- Parser implementation: `{SPD_ACTIVITY_PARSER_VERSION}`",
                f"- Binding policy: `{BINDING_LABEL_POLICY_VERSION}` (unchanged)",
                f"- Exposure policy: `{SPD_LABEL_POLICY_VERSION}` (unchanged)",
                f"- Boundary cases: {int(cases['passed'].sum())}/{len(cases)} passed",
                f"- Raw rows: {len(raw_comparison):,}",
                f"- Selected drug-target pairs: {len(pair_comparison):,}",
                f"- Prompt 2C label/status changes: `{changed}`",
                f"- Preserved legacy-to-current exposure correction: `{transitions}`",
                f"- Current counts: `{counts}`",
                (
                    "- Mapping invariants: 9I52 retained as strict-eligible "
                    "DRD1/P21728; 6HUJ and 8HCQ remain strict-ineligible; "
                    "mapping manifest reports no row multiplication."
                ),
                "- Frozen Phase 1 and external four-state tables were not regenerated.",
                f"- Failures: `{failures}`",
                "",
            ]
        ),
        encoding="utf-8",
    )
    artifacts = {
        path.name: _sha256(path)
        for path in (cases_path, comparison_path, summary_path)
    }
    source_paths = {
        "activity_policy": args.repo_root / "analysis/spd_activity_policy.py",
        "production_adapter": args.repo_root / "analysis/external/spd.py",
        "label_enrichment_adapter": (
            args.repo_root / "analysis/ml/spd_label_enrichment.py"
        ),
        "decimal_audit": Path(__file__).resolve(),
    }
    manifest = {
        "status": "passed" if not failures else "failed",
        "failures": failures,
        "git_sha": _git(args.repo_root, "rev-parse", "HEAD"),
        "branch": _git(args.repo_root, "branch", "--show-current"),
        "parser_version": SPD_ACTIVITY_PARSER_VERSION,
        "policy_versions": {
            "binding": BINDING_LABEL_POLICY_VERSION,
            "exposure": SPD_LABEL_POLICY_VERSION,
        },
        "input_paths": {name: str(path.resolve()) for name, path in inputs.items()},
        "input_hashes": input_hashes,
        "source_paths": {
            name: str(path.resolve()) for name, path in source_paths.items()
        },
        "source_hashes": {
            name: _sha256(path) for name, path in source_paths.items()
        },
        "row_counts": {
            "raw": len(raw_comparison),
            "selected_pairs": len(pair_comparison),
        },
        "before_prompt2c_label_counts": EXPECTED_COUNTS,
        "after_prompt2c_label_counts": counts,
        "prompt2c_changed_counts": changed,
        "preserved_phase_b_transitions": transitions,
        "boundary_case_ids": cases["case_id"].tolist(),
        "mapping_checks": {
            key: mapping_checks.get(key)
            for key in (
                "strict_model_ready_excludes_6huj_8hcq",
                "strict_model_ready_retains_eligible_9i52",
                "legacy_no_row_multiplication",
                "strict_no_row_multiplication",
            )
        },
        "commands": [
            shlex.join(
                [
                    sys.executable,
                    "-m",
                    "analysis.cli.audit_spd_decimal_boundary",
                    *sys.argv[1:],
                ]
            ),
            *args.command_record,
        ],
        "tests": args.test_record,
        "environment": {
            "python": platform.python_version(),
            "implementation": platform.python_implementation(),
            "executable": sys.executable,
            "pandas": pd.__version__,
        },
        "artifact_hashes": artifacts,
        "timestamp_started_utc": started.isoformat(),
        "timestamp_completed_utc": datetime.now(timezone.utc).isoformat(),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    if failures:
        raise RuntimeError("; ".join(failures))
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Audit exact-decimal SPD threshold behavior against Phase B."
    )
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--spd-workbook", type=Path, required=True)
    parser.add_argument("--phase-b-raw", type=Path, required=True)
    parser.add_argument("--phase-b-pairs", type=Path, required=True)
    parser.add_argument("--phase-b-manifest", type=Path, required=True)
    parser.add_argument("--mapping-manifest", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--command-record", action="append", default=[])
    parser.add_argument("--test-record", action="append", default=[])
    return parser


def main(argv: list[str] | None = None) -> int:
    manifest = run_audit(build_parser().parse_args(argv))
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
