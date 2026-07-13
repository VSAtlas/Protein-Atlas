"""One-command, provenance-checked SPD add-on merge pipeline."""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.cli.merge_spd_addon_run import merge_spd_addon_tables
from analysis.ml.audit_suite import run_ml_audit_suite
from analysis.ml.feature_metadata import refresh_tables
from analysis.ml.labels import (
    binary_label_series,
    training_eligibility_mask,
    truthy_series,
)
from analysis.ml.spd_identity_reconciliation import (
    default_fda_mapping_path,
    reconcile_spd_phase1_identity,
)


SCHEMA = "atlas.spd-addon-merge-pipeline.v1"
DEFAULT_AUDITS = {
    "spd_binding_label": "spd_binding_nonleaky",
    "spd_exposure_label": "spd_exposure_nonleaky",
    "combined_activity_ml_label": "spd_binding_nonleaky",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected a JSON object: {path}")
    return payload


def _resolve_recorded_path(value: object, *, manifest: Path) -> Path:
    text = str(value or "").strip()
    if not text:
        raise ValueError(f"missing recorded path in {manifest}")
    path = Path(text).expanduser()
    if path.is_absolute():
        return path.resolve()
    repo_candidate = (Path.cwd() / path).resolve()
    if repo_candidate.exists():
        return repo_candidate
    return (manifest.parent / path).resolve()


def _reference_manifest(addon_dir: Path, explicit: Path | None) -> Path:
    if explicit is not None:
        path = explicit.resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        return path
    preferred = addon_dir / "reference_vina_full" / "reference_vina_manifest.json"
    if preferred.is_file():
        return preferred.resolve()
    candidates = sorted(
        path.resolve()
        for path in addon_dir.glob("reference_vina*/reference_vina_manifest.json")
        if "smoke" not in path.as_posix().casefold()
    )
    if len(candidates) != 1:
        raise ValueError(
            "expected exactly one non-smoke reference Vina manifest under "
            f"{addon_dir}; found {len(candidates)}"
        )
    return candidates[0]


def _recorded_output(payload: dict[str, Any], key: str) -> object:
    value = (payload.get("outputs") or {}).get(key)
    if isinstance(value, dict):
        return value.get("path")
    return value


def _resolve_inputs(
    addon_dir: Path,
    *,
    reference_manifest: Path | None,
    selected_pairs: Path | None,
    score_table: Path | None,
) -> tuple[Path, Path, Path, dict[str, Any]]:
    manifest_path = _reference_manifest(addon_dir, reference_manifest)
    manifest = _load_json(manifest_path)
    selected = (
        selected_pairs.resolve()
        if selected_pairs is not None
        else _resolve_recorded_path(manifest.get("selected_pairs"), manifest=manifest_path)
    )
    scores = (
        score_table.resolve()
        if score_table is not None
        else _resolve_recorded_path(
            _recorded_output(manifest, "scores"), manifest=manifest_path
        )
    )
    for path, label in ((selected, "selected pairs"), (scores, "score table")):
        if not path.is_file():
            raise FileNotFoundError(f"{label} does not exist: {path}")
    expected_selected_sha = str(manifest.get("selected_pairs_sha256") or "").strip()
    if expected_selected_sha and _sha256(selected) != expected_selected_sha:
        raise ValueError(
            "selected-pair checksum differs from the reference-score manifest"
        )
    return manifest_path, selected, scores, manifest


def discover_latest_scored_addon(root: Path) -> tuple[Path, dict[str, Any]]:
    """Select the newest timestamped add-on with a complete score contract."""

    candidates = sorted(
        (
            path.resolve()
            for path in root.glob("target_positive_addon_*")
            if path.is_dir()
        ),
        key=lambda path: (path.name, path.stat().st_mtime_ns),
        reverse=True,
    )
    rejected: list[dict[str, str]] = []
    for candidate in candidates:
        try:
            manifest, selected, scores, _ = _resolve_inputs(
                candidate,
                reference_manifest=None,
                selected_pairs=None,
                score_table=None,
            )
        except (FileNotFoundError, KeyError, TypeError, ValueError) as exc:
            rejected.append(
                {
                    "path": str(candidate),
                    "reason": f"{type(exc).__name__}: {exc}",
                }
            )
            continue
        return candidate, {
            "mode": "latest_completed_score_ready",
            "root": str(root.resolve()),
            "selected": str(candidate),
            "reference_manifest": str(manifest),
            "selected_pairs": str(selected),
            "score_table": str(scores),
            "rejected_newer_candidates": rejected,
        }
    raise FileNotFoundError(
        f"no completed score-ready target_positive_addon_* directory under {root}"
    )


def _pair_key(frame: pd.DataFrame) -> pd.Series:
    required = {"pdb_id", "ligand_base"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"feature-ready table lacks pair keys: {sorted(missing)}")
    return (
        frame["pdb_id"].fillna("").astype(str).str.strip().str.upper()
        + "|"
        + frame["ligand_base"].fillna("").astype(str).str.strip().str.lower()
    )


def _validate_merge(
    feature_table: Path,
    addon_rows: Path,
    identity_summary: dict[str, Any],
    mapping: Path,
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    frame = pd.read_csv(feature_table, low_memory=False)
    addon = pd.read_csv(addon_rows, low_memory=False)
    mapping_record = identity_summary.get("mapping") or {}
    eligibility, eligibility_summary = training_eligibility_mask(frame)
    checks: list[dict[str, Any]] = []

    def add(name: str, passed: bool, detail: object) -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    duplicate_pairs = int(_pair_key(frame).duplicated(keep=False).sum())
    add("no_duplicate_pdb_ligand_pairs", duplicate_pairs == 0, duplicate_pairs)
    add(
        "canonical_mapping_path",
        Path(str(mapping_record.get("path") or "")).resolve() == mapping.resolve(),
        mapping_record.get("path"),
    )
    live_mapping_sha = _sha256(mapping)
    add(
        "canonical_mapping_sha256",
        str(mapping_record.get("sha256") or "") == live_mapping_sha,
        live_mapping_sha,
    )
    spd_columns = [
        column
        for column in (
            "spd_binding_label",
            "spd_exposure_label",
            "spd_ac50_uM",
            "spd_exposure_margin",
        )
        if column in addon.columns
    ]
    addon_spd_values = int(
        sum(addon[column].notna().sum() for column in spd_columns)
    )
    add(
        "external_addons_do_not_create_spd_truth",
        addon_spd_values == 0,
        addon_spd_values,
    )
    add(
        "training_eligibility_nonempty",
        bool(eligibility.any()),
        eligibility_summary,
    )
    identity_blocked = (
        ~truthy_series(frame["identity_training_allowed"])
        if "identity_training_allowed" in frame.columns
        else pd.Series(False, index=frame.index)
    )
    add(
        "identity_blockers_quarantined",
        not eligibility.loc[identity_blocked].any(),
        int(identity_blocked.sum()),
    )
    return checks, eligibility_summary


def run_spd_addon_merge_pipeline(
    *,
    base_table: Path,
    addon_dir: Path,
    out_dir: Path,
    mapping: Path | None = None,
    reference_manifest: Path | None = None,
    selected_pairs: Path | None = None,
    score_table: Path | None = None,
    addon_name: str = "external_activity_addon",
    run_audits: bool = True,
    audit_labels: list[str] | None = None,
    addon_selection: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Merge, identity-reconcile, enrich, validate, and audit one add-on pass."""

    base = base_table.resolve()
    addon = addon_dir.resolve()
    destination = out_dir.resolve()
    canonical_mapping = (mapping or default_fda_mapping_path()).resolve()
    for path, label in (
        (base, "base table"),
        (addon, "add-on directory"),
        (canonical_mapping, "canonical FDA mapping"),
    ):
        if not path.exists():
            raise FileNotFoundError(f"{label} does not exist: {path}")
    manifest_path, selected, scores, reference = _resolve_inputs(
        addon,
        reference_manifest=reference_manifest,
        selected_pairs=selected_pairs,
        score_table=score_table,
    )
    destination.mkdir(parents=True, exist_ok=True)
    raw = destination / "merged_raw.csv"
    addon_rows = destination / "merged_addon_rows.csv"
    merged_manifest = merge_spd_addon_tables(
        base_table=base,
        master_rows=scores,
        selected_pairs=selected,
        out=raw,
        addon_name=addon_name,
        addon_rows_out=addon_rows,
        strict_selected_coverage=True,
        backfill_existing_scores=True,
    )
    reconciled = destination / "identity_reconciled.csv"
    identity_summary = reconcile_spd_phase1_identity(
        raw,
        reconciled,
        mapping_csv=canonical_mapping,
        fail_closed=False,
    )
    feature_dir = destination / "feature_ready"
    refresh = refresh_tables(
        [reconciled],
        out_dir=feature_dir,
        run_dir=None,
        chemical_cluster="auto",
        target_family="auto",
        source_lineage="auto",
    )
    feature_table = Path(str(refresh["outputs"][0])).resolve()
    checks, eligibility = _validate_merge(
        feature_table, addon_rows, identity_summary, canonical_mapping
    )
    pd.DataFrame(checks).to_csv(destination / "merge_validation.csv", index=False)

    audits: dict[str, Any] = {}
    audit_progress_path = destination / "audit_progress.json"
    if run_audits:
        frame = pd.read_csv(feature_table, low_memory=False)
        requested = audit_labels or list(DEFAULT_AUDITS)
        available_labels = set(frame.columns)
        two_class_labels = {
            label: bool(
                binary_label_series(frame[label]).eq(0).any()
                and binary_label_series(frame[label]).eq(1).any()
            )
            for label in requested
            if label in available_labels
        }
        del frame
        gc.collect()
        for label in requested:
            if label not in available_labels:
                audits[label] = {"status": "skipped", "reason": "missing label column"}
                continue
            if not two_class_labels[label]:
                audits[label] = {
                    "status": "skipped",
                    "reason": "label lacks both positive and negative rows",
                }
                continue
            feature_set = DEFAULT_AUDITS.get(label, "spd_binding_nonleaky")
            audit_out = destination / "audits" / label
            audits[label] = {"status": "running", "out_dir": str(audit_out)}
            audit_progress_path.write_text(
                json.dumps(audits, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            try:
                result = run_ml_audit_suite(
                    feature_table,
                    label,
                    audit_out,
                    feature_set=feature_set,
                    claim_mode="exploratory",
                )
                audits[label] = {
                    "status": "completed",
                    "out_dir": str(audit_out),
                    "claim_readiness": result.get("claim_readiness"),
                    "manifest": str(audit_out / "ml_audit_suite_manifest.json"),
                }
                del result
                gc.collect()
            except Exception as exc:
                audits[label] = {
                    "status": "error",
                    "out_dir": str(audit_out),
                    "error": f"{type(exc).__name__}: {exc}",
                }
            audit_progress_path.write_text(
                json.dumps(audits, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )

    failed_checks = [row["check"] for row in checks if not row["passed"]]
    failed_checks.extend(
        f"audit:{label}"
        for label, result in audits.items()
        if result.get("status") == "error"
    )
    if failed_checks:
        status = "blocked"
    elif eligibility["excluded_rows"]:
        status = "ready_with_quarantined_rows"
    else:
        status = "ready"
    manifest = {
        "schema": SCHEMA,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "failed_checks": failed_checks,
        "inputs": {
            "base_table": {"path": str(base), "sha256": _sha256(base)},
            "addon_dir": str(addon),
            "addon_selection": addon_selection or {
                "mode": "explicit",
                "selected": str(addon),
            },
            "reference_manifest": {
                "path": str(manifest_path),
                "sha256": _sha256(manifest_path),
            },
            "selected_pairs": {"path": str(selected), "sha256": _sha256(selected)},
            "score_table": {"path": str(scores), "sha256": _sha256(scores)},
            "canonical_mapping": {
                "path": str(canonical_mapping),
                "sha256": _sha256(canonical_mapping),
            },
        },
        "reference_score_context": {
            "comparison_run_id": reference.get("comparison_run_id"),
            "schema": reference.get("schema"),
        },
        "merge": merged_manifest,
        "identity": {
            "release_status": identity_summary.get("release_status"),
            "counts": identity_summary.get("counts"),
        },
        "training_eligibility": eligibility,
        "validation": checks,
        "audits": audits,
        "outputs": {
            "raw": {"path": str(raw), "sha256": _sha256(raw)},
            "addon_rows": {"path": str(addon_rows), "sha256": _sha256(addon_rows)},
            "identity_reconciled": {
                "path": str(reconciled),
                "sha256": _sha256(reconciled),
            },
            "feature_ready": {
                "path": str(feature_table),
                "sha256": _sha256(feature_table),
            },
            "validation": str(destination / "merge_validation.csv"),
            "audit_progress": str(audit_progress_path) if run_audits else None,
        },
    }
    (destination / "spd_addon_merge_pipeline_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Reliably merge a scored SPD add-on pass, reconcile FDA identities, "
            "refresh ML metadata, and run post-ingestion audits."
        )
    )
    parser.add_argument("--base-table", required=True, type=Path)
    parser.add_argument(
        "--addon-dir",
        type=Path,
        default=None,
        help=(
            "Explicit scored add-on directory. When omitted, select the newest "
            "completed score-ready timestamped directory under --addon-root."
        ),
    )
    parser.add_argument(
        "--addon-root",
        type=Path,
        default=Path("data/AtlasSPD_phase1"),
    )
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--mapping", type=Path, default=None)
    parser.add_argument("--reference-manifest", type=Path, default=None)
    parser.add_argument("--selected-pairs", type=Path, default=None)
    parser.add_argument("--score-table", type=Path, default=None)
    parser.add_argument("--addon-name", default="external_activity_addon")
    parser.add_argument("--audit-label", action="append", default=None)
    parser.add_argument("--skip-audits", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.addon_dir is None:
        addon_dir, addon_selection = discover_latest_scored_addon(
            args.addon_root
        )
    else:
        addon_dir = args.addon_dir
        addon_selection = {
            "mode": "explicit",
            "selected": str(addon_dir.resolve()),
        }
    out_dir = args.out_dir or addon_dir / "merged_pipeline"
    manifest = run_spd_addon_merge_pipeline(
        base_table=args.base_table,
        addon_dir=addon_dir,
        out_dir=out_dir,
        mapping=args.mapping,
        reference_manifest=args.reference_manifest,
        selected_pairs=args.selected_pairs,
        score_table=args.score_table,
        addon_name=args.addon_name,
        run_audits=not args.skip_audits,
        audit_labels=args.audit_label,
        addon_selection=addon_selection,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 2 if manifest["status"] == "blocked" else 0


if __name__ == "__main__":
    raise SystemExit(main())
