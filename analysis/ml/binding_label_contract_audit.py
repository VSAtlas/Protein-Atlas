from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.spd_four_expert_tables import BINDING_LABEL_POLICY_VERSION, _binding_label


def audit_binding_label_contract(
    dataset_path: str | Path,
    out_dir: str | Path,
    *,
    active_um: float = 1.0,
    inactive_um: float = 10.0,
) -> dict[str, Any]:
    dataset = Path(dataset_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    header = pd.read_csv(dataset, nrows=0)
    required = {"spd_binding_label", "spd_ac50_uM", "spd_activity_relation"}
    missing = sorted(required - set(header.columns))
    manifest: dict[str, Any]
    if missing:
        manifest = {
            "status": "skipped",
            "reason": f"missing required columns: {', '.join(missing)}",
            "dataset": str(dataset),
        }
        (out / "binding_label_contract_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
        )
        return manifest

    optional = [
        column
        for column in (
            "drug_id",
            "target_id",
            "pdb_id",
            "spd_binding_label_policy_version",
        )
        if column in header.columns
    ]
    frame = pd.read_csv(dataset, usecols=[*required, *optional], low_memory=False)
    derivation = pd.DataFrame(
        {
            "spd_binding_label": pd.Series(pd.NA, index=frame.index, dtype="Int64"),
            "spd_ac50_uM": frame["spd_ac50_uM"],
            "spd_activity_relation": frame["spd_activity_relation"],
        }
    )
    observed = pd.to_numeric(frame["spd_binding_label"], errors="coerce").astype("Int64")
    version = frame.get(
        "spd_binding_label_policy_version",
        pd.Series("", index=frame.index),
    ).fillna("").astype(str)
    spd_measurement = pd.to_numeric(frame["spd_ac50_uM"], errors="coerce").notna()
    in_scope = spd_measurement | observed.notna() | version.ne("")

    expected = _binding_label(
        derivation,
        active_um=active_um,
        inactive_um=inactive_um,
    )
    relation = frame["spd_activity_relation"].astype("string").str.strip().replace(
        {"\u2264": "<=", "\u2265": ">=", "==": "=", "eq": "=", "lt": "<", "gt": ">"}
    )
    interpretable = in_scope & spd_measurement & relation.isin(
        ["=", "<", "<=", ">", ">="]
    )
    contradiction = interpretable & (
        (expected.notna() & ~observed.eq(expected).fillna(False))
        | (expected.isna() & observed.notna())
    )
    policy_mismatch = in_scope & version.ne(BINDING_LABEL_POLICY_VERSION)

    details = frame.loc[contradiction | policy_mismatch, optional + list(required)].copy()
    details["expected_spd_binding_label"] = expected.loc[details.index]
    details["binding_label_contradiction"] = contradiction.loc[details.index]
    details["binding_policy_version_mismatch"] = policy_mismatch.loc[details.index]
    details.to_csv(out / "binding_label_contract_violations.csv", index=False)
    n_contradictions = int(contradiction.sum())
    n_policy_mismatch = int(policy_mismatch.sum())
    manifest = {
        "status": "failed" if n_contradictions or n_policy_mismatch else "passed",
        "dataset": str(dataset),
        "policy_version_expected": BINDING_LABEL_POLICY_VERSION,
        "active_um": float(active_um),
        "inactive_um": float(inactive_um),
        "n_rows": int(len(frame)),
        "n_rows_audited": int(in_scope.sum()),
        "n_rows_out_of_scope_external": int((~in_scope).sum()),
        "n_interpretable_measurements": int(interpretable.sum()),
        "n_label_contradictions": n_contradictions,
        "n_policy_version_mismatch": n_policy_mismatch,
        "outputs": {
            "violations": str(out / "binding_label_contract_violations.csv"),
            "manifest": str(out / "binding_label_contract_manifest.json"),
        },
    }
    (out / "binding_label_contract_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    return manifest
