from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.spd_external_four_state_merge import _combine_activity_labels
from analysis.ml.spd_four_expert_tables import BINDING_LABEL_POLICY_VERSION, _binding_label


def relabel_spd_binding_censor_aware(
    dataset_path: str | Path,
    out_path: str | Path,
    *,
    active_um: float = 1.0,
    inactive_um: float = 10.0,
) -> dict[str, Any]:
    source = Path(dataset_path)
    out = Path(out_path)
    frame = pd.read_csv(source, low_memory=False)
    required = {"spd_binding_label", "spd_ac50_uM", "spd_activity_relation"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"dataset missing required columns: {', '.join(missing)}")

    old_spd = pd.to_numeric(frame["spd_binding_label"], errors="coerce")
    new_spd = _binding_label(frame, active_um=active_um, inactive_um=inactive_um)
    frame["spd_binding_label_pre_censor_repair"] = frame["spd_binding_label"]
    frame["spd_binding_label"] = new_spd
    frame["spd_binding_label_policy_version"] = BINDING_LABEL_POLICY_VERSION
    transition = pd.Series("unchanged", index=frame.index, dtype="object")
    transition.loc[old_spd.notna() & new_spd.isna()] = "labeled_to_unknown"
    transition.loc[old_spd.isna() & new_spd.notna()] = "unknown_to_labeled"
    transition.loc[old_spd.notna() & new_spd.notna() & old_spd.ne(new_spd)] = (
        "label_changed"
    )
    frame["spd_binding_label_censor_repair_status"] = transition

    combined_recomputed = "combined_activity_ml_label" in frame.columns
    if combined_recomputed:
        old_combined = pd.to_numeric(
            frame["combined_activity_ml_label"], errors="coerce"
        )
    else:
        old_combined = pd.Series(float("nan"), index=frame.index, dtype=float)
    if combined_recomputed:
        for column in (
            "combined_activity_label",
            "combined_activity_ml_label",
            "combined_activity_label_status",
            "combined_activity_label_policy",
        ):
            if column in frame.columns:
                frame[f"{column}_pre_censor_repair"] = frame[column]
        frame = _combine_activity_labels(frame)
    if "combined_activity_ml_label" in frame.columns:
        new_combined = pd.to_numeric(
            frame["combined_activity_ml_label"], errors="coerce"
        )
    else:
        new_combined = pd.Series(float("nan"), index=frame.index, dtype=float)

    changed = old_spd.fillna(-99).ne(pd.to_numeric(new_spd, errors="coerce").fillna(-99))
    audit_columns = [
        column
        for column in (
            "drug_id",
            "target_id",
            "pdb_id",
            "spd_ac50_uM",
            "spd_activity_relation",
            "spd_binding_label_pre_censor_repair",
            "spd_binding_label",
            "spd_binding_label_censor_repair_status",
            "combined_activity_ml_label_pre_censor_repair",
            "combined_activity_ml_label",
        )
        if column in frame.columns
    ]
    changed_rows = frame.loc[changed, audit_columns].copy()
    transitions = (
        pd.DataFrame(
            {
                "old_spd_binding_label": old_spd,
                "new_spd_binding_label": new_spd,
                "activity_relation": frame["spd_activity_relation"],
            }
        )
        .value_counts(dropna=False)
        .rename("n_rows")
        .reset_index()
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out, index=False)
    changed_path = out.with_name(f"{out.stem}.censor_label_changes.csv")
    transitions_path = out.with_name(f"{out.stem}.censor_label_transitions.csv")
    changed_rows.to_csv(changed_path, index=False)
    transitions.to_csv(transitions_path, index=False)
    manifest: dict[str, Any] = {
        "dataset": str(source),
        "output": str(out),
        "policy_version": BINDING_LABEL_POLICY_VERSION,
        "active_um": active_um,
        "inactive_um": inactive_um,
        "n_rows": int(len(frame)),
        "n_old_positive": int(old_spd.eq(1).sum()),
        "n_new_positive": int(new_spd.eq(1).sum()),
        "n_old_negative": int(old_spd.eq(0).sum()),
        "n_new_negative": int(new_spd.eq(0).sum()),
        "n_old_unknown": int(old_spd.isna().sum()),
        "n_new_unknown": int(new_spd.isna().sum()),
        "n_spd_labels_changed": int(changed.sum()),
        "combined_activity_recomputed": combined_recomputed,
        "n_combined_labels_changed": int(
            old_combined.fillna(-99).ne(new_combined.fillna(-99)).sum()
        ),
        "policy": (
            "Exact or upper-bound AC50 at/below the active threshold is positive; "
            "exact or lower-bound AC50 at/above the inactive threshold is negative; "
            "bounds crossing a threshold remain unknown."
        ),
        "outputs": {
            "dataset": str(out),
            "changed_rows": str(changed_path),
            "transitions": str(transitions_path),
        },
    }
    manifest_path = out.with_name(f"{out.stem}.censor_relabel_manifest.json")
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest
