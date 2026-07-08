from __future__ import annotations

import pandas as pd

from analysis.ml.audit_suite import run_ml_audit_suite


def _audit_rows(n: int = 72) -> list[dict[str, object]]:
    rows = []
    for idx in range(n):
        rows.append(
            {
                "drug_id": f"D{idx % 18:02d}",
                "target_id": f"T{idx % 9:02d}",
                "label_source": ["chembl", "papyrus", "toxcast"][idx % 3],
                "source_family": "curated" if idx % 3 != 2 else "hts",
                "scaffold_key": f"S{idx % 12:02d}",
                "chemical_cluster": f"C{idx % 11:02d}",
                "target_family": f"F{idx % 5:02d}",
                "protein_class": f"P{idx % 4:02d}",
                "activity_publication_year": 2000 + (idx % 20),
                "atlas_score": float(idx % 7),
                "consensus_score": float((idx * 5) % 9),
                "y": int(idx % 2 == 0),
            }
        )
    return rows


def test_run_ml_audit_suite_writes_contract_outputs(tmp_path):
    dataset = tmp_path / "audit.csv"
    pd.DataFrame(_audit_rows()).to_csv(dataset, index=False)
    out = tmp_path / "audit"

    manifest = run_ml_audit_suite(
        dataset,
        "y",
        out,
        feature_set="pilot_nonleaky",
        splits=["random", "drug_holdout", "target_holdout", "scaffold_holdout", "chemical_cluster_holdout", "target_family_holdout", "temporal_holdout", "source_holdout"],
        source_col="label_source",
    )

    assert manifest["outputs"]["data_card"].endswith("ml_data_card.json")
    for rel in [
        "ml_audit_suite_manifest.json",
        "ml_data_card.json",
        "ml_audit_claim_readiness.json",
        "leakage/split_overlap_summary.csv",
        "leakage/field_overlap_by_split.csv",
        "leakage/label_balance_by_split.csv",
        "leakage/source_label_balance.csv",
        "leakage/high_risk_columns.csv",
        "leakage/leakage_overlap_audit.json",
        "decoy_bias/decoy_bias_summary.json",
    ]:
        assert (out / rel).exists(), rel

