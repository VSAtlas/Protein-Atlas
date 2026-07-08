from __future__ import annotations

from analysis.run_enrichment import run_enrichment


def test_enrichment_reports_required_metrics_and_shuffle_controls() -> None:
    rows = [
        {"drug_id": "d1", "target_id": "t1", "atlas_score": "4", "literature_supported_label": "1"},
        {"drug_id": "d2", "target_id": "t2", "atlas_score": "3", "literature_supported_label": "0"},
        {"drug_id": "d3", "target_id": "t3", "atlas_score": "2", "literature_supported_label": "0"},
        {"drug_id": "d4", "target_id": "t4", "atlas_score": "1", "literature_supported_label": "0"},
    ]

    summary, permutations = run_enrichment(
        rows,
        score="atlas_score",
        label="literature_supported_label",
        n_permutations=5,
        n_bootstraps=5,
        top_ks=[2],
        seed=1,
    )

    metrics = {row["metric"] for row in summary}
    assert {"EF@1%", "EF@5%", "EF@10%", "AUPRC", "Precision@2", "Recall@2"} <= metrics
    assert {row["null_model"] for row in permutations} == {"drug_label_shuffle", "target_label_shuffle"}

