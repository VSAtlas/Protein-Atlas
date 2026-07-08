from __future__ import annotations

from analysis.run_pair_significance import compute_pair_significance


def test_pair_significance_uses_higher_is_better_empirical_p_values() -> None:
    rows = [
        {"drug_id": "d1", "target_id": "t1", "pdb_id": "p1", "atlas_score": "3"},
        {"drug_id": "d1", "target_id": "t2", "pdb_id": "p2", "atlas_score": "1"},
        {"drug_id": "d2", "target_id": "t1", "pdb_id": "p1", "atlas_score": "2"},
    ]

    out = compute_pair_significance(rows, score_field="atlas_score", direction="higher")

    top = next(row for row in out if row["drug_id"] == "d1" and row["target_id"] == "t1")
    assert top["ligand_empirical_p"] == "0.666667"
    assert top["target_empirical_p"] == "0.666667"
    assert top["global_empirical_p"] == "0.5"
    assert "fdr_q_value" in top


def test_pair_significance_supports_lower_is_better_scores() -> None:
    rows = [
        {"drug_id": "d1", "target_id": "t1", "pdb_id": "p1", "mmgbsa_score": "-10"},
        {"drug_id": "d1", "target_id": "t2", "pdb_id": "p2", "mmgbsa_score": "-3"},
    ]

    out = compute_pair_significance(rows, score_field="mmgbsa_score", direction="lower")

    assert out[0]["global_empirical_p"] == "0.666667"

