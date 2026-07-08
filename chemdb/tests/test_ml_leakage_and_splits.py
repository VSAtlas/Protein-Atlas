from __future__ import annotations

import pandas as pd
import pytest

from analysis.ml.audit_utils import source_holdout
from analysis.ml.leakage_checks import find_leaky_features
from analysis.ml.splits import make_split, split_overlap_summary
from analysis.ml.train_classifier_splits import _exclude_unlabeled_holdout_entities


def _split_frame() -> pd.DataFrame:
    rows = []
    for idx in range(60):
        rows.append(
            {
                "drug_id": f"D{idx % 15:02d}",
                "target_id": f"T{idx % 10:02d}",
                "scaffold_key": f"S{idx % 12:02d}",
                "chemical_cluster": f"C{idx % 11:02d}",
                "target_family": f"F{idx % 5:02d}",
                "protein_class": f"P{idx % 4:02d}",
                "label_source": "chembl;papyrus" if idx % 4 == 0 else "toxcast",
                "activity_publication_year": 2000 + idx,
                "y": idx % 2,
            }
        )
    return pd.DataFrame(rows)


def test_find_leaky_features_blocks_source_status_rescore_and_label_columns():
    features = [
        "label_source",
        "source_family",
        "assay_type",
        "activity_type",
        "tox_label",
        "mechanism_label_status",
        "activity_threshold_nM",
        "pubchem_assay_ids",
        "SCORCH_score_used",
        "rescore_available",
        "uniprot_accession",
    ]
    assert set(features).issubset(set(find_leaky_features(features)))


@pytest.mark.parametrize(
    "split_mode",
    ["drug_holdout", "target_holdout", "scaffold_holdout", "chemical_cluster_holdout", "target_family_holdout"],
)
def test_entity_holdout_splits_have_zero_declared_overlap(split_mode):
    frame = _split_frame()
    train_idx, test_idx = make_split(frame, split_mode=split_mode, seed=11)
    summary = split_overlap_summary(frame.loc[train_idx], frame.loc[test_idx], split_mode)
    assert summary["passes_holdout"] is True


def test_temporal_holdout_orders_train_before_test_and_source_holdout_tokenizes_sources():
    frame = _split_frame()
    train_idx, test_idx = make_split(frame, split_mode="temporal_holdout", seed=11)
    assert frame.loc[train_idx, "activity_publication_year"].max() <= frame.loc[test_idx, "activity_publication_year"].min()

    train_idx, test_idx, summary = source_holdout(frame, "label_source", seed=3, test_fraction=0.5)
    assert summary["passes_holdout"] is True
    held_out = set(summary["held_out_sources"])
    train_sources = set(";".join(frame.loc[train_idx, "label_source"]).split(";"))
    assert not (held_out & train_sources)


def test_exclude_unlabeled_pool_removes_heldout_entities_by_split_mode():
    unlabeled = pd.DataFrame({"drug_id": ["D1", "D2", "D3"], "target_id": ["T1", "T2", "T3"]})
    test = pd.DataFrame({"drug_id": ["D2"], "target_id": ["T9"]})
    filtered = _exclude_unlabeled_holdout_entities(unlabeled, test, "drug_holdout")
    assert filtered["drug_id"].tolist() == ["D1", "D3"]

