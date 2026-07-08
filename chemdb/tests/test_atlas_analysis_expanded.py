from __future__ import annotations

import os
import uuid
from pathlib import Path

import pandas as pd
import pytest

from analysis.controls.matched_controls import run_matched_controls
from analysis.external.spd import aggregate_spd_assays, build_spd_benchmark
from analysis.external.toxcast import build_toxcast_benchmark
from analysis.ml.spd_panel import _extract_target_id
from analysis.graph.graph_scores import compute_pair_mechanism_scores
from analysis.graph.mechanism_graph import build_mechanism_graph
from analysis.ml.build_ml_dataset import build_ml_dataset
from analysis.ml.leakage_checks import assert_no_leakage, find_leaky_features
from analysis.ml.splits import make_split
from analysis.priority_score import compute_priority_score
from analysis.schemas import summarize_pair_table, validate_pair_table


def _work_dir() -> Path:
    base = Path(os.environ.get("ATLAS_TEST_TMPDIR", "/stor/home/mpg2352/atlas_pytest_tmp"))
    path = base / f"analysis_{uuid.uuid4().hex}"
    path.mkdir(parents=True, exist_ok=False)
    return path


def _pair_df() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "drug_id": "D1",
                "target_id": "T1",
                "pdb_id": "P1",
                "atlas_score": 3.0,
                "mmgbsa_score": -20.0,
                "free_cmax": 10.0,
                "exposure_plausibility": 1.0,
                "tissue_expression": 5.0,
                "target_adr_evidence": 1,
                "pathway_evidence": 1,
                "structure_quality": 0.9,
                "protein_class": "kinase",
                "ligand_chemotype": "aryl",
                "literature_supported_label": 1,
            },
            {
                "drug_id": "D2",
                "target_id": "T1",
                "pdb_id": "P1",
                "atlas_score": 1.0,
                "mmgbsa_score": -8.0,
                "free_cmax": None,
                "exposure_plausibility": None,
                "tissue_expression": None,
                "target_adr_evidence": 0,
                "pathway_evidence": 0,
                "structure_quality": 0.7,
                "protein_class": "kinase",
                "ligand_chemotype": "alkyl",
                "literature_supported_label": 0,
            },
            {
                "drug_id": "D1",
                "target_id": "T2",
                "pdb_id": "P2",
                "atlas_score": 2.0,
                "mmgbsa_score": -12.0,
                "free_cmax": 5.0,
                "exposure_plausibility": 0.5,
                "tissue_expression": 2.0,
                "target_adr_evidence": 0,
                "pathway_evidence": 1,
                "structure_quality": 0.5,
                "protein_class": "gpcr",
                "ligand_chemotype": "aryl",
                "literature_supported_label": 0,
            },
            {
                "drug_id": "D3",
                "target_id": "T2",
                "pdb_id": "P2",
                "atlas_score": 0.2,
                "mmgbsa_score": -5.0,
                "free_cmax": 2.0,
                "exposure_plausibility": 0.0,
                "tissue_expression": 1.0,
                "target_adr_evidence": 0,
                "pathway_evidence": 0,
                "structure_quality": 0.4,
                "protein_class": "gpcr",
                "ligand_chemotype": "alkyl",
                "literature_supported_label": 0,
            },
        ]
    )


def test_priority_score_preserves_missingness_and_rescales() -> None:
    df = _pair_df()
    scored = compute_priority_score(df)

    missing_row = scored.loc[scored["drug_id"] == "D2"].iloc[0]
    assert "tissue_expression" in missing_row["priority_score_missing_components"]
    assert "exposure_plausibility" in missing_row["priority_score_missing_components"]
    assert 0.0 <= missing_row["priority_score"] <= 1.0
    assert pd.isna(df.loc[df["drug_id"] == "D2", "tissue_expression"]).iloc[0]


def test_pair_schema_summary_and_required_columns() -> None:
    rows = _pair_df().to_dict("records")

    assert validate_pair_table(rows) == []
    summary = summarize_pair_table(rows)
    assert summary["n_rows"] == 4
    assert summary["n_drugs"] == 3
    assert summary["missingness"]["free_cmax"] == 0.25


def test_spd_margin_labels_and_multiple_assay_aggregation() -> None:
    spd = pd.DataFrame(
        [
            {
                "drug_id": "D1",
                "target_id": "T1",
                "assay_id": "A_slow",
                "assay_name": "cell",
                "ac50_nM": 500.0,
                "free_cmax_nM": 10.0,
                "total_cmax_nM": 100.0,
                "source": "SPD",
                "assay_type": "cell",
            },
            {
                "drug_id": "D1",
                "target_id": "T1",
                "assay_id": "A_fast",
                "assay_name": "target bind",
                "ac50_nM": 50.0,
                "free_cmax_nM": 10.0,
                "total_cmax_nM": 100.0,
                "source": "SPD",
                "assay_type": "biochemical",
            },
        ]
    )

    out = aggregate_spd_assays(spd)

    assert len(out) == 1
    row = out.iloc[0]
    assert row["assay_id"] == "A_fast"
    assert row["assay_count"] == 2
    assert row["exposure_margin"] == 5.0
    assert bool(row["spd_exposure_relevant"]) is True


def test_spd_and_toxcast_benchmarks_join_by_drug_target() -> None:
    root = _work_dir()
    pair_path = root / "pairs.csv"
    spd_path = root / "spd.csv"
    toxcast_path = root / "toxcast.csv"
    _pair_df().to_csv(pair_path, index=False)
    pd.DataFrame(
        [
            {
                "drug_id": "D1",
                "target_id": "T1",
                "assay_id": "A1",
                "assay_name": "bind",
                "ac50_nM": 50.0,
                "free_cmax_nM": 10.0,
                "total_cmax_nM": 100.0,
                "source": "SPD",
            }
        ]
    ).to_csv(spd_path, index=False)
    pd.DataFrame(
        [
            {
                "drug_id": "D1",
                "target_id": "T1",
                "assay_id": "TX1",
                "assay_component": "component",
                "ac50_nM": 1000.0,
                "hit_call": "active",
                "modl_ga": 1.0,
                "source": "ToxCast",
            }
        ]
    ).to_csv(toxcast_path, index=False)

    spd = build_spd_benchmark(pair_path, spd_path, None, root / "spd_out.csv")
    toxcast = build_toxcast_benchmark(pair_path, toxcast_path, None, root / "toxcast_out.csv")

    joined_spd = spd.loc[(spd["drug_id"] == "D1") & (spd["target_id"] == "T1")].iloc[0]
    joined_toxcast = toxcast.loc[(toxcast["drug_id"] == "D1") & (toxcast["target_id"] == "T1")].iloc[0]
    assert bool(joined_spd["spd_exposure_relevant"]) is True
    assert bool(joined_toxcast["toxcast_active"]) is True
    assert bool(joined_toxcast["toxcast_potent"]) is True


def test_spd_assay_text_normalizes_target_id() -> None:
    raw = pd.DataFrame(
        {
            "HumanEntrezGeneSymbol(representative)": [None, None, None, None],
            "assay_group_name": [
                "ABCB11 inhibition",
                "ADRA1A BINDING (NIBR assay)",
                "AR (RAT) AGONIST",
                "ABCD (RAT) INHIBITION",
            ],
            "target_id": [
                "ABCB11 INHIBITION",
                "ADRA1A AGONIST (CRO ASSAY)",
                "AR (RAT) AGONIST",
                "ABCD (RAT) ANTAGONIST",
            ],
        }
    )
    target_ids = _extract_target_id(raw, use_assay_fallback=True)
    assert target_ids.tolist() == ["ABCB11", "ADRA1A", "AR", "ABCD"]


def test_mechanism_graph_scores_triad_completion() -> None:
    root = _work_dir()
    pair_path = root / "pairs.csv"
    sider_path = root / "sider.tsv"
    safety_path = root / "safety.tsv"
    reactome_path = root / "reactome.tsv"
    _pair_df().to_csv(pair_path, index=False)
    pd.DataFrame([{"drug_id": "D1", "adr": "QT prolongation", "source": "SIDER"}]).to_csv(
        sider_path, sep="\t", index=False
    )
    pd.DataFrame(
        [
            {
                "target_id": "T1",
                "adr": "QT prolongation",
                "confidence": 0.9,
                "pubmed_ids": "1",
                "source": "OpenTargets",
            }
        ]
    ).to_csv(safety_path, sep="\t", index=False)
    pd.DataFrame(
        [{"target_id": "T1", "pathway": "cardiac repolarization", "adr": "QT prolongation", "source": "Reactome"}]
    ).to_csv(reactome_path, sep="\t", index=False)

    nodes, edges = build_mechanism_graph(
        pair_path,
        sider_path,
        None,
        safety_path,
        reactome_path,
        None,
        root / "edges.csv",
        root / "nodes.csv",
    )
    scores = compute_pair_mechanism_scores(pair_path, root / "edges.csv", root / "scores.csv")

    assert {"drug", "target", "adr", "pathway"}.issubset(set(nodes["node_type"]))
    assert "target_has_safety_evidence" in set(edges["edge_type"])
    row = scores.loc[(scores["drug_id"] == "D1") & (scores["target_id"] == "T1")].iloc[0]
    assert row["triad_complete"] == 1
    assert row["mechanism_graph_score"] >= 2.0


def test_matched_controls_are_seed_reproducible_and_same_class() -> None:
    root = _work_dir()
    pair_path = root / "pairs.csv"
    _pair_df().to_csv(pair_path, index=False)

    first = run_matched_controls(pair_path, "atlas_score", 2, 2, ["protein_class"], root / "matched_a.csv", seed=7)
    second = run_matched_controls(pair_path, "atlas_score", 2, 2, ["protein_class"], root / "matched_b.csv", seed=7)

    assert first["matched_empirical_p"].tolist() == second["matched_empirical_p"].tolist()
    controls = pd.read_csv(root / "matched_controls.csv")
    assert not controls.empty
    assert controls["match_on"].eq("protein_class").all()


def test_ml_feature_set_excludes_ids_and_holdouts_do_not_overlap() -> None:
    root = _work_dir()
    path = root / "pairs.csv"
    out_path = root / "ml.csv"
    _pair_df().to_csv(path, index=False)

    dataset = build_ml_dataset(path, "literature_supported_label", "full_nonleaky", out_path)
    train_idx, test_idx = make_split(dataset, "drug_holdout", seed=1, test_fraction=0.5)
    train_drugs = set(dataset.loc[train_idx, "drug_id"])
    test_drugs = set(dataset.loc[test_idx, "drug_id"])

    assert "drug_id" in dataset.columns
    assert "atlas_score" in dataset.columns
    assert train_drugs.isdisjoint(test_drugs)
    assert find_leaky_features(["atlas_score", "drug_id", "literature_supported_label"]) == [
        "drug_id",
        "literature_supported_label",
    ]
    with pytest.raises(ValueError):
        assert_no_leakage(["target_id", "atlas_score"])
