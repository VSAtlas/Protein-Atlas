from pathlib import Path

import pytest

from analysis.reporting import target_annotation_groups as tag
from analysis.reporting.target_annotation_groups import resolve_target_annotations


def test_resolve_target_annotations_matches_exact_uniprot() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    adme_meta = resolve_target_annotations(
        repo_root,
        target_name="Pregnane X receptor",
        uniprots=["O75469"],
    )
    assert adme_meta["adme_category"] == "Modifier/regulator"
    assert adme_meta["primary_display_safety"] == "Unassigned"
    assert adme_meta["safety_confidence"] == "low"
    assert adme_meta["safety_sources"] == ["Local heuristic catalog"]


def test_resolve_target_annotations_matches_keywords() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    safety_meta = resolve_target_annotations(
        repo_root,
        target_name="Epidermal growth factor receptor",
        uniprots=[],
    )
    assert safety_meta["primary_display_safety"] == "Dermatologic"
    assert safety_meta["safety_buckets"] == ["Dermatologic"]
    assert safety_meta["safety_confidence"] == "low"
    assert safety_meta["safety_sources"] == ["Local heuristic catalog"]


def test_resolve_target_annotations_uses_biology_prior_when_cached_profile_is_drug_only(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    repo_root = Path(__file__).resolve().parents[2]

    monkeypatch.setattr(
        tag,
        "load_target_safety_profile",
        lambda *_args, **_kwargs: {
            "safety_buckets": ["Hepatotoxicity", "GI"],
            "primary_display_safety": "Hepatotoxicity",
            "secondary_safety_buckets": ["GI"],
            "direct_safety_buckets": [],
            "drug_ae_buckets": ["Hepatotoxicity", "GI"],
            "safety_confidence": "low",
            "safety_sources": ["Open Targets drug adverse events"],
            "safety_bucket_scores": {"Hepatotoxicity": "5.000", "GI": "2.000"},
            "safety_evidence_summary": "Open Targets drug adverse events: 2",
            "raw_adverse_event_examples": ["hepatic function abnormal"],
            "direct_liability_examples": [],
        },
    )

    safety_meta = resolve_target_annotations(
        repo_root,
        target_name="progesterone receptor",
        uniprots=["P06401"],
        gene_symbols=["PGR"],
    )
    assert safety_meta["primary_display_safety"] == "Endocrine"
    assert safety_meta["safety_buckets"][0] == "Endocrine"
    assert "Hepatotoxicity" in safety_meta["safety_buckets"]
