import json
from pathlib import Path

from analysis.reporting import target_safety_evidence as tse


def test_classify_safety_buckets_uses_event_and_tissue_keywords() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    buckets = tse.classify_safety_buckets(
        repo_root,
        event_text="QT prolongation and ventricular arrhythmia",
        tissue_text="CV",
    )
    assert buckets[0] == "Cardiotoxicity / QT"


def test_classify_safety_buckets_requires_specific_drug_ae_terms() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    buckets = tse.classify_safety_buckets(
        repo_root,
        event_text="liver panel review",
        evidence_source="drug_ae",
    )
    assert buckets == []


def test_load_target_safety_profile_prefers_cached_uniprot(tmp_path: Path) -> None:
    repo_root = tmp_path
    cache_dir = repo_root / "pathways" / "cache"
    cache_dir.mkdir(parents=True)
    payload = {
        "version": 1,
        "entries": [],
        "by_uniprot": {
            "P12931": {
                "approved_symbol": "SRC",
                "approved_name": "SRC proto-oncogene",
                "safety_buckets": ["CNS", "Hematologic"],
                "primary_display_safety": "CNS",
                "safety_confidence": "high",
                "safety_sources": ["Open Targets safety", "Open Targets drug adverse events"],
                "safety_bucket_scores": {"CNS": "10.0", "Hematologic": "2.0"},
                "safety_evidence_summary": "Open Targets safety: 2; Open Targets drug adverse events: 1",
            }
        },
        "by_gene": {},
    }
    (cache_dir / "target_safety_aggregated.json").write_text(
        json.dumps(payload), encoding="utf-8"
    )
    tse.clear_target_safety_cache()

    profile = tse.load_target_safety_profile(
        repo_root,
        target_name="SRC proto-oncogene",
        uniprots=["P12931"],
        gene_symbols=[],
    )
    assert profile is not None
    assert profile["primary_display_safety"] == "CNS"
    assert profile["safety_confidence"] == "high"
    assert profile["safety_sources"] == [
        "Open Targets drug adverse events",
        "Open Targets safety",
    ]


def test_build_target_safety_entry_aggregates_direct_and_drug_sources() -> None:
    repo_root = Path(__file__).resolve().parents[2]

    class FakeClient:
        def search_target(self, term: str):
            assert term == "P12931"
            return [{"id": "ENSG_TEST", "name": "SRC", "entity": "target"}]

        def fetch_target(self, ensembl_id: str, known_drug_size: int = 40):
            assert ensembl_id == "ENSG_TEST"
            assert known_drug_size >= 3
            return {
                "approvedSymbol": "SRC",
                "approvedName": "SRC proto-oncogene",
                "proteinIds": [{"id": "P12931", "source": "uniprot_swissprot"}],
                "safetyLiabilities": [
                    {
                        "event": "QT prolongation",
                        "datasource": "Brennan et al. (2024)",
                        "biosamples": [{"tissueLabel": "CV"}],
                        "studies": [{"type": "clinical"}],
                    }
                ],
                "knownDrugs": {
                    "rows": [
                        {
                            "drugId": "CHEMBL1",
                            "prefName": "Drug One",
                            "phase": 4,
                            "mechanismOfAction": "SRC inhibitor",
                        }
                    ]
                },
            }

        def fetch_drug(self, chembl_id: str):
            assert chembl_id == "CHEMBL1"
            return {
                "id": "CHEMBL1",
                "name": "Drug One",
                "adverseEvents": {
                    "rows": [
                        {"name": "thrombocytopenia", "count": 120, "logLR": 2.5}
                    ]
                },
            }

    entry = tse.build_target_safety_entry(
        repo_root,
        tse.TargetSpec(
            pdb_id="1O43",
            target_name="SRC proto-oncogene",
            uniprots=("P12931",),
            gene_symbols=(),
        ),
        client=FakeClient(),
        include_drug_evidence=True,
        max_approved_drugs=3,
    )
    assert entry is not None
    assert entry["primary_display_safety"] == "Cardiotoxicity / QT"
    assert "Hematologic" in entry["safety_buckets"]
    assert entry["secondary_safety_buckets"] == ["Hematologic"]
    assert entry["direct_safety_buckets"] == ["Cardiotoxicity / QT"]
    assert entry["drug_ae_buckets"] == ["Hematologic"]
    assert entry["safety_confidence"] == "high"
    assert entry["safety_sources"] == [
        "Open Targets drug adverse events",
        "Open Targets safety",
    ]
    assert entry["direct_liability_examples"] == ["QT prolongation"]
    assert entry["raw_adverse_event_examples"] == ["thrombocytopenia"]


def test_build_target_safety_entry_prefers_direct_endocrine_over_drug_ae_flood() -> None:
    repo_root = Path(__file__).resolve().parents[2]

    class FakeClient:
        def search_target(self, term: str):
            assert term == "P06401"
            return [{"id": "ENSG_PGR", "name": "PGR", "entity": "target"}]

        def fetch_target(self, ensembl_id: str, known_drug_size: int = 40):
            assert ensembl_id == "ENSG_PGR"
            return {
                "approvedSymbol": "PGR",
                "approvedName": "progesterone receptor",
                "proteinIds": [{"id": "P06401", "source": "uniprot_swissprot"}],
                "safetyLiabilities": [
                    {
                        "event": "progesterone receptor assay positive",
                        "datasource": "Open Targets",
                        "biosamples": [{"tissueLabel": "endocrine"}],
                        "studies": [{"type": "reproductive toxicology"}],
                    }
                ],
                "knownDrugs": {
                    "rows": [
                        {"drugId": "CHEMBL1", "prefName": "Drug One", "phase": 4},
                        {"drugId": "CHEMBL2", "prefName": "Drug Two", "phase": 4},
                    ]
                },
            }

        def fetch_drug(self, chembl_id: str):
            events = {
                "CHEMBL1": [
                    {"name": "hepatic function abnormal", "count": 240, "logLR": 2.4},
                    {"name": "drug-induced liver injury", "count": 180, "logLR": 2.2},
                    {"name": "nausea", "count": 150, "logLR": 1.4},
                ],
                "CHEMBL2": [
                    {"name": "hepatic failure", "count": 120, "logLR": 2.0},
                    {"name": "vomiting", "count": 110, "logLR": 1.3},
                ],
            }
            return {"id": chembl_id, "name": chembl_id, "adverseEvents": {"rows": events[chembl_id]}}

    entry = tse.build_target_safety_entry(
        repo_root,
        tse.TargetSpec(
            pdb_id="1A28",
            target_name="progesterone receptor",
            uniprots=("P06401",),
            gene_symbols=("PGR",),
        ),
        client=FakeClient(),
        include_drug_evidence=True,
        max_approved_drugs=4,
    )
    assert entry is not None
    assert entry["primary_display_safety"] == "Endocrine"
    assert "Hepatotoxicity" in entry["secondary_safety_buckets"]
    assert "GI" in entry["secondary_safety_buckets"]
    assert entry["direct_safety_buckets"] == ["Endocrine"]
    assert "Hepatotoxicity" in entry["drug_ae_buckets"]
    assert entry["direct_liability_examples"] == ["progesterone receptor assay positive"]
    assert "hepatic function abnormal" in entry["raw_adverse_event_examples"]
