from __future__ import annotations

# ruff: noqa: E402

import sys
from pathlib import Path
from unittest.mock import Mock

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools import pdb_query


def _mock_response(status_code: int, payload: dict) -> Mock:
    response = Mock()
    response.status_code = status_code
    response.json.return_value = payload
    return response


def test_retry_backoff_on_transient_search_failure(monkeypatch: pytest.MonkeyPatch) -> None:
    session = Mock()
    session.post.side_effect = [
        _mock_response(500, {}),
        _mock_response(200, {"result_set": [{"identifier": "1ABC"}]}),
    ]

    sleep_calls: list[float] = []
    monkeypatch.setattr(pdb_query.time, "sleep", lambda value: sleep_calls.append(value))

    client = pdb_query.PDBQueryClient(
        session=session,
        rate_limit_rps=0.0,
        use_cache=False,
        max_retries=2,
    )
    payload = {"query": {"type": "terminal"}}
    result = client.search(payload)

    assert result == {"result_set": [{"identifier": "1ABC"}]}
    assert session.post.call_count == 2
    assert sleep_calls, "Expected backoff sleep between retry attempts"


def test_cache_hit_prevents_second_entry_fetch(tmp_path: Path) -> None:
    session = Mock()
    session.get.return_value = _mock_response(200, {"entry": {"id": "1ABC"}})

    client = pdb_query.PDBQueryClient(
        session=session,
        cache_dir=tmp_path / "pdb_query_cache",
        use_cache=True,
        cache_ttl_seconds=3600,
        rate_limit_rps=0.0,
    )

    first = pdb_query._fetch_entry("1ABC", client=client)
    second = pdb_query._fetch_entry("1ABC", client=client)

    assert first == {"entry": {"id": "1ABC"}}
    assert second == {"entry": {"id": "1ABC"}}
    assert session.get.call_count == 1


def test_ligand_filtering_excludes_trivial_only_and_keeps_nontrivial() -> None:
    all_trivial, nontrivial = pdb_query._classify_nonpolymer_comp_ids(["HOH", "SO4"])
    assert all_trivial == ["HOH", "SO4"]
    assert nontrivial == []

    all_ids, nontrivial_ids = pdb_query._classify_nonpolymer_comp_ids(["HOH", "LSD"])
    assert all_ids == ["HOH", "LSD"]
    assert nontrivial_ids == ["LSD"]


def test_ligand_filter_mode_relaxed_allows_glycan_components() -> None:
    all_ids, strict_nontrivial = pdb_query._classify_nonpolymer_comp_ids(
        ["NAG", "HOH"],
        ligand_filter_mode="strict",
    )
    assert all_ids == ["NAG", "HOH"]
    assert strict_nontrivial == []

    _, relaxed_nontrivial = pdb_query._classify_nonpolymer_comp_ids(
        ["NAG", "HOH"],
        ligand_filter_mode="relaxed",
    )
    assert relaxed_nontrivial == ["NAG"]


def test_text_ranked_entries_accepts_ligand_filter_mode(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(pdb_query, "_search_rcsb", lambda *args, **kwargs: ["1ABC"])
    seen: dict[str, str] = {}

    def _fake_eval(**kwargs):
        seen["ligand_filter_mode"] = kwargs["ligand_filter_mode"]
        return pdb_query.RankedEntry(
            gene=kwargs["gene_name"],
            pdb_id=kwargs["entry_id"],
            resolution=2.0,
            method="X-RAY",
            nonpolymer_comp_ids=["NAG"],
            nontrivial_comp_ids=["NAG"],
            rank_score=1.0,
            rank_reason={},
        )

    monkeypatch.setattr(pdb_query, "_evaluate_entry", _fake_eval)

    result = pdb_query.get_ranked_entries_for_text(
        "BRCA DNA repair",
        max_candidates=1,
        max_return=1,
        ligand_filter_mode="relaxed",
    )

    assert [entry.pdb_id for entry in result.entries] == ["1ABC"]
    assert seen["ligand_filter_mode"] == "relaxed"


def test_forced_inclusion_keeps_htr2b_lsd_structures(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(pdb_query, "_search_rcsb", lambda *args, **kwargs: ["1ABC"])

    def _fake_eval(**kwargs):
        return pdb_query.RankedEntry(
            gene=kwargs["gene_name"],
            pdb_id=kwargs["entry_id"],
            resolution=2.0,
            method="X-RAY DIFFRACTION",
            nonpolymer_comp_ids=["LIG"],
            nontrivial_comp_ids=["LIG"],
            rank_score=5.0,
            rank_reason={"score": 5.0},
        )

    monkeypatch.setattr(pdb_query, "_evaluate_entry", _fake_eval)
    monkeypatch.setattr(
        pdb_query,
        "_fetch_entry",
        lambda entry_id, client=None: {
            "rcsb_entry_info": {"nonpolymer_entity_ids": ["1"]},
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
        },
    )
    monkeypatch.setattr(
        pdb_query,
        "_list_nonpolymer_comp_ids",
        lambda entry_id, entry_json, client=None: ["LSD"],
    )

    pdb_ids = pdb_query.get_pdbs_for_gene("HTR2B", max_candidates=10, max_return=1)

    assert "5TVN" in pdb_ids
    assert "7SRS" in pdb_ids
    assert "1ABC" in pdb_ids


def test_ranking_deterministic_tie_break_by_pdb_id() -> None:
    rows = [
        pdb_query.RankedEntry(
            gene="GENE",
            pdb_id="2BBB",
            resolution=2.0,
            method="X-RAY",
            rank_score=10.0,
            rank_reason={},
        ),
        pdb_query.RankedEntry(
            gene="GENE",
            pdb_id="1AAA",
            resolution=2.0,
            method="X-RAY",
            rank_score=10.0,
            rank_reason={},
        ),
    ]

    ordered = pdb_query._rank_entries(rows, max_return=2)

    assert [row.pdb_id for row in ordered] == ["1AAA", "2BBB"]
