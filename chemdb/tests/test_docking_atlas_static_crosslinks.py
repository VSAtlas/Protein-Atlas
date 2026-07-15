from __future__ import annotations

import json
import re
from pathlib import Path

from analysis.reporting.docking_atlas_static import generate_static_explorer


def test_static_explorer_links_protein_drug_and_pair_pages(tmp_path: Path) -> None:
    payload = {
        "release": {"id": "atlas-v0.1", "title": "Atlas conference release"},
        "proteins": [
            {
                "id": "P53/Alpha",
                "display_name": "Protein Alpha",
                "native_redock_status": "qualified",
            }
        ],
        "drugs": [{"id": "Drug X", "display_name": "Drug Example"}],
        "pairs": [
            {
                "pair_cell_id": 11,
                "protein_id": "P53/Alpha",
                "drug_id": "Drug X",
                "final_status": "valid",
                "final_score": 2.75,
                "final_score_source_effective": "declared:consensus_vs_decoy_z",
                "final_score_source_family": "consensus_decoy_standardized",
                "rank_eligible": 1,
                "ranking_track": "qualified_holo",
                "protein_rank": 1,
                "drug_rank": 1,
            },
            {
                "pair_cell_id": 12,
                "protein_id": "P53/Alpha",
                "drug_id": "Drug X",
                "final_status": "invalid",
                "final_score": 1.25,
                "protein_rank": 2,
                "drug_rank": 2,
                "failure_reason": "steric clash",
            },
        ],
        "artifacts": [
            {
                "artifact_id": 7,
                "pair_cell_id": 11,
                "member_name": "pose.sdf",
                "sha256": "a" * 64,
            }
        ],
    }

    manifest = generate_static_explorer(payload, tmp_path)

    assert manifest == {
        "release_id": "atlas-v0.1",
        "counts": {"proteins": 1, "drugs": 1, "pairs": 2},
        "entrypoint": "index.html",
    }
    assert (
        json.loads((tmp_path / "site_manifest.json").read_text(encoding="utf-8"))
        == manifest
    )

    index_html = (tmp_path / "index.html").read_text(encoding="utf-8")
    protein_html = (tmp_path / "proteins" / "p53-alpha.html").read_text(
        encoding="utf-8"
    )
    drug_html = (tmp_path / "drugs" / "drug-x.html").read_text(encoding="utf-8")
    pair_one_path = tmp_path / "pairs" / "p53-alpha-drug-x.html"
    pair_two_path = tmp_path / "pairs" / "p53-alpha-drug-x-2.html"
    pair_html = pair_one_path.read_text(encoding="utf-8")

    assert 'href="proteins/p53-alpha.html"' in index_html
    assert 'href="drugs/drug-x.html"' in index_html
    assert 'href="../drugs/drug-x.html"' in protein_html
    first_pair_link = 'href="../pairs/p53-alpha-drug-x.html"'
    second_pair_link = 'href="../pairs/p53-alpha-drug-x-2.html"'
    assert protein_html.count(first_pair_link) == 1
    assert protein_html.count(second_pair_link) == 1
    assert protein_html.index(first_pair_link) < protein_html.index(second_pair_link)
    assert 'href="../proteins/p53-alpha.html"' in drug_html
    assert drug_html.count(first_pair_link) == 1
    assert drug_html.count(second_pair_link) == 1
    assert drug_html.index(first_pair_link) < drug_html.index(second_pair_link)
    assert 'href="../proteins/p53-alpha.html"' in pair_html
    assert "pose.sdf" in pair_html
    assert "a" * 64 in pair_html
    assert 'href="../drugs/drug-x.html"' in pair_html
    assert pair_two_path.is_file()
    assert "2.75" in protein_html
    assert "1.25" in protein_html
    assert "Qualified HOLO rankings" in protein_html
    assert "Exploratory APO rankings" in protein_html
    assert "Unranked and excluded pair cells" in protein_html
    assert "declared:consensus_vs_decoy_z" in protein_html
    assert "consensus_decoy_standardized" in protein_html

    for relative_path in (
        "analysis.html",
        "assets/atlas.css",
        "assets/atlas.js",
        "assets/analysis_data.json",
    ):
        assert (tmp_path / relative_path).is_file()


def test_static_explorer_progressively_renders_more_than_one_page(
    tmp_path: Path,
) -> None:
    pair_count = 101
    payload = {
        "release": {"id": "atlas-paged", "title": "Paged Atlas"},
        "proteins": [{"id": "target-a", "display_name": "Target A"}],
        "drugs": [
            {"id": f"drug-{index:03d}", "display_name": f"Drug {index:03d}"}
            for index in range(1, pair_count + 1)
        ],
        "pairs": [
            {
                "pair_cell_id": index,
                "protein_id": "target-a",
                "drug_id": f"drug-{index:03d}",
                "final_status": "valid",
                "final_score": float(pair_count - index),
                "rank_eligible": 1,
                "ranking_track": "qualified_holo",
                "rank_within_receptor": index,
                "rank_across_receptors": 1,
            }
            for index in range(1, pair_count + 1)
        ],
    }

    generate_static_explorer(payload, tmp_path)

    index_html = (tmp_path / "index.html").read_text(encoding="utf-8")
    active_cards = re.search(
        r'<div id="drug-list" class="card-grid">(.*?)</div><div class="pager"',
        index_html,
        re.DOTALL,
    )
    assert active_cards is not None
    assert active_cards.group(1).count('class="entity-card"') == 100
    card_source = re.search(
        r'<script type="application/json" data-card-source="drug-list">(.*?)</script>',
        index_html,
        re.DOTALL,
    )
    assert card_source is not None
    card_rows = json.loads(card_source.group(1))
    assert len(card_rows) == pair_count
    assert card_rows[-1]["search"] == "drug-101 drug 101"
    assert 'data-card-pager="drug-list"' in index_html
    assert "data-page-prev" in index_html
    assert "data-page-next" in index_html
    assert "data-page-status" in index_html

    protein_html = (tmp_path / "proteins" / "target-a.html").read_text(encoding="utf-8")
    active_pair_rows = re.search(r"<tbody>(.*?)</tbody>", protein_html, re.DOTALL)
    assert active_pair_rows is not None
    assert active_pair_rows.group(1).count("<tr ") == 100
    pair_source = re.search(
        r'<script type="application/json" data-pair-source="qualified-holo-table">(.*?)</script>',
        protein_html,
        re.DOTALL,
    )
    assert pair_source is not None
    pair_rows = json.loads(pair_source.group(1))
    assert len(pair_rows) == pair_count
    assert all(row["eligible"] is True for row in pair_rows)
    assert "drug 101 drug-101 valid" in pair_rows[-1]["search"]
    assert 'data-pair-pager="qualified-holo-table"' in protein_html


def test_static_explorer_separates_holo_apo_and_unqualified_rank_tracks(
    tmp_path: Path,
) -> None:
    payload = {
        "release": {"id": "atlas-tracks", "title": "Rank tracks"},
        "proteins": [{"id": "target-a", "display_name": "Target A"}],
        "drugs": [
            {"id": "holo", "display_name": "Drug Holo"},
            {"id": "apo", "display_name": "Drug APO"},
            {"id": "untracked", "display_name": "Drug Untracked"},
        ],
        "pairs": [
            {
                "pair_cell_id": 1,
                "protein_id": "target-a",
                "drug_id": "holo",
                "final_status": "valid",
                "final_score": 3.0,
                "final_score_source_effective": "declared:consensus_vs_decoy_z",
                "final_score_source_family": "consensus_decoy_standardized",
                "ranking_track": "qualified_holo",
                "rank_eligible": 1,
                "rank_within_receptor": 1,
            },
            {
                "pair_cell_id": 2,
                "protein_id": "target-a",
                "drug_id": "apo",
                "final_status": "valid",
                "final_score": 2.0,
                "final_score_source_effective": "reconstructed:scorch_vs_decoy_z",
                "final_score_source_family": "scorch_decoy_standardized",
                "ranking_track": "exploratory_apo",
                "apo_exploratory_rank_eligible": 1,
                "apo_rank_within_receptor": 1,
            },
            {
                "pair_cell_id": 3,
                "protein_id": "target-a",
                "drug_id": "untracked",
                "final_status": "valid",
                "final_score": 1.0,
                "rank_eligible": 1,
                "rank_within_receptor": 777,
            },
        ],
    }

    generate_static_explorer(payload, tmp_path)
    protein_html = (tmp_path / "proteins" / "target-a.html").read_text(encoding="utf-8")

    def rows(table_id: str) -> list[dict[str, object]]:
        match = re.search(
            rf'<script type="application/json" data-pair-source="{table_id}">(.*?)</script>',
            protein_html,
            re.DOTALL,
        )
        assert match is not None
        return json.loads(match.group(1))

    holo, apo, excluded = (
        rows("qualified-holo-table"),
        rows("exploratory-apo-table"),
        rows("unranked-pair-table"),
    )
    assert len(holo) == len(apo) == len(excluded) == 1
    assert "Drug Holo" in str(holo[0]["html"])
    assert "declared:consensus_vs_decoy_z" in str(holo[0]["html"])
    assert "consensus_decoy_standardized" in str(holo[0]["html"])
    assert "Drug APO" in str(apo[0]["html"])
    assert "reconstructed:scorch_vs_decoy_z" in str(apo[0]["html"])
    assert "Drug Untracked" in str(excluded[0]["html"])
    assert ">777<" not in str(excluded[0]["html"])
