from __future__ import annotations

import json
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
                "protein_id": "P53/Alpha",
                "drug_id": "Drug X",
                "final_status": "valid",
                "final_score": 2.75,
                "protein_rank": 1,
                "drug_rank": 1,
            },
            {
                "protein_id": "P53/Alpha",
                "drug_id": "Drug X",
                "final_status": "invalid",
                "final_score": 1.25,
                "protein_rank": 2,
                "drug_rank": 2,
                "failure_reason": "steric clash",
            },
        ],
    }

    manifest = generate_static_explorer(payload, tmp_path)

    assert manifest == {
        "release_id": "atlas-v0.1",
        "counts": {"proteins": 1, "drugs": 1, "pairs": 2},
        "entrypoint": "index.html",
    }
    assert json.loads((tmp_path / "site_manifest.json").read_text(encoding="utf-8")) == manifest

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
    assert 'href="../drugs/drug-x.html"' in pair_html
    assert pair_two_path.is_file()
    assert "2.75" in protein_html
    assert "1.25" in protein_html

    for relative_path in (
        "analysis.html",
        "assets/atlas.css",
        "assets/atlas.js",
        "assets/analysis_data.json",
    ):
        assert (tmp_path / relative_path).is_file()
