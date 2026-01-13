import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEEPCOY_DIR = REPO_ROOT / "DeepCoy_duds"
if str(DEEPCOY_DIR) not in sys.path:
    sys.path.insert(0, str(DEEPCOY_DIR))

import external_sources as ext  # noqa: E402


def test_get_json_paged_follows_next_and_audits(monkeypatch):
    responses = [
        {
            "activities": [1, 2],
            "page_meta": {"next": "http://next", "total_count": 4, "page": 1},
        },
        {"activities": [3, 4], "page_meta": {}},
    ]
    calls = []

    def fake_get_json(url, params, headers, timeout, retries, session=None):
        calls.append((url, params))
        return responses[len(calls) - 1]

    monkeypatch.setattr(ext, "_get_json", fake_get_json)
    audit = ext.SourceAudit(enabled=True, max_lines=10)

    items, metas = ext._get_json_paged(
        "http://first",
        {"q": "x"},
        {},
        timeout=5,
        retries=0,
        item_key="activities",
        max_pages=5,
        audit=audit,
        source="chembl",
        purpose="chembl.test",
    )

    assert items == [1, 2, 3, 4]
    assert len(calls) == 2
    assert metas and metas[0].get("next") == "http://next"
    assert len(audit.records) == 2
