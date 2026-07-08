from pathlib import Path
from unittest import mock

from chemdb.tools import fill_fda_mapping_names as tool


class FakeResponse:
    def __init__(self, payload, status_code=200):
        self._payload = payload
        self.status_code = status_code
        self.headers = {}

    def json(self):
        return self._payload


def _make_client(cache_dir: Path) -> tool.CachedJsonClient:
    return tool.CachedJsonClient(
        cache_dir,
        sleep=0.0,
        rate_limiter=None,
        cache_mode="off",
        max_retries=0,
        backoff_base=0.0,
        backoff_max=0.0,
    )


def test_pugview_drug_sections_override_identifier_record_title(tmp_path):
    row = {
        "pubchem_record_title": "",
        "pubchem_name": "",
        "pubchem_unii_list": "",
        "rxnorm_generic_name": "",
    }
    payload = {
        "Record": {
            "RecordTitle": "SCHEMBL3648819",
            "Section": [
                {
                    "TOCHeading": "Drug and Medication Information",
                    "Section": [
                        {
                            "TOCHeading": "Generic Name",
                            "Information": [
                                {
                                    "Value": {
                                        "StringWithMarkup": [
                                            {"String": "Artesunate"},
                                        ]
                                    }
                                }
                            ],
                        }
                    ],
                }
            ],
        }
    }

    def fake_get(url: str, timeout: int = 10):
        if "/pug_view/data/compound/101/JSON/" in url:
            return FakeResponse(payload)
        if "fastidentity/cid/101/cids/JSON" in url:
            return FakeResponse({"IdentifierList": {"CID": []}})
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        result = tool._resolve_pubchem_last_mile(
            row=row,
            pubchem_cid="101",
            pubchem_iupac_name="",
            pubchem_client=_make_client(tmp_path / "cache"),
        )

    assert result.name == "Artesunate"
    assert not tool.is_identifier_like_name(result.name)
    assert row["pubchem_name"] == "Artesunate"


def test_pugview_returns_empty_when_only_identifier_or_iupac_candidates(tmp_path):
    row = {
        "pubchem_record_title": "",
        "pubchem_name": "",
        "pubchem_unii_list": "",
        "rxnorm_generic_name": "",
    }
    payload = {
        "Record": {
            "RecordTitle": "SCHEMBL3648819",
            "Section": [
                {
                    "TOCHeading": "Names and Identifiers",
                    "Section": [
                        {
                            "TOCHeading": "Computed Descriptors",
                            "Section": [
                                {
                                    "TOCHeading": "IUPAC Name",
                                    "Information": [
                                        {
                                            "Value": {
                                                "StringWithMarkup": [
                                                    {
                                                        "String": "2-acetoxybenzoic acid"
                                                    }
                                                ]
                                            }
                                        }
                                    ],
                                }
                            ],
                        },
                        {
                            "TOCHeading": "Other Identifiers",
                            "Information": [
                                {
                                    "Value": {
                                        "StringWithMarkup": [
                                            {"String": "CHEMBL12345"},
                                            {"String": "CAS 50-78-2"},
                                        ]
                                    }
                                }
                            ],
                        },
                    ],
                }
            ],
        }
    }

    def fake_get(url: str, timeout: int = 10):
        if "/pug_view/data/compound/202/JSON/" in url:
            return FakeResponse(payload)
        if "fastidentity/cid/202/cids/JSON" in url:
            return FakeResponse({"IdentifierList": {"CID": []}})
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        result = tool._resolve_pubchem_last_mile(
            row=row,
            pubchem_cid="202",
            pubchem_iupac_name="",
            pubchem_client=_make_client(tmp_path / "cache"),
        )

    assert result.name == ""
