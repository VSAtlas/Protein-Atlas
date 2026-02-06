import csv
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


def _write_csv(path: Path, headers: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)


def test_mapping_quality_gate_replaces_identifier_like_and_placeholders(tmp_path):
    in_csv = tmp_path / "in.csv"
    out_csv = tmp_path / "out.csv"
    cache_dir = tmp_path / "cache"
    headers = [
        "file_num",
        "display_name",
        "pubchem_cid_resolved",
        "pubchem_iupac_name",
        "pubchem_record_title",
        "pubchem_name",
        "pubchem_synonyms",
        "rxnorm_generic_name",
    ]
    _write_csv(
        in_csv,
        headers,
        [
            {
                "file_num": "1",
                "display_name": "UNK_1",
                "pubchem_cid_resolved": "111",
                "pubchem_iupac_name": "",
                "pubchem_record_title": "",
                "pubchem_name": "",
                "pubchem_synonyms": "",
                "rxnorm_generic_name": "",
            },
            {
                "file_num": "2",
                "display_name": "SCHEMBL123",
                "pubchem_cid_resolved": "222",
                "pubchem_iupac_name": "",
                "pubchem_record_title": "SCHEMBL123",
                "pubchem_name": "SCHEMBL123",
                "pubchem_synonyms": "",
                "rxnorm_generic_name": "",
            },
        ],
    )

    def fake_get(url: str, timeout: int = 10):
        if "/cid/111/property/Title,IUPACName/JSON" in url:
            return FakeResponse(
                {
                    "PropertyTable": {
                        "Properties": [
                            {
                                "CID": 111,
                                "Title": "",
                                "IUPACName": "2-acetoxybenzoic acid",
                            }
                        ]
                    }
                }
            )
        if "/cid/222/property/Title,IUPACName/JSON" in url:
            return FakeResponse(
                {
                    "PropertyTable": {
                        "Properties": [
                            {
                                "CID": 222,
                                "Title": "",
                                "IUPACName": "N-(4-hydroxyphenyl)acetamide",
                            }
                        ]
                    }
                }
            )
        if "/cid/111/synonyms/JSON" in url or "/cid/222/synonyms/JSON" in url:
            return FakeResponse({}, status_code=404)
        if "/pug_view/data/compound/111/JSON/" in url:
            return FakeResponse(
                {
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
                                                        {"String": "Aspirin"}
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
            )
        if "/pug_view/data/compound/222/JSON/" in url:
            return FakeResponse(
                {
                    "Record": {
                        "RecordTitle": "CHEMBL222",
                        "Section": [
                            {
                                "TOCHeading": "Drug and Medication Information",
                                "Section": [
                                    {
                                        "TOCHeading": "DrugBank",
                                        "Information": [
                                            {
                                                "Value": {
                                                    "StringWithMarkup": [
                                                        {"String": "Metformin"}
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
            )
        if "rxnav.nlm.nih.gov/REST/approximateTerm.json" in url:
            return FakeResponse({"approximateGroup": {"candidate": []}})
        if "fastidentity/cid/" in url:
            return FakeResponse({"IdentifierList": {"CID": []}})
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        exit_code = tool.main(
            [
                "--in_csv",
                str(in_csv),
                "--out_csv",
                str(out_csv),
                "--enable_unichem",
                "false",
                "--enable_pubchem_rerank",
                "true",
                "--enable_gsrs",
                "false",
                "--cache_dir",
                str(cache_dir),
                "--sleep",
                "0",
                "--max_workers",
                "1",
            ]
        )
    assert exit_code == 0

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert rows[0]["display_name"] == "Aspirin"
    assert rows[1]["display_name"] == "Metformin"

    for row in rows:
        display_name = (row.get("display_name") or "").strip()
        assert display_name
        assert not display_name.upper().startswith("UNK_")
        assert not tool.FDA_PLACEHOLDER_RE.match(display_name)
        assert not tool.is_identifier_like_name(display_name)
        assert not tool.is_bad_display_name(display_name, row.get("pubchem_iupac_name"))

