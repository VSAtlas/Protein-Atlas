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


def _run_tool(in_csv: Path, out_csv: Path, cache_dir: Path) -> int:
    drugcentral_tsv = in_csv.parent / "drugcentral_empty.tsv"
    with drugcentral_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["INN", "INCHIKEY"], delimiter="\t")
        writer.writeheader()

    return tool.main(
        [
            "--in_csv",
            str(in_csv),
            "--out_csv",
            str(out_csv),
            "--drugcentral_structures_tsv",
            str(drugcentral_tsv),
            "--enable_unichem",
            "false",
            "--enable_pubchem_rerank",
            "true",
            "--cache_dir",
            str(cache_dir),
            "--sleep",
            "0",
            "--max_workers",
            "1",
        ]
    )


def test_pubchem_cas_last_mile_resolves_hard_cases(tmp_path):
    in_csv = tmp_path / "in.csv"
    out_csv = tmp_path / "out.csv"
    cache_dir = tmp_path / "cache"
    headers = [
        "file_num",
        "sdf_title",
        "display_name",
        "cas",
        "pubchem_cid_resolved",
        "pubchem_iupac_name",
        "pubchem_record_title",
        "pubchem_name",
        "pubchem_synonyms",
        "pubchem_unii_list",
        "rxnorm_generic_name",
    ]
    _write_csv(
        in_csv,
        headers,
        [
            {
                "file_num": "6710",
                "sdf_title": "fda_6710",
                "display_name": "C19H17N5O7S3",
                "cas": "104010-37-9",
                "pubchem_cid_resolved": "",
                "pubchem_iupac_name": "",
                "pubchem_record_title": "",
                "pubchem_name": "",
                "pubchem_synonyms": "",
                "pubchem_unii_list": "",
                "rxnorm_generic_name": "",
            },
            {
                "file_num": "7686",
                "sdf_title": "fda_7686",
                "display_name": "UNK_7686",
                "cas": "2209104-86-7",
                "pubchem_cid_resolved": "",
                "pubchem_iupac_name": "",
                "pubchem_record_title": "",
                "pubchem_name": "",
                "pubchem_synonyms": "",
                "pubchem_unii_list": "",
                "rxnorm_generic_name": "",
            },
        ],
    )

    def fake_get(url: str, timeout: int = 10):
        if "/compound/name/104010-37-9/cids/JSON" in url:
            return FakeResponse({"IdentifierList": {"CID": [104010]}})
        if "/compound/name/2209104-86-7/cids/JSON" in url:
            return FakeResponse({"IdentifierList": {"CID": [2209104]}})
        if "/pug_view/data/compound/104010/JSON/" in url:
            return FakeResponse(
                {"Record": {"RecordTitle": "ceftiofur sodium", "Section": []}}
            )
        if "/pug_view/data/compound/2209104/JSON/" in url:
            return FakeResponse(
                {"Record": {"RecordTitle": "crisugabalin besylate", "Section": []}}
            )
        if "/compound/cid/104010/property/Title,IUPACName/JSON" in url:
            return FakeResponse(
                {
                    "PropertyTable": {
                        "Properties": [
                            {
                                "CID": 104010,
                                "Title": "ceftiofur sodium",
                                "IUPACName": "",
                            }
                        ]
                    }
                }
            )
        if "/compound/cid/2209104/property/Title,IUPACName/JSON" in url:
            return FakeResponse(
                {
                    "PropertyTable": {
                        "Properties": [
                            {
                                "CID": 2209104,
                                "Title": "crisugabalin besylate",
                                "IUPACName": "",
                            }
                        ]
                    }
                }
            )
        if "/compound/cid/104010/synonyms/JSON" in url:
            return FakeResponse({"InformationList": {"Information": [{"Synonym": []}]}})
        if "/compound/cid/2209104/synonyms/JSON" in url:
            return FakeResponse({"InformationList": {"Information": [{"Synonym": []}]}})
        if "rxnav.nlm.nih.gov/REST/approximateTerm.json" in url:
            return FakeResponse({"approximateGroup": {"candidate": []}})
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        exit_code = _run_tool(in_csv, out_csv, cache_dir)
    assert exit_code == 0

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert rows[0]["display_name"] == "ceftiofur sodium"
    assert rows[0]["pubchem_cid_resolved"] == "104010"
    assert rows[1]["display_name"] == "crisugabalin besylate"
    assert rows[1]["pubchem_cid_resolved"] == "2209104"
