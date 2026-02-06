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


def _fake_get(url: str, timeout: int = 10):
    if "/cid/200/property/Title,IUPACName/JSON" in url:
        return FakeResponse(
            {
                "PropertyTable": {
                    "Properties": [
                        {
                            "CID": 200,
                            "Title": "C17H12Cl2F3N7O2S",
                            "IUPACName": "N,N-bis(2-chloroethyl)-1,3,2-oxazaphosphinan-2-amine",
                        }
                    ]
                }
            }
        )
    if "/cid/201/property/Title,IUPACName/JSON" in url:
        return FakeResponse(
            {
                "PropertyTable": {
                    "Properties": [
                        {
                            "CID": 201,
                            "Title": "3 hydroxy 2 phenyl 1h inden 1 one",
                            "IUPACName": "3-hydroxy-2-phenyl-1H-inden-1-one",
                        }
                    ]
                }
            }
        )
    if "/cid/200/synonyms/JSON" in url:
        return FakeResponse(
            {"InformationList": {"Information": [{"CID": 200, "Synonym": ["Cinnarizine"]}]}}
        )
    if "/cid/201/synonyms/JSON" in url:
        return FakeResponse(
            {"InformationList": {"Information": [{"CID": 201, "Synonym": ["Phenindione"]}]}}
        )
    if "rxnav.nlm.nih.gov/REST/approximateTerm.json" in url:
        return FakeResponse({"approximateGroup": {"candidate": []}})
    if "rxnav.nlm.nih.gov/REST/rxcui/" in url:
        return FakeResponse({"properties": {}})
    return FakeResponse({}, status_code=404)


def test_quality_fix_mode_strict_updates_only_flagged_rows(tmp_path):
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
        "generic_name",
        "rxnorm_generic_name",
        "rxnorm_rxcui",
    ]
    _write_csv(
        in_csv,
        headers,
        [
            {
                "file_num": "1",
                "display_name": "Aspirin",
                "pubchem_cid_resolved": "",
                "pubchem_iupac_name": "",
                "pubchem_record_title": "",
                "pubchem_name": "",
                "pubchem_synonyms": "",
                "generic_name": "",
                "rxnorm_generic_name": "",
                "rxnorm_rxcui": "",
            },
            {
                "file_num": "2",
                "display_name": "C17H12Cl2F3N7O2S",
                "pubchem_cid_resolved": "200",
                "pubchem_iupac_name": "",
                "pubchem_record_title": "",
                "pubchem_name": "",
                "pubchem_synonyms": "",
                "generic_name": "",
                "rxnorm_generic_name": "",
                "rxnorm_rxcui": "",
            },
            {
                "file_num": "3",
                "display_name": "Metformin",
                "pubchem_cid_resolved": "",
                "pubchem_iupac_name": "",
                "pubchem_record_title": "",
                "pubchem_name": "",
                "pubchem_synonyms": "",
                "generic_name": "",
                "rxnorm_generic_name": "",
                "rxnorm_rxcui": "",
            },
            {
                "file_num": "4",
                "display_name": "3 hydroxy 2 phenyl 1h inden 1 one",
                "pubchem_cid_resolved": "201",
                "pubchem_iupac_name": "",
                "pubchem_record_title": "",
                "pubchem_name": "",
                "pubchem_synonyms": "",
                "generic_name": "",
                "rxnorm_generic_name": "",
                "rxnorm_rxcui": "",
            },
            {
                "file_num": "5",
                "display_name": "Ibuprofen",
                "pubchem_cid_resolved": "",
                "pubchem_iupac_name": "",
                "pubchem_record_title": "",
                "pubchem_name": "",
                "pubchem_synonyms": "",
                "generic_name": "",
                "rxnorm_generic_name": "",
                "rxnorm_rxcui": "",
            },
        ],
    )

    with mock.patch.object(tool.requests, "get", side_effect=_fake_get):
        exit_code = tool.main(
            [
                "--in_csv",
                str(in_csv),
                "--out_csv",
                str(out_csv),
                "--quality_fix_mode",
                "strict",
                "--enable_unichem",
                "false",
                "--enable_gsrs",
                "false",
                "--sleep",
                "0",
                "--cache_dir",
                str(cache_dir),
                "--max_workers",
                "1",
            ]
        )
    assert exit_code == 0

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        out_rows = {row["file_num"]: row for row in csv.DictReader(handle)}

    assert out_rows["1"]["display_name"] == "Aspirin"
    assert out_rows["3"]["display_name"] == "Metformin"
    assert out_rows["5"]["display_name"] == "Ibuprofen"
    assert out_rows["2"]["display_name"] == "Cinnarizine"
    assert out_rows["4"]["display_name"] == "Phenindione"
