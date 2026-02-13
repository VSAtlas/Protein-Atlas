import csv
from pathlib import Path
from unittest import mock

from chemdb.tools import fill_fda_mapping_names as tool

HEADER = [
    "file_num",
    "display_name",
    "inchikey",
    "pubchem_cid_resolved",
    "pubchem_iupac_name",
    "generic_name",
]


class FakeResponse:
    def __init__(self, payload, status_code=200):
        self._payload = payload
        self.status_code = status_code
        self.headers = {}

    def json(self):
        return self._payload


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER)
        writer.writeheader()
        writer.writerows(rows)


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["INN", "INCHIKEY"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _fake_post(url: str, json=None, timeout: int = 10):
    if "ebi.ac.uk/unichem/api/v1/compounds" in url:
        inchikey = (json or {}).get("compound", "")
        if inchikey == "INCHIKEY_MET":
            return FakeResponse([{"src_id": 1, "src_compound_id": "CHEMBL_MET"}])
        if inchikey == "INCHIKEY_CHM":
            return FakeResponse([{"src_id": 1, "src_compound_id": "CHEMBL_CHM"}])
        return FakeResponse([])
    return FakeResponse({}, status_code=404)


def _fake_get(url: str, timeout: int = 10):
    if "ebi.ac.uk/chembl/api/data/molecule/CHEMBL_MET.json" in url:
        return FakeResponse({"pref_name": "WrongName", "max_phase": 4})
    if "ebi.ac.uk/chembl/api/data/molecule/CHEMBL_CHM.json" in url:
        return FakeResponse({"pref_name": "Candesartan", "max_phase": "4"})
    if "pubchem.ncbi.nlm.nih.gov" in url and "/property/Title,IUPACName/JSON" in url:
        if "/cid/111/" in url:
            return FakeResponse(
                {
                    "PropertyTable": {
                        "Properties": [
                            {
                                "CID": 111,
                                "Title": "Metformin",
                                "IUPACName": "N,N-dimethylbiguanide",
                            }
                        ]
                    }
                }
            )
        if "/cid/222/" in url:
            return FakeResponse(
                {
                    "PropertyTable": {
                        "Properties": [
                            {
                                "CID": 222,
                                "Title": "Candesartan",
                                "IUPACName": "3-(benzyl)-1-oxa-2-azabicyclo",
                            }
                        ]
                    }
                }
            )
        if "/cid/333/" in url:
            return FakeResponse(
                {
                    "PropertyTable": {
                        "Properties": [
                            {
                                "CID": 333,
                                "Title": "3-(1H-imidazol-1-yl)propanoic acid",
                                "IUPACName": "3-(1H-imidazol-1-yl)propanoic acid",
                            }
                        ]
                    }
                }
            )
    if "pubchem.ncbi.nlm.nih.gov" in url and "/synonyms/JSON" in url:
        if "/cid/111/" in url:
            return FakeResponse(
                {
                    "InformationList": {
                        "Information": [{"CID": 111, "Synonym": ["Metformin"]}]
                    }
                }
            )
        if "/cid/222/" in url:
            return FakeResponse(
                {
                    "InformationList": {
                        "Information": [{"CID": 222, "Synonym": ["Candesartan"]}]
                    }
                }
            )
        if "/cid/333/" in url:
            return FakeResponse(
                {
                    "InformationList": {
                        "Information": [
                            {
                                "CID": 333,
                                "Synonym": [
                                    "3-(1H-imidazol-1-yl)propanoic acid",
                                    "Aspirin",
                                    "ABC-123",
                                ],
                            }
                        ]
                    }
                }
            )
    if "rxnav.nlm.nih.gov/REST/approximateTerm.json" in url:
        return FakeResponse({"approximateGroup": {"candidate": []}})
    if "rxnav.nlm.nih.gov/REST/rxcui/" in url:
        return FakeResponse({"properties": {}})
    return FakeResponse({}, status_code=404)


def test_resolution_sources_order(tmp_path):
    in_csv = tmp_path / "fda_mapping.csv"
    out_csv = tmp_path / "fda_mapping_filled.csv"
    cache_dir = tmp_path / "cache"
    drugcentral_tsv = tmp_path / "drugcentral.tsv"

    _write_tsv(
        drugcentral_tsv,
        [
            {
                "INN": "Metformin",
                "INCHIKEY": "INCHIKEY_MET",
            }
        ],
    )

    rows = [
        {
            "file_num": "1",
            "display_name": "fda_0001",
            "inchikey": "INCHIKEY_MET",
            "pubchem_cid_resolved": "111",
            "pubchem_iupac_name": "",
            "generic_name": "",
        },
        {
            "file_num": "2",
            "display_name": "UNK_2",
            "inchikey": "INCHIKEY_CHM",
            "pubchem_cid_resolved": "222",
            "pubchem_iupac_name": "",
            "generic_name": "",
        },
        {
            "file_num": "3",
            "display_name": "",
            "inchikey": "INCHIKEY_PUB",
            "pubchem_cid_resolved": "333",
            "pubchem_iupac_name": "",
            "generic_name": "",
        },
    ]
    _write_csv(in_csv, rows)

    with mock.patch.object(
        tool.requests, "get", side_effect=_fake_get
    ), mock.patch.object(tool.requests, "post", side_effect=_fake_post):
        exit_code = tool.main(
            [
                "--in_csv",
                str(in_csv),
                "--out_csv",
                str(out_csv),
                "--drugcentral_structures_tsv",
                str(drugcentral_tsv),
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
        reader = csv.DictReader(handle)
        out_rows = list(reader)

    by_file_num = {row["file_num"]: row for row in out_rows}
    assert by_file_num["1"]["display_name"] == "Metformin"
    assert by_file_num["1"]["generic_name"] == "Metformin"
    assert by_file_num["2"]["display_name"] == "Candesartan"
    assert by_file_num["3"]["display_name"] == "Aspirin"
    assert not tool.is_bad_display_name(
        by_file_num["1"]["display_name"], by_file_num["1"].get("pubchem_iupac_name")
    )
    assert not tool.is_bad_display_name(
        by_file_num["2"]["display_name"], by_file_num["2"].get("pubchem_iupac_name")
    )
    assert not tool.is_bad_display_name(
        by_file_num["3"]["display_name"], by_file_num["3"].get("pubchem_iupac_name")
    )
