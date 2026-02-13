import csv
from pathlib import Path
from unittest import mock
from urllib.parse import parse_qs, urlparse

from chemdb.tools import fill_fda_mapping_names as tool

HEADER = [
    "scheme",
    "file_num",
    "stage",
    "path",
    "stage_rank",
    "remark_smiles",
    "remark_inchikey",
    "remark_name",
    "pdbqt_formula",
    "pdbqt_heavy_atoms",
    "sdf_index",
    "match_method",
    "sdf_title",
    "smiles",
    "inchikey",
    "smiles_neutral",
    "cas",
    "id",
    "sdf_formula",
    "sdf_heavy_atoms",
    "pubchem_name",
    "pubchem_cid_resolved",
    "pubchem_record_title",
    "generic_name",
    "brand_names",
    "pubchem_iupac_name",
    "pubchem_synonyms",
    "pubchem_unii_list",
    "display_name",
    "rxnorm_rxcui",
    "rxnorm_generic_name",
    "rxnorm_brand_names",
    "drugcentral_id",
    "drugcentral_generic_name",
    "drugcentral_brand_names",
]


class FakeResponse:
    def __init__(self, payload, status_code=200):
        self._payload = payload
        self.status_code = status_code

    def json(self):
        return self._payload


def _write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER)
        writer.writeheader()
        writer.writerows(rows)


def _make_row(**overrides: str) -> dict[str, str]:
    row = {key: "" for key in HEADER}
    row.update(overrides)
    return row


def _fake_get(url: str, timeout: int = 10):
    if "pubchem.ncbi.nlm.nih.gov" in url and "/inchikey/" in url:
        inchikey = url.split("/inchikey/")[1].split("/")[0]
        if inchikey == "INCHIKEYA":
            return FakeResponse({"IdentifierList": {"CID": [111]}})
        if inchikey == "INCHIKEYB":
            return FakeResponse({"IdentifierList": {"CID": [222]}})
        return FakeResponse({"IdentifierList": {"CID": []}})
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
                                "IUPACName": (
                                    "3-[[2'-(1H-tetrazol-5-yl)[1,1'-biphenyl]-4-yl]"
                                    "methyl]-2-ethoxy-1,3-diazaspiro[4.4]non-1-en-4-one"
                                ),
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
                        "Information": [
                            {"CID": 222, "Synonym": ["Candesartan", "AT1 antagonist"]}
                        ]
                    }
                }
            )
    if "rxnav.nlm.nih.gov/REST/approximateTerm.json" in url:
        term = parse_qs(urlparse(url).query).get("term", [""])[0]
        if "Metformin" in term:
            return FakeResponse(
                {
                    "approximateGroup": {
                        "candidate": [{"rxcui": "860975", "score": "100"}]
                    }
                }
            )
        if "Candesartan" in term:
            return FakeResponse(
                {"approximateGroup": {"candidate": [{"rxcui": "1816", "score": "100"}]}}
            )
        return FakeResponse({"approximateGroup": {"candidate": []}})
    if "rxnav.nlm.nih.gov/REST/rxcui/860975/properties.json" in url:
        return FakeResponse({"properties": {"name": "Metformin", "rxcui": "860975"}})
    if "rxnav.nlm.nih.gov/REST/rxcui/1816/properties.json" in url:
        return FakeResponse({"properties": {"name": "Candesartan", "rxcui": "1816"}})
    return FakeResponse({}, status_code=404)


def test_fill_fda_mapping_names(tmp_path):
    in_csv = tmp_path / "fda_mapping.csv"
    out_csv = tmp_path / "fda_mapping_filled.csv"
    cache_dir = tmp_path / "cache"

    rows = [
        _make_row(
            scheme="rdk",
            file_num="1",
            path="/x/rdk_0000001.pdbqt",
            sdf_title="RowA",
            inchikey="INCHIKEYA",
            display_name="fda_1234",
        ),
        _make_row(
            scheme="rdk",
            file_num="2",
            path="/x/rdk_0000002.pdbqt",
            sdf_title="RowB",
            inchikey="INCHIKEYB",
            display_name=(
                "6,7-dimethoxy-2-[4-(phenylmethyl)piperazin-1-yl]" "quinazolin-4-amine"
            ),
        ),
        _make_row(
            scheme="rdk",
            file_num="3",
            path="/x/rdk_0000003.pdbqt",
            sdf_title="RowC",
            display_name="",
        ),
    ]
    _write_csv(in_csv, rows)
    original_bytes = in_csv.read_bytes()

    with mock.patch.object(tool.requests, "get", side_effect=_fake_get):
        exit_code = tool.main(
            [
                "--in_csv",
                str(in_csv),
                "--out_csv",
                str(out_csv),
                "--sleep",
                "0",
                "--cache_dir",
                str(cache_dir),
                "--max_workers",
                "4",
                "--fail_if_unresolved",
            ]
        )

    assert exit_code == 2
    assert in_csv.read_bytes() == original_bytes

    with in_csv.open("r", encoding="utf-8") as handle:
        in_header = handle.readline().strip()
    with out_csv.open("r", encoding="utf-8") as handle:
        out_header = handle.readline().strip()
    assert out_header == in_header

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        out_rows = list(reader)

    by_file_num = {row["file_num"]: row for row in out_rows}
    assert by_file_num["1"]["display_name"] == "Metformin"
    assert by_file_num["2"]["display_name"] == "Candesartan"
    assert not tool.is_bad_display_name(
        by_file_num["2"]["display_name"], by_file_num["2"].get("pubchem_iupac_name")
    )
    assert by_file_num["3"]["display_name"] == "UNK_3"
    assert tool.is_bad_display_name(
        by_file_num["3"]["display_name"], by_file_num["3"].get("pubchem_iupac_name")
    )
