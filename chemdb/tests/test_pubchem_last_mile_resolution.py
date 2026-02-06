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


def test_pubchem_pugview_recordtitle_resolves_and_counts(tmp_path, capsys):
    in_csv = tmp_path / "in.csv"
    out_csv = tmp_path / "out.csv"
    cache_dir = tmp_path / "cache"
    headers = [
        "file_num",
        "display_name",
        "pubchem_cid_resolved",
        "pubchem_iupac_name",
        "pubchem_record_title",
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
                "pubchem_synonyms": "",
                "rxnorm_generic_name": "",
            }
        ],
    )

    def fake_get(url: str, timeout: int = 10):
        if "/cid/111/property/Title,IUPACName/JSON" in url:
            return FakeResponse(
                {
                    "PropertyTable": {
                        "Properties": [{"CID": 111, "Title": "", "IUPACName": "2-acetoxybenzoic acid"}]
                    }
                }
            )
        if "/cid/111/synonyms/JSON" in url:
            return FakeResponse({}, status_code=404)
        if "/pug_view/data/compound/111/JSON/" in url:
            return FakeResponse({"Record": {"RecordTitle": "Aspirin", "Section": []}})
        if "rxnav.nlm.nih.gov/REST/approximateTerm.json" in url:
            return FakeResponse({"approximateGroup": {"candidate": []}})
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        exit_code = _run_tool(in_csv, out_csv, cache_dir)
    assert exit_code == 0

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["display_name"] == "Aspirin"

    output = capsys.readouterr().out
    assert "pubchem pugview title resolved: 1" in output
    assert "pubchem pugview calls: 1" in output


def test_pubchem_synonyms_404_but_pugview_still_resolves(tmp_path):
    in_csv = tmp_path / "in.csv"
    out_csv = tmp_path / "out.csv"
    cache_dir = tmp_path / "cache"
    headers = [
        "file_num",
        "display_name",
        "pubchem_cid_resolved",
        "pubchem_iupac_name",
        "pubchem_record_title",
        "pubchem_synonyms",
        "rxnorm_generic_name",
    ]
    _write_csv(
        in_csv,
        headers,
        [
            {
                "file_num": "2",
                "display_name": "UNK_2",
                "pubchem_cid_resolved": "222",
                "pubchem_iupac_name": "",
                "pubchem_record_title": "",
                "pubchem_synonyms": "",
                "rxnorm_generic_name": "",
            }
        ],
    )

    def fake_get(url: str, timeout: int = 10):
        if "/cid/222/property/Title,IUPACName/JSON" in url:
            return FakeResponse(
                {
                    "PropertyTable": {
                        "Properties": [{"CID": 222, "Title": "", "IUPACName": "2,3-dimethyl-1H"}]
                    }
                }
            )
        if "/cid/222/synonyms/JSON" in url:
            return FakeResponse({}, status_code=404)
        if "/pug_view/data/compound/222/JSON/" in url:
            return FakeResponse({"Record": {"RecordTitle": "Candesartan", "Section": []}})
        if "rxnav.nlm.nih.gov/REST/approximateTerm.json" in url:
            return FakeResponse({"approximateGroup": {"candidate": []}})
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        exit_code = _run_tool(in_csv, out_csv, cache_dir)
    assert exit_code == 0

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["display_name"] == "Candesartan"


def test_pubchem_smiles_400_uses_inchi_cid_fallback(tmp_path, capsys):
    in_csv = tmp_path / "in.csv"
    out_csv = tmp_path / "out.csv"
    cache_dir = tmp_path / "cache"
    headers = [
        "file_num",
        "display_name",
        "smiles",
        "inchi",
        "pubchem_cid_resolved",
        "pubchem_iupac_name",
        "pubchem_record_title",
        "pubchem_synonyms",
        "rxnorm_generic_name",
    ]
    _write_csv(
        in_csv,
        headers,
        [
            {
                "file_num": "3",
                "display_name": "UNK_3",
                "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
                "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
                "pubchem_cid_resolved": "",
                "pubchem_iupac_name": "",
                "pubchem_record_title": "",
                "pubchem_synonyms": "",
                "rxnorm_generic_name": "",
            }
        ],
    )

    def fake_get(url: str, timeout: int = 10):
        if "/compound/smiles/" in url and "/cids/JSON" in url:
            return FakeResponse({}, status_code=400)
        if "/compound/inchi/" in url and "/cids/JSON" in url:
            return FakeResponse({"IdentifierList": {"CID": [555]}})
        if "/cid/555/property/Title,IUPACName/JSON" in url:
            return FakeResponse(
                {"PropertyTable": {"Properties": [{"CID": 555, "Title": "", "IUPACName": "3,4-dimethyl"}]}}
            )
        if "/cid/555/synonyms/JSON" in url:
            return FakeResponse({}, status_code=404)
        if "/pug_view/data/compound/555/JSON/" in url:
            return FakeResponse({"Record": {"RecordTitle": "Atenolol", "Section": []}})
        if "rxnav.nlm.nih.gov/REST/approximateTerm.json" in url:
            return FakeResponse({"approximateGroup": {"candidate": []}})
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        exit_code = _run_tool(in_csv, out_csv, cache_dir)
    assert exit_code == 0

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["display_name"] == "Atenolol"
    assert rows[0]["pubchem_cid_resolved"] == "555"

    output = capsys.readouterr().out
    assert "pubchem smiles 400 inchi fallback attempted: 1" in output
    assert "pubchem inchi cid resolved: 1" in output


def test_pubchem_parent_connectivity_rescue(tmp_path):
    in_csv = tmp_path / "in.csv"
    out_csv = tmp_path / "out.csv"
    cache_dir = tmp_path / "cache"
    headers = [
        "file_num",
        "display_name",
        "pubchem_cid_resolved",
        "pubchem_iupac_name",
        "pubchem_record_title",
        "pubchem_parent_cid",
        "pubchem_synonyms",
        "rxnorm_generic_name",
    ]
    _write_csv(
        in_csv,
        headers,
        [
            {
                "file_num": "4",
                "display_name": "UNK_4",
                "pubchem_cid_resolved": "777",
                "pubchem_iupac_name": "",
                "pubchem_record_title": "",
                "pubchem_parent_cid": "",
                "pubchem_synonyms": "",
                "rxnorm_generic_name": "",
            }
        ],
    )

    def fake_get(url: str, timeout: int = 10):
        if "/cid/777/property/Title,IUPACName/JSON" in url:
            return FakeResponse(
                {"PropertyTable": {"Properties": [{"CID": 777, "Title": "", "IUPACName": "2,4-dimethyl"}]}}
            )
        if "/cid/777/synonyms/JSON" in url:
            return FakeResponse({}, status_code=404)
        if "fastidentity/cid/777/cids/JSON?identity_type=same_parent_connectivity" in url:
            return FakeResponse({"IdentifierList": {"CID": [888]}})
        if "/pug_view/data/compound/777/JSON/" in url:
            return FakeResponse({"Record": {"RecordTitle": "2,3-dimethyl-1H", "Section": []}})
        if "/pug_view/data/compound/888/JSON/" in url:
            return FakeResponse({"Record": {"RecordTitle": "Clopidogrel", "Section": []}})
        if "rxnav.nlm.nih.gov/REST/approximateTerm.json" in url:
            return FakeResponse({"approximateGroup": {"candidate": []}})
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        exit_code = _run_tool(in_csv, out_csv, cache_dir)
    assert exit_code == 0

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["display_name"] == "Clopidogrel"
    assert rows[0]["pubchem_parent_cid"] == "888"


def test_pubchem_inchikey_2d_cid_fallback(tmp_path):
    in_csv = tmp_path / "in.csv"
    out_csv = tmp_path / "out.csv"
    cache_dir = tmp_path / "cache"
    headers = [
        "file_num",
        "display_name",
        "inchikey",
        "pubchem_cid_resolved",
        "pubchem_iupac_name",
        "pubchem_record_title",
        "pubchem_synonyms",
        "rxnorm_generic_name",
    ]
    _write_csv(
        in_csv,
        headers,
        [
            {
                "file_num": "5",
                "display_name": "UNK_5",
                "inchikey": "APRTUXPAFSNIOL-MSOLQXFVSA-O",
                "pubchem_cid_resolved": "",
                "pubchem_iupac_name": "",
                "pubchem_record_title": "",
                "pubchem_synonyms": "",
                "rxnorm_generic_name": "",
            }
        ],
    )

    def fake_get(url: str, timeout: int = 10):
        if "/inchikey/APRTUXPAFSNIOL-MSOLQXFVSA-O/cids/JSON" in url:
            return FakeResponse({}, status_code=404)
        if "/inchikey/APRTUXPAFSNIOL/cids/JSON" in url:
            return FakeResponse({"IdentifierList": {"CID": [999]}})
        if "/cid/999/property/Title,IUPACName/JSON" in url:
            return FakeResponse(
                {"PropertyTable": {"Properties": [{"CID": 999, "Title": "", "IUPACName": "1,2-dimethyl"}]}}
            )
        if "/cid/999/synonyms/JSON" in url:
            return FakeResponse({}, status_code=404)
        if "/pug_view/data/compound/999/JSON/" in url:
            return FakeResponse({"Record": {"RecordTitle": "Bisoprolol", "Section": []}})
        if "rxnav.nlm.nih.gov/REST/approximateTerm.json" in url:
            return FakeResponse({"approximateGroup": {"candidate": []}})
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        exit_code = _run_tool(in_csv, out_csv, cache_dir)
    assert exit_code == 0

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["display_name"] == "Bisoprolol"
    assert rows[0]["pubchem_cid_resolved"] == "999"


def test_pubchem_inchikey_fastidentity_cid_fallback(tmp_path):
    in_csv = tmp_path / "in.csv"
    out_csv = tmp_path / "out.csv"
    cache_dir = tmp_path / "cache"
    headers = [
        "file_num",
        "display_name",
        "inchikey",
        "pubchem_cid_resolved",
        "pubchem_iupac_name",
        "pubchem_record_title",
        "pubchem_synonyms",
        "rxnorm_generic_name",
    ]
    _write_csv(
        in_csv,
        headers,
        [
            {
                "file_num": "6",
                "display_name": "UNK_6",
                "inchikey": "ZZZZZZZZZZZZZZ-AAAAAAAAAA-Z",
                "pubchem_cid_resolved": "",
                "pubchem_iupac_name": "",
                "pubchem_record_title": "",
                "pubchem_synonyms": "",
                "rxnorm_generic_name": "",
            }
        ],
    )

    def fake_get(url: str, timeout: int = 10):
        if (
            "fastidentity/inchikey/ZZZZZZZZZZZZZZ-AAAAAAAAAA-Z/cids/JSON"
            "?identity_type=same_parent_connectivity" in url
        ):
            return FakeResponse({"IdentifierList": {"CID": [1001]}})
        if "/rest/pug/compound/inchikey/ZZZZZZZZZZZZZZ-AAAAAAAAAA-Z/cids/JSON" in url:
            return FakeResponse({}, status_code=404)
        if "/rest/pug/compound/inchikey/ZZZZZZZZZZZZZZ/cids/JSON" in url:
            return FakeResponse({}, status_code=404)
        if "/cid/1001/property/Title,IUPACName/JSON" in url:
            return FakeResponse(
                {"PropertyTable": {"Properties": [{"CID": 1001, "Title": "", "IUPACName": "2,4-dimethyl"}]}}
            )
        if "/cid/1001/synonyms/JSON" in url:
            return FakeResponse({}, status_code=404)
        if "/pug_view/data/compound/1001/JSON/" in url:
            return FakeResponse({"Record": {"RecordTitle": "Warfarin", "Section": []}})
        if "rxnav.nlm.nih.gov/REST/approximateTerm.json" in url:
            return FakeResponse({"approximateGroup": {"candidate": []}})
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        exit_code = _run_tool(in_csv, out_csv, cache_dir)
    assert exit_code == 0

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["display_name"] == "Warfarin"
    assert rows[0]["pubchem_cid_resolved"] == "1001"
