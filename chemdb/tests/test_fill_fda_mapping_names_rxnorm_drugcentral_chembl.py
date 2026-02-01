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


def _write_tsv(path: Path, headers: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_rxnorm_promotion_sets_display_name(tmp_path):
    in_csv = tmp_path / "fda_mapping.csv"
    out_csv = tmp_path / "fda_mapping_filled.csv"
    cache_dir = tmp_path / "cache"

    headers = [
        "file_num",
        "display_name",
        "rxnorm_generic_name",
        "pubchem_iupac_name",
    ]
    rows = [
        {
            "file_num": "1",
            "display_name": "UNK_1234",
            "rxnorm_generic_name": "Metformin",
            "pubchem_iupac_name": "",
        }
    ]
    _write_csv(in_csv, headers, rows)

    exit_code = tool.main(
        [
            "--in_csv",
            str(in_csv),
            "--out_csv",
            str(out_csv),
            "--enable_unichem",
            "false",
            "--enable_pubchem_rerank",
            "false",
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

    assert out_rows[0]["display_name"] == "Metformin"


def test_drugcentral_smiles_inchikey_derivation(tmp_path, monkeypatch):
    in_csv = tmp_path / "fda_mapping.csv"
    out_csv = tmp_path / "fda_mapping_filled.csv"
    cache_dir = tmp_path / "cache"
    drugcentral_tsv = tmp_path / "drugcentral.tsv"

    headers = [
        "file_num",
        "display_name",
        "inchikey",
        "generic_name",
        "pubchem_iupac_name",
    ]
    rows = [
        {
            "file_num": "1",
            "display_name": "UNK_1",
            "inchikey": "TESTINCHIKEY1234567890ABCDEF",
            "generic_name": "",
            "pubchem_iupac_name": "",
        }
    ]
    _write_csv(in_csv, headers, rows)

    _write_tsv(
        drugcentral_tsv,
        ["INN", "SMILES"],
        [
            {
                "INN": "Metformin",
                "SMILES": "CN(C)NC(=N)N=C(N)N",
            }
        ],
    )

    monkeypatch.setattr(tool, "_rdkit_available", lambda: True)
    monkeypatch.setattr(
        tool, "rdkit_inchikey_from_smiles", lambda value: "TESTINCHIKEY1234567890ABCDEF"
    )

    exit_code = tool.main(
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
            "false",
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

    assert out_rows[0]["display_name"] == "Metformin"
    assert out_rows[0]["generic_name"] == "Metformin"


def _fake_post(url: str, json=None, timeout: int = 10):
    if "ebi.ac.uk/unichem/api/v1/compounds" in url:
        return FakeResponse(
            {
                "compounds": [
                    {
                        "sources": [
                            {"id": 1, "compoundId": "CHEMBL_FIRST_APPROVAL"},
                        ]
                    }
                ]
            }
        )
    return FakeResponse({}, status_code=404)


def _fake_get(url: str, timeout: int = 10):
    if "ebi.ac.uk/chembl/api/data/molecule/CHEMBL_FIRST_APPROVAL.json" in url:
        return FakeResponse({"pref_name": "Candesartan", "first_approval": 2010})
    return FakeResponse({}, status_code=404)


def test_chembl_first_approval_acceptance(tmp_path):
    in_csv = tmp_path / "fda_mapping.csv"
    out_csv = tmp_path / "fda_mapping_filled.csv"
    cache_dir = tmp_path / "cache"

    headers = [
        "file_num",
        "display_name",
        "inchikey",
        "pubchem_iupac_name",
    ]
    rows = [
        {
            "file_num": "1",
            "display_name": "UNK_1",
            "inchikey": "INCHIKEY_CHEMBL",
            "pubchem_iupac_name": "",
        }
    ]
    _write_csv(in_csv, headers, rows)

    with mock.patch.object(
        tool.requests, "get", side_effect=_fake_get
    ), mock.patch.object(tool.requests, "post", side_effect=_fake_post):
        exit_code = tool.main(
            [
                "--in_csv",
                str(in_csv),
                "--out_csv",
                str(out_csv),
                "--enable_unichem",
                "true",
                "--enable_pubchem_rerank",
                "false",
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

    assert out_rows[0]["display_name"] == "Candesartan"
