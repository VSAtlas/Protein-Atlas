import csv
import logging
from pathlib import Path

from chemdb.tools import fill_fda_mapping_names as tool


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


def test_drugcentral_2d_key_matches_full_inchikey_row(tmp_path):
    in_csv = tmp_path / "mapping.csv"
    out_csv = tmp_path / "mapping_out.csv"
    cache_dir = tmp_path / "cache"
    drugcentral_tsv = tmp_path / "drugcentral.tsv"

    headers = [
        "file_num",
        "display_name",
        "inchikey",
        "generic_name",
        "pubchem_iupac_name",
    ]
    _write_csv(
        in_csv,
        headers,
        [
            {
                "file_num": "1",
                "display_name": "UNK_1",
                "inchikey": "ABCDEFGHIJKLMN-QRSTUVWXYY-Z",
                "generic_name": "",
                "pubchem_iupac_name": "",
            }
        ],
    )
    _write_tsv(
        drugcentral_tsv,
        ["INN", "inchi_key_2d"],
        [{"INN": "Metformin", "inchi_key_2d": "ABCDEFGHIJKLMN"}],
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
        rows = list(csv.DictReader(handle))
    assert rows[0]["display_name"] == "Metformin"
    assert rows[0]["generic_name"] == "Metformin"


def test_drugcentral_rdkit_gating_without_inchikey_columns(tmp_path, monkeypatch, caplog):
    drugcentral_tsv = tmp_path / "drugcentral.tsv"
    _write_tsv(
        drugcentral_tsv,
        ["INN", "SMILES"],
        [{"INN": "Metformin", "SMILES": "CN(C)NC(=N)N=C(N)N"}],
    )

    monkeypatch.setattr(tool, "_rdkit_available", lambda: False)
    caplog.set_level(logging.WARNING)

    index = tool.load_drugcentral_index(drugcentral_tsv)

    assert len(index.full) == 0
    assert len(index.two_d) == 0
    warnings = [
        rec.message
        for rec in caplog.records
        if "DrugCentral TSV missing INCHIKEY column" in rec.message
    ]
    assert len(warnings) == 1
