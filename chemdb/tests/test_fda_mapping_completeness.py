import csv
import os
from pathlib import Path

import pytest

from chemdb.tools import fill_fda_mapping_names as tool


def test_fda_mapping_completeness(tmp_path):
    default_path = Path("chemdb/data/fda_mapping_from_pdbqt.csv")
    in_csv = Path(os.environ.get("FDA_MAPPING_IN_CSV", default_path))
    drugcentral_env = os.environ.get("DRUGCENTRAL_STRUCTURES_TSV", "")
    drugcentral_path = Path(drugcentral_env) if drugcentral_env else None

    if not in_csv.exists():
        pytest.skip(f"FDA mapping CSV not found: {in_csv}")
    if drugcentral_path is None or not drugcentral_path.exists():
        pytest.skip("DRUGCENTRAL_STRUCTURES_TSV not set or file missing.")

    out_csv = tmp_path / "fda_mapping_filled.csv"
    report_path = tmp_path / "fda_mapping.unresolved.txt"

    exit_code = tool.main(
        [
            "--in_csv",
            str(in_csv),
            "--out_csv",
            str(out_csv),
            "--report_path",
            str(report_path),
            "--drugcentral_structures_tsv",
            str(drugcentral_path),
            "--enable_unichem",
            "false",
            "--enable_pubchem_rerank",
            "false",
            "--require_all_named",
            "--sleep",
            "0",
            "--cache_mode",
            "off",
            "--max_workers",
            "4",
        ]
    )

    assert exit_code == 0

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            assert not tool.is_bad_display_name(
                row.get("display_name", ""), row.get("pubchem_iupac_name", "")
            )
