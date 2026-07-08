from __future__ import annotations

import csv
from pathlib import Path

import pytest

from docking.docking_consensus_score import compute_consensus_for_variant_ph

REPO_ROOT = Path(__file__).resolve().parents[2]
PREPPED_TEST_LIB = (
    REPO_ROOT
    / "chemdb"
    / "tests"
    / "fixtures"
    / "prepped_ligands"
    / "fda_test_library_10"
)


class DummyPaths:
    def __init__(self, pdb_id: str, variant_root: Path) -> None:
        self.pdb_id = pdb_id
        self._variant_root = variant_root

    def docked_variant_root(
        self, variant_env: str | None, ph_label: str | None
    ) -> Path:
        return self._variant_root


def _load_test_library_bases() -> list[str]:
    if not PREPPED_TEST_LIB.exists():
        return []
    ph_dir = PREPPED_TEST_LIB / "pH7_0"
    bases = []
    if ph_dir.exists():
        bases = sorted({path.stem for path in ph_dir.glob("*.pdbqt")})
    if not bases:
        bases = sorted({path.stem for path in PREPPED_TEST_LIB.glob("*.pdbqt")})
    if len(bases) < 10:
        bases = sorted({path.stem for path in PREPPED_TEST_LIB.rglob("*.pdbqt")})
    return bases


def _write_minimal_score_csv(path: Path, ligands: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["ligand", "score", "valid"])
        writer.writeheader()
        for idx, ligand in enumerate(ligands, start=1):
            writer.writerow(
                {"ligand": ligand, "score": f"{-1.0 * idx:.2f}", "valid": "true"}
            )


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def _assert_consensus_desc(rows: list[dict[str, str]]) -> None:
    scores = [float(row["consensus_score"]) for row in rows]
    assert scores == sorted(scores, reverse=True)


def test_consensus_prefix_outputs_and_no_overwrite(tmp_path: Path) -> None:
    bases = _load_test_library_bases()
    if len(bases) < 10:
        pytest.skip("bundled prepped ligand fixture missing or incomplete.")

    bases = bases[:10]
    base_ligands = [f"{base}.pdbqt" for base in bases]
    prefixed_ligands = [f"fda_dud_{base}.pdbqt" for base in bases]

    run_id = "unit_test_run"
    pdb_id = "TEST"
    variant_root = tmp_path / "docked" / run_id / pdb_id / "HOLO" / "pH7_0"
    paths = DummyPaths(pdb_id=pdb_id, variant_root=variant_root)
    cfg = {"RUN_ID": run_id}

    _write_minimal_score_csv(variant_root / "docking_score_long.csv", base_ligands)
    _write_minimal_score_csv(
        variant_root / "fda_dud_docking_score_long.csv", prefixed_ligands
    )
    _write_minimal_score_csv(
        variant_root / "dud_docking_score_long.csv",
        [f"dud_{base}.pdbqt" for base in bases],
    )

    compute_consensus_for_variant_ph(
        cfg=cfg,
        paths=paths,
        ph_label="pH7_0",
        variant_env="HOLO",
        variant_label="HOLO",
        csv_prefix="",
        run_mode="fda",
    )
    compute_consensus_for_variant_ph(
        cfg=cfg,
        paths=paths,
        ph_label="pH7_0",
        variant_env="HOLO",
        variant_label="HOLO",
        csv_prefix="fda_dud_",
        run_mode="fda_dud",
    )
    compute_consensus_for_variant_ph(
        cfg=cfg,
        paths=paths,
        ph_label="pH7_0",
        variant_env="HOLO",
        variant_label="HOLO",
        csv_prefix="dud_",
        run_mode="dud",
    )

    base_out = variant_root / "consensus_docking_scores.csv"
    prefixed_out = variant_root / "fda_dud_consensus_docking_scores.csv"
    dud_out = variant_root / "dud_consensus_docking_scores.csv"

    assert base_out.exists()
    assert prefixed_out.exists()
    assert dud_out.exists()

    base_rows = _read_csv_rows(base_out)
    prefixed_rows = _read_csv_rows(prefixed_out)
    dud_rows = _read_csv_rows(dud_out)

    assert len(base_rows) == 10
    assert len(prefixed_rows) == 10
    assert len(dud_rows) == 10
    _assert_consensus_desc(base_rows)
    _assert_consensus_desc(prefixed_rows)
    _assert_consensus_desc(dud_rows)

    base_set = {row.get("ligand", "") for row in base_rows}
    prefixed_set = {row.get("ligand", "") for row in prefixed_rows}
    assert base_set == set(base_ligands)
    assert prefixed_set == set(prefixed_ligands)
    assert base_set.isdisjoint(prefixed_set)

    assert all(row.get("run_mode") == "fda" for row in base_rows)
    assert all(row.get("run_mode") == "fda_dud" for row in prefixed_rows)
    assert all(row.get("run_mode") == "dud" for row in dud_rows)
    assert all(row.get("library") == "fda" for row in base_rows)
    assert all(row.get("library") == "fda_dud" for row in prefixed_rows)
    assert all(row.get("library") == "dud" for row in dud_rows)
