from __future__ import annotations

import json
import shutil
from pathlib import Path

from config.runtime_config import load_config
from docking.pockets import active_site_runtime as activesite
from path_router import make_paths


def _hetatm_counts(pdb_path: Path) -> dict[tuple[str, str, str], int]:
    counts: dict[tuple[str, str, str], int] = {}
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip().upper()
            chain = (line[21:22] or "-").strip() or "-"
            resnum = (line[22:26] or "").strip() or "?"
            key = (resname, chain, resnum)
            counts[key] = counts.get(key, 0) + 1
    return counts


def _write_stub_6lu7(path: Path) -> None:
    lines = [
        "ATOM      1  N   ALA A   1       1.000   1.000   1.000  1.00 20.00           N",
        "HETATM    2  C1  PJE A 401       2.000   2.000   2.000  1.00 20.00           C",
        "HETATM    3  C2  PJE A 401       2.500   2.500   2.500  1.00 20.00           C",
        "HETATM    4  C1  02J B 501       3.000   3.000   3.000  1.00 20.00           C",
        "HETATM    5  N1  02J B 501       3.500   3.500   3.500  1.00 20.00           N",
        "HETATM    6  O   010 C 601       4.000   4.000   4.000  1.00 20.00           O",
        "HETATM    7  O   HOH A 700       5.000   5.000   5.000  1.00 20.00           O",
        "END",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_pockets_json_6lu7(tmp_path: Path, monkeypatch) -> None:
    repo_root = Path(__file__).resolve().parents[2]
    source_pdb = repo_root / "input_pdbs" / "6LU7.pdb"
    using_real_fixture = source_pdb.exists()
    input_dir = tmp_path / "input_pdbs"
    input_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = input_dir / "6LU7.pdb"
    if using_real_fixture:
        shutil.copyfile(source_pdb, pdb_path)
    else:
        _write_stub_6lu7(pdb_path)

    cfg = dict(load_config())
    cfg.update(
        {
            "OVERALL_DIR": str(tmp_path),
            "INPUT_DIR": str(input_dir),
            "OUTPUT_DIR": str(tmp_path / "processed_pdbs"),
            "DOCKED_DIR": str(tmp_path / "docked"),
            "LIGANDS_MOL2_DIR": str(tmp_path / "ligands_mol2"),
            "OUTPUT_LIGANDS_DIR": str(tmp_path / "prepped_ligands"),
            "RUN_ID": "POCKETS_TEST",
        }
    )

    for key in ("OUTPUT_DIR", "DOCKED_DIR", "LIGANDS_MOL2_DIR", "OUTPUT_LIGANDS_DIR"):
        Path(cfg[key]).mkdir(parents=True, exist_ok=True)

    # Refresh router roots to honor tmp_path overrides.
    import path_router.path_router as pr

    pr._ROUTER_ROOTS = None  # type: ignore[attr-defined]
    monkeypatch.setenv("ATLAS_RUN_ID", cfg["RUN_ID"])

    paths = make_paths(cfg, base_id="6LU7", pdb_file="6LU7.pdb")

    for key, value in cfg.items():
        if isinstance(value, (str, int, float)):
            monkeypatch.setenv(str(key), str(value))
    activesite.main(str(pdb_path))

    pockets_json = Path(paths.work_dir).parent / "pockets" / "pockets.json"
    assert pockets_json.exists(), f"pockets.json not found at {pockets_json}"

    data = json.loads(pockets_json.read_text(encoding="utf-8"))
    assert data["pdb_id"] == "6LU7"
    pockets = data.get("pockets", [])
    excluded = data.get("excluded", [])
    assert pockets, "Expected at least one pocket entry"

    included_resnames = {p["ligand_resname"] for p in pockets}
    assert {"PJE", "02J"}.issubset(included_resnames)

    excluded_resnames = {e["ligand_resname"] for e in excluded}
    assert "HOH" not in excluded_resnames
    if using_real_fixture:
        assert any(
            e.get("ligand_resname") == "010"
            and "METHANOL" in str(e.get("reason", "")).upper()
            for e in excluded
        )

    het_counts = _hetatm_counts(pdb_path)
    blocked_resnames = {"HOH", "WAT", "010"} if using_real_fixture else {"HOH", "WAT"}
    for pocket in pockets:
        key = (
            pocket["ligand_resname"],
            pocket["ligand_chain"],
            pocket["ligand_resnum"],
        )
        assert key in het_counts, f"Pocket key {key} missing from HETATM map"
        assert pocket["atom_count"] == het_counts[key]
        assert pocket["ligand_resname"] not in blocked_resnames
