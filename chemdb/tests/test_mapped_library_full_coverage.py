import json
import logging
from pathlib import Path

from docking.docking_ligands import prepare_and_filter_ligands
from path_router import make_paths


def _write_fake_lig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["REMARK fake ligand for mapped-library coverage test"]
    for i in range(1, 9):
        lines.append(
            f"ATOM  {i:5d}  C   LIG A   1       0.000   0.000   0.000  1.00  0.00           C"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_prepare_and_filter_ligands_covers_full_mapped_library_with_partial_entries_manifest(
    tmp_path,
):
    prepped_root = tmp_path / "prepped_ligands"
    lib_root = prepped_root / "bench_fabp4"
    (prepped_root / "BNNQ").mkdir(parents=True, exist_ok=True)

    lig_paths = [
        lib_root / "actives" / "actives_final_00001.pdbqt",
        lib_root / "actives" / "actives_final_00002.pdbqt",
        lib_root / "decoys" / "decoys_final_00001.pdbqt",
        lib_root / "decoys" / "decoys_final_00002.pdbqt",
        lib_root / "decoys" / "decoys_final_00003.pdbqt",
        lib_root / "decoys" / "decoys_final_00004.pdbqt",
        lib_root / "decoys" / "decoys_final_00005.pdbqt",
        lib_root / "decoys" / "decoys_final_00006.pdbqt",
    ]
    for lig in lig_paths:
        _write_fake_lig(lig)

    all_rel = [str(p.relative_to(lib_root).as_posix()) for p in lig_paths]
    all_rel_sorted = sorted(all_rel)
    manifest_entries = {
        Path(rel).stem.lower(): rel for rel in all_rel_sorted[:3]
    }  # intentionally undercounted
    manifest_filenames = {
        Path(rel).name.lower(): rel for rel in all_rel_sorted
    }  # complete inventory
    manifest_payload = {
        "mtime": lib_root.stat().st_mtime,
        "entries": manifest_entries,
        "filenames": manifest_filenames,
    }
    (lib_root / "_manifest.json").write_text(
        json.dumps(manifest_payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "INPUT_DIR": str(tmp_path / "input_pdbs"),
        "OUTPUT_DIR": str(tmp_path / "processed_pdbs"),
        "DOCKED_DIR": str(tmp_path / "docked"),
        "PREPPED_LIGANDS_DIR": str(prepped_root),
        "OUTPUT_LIGANDS_DIR": str(prepped_root),
        "LIBRARY_SUBDIR_DEFAULT": "fda_library",
        "TEST_MODE_ENABLE": "dud",
        "TEST_LIBRARY_MAP": {"BNNQ": "bench_fabp4"},
        "LIBRARY_MANIFEST_BUILD_ON_SCAN": False,
        "PH_LIGAND_MODE": "off",
        "RUN_ID": "mapped_lib_cov",
    }
    for key in ("INPUT_DIR", "OUTPUT_DIR", "DOCKED_DIR"):
        Path(cfg[key]).mkdir(parents=True, exist_ok=True)

    paths = make_paths(cfg, base_id="BNNQ", pdb_file="BNNQ.pdb")
    logger = logging.getLogger("test.mapped_library_full_coverage")

    ligands_1, _, _ = prepare_and_filter_ligands(cfg, paths, logger)
    ligands_2, _, _ = prepare_and_filter_ligands(cfg, paths, logger)

    expected = {str(p.resolve()) for p in lig_paths}
    selected_1 = sorted(
        str(Path(p).resolve())
        for p in ligands_1
        if str(Path(p).resolve()).startswith(str(lib_root.resolve()))
    )
    selected_2 = sorted(
        str(Path(p).resolve())
        for p in ligands_2
        if str(Path(p).resolve()).startswith(str(lib_root.resolve()))
    )

    assert set(selected_1) == expected
    assert len(selected_1) == len(expected)
    assert selected_1 == selected_2
