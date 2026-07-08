from __future__ import annotations

from pathlib import Path

import main as main_mod


def test_build_combo_work_items_uses_manifest_tags(tmp_path: Path, monkeypatch) -> None:
    cfg = {"PH_ENSEMBLE": True, "INPUT_DIR": str(tmp_path)}
    (tmp_path / "ABCD.pdb").write_text("PDB", encoding="utf-8")

    monkeypatch.setattr(
        main_mod,
        "load_ph_tags",
        lambda pdb_id, variant=None: ["pH7_0", "pH8_0"],
    )
    monkeypatch.setattr(
        main_mod,
        "select_ph_values_for_protonation",
        lambda _path: [7.2],
    )

    items = main_mod._build_combo_work_items(
        cfg, ["ABCD.pdb"], variant_token="HOLO"
    )
    assert items == [("ABCD.pdb", "pH7_0"), ("ABCD.pdb", "pH8_0")]


def test_build_combo_work_items_falls_back_to_context_ph(
    tmp_path: Path, monkeypatch
) -> None:
    cfg = {"PH_ENSEMBLE": True, "INPUT_DIR": str(tmp_path)}
    (tmp_path / "EFGH.pdb").write_text("PDB", encoding="utf-8")

    monkeypatch.setattr(main_mod, "load_ph_tags", lambda pdb_id, variant=None: [])
    monkeypatch.setattr(
        main_mod,
        "select_ph_values_for_protonation",
        lambda _path: [7.2, 8.4],
    )

    items = main_mod._build_combo_work_items(
        cfg, ["EFGH.pdb"], variant_token="HOLO"
    )
    assert items == [("EFGH.pdb", "pH7_2"), ("EFGH.pdb", "pH8_4")]


def test_build_combo_work_items_non_ph_mode(tmp_path: Path) -> None:
    cfg = {"PH_ENSEMBLE": False, "INPUT_DIR": str(tmp_path)}
    items = main_mod._build_combo_work_items(
        cfg, ["IJKL.pdb"], variant_token="HOLO"
    )
    assert items == [("IJKL.pdb", None)]


def test_build_combo_work_items_string_false_is_non_ph_mode(tmp_path: Path) -> None:
    cfg = {"PH_ENSEMBLE": "false", "INPUT_DIR": str(tmp_path)}
    items = main_mod._build_combo_work_items(
        cfg, ["MNOP.pdb"], variant_token="HOLO"
    )
    assert items == [("MNOP.pdb", None)]
