from pathlib import Path

import protein_prep.prep_utils as prep_utils
from protein_prep import layout


def test_fold_legacy_layout_moves_dirs_and_nolig_file(tmp_path: Path) -> None:
    root = tmp_path
    legacy_nolig = root / "TEST_nolig"
    legacy_ligs = root / "TEST_CLEANED_LIGANDS"
    legacy_nolig.mkdir(parents=True, exist_ok=True)
    legacy_ligs.mkdir(parents=True, exist_ok=True)
    (legacy_nolig / "a.txt").write_text("nolig", encoding="utf-8")
    (legacy_ligs / "b.txt").write_text("lig", encoding="utf-8")
    (root / "TEST_nolig.pdb").write_text("MODEL\nEND\n", encoding="utf-8")

    layout.fold_legacy_layout("TEST", root)

    base = root / "TEST"
    assert (base / "nolig" / "a.txt").exists()
    assert (base / "ligands_raw" / "b.txt").exists()
    moved_nolig = base / "nolig" / "TEST_nolig_phenix_clean.pdb"
    assert moved_nolig.exists()
    assert "MODEL" in moved_nolig.read_text(encoding="utf-8")
    assert not legacy_nolig.exists()
    assert not legacy_ligs.exists()


def test_canon_paths_respects_variant_and_default(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(prep_utils, "config", {"_CURRENT_VARIANT": "APO"}, raising=False)

    with_variant = layout.canon_paths("1abc", tmp_path, variant="HOLO")
    assert with_variant["protein_root"] == (tmp_path / "1ABC" / "HOLO").resolve()
    assert with_variant["variant"] == "HOLO"
    assert with_variant["receptor"] == (tmp_path / "1ABC" / "HOLO" / "receptor").resolve()

    no_variant = layout.canon_paths("1abc", tmp_path, variant=None)
    assert no_variant["protein_root"] == (tmp_path / "1ABC" / "APO").resolve()
