from types import SimpleNamespace

import protein_prep.holo_restore as holo_restore


def test_holo_restore_noop_when_mode_not_holo(tmp_path, monkeypatch):
    holo_restore._hydrate_legacy_globals()

    in_pdb = tmp_path / "TEST.pdb"
    in_pdb.write_text("HETATM    1 ZN   ZN A   1      0.000   0.000   0.000  1.00 20.00          ZN\nEND\n", encoding="utf-8")
    cleaned = tmp_path / "cleaned.pdb"
    cleaned.write_text("ATOM      1  CA  ALA A   1      0.000   0.000   0.000  1.00 20.00           C\nEND\n", encoding="utf-8")

    cfg = {"APO_HOLO_MODE": "apo", "INPUT_DIR": str(tmp_path)}

    out = holo_restore._holo_restore_from_input_if_needed(
        pdb_id="TEST",
        cleaned_pdb=str(cleaned),
        output_pdbqt=str(tmp_path / "out.pdbqt"),
        config=cfg,
    )

    assert out == (0, 0, False)
    text = cleaned.read_text(encoding="utf-8")
    assert "ZN" not in text


def test_holo_restore_appends_missing_holo_candidates(tmp_path, monkeypatch):
    holo_restore._hydrate_legacy_globals()

    in_pdb = tmp_path / "TEST.pdb"
    in_pdb.write_text(
        "HETATM    1 ZN   ZN A   1      1.000   2.000   3.000  1.00 20.00          ZN\n"
        "END\n",
        encoding="utf-8",
    )
    cleaned = tmp_path / "cleaned.pdb"
    cleaned.write_text(
        "ATOM      1  CA  ALA A   1      0.000   0.000   0.000  1.00 20.00           C\n"
        "END\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(
        holo_restore, "_resolve_variant_token", lambda _cfg: "HOLO", raising=False
    )
    monkeypatch.setattr(
        holo_restore,
        "make_paths",
        lambda config, base_id, pdb_file: SimpleNamespace(input_pdb_path=in_pdb),
    )
    monkeypatch.setattr(holo_restore, "load_canonical_metals", lambda _x: {"ZN"})
    monkeypatch.setattr(holo_restore, "load_canonical_cofactors", lambda _x: set())
    monkeypatch.setattr(holo_restore, "load_canonical_waters", lambda _x: {"HOH"})
    monkeypatch.setattr(
        holo_restore,
        "_parse_atoms_from_pdb_like_lines",
        lambda input_lines, canonical_metals: ([], []),
    )
    monkeypatch.setattr(holo_restore, "_find_metal_donors", lambda *a, **k: [])
    monkeypatch.setattr(holo_restore, "fix_element_columns_in_file", lambda *a, **k: None)

    cfg = {"APO_HOLO_MODE": "holo", "INPUT_DIR": str(tmp_path)}

    out = holo_restore._holo_restore_from_input_if_needed(
        pdb_id="TEST",
        cleaned_pdb=str(cleaned),
        output_pdbqt=str(tmp_path / "out.pdbqt"),
        config=cfg,
    )

    assert out == (1, 0, True)
    text = cleaned.read_text(encoding="utf-8")
    assert "ZN" in text
