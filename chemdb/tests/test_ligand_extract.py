from pathlib import Path

from protein_prep import ligand_extract


def _pdb_line(
    record: str,
    serial: int,
    atom: str,
    resname: str,
    chain: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    return (
        f"{record:<6}{serial:>5} {atom:<4} {resname:>3} {chain}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{20.00:>6.2f}          {element:>2}\n"
    )


def _patch_ligand_extract_globals(monkeypatch) -> None:
    monkeypatch.setattr(
        ligand_extract,
        "_hydrate_legacy_globals",
        lambda: None,
        raising=False,
    )
    monkeypatch.setattr(
        ligand_extract,
        "fix_element_columns_in_file",
        lambda *_args, **_kwargs: None,
        raising=False,
    )
    monkeypatch.setattr(
        ligand_extract,
        "ALIASES",
        type("Aliases", (), {"retain_resnames": []})(),
        raising=False,
    )
    monkeypatch.setattr(
        ligand_extract,
        "_normalize_resname",
        lambda token: (token or "").strip().upper(),
        raising=False,
    )
    monkeypatch.setattr(ligand_extract, "_WATER_NAMES", {"HOH"}, raising=False)
    monkeypatch.setattr(ligand_extract, "_RETAIN_VARIANT", set(), raising=False)
    monkeypatch.setattr(
        ligand_extract,
        "_RETAIN_VARIANT_CANONICAL",
        set(),
        raising=False,
    )


def test_extract_ligands_from_filtered_creates_ligand_files(monkeypatch, tmp_path: Path) -> None:
    _patch_ligand_extract_globals(monkeypatch)

    filtered = tmp_path / "filtered.pdb"
    out_dir = tmp_path / "ligands_raw"
    filtered.write_text(
        "".join(
            [
                _pdb_line("ATOM", 1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"),
                _pdb_line("HETATM", 2, "C1", "LIG", "A", 10, 1.0, 0.0, 0.0, "C"),
                _pdb_line("HETATM", 3, "N1", "LIG", "A", 10, 2.0, 0.0, 0.0, "N"),
                _pdb_line("HETATM", 4, "O", "HOH", "A", 11, 3.0, 0.0, 0.0, "O"),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    files = ligand_extract.extract_ligands_from_filtered(filtered, out_dir)

    assert len(files) == 1
    lig_file = files[0]
    assert lig_file.exists()
    text = lig_file.read_text(encoding="utf-8")
    assert "REMARK Extracted LIG_A10" in text
    assert " LIG " in text
    assert " HOH " not in text


def test_write_pristine_reference_creates_stub(monkeypatch, tmp_path: Path) -> None:
    _patch_ligand_extract_globals(monkeypatch)
    monkeypatch.setattr(ligand_extract, "_HAS_RDKIT", False, raising=False)

    lig_dir = tmp_path / "ligands_raw"
    lig_dir.mkdir(parents=True, exist_ok=True)
    lig_file = lig_dir / "LIG_A10.pdb"
    lig_file.write_text(
        _pdb_line("HETATM", 1, "C1", "LIG", "A", 10, 0.0, 0.0, 0.0, "C"),
        encoding="utf-8",
    )

    ligand_extract._write_pristine_reference(lig_file)

    ref_sdf = tmp_path / "reference" / "LIG_A10.sdf"
    assert ref_sdf.exists()
    assert "-Pristine-" in ref_sdf.read_text(encoding="utf-8")


def test_expose_ligand_intermediates_for_debug_noop_when_missing_source(
    monkeypatch,
    tmp_path: Path,
) -> None:
    _patch_ligand_extract_globals(monkeypatch)

    missing_src = tmp_path / "does_not_exist"
    link_dir = tmp_path / "debug_link"

    ligand_extract.expose_ligand_intermediates_for_debug(missing_src, link_dir)

    assert not link_dir.exists()
