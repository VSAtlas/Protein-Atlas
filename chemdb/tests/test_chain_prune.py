from pathlib import Path

from protein_prep import chain_prune


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


def _patch_chain_cfg(monkeypatch) -> None:
    monkeypatch.setattr(
        chain_prune, "_hydrate_legacy_globals", lambda: None, raising=False
    )
    monkeypatch.setattr(
        chain_prune, "_cfg_float", lambda key, default=0.0: 6.0, raising=False
    )
    monkeypatch.setattr(
        chain_prune, "_cfg_int", lambda key, default=0: 200, raising=False
    )
    monkeypatch.setattr(chain_prune, "_cfg_chain_keep_list", lambda: [], raising=False)
    monkeypatch.setattr(
        chain_prune, "_cfg_bool", lambda key, default=False: False, raising=False
    )
    monkeypatch.setattr(
        chain_prune, "_post_write_element_guard", lambda *a, **k: None, raising=False
    )
    monkeypatch.setattr(
        chain_prune,
        "_normalize_resname",
        lambda token: (token or "").strip().upper(),
        raising=False,
    )
    monkeypatch.setattr(chain_prune, "_RETAIN_VARIANT", set(), raising=False)
    monkeypatch.setattr(chain_prune, "_RETAIN_VARIANT_CANONICAL", set(), raising=False)


def test_select_chains_to_keep_and_prune(monkeypatch, tmp_path: Path) -> None:
    _patch_chain_cfg(monkeypatch)

    pdb = tmp_path / "filtered.pdb"
    ligands_dir = tmp_path / "ligands_raw"
    ligands_dir.mkdir(parents=True, exist_ok=True)
    ligand = ligands_dir / "lig_a.pdb"

    pdb.write_text(
        "".join(
            [
                _pdb_line("ATOM", 1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"),
                _pdb_line("ATOM", 2, "N", "ALA", "A", 1, 1.0, 0.0, 0.0, "N"),
                _pdb_line("ATOM", 3, "CA", "GLY", "B", 1, 50.0, 50.0, 50.0, "C"),
                _pdb_line("ATOM", 4, "N", "GLY", "B", 1, 51.0, 50.0, 50.0, "N"),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    ligand.write_text(
        "".join(
            [
                _pdb_line("HETATM", 10, "C1", "LIG", "A", 10, 0.5, 0.0, 0.0, "C"),
                _pdb_line("HETATM", 11, "O1", "LIG", "A", 10, 0.0, 0.5, 0.0, "O"),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    kept = chain_prune.select_chains_to_keep(pdb, ligands_dir)
    assert kept == {"A"}

    pruned = tmp_path / "pruned.pdb"
    chain_prune.prune_to_chains(pdb, kept, pruned)

    text = pruned.read_text(encoding="utf-8")
    assert " A   1" in text
    assert " B   1" not in text
