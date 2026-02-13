from pathlib import Path

from protein_prep import strip_nsr


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
    bfactor: float = 20.0,
) -> str:
    return (
        f"{record:<6}{serial:>5} {atom:<4} {resname:>3} {chain}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{bfactor:>6.2f}          {element:>2}\n"
    )


def _patch_strip_globals(monkeypatch, *, cfg: dict, retain: set[str]) -> None:
    monkeypatch.setattr(strip_nsr, "config", cfg, raising=False)
    monkeypatch.setattr(
        strip_nsr,
        "_resolve_variant_token",
        lambda _cfg, variant=None: variant or "HOLO",
        raising=False,
    )
    monkeypatch.setattr(
        strip_nsr,
        "_normalize_resname",
        lambda token: (token or "").strip().upper(),
        raising=False,
    )
    monkeypatch.setattr(strip_nsr, "_WATER_NAMES", {"HOH"}, raising=False)
    monkeypatch.setattr(strip_nsr, "_ELEM_CANON", {"NA", "ZN", "MG"}, raising=False)
    monkeypatch.setattr(strip_nsr, "_COFACTOR_NAMES", set(), raising=False)
    monkeypatch.setattr(strip_nsr, "_COFACTOR_CANONICAL", set(), raising=False)
    monkeypatch.setattr(strip_nsr, "_COFACTOR_RAW_ALL", set(), raising=False)
    monkeypatch.setattr(strip_nsr, "_RETAIN_VARIANT", set(retain), raising=False)
    monkeypatch.setattr(
        strip_nsr,
        "_RETAIN_VARIANT_CANONICAL",
        {r.upper() for r in retain},
        raising=False,
    )
    monkeypatch.setattr(strip_nsr, "_POLICY_MODE", "LEGACY", raising=False)
    monkeypatch.setattr(strip_nsr, "_ALIASES_BIND_LOGGED", False, raising=False)
    monkeypatch.setattr(strip_nsr, "_ION_KEEP_LOGGED", set(), raising=False)
    monkeypatch.setattr(strip_nsr, "_COFACTOR_DROP_LOGGED", set(), raising=False)
    monkeypatch.setattr(strip_nsr, "_ION_AUDIT_METALS", {"ZN", "MG"}, raising=False)
    monkeypatch.setattr(strip_nsr, "_ION_AUDIT_SIMPLE_IONS", {"NA", "CL"}, raising=False)
    monkeypatch.setattr(
        strip_nsr,
        "_parse_xyz",
        lambda line: (
            float(line[30:38]),
            float(line[38:46]),
            float(line[46:54]),
        ),
        raising=False,
    )
    monkeypatch.setattr(strip_nsr, "_post_write_element_guard", lambda *a, **k: None, raising=False)
    monkeypatch.setattr(
        strip_nsr,
        "_scan_metal_map",
        lambda _p, _rules=None: {},
        raising=False,
    )
    monkeypatch.setattr(
        strip_nsr,
        "_summarize_ions_file",
        lambda _p: {"res_hist": "none", "elem_hist": "none"},
        raising=False,
    )
    monkeypatch.setattr(
        strip_nsr,
        "ALIASES",
        type("Aliases", (), {"retain_resnames": []})(),
        raising=False,
    )


def test_strip_nonstandard_residues_remove_all_waters(monkeypatch, tmp_path: Path) -> None:
    _patch_strip_globals(
        monkeypatch,
        cfg={"WATER_POLICY": "remove_all", "WATER_SITE_RADIUS_ANG": 6.0, "WATER_MAX_BFACTOR": 60.0},
        retain=set(),
    )

    in_pdb = tmp_path / "input.pdb"
    out_pdb = tmp_path / "out.pdb"
    in_pdb.write_text(
        "".join(
            [
                _pdb_line("ATOM", 1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"),
                _pdb_line("HETATM", 2, "O", "HOH", "A", 10, 1.0, 0.0, 0.0, "O"),
                _pdb_line("HETATM", 3, "NA", "NA", "A", 20, 2.0, 0.0, 0.0, "NA"),
                _pdb_line("HETATM", 4, "C1", "LIG", "A", 30, 3.0, 0.0, 0.0, "C"),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    removed, out_path = strip_nsr.strip_nonstandard_residues(in_pdb, out_pdb)
    text = Path(out_path).read_text(encoding="utf-8")

    assert removed == 2
    assert " ALA " in text
    assert " NA " in text
    assert " HOH " not in text
    assert " LIG " not in text


def test_strip_nonstandard_residues_site_only_water_radius(monkeypatch, tmp_path: Path) -> None:
    _patch_strip_globals(
        monkeypatch,
        cfg={"WATER_POLICY": "site_only", "WATER_SITE_RADIUS_ANG": 2.0, "WATER_MAX_BFACTOR": 60.0},
        retain={"COF"},
    )

    in_pdb = tmp_path / "input_site.pdb"
    out_pdb = tmp_path / "out_site.pdb"
    in_pdb.write_text(
        "".join(
            [
                _pdb_line("ATOM", 1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"),
                _pdb_line("HETATM", 2, "ZN", "COF", "A", 5, 0.0, 0.0, 0.0, "ZN"),
                _pdb_line("HETATM", 3, "O", "HOH", "A", 11, 1.0, 0.0, 0.0, "O", bfactor=10.0),
                _pdb_line("HETATM", 4, "O", "HOH", "A", 12, 10.0, 0.0, 0.0, "O", bfactor=10.0),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    strip_nsr.strip_nonstandard_residues(in_pdb, out_pdb)
    text = out_pdb.read_text(encoding="utf-8")

    assert " A  11" in text
    assert " A  12" not in text
    assert " COF " in text
