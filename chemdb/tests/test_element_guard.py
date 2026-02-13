from pathlib import Path

import pytest

import protein_prep.element_guard as element_guard


def _atom_line_with_element(element: str) -> str:
    return (
        "ATOM      1  CA  GLY A   1      11.000  12.000  13.000"
        f"  1.00 20.00          {element:>2}\n"
    )


def test_meeko_preflight_raises_for_helium_and_writes_prefile(tmp_path: Path, monkeypatch):
    input_pdb = tmp_path / "input_he.pdb"
    work_dir = tmp_path / "work"
    input_pdb.write_text(_atom_line_with_element("HE"), encoding="utf-8")

    monkeypatch.setattr(
        element_guard, "fix_element_columns_in_file", lambda *_a, **_k: 0
    )
    monkeypatch.setattr(
        element_guard,
        "scan_helium_counts_with_hits",
        lambda _txt, max_hits=5: (1, ["ATOM:1:HE"]),
    )
    monkeypatch.setattr(element_guard, "rules_version", lambda: "unit-test-rules")

    with pytest.raises(RuntimeError, match="helium_preflight_failed"):
        element_guard._meeko_preflight_or_fail(input_pdb, work_dir)

    prefile = work_dir / "meeko_input_pre_sanitize.pdb"
    assert prefile.exists()
    assert "HE" in prefile.read_text(encoding="utf-8")


def test_post_write_element_guard_applies_fixer(tmp_path: Path, monkeypatch):
    pdb_path = tmp_path / "post_write.pdb"
    pdb_path.write_text(_atom_line_with_element("HE"), encoding="utf-8")

    def fake_fix(src, dst_path=None, rewrite_atoms=False):
        target = Path(dst_path or src)
        txt = Path(src).read_text(encoding="utf-8")
        target.write_text(txt.replace(" HE\n", " H \n"), encoding="utf-8")
        return 1

    monkeypatch.setattr(element_guard, "fix_element_columns_in_file", fake_fix)
    monkeypatch.setattr(
        element_guard, "scan_helium_counts", lambda txt: 1 if " HE\n" in txt else 0
    )

    element_guard._post_write_element_guard("unit_step", pdb_path)

    output = pdb_path.read_text(encoding="utf-8")
    assert " HE\n" not in output
    assert " H \n" in output


def test_meeko_preflight_returns_resolved_path_when_no_helium(
    tmp_path: Path, monkeypatch
):
    input_pdb = tmp_path / "input_ok.pdb"
    work_dir = tmp_path / "work_ok"
    input_pdb.write_text(_atom_line_with_element("C"), encoding="utf-8")

    monkeypatch.setattr(
        element_guard, "fix_element_columns_in_file", lambda *_a, **_k: 0
    )
    monkeypatch.setattr(
        element_guard, "scan_helium_counts_with_hits", lambda _txt, max_hits=5: (0, [])
    )

    resolved = element_guard._meeko_preflight_or_fail(input_pdb, work_dir)

    assert resolved == input_pdb.resolve()
    assert (work_dir / "meeko_input_pre_sanitize.pdb").exists()
