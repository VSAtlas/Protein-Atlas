import json
from pathlib import Path

from protein_prep.metal_site_audit import run_metal_site_audit


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


class _RouterPaths:
    def __init__(self, root: Path):
        self.root = root

    def docked_variant_root(self, variant: str, ph_label: str | None = None) -> Path:
        suffix = ph_label or "none"
        return self.root / variant / suffix


def test_run_metal_site_audit_writes_expected_donor_payload(tmp_path: Path) -> None:
    input_pdb = tmp_path / "input.pdb"
    receptor_pdbqt = tmp_path / "receptor.pdbqt"

    # ZN with two nearby donors in input; receptor keeps one donor to force diff.
    input_pdb.write_text(
        "".join(
            [
                _pdb_line("HETATM", 1, "ZN", "ZN", "A", 100, 0.0, 0.0, 0.0, "ZN"),
                _pdb_line("ATOM", 2, "ND1", "HIS", "A", 10, 1.8, 0.0, 0.0, "N"),
                _pdb_line("ATOM", 3, "SG", "CYS", "A", 11, 0.0, 2.2, 0.0, "S"),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )
    receptor_pdbqt.write_text(
        "".join(
            [
                _pdb_line("HETATM", 1, "ZN", "ZN", "A", 100, 0.0, 0.0, 0.0, "ZN"),
                _pdb_line("ATOM", 2, "ND1", "HIS", "A", 10, 1.8, 0.0, 0.0, "N"),
                "END\n",
            ]
        ),
        encoding="utf-8",
    )

    router = _RouterPaths(tmp_path / "docked")
    run_metal_site_audit(
        pdb_id="TEST",
        router_paths=router,
        input_pdb_path=str(input_pdb),
        receptor_pdb_path=None,
        receptor_pdbqt_path=str(receptor_pdbqt),
        variant_label="HOLO",
        ph_label="pH7_0",
    )

    out_json = router.docked_variant_root("HOLO", "pH7_0") / "metal_site_audit.json"
    assert out_json.exists()

    payload = json.loads(out_json.read_text(encoding="utf-8"))
    assert payload["pdb_id"] == "TEST"
    assert payload["metals"]

    metal_entry = payload["metals"][0]
    donor_atoms = {d["atom_name"] for d in metal_entry["donors_input"]}
    assert {"ND1", "SG"}.issubset(donor_atoms)

    lost_atoms = {d["atom_name"] for d in metal_entry["lost_donors"]}
    assert "SG" in lost_atoms
