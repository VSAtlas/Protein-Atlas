
from protein_prep import automate_protein_prep
import protein_prep.clean_pipeline as clean_pipeline


def test_clean_pdb_delegates_to_clean_pipeline(monkeypatch, tmp_path):
    calls = []

    sentinel = str(tmp_path / "sentinel_cleaned.pdb")

    def fake_clean_pdb(pdb_file, output_root, logger=None):
        calls.append((str(pdb_file), str(output_root), logger))
        return sentinel

    monkeypatch.setattr(clean_pipeline, "clean_pdb", fake_clean_pdb)

    out = automate_protein_prep.clean_pdb("in.pdb", "out_root")

    assert out == sentinel
    assert calls == [("in.pdb", "out_root", None)]
