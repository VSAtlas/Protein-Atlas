import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEEPCOY_DIR = REPO_ROOT / "DeepCoy_duds"
if str(DEEPCOY_DIR) not in sys.path:
    sys.path.insert(0, str(DEEPCOY_DIR))

from DeepCoy_duds.generate_dud_library import query_external_sources  # noqa: E402


def make_opts():
    return {
        "timeout": 1,
        "retries": 0,
        "max_actives_per_source": 5,
        "affinity_cutoff_nm": 1000,
        "iuphar_approved_only": True,
        "iuphar_primary_only": True,
        "chembl_activity_types": ["Ki", "Kd"],
        "chembl_max_phase": None,
        "chembl_max_pages": 2,
        "drugbank_data_dir": None,
    }


def test_query_external_sources_skips_chembl(monkeypatch, tmp_path):
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()

    calls = {"chembl": 0}

    def stub_chembl(*args, **kwargs):
        calls["chembl"] += 1
        return [], {"counts": {}, "cached": False}

    monkeypatch.setattr("generate_dud_library.fetch_chembl_smiles", stub_chembl)
    monkeypatch.setattr(
        "generate_dud_library.fetch_bindingdb_smiles",
        lambda *a, **k: ([], {"counts": {}, "cached": False}),
    )
    monkeypatch.setattr(
        "generate_dud_library.fetch_iuphar_smiles",
        lambda *a, **k: ([], {"counts": {}, "cached": False}),
    )

    actives, provenance, details = query_external_sources(
        "U", [], "PDB", ["bindingdb"], cache_dir, make_opts(), max_workers=1
    )
    assert actives == []
    assert provenance == {}
    assert details == {}
    assert calls["chembl"] == 0
