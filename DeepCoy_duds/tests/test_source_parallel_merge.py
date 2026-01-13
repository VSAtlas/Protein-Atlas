from generate_dud_library import query_external_sources


def make_opts():
    return {
        "timeout": 1,
        "retries": 0,
        "max_actives_per_source": 10,
        "affinity_cutoff_nm": 1000,
        "iuphar_approved_only": True,
        "iuphar_primary_only": True,
        "chembl_activity_types": ["Ki", "Kd"],
        "chembl_max_phase": None,
        "drugbank_data_dir": None,
    }


def test_deterministic_merge_parallel_vs_sequential(monkeypatch, tmp_path):
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()

    def stub_chembl(*args, **kwargs):
        return ["A", "B"], {
            "counts": {"targets": 1, "assays": 1, "activities": 2},
            "cached": False,
        }

    def stub_bindingdb(*args, **kwargs):
        return ["B", "C"], {
            "counts": {"pdb_hits": 1, "uniprot_hits": 1},
            "cached": False,
        }

    def stub_iuphar(*args, **kwargs):
        return ["C", "D"], {
            "counts": {"targets": 1, "interactions": 2},
            "cached": False,
        }

    monkeypatch.setattr("generate_dud_library.fetch_chembl_smiles", stub_chembl)
    monkeypatch.setattr("generate_dud_library.fetch_bindingdb_smiles", stub_bindingdb)
    monkeypatch.setattr("generate_dud_library.fetch_iuphar_smiles", stub_iuphar)

    sources = ["chembl", "bindingdb", "iuphar"]
    opts = make_opts()

    seq, seq_srcs, _ = query_external_sources(
        "U", [], "PDB", sources, cache_dir, opts, max_workers=1
    )
    par, par_srcs, _ = query_external_sources(
        "U", [], "PDB", sources, cache_dir, opts, max_workers=4
    )
    assert seq == par
    assert seq_srcs == par_srcs
    assert [s for s, _ in seq] == ["A", "B", "C", "D"]


def test_cached_source_merge_order(monkeypatch, tmp_path):
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()

    def stub_chembl(*args, **kwargs):
        return ["X"], {"counts": {"targets": 1}, "cached": True}

    def stub_bindingdb(*args, **kwargs):
        return ["Y"], {"counts": {"pdb_hits": 1}, "cached": False}

    monkeypatch.setattr("generate_dud_library.fetch_chembl_smiles", stub_chembl)
    monkeypatch.setattr("generate_dud_library.fetch_bindingdb_smiles", stub_bindingdb)
    monkeypatch.setattr(
        "generate_dud_library.fetch_iuphar_smiles",
        lambda *a, **k: ([], {"counts": {}, "cached": False}),
    )

    sources = ["chembl", "bindingdb", "iuphar"]
    opts = make_opts()

    seq, seq_srcs, _ = query_external_sources(
        "U", [], "PDB", sources, cache_dir, opts, max_workers=1
    )
    par, par_srcs, _ = query_external_sources(
        "U", [], "PDB", sources, cache_dir, opts, max_workers=2
    )
    assert seq == par
    assert seq_srcs == par_srcs
    assert [s for s, _ in seq] == ["X", "Y"]
