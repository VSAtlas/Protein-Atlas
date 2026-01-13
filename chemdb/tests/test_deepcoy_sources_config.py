import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEEPCOY_DIR = REPO_ROOT / "DeepCoy_duds"
if str(DEEPCOY_DIR) not in sys.path:
    sys.path.insert(0, str(DEEPCOY_DIR))

from generate_dud_library import (  # noqa: E402
    DEFAULT_SOURCES,
    build_source_audit_lines,
    _filter_known_sources,
    resolve_sources,
)


def _make_source_opts():
    return {
        "affinity_cutoff_nm": 10000,
        "chembl_activity_types": ["Ki", "Kd", "IC50"],
        "iuphar_approved_only": True,
        "iuphar_primary_only": True,
        "drugbank_data_dir": None,
    }


def test_sources_resolution_precedence():
    cfg = {"DEEPCOY_SOURCES": "bindingdb,iuphar"}
    resolved, origin, raw_value = resolve_sources(None, cfg)
    assert resolved == ["bindingdb", "iuphar"]
    assert origin == "config"
    assert raw_value == "bindingdb,iuphar"

    resolved_cli, origin_cli, _ = resolve_sources("chembl,rcsb", cfg)
    assert resolved_cli == ["chembl", "rcsb"]
    assert origin_cli == "cli"

    resolved_default, origin_default, _ = resolve_sources(None, {})
    assert resolved_default == DEFAULT_SOURCES
    assert origin_default == "default"


def test_active_sources_config_parsing_alias_and_unknown():
    cfg = {"DEEPCOY_ACTIVE_SOURCES": " ChEMBL , pdb_ligands , mystery "}
    resolved, origin, raw_value = resolve_sources(None, cfg)
    assert raw_value.strip().lower().startswith("chembl")
    assert origin == "config_active"
    filtered, unknown = _filter_known_sources(resolved)
    assert filtered == ["chembl", "rcsb"]
    assert unknown == ["mystery"]


def test_audit_url_rendering():
    source_opts = _make_source_opts()
    uniprot = "P00520"
    pdb_id = "1IEP"
    ec_numbers = ["2.7.10.2"]

    lines = build_source_audit_lines(
        uniprot,
        ec_numbers,
        pdb_id,
        source_opts,
        ["chembl", "bindingdb", "iuphar", "rcsb", "rhea"],
    )

    joined = "\n".join(lines)
    assert "target_components__accession=P00520" in joined
    assert "pdb=1IEP" in joined and "cutoff=10000" in joined
    assert "accession=P00520" in joined
    assert "ecNumber=2.7.10.2" in joined
    assert "/core/entry/1IEP" in joined

    lines_no_ec = build_source_audit_lines(
        uniprot,
        [],
        pdb_id,
        source_opts,
        ["iuphar"],
    )
    assert not any("ecNumber=" in line for line in lines_no_ec)
