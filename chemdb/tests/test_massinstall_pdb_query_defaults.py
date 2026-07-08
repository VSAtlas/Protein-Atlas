from __future__ import annotations

# ruff: noqa: E402

import argparse
import csv
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools import massinstall
from tools import pdb_query


def _write_genes_csv(path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "category"])
        writer.writeheader()
        writer.writerow({"gene": "HTR2B", "category": "GPCR"})


def _build_args(genes_file: Path) -> argparse.Namespace:
    return argparse.Namespace(
        pdb_query=True,
        genes_file=str(genes_file),
        gene_uniprot_file=None,
        forced_csv=None,
        trivial_ligands_file=None,
        keep_ligands_file=None,
        max_candidates=None,
        max_return=None,
        resolution_max=None,
        include_forced=None,
        dry_run=True,
        install=False,
        max_install=None,
        out_dir=None,
        selected_out=None,
        missing_out=None,
        category_limits_json=None,
        cache_dir=None,
        cache_ttl_days=None,
        no_cache=False,
        rate_limit=None,
        timeout=None,
        max_retries=None,
        ligand_filter_mode="strict",
    )


def test_select_from_pdb_query_uses_pdb_query_defaults_not_config(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    genes_file = tmp_path / "genes.csv"
    _write_genes_csv(genes_file)

    client_kwargs: dict[str, object] = {}
    call_kwargs: dict[str, object] = {}

    class FakeClient:
        def __init__(self, **kwargs: object) -> None:
            client_kwargs.update(kwargs)

    def fake_get_ranked_entries_for_gene(gene: str, **kwargs: object) -> pdb_query.SelectionResult:
        call_kwargs["gene"] = gene
        call_kwargs.update(kwargs)
        return pdb_query.SelectionResult(
            entries=[
                pdb_query.RankedEntry(
                    gene=gene,
                    pdb_id="1ABC",
                    resolution=2.0,
                    method="X-RAY",
                    nontrivial_comp_ids=["LIG"],
                    rank_score=7.0,
                    forced_included=False,
                    valid_entry=True,
                )
            ],
            stats={},
        )

    monkeypatch.setattr(massinstall, "PDBQueryClient", FakeClient)
    monkeypatch.setattr(
        massinstall,
        "get_ranked_entries_for_gene",
        fake_get_ranked_entries_for_gene,
    )

    cfg = {
        "PDB_QUERY_GENES_FILE": "ignored.csv",
        "PDB_QUERY_CATEGORY_LIMITS": '{"GPCR":0}',
        "PDB_QUERY_GENE_UNIPROT_FILE": "ignored_uniprots.csv",
        "PDB_QUERY_FORCED_CSV": "ignored_forced.csv",
        "PDB_QUERY_TRIVIAL_LIGANDS_FILE": "ignored_trivial.txt",
        "PDB_QUERY_KEEP_LIGANDS_FILE": "ignored_keep.txt",
        "PDB_QUERY_MAX_CANDIDATES": "9",
        "PDB_QUERY_MAX_RETURN": "1",
        "PDB_QUERY_RESOLUTION_MAX": "1.2",
        "PDB_QUERY_INCLUDE_FORCED": "false",
        "PDB_QUERY_RATE_LIMIT": "1.5",
        "PDB_QUERY_TIMEOUT": "5.0",
        "PDB_QUERY_MAX_RETRIES": "1",
        "PDB_QUERY_CACHE_DIR": "ignored-cache",
        "PDB_QUERY_CACHE_TTL_DAYS": "2",
        "PDB_QUERY_NO_CACHE": "true",
    }

    rows, pdb_ids = massinstall._select_from_pdb_query(cfg=cfg, args=_build_args(genes_file))

    assert pdb_ids == ["1ABC"]
    assert rows[0]["category"] == "GPCR"
    assert call_kwargs["gene"] == "HTR2B"
    assert call_kwargs["max_candidates"] == 200
    assert call_kwargs["max_return"] == 10
    assert call_kwargs["resolution_max"] == pdb_query.DEFAULT_RESOLUTION_MAX
    assert call_kwargs["include_forced"] is True
    assert call_kwargs["gene_to_uniprots"] == {}
    assert call_kwargs["trivial_ligands"] == set()
    assert call_kwargs["keep_ligands"] == set()
    assert call_kwargs["ligand_filter_mode"] == "strict"
    assert client_kwargs["cache_dir"] == pdb_query.DEFAULT_CACHE_DIR
    assert client_kwargs["use_cache"] is True
    assert client_kwargs["cache_ttl_seconds"] == pdb_query.DEFAULT_CACHE_TTL_SECONDS
    assert client_kwargs["rate_limit_rps"] == pdb_query.DEFAULT_RATE_LIMIT_RPS
    assert client_kwargs["timeout_seconds"] == pdb_query.DEFAULT_TIMEOUT_SECONDS
    assert client_kwargs["max_retries"] == pdb_query.DEFAULT_MAX_RETRIES


def test_select_from_pdb_query_uses_config_genes_file_when_cli_missing(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    genes_file = tmp_path / "genes.csv"
    _write_genes_csv(genes_file)
    args = _build_args(tmp_path / "unused.csv")
    args.genes_file = None

    class FakeClient:
        def __init__(self, **kwargs: object) -> None:
            del kwargs

    def fake_get_ranked_entries_for_gene(gene: str, **kwargs: object) -> pdb_query.SelectionResult:
        del kwargs
        return pdb_query.SelectionResult(
            entries=[
                pdb_query.RankedEntry(
                    gene=gene,
                    pdb_id="1ABC",
                    resolution=2.0,
                    method="X-RAY",
                    nontrivial_comp_ids=["LIG"],
                    rank_score=7.0,
                    forced_included=False,
                    valid_entry=True,
                )
            ],
            stats={},
        )

    monkeypatch.setattr(massinstall, "PDBQueryClient", FakeClient)
    monkeypatch.setattr(
        massinstall,
        "get_ranked_entries_for_gene",
        fake_get_ranked_entries_for_gene,
    )

    rows, pdb_ids = massinstall._select_from_pdb_query(
        cfg={"PDB_QUERY_GENES_FILE": str(genes_file)},
        args=args,
    )

    assert pdb_ids == ["1ABC"]
    assert rows[0]["gene"] == "HTR2B"
    assert rows[0]["category"] == "GPCR"


def test_select_from_pdb_query_errors_when_no_gene_source(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    args = _build_args(tmp_path / "unused.csv")
    args.genes_file = None
    monkeypatch.chdir(tmp_path)

    with pytest.raises(ValueError, match="needs a genes source"):
        massinstall._select_from_pdb_query(cfg={}, args=args)


def test_massinstall_help_runs_without_pythonpath() -> None:
    result = subprocess.run(
        [sys.executable, "tools/massinstall.py", "--help"],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    assert "Install PDB files from static list or PDB query" in result.stdout
