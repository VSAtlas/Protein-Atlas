from __future__ import annotations

import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import pathway_resolver as pr


class FakeHttpClient:
    def __init__(self, raise_on_call: bool = False) -> None:
        self._responses: dict[tuple[str, str, tuple[tuple[str, str], ...]], object] = {}
        self.raise_on_call = raise_on_call

    def _params_key(self, params: dict[str, str] | None) -> tuple[tuple[str, str], ...]:
        if not params:
            return ()
        return tuple(sorted((str(k), str(v)) for k, v in params.items()))

    def add_json(
        self, url: str, response: object, params: dict[str, str] | None = None
    ) -> None:
        key = ("json", url, self._params_key(params))
        self._responses[key] = response

    def add_text(
        self, url: str, response: object, params: dict[str, str] | None = None
    ) -> None:
        key = ("text", url, self._params_key(params))
        self._responses[key] = response

    def get_json(self, url: str, params: dict[str, str] | None = None) -> object:
        if self.raise_on_call:
            raise RuntimeError("Unexpected HTTP call")
        key = ("json", url, self._params_key(params))
        if key not in self._responses:
            raise KeyError(f"Missing JSON response for {key}")
        resp = self._responses[key]
        if isinstance(resp, Exception):
            raise resp
        return resp

    def get_text(self, url: str, params: dict[str, str] | None = None) -> str:
        if self.raise_on_call:
            raise RuntimeError("Unexpected HTTP call")
        key = ("text", url, self._params_key(params))
        if key not in self._responses:
            raise KeyError(f"Missing text response for {key}")
        resp = self._responses[key]
        if isinstance(resp, Exception):
            raise resp
        return str(resp)


def _logger() -> logging.Logger:
    logger = logging.getLogger("pathway_resolver_test")
    logger.handlers.clear()
    logger.addHandler(logging.NullHandler())
    return logger


def test_pathway_cache_hit(tmp_path: Path) -> None:
    cache_dir = tmp_path / "cache"
    cache = pr.Cache(cache_dir=cache_dir, refresh=False, logger=_logger())
    http = FakeHttpClient()

    search_params = {"query": "glycolysis", "species": "Homo sapiens", "types": "Pathway"}
    http.add_json(
        pr.REACTOME_SEARCH_URL,
        {
            "results": [
                {
                    "stId": "R-HSA-1",
                    "name": "Glycolysis",
                    "score": 5.0,
                    "species": [{"name": "Homo sapiens"}],
                }
            ]
        },
        params=search_params,
    )
    http.add_json(
        pr.REACTOME_PARTICIPANTS_URL.format(pathway_id="R-HSA-1"),
        [
            {"databaseName": "UniProt", "identifier": "Q22222"},
            {"databaseName": "UniProt", "identifier": "P11111"},
        ],
    )

    uniprots = pr.resolve_uniprots(
        "glycolysis", "reactome", "Homo sapiens", cache, http
    )
    assert uniprots == ["P11111", "Q22222"]
    cache_path = cache.pathway_cache_path("reactome", "Homo sapiens", "glycolysis")
    assert cache_path.exists()

    cache_reuse = pr.Cache(cache_dir=cache_dir, refresh=False, logger=_logger())
    http_fail = FakeHttpClient(raise_on_call=True)
    uniprots_cached = pr.resolve_uniprots(
        "glycolysis", "reactome", "Homo sapiens", cache_reuse, http_fail
    )
    assert uniprots_cached == ["P11111", "Q22222"]


def test_uniprot_cache_hit(tmp_path: Path) -> None:
    cache_dir = tmp_path / "cache"
    cache = pr.Cache(cache_dir=cache_dir, refresh=False, logger=_logger())
    http = FakeHttpClient()

    uniprot = "P12345"
    pdbe_url = f"{pr.PDBE_MAPPING_BASES[0]}/{uniprot}"
    http.add_json(
        pdbe_url,
        {
            "P12345": [
                {
                    "pdb_id": "1abc",
                    "resolution": 2.0,
                    "coverage": 0.7,
                    "method": "X-ray",
                    "mutation_count": 0,
                    "has_ligand": False,
                }
            ]
        },
    )

    candidates = pr.map_uniprot_to_pdb(uniprot, cache, http)
    assert candidates and candidates[0]["pdb_id"] == "1ABC"
    cache_path = cache.uniprot_cache_path(uniprot)
    assert cache_path.exists()

    cache_reuse = pr.Cache(cache_dir=cache_dir, refresh=False, logger=_logger())
    http_fail = FakeHttpClient(raise_on_call=True)
    candidates_cached = pr.map_uniprot_to_pdb(uniprot, cache_reuse, http_fail)
    assert candidates_cached == candidates


def test_selection_heuristics() -> None:
    candidates = [
        {"pdb_id": "1aaa", "has_ligand": False, "coverage": 0.9, "resolution": 1.5},
        {"pdb_id": "2bbb", "has_ligand": True, "coverage": 0.9, "resolution": 2.0},
    ]
    selected = pr.select_representatives("P11111", candidates, 1)
    assert selected == ["2BBB"]

    candidates = [
        {"pdb_id": "3ccc", "has_ligand": False, "coverage": 0.4, "resolution": 1.0},
        {"pdb_id": "4ddd", "has_ligand": False, "coverage": 0.8, "resolution": 2.5},
    ]
    selected = pr.select_representatives("P11111", candidates, 1)
    assert selected == ["4DDD"]

    candidates = [
        {"pdb_id": "2bbb", "has_ligand": False, "coverage": 0.5, "resolution": 2.0},
        {"pdb_id": "1aaa", "has_ligand": False, "coverage": 0.5, "resolution": 2.0},
    ]
    selected = pr.select_representatives("P11111", candidates, 2)
    assert selected == ["1AAA", "2BBB"]


def test_end_to_end_output(tmp_path: Path) -> None:
    cache_dir = tmp_path / "cache"
    output_path = tmp_path / "resolved_pdbs.txt"
    config_path = tmp_path / "config.txt"
    config_path.write_text(f"OVERALL_DIR={tmp_path}\n", encoding="utf-8")

    http = FakeHttpClient()
    search_params = {"query": "glycolysis", "species": "Homo sapiens", "types": "Pathway"}
    http.add_json(
        pr.REACTOME_SEARCH_URL,
        {
            "results": [
                {
                    "stId": "R-HSA-1",
                    "name": "Glycolysis",
                    "score": 5.0,
                    "species": [{"name": "Homo sapiens"}],
                }
            ]
        },
        params=search_params,
    )
    http.add_json(
        pr.REACTOME_PARTICIPANTS_URL.format(pathway_id="R-HSA-1"),
        [
            {"databaseName": "UniProt", "identifier": "P11111"},
            {"databaseName": "UniProt", "identifier": "Q22222"},
        ],
    )
    http.add_json(
        f"{pr.PDBE_MAPPING_BASES[0]}/P11111",
        {
            "P11111": [
                {
                    "pdb_id": "1abc",
                    "coverage": 0.5,
                    "resolution": 1.2,
                    "has_ligand": False,
                },
                {
                    "pdb_id": "2bcd",
                    "coverage": 0.7,
                    "resolution": 2.0,
                    "has_ligand": True,
                },
            ]
        },
    )
    http.add_json(
        f"{pr.PDBE_MAPPING_BASES[0]}/Q22222",
        {
            "Q22222": [
                {
                    "pdb_id": "3cde",
                    "coverage": 0.9,
                    "resolution": 2.5,
                    "has_ligand": False,
                }
            ]
        },
    )

    exit_code = pr.run(
        [
            "glycolysis",
            "--source",
            "reactome",
            "--organism",
            "Homo sapiens",
            "--output",
            str(output_path),
            "--cache-dir",
            str(cache_dir),
            "--config",
            str(config_path),
        ],
        http_client=http,
    )
    assert exit_code == 0
    lines = output_path.read_text(encoding="utf-8").strip().splitlines()
    assert lines == ["2BCD", "3CDE"]
