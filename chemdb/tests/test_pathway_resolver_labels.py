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
    logger = logging.getLogger("pathway_resolver_label_test")
    logger.handlers.clear()
    logger.addHandler(logging.NullHandler())
    return logger


def test_label_from_title_variants() -> None:
    assert (
        pr.label_from_title(
            "Crystal structure of C- Kit tyrosine kinase bound to inhibitor"
        )
        == "C-KIT TYROSINE KINASE"
    )
    assert (
        pr.label_from_title(
            "Structural basis of activation of EGFR in complex with inhibitor"
        )
        == "EGFR"
    )
    assert pr.label_from_title("Solution structure of P53 with DNA") == "P53"
    assert (
        pr.label_from_title("Structure of the human enzyme ABC-1.")
        == "THE HUMAN ENZYME ABC-1"
    )


def test_run_writes_labels_and_cache(tmp_path: Path) -> None:
    config_path = tmp_path / "config.txt"
    config_path.write_text(
        f"OVERALL_DIR={tmp_path}\npathway_organism=Homo sapiens\n", encoding="utf-8"
    )

    http = FakeHttpClient()
    search_params = {
        "query": "glycolysis",
        "species": "Homo sapiens",
        "types": "Pathway",
    }
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
        ],
    )
    http.add_json(
        f"{pr.PDBE_MAPPING_BASES[0]}/P11111",
        {
            "P11111": [
                {
                    "pdb_id": "1abc",
                    "coverage": 0.7,
                    "resolution": 1.8,
                }
            ]
        },
    )
    http.add_json(
        f"{pr.PDBE_LIGAND_URL}/1abc",
        {"1abc": [{"chem_comp_id": "ATP"}]},
    )
    http.add_json(
        f"{pr.PDBE_SUMMARY_URL}/1abc",
        {
            "1abc": [
                {
                    "title": "Crystal structure of C- Kit tyrosine kinase bound to inhibitor"
                }
            ]
        },
    )

    exit_code = pr.run(
        [
            "glycolysis",
            "--source",
            "reactome",
            "--no-catalyst-only",
            "--config",
            str(config_path),
        ],
        http_client=http,
    )
    assert exit_code == 0

    output_path = tmp_path / "pathways" / "resolved_pdbs.txt"
    ids_path = tmp_path / "pathways" / "resolved_pdbs_ids.txt"
    label_lines = output_path.read_text(encoding="utf-8").strip().splitlines()
    assert label_lines == ["1ABC\tC-KIT TYROSINE KINASE"]
    id_lines = ids_path.read_text(encoding="utf-8").strip().splitlines()
    assert id_lines == ["1ABC"]

    cache_dir = tmp_path / "pathways" / "cache"
    cache = pr.Cache(cache_dir=cache_dir, refresh=False, logger=_logger())
    assert cache.pathway_cache_path(
        "reactome", "Homo sapiens", "glycolysis", False
    ).exists()
    assert cache.uniprot_cache_path("P11111").exists()
    assert cache.pdb_cache_path("1ABC").exists()
