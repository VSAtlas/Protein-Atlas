import csv
import logging
from pathlib import Path
from unittest import mock
from urllib.parse import parse_qs, urlparse

from chemdb.tools import fill_fda_mapping_names as tool


class FakeResponse:
    def __init__(self, payload, status_code=200):
        self._payload = payload
        self.status_code = status_code
        self.headers = {}

    def json(self):
        return self._payload


def _write_csv(path: Path, headers: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)


def _run_tool(
    in_csv: Path,
    out_csv: Path,
    cache_dir: Path,
    *,
    max_retries: int = 1,
    gsrs_base_url: str = "https://gsrs.example/api/v1",
) -> int:
    return tool.main(
        [
            "--in_csv",
            str(in_csv),
            "--out_csv",
            str(out_csv),
            "--enable_unichem",
            "false",
            "--enable_pubchem_rerank",
            "false",
            "--enable_gsrs",
            "true",
            "--gsrs_base_url",
            gsrs_base_url,
            "--cache_dir",
            str(cache_dir),
            "--sleep",
            "0",
            "--qps",
            "50",
            "--max_retries",
            str(max_retries),
            "--backoff_base",
            "0",
            "--backoff_max",
            "0",
            "--max_workers",
            "1",
        ]
    )


def _query_param(url: str, key: str) -> str:
    parsed = urlparse(url)
    values = parse_qs(parsed.query).get(key, [])
    return values[0] if values else ""


def test_gsrs_resolves_by_unii(tmp_path, capsys):
    in_csv = tmp_path / "in.csv"
    out_csv = tmp_path / "out.csv"
    cache_dir = tmp_path / "cache"
    headers = ["file_num", "display_name", "pubchem_unii_list", "pubchem_iupac_name"]
    _write_csv(
        in_csv,
        headers,
        [
            {
                "file_num": "1",
                "display_name": "UNK_1",
                "pubchem_unii_list": "ABCD123456",
                "pubchem_iupac_name": "",
            }
        ],
    )

    def fake_get(url: str, timeout: int = 10):
        query = _query_param(url, "q")
        if query == "approvalID:ABCD123456":
            return FakeResponse(
                {
                    "content": [
                        {"names": [{"name": "Metformin", "preferred": True}]}
                    ]
                }
            )
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        exit_code = _run_tool(in_csv, out_csv, cache_dir)
    assert exit_code == 0

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["display_name"] == "Metformin"
    output = capsys.readouterr().out
    assert "gsrs resolved total: 1" in output
    assert "gsrs resolved by unii: 1" in output


def test_gsrs_inchikey_2d_fallback(tmp_path, capsys):
    in_csv = tmp_path / "in.csv"
    out_csv = tmp_path / "out.csv"
    cache_dir = tmp_path / "cache"
    headers = ["file_num", "display_name", "inchikey", "pubchem_iupac_name"]
    _write_csv(
        in_csv,
        headers,
        [
            {
                "file_num": "2",
                "display_name": "UNK_2",
                "inchikey": "APRTUXPAFSNIOL-MSOLQXFVSA-O",
                "pubchem_iupac_name": "",
            }
        ],
    )

    def fake_get(url: str, timeout: int = 10):
        query = _query_param(url, "q")
        if "APRTUXPAFSNIOL-MSOLQXFVSA-O" in query:
            return FakeResponse({"content": []})
        if "APRTUXPAFSNIOL" in query:
            return FakeResponse({"content": [{"preferredName": "Amlodipine"}]})
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        exit_code = _run_tool(in_csv, out_csv, cache_dir)
    assert exit_code == 0

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["display_name"] == "Amlodipine"
    output = capsys.readouterr().out
    assert "gsrs resolved by inchikey 2d: 1" in output


def test_gsrs_iupac_only_name_is_rejected(tmp_path):
    in_csv = tmp_path / "in.csv"
    out_csv = tmp_path / "out.csv"
    cache_dir = tmp_path / "cache"
    headers = ["file_num", "display_name", "pubchem_unii_list", "pubchem_iupac_name"]
    _write_csv(
        in_csv,
        headers,
        [
            {
                "file_num": "3",
                "display_name": "UNK_3",
                "pubchem_unii_list": "WXYZ987654",
                "pubchem_iupac_name": "2,3-dimethyl-1H-indole",
            }
        ],
    )

    def fake_get(url: str, timeout: int = 10):
        query = _query_param(url, "q")
        if query == "approvalID:WXYZ987654":
            return FakeResponse(
                {
                    "content": [
                        {
                            "names": [
                                {"name": "2,3-dimethyl-1H-indole", "preferred": True}
                            ]
                        }
                    ]
                }
            )
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        exit_code = _run_tool(in_csv, out_csv, cache_dir)
    assert exit_code == 0

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["display_name"] == "UNK_3"


def test_gsrs_cache_only_success_payloads(tmp_path):
    cache_dir = tmp_path / "cache"
    client = tool.CachedJsonClient(
        cache_dir=cache_dir,
        sleep=0.0,
        cache_mode="use",
        rate_limiter=None,
        max_retries=0,
        backoff_base=0.0,
        backoff_max=0.0,
    )
    gsrs = tool.GSRSClient("https://gsrs.example/api/v1", client)

    def fake_get_fail(url: str, timeout: int = 10):
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get_fail):
        assert gsrs.search_by_unii("ABCD123456") == []
    assert list(cache_dir.rglob("*.json")) == []

    def fake_get_ok(url: str, timeout: int = 10):
        query = _query_param(url, "q")
        if query == "approvalID:ABCD123456":
            return FakeResponse({"content": [{"name": "Dapagliflozin"}]})
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get_ok):
        records = gsrs.search_by_unii("ABCD123456")
    assert records
    assert list(cache_dir.rglob("*.json"))


def test_gsrs_retry_429_then_200_resolves_and_logs(tmp_path, caplog):
    in_csv = tmp_path / "in.csv"
    out_csv = tmp_path / "out.csv"
    cache_dir = tmp_path / "cache"
    headers = ["file_num", "display_name", "pubchem_unii_list", "pubchem_iupac_name"]
    _write_csv(
        in_csv,
        headers,
        [
            {
                "file_num": "4",
                "display_name": "UNK_4",
                "pubchem_unii_list": "QWER123456",
                "pubchem_iupac_name": "",
            }
        ],
    )
    call_state = {"count": 0}

    def fake_get(url: str, timeout: int = 10):
        query = _query_param(url, "q")
        if query == "approvalID:QWER123456":
            call_state["count"] += 1
            if call_state["count"] == 1:
                return FakeResponse({}, status_code=429)
            return FakeResponse({"content": [{"name": "Losartan"}]})
        return FakeResponse({}, status_code=404)

    caplog.set_level(logging.WARNING)
    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        exit_code = _run_tool(in_csv, out_csv, cache_dir, max_retries=1)
    assert exit_code == 0

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0]["display_name"] == "Losartan"
    assert any("status=429" in record.getMessage() for record in caplog.records)


def test_gsrs_multirow_fixture_resolves_all(tmp_path):
    in_csv = tmp_path / "in.csv"
    out_csv = tmp_path / "out.csv"
    cache_dir = tmp_path / "cache"
    headers = [
        "file_num",
        "display_name",
        "pubchem_unii_list",
        "inchikey",
        "smiles",
        "pubchem_iupac_name",
    ]
    _write_csv(
        in_csv,
        headers,
        [
            {
                "file_num": "11",
                "display_name": "UNK_11",
                "pubchem_unii_list": "AAAA111111",
                "inchikey": "",
                "smiles": "",
                "pubchem_iupac_name": "",
            },
            {
                "file_num": "12",
                "display_name": "UNK_12",
                "pubchem_unii_list": "",
                "inchikey": "BBBBBBBBBBBBBB-CCCCCCCCCC-Z",
                "smiles": "",
                "pubchem_iupac_name": "",
            },
            {
                "file_num": "13",
                "display_name": "UNK_13",
                "pubchem_unii_list": "",
                "inchikey": "",
                "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
                "pubchem_iupac_name": "",
            },
        ],
    )

    def fake_get(url: str, timeout: int = 10):
        query = _query_param(url, "q")
        if query == "approvalID:AAAA111111":
            return FakeResponse({"content": [{"name": "Empagliflozin"}]})
        if "BBBBBBBBBBBBBB-CCCCCCCCCC-Z" in query:
            return FakeResponse({"content": []})
        if query == "inchikey:BBBBBBBBBBBBBB":
            return FakeResponse({"content": [{"preferredName": "Rivaroxaban"}]})
        if urlparse(url).path.endswith("/substances/structureSearch"):
            return FakeResponse({"content": [{"displayName": "Aspirin"}]})
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        exit_code = _run_tool(in_csv, out_csv, cache_dir)
    assert exit_code == 0

    with out_csv.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows and all(not row["display_name"].upper().startswith("UNK_") for row in rows)
