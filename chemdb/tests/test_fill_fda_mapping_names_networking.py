import json
import logging
from pathlib import Path
from unittest import mock


from chemdb.tools import fill_fda_mapping_names as tool


class FakeResponse:
    def __init__(self, payload, status_code=200, headers=None):
        self._payload = payload
        self.status_code = status_code
        self.headers = headers or {}

    def json(self):
        return self._payload


def _write_cache(cache_path: Path, url: str, payload: dict) -> None:
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(
        json.dumps({"url": url, "payload": payload}, ensure_ascii=True),
        encoding="utf-8",
    )


def _make_client(
    cache_dir: Path,
    *,
    cache_mode: str,
    max_retries: int = 1,
    backoff_base: float = 0.0,
    backoff_max: float = 0.0,
) -> tool.CachedJsonClient:
    return tool.CachedJsonClient(
        cache_dir,
        sleep=0,
        cache_mode=cache_mode,
        rate_limiter=None,
        max_retries=max_retries,
        backoff_base=backoff_base,
        backoff_max=backoff_max,
    )


def test_do_not_cache_failures(tmp_path):
    cache_dir = tmp_path / "cache"
    client = _make_client(cache_dir, cache_mode="use")
    url = "https://example.test/resource"
    cache_path = client._cache_path(url)
    calls = {"count": 0}

    def fake_get(url_arg, timeout=10):
        calls["count"] += 1
        if calls["count"] == 1:
            assert not cache_path.exists()
            return FakeResponse({}, status_code=429)
        return FakeResponse({"ok": True}, status_code=200)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        payload = client.get_json(url)

    assert payload == {"ok": True}
    assert cache_path.exists()

    existing_payload = {"cached": "good"}
    _write_cache(cache_path, url, existing_payload)
    refresh_client = _make_client(cache_dir, cache_mode="refresh", max_retries=0)

    with mock.patch.object(
        tool.requests, "get", return_value=FakeResponse({}, status_code=429)
    ):
        payload = refresh_client.get_json(url)

    assert payload is None
    cached = json.loads(cache_path.read_text(encoding="utf-8"))
    assert cached["payload"] == existing_payload


def test_retry_backoff_and_logging(tmp_path, caplog):
    cache_dir = tmp_path / "cache"
    client = tool.CachedJsonClient(
        cache_dir,
        sleep=0,
        cache_mode="off",
        rate_limiter=None,
        max_retries=5,
        backoff_base=0.5,
        backoff_max=8.0,
    )
    responses = [
        FakeResponse({}, status_code=503),
        FakeResponse({}, status_code=503),
        FakeResponse({"ok": True}, status_code=200),
    ]
    waits = []

    def fake_get(url_arg, timeout=10):
        return responses.pop(0)

    with mock.patch.object(
        tool.requests, "get", side_effect=fake_get
    ) as mocked_get, mock.patch.object(
        tool.time, "sleep", side_effect=lambda value: waits.append(value)
    ):
        with caplog.at_level(logging.WARNING):
            payload = client.get_json("https://example.test/retry")

    assert payload == {"ok": True}
    assert mocked_get.call_count == 3
    assert len(waits) == 2
    assert waits[1] > waits[0]
    messages = [record.message for record in caplog.records]
    assert any("status=503" in message and "retry=1" in message for message in messages)
    assert any("status=503" in message and "retry=2" in message for message in messages)


def test_cache_modes_use_refresh_off(tmp_path):
    cache_dir = tmp_path / "cache"
    url = "https://example.test/cache"
    seeded_payload = {"seeded": True}

    use_client = _make_client(cache_dir, cache_mode="use")
    cache_path = use_client._cache_path(url)
    _write_cache(cache_path, url, seeded_payload)

    with mock.patch.object(tool.requests, "get") as mocked_get:
        payload = use_client.get_json(url)
    assert payload == seeded_payload
    mocked_get.assert_not_called()

    refresh_client = _make_client(cache_dir, cache_mode="refresh")
    new_payload = {"fresh": True}
    with mock.patch.object(
        tool.requests, "get", return_value=FakeResponse(new_payload, status_code=200)
    ) as mocked_get:
        payload = refresh_client.get_json(url)
    assert payload == new_payload
    assert mocked_get.call_count == 1
    cached = json.loads(cache_path.read_text(encoding="utf-8"))
    assert cached["payload"] == new_payload

    off_client = _make_client(cache_dir, cache_mode="off")
    with mock.patch.object(
        tool.requests, "get", return_value=FakeResponse({"off": True}, status_code=200)
    ) as mocked_get:
        payload = off_client.get_json(url)
    assert payload == {"off": True}
    assert mocked_get.call_count == 1
    cached = json.loads(cache_path.read_text(encoding="utf-8"))
    assert cached["payload"] == new_payload


def test_skip_cid_lookup_when_resolved(tmp_path):
    cache_dir = tmp_path / "cache"
    client = tool.CachedJsonClient(
        cache_dir,
        sleep=0,
        cache_mode="off",
        rate_limiter=None,
        max_retries=0,
        backoff_base=0.0,
        backoff_max=0.0,
    )
    row = {"pubchem_cid_resolved": "2244", "inchikey": "INCHIKEYX"}
    urls = []

    def fake_get(url_arg, timeout=10):
        urls.append(url_arg)
        if "/property/Title,IUPACName/JSON" in url_arg:
            return FakeResponse(
                {
                    "PropertyTable": {
                        "Properties": [{"Title": "Name", "IUPACName": "IUPAC"}]
                    }
                }
            )
        if "/synonyms/JSON" in url_arg:
            return FakeResponse(
                {"InformationList": {"Information": [{"Synonym": ["Name", "Alias"]}]}}
            )
        return FakeResponse({}, status_code=404)

    with mock.patch.object(tool.requests, "get", side_effect=fake_get):
        result = tool.resolve_pubchem(row, client)

    assert result.cid == "2244"
    assert any("/cid/2244/property/Title,IUPACName/JSON" in url for url in urls)
    assert any("/cid/2244/synonyms/JSON" in url for url in urls)
    assert not any("/inchikey/" in url for url in urls)
    assert not any("/smiles/" in url for url in urls)
