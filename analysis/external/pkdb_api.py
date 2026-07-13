from __future__ import annotations

from datetime import datetime, timezone
import io
import json
from pathlib import Path
from typing import Any
import zipfile

import requests  # type: ignore[import-untyped]


PKDB_API_BASE = "https://pk-db.com/api/v1"
PKDB_DOCUMENTED_EXAMPLE = "Abernethy1982"


def _response_count(response: requests.Response) -> int:
    try:
        payload = response.json()
    except ValueError:
        return 0
    data = payload.get("data") if isinstance(payload, dict) else None
    if not isinstance(data, dict):
        return 0
    try:
        return int(data.get("count") or 0)
    except (TypeError, ValueError):
        return 0

def _response_rows(response: requests.Response) -> list[dict[str, Any]]:
    try:
        payload = response.json()
    except ValueError:
        return []
    data = payload.get("data") if isinstance(payload, dict) else None
    rows = data.get("data") if isinstance(data, dict) else None
    return [row for row in rows or [] if isinstance(row, dict)]


def _zip_member_size(content: bytes, member: str) -> int:
    try:
        with zipfile.ZipFile(io.BytesIO(content)) as archive:
            return archive.getinfo(member).file_size
    except (KeyError, OSError, zipfile.BadZipFile):
        return 0


def probe_pkdb_api(
    *,
    out_dir: str | Path,
    timeout: int = 180,
    session: requests.Session | None = None,
) -> dict[str, Any]:
    """Check whether PK-DB's advertised output rows are actually retrievable.

    PK-DB's filter endpoint can report nonzero output counts while both the
    outputs endpoint and exported outputs.csv are empty. That state is a
    source outage, not evidence that the queried study has no PK measurements.
    """

    client = session or requests.Session()
    payload: dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "api_base": PKDB_API_BASE,
        "documented_example": PKDB_DOCUMENTED_EXAMPLE,
        "documentation_url": "https://pk-db.com/api/v1/swagger/",
        "example_notebook_url": (
            "https://github.com/matthiaskoenig/pkdb/blob/develop/docs/"
            "pkdb_api.ipynb"
        ),
        "status": "unavailable",
    }
    try:
        statistics = client.get(f"{PKDB_API_BASE}/statistics/", timeout=timeout)
        statistics.raise_for_status()
        statistics_payload = statistics.json()
        payload["statistics_http_status"] = statistics.status_code
        payload["statistics_output_count"] = int(
            statistics_payload.get("output_count") or 0
        )

        filter_response = client.get(
            f"{PKDB_API_BASE}/filter/",
            params={"studies__name": PKDB_DOCUMENTED_EXAMPLE},
            timeout=timeout,
        )
        filter_response.raise_for_status()
        filter_payload = filter_response.json()
        query_uuid = str(filter_payload.get("uuid") or "")
        advertised_outputs = int(filter_payload.get("outputs") or 0)
        payload["filter_http_status"] = filter_response.status_code
        payload["filter_uuid"] = query_uuid
        payload["filter_advertised_outputs"] = advertised_outputs
        study_response = client.get(
            f"{PKDB_API_BASE}/studies/",
            params={"name": PKDB_DOCUMENTED_EXAMPLE, "format": "json"},
            timeout=timeout,
        )
        study_response.raise_for_status()
        study_rows = _response_rows(study_response)
        study_row = study_rows[0] if study_rows else {}
        outputset = study_row.get("outputset")
        study_output_ids = (
            outputset.get("outputs") if isinstance(outputset, dict) else []
        ) or []
        payload["study_endpoint_count"] = _response_count(study_response)
        payload["study_reported_output_count"] = int(
            study_row.get("output_count") or 0
        )
        payload["study_embedded_output_id_count"] = len(study_output_ids)

        alpha_response = client.get(
            "https://alpha.pk-db.com/api/v1/outputs/",
            params={"page_size": 1, "format": "json"},
            timeout=timeout,
        )
        alpha_response.raise_for_status()
        payload["alpha_outputs_endpoint_count"] = _response_count(alpha_response)

        outputs_response = client.get(
            f"{PKDB_API_BASE}/outputs/",
            params={"uuid": query_uuid, "page_size": 5, "format": "json"},
            timeout=timeout,
        )
        outputs_response.raise_for_status()
        endpoint_outputs = _response_count(outputs_response)
        payload["outputs_http_status"] = outputs_response.status_code
        payload["outputs_endpoint_count"] = endpoint_outputs

        archive_response = client.get(
            f"{PKDB_API_BASE}/filter/",
            params={
                "studies__name": PKDB_DOCUMENTED_EXAMPLE,
                "download": "true",
                "concise": "true",
            },
            timeout=timeout,
        )
        archive_response.raise_for_status()
        archive_output_bytes = _zip_member_size(
            archive_response.content,
            "outputs.csv",
        )
        payload["archive_http_status"] = archive_response.status_code
        payload["archive_bytes"] = len(archive_response.content)
        payload["archive_outputs_csv_bytes"] = archive_output_bytes

        outputs_retrievable = endpoint_outputs > 0 or archive_output_bytes > 3
        if advertised_outputs > 0 and not outputs_retrievable:
            payload["status"] = "server_inconsistent"
            payload["reason"] = (
                "the PostgreSQL-backed study endpoint exposes output counts and "
                "IDs, but Elasticsearch-backed /outputs/ and generated "
                "outputs.csv are empty on both production and alpha; numeric "
                "PK values are not publicly recoverable until the source index is rebuilt"
            )
        elif outputs_retrievable:
            payload["status"] = "available"
            payload["reason"] = "output measurements are retrievable"
        else:
            payload["status"] = "no_documented_example_outputs"
            payload["reason"] = "documented example unexpectedly has no outputs"
    except (OSError, requests.RequestException, TypeError, ValueError) as exc:
        payload["status"] = "request_failed"
        payload["reason"] = f"{type(exc).__name__}: {exc}"

    destination = Path(out_dir)
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "pkdb_api_health.json").write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    return payload
