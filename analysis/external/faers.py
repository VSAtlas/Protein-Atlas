from __future__ import annotations

import hashlib
import json
import math
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.mapping import normalize_columns
from analysis.external.source_tables import read_source_table


OPENFDA_EVENT_URL = "https://api.fda.gov/drug/event.json"

FAERS_PAIR_ALIASES = {
    "drug_id": [
        "drug_id",
        "drug",
        "drug_name",
        "exposureName",
        "exposure_name",
        "targetName",
        "target_name",
    ],
    "adr_id": ["adr_id", "adr", "event", "outcomeName", "outcome_name", "meddra_term", "pt"],
    "ground_truth": ["groundTruth", "ground_truth", "label", "negative_control", "is_negative"],
    "source": ["source", "dataset"],
}


def _clean_term(value: Any) -> str:
    text = "" if pd.isna(value) else str(value)
    text = text.replace("_", " ").replace("#", " ")
    text = " ".join(text.split()).strip()
    return text


def _clean_ohdsi_outcome(value: Any) -> str:
    text = _clean_term(value)
    for prefix in ("OMOP ", "HOI "):
        if text.upper().startswith(prefix.strip()):
            text = text[len(prefix) :]
    parts = [part for part in text.split() if not part.isdigit()]
    return " ".join(parts).strip()


def _openfda_phrase(field: str, value: Any) -> str:
    text = _clean_ohdsi_outcome(value)
    text = "".join(ch if ch.isalnum() or ch.isspace() or ch in "-/" else " " for ch in text)
    text = " ".join(text.split()).strip(" -/")
    if not text:
        return ""
    return f'{field}:"{text}"'


def _drug_search(drug: str) -> str:
    terms = [
        _openfda_phrase("patient.drug.openfda.generic_name", drug),
        _openfda_phrase("patient.drug.openfda.brand_name", drug),
        _openfda_phrase("patient.drug.medicinalproduct", drug),
    ]
    terms = [term for term in terms if term]
    return "(" + " OR ".join(terms) + ")" if terms else ""


def _event_search(event: str) -> str:
    return _openfda_phrase("patient.reaction.reactionmeddrapt", event)


def _cache_key(search: str) -> str:
    return hashlib.sha256(search.encode("utf-8")).hexdigest()[:24]


def _request_total(
    search: str,
    *,
    cache_dir: Path,
    api_key: str | None = None,
    sleep_sec: float = 0.2,
    timeout: int = 60,
    retries: int = 4,
) -> tuple[int, str]:
    if not search:
        return 0, "empty_query"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{_cache_key(search)}.json"
    if cache_path.exists():
        payload = json.loads(cache_path.read_text(encoding="utf-8"))
        return int(payload.get("total") or 0), str(payload.get("status") or "cached")
    params = {"search": search, "limit": "1"}
    if api_key:
        params["api_key"] = api_key
    url = f"{OPENFDA_EVENT_URL}?{urllib.parse.urlencode(params)}"
    status = "failed"
    total = 0
    last_error = ""
    for attempt in range(max(1, retries)):
        if sleep_sec > 0:
            time.sleep(sleep_sec)
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "AtlasAnalysis/1.0"})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                payload = json.loads(response.read().decode("utf-8"))
            total = int(((payload.get("meta") or {}).get("results") or {}).get("total") or 0)
            status = "ok"
            break
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                status = "no_results"
                total = 0
                break
            last_error = f"HTTP {exc.code}"
            if exc.code not in {429, 500, 502, 503, 504}:
                status = "error"
                break
            retry_after = exc.headers.get("Retry-After")
            try:
                delay = float(retry_after) if retry_after else 0.0
            except ValueError:
                delay = 0.0
            time.sleep(max(delay, min(60.0, 2.0 * float(attempt + 1))))
        except Exception as exc:
            last_error = exc.__class__.__name__
            time.sleep(min(30.0, 1.0 + float(attempt)))
    cache_path.write_text(
        json.dumps(
            {
                "search": search,
                "url": url,
                "total": total,
                "status": status,
                "last_error": last_error,
                "source": OPENFDA_EVENT_URL,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return total, status


def _ror_stats(pair_reports: int, drug_reports: int, event_reports: int, total_reports: int) -> dict[str, float]:
    a = max(0.0, float(pair_reports))
    b = max(0.0, float(drug_reports) - a)
    c = max(0.0, float(event_reports) - a)
    d = max(0.0, float(total_reports) - a - b - c)
    ac, bc, cc, dc = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    ror = (ac / bc) / (cc / dc)
    se = math.sqrt((1.0 / ac) + (1.0 / bc) + (1.0 / cc) + (1.0 / dc))
    ror_lower = math.exp(math.log(ror) - 1.96 * se)
    prr = (ac / (ac + bc)) / (cc / (cc + dc))
    expected = ((ac + bc) * (ac + cc)) / max(1.0, ac + bc + cc + dc)
    ebgm = ac / max(expected, 1e-12)
    ebgm05 = math.exp(math.log(ebgm) - 1.645 * se)
    return {
        "ror": ror,
        "ror_lower": ror_lower,
        "prr": prr,
        "ebgm": ebgm,
        "ebgm05": ebgm05,
    }


def load_drug_adr_pairs(paths: list[str | Path], *, negative_only: bool = False) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in paths:
        raw = read_source_table(path)
        df = normalize_columns(raw, FAERS_PAIR_ALIASES)
        if negative_only:
            truth = pd.to_numeric(df["ground_truth"], errors="coerce")
            if truth.notna().any():
                df = df[truth.eq(0)].copy()
            elif "negative" in str(path).lower():
                df = df.copy()
            else:
                df = df.iloc[0:0].copy()
        df["drug_id"] = df["drug_id"].map(_clean_term)
        df["adr_id"] = df["adr_id"].map(_clean_ohdsi_outcome)
        df = df[df["drug_id"].ne("") & df["adr_id"].ne("")].copy()
        df["source"] = df["source"].fillna(Path(path).stem)
        frames.append(df[["drug_id", "adr_id", "source"]])
    if not frames:
        return pd.DataFrame(columns=["drug_id", "adr_id", "source"])
    return pd.concat(frames, ignore_index=True).drop_duplicates(["drug_id", "adr_id"])


def build_faers_disproportionality_table(
    pair_paths: list[str | Path],
    out_path: str | Path,
    *,
    cache_dir: str | Path,
    negative_only: bool = False,
    api_key: str | None = None,
    sleep_sec: float = 0.2,
    max_pairs: int | None = None,
) -> pd.DataFrame:
    pairs = load_drug_adr_pairs(pair_paths, negative_only=negative_only)
    if max_pairs is not None:
        pairs = pairs.head(max_pairs).copy()
    cache = Path(cache_dir)
    total_reports, total_status = _request_total(
        "receivedate:[20040101 TO 30000101]",
        cache_dir=cache,
        api_key=api_key,
        sleep_sec=sleep_sec,
    )
    drug_cache: dict[str, tuple[int, str]] = {}
    event_cache: dict[str, tuple[int, str]] = {}
    rows: list[dict[str, Any]] = []
    for row in pairs.itertuples(index=False):
        drug = str(row.drug_id)
        event = str(row.adr_id)
        drug_query = _drug_search(drug)
        event_query = _event_search(event)
        pair_query = f"{drug_query} AND {event_query}" if drug_query and event_query else ""
        if drug not in drug_cache:
            drug_cache[drug] = _request_total(
                drug_query, cache_dir=cache, api_key=api_key, sleep_sec=sleep_sec
            )
        if event not in event_cache:
            event_cache[event] = _request_total(
                event_query, cache_dir=cache, api_key=api_key, sleep_sec=sleep_sec
            )
        pair_reports, pair_status = _request_total(
            pair_query, cache_dir=cache, api_key=api_key, sleep_sec=sleep_sec
        )
        drug_reports, drug_status = drug_cache[drug]
        event_reports, event_status = event_cache[event]
        stats = _ror_stats(pair_reports, drug_reports, event_reports, total_reports)
        rows.append(
            {
                "drug_id": drug,
                "adr_id": event,
                "source": "openFDA_FAERS",
                "pair_source": row.source,
                "pair_reports": pair_reports,
                "drug_reports": drug_reports,
                "event_reports": event_reports,
                "total_reports": total_reports,
                "pair_query_status": pair_status,
                "drug_query_status": drug_status,
                "event_query_status": event_status,
                "total_query_status": total_status,
                **stats,
            }
        )
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(rows)
    result.to_csv(out, index=False)
    out.with_suffix(".manifest.json").write_text(
        json.dumps(
            {
                "source": OPENFDA_EVENT_URL,
                "pair_paths": [str(path) for path in pair_paths],
                "negative_only": negative_only,
                "rows": int(len(result)),
                "cache_dir": str(cache),
                "total_reports": int(total_reports),
                "sleep_sec": sleep_sec,
                "max_pairs": max_pairs,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return result
