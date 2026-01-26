#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import re
import shutil
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence
from urllib.parse import quote, urlencode

import requests  # type: ignore[import-untyped]

FDA_PLACEHOLDER_RE = re.compile(r"^fda_\d+$", re.IGNORECASE)
IUPAC_PAREN_DIGIT_RE = re.compile(r"\(\d")
IUPAC_DIGIT_LOCANT_RE = re.compile(r"\b\d{1,3}[a-z]?\s*-\s*[A-Za-z]")
IUPAC_NOS_LOCANT_RE = re.compile(r"\bN-\b|\bO-\b|\bS-\b")
IUPAC_COMMA_LOCANT_RE = re.compile(r"\b\d,\d")
IUPAC_TOKEN_RE = re.compile(
    r"(azanium|ylazanium|methoxy|propyl|butyl|penta|hexa|cyclo|carboxylate|benzamide)",
    re.IGNORECASE,
)
IUPAC_LONG_PUNCT_RE = re.compile(r"[,;()]")


def _clean_text(value: Any) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def is_bad_display_name(name: str, pubchem_iupac_name: str | None = None) -> bool:
    text = _clean_text(name)
    if not text:
        return True
    if FDA_PLACEHOLDER_RE.match(text):
        return True
    if pubchem_iupac_name:
        if text.lower() == _clean_text(pubchem_iupac_name).lower():
            return True
    if IUPAC_PAREN_DIGIT_RE.search(text):
        return True
    if IUPAC_DIGIT_LOCANT_RE.search(text):
        return True
    if IUPAC_NOS_LOCANT_RE.search(text):
        return True
    if IUPAC_COMMA_LOCANT_RE.search(text):
        return True
    if IUPAC_TOKEN_RE.search(text):
        return True
    if len(text) > 60 and IUPAC_LONG_PUNCT_RE.search(text):
        return True
    return False


def is_good_display_name(name: str, pubchem_iupac_name: str | None = None) -> bool:
    return not is_bad_display_name(name, pubchem_iupac_name)


@dataclass(frozen=True)
class PubChemResult:
    cid: str | None = None
    title: str | None = None
    iupac_name: str | None = None
    synonyms: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class RxNormResult:
    rxcui: str | None = None
    name: str | None = None


class CachedJsonClient:
    def __init__(self, cache_dir: Path, sleep: float) -> None:
        self.cache_dir = cache_dir
        self.sleep = max(float(sleep), 0.0)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()

    def _cache_path(self, url: str) -> Path:
        digest = hashlib.sha256(url.encode("utf-8")).hexdigest()
        return self.cache_dir / f"{digest}.json"

    def get_json(self, url: str) -> Optional[dict[str, Any]]:
        cache_path = self._cache_path(url)
        with self._lock:
            if cache_path.exists():
                try:
                    cached_payload = json.loads(cache_path.read_text(encoding="utf-8"))
                    if isinstance(cached_payload, dict) and "payload" in cached_payload:
                        return cached_payload.get("payload")
                    return cached_payload
                except Exception:
                    pass

        if self.sleep:
            time.sleep(self.sleep)

        response_payload: Optional[dict[str, Any]] = None
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                response_payload = response.json()
        except Exception:
            response_payload = None

        with self._lock:
            if not cache_path.exists():
                try:
                    tmp_path = cache_path.with_suffix(cache_path.suffix + ".tmp")
                    tmp_path.write_text(
                        json.dumps(
                            {"url": url, "payload": response_payload}, ensure_ascii=True
                        ),
                        encoding="utf-8",
                    )
                    tmp_path.replace(cache_path)
                except Exception:
                    pass

        return response_payload


def _resolution_key(row: dict[str, Any]) -> tuple[str | None, str | None]:
    for key in ("inchikey", "remark_inchikey", "smiles_neutral", "smiles", "remark_smiles"):
        value = _clean_text(row.get(key))
        if value:
            return key, value
    return None, None


def _parse_pubchem_cid(payload: Optional[dict[str, Any]]) -> Optional[str]:
    if not payload:
        return None
    try:
        cids = payload["IdentifierList"]["CID"]
    except Exception:
        return None
    if not cids:
        return None
    return str(cids[0])


def _parse_pubchem_properties(
    payload: Optional[dict[str, Any]],
) -> tuple[Optional[str], Optional[str]]:
    if not payload:
        return None, None
    try:
        props = payload["PropertyTable"]["Properties"]
        if not props:
            return None, None
        entry = props[0]
        return entry.get("Title"), entry.get("IUPACName")
    except Exception:
        return None, None


def _parse_pubchem_synonyms(payload: Optional[dict[str, Any]]) -> list[str]:
    if not payload:
        return []
    try:
        info = payload["InformationList"]["Information"]
        if not info:
            return []
        entry = info[0]
        syns = entry.get("Synonym", [])
        return [str(s) for s in syns if str(s).strip()]
    except Exception:
        return []


def resolve_pubchem(row: dict[str, Any], client: CachedJsonClient) -> PubChemResult:
    key_type, key_value = _resolution_key(row)
    if not key_type or not key_value:
        return PubChemResult()

    if key_type in ("inchikey", "remark_inchikey"):
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
            f"{quote(key_value, safe='')}/cids/JSON"
        )
    else:
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
            f"{quote(key_value, safe='')}/cids/JSON"
        )

    cid = _parse_pubchem_cid(client.get_json(url))
    if not cid:
        return PubChemResult()

    props_url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{cid}/property/Title,IUPACName/JSON"
    )
    syn_url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{cid}/synonyms/JSON"
    )
    title, iupac_name = _parse_pubchem_properties(client.get_json(props_url))
    synonyms = _parse_pubchem_synonyms(client.get_json(syn_url))
    return PubChemResult(cid=cid, title=title, iupac_name=iupac_name, synonyms=synonyms)


def _parse_rxnorm_candidate(payload: Optional[dict[str, Any]]) -> Optional[str]:
    if not payload:
        return None
    try:
        candidates = payload["approximateGroup"]["candidate"]
    except Exception:
        return None
    if not candidates:
        return None

    def score_of(item: dict[str, Any]) -> int:
        try:
            return int(item.get("score", 0))
        except Exception:
            return 0

    best = max(candidates, key=score_of)
    return best.get("rxcui")


def _parse_rxnorm_name(payload: Optional[dict[str, Any]]) -> Optional[str]:
    if not payload:
        return None
    try:
        props = payload["properties"]
        return props.get("name")
    except Exception:
        return None


def resolve_rxnorm(term: str, client: CachedJsonClient) -> RxNormResult:
    query = _clean_text(term)
    if not query:
        return RxNormResult()
    url = f"https://rxnav.nlm.nih.gov/REST/approximateTerm.json?{urlencode({'term': query, 'maxEntries': 1})}"
    rxcui = _parse_rxnorm_candidate(client.get_json(url))
    if not rxcui:
        return RxNormResult()
    props_url = f"https://rxnav.nlm.nih.gov/REST/rxcui/{rxcui}/properties.json"
    name = _parse_rxnorm_name(client.get_json(props_url))
    return RxNormResult(rxcui=rxcui, name=name)


def _maybe_set_field(row: dict[str, Any], field: str, value: Optional[str]) -> bool:
    if field not in row:
        return False
    if _clean_text(row.get(field)):
        return False
    if not _clean_text(value):
        return False
    row[field] = value
    return True


def _pick_first(*values: Optional[str]) -> str:
    for value in values:
        text = _clean_text(value)
        if text:
            return text
    return ""


def _split_synonyms(raw: str) -> list[str]:
    if not raw:
        return []
    parts = [p.strip() for p in re.split(r"[;|]", raw) if p.strip()]
    return parts


def _choose_pubchem_synonym(
    candidates: Iterable[str], pubchem_iupac_name: str | None
) -> str:
    for candidate in candidates:
        if is_good_display_name(candidate, pubchem_iupac_name):
            return candidate
    return ""


def _normalize_numeric(text: str) -> str:
    raw = _clean_text(text)
    if not raw:
        return ""
    try:
        return str(int(float(raw)))
    except Exception:
        return ""


def _fallback_display_name(row: dict[str, Any], row_index: int) -> str:
    file_num = _normalize_numeric(_clean_text(row.get("file_num")))
    if file_num:
        return f"UNK_{file_num}"
    sdf_index = _normalize_numeric(_clean_text(row.get("sdf_index")))
    if sdf_index:
        return f"UNK_{sdf_index}"
    inchikey = _pick_first(row.get("inchikey"), row.get("remark_inchikey"))
    if inchikey:
        return f"UNK_{inchikey[:8]}"
    return f"UNK_ROW{row_index}"


def _row_identifier(row: dict[str, Any], row_index: int) -> str:
    return _pick_first(
        row.get("sdf_title"),
        row.get("path"),
        row.get("file_num"),
        row.get("sdf_index"),
        f"row_{row_index}",
    )


def _select_display_name(
    row: dict[str, Any],
    pubchem: PubChemResult,
    row_index: int,
) -> tuple[str, bool]:
    pubchem_iupac = _pick_first(row.get("pubchem_iupac_name"), pubchem.iupac_name)
    synonyms = pubchem.synonyms
    if not synonyms:
        synonyms = _split_synonyms(_clean_text(row.get("pubchem_synonyms")))
    best_synonym = _choose_pubchem_synonym(synonyms, pubchem_iupac)

    candidates = [
        _clean_text(row.get("rxnorm_generic_name")),
        _clean_text(row.get("drugcentral_generic_name")),
        _clean_text(row.get("generic_name")),
        _clean_text(row.get("pubchem_record_title")),
        _clean_text(row.get("pubchem_name")),
        best_synonym,
    ]
    for candidate in candidates:
        if candidate and is_good_display_name(candidate, pubchem_iupac):
            return candidate, False

    return _fallback_display_name(row, row_index), True


@dataclass
class Summary:
    total_rows: int = 0
    rows_evaluated: int = 0
    rows_updated: int = 0
    pubchem_resolved: int = 0
    rxnorm_resolved: int = 0
    fallbacks_used: int = 0


@dataclass
class RowResult:
    row_index: int
    row: dict[str, Any]
    row_evaluated: bool
    row_updated: bool
    pubchem_resolved: bool
    rxnorm_resolved: bool
    used_fallback: bool
    final_bad: bool
    row_identifier: str


def _process_row(
    row_index: int,
    row: dict[str, Any],
    only_fix_bad_display_names: bool,
    pubchem_client: CachedJsonClient,
    rxnorm_client: CachedJsonClient,
) -> RowResult:
    current_display = _clean_text(row.get("display_name"))
    pubchem_iupac = _clean_text(row.get("pubchem_iupac_name"))
    display_bad = is_bad_display_name(current_display, pubchem_iupac)
    evaluate = (not only_fix_bad_display_names) or display_bad
    row_updated = False
    used_fallback = False
    row_evaluated = False
    pubchem_resolved = False
    rxnorm_resolved = False

    pubchem = PubChemResult()
    if evaluate:
        row_evaluated = True
        pubchem = resolve_pubchem(row, pubchem_client)
        if pubchem.cid:
            pubchem_resolved = True
        row_updated |= _maybe_set_field(row, "pubchem_cid_resolved", pubchem.cid)
        row_updated |= _maybe_set_field(row, "pubchem_record_title", pubchem.title)
        row_updated |= _maybe_set_field(row, "pubchem_name", pubchem.title)
        row_updated |= _maybe_set_field(row, "pubchem_iupac_name", pubchem.iupac_name)
        if pubchem.synonyms:
            syn_value = "; ".join(pubchem.synonyms)
            row_updated |= _maybe_set_field(row, "pubchem_synonyms", syn_value)

        rxnorm_name = _clean_text(row.get("rxnorm_generic_name"))
        if not rxnorm_name:
            candidate = _pick_first(
                pubchem.title,
                row.get("pubchem_record_title"),
                row.get("pubchem_name"),
                row.get("generic_name"),
                current_display,
            )
            rxnorm = resolve_rxnorm(candidate, rxnorm_client)
            if rxnorm.rxcui or rxnorm.name:
                rxnorm_resolved = True
            row_updated |= _maybe_set_field(row, "rxnorm_generic_name", rxnorm.name)
            row_updated |= _maybe_set_field(row, "rxnorm_rxcui", rxnorm.rxcui)

        pubchem_iupac = _clean_text(row.get("pubchem_iupac_name"))
        if is_bad_display_name(current_display, pubchem_iupac):
            new_display, used_fallback = _select_display_name(row, pubchem, row_index)
            if new_display and new_display != current_display:
                row["display_name"] = new_display
                row_updated = True

    final_display = _clean_text(row.get("display_name"))
    final_bad = is_bad_display_name(
        final_display, _clean_text(row.get("pubchem_iupac_name"))
    )
    row_identifier = _row_identifier(row, row_index) if (used_fallback or final_bad) else ""
    return RowResult(
        row_index=row_index,
        row=row,
        row_evaluated=row_evaluated,
        row_updated=row_updated,
        pubchem_resolved=pubchem_resolved,
        rxnorm_resolved=rxnorm_resolved,
        used_fallback=used_fallback,
        final_bad=final_bad,
        row_identifier=row_identifier,
    )


def _process_row_item(
    item: tuple[int, dict[str, Any]],
    only_fix_bad_display_names: bool,
    pubchem_client: CachedJsonClient,
    rxnorm_client: CachedJsonClient,
) -> RowResult:
    row_index, row = item
    return _process_row(
        row_index,
        row,
        only_fix_bad_display_names=only_fix_bad_display_names,
        pubchem_client=pubchem_client,
        rxnorm_client=rxnorm_client,
    )


def fill_mapping_names(
    in_csv: Path,
    out_csv: Path,
    *,
    only_fix_bad_display_names: bool,
    fail_if_unresolved: bool,
    report_path: Path,
    cache_dir: Path,
    sleep: float,
    max_workers: int,
) -> int:
    pubchem_client = CachedJsonClient(cache_dir / "pubchem", sleep)
    rxnorm_client = CachedJsonClient(cache_dir / "rxnorm", sleep)

    unresolved: list[str] = []
    has_bad = False
    summary = Summary()

    in_place = in_csv.resolve() == out_csv.resolve()
    temp_path = out_csv.with_suffix(out_csv.suffix + ".tmp") if in_place else None
    write_path = temp_path if temp_path else out_csv
    write_path.parent.mkdir(parents=True, exist_ok=True)

    if max_workers < 1:
        raise ValueError("--max_workers must be >= 1.")

    with in_csv.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        if not fieldnames:
            raise ValueError("Input CSV has no headers.")

        with write_path.open("w", newline="", encoding="utf-8") as out_handle:
            writer = csv.DictWriter(out_handle, fieldnames=fieldnames)
            writer.writeheader()
            process_row = partial(
                _process_row_item,
                only_fix_bad_display_names=only_fix_bad_display_names,
                pubchem_client=pubchem_client,
                rxnorm_client=rxnorm_client,
            )

            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                for result in executor.map(process_row, enumerate(reader, start=1)):
                    summary.total_rows += 1
                    if result.row_evaluated:
                        summary.rows_evaluated += 1
                    if result.pubchem_resolved:
                        summary.pubchem_resolved += 1
                    if result.rxnorm_resolved:
                        summary.rxnorm_resolved += 1
                    if result.used_fallback:
                        summary.fallbacks_used += 1
                    if result.row_updated:
                        summary.rows_updated += 1
                    if result.row_identifier:
                        unresolved.append(result.row_identifier)
                    if result.final_bad:
                        has_bad = True

                    writer.writerow(
                        {name: result.row.get(name, "") for name in fieldnames}
                    )

    if temp_path:
        temp_path.replace(out_csv)

    report_path.parent.mkdir(parents=True, exist_ok=True)
    ordered_unresolved = list(dict.fromkeys(unresolved))
    report_path.write_text("\n".join(ordered_unresolved) + ("\n" if ordered_unresolved else ""), encoding="utf-8")

    print(f"total rows: {summary.total_rows}")
    print(f"rows evaluated: {summary.rows_evaluated}")
    print(f"rows updated: {summary.rows_updated}")
    print(f"pubchem resolved: {summary.pubchem_resolved}")
    print(f"rxnorm resolved: {summary.rxnorm_resolved}")
    print(f"fallbacks used: {summary.fallbacks_used}")

    if fail_if_unresolved and has_bad:
        return 2
    return 0


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fill FDA mapping display names from PubChem/RxNorm without overwriting inputs."
    )
    parser.add_argument("--in_csv", required=True, help="Input FDA mapping CSV.")
    parser.add_argument("--out_csv", help="Output CSV path (required unless --inplace).")
    parser.add_argument(
        "--inplace",
        action="store_true",
        help="Overwrite --in_csv after creating a timestamped .bak backup.",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.25,
        help="Sleep between external requests (only when uncached).",
    )
    parser.add_argument(
        "--only_fix_bad_display_names",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Only update rows with placeholder/IUPAC-like display names (default: true).",
    )
    parser.add_argument(
        "--fail_if_unresolved",
        action="store_true",
        help="Exit with code 2 if any row still has an invalid display_name.",
    )
    parser.add_argument(
        "--report_path",
        help="Path for unresolved report (default: <out_csv>.unresolved.txt).",
    )
    parser.add_argument(
        "--cache_dir",
        default=".cache/fda_name_fill/",
        help="Directory for JSON response caches.",
    )
    parser.add_argument(
        "--max_workers",
        type=int,
        default=min(32, os.cpu_count() or 4),
        help="Maximum worker threads for parallel resolution (default: min(32, cpu_count)).",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_arg_parser()
    args = parser.parse_args(argv)

    in_path = Path(args.in_csv)
    if not in_path.exists():
        parser.error(f"--in_csv not found: {in_path}")

    if args.inplace:
        if args.out_csv:
            out_path = Path(args.out_csv)
            if out_path.resolve() != in_path.resolve():
                parser.error("--inplace requires --out_csv to match --in_csv when provided.")
        out_path = in_path
        timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_path = in_path.with_name(f"{in_path.name}.{timestamp}.bak")
        shutil.copy2(in_path, backup_path)
    else:
        if not args.out_csv:
            parser.error("--out_csv is required unless --inplace is set.")
        out_path = Path(args.out_csv)
        if out_path.resolve() == in_path.resolve():
            parser.error("--out_csv must differ from --in_csv unless --inplace is set.")

    report_path = Path(args.report_path) if args.report_path else Path(
        f"{out_path}.unresolved.txt"
    )
    cache_dir = Path(args.cache_dir)

    try:
        return fill_mapping_names(
            in_path,
            out_path,
            only_fix_bad_display_names=bool(args.only_fix_bad_display_names),
            fail_if_unresolved=bool(args.fail_if_unresolved),
            report_path=report_path,
            cache_dir=cache_dir,
            sleep=float(args.sleep),
            max_workers=int(args.max_workers),
        )
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
