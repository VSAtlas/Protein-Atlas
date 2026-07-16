from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import re
import time
from typing import Any, Iterable, Mapping
from urllib.parse import urlparse

import pandas as pd
import requests  # type: ignore[import-untyped]

from analysis.external.pk_context import (
    PK_CONTEXT_COLUMNS,
    finalize_context,
    normalize_key,
)


PKDB_API_BASE = "https://pk-db.com/api/v1"
PKDB_TERMS_URL = "https://github.com/matthiaskoenig/pkdb/blob/develop/TERMS_OF_USE.md"
SCHEMA_VERSION = "atlas_pkdb_source_tsv_recovery_v3"

IV_CLEARANCE_CONTEXTS = {
    "systemic_clearance",
    "plasma_clearance",
    "blood_clearance",
    "total_body_clearance",
}
APPARENT_CLEARANCE_CONTEXTS = {"apparent_oral_clearance"}
SOURCE_RIGHTS_COLUMNS = [
    "study_sid",
    "study_reference",
    "source_file_url",
    "training_allowed",
    "rights_reference",
    "rights_basis",
]

DETAIL_COLUMNS = [
    "pkdb_context_id",
    "drug_id",
    "drug_name",
    "inchikey",
    "rdkit_mol_wt",
    "study_sid",
    "study_name",
    "study_reference",
    "study_reference_type",
    "study_reference_id",
    "publication_year",
    "study_date",
    "study_licence",
    "source_file_url",
    "source_file_sha256",
    "source_row_number",
    "measurement_type",
    "measurement_context",
    "statistic_type",
    "raw_value",
    "raw_unit",
    "normalized_value",
    "normalized_unit",
    "normalization_method",
    "matrix",
    "method",
    "group_name",
    "population",
    "species",
    "reported_interventions",
    "intervention_name",
    "comparator_intervention_name",
    "dose_value",
    "dose_unit",
    "route",
    "reference_route",
    "comparator_route",
    "formulation",
    "regimen",
    "endpoint_link_status",
    "endpoint_evidence_ready",
    "context_status",
    "semantic_ready",
    "license_allows_ml_training",
    "source_rights_status",
    "source_rights_reference",
    "source_rights_basis",
    "source_rights_allow_ml_training",
    "training_allowed",
    "model_ready",
]


@dataclass(frozen=True)
class DrugIdentity:
    drug_id: str
    drug_name: str
    inchikey: str
    rdkit_mol_wt: float | None = None


def _text(value: Any) -> str:
    if value is None or (not isinstance(value, (list, tuple, dict)) and pd.isna(value)):
        return ""
    value_text = str(value).strip()
    return (
        "" if value_text.casefold() in {"", "na", "nan", "none", "null"} else value_text
    )


def _numeric(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if pd.notna(number) else None


def _truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return normalize_key(value) in {"1", "true", "yes", "y"}


def _first_text(row: Mapping[str, Any], columns: Iterable[str]) -> str:
    for column in columns:
        value = _text(row.get(column))
        if value:
            return value
    return ""


def _positive_finite_number(value: Any) -> float | None:
    number = _numeric(value)
    if number is None or not math.isfinite(number) or number <= 0:
        return None
    return number


def _valid_publication_year(value: Any) -> int | None:
    text = _text(value)
    match = re.search(r"(?<!\d)(1[6-9]\d{2}|20\d{2})(?!\d)", text)
    if match is None:
        return None
    year = int(match.group(1))
    return year if year <= datetime.now(timezone.utc).year + 1 else None


def _study_reference_metadata(
    study: Mapping[str, Any],
) -> tuple[str, str, str, int | None]:
    """Return a compact reference plus one structured identifier and publication year."""

    reference = study.get("reference")
    reference_type = ""
    reference_id = ""
    reference_text = ""
    publication_year: int | None = None

    if isinstance(reference, Mapping):
        pmid = _text(reference.get("pmid"))
        doi = _text(reference.get("doi"))
        if doi:
            doi = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", doi, flags=re.IGNORECASE)
        sid = _text(reference.get("sid"))
        name = _text(reference.get("name"))
        if pmid:
            reference_type, reference_id = "pmid", pmid
            reference_text = f"PMID:{pmid}"
        elif doi:
            reference_type, reference_id = "doi", doi
            reference_text = f"DOI:{doi}"
        elif sid:
            reference_type, reference_id = "sid", sid
            reference_text = f"SID:{sid}"
        elif name:
            reference_type, reference_id = "name", name
            reference_text = f"NAME:{name}"
        publication_year = _valid_publication_year(
            _first_text(reference, ("date", "publication_date", "year"))
        )
    else:
        reference_text = _text(reference)
        typed = re.fullmatch(
            r"\s*(pubmed|pmid|doi|sid)\s*[|:]\s*(\S.*?)\s*",
            reference_text,
            flags=re.IGNORECASE,
        )
        if typed:
            raw_type, reference_id = typed.groups()
            reference_type = (
                "pmid" if raw_type.casefold() == "pubmed" else raw_type.casefold()
            )
            reference_text = f"{reference_type.upper()}:{reference_id}"
        elif reference_text:
            reference_id = reference_text

    if publication_year is None:
        publication_year = _valid_publication_year(
            _first_text(
                study,
                ("reference_date", "publication_date", "publication_year"),
            )
        )
    return reference_text, reference_type, reference_id, publication_year


def normalize_source_rights(source_rights: pd.DataFrame | None) -> pd.DataFrame:
    """Validate explicit original-source rights evidence for fail-closed recovery."""

    if source_rights is None or source_rights.empty:
        return pd.DataFrame(columns=SOURCE_RIGHTS_COLUMNS)
    rights = source_rights.copy()
    rights.columns = [str(column).strip().casefold() for column in rights.columns]
    required = {"study_sid", "training_allowed", "rights_reference"}
    missing = sorted(required.difference(rights.columns))
    if missing:
        raise ValueError(
            "source-rights table is missing required columns: " + ", ".join(missing)
        )
    for column in SOURCE_RIGHTS_COLUMNS:
        if column not in rights:
            rights[column] = ""
    return rights.loc[:, SOURCE_RIGHTS_COLUMNS].copy()


def _source_rights_decision(
    *,
    study: Mapping[str, Any],
    source_url: str,
    source_rights: pd.DataFrame,
) -> tuple[bool, str, str, str]:
    if source_rights.empty:
        return False, "source_level_training_rights_not_provided", "", ""

    sid = _text(study.get("sid"))
    matches = source_rights.loc[
        source_rights["study_sid"].fillna("").astype(str).str.strip().eq(sid)
    ].copy()
    if matches.empty:
        return False, "source_level_training_rights_not_verified", "", ""

    rights_urls = matches["source_file_url"].fillna("").astype(str).str.strip()
    matches = matches.loc[rights_urls.eq("") | rights_urls.eq(source_url)].copy()
    study_reference = normalize_key(_study_reference_metadata(study)[0])
    rights_references = matches["study_reference"].map(normalize_key)
    matches = matches.loc[
        rights_references.eq("") | rights_references.eq(study_reference)
    ].copy()
    if matches.empty:
        return False, "source_level_training_rights_not_verified_for_source", "", ""

    allowed = matches["training_allowed"].map(_truthy)
    if allowed.any() and not allowed.all():
        return False, "source_level_training_rights_conflict", "", ""
    if not allowed.any():
        return False, "source_level_training_rights_denied", "", ""

    evidence = [
        _text(value)
        for value in matches.loc[allowed, "rights_reference"]
        if _text(value)
    ]
    if not evidence:
        return False, "source_level_training_rights_evidence_missing", "", ""
    bases = [
        _text(value) for value in matches.loc[allowed, "rights_basis"] if _text(value)
    ]
    return (
        True,
        "source_level_training_rights_verified",
        ";".join(dict.fromkeys(evidence)),
        ";".join(dict.fromkeys(bases)),
    )


def build_drug_alias_map(model_table: pd.DataFrame) -> dict[str, DrugIdentity]:
    """Build an exact, order-independent Phase 1 drug-name lookup."""

    name_columns = [
        "generic_name",
        "display_name",
        "generic_name_fda_map",
        "_join_drug_name",
        "drug_name",
    ]
    aliases: dict[str, set[str]] = {}
    names_by_drug: dict[str, set[str]] = {}
    inchikeys_by_drug: dict[str, set[str]] = {}
    mol_weights_by_drug: dict[str, set[float]] = {}
    has_rdkit_mol_wt = "rdkit_mol_wt" in model_table.columns
    for _, row in model_table.iterrows():
        drug_id = _text(row.get("drug_id"))
        drug_name = _first_text(row, name_columns)
        if not drug_id or not drug_name:
            continue
        names_by_drug.setdefault(drug_id, set()).add(drug_name)
        inchikey = _text(row.get("inchikey")).upper()
        if inchikey:
            inchikeys_by_drug.setdefault(drug_id, set()).add(inchikey)
        if has_rdkit_mol_wt:
            molecular_weight = _positive_finite_number(row.get("rdkit_mol_wt"))
            if molecular_weight is not None:
                mol_weights_by_drug.setdefault(drug_id, set()).add(molecular_weight)
        for column in name_columns:
            alias = normalize_key(row.get(column))
            if alias:
                aliases.setdefault(alias, set()).add(drug_id)

    resolved: dict[str, DrugIdentity] = {}
    for alias in sorted(aliases):
        drug_ids = aliases[alias]
        if len(drug_ids) != 1:
            continue
        drug_id = next(iter(drug_ids))
        inchikeys = inchikeys_by_drug.get(drug_id, set())
        if len(inchikeys) > 1:
            continue
        names = names_by_drug[drug_id]
        drug_name = min(names, key=lambda name: (normalize_key(name), name))
        mol_weights = mol_weights_by_drug.get(drug_id, set())
        resolved[alias] = DrugIdentity(
            drug_id=drug_id,
            drug_name=drug_name,
            inchikey=next(iter(inchikeys)) if inchikeys else "",
            rdkit_mol_wt=next(iter(mol_weights)) if len(mol_weights) == 1 else None,
        )
    return resolved


def build_canonical_pk_identity_table(
    model_table: pd.DataFrame,
    fda_mapping: pd.DataFrame | None,
) -> pd.DataFrame:
    """Replace evidence-derived identities with validated ligand-base identities."""

    if fda_mapping is None or fda_mapping.empty or "ligand_base" not in model_table:
        return model_table.copy()

    required = {
        "rdk_id",
        "preferred_identity",
        "identity_parent_inchikey",
        "identity_structure_validated",
    }
    missing = sorted(required.difference(fda_mapping.columns))
    if missing:
        raise ValueError(
            "FDA identity mapping is missing required columns: " + ", ".join(missing)
        )

    model = model_table.copy()
    mapping = fda_mapping.copy()
    model["_pk_rdk_id"] = model["ligand_base"].map(
        lambda value: _text(value).casefold()
    )
    mapping["_pk_rdk_id"] = mapping["rdk_id"].map(lambda value: _text(value).casefold())
    mapping["_pk_name"] = mapping.apply(
        lambda row: _first_text(
            row,
            (
                "preferred_identity",
                "resolved_preferred_name",
                "generic_name",
                "display_name",
            ),
        ),
        axis=1,
    )
    mapping["_pk_parent_inchikey"] = (
        mapping["identity_parent_inchikey"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
    )
    model_rdk_ids = set(model["_pk_rdk_id"]) - {""}
    mapping = mapping.loc[
        mapping["_pk_rdk_id"].isin(model_rdk_ids)
        & mapping["identity_structure_validated"].map(_truthy)
        & mapping["_pk_rdk_id"].ne("")
        & mapping["_pk_name"].ne("")
        & mapping["_pk_parent_inchikey"].ne("")
    ].copy()

    conflicting_rdk_ids = set(
        mapping.groupby("_pk_rdk_id", dropna=False).filter(
            lambda group: (
                group["_pk_parent_inchikey"].nunique() > 1
                or group["_pk_name"].map(normalize_key).nunique() > 1
            )
        )["_pk_rdk_id"]
    )
    mapping = mapping.loc[~mapping["_pk_rdk_id"].isin(conflicting_rdk_ids)]
    mapping = mapping.drop_duplicates("_pk_rdk_id", keep="first")

    def unique_molecular_weight(values: pd.Series) -> float | None:
        numeric = pd.to_numeric(values, errors="coerce").dropna()
        unique = sorted({round(float(value), 8) for value in numeric if value > 0})
        return unique[0] if len(unique) == 1 else None

    molecular_weights = (
        model.groupby("_pk_rdk_id", dropna=False)["rdkit_mol_wt"]
        .agg(unique_molecular_weight)
        .to_dict()
        if "rdkit_mol_wt" in model
        else {}
    )
    mapping["rdkit_mol_wt"] = mapping["_pk_rdk_id"].map(molecular_weights)

    canonical = pd.DataFrame(
        {
            "drug_id": mapping["_pk_parent_inchikey"],
            "generic_name": mapping["_pk_name"],
            "display_name": mapping["_pk_name"],
            "drug_name": mapping["_pk_name"],
            "inchikey": mapping["_pk_parent_inchikey"],
            "rdkit_mol_wt": mapping["rdkit_mol_wt"],
            "ligand_base": mapping["rdk_id"],
            "pk_identity_source": "canonical_fda_mapping_validated_parent",
        }
    )
    admitted = set(mapping["_pk_rdk_id"])
    fallback = model.loc[~model["_pk_rdk_id"].isin(admitted)].drop(
        columns=["_pk_rdk_id"]
    )
    return pd.concat([canonical, fallback], ignore_index=True, sort=False)


def _request_bytes(url: str, *, timeout: int, attempts: int = 3) -> bytes:
    if urlparse(url).netloc != "pk-db.com":
        raise ValueError(f"refusing non-PK-DB URL: {url}")
    error: Exception | None = None
    for attempt in range(attempts):
        try:
            response = requests.get(
                url,
                timeout=timeout,
                headers={"User-Agent": "Atlas-PK-context/1.0"},
            )
            response.raise_for_status()
            return response.content
        except requests.RequestException as exc:
            error = exc
            if attempt + 1 < attempts:
                time.sleep(2**attempt)
    raise RuntimeError(f"PK-DB request failed for {url}: {error}")


def _cached_bytes(
    url: str,
    *,
    cache_path: Path,
    timeout: int,
    refresh: bool,
) -> bytes:
    if cache_path.is_file() and not refresh:
        return cache_path.read_bytes()
    content = _request_bytes(url, timeout=timeout)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = cache_path.with_suffix(cache_path.suffix + ".partial")
    temporary.write_bytes(content)
    temporary.replace(cache_path)
    return content


def _json_bytes(content: bytes, *, source: str) -> dict[str, Any]:
    try:
        payload = json.loads(content)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid JSON from {source}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object from {source}")
    return payload


def fetch_study_index(
    *,
    cache_dir: Path,
    timeout: int = 120,
    refresh: bool = False,
    page_size: int = 250,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    page = 1
    while True:
        url = f"{PKDB_API_BASE}/studies/?format=json&page_size={page_size}&page={page}"
        path = cache_dir / "study_index" / f"page_{page:04d}.json"
        payload = _json_bytes(
            _cached_bytes(url, cache_path=path, timeout=timeout, refresh=refresh),
            source=url,
        )
        data = payload.get("data")
        page_rows = data.get("data") if isinstance(data, dict) else []
        rows.extend(item for item in page_rows or [] if isinstance(item, dict))
        last_page = int(payload.get("last_page") or page)
        if page >= last_page:
            break
        page += 1
    return rows


def _study_drug_matches(
    study: Mapping[str, Any],
    aliases: Mapping[str, DrugIdentity],
) -> dict[str, DrugIdentity]:
    matches: dict[str, DrugIdentity] = {}
    for substance in study.get("substances") or []:
        if not isinstance(substance, Mapping):
            continue
        for field in ("name", "label", "sid"):
            key = normalize_key(substance.get(field))
            if key in aliases:
                matches[key] = aliases[key]
    return matches


def _safe_sid(value: Any) -> str:
    sid = _text(value)
    if not sid or not re.fullmatch(r"[A-Za-z0-9]+", sid):
        raise ValueError(f"invalid PK-DB study SID: {sid!r}")
    return sid


def fetch_study_detail(
    sid: str,
    *,
    cache_dir: Path,
    timeout: int,
    refresh: bool,
) -> dict[str, Any]:
    safe_sid = _safe_sid(sid)
    url = f"{PKDB_API_BASE}/studies/{safe_sid}/?format=json"
    path = cache_dir / "study_detail" / f"{safe_sid}.json"
    return _json_bytes(
        _cached_bytes(url, cache_path=path, timeout=timeout, refresh=refresh),
        source=url,
    )


def _iter_values(value: Any) -> Iterable[Any]:
    if isinstance(value, Mapping):
        for nested in value.values():
            yield from _iter_values(nested)
    elif isinstance(value, list):
        for nested in value:
            yield from _iter_values(nested)
    else:
        yield value


def study_tsv_urls(detail: Mapping[str, Any]) -> list[str]:
    urls: set[str] = set()
    for value in _iter_values(detail):
        candidate = _text(value)
        if candidate.startswith("/media/data/"):
            candidate = f"https://pk-db.com{candidate}"
        if candidate.startswith(
            "https://pk-db.com/media/data/"
        ) and candidate.casefold().endswith(".tsv"):
            urls.add(candidate)
    return sorted(urls)


def _tsv_cache_path(cache_dir: Path, url: str) -> Path:
    basename = Path(urlparse(url).path).name.lstrip(".") or "source.tsv"
    digest = hashlib.sha256(url.encode("utf-8")).hexdigest()[:16]
    return cache_dir / "source_tsv" / f"{digest}_{basename}"


def fetch_tsv(
    url: str,
    *,
    cache_dir: Path,
    timeout: int,
    refresh: bool,
) -> tuple[pd.DataFrame, Path, str]:
    path = _tsv_cache_path(cache_dir, url)
    content = _cached_bytes(url, cache_path=path, timeout=timeout, refresh=refresh)
    digest = hashlib.sha256(content).hexdigest()
    frame = pd.read_csv(io.BytesIO(content), sep="\t", dtype="object")
    frame.columns = [str(column).strip().casefold() for column in frame.columns]
    return frame, path, digest


def _source_kind(frame: pd.DataFrame) -> str:
    columns = set(frame.columns)
    if {"name", "route", "application", "substance", "value", "unit"}.issubset(columns):
        return "interventions"
    if (
        {"measurement_type", "substance"}.issubset(columns)
        and ("intervention" in columns or "interventions" in columns)
        and columns.intersection({"value", "mean", "median"})
    ):
        return "outputs"
    if {"name", "measurement_type"}.issubset(columns) and columns.intersection(
        {"choice", "mean", "median", "min", "max"}
    ):
        return "groups"
    return "other"


def _measurement_context(value: Any) -> str:
    text = normalize_key(value)
    if re.search(r"\bc\s*max\b|\bcmax\b|maximum concentration", text):
        return "cmax"
    if "absolute bioavailability" in text or text in {"fabs", "absolute f"}:
        return "absolute_bioavailability"
    if "bioavailability" in text:
        return "bioavailability_unspecified"
    if "clearance" not in text and not re.fullmatch(r"cl(?:\s+f)?", text):
        return "other"
    if "creatinine" in text:
        return "subject_creatinine_clearance"
    if "renal" in text or text in {"clr", "clearance renal"}:
        return "renal_clearance"
    if "oral" in text or "apparent" in text or text in {"cl f", "clf"}:
        return "apparent_oral_clearance"
    if "total body" in text or text == "total clearance":
        return "total_body_clearance"
    if "systemic" in text:
        return "systemic_clearance"
    if "plasma" in text:
        return "plasma_clearance"
    if "blood" in text:
        return "blood_clearance"
    return "clearance_unspecified"


def _reported_statistic(row: Mapping[str, Any]) -> tuple[str, float | None]:
    for column in ("value", "mean", "median"):
        number = _numeric(row.get(column))
        if number is not None:
            return column, number
    return "", None


def _normalize_clearance(value: float, unit: str) -> tuple[float | None, str]:
    compact = (
        unit.casefold()
        .replace("litres", "l")
        .replace("liters", "l")
        .replace("liter", "l")
        .replace("millilitres", "ml")
        .replace("milliliters", "ml")
        .replace("milliliter", "ml")
        .replace("hours", "h")
        .replace("hour", "h")
        .replace("hrs", "h")
        .replace("hr", "h")
        .replace("minutes", "min")
        .replace("minute", "min")
        .replace("²", "2")
        .replace("^", "")
        .replace(" ", "")
    )
    factors = {
        "l/h": (1.0, "L/h"),
        "ml/h": (0.001, "L/h"),
        "l/min": (60.0, "L/h"),
        "ml/min": (0.06, "L/h"),
        "l/h/kg": (1.0, "L/h/kg"),
        "ml/h/kg": (0.001, "L/h/kg"),
        "l/min/kg": (60.0, "L/h/kg"),
        "ml/min/kg": (0.06, "L/h/kg"),
        "l/h/m2": (1.0, "L/h/m2"),
        "ml/h/m2": (0.001, "L/h/m2"),
        "l/min/m2": (60.0, "L/h/m2"),
        "ml/min/m2": (0.06, "L/h/m2"),
    }
    normalized = factors.get(compact)
    if normalized is None:
        return None, ""
    factor, canonical_unit = normalized
    return value * factor, canonical_unit


def _compact_concentration_unit(unit: str) -> str:
    return (
        unit.casefold()
        .replace("µ", "u")
        .replace("μ", "u")
        .replace("micromoles", "umol")
        .replace("micromole", "umol")
        .replace("nanomoles", "nmol")
        .replace("nanomole", "nmol")
        .replace("picomoles", "pmol")
        .replace("picomole", "pmol")
        .replace("millimoles", "mmol")
        .replace("millimole", "mmol")
        .replace("litres", "l")
        .replace("liters", "l")
        .replace("liter", "l")
        .replace("mcg", "ug")
        .replace(" ", "")
    )


def _normalize_cmax(
    value: float,
    unit: str,
    *,
    molecular_weight: float | None = None,
) -> tuple[float | None, str, str]:
    compact = _compact_concentration_unit(unit)
    molar_factors = {
        "pm": 1e-6,
        "pmol/l": 1e-6,
        "nm": 1e-3,
        "nmol/l": 1e-3,
        "um": 1.0,
        "umol/l": 1.0,
        "mm": 1e3,
        "mmol/l": 1e3,
        "m": 1e6,
        "mol/l": 1e6,
        "pmol/ml": 1e-3,
        "nmol/ml": 1.0,
        "umol/ml": 1e3,
    }
    molar_factor = molar_factors.get(compact)
    if molar_factor is not None:
        return value * molar_factor, "uM", "molar_concentration_unit_scale"

    mass_to_mg_per_l = {
        "pg/l": 1e-9,
        "ng/l": 1e-6,
        "ug/l": 1e-3,
        "mg/l": 1.0,
        "g/l": 1e3,
        "pg/ml": 1e-6,
        "ng/ml": 1e-3,
        "ug/ml": 1.0,
        "mg/ml": 1e3,
        "ng/dl": 1e-5,
        "ug/dl": 1e-2,
        "mg/dl": 10.0,
    }
    mass_factor = mass_to_mg_per_l.get(compact)
    if mass_factor is None:
        return None, "", ""
    if molecular_weight is None:
        return None, "", "mass_concentration_requires_rdkit_mol_wt"
    normalized = value * mass_factor * 1000.0 / molecular_weight
    return normalized, "uM", "mass_concentration_via_rdkit_mol_wt"


def _group_context(groups: pd.DataFrame) -> dict[str, dict[str, str]]:
    context: dict[str, dict[str, str]] = {}
    if groups.empty:
        return context
    for _, row in groups.iterrows():
        name = _text(row.get("name"))
        if not name:
            continue
        info = context.setdefault(name, {"parent": _text(row.get("parent"))})
        measurement = normalize_key(row.get("measurement_type"))
        if measurement:
            info[measurement] = _first_text(
                row, ("choice", "mean", "median", "min", "max")
            )

    for name, info in context.items():
        inherited = dict(info)
        parent = info.get("parent", "")
        visited = {name}
        while parent and parent not in visited and parent in context:
            visited.add(parent)
            parent_info = context[parent]
            for key, value in parent_info.items():
                inherited.setdefault(key, value)
            parent = parent_info.get("parent", "")
        context[name] = inherited
    return context


def _species(group_info: Mapping[str, str]) -> str:
    value = normalize_key(group_info.get("species"))
    if value in {"homo sapiens", "human", "humans"}:
        return "Homo sapiens"
    return _text(group_info.get("species"))


def _population(group_name: str, group_info: Mapping[str, str]) -> str:
    values = [group_name]
    for key in ("healthy", "disease", "age", "sex"):
        value = _text(group_info.get(key))
        if value:
            values.append(f"{key}={value}")
    return "; ".join(dict.fromkeys(values))


def _intervention_map(frame: pd.DataFrame) -> dict[str, dict[str, Any]]:
    if frame.empty:
        return {}
    return {
        _text(row.get("name")): row.to_dict()
        for _, row in frame.iterrows()
        if _text(row.get("name"))
    }


def _route_context(value: Any) -> str:
    route = normalize_key(value)
    if not route or route in {"other", "unknown", "not reported", "not specified"}:
        return ""
    intravenous = bool(re.search(r"\b(?:intravenous|intra venous|iv|i v)\b", route))
    extravascular_tokens = {
        "oral",
        "subcutaneous",
        "intramuscular",
        "inhaled",
        "inhalation",
        "intranasal",
        "nasal",
        "rectal",
        "transdermal",
        "topical",
        "buccal",
        "sublingual",
    }
    extravascular = any(token in route.split() for token in extravascular_tokens)
    if intravenous and extravascular:
        return ""
    return "intravenous" if intravenous else "extravascular"


def _linked_interventions(
    row: Mapping[str, Any],
    interventions: Mapping[str, Mapping[str, Any]],
) -> tuple[str, list[tuple[str, Mapping[str, Any]]], str]:
    raw = _first_text(row, ("intervention", "interventions"))
    names = list(
        dict.fromkeys(item.strip() for item in re.split(r"\|\||;", raw) if item.strip())
    )
    if not names:
        return raw, [], "intervention_not_linked"
    linked = [(name, interventions[name]) for name in names if name in interventions]
    if len(linked) != len(names):
        return raw, linked, "intervention_reference_not_found"
    status = (
        "exact_intervention_link" if len(linked) == 1 else "multiple_intervention_links"
    )
    return raw, linked, status


def _compatible_interventions(
    row: Mapping[str, Any],
    linked: Iterable[tuple[str, Mapping[str, Any]]],
) -> list[tuple[str, Mapping[str, Any]]]:
    analyte = normalize_key(row.get("substance"))
    return [
        (name, intervention)
        for name, intervention in linked
        if analyte and normalize_key(intervention.get("substance")) == analyte
    ]


def _endpoint_intervention_evidence(
    *,
    context: str,
    row: Mapping[str, Any],
    interventions: Mapping[str, Mapping[str, Any]],
) -> tuple[
    str,
    str,
    Mapping[str, Any],
    str,
    Mapping[str, Any],
    str,
    bool,
]:
    reported, linked, raw_link_status = _linked_interventions(row, interventions)
    compatible = _compatible_interventions(row, linked)
    intravenous = [
        item
        for item in compatible
        if _route_context(item[1].get("route")) == "intravenous"
    ]
    extravascular = [
        item
        for item in compatible
        if _route_context(item[1].get("route")) == "extravascular"
    ]

    primary_name = ""
    primary: Mapping[str, Any] = {}
    comparator_name = ""
    comparator: Mapping[str, Any] = {}
    if context == "absolute_bioavailability":
        if len(extravascular) == 1:
            primary_name, primary = extravascular[0]
        if len(intravenous) == 1:
            comparator_name, comparator = intravenous[0]
    elif context in IV_CLEARANCE_CONTEXTS:
        if len(intravenous) == 1:
            primary_name, primary = intravenous[0]
        elif len(compatible) == 1:
            primary_name, primary = compatible[0]
    elif context in APPARENT_CLEARANCE_CONTEXTS:
        if len(extravascular) == 1:
            primary_name, primary = extravascular[0]
        elif len(compatible) == 1:
            primary_name, primary = compatible[0]
    elif len(compatible) == 1:
        primary_name, primary = compatible[0]

    if raw_link_status not in {
        "exact_intervention_link",
        "multiple_intervention_links",
    }:
        endpoint_status = raw_link_status
        endpoint_ready = False
    elif not compatible:
        endpoint_status = "endpoint_intervention_substance_mismatch"
        endpoint_ready = False
    elif context == "absolute_bioavailability":
        if not extravascular:
            endpoint_status = (
                "absolute_bioavailability_endpoint_linked_extravascular_"
                "intervention_not_verified"
            )
            endpoint_ready = False
        elif not intravenous:
            endpoint_status = (
                "absolute_bioavailability_endpoint_linked_iv_comparator_not_verified"
            )
            endpoint_ready = False
        elif len(extravascular) != 1 or len(intravenous) != 1 or len(compatible) != 2:
            endpoint_status = "absolute_bioavailability_endpoint_links_ambiguous"
            endpoint_ready = False
        else:
            endpoint_status = (
                "absolute_bioavailability_endpoint_linked_comparator_verified"
            )
            endpoint_ready = True
    elif context in IV_CLEARANCE_CONTEXTS:
        if not intravenous:
            endpoint_status = "iv_clearance_context_not_verified"
            endpoint_ready = False
        elif len(intravenous) != 1 or len(compatible) != 1:
            endpoint_status = "iv_clearance_context_ambiguous"
            endpoint_ready = False
        else:
            endpoint_status = "iv_clearance_context_verified"
            endpoint_ready = True
    elif context in APPARENT_CLEARANCE_CONTEXTS:
        if not extravascular:
            endpoint_status = "apparent_clearance_extravascular_route_not_verified"
            endpoint_ready = False
        elif len(extravascular) != 1 or len(compatible) != 1:
            endpoint_status = "apparent_clearance_extravascular_route_ambiguous"
            endpoint_ready = False
        else:
            endpoint_status = "apparent_clearance_extravascular_route_verified"
            endpoint_ready = True
    elif len(compatible) != 1:
        endpoint_status = "endpoint_intervention_link_ambiguous"
        endpoint_ready = False
    else:
        endpoint_status = "endpoint_intervention_link_verified"
        endpoint_ready = True

    return (
        reported,
        primary_name,
        primary,
        comparator_name,
        comparator,
        endpoint_status,
        endpoint_ready,
    )


def collect_endpoint_inventory_rows(
    *,
    study: Mapping[str, Any],
    source_frames: Iterable[tuple[str, str, pd.DataFrame]],
    aliases: Mapping[str, DrugIdentity],
) -> pd.DataFrame:
    """Collect non-value endpoint metadata for coverage planning only."""

    output_parts: list[pd.DataFrame] = []
    group_parts: list[pd.DataFrame] = []
    for _, _, frame in source_frames:
        kind = _source_kind(frame)
        if kind == "outputs":
            output_parts.append(frame)
        elif kind == "groups":
            group_parts.append(frame)
    groups = _group_context(
        pd.concat(group_parts, ignore_index=True) if group_parts else pd.DataFrame()
    )
    rows: list[dict[str, Any]] = []
    for frame in output_parts:
        for _, row in frame.iterrows():
            identity = aliases.get(normalize_key(row.get("substance")))
            if identity is None:
                continue
            statistic_type, value = _reported_statistic(row)
            measurement_type = _text(row.get("measurement_type"))
            group_name = _text(row.get("group"))
            rows.append(
                {
                    "study_sid": _text(study.get("sid")),
                    "drug_id": identity.drug_id,
                    "measurement_type": measurement_type,
                    "measurement_context": _measurement_context(measurement_type),
                    "species": _species(groups.get(group_name, {})),
                    "has_numeric_value": value is not None,
                    "has_unit": bool(_text(row.get("unit"))),
                    "statistic_type": statistic_type,
                }
            )
    return pd.DataFrame(
        rows,
        columns=[
            "study_sid",
            "drug_id",
            "measurement_type",
            "measurement_context",
            "species",
            "has_numeric_value",
            "has_unit",
            "statistic_type",
        ],
    )


def parse_study_sources(
    *,
    study: Mapping[str, Any],
    source_frames: Iterable[tuple[str, str, pd.DataFrame]],
    aliases: Mapping[str, DrugIdentity],
    source_rights: pd.DataFrame | None = None,
) -> pd.DataFrame:
    outputs: list[tuple[str, str, pd.DataFrame]] = []
    intervention_parts: list[pd.DataFrame] = []
    group_parts: list[pd.DataFrame] = []
    for url, digest, frame in source_frames:
        kind = _source_kind(frame)
        if kind == "outputs":
            outputs.append((url, digest, frame))
        elif kind == "interventions":
            intervention_parts.append(frame)
        elif kind == "groups":
            group_parts.append(frame)

    interventions = _intervention_map(
        pd.concat(intervention_parts, ignore_index=True)
        if intervention_parts
        else pd.DataFrame()
    )
    groups = _group_context(
        pd.concat(group_parts, ignore_index=True) if group_parts else pd.DataFrame()
    )
    rights = normalize_source_rights(source_rights)
    (
        reference_text,
        reference_type,
        reference_id,
        publication_year,
    ) = _study_reference_metadata(study)
    rows: list[dict[str, Any]] = []
    for source_url, source_digest, frame in outputs:
        for row_number, (_, row) in enumerate(frame.iterrows(), start=2):
            identity = aliases.get(normalize_key(row.get("substance")))
            if identity is None:
                continue
            measurement_type = _text(row.get("measurement_type"))
            context = _measurement_context(measurement_type)
            if context == "other" or context == "subject_creatinine_clearance":
                continue
            statistic_type, value = _reported_statistic(row)
            if value is None:
                continue
            raw_unit = _text(row.get("unit"))
            normalized_value: float | None = None
            normalized_unit = ""
            normalization_method = ""
            if "clearance" in context:
                normalized_value, normalized_unit = _normalize_clearance(
                    value, raw_unit
                )
                if normalized_value is not None:
                    normalization_method = "clearance_unit_scale"
            elif context == "cmax":
                (
                    normalized_value,
                    normalized_unit,
                    normalization_method,
                ) = _normalize_cmax(
                    value,
                    raw_unit,
                    molecular_weight=identity.rdkit_mol_wt,
                )
            elif context == "absolute_bioavailability" and raw_unit.casefold() in {
                "%",
                "percent",
            }:
                normalized_value, normalized_unit = value, "%"
                normalization_method = "reported_percent"

            (
                reported_interventions,
                intervention_name,
                intervention,
                comparator_intervention_name,
                comparator_intervention,
                endpoint_link_status,
                endpoint_evidence_ready,
            ) = _endpoint_intervention_evidence(
                context=context,
                row=row,
                interventions=interventions,
            )
            (
                source_rights_allow_ml_training,
                source_rights_status,
                source_rights_reference,
                source_rights_basis,
            ) = _source_rights_decision(
                study=study,
                source_url=source_url,
                source_rights=rights,
            )
            group_name = _text(row.get("group"))
            group_info = groups.get(group_name, {})
            species = _species(group_info)
            statuses = [endpoint_link_status]
            if not species:
                statuses.append("species_unknown")
            elif species != "Homo sapiens":
                statuses.append("non_human")
            if normalized_value is None:
                statuses.append("unit_not_safely_normalized")
            if normalization_method == "mass_concentration_requires_rdkit_mol_wt":
                statuses.append(normalization_method)
            if context in {"clearance_unspecified", "bioavailability_unspecified"}:
                statuses.append("endpoint_semantics_unspecified")
            if normalize_key(study.get("licence")) == "open":
                statuses.append(
                    "pkdb_study_licence_open_is_not_source_training_permission"
                )
            statuses.append(source_rights_status)
            semantic_ready = (
                species == "Homo sapiens"
                and normalized_value is not None
                and context
                not in {"clearance_unspecified", "bioavailability_unspecified"}
                and endpoint_evidence_ready
            )
            license_allows_ml_training = False
            training_allowed = semantic_ready and source_rights_allow_ml_training
            model_ready = training_allowed
            source_row_id = f"{study.get('sid')}:{source_digest}:{row_number}"
            rows.append(
                {
                    "pkdb_context_id": hashlib.sha256(
                        source_row_id.encode("utf-8")
                    ).hexdigest()[:20],
                    "drug_id": identity.drug_id,
                    "drug_name": identity.drug_name,
                    "inchikey": identity.inchikey,
                    "rdkit_mol_wt": identity.rdkit_mol_wt,
                    "study_sid": _text(study.get("sid")),
                    "study_name": _text(study.get("name")),
                    "study_reference": reference_text,
                    "study_reference_type": reference_type,
                    "study_reference_id": reference_id,
                    "publication_year": publication_year,
                    "study_date": _text(study.get("date")),
                    "study_licence": _text(study.get("licence")),
                    "source_file_url": source_url,
                    "source_file_sha256": source_digest,
                    "source_row_number": row_number,
                    "measurement_type": measurement_type,
                    "measurement_context": context,
                    "statistic_type": statistic_type,
                    "raw_value": value,
                    "raw_unit": raw_unit,
                    "normalized_value": normalized_value,
                    "normalized_unit": normalized_unit,
                    "normalization_method": normalization_method,
                    "matrix": _text(row.get("tissue")),
                    "method": _text(row.get("method")),
                    "group_name": group_name,
                    "population": _population(group_name, group_info),
                    "species": species,
                    "reported_interventions": reported_interventions,
                    "intervention_name": intervention_name,
                    "comparator_intervention_name": comparator_intervention_name,
                    "dose_value": _numeric(intervention.get("value")),
                    "dose_unit": _text(intervention.get("unit")),
                    "route": _text(intervention.get("route")),
                    "reference_route": (
                        "intravenous"
                        if _route_context(comparator_intervention.get("route"))
                        == "intravenous"
                        else ""
                    ),
                    "comparator_route": _text(comparator_intervention.get("route")),
                    "formulation": _text(intervention.get("form")),
                    "regimen": _text(intervention.get("application")),
                    "endpoint_link_status": endpoint_link_status,
                    "endpoint_evidence_ready": endpoint_evidence_ready,
                    "context_status": ";".join(dict.fromkeys(statuses)),
                    "semantic_ready": semantic_ready,
                    "license_allows_ml_training": license_allows_ml_training,
                    "source_rights_status": source_rights_status,
                    "source_rights_reference": source_rights_reference,
                    "source_rights_basis": source_rights_basis,
                    "source_rights_allow_ml_training": (
                        source_rights_allow_ml_training
                    ),
                    "training_allowed": training_allowed,
                    "model_ready": model_ready,
                }
            )
    return pd.DataFrame(rows, columns=DETAIL_COLUMNS)


def to_pk_context(recovered: pd.DataFrame) -> pd.DataFrame:
    if recovered.empty:
        return pd.DataFrame(columns=["pk_context_id", *PK_CONTEXT_COLUMNS])
    rows: list[dict[str, Any]] = []
    for _, source in recovered.iterrows():
        row = {column: pd.NA for column in PK_CONTEXT_COLUMNS}
        context = _text(source.get("measurement_context"))
        row.update(
            {
                "drug_id": source.get("drug_id"),
                "drug_name": source.get("drug_name"),
                "inchikey": source.get("inchikey"),
                "source_name": "PK-DB_recovered_source_TSV",
                "source_version": (
                    str(int(source.get("publication_year")))
                    if pd.notna(source.get("publication_year"))
                    else source.get("study_date")
                ),
                "source_record_id": source.get("pkdb_context_id"),
                "source_url": source.get("source_file_url"),
                "study_id": source.get("study_sid"),
                "reference": source.get("study_reference"),
                "population": source.get("population"),
                "species": source.get("species"),
                "dose_value": source.get("dose_value"),
                "dose_unit": source.get("dose_unit"),
                "route": source.get("route"),
                "regimen": source.get("regimen"),
                "formulation": source.get("formulation"),
                "dose_context_type": (
                    "cmax_study_matched"
                    if context == "cmax"
                    else "pk_endpoint_study_matched"
                ),
                "parent_or_metabolite": "parent",
                "measurement_context": context,
                "reference_route": source.get("reference_route"),
                "training_allowed": source.get("training_allowed"),
                "extraction_method": SCHEMA_VERSION,
                "source_confidence": (
                    "high" if _truthy(source.get("semantic_ready")) else "medium"
                ),
                "context_status": source.get("context_status"),
                "license_note": (
                    f"PK-DB study licence={_text(source.get('study_licence'))}; "
                    "access metadata does not grant third-party ML-training rights; "
                    "source rights status="
                    f"{_text(source.get('source_rights_status'))}; evidence="
                    f"{_text(source.get('source_rights_reference'))}"
                ),
                "missing_reason": (
                    ""
                    if _truthy(source.get("model_ready"))
                    else source.get("context_status")
                ),
            }
        )
        value = source.get("normalized_value")
        unit = source.get("normalized_unit")
        if context == "cmax":
            row["cmax_value_raw"] = source.get("raw_value")
            row["cmax_unit_raw"] = source.get("raw_unit")
            if unit == "uM":
                row["cmax_um"] = value
        elif "clearance" in context:
            row["clearance_value"] = value
            row["clearance_unit"] = unit
        elif context == "absolute_bioavailability":
            row["bioavailability_value"] = value
            row["bioavailability_unit"] = unit
        rows.append(row)
    return finalize_context(pd.DataFrame(rows, columns=PK_CONTEXT_COLUMNS))


def recover_pkdb_context(
    *,
    model_table: pd.DataFrame,
    out_dir: str | Path,
    timeout: int = 120,
    workers: int = 4,
    refresh: bool = False,
    max_studies: int = 0,
    source_rights: pd.DataFrame | None = None,
) -> dict[str, Any]:
    output = Path(out_dir)
    cache_dir = output / "cache"
    output.mkdir(parents=True, exist_ok=True)
    aliases = build_drug_alias_map(model_table)
    rights = normalize_source_rights(source_rights)
    study_index = fetch_study_index(
        cache_dir=cache_dir,
        timeout=timeout,
        refresh=refresh,
    )
    matched = [study for study in study_index if _study_drug_matches(study, aliases)]
    open_access = [
        study for study in matched if normalize_key(study.get("licence")) == "open"
    ]
    selected = open_access
    if max_studies > 0:
        selected = selected[:max_studies]

    details: dict[str, dict[str, Any]] = {}
    failures: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=max(1, min(workers, 8))) as pool:
        detail_futures = {
            pool.submit(
                fetch_study_detail,
                _safe_sid(study.get("sid")),
                cache_dir=cache_dir,
                timeout=timeout,
                refresh=refresh,
            ): _safe_sid(study.get("sid"))
            for study in selected
        }
        for detail_future in as_completed(detail_futures):
            sid = detail_futures[detail_future]
            try:
                details[sid] = detail_future.result()
            except (OSError, RuntimeError, ValueError) as exc:
                failures.append(
                    {"study_sid": sid, "stage": "detail", "reason": str(exc)}
                )

    source_jobs: dict[str, str] = {}
    for sid, detail in details.items():
        for url in study_tsv_urls(detail):
            source_jobs[url] = sid
    source_results: dict[str, tuple[str, pd.DataFrame]] = {}
    with ThreadPoolExecutor(max_workers=max(1, min(workers, 8))) as pool:
        tsv_futures = {
            pool.submit(
                fetch_tsv,
                url,
                cache_dir=cache_dir,
                timeout=timeout,
                refresh=refresh,
            ): url
            for url in source_jobs
        }
        for tsv_future in as_completed(tsv_futures):
            url = tsv_futures[tsv_future]
            try:
                frame, _, digest = tsv_future.result()
                source_results[url] = (digest, frame)
            except (OSError, RuntimeError, ValueError, pd.errors.ParserError) as exc:
                failures.append(
                    {
                        "study_sid": source_jobs[url],
                        "stage": "source_tsv",
                        "reason": str(exc),
                    }
                )

    parts: list[pd.DataFrame] = []
    inventory_parts: list[pd.DataFrame] = []
    for sid, detail in details.items():
        frames = [
            (url, source_results[url][0], source_results[url][1])
            for url in study_tsv_urls(detail)
            if url in source_results
        ]
        parts.append(
            parse_study_sources(
                study=detail,
                source_frames=frames,
                aliases=aliases,
                source_rights=rights,
            )
        )
        inventory_parts.append(
            collect_endpoint_inventory_rows(
                study=detail,
                source_frames=frames,
                aliases=aliases,
            )
        )
    nonempty_parts = [
        part.dropna(axis="columns", how="all") for part in parts if not part.empty
    ]
    recovered = (
        pd.concat(nonempty_parts, ignore_index=True).reindex(columns=DETAIL_COLUMNS)
        if nonempty_parts
        else pd.DataFrame(columns=DETAIL_COLUMNS)
    )
    recovered = recovered.drop_duplicates("pkdb_context_id", keep="first")
    context = to_pk_context(recovered)
    nonempty_inventory = [part for part in inventory_parts if not part.empty]
    inventory_rows = (
        pd.concat(nonempty_inventory, ignore_index=True)
        if nonempty_inventory
        else pd.DataFrame(
            columns=[
                "study_sid",
                "drug_id",
                "measurement_type",
                "measurement_context",
                "species",
                "has_numeric_value",
                "has_unit",
                "statistic_type",
            ]
        )
    )
    endpoint_inventory = (
        inventory_rows.groupby(
            ["measurement_type", "measurement_context"],
            dropna=False,
        )
        .agg(
            rows=("drug_id", "size"),
            drugs=("drug_id", "nunique"),
            studies=("study_sid", "nunique"),
            human_rows=(
                "species",
                lambda values: int(values.eq("Homo sapiens").sum()),
            ),
            numeric_rows=(
                "has_numeric_value",
                lambda values: int(values.fillna(False).astype(bool).sum()),
            ),
            unit_rows=(
                "has_unit",
                lambda values: int(values.fillna(False).astype(bool).sum()),
            ),
        )
        .reset_index()
        .sort_values(["rows", "measurement_type"], ascending=[False, True])
        if not inventory_rows.empty
        else pd.DataFrame(
            columns=[
                "measurement_type",
                "measurement_context",
                "rows",
                "drugs",
                "studies",
                "human_rows",
                "numeric_rows",
                "unit_rows",
            ]
        )
    )
    endpoint_inventory.to_csv(output / "pkdb_endpoint_inventory.csv", index=False)
    recovered.to_csv(output / "pkdb_recovered_context.csv", index=False)
    context.to_csv(output / "pkdb_pk_context.csv", index=False)
    pd.DataFrame(failures, columns=["study_sid", "stage", "reason"]).to_csv(
        output / "pkdb_recovery_failures.csv", index=False
    )

    def count_flag(column: str) -> int:
        values = recovered.get(column, pd.Series(False, index=recovered.index))
        return int(values.fillna(False).astype(bool).sum())

    summary_columns = [
        "measurement_context",
        "rows",
        "drugs",
        "studies",
        "human_rows",
        "semantic_ready_rows",
        "endpoint_evidence_ready_rows",
        "source_rights_verified_rows",
        "training_allowed_rows",
        "model_ready_rows",
    ]
    summary = (
        recovered.groupby("measurement_context", dropna=False)
        .agg(
            rows=("pkdb_context_id", "size"),
            drugs=("drug_id", "nunique"),
            studies=("study_sid", "nunique"),
            human_rows=(
                "species",
                lambda values: int(values.eq("Homo sapiens").sum()),
            ),
            semantic_ready_rows=(
                "semantic_ready",
                lambda values: int(values.fillna(False).astype(bool).sum()),
            ),
            endpoint_evidence_ready_rows=(
                "endpoint_evidence_ready",
                lambda values: int(values.fillna(False).astype(bool).sum()),
            ),
            source_rights_verified_rows=(
                "source_rights_allow_ml_training",
                lambda values: int(values.fillna(False).astype(bool).sum()),
            ),
            training_allowed_rows=(
                "training_allowed",
                lambda values: int(values.fillna(False).astype(bool).sum()),
            ),
            model_ready_rows=(
                "model_ready",
                lambda values: int(values.fillna(False).astype(bool).sum()),
            ),
        )
        .reset_index()
        if not recovered.empty
        else pd.DataFrame(columns=summary_columns)
    )
    summary.to_csv(output / "pkdb_recovery_summary.csv", index=False)

    verified_study_ids = set(
        recovered.loc[
            recovered.get(
                "source_rights_allow_ml_training",
                pd.Series(False, index=recovered.index),
            )
            .fillna(False)
            .astype(bool),
            "study_sid",
        ]
        .fillna("")
        .astype(str)
    )
    selected_study_ids = {_text(study.get("sid")) for study in selected}
    semantic_ready_rows = count_flag("semantic_ready")
    endpoint_evidence_ready_rows = count_flag("endpoint_evidence_ready")
    source_rights_verified_rows = count_flag("source_rights_allow_ml_training")
    training_allowed_rows = count_flag("training_allowed")
    model_ready_rows = count_flag("model_ready")
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "api_base": PKDB_API_BASE,
        "recovery_path": (
            "documented study endpoint plus publicly accessible source TSVs"
        ),
        "reason": (
            "documented Elasticsearch /outputs/ and bulk outputs.csv were empty; "
            "source TSV rows are retained as contextual evidence but remain "
            "quarantined from ML unless original-source rights and endpoint-linked "
            "route/comparator evidence are explicit"
        ),
        "rights_policy": (
            "PK-DB study licence=open indicates public access only and is never "
            "treated as third-party ML-training permission; a source-rights row "
            "with an affirmative training_allowed value and rights_reference is required"
        ),
        "endpoint_policy": (
            "systemic/plasma/blood/total clearance requires a same-analyte linked "
            "intravenous intervention; apparent CL/F requires a same-analyte linked "
            "extravascular intervention; absolute bioavailability requires linked "
            "same-analyte extravascular and intravenous comparator interventions"
        ),
        "study_index_rows": len(study_index),
        "phase1_aliases": len(aliases),
        "matched_studies": len(matched),
        "open_access_studies": len(open_access),
        "public_access_excluded_studies": len(matched) - len(open_access),
        "selection_limited_studies": len(open_access) - len(selected),
        "selected_studies": len(selected),
        "source_rights_manifest_rows": len(rights),
        "source_rights_verified_studies": len(verified_study_ids),
        "rights_excluded_studies": len(selected_study_ids - verified_study_ids),
        "retrieved_study_details": len(details),
        "source_tsv_urls": len(source_jobs),
        "recovered_rows": len(recovered),
        "pk_context_rows": len(context),
        "endpoint_inventory_rows": int(len(inventory_rows)),
        "endpoint_inventory_types": int(len(endpoint_inventory)),
        "endpoint_inventory_human_rows": int(
            inventory_rows["species"].eq("Homo sapiens").sum()
        ),
        "semantic_ready_rows": semantic_ready_rows,
        "endpoint_evidence_ready_rows": endpoint_evidence_ready_rows,
        "source_rights_verified_rows": source_rights_verified_rows,
        "training_allowed_rows": training_allowed_rows,
        "model_ready_rows": model_ready_rows,
        "quarantined_recovered_rows": len(recovered) - model_ready_rows,
        "failure_count": len(failures),
        "summary": summary.to_dict("records"),
        "terms_url": PKDB_TERMS_URL,
    }
    (output / "pkdb_recovery_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest
