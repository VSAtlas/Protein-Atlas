from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from typing import Any

import pandas as pd

from analysis.external.source_tables import download_to_cache


SCHEMA_VERSION = "atlas_spd_cmax_source_recovery_v2"
SPD_WORKBOOK_SHA256 = (
    "6b7681c4cc51670740ca8a20a6317927e47fe3672ae768b5ac4ab2512f06c640"
)
SMIT_COMMIT = "b86bedfdfd13243e25e928f48877f8c026160a92"
SMIT_BASE_URL = (
    "https://raw.githubusercontent.com/inessmit/adverse_event_analysis/"
    f"{SMIT_COMMIT}/plasma_concentrations/results"
)
SMIT_FILES = {
    "total_plasma_concentrations_approved_drugs_with_refs.txt": (
        "a9d98d311be8dc346a4e3eefa1ebd21eab785606e5a941413e5e99a636096b47"
    ),
    "molregno2median_plasma_total_unbound.txt": (
        "2d33f675310fd3db7d915ec9bd93bc14a279e965643d97f952bd17a03190aa90"
    ),
}

REVIEW_CONTEXT_REQUIRED_COLUMNS = frozenset(
    {
        "drug_id",
        "source_name",
        "source_record_id",
        "original_source_record_id",
        "inchikey",
        "reference",
        "citation_status",
        "context_status",
        "dose_context_type",
        "study_id",
        "cmax_um",
        "training_allowed",
    }
)

# Reviewed synonym mappings. Every mapping is also required to reproduce the
# published SPD/Smit Cmax value within the configured numeric tolerance.
SPD_TO_SMIT_NAME = {
    "adipiodone": "iodipamide",
    "alimemazine": "methylpromazine",
    "aminophenazone": "aminopyrine",
    "betanidine": "bethanidine",
    "cefalotin": "cephalothin",
    "colecalciferol": "cholecalciferol",
    "dicoumarol": "dicumarol",
    "dimetindene": "dimethindene",
    "dosulepin": "dothiepin",
    "metamizole": "dipyrone",
    "phenazone": "antipyrine",
    "secbutabarbital": "butabarbital",
    "sodiumaurothiomalate": "goldsodiumthiomalate",
    "sulfafurazole": "sulfisoxazole",
    "sultiame": "sulthiame",
}

_DOSE_RE = re.compile(
    r"(?<![\w.])(?P<value>\d+(?:\.\d+)?)\s*"
    r"(?P<unit>ng|mcg|ug|mg|g)(?:\s*/\s*(?P<denom>kg|m2))?(?!\w)",
    re.IGNORECASE,
)
_ROUTES = {
    "intravenous": (r"\bintravenous(?:ly)?\b", r"\bi\.?v\.?\b"),
    "oral": (r"\boral(?:ly)?\b", r"\bper os\b"),
    "intramuscular": (r"\bintramuscular(?:ly)?\b", r"\bi\.?m\.?\b"),
    "subcutaneous": (r"\bsubcutaneous(?:ly)?\b", r"\bs\.?c\.?\b"),
    "inhaled": (r"\binhal(?:ed|ation)\b",),
    "transdermal": (r"\btransdermal\b",),
    "topical": (r"\btopical(?:ly)?\b",),
}
_FORMULATIONS = {
    "capsule": r"\bcapsules?\b",
    "injection": r"\binjections?\b",
    "infusion": r"\binfusions?\b",
    "solution": r"\bsolutions?\b",
    "suspension": r"\bsuspensions?\b",
    "tablet": r"\btablets?\b",
}
_REGIMENS = (
    ("once_daily", r"\b(?:once daily|once a day|q\.?d\.?)\b"),
    ("twice_daily", r"\b(?:twice daily|twice a day|b\.?i\.?d\.?)\b"),
    ("three_times_daily", r"\b(?:three times daily|t\.?i\.?d\.?)\b"),
    ("four_times_daily", r"\b(?:four times daily|q\.?i\.?d\.?)\b"),
    ("single_dose", r"\b(?:single[- ]dose|single dose)\b"),
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_spd_workbook(
    workbook: str | Path,
    *,
    expected_sha256: str = SPD_WORKBOOK_SHA256,
) -> str:
    path = Path(workbook)
    if not path.is_file():
        raise FileNotFoundError(f"SPD workbook not found: {path}")
    actual = _sha256(path)
    if actual != expected_sha256:
        raise ValueError(
            "checksum mismatch for SPD workbook "
            f"{path}: expected {expected_sha256}, got {actual}"
        )
    return actual


def validate_spd_recovery_cache(
    recovery_dir: str | Path,
    *,
    expected_workbook_sha256: str = SPD_WORKBOOK_SHA256,
) -> bool:
    root = Path(recovery_dir)
    manifest_path = root / "spd_cmax_source_recovery_manifest.json"
    context_path = root / "spd_cmax_administration_context.csv"
    if not manifest_path.is_file() or not context_path.is_file():
        return False
    try:
        manifest = json.loads(manifest_path.read_text())
        columns = set(pd.read_csv(context_path, nrows=0).columns)
    except (OSError, ValueError, pd.errors.ParserError):
        return False
    inputs = manifest.get("inputs", {})
    output_hashes = manifest.get("output_sha256", {})
    smit_files = inputs.get("smit_files", {})
    if (
        manifest.get("schema_version") != SCHEMA_VERSION
        or not REVIEW_CONTEXT_REQUIRED_COLUMNS.issubset(columns)
        or inputs.get("spd_workbook_sha256") != expected_workbook_sha256
        or inputs.get("spd_workbook_expected_sha256") != expected_workbook_sha256
        or inputs.get("smit_commit") != SMIT_COMMIT
        or output_hashes.get("administration_context") != _sha256(context_path)
    ):
        return False
    return all(
        smit_files.get(filename, {}).get("sha256") == expected_hash
        for filename, expected_hash in SMIT_FILES.items()
    )


def _normalize_name(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return re.sub(r"[^a-z0-9]+", "", str(value).casefold())


def _clean(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def _float(value: Any) -> float | None:
    parsed = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    return None if pd.isna(parsed) else float(parsed)


def _relative_error(left: float, right: float) -> float:
    scale = max(abs(left), abs(right), 1e-15)
    return abs(left - right) / scale


def ensure_smit_source(source_dir: str | Path, *, download: bool = True) -> dict[str, Path]:
    root = Path(source_dir)
    root.mkdir(parents=True, exist_ok=True)
    resolved: dict[str, Path] = {}
    for filename, expected_sha256 in SMIT_FILES.items():
        path = root / filename
        if not path.is_file() and download:
            download_to_cache(f"{SMIT_BASE_URL}/{filename}", path, retries=3)
        if not path.is_file():
            raise FileNotFoundError(f"missing pinned Smit source file: {path}")
        actual = _sha256(path)
        if actual != expected_sha256:
            raise ValueError(
                f"checksum mismatch for {path}: expected {expected_sha256}, got {actual}"
            )
        resolved[filename] = path
    return resolved


def _source_kind(source: str) -> str:
    return {
        "NIBR curation": "primary_literature_curation_without_row_citation",
        "Pharmapendium": "licensed_heterogeneous_q3_aggregate",
        "Smit et al.": "public_multi_source_median_aggregate",
    }.get(source, "unresolved")


def load_spd_cmax_inventory(workbook: str | Path) -> pd.DataFrame:
    raw = pd.read_excel(workbook, sheet_name="S Data 2", header=2)
    ids = pd.to_numeric(raw["drugcentral struct id"], errors="coerce").astype("Int64")
    source_by_id = dict(zip(ids.astype("string"), raw["Cmax source"], strict=False))
    rows: list[dict[str, Any]] = []
    for position, row in raw.iterrows():
        cmax = _float(row.get("Cmax tot (uM)"))
        if cmax is None:
            continue
        drug_id = str(int(row["drugcentral struct id"]))
        source = _clean(row.get("Cmax source"))
        parent_value = _float(row.get("drugcentral prescribed id (when prodrug)"))
        parent_id = str(int(parent_value)) if parent_value is not None else ""
        inherited = not source and _clean(row.get("Cmax from parent?")).casefold() == "yes"
        if inherited and parent_id:
            source = _clean(source_by_id.get(parent_id))
        observation_id = f"spd:sutherland2023:sdata2:drugcentral:{drug_id}:cmax_total_um"
        rows.append(
            {
                "observation_id": observation_id,
                "excel_row": int(position) + 4,
                "drugcentral_struct_id": drug_id,
                "drug_name": _clean(row.get("drugcentral name")),
                "inchikey": _clean(row.get("drugcentral inchikey")).upper(),
                "cmax_total_um": cmax,
                "free_cmax_um": _float(row.get("Cmax, free (uM)")),
                "ppb_percent": _float(row.get("PPB %")),
                "cmax_source_raw": source,
                "cmax_source_kind": _source_kind(source),
                "parent_drugcentral_struct_id": parent_id,
                "source_inherited_from_parent": inherited,
                "derived_from_observation_id": (
                    f"spd:sutherland2023:sdata2:drugcentral:{parent_id}:cmax_total_um"
                    if inherited and parent_id
                    else ""
                ),
                "citation_status": "provider_label_only" if source else "unresolved",
                "context_status": "not_reported_in_spd_workbook",
            }
        )
    return pd.DataFrame(rows)


def _parse_context(text: str) -> dict[str, Any]:
    normalized = text.replace("microgram", "ug").replace("micrograms", "ug")
    doses: list[tuple[float, str]] = []
    for match in _DOSE_RE.finditer(normalized):
        unit = match.group("unit").casefold().replace("mcg", "ug")
        if match.group("denom"):
            unit += "/" + match.group("denom").casefold()
        doses.append((float(match.group("value")), unit))
    distinct_doses = list(dict.fromkeys(doses))
    routes = [
        route
        for route, patterns in _ROUTES.items()
        if any(re.search(pattern, normalized, re.IGNORECASE) for pattern in patterns)
    ]
    formulations = [
        name
        for name, pattern in _FORMULATIONS.items()
        if re.search(pattern, normalized, re.IGNORECASE)
    ]
    regimens = [
        name
        for name, pattern in _REGIMENS
        if re.search(pattern, normalized, re.IGNORECASE)
    ]
    steady_state = bool(re.search(r"\bsteady[- ]state\b", normalized, re.IGNORECASE))
    return {
        "candidate_dose_value": distinct_doses[0][0] if len(distinct_doses) == 1 else None,
        "candidate_dose_unit": distinct_doses[0][1] if len(distinct_doses) == 1 else "",
        "candidate_route": routes[0] if len(routes) == 1 else "",
        "candidate_regimen": regimens[0] if len(regimens) == 1 else "",
        "candidate_formulation": formulations[0] if len(formulations) == 1 else "",
        "candidate_steady_state": True if steady_state else None,
        "context_parse_status": (
            "ambiguous_multiple_values"
            if len(distinct_doses) > 1 or len(routes) > 1 or len(regimens) > 1
            else "candidate_extracted"
            if distinct_doses or routes or regimens or formulations or steady_state
            else "no_context_in_source_record"
        ),
    }


def _median_positions(group: pd.DataFrame) -> set[int]:
    ordered = group.loc[group["Molar_value"].notna()].sort_values(
        ["Molar_value", "_source_row"], kind="stable"
    )
    size = len(ordered)
    if size == 0:
        return set()
    positions = [size // 2] if size % 2 else [size // 2 - 1, size // 2]
    return {int(ordered.iloc[position]["_source_row"]) for position in positions}


def map_smit_cmax_sources(
    inventory: pd.DataFrame,
    source_files: dict[str, Path],
    *,
    value_tolerance: float = 0.01,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    medians = pd.read_csv(
        source_files["molregno2median_plasma_total_unbound.txt"], sep="\t"
    )
    observations = pd.read_csv(
        source_files["total_plasma_concentrations_approved_drugs_with_refs.txt"],
        sep="\t",
        low_memory=False,
    )
    medians["_name_key"] = medians["pref_name"].map(_normalize_name)
    median_by_name = medians.drop_duplicates("_name_key").set_index("_name_key")
    observations["_source_row"] = range(len(observations))
    context = observations["description"].fillna("").astype(str).map(_parse_context)
    context_frame = pd.DataFrame(context.tolist(), index=observations.index)
    observations = pd.concat([observations, context_frame], axis=1)

    mapped_rows: list[pd.DataFrame] = []
    mapping_rows: list[dict[str, Any]] = []
    smit_inventory = inventory.loc[inventory["cmax_source_raw"].eq("Smit et al.")]
    for _, spd in smit_inventory.iterrows():
        raw_key = _normalize_name(spd["drug_name"])
        lookup_key = SPD_TO_SMIT_NAME.get(raw_key, raw_key)
        match_method = "curated_synonym" if lookup_key != raw_key else "normalized_name"
        if lookup_key not in median_by_name.index:
            mapping_rows.append(
                {
                    **spd.to_dict(),
                    "smit_match_status": "unmatched_name",
                    "smit_match_method": match_method,
                }
            )
            continue
        median = median_by_name.loc[lookup_key]
        smit_um = float(median["median Molar total plasma concentration"]) * 1e6
        relative_error = _relative_error(float(spd["cmax_total_um"]), smit_um)
        if relative_error > value_tolerance:
            mapping_rows.append(
                {
                    **spd.to_dict(),
                    "smit_match_status": "cmax_value_mismatch",
                    "smit_match_method": match_method,
                    "smit_chembl_id": _clean(median.get("chembl_id")),
                    "smit_molregno": int(median["molregno"]),
                    "smit_median_cmax_um": smit_um,
                    "smit_relative_error": relative_error,
                }
            )
            continue
        group = observations.loc[
            pd.to_numeric(observations["molregno"], errors="coerce").eq(
                float(median["molregno"])
            )
        ].copy()
        contributors = _median_positions(group)
        group["is_median_contributor"] = group["_source_row"].isin(contributors)
        group.insert(0, "spd_observation_id", spd["observation_id"])
        group.insert(1, "spd_drugcentral_struct_id", spd["drugcentral_struct_id"])
        group.insert(2, "spd_drug_name", spd["drug_name"])
        group.insert(3, "spd_cmax_total_um", spd["cmax_total_um"])
        mapped_rows.append(group)
        recomputed = float(pd.to_numeric(group["Molar_value"], errors="coerce").median())
        mapping_rows.append(
            {
                **spd.to_dict(),
                "smit_match_status": "verified_aggregate_match",
                "smit_match_method": match_method,
                "smit_chembl_id": _clean(median.get("chembl_id")),
                "smit_molregno": int(median["molregno"]),
                "smit_median_cmax_um": smit_um,
                "smit_relative_error": relative_error,
                "smit_source_observation_count": int(len(group)),
                "smit_median_contributor_count": int(len(contributors)),
                "smit_median_recomputed_relative_error": _relative_error(
                    recomputed, float(median["median Molar total plasma concentration"])
                ),
            }
        )
    mapped = pd.concat(mapped_rows, ignore_index=True) if mapped_rows else pd.DataFrame()
    return pd.DataFrame(mapping_rows), mapped


def _consensus(group: pd.DataFrame, column: str) -> Any:
    values = [value for value in group[column].tolist() if not pd.isna(value) and value != ""]
    unique = list(dict.fromkeys(values))
    return unique[0] if len(unique) == 1 else pd.NA


def build_administration_context(
    inventory: pd.DataFrame,
    mappings: pd.DataFrame,
    source_observations: pd.DataFrame,
) -> pd.DataFrame:
    mapping_by_id = mappings.set_index("observation_id", drop=False)
    rows: list[dict[str, Any]] = []
    for _, item in inventory.iterrows():
        source = item["cmax_source_raw"]
        row: dict[str, Any] = {
            "observation_id": item["observation_id"],
            "drug_id": item["drugcentral_struct_id"],
            "drug_name": item["drug_name"],
            "inchikey": item["inchikey"],
            "source_name": "SPD",
            "source_version": "Sutherland_2023_Supplementary_Data_2",
            "source_record_id": item["drugcentral_struct_id"],
            "original_source_record_id": "",
            "source_url": "https://doi.org/10.1038/s41467-023-40064-9",
            "cmax_source_raw": source,
            "cmax_um": item["cmax_total_um"],
            "free_cmax_um": item["free_cmax_um"],
            "protein_binding_percent": item["ppb_percent"],
            "measurement_context": "observed_cmax",
            "dose_context_type": "",
            "study_id": "",
            "reference": source,
            "dose_value": pd.NA,
            "dose_unit": "",
            "route": "",
            "regimen": "",
            "formulation": "",
            "steady_state": pd.NA,
            "candidate_dose_value": pd.NA,
            "candidate_dose_unit": "",
            "candidate_route": "",
            "candidate_regimen": "",
            "candidate_formulation": "",
            "candidate_steady_state": pd.NA,
            "citation_status": "provider_label_only",
            "context_status": "not_reported_in_spd_workbook",
            "training_allowed": False,
            "license_note": "SPD supplement CC-BY-4.0; upstream source restrictions may apply",
        }
        if source == "Pharmapendium":
            row["context_status"] = "aggregate_only_proprietary_q3"
        elif source == "NIBR curation":
            row["context_status"] = "source_citation_manifest_not_published"
        elif source == "Smit et al." and item["observation_id"] in mapping_by_id.index:
            match = mapping_by_id.loc[item["observation_id"]]
            row["citation_status"] = _clean(match.get("smit_match_status"))
            row["context_status"] = "aggregate_only_smit_median"
            row["license_note"] = (
                "Smit repository states MIT except SIDER; upstream ChEMBL/publication attribution applies"
            )
            contributors = source_observations.loc[
                source_observations.get("spd_observation_id", pd.Series(dtype="object"))
                .astype(str)
                .eq(str(item["observation_id"]))
                & source_observations.get(
                    "is_median_contributor", pd.Series(False, index=source_observations.index)
                ).fillna(False)
            ]
            if not contributors.empty:
                for field in (
                    "candidate_dose_value",
                    "candidate_dose_unit",
                    "candidate_route",
                    "candidate_regimen",
                    "candidate_formulation",
                    "candidate_steady_state",
                ):
                    row[field] = _consensus(contributors, field)
                pmids = [
                    str(int(value))
                    for value in pd.to_numeric(contributors["pubmed_id"], errors="coerce")
                    .dropna()
                    .unique()
                ]
                titles = [
                    value
                    for value in contributors["doc_title"].fillna("").astype(str).unique()
                    if value
                ]
                if len(pmids) == 1:
                    row["study_id"] = f"PMID:{pmids[0]}"
                if len(titles) == 1:
                    row["reference"] = titles[0]
                source_count = int(match.get("smit_source_observation_count") or 0)
                complete_candidate = any(
                    not pd.isna(row[field]) and row[field] != ""
                    for field in (
                        "candidate_dose_value",
                        "candidate_route",
                        "candidate_regimen",
                        "candidate_formulation",
                    )
                )
                if source_count == 1 and len(pmids) == 1 and complete_candidate:
                    # The Smit table stores a derived per-drug median. A sole
                    # contributor with parseable prose is useful for review,
                    # but it does not prove that every parsed field describes
                    # the exact Cmax scenario. Keep it out of training until a
                    # researcher explicitly adjudicates the source record.
                    row["context_status"] = "single_source_context_candidate"
                    row["citation_status"] = "original_source_record_candidate"
        rows.append(row)
    return pd.DataFrame(rows)


def recover_spd_cmax_context(
    *,
    workbook: str | Path,
    out_dir: str | Path,
    source_dir: str | Path,
    download: bool = True,
    value_tolerance: float = 0.01,
    expected_workbook_sha256: str | None = SPD_WORKBOOK_SHA256,
) -> dict[str, Any]:
    output = Path(out_dir)
    output.mkdir(parents=True, exist_ok=True)
    workbook_path = Path(workbook)
    workbook_sha256 = (
        verify_spd_workbook(
            workbook_path,
            expected_sha256=expected_workbook_sha256,
        )
        if expected_workbook_sha256 is not None
        else _sha256(workbook_path)
    )
    source_files = ensure_smit_source(source_dir, download=download)
    inventory = load_spd_cmax_inventory(workbook_path)
    mappings, observations = map_smit_cmax_sources(
        inventory, source_files, value_tolerance=value_tolerance
    )
    context = build_administration_context(inventory, mappings, observations)

    inventory_path = output / "spd_cmax_source_inventory.csv"
    mapping_path = output / "spd_smit_cmax_mapping.csv"
    observation_path = output / "spd_smit_cmax_source_observations.csv"
    context_path = output / "spd_cmax_administration_context.csv"
    inventory.to_csv(inventory_path, index=False)
    mappings.to_csv(mapping_path, index=False)
    observations.to_csv(observation_path, index=False)
    context.to_csv(context_path, index=False)

    status_counts = context["context_status"].value_counts(dropna=False).to_dict()
    manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "spd_workbook": str(workbook_path),
            "spd_workbook_sha256": workbook_sha256,
            "spd_workbook_expected_sha256": expected_workbook_sha256,
            "smit_commit": SMIT_COMMIT,
            "smit_files": {
                name: {"path": str(path), "sha256": _sha256(path)}
                for name, path in source_files.items()
            },
        },
        "counts": {
            "spd_cmax_rows": int(len(inventory)),
            "spd_cmax_source_counts": inventory["cmax_source_raw"]
            .value_counts(dropna=False)
            .to_dict(),
            "smit_rows": int(inventory["cmax_source_raw"].eq("Smit et al.").sum()),
            "smit_verified_aggregate_matches": int(
                mappings["smit_match_status"].eq("verified_aggregate_match").sum()
            ),
            "training_allowed_context_rows": int(context["training_allowed"].fillna(False).sum()),
            "context_status_counts": status_counts,
        },
        "policy": {
            "spd_source_field": "provider/aggregate label, not row-level citation",
            "nibr": "citation manifest absent; administration remains unavailable",
            "pharmapendium": "heterogeneous third-quartile aggregate; no synthetic study context",
            "smit": "median aggregate; contributor context remains candidate-only until explicit source adjudication",
            "missing_context": "null, never inferred from a label or another PK scenario",
        },
        "outputs": {
            "inventory": str(inventory_path),
            "smit_mapping": str(mapping_path),
            "smit_source_observations": str(observation_path),
            "administration_context": str(context_path),
        },
        "output_sha256": {
            "inventory": _sha256(inventory_path),
            "smit_mapping": _sha256(mapping_path),
            "smit_source_observations": _sha256(observation_path),
            "administration_context": _sha256(context_path),
        },
    }
    manifest_path = output / "spd_cmax_source_recovery_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    manifest["manifest"] = str(manifest_path)
    return manifest


def merge_verified_spd_administration_context(
    spd_context: pd.DataFrame,
    recovered: str | Path | pd.DataFrame | None,
) -> pd.DataFrame:
    if recovered is None:
        return spd_context
    source = (
        recovered.copy()
        if isinstance(recovered, pd.DataFrame)
        else pd.read_csv(recovered, low_memory=False)
    )
    required = REVIEW_CONTEXT_REQUIRED_COLUMNS
    if not required.issubset(source.columns):
        missing = sorted(required.difference(source.columns))
        raise ValueError(f"SPD administration context missing columns: {missing}")

    def identifier(values: pd.Series) -> pd.Series:
        return (
            values.astype("string")
            .str.strip()
            .str.replace(r"\.0+$", "", regex=True)
            .str.casefold()
        )

    requested = source["training_allowed"].astype("string").str.casefold().isin(
        {"1", "true", "yes", "y"}
    )
    administration_fields = [
        column
        for column in ("dose_value", "route", "regimen", "formulation")
        if column in source.columns
    ]
    has_administration = (
        source[administration_fields]
        .apply(
            lambda column: column.notna()
            & column.astype("string").str.strip().ne("")
        )
        .any(axis=1)
        if administration_fields
        else pd.Series(False, index=source.index)
    )
    provenance_verified = (
        source["source_name"].astype("string").str.strip().eq("SPD")
        & source["citation_status"].eq("verified_original_source_record")
        & source["context_status"].eq("verified_single_source_observation")
        & source["dose_context_type"].eq("cmax_study_matched")
        & source["study_id"].notna()
        & source["study_id"].astype("string").str.strip().ne("")
        & source["source_record_id"].notna()
        & source["source_record_id"].astype("string").str.strip().ne("")
        & source["original_source_record_id"].notna()
        & source["original_source_record_id"].astype("string").str.strip().ne("")
        & source["reference"].notna()
        & source["reference"].astype("string").str.strip().ne("")
        & has_administration
    ).fillna(False)

    context = spd_context.copy()
    context["_drug_key"] = identifier(context["drug_id"])
    if context["_drug_key"].duplicated().any():
        raise ValueError("SPD PK context contains duplicate drug_id rows")
    expected = context.set_index("_drug_key")
    source["_drug_key"] = identifier(source["drug_id"])
    known_drug = source["_drug_key"].isin(expected.index)

    expected_cmax = source["_drug_key"].map(
        pd.to_numeric(expected["cmax_um"], errors="coerce")
    )
    supplied_cmax = pd.to_numeric(source["cmax_um"], errors="coerce")
    same_cmax = (
        supplied_cmax.notna()
        & expected_cmax.notna()
        & (supplied_cmax - expected_cmax)
        .abs()
        .le(1e-9 + 1e-6 * expected_cmax.abs())
    )

    supplied_record = identifier(source["source_record_id"])
    expected_record = source["_drug_key"].map(expected["source_record_id"])
    same_source_record = supplied_record.eq(identifier(expected_record)).fillna(False)

    supplied_inchikey = source["inchikey"].astype("string").str.strip().str.upper()
    expected_inchikey = (
        source["_drug_key"]
        .map(expected["inchikey"])
        .astype("string")
        .str.strip()
        .str.upper()
    )
    same_inchikey = (
        supplied_inchikey.notna()
        & supplied_inchikey.ne("")
        & expected_inchikey.notna()
        & expected_inchikey.ne("")
        & supplied_inchikey.eq(expected_inchikey)
    ).fillna(False)

    approved = (
        requested
        & provenance_verified
        & known_drug
        & same_cmax
        & same_source_record
        & same_inchikey
    )
    invalid_requested = requested & ~approved
    if invalid_requested.any():
        bad_rows = source.loc[invalid_requested, "drug_id"].astype(str).tolist()
        raise ValueError(
            "training-approved SPD administration rows failed provenance or exact "
            f"Cmax/identity validation: {bad_rows[:10]}"
        )

    allowed = source.loc[approved].copy()
    if allowed.empty:
        return spd_context
    if allowed["_drug_key"].duplicated().any():
        raise ValueError("verified SPD administration context contains duplicate drug_id rows")
    updates = allowed.set_index("_drug_key")
    out = spd_context.copy()
    keys = identifier(out["drug_id"])
    original_records = keys.map(updates["original_source_record_id"])
    out["source_record_id"] = original_records.where(
        original_records.notna(), out["source_record_id"]
    )
    fields = (
        "study_id",
        "reference",
        "dose_value",
        "dose_unit",
        "route",
        "regimen",
        "formulation",
        "steady_state",
        "dose_context_type",
        "measurement_context",
        "context_status",
        "training_allowed",
        "license_note",
    )
    for field in fields:
        if field in updates:
            values = keys.map(updates[field])
            out[field] = values.where(values.notna(), out[field])
    return out


__all__ = [
    "REVIEW_CONTEXT_REQUIRED_COLUMNS",
    "SMIT_COMMIT",
    "SPD_WORKBOOK_SHA256",
    "build_administration_context",
    "ensure_smit_source",
    "load_spd_cmax_inventory",
    "map_smit_cmax_sources",
    "merge_verified_spd_administration_context",
    "recover_spd_cmax_context",
    "validate_spd_recovery_cache",
    "verify_spd_workbook",
]
