from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from collections.abc import Callable, Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


AUDITED_PHASE1_AFFECTED_COUNTS = {
    "9I52": 388,
    "6HUJ": 570,
    "8HCQ": 278,
}
AUDITED_PHASE1_9I52_LABEL_COUNTS = {
    "legacy": {
        "spd_binding_label": {
            "positive": 16,
            "negative": 372,
            "excluded": 0,
            "unknown": 0,
        },
        "spd_exposure_label": {
            "positive": 9,
            "negative": 166,
            "excluded": 0,
            "unknown": 213,
        },
    },
    "strict": {
        "spd_binding_label": {
            "positive": 15,
            "negative": 350,
            "excluded": 0,
            "unknown": 23,
        },
        "spd_exposure_label": {
            "positive": 7,
            "negative": 165,
            "excluded": 0,
            "unknown": 216,
        },
    },
}
AUDITED_PHASE1_9I52_BINARY_FLIPS = {
    "spd_binding_label": 28,
    "spd_exposure_label": 15,
}
AUDITED_PHASE1_9I52_PANEL_MATCHES = 385
AUDITED_PHASE1_9I52_UNMATCHED_DRUGS = frozenset(
    {"arformoterol", "isometheptene", "terbutaline"}
)
EXPECTED_AFFECTED_PDB_SCOPE = frozenset(AUDITED_PHASE1_AFFECTED_COUNTS)
LABEL_COLUMNS = (
    "spd_binding_label",
    "spd_exposure_label",
    "tissue_site_label",
    "mechanism_ml_label",
)
TARGET_COLUMNS = ("target_id", "target_gene", "target_uniprot")
MAPPING_AUDIT_COLUMNS = (
    "spd_receptor_mapping_policy_version",
    "spd_receptor_mapping_status",
    "spd_receptor_mapping_strict_eligible",
    "spd_receptor_mapping_exclusion_reason",
    "spd_receptor_mapping_contract_source",
    "spd_receptor_mapping_contract_sha256",
    "spd_receptor_mapping_changed",
    "spd_receptor_mapping_requires_target_derived_rebuild",
    "spd_receptor_mapping_legacy_target_id",
    "spd_receptor_mapping_legacy_target_gene",
    "spd_receptor_mapping_legacy_target_uniprot",
    "spd_receptor_mapping_revised_target_id",
    "spd_receptor_mapping_revised_target_gene",
    "spd_receptor_mapping_revised_target_uniprot",
    "spd_receptor_mapping_selected_chain",
    "spd_receptor_mapping_site_chain",
)
EXPERT_TABLES = {
    "binding": "ml_spd_binding_table.csv",
    "exposure": "ml_spd_exposure_table.csv",
    "tissue": "ml_spd_tissue_table.csv",
    "mechanism": "ml_spd_mechanism_table.csv",
}


class ContractViolation(RuntimeError):
    """Raised after comparison artifacts record a failed mapping contract."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _resolve_contract_path(path: str | Path | None) -> Path:
    if path is not None:
        resolved = Path(path).resolve()
    else:
        from analysis.ml.spd_receptor_mapping import DEFAULT_CONTRACT_PATH

        resolved = Path(DEFAULT_CONTRACT_PATH).resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"Receptor mapping contract does not exist: {resolved}")
    return resolved


def _required_file(path: str | Path, name: str) -> Path:
    resolved = Path(path).resolve()
    if not resolved.is_file():
        raise FileNotFoundError(f"{name} does not exist: {resolved}")
    return resolved


def _optional_file(path: str | Path | None, name: str) -> Path | None:
    if path is None:
        return None
    return _required_file(path, name)


def _prepare_out_dir(path: str | Path) -> Path:
    out = Path(path).resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Output directory must be absent or empty: {out}")
    out.mkdir(parents=True, exist_ok=True)
    return out


def _git_state(repo_root: Path) -> dict[str, Any]:
    def run(*args: str) -> str:
        completed = subprocess.run(
            ["git", *args],
            cwd=repo_root,
            check=True,
            capture_output=True,
            text=True,
        )
        return completed.stdout.rstrip()

    status = run("status", "--short")
    return {
        "sha": run("rev-parse", "HEAD"),
        "branch": run("branch", "--show-current"),
        "dirty": bool(status),
        "status_short": status.splitlines(),
    }


def _code_hashes(repo_root: Path) -> dict[str, str | None]:
    paths = {
        "orchestration_cli": Path(__file__).resolve(),
        "four_expert_builder": repo_root / "analysis/ml/spd_four_expert_tables.py",
        "label_enrichment": repo_root / "analysis/ml/spd_label_enrichment.py",
        "receptor_mapping_module": repo_root / "analysis/ml/spd_receptor_mapping.py",
        "model_ready_deduplication": repo_root / "analysis/ml/build_ml_dataset.py",
        "mechanism_projection": (
            repo_root / "analysis/ml/mechanism_label_projection.py"
        ),
        "tissue_projection": (
            repo_root / "analysis/ml/tissue_expression_projection.py"
        ),
    }
    return {
        name: _sha256(path) if path.is_file() else None for name, path in paths.items()
    }


def _file_record(path: Path) -> dict[str, Any]:
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": _sha256(path),
    }


def _resolve_evidence_path(root: Path, relative_path: object) -> Path:
    if not isinstance(relative_path, str) or not relative_path.strip():
        raise ContractViolation("Contract evidence path must be a nonempty string")
    declared = Path(relative_path)
    if declared.is_absolute():
        raise ContractViolation(
            f"Contract evidence paths must be relative to --evidence-root: {declared}"
        )
    resolved = (root / declared).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ContractViolation(
            f"Contract evidence escapes --evidence-root: {declared}"
        ) from exc
    return resolved


def _prepared_receptor_chains(path: Path) -> set[str]:
    chains: set[str] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21:
                chain = line[21].strip()
                if chain:
                    chains.add(chain)
    return chains


def _validate_9i52_ligand_site_evidence(path: Path) -> dict[str, Any]:
    review = pd.read_csv(path, low_memory=False)
    required = {
        "pdb_id",
        "author_chain_id",
        "candidate_primary_accession",
        "candidate_gene_symbols",
        "chain_level_classification",
        "docking_site_evidence_json",
    }
    missing = sorted(required - set(review.columns))
    if missing:
        raise ContractViolation(
            "9I52 audit mapping review lacks ligand-site columns: " + ", ".join(missing)
        )
    pdb = _normalized_pdb(review["pdb_id"])
    chain = review["author_chain_id"].fillna("").astype(str).str.strip().str.upper()
    accession = (
        review["candidate_primary_accession"]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
    )
    genes = (
        review["candidate_gene_symbols"]
        .fillna("")
        .astype(str)
        .str.upper()
        .str.split(r"[;,| ]+", regex=True)
    )
    classification = (
        review["chain_level_classification"].fillna("").astype(str).str.strip()
    )
    candidate = (
        pdb.eq("9I52")
        & chain.eq("R")
        & accession.eq("P21728")
        & genes.map(lambda values: "DRD1" in values)
        & classification.eq("probable_true_target_site_segment")
    )
    if not bool(candidate.any()):
        raise ContractViolation(
            "9I52 ligand-site evidence does not identify chain R as DRD1/P21728"
        )

    contact_confirmed = False
    for value in review.loc[candidate, "docking_site_evidence_json"].dropna():
        try:
            payload = json.loads(str(value))
        except json.JSONDecodeError:
            continue
        for contact in payload.get("segment_contacts", []):
            if (
                str(contact.get("chain_id", "")).strip().upper() == "R"
                and str(contact.get("reported_accession", "")).strip().upper()
                == "P21728"
                and int(contact.get("ligand_atoms_within_5A", 0)) > 0
            ):
                contact_confirmed = True
                break
    if not contact_confirmed:
        raise ContractViolation(
            "9I52 audit row lacks ligand-site contact evidence for chain R/P21728"
        )
    return {
        "path": str(path),
        "matching_rows": int(candidate.sum()),
        "chain": "R",
        "gene": "DRD1",
        "uniprot": "P21728",
        "ligand_site_contact_confirmed": True,
    }


def _validate_contract_evidence(
    contract_path: Path,
    evidence_root: Path,
) -> tuple[dict[str, Any], list[Path]]:
    try:
        payload = json.loads(contract_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ContractViolation(f"Contract is not valid JSON: {contract_path}") from exc
    mappings = payload.get("mappings")
    if not isinstance(mappings, list):
        raise ContractViolation("Contract mappings must be a list")

    records: list[dict[str, Any]] = []
    evidence_paths: list[Path] = []
    review_9i52: Path | None = None
    prepared_9i52: Path | None = None
    mapping_9i52: Mapping[str, Any] | None = None
    for mapping in mappings:
        if not isinstance(mapping, Mapping):
            raise ContractViolation("Every contract mapping must be an object")
        pdb_id = str(mapping.get("pdb_id", "")).strip().upper()
        if pdb_id == "9I52":
            mapping_9i52 = mapping
        evidence = mapping.get("evidence")
        if not isinstance(evidence, list):
            raise ContractViolation(f"{pdb_id}: contract evidence must be a list")
        kinds_seen: set[str] = set()
        for item in evidence:
            if not isinstance(item, Mapping):
                raise ContractViolation(f"{pdb_id}: evidence entry must be an object")
            kind = str(item.get("kind", "")).strip()
            if kind in kinds_seen:
                raise ContractViolation(f"{pdb_id}: duplicate evidence kind {kind}")
            kinds_seen.add(kind)
            path = _resolve_evidence_path(evidence_root, item.get("path"))
            if not path.is_file():
                raise ContractViolation(f"{pdb_id} {kind} evidence is missing: {path}")
            expected_sha = str(item.get("sha256", "")).strip().lower()
            observed_sha = _sha256(path)
            if observed_sha != expected_sha:
                raise ContractViolation(
                    f"{pdb_id} {kind} SHA-256 mismatch: "
                    f"expected {expected_sha}, observed {observed_sha}"
                )
            evidence_paths.append(path)
            record: dict[str, Any] = {
                "pdb_id": pdb_id,
                "kind": kind,
                "declared_path": str(item.get("path")),
                "resolved_path": str(path),
                "sha256": observed_sha,
                "size_bytes": path.stat().st_size,
            }
            if kind == "prepared_receptor":
                chains = sorted(_prepared_receptor_chains(path))
                record["prepared_chains"] = chains
                if pdb_id == "9I52":
                    prepared_9i52 = path
            if kind == "audit_mapping_review" and pdb_id == "9I52":
                review_9i52 = path
            records.append(record)

    if mapping_9i52 is None or review_9i52 is None or prepared_9i52 is None:
        raise ContractViolation(
            "Contract lacks complete 9I52 mapping-review/prepared-receptor evidence"
        )
    expected_9i52 = {
        "selected_chain": "R",
        "site_chain": "R",
        "revised_target_id": "DRD1",
        "revised_target_gene": "DRD1",
        "revised_target_uniprot": "P21728",
    }
    mismatches = {
        field: mapping_9i52.get(field)
        for field, expected in expected_9i52.items()
        if mapping_9i52.get(field) != expected
    }
    if mismatches:
        raise ContractViolation(f"9I52 contract identity/chain mismatch: {mismatches}")
    chains_9i52 = _prepared_receptor_chains(prepared_9i52)
    if "R" not in chains_9i52:
        raise ContractViolation(
            f"9I52 prepared receptor lacks selected/site chain R: {prepared_9i52}"
        )
    ligand_site = _validate_9i52_ligand_site_evidence(review_9i52)
    return (
        {
            "evidence_root": str(evidence_root),
            "files": records,
            "all_declared_hashes_match": True,
            "9i52_prepared_chain_r_present": True,
            "9i52_ligand_site": ligand_site,
        },
        list(dict.fromkeys(evidence_paths)),
    )


def _bool_series(series: pd.Series, *, column: str) -> pd.Series:
    if pd.api.types.is_bool_dtype(series.dtype):
        return series.fillna(False).astype(bool)
    normalized = series.fillna("").astype(str).str.strip().str.lower()
    valid = normalized.isin({"true", "false", "1", "0"})
    if not bool(valid.all()):
        examples = sorted(normalized.loc[~valid].unique().tolist())[:5]
        raise ContractViolation(
            f"{column} must contain only booleans; invalid values: {examples}"
        )
    return normalized.isin({"true", "1"})


def _normalized_pdb(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.upper()


def _value_counts(series: pd.Series) -> dict[str, int]:
    values = series.astype("object").where(series.notna(), "<missing>").astype(str)
    return {
        key: int(value) for key, value in values.value_counts().sort_index().items()
    }


def _semantic_equal(left: pd.Series, right: pd.Series) -> pd.Series:
    left_missing = left.isna()
    right_missing = right.isna()
    left_text = left.fillna("").astype(str).str.strip()
    right_text = right.fillna("").astype(str).str.strip()
    left_numeric = pd.to_numeric(left_text, errors="coerce")
    right_numeric = pd.to_numeric(right_text, errors="coerce")
    numeric_equal = (
        left_numeric.notna() & right_numeric.notna() & left_numeric.eq(right_numeric)
    )
    return (left_missing & right_missing) | (
        ~left_missing & ~right_missing & (left_text.eq(right_text) | numeric_equal)
    )


def _label_counts(frame: pd.DataFrame) -> dict[str, dict[str, int]]:
    counts: dict[str, dict[str, int]] = {}
    for column in LABEL_COLUMNS:
        if column not in frame:
            continue
        numeric = pd.to_numeric(frame[column], errors="coerce")
        counts[column] = {
            "positive": int(numeric.eq(1).sum()),
            "negative": int(numeric.eq(0).sum()),
            "excluded": int(numeric.eq(-1).sum()),
            "unknown": int(numeric.isna().sum()),
        }
    return counts


def _label_change_counts(
    legacy: pd.DataFrame,
    strict: pd.DataFrame,
    mask: pd.Series,
) -> dict[str, int]:
    counts: dict[str, int] = {}
    for column in LABEL_COLUMNS:
        if column not in legacy and column not in strict:
            continue
        before = (
            legacy[column] if column in legacy else pd.Series(pd.NA, index=legacy.index)
        )
        after = (
            strict[column] if column in strict else pd.Series(pd.NA, index=strict.index)
        )
        before_numeric = pd.to_numeric(before.loc[mask], errors="coerce")
        after_numeric = pd.to_numeric(after.loc[mask], errors="coerce")
        counts[column] = int(
            (
                (before_numeric.eq(0) & after_numeric.eq(1))
                | (before_numeric.eq(1) & after_numeric.eq(0))
            ).sum()
        )
    return counts


def _label_transition_counts(
    legacy: pd.DataFrame,
    strict: pd.DataFrame,
    mask: pd.Series,
) -> dict[str, dict[str, int]]:
    """Separate binary flips from known/unknown label-state changes."""

    counts: dict[str, dict[str, int]] = {}
    for column in LABEL_COLUMNS:
        if column not in legacy and column not in strict:
            continue
        before = pd.to_numeric(
            legacy[column].loc[mask]
            if column in legacy
            else pd.Series(pd.NA, index=legacy.index[mask]),
            errors="coerce",
        )
        after = pd.to_numeric(
            strict[column].loc[mask]
            if column in strict
            else pd.Series(pd.NA, index=strict.index[mask]),
            errors="coerce",
        )
        binary_flips = (before.eq(0) & after.eq(1)) | (before.eq(1) & after.eq(0))
        known_to_unknown = before.isin([0, 1]) & after.isna()
        unknown_to_known = before.isna() & after.isin([0, 1])
        counts[column] = {
            "categorical_changes": int((~_semantic_equal(before, after)).sum()),
            "binary_flips": int(binary_flips.sum()),
            "to_unknown": int(known_to_unknown.sum()),
            "from_unknown": int(unknown_to_known.sum()),
        }
    return counts


def _read_expert_tables(out_dir: Path) -> dict[str, pd.DataFrame]:
    tables: dict[str, pd.DataFrame] = {}
    for expert, filename in EXPERT_TABLES.items():
        path = out_dir / filename
        if not path.is_file():
            raise ContractViolation(
                f"Builder did not write expected {expert} table: {path}"
            )
        tables[expert] = pd.read_csv(path, low_memory=False)
    return tables


def _read_model_ready_tables(out_dir: Path) -> dict[str, pd.DataFrame | None]:
    tables: dict[str, pd.DataFrame | None] = {}
    for expert in EXPERT_TABLES:
        path = out_dir / "model_ready" / f"ml_spd_{expert}_model_ready.csv"
        tables[expert] = pd.read_csv(path, low_memory=False) if path.is_file() else None
    return tables


def _comparison_source(
    tables: Mapping[str, pd.DataFrame],
) -> pd.DataFrame:
    source = tables["binding"].copy()
    for expert, label_column in zip(EXPERT_TABLES, LABEL_COLUMNS, strict=True):
        expert_table = tables[expert]
        if label_column in expert_table.columns:
            source[label_column] = expert_table[label_column]
    return source


def _model_ready_pdb_counts(
    tables: Mapping[str, pd.DataFrame | None],
) -> dict[str, dict[str, int]]:
    counts: dict[str, dict[str, int]] = {}
    for expert, frame in tables.items():
        if frame is None:
            counts[expert] = {}
            continue
        if "pdb_id" not in frame.columns:
            raise ContractViolation(f"{expert} model-ready output lacks pdb_id")
        counts[expert] = {
            pdb: int(count)
            for pdb, count in _normalized_pdb(frame["pdb_id"])
            .value_counts()
            .sort_index()
            .items()
        }
    return counts


def _strict_model_ready_checks(
    strict_raw: Mapping[str, pd.DataFrame],
    strict_model_ready: Mapping[str, pd.DataFrame | None],
) -> tuple[bool, bool, dict[str, dict[str, int]]]:
    from analysis.ml.build_ml_dataset import (
        add_standard_ml_columns,
        deduplicate_ml_rows,
    )
    from analysis.ml.labels import binary_label_series

    counts = _model_ready_pdb_counts(strict_model_ready)
    unresolved_excluded = all(
        pdb not in expert_counts
        for expert_counts in counts.values()
        for pdb in ("6HUJ", "8HCQ")
    )
    eligible_retained = True
    label_by_expert = dict(zip(EXPERT_TABLES, LABEL_COLUMNS, strict=True))
    for expert, label_column in label_by_expert.items():
        raw = strict_raw[expert]
        label = binary_label_series(raw[label_column])
        eligible = _bool_series(
            raw["spd_receptor_mapping_strict_eligible"],
            column="spd_receptor_mapping_strict_eligible",
        )
        expected_frame = raw.loc[label.notna() & eligible].copy()
        if not expected_frame.empty:
            expected_frame[label_column] = label.loc[expected_frame.index].astype(int)
            expected_frame = deduplicate_ml_rows(
                add_standard_ml_columns(expected_frame),
                label_column,
            )
        expected_counts = {
            pdb: int(count)
            for pdb, count in _normalized_pdb(
                expected_frame.get("pdb_id", pd.Series(dtype="object"))
            )
            .value_counts()
            .sort_index()
            .items()
        }
        if counts[expert] != expected_counts:
            eligible_retained = False
    return unresolved_excluded, eligible_retained, counts


def _assert_expert_alignment(
    tables: Mapping[str, pd.DataFrame],
    *,
    expected_rows: int,
    mode: str,
) -> None:
    for expert, frame in tables.items():
        if len(frame) != expected_rows:
            raise ContractViolation(
                f"{mode} {expert} table multiplied or dropped rows: "
                f"expected {expected_rows}, observed {len(frame)}"
            )
    reference = tables["binding"]
    identity_cols = [
        column
        for column in ("pdb_id", "ligand_base", "drug_id")
        if column in reference.columns
    ]
    if not identity_cols:
        raise ContractViolation("Expert outputs lack a stable identity column")
    for expert, frame in tables.items():
        if expert == "binding":
            continue
        for column in identity_cols:
            if column not in frame.columns or not reference[column].equals(
                frame[column]
            ):
                raise ContractViolation(
                    f"{mode} {expert} row order differs from binding for {column}"
                )


def _comparison_frame(
    legacy: pd.DataFrame,
    strict: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.Series]:
    if len(legacy) != len(strict):
        raise ContractViolation(
            f"Legacy/strict row counts differ: {len(legacy)} versus {len(strict)}"
        )
    missing = [column for column in MAPPING_AUDIT_COLUMNS if column not in strict]
    if missing:
        raise ContractViolation(
            "Strict output lacks required receptor mapping audit columns: "
            + ", ".join(missing)
        )
    if "pdb_id" not in strict:
        raise ContractViolation("Strict output lacks pdb_id")

    eligible = _bool_series(
        strict["spd_receptor_mapping_strict_eligible"],
        column="spd_receptor_mapping_strict_eligible",
    )
    changed = _bool_series(
        strict["spd_receptor_mapping_changed"],
        column="spd_receptor_mapping_changed",
    )
    status = strict["spd_receptor_mapping_status"].fillna("").astype(str).str.strip()
    exclusion = (
        strict["spd_receptor_mapping_exclusion_reason"]
        .fillna("")
        .astype(str)
        .str.strip()
    )
    affected = changed | ~eligible | exclusion.ne("")

    comparison = pd.DataFrame(
        {
            "source_row_number": pd.RangeIndex(start=0, stop=len(strict), step=1),
            "pdb_id": _normalized_pdb(strict["pdb_id"]),
        }
    )
    for column in (
        "ligand_base",
        "drug_id",
        "generic_name",
        "display_name",
    ):
        if column in strict:
            comparison[column] = strict[column]
    for column in TARGET_COLUMNS:
        legacy_fallback = f"spd_receptor_mapping_legacy_{column}"
        revised_fallback = f"spd_receptor_mapping_revised_{column}"
        if legacy_fallback in strict:
            comparison[f"before_{column}"] = strict[legacy_fallback]
        elif column in legacy:
            comparison[f"before_{column}"] = legacy[column]
        else:
            comparison[f"before_{column}"] = pd.NA
        if revised_fallback in strict:
            comparison[f"after_{column}"] = strict[revised_fallback]
        elif column in strict:
            comparison[f"after_{column}"] = strict[column]
        else:
            comparison[f"after_{column}"] = pd.NA
    for column in (
        "spd_drug_id",
        "spd_target_id",
        "spd_ac50_uM",
        "spd_activity_relation",
        "spd_label_status",
        "spd_missing_reason",
        "spd_assay_ids",
        "spd_source",
        "free_cmax_um",
        "exposure_margin",
    ):
        if column in legacy or column in strict:
            comparison[f"before_{column}"] = (
                legacy[column]
                if column in legacy
                else pd.Series(pd.NA, index=legacy.index)
            )
            comparison[f"after_{column}"] = (
                strict[column]
                if column in strict
                else pd.Series(pd.NA, index=strict.index)
            )
    for column in MAPPING_AUDIT_COLUMNS:
        comparison[column] = strict[column]
    for column in LABEL_COLUMNS:
        if column in legacy or column in strict:
            comparison[f"before_{column}"] = (
                legacy[column]
                if column in legacy
                else pd.Series(pd.NA, index=legacy.index)
            )
            comparison[f"after_{column}"] = (
                strict[column]
                if column in strict
                else pd.Series(pd.NA, index=strict.index)
            )
            comparison[f"{column}_changed"] = ~_semantic_equal(
                comparison[f"before_{column}"],
                comparison[f"after_{column}"],
            )
    comparison["strict_eligible"] = eligible
    comparison["mapping_affected"] = affected
    comparison["spd_receptor_mapping_status"] = status
    return comparison.loc[affected].reset_index(drop=True), affected


def _unrelated_rows_unchanged(
    legacy: pd.DataFrame,
    strict: pd.DataFrame,
    affected: pd.Series,
) -> tuple[bool, list[str]]:
    outside = ~affected
    ignored = set(MAPPING_AUDIT_COLUMNS)
    common = [
        column
        for column in legacy.columns
        if column in strict.columns and column not in ignored
    ]
    changed_columns: list[str] = []
    for column in common:
        before = legacy.loc[outside, column].reset_index(drop=True)
        after = strict.loc[outside, column].reset_index(drop=True)
        if not bool(_semantic_equal(before, after).all()):
            changed_columns.append(column)
    return not changed_columns, changed_columns


def _write_comparison_markdown(
    path: Path,
    *,
    summary: pd.DataFrame,
    checks: Mapping[str, Any],
    label_counts: Mapping[str, Any],
    eligibility: Mapping[str, Any],
) -> None:
    lines = [
        "# SPD receptor mapping before/after",
        "",
        f"Contract status: **{'PASS' if all(checks.values()) else 'FAIL'}**",
        "",
        "## Affected receptor scope",
        "",
        "| PDB | Legacy target | Revised status | Rows affected | Binding label changes | Exposure label changes | Strict eligibility |",
        "|---|---|---|---:|---|---|---|",
    ]
    lines.extend(
        "| "
        + " | ".join(
            [
                str(row.pdb_id),
                str(row.legacy_target),
                str(row.revised_status),
                str(int(row.rows_affected)),
                str(row.binding_label_changes),
                str(row.exposure_label_changes),
                str(row.strict_eligibility),
            ]
        )
        + " |"
        for row in summary.itertuples(index=False)
    )
    lines.extend(
        [
            "",
            "## Strict eligibility",
            "",
            f"- Eligible rows: {eligibility['eligible_rows']}",
            f"- Ineligible rows: {eligibility['ineligible_rows']}",
            f"- Status counts: `{json.dumps(eligibility['status_counts'], sort_keys=True)}`",
            (
                "- Exclusion reasons: "
                f"`{json.dumps(eligibility['exclusion_reason_counts'], sort_keys=True)}`"
            ),
            "",
            "## Label counts",
            "",
            "| Expert label | Mode | Positive | Negative | Excluded | Unknown |",
            "|---|---|---:|---:|---:|---:|",
        ]
    )
    for mode in ("legacy", "strict"):
        for label, values in label_counts[mode].items():
            lines.append(
                f"| {label} | {mode} | {values['positive']} | "
                f"{values['negative']} | {values['excluded']} | {values['unknown']} |"
            )
    lines.extend(
        [
            "",
            "## Contract checks",
            "",
        ]
    )
    lines.extend(
        f"- {'PASS' if passed else 'FAIL'}: `{name}`" for name, passed in checks.items()
    )
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def _format_label_transition(
    before: Mapping[str, int],
    after: Mapping[str, int],
    transitions: Mapping[str, int],
) -> str:
    return (
        f"{before['positive']}+/{before['negative']}-/{before['unknown']} unknown"
        f" -> {after['positive']}+/{after['negative']}-/{after['unknown']} unknown"
        f"; {transitions['categorical_changes']} categorical changes"
        f" ({transitions['binary_flips']} binary flips, "
        f"{transitions['to_unknown']} to unknown, "
        f"{transitions['from_unknown']} from unknown)"
    )


def _mapping_summary_frame(
    comparison: pd.DataFrame,
    legacy: pd.DataFrame,
    strict: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    legacy_pdb = _normalized_pdb(legacy["pdb_id"])
    strict_pdb = _normalized_pdb(strict["pdb_id"])
    for pdb_id in sorted(EXPECTED_AFFECTED_PDB_SCOPE):
        legacy_mask = legacy_pdb.eq(pdb_id)
        strict_mask = strict_pdb.eq(pdb_id)
        delta = comparison.loc[comparison["pdb_id"].eq(pdb_id)]
        if delta.empty:
            raise ContractViolation(f"Missing affected mapping summary for {pdb_id}")
        binding_before = _label_counts(legacy.loc[legacy_mask]).get(
            "spd_binding_label", {}
        )
        binding_after = _label_counts(strict.loc[strict_mask]).get(
            "spd_binding_label", {}
        )
        exposure_before = _label_counts(legacy.loc[legacy_mask]).get(
            "spd_exposure_label", {}
        )
        exposure_after = _label_counts(strict.loc[strict_mask]).get(
            "spd_exposure_label", {}
        )
        transitions = _label_transition_counts(legacy, strict, legacy_mask)
        status = str(delta["spd_receptor_mapping_status"].iloc[0])
        revised = delta["after_target_gene"].dropna().astype(str).unique().tolist()
        revised_status = f"{status}:{revised[0]}" if revised else status
        eligible = _bool_series(
            delta["spd_receptor_mapping_strict_eligible"],
            column="spd_receptor_mapping_strict_eligible",
        )
        reasons = sorted(
            {
                value
                for value in delta["spd_receptor_mapping_exclusion_reason"]
                .fillna("")
                .astype(str)
                .str.strip()
                if value
            }
        )
        eligibility = (
            "eligible"
            if bool(eligible.all())
            else "ineligible:" + ";".join(reasons or ["unspecified"])
        )
        rows.append(
            {
                "pdb_id": pdb_id,
                "legacy_target": str(delta["before_target_gene"].iloc[0]),
                "revised_status": revised_status,
                "rows_affected": int(len(delta)),
                "binding_label_changes": _format_label_transition(
                    binding_before,
                    binding_after,
                    transitions["spd_binding_label"],
                ),
                "exposure_label_changes": _format_label_transition(
                    exposure_before,
                    exposure_after,
                    transitions["spd_exposure_label"],
                ),
                "strict_eligibility": eligibility,
            }
        )
    return pd.DataFrame(rows)


def run_receptor_mapping_comparison(
    *,
    spd_table: str | Path,
    spd_panel: str | Path,
    target_map: str | Path,
    ligand_map: str | Path,
    out_dir: str | Path,
    contract: str | Path | None = None,
    evidence_root: str | Path,
    target_metadata: str | Path | None = None,
    mechanism_labels: str | Path | None = None,
    target_expression: str | Path | None = None,
    run_dir: str | Path | None = None,
    expected_affected_counts: Mapping[str, int] | None = AUDITED_PHASE1_AFFECTED_COUNTS,
    builder: Callable[..., dict[str, Any]] | None = None,
) -> dict[str, Any]:
    source_path = _required_file(spd_table, "SPD table")
    panel_path = _required_file(spd_panel, "SPD panel")
    target_map_path = _required_file(target_map, "target map")
    ligand_map_path = _required_file(ligand_map, "ligand map")
    contract_path = _resolve_contract_path(contract)
    evidence_root_path = Path(evidence_root).resolve()
    if not evidence_root_path.is_dir():
        raise FileNotFoundError(f"Evidence root does not exist: {evidence_root_path}")
    evidence_validation, evidence_paths = _validate_contract_evidence(
        contract_path,
        evidence_root_path,
    )
    metadata_path = _optional_file(target_metadata, "target metadata")
    mechanism_path = _optional_file(mechanism_labels, "mechanism labels")
    expression_path = _optional_file(target_expression, "target expression")
    run_path = Path(run_dir).resolve() if run_dir is not None else None
    out = _prepare_out_dir(out_dir)

    input_paths = {
        "spd_table": source_path,
        "spd_panel": panel_path,
        "target_map": target_map_path,
        "ligand_map": ligand_map_path,
        "contract": contract_path,
    }
    for name, path in (
        ("target_metadata", metadata_path),
        ("mechanism_labels", mechanism_path),
        ("target_expression", expression_path),
    ):
        if path is not None:
            input_paths[name] = path
    input_before = {name: _file_record(path) for name, path in input_paths.items()}

    source_rows = len(pd.read_csv(source_path, usecols=lambda _column: True))
    build_kwargs: dict[str, Any] = {
        "spd_panel_path": panel_path,
        "target_map_path": target_map_path,
        "ligand_map_path": ligand_map_path,
        "target_metadata_path": metadata_path,
        "mechanism_label_path": mechanism_path,
        "target_expression_path": expression_path,
        "run_dir": run_path,
        "receptor_mapping_contract_path": contract_path,
    }
    build_manifests: dict[str, Any] = {}
    mode_tables: dict[str, dict[str, pd.DataFrame]] = {}
    model_ready_tables: dict[str, dict[str, pd.DataFrame | None]] = {}
    if builder is None:
        from analysis.ml.spd_four_expert_tables import (
            build_spd_four_expert_tables,
        )

        builder = build_spd_four_expert_tables
    for mode in ("legacy", "strict"):
        mode_out = out / mode
        build_manifests[mode] = builder(
            source_path,
            mode_out,
            receptor_mapping_mode=mode,
            **build_kwargs,
        )
        if _sha256(source_path) != input_before["spd_table"]["sha256"]:
            raise ContractViolation(
                f"Builder modified serialized SPD input in {mode} mode"
            )
        mode_tables[mode] = _read_expert_tables(mode_out)
        model_ready_tables[mode] = _read_model_ready_tables(mode_out)
        _assert_expert_alignment(
            mode_tables[mode],
            expected_rows=source_rows,
            mode=mode,
        )

    legacy_binding = mode_tables["legacy"]["binding"]
    strict_binding = mode_tables["strict"]["binding"]
    legacy_comparison_source = _comparison_source(mode_tables["legacy"])
    strict_comparison_source = _comparison_source(mode_tables["strict"])
    comparison, affected = _comparison_frame(
        legacy_comparison_source,
        strict_comparison_source,
    )
    affected_scope = set(comparison["pdb_id"].dropna().astype(str))
    affected_counts = {
        pdb: int(count)
        for pdb, count in comparison["pdb_id"].value_counts().sort_index().items()
    }
    unrelated_unchanged, unrelated_changed_columns = _unrelated_rows_unchanged(
        legacy_binding,
        strict_binding,
        affected,
    )

    strict_eligible = _bool_series(
        strict_binding["spd_receptor_mapping_strict_eligible"],
        column="spd_receptor_mapping_strict_eligible",
    )
    eligibility = {
        "eligible_rows": int(strict_eligible.sum()),
        "ineligible_rows": int((~strict_eligible).sum()),
        "status_counts": _value_counts(strict_binding["spd_receptor_mapping_status"]),
        "exclusion_reason_counts": _value_counts(
            strict_binding["spd_receptor_mapping_exclusion_reason"]
        ),
    }
    label_counts = {
        "legacy": _label_counts(legacy_comparison_source),
        "strict": _label_counts(strict_comparison_source),
    }
    pdb_9i52 = _normalized_pdb(strict_binding["pdb_id"]).eq("9I52")
    label_counts_9i52 = {
        "legacy": _label_counts(legacy_comparison_source.loc[pdb_9i52]),
        "strict": _label_counts(strict_comparison_source.loc[pdb_9i52]),
    }
    label_changes_9i52 = _label_change_counts(
        legacy_comparison_source,
        strict_comparison_source,
        pdb_9i52,
    )
    label_transitions_9i52 = _label_transition_counts(
        legacy_comparison_source,
        strict_comparison_source,
        pdb_9i52,
    )
    strict_9i52 = strict_binding.loc[pdb_9i52]
    panel_match_mask = pd.to_numeric(
        strict_9i52.get(
            "spd_ac50_uM",
            pd.Series(pd.NA, index=strict_9i52.index),
        ),
        errors="coerce",
    ).notna()
    panel_match_count_9i52 = int(panel_match_mask.sum())
    unmatched_drugs_9i52 = sorted(
        {
            str(value).strip().lower()
            for value in strict_9i52.loc[~panel_match_mask, "drug_id"].dropna()
            if str(value).strip()
        }
    )
    unmatched_rows_9i52: list[dict[str, Any]] = []
    for source_index in strict_9i52.index[~panel_match_mask]:
        before_binding = pd.to_numeric(
            pd.Series(
                [legacy_comparison_source.at[source_index, "spd_binding_label"]]
            ),
            errors="coerce",
        ).iloc[0]
        after_binding = pd.to_numeric(
            pd.Series(
                [strict_comparison_source.at[source_index, "spd_binding_label"]]
            ),
            errors="coerce",
        ).iloc[0]
        unmatched_rows_9i52.append(
            {
                "source_row_number": int(source_index),
                "ligand_base": str(strict_binding.at[source_index, "ligand_base"]),
                "drug_id": str(strict_binding.at[source_index, "drug_id"]),
                "legacy_target_gene": str(
                    legacy_binding.at[source_index, "target_gene"]
                ),
                "strict_target_gene": str(
                    strict_binding.at[source_index, "target_gene"]
                ),
                "legacy_binding_label": (
                    int(before_binding) if pd.notna(before_binding) else None
                ),
                "strict_binding_label": (
                    int(after_binding) if pd.notna(after_binding) else None
                ),
                "strict_spd_missing_reason": (
                    None
                    if "spd_missing_reason" not in strict_binding
                    or pd.isna(
                        strict_binding.at[source_index, "spd_missing_reason"]
                    )
                    else str(
                        strict_binding.at[source_index, "spd_missing_reason"]
                    )
                ),
            }
        )
    (
        unresolved_model_ready_excluded,
        eligible_9i52_model_ready_retained,
        strict_model_ready_pdb_counts,
    ) = _strict_model_ready_checks(
        mode_tables["strict"],
        model_ready_tables["strict"],
    )
    legacy_model_ready_pdb_counts = _model_ready_pdb_counts(
        model_ready_tables["legacy"]
    )
    outside_contract = ~_normalized_pdb(strict_binding["pdb_id"]).isin(
        EXPECTED_AFFECTED_PDB_SCOPE
    )

    checks: dict[str, bool] = {
        "identical_serialized_input_path": True,
        "all_inputs_preserved": all(
            _sha256(path) == input_before[name]["sha256"]
            for name, path in input_paths.items()
        ),
        "all_contract_evidence_preserved": all(
            _sha256(path)
            == next(
                record["sha256"]
                for record in evidence_validation["files"]
                if record["resolved_path"] == str(path)
            )
            for path in evidence_paths
        ),
        "legacy_no_row_multiplication": len(legacy_binding) == source_rows,
        "strict_no_row_multiplication": len(strict_binding) == source_rows,
        "legacy_strict_row_count_equal": len(legacy_binding) == len(strict_binding),
        "affected_pdb_scope_exact": affected_scope == EXPECTED_AFFECTED_PDB_SCOPE,
        "unrelated_pdb_rows_unchanged": unrelated_unchanged,
        "strict_eligibility_populated": bool(strict_eligible.notna().all()),
        "out_of_contract_rows_strict_eligible": bool(
            strict_eligible.loc[outside_contract].all()
        ),
        "strict_model_ready_excludes_6huj_8hcq": unresolved_model_ready_excluded,
        "strict_model_ready_retains_eligible_9i52": (
            eligible_9i52_model_ready_retained
        ),
    }
    if expected_affected_counts is not None:
        normalized_expected = {
            str(pdb).strip().upper(): int(count)
            for pdb, count in expected_affected_counts.items()
        }
        checks["audited_phase1_affected_counts_exact"] = (
            affected_counts == normalized_expected
        )
        if normalized_expected == AUDITED_PHASE1_AFFECTED_COUNTS:
            checks["audited_phase1_9i52_label_counts_exact"] = all(
                label_counts_9i52[mode].get(label) == expected
                for mode, mode_counts in AUDITED_PHASE1_9I52_LABEL_COUNTS.items()
                for label, expected in mode_counts.items()
            )
            checks["audited_phase1_9i52_binary_flips_exact"] = all(
                label_changes_9i52.get(label) == expected
                for label, expected in AUDITED_PHASE1_9I52_BINARY_FLIPS.items()
            )
            checks["audited_phase1_9i52_panel_matches_exact"] = (
                panel_match_count_9i52 == AUDITED_PHASE1_9I52_PANEL_MATCHES
            )
            checks["audited_phase1_9i52_unmatched_drugs_exact"] = (
                set(unmatched_drugs_9i52) == AUDITED_PHASE1_9I52_UNMATCHED_DRUGS
            )
            checks["audited_phase1_9i52_unmatched_labels_not_transferred"] = bool(
                pd.to_numeric(
                    strict_comparison_source.loc[
                        strict_9i52.index[~panel_match_mask],
                        "spd_binding_label",
                    ],
                    errors="coerce",
                ).isna().all()
            )
    checks = {name: bool(passed) for name, passed in checks.items()}

    row_delta_path = out / "receptor_mapping_row_deltas.csv"
    comparison.to_csv(row_delta_path, index=False)
    summary_frame = _mapping_summary_frame(
        comparison,
        legacy_comparison_source,
        strict_comparison_source,
    )
    comparison_path = out / "receptor_mapping_before_after.csv"
    summary_frame.to_csv(comparison_path, index=False)
    markdown_path = out / "receptor_mapping_before_after.md"
    _write_comparison_markdown(
        markdown_path,
        summary=summary_frame,
        checks=checks,
        label_counts=label_counts,
        eligibility=eligibility,
    )

    repo_root = _repo_root()
    manifest: dict[str, Any] = {
        "status": "pass" if all(checks.values()) else "fail",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": {
            "executable": sys.executable,
            "version": platform.python_version(),
            "platform": platform.platform(),
            "pandas_version": pd.__version__,
        },
        "git": _git_state(repo_root),
        "code_sha256": _code_hashes(repo_root),
        "inputs": input_before,
        "contract_evidence_validation": evidence_validation,
        "serialized_input_rows": source_rows,
        "modes": {
            "legacy": {
                "out_dir": str(out / "legacy"),
                "builder_manifest": build_manifests["legacy"],
            },
            "strict": {
                "out_dir": str(out / "strict"),
                "builder_manifest": build_manifests["strict"],
            },
        },
        "affected_pdb_scope": sorted(affected_scope),
        "expected_affected_pdb_scope": sorted(EXPECTED_AFFECTED_PDB_SCOPE),
        "affected_counts": affected_counts,
        "expected_affected_counts": (
            {
                str(pdb).strip().upper(): int(count)
                for pdb, count in expected_affected_counts.items()
            }
            if expected_affected_counts is not None
            else None
        ),
        "strict_eligibility": eligibility,
        "label_counts": label_counts,
        "model_ready_pdb_counts": {
            "legacy": legacy_model_ready_pdb_counts,
            "strict": strict_model_ready_pdb_counts,
        },
        "audited_phase1_9i52": {
            "label_counts": label_counts_9i52,
            "binary_flip_counts": label_changes_9i52,
            "label_transition_counts": label_transitions_9i52,
            "panel_match_rows": panel_match_count_9i52,
            "unmatched_drugs": unmatched_drugs_9i52,
            "unmatched_panel_rows": unmatched_rows_9i52,
        },
        "unrelated_changed_columns": unrelated_changed_columns,
        "checks": checks,
        "artifacts": {
            "comparison_csv": str(comparison_path),
            "comparison_markdown": str(markdown_path),
            "row_deltas_csv": str(row_delta_path),
        },
        "preservation_note": (
            "The same serialized SPD input path was passed to both builds; "
            "all input hashes were rechecked after both builds."
        ),
        "scope_note": (
            "This command rebuilds only legacy/strict SPD expert-table comparisons. "
            "It does not merge external four-state evidence, mutate frozen Phase1 "
            "tables, or train models."
        ),
    }
    manifest_path = out / "mapping_provenance_manifest.json"
    manifest["artifacts"]["manifest"] = str(manifest_path)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    if not all(checks.values()):
        failures = [name for name, passed in checks.items() if not passed]
        raise ContractViolation(
            "Receptor mapping comparison contract failed: " + ", ".join(failures)
        )
    return manifest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Rebuild identical SPD input under legacy and strict receptor mapping "
            "and write a fail-closed before/after contract."
        )
    )
    parser.add_argument("--spd-table", required=True, type=Path)
    parser.add_argument("--spd-panel", required=True, type=Path)
    parser.add_argument("--target-map", required=True, type=Path)
    parser.add_argument("--ligand-map", required=True, type=Path)
    parser.add_argument("--target-metadata", type=Path)
    parser.add_argument("--mechanism-labels", type=Path)
    parser.add_argument("--target-expression", type=Path)
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--contract", type=Path)
    parser.add_argument(
        "--evidence-root",
        required=True,
        type=Path,
        help=(
            "Root used to resolve and hash every contract evidence path; "
            "use the original repository root for the audited Phase1 contract."
        ),
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument(
        "--expected-profile",
        choices=("audited-phase1", "scope-only"),
        default="audited-phase1",
        help=(
            "audited-phase1 enforces 9I52=388, 6HUJ=570, 8HCQ=278; "
            "scope-only still enforces the exact three-PDB affected scope."
        ),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    expected_counts = (
        AUDITED_PHASE1_AFFECTED_COUNTS
        if args.expected_profile == "audited-phase1"
        else None
    )
    manifest = run_receptor_mapping_comparison(
        spd_table=args.spd_table,
        spd_panel=args.spd_panel,
        target_map=args.target_map,
        ligand_map=args.ligand_map,
        target_metadata=args.target_metadata,
        mechanism_labels=args.mechanism_labels,
        target_expression=args.target_expression,
        run_dir=args.run_dir,
        contract=args.contract,
        evidence_root=args.evidence_root,
        out_dir=args.out_dir,
        expected_affected_counts=expected_counts,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
