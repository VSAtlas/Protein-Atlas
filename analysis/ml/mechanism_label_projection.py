from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_MECHANISM_LABEL_CANDIDATES = [
    Path("negative_evidence/mechanism_four_state_labels_deep_sources_v7_adrecs_positive_source_balanced.csv"),
    Path("negative_evidence/mechanism_four_state_labels_paper_sources_v9_source_balanced.csv"),
    Path("negative_evidence/mechanism_four_state_labels_deep_sources_source_balanced.csv"),
    Path("negative_evidence/mechanism_four_state_labels.csv"),
    Path("data/pilotstudy/negative_evidence/mechanism_four_state_labels_deep_sources_v7_adrecs_positive_source_balanced.csv"),
    Path("data/pilotstudy/negative_evidence/mechanism_four_state_labels_paper_sources_v9_source_balanced.csv"),
    Path("data/pilotstudy/negative_evidence/mechanism_four_state_labels_deep_sources_source_balanced.csv"),
    Path("data/pilotstudy/negative_evidence/mechanism_four_state_labels.csv"),
]

EVIDENCE_NUMERIC_COLS = [
    "mechanism_graph_score",
    "mechanism_path_count",
    "drug_adr_known",
    "target_adr_known",
    "target_pathway_adr_link",
    "triad_complete",
    "target_adr_evidence",
    "pathway_evidence",
]
PROVENANCE_COLS = [
    "label_source",
    "source_family",
    "upstream_source",
    "source_label_policy",
    "four_state_label_status",
    "mechanism_label_status",
    "negative_evidence_type",
    "negative_source",
    "negative_confidence",
]
LABEL_COLS = ["mechanism_ml_label", "four_state_ml_label", "mechanism_label"]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _lower(value: Any) -> str:
    return _clean(value).lower()


def _upper(value: Any) -> str:
    return _clean(value).upper()


def _resolve_existing(explicit: str | Path | None, candidates: list[Path], run_dir: Path | None = None) -> Path | None:
    if explicit is not None:
        path = Path(explicit)
        return path if path.exists() else None
    roots = []
    if run_dir is not None:
        roots.append(run_dir)
    roots.extend([Path.cwd(), _repo_root()])
    for root in roots:
        for rel in candidates:
            path = rel if rel.is_absolute() else root / rel
            if path.exists():
                return path
    return None


def _first_nonempty_frame(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    out = pd.Series("", index=df.index, dtype="object")
    for col in cols:
        if col not in df.columns:
            continue
        values = df[col].map(_clean)
        out = out.where(out.astype(str).str.len().gt(0), values)
    return out


def _join_unique(values: pd.Series) -> str:
    seen: list[str] = []
    for value in values:
        text = _clean(value)
        if text and text not in seen:
            seen.append(text)
    return ";".join(seen)


def _label_from_row(row: pd.Series) -> float | pd.NA:
    for col in LABEL_COLS:
        if col not in row.index:
            continue
        value = pd.to_numeric(row.get(col), errors="coerce")
        if pd.notna(value) and float(value) in {-1.0, 0.0, 1.0}:
            return float(value)
    return pd.NA


def _collapse_labels(label_table: pd.DataFrame) -> pd.DataFrame:
    work = label_table.copy()
    work["_mechanism_ligand_key"] = _first_nonempty_frame(
        work, ["ligand_base", "drug_id", "generic_name", "display_name"]
    ).map(_lower)
    work["_mechanism_target_key"] = _first_nonempty_frame(
        work, ["target_gene", "gene_symbol", "target_id", "target_uniprot", "pdb_id"]
    ).map(_upper)
    work = work[
        work["_mechanism_ligand_key"].astype(str).str.len().gt(0)
        & work["_mechanism_target_key"].astype(str).str.len().gt(0)
    ].copy()
    if work.empty:
        return pd.DataFrame(columns=["_mechanism_ligand_key", "_mechanism_target_key"])
    work["_resolved_mechanism_label"] = work.apply(_label_from_row, axis=1)

    rows: list[dict[str, Any]] = []
    for (ligand_key, target_key), group in work.groupby(["_mechanism_ligand_key", "_mechanism_target_key"], dropna=False):
        labels = sorted(
            {float(value) for value in pd.to_numeric(group["_resolved_mechanism_label"], errors="coerce").dropna()}
        )
        if not labels:
            mechanism_label = pd.NA
            mechanism_ml_label = pd.NA
            status = "unknown_projected_no_resolved_label"
        elif -1.0 in labels or (0.0 in labels and 1.0 in labels):
            mechanism_label = -1
            mechanism_ml_label = pd.NA
            status = "excluded_projected_conflicting_mechanism_sources"
        else:
            mechanism_label = int(labels[-1])
            mechanism_ml_label = mechanism_label
            status = "projected_strict_positive" if mechanism_label == 1 else "projected_measured_or_reliable_negative"
        row: dict[str, Any] = {
            "_mechanism_ligand_key": ligand_key,
            "_mechanism_target_key": target_key,
            "mechanism_label": mechanism_label,
            "mechanism_ml_label": mechanism_ml_label,
            "mechanism_label_status": status,
            "mechanism_label_projection_key": "ligand_base_or_drug_name+target_gene",
            "mechanism_label_source": _join_unique(group.get("label_source", pd.Series(dtype=object))),
            "mechanism_label_projection_source_rows": int(len(group)),
        }
        for col in PROVENANCE_COLS:
            if col in group.columns and col not in {"mechanism_label_status"}:
                row[f"projected_{col}"] = _join_unique(group[col])
        for col in EVIDENCE_NUMERIC_COLS:
            if col in group.columns:
                numeric = pd.to_numeric(group[col], errors="coerce")
                row[col] = numeric.max() if numeric.notna().any() else pd.NA
        rows.append(row)
    return pd.DataFrame(rows)


def project_mechanism_labels_for_run_master(
    source: pd.DataFrame,
    *,
    mechanism_label_path: str | Path | None = None,
    run_dir: str | Path | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Project curated mechanism labels onto a run table with stable ligand/target keys.

    This is a convenience join, not label invention. It only transfers labels
    from existing mechanism evidence tables and leaves non-overlapping rows
    unknown.
    """

    summary: dict[str, Any] = {
        "status": "skipped",
        "reason": "not_needed",
        "input_rows": int(len(source)),
    }
    if "mechanism_ml_label" in source.columns and pd.to_numeric(source["mechanism_ml_label"], errors="coerce").notna().any():
        return source, summary
    if not ({"ligand_base", "drug_id", "generic_name", "display_name"} & set(source.columns)):
        summary.update({"reason": "missing_ligand_key_columns"})
        return source, summary
    if not ({"target_gene", "gene_symbol", "target_id", "target_uniprot", "pdb_id"} & set(source.columns)):
        summary.update({"reason": "missing_target_key_columns"})
        return source, summary

    label_path = _resolve_existing(mechanism_label_path, DEFAULT_MECHANISM_LABEL_CANDIDATES, Path(run_dir) if run_dir else None)
    if label_path is None:
        summary.update({"reason": "missing_mechanism_label_source"})
        return source, summary

    label_table = pd.read_csv(label_path, low_memory=False)
    collapsed = _collapse_labels(label_table)
    if collapsed.empty:
        summary.update({"reason": "mechanism_label_source_has_no_projectable_keys", "mechanism_label_path": str(label_path)})
        return source, summary

    out = source.copy()
    out["_mechanism_ligand_key"] = _first_nonempty_frame(
        out, ["ligand_base", "drug_id", "generic_name", "display_name"]
    ).map(_lower)
    out["_mechanism_target_key"] = _first_nonempty_frame(
        out, ["target_gene", "gene_symbol", "target_id", "target_uniprot", "pdb_id"]
    ).map(_upper)
    out = out.merge(collapsed, on=["_mechanism_ligand_key", "_mechanism_target_key"], how="left", suffixes=("", "_projected"))
    if "mechanism_ml_label_projected" in out.columns:
        out["mechanism_ml_label"] = out.get("mechanism_ml_label", pd.Series(pd.NA, index=out.index))
        out["mechanism_ml_label"] = out["mechanism_ml_label"].where(out["mechanism_ml_label"].notna(), out["mechanism_ml_label_projected"])
        out = out.drop(columns=["mechanism_ml_label_projected"])
    out["mechanism_label_projected_available"] = pd.to_numeric(out.get("mechanism_ml_label"), errors="coerce").notna()
    out["mechanism_label_projection_source"] = str(label_path)
    out = out.loc[:, ~out.columns.duplicated()].copy()

    projected = pd.to_numeric(out.get("mechanism_ml_label"), errors="coerce")
    summary.update(
        {
            "status": "projected",
            "reason": "joined_existing_mechanism_label_table",
            "mechanism_label_path": str(label_path),
            "rows": int(len(out)),
            "source_label_rows": int(len(label_table)),
            "source_projectable_pairs": int(len(collapsed)),
            "projected_label_rows": int(projected.notna().sum()),
            "projected_positive_rows": int(projected.eq(1).sum()),
            "projected_negative_rows": int(projected.eq(0).sum()),
            "unique_projected_pairs": int(
                out.loc[projected.notna(), ["_mechanism_ligand_key", "_mechanism_target_key"]]
                .drop_duplicates()
                .shape[0]
            ),
        }
    )
    return out, summary
