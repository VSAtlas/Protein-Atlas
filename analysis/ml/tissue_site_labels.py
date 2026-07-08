from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.adr_site_mapping import add_adr_site_mapping

DEFAULT_TISSUE_LABEL_SOURCES = [
    Path("data/external/adrecs_target/drug_target_adr_evidence.tsv"),
    Path("data/external/adrecs_target/target_adr_evidence.tsv"),
    Path("data/external/sider/drug_adr.tsv"),
    Path("data/external/ohdsi/omopReferenceSet.csv"),
    Path("data/external/ohdsi/ohdsiNegativeControls.csv"),
    Path("data/external/ohdsi/ohdsiDevelopmentNegativeControls.csv"),
    Path("data/external/faers/ohdsi_omop_faers_disproportionality.csv"),
]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _key(value: Any) -> str:
    return _clean(value).upper()


def _read(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def _existing_sources(run_dir: Path | None = None) -> list[Path]:
    roots = []
    if run_dir is not None:
        roots.append(run_dir)
    roots.extend([Path.cwd(), _repo_root()])
    found: list[Path] = []
    for rel in DEFAULT_TISSUE_LABEL_SOURCES:
        for root in roots:
            path = rel if rel.is_absolute() else root / rel
            if path.exists():
                found.append(path)
                break
    return found


def _first_nonempty(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    out = pd.Series("", index=df.index, dtype="object")
    for col in cols:
        if col not in df.columns:
            continue
        values = df[col].map(_clean)
        out = out.where(out.astype(str).str.len().gt(0), values)
    return out


def _ground_truth_state(frame: pd.DataFrame) -> pd.Series:
    for col in ["groundTruth", "ground_truth", "label_state", "label", "truth"]:
        if col in frame.columns:
            raw = frame[col]
            lowered = raw.fillna("").astype(str).str.strip().str.lower()
            out = pd.Series(pd.NA, index=frame.index, dtype="Float64")
            out = out.mask(lowered.isin({"1", "true", "positive", "pos"}), 1.0)
            out = out.mask(lowered.isin({"0", "false", "negative", "neg"}), 0.0)
            numeric = pd.to_numeric(raw, errors="coerce")
            out = out.where(out.notna(), numeric.where(numeric.isin([0, 1])))
            return out
    return pd.Series(pd.NA, index=frame.index, dtype="Float64")


def _source_name(path: Path) -> str:
    parts = [part.lower() for part in path.parts]
    for token in ["adrecs", "sider", "ohdsi", "omop", "faers", "dilirank", "toxrefdb", "cipa"]:
        if any(token in part for part in parts):
            return token
    return path.stem.lower()


def _evidence_from_file(path: Path) -> pd.DataFrame:
    raw = _read(path)
    if raw.empty:
        return pd.DataFrame()
    source = _source_name(path)
    work = raw.copy()
    work["drug_id"] = _first_nonempty(work, ["drug_id", "drug", "drug_name", "compound_id", "stitch_id", "rxnorm", "ingredient"])
    work["target_id"] = _first_nonempty(work, ["target_id", "target", "uniprot", "target_uniprot", "gene_symbol", "target_gene"])
    work["adr_term"] = _first_nonempty(work, ["adr_term", "adr", "event", "event_name", "outcome", "safety_event", "meddra_pt"])
    work = add_adr_site_mapping(work)
    label_state = pd.Series(pd.NA, index=work.index, dtype="Float64")
    confidence = pd.Series(0.65, index=work.index, dtype="Float64")
    evidence_type = pd.Series("independent_tissue_site_evidence", index=work.index, dtype="object")

    if source in {"adrecs", "sider"}:
        label_state = label_state.mask(work["adr_site_group"].notna(), 1.0)
        confidence = confidence.mask(label_state.eq(1), 0.85 if source == "adrecs" else 0.75)
        evidence_type = evidence_type.mask(label_state.eq(1), f"{source}_positive_drug_or_target_site")
    elif source in {"ohdsi", "omop"}:
        label_state = _ground_truth_state(work)
        confidence = confidence.mask(label_state.notna(), 0.80)
        evidence_type = evidence_type.mask(label_state.notna(), f"{source}_clinical_control_site")
    elif source == "faers":
        # Local FAERS table is generated as a non-signal control table. It is weak evidence only.
        label_state = label_state.mask(work["adr_site_group"].notna(), 0.0)
        confidence = confidence.mask(label_state.eq(0), 0.40)
        evidence_type = evidence_type.mask(label_state.eq(0), "faers_non_signal_weak_site_control")

    out = pd.DataFrame(
        {
            "drug_key": work["drug_id"].map(_key),
            "target_key": work["target_id"].map(_key),
            "adr_site_group": work["adr_site_group"].map(lambda value: _clean(value).lower().replace(" ", "_")),
            "tissue_site_label_evidence_state": label_state,
            "tissue_site_label_source": source,
            "tissue_site_label_evidence_type": evidence_type,
            "tissue_site_label_confidence": confidence,
        }
    )
    out = out[out["adr_site_group"].astype(str).str.len().gt(0)]
    out = out[out["tissue_site_label_evidence_state"].notna()]
    return out


def build_independent_tissue_site_evidence(run_dir: str | Path | None = None) -> pd.DataFrame:
    rows = []
    for path in _existing_sources(Path(run_dir) if run_dir else None):
        try:
            evidence = _evidence_from_file(path)
        except Exception:
            continue
        if not evidence.empty:
            evidence["tissue_site_label_source_path"] = str(path)
            rows.append(evidence)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True).drop_duplicates()


def _collapse(states: pd.Series) -> float | pd.NA:
    numeric = pd.to_numeric(states, errors="coerce").dropna()
    if numeric.empty:
        return pd.NA
    has_pos = numeric.eq(1).any()
    has_neg = numeric.eq(0).any()
    if has_pos and has_neg:
        return -1.0
    if has_pos:
        return 1.0
    if has_neg:
        return 0.0
    return pd.NA


def add_independent_tissue_site_labels(
    df: pd.DataFrame,
    *,
    run_dir: str | Path | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Project independent drug/site, target/site, and drug-target/site labels onto rows.

    This never derives labels from expression scores. It keeps ambiguous/conflicting
    projections as -1 so supervised tissue training can exclude them.
    """

    out = df.copy()
    evidence = build_independent_tissue_site_evidence(run_dir)
    if evidence.empty or "adr_site_group" not in out.columns:
        return out, {"status": "no_independent_tissue_label_evidence", "evidence_rows": int(len(evidence))}

    out["_drug_key_for_tissue_label"] = _first_nonempty(
        out, ["drug_id", "ligand_base", "generic_name", "display_name", "mapped_drug_name"]
    ).map(_key)
    out["_target_key_for_tissue_label"] = _first_nonempty(
        out, ["target_uniprot", "uniprot", "target_id", "target_gene", "gene_symbol"]
    ).map(_key)
    out["_site_key_for_tissue_label"] = out["adr_site_group"].map(lambda value: _clean(value).lower().replace(" ", "_"))

    projections: list[pd.DataFrame] = []
    key_specs = [
        ("drug_target_site", ["drug_key", "target_key", "adr_site_group"], ["_drug_key_for_tissue_label", "_target_key_for_tissue_label", "_site_key_for_tissue_label"]),
        ("target_site", ["target_key", "adr_site_group"], ["_target_key_for_tissue_label", "_site_key_for_tissue_label"]),
        ("drug_site", ["drug_key", "adr_site_group"], ["_drug_key_for_tissue_label", "_site_key_for_tissue_label"]),
    ]
    base = out.reset_index(names="_atlas_row_id")
    for projection_type, evidence_keys, row_keys in key_specs:
        ev = evidence[evidence[evidence_keys[0]].astype(str).str.len().gt(0)].copy()
        for key in evidence_keys[1:]:
            ev = ev[ev[key].astype(str).str.len().gt(0)]
        if ev.empty:
            continue
        ev = ev.rename(columns={e: r for e, r in zip(evidence_keys, row_keys, strict=False)})
        merged = base[["_atlas_row_id", *row_keys]].merge(ev, on=row_keys, how="inner")
        if merged.empty:
            continue
        merged["tissue_site_label_projection_type"] = projection_type
        projections.append(merged)

    if not projections:
        return out.drop(columns=[c for c in out.columns if c.startswith("_") and c.endswith("for_tissue_label")], errors="ignore"), {
            "status": "no_projected_tissue_labels",
            "evidence_rows": int(len(evidence)),
        }
    projected = pd.concat(projections, ignore_index=True)
    grouped = projected.groupby("_atlas_row_id", dropna=False).agg(
        tissue_site_label=("tissue_site_label_evidence_state", _collapse),
        tissue_site_label_source=("tissue_site_label_source", lambda s: ";".join(sorted(set(map(str, s))))),
        tissue_site_label_evidence_type=("tissue_site_label_evidence_type", lambda s: ";".join(sorted(set(map(str, s))))),
        tissue_site_label_confidence=("tissue_site_label_confidence", "max"),
        tissue_site_label_projection_type=("tissue_site_label_projection_type", lambda s: ";".join(sorted(set(map(str, s))))),
    )
    for col in grouped.columns:
        if col not in out.columns:
            out[col] = pd.NA
        out.loc[grouped.index, col] = grouped[col]
    out["tissue_site_label_policy"] = "independent_clinical_or_curated_site_label_not_expression_derived"
    out = out.drop(columns=[c for c in out.columns if c.startswith("_") and c.endswith("for_tissue_label")], errors="ignore")
    labels = pd.to_numeric(out.get("tissue_site_label", pd.Series(pd.NA, index=out.index)), errors="coerce")
    return out, {
        "status": "projected",
        "evidence_rows": int(len(evidence)),
        "projected_rows": int(labels.notna().sum()),
        "positive_rows": int(labels.eq(1).sum()),
        "negative_rows": int(labels.eq(0).sum()),
        "conflicting_rows": int(labels.eq(-1).sum()),
        "policy": "Tissue labels are projected from independent ADR/clinical site evidence, not from expression features.",
    }
