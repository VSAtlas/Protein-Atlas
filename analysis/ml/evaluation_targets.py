from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.bigbind import build_bigbind_benchmark
from analysis.ml.labels import binary_label_series, truthy_series
from analysis.ml.mechanism_targets import build_mechanism_scores_from_edges, mechanism_label_status
from analysis.ml.spd_panel import stage_spd_full_panel


COMMON_ID_COLS = ["drug_id", "target_id", "pdb_id", "ligand_chemotype", "protein_class"]
PK_FEATURES = ["free_cmax_um", "cmax_um", "fraction_unbound_plasma"]


def _read_existing(path: str | Path | None) -> pd.DataFrame:
    if path is None or not Path(path).exists():
        return pd.DataFrame()
    return pd.read_csv(path, low_memory=False)


def _write_manifest(out_dir: Path, manifest: dict[str, Any]) -> None:
    (out_dir / "evaluation_target_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def _keep_columns(df: pd.DataFrame, extra: list[str]) -> pd.DataFrame:
    keep = [col for col in [*COMMON_ID_COLS, *extra] if col in df.columns]
    return df[keep].copy()


def _write_target(df: pd.DataFrame, out_path: Path, label_col: str) -> dict[str, Any]:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    labels = binary_label_series(df[label_col]) if label_col in df.columns else pd.Series(pd.NA, index=df.index)
    return {
        "path": str(out_path),
        "rows": int(len(df)),
        "label_col": label_col,
        "labelable_rows": int(labels.notna().sum()),
        "positive_rows": int(labels.eq(1).sum()),
        "negative_rows": int(labels.eq(0).sum()),
    }


def _bioactivity_target(df: pd.DataFrame, out_dir: Path) -> dict[str, Any]:
    if df.empty or "bioactivity_ml_label" not in df.columns:
        return {"status": "missing", "reason": "bioactivity_ml_label not available"}
    out = _keep_columns(
        df,
        [
            "bioactivity_ml_label",
            "label_source",
            "label_source_count",
            "bioactivity_ml_cross_source_disagreement",
            "chembl_label_available",
            "papyrus_label_available",
            "toxcast_label_available",
            "activity_publication_year",
            "database_release_year",
            "atlas_score",
            "consensus_score",
            *PK_FEATURES,
        ],
    )
    out["evaluation_target"] = "bioactivity_assay_activity"
    summary = _write_target(out, out_dir / "bioactivity_assay_target.csv", "bioactivity_ml_label")
    summary.update(
        {
            "status": "written",
            "primary_claim": "Atlas recovers measured in vitro activity where assay evidence exists.",
            "primary_feature_set": "pilot_binding_only",
            "pk_policy": "not primary for in-vitro assay activity; use only as sensitivity/provenance audit",
            "recommended_splits": ["drug_holdout", "target_holdout", "scaffold_holdout", "source_holdout"],
            "major_gap": "ChEMBL/Papyrus to ToxCast source transfer is currently near-random in the pilot.",
        }
    )
    return summary


def _exposure_target(df: pd.DataFrame, out_dir: Path) -> dict[str, Any]:
    if df.empty:
        return {"status": "missing", "reason": "SPD exposure table not available"}
    label_col = "ml_binary_label" if "ml_binary_label" in df.columns else "spd_exposure_relevant"
    if label_col not in df.columns:
        return {"status": "missing", "reason": "SPD binary exposure label not available"}
    out = _keep_columns(
        df,
        [
            label_col,
            "spd_label_status",
            "ml_supervised_eligible",
            "ml_exclude_reason",
            "exposure_margin",
            "ac50_nM",
            "free_cmax_nM",
            "total_cmax_nM",
            "atlas_score",
            "consensus_score",
            *PK_FEATURES,
        ],
    )
    out["evaluation_target"] = "spd_exposure_relevance"
    summary = _write_target(out, out_dir / "spd_exposure_relevance_target.csv", label_col)
    summary.update(
        {
            "status": "written",
            "primary_claim": "Atlas prioritizes exposure-relevant secondary pharmacology where SPD assays exist.",
            "primary_feature_set": "pilot_binding_only",
            "forbidden_predictive_features": ["exposure_margin", "ac50_nM", "free_cmax_nM", "total_cmax_nM"],
            "pk_policy": "PK/free-Cmax defines the SPD label; use PK models only as exposure-context sensitivity, not as noncircular prediction.",
            "recommended_splits": ["drug_holdout", "target_holdout"],
            "major_gap": "Pilot SPD target is small until the full SPD protein panel is run through Atlas.",
        }
    )
    return summary


def _mechanism_target(pair_df: pd.DataFrame, mechanism_df: pd.DataFrame, out_dir: Path) -> dict[str, Any]:
    source = mechanism_df.copy() if not mechanism_df.empty else pair_df.copy()
    if source.empty:
        return {"status": "missing", "reason": "no pair or mechanism score table available"}
    out = _keep_columns(
        source,
        [
            "literature_supported_label",
            "drug_adr_known",
            "target_adr_known",
            "target_pathway_adr_link",
            "drug_target_known",
            "triad_complete",
            "mechanism_path_count",
            "mechanism_graph_score",
            "target_adr_evidence",
            "pathway_evidence",
            "atlas_score",
            "consensus_score",
            *PK_FEATURES,
        ],
    )
    label = pd.Series(pd.NA, index=out.index, dtype="Int64")
    if "literature_supported_label" in out.columns:
        parsed = binary_label_series(out["literature_supported_label"])
        label = label.mask(parsed.notna(), parsed)
    if "triad_complete" in out.columns:
        triad = truthy_series(out["triad_complete"])
        label = label.mask(triad, 1)
    elif "mechanism_graph_score" in out.columns:
        score = pd.to_numeric(out["mechanism_graph_score"], errors="coerce")
        label = label.mask(score.gt(0), 1)
    out["mechanism_ml_label"] = label
    out["mechanism_label_status"] = mechanism_label_status(out)
    out["evaluation_target"] = "adr_mechanism_support"
    summary = _write_target(out, out_dir / "mechanism_adr_support_target.csv", "mechanism_ml_label")
    summary.update(
        {
            "status": "written",
            "primary_claim": "Atlas enriches for plausible ADR mechanism support rather than direct clinical ADR prediction.",
            "primary_feature_set": "pilot_binding_only plus nonlabel biological context only in ablations",
            "forbidden_predictive_features": [
                "literature_supported_label",
                "mechanism_ml_label",
                "drug_adr_known",
                "target_adr_known",
                "target_pathway_adr_link",
                "triad_complete",
            ],
            "pk_policy": "exposure can stratify mechanism plausibility but is not sufficient for mechanism support.",
            "recommended_splits": ["drug_holdout", "target_holdout", "adr_panel_holdout"],
            "major_gap": "Requires local SIDER/OpenTargets/CTD/Reactome graph coverage; missing mechanism evidence remains unknown.",
        }
    )
    return summary


def build_evaluation_targets(
    out_dir: str | Path,
    *,
    pair_table: str | Path | None = None,
    bioactivity_table: str | Path | None = None,
    spd_table: str | Path | None = None,
    spd_full_panel: str | Path | None = None,
    spd_mapping: str | Path | None = None,
    mechanism_scores: str | Path | None = None,
    mechanism_edges: str | Path | None = None,
    bigbind_table: str | Path | None = None,
    bigbind_mapping: str | Path | None = None,
) -> dict[str, Any]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    pair_df = _read_existing(pair_table)
    bioactivity_df = _read_existing(bioactivity_table)
    spd_df = _read_existing(spd_table)
    mechanism_path = Path(mechanism_scores) if mechanism_scores is not None else None
    if (
        (mechanism_path is None or not mechanism_path.exists())
        and mechanism_edges is not None
        and pair_table is not None
        and Path(mechanism_edges).exists()
    ):
        mechanism_path = out / "pair_mechanism_scores.csv"
        build_mechanism_scores_from_edges(pair_table, mechanism_edges, mechanism_path)
    mechanism_df = _read_existing(mechanism_path)
    manifest: dict[str, Any] = {
        "purpose": "Separate Atlas evaluation targets so bioactivity, exposure relevance, and ADR mechanism support are not mixed.",
        "bigbind_note": (
            "BigBind is appropriate for docking/bioactivity benchmarking. Existing runtime BigBind utilities live in "
            "src/ml/data/bigbind.py; the Atlas importer is only a local pair-table join shim."
        ),
        "targets": {
            "bioactivity": _bioactivity_target(bioactivity_df, out),
            "exposure_relevance": _exposure_target(spd_df, out),
            "mechanism_adr_support": _mechanism_target(pair_df, mechanism_df, out),
        },
        "remaining_publication_gaps": [
            "Run full SPD target/protein panel through Atlas before making exposure-relevance ML claims.",
            "Stage local SIDER/OpenTargets Safety/CTD/Reactome graph files for mechanism labels and ADR panels.",
            "Do not use PK/free-Cmax as a primary predictor for generic in-vitro bioactivity.",
            "Treat BigBind as docking/bioactivity benchmark support, not exposure or ADR validation.",
            "Report source-holdout results; current ChEMBL/Papyrus to ToxCast transfer is not a positive claim.",
        ],
    }
    if bigbind_table is not None and pair_table is not None and Path(bigbind_table).exists():
        bigbind_out = out / "bigbind_bioactivity_target.csv"
        build_bigbind_benchmark(pair_table, bigbind_table, bigbind_mapping, bigbind_out)
        manifest["targets"]["bigbind_bioactivity"] = {
            "status": "written",
            "path": str(bigbind_out),
            "primary_claim": "Atlas docking scores recover BigBind-style active/inactive benchmark labels.",
            "primary_feature_set": "binding/docking features only",
            "pk_policy": "not relevant to BigBind in-vitro docking benchmark labels",
        }
    if spd_full_panel is not None and pair_table is not None and Path(spd_full_panel).exists():
        manifest["targets"]["spd_full_panel_coverage"] = stage_spd_full_panel(
            spd_full_panel,
            pair_table,
            out / "spd_full_panel",
            mapping_path=spd_mapping,
        )
    _write_manifest(out, manifest)
    return manifest
