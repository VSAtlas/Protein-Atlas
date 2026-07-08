from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import pandas as pd

from analysis.graph.graph_scores import compute_pair_mechanism_scores
from analysis.graph.mechanism_graph import build_mechanism_graph
from analysis.external.papyrus_chembl import build_bioactivity_benchmark
from analysis.external.source_registry import source_registry_from_config, write_registry_manifest
from analysis.external.spd import build_spd_benchmark
from analysis.external.toxcast import build_toxcast_benchmark
from analysis.io import output_dirs


def _existing_path(value: Any) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(str(value))
    return path if path.exists() else None


def _inputs_mapping(config: Mapping[str, Any]) -> Mapping[str, Any]:
    raw = config.get("inputs")
    return raw if isinstance(raw, Mapping) else {}


def _mapping_path(config: Mapping[str, Any], spec: Mapping[str, Any]) -> str | Path | None:
    if spec.get("mapping"):
        return str(spec["mapping"])
    return _inputs_mapping(config).get("mapping_file")


def _source_path(config: Mapping[str, Any], spec: Mapping[str, Any], legacy_key: str) -> Path | None:
    path = _existing_path(spec.get("path"))
    if path is not None:
        return path
    return _existing_path(_inputs_mapping(config).get(legacy_key))


def _coverage_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for col in sorted(c for c in df.columns if c.endswith("_label_status")):
        prefix = col[: -len("_label_status")]
        counts = df[col].value_counts(dropna=False)
        rows.append(
            {
                "source": prefix,
                "status_column": col,
                "n_rows": len(df),
                "n_drugs": int(df["drug_id"].nunique()) if "drug_id" in df else 0,
                "n_targets": int(df["target_id"].nunique()) if "target_id" in df else 0,
                "n_drug_target_pairs": int(df[["drug_id", "target_id"]].drop_duplicates().shape[0])
                if {"drug_id", "target_id"}.issubset(df.columns)
                else 0,
                "status_counts": ";".join(f"{key}:{value}" for key, value in counts.items()),
            }
        )
    return pd.DataFrame(rows)


def _missing_coverage_row(source_name: str, status: str) -> dict[str, object]:
    return {
        "source": source_name,
        "status_column": f"{source_name}_label_status",
        "n_rows": 0,
        "n_drugs": 0,
        "n_targets": 0,
        "n_drug_target_pairs": 0,
        "status_counts": f"{status}:0",
    }


def _stage_mechanism_sources(
    pair_table_path: str | Path,
    config: Mapping[str, Any],
    registry: Mapping[str, Mapping[str, Any]],
    base_out: Path,
) -> list[dict[str, object]]:
    mechanism_specs = [
        ("sider", "sider_file"),
        ("ctd", "ctd_file"),
        ("opentargets_safety", "opentargets_safety_file"),
        ("reactome", "reactome_file"),
    ]
    rows: list[dict[str, object]] = []
    paths: dict[str, Path | None] = {}
    for source_name, legacy_key in mechanism_specs:
        spec = registry.get(source_name, {})
        if spec.get("enabled") is False:
            rows.append(
                {
                    "source": source_name,
                    "status": "disabled",
                    "path": "",
                    "exists": False,
                    "role": spec.get("role", ""),
                    "label_family": spec.get("label_family", ""),
                }
            )
            paths[source_name] = None
            continue
        source_path = _source_path(config, spec, legacy_key)
        paths[source_name] = source_path
        rows.append(
            {
                "source": source_name,
                "status": "loaded" if source_path is not None else "missing_local_file",
                "path": str(source_path or spec.get("path") or ""),
                "exists": source_path is not None,
                "role": spec.get("role", ""),
                "label_family": spec.get("label_family", ""),
            }
        )
    if any(path is not None for path in paths.values()):
        graph_dir = base_out / "mechanism_graph"
        graph_dir.mkdir(parents=True, exist_ok=True)
        edges_path = graph_dir / "mechanism_graph_edges.csv"
        nodes_path = graph_dir / "mechanism_graph_nodes.csv"
        scores_path = graph_dir / "pair_mechanism_scores.csv"
        build_mechanism_graph(
            pair_table_path,
            paths.get("sider"),
            paths.get("ctd"),
            paths.get("opentargets_safety"),
            paths.get("reactome"),
            None,
            edges_path,
            nodes_path,
        )
        compute_pair_mechanism_scores(pair_table_path, edges_path, scores_path)
        rows.append(
            {
                "source": "mechanism_graph",
                "status": "written",
                "path": str(scores_path),
                "exists": True,
                "role": "pair_level_mechanism_scoring",
                "label_family": "mechanism_evidence",
            }
        )
    return rows


def deduplicate_drug_target_labels(df: pd.DataFrame, out_path: str | Path) -> pd.DataFrame:
    if not {"drug_id", "target_id"}.issubset(df.columns):
        out = df.copy()
    else:
        aggregations: dict[str, str] = {}
        for col in df.columns:
            if col in {"drug_id", "target_id"}:
                continue
            if col.endswith(("_active", "_inactive", "_relevant", "_weak", "_unlikely", "_assayed")):
                aggregations[col] = "max"
            elif col.endswith(("_activity_nM", "_threshold_nM", "exposure_margin", "ac50_nM")):
                aggregations[col] = "min"
            else:
                aggregations[col] = "first"
        out = df.groupby(["drug_id", "target_id"], dropna=False).agg(aggregations).reset_index()
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    return out


def build_external_label_tables(
    pair_table_path: str | Path,
    config: Mapping[str, Any],
    out_dir: str | Path | None = None,
) -> pd.DataFrame:
    analysis_dir, _figures_dir, _models_dir = output_dirs(config)
    base_out = Path(out_dir or (analysis_dir / "external_labels"))
    base_out.mkdir(parents=True, exist_ok=True)
    registry = source_registry_from_config(config)
    write_registry_manifest(config, base_out / "source_registry_manifest.json")

    current_path = Path(pair_table_path)
    current_df = pd.read_csv(current_path)
    source_specs = [
        ("spd", "spd_file", "spd_atlas_labels.csv"),
        ("toxcast", "toxcast_file", "toxcast_atlas_labels.csv"),
        ("papyrus", "papyrus_file", "papyrus_atlas_labels.csv"),
        ("chembl", "chembl_file", "chembl_atlas_labels.csv"),
    ]
    stage_rows: list[dict[str, object]] = []
    missing_coverage_rows: list[dict[str, object]] = []
    for source_name, legacy_key, filename in source_specs:
        spec = registry.get(source_name, {})
        if spec.get("enabled") is False:
            stage_rows.append(
                {
                    "source": source_name,
                    "status": "disabled",
                    "path": "",
                    "exists": False,
                    "role": spec.get("role", ""),
                    "label_family": spec.get("label_family", ""),
                }
            )
            continue
        source_path = _source_path(config, spec, legacy_key)
        if source_path is None:
            stage_rows.append(
                {
                    "source": source_name,
                    "status": "missing_local_file",
                    "path": str(spec.get("path") or ""),
                    "exists": False,
                    "role": spec.get("role", ""),
                    "label_family": spec.get("label_family", ""),
                }
            )
            missing_coverage_rows.append(_missing_coverage_row(source_name, "missing_local_file"))
            continue
        out_path = base_out / filename
        mapping = _mapping_path(config, spec)
        stage_rows.append(
            {
                "source": source_name,
                "status": "loaded",
                "path": str(source_path),
                "exists": True,
                "role": spec.get("role", ""),
                "label_family": spec.get("label_family", ""),
                "output": str(out_path),
            }
        )
        if source_name == "spd":
            thresholds = config.get("exposure", {}).get("margin_thresholds", {}) if isinstance(config.get("exposure"), Mapping) else {}
            current_df = build_spd_benchmark(
                current_path,
                source_path,
                mapping,
                out_path,
                float(thresholds.get("strong", 10.0)),
                float(thresholds.get("weak", 100.0)),
            )
        elif source_name == "toxcast":
            threshold = float(spec.get("potent_threshold_nM", config.get("toxcast", {}).get("potent_threshold_nM", 10000.0) if isinstance(config.get("toxcast"), Mapping) else 10000.0))
            current_df = build_toxcast_benchmark(current_path, source_path, mapping, out_path, threshold)
        else:
            threshold = float(spec.get("active_threshold_nM", config.get("bioactivity", {}).get("active_threshold_nM", 10000.0) if isinstance(config.get("bioactivity"), Mapping) else 10000.0))
            current_df = build_bioactivity_benchmark(
                current_path,
                source_path,
                mapping,
                out_path,
                active_threshold_nM=threshold,
                label_prefix=source_name,
            )
        current_path = out_path

    mechanism_stage_rows = _stage_mechanism_sources(pair_table_path, config, registry, base_out)
    joined_path = base_out / "atlas_external_labels.csv"
    current_df.to_csv(joined_path, index=False)
    deduplicate_drug_target_labels(current_df, base_out / "atlas_external_labels_dedup_drug_target.csv")
    coverage = _coverage_summary(current_df)
    if missing_coverage_rows:
        coverage = pd.concat([coverage, pd.DataFrame(missing_coverage_rows)], ignore_index=True)
    coverage.to_csv(base_out / "source_coverage_summary.csv", index=False)
    pd.DataFrame(stage_rows).to_csv(base_out / "external_source_stage_summary.csv", index=False)
    pd.DataFrame(mechanism_stage_rows).to_csv(base_out / "mechanism_source_stage_summary.csv", index=False)
    if base_out.name == "external_labels":
        coverage.to_csv(base_out.parent / "source_coverage_summary.csv", index=False)
        pd.DataFrame(stage_rows).to_csv(base_out.parent / "external_source_stage_summary.csv", index=False)
        pd.DataFrame(mechanism_stage_rows).to_csv(base_out.parent / "mechanism_source_stage_summary.csv", index=False)
    return current_df
