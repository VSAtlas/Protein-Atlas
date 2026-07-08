from __future__ import annotations

import importlib
import json
import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pandas.api import types as ptypes


DATASET_SUFFIXES = (".csv", ".csv.gz", ".tsv", ".tsv.gz", ".parquet", ".jsonl", ".json")
DEFAULT_TOOLS = ("profile", "phik", "dython", "mutual-info", "network")
PROTECTED_PATH_PARTS = {
    "ci",
    "docked",
    "input_pdbs",
    "logs",
    "post_docked",
    "prepped_ligands",
    "processed_pdbs",
}
ID_LIKE_NAMES = {
    "drug_id",
    "target_id",
    "pdb_id",
    "ligand_id",
    "compound_id",
    "uniprot_id",
    "smiles",
    "inchi",
    "inchikey",
}


@dataclass(frozen=True)
class DatasetResolution:
    path: Path
    matches: tuple[Path, ...] = ()


def resolve_dataset_path(
    dataset: str | Path | None,
    *,
    repo_root: str | Path = ".",
    dataset_name: str | None = None,
    max_matches: int = 25,
) -> DatasetResolution:
    """Resolve an explicit dataset path or search by a filename/path substring."""

    root = Path(repo_root).resolve()
    query = dataset_name or ""
    if dataset:
        raw = Path(dataset).expanduser()
        candidate = raw if raw.is_absolute() else root / raw
        if candidate.exists():
            _reject_protected_path(candidate)
            return DatasetResolution(path=candidate.resolve())
        if not raw.suffix and len(raw.parts) == 1:
            query = str(dataset)
        else:
            raise FileNotFoundError(f"dataset not found: {candidate}")

    query = query.strip()
    if not query:
        raise ValueError("provide --dataset PATH or --dataset-name TEXT")

    matches = _search_dataset_candidates(root, query, max_matches=max_matches)
    if not matches:
        raise FileNotFoundError(f"no dataset files matched query: {query!r}")
    if len(matches) > 1:
        rendered = "\n".join(f"  - {path}" for path in matches[:max_matches])
        raise ValueError(
            f"dataset query {query!r} matched multiple files; pass --dataset explicitly:\n{rendered}"
        )
    return DatasetResolution(path=matches[0].resolve(), matches=tuple(matches))


def find_dataset_candidates(
    query: str,
    *,
    repo_root: str | Path = ".",
    max_matches: int = 25,
) -> list[Path]:
    return _search_dataset_candidates(Path(repo_root).resolve(), query, max_matches=max_matches)


def run_dataset_eda(
    dataset: str | Path | None,
    out_dir: str | Path,
    *,
    repo_root: str | Path = ".",
    dataset_name: str | None = None,
    tools: Sequence[str] | None = None,
    label_col: str | None = None,
    read_rows: int | None = None,
    sample_rows: int = 10_000,
    random_state: int = 13,
    max_association_columns: int = 80,
    max_category_levels: int = 200,
    max_mi_features: int = 300,
    network_threshold: float = 0.35,
    network_max_edges: int = 300,
    ydata_minimal: bool = True,
    dython_nominal_assoc: str = "cramer",
    include_id_like: bool = False,
    include_high_cardinality: bool = False,
    fail_on_missing: bool = False,
) -> dict[str, Any]:
    selected_tools = _normalize_tools(tools)
    root = Path(repo_root).resolve()
    resolved = resolve_dataset_path(dataset, repo_root=root, dataset_name=dataset_name)
    out = Path(out_dir)
    if not out.is_absolute():
        out = root / out
    out.mkdir(parents=True, exist_ok=True)

    frame = _read_dataset(resolved.path, read_rows=read_rows)
    sampled = _sample_frame(frame, sample_rows=sample_rows, random_state=random_state)
    association_frame, association_skips = _select_association_columns(
        sampled,
        label_col=label_col,
        max_columns=max_association_columns,
        max_category_levels=max_category_levels,
        include_high_cardinality=include_high_cardinality,
    )

    outputs: dict[str, str] = {}
    skipped: dict[str, str] = {}
    errors: dict[str, str] = {}
    parameters: dict[str, Any] = {
        "tools": selected_tools,
        "label_col": label_col,
        "read_rows": read_rows,
        "sample_rows": sample_rows,
        "sampled_rows": int(len(sampled)),
        "random_state": random_state,
        "max_association_columns": max_association_columns,
        "max_category_levels": max_category_levels,
        "max_mi_features": max_mi_features,
        "network_threshold": network_threshold,
        "network_max_edges": network_max_edges,
        "ydata_minimal": ydata_minimal,
        "dython_nominal_assoc": dython_nominal_assoc,
        "include_id_like": include_id_like,
        "include_high_cardinality": include_high_cardinality,
    }

    outputs["dataset_shape"] = _write_json(
        out / "dataset_shape.json",
        {
            "dataset": str(resolved.path),
            "rows_loaded": int(len(frame)),
            "columns_loaded": int(frame.shape[1]),
            "rows_analyzed": int(len(sampled)),
            "sampled": bool(len(sampled) != len(frame)),
        },
    )
    outputs["column_summary"] = _write_column_summary(frame, out / "column_summary.csv")
    outputs["association_columns"] = _write_association_column_report(
        association_frame,
        association_skips,
        out,
    )
    if label_col:
        if label_col in frame.columns:
            outputs["label_summary"] = _write_label_summary(frame, label_col, out / "label_summary.csv")
        else:
            skipped["label_summary"] = f"label column not found: {label_col}"

    matrices: dict[str, pd.DataFrame] = {}
    if "profile" in selected_tools:
        _run_step(
            "profile",
            outputs,
            skipped,
            errors,
            fail_on_missing,
            lambda: _run_ydata_profile(sampled, out, minimal=ydata_minimal),
        )
    if "phik" in selected_tools:
        _run_step(
            "phik",
            outputs,
            skipped,
            errors,
            fail_on_missing,
            lambda: _run_phik(association_frame, out, matrices),
        )
    if "dython" in selected_tools:
        _run_step(
            "dython",
            outputs,
            skipped,
            errors,
            fail_on_missing,
            lambda: _run_dython(
                association_frame,
                out,
                matrices,
                nominal_assoc=dython_nominal_assoc,
            ),
        )
    if "mutual-info" in selected_tools:
        if label_col:
            _run_step(
                "mutual-info",
                outputs,
                skipped,
                errors,
                fail_on_missing,
                lambda: _run_mutual_info(
                    sampled,
                    label_col,
                    out,
                    max_features=max_mi_features,
                    max_category_levels=max_category_levels,
                    random_state=random_state,
                    include_id_like=include_id_like,
                    include_high_cardinality=include_high_cardinality,
                ),
            )
        else:
            skipped["mutual-info"] = "requires --label-col"
    if "network" in selected_tools:
        _run_step(
            "network",
            outputs,
            skipped,
            errors,
            fail_on_missing,
            lambda: _run_network(
                association_frame,
                out,
                matrices,
                threshold=network_threshold,
                max_edges=network_max_edges,
            ),
        )

    manifest: dict[str, Any] = {
        "dataset": str(resolved.path),
        "out_dir": str(out),
        "shape": {
            "rows_loaded": int(len(frame)),
            "columns_loaded": int(frame.shape[1]),
            "rows_analyzed": int(len(sampled)),
            "columns_analyzed_for_association": int(association_frame.shape[1]),
        },
        "parameters": parameters,
        "outputs": outputs,
        "skipped": skipped,
        "errors": errors,
    }
    manifest_path = out / "dataset_eda_manifest.json"
    outputs["manifest"] = str(manifest_path)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    _write_index_html(out, manifest)
    return manifest


def _search_dataset_candidates(root: Path, query: str, *, max_matches: int) -> list[Path]:
    direct_matches = _direct_dataset_matches(root, query, max_matches=max_matches)
    if direct_matches:
        return direct_matches

    query_lower = query.lower()
    roots = [root / "data", root / "outputs" / "data", root / "analysis", root / "outputs" / "analysis"]
    matches: list[Path] = []
    for search_root in roots:
        if not search_root.exists():
            continue
        for path in _iter_dataset_files(search_root):
            text = str(path.relative_to(root) if path.is_relative_to(root) else path).lower()
            if query_lower in text:
                matches.append(path)
                if len(matches) >= max_matches:
                    return sorted(matches, key=lambda item: (len(str(item)), str(item)))
    return sorted(matches, key=lambda item: (len(str(item)), str(item)))



def _direct_dataset_matches(root: Path, query: str, *, max_matches: int) -> list[Path]:
    raw = Path(query).expanduser()
    candidates = [raw if raw.is_absolute() else root / raw]
    if not raw.is_absolute():
        candidates.extend(
            [
                root / "data" / raw,
                root / "outputs" / "data" / raw,
                root / "analysis" / raw,
                root / "outputs" / "analysis" / raw,
            ]
        )
    for candidate in candidates:
        if not candidate.exists():
            continue
        _reject_protected_path(candidate)
        if candidate.is_file() and _dataset_suffix(candidate):
            return [candidate.resolve()]
        if candidate.is_dir():
            matches: list[Path] = []
            for path in _iter_dataset_files(candidate):
                matches.append(path.resolve())
                if len(matches) >= max_matches:
                    break
            if matches:
                return sorted(matches, key=lambda item: (len(str(item)), str(item)))
    return []


def _iter_dataset_files(root: Path) -> Any:
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if PROTECTED_PATH_PARTS.intersection(path.parts):
            continue
        if _dataset_suffix(path):
            yield path


def _reject_protected_path(path: Path) -> None:
    if PROTECTED_PATH_PARTS.intersection(path.parts):
        raise ValueError(f"refusing to read protected Atlas path: {path}")


def _dataset_suffix(path: Path) -> str:
    suffixes = "".join(path.suffixes[-2:]).lower()
    if suffixes in DATASET_SUFFIXES:
        return suffixes
    suffix = path.suffix.lower()
    return suffix if suffix in DATASET_SUFFIXES else ""


def _normalize_tools(tools: Sequence[str] | None) -> list[str]:
    if not tools:
        return list(DEFAULT_TOOLS)
    normalized: list[str] = []
    aliases = {
        "all": list(DEFAULT_TOOLS),
        "ydata": ["profile"],
        "ydata-profiling": ["profile"],
        "mi": ["mutual-info"],
        "mutual_info": ["mutual-info"],
        "mutual-information": ["mutual-info"],
        "pyvis": ["network"],
        "networkx": ["network"],
    }
    for raw in tools:
        for item in str(raw).split(","):
            name = item.strip().lower()
            if not name:
                continue
            expanded = aliases.get(name, [name])
            for value in expanded:
                if value not in DEFAULT_TOOLS:
                    raise ValueError(f"unknown EDA tool: {value}")
                if value not in normalized:
                    normalized.append(value)
    return normalized


def _read_dataset(path: Path, *, read_rows: int | None) -> pd.DataFrame:
    suffix = _dataset_suffix(path)
    if suffix in {".csv", ".csv.gz"}:
        return pd.read_csv(path, low_memory=False, nrows=read_rows)
    if suffix in {".tsv", ".tsv.gz"}:
        return pd.read_csv(path, sep="\t", low_memory=False, nrows=read_rows)
    if suffix == ".parquet":
        frame = pd.read_parquet(path)
        return frame.head(read_rows) if read_rows else frame
    if suffix == ".jsonl":
        return pd.read_json(path, lines=True, nrows=read_rows)
    if suffix == ".json":
        frame = pd.read_json(path)
        return frame.head(read_rows) if read_rows else frame
    raise ValueError(f"unsupported dataset suffix for {path}")


def _sample_frame(frame: pd.DataFrame, *, sample_rows: int, random_state: int) -> pd.DataFrame:
    if sample_rows <= 0 or len(frame) <= sample_rows:
        return frame.copy()
    return frame.sample(n=sample_rows, random_state=random_state).reset_index(drop=True)


def _write_column_summary(frame: pd.DataFrame, path: Path) -> str:
    rows: list[dict[str, Any]] = []
    total = max(len(frame), 1)
    for col in frame.columns:
        series = frame[col]
        non_null = int(series.notna().sum())
        row: dict[str, Any] = {
            "column": str(col),
            "dtype": str(series.dtype),
            "non_null": non_null,
            "missing": int(series.isna().sum()),
            "missing_fraction": float(series.isna().sum() / total),
            "n_unique": int(series.nunique(dropna=True)),
            "is_numeric": bool(ptypes.is_numeric_dtype(series)),
            "is_bool": bool(ptypes.is_bool_dtype(series)),
            "is_datetime": bool(ptypes.is_datetime64_any_dtype(series)),
        }
        if ptypes.is_numeric_dtype(series):
            numeric = pd.to_numeric(series, errors="coerce")
            row.update(
                {
                    "min": _safe_float(numeric.min()),
                    "max": _safe_float(numeric.max()),
                    "mean": _safe_float(numeric.mean()),
                    "std": _safe_float(numeric.std()),
                }
            )
        else:
            counts = series.dropna().astype(str).value_counts().head(1)
            if not counts.empty:
                row["top_value"] = str(counts.index[0])[:200]
                row["top_count"] = int(counts.iloc[0])
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


def _write_label_summary(frame: pd.DataFrame, label_col: str, path: Path) -> str:
    counts = frame[label_col].fillna("<NA>").astype(str).value_counts(dropna=False)
    summary = counts.rename_axis("label").reset_index(name="count")
    summary["fraction"] = summary["count"] / max(int(summary["count"].sum()), 1)
    summary.to_csv(path, index=False)
    return str(path)


def _select_association_columns(
    frame: pd.DataFrame,
    *,
    label_col: str | None,
    max_columns: int,
    max_category_levels: int,
    include_high_cardinality: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected: list[str] = []
    skipped: list[dict[str, Any]] = []
    scored: list[tuple[tuple[int, float, int, str], str]] = []
    for col in frame.columns:
        series = frame[col]
        n_unique = int(series.nunique(dropna=True))
        missing_fraction = float(series.isna().mean())
        if n_unique <= 1:
            skipped.append({"column": col, "reason": "constant_or_empty", "n_unique": n_unique})
            continue
        if not include_high_cardinality and not ptypes.is_numeric_dtype(series) and n_unique > max_category_levels:
            skipped.append({"column": col, "reason": "high_cardinality", "n_unique": n_unique})
            continue
        if label_col and col == label_col:
            scored.append(((0, missing_fraction, n_unique, str(col)), col))
        elif ptypes.is_numeric_dtype(series):
            scored.append(((1, missing_fraction, n_unique, str(col)), col))
        else:
            scored.append(((2, missing_fraction, n_unique, str(col)), col))
    for _, col in sorted(scored, key=lambda item: item[0]):
        if len(selected) >= max_columns:
            skipped.append(
                {
                    "column": col,
                    "reason": "over_max_association_columns",
                    "n_unique": int(frame[col].nunique(dropna=True)),
                }
            )
            continue
        selected.append(col)
    return frame[selected].copy(), pd.DataFrame(skipped)


def _write_association_column_report(
    association_frame: pd.DataFrame,
    skipped: pd.DataFrame,
    out: Path,
) -> str:
    selected_path = out / "association_input_columns.txt"
    selected_path.write_text("\n".join(map(str, association_frame.columns)) + "\n", encoding="utf-8")
    skipped_path = out / "association_skipped_columns.csv"
    skipped.to_csv(skipped_path, index=False)
    return str(selected_path)


def _run_step(
    name: str,
    outputs: dict[str, str],
    skipped: dict[str, str],
    errors: dict[str, str],
    fail_on_missing: bool,
    callback: Any,
) -> None:
    try:
        produced = callback()
        if produced:
            outputs.update(produced)
    except OptionalDependencyMissing as exc:
        skipped[name] = str(exc)
        if fail_on_missing:
            raise
    except Exception as exc:
        errors[name] = f"{type(exc).__name__}: {exc}"
        if fail_on_missing:
            raise


class OptionalDependencyMissing(RuntimeError):
    pass


def _optional_import(module: str, package: str | None = None) -> Any:
    try:
        return importlib.import_module(module)
    except Exception as exc:  # pragma: no cover - depends on local optional env
        raise OptionalDependencyMissing(f"install optional package {package or module}: {exc}") from exc


def _run_ydata_profile(frame: pd.DataFrame, out: Path, *, minimal: bool) -> dict[str, str]:
    ydata = _optional_import("ydata_profiling", "ydata-profiling")
    profile_cls = getattr(ydata, "ProfileReport")
    profile = profile_cls(
        frame,
        title="Atlas dataset EDA profile",
        minimal=minimal,
        explorative=not minimal,
    )
    html_path = out / "ydata_profile.html"
    json_path = out / "ydata_profile.json"
    profile.to_file(str(html_path))
    json_text = profile.to_json()
    json_path.write_text(json_text, encoding="utf-8")
    return {"ydata_profile_html": str(html_path), "ydata_profile_json": str(json_path)}


def _run_phik(
    frame: pd.DataFrame,
    out: Path,
    matrices: dict[str, pd.DataFrame],
) -> dict[str, str]:
    if frame.shape[1] < 2:
        raise ValueError("need at least two association columns for Phi_K")
    _optional_import("phik", "phik")
    interval_cols = [str(col) for col in frame.columns if ptypes.is_numeric_dtype(frame[col])]
    matrix = frame.phik_matrix(interval_cols=interval_cols or None)  # type: ignore[attr-defined]
    matrix = _clean_matrix(matrix)
    matrices["phik"] = matrix
    csv_path = out / "phik_matrix.csv"
    json_path = out / "phik_matrix.json"
    matrix.to_csv(csv_path)
    json_path.write_text(matrix.to_json(orient="index", indent=2), encoding="utf-8")
    return {"phik_matrix": str(csv_path), "phik_matrix_json": str(json_path)}


def _run_dython(
    frame: pd.DataFrame,
    out: Path,
    matrices: dict[str, pd.DataFrame],
    *,
    nominal_assoc: str,
) -> dict[str, str]:
    if frame.shape[1] < 2:
        raise ValueError("need at least two association columns for dython")
    nominal = _optional_import("dython.nominal", "dython")
    result = nominal.associations(
        frame,
        nominal_columns="auto",
        nom_nom_assoc=nominal_assoc,
        num_num_assoc="pearson",
        compute_only=True,
        plot=False,
    )
    matrix = _extract_dython_matrix(result)
    matrix = _clean_matrix(matrix)
    matrices["dython"] = matrix
    csv_path = out / "dython_associations.csv"
    json_path = out / "dython_associations.json"
    matrix.to_csv(csv_path)
    json_path.write_text(matrix.to_json(orient="index", indent=2), encoding="utf-8")
    return {"dython_associations": str(csv_path), "dython_associations_json": str(json_path)}


def _extract_dython_matrix(result: Any) -> pd.DataFrame:
    if isinstance(result, pd.DataFrame):
        return result
    if isinstance(result, dict):
        for key in ("corr", "correlation", "associations"):
            value = result.get(key)
            if isinstance(value, pd.DataFrame):
                return value
        for value in result.values():
            if isinstance(value, pd.DataFrame):
                return value
    if isinstance(result, tuple):
        for value in result:
            if isinstance(value, pd.DataFrame):
                return value
    return pd.DataFrame(result)


def _run_mutual_info(
    frame: pd.DataFrame,
    label_col: str | None,
    out: Path,
    *,
    max_features: int,
    max_category_levels: int,
    random_state: int,
    include_id_like: bool,
    include_high_cardinality: bool,
) -> dict[str, str]:
    if not label_col:
        raise ValueError("mutual information requires --label-col")
    if label_col not in frame.columns:
        raise ValueError(f"label column not found: {label_col}")
    sklearn_feature_selection = _optional_import("sklearn.feature_selection", "scikit-learn")
    y_raw = frame[label_col]
    keep_rows = y_raw.notna()
    y_values = y_raw[keep_rows].astype(str)
    y_encoded, labels = pd.factorize(y_values, sort=True)
    if len(labels) < 2:
        raise ValueError(f"label column needs at least two classes: {label_col}")

    feature_cols, excluded = _select_mutual_info_features(
        frame.loc[keep_rows],
        label_col=label_col,
        max_features=max_features,
        max_category_levels=max_category_levels,
        include_id_like=include_id_like,
        include_high_cardinality=include_high_cardinality,
    )
    if not feature_cols:
        raise ValueError("no eligible features for mutual information")

    x_values, discrete_mask = _encode_mi_features(frame.loc[keep_rows, feature_cols])
    scores = sklearn_feature_selection.mutual_info_classif(
        x_values,
        y_encoded,
        discrete_features=discrete_mask,
        random_state=random_state,
    )
    rows: list[dict[str, Any]] = []
    for col, score, discrete in zip(feature_cols, scores, discrete_mask):
        series = frame[col]
        rows.append(
            {
                "feature": col,
                "mutual_information": float(score),
                "discrete": bool(discrete),
                "dtype": str(series.dtype),
                "n_unique": int(series.nunique(dropna=True)),
                "missing_fraction": float(series.isna().mean()),
            }
        )
    ranking = pd.DataFrame(rows).sort_values("mutual_information", ascending=False)
    ranking_path = out / "mutual_information.csv"
    excluded_path = out / "mutual_information_excluded_features.csv"
    ranking.to_csv(ranking_path, index=False)
    pd.DataFrame(excluded).to_csv(excluded_path, index=False)
    return {
        "mutual_information": str(ranking_path),
        "mutual_information_excluded_features": str(excluded_path),
    }


def _select_mutual_info_features(
    frame: pd.DataFrame,
    *,
    label_col: str,
    max_features: int,
    max_category_levels: int,
    include_id_like: bool,
    include_high_cardinality: bool,
) -> tuple[list[str], list[dict[str, Any]]]:
    scored: list[tuple[tuple[int, float, int, str], str]] = []
    excluded: list[dict[str, Any]] = []
    for col in frame.columns:
        if col == label_col:
            continue
        series = frame[col]
        n_unique = int(series.nunique(dropna=True))
        if n_unique <= 1:
            excluded.append({"feature": col, "reason": "constant_or_empty", "n_unique": n_unique})
            continue
        if not include_id_like and _is_id_like(str(col)):
            excluded.append({"feature": col, "reason": "id_like_metadata", "n_unique": n_unique})
            continue
        if not include_high_cardinality and not ptypes.is_numeric_dtype(series) and n_unique > max_category_levels:
            excluded.append({"feature": col, "reason": "high_cardinality", "n_unique": n_unique})
            continue
        type_rank = 0 if ptypes.is_numeric_dtype(series) else 1
        scored.append(((type_rank, float(series.isna().mean()), n_unique, str(col)), col))
    selected = [col for _, col in sorted(scored, key=lambda item: item[0])[:max_features]]
    for _, col in sorted(scored, key=lambda item: item[0])[max_features:]:
        excluded.append(
            {
                "feature": col,
                "reason": "over_max_mi_features",
                "n_unique": int(frame[col].nunique(dropna=True)),
            }
        )
    return selected, excluded


def _encode_mi_features(frame: pd.DataFrame) -> tuple[np.ndarray, list[bool]]:
    arrays: list[np.ndarray] = []
    discrete_mask: list[bool] = []
    for col in frame.columns:
        series = frame[col]
        if ptypes.is_numeric_dtype(series):
            numeric = pd.to_numeric(series, errors="coerce")
            fill_value = numeric.median()
            if pd.isna(fill_value):
                fill_value = 0.0
            arrays.append(numeric.fillna(float(fill_value)).to_numpy(dtype=float))
            discrete_mask.append(False)
        else:
            encoded, _ = pd.factorize(series.fillna("<NA>").astype(str), sort=True)
            arrays.append(encoded.astype(float))
            discrete_mask.append(True)
    return np.column_stack(arrays), discrete_mask


def _run_network(
    frame: pd.DataFrame,
    out: Path,
    matrices: dict[str, pd.DataFrame],
    *,
    threshold: float,
    max_edges: int,
) -> dict[str, str]:
    networkx = _optional_import("networkx", "networkx")
    source, matrix = _choose_network_matrix(frame, matrices)
    graph = networkx.Graph(source_matrix=source, threshold=threshold)
    for col in matrix.columns:
        graph.add_node(str(col))

    edge_rows: list[dict[str, Any]] = []
    cols = list(matrix.columns)
    for i, left in enumerate(cols):
        for right in cols[i + 1 :]:
            weight = _safe_float(matrix.loc[left, right])
            if weight is None:
                continue
            abs_weight = abs(weight)
            if abs_weight < threshold:
                continue
            edge_rows.append(
                {
                    "source": str(left),
                    "target": str(right),
                    "weight": float(weight),
                    "abs_weight": float(abs_weight),
                }
            )
    edge_rows = sorted(edge_rows, key=lambda item: item["abs_weight"], reverse=True)[:max_edges]
    for row in edge_rows:
        graph.add_edge(
            row["source"],
            row["target"],
            weight=row["weight"],
            abs_weight=row["abs_weight"],
        )

    nodes_path = out / "association_network_nodes.csv"
    edges_path = out / "association_network_edges.csv"
    graphml_path = out / "association_network.graphml"
    summary_path = out / "association_network_summary.json"
    pd.DataFrame({"node": [str(node) for node in graph.nodes]}).to_csv(nodes_path, index=False)
    pd.DataFrame(edge_rows).to_csv(edges_path, index=False)
    networkx.write_graphml(graph, graphml_path)
    summary = {
        "source_matrix": source,
        "threshold": threshold,
        "max_edges": max_edges,
        "nodes": int(graph.number_of_nodes()),
        "edges": int(graph.number_of_edges()),
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    outputs = {
        "association_network_nodes": str(nodes_path),
        "association_network_edges": str(edges_path),
        "association_network_graphml": str(graphml_path),
        "association_network_summary": str(summary_path),
    }
    try:
        outputs["association_network_html"] = _write_pyvis_network(graph, out / "association_network.html")
    except OptionalDependencyMissing:
        pass
    return outputs


def _choose_network_matrix(
    frame: pd.DataFrame,
    matrices: dict[str, pd.DataFrame],
) -> tuple[str, pd.DataFrame]:
    for name in ("phik", "dython"):
        matrix = matrices.get(name)
        if matrix is not None and matrix.shape[0] >= 2:
            return name, matrix
    numeric = frame.select_dtypes(include=[np.number])
    if numeric.shape[1] < 2:
        raise ValueError("network requires a Phi_K/dython matrix or at least two numeric columns")
    return "pearson_numeric_fallback", _clean_matrix(numeric.corr())


def _write_pyvis_network(graph: Any, path: Path) -> str:
    pyvis_network = _optional_import("pyvis.network", "pyvis")
    net = pyvis_network.Network(
        height="820px",
        width="100%",
        bgcolor="#ffffff",
        font_color="#111827",
        cdn_resources="in_line",
    )
    net.barnes_hut(gravity=-35_000, central_gravity=0.25, spring_length=150)
    for node in graph.nodes:
        degree = graph.degree[node]
        net.add_node(str(node), label=str(node), title=f"{node}<br>degree: {degree}", value=max(degree, 1))
    for left, right, data in graph.edges(data=True):
        weight = float(data.get("weight", 0.0))
        abs_weight = float(data.get("abs_weight", abs(weight)))
        net.add_edge(
            str(left),
            str(right),
            value=max(abs_weight * 10.0, 1.0),
            title=f"association: {weight:.3f}",
        )
    net.write_html(str(path), notebook=False)
    return str(path)


def _clean_matrix(matrix: pd.DataFrame) -> pd.DataFrame:
    matrix = pd.DataFrame(matrix).copy()
    matrix.index = matrix.index.map(str)
    matrix.columns = matrix.columns.map(str)
    return matrix.replace([np.inf, -np.inf], np.nan)


def _safe_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(result) or math.isinf(result):
        return None
    return result


def _is_id_like(name: str) -> bool:
    lowered = name.lower()
    return lowered in ID_LIKE_NAMES or lowered.endswith("_id") or lowered.endswith("_name")


def _write_json(path: Path, payload: dict[str, Any]) -> str:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return str(path)


def _write_index_html(out: Path, manifest: dict[str, Any]) -> None:
    rows = []
    for key, value in sorted(manifest.get("outputs", {}).items()):
        target = Path(value)
        label = target.name
        rows.append(f'<li><a href="{label}">{key}</a></li>')
    skipped = "".join(
        f"<li><strong>{key}</strong>: {value}</li>"
        for key, value in sorted(manifest.get("skipped", {}).items())
    )
    errors = "".join(
        f"<li><strong>{key}</strong>: {value}</li>"
        for key, value in sorted(manifest.get("errors", {}).items())
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Atlas Dataset EDA</title>
  <style>
    body {{ font-family: system-ui, sans-serif; margin: 32px; color: #111827; }}
    code {{ background: #f3f4f6; padding: 2px 4px; border-radius: 4px; }}
    li {{ margin: 6px 0; }}
  </style>
</head>
<body>
  <h1>Atlas Dataset EDA</h1>
  <p><strong>Dataset:</strong> <code>{manifest.get("dataset", "")}</code></p>
  <p><strong>Rows analyzed:</strong> {manifest.get("shape", {}).get("rows_analyzed", "")}</p>
  <h2>Outputs</h2>
  <ul>
    {''.join(rows)}
  </ul>
  <h2>Skipped</h2>
  <ul>{skipped or "<li>None</li>"}</ul>
  <h2>Errors</h2>
  <ul>{errors or "<li>None</li>"}</ul>
</body>
</html>
"""
    (out / "index.html").write_text(html, encoding="utf-8")
