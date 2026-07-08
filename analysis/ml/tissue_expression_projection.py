from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.tissue_site_labels import add_independent_tissue_site_labels
from analysis.ml.tissue_site_scoring import add_tissue_site_scores


DEFAULT_TARGET_EXPRESSION_CANDIDATES = [
    Path("data/external/hpa_gtex/tissue_expression.tsv"),
    Path("data/external/hpa_gtex/tissue_expression.csv"),
    Path("data/external/hpa/tissue_expression.tsv"),
    Path("data/external/hpa/tissue_expression.csv"),
    Path("data/external/gtex/tissue_expression.tsv"),
    Path("data/external/gtex/tissue_expression.csv"),
    Path("data/external/bgee/tissue_expression.tsv"),
    Path("data/external/bgee/tissue_expression.csv"),
    Path("data/external/opentargets/expression.tsv"),
    Path("data/external/opentargets/expression.csv"),
    Path("data/tissue_expression/tissue_expression.tsv"),
    Path("data/tissue_expression/tissue_expression.csv"),
]

DEFAULT_TARGET_ADR_CONTEXT_CANDIDATES = [
    Path("data/external/adrecs_target/target_adr_evidence.tsv"),
    Path("data/external/paper_sources/normalized_v2/combined_target_adr_evidence.tsv"),
    Path("data/external/opentargets/safety.tsv"),
]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


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


def _read_table(path: Path) -> pd.DataFrame:
    sep = "	" if path.suffix.lower() in {".tsv", ".tab", ".txt"} else ","
    return pd.read_csv(path, sep=sep, low_memory=False)


def _first_nonempty_frame(df: pd.DataFrame, cols: list[str]) -> pd.Series:
    out = pd.Series("", index=df.index, dtype="object")
    for col in cols:
        if col not in df.columns:
            continue
        values = df[col].map(_clean)
        out = out.where(out.astype(str).str.len().gt(0), values)
    return out


def _merge_expression(pair_table: pd.DataFrame, expression: pd.DataFrame) -> pd.DataFrame:
    if expression.empty:
        return pair_table.copy()
    left = pair_table.copy()
    right = expression.copy()
    if "target_gene" not in left.columns:
        left["target_gene"] = _first_nonempty_frame(left, ["gene_symbol", "target_id", "target_uniprot"]).map(_upper)
    if "target_gene" not in right.columns:
        for col in ["gene_symbol", "gene", "target_id", "uniprot", "target_uniprot"]:
            if col in right.columns:
                right["target_gene"] = right[col].map(_upper)
                break
    if "target_gene" not in right.columns:
        return left
    for key_set in [["target_gene", "adr_site_group"], ["target_gene", "site_name"], ["target_gene"]]:
        if set(key_set).issubset(left.columns) and set(key_set).issubset(right.columns):
            right_subset = right.drop_duplicates(key_set, keep="first")
            return left.merge(right_subset, on=key_set, how="left", suffixes=("", "_expression"))
    right_subset = right.drop_duplicates(["target_gene"], keep="first")
    return left.merge(right_subset, on="target_gene", how="left", suffixes=("", "_expression"))


def _read_existing_tables(candidates: list[Path], run_dir: Path | None = None) -> list[pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    roots = []
    if run_dir is not None:
        roots.append(run_dir)
    roots.extend([Path.cwd(), _repo_root()])
    for rel in candidates:
        path = rel if rel.is_absolute() else None
        paths = [path] if path is not None else [root / rel for root in roots]
        for candidate in paths:
            if candidate is not None and candidate.exists():
                frames.append(_read_table(candidate))
                break
    return frames


def _target_adr_context(run_dir: Path | None = None) -> pd.DataFrame:
    frames = _read_existing_tables(DEFAULT_TARGET_ADR_CONTEXT_CANDIDATES, run_dir)
    rows: list[pd.DataFrame] = []
    for frame in frames:
        if frame.empty:
            continue
        work = frame.copy()
        if "benchmark_only" in work.columns:
            benchmark = work["benchmark_only"].fillna(False).astype(str).str.lower().isin({"1", "true", "yes"})
            work = work[~benchmark].copy()
        if "label_state" in work.columns:
            label_state = pd.to_numeric(work["label_state"], errors="coerce")
            work = work[label_state.eq(1) | label_state.isna()].copy()
        if work.empty:
            continue
        target_gene = _first_nonempty_frame(work, ["gene_symbol", "target_gene", "approved_symbol", "_source_approved_symbol"]).map(_upper)
        target_id = _first_nonempty_frame(work, ["target_id", "uniprot", "target_uniprot"]).map(_upper)
        adr_term = _first_nonempty_frame(work, ["adr_term", "adr", "disease_id", "event", "safety_event"])
        source = _first_nonempty_frame(work, ["source_name", "source", "upstream_source"])
        context = pd.DataFrame(
            {
                "target_gene": target_gene,
                "target_id": target_id,
                "_target_adr_term": adr_term,
                "_target_adr_source": source,
            }
        )
        context = context[context["_target_adr_term"].astype(str).str.len().gt(0)]
        rows.append(context)
    if not rows:
        return pd.DataFrame()
    combined = pd.concat(rows, ignore_index=True)
    combined = combined[combined["target_gene"].astype(str).str.len().gt(0) | combined["target_id"].astype(str).str.len().gt(0)]
    if combined.empty:
        return pd.DataFrame()

    def join_limited(values: pd.Series, limit: int = 40) -> str:
        seen: list[str] = []
        for value in values.dropna().astype(str):
            for token in value.split(";"):
                cleaned = _clean(token)
                if cleaned and cleaned not in seen:
                    seen.append(cleaned)
                if len(seen) >= limit:
                    return ";".join(seen)
        return ";".join(seen)

    aggregations = {
        "target_adr_terms_aggregated": ("_target_adr_term", join_limited),
        "target_adr_sources_aggregated": ("_target_adr_source", join_limited),
        "target_adr_evidence_count": ("_target_adr_term", "nunique"),
    }
    by_gene = combined[combined["target_gene"].astype(str).str.len().gt(0)].groupby("target_gene", dropna=False).agg(**aggregations).reset_index()
    by_id = combined[combined["target_id"].astype(str).str.len().gt(0)].groupby("target_id", dropna=False).agg(**aggregations).reset_index()
    return pd.concat([by_gene, by_id], ignore_index=True, sort=False).drop_duplicates()


def _merge_target_adr_context(source: pd.DataFrame, run_dir: Path | None = None) -> tuple[pd.DataFrame, dict[str, Any]]:
    context = _target_adr_context(run_dir)
    if context.empty:
        return source.copy(), {"status": "no_target_adr_context", "target_adr_context_rows": 0}
    left = source.copy()
    if "target_gene" not in left.columns:
        left["target_gene"] = _first_nonempty_frame(left, ["gene_symbol", "target_id", "target_uniprot"]).map(_upper)
    merged = left
    if "target_gene" in context.columns:
        gene_context = context.dropna(subset=["target_gene"]).drop_duplicates("target_gene", keep="first")
        if not gene_context.empty:
            merged = merged.merge(gene_context, on="target_gene", how="left", suffixes=("", "_target_adr_gene"))
    if "target_id" in context.columns and "target_id" in merged.columns:
        id_context = context.dropna(subset=["target_id"]).drop_duplicates("target_id", keep="first")
        if not id_context.empty:
            merged = merged.merge(id_context, on="target_id", how="left", suffixes=("", "_target_adr_id"))
    for col in ["target_adr_terms_aggregated", "target_adr_sources_aggregated", "target_adr_evidence_count"]:
        id_col = f"{col}_target_adr_id"
        gene_col = col
        if id_col in merged.columns:
            if gene_col in merged.columns:
                merged[gene_col] = merged[gene_col].where(merged[gene_col].notna(), merged[id_col])
                merged = merged.drop(columns=[id_col])
            else:
                merged = merged.rename(columns={id_col: gene_col})
    terms = merged.get("target_adr_terms_aggregated", pd.Series(pd.NA, index=merged.index))
    if "adr_term" not in merged.columns:
        merged["adr_term"] = terms
    else:
        empty_adr = merged["adr_term"].isna() | merged["adr_term"].astype(str).str.strip().eq("")
        merged.loc[empty_adr, "adr_term"] = terms.loc[empty_adr]
    source_values = merged.get("target_adr_sources_aggregated", pd.Series(pd.NA, index=merged.index))
    if "adr_site_mapping_source" not in merged.columns:
        merged["adr_site_mapping_source"] = pd.NA
    has_terms = terms.notna() & terms.astype(str).str.strip().ne("")
    empty_source = merged["adr_site_mapping_source"].isna() | merged["adr_site_mapping_source"].astype(str).str.strip().eq("")
    merged.loc[has_terms & empty_source, "adr_site_mapping_source"] = (
        "target_adr_evidence:" + source_values.fillna("local_target_adr_sources").astype(str)
    )
    return merged, {
        "status": "joined",
        "target_adr_context_rows": int(len(context)),
        "rows_with_target_adr_terms": int(has_terms.sum()),
    }


def add_tissue_expression_context(
    source: pd.DataFrame,
    *,
    target_expression_path: str | Path | None = None,
    run_dir: str | Path | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Add deterministic tissue-expression score columns when local expression data exists.

    Missing source files or expression rows remain missing. This function does
    not create supervised tissue labels from expression rules.
    """

    expression_path = _resolve_existing(
        target_expression_path,
        DEFAULT_TARGET_EXPRESSION_CANDIDATES,
        Path(run_dir) if run_dir else None,
    )
    expression = pd.DataFrame()
    if expression_path is not None:
        expression = _read_table(expression_path)
    with_adr_context, adr_context_summary = _merge_target_adr_context(source, Path(run_dir) if run_dir else None)
    merged = _merge_expression(with_adr_context, expression) if not expression.empty else with_adr_context.copy()
    scored, score_summary = add_tissue_site_scores(merged)
    scored, label_summary = add_independent_tissue_site_labels(scored, run_dir=run_dir)
    summary: dict[str, Any] = {
        "status": "scored" if expression_path is not None else "no_expression_source",
        "target_expression_path": str(expression_path) if expression_path is not None else None,
        "input_rows": int(len(source)),
        "rows": int(len(scored)),
        "expression_source_rows": int(len(expression)),
        "site_relevance_nonmissing": int(scored.get("site_relevance_score", pd.Series(index=scored.index)).notna().sum()),
        "expression_presence_nonmissing": int(scored.get("expression_presence_score", pd.Series(index=scored.index)).notna().sum()),
        "site_specificity_nonmissing": int(scored.get("site_specificity_score", pd.Series(index=scored.index)).notna().sum()),
        "tissue_site_label_nonmissing": int(scored.get("tissue_site_label", pd.Series(index=scored.index)).notna().sum()) if "tissue_site_label" in scored.columns else 0,
        "target_adr_context": adr_context_summary,
        "independent_tissue_site_labels": label_summary,
        "policy": "Expression scores are deterministic context features only. Supervised tissue labels, when present, are projected from independent ADR/clinical site evidence.",
    }
    summary.update({f"scoring_{key}": value for key, value in score_summary.items()})
    return scored.loc[:, ~scored.columns.duplicated()].copy(), summary
