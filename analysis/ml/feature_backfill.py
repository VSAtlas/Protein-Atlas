from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd

MISSING_TEXT = {"", "nan", "none", "null", "na", "n/a", "other / unassigned", "other-unassigned"}
SENSITIVE_SITE_TOKENS = ("heart", "liver", "kidney", "brain", "cns", "lung", "reproductive", "blood", "immune")
SCORE_COLUMNS = ["atlas_score", "consensus_score", "SCORCH_score_used", "final_score", "banana_score"]
PAIR_FEATURE_COLUMNS = [
    "atlas_score",
    "consensus_score",
    "SCORCH_score_used",
    "final_score",
    "banana_score",
    "banana_binding_probability",
    "banana_score_normalized",
    "banana_atlas_blend_score",
    "binding_expert_score",
    "z_selected",
    "z_selected_source",
]
TARGET_FEATURE_COLUMNS = [
    "protein_class",
    "target_family",
    "target_gene",
    "target_uniprot",
    "structure_quality",
    "structure_quality_source",
    "pdb_resolution",
    "pdb_rank_score",
    "pdb_selection_method",
]
LIGAND_FEATURE_COLUMNS = [
    "ligand_chemotype",
    "ligand_chemotype_source",
    "scaffold_key",
    "scaffold_source",
    "smiles",
    "canonical_smiles",
    "display_name",
    "generic_name",
]


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _lower(value: Any) -> str:
    return _clean(value).lower()


def _is_missing(value: Any) -> bool:
    return _lower(value) in MISSING_TEXT


def _missing_mask(series: pd.Series) -> pd.Series:
    text = series.astype("string").str.strip().str.lower()
    return series.isna() | text.isin(MISSING_TEXT)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _first_present(row: pd.Series, columns: list[str]) -> str:
    for col in columns:
        if col in row.index:
            value = _clean(row.get(col))
            if value:
                return value
    return ""


def _ensure_aliases(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "drug_id" not in out.columns and "ligand_base" in out.columns:
        out["drug_id"] = out["ligand_base"]
        out["drug_id_source"] = "ligand_base_alias"
    if "target_id" not in out.columns:
        for col in ("target_gene", "gene_symbol", "target_uniprot", "pdb_id"):
            if col in out.columns:
                out["target_id"] = out[col]
                out["target_id_source"] = f"{col}_alias"
                break
    if "atlas_score" not in out.columns:
        for col in ("z_selected", "final_score", "consensus_score"):
            if col in out.columns:
                out["atlas_score"] = pd.to_numeric(out[col], errors="coerce")
                out["atlas_score_source_for_ml"] = col
                break
    if "consensus_score" not in out.columns:
        for col in ("final_score", "z_selected"):
            if col in out.columns:
                out["consensus_score"] = pd.to_numeric(out[col], errors="coerce")
                out["consensus_score_source_for_ml"] = col
                break
    return out


def _normalize(series: pd.Series) -> pd.Series:
    vals = pd.to_numeric(series, errors="coerce")
    lo = vals.min(skipna=True)
    hi = vals.max(skipna=True)
    if pd.isna(lo) or pd.isna(hi) or hi == lo:
        return pd.Series(pd.NA, index=series.index, dtype="Float64")
    return ((vals - lo) / (hi - lo)).astype("Float64")


def _score_to_probability(value: Any) -> float | None:
    score = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(score):
        return None
    # BANANA reports a binding-like score; this monotone transform is only a
    # compact prior used when a calibrated BANANA probability is absent.
    import math

    return float(1.0 / (1.0 + math.exp(-float(score))))


def _weighted_average(df: pd.DataFrame, weights: dict[str, float]) -> pd.Series:
    numer = pd.Series(0.0, index=df.index)
    denom = pd.Series(0.0, index=df.index)
    for col, weight in weights.items():
        if col not in df.columns:
            continue
        vals = pd.to_numeric(df[col], errors="coerce")
        present = vals.notna()
        numer = numer.add(vals.fillna(0.0) * float(weight))
        denom = denom.add(present.astype(float) * float(weight))
    return (numer / denom).where(denom.gt(0), pd.NA).astype("Float64")


def _existing_candidate_tables(run_dir: Path | None, extra_tables: list[Path] | None) -> list[Path]:
    candidates: list[Path] = []
    if run_dir is not None:
        candidates.extend(
            [
                run_dir / "master_rows.csv",
                run_dir / "heatmap_input.csv",
                run_dir / "publication_heatmap_input.csv",
                run_dir / "banana_scores.csv",
                run_dir / "banana_binding_expert.csv",
                run_dir / "binding_expert.csv",
                run_dir / "ml" / "banana_scores.csv",
                run_dir / "ml" / "banana_binding_expert.csv",
            ]
        )
        for folder in (run_dir / "banana", run_dir / "banana_inputs", run_dir / "ml"):
            if folder.exists():
                candidates.extend(sorted(folder.glob("*banana*.csv"))[:50])
    candidates.extend(extra_tables or [])
    out: list[Path] = []
    seen: set[Path] = set()
    for path in candidates:
        try:
            resolved = path.resolve()
        except Exception:
            resolved = path
        if resolved in seen or not path.exists() or not path.is_file():
            continue
        seen.add(resolved)
        out.append(path)
    return out


def _read_candidate(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, low_memory=False)
    except Exception:
        return pd.DataFrame()
    if df.empty:
        return df
    df = _ensure_aliases(df)
    if "pdb_id" in df.columns:
        df["pdb_id"] = df["pdb_id"].map(_clean)
    if "ligand_base" in df.columns:
        df["ligand_base"] = df["ligand_base"].map(_clean)
    if "drug_id" in df.columns:
        df["drug_id"] = df["drug_id"].map(_clean)
    if "target_id" in df.columns:
        df["target_id"] = df["target_id"].map(_clean)
    return df


def _fill_from_source(
    out: pd.DataFrame,
    source: pd.DataFrame,
    *,
    keys: list[str],
    columns: list[str],
    source_name: str,
) -> tuple[pd.DataFrame, dict[str, int]]:
    if source.empty or not set(keys).issubset(out.columns) or not set(keys).issubset(source.columns):
        return out, {}
    fill_cols = [col for col in columns if col in source.columns]
    if not fill_cols:
        return out, {}
    left = out.copy()
    for col in fill_cols:
        if col not in left.columns:
            left[col] = pd.NA
    right = source[[*keys, *fill_cols]].copy()
    for key in keys:
        right[key] = right[key].map(_clean)
    right = right.drop_duplicates(keys, keep="first")
    merged = left.merge(right, on=keys, how="left", suffixes=("", "__fill"))
    counts: dict[str, int] = {}
    for col in fill_cols:
        fill_col = f"{col}__fill"
        if fill_col not in merged.columns:
            continue
        if col not in merged.columns:
            merged[col] = pd.NA
        before = int((~_missing_mask(merged[col])).sum())
        missing = _missing_mask(merged[col])
        fill_present = ~_missing_mask(merged[fill_col])
        merged.loc[missing & fill_present, col] = merged.loc[missing & fill_present, fill_col]
        after = int((~_missing_mask(merged[col])).sum())
        added = int(after - before)
        if added > 0:
            counts[col] = added
            source_col = f"{col}_feature_source"
            if source_col not in merged.columns:
                merged[source_col] = ""
            merged.loc[missing & fill_present, source_col] = source_name
    return merged.drop(columns=[col for col in merged.columns if col.endswith("__fill")]), counts


def _fill_pair_features(df: pd.DataFrame, run_dir: Path | None, extra_tables: list[Path] | None) -> tuple[pd.DataFrame, dict[str, Any]]:
    out = df.copy()
    summary: dict[str, Any] = {"candidate_tables": [], "filled": {}}
    key_sets = [
        ["drug_id", "target_id", "pdb_id"],
        ["drug_id", "pdb_id"],
        ["ligand_base", "target_id", "pdb_id"],
        ["ligand_base", "pdb_id"],
    ]
    columns = [*PAIR_FEATURE_COLUMNS, *TARGET_FEATURE_COLUMNS, *LIGAND_FEATURE_COLUMNS]
    for path in _existing_candidate_tables(run_dir, extra_tables):
        source = _read_candidate(path)
        if source.empty:
            continue
        source_summary: dict[str, Any] = {"path": str(path), "columns": int(len(source.columns)), "filled": {}}
        for keys in key_sets:
            out, counts = _fill_from_source(out, source, keys=keys, columns=columns, source_name=str(path))
            for col, value in counts.items():
                source_summary["filled"][col] = int(source_summary["filled"].get(col, 0)) + int(value)
                summary["filled"][col] = int(summary["filled"].get(col, 0)) + int(value)
        if source_summary["filled"]:
            summary["candidate_tables"].append(source_summary)
    return out, summary


def _murcko_scaffold_key(smiles: str) -> str:
    if not smiles:
        return ""
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        if not scaffold:
            scaffold = Chem.MolToSmiles(mol, isomericSmiles=False)
        digest = hashlib.sha1(scaffold.encode("utf-8")).hexdigest()[:12]
        return f"murcko:{digest}"
    except Exception:
        return ""


def _load_ligand_mapping(repo_root: Path) -> dict[str, dict[str, str]]:
    path = repo_root / "chemdb/data/fda_mapping_from_pdbqt.csv"
    if not path.exists():
        return {}
    wanted = {
        "path",
        "generic_name",
        "display_name",
        "pubchem_name",
        "drugcentral_generic_name",
        "smiles_neutral",
        "smiles",
        "canonical_smiles",
    }
    try:
        mapping = pd.read_csv(path, low_memory=False, usecols=lambda col: col in wanted)
    except Exception:
        return {}
    lookup: dict[str, dict[str, str]] = {}
    for row in mapping.to_dict("records"):
        smiles = _clean(row.get("smiles_neutral") or row.get("canonical_smiles") or row.get("smiles"))
        record = {
            "smiles": smiles,
            "canonical_smiles": smiles,
            "scaffold_key": _murcko_scaffold_key(smiles),
            "scaffold_source": "fda_mapping_from_pdbqt_murcko" if smiles else "",
            "generic_name": _clean(row.get("generic_name") or row.get("drugcentral_generic_name")),
            "display_name": _clean(row.get("display_name") or row.get("pubchem_name")),
        }
        path_text = _clean(row.get("path"))
        if path_text:
            lookup[f"base:{Path(path_text).stem.lower()}"] = record
        for field in ("generic_name", "display_name", "pubchem_name", "drugcentral_generic_name"):
            key = _lower(row.get(field))
            if key:
                lookup[f"name:{key}"] = record
    return lookup


def _ligand_mapping_for_row(row: pd.Series, lookup: dict[str, dict[str, str]]) -> dict[str, str]:
    for field, key_type in (
        ("ligand_base", "base"),
        ("drug_id", "name"),
        ("display_name", "name"),
        ("generic_name", "name"),
    ):
        value = _lower(row.get(field)) if field in row.index else ""
        if not value:
            continue
        found = lookup.get(f"{key_type}:{value}")
        if found:
            return found
    return {}


def _lookup_records_for_rows(out: pd.DataFrame, lookup: dict[str, dict[str, str]]) -> pd.Series:
    keys = pd.Series("", index=out.index, dtype="object")
    if "ligand_base" in out.columns:
        base_keys = "base:" + out["ligand_base"].fillna("").astype(str).str.strip().str.lower()
        keys = keys.where(base_keys.str.len().le(5), base_keys)
    for col in ("drug_id", "display_name", "generic_name"):
        if col not in out.columns:
            continue
        name_keys = "name:" + out[col].fillna("").astype(str).str.strip().str.lower()
        keys = keys.where(keys.astype(str).str.len().gt(0), name_keys.where(name_keys.str.len().gt(5), ""))
    return keys.map(lambda key: lookup.get(str(key), {}))


def _record_values(records: pd.Series, field: str) -> pd.Series:
    return records.map(lambda rec: rec.get(field, "") if isinstance(rec, dict) else "")


def _fill_ligand_features(df: pd.DataFrame, repo_root: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    out = df.copy()
    for col in ["ligand_chemotype", "ligand_chemotype_source", "scaffold_key", "scaffold_source"]:
        if col not in out.columns:
            out[col] = ""
    for col in ("smiles", "canonical_smiles", "display_name", "generic_name"):
        if col not in out.columns:
            out[col] = pd.NA
    lookup = _load_ligand_mapping(repo_root)
    records = _lookup_records_for_rows(out, lookup)
    filled = {"ligand_chemotype": 0, "scaffold_key": 0, "smiles": 0}

    for col in ("smiles", "canonical_smiles", "display_name", "generic_name"):
        values = _record_values(records, col)
        missing = _missing_mask(out[col])
        present = ~_missing_mask(values)
        out.loc[missing & present, col] = values.loc[missing & present]
        if col in {"smiles", "canonical_smiles"}:
            filled["smiles"] += int((missing & present).sum())

    scaffold_values = _record_values(records, "scaffold_key")
    scaffold_missing = _missing_mask(out["scaffold_key"])
    scaffold_present = ~_missing_mask(scaffold_values)
    out.loc[scaffold_missing & scaffold_present, "scaffold_key"] = scaffold_values.loc[scaffold_missing & scaffold_present]
    out.loc[scaffold_missing & scaffold_present, "scaffold_source"] = "fda_mapping_from_pdbqt_murcko"
    filled["scaffold_key"] = int((scaffold_missing & scaffold_present).sum())

    chemotype_missing = _missing_mask(out["ligand_chemotype"])
    scaffold_available = ~_missing_mask(out["scaffold_key"])
    chemotype_fill = chemotype_missing & scaffold_available
    out.loc[chemotype_fill, "ligand_chemotype"] = "Other / Structural scaffold assigned"
    out.loc[chemotype_fill, "ligand_chemotype_source"] = "fda_mapping_murcko"
    filled["ligand_chemotype"] = int(chemotype_fill.sum())
    return out, {"ligand_mapping_rows": int(len(lookup)), "filled": filled}


def _target_metadata_files(repo_root: Path, run_dir: Path | None) -> list[Path]:
    candidates: list[Path] = []
    if run_dir is not None:
        candidates.extend(
            [
                run_dir / "target_install" / "spd_targets_selected.csv",
                run_dir / "target_install" / "selected.csv",
                run_dir / "spd_targets_selected.csv",
            ]
        )
    gene_list = repo_root / "analysis/gene_list"
    if gene_list.exists():
        candidates.extend(sorted(gene_list.glob("spd*selected*.csv"))[:100])
        candidates.extend(sorted(gene_list.glob("spd*allpass*.csv"))[:50])
    out: list[Path] = []
    seen: set[Path] = set()
    for path in candidates:
        if path.exists() and path.is_file() and path.resolve() not in seen:
            seen.add(path.resolve())
            out.append(path)
    return out


def _load_target_metadata(repo_root: Path, run_dir: Path | None) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for path in _target_metadata_files(repo_root, run_dir):
        try:
            df = pd.read_csv(path, low_memory=False)
        except Exception:
            continue
        if "pdb_id" not in df.columns:
            continue
        keep = [col for col in ["gene", "target_gene", "pdb_id", "rank_score", "resolution", "method"] if col in df.columns]
        part = df[keep].copy()
        part["_target_metadata_source"] = str(path)
        rows.append(part)
    if not rows:
        return pd.DataFrame()
    meta = pd.concat(rows, ignore_index=True)
    meta["pdb_id"] = meta["pdb_id"].map(_clean)
    meta["pdb_rank_score"] = pd.to_numeric(meta.get("rank_score"), errors="coerce")
    meta["pdb_resolution"] = pd.to_numeric(meta.get("resolution"), errors="coerce")
    rank_quality = (meta["pdb_rank_score"] / 100.0).clip(lower=0.0, upper=1.0)
    resolution_quality = (1.0 - ((meta["pdb_resolution"] - 1.0) / 4.0)).clip(lower=0.0, upper=1.0)
    meta["structure_quality"] = rank_quality.where(rank_quality.notna(), resolution_quality)
    meta["structure_quality_source"] = meta["_target_metadata_source"]
    if "target_gene" not in meta.columns and "gene" in meta.columns:
        meta["target_gene"] = meta["gene"]
    if "pdb_selection_method" not in meta.columns and "method" in meta.columns:
        meta["pdb_selection_method"] = meta["method"]
    sort_cols = ["pdb_id", "structure_quality"]
    return meta.sort_values(sort_cols, ascending=[True, False]).drop_duplicates("pdb_id", keep="first")


def _fill_target_features(df: pd.DataFrame, repo_root: Path, run_dir: Path | None) -> tuple[pd.DataFrame, dict[str, Any]]:
    meta = _load_target_metadata(repo_root, run_dir)
    if meta.empty:
        return df.copy(), {"target_metadata_files": 0, "filled": {}}
    out, counts = _fill_from_source(
        df.copy(),
        meta,
        keys=["pdb_id"],
        columns=["target_gene", "structure_quality", "structure_quality_source", "pdb_resolution", "pdb_rank_score", "pdb_selection_method"],
        source_name="target_install_metadata",
    )
    return out, {"target_metadata_rows": int(len(meta)), "filled": counts}


def _fill_sensitive_offsite(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    out = df.copy()
    if "sensitive_offsite_expression_score" in out.columns and out["sensitive_offsite_expression_score"].notna().any():
        return out, {"sensitive_offsite_expression_score": 0}
    site = pd.Series("", index=out.index, dtype="object")
    for col in ("site_name", "adr_site_group", "adr_soc"):
        if col in out.columns:
            site = site.str.cat(out[col].fillna("").astype(str), sep=" ")
    values = site.str.lower().map(lambda text: any(token in text for token in SENSITIVE_SITE_TOKENS) if text.strip() else pd.NA)
    out["sensitive_offsite_expression_score"] = values.astype("Float64")
    return out, {"sensitive_offsite_expression_score": int(out["sensitive_offsite_expression_score"].notna().sum())}


def _finalize_score_priors(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    out = df.copy()
    filled: dict[str, int] = {}
    for col in SCORE_COLUMNS:
        if col in out.columns:
            norm_col = f"{col}_normalized"
            if norm_col not in out.columns or out[norm_col].isna().all():
                out[norm_col] = _normalize(out[col])
                filled[norm_col] = int(out[norm_col].notna().sum())
    if "banana_binding_probability" not in out.columns and "banana_score" in out.columns:
        out["banana_binding_probability"] = pd.to_numeric(out["banana_score"], errors="coerce").map(_score_to_probability)
        filled["banana_binding_probability"] = int(out["banana_binding_probability"].notna().sum())
    if "banana_atlas_blend_score" not in out.columns:
        out["banana_atlas_blend_score"] = _weighted_average(
            out,
            {
                "banana_binding_probability": 0.50,
                "banana_score_normalized": 0.20,
                "atlas_score_normalized": 0.15,
                "consensus_score_normalized": 0.10,
                "SCORCH_score_used_normalized": 0.20,
                "final_score_normalized": 0.10,
            },
        )
        filled["banana_atlas_blend_score"] = int(out["banana_atlas_blend_score"].notna().sum())
    if "binding_expert_score" not in out.columns:
        out["binding_expert_score"] = _weighted_average(
            out,
            {
                "banana_atlas_blend_score": 1.0,
                "atlas_score_normalized": 0.25,
                "consensus_score_normalized": 0.20,
            },
        )
        filled["binding_expert_score"] = int(out["binding_expert_score"].notna().sum())
    if "binding_expert_source" not in out.columns and "binding_expert_score" in out.columns:
        out["binding_expert_source"] = "feature_backfill_score_prior"
    return out, {"filled": filled}


def backfill_training_features(
    df: pd.DataFrame,
    *,
    run_dir: str | Path | None = None,
    repo_root: str | Path | None = None,
    extra_feature_tables: list[str | Path] | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Backfill reusable ML feature columns from local Atlas run artifacts.

    This function only joins existing local evidence. It never launches expensive
    tools such as BANANA inference and never converts missing data to negative
    evidence. Missing columns remain missing and are reported in the summary.
    """

    root = Path(repo_root).resolve() if repo_root is not None else _repo_root()
    run = Path(run_dir).resolve() if run_dir is not None else None
    out = _ensure_aliases(df)
    if "pdb_id" in out.columns:
        out["pdb_id"] = out["pdb_id"].map(_clean)
    if "ligand_base" in out.columns:
        out["ligand_base"] = out["ligand_base"].map(_clean)
    if "drug_id" in out.columns:
        out["drug_id"] = out["drug_id"].map(_clean)
    summary: dict[str, Any] = {
        "input_rows": int(len(out)),
        "run_dir": str(run) if run else None,
        "repo_root": str(root),
        "policy": "Backfills local feature metadata only; missing data remain unknown and heavy inference is not auto-launched.",
    }
    extras = [Path(path) for path in extra_feature_tables or []]
    out, pair_summary = _fill_pair_features(out, run, extras)
    out, ligand_summary = _fill_ligand_features(out, root)
    out, target_summary = _fill_target_features(out, root, run)
    out, sensitive_summary = _fill_sensitive_offsite(out)
    out, score_summary = _finalize_score_priors(out)
    summary.update(
        {
            "pair_feature_backfill": pair_summary,
            "ligand_feature_backfill": ligand_summary,
            "target_feature_backfill": target_summary,
            "sensitive_offsite_backfill": sensitive_summary,
            "score_prior_backfill": score_summary,
        }
    )
    tracked = [
        "banana_score_normalized",
        "binding_expert_score",
        "structure_quality",
        "target_family",
        "ligand_chemotype",
        "rdkit_mol_wt",
        "sensitive_offsite_expression_score",
    ]
    summary["tracked_feature_nonmissing"] = {
        col: int(out[col].notna().sum()) if col in out.columns else 0 for col in tracked
    }
    return out, summary


def write_full_feature_table(
    input_path: str | Path,
    out_path: str | Path,
    *,
    run_dir: str | Path | None = None,
    repo_root: str | Path | None = None,
    extra_feature_tables: list[str | Path] | None = None,
) -> pd.DataFrame:
    source = pd.read_csv(input_path, low_memory=False)
    enriched, summary = backfill_training_features(
        source,
        run_dir=run_dir,
        repo_root=repo_root,
        extra_feature_tables=extra_feature_tables,
    )
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    enriched.to_csv(out, index=False)
    summary.update({"input_path": str(input_path), "out_path": str(out), "output_columns": int(len(enriched.columns))})
    out.with_suffix(".manifest.json").write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    return enriched
