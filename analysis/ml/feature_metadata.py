from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.feature_backfill import backfill_training_features
from analysis.ml.chemotype_backfill import add_structural_ligand_chemotypes
from analysis.ml.chemical_clusters import add_butina_chemical_clusters
from analysis.ml.ligand_descriptors import add_ligand_physchem_descriptors


UNASSIGNED = {"", "nan", "none", "null", "other / unassigned", "other-unassigned"}

GENE_FAMILY_RULES: list[tuple[str, tuple[str, ...]]] = [
    ("GPCR", ("ADRA", "ADRB", "CHRM", "DRD", "HTR", "HRH", "OPR", "SSTR", "CCR", "CXCR", "PTGER", "TACR")),
    ("ION_CHANNEL", ("KCN", "SCN", "CACN", "GABR", "GRIA", "GRIN", "CHRNA", "CHRNB", "HTR3")),
    ("KINASE", ("AKT", "ALK", "ABL", "BRAF", "CDK", "EGFR", "ERBB", "FGFR", "FLT", "JAK", "KIT", "MAPK", "MET", "MTOR", "PIK3", "SRC", "SYK", "VEGFR")),
    ("NUCLEAR_RECEPTOR", ("AR", "ESR", "NR", "PPAR", "RARA", "RXR", "THRA", "THRB", "VDR")),
    ("TRANSPORTER", ("ABCB", "ABCC", "ABCG", "SLC")),
    ("CYP_ENZYME", ("CYP",)),
    ("PROTEASE", ("ACE", "CASP", "F", "MMP", "PLG", "REN", "TMPRSS", "TPS")),
    ("OXIDOREDUCTASE", ("COX", "MAO", "PTGS")),
    ("PHOSPHATASE", ("PTP", "PPP", "DUSP")),
]

SOURCE_FAMILY_RULES: list[tuple[str, tuple[str, ...]]] = [
    ("spd", ("spd", "secondary pharmacology")),
    ("toxcast", ("toxcast", "invitrodb", "epa")),
    ("chembl", ("chembl",)),
    ("papyrus", ("papyrus",)),
    ("bindingdb", ("bindingdb",)),
    ("pubchem", ("pubchem", "bioassay")),
    ("adrecs", ("adrecs",)),
    ("opentargets", ("open targets", "opentargets")),
    ("ctd", ("ctd",)),
    ("sider", ("sider",)),
]

TARGET_FAMILY_CANONICAL = {
    "gpcr": "GPCR",
    "ion channel": "Ion Channel",
    "ion_channel": "Ion Channel",
    "kinase": "Kinase",
    "nuclear receptor": "Nuclear Hormone Receptor",
    "nuclear hormone receptor": "Nuclear Hormone Receptor",
    "nuclear_receptor": "Nuclear Hormone Receptor",
    "transporter": "Transporter",
    "enzyme": "Enzyme",
    "cyp enzyme": "Enzyme",
    "cyp_enzyme": "Enzyme",
    "oxidoreductase": "Enzyme",
    "phosphatase": "Enzyme",
    "protease": "Protease",
}


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _clean_lower(value: Any) -> str:
    return _clean(value).lower()


def _is_missing(value: Any) -> bool:
    return _clean_lower(value) in UNASSIGNED


def _canonical_target_family(value: Any) -> str:
    cleaned = _clean(value)
    key = cleaned.casefold().replace("-", " ")
    return TARGET_FAMILY_CANONICAL.get(key, cleaned)


def _first_present(row: pd.Series, fields: list[str]) -> tuple[str, str]:
    for field in fields:
        if field in row.index and not _is_missing(row.get(field)):
            return _clean(row.get(field)), field
    return "", ""


def _murcko_from_smiles(smiles: str) -> str:
    if not _clean(smiles):
        return ""
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold

        mol = Chem.MolFromSmiles(str(smiles).strip())
        if mol is None:
            return ""
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        return scaffold or Chem.MolToSmiles(mol, isomericSmiles=False)
    except Exception:
        return ""


def _target_family_from_gene(value: Any) -> tuple[str, str]:
    gene = _clean(value).upper()
    if not gene:
        return "", ""
    for family, prefixes in GENE_FAMILY_RULES:
        if any(gene.startswith(prefix) for prefix in prefixes):
            return family, "gene_symbol_heuristic"
    return "OTHER_TARGET_FAMILY", "gene_symbol_fallback"


def _source_family_from_row(row: pd.Series) -> tuple[str, str]:
    text = " ".join(
        _clean_lower(row.get(col))
        for col in ["source_family", "upstream_source", "label_source", "source_objective", "source"]
        if col in row.index
    )
    for family, tokens in SOURCE_FAMILY_RULES:
        if any(token in text for token in tokens):
            return family, "source_text_rule"
    return "", ""


def enrich_ml_feature_metadata(
    df: pd.DataFrame,
    *,
    chemical_cluster: str = "auto",
    target_family: str = "auto",
    source_lineage: str = "auto",
    drop_columns: list[str] | None = None,
    run_dir: Path | None = None,
    repo_root: Path | None = None,
    extra_feature_tables: list[str | Path] | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    out = df.copy()
    summary: dict[str, Any] = {
        "input_rows": int(len(out)),
        "chemical_cluster_mode": chemical_cluster,
        "target_family_mode": target_family,
        "source_lineage_mode": source_lineage,
    }

    if drop_columns:
        removable = [col for col in drop_columns if col in out.columns]
        out = out.drop(columns=removable)
        summary["dropped_columns"] = removable

    out, backfill_summary = backfill_training_features(
        out,
        run_dir=run_dir,
        repo_root=repo_root,
        extra_feature_tables=extra_feature_tables,
    )
    summary["feature_backfill"] = backfill_summary

    out, descriptor_summary = add_ligand_physchem_descriptors(out, repo_root=repo_root)
    summary.update(descriptor_summary)

    out, chemotype_summary = add_structural_ligand_chemotypes(out, repo_root=repo_root)
    summary.update(chemotype_summary)

    if chemical_cluster != "none":
        out = _add_chemical_cluster(out, mode=chemical_cluster)
        if out.attrs.get("chemical_cluster_summary"):
            summary.update(out.attrs["chemical_cluster_summary"])
    if target_family != "none":
        out = _add_target_family(out, mode=target_family)
    if source_lineage != "none":
        out = _add_source_lineage(out, mode=source_lineage)

    for col in ["chemical_cluster", "target_family", "source_family", "upstream_source"]:
        if col in out.columns:
            summary[f"{col}_nonmissing"] = int(out[col].map(lambda value: not _is_missing(value)).sum())
            summary[f"{col}_unique"] = int(out[col].nunique(dropna=True))
    return out, summary


def _add_chemical_cluster(df: pd.DataFrame, *, mode: str) -> pd.DataFrame:
    if mode in {"butina", "ecfp"}:
        clustered, summary = add_butina_chemical_clusters(df)
        clustered.attrs["chemical_cluster_summary"] = summary
        return clustered

    out = df.copy()
    existing = out["chemical_cluster"].map(_clean) if "chemical_cluster" in out.columns else pd.Series("", index=out.index)
    assigned = existing.map(lambda value: not _is_missing(value))
    values = existing.where(assigned, "")
    sources = pd.Series("existing", index=out.index, dtype="object").where(assigned, "")

    unresolved = ~assigned
    if mode in {"auto", "scaffold"} and unresolved.any():
        scaffold = _first_available_column(out, ["scaffold_key", "murcko_scaffold", "chemical_scaffold"]).map(_clean)
        mask = unresolved & scaffold.map(lambda value: not _is_missing(value))
        values.loc[mask] = scaffold.loc[mask]
        sources.loc[mask] = "scaffold_column"
        unresolved = values.map(_is_missing)

    if mode in {"auto", "smiles"} and unresolved.any():
        smiles = _first_available_column(out, ["smiles", "canonical_smiles", "ligand_smiles"]).map(_clean)
        unique_smiles = smiles.loc[unresolved & smiles.astype(str).str.len().gt(0)].drop_duplicates()
        scaffold_lookup = {value: _murcko_from_smiles(value) for value in unique_smiles.tolist()}
        scaffold = smiles.map(lambda value: scaffold_lookup.get(value, ""))
        mask = unresolved & scaffold.astype(str).str.len().gt(0)
        values.loc[mask] = "murcko_smiles:" + scaffold.loc[mask].astype(str)
        sources.loc[mask] = "smiles"
        unresolved = values.map(_is_missing)

    if mode in {"auto", "chemotype"} and unresolved.any():
        chemotype = _first_available_column(out, ["ligand_chemotype", "chemotype"]).map(_clean)
        mask = unresolved & chemotype.map(lambda value: not _is_missing(value))
        values.loc[mask] = "chemotype:" + chemotype.loc[mask].astype(str)
        sources.loc[mask] = "ligand_chemotype"
        unresolved = values.map(_is_missing)

    values = values.where(values.map(lambda value: not _is_missing(value)), "chemical_cluster_unassigned")
    sources = sources.where(sources.astype(str).str.len().gt(0), "unassigned")
    out["chemical_cluster"] = values
    out["chemical_cluster_source"] = sources
    return out


def _first_available_column(df: pd.DataFrame, fields: list[str]) -> pd.Series:
    out = pd.Series("", index=df.index, dtype="object")
    for field in fields:
        if field not in df.columns:
            continue
        values = df[field].map(_clean)
        out = out.where(out.astype(str).str.len().gt(0), values)
    return out

def _add_target_family(df: pd.DataFrame, *, mode: str) -> pd.DataFrame:
    out = df.copy()
    existing = out["target_family"].map(_clean) if "target_family" in out.columns else pd.Series("", index=out.index)
    assigned = existing.map(lambda value: not _is_missing(value))
    values = existing.where(assigned, "")
    sources = pd.Series("existing", index=out.index, dtype="object").where(assigned, "")

    unresolved = ~assigned
    if mode in {"auto", "protein_class"} and unresolved.any():
        for field in ["protein_family", "target_class", "protein_class"]:
            if field not in out.columns:
                continue
            field_values = out[field].map(_clean)
            mask = unresolved & field_values.map(lambda value: not _is_missing(value))
            values.loc[mask] = field_values.loc[mask]
            sources.loc[mask] = f"{field}_column"
            unresolved = values.map(_is_missing)

    if mode in {"auto", "gene_heuristic"} and unresolved.any():
        gene_values = _first_available_column(out, ["target_gene", "gene_symbol", "target_id"]).map(_clean)
        unique_genes = gene_values.loc[unresolved & gene_values.astype(str).str.len().gt(0)].drop_duplicates()
        family_lookup = {gene: _target_family_from_gene(gene) for gene in unique_genes.tolist()}
        families = gene_values.map(lambda gene: family_lookup.get(gene, ("", ""))[0])
        family_sources = gene_values.map(lambda gene: family_lookup.get(gene, ("", ""))[1])
        mask = unresolved & families.astype(str).str.len().gt(0)
        values.loc[mask] = families.loc[mask]
        sources.loc[mask] = family_sources.loc[mask]

    values = values.where(values.map(lambda value: not _is_missing(value)), "TARGET_FAMILY_UNASSIGNED")
    values = values.map(_canonical_target_family)
    sources = sources.where(sources.astype(str).str.len().gt(0), "unassigned")
    out["target_family"] = values
    out["target_family_source"] = sources
    return out


def _add_source_lineage(df: pd.DataFrame, *, mode: str) -> pd.DataFrame:
    out = df.copy()
    if "source_family" not in out.columns:
        evidence_text = _combine_lower_text(
            out,
            ["evidence_sources", "parent_sources", "label_source", "source", "source_name", "upstream_source"],
        )
        families, family_sources = _source_families_from_text(evidence_text, source_label="evidence_source_text_rule")
        unresolved = families.eq("source_family_unassigned")
        if unresolved.any():
            objective_text = _combine_lower_text(out, ["source_objective"])
            fallback_families, fallback_sources = _source_families_from_text(
                objective_text,
                source_label="source_objective_text_rule",
            )
            fallback_mask = unresolved & fallback_families.ne("source_family_unassigned")
            families.loc[fallback_mask] = fallback_families.loc[fallback_mask]
            family_sources.loc[fallback_mask] = fallback_sources.loc[fallback_mask]
        out["source_family"] = families
        out["source_family_source"] = family_sources
    if "upstream_source" not in out.columns:
        for source_col in ["evidence_sources", "parent_sources", "label_source", "source", "source_name", "source_objective"]:
            if source_col in out.columns:
                out["upstream_source"] = out[source_col]
                out["upstream_source_source"] = source_col
                break
    return out


def _combine_lower_text(df: pd.DataFrame, fields: list[str]) -> pd.Series:
    text = pd.Series("", index=df.index, dtype="object")
    for field in fields:
        if field in df.columns:
            text = text.str.cat(df[field].map(_clean_lower), sep=" ")
    return text.str.strip()


def _source_families_from_text(text: pd.Series, *, source_label: str) -> tuple[pd.Series, pd.Series]:
    families = pd.Series("source_family_unassigned", index=text.index, dtype="object")
    sources = pd.Series("unassigned", index=text.index, dtype="object")
    match_count = pd.Series(0, index=text.index, dtype="int64")
    first_family = pd.Series("", index=text.index, dtype="object")
    for family, tokens in SOURCE_FAMILY_RULES:
        mask = pd.Series(False, index=text.index)
        for token in tokens:
            mask = mask | text.str.contains(token, regex=False, na=False)
        first_family = first_family.where(match_count.ne(0), family)
        match_count = match_count + mask.astype(int)
    single = match_count.eq(1)
    multiple = match_count.gt(1)
    families.loc[single] = first_family.loc[single]
    sources.loc[single] = source_label
    families.loc[multiple] = "multi_source"
    sources.loc[multiple] = f"{source_label}_multi"
    return families, sources


def refresh_tables(
    paths: list[Path],
    *,
    out_dir: Path | None = None,
    in_place: bool = False,
    run_dir: Path | None = None,
    chemical_cluster: str = "auto",
    target_family: str = "auto",
    source_lineage: str = "auto",
    drop_columns: list[str] | None = None,
    extra_feature_tables: list[str | Path] | None = None,
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    outputs: list[str] = []
    for path in paths:
        df = pd.read_csv(path, low_memory=False)
        enriched, summary = enrich_ml_feature_metadata(
            df,
            chemical_cluster=chemical_cluster,
            target_family=target_family,
            source_lineage=source_lineage,
            drop_columns=drop_columns,
            run_dir=run_dir,
            extra_feature_tables=extra_feature_tables,
        )
        if in_place:
            out_path = path
        else:
            base = out_dir or path.parent
            if run_dir is not None:
                try:
                    rel = path.resolve().relative_to(run_dir.resolve())
                except ValueError:
                    rel = Path(path.name)
            else:
                rel = Path(path.name)
            out_path = base / rel
        out_path.parent.mkdir(parents=True, exist_ok=True)
        enriched.to_csv(out_path, index=False)
        summary.update({"input": str(path), "output": str(out_path), "columns": int(len(enriched.columns))})
        rows.append(summary)
        outputs.append(str(out_path))
    manifest = {
        "tables": rows,
        "outputs": outputs,
        "n_tables": len(paths),
        "chemical_cluster": chemical_cluster,
        "target_family": target_family,
        "source_lineage": source_lineage,
        "in_place": in_place,
    }
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "ml_feature_metadata_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        pd.DataFrame(rows).to_csv(out_dir / "ml_feature_metadata_summary.csv", index=False)
    return manifest
