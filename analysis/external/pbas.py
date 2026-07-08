from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.external.mapping import apply_pair_mapping, normalize_columns, read_optional_mapping
from analysis.external.source_tables import read_source_table


PBAS_COMMENT = (
    "Sawada et al. introduced PBAS-style proteome-wide docking profiles: "
    "drug-level vectors of predicted binding affinities across human proteins, "
    "used for therapeutic indication and side-effect prediction. Atlas uses "
    "PBAS/Sawada-style profiles as a baseline and adds pair-level statistical "
    "significance, Vina/SCORCH consensus, MM/GBSA refinement, exposure "
    "plausibility, tissue expression, target-ADR evidence, pathway evidence, "
    "and mechanism-graph interpretation."
)


def _read_matrix(path: str | Path) -> pd.DataFrame:
    df = read_source_table(path)
    if {"drug_id", "target_id"}.issubset(df.columns):
        return df
    first = df.columns[0]
    return df.melt(id_vars=[first], var_name="target_id", value_name="pbas_score").rename(columns={first: "drug_id"})


def _load_metadata(path: str | Path | None, aliases: dict[str, list[str]]) -> pd.DataFrame:
    if path is None or not Path(path).exists():
        return pd.DataFrame()
    return normalize_columns(read_source_table(path), aliases)


def build_pbas_pair_scores(
    pbas_matrix_path: str | Path,
    out_path: str | Path,
    *,
    drug_metadata_path: str | Path | None = None,
    target_metadata_path: str | Path | None = None,
    pocket_info_path: str | Path | None = None,
    mapping_path: str | Path | None = None,
) -> pd.DataFrame:
    pbas = normalize_columns(
        _read_matrix(pbas_matrix_path),
        {
            "drug_id": ["drug_id", "drug", "compound_id", "chembl_id", "pbas_drug_id"],
            "target_id": ["target_id", "target", "uniprot", "gene", "protein_id", "pbas_target_id"],
            "pbas_score": ["pbas_score", "score", "vina_score", "binding_score"],
        },
    )
    pbas = apply_pair_mapping(pbas, read_optional_mapping(mapping_path))
    pbas["pbas_score"] = pd.to_numeric(pbas["pbas_score"], errors="coerce")
    pbas["pbas_available"] = pbas["pbas_score"].notna()
    pbas["pbas_source"] = "Sawada_PBAS"
    pbas["pbas_structure_source"] = "AlphaFold/PBAS"

    drug_meta = _load_metadata(
        drug_metadata_path,
        {"drug_id": ["drug_id", "drug", "compound_id", "chembl_id", "pbas_drug_id"], "pbas_drug_name": ["name", "drug_name"]},
    )
    if not drug_meta.empty:
        pbas = pbas.merge(drug_meta.drop_duplicates("drug_id"), on="drug_id", how="left")

    target_meta = _load_metadata(
        target_metadata_path,
        {"target_id": ["target_id", "target", "uniprot", "gene", "protein_id", "pbas_target_id"], "pbas_target_name": ["name", "gene", "protein_name"]},
    )
    if not target_meta.empty:
        pbas = pbas.merge(target_meta.drop_duplicates("target_id"), on="target_id", how="left")

    pocket = _load_metadata(
        pocket_info_path,
        {"target_id": ["target_id", "target", "uniprot", "protein_id"], "pbas_pocket_id": ["pocket_id", "pocket", "site_id"]},
    )
    if not pocket.empty:
        pbas = pbas.merge(pocket.drop_duplicates("target_id"), on="target_id", how="left")
    if "pbas_pocket_id" not in pbas:
        pbas["pbas_pocket_id"] = pd.NA

    scored = pbas["pbas_score"].notna()
    pbas["pbas_rank_percentile"] = pd.NA
    if scored.any():
        pbas.loc[scored, "pbas_rank_percentile"] = pbas.loc[scored, "pbas_score"].rank(pct=True, ascending=False)

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "drug_id",
        "target_id",
        "pbas_score",
        "pbas_rank_percentile",
        "pbas_source",
        "pbas_structure_source",
        "pbas_pocket_id",
        "pbas_available",
    ]
    extra = [col for col in pbas.columns if col not in columns]
    pbas[[*columns, *extra]].to_csv(out, index=False)
    _coverage_summary(pbas).to_csv(out.with_name("pbas_coverage_summary.csv"), index=False)
    return pbas


def _coverage_summary(df: pd.DataFrame) -> pd.DataFrame:
    unmapped_drug = df["drug_id"].isna() | df["drug_id"].astype(str).str.strip().eq("")
    unmapped_target = df["target_id"].isna() | df["target_id"].astype(str).str.strip().eq("")
    return pd.DataFrame(
        [
            {
                "source": "Sawada_PBAS",
                "n_rows": len(df),
                "n_drugs": int(df["drug_id"].nunique(dropna=True)),
                "n_targets": int(df["target_id"].nunique(dropna=True)),
                "n_available_scores": int(df["pbas_score"].notna().sum()),
                "n_unmapped_drugs": int(unmapped_drug.sum()),
                "n_unmapped_targets": int(unmapped_target.sum()),
                "manifest_note": PBAS_COMMENT,
            }
        ]
    )
