from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.external.mapping import normalize_columns
from analysis.external.source_tables import read_source_table


PK_ALIASES = {
    "drug_id": ["drug_id", "atlas_drug_id", "ligand_base", "compound_id", "drugbank_id", "substance", "name"],
    "drug_name": ["drug_name", "display_name", "generic_name", "compound_name", "name", "substance"],
    "cmax_um": ["cmax_um", "cmax_uM", "free_total_cmax_um", "total_cmax_um", "combined_cmax_um"],
    "free_cmax_um": ["free_cmax_um", "free_cmax_uM", "cmax_free_um", "combined_free_cmax_um"],
    "fraction_unbound_plasma": ["fraction_unbound_plasma", "fraction_unbound", "fu", "fu_plasma"],
    "protein_binding_percent": ["protein_binding_percent", "protein_binding", "bound_percent", "value"],
    "pk_parameter": ["pk_parameter", "parameter", "parameter_name", "endpoint", "metric", "measurement", "test"],
    "unit": ["unit", "normalized_unit", "cmax_unit"],
    "normalized_value": ["normalized_value", "cmax_normalized_value", "value"],
    "source": ["source", "pk_source", "exposure_source"],
}


def _text_series(df: pd.DataFrame, col: str) -> pd.Series:
    return df.get(col, pd.Series("", index=df.index)).fillna("").astype(str).str.strip()


def _unit_to_cmax_um(value: pd.Series, unit: pd.Series) -> pd.Series:
    raw = pd.to_numeric(value, errors="coerce")
    normalized_unit = unit.fillna("").astype(str).str.strip().str.lower()
    out = pd.Series(pd.NA, index=value.index, dtype="Float64")
    out = out.where(~normalized_unit.isin({"um", "µm", "umol/l", "µmol/l", "micromolar"}), raw)
    out = out.where(~normalized_unit.isin({"nm", "nmol/l", "nanomolar"}), raw / 1000.0)
    out = out.where(~normalized_unit.isin({"mm", "mmol/l", "millimolar"}), raw * 1000.0)
    out = out.where(~normalized_unit.isin({"m", "mol/l", "molar"}), raw * 1_000_000.0)
    return out


def _fraction_unbound_from_binding(series: pd.Series) -> pd.Series:
    binding = pd.to_numeric(series, errors="coerce")
    fraction = 1.0 - (binding / 100.0)
    return fraction.where(fraction.between(0.0, 1.0))


def normalize_pk_table(path: str | Path, *, source_name: str) -> pd.DataFrame:
    raw = read_source_table(path)
    df = normalize_columns(raw, PK_ALIASES)
    raw_columns = {str(col).strip().lower() for col in raw.columns}
    out = pd.DataFrame(index=df.index)
    out["drug_id"] = _text_series(df, "drug_id")
    out["drug_name"] = _text_series(df, "drug_name")
    out["cmax_um"] = pd.to_numeric(df.get("cmax_um"), errors="coerce")
    converted = _unit_to_cmax_um(df.get("normalized_value", pd.Series(pd.NA, index=df.index)), _text_series(df, "unit"))
    parameter = _text_series(df, "pk_parameter").str.lower()
    has_parameter = parameter.astype(str).str.len().gt(0).any()
    if has_parameter:
        cmax_like = parameter.str.contains(
            r"\bc[\s_-]*max\b|max(?:imum)? concentration|peak plasma",
            regex=True,
            na=False,
        )
        free_cmax_like = cmax_like & parameter.str.contains(r"free|unbound", regex=True, na=False)
        out["cmax_um"] = out["cmax_um"].fillna(converted.where(cmax_like & ~free_cmax_like))
    else:
        out["cmax_um"] = out["cmax_um"].fillna(converted)
    out["free_cmax_um"] = pd.to_numeric(df.get("free_cmax_um"), errors="coerce")
    if has_parameter:
        out["free_cmax_um"] = out["free_cmax_um"].fillna(converted.where(free_cmax_like))
    out["fraction_unbound_plasma"] = pd.to_numeric(df.get("fraction_unbound_plasma"), errors="coerce")
    has_binding_field = bool(
        {"protein_binding_percent", "protein binding percent", "protein_binding", "bound_percent"} & raw_columns
    )
    if has_binding_field or "protein_binding" in source_name.lower():
        out["fraction_unbound_plasma"] = out["fraction_unbound_plasma"].fillna(
            _fraction_unbound_from_binding(df.get("protein_binding_percent", pd.Series(pd.NA, index=df.index)))
        )
    out["free_cmax_um"] = out["free_cmax_um"].fillna(out["cmax_um"] * out["fraction_unbound_plasma"])
    source_text = _text_series(df, "source")
    out["exposure_source"] = source_text.where(source_text.astype(str).str.len() > 0, source_name)
    out["pk_source"] = source_name
    out = out[(out["drug_id"].astype(str).str.len() > 0) | (out["drug_name"].astype(str).str.len() > 0)].copy()
    return out


def build_combined_pk_table(
    out_path: str | Path,
    *,
    existing_pk: str | Path | None = None,
    spd_pk: str | Path | None = None,
    openfda_pk: str | Path | None = None,
    pkdb: str | Path | None = None,
    ncats_inxight_pk: str | Path | None = None,
    drugbank_cmax: str | Path | None = None,
    drugbank_protein_binding: str | Path | None = None,
) -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    for source_name, path in (
        ("existing_pk", existing_pk),
        ("SPD", spd_pk),
        ("openFDA_SPL", openfda_pk),
        ("PK-DB", pkdb),
        ("NCATS_Inxight_FRDB", ncats_inxight_pk),
        ("DrugBank_Cmax", drugbank_cmax),
        ("DrugBank_protein_binding", drugbank_protein_binding),
    ):
        if path is not None and Path(path).exists():
            parts.append(normalize_pk_table(path, source_name=source_name))
    if not parts:
        combined = pd.DataFrame(
            columns=[
                "drug_id",
                "drug_name",
                "free_cmax_um",
                "cmax_um",
                "fraction_unbound_plasma",
                "exposure_source",
                "pk_source",
            ]
        )
    else:
        combined = pd.concat(parts, ignore_index=True)
        combined["_pk_key"] = combined["drug_id"].where(
            combined["drug_id"].astype(str).str.len() > 0,
            combined["drug_name"].str.lower(),
        )
        combined["_has_free"] = pd.to_numeric(combined["free_cmax_um"], errors="coerce").notna()
        combined["_has_cmax_fu"] = pd.to_numeric(combined["cmax_um"], errors="coerce").notna() & pd.to_numeric(
            combined["fraction_unbound_plasma"], errors="coerce"
        ).notna()
        combined = combined.sort_values(["_pk_key", "_has_free", "_has_cmax_fu"], ascending=[True, False, False])
        combined = combined.drop_duplicates("_pk_key", keep="first").drop(columns=["_pk_key", "_has_free", "_has_cmax_fu"])
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out, index=False)
    summary = pd.DataFrame(
        [
            {
                "n_rows": len(combined),
                "n_free_cmax": int(pd.to_numeric(combined.get("free_cmax_um"), errors="coerce").notna().sum()),
                "n_cmax": int(pd.to_numeric(combined.get("cmax_um"), errors="coerce").notna().sum()),
                "n_fraction_unbound": int(
                    pd.to_numeric(combined.get("fraction_unbound_plasma"), errors="coerce").notna().sum()
                ),
                "sources": ";".join(sorted({str(v) for v in combined.get("pk_source", pd.Series(dtype=object)).dropna()})),
            }
        ]
    )
    summary.to_csv(out.with_name(out.stem + ".summary.csv"), index=False)
    return combined
