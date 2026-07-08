from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.external.compound_index import load_atlas_compound_index, map_by_inchikey


def _smiles_inchikey(smiles: str) -> str:
    try:
        from rdkit import Chem
    except Exception:
        return ""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return ""
    return str(Chem.MolToInchiKey(mol) or "")


def _stable_compound_id(smiles: str) -> str:
    digest = hashlib.sha1(str(smiles).encode("utf-8")).hexdigest()[:16]
    return f"muv_smiles:{digest}"


def stage_muv_benchmark(
    muv_csv_path: str | Path,
    out_path: str | Path,
    *,
    mapping_path: str | Path | None = None,
    map_to_atlas: bool = True,
) -> pd.DataFrame:
    """Normalize MoleculeNet MUV into long-form benchmark-only source rows."""

    raw = pd.read_csv(muv_csv_path, low_memory=False)
    if "SMILES" not in raw.columns:
        raise ValueError("MUV CSV must contain a SMILES column")
    task_cols = [col for col in raw.columns if str(col).startswith("MUV-")]
    if not task_cols:
        raise ValueError("MUV CSV must contain MUV-* task columns")
    compounds = raw[["SMILES"]].drop_duplicates().copy()
    compounds["inchikey"] = compounds["SMILES"].map(_smiles_inchikey)
    mapped = compounds.copy()
    if mapping_path is not None and map_to_atlas:
        atlas_index = load_atlas_compound_index(mapping_path)
        mapped = map_by_inchikey(compounds, atlas_index, source_inchikey_col="inchikey")
    mapped["drug_id"] = mapped.get("drug_id", pd.Series(pd.NA, index=mapped.index))
    mapped["drug_id"] = mapped["drug_id"].fillna(mapped["SMILES"].map(_stable_compound_id))
    mapped["drug_name"] = mapped.get("drug_name", pd.Series(pd.NA, index=mapped.index))
    mapped["drug_name"] = mapped["drug_name"].fillna(mapped["drug_id"])
    mapped["atlas_mapped"] = ~mapped["drug_id"].astype(str).str.startswith("muv_smiles:")
    id_map = mapped[["SMILES", "drug_id", "drug_name", "inchikey", "atlas_mapped"]].drop_duplicates("SMILES")
    long = raw.melt(id_vars=["SMILES"], value_vars=task_cols, var_name="target_id", value_name="activity_outcome")
    long = long[long["activity_outcome"].notna()].copy()
    long["activity_outcome"] = pd.to_numeric(long["activity_outcome"], errors="coerce")
    long = long[long["activity_outcome"].isin([0, 1])].copy()
    long = long.merge(id_map, on="SMILES", how="left")
    out = pd.DataFrame(
        {
            "drug_id": long["drug_id"],
            "drug_name": long["drug_name"],
            "inchikey": long["inchikey"],
            "smiles": long["SMILES"],
            "target_id": long["target_id"],
            "assay_id": long["target_id"],
            "activity_outcome": long["activity_outcome"].astype(int),
            "activity_type": "MUV benchmark binary task",
            "activity_relation": "",
            "activity_units": "",
            "source": "MUV",
            "atlas_mapped": long["atlas_mapped"].fillna(False).astype(bool),
        }
    ).drop_duplicates()
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep="\t", index=False)
    manifest: dict[str, Any] = {
        "source": "MoleculeNet MUV",
        "input": str(muv_csv_path),
        "output": str(path),
        "rows": int(len(out)),
        "wide_rows": int(len(raw)),
        "tasks": task_cols,
        "atlas_mapped_rows": int(out["atlas_mapped"].sum()),
        "benchmark_only": True,
        "label_policy": "MUV labels are virtual-screening benchmark labels and must not become production Atlas truth labels.",
    }
    path.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return out
