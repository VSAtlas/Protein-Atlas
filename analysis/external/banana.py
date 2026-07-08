from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd


FDA_MAPPING_COLUMNS = {
    "path",
    "generic_name",
    "display_name",
    "pubchem_name",
    "drugcentral_generic_name",
    "smiles_neutral",
    "smiles",
    "inchikey",
}


def _atlas_tools_dir() -> Path:
    repo_root = Path(__file__).resolve().parents[2]
    default_tools = repo_root.parent.parent / "tools"
    return Path(os.environ.get("ATLAS_TOOLS_DIR") or default_tools).expanduser()


def default_banana_root() -> Path:
    return Path(os.environ.get("BANANA_ROOT") or (_atlas_tools_dir() / "banana")).expanduser()


def default_banana_python() -> Path | None:
    explicit = os.environ.get("BANANA_PYTHON")
    if explicit:
        return Path(explicit).expanduser()
    candidate = _atlas_tools_dir() / "envs" / "banana" / "bin" / "python"
    return candidate if candidate.exists() else None


def _clean(value: Any) -> str:
    return " ".join(str(value or "").strip().split())


def _key(value: Any) -> str:
    return _clean(value).lower()


def _read_table(path: str | Path | None) -> pd.DataFrame:
    if path is None or not Path(path).exists():
        return pd.DataFrame()
    return pd.read_csv(path, low_memory=False)


def _load_fda_smiles(mapping_path: str | Path | None) -> pd.DataFrame:
    if mapping_path is None or not Path(mapping_path).exists():
        return pd.DataFrame(columns=["ligand_base", "drug_key", "smiles", "inchikey"])
    mapping = pd.read_csv(
        mapping_path,
        low_memory=False,
        usecols=lambda col: col in FDA_MAPPING_COLUMNS,
    )
    rows: list[dict[str, str]] = []
    for row in mapping.to_dict("records"):
        smiles = _clean(row.get("smiles_neutral")) or _clean(row.get("smiles"))
        if not smiles:
            continue
        ligand_base = Path(_clean(row.get("path"))).stem if _clean(row.get("path")) else ""
        names = {
            _key(row.get("generic_name")),
            _key(row.get("display_name")),
            _key(row.get("pubchem_name")),
            _key(row.get("drugcentral_generic_name")),
        }
        for name in sorted(name for name in names if name):
            rows.append(
                {
                    "ligand_base": ligand_base,
                    "drug_key": name,
                    "smiles": smiles,
                    "inchikey": _clean(row.get("inchikey")),
                }
            )
        if ligand_base:
            rows.append(
                {
                    "ligand_base": ligand_base,
                    "drug_key": "",
                    "smiles": smiles,
                    "inchikey": _clean(row.get("inchikey")),
                }
            )
    if not rows:
        return pd.DataFrame(columns=["ligand_base", "drug_key", "smiles", "inchikey"])
    return pd.DataFrame(rows).drop_duplicates(["ligand_base", "drug_key", "smiles"])


def _merge_smiles(pair: pd.DataFrame, smiles: pd.DataFrame) -> pd.DataFrame:
    out = pair.copy()
    if out.empty or smiles.empty:
        out["banana_smiles"] = pd.NA
        out["banana_input_status"] = "missing_fda_mapping_smiles"
        return out
    drug_names = pd.Series("", index=out.index, dtype="object")
    for col in ["generic_name", "mapped_drug_name", "display_name", "drug_id"]:
        if col not in out.columns:
            continue
        values = out[col].fillna("").astype(str).str.strip()
        drug_names = drug_names.where(drug_names.astype(str).str.len().gt(0), values)
    out["_banana_drug_key"] = drug_names.map(_key)
    out["_banana_ligand_base"] = (
        out.get("ligand_base", pd.Series("", index=out.index)).fillna("").astype(str)
    )

    by_base = smiles[smiles["ligand_base"].astype(str).str.len().gt(0)].drop_duplicates("ligand_base")
    merged = out.merge(
        by_base[["ligand_base", "smiles", "inchikey"]].rename(
            columns={
                "ligand_base": "_banana_ligand_base",
                "smiles": "banana_smiles",
                "inchikey": "banana_inchikey",
            }
        ),
        on="_banana_ligand_base",
        how="left",
    )
    by_name = smiles[smiles["drug_key"].astype(str).str.len().gt(0)].drop_duplicates("drug_key")
    missing = merged["banana_smiles"].isna()
    if missing.any():
        name_join = merged.loc[missing, ["_banana_drug_key"]].merge(
            by_name[["drug_key", "smiles", "inchikey"]].rename(
                columns={
                    "drug_key": "_banana_drug_key",
                    "smiles": "_banana_name_smiles",
                    "inchikey": "_banana_name_inchikey",
                }
            ),
            on="_banana_drug_key",
            how="left",
        )
        merged.loc[missing, "banana_smiles"] = name_join["_banana_name_smiles"].to_numpy()
        merged.loc[missing, "banana_inchikey"] = name_join["_banana_name_inchikey"].to_numpy()
    merged["banana_input_status"] = "ready_missing_pocket"
    merged.loc[merged["banana_smiles"].isna(), "banana_input_status"] = "missing_smiles"
    return merged.drop(columns=["_banana_drug_key", "_banana_ligand_base"], errors="ignore")


def build_banana_input_table(
    pair_table_path: str | Path,
    out_path: str | Path,
    *,
    fda_mapping_path: str | Path | None = "chemdb/data/fda_mapping_from_pdbqt.csv",
    pocket_map_path: str | Path | None = None,
) -> pd.DataFrame:
    """Build a BANANA-ready table with SMILES and optional pocket PDB paths.

    BANANA needs a SMILES string and a pocket-residue PDB file. This function does
    not infer missing pocket files; missing pockets remain explicit.
    """

    pair = _read_table(pair_table_path)
    if pair.empty:
        raise ValueError(f"pair table is empty or missing: {pair_table_path}")
    out = _merge_smiles(pair, _load_fda_smiles(fda_mapping_path))
    pocket = _read_table(pocket_map_path)
    if not pocket.empty and {"pdb_id", "pocket_pdb"}.issubset(pocket.columns):
        out = out.merge(
            pocket[["pdb_id", "pocket_pdb"]].drop_duplicates("pdb_id"),
            on="pdb_id",
            how="left",
        )
        has_smiles = out["banana_smiles"].notna()
        has_pocket = out["pocket_pdb"].fillna("").astype(str).str.len().gt(0)
        out.loc[has_smiles & has_pocket, "banana_input_status"] = "ready"
        out.loc[has_smiles & ~has_pocket, "banana_input_status"] = "missing_pocket_pdb"
    else:
        out["pocket_pdb"] = pd.NA
    keep = [
        col
        for col in [
            "drug_id",
            "target_id",
            "pdb_id",
            "ligand_base",
            "banana_smiles",
            "banana_inchikey",
            "pocket_pdb",
            "banana_input_status",
        ]
        if col in out.columns
    ]
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out[keep].to_csv(path, index=False)
    _write_summary(path.with_suffix(".summary.json"), out, "banana_input_status")
    return out[keep]


def _write_summary(path: Path, df: pd.DataFrame, status_col: str) -> None:
    status = df[status_col].fillna("missing").astype(str).value_counts().to_dict() if status_col in df.columns else {}
    payload = {
        "rows": int(len(df)),
        "status_counts": {str(k): int(v) for k, v in status.items()},
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _score_to_probability(value: float) -> float:
    if not math.isfinite(value):
        return float("nan")
    if value <= -99:
        return float("nan")
    return 1.0 / (1.0 + math.exp(-max(min(value, 40.0), -40.0)))


def run_banana_inference(
    banana_input_path: str | Path,
    out_path: str | Path,
    *,
    banana_root: str | Path | None = None,
    python_executable: str | Path | None = None,
    batch_size: int = 64,
    num_workers: int = 4,
    no_gpu: bool = True,
) -> pd.DataFrame:
    """Run BANANA inference for rows with ready SMILES and pocket PDB paths."""

    data = _read_table(banana_input_path)
    status = (
        data["banana_input_status"]
        if "banana_input_status" in data.columns
        else pd.Series("", index=data.index)
    )
    ready = data[status.astype(str).eq("ready")].copy()
    if ready.empty:
        raise ValueError("BANANA input table has no rows with banana_input_status=ready")
    root = Path(banana_root or default_banana_root()).resolve()
    if not (root / "inference.py").exists():
        raise FileNotFoundError(f"BANANA inference.py not found under {root}")
    python = str(python_executable or default_banana_python() or sys.executable)
    rows: list[pd.DataFrame] = []
    tmp_dir = Path(out_path).resolve().parent / "_banana_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    for pocket_pdb, group in ready.groupby("pocket_pdb"):
        pocket_path = Path(str(pocket_pdb))
        if not pocket_path.exists():
            continue
        safe = str(group["pdb_id"].iloc[0] if "pdb_id" in group.columns else pocket_path.stem)
        smi_path = tmp_dir / f"{safe}.smi"
        raw_out = tmp_dir / f"{safe}.banana_raw.txt"
        smi_path.write_text(
            "\n".join(group["banana_smiles"].astype(str).tolist()) + "\n",
            encoding="utf-8",
        )
        cmd = [
            python,
            "inference.py",
            str(smi_path.resolve()),
            str(pocket_path.resolve()),
            "--out_file",
            str(raw_out.resolve()),
            "--batch_size",
            str(batch_size),
            "--num_workers",
            str(num_workers),
        ]
        if no_gpu:
            cmd.append("--no_gpu")
        subprocess.run(cmd, cwd=root, check=True)
        scores = pd.to_numeric(pd.read_csv(raw_out, header=None)[0], errors="coerce")
        scored = group.reset_index(drop=True).copy()
        scored["banana_score"] = scores.reindex(scored.index).to_numpy()
        rows.append(scored)
    if not rows:
        raise RuntimeError("BANANA inference produced no scored rows")
    out = pd.concat(rows, ignore_index=True)
    out["banana_binding_probability"] = pd.to_numeric(out["banana_score"], errors="coerce").map(_score_to_probability)
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    _write_summary(path.with_suffix(".summary.json"), out, "banana_input_status")
    return out


def _normalize(series: pd.Series) -> pd.Series:
    vals = pd.to_numeric(series, errors="coerce")
    lo = vals.min(skipna=True)
    hi = vals.max(skipna=True)
    if pd.isna(lo) or pd.isna(hi) or hi == lo:
        return pd.Series(pd.NA, index=series.index, dtype="Float64")
    return ((vals - lo) / (hi - lo)).astype("Float64")


def _weighted_available_average(df: pd.DataFrame, weights: dict[str, float]) -> pd.Series:
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


def build_banana_binding_expert(
    pair_table_path: str | Path,
    out_path: str | Path,
    *,
    banana_scores_path: str | Path | None = None,
) -> pd.DataFrame:
    """Join BANANA scores and Atlas scores into a binding-expert table."""

    pair = _read_table(pair_table_path)
    if pair.empty:
        raise ValueError(f"pair table is empty or missing: {pair_table_path}")
    keep = [
        col
        for col in [
            "drug_id",
            "target_id",
            "pdb_id",
            "ligand_base",
            "atlas_score",
            "consensus_score",
            "SCORCH_score_used",
            "final_score",
        ]
        if col in pair.columns
    ]
    out = pair[keep].copy()
    scores = _read_table(banana_scores_path)
    if not scores.empty:
        score_cols = [
            col
            for col in [
                "drug_id",
                "target_id",
                "pdb_id",
                "banana_score",
                "banana_binding_probability",
                "banana_input_status",
            ]
            if col in scores.columns
        ]
        keys = [col for col in ["drug_id", "target_id", "pdb_id"] if col in out.columns and col in scores.columns]
        out = out.merge(scores[score_cols].drop_duplicates(keys), on=keys, how="left") if keys else out
    else:
        out["banana_input_status"] = "missing_banana_scores"
    if "banana_score" in out.columns:
        banana_values = pd.to_numeric(out["banana_score"], errors="coerce")
        out["banana_score_normalized"] = _normalize(out["banana_score"].mask(banana_values.le(-99)))
    if "banana_binding_probability" not in out.columns and "banana_score" in out.columns:
        out["banana_binding_probability"] = pd.to_numeric(out["banana_score"], errors="coerce").map(_score_to_probability)
    for col in ["atlas_score", "consensus_score", "SCORCH_score_used", "final_score"]:
        if col in out.columns:
            out[f"{col}_normalized"] = _normalize(out[col])
    out["banana_atlas_blend_score"] = _weighted_available_average(
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
    out["binding_expert_score"] = _weighted_available_average(
        out,
        {
            "banana_atlas_blend_score": 1.0,
            "atlas_score_normalized": 0.25,
            "consensus_score_normalized": 0.20,
        },
    )
    out["binding_expert_source"] = "banana_atlas_blend"
    out.loc[out.get("banana_score", pd.Series(pd.NA, index=out.index)).isna(), "binding_expert_source"] = "atlas_consensus_fallback"
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    manifest = {
        "rows": int(len(out)),
        "n_with_banana_score": int(out.get("banana_score", pd.Series(pd.NA, index=out.index)).notna().sum()),
        "n_with_binding_expert_score": int(out["binding_expert_score"].notna().sum()),
        "interpretation": (
            "BANANA/BigBind is used as a binding/activity expert. Atlas, consensus, "
            "SCORCH, and final scores are complementary calibration/prioritization layers."
        ),
    }
    path.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return out
