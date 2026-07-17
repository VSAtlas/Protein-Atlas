from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.calibration.metrics import calibration_metrics


def _first_col(df: pd.DataFrame, names: list[str]) -> str | None:
    for name in names:
        if name in df.columns and df[name].notna().any():
            return name
    return None


def _fingerprint(smiles: str):
    try:
        from rdkit import Chem, DataStructs
        from rdkit.Chem import rdFingerprintGenerator

        mol = Chem.MolFromSmiles(str(smiles or "").strip())
        if mol is None:
            return None, DataStructs
        generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        return generator.GetFingerprint(mol), DataStructs
    except Exception:
        return None, None


def append_chemical_fingerprint_ad(
    pred: pd.DataFrame,
    train: pd.DataFrame,
    test: pd.DataFrame,
    *,
    label_col: str,
    out_path: str | Path,
) -> pd.DataFrame:
    smiles_col = _first_col(train, ["smiles", "canonical_smiles", "ligand_smiles"])
    if smiles_col is None or smiles_col not in test.columns:
        return pred
    train_fps = []
    data_structs = None
    for smiles in train[smiles_col].dropna().astype(str).unique():
        fp, data_structs = _fingerprint(smiles)
        if fp is not None:
            train_fps.append(fp)
    if not train_fps or data_structs is None:
        return pred
    similarities: list[float | None] = []
    for smiles in test[smiles_col].astype(str).tolist():
        fp, _ = _fingerprint(smiles)
        if fp is None:
            similarities.append(None)
            continue
        sims = data_structs.BulkTanimotoSimilarity(fp, train_fps)
        similarities.append(float(max(sims)) if sims else None)
    out = pred.copy()
    out["chemical_fingerprint_train_max_tanimoto"] = similarities
    out["chemical_fingerprint_ad"] = [
        "in_domain" if value is not None and pd.notna(value) and float(value) >= 0.4 else "out_of_domain"
        for value in similarities
    ]
    _write_group_summary(
        out,
        label_col,
        "chemical_fingerprint_ad",
        Path(out_path) / "chemical_fingerprint_ad_summary.csv",
    )
    return out


def append_target_family_ad(
    pred: pd.DataFrame,
    train: pd.DataFrame,
    test: pd.DataFrame,
    *,
    label_col: str,
    out_path: str | Path,
) -> pd.DataFrame:
    family_col = _first_col(train, ["target_family", "protein_family", "target_class", "protein_class"])
    if family_col is None or family_col not in test.columns:
        return pred
    train_families = set(train[family_col].dropna().astype(str))
    out = pred.copy()
    if family_col not in out.columns:
        out[family_col] = test[family_col].values
    out["target_family_seen_in_train"] = test[family_col].astype(str).isin(train_families).values
    out["target_family_ad"] = out["target_family_seen_in_train"].map(lambda value: "in_domain" if value else "out_of_domain")
    _write_group_summary(
        out,
        label_col,
        "target_family_ad",
        Path(out_path) / "target_similarity_ad_summary.csv",
    )
    return out


def _write_group_summary(df: pd.DataFrame, label_col: str, group_col: str, out_path: Path) -> None:
    if group_col not in df.columns or label_col not in df.columns:
        return
    rows = []
    for value, group in df.groupby(group_col, dropna=False):
        labels = pd.to_numeric(group[label_col], errors="coerce")
        scores = pd.to_numeric(group["ml_prediction_score"], errors="coerce") if "ml_prediction_score" in group else pd.Series(dtype=float)
        valid = labels.notna() & scores.notna()
        metrics = (
            calibration_metrics(scores.loc[valid].tolist(), labels.loc[valid].astype(int).tolist())
            if valid.sum() and labels.loc[valid].nunique() > 1
            else {}
        )
        rows.append(
            {
                group_col: value,
                "n": int(len(group)),
                "positive_rate": float(labels.mean()) if labels.notna().any() else None,
                **metrics,
            }
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_path, index=False)


def append_abstention_flags(pred: pd.DataFrame, *, out_path: str | Path | None = None) -> pd.DataFrame:
    """Add conservative review flags for predictions outside training support.

    The score is retained. The flag says a prediction should be interpreted as
    out-of-domain until a closer chemical/target/source neighborhood exists.
    """

    out = pred.copy()
    reason_cols = [
        "applicability_domain",
        "chemical_fingerprint_ad",
        "target_family_ad",
    ]
    reasons: list[str] = []
    for _, row in out.iterrows():
        row_reasons = []
        for col in reason_cols:
            if col in out.columns and str(row.get(col, "")).strip().lower() == "out_of_domain":
                row_reasons.append(col)
        reasons.append(";".join(row_reasons))
    out["abstention_reasons"] = reasons
    out["abstention_flag"] = ["review_out_of_domain" if reason else "score_in_domain" for reason in reasons]
    if out_path is not None:
        rows = []
        for flag, group in out.groupby("abstention_flag", dropna=False):
            rows.append({"abstention_flag": flag, "n": int(len(group))})
        Path(out_path).mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(Path(out_path) / "model_abstention_summary.csv", index=False)
    return out
