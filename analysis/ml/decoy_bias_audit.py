from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.calibration.metrics import calibration_metrics


DECOY_MARKERS = ("decoy", "dud-e", "dude", "deepcoy", "dekois", "muv")
PROPERTY_COLS = [
    "molecular_weight",
    "mol_weight",
    "mw",
    "logp",
    "clogp",
    "tpsa",
    "hbd",
    "hba",
    "rotatable_bonds",
    "num_rotatable_bonds",
    "formal_charge",
    "heavy_atom_count",
    "ring_count",
    "qed",
]
LIGAND_ONLY_COLS = [
    "ligand_chemotype",
    "scaffold_key",
    "chemical_cluster",
    "ligand_cluster",
    "molecular_weight",
    "mol_weight",
    "mw",
    "logp",
    "clogp",
    "tpsa",
    "hbd",
    "hba",
    "rotatable_bonds",
    "num_rotatable_bonds",
    "formal_charge",
    "heavy_atom_count",
    "ring_count",
    "qed",
]


def _truthy_text(value: object) -> bool:
    text = str(value or "").strip().lower()
    return text in {"1", "true", "yes", "y", "decoy"} or any(marker in text for marker in DECOY_MARKERS)


def _decoy_labels(df: pd.DataFrame) -> pd.Series:
    label = pd.Series(False, index=df.index)
    for col in ["is_decoy", "benchmark_only", "source_role", "label_source", "negative_evidence_type", "source_name"]:
        if col in df.columns:
            label |= df[col].map(_truthy_text)
    return label.astype(int)


def _standardized_mean_difference(a: pd.Series, b: pd.Series) -> float | None:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if a.empty or b.empty:
        return None
    pooled = ((float(a.var()) + float(b.var())) / 2.0) ** 0.5
    if pooled == 0:
        pooled = 1.0
    return float((a.mean() - b.mean()) / pooled)


def _property_matching(df: pd.DataFrame, decoy_col: str, out: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_cols = [col for col in ["target_id", "pdb_id", "decoy_source"] if col in df.columns]
    groups = df.groupby(group_cols, dropna=False) if group_cols else [("all", df)]
    for group_key, group in groups:
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        base = {col: value for col, value in zip(group_cols, group_key, strict=False)}
        active = group[group[decoy_col].eq(0)]
        decoy = group[group[decoy_col].eq(1)]
        for prop in PROPERTY_COLS:
            if prop not in group.columns:
                continue
            rows.append(
                {
                    **base,
                    "property": prop,
                    "n_active": int(len(active)),
                    "n_decoy": int(len(decoy)),
                    "active_mean": float(pd.to_numeric(active[prop], errors="coerce").mean())
                    if pd.to_numeric(active[prop], errors="coerce").notna().any()
                    else None,
                    "decoy_mean": float(pd.to_numeric(decoy[prop], errors="coerce").mean())
                    if pd.to_numeric(decoy[prop], errors="coerce").notna().any()
                    else None,
                    "standardized_mean_difference": _standardized_mean_difference(active[prop], decoy[prop]),
                }
            )
    table = pd.DataFrame(rows)
    table.to_csv(out / "decoy_property_matching.csv", index=False)
    return table


def _ligand_only_separability(df: pd.DataFrame, decoy_col: str, out: Path, seed: int) -> pd.DataFrame:
    feature_cols = [col for col in LIGAND_ONLY_COLS if col in df.columns]
    if len(feature_cols) == 0 or df[decoy_col].nunique() < 2 or len(df) < 20:
        table = pd.DataFrame(
            [{"status": "skipped", "reason": "insufficient_ligand_features_or_classes"}]
        )
        table.to_csv(out / "decoy_ligand_only_separability.csv", index=False)
        return table
    try:
        from sklearn.linear_model import LogisticRegression
        from sklearn.model_selection import train_test_split
    except ImportError:
        table = pd.DataFrame([{"status": "skipped", "reason": "sklearn_unavailable"}])
        table.to_csv(out / "decoy_ligand_only_separability.csv", index=False)
        return table
    x = df[feature_cols].copy()
    for col in x.columns:
        if x[col].dtype == object:
            x[col] = x[col].fillna("missing")
        else:
            vals = pd.to_numeric(x[col], errors="coerce")
            x[col] = vals.fillna(vals.median() if vals.notna().any() else 0.0)
    x = pd.get_dummies(x, dummy_na=True).astype(float)
    y = df[decoy_col].astype(int)
    try:
        x_train, x_test, y_train, y_test = train_test_split(
            x,
            y,
            test_size=0.25,
            random_state=seed,
            stratify=y if y.nunique() > 1 and y.value_counts().min() >= 2 else None,
        )
        model = LogisticRegression(max_iter=1000, class_weight="balanced")
        model.fit(x_train, y_train)
        probs = model.predict_proba(x_test)[:, 1]
        metrics = calibration_metrics([float(v) for v in probs], y_test.astype(int).tolist())
        table = pd.DataFrame(
            [
                {
                    "status": "ok",
                    "n_train": int(len(x_train)),
                    "n_test": int(len(x_test)),
                    "feature_count": int(x.shape[1]),
                    **metrics,
                }
            ]
        )
    except Exception as exc:
        table = pd.DataFrame([{"status": "error", "reason": str(exc)}])
    table.to_csv(out / "decoy_ligand_only_separability.csv", index=False)
    return table


def _scaffold_bias(df: pd.DataFrame, decoy_col: str, out: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if "scaffold_key" not in df.columns:
        table = pd.DataFrame([{"status": "skipped", "reason": "missing_scaffold_key"}])
        table.to_csv(out / "decoy_scaffold_bias.csv", index=False)
        return table
    group_cols = [col for col in ["target_id", "pdb_id"] if col in df.columns]
    groups = df.groupby(group_cols, dropna=False) if group_cols else [("all", df)]
    for group_key, group in groups:
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        active_scaffolds = set(group.loc[group[decoy_col].eq(0), "scaffold_key"].dropna().astype(str))
        decoy_scaffolds = set(group.loc[group[decoy_col].eq(1), "scaffold_key"].dropna().astype(str))
        overlap = active_scaffolds & decoy_scaffolds
        row = {col: value for col, value in zip(group_cols, group_key, strict=False)}
        row.update(
            {
                "status": "ok",
                "active_scaffolds": int(len(active_scaffolds)),
                "decoy_scaffolds": int(len(decoy_scaffolds)),
                "overlap_scaffolds": int(len(overlap)),
                "decoy_overlap_fraction": len(overlap) / len(decoy_scaffolds) if decoy_scaffolds else None,
                "top_active_scaffold_fraction": float(
                    group.loc[group[decoy_col].eq(0), "scaffold_key"].value_counts(normalize=True).iloc[0]
                )
                if active_scaffolds
                else None,
            }
        )
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(out / "decoy_scaffold_bias.csv", index=False)
    return table


def audit_decoy_bias(
    dataset_path: str | Path,
    out_dir: str | Path,
    *,
    seed: int = 42,
) -> dict[str, Any]:
    df = pd.read_csv(dataset_path, low_memory=False)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    work = df.copy()
    work["_decoy_audit_label"] = _decoy_labels(work)
    if "decoy_source" not in work.columns:
        source = pd.Series("unknown", index=work.index)
        for col in ["source_name", "label_source", "negative_source"]:
            if col in work.columns:
                source = source.mask(source.eq("unknown"), work[col].fillna("").astype(str))
        work["decoy_source"] = source
    n_decoy = int(work["_decoy_audit_label"].sum())
    n_active_or_nondecoy = int((work["_decoy_audit_label"] == 0).sum())
    if n_decoy == 0:
        manifest = {
            "status": "skipped",
            "reason": "no_decoy_or_benchmark_only_rows_detected",
            "n_rows": int(len(work)),
        }
        (out / "decoy_bias_summary.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
        return manifest
    property_table = _property_matching(work, "_decoy_audit_label", out)
    separability = _ligand_only_separability(work, "_decoy_audit_label", out, seed)
    scaffold = _scaffold_bias(work, "_decoy_audit_label", out)
    flags: list[str] = []
    if not separability.empty and "AUROC" in separability.columns:
        auc = pd.to_numeric(separability["AUROC"], errors="coerce").dropna()
        if not auc.empty and float(auc.iloc[0]) >= 0.65:
            flags.append("ligand_only_decoy_separability_high")
    if "standardized_mean_difference" in property_table.columns:
        smd = pd.to_numeric(property_table["standardized_mean_difference"], errors="coerce").abs()
        if smd.dropna().ge(0.5).any():
            flags.append("active_decoy_property_shift_detected")
    if "top_active_scaffold_fraction" in scaffold.columns:
        concentration = pd.to_numeric(scaffold["top_active_scaffold_fraction"], errors="coerce")
        if concentration.dropna().ge(0.5).any():
            flags.append("active_scaffold_concentration_high")
    manifest = {
        "status": "ok",
        "n_rows": int(len(work)),
        "n_decoy": n_decoy,
        "n_active_or_nondecoy": n_active_or_nondecoy,
        "flags": flags,
        "claim_policy": (
            "Decoy rows are benchmark-only. If ligand-only or property-only models "
            "separate decoys well, docking/ML benchmark claims require caveats or "
            "additional experimental inactive validation."
        ),
        "outputs": {
            "property_matching": str(out / "decoy_property_matching.csv"),
            "ligand_only_separability": str(out / "decoy_ligand_only_separability.csv"),
            "scaffold_bias": str(out / "decoy_scaffold_bias.csv"),
            "summary": str(out / "decoy_bias_summary.json"),
        },
    }
    (out / "decoy_bias_summary.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest
