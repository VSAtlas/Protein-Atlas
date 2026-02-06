from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors


_PROPERTY_COLUMNS = ("MW", "LogP", "TPSA", "HBD", "HBA")


def _smiles_to_mol(smiles: str) -> Chem.Mol | None:
    text = str(smiles or "").strip()
    if not text:
        return None
    return Chem.MolFromSmiles(text)


def _compute_properties_for_df(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float]] = []
    for smiles in df.get("lig_smiles", pd.Series(dtype=object)).tolist():
        mol = _smiles_to_mol(str(smiles))
        if mol is None:
            rows.append({k: float("nan") for k in _PROPERTY_COLUMNS})
            continue
        rows.append(
            {
                "MW": float(Descriptors.MolWt(mol)),
                "LogP": float(Crippen.MolLogP(mol)),
                "TPSA": float(rdMolDescriptors.CalcTPSA(mol)),
                "HBD": float(Lipinski.NumHDonors(mol)),
                "HBA": float(Lipinski.NumHAcceptors(mol)),
            }
        )
    return pd.DataFrame(rows, index=df.index)


def compute_property_bins(
    df: pd.DataFrame,
    cols: list[str] | tuple[str, ...] | None = None,
    n_bins: int = 10,
) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    cols_used = list(cols or _PROPERTY_COLUMNS)
    out = df.copy()
    props = _compute_properties_for_df(out)
    for col in cols_used:
        if col not in props.columns:
            continue
        series = props[col]
        out[f"{col}_bin"] = pd.qcut(
            series.rank(method="first"),
            q=max(2, int(n_bins)),
            labels=False,
            duplicates="drop",
        )
    return out


def _group_value(row: pd.Series, *, within_group: bool, group_key: str) -> str:
    if not within_group:
        return "__global__"
    key = str(group_key).strip().lower()
    if key == "target":
        return str(row.get("ex_rec_pdb", "")).strip()
    if key == "pocket":
        return str(row.get("pocket", "")).strip()
    if key == "target_pocket":
        return f"{str(row.get('ex_rec_pdb', '')).strip()}::{str(row.get('pocket', '')).strip()}"
    return "__global__"


def mine_property_matched_negatives(
    df: pd.DataFrame,
    per_active_k: int,
    seed: int,
    *,
    within_group: bool = True,
    group_key: str = "target",
) -> pd.DataFrame:
    if df.empty:
        return df.iloc[0:0].copy()
    working = compute_property_bins(df, n_bins=10).copy()
    rng = np.random.default_rng(int(seed))
    negatives = working.loc[working["active"].astype(int) <= 0].copy()
    actives = working.loc[working["active"].astype(int) > 0].copy()

    selected_idx: list[int] = []
    bin_cols = [c for c in working.columns if c.endswith("_bin")]
    if not bin_cols:
        return negatives.iloc[0:0].copy()

    for active_idx, active_row in actives.iterrows():
        group_val = _group_value(active_row, within_group=within_group, group_key=group_key)
        mask = np.ones(len(negatives), dtype=bool)
        if within_group:
            neg_groups = negatives.apply(
                lambda row: _group_value(row, within_group=True, group_key=group_key), axis=1
            )
            mask &= neg_groups.to_numpy(dtype=object) == group_val
        for col in bin_cols:
            mask &= negatives[col].to_numpy() == active_row[col]
        candidates = negatives.loc[mask]
        if candidates.empty:
            continue
        choose_n = min(int(per_active_k), len(candidates))
        picked = rng.choice(candidates.index.to_numpy(dtype=int), size=choose_n, replace=False)
        selected_idx.extend(int(x) for x in picked.tolist())
    if not selected_idx:
        return negatives.iloc[0:0].copy()
    return negatives.loc[sorted(set(selected_idx))].copy()


def _compute_morgan_fps(smiles: pd.Series, n_bits: int = 2048) -> dict[int, Any]:
    fps: dict[int, Any] = {}
    for idx, smi in smiles.items():
        mol = _smiles_to_mol(str(smi))
        if mol is None:
            continue
        fps[int(idx)] = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=int(n_bits))
    return fps


def mine_similarity_negatives(
    df: pd.DataFrame,
    per_active_k: int,
    seed: int,
    *,
    within_group: bool = True,
    group_key: str = "target",
    max_candidates: int = 2000,
) -> pd.DataFrame:
    if df.empty:
        return df.iloc[0:0].copy()

    rng = np.random.default_rng(int(seed))
    work = df.copy()
    active_mask = work["active"].astype(int) > 0
    negatives = work.loc[~active_mask].copy()
    actives = work.loc[active_mask].copy()
    fps = _compute_morgan_fps(work["lig_smiles"])

    selected: list[int] = []
    for active_idx, active_row in actives.iterrows():
        fp_a = fps.get(int(active_idx))
        if fp_a is None:
            continue

        if within_group:
            group_val = _group_value(active_row, within_group=True, group_key=group_key)
            cand = negatives.loc[
                negatives.apply(
                    lambda row: _group_value(row, within_group=True, group_key=group_key) == group_val,
                    axis=1,
                )
            ].copy()
        else:
            cand = negatives.copy()
        if cand.empty:
            continue
        cand = cand.loc[
            cand.get("murcko_scaffold", "").astype(str)
            != str(active_row.get("murcko_scaffold", ""))
        ].copy()
        if cand.empty:
            continue

        if len(cand) > int(max_candidates):
            sampled_idx = rng.choice(
                cand.index.to_numpy(dtype=int),
                size=int(max_candidates),
                replace=False,
            )
            cand = cand.loc[sampled_idx]

        sims: list[tuple[int, float]] = []
        for neg_idx in cand.index.tolist():
            fp_n = fps.get(int(neg_idx))
            if fp_n is None:
                continue
            sim = float(DataStructs.TanimotoSimilarity(fp_a, fp_n))
            sims.append((int(neg_idx), sim))
        sims.sort(key=lambda x: (x[1], -x[0]), reverse=True)
        selected.extend(idx for idx, _ in sims[: max(0, int(per_active_k))])

    if not selected:
        return negatives.iloc[0:0].copy()
    return negatives.loc[sorted(set(selected))].copy()


def merge_hard_negatives(
    df: pd.DataFrame,
    *,
    policy: str = "property",
    ratio: float = 1.0,
    per_active_k: int = 5,
    within_group: bool = True,
    group_key: str = "target",
    max_candidates: int = 2000,
    seed: int = 42,
) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    policy_norm = str(policy).strip().lower()
    if ratio <= 0:
        return df.copy()

    if policy_norm == "property":
        mined = mine_property_matched_negatives(
            df,
            per_active_k=per_active_k,
            seed=seed,
            within_group=within_group,
            group_key=group_key,
        )
    elif policy_norm == "similarity":
        mined = mine_similarity_negatives(
            df,
            per_active_k=per_active_k,
            seed=seed,
            within_group=within_group,
            group_key=group_key,
            max_candidates=max_candidates,
        )
    elif policy_norm == "hybrid":
        p = mine_property_matched_negatives(
            df,
            per_active_k=per_active_k,
            seed=seed,
            within_group=within_group,
            group_key=group_key,
        )
        s = mine_similarity_negatives(
            df,
            per_active_k=per_active_k,
            seed=seed,
            within_group=within_group,
            group_key=group_key,
            max_candidates=max_candidates,
        )
        mined = pd.concat([p, s], ignore_index=False).drop_duplicates()
    else:
        raise ValueError("policy must be one of: property, similarity, hybrid")

    if mined.empty:
        return df.copy()

    max_add = max(1, int(round(float(ratio) * len(df))))
    mined = mined.head(max_add)
    out = pd.concat([df, mined], ignore_index=True)
    out["_hard_negative_augmented"] = 0
    out.loc[len(df) :, "_hard_negative_augmented"] = 1
    return out
