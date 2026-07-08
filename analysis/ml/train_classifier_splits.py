from __future__ import annotations

import pandas as pd


def _custom_split(
    data: pd.DataFrame,
    split_column: str | None,
    train_values: list[str] | None,
    test_values: list[str] | None,
) -> tuple[pd.Index, pd.Index, dict[str, object]] | None:
    if not split_column:
        return None
    if split_column not in data.columns:
        raise ValueError(f"custom split column {split_column!r} not found")
    values = data[split_column].fillna("").astype(str)

    def matches(allowed: list[str] | None) -> pd.Series:
        if not allowed:
            return pd.Series(False, index=data.index)
        mask = pd.Series(False, index=data.index)
        for value in allowed:
            token = str(value).strip()
            if not token:
                continue
            mask |= values.eq(token) | values.str.split(";").map(lambda parts: token in parts)
        return mask

    test_mask = matches(test_values)
    train_mask = matches(train_values) if train_values else ~test_mask
    train_mask &= ~test_mask
    if not train_mask.any() or not test_mask.any():
        raise ValueError(
            f"custom split produced empty train/test: train={int(train_mask.sum())} test={int(test_mask.sum())}"
        )
    summary = {
        "split_mode": "custom_source_holdout",
        "split_column": split_column,
        "train_values": train_values or ["not_test_values"],
        "test_values": test_values or [],
        "n_train": int(train_mask.sum()),
        "n_test": int(test_mask.sum()),
        "passes_holdout": True,
        "overlaps": {},
    }
    return data.index[train_mask], data.index[test_mask], summary


def _temporal_split(
    data: pd.DataFrame,
    year_col: str | None,
    cutoff_year: int | None,
    *,
    max_missing_fraction: float = 0.20,
) -> tuple[pd.Index, pd.Index, dict[str, object]] | None:
    if not year_col and cutoff_year is None:
        return None
    candidates = [year_col] if year_col else [
        "activity_publication_year",
        "database_release_year",
        "label_publication_year",
        "evidence_publication_year",
    ]
    field = next((col for col in candidates if col and col in data.columns), None)
    if not field:
        raise ValueError("temporal validation requires a year column")
    years = pd.to_numeric(data[field], errors="coerce")
    if cutoff_year is None:
        cutoff = int(years.dropna().quantile(0.8))
    else:
        cutoff = int(cutoff_year)
    observed_year = years.notna()
    missing_fraction = float((~observed_year).mean()) if len(years) else 1.0
    if missing_fraction > max_missing_fraction:
        raise ValueError(
            f"temporal validation year coverage too low for {field}: "
            f"missing_fraction={missing_fraction:.3f} > {max_missing_fraction:.3f}"
        )
    train_mask = observed_year & years.le(cutoff)
    test_mask = observed_year & years.gt(cutoff)
    if not train_mask.any() or not test_mask.any():
        raise ValueError(
            f"temporal split produced empty train/test using {field} cutoff={cutoff}"
        )
    return data.index[train_mask], data.index[test_mask], {
        "split_mode": "temporal_holdout",
        "year_col": field,
        "cutoff_year": cutoff,
        "n_train": int(train_mask.sum()),
        "n_test": int(test_mask.sum()),
        "n_missing_year_excluded": int((~observed_year).sum()),
        "missing_year_fraction": missing_fraction,
        "passes_holdout": True,
        "overlaps": {},
    }
def _fixed_validation_split(
    train: pd.DataFrame,
    *,
    validation_fold_col: str | None,
    validation_fold_value: str | None,
    seed: int,
    validation_fraction: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Index]:
    if validation_fold_col:
        if validation_fold_col not in train.columns:
            raise ValueError(f"validation fold column {validation_fold_col!r} not found")
        value = "" if validation_fold_value is None else str(validation_fold_value)
        mask = train[validation_fold_col].fillna("").astype(str).eq(value)
        if not mask.any() or mask.all():
            raise ValueError(
                f"validation fold produced invalid train/validation sizes using "
                f"{validation_fold_col}={value!r}"
            )
        return train.loc[~mask].copy(), train.loc[mask].copy(), train.index[mask]
    if validation_fraction <= 0 or len(train) < 5:
        return train, train.iloc[0:0].copy(), pd.Index([])
    val_idx = train.sample(frac=validation_fraction, random_state=seed).index
    return train.drop(index=val_idx).copy(), train.loc[val_idx].copy(), val_idx


def _clip_weight(value: object, default: float = 1.0) -> float:
    try:
        parsed = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return default
    if pd.isna(parsed):
        return default
    return min(1.0, max(0.05, parsed))


def _sample_weights(frame: pd.DataFrame, label_col: str, *, pu_mode: str) -> pd.Series | None:
    weights = pd.Series(1.0, index=frame.index)
    used = False
    if pu_mode in {"positive_unlabeled_weighted", "case_control_matched"} and "_sample_weight" in frame.columns:
        weights = pd.to_numeric(frame["_sample_weight"], errors="coerce").fillna(1.0).clip(lower=0.05, upper=1.0)
        used = True
    labels = pd.to_numeric(frame[label_col], errors="coerce")
    if "negative_confidence" in frame.columns:
        confidence = frame["negative_confidence"].map(_clip_weight)
        neg_mask = labels.eq(0)
        if neg_mask.any():
            weights.loc[neg_mask] = weights.loc[neg_mask] * confidence.loc[neg_mask]
            used = True
    if "negative_evidence_type" in frame.columns:
        evidence_weight = {
            "fold_specific_reliable_negative": 0.35,
            "reliable_negative": 0.35,
            "pu_temporary_unlabeled_negative": 0.25,
            "faers_nonsignal": 0.25,
            "offsides_nonsignal": 0.25,
            "omop_negative_control": 0.75,
            "measured_inactive": 1.0,
            "measured_or_reliable_negative": 0.75,
        }
        mapped = frame["negative_evidence_type"].fillna("").astype(str).map(evidence_weight).dropna()
        if not mapped.empty:
            weights.loc[mapped.index] = weights.loc[mapped.index] * mapped
            used = True
    return weights if used else None


def _exclude_unlabeled_holdout_entities(unlabeled: pd.DataFrame, test: pd.DataFrame, split_mode: str) -> pd.DataFrame:
    if unlabeled.empty:
        return unlabeled
    fields_by_mode = {
        "drug_holdout": ["drug_id"],
        "target_holdout": ["target_id"],
        "scaffold_holdout": ["scaffold_key", "ligand_chemotype"],
        "chemical_cluster_holdout": ["chemical_cluster", "scaffold_key", "ligand_chemotype"],
        "target_family_holdout": ["target_family", "protein_class"],
        "protein_class_holdout": ["protein_class"],
    }
    out = unlabeled.copy()
    for field in fields_by_mode.get(split_mode, []):
        if field not in out.columns or field not in test.columns:
            continue
        held_out = set(test[field].dropna().astype(str))
        if held_out:
            return out[~out[field].fillna("").astype(str).isin(held_out)].copy()
    return out
