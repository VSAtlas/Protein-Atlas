from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import pandas as pd


REQUIRED_COLUMNS = (
    "lig_smiles",
    "active",
    "pocket",
    "ex_rec_pdb",
)


@dataclass(frozen=True)
class BigBindDataset:
    bigbind_root: Path
    requested_splits: tuple[str, ...]
    holdout_match_column: str
    full_df: pd.DataFrame
    train_df: pd.DataFrame
    holdout_df: pd.DataFrame


def resolve_bigbind_root(bigbind_dir: Path) -> Path:
    expanded = bigbind_dir.expanduser().resolve()
    direct = expanded
    nested = expanded / "BigBindV1.5"
    if (direct / "activities_train.csv").exists():
        return direct
    if (nested / "activities_train.csv").exists():
        return nested
    raise FileNotFoundError(
        f"Could not find BigBind activities CSV files under {expanded}. "
        "Expected activities_{train,val,test}.csv in that directory "
        "or in a BigBindV1.5 subdirectory."
    )


def resolve_pocket_path(bigbind_root: Path, raw_pocket_path: str) -> Path:
    candidate = Path(raw_pocket_path)
    if candidate.is_absolute():
        return candidate
    direct = bigbind_root / candidate
    if direct.exists():
        return direct
    nested = bigbind_root / "BigBindV1.5" / candidate
    if nested.exists():
        return nested
    return direct


def _normalize_text_column(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.lower()


def _to_binary_active(value: object) -> int:
    if isinstance(value, bool):
        return int(value)
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y"}:
        return 1
    if text in {"0", "false", "f", "no", "n"}:
        return 0
    return int(float(text)) if text else 0


def _dataset_stats(df: pd.DataFrame) -> tuple[int, int, float]:
    n_rows = int(len(df))
    n_actives = int(df["active"].sum()) if n_rows else 0
    active_fraction = float(n_actives / n_rows) if n_rows else float("nan")
    return n_rows, n_actives, active_fraction


def _log_dataset_stats(logger: logging.Logger | None, label: str, df: pd.DataFrame) -> None:
    if logger is None:
        return
    n_rows, n_actives, active_fraction = _dataset_stats(df)
    logger.info(
        "[bigbind] %s: N=%d n_actives=%d actives_fraction=%.6f",
        label,
        n_rows,
        n_actives,
        active_fraction,
    )


def _load_split_csv(bigbind_root: Path, split: str) -> pd.DataFrame:
    split_name = split.strip().lower()
    split_path = bigbind_root / f"activities_{split_name}.csv"
    if not split_path.exists():
        raise FileNotFoundError(f"Missing BigBind split file: {split_path}")
    frame = pd.read_csv(split_path, low_memory=False)
    missing = [column for column in REQUIRED_COLUMNS if column not in frame.columns]
    if missing:
        raise KeyError(f"{split_path} missing required columns: {missing}")
    frame = frame.copy()
    frame["split"] = split_name
    return frame


def load_train_holdout_from_bigbind(
    *,
    bigbind_dir: Path,
    train_pdb: str,
    splits: Sequence[str],
    max_rows: int | None,
    random_seed: int,
    logger: logging.Logger | None = None,
) -> BigBindDataset:
    bigbind_root = resolve_bigbind_root(bigbind_dir)
    requested_splits = tuple(split.strip().lower() for split in splits if split.strip())
    if not requested_splits:
        raise ValueError("At least one split is required.")

    split_frames = [_load_split_csv(bigbind_root, split) for split in requested_splits]
    full_df = pd.concat(split_frames, ignore_index=True)
    full_df["active"] = full_df["active"].map(_to_binary_active).astype(int)

    if max_rows is not None and max_rows > 0 and len(full_df) > max_rows:
        full_df = full_df.sample(n=max_rows, random_state=random_seed).reset_index(drop=True)
        if logger:
            logger.info("[bigbind] Applied max_rows=%d after loading requested splits.", max_rows)

    target = train_pdb.strip().lower()
    if not target:
        raise ValueError("Config key 'train_pdb' cannot be empty.")

    pdb_mask = _normalize_text_column(full_df["ex_rec_pdb"]) == target
    if pdb_mask.any():
        holdout_mask = pdb_mask
        holdout_match_column = "ex_rec_pdb"
    else:
        pocket_mask = _normalize_text_column(full_df["pocket"]) == target
        holdout_mask = pocket_mask
        holdout_match_column = "pocket"

    holdout_df = full_df.loc[holdout_mask].copy()
    train_df = full_df.loc[~holdout_mask].copy()

    if holdout_df.empty:
        raise ValueError(
            "No holdout rows matched train_pdb. "
            f"Tried ex_rec_pdb and pocket matching for '{train_pdb}'."
        )
    if train_df.empty:
        raise ValueError(
            f"All rows matched holdout target '{train_pdb}', leaving no training rows."
        )

    _log_dataset_stats(logger, "all_rows", full_df)
    _log_dataset_stats(logger, "train_rows", train_df)
    _log_dataset_stats(logger, "holdout_rows", holdout_df)
    if logger:
        logger.info(
            "[bigbind] Holdout matching column for train_pdb='%s': %s",
            train_pdb,
            holdout_match_column,
        )

    return BigBindDataset(
        bigbind_root=bigbind_root,
        requested_splits=requested_splits,
        holdout_match_column=holdout_match_column,
        full_df=full_df,
        train_df=train_df,
        holdout_df=holdout_df,
    )
