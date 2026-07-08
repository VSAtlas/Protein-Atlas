from __future__ import annotations

import argparse
import logging
import re
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from config.output_paths import output_root, run_output_dir

COMPONENT = "[interactions-parquet]"
CHUNK_SIZE = 200_000
PART_FILE_RE = re.compile(r"^part-(\d+)\.parquet$")
TRUE_VALUES = {"1", "true", "yes", "y", "on", "t"}
FALSE_VALUES = {"0", "false", "no", "n", "off", "f"}

REQUIRED_COLUMNS = [
    "run_id",
    "target_id",
    "pdb_id",
    "variant",
    "ph_label",
    "ligand_base",
    "ligand_display",
    "z_selected",
]
STRING_COLUMNS = [
    "z_selected_source",
    "library",
    "run_mode",
    "ligand_file",
    "pose_invalid_reason_top",
    "pocket_method",
]
BOOL_COLUMNS = [
    "pose_valid_any",
    "is_decoy",
    "is_control",
]
INT_COLUMNS = ["rank"]
FLOAT32_COLUMNS = [
    "pct_rank",
    "center_x",
    "center_y",
    "center_z",
    "box_x",
    "box_y",
    "box_z",
    "ef1",
    "roc_auc",
    "roc_auc_adj",
]


def _configure_logging(verbose: bool) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")
    return logging.getLogger("interactions-parquet")


def _normalize_text_series(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _coerce_bool_series(series: pd.Series) -> pd.Series:
    normalized = _normalize_text_series(series).str.lower()
    out = pd.Series(pd.NA, index=series.index, dtype="boolean")
    out.loc[normalized.isin(TRUE_VALUES)] = True
    out.loc[normalized.isin(FALSE_VALUES)] = False
    return out


def _coerce_float_series(series: pd.Series, dtype: str) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    finite_mask = np.isfinite(numeric.to_numpy(dtype="float64", na_value=np.nan))
    numeric = numeric.where(finite_mask, np.nan)
    return numeric.astype(dtype)


def _coerce_int_series(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    values = numeric.to_numpy(dtype="float64", na_value=np.nan)
    finite_mask = np.isfinite(values)
    integer_mask = finite_mask & np.isclose(values, np.round(values))
    numeric = numeric.where(integer_mask, np.nan)
    return numeric.astype("Int32")


def _split_target_id(target_id: pd.Series) -> Tuple[pd.Series, pd.Series, pd.Series]:
    parts = target_id.str.split("|", n=2, expand=True, regex=False)
    for idx in range(3):
        if idx not in parts.columns:
            parts[idx] = ""
    pdb = _normalize_text_series(parts[0])
    variant = _normalize_text_series(parts[1])
    ph = _normalize_text_series(parts[2])
    return pdb, variant, ph


def _parse_targets(raw: Optional[str]) -> Optional[Set[str]]:
    if raw is None:
        return None
    parsed = {token.strip() for token in raw.split(",") if token.strip()}
    return parsed or None


def _parse_partition_cols(raw: str) -> List[str]:
    cols = [token.strip() for token in str(raw or "").split(",") if token.strip()]
    if not cols:
        raise ValueError("partition columns cannot be empty")
    return cols


def _resolve_input_path(repo_root: Path, run_id: str, raw_input: Optional[str]) -> Path:
    if raw_input:
        input_path = Path(raw_input)
        if not input_path.is_absolute():
            input_path = repo_root / input_path
        return input_path
    candidates = [
        output_root(repo_root, "data") / run_id / "master_rows.csv",
        output_root(repo_root, "data") / run_id / "heatmap_input.csv",
        repo_root / "data" / run_id / "master_rows.csv",
        repo_root / "data" / run_id / "heatmap_input.csv",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "No default input CSV found; expected outputs/data/<run_id>/master_rows.csv "
        "or outputs/data/<run_id>/heatmap_input.csv"
    )


def _resolve_output_dir(repo_root: Path, run_id: str, raw_out_dir: Optional[str]) -> Path:
    if raw_out_dir:
        out_dir = Path(raw_out_dir)
        if not out_dir.is_absolute():
            out_dir = repo_root / out_dir
        return out_dir
    data_dir = run_output_dir(repo_root, "data", run_id)
    return data_dir / "dataset" / "interactions"


def _validate_input_columns(columns: Sequence[str]) -> None:
    colset = set(columns)
    if "z_selected" not in colset:
        raise ValueError("Missing required column: z_selected")
    has_target_id = "target_id" in colset
    has_target_components = {"pdb_id", "variant", "ph_label"}.issubset(colset)
    if not has_target_id and not has_target_components:
        raise ValueError(
            "Missing required target columns: expected target_id or "
            "all of pdb_id, variant, ph_label"
        )
    if "ligand_display" not in colset and "ligand_base" not in colset:
        raise ValueError("Missing required columns: ligand_display or ligand_base")


def _normalize_chunk(
    chunk: pd.DataFrame, run_id: str, targets: Optional[Set[str]]
) -> pd.DataFrame:
    out = pd.DataFrame(index=chunk.index)

    z_selected = _coerce_float_series(chunk["z_selected"], "float64")
    valid_t = z_selected.notna()
    if not bool(valid_t.any()):
        return out

    chunk = chunk.loc[valid_t].copy()
    out = out.loc[valid_t].copy()
    z_selected = z_selected.loc[valid_t]
    out["z_selected"] = z_selected

    if "target_id" in chunk.columns:
        target_id = _normalize_text_series(chunk["target_id"])
    else:
        target_id = pd.Series("", index=chunk.index, dtype="string")

    has_target_components = {"pdb_id", "variant", "ph_label"}.issubset(chunk.columns)
    if has_target_components:
        pdb_raw = _normalize_text_series(chunk["pdb_id"])
        variant_raw = _normalize_text_series(chunk["variant"])
        ph_raw = _normalize_text_series(chunk["ph_label"])
        derived_target_id = pdb_raw + "|" + variant_raw + "|" + ph_raw
        target_id = target_id.where(target_id != "", derived_target_id)

    target_id = target_id.replace("nan", "", regex=False)
    pdb_from_target, variant_from_target, ph_from_target = _split_target_id(target_id)

    if "pdb_id" in chunk.columns:
        pdb_id = _normalize_text_series(chunk["pdb_id"]).where(
            _normalize_text_series(chunk["pdb_id"]) != "", pdb_from_target
        )
    else:
        pdb_id = pdb_from_target

    if "variant" in chunk.columns:
        variant = _normalize_text_series(chunk["variant"]).where(
            _normalize_text_series(chunk["variant"]) != "", variant_from_target
        )
    else:
        variant = variant_from_target

    if "ph_label" in chunk.columns:
        ph_label = _normalize_text_series(chunk["ph_label"]).where(
            _normalize_text_series(chunk["ph_label"]) != "", ph_from_target
        )
    else:
        ph_label = ph_from_target

    target_id = pdb_id + "|" + variant + "|" + ph_label
    target_id_clean = target_id.str.strip()
    valid_target = (target_id_clean != "") & (target_id_clean != "||")
    if not bool(valid_target.any()):
        return pd.DataFrame(columns=REQUIRED_COLUMNS)

    out = out.loc[valid_target].copy()
    chunk = chunk.loc[valid_target].copy()
    target_id = target_id.loc[valid_target]
    pdb_id = pdb_id.loc[valid_target]
    variant = variant.loc[valid_target]
    ph_label = ph_label.loc[valid_target]

    if "ligand_display" in chunk.columns:
        ligand_display = _normalize_text_series(chunk["ligand_display"])
    else:
        ligand_display = pd.Series("", index=chunk.index, dtype="string")
    if "ligand_base" in chunk.columns:
        ligand_base = _normalize_text_series(chunk["ligand_base"])
    else:
        ligand_base = pd.Series("", index=chunk.index, dtype="string")

    ligand_display = ligand_display.mask(
        ligand_display.str.lower() == "nan", ""
    )
    ligand_base = ligand_base.mask(ligand_base.str.lower() == "nan", "")
    ligand_display = ligand_display.where(ligand_display != "", ligand_base)
    ligand_base = ligand_base.where(ligand_base != "", ligand_display)

    valid_ligand = (ligand_display.str.strip() != "") | (ligand_base.str.strip() != "")
    if not bool(valid_ligand.any()):
        return pd.DataFrame(columns=REQUIRED_COLUMNS)

    out = out.loc[valid_ligand].copy()
    chunk = chunk.loc[valid_ligand].copy()
    target_id = target_id.loc[valid_ligand]
    pdb_id = pdb_id.loc[valid_ligand]
    variant = variant.loc[valid_ligand]
    ph_label = ph_label.loc[valid_ligand]
    ligand_display = ligand_display.loc[valid_ligand]
    ligand_base = ligand_base.loc[valid_ligand]

    out["run_id"] = run_id
    out["target_id"] = target_id
    out["pdb_id"] = pdb_id
    out["variant"] = variant
    out["ph_label"] = ph_label
    out["ligand_base"] = ligand_base
    out["ligand_display"] = ligand_display

    for col in STRING_COLUMNS:
        if col in chunk.columns:
            out[col] = _normalize_text_series(chunk[col])

    for col in BOOL_COLUMNS:
        if col in chunk.columns:
            out[col] = _coerce_bool_series(chunk[col])

    for col in INT_COLUMNS:
        if col in chunk.columns:
            out[col] = _coerce_int_series(chunk[col])

    for col in FLOAT32_COLUMNS:
        if col in chunk.columns:
            out[col] = _coerce_float_series(chunk[col], "float32")

    fdr_cols = sorted(col for col in chunk.columns if col.startswith("fdr_"))
    for col in fdr_cols:
        if col in out.columns:
            continue
        if col == "fdr_reliable" or col.startswith("fdr_hit_"):
            out[col] = _coerce_bool_series(chunk[col])
        elif col.endswith("_field"):
            out[col] = _normalize_text_series(chunk[col])
        else:
            out[col] = _coerce_float_series(chunk[col], "float32")

    if targets is not None:
        out = out[out["target_id"].isin(targets)].copy()

    ordered = REQUIRED_COLUMNS + [
        col
        for col in (
            STRING_COLUMNS
            + BOOL_COLUMNS
            + INT_COLUMNS
            + FLOAT32_COLUMNS
            + fdr_cols
        )
        if col in out.columns and col not in REQUIRED_COLUMNS
    ]
    return out[ordered]


def _init_part_counters(
    out_dir: Path, partition_cols: Sequence[str]
) -> Dict[Tuple[str, ...], int]:
    counters: Dict[Tuple[str, ...], int] = {}
    if not out_dir.exists():
        return counters
    for part_path in out_dir.rglob("part-*.parquet"):
        match = PART_FILE_RE.match(part_path.name)
        if match is None:
            continue
        idx = int(match.group(1))
        rel_parent = part_path.parent.relative_to(out_dir)
        values: Dict[str, str] = {}
        for segment in rel_parent.parts:
            if "=" not in segment:
                continue
            col, value = segment.split("=", 1)
            values[col] = value
        if not all(col in values for col in partition_cols):
            continue
        key = tuple(values[col] for col in partition_cols)
        counters[key] = max(counters.get(key, 0), idx + 1)
    return counters


def _partition_path(out_dir: Path, partition_cols: Sequence[str], key: Sequence[str]) -> Path:
    partition_dir = out_dir
    for col, value in zip(partition_cols, key):
        partition_dir = partition_dir / f"{col}={value}"
    return partition_dir


def _write_chunk(
    chunk: pd.DataFrame,
    out_dir: Path,
    partition_cols: Sequence[str],
    counters: Dict[Tuple[str, ...], int],
) -> Set[str]:
    targets_written: Set[str] = set()
    grouped = chunk.groupby(list(partition_cols), sort=False, dropna=False)
    for raw_key, partition_df in grouped:
        if not isinstance(raw_key, tuple):
            raw_key = (raw_key,)
        key = tuple(str(value) for value in raw_key)
        partition_dir = _partition_path(out_dir, partition_cols, key)
        partition_dir.mkdir(parents=True, exist_ok=True)
        part_idx = counters.get(key, 0)
        out_path = partition_dir / f"part-{part_idx:05d}.parquet"
        table = pa.Table.from_pandas(partition_df, preserve_index=False)
        pq.write_table(table, out_path, compression="snappy")
        counters[key] = part_idx + 1
        targets_written.update(partition_df["target_id"].dropna().astype(str).tolist())
    return targets_written


def convert_interactions_to_parquet(
    run_id: str,
    repo_root: Path,
    input_path: Optional[Path] = None,
    out_dir: Optional[Path] = None,
    partition_cols: Optional[Sequence[str]] = None,
    overwrite: bool = False,
    targets: Optional[Set[str]] = None,
    logger: Optional[logging.Logger] = None,
) -> Dict[str, Any]:
    log = logger or logging.getLogger("interactions-parquet")
    resolved_input = input_path or _resolve_input_path(repo_root, run_id, None)
    resolved_out_dir = out_dir or _resolve_output_dir(repo_root, run_id, None)
    resolved_partition_cols = list(partition_cols or ["run_id", "target_id"])

    if not resolved_input.exists():
        raise FileNotFoundError(f"Input CSV not found: {resolved_input}")
    if resolved_input.suffix.lower() != ".csv":
        raise ValueError(f"Expected CSV input, got: {resolved_input}")

    header_df = pd.read_csv(resolved_input, nrows=0, dtype=str)
    _validate_input_columns(list(header_df.columns))

    if overwrite and resolved_out_dir.exists():
        shutil.rmtree(resolved_out_dir)
    resolved_out_dir.mkdir(parents=True, exist_ok=True)

    counters = _init_part_counters(resolved_out_dir, resolved_partition_cols)
    total_rows_in = 0
    total_rows_out = 0
    written_targets: Set[str] = set()

    for chunk in pd.read_csv(
        resolved_input,
        dtype=str,
        chunksize=CHUNK_SIZE,
        keep_default_na=False,
    ):
        total_rows_in += len(chunk)
        normalized = _normalize_chunk(chunk, run_id=run_id, targets=targets)
        if normalized.empty:
            continue
        missing_partition_cols = [
            col for col in resolved_partition_cols if col not in normalized.columns
        ]
        if missing_partition_cols:
            raise ValueError(
                "Missing partition columns in normalized dataset: "
                + ", ".join(missing_partition_cols)
            )
        total_rows_out += len(normalized)
        written_targets.update(
            _write_chunk(
                normalized,
                resolved_out_dir,
                resolved_partition_cols,
                counters,
            )
        )

    stats = {
        "input_path": resolved_input,
        "output_dir": resolved_out_dir,
        "rows_in": total_rows_in,
        "rows_out": total_rows_out,
        "targets_written": len(written_targets),
        "partition_cols": resolved_partition_cols,
    }
    log.info(
        (
            "%s action=convert status=ok input=%s rows_in=%d rows_out=%d "
            "targets_written=%d out_dir=%s partition_cols=%s"
        ),
        COMPONENT,
        resolved_input,
        total_rows_in,
        total_rows_out,
        len(written_targets),
        resolved_out_dir,
        ",".join(resolved_partition_cols),
    )
    return stats


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert master_rows.csv/heatmap_input.csv into a partitioned "
            "Parquet interactions dataset."
        )
    )
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--input", default=None)
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--partition-cols", default="run_id,target_id")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--targets",
        default=None,
        help="Comma-separated target_id values to build.",
    )
    parser.add_argument("--verbose", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _build_arg_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    logger = _configure_logging(args.verbose)
    repo_root = Path(args.repo_root).resolve()
    partition_cols = _parse_partition_cols(args.partition_cols)
    targets = _parse_targets(args.targets)
    input_path = _resolve_input_path(repo_root, args.run_id, args.input)
    out_dir = _resolve_output_dir(repo_root, args.run_id, args.out_dir)

    logger.info(
        "%s action=start run_id=%s input=%s out_dir=%s targets=%s",
        COMPONENT,
        args.run_id,
        input_path,
        out_dir,
        ",".join(sorted(targets)) if targets else "*",
    )
    try:
        convert_interactions_to_parquet(
            run_id=args.run_id,
            repo_root=repo_root,
            input_path=input_path,
            out_dir=out_dir,
            partition_cols=partition_cols,
            overwrite=args.overwrite,
            targets=targets,
            logger=logger,
        )
    except Exception as exc:
        logger.error("%s action=fail error=%s", COMPONENT, exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
