from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd


LABEL_CANDIDATES = ("Y", "Label", "label", "activity_outcome", "Outcome", "score", "Score")
SMILES_CANDIDATES = ("Drug", "Drug1", "SMILES", "smiles", "canonical_smiles", "Drug_SMILES")
DRUG_ID_CANDIDATES = ("Drug_ID", "Drug_ID1", "drug_id", "compound_id", "mol_id", "Name")
DRUG_NAME_CANDIDATES = ("Drug_Name", "drug_name", "Name", "Compound", "compound_name")
TARGET_ID_CANDIDATES = ("Target_ID", "target_id", "Protein_ID", "Gene", "Target_Name", "target_name")
TARGET_SEQUENCE_CANDIDATES = ("Target", "Target_Sequence", "target_sequence", "Protein", "protein_sequence")


def _safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")
    return token or "tdc_dataset"


def _stable_id(prefix: str, value: Any) -> str:
    digest = hashlib.sha1(str(value).encode("utf-8")).hexdigest()[:16]
    return f"{prefix}:{digest}"


def _separator(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def _resolve_column(
    frame: pd.DataFrame,
    requested: str | None,
    candidates: tuple[str, ...],
    *,
    role: str,
    required: bool = True,
) -> str | None:
    if requested:
        if requested not in frame.columns:
            raise ValueError(f"{role} column {requested!r} not found")
        return requested
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate
    if required:
        raise ValueError(f"could not infer {role} column; pass --{role}-col")
    return None


def _frame_from_tdc_data(data: Any) -> pd.DataFrame:
    if isinstance(data, pd.DataFrame):
        return data.copy()
    if isinstance(data, dict):
        frames: list[pd.DataFrame] = []
        for split_name, split_data in data.items():
            split_frame = pd.DataFrame(split_data).copy()
            split_frame["tdc_split"] = str(split_name)
            frames.append(split_frame)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return pd.DataFrame(data)


def _load_tdc_frame(
    *,
    tdc_module: str,
    tdc_class: str,
    name: str,
    cache_dir: Path | None,
    use_tdc_split: bool,
    split_method: str,
) -> pd.DataFrame:
    try:
        module = importlib.import_module(f"tdc.{tdc_module}")
    except ImportError as exc:
        raise RuntimeError(
            "PyTDC is not installed in this environment. Install the optional molecular tools "
            "or pass --input-csv with an already downloaded dataset."
        ) from exc
    try:
        dataset_cls = getattr(module, tdc_class)
    except AttributeError as exc:
        raise ValueError(f"tdc.{tdc_module} has no dataset class {tdc_class!r}") from exc
    kwargs: dict[str, Any] = {"name": name}
    if cache_dir is not None:
        kwargs["path"] = str(cache_dir)
    try:
        dataset = dataset_cls(**kwargs)
    except TypeError:
        if "path" not in kwargs:
            raise
        kwargs.pop("path")
        dataset = dataset_cls(**kwargs)
    if use_tdc_split:
        try:
            split = dataset.get_split(method=split_method)
        except TypeError:
            split = dataset.get_split()
        return _frame_from_tdc_data(split)
    return _frame_from_tdc_data(dataset.get_data())


def _coerce_label(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    if int(numeric.notna().sum()) == int(values.notna().sum()):
        return numeric
    return values.astype("string")


def _series_or_missing(frame: pd.DataFrame, col: str | None) -> pd.Series:
    if col and col in frame.columns:
        return frame[col]
    return pd.Series(pd.NA, index=frame.index)


def stage_tdc_frame(
    frame: pd.DataFrame,
    out_path: str | Path,
    *,
    dataset_name: str,
    tdc_module: str,
    tdc_class: str,
    input_source: str,
    label_col: str | None = None,
    smiles_col: str | None = None,
    drug_id_col: str | None = None,
    drug_name_col: str | None = None,
    target_id_col: str | None = None,
    target_sequence_col: str | None = None,
    single_target_id: str | None = None,
    activity_type: str | None = None,
    activity_units: str = "",
    training_allowed: bool = False,
) -> dict[str, Any]:
    if frame.empty:
        raise ValueError("TDC input is empty")
    resolved_label = _resolve_column(frame, label_col, LABEL_CANDIDATES, role="label")
    resolved_smiles = _resolve_column(frame, smiles_col, SMILES_CANDIDATES, role="smiles")
    resolved_drug_id = _resolve_column(
        frame,
        drug_id_col,
        DRUG_ID_CANDIDATES,
        role="drug-id",
        required=False,
    )
    resolved_drug_name = _resolve_column(
        frame,
        drug_name_col,
        DRUG_NAME_CANDIDATES,
        role="drug-name",
        required=False,
    )
    resolved_target_id = _resolve_column(
        frame,
        target_id_col,
        TARGET_ID_CANDIDATES,
        role="target-id",
        required=False,
    )
    resolved_target_sequence = _resolve_column(
        frame,
        target_sequence_col,
        TARGET_SEQUENCE_CANDIDATES,
        role="target-sequence",
        required=False,
    )
    source_name = f"TDC:{dataset_name}"
    smiles = frame[resolved_smiles].astype(str).str.strip()
    labels = _coerce_label(frame[resolved_label])
    staged = pd.DataFrame(
        {
            "drug_id": _series_or_missing(frame, resolved_drug_id),
            "drug_name": _series_or_missing(frame, resolved_drug_name),
            "smiles": smiles,
            "target_id": _series_or_missing(frame, resolved_target_id),
            "target_sequence": _series_or_missing(frame, resolved_target_sequence),
            "activity_outcome": labels,
            "activity_type": activity_type or f"TDC {tdc_class} {dataset_name}",
            "activity_relation": "",
            "activity_units": activity_units,
            "assay_id": dataset_name,
            "source": source_name,
            "label_source": source_name,
            "source_family": "TDC",
            "upstream_source": source_name,
            "tdc_module": tdc_module,
            "tdc_class": tdc_class,
            "tdc_dataset": dataset_name,
            "benchmark_only": not training_allowed,
            "training_allowed": training_allowed,
        }
    )
    if "tdc_split" in frame.columns:
        staged["tdc_split"] = frame["tdc_split"]
    staged = staged[staged["smiles"].ne("") & staged["smiles"].ne("nan")]
    staged = staged[staged["activity_outcome"].notna()].copy()
    if staged.empty:
        raise ValueError("no rows remain after dropping missing SMILES and labels")
    missing_drug = staged["drug_id"].isna() | staged["drug_id"].astype(str).str.strip().eq("")
    staged.loc[missing_drug, "drug_id"] = staged.loc[missing_drug, "smiles"].map(
        lambda value: _stable_id("tdc_smiles", value)
    )
    missing_name = staged["drug_name"].isna() | staged["drug_name"].astype(str).str.strip().eq("")
    staged.loc[missing_name, "drug_name"] = staged.loc[missing_name, "drug_id"]
    missing_target = staged["target_id"].isna() | staged["target_id"].astype(str).str.strip().eq("")
    if single_target_id:
        staged.loc[missing_target, "target_id"] = single_target_id
    elif resolved_target_sequence:
        staged.loc[missing_target, "target_id"] = staged.loc[missing_target, "target_sequence"].map(
            lambda value: _stable_id("tdc_target", value)
        )
    else:
        staged.loc[missing_target, "target_id"] = f"tdc:{_safe_token(dataset_name)}"
    staged = staged.drop_duplicates().reset_index(drop=True)

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    sep = _separator(out)
    staged.to_csv(out, sep=sep, index=False)
    manifest_path = out.with_suffix(".manifest.json")
    manifest: dict[str, Any] = {
        "adapter": "tdc",
        "input": input_source,
        "output": str(out),
        "rows": int(len(staged)),
        "input_rows": int(len(frame)),
        "dataset": dataset_name,
        "tdc_module": tdc_module,
        "tdc_class": tdc_class,
        "columns": {
            "label": resolved_label,
            "smiles": resolved_smiles,
            "drug_id": resolved_drug_id,
            "drug_name": resolved_drug_name,
            "target_id": resolved_target_id,
            "target_sequence": resolved_target_sequence,
        },
        "label_counts": {
            str(key): int(value)
            for key, value in staged["activity_outcome"].value_counts(dropna=False).items()
        },
        "benchmark_only": not training_allowed,
        "training_allowed": training_allowed,
        "policy": (
            "TDC rows are staged as external benchmark evidence by default. Promote to training "
            "only after source lineage, licensing, overlap, leakage, and split-policy review."
        ),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Stage a PyTDC or local TDC-style dataset for Atlas benchmarks.")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--input-csv", type=Path, default=None)
    source.add_argument("--name", default=None, help="PyTDC dataset name, for example hERG.")
    parser.add_argument("--dataset-name", default=None, help="Dataset label for --input-csv staging.")
    parser.add_argument("--tdc-module", default="single_pred", help="PyTDC module under tdc.*, e.g. single_pred.")
    parser.add_argument("--tdc-class", default="Tox", help="PyTDC dataset class, e.g. Tox, ADME, DTI.")
    parser.add_argument("--tdc-cache-dir", type=Path, default=None)
    parser.add_argument("--use-tdc-split", action="store_true")
    parser.add_argument("--split-method", default="random")
    parser.add_argument("--label-col", default=None)
    parser.add_argument("--smiles-col", default=None)
    parser.add_argument("--drug-id-col", default=None)
    parser.add_argument("--drug-name-col", default=None)
    parser.add_argument("--target-id-col", default=None)
    parser.add_argument("--target-sequence-col", default=None)
    parser.add_argument("--single-target-id", default=None)
    parser.add_argument("--activity-type", default=None)
    parser.add_argument("--activity-units", default="")
    parser.add_argument(
        "--training-allowed",
        action="store_true",
        help="Mark staged rows as training-eligible after external source review; default is benchmark-only.",
    )
    parser.add_argument("--out", required=True, type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.input_csv is not None:
        frame = pd.read_csv(args.input_csv, low_memory=False)
        dataset_name = args.dataset_name or args.input_csv.stem
        input_source = str(args.input_csv)
    else:
        dataset_name = str(args.name)
        frame = _load_tdc_frame(
            tdc_module=args.tdc_module,
            tdc_class=args.tdc_class,
            name=dataset_name,
            cache_dir=args.tdc_cache_dir,
            use_tdc_split=args.use_tdc_split,
            split_method=args.split_method,
        )
        input_source = f"PyTDC:{args.tdc_module}.{args.tdc_class}:{dataset_name}"
    stage_tdc_frame(
        frame,
        args.out,
        dataset_name=dataset_name,
        tdc_module=args.tdc_module,
        tdc_class=args.tdc_class,
        input_source=input_source,
        label_col=args.label_col,
        smiles_col=args.smiles_col,
        drug_id_col=args.drug_id_col,
        drug_name_col=args.drug_name_col,
        target_id_col=args.target_id_col,
        target_sequence_col=args.target_sequence_col,
        single_target_id=args.single_target_id,
        activity_type=args.activity_type,
        activity_units=args.activity_units,
        training_allowed=args.training_allowed,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
