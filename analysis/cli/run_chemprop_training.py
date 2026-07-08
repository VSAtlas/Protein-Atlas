from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd


SMILES_CANDIDATES = (
    "canonical_smiles",
    "smiles",
    "SMILES",
    "standardized_smiles",
    "Drug",
)

DEFAULT_METADATA_COLS = [
    "drug_id",
    "target_id",
    "pdb_id",
    "scaffold_key",
    "chemical_cluster",
    "ligand_chemotype",
    "target_family",
    "protein_class",
    "label_source",
    "source_family",
    "upstream_source",
    "assay_type",
    "endpoint_type",
    "activity_publication_year",
    "database_release_year",
]


def _resolve_column(
    frame: pd.DataFrame,
    requested: str | None,
    candidates: tuple[str, ...],
    *,
    role: str,
) -> str:
    if requested:
        if requested not in frame.columns:
            raise ValueError(f"{role} column {requested!r} not found")
        return requested
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate
    raise ValueError(f"could not infer {role} column; pass --{role}-col")


def _infer_task_type(labels: pd.Series, requested: str) -> str:
    if requested != "auto":
        return requested
    numeric = pd.to_numeric(labels, errors="coerce").dropna()
    observed = set(numeric.astype(float).unique().tolist())
    if observed and observed.issubset({0.0, 1.0}):
        return "classification"
    return "regression"


def _prepare_label(labels: pd.Series, task_type: str) -> pd.Series:
    numeric = pd.to_numeric(labels, errors="coerce")
    if task_type == "classification":
        invalid = sorted(set(numeric.dropna().astype(float).unique().tolist()) - {0.0, 1.0})
        if invalid:
            raise ValueError(
                "Chemprop classification labels must be binary 0/1 after numeric coercion; "
                f"found non-binary values {invalid[:5]}"
            )
        return numeric.astype(int)
    return numeric


def _selected_metadata_cols(
    frame: pd.DataFrame,
    requested: list[str] | None,
    *,
    smiles_col: str,
    label_col: str,
) -> list[str]:
    candidates = requested if requested is not None else DEFAULT_METADATA_COLS
    seen: set[str] = set()
    selected: list[str] = []
    for col in candidates:
        if col in {smiles_col, label_col} or col in seen or col not in frame.columns:
            continue
        selected.append(col)
        seen.add(col)
    return selected


def _label_counts(values: pd.Series) -> dict[str, int]:
    return {str(key): int(value) for key, value in values.value_counts(dropna=False).items()}


def stage_chemprop_input(
    dataset_path: str | Path,
    label_col: str,
    out_dir: str | Path,
    *,
    smiles_col: str | None = None,
    chemprop_label_col: str | None = None,
    task_type: str = "auto",
    metadata_cols: list[str] | None = None,
    deduplicate_smiles: bool = False,
) -> dict[str, Any]:
    dataset = Path(dataset_path)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(dataset, low_memory=False)
    if label_col not in frame.columns:
        raise ValueError(f"label column {label_col!r} not found")
    resolved_smiles = _resolve_column(frame, smiles_col, SMILES_CANDIDATES, role="smiles")
    work_cols = [resolved_smiles, label_col]
    selected_metadata = _selected_metadata_cols(
        frame,
        metadata_cols,
        smiles_col=resolved_smiles,
        label_col=label_col,
    )
    work = frame[[*work_cols, *selected_metadata]].copy()
    work["_atlas_row_id"] = frame.index.astype(int)
    work[resolved_smiles] = work[resolved_smiles].astype(str).str.strip()
    nonempty_smiles = work[resolved_smiles].ne("") & work[resolved_smiles].ne("nan")
    nonmissing_label = pd.to_numeric(work[label_col], errors="coerce").notna()
    staged = work.loc[nonempty_smiles & nonmissing_label].copy()
    if staged.empty:
        raise ValueError("no rows remain after dropping missing SMILES and labels")
    resolved_task_type = _infer_task_type(staged[label_col], task_type)
    output_label_col = chemprop_label_col or label_col
    staged_label = _prepare_label(staged[label_col], resolved_task_type)
    if deduplicate_smiles:
        duplicate_conflicts = (
            pd.DataFrame({"smiles": staged[resolved_smiles], "label": staged_label})
            .groupby("smiles")["label"]
            .nunique(dropna=True)
        )
        conflicting = duplicate_conflicts[duplicate_conflicts > 1]
        if not conflicting.empty:
            examples = conflicting.index.astype(str).tolist()[:5]
            raise ValueError(f"cannot deduplicate SMILES with conflicting labels; examples: {examples}")
        staged = staged.assign(_chemprop_label=staged_label).drop_duplicates(resolved_smiles, keep="first")
        staged_label = staged["_chemprop_label"]
    staged = staged.reset_index(drop=True)
    staged["_chemprop_row_id"] = staged.index.astype(int)

    chemprop_frame = pd.DataFrame(
        {
            "smiles": staged[resolved_smiles],
            output_label_col: staged_label.reset_index(drop=True),
        }
    )
    chemprop_path = out / "chemprop_input.csv"
    metadata_path = out / "chemprop_metadata.csv"
    manifest_path = out / "chemprop_training_manifest.json"
    chemprop_frame.to_csv(chemprop_path, index=False)
    metadata_frame = staged[
        ["_chemprop_row_id", "_atlas_row_id", resolved_smiles, label_col, *selected_metadata]
    ].copy()
    metadata_frame.to_csv(metadata_path, index=False)
    manifest: dict[str, Any] = {
        "adapter": "chemprop",
        "dataset": str(dataset),
        "label_col": label_col,
        "chemprop_label_col": output_label_col,
        "smiles_col": resolved_smiles,
        "task_type": resolved_task_type,
        "input_rows": int(len(frame)),
        "staged_rows": int(len(chemprop_frame)),
        "dropped_missing_smiles": int((~nonempty_smiles).sum()),
        "dropped_missing_label": int((nonempty_smiles & ~nonmissing_label).sum()),
        "deduplicate_smiles": bool(deduplicate_smiles),
        "label_counts": _label_counts(chemprop_frame[output_label_col]),
        "outputs": {
            "chemprop_input": str(chemprop_path),
            "metadata": str(metadata_path),
            "manifest": str(manifest_path),
        },
        "policy": (
            "Chemprop is staged as an external molecular baseline. Atlas leakage, split, "
            "source-transfer, and calibration audits remain required before interpreting claims."
        ),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def _chemprop_command(
    *,
    chemprop_bin: str | None,
    command_style: str,
    data_path: Path,
    task_type: str,
    model_dir: Path,
    extra_args: list[str] | None,
) -> list[str]:
    if command_style == "v1":
        executable = chemprop_bin or "chemprop_train"
        command = [
            executable,
            "--data_path",
            str(data_path),
            "--dataset_type",
            task_type,
            "--save_dir",
            str(model_dir),
        ]
    else:
        executable = chemprop_bin or "chemprop"
        command = [
            executable,
            "train",
            "--data-path",
            str(data_path),
            "--task-type",
            task_type,
            "--output-dir",
            str(model_dir),
        ]
    command.extend(extra_args or [])
    return command


def run_chemprop(
    manifest: dict[str, Any],
    *,
    chemprop_bin: str | None = None,
    command_style: str = "v2",
    extra_args: list[str] | None = None,
) -> int:
    outputs = manifest["outputs"]
    data_path = Path(outputs["chemprop_input"])
    manifest_path = Path(outputs["manifest"])
    model_dir = manifest_path.parent / "chemprop_model"
    model_dir.mkdir(parents=True, exist_ok=True)
    command = _chemprop_command(
        chemprop_bin=chemprop_bin,
        command_style=command_style,
        data_path=data_path,
        task_type=str(manifest["task_type"]),
        model_dir=model_dir,
        extra_args=extra_args,
    )
    executable = command[0]
    if Path(executable).name == executable and shutil.which(executable) is None:
        raise FileNotFoundError(
            f"Chemprop executable {executable!r} not found on PATH; install optional molecular tools "
            "or pass --chemprop-bin"
        )
    stdout_path = manifest_path.parent / "chemprop_stdout.log"
    stderr_path = manifest_path.parent / "chemprop_stderr.log"
    with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open("w", encoding="utf-8") as stderr:
        completed = subprocess.run(command, check=False, stdout=stdout, stderr=stderr)
    manifest["chemprop_run"] = {
        "command": command,
        "returncode": int(completed.returncode),
        "model_dir": str(model_dir),
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return int(completed.returncode)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Stage an Atlas model-ready CSV for Chemprop and optionally run Chemprop."
    )
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--smiles-col", default=None)
    parser.add_argument("--chemprop-label-col", default=None)
    parser.add_argument("--task-type", choices=["auto", "classification", "regression"], default="auto")
    parser.add_argument("--metadata-cols", nargs="*", default=None)
    parser.add_argument("--deduplicate-smiles", action="store_true")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--run", action="store_true", help="Run Chemprop after staging the input CSV.")
    parser.add_argument("--chemprop-bin", default=None)
    parser.add_argument(
        "--command-style",
        choices=["v2", "v1"],
        default="v2",
        help="Chemprop CLI shape: v2 uses `chemprop train`; v1 uses `chemprop_train`.",
    )
    parser.add_argument(
        "--chemprop-extra-arg",
        action="append",
        default=None,
        help="Extra argument appended to the Chemprop command. Repeat for multiple values.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    manifest = stage_chemprop_input(
        args.dataset,
        args.label,
        args.out_dir,
        smiles_col=args.smiles_col,
        chemprop_label_col=args.chemprop_label_col,
        task_type=args.task_type,
        metadata_cols=args.metadata_cols,
        deduplicate_smiles=args.deduplicate_smiles,
    )
    if not args.run:
        return 0
    return run_chemprop(
        manifest,
        chemprop_bin=args.chemprop_bin,
        command_style=args.command_style,
        extra_args=args.chemprop_extra_arg,
    )


if __name__ == "__main__":
    raise SystemExit(main())
