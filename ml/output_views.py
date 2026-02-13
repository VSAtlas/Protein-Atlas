from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd


_PREFERRED_COLUMNS = (
    "trial_id",
    "rank_by_inner",
    "status",
    "mode",
    "feature_variant",
    "model_family",
    "calibration_enabled",
    "calibration_method",
    "inner_mean_EF@1%",
    "inner_mean_PR_AUC",
    "holdout_PR_AUC",
    "holdout_EF@1%",
    "active",
    "p_active",
    "p_active_uncalibrated",
    "prob_calibrated",
    "prob_uncalibrated",
    "lig_smiles",
    "ex_rec_pdb",
    "pocket",
    "scaffold",
    "murcko_scaffold",
    "error",
)


def _truncate_text(value: object, max_chars: int) -> str:
    text = str(value)
    if len(text) <= max_chars:
        return text
    return text[: max(0, max_chars - 3)] + "..."


def _choose_columns(all_columns: list[str], max_cols: int) -> list[str]:
    selected: list[str] = []
    for name in _PREFERRED_COLUMNS:
        if name in all_columns and name not in selected:
            selected.append(name)
        if len(selected) >= max_cols:
            return selected
    for name in all_columns:
        if name in selected:
            continue
        selected.append(name)
        if len(selected) >= max_cols:
            break
    return selected


def _preview_dataframe(
    csv_path: Path,
    *,
    max_rows: int,
    max_cols: int,
    cell_chars: int,
) -> tuple[pd.DataFrame, list[str]]:
    header = pd.read_csv(csv_path, nrows=0)
    all_columns = [str(col) for col in header.columns.tolist()]
    selected_columns = _choose_columns(all_columns, max_cols=max_cols)
    preview = pd.read_csv(
        csv_path,
        usecols=selected_columns,
        nrows=max_rows,
        low_memory=False,
    )
    for col in preview.columns:
        preview[col] = preview[col].map(lambda value: _truncate_text(value, cell_chars))
    return preview, all_columns


def write_truncated_view(
    csv_path: Path,
    *,
    max_rows: int = 60,
    max_cols: int = 14,
    cell_chars: int = 44,
) -> Path:
    csv_path = Path(csv_path)
    preview_df, all_columns = _preview_dataframe(
        csv_path,
        max_rows=max_rows,
        max_cols=max_cols,
        cell_chars=cell_chars,
    )
    out_path = csv_path.with_name(f"{csv_path.stem}_truncated.txt")
    summary_lines = [
        f"source_csv: {csv_path.name}",
        f"columns_total: {len(all_columns)}",
        f"columns_shown: {len(preview_df.columns)}",
        f"rows_shown: {len(preview_df)}",
        f"limits: max_rows={max_rows} max_cols={max_cols} cell_chars={cell_chars}",
        "",
    ]
    table_text = preview_df.to_string(index=False) if not preview_df.empty else "<no rows>"
    out_path.write_text("\n".join(summary_lines) + table_text + "\n", encoding="utf-8")
    return out_path


def write_truncated_views_for_run_dir(
    run_dir: Path,
    *,
    max_rows: int = 60,
    max_cols: int = 14,
    cell_chars: int = 44,
) -> list[Path]:
    run_dir = Path(run_dir)
    written: list[Path] = []
    for csv_path in sorted(run_dir.glob("*.csv")):
        try:
            out = write_truncated_view(
                csv_path,
                max_rows=max_rows,
                max_cols=max_cols,
                cell_chars=cell_chars,
            )
            written.append(out)
        except Exception:
            # Keep runner execution robust even if one CSV is malformed.
            continue
    return written


def _iter_run_dirs(outputs_root: Path) -> Iterable[Path]:
    for path in sorted(outputs_root.iterdir()):
        if path.is_dir():
            yield path


def convert_existing_outputs(
    *,
    outputs_root: Path,
    include_outputs_root_csvs: bool = True,
    max_rows: int = 60,
    max_cols: int = 14,
    cell_chars: int = 44,
) -> list[Path]:
    outputs_root = Path(outputs_root)
    written: list[Path] = []
    if include_outputs_root_csvs:
        for csv_path in sorted(outputs_root.glob("*.csv")):
            try:
                written.append(
                    write_truncated_view(
                        csv_path,
                        max_rows=max_rows,
                        max_cols=max_cols,
                        cell_chars=cell_chars,
                    )
                )
            except Exception:
                continue
    for run_dir in _iter_run_dirs(outputs_root):
        written.extend(
            write_truncated_views_for_run_dir(
                run_dir,
                max_rows=max_rows,
                max_cols=max_cols,
                cell_chars=cell_chars,
            )
        )
    return written


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate truncated text views for ML CSV outputs."
    )
    parser.add_argument(
        "--run-dir",
        default=None,
        help="Single run directory to convert. If omitted, convert ml/outputs.",
    )
    parser.add_argument(
        "--outputs-root",
        default="ml/outputs",
        help="Root outputs directory when --run-dir is omitted.",
    )
    parser.add_argument("--max-rows", type=int, default=60)
    parser.add_argument("--max-cols", type=int, default=14)
    parser.add_argument("--cell-chars", type=int, default=44)
    return parser


def main() -> None:
    args = _build_arg_parser().parse_args()
    if args.run_dir:
        written = write_truncated_views_for_run_dir(
            Path(args.run_dir),
            max_rows=int(args.max_rows),
            max_cols=int(args.max_cols),
            cell_chars=int(args.cell_chars),
        )
    else:
        written = convert_existing_outputs(
            outputs_root=Path(args.outputs_root),
            max_rows=int(args.max_rows),
            max_cols=int(args.max_cols),
            cell_chars=int(args.cell_chars),
        )
    print(f"[ml.output_views] written={len(written)}")


if __name__ == "__main__":
    main()
