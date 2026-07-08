from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd  # type: ignore[import-untyped]

from analysis.dud_eval_core.discovery import make_target_key
from analysis.dud_eval_core.labels import _make_placeholder_row
from analysis.dud_eval_core.log import dbg
from analysis.dud_eval_core.types import TargetEvaluation, TargetSpec

DUD_E_PROTEIN_CLASSES: Dict[str, Tuple[str, ...]] = {
    "Protease": (
        "1SQT",
        "1W7X",
        "1XL2",
        "1YPE",
        "2AYW",
        "2CNK",
        "2I78",
        "2OI0",
        "2ZEC",
        "3BKL",
        "3CHP",
        "3G6Z",
        "3KL6",
        "3L5D",
        "830C",
    ),
    "Nuclear receptor": (
        "1MV9",
        "1Q4X",
        "1SJ0",
        "2AA2",
        "2AM9",
        "2FSZ",
        "2GTK",
        "2P54",
        "2ZNP",
        "3BQD",
        "3KBA",
    ),
    "GPCR": ("2VT4", "3EML", "3NY8", "3ODU", "3PBL"),
    "Ion channel": ("1VSO", "3KGC"),
    "Cytochrome P450": ("1R9O", "3NXU"),
    "Kinase": (
        "1H00",
        "2ETR",
        "2HZI",
        "2OWB",
        "2P2I",
        "2QD9",
        "2RGP",
        "3BIZ",
        "3BZ3",
        "3C4F",
        "3CQW",
        "3D0E",
        "3D4Q",
        "3EL8",
        "3EQH",
        "3HMM",
        "3KRJ",
    ),
    "Miscellaneous": ("1UYG", "3CJO"),
    "Other enzyme": (
        "1B9V",
        "1BCD",
        "1C8K",
        "1D3G",
        "1E66",
        "1J4H",
        "1KVO",
        "1L2S",
        "1LI4",
        "1LRU",
        "1NJS",
        "1QW6",
        "1S3B",
        "1SYN",
        "1UDT",
        "1ZW5",
        "2AZR",
        "2B8T",
        "2E1W",
        "2H7L",
        "2HV5",
        "2OYU",
        "2V3F",
        "3BGS",
        "3BWM",
        "3CCW",
        "3E37",
        "3F07",
        "3F9M",
        "3FRJ",
        "3L3M",
        "3LAN",
        "3MAX",
        "3NF7",
        "3NXO",
    ),
}

def _normalize_pdb_token(value: Any) -> str:
    token = str(value or "").strip().upper()
    if token.endswith(".PDB"):
        token = token[:-4]
    return token


def _write_placeholder_metrics(out_dir: Path, placeholder: Any) -> None:
    """
    Ensure per-target output layout exists even when metrics cannot be computed.
    """
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        payload = (
            placeholder.to_dict()
            if hasattr(placeholder, "to_dict")
            else dict(placeholder)
        )
        if not isinstance(payload, dict) or not payload:
            payload = {"status_reason": "missing_metrics"}
        pd.DataFrame([payload]).to_csv(out_dir / "metrics.tsv", sep="\t", index=False)
    except Exception as exc:
        dbg("WARN", "metrics", f"placeholder_metrics_write_failed dir={out_dir} err={exc}")


def _class_slug(class_name: str) -> str:
    return class_name.lower().replace(" ", "_")


_CLASS_AGGREGATE_PLOT_FILENAME = "metric_boxplots.png"
_CLASS_OVERALL_PLOT_DIRNAME = "overall_by_metric"


def _class_metric_base_columns(bedroc_alpha: float) -> List[str]:
    return [
        "N",
        "n_actives",
        "actives_fraction",
        "ROC_AUC",
        "PR_AUC",
        "logAUC",
        "logAUC_adj",
        f"BEDROC_alpha_{bedroc_alpha:g}",
        "EF@1%",
        "EF@2%",
        "EF@5%",
        "EF@10%",
    ]


def _class_metric_stat_columns(bedroc_alpha: float) -> List[str]:
    out: List[str] = []
    for metric in _class_metric_base_columns(bedroc_alpha):
        out.extend([f"{metric}_std", f"{metric}_median", f"{metric}_iqr"])
    return out


def _class_plot_metric_columns(bedroc_alpha: float) -> List[str]:
    return [
        "ROC_AUC",
        "PR_AUC",
        "logAUC_adj",
        f"BEDROC_alpha_{bedroc_alpha:g}",
        "EF@1%",
        "EF@2%",
        "EF@5%",
        "EF@10%",
    ]


def _class_metric_plot_slug(metric_name: str) -> str:
    token = str(metric_name).strip().lower()
    token = token.replace("@", "_at_")
    token = token.replace("%", "pct")
    token = re.sub(r"[^a-z0-9]+", "_", token)
    return token.strip("_") or "metric"


def _finite_metric_values(rows: List[Dict[str, Any]], metric_name: str) -> np.ndarray:
    values: List[float] = []
    for row in rows:
        raw = row.get(metric_name, float("nan"))
        try:
            value = float(raw)
        except (TypeError, ValueError):
            continue
        if np.isfinite(value):
            values.append(value)
    return np.asarray(values, dtype=float)


def _class_metric_summary(rows: List[Dict[str, Any]], metric_name: str) -> Dict[str, float]:
    values = _finite_metric_values(rows, metric_name)
    if values.size == 0:
        return {
            metric_name: float("nan"),
            f"{metric_name}_std": float("nan"),
            f"{metric_name}_median": float("nan"),
            f"{metric_name}_iqr": float("nan"),
        }
    q1 = float(np.percentile(values, 25))
    q3 = float(np.percentile(values, 75))
    std = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
    return {
        metric_name: float(np.mean(values)),
        f"{metric_name}_std": std,
        f"{metric_name}_median": float(np.median(values)),
        f"{metric_name}_iqr": float(q3 - q1),
    }


def _build_class_metric_row(
    *,
    class_name: str,
    mode_slug: str,
    run_id: Optional[str],
    bedroc_alpha: float,
    status_reason: str,
    targets_total_in_class: int,
    targets_present_in_run: int,
    targets_evaluable: int,
    pooled_rows: int,
    pooled_actives: int,
    metric_rows: List[Dict[str, Any]],
) -> Dict[str, Any]:
    row = _make_placeholder_row(
        pdb_id=f"class:{class_name}",
        variant=None,
        ph_tag=None,
        run_id=run_id,
        bedroc_alpha=bedroc_alpha,
        status_reason=status_reason,
    ).to_dict()
    row["mode"] = mode_slug
    row["aggregation_method"] = "macro_per_target"
    row["targets_total_in_class"] = int(max(targets_total_in_class, 0))
    row["targets_present_in_run"] = int(max(targets_present_in_run, 0))
    row["targets_evaluable"] = int(max(targets_evaluable, 0))
    row["pooled_rows"] = int(max(pooled_rows, 0))
    row["pooled_actives"] = int(max(min(pooled_actives, pooled_rows), 0))
    for metric_name in _class_metric_base_columns(bedroc_alpha):
        row.update(_class_metric_summary(metric_rows, metric_name))
    row["status_reason"] = status_reason
    return row


def _make_class_placeholder_metrics_row(
    *,
    class_name: str,
    mode_slug: str,
    run_id: Optional[str],
    bedroc_alpha: float,
    status_reason: str,
    targets_total_in_class: int,
    targets_present_in_run: int,
    targets_evaluable: int,
    pooled_rows: int,
    pooled_actives: int,
) -> Dict[str, Any]:
    return _build_class_metric_row(
        class_name=class_name,
        mode_slug=mode_slug,
        run_id=run_id,
        bedroc_alpha=bedroc_alpha,
        status_reason=status_reason,
        targets_total_in_class=targets_total_in_class,
        targets_present_in_run=targets_present_in_run,
        targets_evaluable=targets_evaluable,
        pooled_rows=pooled_rows,
        pooled_actives=pooled_actives,
        metric_rows=[],
    )


def _write_class_metric_boxplots(
    *,
    out_dir: Path,
    class_name: str,
    mode_slug: str,
    bedroc_alpha: float,
    metric_rows: List[Dict[str, Any]],
) -> None:
    try:
        import matplotlib.pyplot as plt  # type: ignore[import-untyped]
    except Exception as exc:
        dbg(
            "WARN",
            "class.aggregate",
            f"boxplot_import_failed dir={out_dir} err={exc}",
        )
        return

    metric_names = [
        metric
        for metric in _class_plot_metric_columns(bedroc_alpha)
        if _finite_metric_values(metric_rows, metric).size > 0
    ]
    if not metric_names:
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    n_cols = min(3, len(metric_names))
    n_rows = int(math.ceil(len(metric_names) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.5 * n_cols, 3.8 * n_rows))
    axes_arr = np.atleast_1d(axes).ravel()
    for idx, metric_name in enumerate(metric_names):
        ax = axes_arr[idx]
        values = _finite_metric_values(metric_rows, metric_name)
        ax.boxplot(
            [values],
            widths=0.5,
            showmeans=True,
            patch_artist=True,
            boxprops={"facecolor": "#dbeafe", "edgecolor": "#1d4ed8"},
            whiskerprops={"color": "#1d4ed8"},
            capprops={"color": "#1d4ed8"},
            medianprops={"color": "#b45309", "linewidth": 2},
            meanprops={
                "marker": "D",
                "markerfacecolor": "#111827",
                "markeredgecolor": "#111827",
                "markersize": 5,
            },
        )
        mean_val = float(np.mean(values))
        std_val = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
        ax.errorbar(
            [1.0],
            [mean_val],
            yerr=[std_val],
            fmt="none",
            ecolor="#dc2626",
            elinewidth=1.5,
            capsize=5,
        )
        ax.set_title(metric_name)
        ax.set_xticks([1])
        ax.set_xticklabels([f"n={values.size}"])
        ax.grid(axis="y", linestyle=":", alpha=0.35)
    for ax in axes_arr[len(metric_names) :]:
        ax.axis("off")
    fig.suptitle(f"{class_name} ({mode_slug}) metric distributions")
    fig.tight_layout()
    fig.savefig(out_dir / _CLASS_AGGREGATE_PLOT_FILENAME, dpi=200)
    plt.close(fig)


def _write_overall_class_metric_boxplots(
    *,
    out_dir: Path,
    mode_slug: str,
    bedroc_alpha: float,
    class_metric_rows: Dict[str, List[Dict[str, Any]]],
    write_placeholder_images: bool,
) -> None:
    try:
        import matplotlib.pyplot as plt  # type: ignore[import-untyped]
    except Exception as exc:
        dbg(
            "WARN",
            "class.aggregate",
            f"overall_plot_import_failed dir={out_dir} err={exc}",
        )
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    for metric_name in _class_plot_metric_columns(bedroc_alpha):
        plot_path = out_dir / f"{_class_metric_plot_slug(metric_name)}_by_class.png"
        class_names: List[str] = []
        class_values: List[np.ndarray] = []
        for class_name in DUD_E_PROTEIN_CLASSES:
            values = _finite_metric_values(class_metric_rows.get(class_name, []), metric_name)
            if values.size == 0:
                continue
            class_names.append(class_name)
            class_values.append(values)

        if not class_values:
            if not write_placeholder_images:
                continue
            try:
                plt.figure(figsize=(8.0, 4.0))
                plt.axis("off")
                plt.text(
                    0.5,
                    0.56,
                    f"No evaluable class data for {metric_name}",
                    ha="center",
                    va="center",
                )
                plt.text(0.5, 0.40, f"mode: {mode_slug}", ha="center", va="center")
                plt.tight_layout()
                plt.savefig(plot_path, dpi=200)
                plt.close()
            except Exception as exc:
                dbg(
                    "WARN",
                    "class.aggregate",
                    (
                        f"overall_placeholder_plot_write_failed mode={mode_slug} "
                        f"metric={metric_name} err={exc}"
                    ),
                )
            continue

        fig, ax = plt.subplots(figsize=(max(8.0, 1.2 * len(class_names) + 4.0), 5.0))
        ax.boxplot(
            class_values,
            tick_labels=class_names,
            showmeans=True,
            patch_artist=True,
            boxprops={"facecolor": "#dbeafe", "edgecolor": "#1d4ed8"},
            whiskerprops={"color": "#1d4ed8"},
            capprops={"color": "#1d4ed8"},
            medianprops={"color": "#b45309", "linewidth": 2},
            meanprops={
                "marker": "D",
                "markerfacecolor": "#111827",
                "markeredgecolor": "#111827",
                "markersize": 5,
            },
        )
        for idx, values in enumerate(class_values, start=1):
            mean_val = float(np.mean(values))
            std_val = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
            ax.errorbar(
                [float(idx)],
                [mean_val],
                yerr=[std_val],
                fmt="none",
                ecolor="#dc2626",
                elinewidth=1.5,
                capsize=5,
            )
        ax.set_title(f"{metric_name} by class ({mode_slug})")
        ax.set_ylabel(metric_name)
        ax.tick_params(axis="x", rotation=35)
        ax.grid(axis="y", linestyle=":", alpha=0.35)
        fig.tight_layout()
        fig.savefig(plot_path, dpi=200)
        plt.close(fig)


def _write_placeholder_class_plots(
    *,
    out_dir: Path,
    class_name: str,
    mode_slug: str,
    status_reason: str,
) -> None:
    try:
        import matplotlib.pyplot as plt  # type: ignore[import-untyped]
    except Exception as exc:
        dbg(
            "WARN",
            "class.aggregate",
            f"placeholder_plot_import_failed dir={out_dir} err={exc}",
        )
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    for filename in (_CLASS_AGGREGATE_PLOT_FILENAME,):
        try:
            plt.figure(figsize=(5.5, 3.5))
            plt.axis("off")
            plt.text(0.5, 0.58, "No evaluable class data", ha="center", va="center")
            plt.text(
                0.5,
                0.42,
                f"status: {status_reason}",
                ha="center",
                va="center",
                fontsize=9,
            )
            plt.title(f"{class_name} ({mode_slug})")
            plt.tight_layout()
            plt.savefig(out_dir / filename, dpi=200)
            plt.close()
        except Exception as exc:
            dbg(
                "WARN",
                "class.aggregate",
                (
                    f"placeholder_plot_write_failed mode={mode_slug} class={class_name} "
                    f"plot={filename} err={exc}"
                ),
            )


def _target_keys_from_summary_df(summary_df: pd.DataFrame) -> Set[str]:
    if summary_df.empty or "pdb_id" not in summary_df.columns:
        return set()

    keys: Set[str] = set()
    has_variant = "variant" in summary_df.columns
    has_ph = "pH" in summary_df.columns
    for _, row in summary_df.iterrows():
        pdb_id = str(row.get("pdb_id", "") or "").strip()
        if not pdb_id:
            continue
        variant_raw = row.get("variant", "") if has_variant else ""
        ph_raw = row.get("pH", "") if has_ph else ""
        variant = ""
        ph_tag = ""
        if not pd.isna(variant_raw):
            variant = str(variant_raw).strip()
        if not pd.isna(ph_raw):
            ph_tag = str(ph_raw).strip()
        keys.add(make_target_key(pdb_id, variant or None, ph_tag or None))
    return keys


def _emit_class_aggregate_plots(
    *,
    analysis_root: Path,
    run_id: Optional[str],
    target_specs_by_key: Dict[str, TargetSpec],
    eval_results: Dict[str, TargetEvaluation],
    included_target_keys: Set[str],
    mode_slug: str,
    score_high_is_better: bool,
    hist_xlabel: str,
    title_suffix: str,
    mode_label: str,
    bedroc_alpha: float,
    logauc_lambda: float,
    class_min_targets: int,
    write_placeholder_images: bool,
) -> Tuple[Dict[str, TargetEvaluation], List[Dict[str, Any]]]:
    min_targets = max(int(class_min_targets), 1)

    class_results: Dict[str, TargetEvaluation] = {}
    generation_rows: List[Dict[str, Any]] = []
    class_metric_rows_by_class: Dict[str, List[Dict[str, Any]]] = {}
    class_root = analysis_root / "class_aggregates" / mode_slug
    class_root.mkdir(parents=True, exist_ok=True)

    for class_name, class_members in DUD_E_PROTEIN_CLASSES.items():
        class_slug = _class_slug(class_name)
        class_dir = class_root / class_slug
        class_dir.mkdir(parents=True, exist_ok=True)
        member_set = {_normalize_pdb_token(token) for token in class_members}
        class_metric_rows_by_class[class_name] = []
        present_target_keys: List[str] = []
        pooled_frames: List[pd.DataFrame] = []
        metric_rows: List[Dict[str, Any]] = []
        evaluable_target_keys: Set[str] = set()
        ligand_basenames: Set[str] = set()
        run_ids_present: Set[str] = set()
        has_run_id_col = False
        status_reason = "ok"
        for target_key in sorted(included_target_keys):
            spec = target_specs_by_key.get(target_key)
            if spec is None:
                continue
            if _normalize_pdb_token(spec.pdb_id) not in member_set:
                continue
            present_target_keys.append(target_key)
            eval_meta = eval_results.get(target_key)
            if eval_meta is None:
                continue
            ligand_basenames.update(eval_meta.ligand_basenames)
            run_ids_present.update(eval_meta.run_ids)
            has_run_id_col = has_run_id_col or eval_meta.has_run_id_column
            if eval_meta.z_selected is not None:
                z_selected = eval_meta.z_selected
                if (
                    isinstance(z_selected, pd.DataFrame)
                    and not z_selected.empty
                    and "is_active" in z_selected.columns
                    and "best_score" in z_selected.columns
                ):
                    frame = z_selected.loc[:, ["is_active", "best_score"]].copy()
                    frame["is_active"] = pd.to_numeric(frame["is_active"], errors="coerce")
                    frame["best_score"] = pd.to_numeric(
                        frame["best_score"], errors="coerce"
                    )
                    frame = frame.dropna(subset=["is_active", "best_score"])
                    if not frame.empty:
                        frame["is_active"] = frame["is_active"].astype(int)
                        pooled_frames.append(frame)
            if eval_meta.metrics is None:
                continue
            metric_row = (
                eval_meta.metrics.to_dict()
                if hasattr(eval_meta.metrics, "to_dict")
                else dict(eval_meta.metrics)
            )
            if not isinstance(metric_row, dict):
                continue
            metric_rows.append(metric_row)
            class_metric_rows_by_class[class_name].append(metric_row)
            evaluable_target_keys.add(target_key)

        pooled_rows = int(sum(len(frame) for frame in pooled_frames))
        pooled_actives = int(
            sum(int(frame["is_active"].sum()) for frame in pooled_frames)
        )
        class_eval = TargetEvaluation(
            metrics=None,
            ligand_basenames=set(),
            has_run_id_column=False,
            run_ids=set(),
            status_reason="no_evaluable_targets",
            z_selected=None,
        )

        if len(present_target_keys) < min_targets:
            if len(present_target_keys) == 0:
                status_reason = "no_targets_in_run"
            else:
                status_reason = f"insufficient_targets(min={min_targets})"
        elif not metric_rows:
            status_reason = "no_evaluable_targets"
        else:
            try:
                class_row = _build_class_metric_row(
                    class_name=class_name,
                    mode_slug=mode_slug,
                    run_id=run_id,
                    bedroc_alpha=bedroc_alpha,
                    status_reason="ok",
                    targets_total_in_class=len(member_set),
                    targets_present_in_run=len(present_target_keys),
                    targets_evaluable=len(evaluable_target_keys),
                    pooled_rows=pooled_rows,
                    pooled_actives=pooled_actives,
                    metric_rows=metric_rows,
                )
                _write_placeholder_metrics(class_dir, class_row)
                _write_class_metric_boxplots(
                    out_dir=class_dir,
                    class_name=class_name,
                    mode_slug=mode_slug,
                    bedroc_alpha=bedroc_alpha,
                    metric_rows=metric_rows,
                )
                class_eval = TargetEvaluation(
                    metrics=pd.Series(class_row),
                    ligand_basenames=ligand_basenames,
                    has_run_id_column=has_run_id_col,
                    run_ids=run_ids_present,
                    status_reason="ok",
                    z_selected=None,
                )
                status_reason = class_eval.status_reason or "ok"
            except Exception as exc:
                status_reason = "compute_error"
                dbg(
                    "WARN",
                    "class.aggregate",
                    f"mode={mode_slug} class={class_name} err={exc}",
                )

        if class_eval.metrics is None:
            placeholder = _make_class_placeholder_metrics_row(
                class_name=class_name,
                mode_slug=mode_slug,
                run_id=run_id,
                bedroc_alpha=bedroc_alpha,
                status_reason=status_reason,
                targets_total_in_class=len(member_set),
                targets_present_in_run=len(present_target_keys),
                targets_evaluable=len(evaluable_target_keys),
                pooled_rows=pooled_rows,
                pooled_actives=pooled_actives,
            )
            _write_placeholder_metrics(class_dir, placeholder)
            if write_placeholder_images:
                _write_placeholder_class_plots(
                    out_dir=class_dir,
                    class_name=class_name,
                    mode_slug=mode_slug,
                    status_reason=status_reason,
                )
            class_eval = TargetEvaluation(
                metrics=pd.Series(placeholder),
                ligand_basenames=class_eval.ligand_basenames,
                has_run_id_column=class_eval.has_run_id_column,
                run_ids=class_eval.run_ids,
                status_reason=status_reason,
                z_selected=class_eval.z_selected,
            )

        class_results[class_name] = class_eval
        source_metrics = (
            class_eval.metrics.to_dict()
            if class_eval.metrics is not None and hasattr(class_eval.metrics, "to_dict")
            else (
                dict(class_eval.metrics)
                if class_eval.metrics is not None
                else {}
            )
        )
        generation_row = {
            "run_id": str(run_id) if run_id else "(none)",
            "mode": mode_slug,
            "class_name": class_name,
            "class_slug": class_slug,
            "status_reason": status_reason,
            "aggregation_method": "macro_per_target",
            "targets_total_in_class": len(member_set),
            "targets_present_in_run": len(present_target_keys),
            "targets_evaluable": len(evaluable_target_keys),
            "pooled_rows": pooled_rows,
            "pooled_actives": pooled_actives,
        }
        for metric_name in _class_metric_base_columns(bedroc_alpha):
            generation_row[metric_name] = source_metrics.get(metric_name, float("nan"))
        for metric_name in _class_metric_stat_columns(bedroc_alpha):
            generation_row[metric_name] = source_metrics.get(metric_name, float("nan"))
        generation_rows.append(generation_row)
        dbg(
            "INFO",
            "class.aggregate",
            (
                f"mode={mode_slug} class={class_name} present={len(present_target_keys)} "
                f"evaluable={len(evaluable_target_keys)} ligands={pooled_rows} "
                f"status={status_reason}"
            ),
        )

    _write_overall_class_metric_boxplots(
        out_dir=class_root / _CLASS_OVERALL_PLOT_DIRNAME,
        mode_slug=mode_slug,
        bedroc_alpha=bedroc_alpha,
        class_metric_rows=class_metric_rows_by_class,
        write_placeholder_images=write_placeholder_images,
    )
    return class_results, generation_rows

