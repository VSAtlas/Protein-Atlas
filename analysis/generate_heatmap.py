from __future__ import annotations

import argparse
import html
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

from analysis._common import clean_text, format_float, parse_float, read_csv_rows, write_csv_rows


VALUE_CHOICES = (
    "atlas_score",
    "mmgbsa_score",
    "exposure_adjusted_score",
    "priority_score",
    "calibrated_activity_probability",
    "ml_prediction_score",
    "mechanism_graph_score",
    "matched_empirical_p",
    "fdr_q_value",
)
ORDER_CHOICES = ("input", "protein_class", "ligand_chemotype", "hierarchical", "mean_score")


def _pair_key(row: Mapping[str, Any]) -> tuple[str, str]:
    return clean_text(row.get("target_id")), clean_text(row.get("drug_id"))


def _significance_marker(row: Mapping[str, Any]) -> str:
    q_value = parse_float(row.get("fdr_q_value") or row.get("q_value"))
    p_value = parse_float(row.get("combined_empirical_p") or row.get("global_empirical_p"))
    if q_value is not None and q_value < 0.05:
        return "***"
    if q_value is not None and q_value < 0.10:
        return "**"
    if p_value is not None and p_value < 0.05:
        return "*"
    return ""


def _load_significance(path: Path | None) -> dict[tuple[str, str], str]:
    if path is None or not path.exists():
        return {}
    markers: dict[tuple[str, str], str] = {}
    for row in read_csv_rows(path):
        target_id = clean_text(row.get("target_id"))
        drug_id = clean_text(row.get("drug_id"))
        marker = _significance_marker(row)
        if target_id and drug_id and marker:
            markers[(target_id, drug_id)] = marker
    return markers


def _mean(values: Sequence[float]) -> float:
    return sum(values) / len(values) if values else float("-inf")


def _ordered_labels(
    labels: list[str],
    matrix: dict[tuple[str, str], float],
    *,
    axis: str,
    order: str,
    meta: dict[str, str],
    other_labels: list[str],
) -> list[str]:
    if order == "input":
        return labels
    if order in {"protein_class", "ligand_chemotype"}:
        return sorted(labels, key=lambda label: (meta.get(label, "Unassigned"), label))
    if order == "mean_score":
        def score(label: str) -> float:
            vals = [
                matrix[(label, other)] if axis == "row" else matrix[(other, label)]
                for other in other_labels
                if (label, other) in matrix if axis == "row"
            ]
            if axis == "col":
                vals = [matrix[(other, label)] for other in other_labels if (other, label) in matrix]
            return _mean(vals)

        return sorted(labels, key=lambda label: (-score(label), label))
    if order == "hierarchical":
        clustered = _hierarchical_order(labels, matrix, axis=axis, other_labels=other_labels)
        if clustered:
            return clustered
    return labels


def _hierarchical_order(
    labels: list[str],
    matrix: dict[tuple[str, str], float],
    *,
    axis: str,
    other_labels: list[str],
) -> list[str]:
    if len(labels) < 3:
        return labels
    try:
        import numpy as np
        from scipy.cluster.hierarchy import leaves_list, linkage
    except Exception:
        return []
    rows: list[list[float]] = []
    all_values = list(matrix.values())
    fill = sum(all_values) / len(all_values) if all_values else 0.0
    for label in labels:
        vector = []
        for other in other_labels:
            key = (label, other) if axis == "row" else (other, label)
            vector.append(matrix.get(key, fill))
        rows.append(vector)
    try:
        order = leaves_list(linkage(np.array(rows), method="average", metric="euclidean"))
    except Exception:
        return []
    return [labels[int(idx)] for idx in order]


def _color(value: float | None, lo: float, hi: float) -> str:
    if value is None or not math.isfinite(value):
        return "#f1f3f4"
    if hi == lo:
        frac = 0.5
    else:
        frac = max(0.0, min(1.0, (value - lo) / (hi - lo)))
    if frac < 0.5:
        t = frac / 0.5
        r = int(49 + (250 - 49) * t)
        g = int(101 + (247 - 101) * t)
        b = int(163 + (247 - 163) * t)
    else:
        t = (frac - 0.5) / 0.5
        r = int(250 + (178 - 250) * t)
        g = int(247 + (24 - 247) * t)
        b = int(247 + (43 - 247) * t)
    return f"#{r:02x}{g:02x}{b:02x}"


def build_heatmap_matrix(
    pair_rows: Sequence[Mapping[str, Any]],
    *,
    value: str,
    row_order: str,
    col_order: str,
) -> tuple[list[str], list[str], dict[tuple[str, str], float], dict[str, str], dict[str, str]]:
    rows: list[str] = []
    cols: list[str] = []
    matrix: dict[tuple[str, str], float] = {}
    row_meta: dict[str, str] = {}
    col_meta: dict[str, str] = {}
    for row in pair_rows:
        target_id, drug_id = _pair_key(row)
        if not target_id or not drug_id:
            continue
        if target_id not in rows:
            rows.append(target_id)
        if drug_id not in cols:
            cols.append(drug_id)
        row_meta.setdefault(target_id, clean_text(row.get("protein_class")) or "Unassigned")
        col_meta.setdefault(drug_id, clean_text(row.get("ligand_chemotype")) or "Unassigned")
        score = parse_float(row.get(value))
        if score is None:
            continue
        key = (target_id, drug_id)
        previous = matrix.get(key)
        if previous is None:
            matrix[key] = score
        elif value == "mmgbsa_score":
            matrix[key] = min(previous, score)
        else:
            matrix[key] = max(previous, score)
    rows = _ordered_labels(rows, matrix, axis="row", order=row_order, meta=row_meta, other_labels=cols)
    cols = _ordered_labels(cols, matrix, axis="col", order=col_order, meta=col_meta, other_labels=rows)
    return rows, cols, matrix, row_meta, col_meta


def write_matrix_csv(path: Path, rows: Sequence[str], cols: Sequence[str], matrix: Mapping[tuple[str, str], float]) -> None:
    out_rows = []
    for target_id in rows:
        row: dict[str, Any] = {"target_id": target_id}
        for drug_id in cols:
            row[drug_id] = format_float(matrix.get((target_id, drug_id)))
        out_rows.append(row)
    write_csv_rows(path, out_rows, ["target_id", *cols])


def write_svg(
    path: Path,
    rows: Sequence[str],
    cols: Sequence[str],
    matrix: Mapping[tuple[str, str], float],
    *,
    value: str,
    markers: Mapping[tuple[str, str], str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cell = 24
    left = 190
    top = 130
    width = left + cell * len(cols) + 30
    height = top + cell * len(rows) + 40
    values = list(matrix.values())
    lo = min(values) if values else 0.0
    hi = max(values) if values else 1.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>text{font-family:Arial,sans-serif;font-size:10px}.label{fill:#202124}.marker{font-size:9px;font-weight:700;fill:#111}</style>",
        f'<text x="10" y="22" font-size="15" font-weight="700">Atlas heatmap: {html.escape(value)}</text>',
        '<text x="10" y="40" fill="#555">Cell titles preserve z-score values for hover in SVG viewers.</text>',
    ]
    for col_idx, drug_id in enumerate(cols):
        x = left + col_idx * cell + cell * 0.5
        parts.append(
            f'<text class="label" transform="translate({x:.1f},{top - 8}) rotate(-55)" text-anchor="start">{html.escape(drug_id[:42])}</text>'
        )
    for row_idx, target_id in enumerate(rows):
        y = top + row_idx * cell + cell * 0.65
        parts.append(f'<text class="label" x="{left - 8}" y="{y:.1f}" text-anchor="end">{html.escape(target_id[:52])}</text>')
    for row_idx, target_id in enumerate(rows):
        for col_idx, drug_id in enumerate(cols):
            score = matrix.get((target_id, drug_id))
            x = left + col_idx * cell
            y = top + row_idx * cell
            fill = _color(score, lo, hi)
            title_score = "missing" if score is None else format_float(score)
            marker = markers.get((target_id, drug_id), "")
            parts.append(f'<g><title>{html.escape(target_id)} x {html.escape(drug_id)}; {html.escape(value)} z score={html.escape(title_score)}; significance={html.escape(marker or "none")}</title>')
            parts.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{fill}" stroke="#ffffff" stroke-width="1"/>')
            if marker:
                parts.append(f'<text class="marker" x="{x + cell / 2:.1f}" y="{y + cell * 0.65:.1f}" text-anchor="middle">{html.escape(marker)}</text>')
            parts.append("</g>")
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def write_png(
    path: Path,
    rows: Sequence[str],
    cols: Sequence[str],
    matrix: Mapping[tuple[str, str], float],
    *,
    value: str,
    markers: Mapping[tuple[str, str], str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    data = np.array(
        [[matrix.get((target_id, drug_id), np.nan) for drug_id in cols] for target_id in rows],
        dtype=float,
    )
    fig_w = max(6.0, min(24.0, 2.0 + 0.28 * len(cols)))
    fig_h = max(4.0, min(24.0, 1.8 + 0.24 * len(rows)))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=160)
    masked = np.ma.masked_invalid(data)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#f1f3f4")
    image = ax.imshow(masked, aspect="auto", cmap=cmap)
    ax.set_title(f"Atlas heatmap: {value}")
    ax.set_xlabel("Ligand")
    ax.set_ylabel("Target")
    if len(cols) <= 80:
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(cols, rotation=65, ha="right", fontsize=6)
    else:
        ax.set_xticks([])
    if len(rows) <= 80:
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(rows, fontsize=6)
    else:
        ax.set_yticks([])
    for row_idx, target_id in enumerate(rows):
        for col_idx, drug_id in enumerate(cols):
            marker = markers.get((target_id, drug_id), "")
            if marker:
                ax.text(col_idx, row_idx, marker, ha="center", va="center", fontsize=6, color="black")
    fig.colorbar(image, ax=ax, label=value)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate Atlas target-ligand heatmaps.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--value", choices=VALUE_CHOICES, default="atlas_score")
    parser.add_argument("--row-order", choices=ORDER_CHOICES, default="protein_class")
    parser.add_argument("--col-order", choices=ORDER_CHOICES, default="ligand_chemotype")
    parser.add_argument("--significance-table", type=Path, default=None)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--matrix-out", type=Path, default=Path("outputs/analysis/heatmap_matrix.csv"))
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--svg-out", type=Path, default=None)
    args = parser.parse_args(argv)

    pair_rows = read_csv_rows(args.pair_table)
    rows, cols, matrix, _row_meta, _col_meta = build_heatmap_matrix(
        pair_rows,
        value=args.value,
        row_order=args.row_order,
        col_order=args.col_order,
    )
    markers = _load_significance(args.significance_table)
    write_matrix_csv(args.matrix_out, rows, cols, matrix)
    write_png(args.out, rows, cols, matrix, value=args.value, markers=markers)
    write_svg(args.svg_out or args.out.with_suffix(".svg"), rows, cols, matrix, value=args.value, markers=markers)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
