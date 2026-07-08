from __future__ import annotations

import argparse
import math
import random
from pathlib import Path
from statistics import quantiles
from typing import Any, Mapping, Sequence

from analysis._common import (
    clean_text,
    format_float,
    parse_binary,
    parse_float,
    read_csv_rows,
    write_csv_rows,
)


SUMMARY_COLUMNS = [
    "score",
    "label",
    "metric",
    "observed",
    "ci_low",
    "ci_high",
    "drug_label_shuffle_p",
    "target_label_shuffle_p",
    "degree_preserving_label_shuffle_p",
    "n_pairs",
    "n_labeled_positive",
    "positive_rate_overall",
]
PERMUTATION_COLUMNS = ["score", "label", "null_model", "permutation", "metric", "value"]


def _is_positive(value: Any) -> int:
    parsed = parse_binary(value)
    if parsed is not None:
        return parsed
    text = clean_text(value).lower()
    return 1 if text in {"literature_supported", "supported", "known", "curated"} else 0


def _ranked(
    rows: Sequence[Mapping[str, Any]],
    score: str,
    labels: Sequence[int],
    *,
    higher_is_better: bool = True,
) -> list[tuple[float, int]]:
    ranked: list[tuple[float, int]] = []
    for row, label in zip(rows, labels):
        score_value = parse_float(row.get(score))
        if score_value is not None:
            ranked.append((score_value, label))
    ranked.sort(key=lambda item: item[0], reverse=higher_is_better)
    return ranked


def _ef_at_fraction(ranked: Sequence[tuple[float, int]], fraction: float) -> float:
    if not ranked:
        return 0.0
    positives = sum(label for _score, label in ranked)
    if positives <= 0:
        return 0.0
    top_n = max(1, int(len(ranked) * fraction + 0.999999))
    top_rate = sum(label for _score, label in ranked[:top_n]) / top_n
    overall_rate = positives / len(ranked)
    return top_rate / overall_rate if overall_rate > 0 else 0.0


def _precision_at_k(ranked: Sequence[tuple[float, int]], k: int) -> float:
    if not ranked:
        return 0.0
    top = ranked[: max(1, min(k, len(ranked)))]
    return sum(label for _score, label in top) / len(top)


def _recall_at_k(ranked: Sequence[tuple[float, int]], k: int) -> float:
    positives = sum(label for _score, label in ranked)
    if positives <= 0:
        return 0.0
    top = ranked[: max(1, min(k, len(ranked)))]
    return sum(label for _score, label in top) / positives


def _average_precision(ranked: Sequence[tuple[float, int]]) -> float:
    positives = sum(label for _score, label in ranked)
    if positives <= 0:
        return 0.0
    hits = 0
    precision_sum = 0.0
    for idx, (_score, label) in enumerate(ranked, start=1):
        if label:
            hits += 1
            precision_sum += hits / idx
    return precision_sum / positives


def _metrics(ranked: Sequence[tuple[float, int]], top_ks: Sequence[int]) -> dict[str, float]:
    out = {
        "EF@1%": _ef_at_fraction(ranked, 0.01),
        "EF@5%": _ef_at_fraction(ranked, 0.05),
        "EF@10%": _ef_at_fraction(ranked, 0.10),
        "AUPRC": _average_precision(ranked),
    }
    for k in top_ks:
        out[f"Precision@{k}"] = _precision_at_k(ranked, k)
        out[f"Recall@{k}"] = _recall_at_k(ranked, k)
    return out


def _shuffle_entity_labels(
    rows: Sequence[Mapping[str, Any]],
    labels: Sequence[int],
    *,
    entity_field: str,
    rng: random.Random,
) -> list[int]:
    entities: list[str] = []
    label_by_entity: dict[str, int] = {}
    for row, label in zip(rows, labels):
        entity = clean_text(row.get(entity_field))
        if entity not in label_by_entity:
            entities.append(entity)
            label_by_entity[entity] = 0
        label_by_entity[entity] = max(label_by_entity[entity], label)
    shuffled_values = [label_by_entity[entity] for entity in entities]
    rng.shuffle(shuffled_values)
    shuffled_by_entity = dict(zip(entities, shuffled_values))
    return [shuffled_by_entity.get(clean_text(row.get(entity_field)), 0) for row in rows]


def _degree_preserving_shuffle(
    rows: Sequence[Mapping[str, Any]],
    labels: Sequence[int],
    *,
    rng: random.Random,
) -> list[int]:
    drug_degree: dict[str, int] = {}
    target_degree: dict[str, int] = {}
    for row in rows:
        drug_degree[clean_text(row.get("drug_id"))] = drug_degree.get(clean_text(row.get("drug_id")), 0) + 1
        target_degree[clean_text(row.get("target_id"))] = target_degree.get(clean_text(row.get("target_id")), 0) + 1
    bins: dict[tuple[int, int], list[int]] = {}
    for idx, row in enumerate(rows):
        key = (
            min(5, int(math.log2(max(1, drug_degree.get(clean_text(row.get("drug_id")), 1))))),
            min(5, int(math.log2(max(1, target_degree.get(clean_text(row.get("target_id")), 1))))),
        )
        bins.setdefault(key, []).append(idx)
    shuffled = list(labels)
    for indices in bins.values():
        values = [labels[idx] for idx in indices]
        rng.shuffle(values)
        for idx, value in zip(indices, values):
            shuffled[idx] = value
    return shuffled


def _ci(values: Sequence[float]) -> tuple[str, str]:
    if len(values) < 2:
        return "", ""
    qs = quantiles(values, n=100, method="inclusive")
    return format_float(qs[1]), format_float(qs[96])


def run_enrichment(
    rows: Sequence[Mapping[str, Any]],
    *,
    score: str,
    label: str,
    n_permutations: int,
    n_bootstraps: int,
    top_ks: Sequence[int],
    seed: int,
    direction: str = "higher",
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rng = random.Random(seed)
    higher_is_better = direction != "lower"
    labels = [_is_positive(row.get(label)) for row in rows]
    ranked = _ranked(rows, score, labels, higher_is_better=higher_is_better)
    observed = _metrics(ranked, top_ks)
    positive_count = sum(label_value for _score, label_value in ranked)
    positive_rate = positive_count / len(ranked) if ranked else 0.0

    null_metrics: dict[str, dict[str, list[float]]] = {
        "drug_label_shuffle": {name: [] for name in observed},
        "target_label_shuffle": {name: [] for name in observed},
        "degree_preserving_label_shuffle": {name: [] for name in observed},
    }
    permutation_rows: list[dict[str, Any]] = []
    for idx in range(n_permutations):
        for null_model, entity_field in (("drug_label_shuffle", "drug_id"), ("target_label_shuffle", "target_id")):
            shuffled = _shuffle_entity_labels(rows, labels, entity_field=entity_field, rng=rng)
            values = _metrics(_ranked(rows, score, shuffled, higher_is_better=higher_is_better), top_ks)
            for metric, value in values.items():
                null_metrics[null_model][metric].append(value)
                permutation_rows.append(
                    {
                        "score": score,
                        "label": label,
                        "null_model": null_model,
                        "permutation": idx + 1,
                        "metric": metric,
                        "value": format_float(value),
                    }
                )
        shuffled = _degree_preserving_shuffle(rows, labels, rng=rng)
        values = _metrics(_ranked(rows, score, shuffled, higher_is_better=higher_is_better), top_ks)
        for metric, value in values.items():
            null_metrics["degree_preserving_label_shuffle"][metric].append(value)
            permutation_rows.append(
                {
                    "score": score,
                    "label": label,
                    "null_model": "degree_preserving_label_shuffle",
                    "permutation": idx + 1,
                    "metric": metric,
                    "value": format_float(value),
                }
            )

    bootstrap_values: dict[str, list[float]] = {name: [] for name in observed}
    if ranked and n_bootstraps > 0:
        for _idx in range(n_bootstraps):
            sample = [ranked[rng.randrange(len(ranked))] for _ in ranked]
            for metric, value in _metrics(sample, top_ks).items():
                bootstrap_values[metric].append(value)

    summary_rows: list[dict[str, Any]] = []
    for metric, value in observed.items():
        ci_low, ci_high = _ci(bootstrap_values.get(metric, []))
        drug_null = null_metrics["drug_label_shuffle"][metric]
        target_null = null_metrics["target_label_shuffle"][metric]
        degree_null = null_metrics["degree_preserving_label_shuffle"][metric]
        drug_p = (sum(1 for v in drug_null if v >= value) + 1) / (len(drug_null) + 1) if drug_null else 1.0
        target_p = (sum(1 for v in target_null if v >= value) + 1) / (len(target_null) + 1) if target_null else 1.0
        degree_p = (sum(1 for v in degree_null if v >= value) + 1) / (len(degree_null) + 1) if degree_null else 1.0
        summary_rows.append(
            {
                "score": score,
                "label": label,
                "metric": metric,
                "observed": format_float(value),
                "ci_low": ci_low,
                "ci_high": ci_high,
                "drug_label_shuffle_p": format_float(drug_p),
                "target_label_shuffle_p": format_float(target_p),
                "degree_preserving_label_shuffle_p": format_float(degree_p),
                "n_pairs": str(len(ranked)),
                "n_labeled_positive": str(positive_count),
                "positive_rate_overall": format_float(positive_rate),
            }
        )
    return summary_rows, permutation_rows


def write_enrichment_curve(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
    *,
    score: str,
    label: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [_is_positive(row.get(label)) for row in rows]
    ranked = _ranked(rows, score, labels)
    positives = sum(label_value for _score, label_value in ranked)
    xs: list[float] = []
    precision: list[float] = []
    recall: list[float] = []
    hits = 0
    for idx, (_score, label_value) in enumerate(ranked, start=1):
        hits += label_value
        xs.append(idx / len(ranked) if ranked else 0.0)
        precision.append(hits / idx)
        recall.append(hits / positives if positives else 0.0)
    fig, ax = plt.subplots(figsize=(6.5, 4.2), dpi=160)
    ax.plot(xs, precision, label="Precision")
    ax.plot(xs, recall, label="Recall")
    ax.set_xlabel("Top-ranked fraction")
    ax.set_ylabel("Literature-supported label recovery")
    ax.set_title(f"Atlas enrichment curve: {score}")
    ax.set_ylim(0, 1.02)
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _parse_top_ks(value: str) -> list[int]:
    return [int(part) for part in value.split(",") if part.strip()]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run Atlas enrichment analysis with label-shuffle controls.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--score", default="atlas_score")
    parser.add_argument("--label", default="literature_supported_label")
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--bootstraps", type=int, default=1000)
    parser.add_argument("--top-k", default="10,25,50")
    parser.add_argument("--seed", type=int, default=1337)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--permutation-out", type=Path, default=None)
    parser.add_argument("--curve-out", type=Path, default=None)
    args = parser.parse_args(argv)

    rows = read_csv_rows(args.pair_table)
    summary, permutations = run_enrichment(
        rows,
        score=args.score,
        label=args.label,
        n_permutations=args.permutations,
        n_bootstraps=args.bootstraps,
        top_ks=_parse_top_ks(args.top_k),
        seed=args.seed,
    )
    write_csv_rows(args.out, summary, SUMMARY_COLUMNS)
    permutation_out = args.permutation_out or args.out.with_name("permutation_results.csv")
    write_csv_rows(permutation_out, permutations, PERMUTATION_COLUMNS)
    if args.curve_out:
        write_enrichment_curve(args.curve_out, rows, score=args.score, label=args.label)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
