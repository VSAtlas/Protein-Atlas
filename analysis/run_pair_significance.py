from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any, Mapping, Sequence

from analysis._common import clean_text, format_float, parse_float, read_csv_rows, write_csv_rows


OUTPUT_COLUMNS = [
    "drug_id",
    "target_id",
    "pdb_id",
    "score",
    "score_field",
    "score_direction",
    "ligand_empirical_p",
    "target_empirical_p",
    "global_empirical_p",
    "combined_empirical_p",
    "fdr_q_value",
    "significant_q_0_10",
    "significant_q_0_05",
]


def _empirical_p(observed: float, null_scores: list[float], *, higher_is_better: bool) -> float:
    if higher_is_better:
        extreme = sum(1 for score in null_scores if score >= observed)
    else:
        extreme = sum(1 for score in null_scores if score <= observed)
    return (extreme + 1.0) / (len(null_scores) + 1.0)


def _bh_q_values(p_values: list[float]) -> list[float]:
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    q_values = [1.0] * len(p_values)
    running = 1.0
    total = len(p_values)
    for rank_from_end, (idx, p_value) in enumerate(reversed(indexed), start=1):
        rank = total - rank_from_end + 1
        running = min(running, p_value * total / rank)
        q_values[idx] = min(1.0, running)
    return q_values


def _score_rows(rows: Sequence[Mapping[str, Any]], score_field: str) -> list[tuple[Mapping[str, Any], float]]:
    scored = []
    for row in rows:
        score = parse_float(row.get(score_field))
        if score is not None:
            scored.append((row, score))
    return scored


def compute_pair_significance(
    rows: Sequence[Mapping[str, Any]],
    *,
    score_field: str,
    direction: str,
) -> list[dict[str, Any]]:
    higher = direction == "higher"
    scored = _score_rows(rows, score_field)
    global_scores = [score for _row, score in scored]
    by_ligand: dict[str, list[float]] = {}
    by_target: dict[str, list[float]] = {}
    for row, score in scored:
        by_ligand.setdefault(clean_text(row.get("drug_id")), []).append(score)
        by_target.setdefault(clean_text(row.get("target_id")), []).append(score)

    out: list[dict[str, Any]] = []
    combined: list[float] = []
    for row, score in scored:
        ligand_scores = by_ligand.get(clean_text(row.get("drug_id")), [])
        target_scores = by_target.get(clean_text(row.get("target_id")), [])
        ligand_p = _empirical_p(score, ligand_scores, higher_is_better=higher)
        target_p = _empirical_p(score, target_scores, higher_is_better=higher)
        global_p = _empirical_p(score, global_scores, higher_is_better=higher)
        combined_p = max(ligand_p, target_p, global_p)
        combined.append(combined_p)
        out.append(
            {
                "drug_id": clean_text(row.get("drug_id")),
                "target_id": clean_text(row.get("target_id")),
                "pdb_id": clean_text(row.get("pdb_id")),
                "score": format_float(score),
                "score_field": score_field,
                "score_direction": direction,
                "ligand_empirical_p": format_float(ligand_p),
                "target_empirical_p": format_float(target_p),
                "global_empirical_p": format_float(global_p),
                "combined_empirical_p": format_float(combined_p),
            }
        )

    q_values = _bh_q_values(combined)
    for row, q_value in zip(out, q_values):
        row["fdr_q_value"] = format_float(q_value)
        row["significant_q_0_10"] = "1" if q_value < 0.10 else "0"
        row["significant_q_0_05"] = "1" if q_value < 0.05 else "0"
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Compute empirical pair-level Atlas significance.")
    parser.add_argument("--pair-table", required=True, type=Path)
    parser.add_argument("--score", default="atlas_score")
    parser.add_argument("--direction", choices=("higher", "lower"), default="higher")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)

    rows = read_csv_rows(args.pair_table)
    out = compute_pair_significance(rows, score_field=args.score, direction=args.direction)
    write_csv_rows(args.out, out, OUTPUT_COLUMNS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
