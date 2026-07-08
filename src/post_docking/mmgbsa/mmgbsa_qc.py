from __future__ import annotations

import csv
import json
import math
import os
import statistics
from pathlib import Path
from typing import Any, Mapping, Sequence


def _cfg_float(cfg: Mapping[str, Any], key: str, default: float) -> float:
    try:
        return float(cfg.get(key, default))
    except Exception:
        return default


def _cfg_int(cfg: Mapping[str, Any], key: str, default: int) -> int:
    try:
        return int(cfg.get(key, default))
    except Exception:
        return default


def _mean(values: Sequence[float]) -> float | None:
    return float(statistics.mean(values)) if values else None


def _sample_sd(values: Sequence[float]) -> float | None:
    if len(values) < 2:
        return None
    return float(statistics.stdev(values))


def summarize_values(values: Sequence[float]) -> dict[str, Any]:
    vals = [float(v) for v in values]
    n = len(vals)
    sd = _sample_sd(vals)
    sem = (sd / math.sqrt(n)) if sd is not None and n > 0 else None
    ci95 = (1.96 * sem) if sem is not None else None
    return {
        "n": n,
        "mean": _mean(vals),
        "median": float(statistics.median(vals)) if vals else None,
        "sd": sd,
        "sem": sem,
        "ci95": ci95,
        "min": min(vals) if vals else None,
        "max": max(vals) if vals else None,
    }


def block_average(values: Sequence[float], cfg: Mapping[str, Any]) -> dict[str, Any]:
    vals = [float(v) for v in values]
    n = len(vals)
    if not vals:
        return {"n_blocks": 0, "block_size": 0, "block_means": []}

    block_size = _cfg_int(cfg, "MMGBSA_QC_BLOCK_SIZE", 0)
    if block_size <= 0:
        target_blocks = max(1, _cfg_int(cfg, "MMGBSA_QC_MIN_BLOCKS", 5))
        block_size = max(1, n // target_blocks)

    blocks = [vals[idx : idx + block_size] for idx in range(0, n, block_size)]
    if len(blocks) > 1 and len(blocks[-1]) < max(1, block_size // 2):
        blocks[-2].extend(blocks[-1])
        blocks = blocks[:-1]
    block_means = [float(statistics.mean(block)) for block in blocks if block]
    stats = summarize_values(block_means)
    return {
        "n_blocks": len(block_means),
        "block_size": block_size,
        "block_means": block_means,
        "block_mean": stats["mean"],
        "block_sd": stats["sd"],
        "block_sem": stats["sem"],
        "block_ci95": stats["ci95"],
        "block_range": (
            max(block_means) - min(block_means) if len(block_means) >= 2 else None
        ),
    }


def detect_outliers(values: Sequence[float], z_threshold: float = 2.5) -> list[int]:
    vals = [float(v) for v in values]
    if len(vals) < 3:
        return []
    center = statistics.median(vals)
    deviations = [abs(v - center) for v in vals]
    mad = statistics.median(deviations)
    if mad <= 0:
        sd = statistics.stdev(vals)
        if sd <= 0:
            return []
        return [idx for idx, val in enumerate(vals) if abs(val - center) / sd > z_threshold]
    return [
        idx
        for idx, val in enumerate(vals)
        if 0.6745 * abs(val - center) / mad > z_threshold
    ]


def frame_qc(values: Sequence[float], cfg: Mapping[str, Any]) -> dict[str, Any]:
    vals = [float(v) for v in values]
    stats = summarize_values(vals)
    blocks = block_average(vals, cfg)
    min_frames = _cfg_int(cfg, "MMGBSA_QC_MIN_FRAMES", 1)
    max_sem = _cfg_float(cfg, "MMGBSA_QC_MAX_SEM_KCAL", 999999.0)
    max_block_range = _cfg_float(cfg, "MMGBSA_QC_MAX_BLOCK_RANGE_KCAL", 999999.0)
    min_blocks = _cfg_int(cfg, "MMGBSA_QC_MIN_BLOCKS", 1)

    sem = stats.get("sem")
    block_range = blocks.get("block_range")
    pass_min_frames = int(stats["n"]) >= min_frames
    pass_sem = sem is None or float(sem) <= max_sem
    pass_blocks = int(blocks.get("n_blocks", 0)) >= min_blocks
    pass_block_range = block_range is None or float(block_range) <= max_block_range
    first_half = vals[: len(vals) // 2]
    second_half = vals[len(vals) // 2 :]
    drift = None
    if first_half and second_half:
        drift = float(statistics.mean(second_half) - statistics.mean(first_half))
    ok = pass_min_frames and pass_sem and pass_blocks and pass_block_range
    return {
        **stats,
        **{f"block_{k}": v for k, v in blocks.items() if k != "block_means"},
        "block_means": blocks.get("block_means", []),
        "drift_second_minus_first": drift,
        "qc_pass": ok,
        "qc_checks": {
            "min_frames": pass_min_frames,
            "sem": pass_sem,
            "min_blocks": pass_blocks,
            "block_range": pass_block_range,
        },
    }


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    tmp_path = path.with_suffix(path.suffix + ".part")
    with tmp_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp_path, path)


def _line_svg(values: Sequence[float], title: str) -> str:
    width = 760
    height = 280
    pad = 36
    vals = [float(v) for v in values]
    if not vals:
        return "<svg xmlns='http://www.w3.org/2000/svg' width='760' height='280'></svg>\n"
    lo = min(vals)
    hi = max(vals)
    if math.isclose(lo, hi):
        lo -= 1.0
        hi += 1.0
    denom_x = max(1, len(vals) - 1)
    points = []
    for idx, val in enumerate(vals):
        x = pad + idx * (width - 2 * pad) / denom_x
        y = height - pad - (val - lo) * (height - 2 * pad) / (hi - lo)
        points.append(f"{x:.2f},{y:.2f}")
    return (
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' "
        f"viewBox='0 0 {width} {height}'>\n"
        "<rect width='100%' height='100%' fill='white'/>\n"
        f"<text x='{pad}' y='22' font-family='sans-serif' font-size='14'>{title}</text>\n"
        f"<line x1='{pad}' y1='{height-pad}' x2='{width-pad}' y2='{height-pad}' stroke='#333'/>\n"
        f"<line x1='{pad}' y1='{pad}' x2='{pad}' y2='{height-pad}' stroke='#333'/>\n"
        f"<polyline fill='none' stroke='#1f77b4' stroke-width='2' points='{' '.join(points)}'/>\n"
        f"<text x='{pad}' y='{height-8}' font-family='sans-serif' font-size='11'>frame</text>\n"
        f"<text x='{width-pad-80}' y='{height-8}' font-family='sans-serif' font-size='11'>n={len(vals)}</text>\n"
        "</svg>\n"
    )


def write_frame_qc_artifacts(
    out_dir: Path,
    frame_values: Sequence[float],
    cfg: Mapping[str, Any],
    *,
    label: str,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    vals = [float(v) for v in frame_values]
    qc = frame_qc(vals, cfg)

    frame_csv = out_dir / "frame_delta_total.csv"
    with frame_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["frame_index", "delta_total"])
        writer.writeheader()
        for idx, val in enumerate(vals, start=1):
            writer.writerow({"frame_index": idx, "delta_total": f"{val:.8g}"})

    block_csv = out_dir / "block_means.csv"
    with block_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["block_index", "delta_total_mean"])
        writer.writeheader()
        for idx, val in enumerate(qc.get("block_means", []), start=1):
            writer.writerow({"block_index": idx, "delta_total_mean": f"{float(val):.8g}"})

    svg_path = out_dir / "delta_total_convergence.svg"
    svg_path.write_text(_line_svg(vals, f"{label} DELTA TOTAL convergence"), encoding="utf-8")

    qc_payload = {
        "label": label,
        "qc": qc,
        "artifacts": {
            "frame_csv": str(frame_csv),
            "block_csv": str(block_csv),
            "convergence_svg": str(svg_path),
        },
    }
    _write_json(out_dir / "mmgbsa_frame_qc.json", qc_payload)
    return qc_payload
