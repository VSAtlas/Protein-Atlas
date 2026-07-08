"""CLI for the protein-prep benchmark."""

from __future__ import annotations

import argparse
import json
import logging
import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Mapping, Sequence

from protein_prep.benchmark.reporting import write_outputs
from protein_prep.prep_benchmark import DEFAULT_PDB_IDS, BenchmarkRow, run_one

MAX_TARGET_WORKERS = 32


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Open receptor-preparation benchmark harness.")
    parser.add_argument("--pdb-ids", help="Comma-separated PDB IDs. Defaults to a 20-target open co-crystal seed set.")
    parser.add_argument("--pdb-id-file", type=Path, help="Text file with one PDB ID per line.")
    parser.add_argument("--out-dir", type=Path, default=Path("outputs/protein_prep_benchmark"))
    parser.add_argument("--limit", type=int, default=0, help="Limit number of PDB IDs for smoke runs.")
    parser.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("PROTEIN_PREP_BENCHMARK_WORKERS", "1")),
        help="PDB-level parallel workers for prep-only benchmark runs. Capped at 32.",
    )
    parser.add_argument(
        "--parallel-backend",
        choices=("process", "thread"),
        default=os.environ.get("PROTEIN_PREP_BENCHMARK_BACKEND", "process"),
        help="Parallel executor backend used when --workers > 1.",
    )
    parser.add_argument(
        "--add-missing-residues",
        action="store_true",
        help="Allow PDBFixer to build missing residue segments. Off by default because blind loop building can hurt docking receptor chemistry.",
    )
    parser.add_argument(
        "--run-minimize",
        action="store_true",
        help="Run the restrained OpenMM protein-only minimization feasibility check. This can be slow on large receptors.",
    )
    parser.add_argument("--run-redock", action="store_true", help="Run optional Vina redocking when obabel, Meeko ligand prep, and vina are available.")
    parser.add_argument(
        "--run-autodock4zn",
        action="store_true",
        help="Run AutoDock4Zn receptor/map generation for eligible Zn sites when required tools are installed.",
    )
    parser.add_argument(
        "--skip-water-analysis",
        action="store_true",
        help="Skip water-policy Meeko exports for receptor-prep-only smoke runs.",
    )
    parser.add_argument(
        "--publication-mode",
        action="store_true",
        help="Require strict publication-readiness gates: redocking and restrained minimization when geometry cleanup was needed.",
    )
    parser.add_argument(
        "--publication-redock-rmsd-max-a",
        type=float,
        default=2.5,
        help="Maximum accepted redocking RMSD in Angstrom for --publication-mode.",
    )
    parser.add_argument("--box-size", type=float, default=22.0, help="Cubic redocking box size in Angstrom.")
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--log-level", default=os.environ.get("LOG_LEVEL", "INFO"))
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, str(args.log_level).upper(), logging.INFO),
        format="%(levelname)s:%(message)s",
    )
    pdb_ids = _parse_pdb_ids(args.pdb_ids, args.pdb_id_file)
    if args.limit and args.limit > 0:
        pdb_ids = pdb_ids[: args.limit]

    workers = _effective_workers(
        requested=int(args.workers),
        target_count=len(pdb_ids),
        docking_enabled=bool(args.run_redock or args.publication_mode or args.run_autodock4zn),
    )
    rows, failures = _run_targets(
        pdb_ids,
        args,
        workers=workers,
        backend=str(args.parallel_backend),
    )
    write_outputs(rows, args.out_dir)
    if failures:
        (args.out_dir / "failures.json").write_text(
            json.dumps(failures, indent=2, sort_keys=True),
            encoding="utf-8",
        )
    logging.info(
        "[benchmark] wrote %d rows to %s failures=%d workers=%d",
        len(rows),
        args.out_dir,
        len(failures),
        workers,
    )
    return 0 if rows else 1


def _parse_pdb_ids(raw: str | None, pdb_file: Path | None) -> list[str]:
    ids: list[str] = []
    if raw:
        ids.extend(part.strip().upper() for part in raw.split(",") if part.strip())
    if pdb_file:
        with pdb_file.open("r", encoding="utf-8") as handle:
            ids.extend(line.strip().upper() for line in handle if line.strip())
    return ids or list(DEFAULT_PDB_IDS)


def _effective_workers(
    *,
    requested: int,
    target_count: int,
    docking_enabled: bool,
) -> int:
    if target_count <= 1:
        return 1
    workers = max(1, min(int(requested), target_count, MAX_TARGET_WORKERS))
    if workers > 1 and docking_enabled:
        logging.warning(
            "[benchmark.parallel] forcing workers=1 for docking/autodock4zn modes to avoid exceeding the 32-core docking cap"
        )
        return 1
    return workers


def _run_targets(
    pdb_ids: Sequence[str],
    args: argparse.Namespace,
    *,
    workers: int,
    backend: str,
) -> tuple[list[BenchmarkRow], dict[str, str]]:
    tasks = [_task_payload(pdb_id, args) for pdb_id in pdb_ids]
    if workers <= 1:
        return _run_targets_serial(tasks)
    return _run_targets_parallel(tasks, workers=workers, backend=backend)


def _run_targets_serial(tasks: Sequence[dict[str, Any]]) -> tuple[list[BenchmarkRow], dict[str, str]]:
    rows: list[BenchmarkRow] = []
    failures: dict[str, str] = {}
    for task in tasks:
        pdb_id, row, error = _run_one_task(task)
        if row is not None:
            rows.append(row)
        else:
            failures[pdb_id] = error
    return rows, failures


def _run_targets_parallel(
    tasks: Sequence[dict[str, Any]],
    *,
    workers: int,
    backend: str,
) -> tuple[list[BenchmarkRow], dict[str, str]]:
    rows_by_index: dict[int, BenchmarkRow] = {}
    failures: dict[str, str] = {}
    executor_cls = ThreadPoolExecutor if backend == "thread" else ProcessPoolExecutor
    logging.info(
        "[benchmark.parallel] targets=%d workers=%d backend=%s",
        len(tasks),
        workers,
        backend,
    )
    with executor_cls(max_workers=workers) as executor:
        futures = {executor.submit(_run_one_task, task): idx for idx, task in enumerate(tasks)}
        for future in as_completed(futures):
            idx = futures[future]
            pdb_id = str(tasks[idx]["pdb_id"])
            try:
                target, row, error = future.result()
            except Exception as exc:
                failures[pdb_id] = str(exc)
                logging.warning("[benchmark] target failed pdb=%s err=%s", pdb_id, exc)
                continue
            if row is not None:
                rows_by_index[idx] = row
                logging.info("[benchmark] target complete pdb=%s", target)
            else:
                failures[target] = error
                logging.warning("[benchmark] target failed pdb=%s err=%s", target, error)
    rows = [rows_by_index[idx] for idx in range(len(tasks)) if idx in rows_by_index]
    return rows, failures


def _run_one_task(task: Mapping[str, Any]) -> tuple[str, BenchmarkRow | None, str]:
    pdb_id = str(task["pdb_id"])
    try:
        logging.info("[benchmark] target=%s", pdb_id)
        row = run_one(
            pdb_id,
            Path(task["out_dir"]),
            add_missing_residues=bool(task["add_missing_residues"]),
            run_minimize=bool(task["run_minimize"]),
            run_redock=bool(task["run_redock"]),
            run_autodock4zn=bool(task["run_autodock4zn"]),
            run_water_analysis=bool(task["run_water_analysis"]),
            publication_mode=bool(task["publication_mode"]),
            publication_redock_rmsd_max=float(task["publication_redock_rmsd_max"]),
            box_size=float(task["box_size"]),
            exhaustiveness=int(task["exhaustiveness"]),
        )
        return pdb_id, row, ""
    except Exception as exc:
        return pdb_id, None, str(exc)


def _task_payload(pdb_id: str, args: argparse.Namespace) -> dict[str, Any]:
    return {
        "pdb_id": pdb_id,
        "out_dir": str(args.out_dir),
        "add_missing_residues": bool(args.add_missing_residues),
        "run_minimize": bool(args.run_minimize),
        "run_redock": bool(args.run_redock or args.publication_mode),
        "run_autodock4zn": bool(args.run_autodock4zn),
        "run_water_analysis": not bool(args.skip_water_analysis),
        "publication_mode": bool(args.publication_mode),
        "publication_redock_rmsd_max": float(args.publication_redock_rmsd_max_a),
        "box_size": float(args.box_size),
        "exhaustiveness": int(args.exhaustiveness),
    }


__all__ = ["build_parser", "main"]
