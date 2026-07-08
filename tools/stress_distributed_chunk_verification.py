#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Synthetic high-scale stress for distributed chunk verification and "
            "coverage snapshot writes. No docking engines are executed."
        )
    )
    parser.add_argument("--run-id", default="stress_chunk_verify")
    parser.add_argument(
        "--root",
        default="/stor/home/mpg2352/stress_chunk_verify",
        help="Scratch root for generated synthetic artifacts.",
    )
    parser.add_argument("--combos", type=int, default=80)
    parser.add_argument("--chunks-per-combo", type=int, default=12)
    parser.add_argument("--ligands-per-chunk", type=int, default=64)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument(
        "--fault-missing-summary-every",
        type=int,
        default=0,
        help="If >0, every Nth combo omits docking_score_summary.csv.",
    )
    parser.add_argument(
        "--fault-nul-summary-every",
        type=int,
        default=0,
        help="If >0, every Nth combo writes a summary containing NUL bytes.",
    )
    parser.add_argument("--keep", action="store_true")
    parser.add_argument("--json-out", default="")
    return parser.parse_args()


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _import_main(repo_root: Path):
    src_root = repo_root / "src"
    if str(src_root) not in sys.path:
        sys.path.insert(0, str(src_root))
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    import main  # noqa: PLC0415

    return main


def _write_summary_csv(path: Path, run_id: str, ligands: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["run_id,Ligand,stage1\n"]
    for lig in ligands:
        lines.append(f"{run_id},{lig}.pdbqt,-7.0\n")
    path.write_text("".join(lines), encoding="utf-8")


def _write_summary_nul_csv(path: Path, run_id: str, ligands: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    blob = f"run_id,Ligand,stage1\n{run_id},{ligands[0]}.pdbqt,\x00\n".encode("utf-8")
    path.write_bytes(blob)


def _combo_id(index: int) -> str:
    return f"P{index:04d}"


def _ph_for_index(index: int) -> str:
    values = ("pH6_5", "pH6_9", "pH7_0", "pH7_4", "pH7_8")
    return values[index % len(values)]


def main() -> int:
    ns = _parse_args()
    repo_root = _repo_root()
    mod = _import_main(repo_root)

    scratch = Path(str(ns.root)).resolve()
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True, exist_ok=True)

    cfg = {
        "OVERALL_DIR": str(scratch),
        "DOCKED_DIR": str(scratch / "docked"),
        "RUN_ID": str(ns.run_id),
    }
    run_id = str(ns.run_id)

    work_items: list[dict[str, Any]] = []
    build_started = time.perf_counter()
    for ci in range(max(1, int(ns.combos))):
        pdb_id = _combo_id(ci)
        variant = "HOLO"
        ph = _ph_for_index(ci)
        combo_dir = Path(cfg["DOCKED_DIR"]) / run_id / pdb_id / variant / ph
        summary_path = combo_dir / "docking_score_summary.csv"

        ligands = [
            f"decoys_final_{ci:04d}_{li:05d}"
            for li in range(max(1, int(ns.chunks_per_combo)) * max(1, int(ns.ligands_per_chunk)))
        ]
        missing_every = int(ns.fault_missing_summary_every or 0)
        nul_every = int(ns.fault_nul_summary_every or 0)
        if missing_every <= 0 or ((ci + 1) % missing_every != 0):
            if nul_every > 0 and ((ci + 1) % nul_every == 0):
                _write_summary_nul_csv(summary_path, run_id, ligands)
            else:
                _write_summary_csv(summary_path, run_id, ligands)

        chunk_size = max(1, int(ns.ligands_per_chunk))
        for chunk_idx in range(max(1, int(ns.chunks_per_combo))):
            start = chunk_idx * chunk_size
            chunk_ligs = ligands[start : start + chunk_size]
            work_items.append(
                {
                    "chunk_id": f"{pdb_id}_{ph}_{chunk_idx:03d}",
                    "pdb_id": pdb_id,
                    "variant": variant,
                    "ph": ph,
                    "chunk_ligs": chunk_ligs,
                }
            )

    counters = {
        "chunks_total": len(work_items),
        "chunks_ok": 0,
        "chunks_degraded": 0,
        "chunks_fatal": 0,
        "exceptions": 0,
    }

    def _run_item(item: dict[str, Any]) -> tuple[str, str]:
        ok, reason = mod._verify_chunk_combo_outputs(  # type: ignore[attr-defined]
            cfg,
            run_id=run_id,
            pdb_id=item["pdb_id"],
            variant_label=item["variant"],
            ph_tag=item["ph"],
            library_name="",
            chunk_ligand_bases=item["chunk_ligs"],
        )
        expected = {mod._chunk_ligand_key(x) for x in item["chunk_ligs"]}  # type: ignore[attr-defined]
        expected = {x for x in expected if x}
        mod._update_combo_coverage_snapshot(  # type: ignore[attr-defined]
            cfg,
            run_id=run_id,
            pdb_id=item["pdb_id"],
            variant_label=item["variant"],
            ph_tag=item["ph"],
            library_name="stress",
            expected_keys=expected,
            docking_keys=(set() if reason != "ok" else set(expected)),
        )
        if not ok:
            return "fatal", reason
        if reason == "ok":
            return "ok", reason
        return "degraded", reason

    with ThreadPoolExecutor(max_workers=max(1, int(ns.workers))) as pool:
        futures = [pool.submit(_run_item, item) for item in work_items]
        for fut in as_completed(futures):
            try:
                status, _reason = fut.result()
            except Exception:
                counters["exceptions"] += 1
                continue
            if status == "ok":
                counters["chunks_ok"] += 1
            elif status == "degraded":
                counters["chunks_degraded"] += 1
            else:
                counters["chunks_fatal"] += 1

    elapsed = time.perf_counter() - build_started
    out = {
        "run_id": run_id,
        "root": str(scratch),
        "elapsed_sec": round(elapsed, 3),
        "combos": int(ns.combos),
        "chunks_per_combo": int(ns.chunks_per_combo),
        "ligands_per_chunk": int(ns.ligands_per_chunk),
        "workers": int(ns.workers),
        "fault_missing_summary_every": int(ns.fault_missing_summary_every),
        "fault_nul_summary_every": int(ns.fault_nul_summary_every),
        "results": counters,
    }
    print(json.dumps(out, indent=2, sort_keys=True))

    if ns.json_out:
        json_out = Path(str(ns.json_out))
        json_out.parent.mkdir(parents=True, exist_ok=True)
        json_out.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")

    if not ns.keep:
        shutil.rmtree(scratch, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
