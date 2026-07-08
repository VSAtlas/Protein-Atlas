from __future__ import annotations

import json
import random
import shutil
import shlex
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class WorkloadTask:
    task_id: str
    cores: int
    command: str
    expected_sec: float | None = None
    tags: list[str] | None = None



def _busy_loop_command(duration_sec: float, cores: int = 1) -> str:
    secs = max(0.01, float(duration_sec))
    workers = max(1, int(cores))
    script = (
        "import time\\n"
        f"end=time.perf_counter()+{secs:.6f}\\n"
        "x=0\\n"
        "while time.perf_counter()<end:\\n"
        "    x+=1\\n"
    )
    script_q = script.replace("\\", "\\\\").replace("'", "\\'")
    return (
        f"{shlex.quote(sys.executable)} -c \"import builtins,os,multiprocessing as mp; "
        f"script='{script_q}'; "
        f"workers=int(os.environ.get('UTIL_BENCH_GRANTED_WORKERS','{workers}') or '{workers}'); "
        "workers=max(1,workers); "
        "runner=getattr(builtins,'exec'); "
        "procs=[mp.Process(target=runner,args=(script,)) for _ in range(workers)]; "
        "[p.start() for p in procs]; [p.join() for p in procs]\""
    )



def _sample_cores(rng: random.Random, core_weights: dict[int, float]) -> int:
    items = sorted((int(k), float(v)) for k, v in core_weights.items() if int(k) > 0)
    if not items:
        return 1
    total = sum(w for _, w in items)
    if total <= 0:
        return items[0][0]
    pick = rng.random() * total
    cur = 0.0
    for cores, weight in items:
        cur += weight
        if pick <= cur:
            return cores
    return items[-1][0]



def generate_synthetic_workload(
    *,
    n_tasks: int,
    seed: int = 0,
    mean_sec: float = 1.0,
    sigma: float = 0.55,
    outlier_fraction: float = 0.1,
    outlier_scale: float = 3.0,
    core_weights: dict[int, float] | None = None,
    run_name: str = "synthetic_util_bench",
    notes: str = "",
) -> dict[str, Any]:
    """Generate a synthetic CPU-bound workload manifest."""
    rng = random.Random(seed)
    weights = core_weights or {1: 0.65, 2: 0.25, 4: 0.1}

    tasks: list[dict[str, Any]] = []
    mu = max(0.01, mean_sec)

    for idx in range(int(n_tasks)):
        dur = rng.lognormvariate(mu=0.0, sigma=max(0.01, sigma)) * mu
        if rng.random() < max(0.0, min(1.0, outlier_fraction)):
            dur *= max(1.0, outlier_scale)
        dur = max(0.05, float(dur))
        cores = _sample_cores(rng, weights)
        task_id = f"t{idx:04d}"
        tasks.append(
            {
                "id": task_id,
                "cores": int(cores),
                "command": _busy_loop_command(dur, cores=cores),
                "expected_sec": round(dur, 6),
                "tags": ["synthetic", f"cores:{cores}"],
            }
        )

    return {
        "run": {
            "name": run_name,
            "seed": int(seed),
            "notes": notes,
        },
        "tasks": tasks,
    }



def generate_stage_collapse_workload(
    *,
    proteins: int,
    stage1_tasks: int,
    stage2_keep_frac: float = 0.25,
    stage3_keep_frac: float = 0.10,
    seed: int = 0,
    core_weights: dict[int, float] | None = None,
    stage1_mean_sec: float = 0.4,
    stage2_mean_sec: float = 0.3,
    stage3_mean_sec: float = 0.25,
    run_name: str = "stage_collapse_util_bench",
) -> dict[str, Any]:
    """Generate a pipeline-like workload with stage cardinality collapse."""
    rng = random.Random(seed)
    weights = core_weights or {1: 0.7, 2: 0.2, 4: 0.1}
    s1 = max(1, int(stage1_tasks))
    s2 = max(1, int(round(s1 * max(0.01, stage2_keep_frac))))
    # Interpret stage3_keep_frac against stage1 for stable shaping regardless of stage2 rounding.
    s3 = max(1, int(round(s1 * max(0.01, stage3_keep_frac))))

    tasks: list[dict[str, Any]] = []
    for p in range(max(1, int(proteins))):
        pid = f"P{p:02d}"
        for stage_name, n_tasks, mean_sec in (
            ("stage1", s1, stage1_mean_sec),
            ("stage2", s2, stage2_mean_sec),
            ("stage3", s3, stage3_mean_sec),
        ):
            for t in range(n_tasks):
                dur = max(0.05, rng.lognormvariate(mu=0.0, sigma=0.35) * mean_sec)
                cores = _sample_cores(rng, weights)
                tid = f"{pid}_{stage_name}_{t:04d}"
                tasks.append(
                    {
                        "id": tid,
                        "cores": int(cores),
                        "command": _busy_loop_command(dur, cores=cores),
                        "expected_sec": round(dur, 6),
                        "tags": [pid, stage_name, f"cores:{cores}"],
                    }
                )

    return {
        "run": {
            "name": run_name,
            "seed": int(seed),
            "notes": (
                f"stage-collapse proteins={proteins} stage1={s1} "
                f"stage2={s2} stage3={s3}"
            ),
        },
        "tasks": tasks,
    }


def validate_workload(workload: dict[str, Any]) -> None:
    if not isinstance(workload, dict):
        raise ValueError("workload must be an object")

    tasks = workload.get("tasks")
    if not isinstance(tasks, list) or not tasks:
        raise ValueError("workload.tasks must be a non-empty list")

    seen: set[str] = set()
    for row in tasks:
        if not isinstance(row, dict):
            raise ValueError("task entry must be an object")
        task_id = str(row.get("id", "")).strip()
        if not task_id:
            raise ValueError("task.id is required")
        if task_id in seen:
            raise ValueError(f"duplicate task id: {task_id}")
        seen.add(task_id)

        cores = int(row.get("cores", 0))
        if cores < 1:
            raise ValueError(f"task {task_id} has invalid cores={cores}")

        cmd = str(row.get("command", "")).strip()
        if not cmd:
            raise ValueError(f"task {task_id} has empty command")



def load_workload(path: str | Path) -> dict[str, Any]:
    p = Path(path)
    obj = json.loads(p.read_text(encoding="utf-8"))
    validate_workload(obj)
    return obj



def write_workload_json(path: str | Path, workload: dict[str, Any]) -> Path:
    validate_workload(workload)
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(workload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return p



def prepare_synthetic_ligand_fixture(
    *,
    repo_root: str | Path,
    out_root: str | Path,
    new_pdb_id: str,
    source_pdb_id: str = "TEST",
    source_library_subdir: str = "test_library_10",
    clone_library_subdir: str | None = None,
    n_ligands: int = 16,
) -> dict[str, Any]:
    """
    Create a synthetic local fixture derived from TEST/test_library_10.

    Returns a mapping payload suitable for TEST_LIBRARY_MAP injection at runtime.
    """
    root = Path(repo_root)
    out = Path(out_root)
    clone_subdir = clone_library_subdir or f"{source_library_subdir}_clone_{new_pdb_id.lower()}"

    src_pdb = root / "input_pdbs" / f"{source_pdb_id}.pdb"
    if not src_pdb.exists():
        fixture_pdb = (
            root
            / "chemdb"
            / "tests"
            / "fixtures"
            / "input_pdbs"
            / f"{source_pdb_id}.pdb"
        )
        if fixture_pdb.exists():
            src_pdb = fixture_pdb
    if not src_pdb.exists():
        raise FileNotFoundError(f"missing source PDB: {src_pdb}")

    src_lib = root / "prepped_ligands" / source_library_subdir
    if not src_lib.is_dir():
        fixture_lib = (
            root
            / "chemdb"
            / "tests"
            / "fixtures"
            / "prepped_ligands"
            / source_library_subdir
        )
        if fixture_lib.is_dir():
            src_lib = fixture_lib
    if not src_lib.is_dir():
        raise FileNotFoundError(f"missing source library: {src_lib}")

    dst_input = out / "input_pdbs"
    dst_prepped = out / "prepped_ligands" / clone_subdir
    dst_input.mkdir(parents=True, exist_ok=True)
    dst_prepped.mkdir(parents=True, exist_ok=True)

    renamed_pdb = dst_input / f"{new_pdb_id.upper()}.pdb"
    shutil.copy2(src_pdb, renamed_pdb)

    copied = 0
    for lig in sorted(src_lib.glob("*.pdbqt")):
        shutil.copy2(lig, dst_prepped / lig.name)
        copied += 1
        if copied >= max(1, int(n_ligands)):
            break
    if copied == 0:
        raise FileNotFoundError(f"no pdbqt ligands found in {src_lib}")

    return {
        "input_dir": str(dst_input),
        "prepped_ligands_dir": str(out / "prepped_ligands"),
        "source_pdb_id": source_pdb_id.upper(),
        "new_pdb_id": new_pdb_id.upper(),
        "source_library_subdir": source_library_subdir,
        "clone_library_subdir": clone_subdir,
        "test_library_map": {
            source_pdb_id.upper(): source_library_subdir,
            new_pdb_id.upper(): clone_subdir,
        },
        "copied_ligands": copied,
    }
