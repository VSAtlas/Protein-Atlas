# HPC Throughput Benchmarking

Atlas throughput work should use a tiered benchmark ladder. The goal is to keep development fast while reserving final performance claims for real Slurm runs.

## Tiers

| Tier | Command | Purpose | Expected claim strength |
| --- | --- | --- | --- |
| Micro | `atlas throughput bench --profile micro --mode local` | Fast scheduler and SCORCH/reporting iteration on a tiny `test_library_20` workload. | Development only. |
| Smoke | `atlas throughput bench --profile smoke --mode local` | Bench2-compatible local correctness and throughput sanity. | Local-machine comparison. |
| Slurm simulation | `atlas throughput bench --profile smoke --mode slurm-sim --array 0-2%3 --cpus-per-task 4` | Array sharding, chunking, CPU caps, and finalizer-compatible artifact flow. | Orchestration only. |
| Slurm planning | `atlas throughput bench --profile medium --mode slurm-dry-run --array 0-31%4 --cpus-per-task 8` | Submission preview, `%` concurrency, exported args, finalizer dependency, and CPU budget visibility. | Launch readiness only. |
| Real cluster | Use the reviewed `slurm_submit_preview.txt`, then launch with explicit benchmark `--main-args`, for example `atlas slurm submit --run-id <RUN_ID> --array 0-31%4 --cpus-per-task 8 --main-args "-bench2 --run-id {run_id}"`. For a low-SU bench2-fast utilization canary, preview with `atlas slurm submit --bench2-canary --run-id <RUN_ID> --dry-run`. | Queue/controller/filesystem/accounting behavior under real Slurm. | Performance evidence. |

## Fixed Seeds

`atlas throughput bench` sets a benchmark-only Vina seed, defaulting to `1337`. Override with `--seed`.

This seed is exported with `ATLAS_THROUGHPUT_BENCH=1` and `ATLAS_THROUGHPUT_BENCH_SEED`; Atlas consumes it only for throughput-benchmark launches. Production docking profiles and ordinary `atlas slurm submit` calls do not receive fixed Vina seeds by default.

Fixed seeds improve candidate-vs-baseline comparability, but exact reproducibility still requires identical receptor/ligand inputs, Vina version, CPU/threading, parameters, and platform.

## Cached Prep Micro

The `micro` profile is intended to exercise scheduling and reporting, not protein prep. By default it requires cached HOLO receptor prep artifacts and fails fast when they are missing.

Use `--prep-cache-run-id <RUN_ID>` to stage receptor prep from a prior bench2-compatible run into the new micro run. Use `--allow-fresh-prep` only when deliberately testing prep behavior.

The harness writes `prep_cache_summary.json` beside the other throughput artifacts.

## Slurm Launch Checklist

Before moving to a Slurm machine:

- run `atlas throughput bench --profile medium --mode slurm-dry-run --array 0-31%4 --cpus-per-task 8`;
- confirm `slurm_submit_preview.txt` has one array submission, `%` concurrency, exported `ATLAS_MAIN_ARGS`, and a finalizer dependency;
- for whole-node utilization checks, use `atlas slurm submit --bench2-canary --run-id <RUN_ID> --dry-run`, then launch the same command without `--dry-run`; this preset uses three exclusive whole-node array tasks for 15 minutes and exports `ATLAS_USE_NODE_CPUS=1` so Atlas prefers `SLURM_CPUS_ON_NODE`;
- for local/non-Slurm testing, keep `concurrent_cpu_budget <= 32`; on the cluster, choose concurrency and `--cpus-per-task` for the actual allocation and queue policy;
- on the cluster, launch with explicit `--main-args`, monitor with `atlas slurm progress <RUN_ID>`, and let the finalizer complete;
- after completion, run `atlas analysis report <RUN_ID> --status-html`, `atlas analysis throughput integrity <RUN_ID>`, and `atlas artifacts measure --run-id <RUN_ID>`;
- check `outputs/manifests/<RUN_ID>/run_efficiency.json` for `sacct_available=true` and keep `outputs/manifests/<RUN_ID>/slurm_accounting/` with the raw `sacct` command/stdout/stderr.

## Slurm Readiness

Each throughput bench invocation writes:

`outputs/data/<RUN_ID>/throughput_bench/slurm_readiness.json`

This summarizes:

- profile, main args, fixed seed, array spec, `%` concurrency, and concurrent CPU budget;
- integrity and acceptance status when available;
- Slurm dry-run preview path;
- `run_efficiency.json` and `sacct` availability when a real Slurm run has completed;
- remaining gaps before cluster performance claims.

## Remaining Cluster-Only Gaps

Local simulation does not validate:

- queue wait or backfill behavior;
- Slurm controller throttles and site array limits;
- node heterogeneity;
- shared filesystem contention;
- real CPU and memory accounting;
- `sacct` efficiency and MaxRSS data.

Before treating a build as HPC-ready, run a real Slurm pass and collect:

- `atlas slurm progress <RUN_ID>`;
- finalizer completion;
- `atlas analysis throughput integrity <RUN_ID>`;
- acceptance against the current fastest bench2 baseline;
- artifact sizing with `atlas artifacts measure --run-id <RUN_ID>`;
- `run_efficiency.json` with `sacct_available=true`.

## References

- Slurm job arrays and `%` concurrency: https://slurm.schedmd.com/job_array.html
- Slurm accounting fields via `sacct`: https://slurm.schedmd.com/sacct.html
- Local executor development before cluster execution: https://docs.seqera.io/nextflow/executor
- Vina seed/reproducibility behavior: https://vina.scripps.edu/manual/
