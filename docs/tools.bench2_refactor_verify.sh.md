# tools/bench2_refactor_verify.sh

## 2026-03-04

- Added a single-command bench2-fast verification loop for refactor work.
- Runs baseline and candidate benchmarks, captures each run, and compares results via:
  - `tools/benchmark_bench2_compare.py capture`
  - `tools/benchmark_bench2_compare.py compare`
- Supports direct mode and local simulated-distributed mode:
  - direct: `python main.py -bench2 -fast --run-id <id>`
  - simulated distributed: `python tools/simulate_slurm_array.py ... --main-args "-bench2 -fast"`
- Added safety guard for simulated distributed runs:
  - requires `tasks * cpus-per-task <= 32`.
- Supports `--skip-runs` to operate on existing run artifacts only.

### Example

```bash
tools/bench2_refactor_verify.sh \
  --baseline-run-id bench2_base_20260304_a \
  --candidate-run-id bench2_opt_20260304_a \
  --simulate-distributed \
  --tasks 3 \
  --cpus-per-task 4
```
