# Parallel Agent Plan: Green CI + QoL CLI

**Goal:** `atlas dev verify --full --fix --smoke` passes on `main`, with regression coverage for primary QoL commands.

**Orchestration:** Headless `agent` processes in separate shells, each in an isolated `--worktree`. This document is the source of truth; the launcher is `tools/launch_parallel_ci_qol_agents.sh`.

**Logs:** `/stor/home/mpg2352/atlas_parallel_agents/<RUN_ID>/`

**Models:** `composer-2.5-fast` (default) or `composer-2.5` for integrator / hard refactors.

---

## Constraints (every agent)

- Edit only under this repo; never touch `tests/`, `ci/`, `input_pdbs/`, `docked/`.
- Prefer `chemdb/tests/` for new regression tests.
- Use `src.path_router` / `import src.path_router` — no hardcoded data roots.
- Do not read raw run logs; use `atlas debug` / `atlas status` patterns in tests only via fixtures.
- Finish with focused verify (see per-stream) and print `AGENT_DONE:<stream_id>` on success or `AGENT_FAILED:<stream_id>` on block.

---

## Current baseline (2026-05-23)

| Gate | Status | Notes |
|------|--------|-------|
| `tach check` | pass | |
| `ruff check .` | **fail** | ~10 issues (mostly unused imports; 6+ auto-fixable) |
| `mypy .` (full) | **fail** | e.g. `tools/finalize_distributed_run.py:334`; broader debt likely |
| `mypy` (scoped gate) | **fail** | Same finalize line in default `MYPY_TARGETS` |
| `vulture .` | **fail** | Many 60% confidence hits in `tools/` |
| `architecture_gate.py` | **fail** | New Xenon C/D/F blocks outside `tools/xenon_baseline.txt` |
| `python main.py --test -fast` | unknown | Run after static gates clean |

---

## Phase 0 — Orchestrator (once)

```bash
export ATLAS_AGENT_RUN_ID="$(date +%Y%m%d_%H%M%S)"
bash tools/launch_parallel_ci_qol_agents.sh launch
bash tools/launch_parallel_ci_qol_agents.sh status   # poll
```

Optional: watch a stream:

```bash
tail -f "/stor/home/mpg2352/atlas_parallel_agents/${ATLAS_AGENT_RUN_ID}/ci-ruff.log"
```

---

## Phase 1 — Parallel streams (isolated worktrees)

Each stream owns disjoint paths to minimize merge pain.

### CI streams (green gate)

| ID | Worktree branch | Owns | Verify before done |
|----|-----------------|------|-------------------|
| `ci-ruff` | `ci-ruff` | Repo-wide Ruff | `ruff check .` |
| `ci-mypy-gate` | `ci-mypy-gate` | `MYPY_TARGETS` in `tools/quality_gate.sh` + `tools/finalize_distributed_run.py` | `mypy --follow-imports skip` on those targets |
| `ci-mypy-analysis` | `ci-mypy-analysis` | `analysis/` (exclude generated templates if needed) | `mypy analysis --follow-imports skip` |
| `ci-vulture` | `ci-vulture` | `tools/*.py` dead code only | `vulture tools --min-confidence 80` or justified noqa + gate |
| `ci-xenon` | `ci-xenon` | `tools/xenon_baseline.txt` + **refactor preferred** for `src/cli/doctor.py`, `src/cli/main_pipeline.py` if small; else baseline file entries | `python tools/architecture_gate.py` |

**ci-xenon new violations (add to baseline only if refactor is out of scope):**

- `src/protein_prep/water_dude_benchmark.py`
- `src/protein_prep/geometry/constructive.py`
- `src/docking/active_preserving_cap.py`
- `src/post_docking/rescoring/scorch_provisional_cache.py`
- `src/post_docking/mmgbsa/mmgbsa_wang2016_production.py`
- `src/cli/planner_chunking.py`
- `src/cli/throughput_bench.py`
- `src/cli/main_pipeline.py`
- `src/cli/doctor.py`

### QoL streams (CLI regression)

| ID | Worktree branch | Owns | Verify before done |
|----|-----------------|------|-------------------|
| `qol-test-init` | `qol-test-init` | New `chemdb/tests/test_cli_init_new_run.py`; minimal `src/cli/qol_cli.py` only if required | `pytest chemdb/tests/test_cli_init_new_run.py -q` |
| `qol-test-status` | `qol-test-status` | `chemdb/tests/test_status_dashboard.py`, fixtures under `chemdb/tests/fixtures/` | `pytest chemdb/tests/test_status_dashboard.py -q` |
| `qol-test-doctor` | `qol-test-doctor` | `chemdb/tests/test_cli_doctor.py`, `src/cli/doctor.py` | `pytest chemdb/tests/test_cli_doctor.py -q` |
| `qol-test-help-dev` | `qol-test-help-dev` | `chemdb/tests/test_cli_help.py`; hermetic `atlas dev verify` wiring test | `pytest chemdb/tests/test_cli_help.py -q` |

**Do not run in parallel:** `qol-split-cli` (deferred Phase 3) — splits `src/cli/qol_cli.py` (~100KB) and collides with all QoL streams.

---

## Phase 2 — Integrator (sequential merges)

After all streams print `AGENT_DONE:*` (or you accept partial success):

```bash
bash tools/merge_parallel_agent_worktrees.sh
cd /stor/work/VDS_Beckham/atlas/code/protein_automation
atlas dev verify --full --fix --smoke
```

Merge order (low conflict → high):

1. `ci-ruff`
2. `ci-mypy-gate`
3. `ci-mypy-analysis`
4. `ci-vulture`
5. `ci-xenon`
6. `qol-test-init`
7. `qol-test-status`
8. `qol-test-doctor`
9. `qol-test-help-dev`

On conflict: prefer QoL test intent + CI fixes; re-run focused pytest for affected streams.

---

## Phase 3 — Deferred (post-green)

| ID | Purpose |
|----|---------|
| `qol-split-cli` | Split `src/cli/qol_cli.py` into `src/cli/qol/` package; keep `atlas` entrypoints stable |
| `qol-status-live-slurm` | Live Slurm polling in `atlas status --watch` (cluster validation) |
| `ci-mypy-repo-wide` | Full `mypy .` debt burn beyond gate targets |

---

## Agent prompt template

Each launcher invocation uses:

```text
You are stream <ID> for Atlas protein_automation.
Read docs/parallel_agent_plan_green_ci_qol.md section for <ID>.
<stream-specific instructions>
When finished, run the stream Verify command and print AGENT_DONE:<ID> or AGENT_FAILED:<ID>.
Append one entry to docs/agent_handoffs.md (append-only).
Do not commit unless changes are complete and tests pass.
```

---

## Success criteria

- [ ] `bash tools/quality_gate.sh --full --fix` exits 0
- [ ] `python tools/architecture_gate.py` exits 0
- [ ] `pytest chemdb/tests/test_cli_*.py chemdb/tests/test_status_dashboard.py -q` passes
- [ ] `atlas dev verify --full --fix --smoke` exits 0
- [ ] `docs/agent_handoffs.md` has one line per completed stream

---

## Monitoring checklist

```bash
# Running agents
pgrep -af 'agent --print' | grep protein_automation

# Done markers
grep -h 'AGENT_DONE\|AGENT_FAILED' "/stor/home/mpg2352/atlas_parallel_agents/${ATLAS_AGENT_RUN_ID}/"*.log

# Worktrees
git worktree list
```

---

## Token budget tips

- Run **8 Phase-1 agents** concurrently (4 CI + 4 QoL).
- Use **integrator** as a single `composer-2.5` agent with merge + full verify prompt.
- Re-launch only failed streams with `bash tools/launch_parallel_ci_qol_agents.sh launch ci-ruff`.
