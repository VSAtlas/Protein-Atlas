# Parallel Agent Plan: Simplify Rounds (6–9+ modules)

**Goal:** Shrink and de-complexify the largest/most-maintained modules without behavior changes. Each round runs in an isolated worktree and may spawn **nested headless subagents** (`agent --print --worktree ...`).

**Orchestration:** `tools/launch_parallel_simplify_agents.sh`

**Logs:** `/stor/home/mpg2352/atlas_parallel_agents/<RUN_ID>/`

**Models:** `composer-2.5-fast` (workers/subagents), `composer-2.5` (integrator merges)

**Gate after each wave:** `atlas dev verify --full --fix --smoke` (or focused verify if wave is narrow)

---

## Rules (every primary agent and subagent)

- Edit only this repo; never `tests/`, `ci/`, `input_pdbs/`, `docked/`.
- **No behavior changes** unless fixing an obvious bug; preserve public function names and CLI flags.
- Files must stay **under 100 KB**; split into submodules when larger.
- Use canonical imports: `cli`, `path_router`, `config`, `post_docking` — never new `src.*` imports in `analysis/`.
- Prefer extract-function / extract-module over clever abstraction.
- Add or extend `chemdb/tests/` only when splitting public entrypoints.
- Subagents: max **3 per primary**; each subagent gets its own `--worktree <stream>-sub<N>`.
- Primary agent merges sub-worktrees into its stream branch before `AGENT_DONE`.
- Print `AGENT_DONE:<stream>` or `AGENT_FAILED:<stream>`; append one line to `docs/agent_handoffs.md`.

### Subagent spawn template (primary agents)

```bash
REPO=/stor/work/VDS_Beckham/atlas/code/protein_automation
STREAM=<stream-id>
SUB=<sub-id>
agent --print --force --trust --yolo --model composer-2.5-fast \
  --workspace "$REPO" --worktree "${STREAM}-${SUB}" \
  "<paste sub-task prompt from table below>"
```

Poll: `tail -f .../logs/${STREAM}-${SUB}.log` or `grep AGENT_DONE ...`

---

## Why these targets

| Module | Lines | Pain |
|--------|------:|------|
| `analysis/reporting/run_report_core.py` | 4436 | Report monolith |
| `analysis/reporting/heatmap_html.py` | 2808 | HTML + client payload |
| `analysis/dud_eval_core/orchestrate.py` | 3246 | DUD eval orchestration |
| `src/cli/qol_cli.py` | 2535 | QoL dispatcher (over 100 KB) |
| `src/cli/run_manifest_runtime.py` | 2125 | Manifest I/O |
| `src/cli/postrun_hooks_runtime.py` | 1844 | Post-run hooks |
| `src/cli/distributed_runtime_claim_loop.py` | 1222 | Distributed claim loop |
| `src/cli/distributed_chunk_planner.py` | 983 | Chunk planning |
| `src/cli/status_dashboard.py` | 1157 | Status / Slurm dashboard |

Many are also listed in `tools/xenon_baseline.txt` — removing a path from that file after refactor is a success metric.

---

## Nine primary rounds (Wave 1 = all parallel-safe)

| ID | Worktree | Owns | Target outcome | Verify |
|----|----------|------|----------------|--------|
| `simp-01-report` | `simp-01-report` | `analysis/reporting/run_report_core.py` → `analysis/reporting/run_report/` package | No file >800 lines; stable imports | `pytest chemdb/tests/test_run_report*.py -q` |
| `simp-02-heatmap` | `simp-02-heatmap` | `heatmap_html.py`, `heatmap_html_clustergrammer_template.py` | Template vs logic split | `pytest chemdb/tests/test_heatmap*.py chemdb/tests/test_report_html*.py -q` |
| `simp-03-dudeval` | `simp-03-dudeval` | `analysis/dud_eval_core/orchestrate.py` | Orchestration vs reporting split | `pytest chemdb/tests/test_dud_eval*.py -q` |
| `simp-04-qol` | `simp-04-qol` | `src/cli/qol_cli.py` → `src/cli/qol/` | Package <100 KB/file; `atlas` entry unchanged | `pytest chemdb/tests/test_cli_*.py -q` |
| `simp-05-manifest` | `simp-05-manifest` | `run_manifest_runtime.py`, `run_manifest_support.py` | Read/write helpers extracted | `pytest chemdb/tests/test_run_manifest.py -q` |
| `simp-06-postrun` | `simp-06-postrun` | `postrun_hooks_runtime.py`, `postrun_hooks_support.py` | Hook phases as submodules | `pytest chemdb/tests/test_scorch_* chemdb/tests/test_run_lifecycle*.py -q` |
| `simp-07-dist-plan` | `simp-07-dist-plan` | `distributed_chunk_planner.py`, `planner_chunking.py`, `planner_manifest.py` | Planner-only package surface | `pytest chemdb/tests/test_distributed_chunk_planning.py -q` |
| `simp-08-dist-run` | `simp-08-dist-run` | `distributed_runtime_claim_loop.py`, `distributed_chunk_runtime.py` | Claim vs completion split | `pytest chemdb/tests/test_distributed_chunk_regressions.py -q` |
| `simp-09-status` | `simp-09-status` | `status_dashboard.py`, `doctor.py` | Panel builders extracted | `pytest chemdb/tests/test_status_dashboard.py chemdb/tests/test_cli_doctor.py -q` |

**Conflict note:** Run **simp-04-qol** alone if other waves touch `src/cli/qol_cli.py` imports. Waves 5–9 all touch `src/cli/` but different files — merge order matters (see below).

---

## Subagent task breakdown (optional per primary)

### simp-01-report
| Sub | Task |
|-----|------|
| `sub-a` | Extract HTML/table writers from `run_report_core.py` |
| `sub-b` | Extract data assembly / manifest readers |
| `sub-c` | Wire `__init__.py` re-exports; delete dead helpers |

### simp-02-heatmap
| Sub | Task |
|-----|------|
| `sub-a` | Move clustergrammer template strings to template module only |
| `sub-b` | Extract payload builders + asset resolution |
| `sub-c` | Tests: ensure report HTML still embeds heatmap block |

### simp-03-dudeval
| Sub | Task |
|-----|------|
| `sub-a` | Extract CLI/main() and argparse from orchestrate |
| `sub-b` | Extract target discovery loop |
| `sub-c` | Extract metrics/report emission |

### simp-04-qol
| Sub | Task |
|-----|------|
| `sub-a` | `qol/init_demo.py` — init, demo, smoke |
| `sub-b` | `qol/status_runs.py` — status, runs, doctor dispatch |
| `sub-c` | `qol/dev_slurm.py` — dev, slurm, throughput; thin `qol_cli.py` dispatcher |

### simp-05-manifest
| Sub | Task |
|-----|------|
| `sub-a` | YAML load/save atomic helpers |
| `sub-b` | Protein/combo mutation helpers |
| `sub-c` | Re-export compatibility aliases for tests |

### simp-06-postrun
| Sub | Task |
|-----|------|
| `sub-a` | SCORCH barrier + queue drain |
| `sub-b` | Report + heatmap hooks |
| `sub-c` | Retention + throughput integrity hooks |

### simp-07-dist-plan
| Sub | Task |
|-----|------|
| `sub-a` | Ligand enumeration + chunk key generation |
| `sub-b` | Manifest/sample planner modes |
| `sub-c` | Validation + telemetry only |

### simp-08-dist-run
| Sub | Task |
|-----|------|
| `sub-a` | Claim loop + steal/rebalance |
| `sub-b` | Chunk result read/write + terminal states |
| `sub-c` | Integration with `distributed_context.py` (read-only unless required) |

### simp-09-status
| Sub | Task |
|-----|------|
| `sub-a` | Slurm parser fixtures path + `parse_squeue`/`parse_sacct` |
| `sub-b` | HTML renderer + ETA cache |
| `sub-c` | Doctor display helpers (no message churn) |

---

## Execution waves

### Wave 1 (9 parallel — reporting + CLI + runtime disjoint files)

```bash
export ATLAS_AGENT_RUN_ID="simplify_wave1_$(date +%Y%m%d)"
bash tools/launch_parallel_simplify_agents.sh launch wave1
bash tools/launch_parallel_simplify_agents.sh status
```

Streams: all nine IDs above.

### Wave 2 (optional follow-ups — second-pass debt)

| ID | Focus |
|----|--------|
| `simp-10-schema` | `master_schema_export.py` (1633 lines) |
| `simp-11-ml` | `train_classifier.py`, `audit_suite.py` |
| `simp-12-scorch` | `rescoring_scorch.py`, `scorch_orchestration.py` |

```bash
bash tools/launch_parallel_simplify_agents.sh launch wave2
```

### Integrator (after each wave)

```bash
bash tools/merge_parallel_simplify_worktrees.sh wave1
cd /stor/work/VDS_Beckham/atlas/code/protein_automation
atlas dev verify --full --fix --smoke
```

**Merge order (wave 1):** 01 → 02 → 03 → 09 → 04 → 05 → 06 → 07 → 08

---

## Success metrics (per stream)

- [ ] Primary module reduced by **≥25% lines** OR split into package with no file >800 lines / 100 KB
- [ ] `python tools/architecture_gate.py` still passes
- [ ] Stream verify pytest passes
- [ ] Removed baselined paths from `tools/xenon_baseline.txt` when Xenon clean
- [ ] `ruff check` on touched paths clean

---

## Token-burn tips

- Launch **wave1** (9 agents) overnight; use `composer-2.5-fast` for primaries/subagents.
- Use **one integrator** `composer-2.5` merge + verify per wave.
- Re-launch single stream: `bash tools/launch_parallel_simplify_agents.sh launch simp-04-qol`
- Skip wave2 until wave1 verify is green.

---

## Monitoring

```bash
grep -h 'AGENT_DONE\|AGENT_FAILED' "/stor/home/mpg2352/atlas_parallel_agents/${ATLAS_AGENT_RUN_ID}/"*.log
git worktree list | grep simp-
```
