# Professionalism Audit (Core Pipeline)

Date: 2026-03-02
Scope: `main.py`, `src/`, `tools/`, `docs/` (core pipeline only)
Goal: Improve maintainability and engineering hygiene without changing pipeline behavior.

## Implemented in this pass

1. CI quality workflow
- Added `.github/workflows/quality.yml`.
- Runs consistent lint/type gates on push and pull request.

2. Standardized quality entrypoint
- Added `tools/quality_gate.sh` as one command for scoped Ruff + MyPy checks.
- Added optional `--smoke` to include `python main.py --test -fast`.

3. Local code simplification
- Simplified duplicated command/log emission logic in `tools/simulate_slurm_array.py` by centralizing command construction/logging helpers.
- No CLI or runtime behavior changed.

4. Editor consistency baseline
- Added `.editorconfig` for line endings, trailing whitespace trimming, and Python/YAML/Markdown indentation consistency.

5. Hotspot visibility tool
- Added `tools/repo_hotspots.py` (+ docs) to rank large Python modules by bytes/lines so decomposition work is prioritized by objective size data.

## High-priority gaps still present

1. Oversized orchestration files
- `main.py` and `benchmark_mode.py` remain very large and hard to reason about.
- Risk: slower onboarding, harder regression isolation, high merge/conflict pressure.
- Recommended next move: continue extracting runtime responsibilities into `src/cli/*` modules while preserving existing function boundaries and outputs.

2. Mixed layout and import surface
- Legacy root modules and `src/*` packages coexist.
- Risk: tool confusion (`mypy` duplicate module names), harder static analysis.
- Recommended next move: define a single canonical import namespace policy and migrate incrementally behind compatibility wrappers.

3. Inconsistent quality scope vs repo breadth
- Full-repo lint/type is noisy due legacy/vendor-style areas.
- Risk: teams ignore checks because baseline appears unactionable.
- Recommended next move: keep scoped quality gate for active code and expand scope by directory as debt is reduced.

4. Missing top-level project orientation document
- A root onboarding guide is still missing.
- Risk: new contributors lack fast setup/run/debug path.
- Status: intentionally not added in this pass per user request.

## Suggested next slice (no behavior change)

1. Split one additional orchestration seam from `main.py` into `src/cli/`:
- candidate: distributed phase dispatch + queue-service bootstrap glue.

2. Add a compatibility map doc:
- one-page “canonical module paths” list for common imports (`path_router`, `post_docking`, `analysis`).

3. Expand quality gate coverage:
- include one additional stable directory at a time (for example `src/docking` subset), then ratchet.

## Benchmark verification plan linkage

The bench2 simulation + capture/compare flow is used as the non-functional regression guardrail:
- run baseline simulation (`-bench2 -fast`) and capture;
- apply professionalism-only changes;
- run identical simulation and compare with `tools/benchmark_bench2_compare.py`.
