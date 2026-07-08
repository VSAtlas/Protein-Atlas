# tools/bench2_aa_variance.py

## Purpose

Provide a repeatable A/A noise guardrail for bench2 quality deltas so ROC_AUC/PR_AUC drift can be classified as likely random variance vs likely regression.

## Inputs

- `analysis/benchmarks/<run_id>/comparison.json` produced by `tools/benchmark_bench2_compare.py compare`.
- A/A pairs should be same-day runs with identical code and invocation.

## Commands

Build variance profile from A/A comparisons:

```bash
python tools/bench2_aa_variance.py fit \
  --aa-comparison bench2_refactor_cycle2_candidate_20260304 \
  --aa-comparison bench2_refactor_cycle3_candidate_20260304 \
  --out analysis/benchmarks/bench2_aa_profile_20260304.json
```

Assess a new candidate comparison against that profile:

```bash
python tools/bench2_aa_variance.py assess \
  --profile analysis/benchmarks/bench2_aa_profile_20260304.json \
  --candidate-comparison bench2_refactor_cycle4_candidate_20260304 \
  --out analysis/benchmarks/bench2_aa_eval_cycle4_20260304.json
```

## Decision Protocol

1. Run at least 2-3 A/A pairs first and fit a profile.
2. Treat `likely_noise_only=true` as noise-compatible drift.
3. If flagged, rerun one more paired same-day baseline/candidate.
4. Treat repeated flags on the same key as likely real regression and block merge until triaged.
