# PR Integration Plan

Date: 2026-05-02

Source: GitHub open PR snapshot for `VSAtlas/protein-stuff`, plus changed-file lists from the GitHub connector. This report is an integration queue plan only; it does not approve merges and does not validate runtime behavior.

## Executive Summary

There are 67 open PRs in the queue, mostly stacked agent-generated branches. Many are based on generated intermediate branches such as `brcf-update-FAST_EDITION-v0.*` or older path-router rollout branches instead of the current `main`. Treat every branch as stale until it is rebased or resliced onto `main`.

The queue should not be merged as a chronological stack. It has overlapping ownership of `main.py`, `prep_ligands.py`, `automate_protein_prep.py`, `path_router.py`, `run_vina.py`, and newly extracted helper modules. Several PRs are duplicates or superseded by later PRs. PR #54 is especially unsafe as an integration unit because its file list includes generated caches, IDE files, ligand artifacts, input files, analysis outputs, and other non-code material.

Recommended policy:

1. Close or supersede contaminated and duplicate branches first.
2. Recreate a small number of clean integration branches from current `main`.
3. Merge behavior-preserving extraction slices before scientific behavior changes.
4. Require the full Phase 11 acceptance suite for each candidate merge, plus focused tests listed below.

## Active PRs By Ownership Area

| Area | Active PRs / branches | Initial disposition |
| --- | --- | --- |
| `main.py` pathing and orchestration | #1 `codex/refactor-path-handling-in-main.py`, #13, #14, #15, #22, #24, #25, #26, #27, #28, #29, #30, #37, #38, #42, #43, #45, #46, #47, #55, #58, #59, #60, #61, #62, #63, #64, #65, #67 | Heavy overlap; reslice by topic before merge. |
| Docking helpers | #4 `codex/update-run_vina.py-for-path-management`, #10, #23, #24, #27, #28, #61, #64, #65, #66 | Merge only after path-router baseline is settled. |
| pH handling and ligand microstates | #5, #22, #23, #24, #29, #30, #49, #51, #52, #53, #54, #55, #56, #57, #61, #62, #63 | Risky scientific behavior area; split refactors from behavior changes. |
| APO/HOLO handling and receptor prep | #1, #9, #30, #36, #37, #38, #39, #40, #42, #43, #44, #45, #46, #47, #48, #60 | Highest scientific-risk cluster; requires HOLO/APO focused validation. |
| Logging | #15, #18, #20, #21, #28, #35, #57, #58, #59, #62, #67 | Mostly low scientific risk, but root logger changes affect all tests. |
| Single-ligand helpers | #25, #64, #65, #67 | Keep #25 behavior separate from #65/#67 extraction/refactor work. |

## PR Classification

Legend: `refactor`, `bug fix`, `risky behavior`, `incomplete`, `superseded`.

| PR | Branch | Files observed | Classification | Merge-readiness note |
| --- | --- | --- | --- | --- |
| #67 | `codex/move-debug-wrappers-to-debug_fs.py` | `debug_fs.py`, `main.py` | refactor, incomplete | Failed smoke with missing `yaml`; validate after #64/#65 because it imports extracted helpers. |
| #66 | `codex/fix-near-miss-recentering-in-docking.py` | `docking.py` | bug fix, incomplete | Good candidate after docking helper module lands; failed smoke due environment. |
| #65 | `codex/fix-importerror-in-docking.py` | `docking.py`, `main.py`, `single_ligand_index.py` | refactor, incomplete | Overlaps #64 and #67; merge only after deciding final single-ligand helper ownership. |
| #64 | `codex/refactor-docking-helpers-into-docking.py` | `checkpoints.py`, `docking.py`, `main.py`, `path_router.py`, `record_data.py`, `single_ligand_index.py` | refactor, incomplete | Broad extraction branch; likely prerequisite for #65-#67 but must be resliced. |
| #63 | `codex/remove-internal-config-keys-for-ph-handling` | `main.py` | refactor, risky behavior, incomplete | Claims behavior preservation but changes pH state propagation; require pH ensemble tests. |
| #62 | `codex/add-ph-aware-debug-logging-to-ligand-pipeline` | `main.py`, `prep_ligands.py` | logging, incomplete | Low risk if log-only; verify no ligand-selection side effects. |
| #61 | `codex/refactor-ph-ensemble-docking-logic-to-module` | `fallback_recenter.py`, `main.py`, `multi_stage_docking.py`, `ph_ensemble_docking.py`, `prep_ligands.py` | refactor, risky behavior, incomplete | Extracts pH orchestration; overlaps #22-#24 and #63. |
| #60 | `codex/refactor-apo/holo-handling-into-new-module` | `apo_holo_mode.py`, `fallback_recenter.py`, `main.py`, `prep_ligands.py` | refactor, risky behavior, incomplete | Valuable extraction, but APO/HOLO semantics need focused tests. |
| #59 | `codex/fix-logging-issue-in-bootstrap_root_logging` | `logging_topics.py`, `main.py`, `prep_ligands.py` | bug fix, incomplete | Low scientific risk; root logger changes require full smoke tests. |
| #58 | `codex/refactor-logging-filter-logic-to-logging_topics.py` | `logging_topics.py`, `main.py` | refactor, logging, incomplete | Merge before #59 if keeping both; otherwise superseded by #59. |
| #57 | `codex/implement-manifest-fallback-and-logging-in-prep_ligands.py` | `prep_ligands.py` | bug fix, risky behavior, incomplete | Fallback behavior can change ligand coverage; require ligand microstate tests. |
| #56 | `codex/update-ph-aware-microstate-prep-logic` | `prep_ligands.py` | bug fix, risky behavior, incomplete | Candidate after #49/#52 direction is chosen. |
| #55 | `codex/honor-root_dir-in-docking-pipeline` | `library_index.py`, `main.py`, `prep_ligands.py` | bug fix, risky behavior, incomplete | Overlaps #56 and #63; verify selected ligand root exactly. |
| #54 | `codex/add-enumerate_ligands_for_docking-function` | thousands of source, data, cache, IDE, CI, ligand, input, and output files | incomplete, superseded/blocker | Close or recreate from scratch. Do not merge this branch. |
| #53 | `codex/add-alias-index-for-microstate-reuse` | not re-fetched in full; summary indicates microstate registry code | bug fix, risky behavior, incomplete | Potentially useful but should be resliced after #56/#57. |
| #52 | `codex/update-pdbqt-filename-for-ph-subdirs` | `prep_ligands.py` | bug fix, risky behavior, incomplete | Supersedes #51. Requires backward-compatible filename tests. |
| #51 | `codex/modify-pdbqt-filenames-for-ph-specificity` | `prep_ligands.py` | superseded, incomplete | Close in favor of #52 unless #52 loses behavior. |
| #50 | `codex/update-pdbqt-naming-based-on-sdf-parent-names` | `prep_ligands.py` | bug fix, incomplete | Independent naming fix; validate with ligand prep fixtures before pH filename work. |
| #49 | `codex/extend-bulk-ligand-prep-to-support-ph-list` | `prep_ligands.py` | risky behavior, incomplete | Core multi-pH ligand prep behavior; requires explicit approval and tests. |
| #48 | `codex/add-coordination-summary-to-metal-audit-output` | `automate_protein_prep.py` | bug fix, incomplete | Output-shape change to audit JSON; merge after #47 if kept. |
| #47 | `codex/add-metal-site-audit-comparison` | `automate_protein_prep.py`, `main.py` | bug fix, risky behavior, incomplete | Audit-only intent, but touches prep path; verify artifact paths. |
| #46 | `codex/handle-holo-regen-rebuild-inline` | `main.py` | bug fix, risky behavior, incomplete | Supersedes #45; changes HOLO receptor regeneration timing. |
| #45 | `codex/add-receptor-regeneration-after-holo-restore` | `main.py` | superseded, risky behavior, incomplete | Close in favor of #46 if #46 is correct. |
| #44 | `codex/gate-metal-coordinating-ligand-restoration-by-docking-box` | `automate_protein_prep.py` | risky behavior, incomplete | Scientific behavior change; validate metal/cofactor retention. |
| #43 | `codex/plumb-docking-context-into-holo-restore` | `automate_protein_prep.py`, `main.py` | refactor, risky behavior, incomplete | Likely prerequisite for #44; can be resliced as plumbing only. |
| #42 | `codex/respect-docking-box-when-restoring-holo-ligands` | `automate_protein_prep.py`, `main.py` | risky behavior, incomplete | Likely superseded by #43/#44. |
| #41 | `codex/handle-variant-and-ph-metadata-in-dud-evaluator` | `analysis/dud_eval.py` | bug fix, incomplete | Analysis-output behavior; outside core docking but relevant to APO/HOLO/pH reporting. |
| #40 | `codex/add-metal-coordination-audit-json` | `automate_protein_prep.py` | incomplete, superseded | Likely superseded by #47/#48. |
| #39 | `codex/guard-pdb2pqr-hetero-retention-and-pass-variant-context` | `activesite.py`, `automate_protein_prep.py`, `chemdb/aliases.yaml`, `propka_wire.py` | bug fix, risky behavior, incomplete | Strongest candidate in #36-#39 chain if tests pass. |
| #38 | `codex/ensure-variant-scoped-prep-keeps-metals` | `automate_protein_prep.py`, `main.py`, `propka_wire.py` | superseded, risky behavior, incomplete | Close if #39 covers it. |
| #37 | `codex/preserve-metals-in-holo-prep-and-scope-variant-outputs` | `automate_protein_prep.py`, `main.py`, `propka_wire.py` | superseded, risky behavior, incomplete | Close if #39 covers it. |
| #36 | `codex/safeguard-metals-through-pdb2pqr-and-reduce-guard` | `automate_protein_prep.py`, `propka_wire.py` | bug fix, risky behavior, incomplete | May be base concept for #39; do not merge independently without comparison. |
| #35 | `codex/instrument-clean_pdb-ion-tracing-checkpoints` | summary indicates `automate_protein_prep.py` | logging, incomplete | Keep only if needed after #39 audit/metal work. |
| #30 | `codex/ensure-ph-ensemble-artifacts-honor-apo-holo-variants` | `activesite.py`, `automate_protein_prep.py`, `main.py`, `path_router.py`, `ph_ensemble.py`, `tests/test_ions_acceptance.py` | risky behavior, incomplete | Cross-cuts pH and APO/HOLO; merge only after both baselines settle. |
| #29 | `codex/fix-emit_vina_config-invocation-and-add-config-logging` | `library_index.py`, `main.py`, `ph_ensemble.py` | bug fix, logging, incomplete | Candidate after path-router pH config API is final. |
| #28 | `codex/add-debug-logging-for-apo-holo-mode-normalization` | `input_and_export_functions.py`, `main.py`, `path_router.py`, `run_vina.py` | logging, incomplete | Low risk if log-only; overlaps path-router files. |
| #27 | `codex/align-config-routing-with-path-helpers-for-variant-and-ph-modes` | `input_and_export_functions.py`, `main.py`, `path_router.py`, `run_vina.py` | refactor, risky behavior, incomplete | Path-router routing change; validate with path tests before pH/HOLO branches. |
| #26 | `codex/add-automatic-run-selection-and-include-run-ids-in-summaries` | `analysis/dud_eval.py`, `main.py`, `path_router.py` | bug fix, incomplete | Analysis/reporting queue; keep separate from docking merges. |
| #25 | `codex/restore-single-ligand-resolution` | `main.py` | bug fix, incomplete | Good standalone candidate before #65 if single-ligand tests pass. |
| #24 | `codex/execute-multi-stage-docking-for-ph-ensemble-members` | `config.txt`, `input_and_export_functions.py`, `main.py`, `path_router.py`, `run_vina.py` | risky behavior, incomplete | Do not merge config mutation from PR; reslice code only. |
| #23 | `codex/support-ph-aware-docking-paths` | `path_router.py` | refactor, risky behavior, incomplete | Likely prerequisite for #24; verify legacy layout. |
| #22 | `codex/add-PH_ENSEMBLE-pipeline-integration` | `main.py`, `ph_ensemble.py` | risky behavior, incomplete | Feature integration branch; depends on #23/#24 decisions. |
| #21 | `codex/enhance-control-log-instrumentation` | `analysis/dud_eval.py` | logging, incomplete | Merge only in analysis queue. |
| #20 | `codex/add-optional-control-redock-report` | `analysis/dud_eval.py` | bug fix/feature, incomplete | Analysis queue; verify output contracts. |
| #18 | `codex/add-debug-logging-to-dud_eval.py` | `analysis/dud_eval.py` | logging, incomplete | Superseded by later DUD evaluator modularization if already on main. |
| #15 | `codex/add-run_id-to-docking-score-csvs` | `input_and_export_functions.py`, `main.py` | bug fix, output interface change, incomplete | Requires CSV schema/output compatibility review. |
| #14 | `codex/update-ligand-library-selection-in-test-mode` | `main.py` | bug fix, incomplete | Candidate if `test_ligand_run_modes.py` covers it. |
| #13 | `codex/fix-ligand-extraction-argument-errors` | `main.py` | bug fix, incomplete | Failed due missing `numpy`; retest in correct environment. |
| #11 | `codex/refactor-path-handling-to-centralized-router-o45rsu` | `input_and_export_functions.py` | refactor, incomplete | Likely superseded by current path-router baseline. |
| #10 | `codex/refactor-path-handling-to-centralized-router` | `fallback_recenter.py` | refactor, incomplete | Likely superseded by later docking/fallback branches. |
| #9 | `codex/update-path-handling-in-activesite.py` | `activesite.py` | refactor, incomplete | Old path-router rollout; compare against current main before keeping. |
| #5 | `codex/add-router-aware-convenience-for-pdb-selection` | `context_ph.py` | refactor, incomplete | Old pH path-router helper; likely superseded by newer pH work. |
| #4 | `codex/update-run_vina.py-for-path-management` | `run_vina.py` | refactor, incomplete | Old path-router rollout; likely superseded by #23/#27. |
| #3 | `codex/update-prep_ligands.py-with-paths-import-and-init` | `prep_ligands.py` | refactor, incomplete | Old path-router rollout; likely superseded by #49-#57. |
| #2 | `codex/update-input_and_export_functions.py` | `input_and_export_functions.py` | refactor, incomplete | Old path-router rollout; likely superseded by current main or #11/#27. |
| #1 | `codex/refactor-path-handling-in-main.py` | `automate_protein_prep.py`, `main.py` | refactor, draft, incomplete | Do not merge draft; use only as historical context. |

## Overlap And Conflict Risks

`main.py` is the central collision point. PRs #13-#15, #22, #24-#30, #37/#38, #42/#43/#45/#46/#47, and #55/#58-#65/#67 all modify it. A direct queue merge will produce repeated conflict resolution and a high risk of silently mixing pH, APO/HOLO, logging, and single-ligand assumptions.

`prep_ligands.py` has conflicting ownership between pH/microstate work (#49-#57, #62), logging work (#59), APO/HOLO extraction (#60), and docking helper extraction (#61). Merge pH behavior only after root selection and filename conventions are chosen.

`automate_protein_prep.py` has overlapping HOLO/metal/cofactor changes (#36-#48 plus #30). These are not simple refactors: they affect retained receptor chemistry and should be treated as scientific behavior changes.

`path_router.py` and `run_vina.py` are modified by older path-router rollout PRs and newer pH/variant branches (#23/#24/#27/#28/#30/#64). Establish the current router contract first, then reslice dependent changes.

The PR bases indicate a generated branch stack, not a clean queue against `main`. Several heads use older bases; mergeability against their base does not imply merge-readiness against `main`.

## Recommended Merge Order

1. Queue cleanup, no code merge:
   - Close #54 or replace with a clean branch.
   - Close clear duplicates/superseded PRs: #1, #2, #3, #4, #5, #10, #11, #18, #37, #38, #40, #42, #45, #51.
   - Mark all remaining PRs as needing rebase/reslice onto current `main`.

2. Low-risk standalone bug fixes:
   - #13 ligand extraction canonical directory arguments.
   - #14 test-mode ligand root union.
   - #25 single-ligand resolution.
   - #66 near-miss recentering centroid lookup, after `docking.py` exists or is resliced without requiring #64.

3. Logging-only queue:
   - #58 logging topic extraction.
   - #59 root handler bootstrap fix.
   - #62 pH ligand debug logging only if diff proves no selection behavior change.
   - #28 APO/HOLO routing debug logging only if path-router changes are already present or removed.

4. Behavior-preserving extraction queue:
   - #64 docking helper extraction, resliced to avoid unrelated `path_router.py` and `record_data.py` changes if possible.
   - #65 single-ligand helper extraction.
   - #67 debug filesystem helper extraction.
   - #60 APO/HOLO helper extraction only if it is strictly mechanical; otherwise move to scientific queue.
   - #61 pH ensemble helper extraction only if it is strictly mechanical; otherwise move to scientific queue.

5. Path-router and output contract queue:
   - #23 pH-aware docking paths.
   - #27 config routing with path helpers.
   - #15 run_id in logs and CSV outputs.
   - #26/41 DUD evaluator metadata branches, separately from runtime docking.

6. pH scientific behavior queue:
   - #49 multi-pH ligand prep.
   - #50 parent ligand naming.
   - #52 pH-specific PDBQT filenames.
   - #53 duplicate microstate alias prevention.
   - #56 root_dir-aware microstate prep.
   - #57 manifest fallback.
   - #55 custom ligand roots.
   - #63 explicit pH context instead of hidden config keys.
   - #22/#24 pH ensemble pipeline and multi-stage docking only after the above are stable.

7. APO/HOLO and receptor chemistry queue:
   - #36 or #39, not both unless manually reconciled.
   - #43 docking context plumbing.
   - #44 coordination-aware HOLO restoration policy.
   - #46 inline HOLO regen rebuild.
   - #47/#48 metal audit output additions.
   - #30 cross-product pH ensemble plus APO/HOLO variants last.

## Required Tests Per PR Class

All candidate merges:

```bash
tools/quality_gate.sh
pytest chemdb/tests/test_main_full_run.py -q
pytest chemdb/tests/test_ligand_run_modes.py -q
pytest chemdb/tests/test_path_router.py -q
pytest chemdb/tests/test_record_data_csv.py -q
```

Single-ligand and ligand-root changes (#13/#14/#25/#50/#55/#65):

```bash
pytest chemdb/tests/test_ligand_run_modes.py -q
pytest chemdb/tests/test_library_mode.py -q
pytest chemdb/tests/test_record_data_csv.py -q
```

pH and microstate changes (#22/#23/#24/#29/#30/#49/#52-#57/#61-#63):

```bash
pytest chemdb/tests/test_path_router.py -q
pytest chemdb/tests/test_control_redock_after_ph_ensemble.py -q
pytest chemdb/tests/test_control_redock_per_ph.py -q
pytest chemdb/tests/test_ph_ensemble_withh_resolution.py -q
pytest chemdb/tests/test_ligand_microstates_ph_debug.py -q
```

APO/HOLO and receptor chemistry changes (#30/#36/#39/#43/#44/#46-#48/#60):

```bash
pytest chemdb/tests/test_holo_restore.py -q
pytest chemdb/tests/test_ion_audit.py -q
pytest chemdb/tests/test_metal_site_audit.py -q
pytest chemdb/tests/test_receptor_prep_unit.py -q
pytest chemdb/tests/test_protein_prep_refactor_receptor_protonation.py -q
```

Logging-only changes (#18/#21/#28/#35/#58/#59/#62/#67):

```bash
pytest chemdb/tests/test_summarize_run_errors_hermetic.py -q
pytest chemdb/tests/test_main_full_run.py -q
```

Analysis/DUD evaluator changes (#18/#20/#21/#26/#41):

```bash
pytest chemdb/tests/test_dud_eval_modularization_smoke.py -q
pytest chemdb/tests/test_dud_eval_runid_layout.py -q
pytest chemdb/tests/test_dud_eval_microstates.py -q
pytest chemdb/tests/test_control_redock_after_ph_ensemble.py -q
```

## Branches To Close Or Supersede

Close/recreate immediately:

- #54 `codex/add-enumerate_ligands_for_docking-function`: contaminated by generated artifacts, caches, IDE files, input/output data, and broad unrelated changes.

Close if current main already contains equivalent path-router work:

- #1, #2, #3, #4, #5, #9, #10, #11.

Close as superseded after comparison:

- #18 in favor of later DUD evaluator logging/modularization work.
- #37 and #38 in favor of #39.
- #40 in favor of #47/#48.
- #42 in favor of #43/#44.
- #45 in favor of #46.
- #51 in favor of #52.
- #58 if #59 includes the intended logging topic extraction and root-handler fix.

## Blockers

- The queue has no reliable `main`-based merge order. Most PRs are based on generated intermediate branches.
- Several PR bodies report only `py_compile`, no tests, or failed runtime smoke tests due missing `yaml`, `Bio`, or `numpy`. That indicates environment drift and blocks merge readiness.
- #54 is not reviewable as a code PR because it includes large generated/runtime artifacts and unrelated files.
- pH and APO/HOLO changes are scientific behavior changes unless proven otherwise by focused tests and diff review.
- The Phase 11 acceptance suite has not been run against these PRs in this audit.

## Verification During This Audit

The documentation-only patch was checked against the current working tree, not against any open PR branch.

Passed:

```bash
/stor/home/mpg2352/micromamba/bin/micromamba run -n docking-env tools/quality_gate.sh
```

Blocked or failing in the current environment:

```bash
/stor/home/mpg2352/micromamba/bin/micromamba run -n docking-env pytest chemdb/tests/test_main_full_run.py -q
/stor/home/mpg2352/micromamba/bin/micromamba run -n docking-env pytest chemdb/tests/test_ligand_run_modes.py -q
/stor/home/mpg2352/micromamba/bin/micromamba run -n docking-env pytest chemdb/tests/test_path_router.py -q
/stor/home/mpg2352/micromamba/bin/micromamba run -n docking-env pytest chemdb/tests/test_record_data_csv.py -q
/stor/home/mpg2352/micromamba/bin/micromamba run -n docking-env python main.py --test -fast
```

Observed blockers:

- Historical blocker resolved: receptor prep no longer requires `PREPARE_RECEPTOR_SCRIPT` in runtime config validation.
- `test_path_router.py` has 25 passing tests and 2 failures where router roots remain relative (`.` and `docked/test_run`) while the tests expect resolved absolute paths.

## Next Actions

1. Close or label the superseded PRs above before opening more refactor slices.
2. Pick one owner for each domain: docking helper extraction, pH ligand prep, APO/HOLO receptor prep, logging, and single-ligand helpers.
3. Recreate clean branches from current `main` with one ownership boundary each.
4. Require a short PR body section: behavior class, changed files, test commands, and whether scientific behavior can change.
5. Add the environment doctor task before relying on agent smoke-test results, because multiple PRs failed for missing Python dependencies that should exist in the intended environment.
