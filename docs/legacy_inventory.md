## Legacy Inventory

This inventory classifies non-canonical scripts before deletion or refactor. It is intentionally conservative: quarantined means "do not add new code here", not "safe to remove today".

Status legend:

- `maintained`: actively maintained and covered by quality gates/tests.
- `compatibility-only`: retained as an import or command bridge; new logic belongs elsewhere.
- `quarantined`: retained for compatibility or historical reasons; not canonical for new work.
- `candidate for deletion`: remove only after the deletion rule below is satisfied.

Deletion rule:

Do not delete a legacy file until all are true:

1. no tests import it,
2. no docs recommend it,
3. no CLI path reaches it,
4. replacement exists or feature is explicitly unsupported.

## Initial inventory

| Path | Status | Keep? | Reason | Replacement / canonical owner |
| --- | --- | --- | --- | --- |
| `main.py` | compatibility-only | yes | developer-compatible root entrypoint for `atlas` behavior | `src/cli/atlas_main_cli.py` and `src/cli/` runtime modules |
| `src/cli/main_compat.py` | compatibility-only | yes | compatibility import bridge for root entry behavior | `src/cli/atlas_main_cli.py` |
| `src/cli/atlas_main_cli.py` | maintained | yes | packaged `atlas` console entrypoint | `src/cli/` |
| `tools/finalize_distributed_run.py` | maintained | yes | operational utility in current quality gate | `src/cli/` distributed lifecycle modules |
| `tools/simulate_slurm_array.py` | maintained | yes | distributed-mode simulation utility in current quality gate | `src/cli/` distributed/chunk modules |
| `tools/benchmark_bench2_compare.py` | maintained | yes | benchmark comparison utility in current quality gate | `analysis/dud_eval_core/` and `chemdb/bench/util_bench/` |
| `tools/benchmark_mode.py` | quarantined | maybe | legacy benchmark entry flow | `src/cli/` benchmark profile/runtime modules |
| `tools/benchmark_auto_analysis.py` | quarantined | maybe | legacy post-run benchmark helper | `analysis/reporting/` + `analysis/dud_eval_core/` |
| `tools/identify_ligand.py` | quarantined | maybe | one-off utility path | `chemdb/` (future `chemdb/tools/`) |
| `tools/data_analysis.py` | quarantined | maybe | older analysis helper outside maintained reporting owners | `analysis/reporting/` |
| `tools/pose_bust.py` | quarantined | maybe | external pose-checking utility not part of the current public CLI | `src/post_docking/` or `analysis/reporting/` after tests |
| `src/protein_prep/propka_wire.py` | quarantined | maybe | integration shim / legacy coupling | `src/protein_prep/` canonical prep interfaces |
| `src/docking/capture_pose.py` | quarantined | maybe | older pose-capture implementation excluded from current gates | `src/docking/capture_pose_*` maintained modules |
| `src/prep_ligands/dedup_sdfs.py` | quarantined | maybe | one-off ligand-library cleanup helper | `src/prep_ligands/library_index.py` or `src/prep_ligands/` owner module |
| `tools/molprobity_coot.py` | quarantined | maybe | external-tool wrapper utility | `src/protein_prep/` tool runner interfaces |

## Notes

- Quarantined does not mean unused; it means non-canonical for new implementation work.
- Promotion from quarantined to maintained requires explicit tests and gate inclusion.
- Candidate deletion requires checking imports, docs, CLI reachability, and replacement/unsupported status in the same patch or handoff.
- Runtime output directories are not legacy code. They should be ignored or documented as fixtures, not added to this inventory as source modules.
