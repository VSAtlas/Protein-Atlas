## Module Ownership

This ownership map defines where new maintained implementation code should live.

Policy:

- Keep root-level scripts as wrappers when possible.
- Do not add new maintained logic to compatibility wrappers.
- Route new functionality to canonical owner modules.

## Canonical ownership map

| Responsibility | Canonical owner |
| --- | --- |
| CLI parsing and run context | `src/cli/` |
| Path resolution | `src/path_router/` |
| Docking orchestration | `src/docking/` |
| Ligand preparation | `src/prep_ligands/` |
| Receptor/protein preparation | `src/protein_prep/` |
| Post-docking rescoring | `src/post_docking/` |
| MM/GBSA | `src/post_docking/mmgbsa/` |
| Reporting and heatmap | `analysis/reporting/` |
| DUD and benchmark evaluation | `analysis/dud_eval_core/` |
| Chemical aliases and mapping | `chemdb/` |

## Ownership details

| Owner | Owns | Does not own |
| --- | --- | --- |
| `src/cli/` | CLI argument normalization, run context, run lifecycle, manifest support, distributed/chunk planning, smoke/test-mode setup, post-run hook orchestration | Docking/scoring algorithms, reporting computations |
| `src/path_router/` | Portable path construction, run-scoped path routing, pH-aware path context, compatibility path shims | Hardcoded host paths, ad hoc path joins in feature modules |
| `src/docking/` | Docking orchestration, engine adapters, receptor/ligand docking runtime phases, control redocking, pH ensemble docking, pose capture, docking score I/O | CLI parsing, reporting tables, ligand preparation chemistry |
| `src/prep_ligands/` | Ligand conversion/preparation, library indexing, microstate handling, test/FDA library prep helpers | Docking engine execution, final report rendering |
| `src/protein_prep/` | Receptor cleaning, protonation, water/ion handling, altloc and chain pruning, external protein-prep tool runners | Docking run orchestration, path policy |
| `src/post_docking/` | Artifact retention, rescoring, reranking, SCORCH/CNN support, MM/GBSA workflow modules | Primary docking execution, heatmap/report UI |
| `analysis/reporting/` | Run reports, heatmap HTML, manifest-derived reporting, master schema export, target/ligand annotations used by reports | Docking/scoring runtime behavior |
| `analysis/dud_eval_core/` | DUD/decoy evaluation metrics, labels, schema guessing, evaluation orchestration | General CLI runtime and docking execution |
| `chemdb/` | Chemical aliases, mapping databases/utilities, target difficulty helpers, benchmark utility package | Runtime docking outputs or public CLI ownership |
| `src/ml/` | Legacy `ml.*` package implementation and BigBind/model helper utilities | Runtime ML datasets, model output directories, or Atlas ADR analysis ML under `atlas/analysis/ml/` |

## Wrapper guidance

Compatibility wrappers should include explicit comments:

- "Compatibility wrapper only."
- "No new logic here."
- "Use <canonical_module> for new code."

Initial wrappers to treat as compatibility boundaries:

- `main.py`
- `src/cli/atlas_main_cli.py`
- `src/cli/main_compat.py`

## Quality gate policy

`tools/quality_gate.sh` is the source of truth for the current maintained quality surface. It intentionally gates owner modules before legacy scripts:

- Ruff is scoped to maintained CLI/config/path/rescoring targets plus selected tools and `main.py`.
- MyPy is scoped to selected maintained runtime/config/tool modules.
- Quarantined scripts, generated/runtime directories, and protected test/runtime data stay excluded until intentionally adopted.

When expanding the gate, add one owner directory at a time and fix only import errors, clear dead code, formatting, and type annotations that improve clarity. Do not use quality-gate expansion as cover for scientific behavior changes.

## Root script policy

Root-level scripts and bridge modules are compatibility surfaces. They may parse or forward arguments only when that behavior is already established. New behavior belongs in the canonical owner modules above.

Every compatibility wrapper should make its status clear near the top of the file:

```text
Compatibility wrapper only.
No new logic here.
Use <canonical_module> for new code.
```
