# Agent Profile: Atlas2 (BRCF)

## Identity & Scope
You are a coding and refactor assistant working inside the Atlas2 repo.
Primary focus: path-portable, variant-aware docking pipeline (APO vs HOLO, pH ensembles), robust logging, and minimal, safe patches.

## Environment Constraints
- Runtime host: UT Austin BRCF HPC (AMD GPU nodes; no internet).
- Python env: `docking-env` (micromamba/conda).
- Assume Linux shell; avoid OS-specific features.
- Do not run network calls, package installs, or system-level changes.

## Repository Layout (canonical paths)
- Input PDBs: `input_pdbs/`
- Processed: `processed_pdbs/<PDB>/{APO,HOLO}/receptor/`
- Prepped ligands: `prepped_ligands/<PDB>/`
- Docked outputs: `docked/<PDB>/{APO,HOLO}/`
- Config: `config.txt` (OVERALL_DIR, feature flags)
- Path router lives in code; do **not** re-invent path joins outside the router.

## Path Rules
- Centralize project paths via the existing path router (or helpers). 
- No ad-hoc `os.path.join` scattered in modules; consume typed returns from the router.
- Respect variant tokens (APO/HOLO/apo_vs_holo) and pH tokens in all outputs.

## Safety Rails
- Never delete or rewrite outside project directories.
- Never run `rm -rf` or destructive ops on: `input_pdbs/`, `raw/`, or global envs.
- Keep secrets out of code/logs. Assume logs are reviewable.

## Element/ion lists
- Use aliases.yaml for any lists of elements or  ions
- Do not remove past lists or edit any past lists in aliases.yaml
- If you find a list that would be appropriate to move to aliases.yaml move the list there and then update activesite.py so that it can read the updated list of aliases

## Code Patch Style
- Provide full replacement snippets (no diff markers like `+`/`-`).
- Include a brief header in replies:
  1) **Sharper coding prompt (reformat)**
  2) **Assumed folder structure**
  3) **Next action**
- Add concise in-code comments marking anchors and rationale.
- Prefer minimal, centralized edits over wide churn.

## Logging Expectations
- By default, show DEBUG?CRITICAL on console and in run logs.
- Always create a per-protein log file, even on early failure.
- Avoid duplicate handlers and double-emission; fix at a central setup.

## Active-Site & Variants
- If center was extracted from a ligand, **do not** overwrite with P2Rank.
- Respect APO vs HOLO semantics; if APO and HOLO are byte-identical post-prep, treat as ligand-free (APO) and deduplicate according to current policy.

## pH Ensemble
- Global mode default: allow whole-protein PROPKA-guided renaming; no center required.
- If radius is specified without a center, emit a clear error.

## Allowed Commands (examples)
- `python main.py --pdbs <PDB> --stages prep,dock`
- `python ph_ensemble.py --pdb <PDB> --ph-list 6.0 8.0 --global`
- `python analysis/dud_eval.py --docked-root docked/<PDB>`

## Forbidden Actions
- Installing packages, accessing the internet, touching user home outside the repo.
- Writing to arbitrary system paths or cluster-level config.
- Editing security policies or scheduler configs.

## Review & Tests
- Provide micro-tests or quick verification commands where possible.
- Log key decisions (e.g., center selection, dedup triggers, pocket IDs).

## Commit Hygiene
- Small, readable commits tied to a single concern (e.g.,  logging router: attach per-protein file handler ).
- Keep feature flags/env overrides documented in comments near the path router.

*(This file serves as durable project guidance. If the runner supports multiple agent files, this one is the root default.)*

Instrumentation & Debug Logging Policy

Goal: Ensure every critical decision and file-side effect is auditable. Add concise log lines at previously silent branches. Use only these levels: DEBUG, INFO, WARNING, ERROR.

Log format (uniform)
[level] [component] key1=val1 key2=val2 ... | msg


Required keys when available: run_id, pdb_id, variant, stage, path= (for file ops), n= (counts), time_ms= (timings).

Level semantics

DEBUG – fine-grained, frequent, cheap to emit; inputs/outputs of key functions, loop counters, branch choices.

INFO – milestones and summaries; start/finish of a stage, output paths, final counts.

WARNING – unexpected but recoverable conditions; fallbacks, partial data, skipped items.

ERROR – unrecoverable for the current unit (pdb/ligand/stage); include why and what was skipped. Do not crash the whole run unless upstream requires.

Critical checkpoints (add logs even if they don’t exist yet)

Path router & config resolution

INFO [path-router] run_id=… pdb_id=… variant=… OVERALL_DIR=…

DEBUG [path-router] resolved receptor_pdbqt=… prepped_dir=… docked_dir=…

Input PDB ingest & receptor prep

INFO [receptor-prep] pdb_id=… stage=load_pdb source=…

WARNING [receptor-prep] missing_altloc_handling=… action=defaulted

ERROR [receptor-prep] reduce_failed rc=… path=…

pH / PROPKA / microstate selection

INFO [microstates] ligand=… pH=7.4 forms_considered=… chosen=…

DEBUG [microstates] rule=major_microstate basis=pKa_tool=propka

Ligand prep & filtering

INFO [lig-prep] pdb_id=… n_in=… n_ok=… n_skipped=… out_dir=…

WARNING [lig-prep] ligand=… reason=malformed wrote=malformed_ligands.txt

Pocketing / center & box definition

INFO [pocket] pdb_id=… method=p2rank n_pockets=… chosen=… center=… radius=…

WARNING [pocket] fallback=ligand_center reason=p2rank_empty

Control redock (crystal vs docked)

INFO [control-centers] pdb_id=… n=… max_spread_A=… policy=…

DEBUG [control-redock] lig=… rmsd_A=… score_kcal=… chosen=(0|1)

Docking invocation & pose validation

INFO [docking] engine=vina pdb_id=… n_ligands=… exhaustiveness=… n_modes=…

DEBUG [pose-validate] ligand=… pose_i=… status=valid|too_far|clash dist_to_center_A=…

Scoring / rescoring / consensus

INFO [rescore] method=gnina_cnn pdb_id=… n_poses=…

DEBUG [consensus] ligand=… vina=… cnn=… combined=…

Evaluation & exports

INFO [eval] run_id=… out=analysis/out/summary_<runid>.tsv targets=…

DEBUG [eval] target=… ROC_AUC=… PR_AUC=… EF@1%=…

Failures & early exits

ERROR [stage] pdb_id=… reason=… action=skipped

Places to retrofit logs now (previously sparse)

Dedup of APO/HOLO or pH states (which kept, which collapsed, why).

Run-ID filtering in dud_eval.py (rows kept vs total).

Reduce/element-fixing branch choice (skip/rerun/nohyd).

Fallback pocket recentering loop (which ligand re-centered, success/fail).

Malformed ligand pipeline (first failure cause per ligand).

GNINA CNN rescoring enablement (on/off, model, thresholds).

Consensus ablation winner (delta metrics per target).

CSV/TSV writes (exact output paths, row counts).

Message style guide

Be structured: prefer key=value pairs over prose.

Be stable: avoid dynamic wording; keep tags grep-friendly (e.g., [control-centers], [rescore]).

No secrets: never log credentials or full file contents.

One fact per field: don’t cram lists; emit multiple DEBUG lines if needed.

Minimal examples (copy/paste patterns)
INFO  [path-router] run_id=v2025_11_07 pdb_id=1T46 variant=HOLO OVERALL_DIR=/…/atlas2
DEBUG [path-router] receptor_pdbqt=/…/processed_pdbs/1T46/HOLO/receptor/receptor.pdbqt

INFO  [control-centers] pdb_id=1T46 n=2 max_spread_A=0.84 policy=best_redock
DEBUG [control-redock] lig=ATP rmsd_A=1.72 score_kcal=-10.1 chosen=1
DEBUG [control-redock] lig=ADP rmsd_A=2.63 score_kcal=-9.2 chosen=0

INFO  [docking] engine=vina pdb_id=1T46 n_ligands=10861 exhaustiveness=8 n_modes=9
DEBUG [pose-validate] ligand=fda_1158 pose_i=0 status=valid dist_to_center_A=3.1

INFO  [eval] run_id=v2025_11_07 out=analysis/out/summary_v2025_11_07.tsv targets=24
DEBUG [eval] target=1T46 ROC_AUC=0.760 PR_AUC=0.073 EF@1%=6.33
WARNING [eval] run_id requested but column missing; proceeding unfiltered
ERROR [receptor-prep] pdb_id=3VO3 reason=reduce_failed rc=1 action=skipped

Reviewer-friendly checkboxes (per PR)

 New code paths emit at least one INFO and one DEBUG line.

 All fallbacks emit a WARNING with reason= and fallback=.

 Any skipped unit emits an ERROR with reason= and action=skipped.

 Outputs log the exact path and row counts.

 No duplicate handlers; logs appear once in console and in per-protein file.

Implementation hint for Codex: Add logs where decisions are made (conditionals, returns, exceptions) and where side effects happen (file writes, deletions, external tool calls). Keep messages ≤ 120 chars, keys in snake_case, and reuse the component tags above.


## Patch Boundaries & Protected Paths (Must Read)

**Goal:** Prevent patch collisions between automation and local/CI setup by strictly scoping what the agent may edit.

### Allowed default patch scope
Focus patches in:
- `code/protein_automation/**` (Python pipeline code)
- `analysis/**` (non-test analysis helpers)
- `docs/**` (text only)

### Protected paths (DO NOT EDIT unless explicitly requested)
- `tests/**`
- `tests/data/**`  ? cached artifacts only; never commit binaries here
- `.github/workflows/**`
- `ci/**`
- `input_pdbs/**`  ? source inputs; never auto-rewrite
- `processed_pdbs/**`, `docked/**`  ? build outputs; hands off
- `agents.md`  ? only edit when asked to change policy

**If your patch touches any protected path, abort** unless the PR title or top of the user�s prompt contains one of:
- `override:tests`, `override:ci`, or `override:agents`

### Binary / data policy
- Never vendor PDBs, SDFs, MOL2, or other binaries into the repo.
- Tests must **download/cache** required artifacts at runtime under `tests/data/` and rely on `.gitignore`.
- If a required file is missing, update tests to fetch or to skip with a clear message; do **not** add the file to git.

### Environment & CI policy
- Use **only** `ci/run_in_env.sh` to create/activate the environment and run tests.
- Do not modify `.github/workflows/**` unless the user asked for CI changes (`override:ci`).
- Do not add new package installs in code; place them in `environment.yml` and rely on the CI/bootstrap script.

### Test policy
- Do not edit `tests/**` unless the user asked for test changes (`override:tests`).
- Acceptance tests (e.g., ions retention) are the source of truth; code should be modified to satisfy them.

### Conflict-avoidance with setup scripts
- If a patch would change files that the setup/bootstrap scripts also touch in the same run, **abort** and emit:
  `ERROR [agent-guard] protected_paths_touched=<list>; action=abort`
- Prefer adding **logging** or **configuration flags** instead of restructuring shared bootstrapping code.

### PR checklist (agent MUST enforce)
- [ ] No changes in protected paths (unless override tag present).
- [ ] No binaries added to the repo.
- [ ] Patch is single-topic and minimal.
- [ ] New/changed code emits clear logs for key decisions (grep-able tags).
- [ ] Tests pass locally via `ci/run_in_env.sh`.

### Examples
**OK:** Add debug lines in `code/protein_automation/automate_protein_prep.py` and `main.py`.
**NOT OK:** Edit `tests/test_ions_acceptance.py` to �make it pass� (unless PR is tagged `override:tests`).

