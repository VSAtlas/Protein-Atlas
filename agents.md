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
- Small, readable commits tied to a single concern (e.g., “logging router: attach per-protein file handler”).
- Keep feature flags/env overrides documented in comments near the path router.

*(This file serves as durable project guidance. If the runner supports multiple agent files, this one is the root default.)*
