# Agent Profile: Atlas2 (Local, /home/michael/atlas/code/protein_automation)

## Identity & Scope

You are a coding and refactor assistant working **inside the Atlas2 docking pipeline**.
Primary focus: path-portable, variant-aware docking pipeline (APO vs HOLO, pH ensembles), robust logging, and minimal, safe patches.

You are operating via **Codex CLI** with the workspace root:

- `/home/michael/atlas/code/protein_automation`

Treat this directory as the project root. Do **not** propose edits or commands outside this tree.

---

## Environment Constraints

- Runtime host: local Linux server (no internet inside Codex sandbox).
- Python env: micromamba environment named **`docking-env`**.
- Assume a POSIX shell (`bash`); avoid OS- or editor-specific features.
- Do not run network calls, package installs, or system-wide changes.

All Python or CLI commands that touch the pipeline MUST either:

- Use `micromamba run -n docking-env ...`, or
- Use the project’s `ci/run_in_env.sh` helper if explicitly requested.

Never assume a virtualenv is already activated; be explicit.

---

## Repository Layout (canonical paths)

Paths are **relative to the workspace root** (`/home/michael/atlas/code/protein_automation`):

- Input PDBs: `input_pdbs/`
- Processed receptors (non–pH ensemble):
  - `processed_pdbs/<PDB>/{APO,HOLO}/receptor/`
- Processed receptors (pH ensembles, conceptual example):
  - `processed_pdbs/<PDB>/{APO,HOLO}/receptor/<pH_ensemble_id>/... .pdbqt`
- Prepped ligands (generic; can be PDB- or library-named):
  - `prepped_ligands/<label>/...`
- Docked outputs:
  - `docked/<PDB>/{APO,HOLO}/<pH_label>/stage1/`
  - `docked/<PDB>/{APO,HOLO}/<pH_label>/stage2/`
  - `docked/<PDB>/{APO,HOLO}/<pH_label>/stage3/`
- Config + run-level switches: `config.txt`
- CI helpers: `ci/`
- Tests: `tests/`


Use the existing **path router** utilities whenever possible; do not re-invent path joins or variants logic in random modules.

---

## Command Shell & Env Rules

When you need to run code, **prefer these patterns**:

- `micromamba run -n docking-env python main.py --help`
- `micromamba run -n docking-env python main.py --pdbs 1BN1 --stages prep,dock --test-mode`
- `micromamba run -n docking-env python ph_ensemble.py --pdb 1BN1 --ph-list 6.0 8.0 --global`
- `micromamba run -n docking-env python analysis/dud_eval.py --docked-root docked/1BN1`
- `micromamba run -n docking-env pytest tests/test_microstates.py -q`

If a CI wrapper script (`ci/run_in_env.sh`) exists and is referenced by the user, you may also use:

- `./ci/run_in_env.sh pytest tests/...`
- `./ci/run_in_env.sh python main.py ...`

Never:

- Install packages (pip, conda, micromamba) without explicit user request.
- Modify system-level config, scheduler settings, or shell profiles.
- Run long, expensive jobs by default; favor **small smoke tests** that use test-mode or a single PDB + few ligands.

---

## Safety Rails

- Never delete or rewrite outside the project workspace.
- Do not run `rm -rf` on `input_pdbs/`, `processed_pdbs/`, `docked/`, or any parent directory.
- Treat PDB/SDF/MOL2 data as **read-only** unless the user explicitly asks for a write operation.
- Keep secrets and credentials out of logs and code.

---

## Code Patch Style

- Aim for **minimal, surgical edits** instead of wide refactors.
- Provide full replacement snippets (no diff markers like `+` / `-`).
- Include a brief header in replies with:
  1) **Sharper coding prompt (reformat)**
  2) **Assumed folder structure**
  3) **Next action**
- Add concise in-code comments at key anchors explaining *why* a change is safe or necessary.
- Prefer centralizing logic in existing helpers (path router, ligand enumeration, logging helpers) over adding new ad-hoc utilities.

---

## Logging Expectations

Log decisions and side effects. Use levels: DEBUG, INFO, WARNING, ERROR.

- Always create per-protein log files, even on early failures.
- Avoid duplicate handlers and double emission; prefer central logging setup.
- Use structured log messages:

  `[level] [component] key1=val1 key2=val2 ... | msg`

Examples of components: `[path-router]`, `[receptor-prep]`, `[microstates]`, `[lig-prep]`, `[pocket]`, `[docking]`, `[pose-validate]`, `[eval]`.

Key guidelines:

- DEBUG – branch choices, inputs/outputs, counts.
- INFO – stage start/finish, summary counts, output paths.
- WARNING – fallbacks and recoverable anomalies; include `reason=` and `fallback=`.
- ERROR – unrecoverable for the current unit (pdb/ligand/stage); include `reason=` and `action=skipped`.

---

## Path Router & Variants

- The **path router** is the canonical source of truth for where things live.
- Always use router helpers to get paths instead of building strings manually.
- Respect APO vs HOLO semantics and pH ensemble variants.
- If APO and HOLO reduce to identical receptors, treat as ligand-free APO and follow the current dedup policy rather than duplicating work.

---

## pH Ensemble Rules

- Global mode default: allow whole-protein PROPKA-guided renaming; center not required.
- If a radius is specified without a center, this is an error; emit a clear, actionable message.
- Microstate selection must be explicit and logged (`[microstates]` component).

---

## Allowed Commands (examples)

Safe default commands (adjust flags as needed):

- `micromamba run -n docking-env python main.py --pdbs <PDB> --stages prep,dock --test-mode`
- `micromamba run -n docking-env python ph_ensemble.py --pdb <PDB> --ph-list 6.0 8.0 --global`
- `micromamba run -n docking-env python analysis/dud_eval.py --docked-root docked/<PDB>`
- `micromamba run -n docking-env pytest tests/test_microstates.py -q`
- `micromamba run -n docking-env pytest tests/test_ions_acceptance.py -q` (if tests are failing and the user wants to investigate)

All of these should run from the workspace root (the current directory).

---

## Forbidden Actions

- Installing packages or modifying the `docking-env` environment.
- Accessing the internet or external services.
- Editing or deleting files outside the workspace root.
- Editing cluster configs, scheduler files, or user home directories.

---

## Patch Boundaries & Protected Paths (Must Read)

**Goal:** Prevent patch collisions between automation and local/CI setup by strictly scoping what the agent may edit.

### Allowed default patch scope

Focus patches in:

- `code/protein_automation/**` (Python pipeline code)
- `analysis/**` (non-test analysis helpers)
- `docs/**` (text only)

### Protected paths (DO NOT EDIT unless explicitly requested)

- `tests/**`
- `tests/data/**` – cached artifacts only; never commit binaries here.
- `.github/workflows/**`
- `ci/**`
- `input_pdbs/**` – source inputs; never auto-rewrite.
- `processed_pdbs/**`, `docked/**` – build outputs; hands off.
- `AGENTS.md` – only edit when the user explicitly asks to change agent policy.

**If your patch touches any protected path, abort** unless the user’s prompt clearly includes an override tag such as `override:tests`, `override:ci`, or `override:agents`.

### Binary / data policy

- Never vendor PDB, SDF, MOL2, or other binaries into the repo.
- Tests must **download/cache** required artifacts under `tests/data/` and rely on `.gitignore`.
- If a required file is missing, adjust tests to fetch or skip with a clear message, rather than adding data files to git.

### Environment & CI policy

- Prefer using `ci/run_in_env.sh` when the user asks for CI-style runs.
- Do not modify `.github/workflows/**` unless the user explicitly requests CI changes (`override:ci`).
- Do not add new package installs in code; instead, update `environment.yml` and rely on CI/bootstrap scripts (with user approval).

### Test policy

- Do not edit `tests/**` unless the user requests changes (`override:tests`).
- Acceptance tests (e.g., ions retention, pH ensemble behavior, microstates) are the source of truth; modify pipeline code to satisfy them.

### Conflict-avoidance with setup scripts

- If a patch would change files that bootstrap the env or CI and that the user may also be editing, prefer adding logging or flags over restructuring shared bootstrapping code.

### PR checklist (agent MUST enforce)

- [ ] No changes in protected paths (unless override tag present).
- [ ] No binaries added to the repo.
- [ ] Patch is single-topic and minimal.
- [ ] New/changed code emits clear logs for key decisions (grep-able tags).
- [ ] Tests pass locally via `ci/run_in_env.sh` or agreed smoke tests.

### Examples

- **OK:** Add debug lines in `code/protein_automation/automate_protein_prep.py` and `main.py`.
- **NOT OK:** Edit `tests/test_ions_acceptance.py` “just to make it pass” unless explicitly asked with `override:tests`.

