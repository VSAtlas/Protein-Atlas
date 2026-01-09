# Agent Profile: Atlas2 (Local Docking Pipeline)

## 1. Critical Token & Safety Rules (READ FIRST)
**You are strictly strictly bound by these token-saving rules:**
1.  **NO MASSIVE READS:** Never use `cat`, `read`, or `type` on data files (`.csv`, `.pdb`, `.sdf`, `.mol2`, `.pdbqt`, `.log`) without checking size first.
    * *Allowed:* `head -n 5 filename` to check headers/structure.
    * *Allowed:* `grep` or `rg` to find specific lines.
2.  **NO TREE DUMPS:** Do not run `tree`, `ls -R`, or `find .` on the root. Use targeted `ls` only.
3.  **SCOPE:** Edit ONLY `code/protein_automation/` and `docs/`.
    * **PROTECTED:** `tests/`, `ci/`, `input_pdbs/`, `docked/`. (Override required: `override:tests`).
4.  **ENVIRONMENT:** Always use: `micromamba run -n docking-env ...`
    * Ignore `libmamba` lockfile errors; assume success if retry happens.

## 2. Repository Architecture
**Root:** `/home/michael/atlas/code/protein_automation`
**Source:** `src/` (Packages: `docking`, `path_router`, `cli`, `prep_ligands`, `post_docking`).
**Inputs:** `input_pdbs/`, `processed_pdbs/` (Read-Only).
**Outputs:** `docked/<RUN_ID>/...`, `post_docked/<RUN_ID>/...` (Engine-scoped).

**Pathing Rule:**
* Do NOT hardcode paths. Import and use `src.path_router.path_router`.
* Respect `APO` vs `HOLO` and `pH` variants.

## 3. Docking Engine Standards
* **Isolation:** Output to engine-specific folders to avoid overwrites.
    * Vina: `docked/.../stage1`
    * GNINA: `docked/.../gnina_stage1` (CSVs: `gnina_docking_score_summary.csv`)
    * LeDock: `docked/.../ledock_stage1` (Requires MOL2 libs)
    * Decoys (DUD): `docked/.../dud_stage1` or `gnina_dud_stage1`.
* **SCORCH (Rescoring):**
    * Must handle FDA and DUD modes separately.
    * DUD runs must strip `_dud` suffixes from basenames before filtering.

## 4. Coding & Logging Standards
* **Style:** Minimal, surgical patches. No wide refactors.
* **Logging:** Use structured tags for grep-ability:
    * `[docking]`, `[path-router]`, `[gnina]`, `[scorch]`.
    * Levels: DEBUG (verbose), INFO (status), WARNING (fallback), ERROR (skip).
* **DeepCoy:** Only run autogen if `TEST_LIBRARY_MAP` fails or path is missing.
Always lint your code with ruff before showing it to me

## 5. Common Commands (Reference)
* **Smoke Test:** `micromamba run -n docking-env python main.py --pdbs 1BN1 --stages prep,dock --test-mode`
* **pH Ensemble:** `micromamba run -n docking-env python ph_ensemble.py --pdb 1BN1 --ph-list 6.0 8.0 --global`
* **Rescoring:** `micromamba run -n docking-env python -m post_docking.rescoring.rescoring_scorch --run-id <ID> --overwrite`
* **Tests:** `micromamba run -n docking-env pytest chemdb/tests/test_path_router.py -q`

## 6. Code Exploration Strategy
**You must minimize token usage when learning the codebase.**
1.  **Map First:** Read `REPO_MAP.txt` to see which files exist.
    * *Do NOT run `ls -R`, `find .`, or `tree`.*
2.  **Skeletons Second:** When exploring a new module, read the **skeleton file** first.
    * *Source:* `src/docking/docking.py` (Expensive, 4000 tokens)
    * *Target:* `_skeletons/src/docking/docking.py` (Cheap, 400 tokens)
3.  **Full Source Last:** Only read the real `src/` file if you are **editing** it or debugging a specific crash.
4.  **Search:** Use `grep -r "pattern" _skeletons/` to find function definitions cheaply.