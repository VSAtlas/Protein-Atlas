# Agent Profile: Atlas2 (Local Docking Pipeline)

## 1. CRITICAL SAFETY & TOKEN RULES
* **NO MASSIVE READS:** Check file size first. Use `head -n 5` or `csvstat`, NEVER `cat` data files (`.pdb`, `.sdf`, `.csv`, `.log`).
* **NO TREE DUMPS:** Do not run `ls -R`, `find .`, or `tree`. Use `REPO_MAP.txt`.
* **SCOPE:** Edit ONLY `code/protein_automation/` and `docs/`. Protected: `tests/`, `ci/`, `input_pdbs/`, `docked/`.
* **ENVIRONMENT:** `direnv` is active. You can run `python` directly. (Ignore lockfile errors; assume success).
    * *Only ask for elevated permissions if you run into repeated failures.*

## 2. REPO ARCHITECTURE
* **Root:** `/home/michael/atlas/code/protein_automation`
* **Source:** `src/` (Packages: `docking`, `path_router`, `cli`, `prep_ligands`, `post_docking`).
* **Outputs:** `docked/<RUN_ID>/...` (Engine-scoped).
* **Pathing:** NEVER hardcode. Import `src.path_router.path_router`. Respect `APO`/`HOLO` and `pH`.

## 3. CODING & EXECUTION STANDARDS
* **Style:** Surgical patches only. No wide refactors.
* **Safety Stack (MANDATORY):**
    1. **Lint:** Run `ruff check --fix .` to clean syntax.
    2. **Types:** Run `mypy <filename>` to check logic errors.
    3. **Test:** Run `pytest --cov=src tests/` to verify logic + coverage.
* **Commit Rule:** The repo uses `pre-commit`. You CANNOT commit code that fails Ruff or MyPy. Run them manually first to ensure success.
* **Verification:** Run `git diff` immediately after editing to verify changes.
* **Isolation:** Output to `docked/.../stage1` (etc) to avoid overwrites.

## 4. CODE EXPLORATION (LOW-TOKEN)
1.  **Map First:** Read `REPO_MAP.txt` to locate files.
2.  **Skeletons Second:** Read `_skeletons/src/...` to learn functions/args cheaply (400 tokens).
3.  **Source Last:** Only read real `src/` files if editing or debugging a specific crash.
4.  **Search:** Use `rg "pattern" _skeletons/` (ripgrep) to find definitions cheaply.

## 5. DEBUGGING TOOLKIT & DATA SAFETY
**Log Rules:** Do NOT read raw logs. Use the summarizer.
* **Crash Scan:** `python analysis/summarize_run_errors.py <RUN_ID> --recent-hours 48 --top 5`
* **Deep Scan:** `python analysis/summarize_run_errors.py <RUN_ID> --deep --top 10`

**Data Rules:** Do NOT `cat` CSVs.
* **Inspect:** `csvstat <file>` (Stats/Columns) or `csvlook <file> | head` (Table View).
* **Filter:** `csvgrep -c "col_name" -r "pattern" <file>` to extract specific rows.

## 6. CHEATSHEET (Common Commands)
* **Smoke Test:** `python main.py --test -fast --test-fda`
* **Rescoring:** `python -m post_docking.rescoring.rescoring_scorch --run-id <ID> --overwrite`
* **Tests:** `pytest chemdb/tests/test_path_router.py -q`