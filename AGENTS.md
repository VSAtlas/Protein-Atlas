# Agent Profile: Atlas2 (Local Docking Pipeline)

## 1. CRITICAL SAFETY & TOKEN RULES
* **NO MASSIVE READS:** Check file size first. Use `head -n 5` or `grep`, NEVER `cat` data files (`.pdb`, `.sdf`, `.log`).
* **NO TREE DUMPS:** Do not run `ls -R`, `find .`, or `tree`. Use `REPO_MAP.txt`.
* **SCOPE:** Edit ONLY `code/protein_automation/` and `docs/`. Protected: `tests/`, `ci/`, `input_pdbs/`, `docked/`.
* **ENVIRONMENT:** Always run: `micromamba run -n docking-env ...` (Ignore lockfile errors; assume success). ... ask for elevated permissions if you run  into repeated  failures

## 2. REPO ARCHITECTURE
* **Root:** `/home/michael/atlas/code/protein_automation`
* **Source:** `src/` (Packages: `docking`, `path_router`, `cli`, `prep_ligands`, `post_docking`).
* **Outputs:** `docked/<RUN_ID>/...` (Engine-scoped).
* **Pathing:** NEVER hardcode. Import `src.path_router.path_router`. Respect `APO`/`HOLO` and `pH`.

## 3. CODING & EXECUTION STANDARDS
* **Style:** Surgical patches only. No wide refactors.
* **Linting:** You MUST run `ruff check` on every file you touch before finishing.
* **Verification:** You MUST run `git diff` (or `cat` changed lines) to verify edits immediately.
* **DeepCoy:** Only run autogen if `TEST_LIBRARY_MAP` fails.
* **Isolation:**
    * Vina: `docked/.../stage1`
    * GNINA: `docked/.../gnina_stage1`
    * LeDock: `docked/.../ledock_stage1`

## 4. CODE EXPLORATION (LOW-TOKEN)
1.  **Map First:** Read `REPO_MAP.txt` to locate files.
2.  **Skeletons Second:** Read `_skeletons/src/...` to learn functions/args cheaply (400 tokens).
3.  **Source Last:** Only read real `src/` files if editing or debugging a specific crash.
4.  **Search:** Use `grep -r "pattern" _skeletons/`.

## 5. DEBUGGING TOOLKIT (LOGS)
**Do NOT read raw logs. Use the summarizer.**

* **Standard Crash:**
    `micromamba run -n docking-env python analysis/summarize_run_errors.py <RUN_ID> --recent-hours 48 --top 5`
* **Deep Scan (Hidden Errors):**
    `... python analysis/summarize_run_errors.py <RUN_ID> --deep --top 10`
* **Reporting:** Paste ONLY the "signature blocks" (count + error line). Never paste full log files into chat.

## 6. CHEATSHEET (Common Commands)
* **Smoke Test:** `... python main.py --test -fast --test-fda`
* **Rescoring:** `... python -m post_docking.rescoring.rescoring_scorch --run-id <ID> --overwrite`
* **Tests:** `... pytest chemdb/tests/test_path_router.py -q`