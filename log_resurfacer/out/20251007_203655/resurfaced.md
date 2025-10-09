# Atlas Log Resurfacer — Summary

*Fatal:* 0  *Hard:* 30  *Soft:* 6  *Info:* 0

**Health score:** 0 / 100

## 1. [HARD] ligprep.rdkit_valence — 2025
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 6  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_203655.log`
    - `[RDKit] ERROR:[20:37:07] Explicit valence for atom # 18 He, 1, is greater than permitted`
    - `rdkit.Chem.rdchem.AtomValenceException: Explicit valence for atom # 18 He, 1, is greater than permitted`
    - `[RDKit] ERROR:[20:37:08] Explicit valence for atom # 18 He, 1, is greater than permitted`

## 2. [HARD] dock.no_score — 1T46 — STI_A3.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-07 20:37:14,234 - WARNING - No score for STI_A3.pdbqt`

## 3. [HARD] dock.no_score — 1T46 — STI_A3.sanitized.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-07 20:37:14,235 - WARNING - No score for STI_A3.sanitized.pdbqt`
- **Joins:**
  - ligand=STI_A3.sanitized.pdbqt prep_status=ok reason=

## 4. [HARD] dock.no_score — 1T46 — rdk_0000010.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-07 20:37:14,273 - WARNING - No score for rdk_0000010.pdbqt`

## 5. [HARD] dock.no_score — 1T46 — rdk_0002911.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-07 20:37:14,273 - WARNING - No score for rdk_0002911.pdbqt`

## 6. [HARD] dock.no_score — 1T46 — rdk_0000043.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-07 20:37:14,273 - WARNING - No score for rdk_0000043.pdbqt`

## 7. [HARD] dock.no_score — 1T46 — rdk_0000067.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-07 20:37:14,273 - WARNING - No score for rdk_0000067.pdbqt`

## 8. [HARD] dock.no_score — 1T46 — rdk_0000076.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-07 20:37:14,275 - WARNING - No score for rdk_0000076.pdbqt`

## 9. [HARD] dock.no_score — 1T46 — rdk_0000119.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-07 20:37:14,277 - WARNING - No score for rdk_0000119.pdbqt`

## 10. [HARD] dock.no_score — 1T46 — rdk_0000123.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-07 20:37:14,281 - WARNING - No score for rdk_0000123.pdbqt`

## 11. [HARD] dock.no_score — 1T46 — rdk_0000160.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-07 20:37:14,281 - WARNING - No score for rdk_0000160.pdbqt`

## 12. [HARD] dock.no_score — 1T46 — rdk_0000125.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-07 20:37:14,283 - WARNING - No score for rdk_0000125.pdbqt`

## 13. [HARD] dock.no_score — 2025 — STI_A3.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_203655.log`
    - `2025-10-07 20:37:14,234 - WARNING - No score for STI_A3.pdbqt`

## 14. [HARD] dock.no_score — 2025 — STI_A3.sanitized.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_203655.log`
    - `2025-10-07 20:37:14,235 - WARNING - No score for STI_A3.sanitized.pdbqt`

## 15. [HARD] dock.no_score — 2025 — rdk_0000010.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_203655.log`
    - `Docking (bench_pocket1_single):   0%|          | 0/9 [00:00<?, ?ligand/s][A2025-10-07 20:37:14,273 - WARNING - No score for rdk_0000010.pdbqt`

## 16. [HARD] dock.no_score — 2025 — rdk_0002911.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_203655.log`
    - `2025-10-07 20:37:14,273 - WARNING - No score for rdk_0002911.pdbqt`

## 17. [HARD] dock.no_score — 2025 — rdk_0000043.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_203655.log`
    - `2025-10-07 20:37:14,273 - WARNING - No score for rdk_0000043.pdbqt`

## 18. [HARD] dock.no_score — 2025 — rdk_0000067.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_203655.log`
    - `2025-10-07 20:37:14,273 - WARNING - No score for rdk_0000067.pdbqt`

## 19. [HARD] dock.no_score — 2025 — rdk_0000076.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_203655.log`
    - `2025-10-07 20:37:14,275 - WARNING - No score for rdk_0000076.pdbqt`

## 20. [HARD] dock.no_score — 2025 — rdk_0000119.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_203655.log`
    - `2025-10-07 20:37:14,277 - WARNING - No score for rdk_0000119.pdbqt`

## 21. [HARD] dock.no_score — 2025 — rdk_0000123.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_203655.log`
    - `2025-10-07 20:37:14,281 - WARNING - No score for rdk_0000123.pdbqt`

## 22. [HARD] dock.no_score — 2025 — rdk_0000160.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_203655.log`
    - `2025-10-07 20:37:14,281 - WARNING - No score for rdk_0000160.pdbqt`

## 23. [HARD] dock.no_score — 2025 — rdk_0000125.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_203655.log`
    - `2025-10-07 20:37:14,283 - WARNING - No score for rdk_0000125.pdbqt`

## 24. [HARD] ligprep.tsv_failure — 6WTN — RXT_A1204.sanitized.pdb
- **Message:** Ligand prep failed per TSV
- **Count:** 1  |  **Recency:** 
- **Hint:** Treat TSV as ground truth; review the 'reason' column and quarantine failing ligand.
  - **File:** `prepped_ligands/6WTN/ligand_prep_status.tsv`
    - `TSV failure for RXT_A1204.sanitized.pdb: status=fail reason=postcheck_fail_or_small`

## 25. [HARD] ligprep.tsv_failure — 6WTN — RXT_A1204.sanitized
- **Message:** Ligand prep failed per TSV
- **Count:** 1  |  **Recency:** 
- **Hint:** Treat TSV as ground truth; review the 'reason' column and quarantine failing ligand.
  - **File:** `prepped_ligands/6WTN/ligand_prep_status.tsv`
    - `TSV failure for RXT_A1204.sanitized: status=fail reason=postcheck_fail_or_small`

## 26. [SOFT] render.enqueue_none — 1T46
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Control selection failed or disabled; check ‘hint_count’ and control whitelist.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-07 20:37:14,087 - INFO - [render-enqueue] 1T46 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46_NOLIG/receptor/1T46_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/1T46.pdb | ctrl=None | rdk=None`

## 27. [SOFT] render.no_ctrl_rdk_paths — 1T46
- **Message:** Screenshots skipped: no ctrl/rdk pose paths
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Why no valid poses? Correlate with ‘No score’/validation; ensure control selection produced a pose.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-07 20:37:14,288 - INFO - [render-skip] 1T46 bench_pocket1_single: no ctrl/rdk pose paths → no screenshots will be generated.`

## 28. [SOFT] render.enqueue_none — 2025
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Control selection failed or disabled; check ‘hint_count’ and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_203655.log`
    - `2025-10-07 20:37:14,087 - INFO - [render-enqueue] 1T46 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46_NOLIG/receptor/1T46_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/1T46.pdb | ctrl=None | rdk=None`

## 29. [SOFT] render.no_ctrl_rdk_paths — 2025
- **Message:** Screenshots skipped: no ctrl/rdk pose paths
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** Why no valid poses? Correlate with ‘No score’/validation; ensure control selection produced a pose.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_203655.log`
    - `2025-10-07 20:37:14,288 - INFO - [render-skip] 1T46 bench_pocket1_single: no ctrl/rdk pose paths → no screenshots will be generated.`

## 30. [SOFT] recentering.skipped_zero — 1T46
- **Message:** Early recenter skipped: evaluated=0
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** No evaluated poses to guide recentering; widen first-pass or ensure at least one pose completes.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-07 20:37:14,284 - INFO - Early recenter skipped: evaluated=0 < threshold.`

## 31. [SOFT] recentering.skipped_zero — 2025
- **Message:** Early recenter skipped: evaluated=0
- **Count:** 1  |  **Recency:** 2025-10-07 20:37:14
- **Hint:** No evaluated poses to guide recentering; widen first-pass or ensure at least one pose completes.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_203655.log`
    - `2025-10-07 20:37:14,284 - INFO - Early recenter skipped: evaluated=0 < threshold.`

