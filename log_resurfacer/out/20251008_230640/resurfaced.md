# Atlas Log Resurfacer - Summary

*Fatal:* 0  *Hard:* 9  *Soft:* 26  *Info:* 0

**Health score:** 0 / 100

## 1. [HARD] ligprep.rdkit_valence - 2025
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 8  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_230640.log`
    - `[23:08:02] Explicit valence for atom # 10 O, 3, is greater than permitted`
    - `[23:08:02] Explicit valence for atom # 10 O, 3, is greater than permitted`
    - `[23:08:03] Explicit valence for atom # 7 O, 3, is greater than permitted`

## 2. [HARD] ligprep.mgl_partial_write - 2025
- **Message:** Partial write during conversion
- **Count:** 1  |  **Recency:** 
- **Hint:** Conversion partial; capture stderr, retry fallback, then quarantine if persistent.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_230640.log`
    - `WARNING: 79 atoms of 84 in CDL_F606.sanitized.sanitized.mol2  were not written`

## 3. [SOFT] render.enqueue_none - 2025
- **Message:** Controls/RDK not enqueued
- **Count:** 12  |  **Recency:** 2025-10-08 23:11:13
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_230640.log`
    - `2025-10-08 23:09:31,307 - INFO - [render-enqueue] 6O0K bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K_NOLIG/receptor/6O0K_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0K.pdb | ctrl=None | rdk=None`
    - `2025-10-08 23:09:33,836 - INFO - [render-enqueue] 6O0L bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L_NOLIG/receptor/6O0L_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0L.pdb | ctrl=None | rdk=None`
    - `2025-10-08 23:09:51,869 - INFO - [render-enqueue] 4XV2 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2_NOLIG/receptor/4XV2_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XV2.pdb | ctrl=None | rdk=None`

## 4. [SOFT] ligprep.obabel_h_charge - 6ADQ
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_230640.log`
    - `[23:08:37] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9YF_M504.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[23:08:44] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_F606.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 5. [SOFT] render.enqueue_none - 3ZOS
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 23:11:13
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3ZOS/protein.log`
    - `2025-10-08 23:11:13,288 - INFO - [render-enqueue] 3ZOS bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS_NOLIG/receptor/3ZOS_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ZOS.pdb | ctrl=None | rdk=None`

## 6. [SOFT] render.enqueue_none - 4RT7
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 23:11:06
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4RT7/protein.log`
    - `2025-10-08 23:11:06,617 - INFO - [render-enqueue] 4RT7 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7_NOLIG/receptor/4RT7_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4RT7.pdb | ctrl=None | rdk=None`

## 7. [SOFT] render.enqueue_none - 4XUF
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 23:10:59
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4XUF/protein.log`
    - `2025-10-08 23:10:59,763 - INFO - [render-enqueue] 4XUF bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF_NOLIG/receptor/4XUF_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XUF.pdb | ctrl=None | rdk=None`

## 8. [SOFT] render.enqueue_none - 2WGJ
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 23:10:57
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2WGJ/protein.log`
    - `2025-10-08 23:10:57,852 - INFO - [render-enqueue] 2WGJ bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ_NOLIG/receptor/2WGJ_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2WGJ.pdb | ctrl=None | rdk=None`

## 9. [SOFT] render.enqueue_none - 5L2I
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 23:10:57
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5L2I/protein.log`
    - `2025-10-08 23:10:57,416 - INFO - [render-enqueue] 5L2I bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L2I_NOLIG/receptor/5L2I_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5L2I.pdb | ctrl=None | rdk=None`

## 10. [SOFT] render.enqueue_none - 2HYY
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 23:10:28
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2HYY/protein.log`
    - `2025-10-08 23:10:28,945 - INFO - [render-enqueue] 2HYY bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY_NOLIG/receptor/2HYY_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2HYY.pdb | ctrl=None | rdk=None`

## 11. [SOFT] render.enqueue_none - 5L7I
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 23:10:11
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5L7I/protein.log`
    - `2025-10-08 23:10:11,246 - INFO - [render-enqueue] 5L7I bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I_NOLIG/receptor/5L7I_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5L7I.pdb | ctrl=None | rdk=None`

## 12. [SOFT] render.enqueue_none - 5MO4
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 23:09:54
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5MO4/protein.log`
    - `2025-10-08 23:09:54,150 - INFO - [render-enqueue] 5MO4 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4_NOLIG/receptor/5MO4_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5MO4.pdb | ctrl=None | rdk=None`

## 13. [SOFT] render.enqueue_none - 2GQG
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 23:09:52
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2GQG/protein.log`
    - `2025-10-08 23:09:52,507 - INFO - [render-enqueue] 2GQG bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG_NOLIG/receptor/2GQG_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2GQG.pdb | ctrl=None | rdk=None`

## 14. [SOFT] render.enqueue_none - 4XV2
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 23:09:51
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4XV2/protein.log`
    - `2025-10-08 23:09:51,869 - INFO - [render-enqueue] 4XV2 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2_NOLIG/receptor/4XV2_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XV2.pdb | ctrl=None | rdk=None`

## 15. [SOFT] render.enqueue_none - 6O0L
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 23:09:33
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6O0L/protein.log`
    - `2025-10-08 23:09:33,836 - INFO - [render-enqueue] 6O0L bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L_NOLIG/receptor/6O0L_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0L.pdb | ctrl=None | rdk=None`

## 16. [SOFT] render.enqueue_none - 6O0K
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 23:09:31
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6O0K/protein.log`
    - `2025-10-08 23:09:31,307 - INFO - [render-enqueue] 6O0K bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K_NOLIG/receptor/6O0K_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0K.pdb | ctrl=None | rdk=None`

