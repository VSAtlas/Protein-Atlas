# Atlas Log Resurfacer - Summary

*Fatal:* 0  *Hard:* 8  *Soft:* 27  *Info:* 0

**Health score:** 0 / 100

## 1. [HARD] ligprep.rdkit_valence - 2025
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 7  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_203904.log`
    - `[20:40:24] Explicit valence for atom # 10 O, 3, is greater than permitted`
    - `[20:40:25] Explicit valence for atom # 10 O, 3, is greater than permitted`
    - `[20:40:25] Explicit valence for atom # 7 O, 3, is greater than permitted`

## 2. [HARD] ligprep.mgl_partial_write - 2025
- **Message:** Partial write during conversion
- **Count:** 1  |  **Recency:** 
- **Hint:** Conversion partial; capture stderr, retry fallback, then quarantine if persistent.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_203904.log`
    - `WARNING: 72 atoms of 77 in CDL_N604.sanitized.mol2  were not written`

## 3. [SOFT] render.enqueue_none - 2025
- **Message:** Controls/RDK not enqueued
- **Count:** 12  |  **Recency:** 2025-10-08 20:43:01
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_203904.log`
    - `2025-10-08 20:41:26,523 - INFO - [render-enqueue] 6O0K bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K_NOLIG/receptor/6O0K_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0K.pdb | ctrl=None | rdk=None`
    - `2025-10-08 20:41:29,905 - INFO - [render-enqueue] 6O0L bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L_NOLIG/receptor/6O0L_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0L.pdb | ctrl=None | rdk=None`
    - `2025-10-08 20:41:51,671 - INFO - [render-enqueue] 4XV2 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2_NOLIG/receptor/4XV2_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XV2.pdb | ctrl=None | rdk=None`

## 4. [SOFT] ligprep.obabel_h_charge - 6ADQ
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 3  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_203904.log`
    - `[20:40:40] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9YF_C303.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:40:59] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/PLM_W203.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:41:05] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_N604.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 5. [SOFT] render.enqueue_none - 4XUF
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 20:43:01
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4XUF/protein.log`
    - `2025-10-08 20:43:01,066 - INFO - [render-enqueue] 4XUF bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF_NOLIG/receptor/4XUF_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XUF.pdb | ctrl=None | rdk=None`

## 6. [SOFT] render.enqueue_none - 3ZOS
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 20:42:54
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3ZOS/protein.log`
    - `2025-10-08 20:42:54,313 - INFO - [render-enqueue] 3ZOS bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS_NOLIG/receptor/3ZOS_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ZOS.pdb | ctrl=None | rdk=None`

## 7. [SOFT] render.enqueue_none - 2HYY
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 20:42:53
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2HYY/protein.log`
    - `2025-10-08 20:42:53,976 - INFO - [render-enqueue] 2HYY bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY_NOLIG/receptor/2HYY_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2HYY.pdb | ctrl=None | rdk=None`

## 8. [SOFT] render.enqueue_none - 4RT7
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 20:42:33
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4RT7/protein.log`
    - `2025-10-08 20:42:33,394 - INFO - [render-enqueue] 4RT7 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7_NOLIG/receptor/4RT7_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4RT7.pdb | ctrl=None | rdk=None`

## 9. [SOFT] render.enqueue_none - 2WGJ
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 20:42:30
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2WGJ/protein.log`
    - `2025-10-08 20:42:30,680 - INFO - [render-enqueue] 2WGJ bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ_NOLIG/receptor/2WGJ_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2WGJ.pdb | ctrl=None | rdk=None`

## 10. [SOFT] render.enqueue_none - 5L7I
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 20:42:23
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5L7I/protein.log`
    - `2025-10-08 20:42:23,880 - INFO - [render-enqueue] 5L7I bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I_NOLIG/receptor/5L7I_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5L7I.pdb | ctrl=None | rdk=None`

## 11. [SOFT] render.enqueue_none - 5L2I
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 20:42:08
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5L2I/protein.log`
    - `2025-10-08 20:42:08,771 - INFO - [render-enqueue] 5L2I bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L2I_NOLIG/receptor/5L2I_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5L2I.pdb | ctrl=None | rdk=None`

## 12. [SOFT] render.enqueue_none - 2GQG
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 20:42:02
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2GQG/protein.log`
    - `2025-10-08 20:42:02,813 - INFO - [render-enqueue] 2GQG bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG_NOLIG/receptor/2GQG_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2GQG.pdb | ctrl=None | rdk=None`

## 13. [SOFT] render.enqueue_none - 4XV2
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 20:41:51
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4XV2/protein.log`
    - `2025-10-08 20:41:51,671 - INFO - [render-enqueue] 4XV2 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2_NOLIG/receptor/4XV2_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XV2.pdb | ctrl=None | rdk=None`

## 14. [SOFT] render.enqueue_none - 5MO4
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 20:41:51
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5MO4/protein.log`
    - `2025-10-08 20:41:51,811 - INFO - [render-enqueue] 5MO4 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4_NOLIG/receptor/5MO4_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5MO4.pdb | ctrl=None | rdk=None`

## 15. [SOFT] render.enqueue_none - 6O0L
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 20:41:29
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6O0L/protein.log`
    - `2025-10-08 20:41:29,905 - INFO - [render-enqueue] 6O0L bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L_NOLIG/receptor/6O0L_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0L.pdb | ctrl=None | rdk=None`

## 16. [SOFT] render.enqueue_none - 6O0K
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 20:41:26
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6O0K/protein.log`
    - `2025-10-08 20:41:26,523 - INFO - [render-enqueue] 6O0K bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K_NOLIG/receptor/6O0K_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0K.pdb | ctrl=None | rdk=None`

