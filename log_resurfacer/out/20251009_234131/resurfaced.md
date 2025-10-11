# Atlas Log Resurfacer - Summary

*Fatal:* 0  *Hard:* 6  *Soft:* 54  *Info:* 0

**Health score:** 0 / 100

## 1. [HARD] ligprep.rdkit_valence - 2025
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 5  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251009_234131.log`
    - `[23:42:58] Explicit valence for atom # 10 O, 3, is greater than permitted`
    - `[23:43:20] Explicit valence for atom # 10 O, 3, is greater than permitted`
    - `[23:43:23] Explicit valence for atom # 7 O, 3, is greater than permitted`

## 2. [HARD] ligprep.mgl_partial_write - 2025
- **Message:** Partial write during conversion
- **Count:** 1  |  **Recency:** 
- **Hint:** Conversion partial; capture stderr, retry fallback, then quarantine if persistent.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251009_234131.log`
    - `WARNING: 88 atoms of 98 in CDL_A502.sanitized.mol2  were not written`

## 3. [SOFT] render.enqueue_none - 2025
- **Message:** Controls/RDK not enqueued
- **Count:** 25  |  **Recency:** 2025-10-09 23:49:47
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251009_234131.log`
    - `2025-10-09 23:45:50,782 - INFO - [render-enqueue] 6O0L bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/receptor/6O0L_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0L.pdb | ctrl=None | rdk=None`
    - `2025-10-09 23:45:56,087 - INFO - [render-enqueue] 3WZE bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/receptor/3WZE_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3WZE.pdb | ctrl=None | rdk=None`
    - `2025-10-09 23:45:59,019 - INFO - [render-enqueue] 6O0K bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/receptor/6O0K_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0K.pdb | ctrl=None | rdk=None`

## 4. [SOFT] ligprep.obabel_h_charge - 6ADQ
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251009_234131.log`
    - `[23:43:00] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_O304.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[23:43:25] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_A502.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 5. [SOFT] render.enqueue_none - 5L2I
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:49:47
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5L2I/protein.log`
    - `2025-10-09 23:49:47,863 - INFO - [render-enqueue] 5L2I bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L2I/receptor/5L2I_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5L2I.pdb | ctrl=None | rdk=None`

## 6. [SOFT] render.enqueue_none - 3G0E
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:49:21
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3G0E/protein.log`
    - `2025-10-09 23:49:21,543 - INFO - [render-enqueue] 3G0E bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3G0E/receptor/3G0E_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3G0E.pdb | ctrl=None | rdk=None`

## 7. [SOFT] render.enqueue_none - 3OXZ
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:49:18
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3OXZ/protein.log`
    - `2025-10-09 23:49:18,165 - INFO - [render-enqueue] 3OXZ bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/receptor/3OXZ_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3OXZ.pdb | ctrl=None | rdk=None`

## 8. [SOFT] render.enqueue_none - 5I96
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:49:18
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5I96/protein.log`
    - `2025-10-09 23:49:18,843 - INFO - [render-enqueue] 5I96 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/receptor/5I96_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5I96.pdb | ctrl=None | rdk=None`

## 9. [SOFT] render.enqueue_none - 2XP2
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:49:14
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2XP2/protein.log`
    - `2025-10-09 23:49:14,633 - INFO - [render-enqueue] 2XP2 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/receptor/2XP2_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2XP2.pdb | ctrl=None | rdk=None`

## 10. [SOFT] render.enqueue_none - 3LXK
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:49:11
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3LXK/protein.log`
    - `2025-10-09 23:49:11,584 - INFO - [render-enqueue] 3LXK bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/receptor/3LXK_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3LXK.pdb | ctrl=None | rdk=None`

## 11. [SOFT] render.enqueue_none - 3ZOS
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:49:11
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3ZOS/protein.log`
    - `2025-10-09 23:49:11,696 - INFO - [render-enqueue] 3ZOS bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/receptor/3ZOS_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ZOS.pdb | ctrl=None | rdk=None`

## 12. [SOFT] render.enqueue_none - 4RT7
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:49:10
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4RT7/protein.log`
    - `2025-10-09 23:49:10,636 - INFO - [render-enqueue] 4RT7 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/receptor/4RT7_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4RT7.pdb | ctrl=None | rdk=None`

## 13. [SOFT] render.enqueue_none - 4XUF
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:49:09
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4XUF/protein.log`
    - `2025-10-09 23:49:09,954 - INFO - [render-enqueue] 4XUF bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/receptor/4XUF_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XUF.pdb | ctrl=None | rdk=None`

## 14. [SOFT] render.enqueue_none - 2WGJ
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:49:08
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2WGJ/protein.log`
    - `2025-10-09 23:49:08,615 - INFO - [render-enqueue] 2WGJ bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/receptor/2WGJ_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2WGJ.pdb | ctrl=None | rdk=None`

## 15. [SOFT] render.enqueue_none - 6JQR
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:49:05
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6JQR/protein.log`
    - `2025-10-09 23:49:05,079 - INFO - [render-enqueue] 6JQR bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/receptor/6JQR_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6JQR.pdb | ctrl=None | rdk=None`

## 16. [SOFT] render.enqueue_none - 4U5J
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:49:00
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4U5J/protein.log`
    - `2025-10-09 23:49:00,724 - INFO - [render-enqueue] 4U5J bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/receptor/4U5J_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4U5J.pdb | ctrl=None | rdk=None`

## 17. [SOFT] render.enqueue_none - 3ERT
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:48:45
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3ERT/protein.log`
    - `2025-10-09 23:48:45,939 - INFO - [render-enqueue] 3ERT bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/receptor/3ERT_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ERT.pdb | ctrl=None | rdk=None`

## 18. [SOFT] render.enqueue_none - 2HYY
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:47:14
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2HYY/protein.log`
    - `2025-10-09 23:47:14,965 - INFO - [render-enqueue] 2HYY bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/receptor/2HYY_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2HYY.pdb | ctrl=None | rdk=None`

## 19. [SOFT] render.enqueue_none - 2HYY
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:47:14
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251009_234131.log`
    - `Docking (bench_pocket1_single):  78%|███████▊  | 7/9 [00:08<00:01,  1.05ligand/s][A2025-10-09 23:47:14,965 - INFO - [render-enqueue] 2HYY bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/receptor/2HYY_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2HYY.pdb | ctrl=None | rdk=None`

## 20. [SOFT] render.enqueue_none - 3CS9
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:47:09
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3CS9/protein.log`
    - `2025-10-09 23:47:09,306 - INFO - [render-enqueue] 3CS9 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/receptor/3CS9_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3CS9.pdb | ctrl=None | rdk=None`

## 21. [SOFT] render.enqueue_none - 5L7I
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:46:55
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5L7I/protein.log`
    - `2025-10-09 23:46:55,630 - INFO - [render-enqueue] 5L7I bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/receptor/5L7I_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5L7I.pdb | ctrl=None | rdk=None`

## 22. [SOFT] render.enqueue_none - 5MO4
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:46:54
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5MO4/protein.log`
    - `2025-10-09 23:46:54,271 - INFO - [render-enqueue] 5MO4 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/receptor/5MO4_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5MO4.pdb | ctrl=None | rdk=None`

## 23. [SOFT] render.enqueue_none - 2GQG
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:46:52
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2GQG/protein.log`
    - `2025-10-09 23:46:52,775 - INFO - [render-enqueue] 2GQG bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/receptor/2GQG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2GQG.pdb | ctrl=None | rdk=None`

## 24. [SOFT] render.enqueue_none - 1T46
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:46:50
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-09 23:46:50,597 - INFO - [render-enqueue] 1T46 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/receptor/1T46_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/1T46.pdb | ctrl=None | rdk=None`

## 25. [SOFT] render.enqueue_none - 3WZD
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:46:20
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3WZD/protein.log`
    - `2025-10-09 23:46:20,607 - INFO - [render-enqueue] 3WZD bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/receptor/3WZD_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3WZD.pdb | ctrl=None | rdk=None`

## 26. [SOFT] render.enqueue_none - 4AG8
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:46:19
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4AG8/protein.log`
    - `2025-10-09 23:46:19,537 - INFO - [render-enqueue] 4AG8 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/receptor/4AG8_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4AG8.pdb | ctrl=None | rdk=None`

## 27. [SOFT] render.enqueue_none - 4XV2
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:46:15
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4XV2/protein.log`
    - `2025-10-09 23:46:15,086 - INFO - [render-enqueue] 4XV2 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/receptor/4XV2_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XV2.pdb | ctrl=None | rdk=None`

## 28. [SOFT] render.enqueue_none - 4ASD
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:46:14
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4ASD/protein.log`
    - `2025-10-09 23:46:14,414 - INFO - [render-enqueue] 4ASD bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/receptor/4ASD_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4ASD.pdb | ctrl=None | rdk=None`

## 29. [SOFT] render.enqueue_none - 6O0K
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:45:59
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6O0K/protein.log`
    - `2025-10-09 23:45:59,019 - INFO - [render-enqueue] 6O0K bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/receptor/6O0K_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0K.pdb | ctrl=None | rdk=None`

## 30. [SOFT] render.enqueue_none - 3WZE
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:45:56
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3WZE/protein.log`
    - `2025-10-09 23:45:56,087 - INFO - [render-enqueue] 3WZE bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/receptor/3WZE_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3WZE.pdb | ctrl=None | rdk=None`

## 31. [SOFT] render.enqueue_none - 6O0L
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 23:45:50
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6O0L/protein.log`
    - `2025-10-09 23:45:50,782 - INFO - [render-enqueue] 6O0L bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/receptor/6O0L_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0L.pdb | ctrl=None | rdk=None`

