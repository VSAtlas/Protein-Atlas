# Atlas Log Resurfacer - Summary

*Fatal:* 0  *Hard:* 45  *Soft:* 53  *Info:* 0

**Health score:** 0 / 100

## 1. [HARD] ligprep.rdkit_valence - 2025
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 23  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:30] Can't kekulize mol.  Unkekulized atoms: 19 20 21 22 23`
    - `[02:05:30] Explicit valence for atom # 14 C, 5, is greater than permitted`
    - `[02:05:31] Explicit valence for atom # 35 C, 5, is greater than permitted`

## 2. [HARD] ligprep.rdkit_valence - 3CS9
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] sanitize [02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.sanitized.sanitized.pdb: WARNING:root:[obabel primary stderr] 1 molecule converted`
    - `[02:05:34] sanitize [02:05:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_C600.sanitized.sanitized.pdb: WARNING:root:[obabel primary stderr] 1 molecule converted`
    - `[02:05:37] sanitize [02:05:37] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.sanitized.sanitized.pdb: Saved fixed PDB to /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.sanitized.pdb`

## 3. [HARD] ligprep.mgl_partial_write - 2025
- **Message:** Partial write during conversion
- **Count:** 3  |  **Recency:** 
- **Hint:** Conversion partial; capture stderr, retry fallback, then quarantine if persistent.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `WARNING: 20 atoms of 23 in RXT_A601.sanitized.sanitized.mol2  were not written`
    - `WARNING: 20 atoms of 23 in RXT_B601.sanitized.sanitized.mol2  were not written`
    - `WARNING: 73 atoms of 82 in CDL_T201.sanitized.mol2  were not written`

## 4. [HARD] ligprep.rdkit_valence - 6U4J
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] sanitize [02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.sanitized.sanitized.pdb: [02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_A1202.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:33] sanitize [02:05:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.sanitized.sanitized.pdb: WARNING:root:[obabel primary stderr] 1 molecule converted`

## 5. [HARD] ligprep.rdkit_valence - 6O0L
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] sanitize [02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.sanitized.sanitized.pdb: Saved fixed PDB to /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.sanitized.pdb`
    - `[02:05:33] sanitize [02:05:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_A301.sanitized.sanitized.pdb: WARNING:root:[obabel primary stderr] 1 molecule converted`

## 6. [HARD] ligprep.rdkit_valence - 4XV2
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] sanitize [02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.sanitized.pdb: WARNING:root:MGLTools failed for RXT_A601.sanitized.sanitized.mol2 (will try fallback):`
    - `[02:05:34] sanitize [02:05:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.sanitized.sanitized.pdb: WARNING:root:MGLTools failed for RXT_B601.sanitized.sanitized.mol2 (will try fallback):`

## 7. [HARD] ligprep.rdkit_valence - 2WGJ
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 2025-10-08 02:05:30
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:30] sanitize [02:05:30] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.sanitized.pdb: 2025-10-08 02:05:30,835 - INFO - [prep] running ligand prep (raw newer than prepped or first run).`

## 8. [HARD] ligprep.rdkit_valence - 3OG7
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 2025-10-08 02:05:30
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:30] sanitize [02:05:30] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.sanitized.pdb: 2025-10-08 02:05:30,860 - INFO - FORCE_REPROCESS=True | cleaned_exists=True receptor_exists=False`

## 9. [HARD] ligprep.rdkit_valence - 6O0K
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] sanitize [02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.sanitized.sanitized.pdb: WARNING:root:[obabel primary stderr] 1 molecule converted`

## 10. [HARD] ligprep.rdkit_valence - 5L2I
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] sanitize [02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L2I/ligands_raw/LQQ_A900.sanitized.sanitized.pdb: Saved fixed PDB to /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.sanitized.sanitized.pdb`

## 11. [HARD] ligprep.rdkit_valence - 3G0E
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] sanitize [02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3G0E/ligands_raw/B49_A9000.sanitized.sanitized.pdb: [02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/ligands_raw/AXI_A2000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 12. [HARD] ligprep.rdkit_valence - 5I96
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:32] sanitize [02:05:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.sanitized.sanitized.pdb: [02:05:32] DEPRECATION WARNING: please use MorganGenerator`

## 13. [HARD] ligprep.rdkit_valence - 3LXK
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:36] sanitize [02:05:36] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.sanitized.sanitized.pdb: Saved fixed PDB to /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD_NOLIG/work/4ASD_NOLIG_elemfix.pdb`

## 14. [HARD] ligprep.rdkit_valence - 2XP2
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:36] sanitize [02:05:36] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.sanitized.sanitized.pdb: Saved fixed PDB to /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.sanitized.pdb`

## 15. [HARD] ligprep.rdkit_valence - 5MO4
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:41] sanitize [02:05:41] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.sanitized.sanitized.pdb: WARNING:root:Phenix not available; skipping pdbtools polish`

## 16. [SOFT] ligprep.obabel_h_charge - 3WZD
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 5  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1203.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:35] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/LEV_A1201.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:38] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1204.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 17. [SOFT] ligprep.obabel_h_charge - 3CS9
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_C600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:37] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 18. [SOFT] ligprep.obabel_h_charge - 3ZOS
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 3  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1004.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_B1000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:36] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 19. [SOFT] ligprep.obabel_h_charge - 6U4J
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 3  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:31] sanitize [02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.sanitized.sanitized.pdb: [02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_A1202.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 20. [SOFT] ligprep.obabel_h_charge - 5I96
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 3  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:35] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_B503.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:38] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_A506.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 21. [SOFT] ligprep.obabel_h_charge - 5L7I
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 3  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:35] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1204.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:38] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_B1202.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:40] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1203.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 22. [SOFT] render.enqueue_none - 4U2P
- **Message:** Controls/RDK not enqueued
- **Count:** 2  |  **Recency:** 2025-10-08 02:15:11
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4U2P/protein.log`
    - `2025-10-08 02:04:57,254 - INFO - [render-enqueue] 4U2P bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U2P/receptor/4U2P_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4U2P.pdb | ctrl=None | rdk=None`
    - `2025-10-08 02:15:11,976 - INFO - [render-enqueue] 4U2P bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U2P_NOLIG/receptor/4U2P_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4U2P.pdb | ctrl=None | rdk=None`

## 23. [SOFT] ligprep.obabel_h_charge - 4XUF
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:30] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_B1001.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_A1001.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 24. [SOFT] ligprep.obabel_h_charge - 6O0K
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 25. [SOFT] ligprep.obabel_h_charge - 4U5J
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 26. [SOFT] ligprep.obabel_h_charge - 3WZE
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/DTT_A1203.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/BAX_A1201.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 27. [SOFT] ligprep.obabel_h_charge - 5MO4
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:41] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 28. [SOFT] ligprep.obabel_h_charge - 6WTN
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 29. [SOFT] ligprep.obabel_h_charge - 3G0E
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3G0E/ligands_raw/B49_A9000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:31] sanitize [02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3G0E/ligands_raw/B49_A9000.sanitized.sanitized.pdb: [02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/ligands_raw/AXI_A2000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 30. [SOFT] ligprep.obabel_h_charge - 6JQR
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/CXS_A1007.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 31. [SOFT] ligprep.obabel_h_charge - 6O0L
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_A301.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 32. [SOFT] ligprep.obabel_h_charge - 4XV2
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[02:05:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 33. [SOFT] render.enqueue_none - 2025
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 02:15:11
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `2025-10-08 02:15:11,976 - INFO - [render-enqueue] 4U2P bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U2P_NOLIG/receptor/4U2P_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4U2P.pdb | ctrl=None | rdk=None`

## 34. [SOFT] ligprep.obabel_h_charge - 2WGJ
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:30] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 35. [SOFT] ligprep.obabel_h_charge - 3OG7
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:30] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 36. [SOFT] ligprep.obabel_h_charge - 4ASD
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/ligands_raw/BAX_A1500.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 37. [SOFT] ligprep.obabel_h_charge - 3OXZ
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 38. [SOFT] ligprep.obabel_h_charge - 3ERT
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 39. [SOFT] ligprep.obabel_h_charge - 4RT7
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 40. [SOFT] ligprep.obabel_h_charge - 3LXK
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:36] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 41. [SOFT] ligprep.obabel_h_charge - 2XP2
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:36] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 42. [SOFT] ligprep.obabel_h_charge - 6ADQ
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020419.log`
    - `[02:05:39] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_T201.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

