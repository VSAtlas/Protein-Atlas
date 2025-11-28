# Atlas Log Resurfacer - Summary

*Fatal:* 0  *Hard:* 2  *Soft:* 2  *Info:* 0

**Health score:** 76 / 100

## 1. [HARD] ligprep.rdkit_valence - 2025
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020128.log`
    - `[02:01:34] Can't kekulize mol.  Unkekulized atoms: 19 20 21 22 23`

## 2. [HARD] ligprep.rdkit_valence - 3ZBF
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020128.log`
    - `[02:01:34] sanitize [02:01:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.sanitized.pdb: Saved fixed PDB to /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.pdb`

## 3. [SOFT] render.enqueue_none - 6WTN
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-08 02:01:46
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6WTN/protein.log`
    - `2025-10-08 02:01:46,191 - INFO - [render-enqueue] 6WTN bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/receptor/6WTN_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6WTN.pdb | ctrl=None | rdk=None`

## 4. [SOFT] ligprep.obabel_h_charge - 3ZBF
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251008_020128.log`
    - `[02:01:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

