# Atlas Log Resurfacer - Summary

*Fatal:* 0  *Hard:* 4  *Soft:* 4  *Info:* 0

**Health score:** 52 / 100

## 1. [HARD] ligprep.rdkit_valence - 2025
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251009_210232.log`
    - `[21:02:48] Explicit valence for atom # 18 C, 5, is greater than permitted`
    - `[21:02:49] Explicit valence for atom # 18 C, 5, is greater than permitted`

## 2. [HARD] ligprep.rdkit_valence - 6WTN
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251009_210232.log`
    - `[21:02:48] sanitize [21:02:48] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.sanitized.pdb: Saved fixed PDB to /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.sanitized.pdb`
    - `[21:02:49] sanitize [21:02:49] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.pdb: /stor/home/mpg2352/micromamba/envs/docking-env/lib/python3.10/site-packages/Bio/PDB/PDBParser.py:384: PDBConstructionWarning: Ignoring unrecognized record 'END' at line 2274`

## 3. [SOFT] ligprep.obabel_h_charge - 6WTN
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251009_210232.log`
    - `[21:02:48] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[21:02:49] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 4. [SOFT] render.enqueue_none - 6WTN
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 21:03:19
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6WTN/protein.log`
    - `2025-10-09 21:03:19,877 - INFO - [render-enqueue] 6WTN bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/receptor/6WTN_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6WTN.pdb | ctrl=None | rdk=None`

## 5. [SOFT] render.enqueue_none - 2025
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-09 21:03:19
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251009_210232.log`
    - `2025-10-09 21:03:19,877 - INFO - [render-enqueue] 6WTN bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/receptor/6WTN_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6WTN.pdb | ctrl=None | rdk=None`

