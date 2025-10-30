# Atlas Log Resurfacer - Summary

*Fatal:* 0  *Hard:* 44  *Soft:* 10  *Info:* 0

**Health score:** 0 / 100

## 1. [HARD] ligprep.rdkit_valence - 
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 44  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251013_210841.log`
    - `[21:08:46] Explicit valence for atom # 27 C, 5, is greater than permitted`
    - `[21:08:46] Explicit valence for atom # 27 C, 5, is greater than permitted`
    - `[21:08:46] Explicit valence for atom # 27 C, 5, is greater than permitted`
- **Step:** ligprep.sanitize

## 2. [SOFT] ligprep.obabel_h_charge - 3CS9
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 8  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251013_210841.log`
    - `[21:08:46] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[21:08:46] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[21:08:48] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_D600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_C600.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_D600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_C600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_A600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_B600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_C600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_D600.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251013_210841.log`

## 3. [SOFT] render.enqueue_none - 3CS9 - 3CS9_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-13 21:09:09
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3CS9/protein.log`
    - `2025-10-13 21:09:09,423 - INFO - [render-enqueue] 3CS9 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/receptor/3CS9_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3CS9.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__CONTROL_front.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/protein.log`

## 4. [SOFT] render.enqueue_none - 3CS9 - 3CS9_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-13 21:09:09
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251013_210841.log`
    - `2025-10-13 21:09:09,423 - INFO - [render-enqueue] 3CS9 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/receptor/3CS9_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3CS9.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/bench_pocket1_single/3CS9_cleaned__CONTROL_front.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251013_210841.log`

