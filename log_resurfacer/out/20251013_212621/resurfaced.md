# Atlas Log Resurfacer - Summary

*Fatal:* 0  *Hard:* 76  *Soft:* 8  *Info:* 0

**Health score:** 0 / 100

## 1. [HARD] ligprep.rdkit_valence - 
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 72  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251013_212621.log`
    - `[21:26:27] Explicit valence for atom # 27 C, 5, is greater than permitted`
    - `[21:26:27] Explicit valence for atom # 27 C, 5, is greater than permitted`
    - `[21:26:27] Explicit valence for atom # 27 C, 5, is greater than permitted`
- **Step:** ligprep.sanitize

## 2. [HARD] ligprep.rdkit_valence - 3CS9
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251013_212621.log`
    - `[21:26:27] sanitize [21:26:27] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.sanitized.sanitized.pdb: WARNING:root:[ligprep] ADT wrote fewer atoms (63 -> 45); invoking OBabel MOL2->PDBQT fallback`
    - `[21:26:28] sanitize [21:26:28] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.sanitized.sanitized.pdb: WARNING:root:[ligprep] ADT wrote fewer atoms (63 -> 45); invoking OBabel MOL2->PDBQT fallback`
    - `[21:26:29] sanitize [21:26:29] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_C600.sanitized.sanitized.pdb: [ligprep-extracted] discovery_root=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw exists=True is_dir=True`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_D600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_C600.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_D600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_C600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_C600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_D600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_A600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_A600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_B600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_C600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_B600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_D600.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251013_212621.log`

## 3. [SOFT] ligprep.obabel_h_charge - 3CS9
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 8  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251013_212621.log`
    - `[21:26:27] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[21:26:27] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[21:26:28] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
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
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_C600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_D600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_A600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_A600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_B600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_C600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_B600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_D600.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251013_212621.log`

