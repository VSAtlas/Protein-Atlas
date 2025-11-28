# Atlas Log Resurfacer - Summary

*Fatal:* 0  *Hard:* 64  *Soft:* 174  *Info:* 0

**Health score:** 0 / 100

## 1. [HARD] ligprep.rdkit_valence - 
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 63  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:55] Explicit valence for atom # 9 C, 5, is greater than permitted`
    - `[18:33:56] Explicit valence for atom # 19 C, 5, is greater than permitted`
    - `[18:33:56] Can't kekulize mol.  Unkekulized atoms: 19 20 21 22 23`
- **Step:** ligprep.sanitize

## 2. [HARD] ligprep.mgl_partial_write - 6ADQ - CDL_N604.sanitized
- **Message:** Partial write during conversion
- **Count:** 1  |  **Recency:** 
- **Hint:** Conversion partial; capture stderr, retry fallback, then quarantine if persistent.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `WARNING: 72 atoms of 77 in CDL_N604.sanitized.mol2  were not written`
- **Step:** ligprep.write_mol2
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_N604.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_N604.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6ADQ/CDL_N604.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6ADQ/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6ADQ/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 3. [SOFT] ligprep.obabel_h_charge - 3WZD
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 10  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/LEV_A1201.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:58] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1203.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:34:00] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1204.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/EDO_A1207.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/EDO_A1209.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/LEV_A1201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1203.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1204.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/SO4_A1211.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1205.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/GOL_A1210.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/EDO_A1208.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/EDO_A1206.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1202.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1205.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1204.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/LEV_A1201.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1203.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1203.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1205.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/LEV_A1201.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1202.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1204.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/LEV_A1201.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1203.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1204.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1205.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1202.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 4. [SOFT] ligprep.obabel_h_charge - 3CS9
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 8  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_D600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:59] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_D600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_D600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_C600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_D600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_C600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_A600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_B600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_C600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_B600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_D600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_D600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_A600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_C600.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 5. [SOFT] ligprep.obabel_h_charge - 2HYY
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 8  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_D600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:58] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:34:00] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_C600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_D600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_C600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_D600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_C600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_C600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_B600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_D600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_A600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_D600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_C600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_A600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_B600.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 6. [SOFT] ligprep.obabel_h_charge - 5L7I
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 8  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1203.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:58] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_A1202.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:34:00] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_B1202.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_A1202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1204.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1203.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_B1202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/NAG_A1201.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1203.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_B1202.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1204.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_A1202.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/MPG_B1204.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/MPG_B1203.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/VIS_B1202.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/MPG_B1204.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/VIS_A1202.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/VIS_B1202.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/VIS_A1202.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/MPG_B1203.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 7. [SOFT] ligprep.obabel_h_charge - 5I96
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 6  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_A506.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_B503.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:34:00] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/ACT_B505.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/ACT_A508.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_B501.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_B503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A505.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A502.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A504.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_A506.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_A506.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_B503.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_B503.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_A506.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/69Q_B502.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_A506.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/69Q_B502.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_B503.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 8. [SOFT] ligprep.obabel_h_charge - 3ZOS
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 6  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1004.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1004.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:59] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_B1000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/EDO_A1001.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1004.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_B1000.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/EDO_B1001.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/EDO_A1002.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1000.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1004.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_B1000.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1000.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/0LI_A1004.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/0LI_B1000.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/0LI_A1000.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/0LI_A1000.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/0LI_A1004.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/0LI_B1000.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 9. [SOFT] ligprep.obabel_h_charge - 4XUF
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_B1001.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_B1001.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:59] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_A1001.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_B1001.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_A1001.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_B1001.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_A1001.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XUF/P30_B1001.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XUF/P30_A1001.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XUF/P30_B1001.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XUF/P30_A1001.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XUF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 10. [SOFT] ligprep.obabel_h_charge - 6U4J
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:59] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/FLC_B503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/FLC_B503.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_B504.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_A503.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_B504.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_A503.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 11. [SOFT] ligprep.obabel_h_charge - 5MO4
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:59] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/NIL_A601.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/NIL_A601.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/AY7_A602.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/AY7_A602.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 12. [SOFT] ligprep.obabel_h_charge - 4XV2
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:59] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_A801.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_B801.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_A801.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_B801.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 13. [SOFT] ligprep.obabel_h_charge - 2GQG
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:58] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/GOL_B1.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/1N1_B502.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/1N1_A501.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/1N1_A501.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/1N1_B502.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 14. [SOFT] ligprep.obabel_h_charge - 4U5J
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:59] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_B601.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_A601.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_A601.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_B601.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 15. [SOFT] ligprep.obabel_h_charge - 3WZE
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/DTT_A1203.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/BAX_A1201.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:58] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/DTT_A1203.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/ACT_A1202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/DTT_A1203.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/BAX_A1201.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/DTT_A1203.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/BAX_A1201.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/BAX_A1201.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/DTT_A1203.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/BAX_A1201.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/DTT_A1203.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 16. [SOFT] ligprep.obabel_h_charge - 6O0K
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:59] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/2PE_A302.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/2PE_A302.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/LBM_A301.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/LBM_A301.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0K/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 17. [SOFT] ligprep.obabel_h_charge - 6JQR
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:58] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/CXS_A1007.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:59] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/GOL_A1006.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/GOL_A1008.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/CXS_A1007.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/SO4_A1004.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/SO4_A1003.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/SO4_A1002.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/GOL_A1005.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/CXS_A1007.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/CXS_A1007.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/C6F_A1001.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/C6F_A1001.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/CXS_A1007.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 18. [SOFT] ligprep.obabel_h_charge - 6O0L
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:58] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:34:01] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_A301.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/PEG_C303.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/PEG_A303.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/PEG_C302.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/PEG_A302.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_A301.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_A301.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/LBM_A301.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/LBM_A301.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/LBM_C301.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/LBM_C301.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0L/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 19. [SOFT] ligprep.obabel_h_charge - 3LXK
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/MI1_A1125.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/MI1_A1125.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 20. [SOFT] ligprep.obabel_h_charge - 2XP2
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/VGH_A9000.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/VGH_A9000.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 21. [SOFT] ligprep.obabel_h_charge - 6WTN
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/EDO_A1201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/EDO_A1203.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/EDO_A1202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/RXT_A1204.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/RXT_A1204.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 22. [SOFT] ligprep.obabel_h_charge - 4RT7
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4RT7/P30_A1001.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4RT7/P30_A1001.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4RT7/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 23. [SOFT] ligprep.obabel_h_charge - 2WGJ
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/VGH_A2346.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/VGH_A2346.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 24. [SOFT] ligprep.obabel_h_charge - 3ERT
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ERT/OHT_A600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ERT/OHT_A600.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ERT/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 25. [SOFT] ligprep.obabel_h_charge - 1T46
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/PO4_A5.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/PO4_A4.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/1T46/STI_A3.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/1T46/STI_A3.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/1T46/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 26. [SOFT] ligprep.obabel_h_charge - 3G0E
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3G0E/ligands_raw/B49_A9000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3G0E/ligands_raw/B49_A9000.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3G0E/ligands_raw/B49_A9000.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3G0E/ligands_raw/B49_A9000.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3G0E/B49_A9000.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3G0E/B49_A9000.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3G0E/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3G0E/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 27. [SOFT] ligprep.obabel_h_charge - 3ZBF
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/VGH_A3000.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/VGH_A3000.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 28. [SOFT] ligprep.obabel_h_charge - 3OG7
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/032_A1.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/032_A1.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 29. [SOFT] ligprep.obabel_h_charge - 3OXZ
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OXZ/0LI_A1.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OXZ/0LI_A1.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OXZ/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 30. [SOFT] ligprep.obabel_h_charge - 4ASD
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/ligands_raw/BAX_A1500.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/ligands_raw/BAX_A1500.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/ligands_raw/BAX_A1500.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/ligands_raw/BAX_A1500.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4ASD/BAX_A1500.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4ASD/BAX_A1500.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4ASD/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 31. [SOFT] ligprep.obabel_h_charge - 4AG8
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:33:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/ligands_raw/AXI_A2000.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:33:58] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/ligands_raw/AXI_A2000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/ligands_raw/AXI_A2000.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/ligands_raw/AXI_A2000.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4AG8/AXI_A2000.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4AG8/AXI_A2000.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4AG8/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 32. [SOFT] ligprep.obabel_h_charge - 6ADQ
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `[18:34:12] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9YF_M504.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:34:19] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_N604.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9YF_M504.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_N604.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_O304.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_B611.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_A502.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_B608.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9Y0_W201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_F606.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_R605.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_T201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9Y0_B612.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_B610.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/PLM_W203.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9YF_C303.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_B609.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/PLM_K203.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9Y0_P202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_B607.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/PLM_Y301.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_N605.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_H201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9Y0_K201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_P201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9Y0_S301.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9XX_Y302.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9XX_Z302.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_D201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_B605.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9XX_K204.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_N609.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9Y0_G301.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_M503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_N607.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9YF_W202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9Y0_D202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9Y0_B606.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/PLM_Z301.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_N601.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_N606.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9YF_M501.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_R606.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_F605.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_N608.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9YF_A504.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_B603.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_B604.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9YF_A503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_C304.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9YF_K202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9YF_O303.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9XX_W204.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_N610.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_H201.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_F606.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9Y0_B612.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_R605.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_N605.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9YF_C303.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_T201.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/PLM_Y301.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_B608.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9Y0_W201.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9Y0_K201.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_P201.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_B607.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9YF_M504.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_O304.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_B609.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9Y0_P202.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/PLM_W203.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_N604.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_B611.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_B610.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_A502.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/PLM_K203.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6ADQ/CDL_N604.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6ADQ/9YF_M504.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6ADQ/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6ADQ/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 33. [SOFT] render.enqueue_none - 5L2I - 5L2I_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:35
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5L2I/protein.log`
    - `2025-10-11 18:40:35,071 - INFO - [render-enqueue] 5L2I bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L2I_NOLIG/receptor/5L2I_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5L2I.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L2I/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/protein.log`

## 34. [SOFT] render.enqueue_none - 5L2I - 5L2I_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:35
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:40:35,071 - INFO - [render-enqueue] 5L2I bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L2I_NOLIG/receptor/5L2I_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5L2I.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L2I/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 35. [SOFT] render.enqueue_none - 6U4J - 6U4J_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:31
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6U4J/protein.log`
    - `2025-10-11 18:40:31,872 - INFO - [render-enqueue] 6U4J bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J_NOLIG/receptor/6U4J_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6U4J.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/protein.log`

## 36. [SOFT] render.enqueue_none - 6U4J - 6U4J_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:31
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:40:31,872 - INFO - [render-enqueue] 6U4J bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J_NOLIG/receptor/6U4J_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6U4J.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 37. [SOFT] render.enqueue_none - 4XUF - 4XUF_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:13
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4XUF/protein.log`
    - `2025-10-11 18:40:13,473 - INFO - [render-enqueue] 4XUF bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF_NOLIG/receptor/4XUF_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XUF.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XUF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/protein.log`

## 38. [SOFT] render.enqueue_none - 4XUF - 4XUF_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:13
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:40:13,473 - INFO - [render-enqueue] 4XUF bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF_NOLIG/receptor/4XUF_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XUF.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XUF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 39. [SOFT] render.enqueue_none - 2XP2 - 2XP2_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:11
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2XP2/protein.log`
    - `2025-10-11 18:40:11,847 - INFO - [render-enqueue] 2XP2 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2_NOLIG/receptor/2XP2_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2XP2.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/protein.log`

## 40. [SOFT] render.enqueue_none - 2XP2 - 2XP2_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:11
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:40:11,847 - INFO - [render-enqueue] 2XP2 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2_NOLIG/receptor/2XP2_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2XP2.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 41. [SOFT] render.enqueue_none - 2WGJ - 2WGJ_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:10
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2WGJ/protein.log`
    - `2025-10-11 18:40:10,259 - INFO - [render-enqueue] 2WGJ bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ_NOLIG/receptor/2WGJ_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2WGJ.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/protein.log`

## 42. [SOFT] render.enqueue_none - 4U5J - 4U5J_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:10
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4U5J/protein.log`
    - `2025-10-11 18:40:10,989 - INFO - [render-enqueue] 4U5J bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J_NOLIG/receptor/4U5J_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4U5J.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/protein.log`

## 43. [SOFT] render.enqueue_none - 2WGJ - 2WGJ_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:10
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:40:10,259 - INFO - [render-enqueue] 2WGJ bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ_NOLIG/receptor/2WGJ_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2WGJ.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 44. [SOFT] render.enqueue_none - 4U5J - 4U5J_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:10
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:40:10,989 - INFO - [render-enqueue] 4U5J bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J_NOLIG/receptor/4U5J_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4U5J.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 45. [SOFT] render.enqueue_none - 3G0E - 3G0E_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:04
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3G0E/protein.log`
    - `2025-10-11 18:40:04,667 - INFO - [render-enqueue] 3G0E bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3G0E_NOLIG/receptor/3G0E_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3G0E.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3G0E/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3G0E/protein.log`

## 46. [SOFT] render.enqueue_none - 3ZBF - 3ZBF_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:04
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3ZBF/protein.log`
    - `2025-10-11 18:40:04,410 - INFO - [render-enqueue] 3ZBF bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF_NOLIG/receptor/3ZBF_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ZBF.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`

## 47. [SOFT] render.enqueue_none - 3ZBF - 3ZBF_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:04
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:40:04,410 - INFO - [render-enqueue] 3ZBF bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF_NOLIG/receptor/3ZBF_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ZBF.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 48. [SOFT] render.enqueue_none - 3G0E - 3G0E_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:04
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:40:04,667 - INFO - [render-enqueue] 3G0E bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3G0E_NOLIG/receptor/3G0E_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3G0E.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3G0E/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3G0E/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 49. [SOFT] render.enqueue_none - 3OXZ - 3OXZ_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:02
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3OXZ/protein.log`
    - `2025-10-11 18:40:02,554 - INFO - [render-enqueue] 3OXZ bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ_NOLIG/receptor/3OXZ_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3OXZ.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OXZ/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/protein.log`

## 50. [SOFT] render.enqueue_none - 3OXZ - 3OXZ_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:02
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:40:02,554 - INFO - [render-enqueue] 3OXZ bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ_NOLIG/receptor/3OXZ_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3OXZ.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OXZ/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 51. [SOFT] render.enqueue_none - 3OG7 - 3OG7_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:00
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3OG7/protein.log`
    - `2025-10-11 18:40:00,951 - INFO - [render-enqueue] 3OG7 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7_NOLIG/receptor/3OG7_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3OG7.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/protein.log`

## 52. [SOFT] render.enqueue_none - 3ZOS - 3ZOS_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:00
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3ZOS/protein.log`
    - `2025-10-11 18:40:00,337 - INFO - [render-enqueue] 3ZOS bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS_NOLIG/receptor/3ZOS_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ZOS.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/protein.log`

## 53. [SOFT] render.enqueue_none - 4RT7 - 4RT7_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:00
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4RT7/protein.log`
    - `2025-10-11 18:40:00,263 - INFO - [render-enqueue] 4RT7 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7_NOLIG/receptor/4RT7_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4RT7.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4RT7/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/protein.log`

## 54. [SOFT] render.enqueue_none - 4RT7 - 4RT7_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:00
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:40:00,263 - INFO - [render-enqueue] 4RT7 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7_NOLIG/receptor/4RT7_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4RT7.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4RT7/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 55. [SOFT] render.enqueue_none - 3ZOS - 3ZOS_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:00
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:40:00,337 - INFO - [render-enqueue] 3ZOS bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS_NOLIG/receptor/3ZOS_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ZOS.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 56. [SOFT] render.enqueue_none - 3OG7 - 3OG7_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:40:00
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:40:00,951 - INFO - [render-enqueue] 3OG7 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7_NOLIG/receptor/3OG7_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3OG7.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 57. [SOFT] render.enqueue_none - 5I96 - 5I96_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:39:58
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5I96/protein.log`
    - `2025-10-11 18:39:58,087 - INFO - [render-enqueue] 5I96 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96_NOLIG/receptor/5I96_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5I96.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/protein.log`

## 58. [SOFT] render.enqueue_none - 5I96 - 5I96_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:39:58
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:39:58,087 - INFO - [render-enqueue] 5I96 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96_NOLIG/receptor/5I96_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5I96.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 59. [SOFT] render.enqueue_none - 3LXK - 3LXK_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:39:54
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3LXK/protein.log`
    - `2025-10-11 18:39:54,631 - INFO - [render-enqueue] 3LXK bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK_NOLIG/receptor/3LXK_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3LXK.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/protein.log`

## 60. [SOFT] render.enqueue_none - 3LXK - 3LXK_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:39:54
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:39:54,631 - INFO - [render-enqueue] 3LXK bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK_NOLIG/receptor/3LXK_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3LXK.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 61. [SOFT] render.enqueue_none - 6WTN - 6WTN_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:39:53
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6WTN/protein.log`
    - `2025-10-11 18:39:53,821 - INFO - [render-enqueue] 6WTN bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN_NOLIG/receptor/6WTN_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6WTN.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/protein.log`

## 62. [SOFT] render.enqueue_none - 6WTN - 6WTN_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:39:53
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:39:53,821 - INFO - [render-enqueue] 6WTN bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN_NOLIG/receptor/6WTN_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6WTN.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 63. [SOFT] render.enqueue_none - 6JQR - 6JQR_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:39:45
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6JQR/protein.log`
    - `2025-10-11 18:39:45,882 - INFO - [render-enqueue] 6JQR bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR_NOLIG/receptor/6JQR_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6JQR.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/protein.log`

## 64. [SOFT] render.enqueue_none - 6JQR - 6JQR_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:39:45
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:39:45,882 - INFO - [render-enqueue] 6JQR bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR_NOLIG/receptor/6JQR_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6JQR.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 65. [SOFT] render.enqueue_none - 3ERT - 3ERT_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:39:35
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3ERT/protein.log`
    - `2025-10-11 18:39:35,601 - INFO - [render-enqueue] 3ERT bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT_NOLIG/receptor/3ERT_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ERT.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ERT/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/protein.log`

## 66. [SOFT] render.enqueue_none - 3ERT - 3ERT_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:39:35
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:39:35,601 - INFO - [render-enqueue] 3ERT bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT_NOLIG/receptor/3ERT_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ERT.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ERT/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 67. [SOFT] render.enqueue_none - 2HYY - 2HYY_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:38:31
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2HYY/protein.log`
    - `2025-10-11 18:38:31,848 - INFO - [render-enqueue] 2HYY bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY_NOLIG/receptor/2HYY_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2HYY.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/protein.log`

## 68. [SOFT] render.enqueue_none - 2HYY - 2HYY_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:38:31
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:38:31,848 - INFO - [render-enqueue] 2HYY bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY_NOLIG/receptor/2HYY_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2HYY.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 69. [SOFT] render.enqueue_none - 5L7I - 5L7I_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:38:21
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5L7I/protein.log`
    - `2025-10-11 18:38:21,548 - INFO - [render-enqueue] 5L7I bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I_NOLIG/receptor/5L7I_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5L7I.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/protein.log`

## 70. [SOFT] render.enqueue_none - 5L7I - 5L7I_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:38:21
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:38:21,548 - INFO - [render-enqueue] 5L7I bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I_NOLIG/receptor/5L7I_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5L7I.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 71. [SOFT] render.enqueue_none - 3CS9 - 3CS9_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:55
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3CS9/protein.log`
    - `2025-10-11 18:37:55,061 - INFO - [render-enqueue] 3CS9 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9_NOLIG/receptor/3CS9_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3CS9.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/protein.log`

## 72. [SOFT] render.enqueue_none - 3CS9 - 3CS9_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:55
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:37:55,061 - INFO - [render-enqueue] 3CS9 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9_NOLIG/receptor/3CS9_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3CS9.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 73. [SOFT] render.enqueue_none - 1T46 - 1T46_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:49
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-11 18:37:49,809 - INFO - [render-enqueue] 1T46 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46_NOLIG/receptor/1T46_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/1T46.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/1T46/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/protein.log`

## 74. [SOFT] render.enqueue_none - 1T46 - 1T46_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:49
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:37:49,809 - INFO - [render-enqueue] 1T46 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46_NOLIG/receptor/1T46_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/1T46.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/1T46/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 75. [SOFT] render.enqueue_none - 5MO4 - 5MO4_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:48
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5MO4/protein.log`
    - `2025-10-11 18:37:48,637 - INFO - [render-enqueue] 5MO4 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4_NOLIG/receptor/5MO4_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5MO4.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/protein.log`

## 76. [SOFT] render.enqueue_none - 5MO4 - 5MO4_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:48
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:37:48,637 - INFO - [render-enqueue] 5MO4 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4_NOLIG/receptor/5MO4_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5MO4.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 77. [SOFT] render.enqueue_none - 2GQG - 2GQG_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:43
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2GQG/protein.log`
    - `2025-10-11 18:37:43,501 - INFO - [render-enqueue] 2GQG bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG_NOLIG/receptor/2GQG_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2GQG.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/protein.log`

## 78. [SOFT] render.enqueue_none - 3WZD - 3WZD_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:43
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3WZD/protein.log`
    - `2025-10-11 18:37:43,234 - INFO - [render-enqueue] 3WZD bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD_NOLIG/receptor/3WZD_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3WZD.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/protein.log`

## 79. [SOFT] render.enqueue_none - 3WZD - 3WZD_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:43
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:37:43,234 - INFO - [render-enqueue] 3WZD bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD_NOLIG/receptor/3WZD_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3WZD.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 80. [SOFT] render.enqueue_none - 2GQG - 2GQG_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:43
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:37:43,501 - INFO - [render-enqueue] 2GQG bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG_NOLIG/receptor/2GQG_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2GQG.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 81. [SOFT] render.enqueue_none - 4AG8 - 4AG8_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:35
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4AG8/protein.log`
    - `2025-10-11 18:37:35,459 - INFO - [render-enqueue] 4AG8 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8_NOLIG/receptor/4AG8_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4AG8.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4AG8/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/protein.log`

## 82. [SOFT] render.enqueue_none - 4AG8 - 4AG8_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:35
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:37:35,459 - INFO - [render-enqueue] 4AG8 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8_NOLIG/receptor/4AG8_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4AG8.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4AG8/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 83. [SOFT] render.enqueue_none - 4XV2 - 4XV2_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:31
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4XV2/protein.log`
    - `2025-10-11 18:37:31,598 - INFO - [render-enqueue] 4XV2 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2_NOLIG/receptor/4XV2_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XV2.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/protein.log`

## 84. [SOFT] render.enqueue_none - 4XV2 - 4XV2_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:31
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:37:31,598 - INFO - [render-enqueue] 4XV2 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2_NOLIG/receptor/4XV2_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XV2.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 85. [SOFT] render.enqueue_none - 3WZE - 3WZE_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:08
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3WZE/protein.log`
    - `2025-10-11 18:37:08,557 - INFO - [render-enqueue] 3WZE bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE_NOLIG/receptor/3WZE_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3WZE.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/protein.log`

## 86. [SOFT] render.enqueue_none - 4ASD - 4ASD_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:08
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4ASD/protein.log`
    - `2025-10-11 18:37:08,225 - INFO - [render-enqueue] 4ASD bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD_NOLIG/receptor/4ASD_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4ASD.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4ASD/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/protein.log`

## 87. [SOFT] render.enqueue_none - 4ASD - 4ASD_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:08
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:37:08,225 - INFO - [render-enqueue] 4ASD bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD_NOLIG/receptor/4ASD_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4ASD.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4ASD/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 88. [SOFT] render.enqueue_none - 3WZE - 3WZE_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:08
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:37:08,557 - INFO - [render-enqueue] 3WZE bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE_NOLIG/receptor/3WZE_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3WZE.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 89. [SOFT] render.enqueue_none - 6O0L - 6O0L_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:03
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6O0L/protein.log`
    - `2025-10-11 18:37:03,228 - INFO - [render-enqueue] 6O0L bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L_NOLIG/receptor/6O0L_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0L.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0L/protein.log`

## 90. [SOFT] render.enqueue_none - 6O0L - 6O0L_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:03
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:37:03,228 - INFO - [render-enqueue] 6O0L bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L_NOLIG/receptor/6O0L_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0L.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0L/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

## 91. [SOFT] render.enqueue_none - 6O0K - 6O0K_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:00
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6O0K/protein.log`
    - `2025-10-11 18:37:00,388 - INFO - [render-enqueue] 6O0K bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K_NOLIG/receptor/6O0K_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0K.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0K/protein.log`

## 92. [SOFT] render.enqueue_none - 6O0K - 6O0K_NOLIG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-11 18:37:00
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`
    - `2025-10-11 18:37:00,388 - INFO - [render-enqueue] 6O0K bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K_NOLIG/receptor/6O0K_NOLIG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0K.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0K/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251011_183302.log`

