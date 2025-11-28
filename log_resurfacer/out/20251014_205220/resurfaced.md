# Atlas Log Resurfacer - Summary

*Fatal:* 0  *Hard:* 625  *Soft:* 110  *Info:* 0

**Health score:** 0 / 100

## 1. [HARD] ligprep.rdkit_valence - 
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 577  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:54] Explicit valence for atom # 19 C, 5, is greater than permitted`
    - `[20:52:54] Explicit valence for atom # 14 C, 5, is greater than permitted`
    - `[20:52:54] Explicit valence for atom # 14 C, 5, is greater than permitted`
- **Step:** ligprep.sanitize

## 2. [HARD] ligprep.rdkit_valence - 2HYY
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:55] sanitize [20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.sanitized.sanitized.pdb: 1 molecule converted`
    - `[20:52:58] sanitize [20:52:58] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.sanitized.sanitized.pdb: [ligprep] writer=mgltools lig=P30_A1001.sanitized.sanitized.protoB in_atoms=74 out_atoms=76 ok=True reason=`
    - `[20:53:03] sanitize [20:53:03] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_C600.sanitized.sanitized.pdb: 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_D600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_C600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_C600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_D600.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_A600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_B600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_D600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_A600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_C600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_C600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_B600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_D600.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 3. [HARD] ligprep.rdkit_valence - 3CS9
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:55] sanitize [20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.sanitized.sanitized.pdb: 1 molecule converted`
    - `[20:52:59] sanitize [20:52:59] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.sanitized.sanitized.pdb: [20:52:59] Explicit valence for atom # 16 C, 5, is greater than permitted`
    - `[20:53:06] sanitize [20:53:06] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_C600.sanitized.sanitized.pdb: [20:53:06] Explicit valence for atom # 1 C, 5, is greater than permitted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_D600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_C600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_C600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_D600.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_B600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_A600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_C600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_A600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_D600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_B600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_D600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_C600.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 4. [HARD] ligprep.rdkit_valence - 5I96
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 3  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] sanitize [20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.sanitized.sanitized.pdb: [20:52:56] Explicit valence for atom # 18 C, 5, is greater than permitted`
    - `[20:53:00] sanitize [20:53:00] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_A506.sanitized.sanitized.pdb: 1 molecule converted`
    - `[20:53:06] sanitize [20:53:06] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_B503.sanitized.sanitized.pdb: [20:53:06] Explicit valence for atom # 27 C, 5, is greater than permitted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_A506.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A504.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A502.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A505.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_B503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_B501.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/ACT_A508.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/ACT_B505.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_B503.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_A506.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/69Q_B502.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/69Q_B502.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_B503.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_B503.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_A506.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_A506.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 5. [HARD] ligprep.rdkit_valence - 3ZOS - STI_D600.sanitized.sanitized.arom
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 3  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:57] sanitize [20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1000.sanitized.sanitized.pdb: 1 molecule converted`
    - `[20:53:03] sanitize [20:53:03] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1004.sanitized.sanitized.pdb: WARNING:root:[ligprep] addHs stderr (primary) 1 molecule converted`
    - `[20:53:09] sanitize [20:53:09] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_B1000.sanitized.sanitized.pdb: Running Open Babel: /stor/home/mpg2352/micromamba/envs/docking-env/bin/obabel -isdf /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_D600.sanitized.sanitized.arom.sdf -omol2 -O /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_D600.sanitized.sanitized.arom.mol2`
- **Step:** ligprep.sanitize
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 6. [HARD] ligprep.rdkit_valence - 6O0L
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] sanitize [20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_A301.sanitized.sanitized.pdb: [20:52:56] Explicit valence for atom # 27 N, 4, is greater than permitted`
    - `[20:53:00] sanitize [20:53:00] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.sanitized.sanitized.pdb: 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_A301.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/PEG_A302.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/PEG_C302.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/PEG_A303.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/PEG_C303.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_A301.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/LBM_C301.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/LBM_A301.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/LBM_C301.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/LBM_A301.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0L/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 7. [HARD] ligprep.rdkit_valence - 4XV2
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] sanitize [20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.sanitized.pdb: 1 molecule converted`
    - `[20:53:00] sanitize [20:53:00] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.sanitized.sanitized.pdb: [20:53:01] Explicit valence for atom # 11 C, 5, is greater than permitted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_A801.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_A801.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_B801.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_B801.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 8. [HARD] ligprep.rdkit_valence - 4U5J
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] sanitize [20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.sanitized.sanitized.pdb: 1 molecule converted`
    - `[20:53:00] sanitize [20:53:00] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.sanitized.sanitized.pdb: WARNING:root:[obabel primary stderr] 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_B601.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_A601.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_A601.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_B601.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 9. [HARD] ligprep.rdkit_valence - 2GQG - BAX_A1201.sanitized.sanitized.protoB.arom
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:57] sanitize [20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.sanitized.sanitized.pdb: /stor/home/mpg2352/micromamba/envs/docking-env/bin/obabel -isdf /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/BAX_A1201.sanitized.sanitized.protoB.arom.sdf -omol2 -O /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/BAX_A1201.sanitized.sanitized.protoB.arom.mol2`
    - `[20:53:01] sanitize [20:53:01] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.sanitized.sanitized.pdb: 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 10. [HARD] ligprep.rdkit_valence - 4XUF
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:57] sanitize [20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_A1001.sanitized.sanitized.pdb: 1 molecule converted`
    - `[20:53:02] sanitize [20:53:02] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_B1001.sanitized.sanitized.pdb: 1 molecule converted`
- **Step:** ligprep.sanitize
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
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 11. [HARD] ligprep.rdkit_valence - 6JQR
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:57] sanitize [20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.sanitized.sanitized.pdb: [20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:53:01] sanitize [20:53:01] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/CXS_A1007.sanitized.sanitized.pdb: 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/SO4_A1004.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/SO4_A1003.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/GOL_A1005.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/SO4_A1002.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/CXS_A1007.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/GOL_A1006.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/GOL_A1008.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/CXS_A1007.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/C6F_A1001.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/CXS_A1007.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/CXS_A1007.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/C6F_A1001.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 12. [HARD] ligprep.rdkit_valence - 6U4J
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:58] sanitize [20:52:58] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.sanitized.sanitized.pdb: 1 molecule converted`
    - `[20:53:02] sanitize [20:53:02] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.sanitized.sanitized.pdb: 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/FLC_B503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/FLC_B503.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_A503.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_B504.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_B504.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_A503.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 13. [HARD] dock.no_score - 3ZBF - rdk_0002898.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-14 21:23:02
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/3ZBF/protein.log`
    - `2025-10-14 21:23:02,479 - WARNING - No score for rdk_0002898.pdbqt`
- **Step:** dock.score
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`

## 14. [HARD] dock.no_score - 3ZBF - rdk_0003480.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-14 21:21:38
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/3ZBF/protein.log`
    - `2025-10-14 21:21:38,297 - WARNING - No score for rdk_0003480.pdbqt`
- **Step:** dock.score
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`

## 15. [HARD] dock.no_score - 3ZBF - rdk_0003045.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-14 21:20:06
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/3ZBF/protein.log`
    - `2025-10-14 21:20:06,487 - WARNING - No score for rdk_0003045.pdbqt`
- **Step:** dock.score
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`

## 16. [HARD] dock.no_score - 3ZBF - rdk_0002965.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-14 21:14:53
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/3ZBF/protein.log`
    - `2025-10-14 21:14:53,963 - WARNING - No score for rdk_0002965.pdbqt`
- **Step:** dock.score
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`

## 17. [HARD] dock.no_score - 3ZBF - rdk_0003776.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-14 21:09:20
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/3ZBF/protein.log`
    - `2025-10-14 21:09:20,568 - WARNING - No score for rdk_0003776.pdbqt`
- **Step:** dock.score
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`

## 18. [HARD] dock.no_score - 3ZBF - rdk_0004318.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-14 21:08:28
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/3ZBF/protein.log`
    - `2025-10-14 21:08:28,207 - WARNING - No score for rdk_0004318.pdbqt`
- **Step:** dock.score
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`

## 19. [HARD] dock.no_score - 3ZBF - rdk_0002968.pdbqt
- **Message:** No score for 1 ligand(s)
- **Count:** 1  |  **Recency:** 2025-10-14 21:08:28
- **Hint:** Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.
  - **File:** `docked/3ZBF/protein.log`
    - `2025-10-14 21:08:28,207 - WARNING - No score for rdk_0002968.pdbqt`
- **Step:** dock.score
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`

## 20. [HARD] ligprep.rdkit_valence - 3LXK
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:55] sanitize [20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.sanitized.sanitized.pdb: 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/MI1_A1125.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/MI1_A1125.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 21. [HARD] ligprep.rdkit_valence - 3ZBF
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:55] sanitize [20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.sanitized.pdb: [20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/VGH_A3000.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/VGH_A3000.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 22. [HARD] ligprep.rdkit_valence - 3OG7
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:55] sanitize [20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.sanitized.pdb: WARNING:root:[ligprep] addHs stderr (primary) 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/032_A1.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/032_A1.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 23. [HARD] ligprep.rdkit_valence - 6WTN
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] sanitize [20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.sanitized.pdb: WARNING:root:[obabel primary stderr] 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/EDO_A1201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/EDO_A1203.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/EDO_A1202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/RXT_A1204.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/RXT_A1204.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 24. [HARD] ligprep.rdkit_valence - 3ERT
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] sanitize [20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.sanitized.sanitized.pdb: [20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ERT/OHT_A600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ERT/OHT_A600.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ERT/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 25. [HARD] ligprep.rdkit_valence - 3OXZ
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] sanitize [20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.sanitized.sanitized.pdb: 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OXZ/0LI_A1.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OXZ/0LI_A1.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OXZ/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 26. [HARD] ligprep.rdkit_valence - 5L2I
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] sanitize [20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L2I/ligands_raw/LQQ_A900.sanitized.sanitized.pdb: 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L2I/ligands_raw/LQQ_A900.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L2I/ligands_raw/LQQ_A900.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L2I/LQQ_A900.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L2I/LQQ_A900.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L2I/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 27. [HARD] ligprep.rdkit_valence - 2WGJ
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] sanitize [20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.sanitized.pdb: 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/VGH_A2346.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/VGH_A2346.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 28. [HARD] ligprep.rdkit_valence - 1T46 - 0LI_A1.sanitized.sanitized.protoB
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] sanitize [20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.sanitized.sanitized.pdb: [ligprep] Calling prepare_ligand4 on 0LI_A1.sanitized.sanitized.protoB.mol2 -> 0LI_A1.sanitized.pdbqt`
- **Step:** ligprep.sanitize
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/1T46/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 29. [HARD] ligprep.rdkit_valence - 2XP2
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:57] sanitize [20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.sanitized.sanitized.pdb: 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/VGH_A9000.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/VGH_A9000.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 30. [HARD] ligprep.rdkit_valence - 4RT7
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:57] sanitize [20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.sanitized.sanitized.pdb: [ligprep] writer=mgltools lig=AXI_A2000.sanitized.sanitized.protoB in_atoms=46 out_atoms=49 ok=True reason=`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4RT7/P30_A1001.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4RT7/P30_A1001.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4RT7/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 31. [HARD] ligprep.rdkit_valence - 6O0K
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:53:00] sanitize [20:53:00] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.sanitized.sanitized.pdb: WARNING:root:[ligprep] addHs stderr (primary) 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/2PE_A302.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/LBM_A301.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/LBM_A301.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/2PE_A302.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0K/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 32. [HARD] ligprep.rdkit_valence - 5MO4
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:53:02] sanitize [20:53:02] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.sanitized.sanitized.pdb: 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/NIL_A601.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/AY7_A602.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/AY7_A602.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/NIL_A601.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 33. [SOFT] ligprep.obabel_h_charge - 3WZD
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 10  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1202.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1202.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:53:01] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1203.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/EDO_A1208.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/EDO_A1206.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/GOL_A1210.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1205.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/SO4_A1211.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/LEV_A1201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1203.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1204.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/EDO_A1207.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/EDO_A1209.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1205.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1202.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1203.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/LEV_A1201.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1204.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1204.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1205.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1203.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1202.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/LEV_A1201.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1205.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/LEV_A1201.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1204.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1203.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1202.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 34. [SOFT] ligprep.obabel_h_charge - 2HYY
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 8  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_D600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_C600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_C600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_D600.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_A600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_B600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_D600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_A600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_C600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_C600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_B600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_D600.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 35. [SOFT] ligprep.obabel_h_charge - 3CS9
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 8  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:58] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_D600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_C600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_A600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_C600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_B600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3CS9/ligands_raw/NIL_D600.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_B600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_A600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_C600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_A600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_D600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_B600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_D600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/NIL_C600.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3CS9/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3CS9/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 36. [SOFT] ligprep.obabel_h_charge - 5L7I
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 8  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1203.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1203.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:53:02] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1204.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_A1202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1203.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1204.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_B1202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/NAG_A1201.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_A1202.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1204.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_B1202.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1203.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/MPG_B1203.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/VIS_A1202.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/VIS_A1202.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/VIS_B1202.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/MPG_B1204.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/VIS_B1202.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/MPG_B1203.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/MPG_B1204.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 37. [SOFT] ligprep.obabel_h_charge - 5I96
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 6  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:58] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_A506.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_A506.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A504.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A502.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A505.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_B503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_B501.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/ACT_A508.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/ACT_B505.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_B503.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_A506.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/69Q_B502.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/69Q_B502.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_B503.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_B503.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_A506.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_A506.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 38. [SOFT] ligprep.obabel_h_charge - 6JQR
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 5  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:57] sanitize [20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.sanitized.sanitized.pdb: [20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/SO4_A1004.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/SO4_A1003.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/GOL_A1005.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/SO4_A1002.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/CXS_A1007.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/GOL_A1006.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/GOL_A1008.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/CXS_A1007.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/C6F_A1001.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/CXS_A1007.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/CXS_A1007.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/C6F_A1001.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 39. [SOFT] ligprep.obabel_h_charge - 3ZOS
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 5  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:53:02] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1004.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:53:02] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1004.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_B1000.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1004.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/EDO_A1001.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1000.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/EDO_A1002.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/EDO_B1001.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1004.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_B1000.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1000.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/0LI_A1000.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/0LI_A1004.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/0LI_B1000.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/0LI_A1000.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/0LI_A1004.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/0LI_B1000.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 40. [SOFT] ligprep.obabel_h_charge - 6O0L
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_A301.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_A301.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:59] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_A301.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/PEG_A302.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/PEG_C302.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/PEG_A303.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/PEG_C303.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_A301.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/LBM_C301.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/LBM_A301.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/LBM_C301.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/LBM_A301.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0L/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 41. [SOFT] ligprep.obabel_h_charge - 6O0K
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:59] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/2PE_A302.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/LBM_A301.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/LBM_A301.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/2PE_A302.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0K/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 42. [SOFT] ligprep.obabel_h_charge - 4XV2
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:53:00] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_A801.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_A801.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_B801.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_B801.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 43. [SOFT] ligprep.obabel_h_charge - 5MO4
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:53:01] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/NIL_A601.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/AY7_A602.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/AY7_A602.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/NIL_A601.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 44. [SOFT] ligprep.obabel_h_charge - 3WZE
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/BAX_A1201.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/BAX_A1201.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:53:00] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/DTT_A1203.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/BAX_A1201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/DTT_A1203.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/ACT_A1202.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/BAX_A1201.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/DTT_A1203.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/BAX_A1201.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/BAX_A1201.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/DTT_A1203.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/DTT_A1203.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 45. [SOFT] ligprep.obabel_h_charge - 2GQG
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:53:00] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/GOL_B1.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/1N1_A501.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/1N1_A501.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/1N1_B502.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/1N1_B502.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 46. [SOFT] ligprep.obabel_h_charge - 4XUF
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_A1001.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_A1001.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:53:02] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_B1001.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
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
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 47. [SOFT] ligprep.obabel_h_charge - 6U4J
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:53:01] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/FLC_B503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/FLC_B503.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_A503.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_B504.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_B504.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_A503.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 48. [SOFT] ligprep.obabel_h_charge - 3ZBF
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 3  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:55] sanitize [20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.sanitized.pdb: [20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/VGH_A3000.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/VGH_A3000.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 49. [SOFT] ligprep.obabel_h_charge - 3ERT
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 3  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:56] sanitize [20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.sanitized.sanitized.pdb: [20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ERT/OHT_A600.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ERT/OHT_A600.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ERT/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 50. [SOFT] ligprep.obabel_h_charge - 4U5J
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 3  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:59] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:59] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_B601.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_A601.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_A601.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_B601.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 51. [SOFT] ligprep.obabel_h_charge - 3LXK
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:54] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:54] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/MI1_A1125.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/MI1_A1125.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 52. [SOFT] ligprep.obabel_h_charge - 3OG7
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/032_A1.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/032_A1.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 53. [SOFT] ligprep.obabel_h_charge - 6WTN
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/EDO_A1201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/EDO_A1203.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/EDO_A1202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/RXT_A1204.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/RXT_A1204.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 54. [SOFT] ligprep.obabel_h_charge - 4AG8
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/ligands_raw/AXI_A2000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/ligands_raw/AXI_A2000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/ligands_raw/AXI_A2000.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/ligands_raw/AXI_A2000.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4AG8/AXI_A2000.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4AG8/AXI_A2000.sanitized.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4AG8/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 55. [SOFT] ligprep.obabel_h_charge - 3OXZ
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:55] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OXZ/0LI_A1.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OXZ/0LI_A1.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OXZ/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 56. [SOFT] ligprep.obabel_h_charge - 4ASD
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/ligands_raw/BAX_A1500.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/ligands_raw/BAX_A1500.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/ligands_raw/BAX_A1500.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/ligands_raw/BAX_A1500.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4ASD/BAX_A1500.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4ASD/BAX_A1500.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4ASD/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 57. [SOFT] ligprep.obabel_h_charge - 1T46
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/PO4_A4.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/PO4_A5.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/1T46/STI_A3.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/1T46/STI_A3.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/1T46/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 58. [SOFT] ligprep.obabel_h_charge - 2XP2
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/VGH_A9000.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/VGH_A9000.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 59. [SOFT] ligprep.obabel_h_charge - 4RT7
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[20:52:57] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4RT7/P30_A1001.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4RT7/P30_A1001.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4RT7/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

## 60. [SOFT] ligprep.obabel_h_charge - 2WGJ
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`
    - `[20:52:56] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/VGH_A2346.sanitized.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/VGH_A2346.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_205220.log`

