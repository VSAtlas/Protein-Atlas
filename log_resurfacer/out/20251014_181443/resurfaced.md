# Atlas Log Resurfacer - Summary

*Fatal:* 0  *Hard:* 384  *Soft:* 158  *Info:* 0

**Health score:** 0 / 100

## 1. [HARD] ligprep.rdkit_valence - 
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 359  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:31] Can't kekulize mol.  Unkekulized atoms: 19 20 21 22 23`
    - `[18:15:31] Can't kekulize mol.  Unkekulized atoms: 19 20 21 22 23`
    - `[18:15:31] Explicit valence for atom # 14 N, 4, is greater than permitted`
- **Step:** ligprep.sanitize

## 2. [HARD] ligprep.rdkit_valence - 2HYY
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] sanitize [18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.sanitized.pdb: 1 molecule converted`
    - `[18:15:35] sanitize [18:15:35] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.sanitized.pdb: 1 molecule converted`
    - `[18:15:41] sanitize [18:15:41] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_C600.sanitized.pdb: [ligprep] post-ADT PDBQT atoms=42 ok=True reason=`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_C600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_D600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_D600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_C600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_A600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_D600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_C600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_B600.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 3. [HARD] ligprep.rdkit_valence - 4XV2
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:33] sanitize [18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.pdb: ==============================`
    - `[18:15:35] sanitize [18:15:35] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.sanitized.pdb: [18:15:35] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_A801.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_B801.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 4. [HARD] ligprep.rdkit_valence - 6O0L - AY7_A602.sanitized.protoB
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:33] sanitize [18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_A301.sanitized.pdb: WARNING:root:[ligprep] OBabel -h failed (or no H added); attempting RDKit AddHs fallback: AY7_A602.sanitized.protoB.mol2`
    - `[18:15:36] sanitize [18:15:36] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.sanitized.pdb: 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0L/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 5. [HARD] ligprep.rdkit_valence - 2GQG
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:33] sanitize [18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.sanitized.pdb: 1 molecule converted`
    - `[18:15:36] sanitize [18:15:36] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.sanitized.pdb: [18:15:36] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1004.sanitized.pdb: warning - aromatic N with 3 aromatic bonds - skipping charge guess for this atom`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/GOL_B1.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/1N1_A501.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/1N1_B502.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 6. [HARD] ligprep.rdkit_valence - 4U5J
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:33] sanitize [18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.sanitized.pdb: WARNING:root:[ligprep] addHs stderr (primary) 1 molecule converted`
    - `[18:15:37] sanitize [18:15:37] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.sanitized.pdb: [18:15:37] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/DTT_A1203.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_B601.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_A601.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 7. [HARD] ligprep.rdkit_valence - 6U4J
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:34] sanitize [18:15:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.sanitized.pdb: [18:15:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:39] sanitize [18:15:39] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.sanitized.pdb: WARNING:root:[ligprep] addHs stderr (primary) 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/FLC_B503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/FLC_B503.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_A503.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_B504.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 8. [HARD] ligprep.rdkit_valence - 3ZBF
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] sanitize [18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.pdb: WARNING:root:[obabel primary stderr] 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/VGH_A3000.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 9. [HARD] ligprep.rdkit_valence - 2XP2
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] sanitize [18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.sanitized.pdb: [18:15:32] Explicit valence for atom # 0 C, 5, is greater than permitted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/VGH_A9000.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 10. [HARD] ligprep.rdkit_valence - 3LXK
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] sanitize [18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.sanitized.pdb: WARNING:root:[obabel primary stderr] 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/MI1_A1125.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 11. [HARD] ligprep.rdkit_valence - 2WGJ
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] sanitize [18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.pdb: [18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1000.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/VGH_A2346.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 12. [HARD] ligprep.rdkit_valence - 6WTN - BAX_A1201.sanitized.protoB
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:33] sanitize [18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.pdb: WARNING:root:[ligprep] obabel_h_charge: H not added (pre=16 post=16); forcing RDKit AddHs fallback: BAX_A1201.sanitized.protoB.mol2`
- **Step:** ligprep.sanitize
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 13. [HARD] ligprep.rdkit_valence - 1T46
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:33] sanitize [18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.sanitized.pdb: 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/PO4_A5.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/PO4_A4.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/1T46/STI_A3.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/1T46/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 14. [HARD] ligprep.rdkit_valence - 5L2I
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:33] sanitize [18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L2I/ligands_raw/LQQ_A900.sanitized.pdb: 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L2I/ligands_raw/LQQ_A900.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L2I/ligands_raw/LQQ_A900.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L2I/LQQ_A900.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L2I/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 15. [HARD] ligprep.rdkit_valence - 3OG7
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:34] sanitize [18:15:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.pdb: WARNING:root:[ligprep] addHs stderr (primary) 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/032_A1.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 16. [HARD] ligprep.rdkit_valence - 5I96
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:34] sanitize [18:15:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.sanitized.pdb: [18:15:34] Explicit valence for atom # 1 C, 5, is greater than permitted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_A506.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A504.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/ACT_A508.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/ACT_B505.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A502.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_B503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A505.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_B501.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_A506.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_B503.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/69Q_B502.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_A506.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_B503.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 17. [HARD] ligprep.rdkit_valence - 5MO4
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:36] sanitize [18:15:36] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.sanitized.pdb: WARNING:root:[obabel primary stderr] 1 molecule converted`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/NIL_A601.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/AY7_A602.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 18. [HARD] ligprep.rdkit_valence - 6O0K
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:38] sanitize [18:15:38] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.sanitized.pdb: WARNING:root:[arom-mismatch] 0LI_A1004.sanitized.protoB.mol2: MOL2_arom=21 > PDBQT_arom=18 (trying rescue)`
- **Step:** ligprep.sanitize
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/LBM_A301.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/2PE_A302.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0K/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 19. [SOFT] ligprep.obabel_h_charge - 3WZD
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 10  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1202.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1202.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:35] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1203.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1205.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/GOL_A1210.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/EDO_A1208.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/EDO_A1206.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/EDO_A1207.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/EDO_A1209.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/LEV_A1201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1203.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1204.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/SO4_A1211.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1203.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1204.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/LEV_A1201.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1205.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/ligands_raw/DTT_A1202.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/LEV_A1201.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1204.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1205.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1202.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/DTT_A1203.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 20. [SOFT] ligprep.obabel_h_charge - 5L7I
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 8  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1203.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1203.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:36] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1204.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_B1202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/NAG_A1201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_A1202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1203.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1204.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1204.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/MPG_B1203.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_B1202.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/ligands_raw/VIS_A1202.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/MPG_B1204.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/VIS_B1202.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/VIS_A1202.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/MPG_B1203.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 21. [SOFT] ligprep.obabel_h_charge - 2HYY
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 7  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_C600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_D600.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_D600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_C600.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_A600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_D600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_C600.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/STI_B600.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 22. [SOFT] ligprep.obabel_h_charge - 5I96
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 6  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:38] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_A506.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_A506.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A504.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/ACT_A508.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/ACT_B505.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A502.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_B503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_A505.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/GOL_B501.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/69Q_B502.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_A506.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/ligands_raw/NDP_B503.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/69Q_B502.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_A506.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/NDP_B503.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 23. [SOFT] ligprep.obabel_h_charge - 4XV2
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 5  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_A801.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/ligands_raw/P06_B801.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_A801.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/P06_B801.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 24. [SOFT] ligprep.obabel_h_charge - 3ZOS
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 5  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1000.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:36] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1004.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:36] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1004.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1004.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/EDO_A1001.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_B1000.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/EDO_B1001.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1000.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/EDO_A1002.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1000.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1004.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_B1000.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/0LI_B1000.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/0LI_A1000.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/0LI_A1004.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 25. [SOFT] ligprep.obabel_h_charge - 4U5J
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 5  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:36] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_A601.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/ligands_raw/RXT_B601.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_B601.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/RXT_A601.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 26. [SOFT] ligprep.obabel_h_charge - 6U4J
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 5  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:34] sanitize [18:15:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.sanitized.pdb: [18:15:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/FLC_B503.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/FLC_B503.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_B504.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/ligands_raw/PWV_A503.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_A503.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/PWV_B504.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 27. [SOFT] ligprep.obabel_h_charge - 4XUF
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_A1001.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_A1001.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:34] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_B1001.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_A1001.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_B1001.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_A1001.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/ligands_raw/P30_B1001.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XUF/P30_A1001.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XUF/P30_B1001.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XUF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 28. [SOFT] ligprep.obabel_h_charge - 5MO4
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:35] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/AY7_A602.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/ligands_raw/NIL_A601.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/NIL_A601.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/AY7_A602.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 29. [SOFT] ligprep.obabel_h_charge - 2GQG
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:35] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/GOL_B1.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/1N1_A501.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/1N1_B502.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 30. [SOFT] ligprep.obabel_h_charge - 6O0K
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:37] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/2PE_A302.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/ligands_raw/LBM_A301.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/LBM_A301.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/2PE_A302.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0K/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 31. [SOFT] ligprep.obabel_h_charge - 6JQR
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:37] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/CXS_A1007.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/SO4_A1002.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/GOL_A1005.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/SO4_A1004.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/SO4_A1003.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/GOL_A1008.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/GOL_A1006.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/CXS_A1007.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/CXS_A1007.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/ligands_raw/C6F_A1001.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/C6F_A1001.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/CXS_A1007.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 32. [SOFT] ligprep.obabel_h_charge - 2WGJ
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 3  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:32] sanitize [18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.pdb: [18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/ligands_raw/0LI_A1000.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/ligands_raw/VGH_A2346.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/VGH_A2346.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 33. [SOFT] ligprep.obabel_h_charge - 6O0L
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 3  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_A301.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_A301.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:35] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/ligands_raw/LBM_C301.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
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
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/LBM_C301.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/LBM_A301.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0L/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 34. [SOFT] ligprep.obabel_h_charge - 3WZE
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 3  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/BAX_A1201.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/BAX_A1201.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:37] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/DTT_A1203.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/BAX_A1201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/ACT_A1202.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/DTT_A1203.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/DTT_A1203.sanitized.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/ligands_raw/BAX_A1201.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/BAX_A1201.pdbqt`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/DTT_A1203.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 35. [SOFT] ligprep.obabel_h_charge - 3ZBF
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/ligands_raw/VGH_A3000.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/VGH_A3000.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 36. [SOFT] ligprep.obabel_h_charge - 4ASD
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/ligands_raw/BAX_A1500.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/ligands_raw/BAX_A1500.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/ligands_raw/BAX_A1500.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/ligands_raw/BAX_A1500.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4ASD/BAX_A1500.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4ASD/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 37. [SOFT] ligprep.obabel_h_charge - 2XP2
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/ligands_raw/VGH_A9000.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/VGH_A9000.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 38. [SOFT] ligprep.obabel_h_charge - 3LXK
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:31] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/ligands_raw/MI1_A1125.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/MI1_A1125.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 39. [SOFT] ligprep.obabel_h_charge - 4RT7
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/ligands_raw/P30_A1001.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4RT7/P30_A1001.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4RT7/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 40. [SOFT] ligprep.obabel_h_charge - 3OXZ
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/ligands_raw/0LI_A1.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OXZ/0LI_A1.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OXZ/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 41. [SOFT] ligprep.obabel_h_charge - 3ERT
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/ligands_raw/OHT_A600.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ERT/OHT_A600.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ERT/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 42. [SOFT] ligprep.obabel_h_charge - 4AG8
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/ligands_raw/AXI_A2000.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/ligands_raw/AXI_A2000.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/ligands_raw/AXI_A2000.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/ligands_raw/AXI_A2000.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4AG8/AXI_A2000.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4AG8/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 43. [SOFT] ligprep.obabel_h_charge - 6WTN
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `Running Open Babel: [18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/EDO_A1201.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/EDO_A1203.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/EDO_A1202.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/ligands_raw/RXT_A1204.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/RXT_A1204.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 44. [SOFT] ligprep.obabel_h_charge - 1T46
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:32] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/PO4_A5.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/PO4_A4.pdb`
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/1T46/STI_A3.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/1T46/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 45. [SOFT] ligprep.obabel_h_charge - 3OG7
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `[18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[18:15:33] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h
- **Paths:**
  - pre_raw_pdb: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.pdb`
  - pre_raw_pdb_sanitized: `/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/ligands_raw/032_A1.sanitized.pdb`
  - post_pdbqt: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/032_A1.pdbqt`
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/ligand_prep_status.tsv`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 46. [SOFT] render.enqueue_none - 5L7I - 5L7I_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:21:08
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5L7I/protein.log`
    - `2025-10-14 18:21:08,685 - INFO - [render-enqueue] 5L7I bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/receptor/5L7I_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5L7I.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__CONTROL_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/protein.log`

## 47. [SOFT] render.enqueue_none - 5L7I - 5L7I_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:21:08
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:21:08,685 - INFO - [render-enqueue] 5L7I bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L7I/receptor/5L7I_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5L7I.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L7I/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/bench_pocket1_single/5L7I_cleaned__CONTROL_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L7I/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 48. [SOFT] render.enqueue_none - 2HYY - 2HYY_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:20:49
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2HYY/protein.log`
    - `2025-10-14 18:20:49,506 - INFO - [render-enqueue] 2HYY bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/receptor/2HYY_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2HYY.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__CONTROL+RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/protein.log`

## 49. [SOFT] render.enqueue_none - 2HYY - 2HYY_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:20:49
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:20:49,506 - INFO - [render-enqueue] 2HYY bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/receptor/2HYY_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2HYY.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2HYY/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/bench_pocket1_single/2HYY_cleaned__CONTROL+RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2HYY/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 50. [SOFT] render.enqueue_none - 3ZOS - 3ZOS_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:20:42
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3ZOS/protein.log`
    - `2025-10-14 18:20:42,136 - INFO - [render-enqueue] 3ZOS bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/receptor/3ZOS_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ZOS.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__CONTROL+RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/protein.log`

## 51. [SOFT] render.enqueue_none - 5I96 - 5I96_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:20:42
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5I96/protein.log`
    - `2025-10-14 18:20:42,515 - INFO - [render-enqueue] 5I96 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/receptor/5I96_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5I96.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__CONTROL+RDKclosest_front.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/protein.log`

## 52. [SOFT] render.enqueue_none - 3ZOS - 3ZOS_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:20:42
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:20:42,136 - INFO - [render-enqueue] 3ZOS bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZOS/receptor/3ZOS_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ZOS.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZOS/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/bench_pocket1_single/3ZOS_cleaned__CONTROL+RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZOS/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 53. [SOFT] render.enqueue_none - 5I96 - 5I96_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:20:42
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:20:42,515 - INFO - [render-enqueue] 5I96 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5I96/receptor/5I96_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5I96.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5I96/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/bench_pocket1_single/5I96_cleaned__CONTROL+RDKclosest_front.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5I96/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 54. [SOFT] render.enqueue_none - 6U4J - 6U4J_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:20:28
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6U4J/protein.log`
    - `2025-10-14 18:20:28,587 - INFO - [render-enqueue] 6U4J bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/receptor/6U4J_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6U4J.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__CONTROL+RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/protein.log`

## 55. [SOFT] render.enqueue_none - 6U4J - 6U4J_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:20:28
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:20:28,587 - INFO - [render-enqueue] 6U4J bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6U4J/receptor/6U4J_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6U4J.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6U4J/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/bench_pocket1_single/6U4J_cleaned__CONTROL+RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6U4J/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 56. [SOFT] render.enqueue_none - 2GQG - 2GQG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:32
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2GQG/protein.log`
    - `2025-10-14 18:19:32,619 - INFO - [render-enqueue] 2GQG bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/receptor/2GQG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2GQG.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__CONTROL+RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/protein.log`

## 57. [SOFT] render.enqueue_none - 2GQG - 2GQG_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:32
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:32,619 - INFO - [render-enqueue] 2GQG bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/receptor/2GQG_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2GQG.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2GQG/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/bench_pocket1_single/2GQG_cleaned__CONTROL+RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2GQG/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 58. [SOFT] render.enqueue_none - 3OG7 - 3OG7_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:31
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3OG7/protein.log`
    - `2025-10-14 18:19:31,499 - INFO - [render-enqueue] 3OG7 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/receptor/3OG7_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3OG7.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/protein.log`

## 59. [SOFT] render.enqueue_none - 3OG7 - 3OG7_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:31
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:31,499 - INFO - [render-enqueue] 3OG7 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OG7/receptor/3OG7_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3OG7.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OG7/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/bench_pocket1_single/3OG7_cleaned__RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OG7/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 60. [SOFT] render.enqueue_none - 4XUF - 4XUF_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:30
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4XUF/protein.log`
    - `2025-10-14 18:19:30,541 - INFO - [render-enqueue] 4XUF bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/receptor/4XUF_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XUF.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XUF/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__CONTROL+RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/protein.log`

## 61. [SOFT] render.enqueue_none - 4XUF - 4XUF_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:30
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:30,541 - INFO - [render-enqueue] 4XUF bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XUF/receptor/4XUF_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XUF.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XUF/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/bench_pocket1_single/4XUF_cleaned__CONTROL+RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XUF/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 62. [SOFT] render.enqueue_none - 5MO4 - 5MO4_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:26
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5MO4/protein.log`
    - `2025-10-14 18:19:26,505 - INFO - [render-enqueue] 5MO4 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/receptor/5MO4_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5MO4.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/protein.log`

## 63. [SOFT] render.enqueue_none - 5MO4 - 5MO4_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:26
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:26,505 - INFO - [render-enqueue] 5MO4 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5MO4/receptor/5MO4_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5MO4.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5MO4/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/bench_pocket1_single/5MO4_cleaned__RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5MO4/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 64. [SOFT] render.enqueue_none - 4U5J - 4U5J_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:25
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4U5J/protein.log`
    - `2025-10-14 18:19:25,978 - INFO - [render-enqueue] 4U5J bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/receptor/4U5J_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4U5J.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/protein.log`

## 65. [SOFT] render.enqueue_none - 4U5J - 4U5J_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:25
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:25,978 - INFO - [render-enqueue] 4U5J bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4U5J/receptor/4U5J_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4U5J.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4U5J/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/bench_pocket1_single/4U5J_cleaned__RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4U5J/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 66. [SOFT] render.enqueue_none - 4RT7 - 4RT7_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:21
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4RT7/protein.log`
    - `2025-10-14 18:19:21,992 - INFO - [render-enqueue] 4RT7 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/receptor/4RT7_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4RT7.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4RT7/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__CONTROL+RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/protein.log`

## 67. [SOFT] render.enqueue_none - 4RT7 - 4RT7_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:21
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:21,992 - INFO - [render-enqueue] 4RT7 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4RT7/receptor/4RT7_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4RT7.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4RT7/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/bench_pocket1_single/4RT7_cleaned__CONTROL+RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4RT7/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 68. [SOFT] render.enqueue_none - 2XP2 - 2XP2_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:20
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2XP2/protein.log`
    - `2025-10-14 18:19:20,442 - INFO - [render-enqueue] 2XP2 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/receptor/2XP2_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2XP2.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/protein.log`

## 69. [SOFT] render.enqueue_none - 3WZD - 3WZD_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:20
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3WZD/protein.log`
    - `2025-10-14 18:19:20,379 - INFO - [render-enqueue] 3WZD bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/receptor/3WZD_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3WZD.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/protein.log`

## 70. [SOFT] render.enqueue_none - 3WZD - 3WZD_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:20
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:20,379 - INFO - [render-enqueue] 3WZD bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZD/receptor/3WZD_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3WZD.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZD/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/bench_pocket1_single/3WZD_cleaned__RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZD/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 71. [SOFT] render.enqueue_none - 2XP2 - 2XP2_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:20
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:20,442 - INFO - [render-enqueue] 2XP2 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2XP2/receptor/2XP2_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2XP2.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2XP2/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/bench_pocket1_single/2XP2_cleaned__RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2XP2/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 72. [SOFT] render.enqueue_none - 6WTN - 6WTN_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:19
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6WTN/protein.log`
    - `2025-10-14 18:19:19,437 - INFO - [render-enqueue] 6WTN bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/receptor/6WTN_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6WTN.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/protein.log`

## 73. [SOFT] render.enqueue_none - 6WTN - 6WTN_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:19
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:19,437 - INFO - [render-enqueue] 6WTN bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6WTN/receptor/6WTN_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6WTN.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6WTN/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/bench_pocket1_single/6WTN_cleaned__RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6WTN/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 74. [SOFT] render.enqueue_none - 2WGJ - 2WGJ_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:18
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/2WGJ/protein.log`
    - `2025-10-14 18:19:18,267 - INFO - [render-enqueue] 2WGJ bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/receptor/2WGJ_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2WGJ.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/protein.log`

## 75. [SOFT] render.enqueue_none - 3ERT - 3ERT_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:18
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3ERT/protein.log`
    - `2025-10-14 18:19:18,273 - INFO - [render-enqueue] 3ERT bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/receptor/3ERT_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ERT.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ERT/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__CONTROL+RDKclosest_front.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/protein.log`

## 76. [SOFT] render.enqueue_none - 2WGJ - 2WGJ_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:18
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:18,267 - INFO - [render-enqueue] 2WGJ bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2WGJ/receptor/2WGJ_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/2WGJ.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/2WGJ/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/bench_pocket1_single/2WGJ_cleaned__RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/2WGJ/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 77. [SOFT] render.enqueue_none - 3ERT - 3ERT_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:18
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:18,273 - INFO - [render-enqueue] 3ERT bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ERT/receptor/3ERT_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ERT.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ERT/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/bench_pocket1_single/3ERT_cleaned__CONTROL+RDKclosest_front.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ERT/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 78. [SOFT] render.enqueue_none - 3LXK - 3LXK_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:17
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3LXK/protein.log`
    - `2025-10-14 18:19:17,358 - INFO - [render-enqueue] 3LXK bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/receptor/3LXK_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3LXK.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__CONTROL_front.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/protein.log`

## 79. [SOFT] render.enqueue_none - 5L2I - 5L2I_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:17
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/5L2I/protein.log`
    - `2025-10-14 18:19:17,921 - INFO - [render-enqueue] 5L2I bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L2I/receptor/5L2I_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5L2I.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L2I/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/protein.log`

## 80. [SOFT] render.enqueue_none - 3LXK - 3LXK_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:17
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:17,358 - INFO - [render-enqueue] 3LXK bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3LXK/receptor/3LXK_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3LXK.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3LXK/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/bench_pocket1_single/3LXK_cleaned__CONTROL_front.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3LXK/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 81. [SOFT] render.enqueue_none - 5L2I - 5L2I_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:17
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:17,921 - INFO - [render-enqueue] 5L2I bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/5L2I/receptor/5L2I_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/5L2I.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/5L2I/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/bench_pocket1_single/5L2I_cleaned__RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/5L2I/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 82. [SOFT] render.enqueue_none - 4AG8 - 4AG8_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:15
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4AG8/protein.log`
    - `2025-10-14 18:19:15,826 - INFO - [render-enqueue] 4AG8 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/receptor/4AG8_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4AG8.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4AG8/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__CONTROL_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/protein.log`

## 83. [SOFT] render.enqueue_none - 4XV2 - 4XV2_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:15
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4XV2/protein.log`
    - `2025-10-14 18:19:15,760 - INFO - [render-enqueue] 4XV2 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/receptor/4XV2_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XV2.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__CONTROL+RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/protein.log`

## 84. [SOFT] render.enqueue_none - 4XV2 - 4XV2_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:15
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:15,760 - INFO - [render-enqueue] 4XV2 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4XV2/receptor/4XV2_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4XV2.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4XV2/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/bench_pocket1_single/4XV2_cleaned__CONTROL+RDKclosest_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4XV2/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 85. [SOFT] render.enqueue_none - 4AG8 - 4AG8_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:15
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:15,826 - INFO - [render-enqueue] 4AG8 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4AG8/receptor/4AG8_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4AG8.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4AG8/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/bench_pocket1_single/4AG8_cleaned__CONTROL_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4AG8/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 86. [SOFT] render.enqueue_none - 1T46 - 1T46_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:14
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/1T46/protein.log`
    - `2025-10-14 18:19:14,957 - INFO - [render-enqueue] 1T46 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/receptor/1T46_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/1T46.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/1T46/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__CONTROL_front.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/protein.log`

## 87. [SOFT] render.enqueue_none - 3ZBF - 3ZBF_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:14
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3ZBF/protein.log`
    - `2025-10-14 18:19:14,428 - INFO - [render-enqueue] 3ZBF bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/receptor/3ZBF_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ZBF.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`

## 88. [SOFT] render.enqueue_none - 6JQR - 6JQR_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:14
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6JQR/protein.log`
    - `2025-10-14 18:19:14,795 - INFO - [render-enqueue] 6JQR bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/receptor/6JQR_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6JQR.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__CONTROL_front.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/protein.log`

## 89. [SOFT] render.enqueue_none - 3ZBF - 3ZBF_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:14
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:14,428 - INFO - [render-enqueue] 3ZBF bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3ZBF/receptor/3ZBF_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3ZBF.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3ZBF/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/bench_pocket1_single/3ZBF_cleaned__RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3ZBF/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 90. [SOFT] render.enqueue_none - 6JQR - 6JQR_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:14
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:14,795 - INFO - [render-enqueue] 6JQR bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6JQR/receptor/6JQR_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6JQR.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6JQR/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/bench_pocket1_single/6JQR_cleaned__CONTROL_front.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6JQR/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 91. [SOFT] render.enqueue_none - 1T46 - 1T46_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:14
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:14,957 - INFO - [render-enqueue] 1T46 bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/receptor/1T46_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/1T46.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/1T46/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/bench_pocket1_single/1T46_cleaned__CONTROL_front.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/1T46/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 92. [SOFT] render.enqueue_none - 3OXZ - 3OXZ_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:13
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3OXZ/protein.log`
    - `2025-10-14 18:19:13,640 - INFO - [render-enqueue] 3OXZ bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/receptor/3OXZ_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3OXZ.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OXZ/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__CONTROL+RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/protein.log`

## 93. [SOFT] render.enqueue_none - 3OXZ - 3OXZ_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:13
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:13,640 - INFO - [render-enqueue] 3OXZ bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3OXZ/receptor/3OXZ_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3OXZ.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3OXZ/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/bench_pocket1_single/3OXZ_cleaned__CONTROL+RDKclosest_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3OXZ/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 94. [SOFT] render.enqueue_none - 3WZE - 3WZE_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:12
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/3WZE/protein.log`
    - `2025-10-14 18:19:12,343 - INFO - [render-enqueue] 3WZE bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/receptor/3WZE_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3WZE.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__CONTROL_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/protein.log`

## 95. [SOFT] render.enqueue_none - 3WZE - 3WZE_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:12
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:12,343 - INFO - [render-enqueue] 3WZE bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/3WZE/receptor/3WZE_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/3WZE.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/3WZE/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__CONTROL+RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/bench_pocket1_single/3WZE_cleaned__CONTROL_side.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/3WZE/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 96. [SOFT] render.enqueue_none - 4ASD - 4ASD_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:10
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/4ASD/protein.log`
    - `2025-10-14 18:19:10,036 - INFO - [render-enqueue] 4ASD bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/receptor/4ASD_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4ASD.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4ASD/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__CONTROL+RDKclosest_front.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/protein.log`

## 97. [SOFT] render.enqueue_none - 4ASD - 4ASD_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:19:10
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:19:10,036 - INFO - [render-enqueue] 4ASD bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/4ASD/receptor/4ASD_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/4ASD.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/4ASD/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__RDKclosest_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__CONTROL_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__CONTROL+RDKclosest_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__CONTROL+RDKclosest_top.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/bench_pocket1_single/4ASD_cleaned__CONTROL+RDKclosest_front.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/4ASD/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 98. [SOFT] render.enqueue_none - 6O0K - 6O0K_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:18:57
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6O0K/protein.log`
    - `2025-10-14 18:18:57,982 - INFO - [render-enqueue] 6O0K bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/receptor/6O0K_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0K.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0K/bench_pocket1_single/6O0K_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0K/bench_pocket1_single/6O0K_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0K/bench_pocket1_single/6O0K_cleaned__CONTROL_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0K/protein.log`

## 99. [SOFT] render.enqueue_none - 6O0L - 6O0L_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:18:57
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `docked/6O0L/protein.log`
    - `2025-10-14 18:18:57,123 - INFO - [render-enqueue] 6O0L bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/receptor/6O0L_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0L.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0L/bench_pocket1_single/6O0L_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0L/bench_pocket1_single/6O0L_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0L/bench_pocket1_single/6O0L_cleaned__CONTROL_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0L/protein.log`

## 100. [SOFT] render.enqueue_none - 6O0L - 6O0L_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:18:57
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:18:57,123 - INFO - [render-enqueue] 6O0L bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0L/receptor/6O0L_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0L.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0L/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0L/bench_pocket1_single/6O0L_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0L/bench_pocket1_single/6O0L_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0L/bench_pocket1_single/6O0L_cleaned__CONTROL_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0L/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

## 101. [SOFT] render.enqueue_none - 6O0K - 6O0K_cleaned
- **Message:** Controls/RDK not enqueued
- **Count:** 1  |  **Recency:** 2025-10-14 18:18:57
- **Hint:** Control selection failed or disabled; check 'hint_count' and control whitelist.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`
    - `2025-10-14 18:18:57,982 - INFO - [render-enqueue] 6O0K bench_pocket1_single: receptor=/stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6O0K/receptor/6O0K_cleaned.pdb | orig=/stor/home/mpg2352/atlas/code/protein_automation/input_pdbs/6O0K.pdb | ctrl=None | rdk=None`
- **Step:** render
- **Paths:**
  - ligand_prep_status_tsv: `/stor/home/mpg2352/atlas/code/protein_automation/prepped_ligands/6O0K/ligand_prep_status.tsv`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0K/bench_pocket1_single/6O0K_cleaned__CONTROL_front.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0K/bench_pocket1_single/6O0K_cleaned__CONTROL_side.png`
  - dock_inputs: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0K/bench_pocket1_single/6O0K_cleaned__CONTROL_top.png`
  - protein_log: `/stor/home/mpg2352/atlas/code/protein_automation/docked/6O0K/protein.log`
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251014_181443.log`

