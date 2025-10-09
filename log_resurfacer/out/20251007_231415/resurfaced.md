# Atlas Log Resurfacer - Summary

*Fatal:* 0  *Hard:* 14  *Soft:* 7  *Info:* 0

**Health score:** 0 / 100

## 1. [HARD] ligprep.rdkit_valence - 2025
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 7  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_231415.log`
    - `[23:14:21] Explicit valence for atom # 10 C, 5, is greater than permitted`
    - `[23:14:21] Explicit valence for atom # 11 C, 5, is greater than permitted`
    - `[23:14:21] Explicit valence for atom # 10 C, 5, is greater than permitted`

## 2. [HARD] ligprep.rdkit_valence - 2HYY
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_231415.log`
    - `[23:14:21] sanitize [23:14:21] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_D600.sanitized.sanitized.pdb: Saved fixed PDB to /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_A600.sanitized.pdb`
    - `[23:14:21] sanitize [23:14:21] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.sanitized.sanitized.pdb: WARNING:root:MODELLER failed on /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46_NOLIG/work/1T46_NOLIG_elemfix.pdb: pdbnam_____E> Filename for PDB code not found: /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46_NOLIG/work/1T46_NOLIG_elemfix.pdb`
    - `[23:14:21] sanitize [23:14:21] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_C600.sanitized.sanitized.pdb: WARNING:root:MODELLER failed on /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG_NOLIG/work/2GQG_NOLIG_elemfix.pdb: pdbnam_____E> Filename for PDB code not found: /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG_NOLIG/work/2GQG_NOLIG_elemfix.pdb`

## 3. [HARD] ligprep.rdkit_valence - 2GQG
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_231415.log`
    - `[23:14:21] sanitize [23:14:21] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.sanitized.sanitized.pdb: WARNING:root:[obabel primary stderr] 1 molecule converted`
    - `[23:14:21] sanitize [23:14:21] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.sanitized.sanitized.pdb: WARNING:root:[obabel primary stderr] 1 molecule converted`

## 4. [HARD] ligprep.rdkit_valence - 1T46
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_231415.log`
    - `[23:14:21] sanitize [23:14:21] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.sanitized.sanitized.pdb: WARNING:root:[obabel primary stderr] 1 molecule converted`

## 5. [SOFT] ligprep.obabel_h_charge - 2HYY
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 4  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_231415.log`
    - `[23:14:21] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_D600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[23:14:21] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_B600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[23:14:21] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2HYY/ligands_raw/STI_C600.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 6. [SOFT] ligprep.obabel_h_charge - 2GQG
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_231415.log`
    - `[23:14:21] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_B502.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[23:14:21] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/2GQG/ligands_raw/1N1_A501.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

## 7. [SOFT] ligprep.obabel_h_charge - 1T46
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 1  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251007_231415.log`
    - `[23:14:21] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/1T46/ligands_raw/STI_A3.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

