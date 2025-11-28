# Atlas Log Resurfacer - Summary

*Fatal:* 0  *Hard:* 10  *Soft:* 8  *Info:* 0

**Health score:** 0 / 100

## 1. [HARD] ligprep.rdkit_valence - 2025
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 9  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `logs/bench_20251010_131013.log`
    - `[13:11:51] Explicit valence for atom # 10 O, 3, is greater than permitted`
    - `[13:11:58] Explicit valence for atom # 10 O, 3, is greater than permitted`
    - `[13:12:01] Explicit valence for atom # 7 O, 3, is greater than permitted`
- **Step:** ligprep.sanitize
- **Paths:**
  - bench_logs: `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251010_131013.log`

## 2. [HARD] ligprep.mgl_partial_write - 2025 - CDL_N605.sanitized
- **Message:** Partial write during conversion
- **Count:** 1  |  **Recency:** 
- **Hint:** Conversion partial; capture stderr, retry fallback, then quarantine if persistent.
  - **File:** `logs/bench_20251010_131013.log`
    - `WARNING: 70 atoms of 80 in CDL_N605.sanitized.mol2  were not written`
- **Step:** ligprep.write_mol2

## 3. [SOFT] ligprep.obabel_h_charge - 6ADQ
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 8  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `logs/bench_20251010_131013.log`
    - `[13:12:05] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_B608.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[13:12:18] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9Y0_W201.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[13:12:51] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/MQ9_B609.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
- **Step:** ligprep.add_h

