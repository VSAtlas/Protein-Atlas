# Atlas Log Resurfacer - Summary

*Fatal:* 0  *Hard:* 10  *Soft:* 2  *Info:* 0

**Health score:** 0 / 100

## 1. [HARD] ligprep.rdkit_valence - 2025
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 9  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251009_194106.log`
    - `[19:42:15] Explicit valence for atom # 10 O, 3, is greater than permitted`
    - `[19:42:15] Explicit valence for atom # 10 O, 3, is greater than permitted`
    - `[19:42:15] Explicit valence for atom # 7 O, 3, is greater than permitted`

## 2. [HARD] ligprep.mgl_partial_write - 2025
- **Message:** Partial write during conversion
- **Count:** 1  |  **Recency:** 
- **Hint:** Conversion partial; capture stderr, retry fallback, then quarantine if persistent.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251009_194106.log`
    - `WARNING: 63 atoms of 69 in CDL_B611.sanitized.mol2  were not written`

## 3. [SOFT] ligprep.obabel_h_charge - 6ADQ
- **Message:** OpenBabel H/charge warning(s)
- **Count:** 2  |  **Recency:** 
- **Hint:** Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251009_194106.log`
    - `[19:42:16] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/9Y0_B612.sanitized.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`
    - `[19:42:17] /stor/home/mpg2352/atlas/code/protein_automation/processed_pdbs/6ADQ/ligands_raw/CDL_B611.sanitized.pdb: Warning - no explicit hydrogens in mol2 file but needed for formal charge estimation.`

