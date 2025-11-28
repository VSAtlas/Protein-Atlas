# Atlas Log Resurfacer - Summary

*Fatal:* 0  *Hard:* 8  *Soft:* 0  *Info:* 0

**Health score:** 20 / 100

## 1. [HARD] ligprep.rdkit_valence - 
- **Message:** RDKit sanitize/valence/kekulize error(s)
- **Count:** 8  |  **Recency:** 
- **Hint:** Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.
  - **File:** `/stor/home/mpg2352/atlas/code/protein_automation/logs/bench_20251013_152831.log`
    - `[15:29:06] Explicit valence for atom # 18 C, 5, is greater than permitted`
    - `[15:29:06] ERROR: Explicit valence for atom # 18 C, 5, is greater than permitted`
    - `[15:29:07] Explicit valence for atom # 18 C, 5, is greater than permitted`
- **Step:** ligprep.sanitize

