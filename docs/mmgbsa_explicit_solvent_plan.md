# MMGBSA Explicit-Solvent Workflow Plan

## Goal

Make `atlas-mmgbsa batch-from-report` able to take selected heatmap/report hits and
run a publication-oriented explicit-solvent preparation, equilibration, production,
stripping, QC, and MMPBSA workflow with minimal per-target manual setup.

## Implementation Milestones

1. Preserve MMGBSA-ready pose artifacts.

   Report-driven MMGBSA should not depend on reconstructing ligand chemistry from
   retained PDBQT files. Docking retention should keep, or be able to restore,
   coordinate-matched SDF or MOL2 poses with authoritative bond orders and charge
   provenance. PDBQT-to-SDF conversion can remain a smoke-test fallback only.
   Status: partially wired in MMGBSA prep through
   `MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE`, chemistry-source provenance, and strict
   production gating. Update: Meeko SDF ligand prep now writes
   `*.ligprep_source.json`, successful Vina poses write coordinate-matched `*.sdf`
   and `*.mmgbsa_pose.json` beside the docked PDBQT, report-driven MMGBSA prefers
   those SDFs, and archive extraction strips retained `post_docked/<run_id>/...`
   or `docked/<run_id>/...` prefixes so restored artifacts land in the searchable
   run layout. Legacy runs without source sidecars still need regeneration from
   source chemistry.

2. Add explicit root and artifact controls to the batch CLI.

   The batch command should accept report CSVs, optional retained artifact archives,
   and explicit `--post-docked-dir`, `--docked-dir`, `--configs-dir`, and
   `--overall-dir` overrides so old runs and staged artifacts can be tested without
   editing global config files.

3. Build solvated systems.

   Add a dedicated explicit-system builder that writes and runs tleap for:

   - receptor/ligand complex topology from the existing MMGBSA receptor and ligand
     prep outputs,
   - `solvateBox` or `solvateOct` with configurable water model and buffer,
   - neutralization and optional salt concentration with `addIonsRand`,
   - separate dry-system topologies for MMPBSA after stripping solvent and ions.
   Status: implemented in `src/post_docking/mmgbsa/mmgbsa_explicit.py` and wired
   into production-mode MD replicates. Defaults use TIP3P, 10 A water buffer,
   neutralization, and 0.15 M salt.

4. Add explicit MD stages.

   Implement explicit-solvent minimization, heating, density equilibration,
   unrestrained equilibration, and production inputs. Settings should be generated
   from `MMGBSA_PROTOCOL=production` defaults but remain overrideable.
   Status: implemented with restrained minimization, unrestrained minimization,
   heating, density equilibration, unrestrained equilibration, and production
   `sander` stages.

5. Strip and window snapshots.

   Use cpptraj to strip water/ions from production trajectories, write stripped
   snapshots, and feed those snapshots into MMPBSA with the existing
   `MMGBSA_ANALYSIS_START_PS`, `MMGBSA_ANALYSIS_END_PS`, and interval controls.
   Status: implemented through explicit-solvent `cpptraj` stripping before the
   existing MMPBSA runner consumes the dry complex/receptor/ligand topologies.

6. Promote QC to first-class output.

   Extend current QC with explicit-solvent metadata: complex RMSD, ligand RMSD,
   protein RMSF, box/density checks, temperature/pressure summaries, per-frame
   energy convergence, block averages, and replicate statistics.

7. Harden chemistry review gates.

   Production mode should require recorded ligand charge/protonation/tautomer
   decisions, metal/water retention decisions, APO/HOLO and pH provenance, and a
   methods JSON that records force fields, water model, salt, frame window, and
   trajectory stripping masks. Status: ligand subgates and methods JSON coverage
   are wired for protonation, tautomer, stereochemistry, net charge, parameter
   review, input chemistry authority, force field, charge method, and precharged
   RESP/QM MOL2 source when supplied. Water/metal retention provenance is present
   in per-target/per-ligand review records under `mmgbsa/review_records/`. Atlas
   consumes precomputed RESP/QM MOL2 files and attempts internal AmberTools RESP
   fitting when needed; non-RESP fallback charge methods are flagged as not
   publication-ready.

8. Add tests and short smoke profiles.

   Add unit tests for tleap script generation, path/root overrides, archive
   extraction, snapshot-window conversion, and QC summaries. Add a no-MD smoke
   that runs ligand prep plus tleap on one restored report hit, and a tiny explicit
   MD smoke that runs only a few steps under a test preset.

## Pilot-Like Smoke Result

No literal `pilotstudy` artifact directory was found under the known workspace roots
during this pass. A retained pilot-like run was available at
`20260427_163422/TEST/HOLO/pH7_0`.

The batch workflow correctly selected:

- run: `20260427_163422`
- target: `TEST`
- variant: `HOLO`
- pH: `pH7_0`
- stage: `stage3`
- ligand: `actives_test_library_10_2`

The retained artifact archive was `artifacts.tar.zst`; extraction now supports this
format. A temporary PDBQT-to-SDF conversion allowed the smoke to verify report
selection, receptor matching, ligand prep, and tleap invocation. LEaP then reported
missing parameters from the converted ligand chemistry, which confirms why retained
SDF/MOL2 pose artifacts are required for publication-grade MMGBSA.
