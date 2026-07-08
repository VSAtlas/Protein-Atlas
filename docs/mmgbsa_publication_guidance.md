# MM/GBSA Follow-Up Workflow

MM/GBSA is now treated as follow-up evidence for selected hits, not as an
automatic post-run pass over every docking output. Normal pipeline finalization
skips MM/GBSA unless `MMGBSA_AUTO_RUN=true`; selected hits can be run with:

```bash
atlas-mmgbsa batch-from-report \
  --report analysis/figures/<run_id>/significant_hits.csv \
  --run-id <run_id> \
  --require-significant \
  --dry-run
```

Remove `--dry-run` after the plan points at the expected PDB, variant, pH,
stage, and ligand IDs. Artifact archives can be expanded first with repeated
`--artifact <archive>` arguments; archives are checked for path traversal before
extraction.

For a publication-oriented follow-up run, start with:

```txt
MMGBSA_PROTOCOL=production
MMGBSA_CHEMISTRY_REVIEW_STATUS=reviewed
MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE=true
MMGBSA_INPUT_CHEMISTRY_SOURCE=curated_sdf_ligprep_pH7.4
MMGBSA_PROTONATION_REVIEW_STATUS=reviewed
MMGBSA_TAUTOMER_REVIEW_STATUS=reviewed
MMGBSA_STEREOCHEMISTRY_REVIEW_STATUS=reviewed
MMGBSA_NET_CHARGE_REVIEW_STATUS=reviewed
MMGBSA_PARAMETER_REVIEW_STATUS=reviewed
```

`production` expands to explicit-solvent ensemble defaults: 3 replicates,
TIP3P solvation, neutralization plus 0.15 M salt, restrained minimization,
unrestrained minimization, heating, density equilibration, unrestrained
equilibration, 10000 ps production, 20 ps frame spacing, stripped snapshots for
MMPBSA after a 2000 ps analysis burn-in, frame/block QC, and replicate
uncertainty reporting. `screening` keeps the lightweight one-frame behavior.

## Publication-Ready Expectations

Ligand preparation should start from chemically meaningful structures: curated
SDF/MOL2, ligand-only PDB with correct hydrogens and bonding, or a crystal
ligand that has been deliberately protonated and assigned a tautomer,
stereochemistry, and net charge. Do not treat a PDBQT-derived SDF as
chemistry-authoritative without review; PDBQT omits bond-order and charge
provenance that MM/GBSA parameterization depends on.

Atlas now records ligand-prep provenance next to each prepared ligand in
`*.mmgbsa_prep.json`: input SDF, chemistry source, GAFF/GAFF2 atom typing,
AM1-BCC versus gas fallback, RDKit charge/radical checks, `parmchk2` warnings,
and whether the ligand prep passed the current publication gate. In production
mode, strict ligand chemistry rejects unreviewed/non-authoritative inputs, gas
fallback charges, unknown atom types, `parmchk2` warnings marked for review, and
unreviewed protonation/tautomer/stereo/net-charge/parameter subgates.

For future Vina runs prepared through the Meeko SDF path, Atlas writes
`*.ligprep_source.json` beside each prepared PDBQT and `*.sdf` plus
`*.mmgbsa_pose.json` beside each successful docked PDBQT. These docked-pose SDFs
reuse the source SDF chemistry and replace only the coordinates with the first
docked pose. Report-driven MM/GBSA prefers these authoritative docked-pose SDFs,
including after retained artifact archive extraction. PDBQT-to-SDF conversion
should remain a smoke-test fallback because it reconstructs bond orders/charges
after the fact.

Report rows can also provide authoritative source chemistry directly with
`sdf_path`/`mol2_path` columns. This is how Wang/PDBbind control staging now
works: `atlas-mmgbsa benchmark-wang2016 --pdbbind-root ...` or
`--pdbbind-zip ...` writes `wang2016_pdbbind_ligands.csv`, and
`atlas-mmgbsa batch-from-report --report wang2016_pdbbind_ligands.csv ...`
materializes those SDF/MOL2 inputs into the normal post-docked stage before
MM/GBSA. If an SDF does not sanitize but the paired MOL2 is present, Atlas
derives the staged SDF from MOL2 and records the source in
`*.mmgbsa_pose.json`.

Each prepared ligand now gets first-class ligand-state provenance in
`*.mmgbsa_prep.json` and the MM/GBSA review record: RDKit formal charge,
radicals, fragments, element set, tautomer enumeration, optional Dimorphite-DL
protomer enumeration near the target pH, and a severity-ranked
`suspicious_chemistry` classifier. Production mode blocks `blocker` and
`critical` chemistry findings rather than silently treating ambiguous ligand
states as publishable.

Ligand-state inference reuses the same pH context stack as receptor prep. If an
MM/GBSA pH label is present, Atlas parses it with the docking pH-ensemble label
parser. Otherwise it uses `protein_prep.openmm_repair.get_target_ph_for_prep`,
which delegates to the shared PDB header/context pH selector. The inference
preserves authoritative crystal/PDBbind input chemistry by default, ranks
Dimorphite-DL protomers and RDKit tautomers with chemistry filters and simple
evidence terms, and records the selected state plus ranking margin. Ambiguous
state ranking is a warning in screening mode and a production blocker.

For RESP/QM-derived charges, place the precharged MOL2 next to the selected SDF
as `<ligand>.resp.mol2`, `<ligand>.qm_resp.mol2`, `<ligand>.qm.mol2`, or
`<ligand>.precharged.mol2`; `resp/`, `qm/`, and `mol2/` subdirectories beside
the SDF are also searched. Set the existing ligand charge method to `resp`,
`qm_resp`, `qm`, or `precharged`; Atlas copies that MOL2, records it in
`*.mmgbsa_prep.json`, and runs `parmchk2` to validate/write the frcmod. If no
precharged MOL2 is found, Atlas attempts AmberTools internal `antechamber -c
resp` and records the command logs. A Mulliken fallback is recorded as
non-publication-ready and is blocked by strict production mode.

For true QM/RESP handoff, use:

```bash
atlas-mmgbsa resp-plan \
  --ligand path/to/authoritative_pose.sdf \
  --out-dir /stor/home/mpg2352/mmgbsa_resp/LIG \
  --net-charge 0 \
  --run-antechamber
```

This writes `resp_workflow.json`, a Gaussian-style AmberTools input
`*.resp.gcrt`, and two scripts: `01_make_qm_input.sh` and
`02_finalize_resp.sh`. Run and review the QM calculation externally, place the
completed Gaussian/GAMESS output at the path recorded in `resp_workflow.json`,
then run:

```bash
atlas-mmgbsa resp-finalize \
  --qm-output /stor/home/mpg2352/mmgbsa_resp/LIG/LIG.resp.log \
  --out-dir /stor/home/mpg2352/mmgbsa_resp/LIG \
  --net-charge 0
```

The finalized `*.resp.mol2` and `*.resp.frcmod` can then be placed beside the
authoritative SDF for production MM/GBSA. This script makes the workflow
auditable, but it does not remove the need to review protonation, tautomer,
stereochemistry, net charge, QM settings, and ESP/RESP quality.

Each ligand now gets a structured review record under
`mmgbsa/review_records/<target>/<ligand>.mmgbsa_review.json`. The record captures
ligand protonation, tautomer, stereochemistry, net charge, parameter review,
input chemistry authority, water policy, metal policy, APO/HOLO, pH, and notes.
Production validation requires these sections to be reviewed/approved/ok and the
input chemistry to be authoritative.

Amber's MMPBSA.py tutorials describe MM/PB(GB)SA as post-processing over
representative trajectory snapshots using complex, receptor, ligand topologies
and trajectory files:
https://ambermd.org/tutorials/advanced/tutorial3/py_script/section1.php

For publication use, treat docking-pose one-frame MM/GBSA as triage only. A
stronger run should document:

- complete AmberTools version/prefix and exact `MMPBSA.py` command inputs;
- ligand protonation/tautomer and net charge handling, including failed
  `antechamber` fallbacks;
- receptor APO/HOLO and pH state used for each hit;
- MD equilibration and production protocol before endpoint analysis;
- number of frames, frame spacing, convergence behavior, block averages,
  replicate uncertainty, and replicate aggregation;
- GB model (`igb`), salt concentration, and whether PB and entropy terms were
  run;
- whether waters/metals were retained or stripped and why.

The current Atlas implementation differs from publication-grade MM/GBSA in
important ways:

- default `ONEFRAME` mode converts a single docked complex coordinate into a
  one-frame trajectory, which is useful for regression testing but not a
  converged ensemble;
- `MMGBSA_PROTOCOL=screening` keeps short throughput-oriented defaults, while
  `MMGBSA_PROTOCOL=production` enables longer MD, frame-window controls, QC
  artifacts, and stricter validation;
- entropy is not enabled by default, so results are endpoint binding-energy
  estimates rather than full binding free energies;
- ligand and receptor preparation are still automated, so production mode
  requires an explicit chemistry-review acknowledgment for protonation, metal
  coordination, strained poses, and active-site waters;
- existing legacy runs do not automatically gain authoritative docked-pose SDFs
  unless their ligand PDBQTs have `*.ligprep_source.json` sidecars or the
  artifacts are regenerated from source SDF libraries;
- RESP/QM-derived ligand charge fitting now has an AmberTools internal attempt,
  but precomputed RESP/QM MOL2 remains preferred for publication work because
  RESP setup details should be reviewed and preserved.
- explicit-solvent system building, equilibration, production, and trajectory
  stripping are implemented as the production default, but long production runs
  still need real-system validation before manuscript use.
- ligand protomer/tautomer inference now chooses a default state and records the
  ranking, but Atlas does not yet generate and run every plausible state through
  redocking or free-energy comparison. Production mode should treat ambiguous
  charge/protomer/tautomer findings as a review queue or multi-state run.

## Benchmarking Against Published Work

A practical benchmark target is Wang et al. 2016, "Calculating protein-ligand
binding affinities with MMPBSA: Method and Error Analysis"
(https://pmc.ncbi.nlm.nih.gov/articles/PMC5018451/). The paper selected
PDBbind complexes with experimental affinities across trypsin beta, thrombin
alpha, CDK/PKA, urokinase-type plasminogen activator, beta-glucosidase A, and
factor Xa. It is useful because it compares computed relative binding affinities
against experiment across receptor families rather than a single cherry-picked
complex.

For Atlas, the benchmark should be implemented as a curated fixture manifest:
PDB ID, ligand identifier, receptor family, experimental affinity, pH/protonation
notes, retained crystallographic waters/ions, and expected family assignment.
Run Atlas with `MMGBSA_PROTOCOL=production`, then compare per-family Spearman or
Pearson correlation and RMSE/MAE against experiment. Exact energy values should
not be expected to match the paper unless force fields, charge derivation,
water/ion retention, trajectory length, PB/GB model, dielectric settings, and
snapshot windows are deliberately matched.

Atlas includes a benchmark harness:

```bash
atlas-mmgbsa benchmark-wang2016 \
  --out-dir /stor/home/mpg2352/mmgbsa_benchmarks/wang2016
```

The command downloads the public supplementary DOCX when reachable, extracts all
tables to CSV/JSON, writes the paper-level Table 1 and Table 8 target metrics,
and emits `wang2016_benchmark_manifest.json`. If you have an Atlas MM/GBSA CSV
with matching PDB IDs and a `delta_total` column, pass it with
`--atlas-results`; if the parsed supplement table uses a different reference
energy column, set `--reference-energy-column`.

To prepare the representative crystal/redocking fixture set from the paper:

```bash
atlas-mmgbsa benchmark-wang2016 \
  --out-dir /stor/home/mpg2352/mmgbsa_benchmarks/wang2016 \
  --download-pdbs
```

This downloads `1Q8W`, `2J75`, `1VZQ`, `1O2Q`, `1C5Z`, and `1LPZ` into
`<out-dir>/input_pdbs/` and writes `wang2016_downloaded_pdbs.csv`. These
structures are useful for two related checks: Atlas MM/GBSA on crystal-like
control poses, and Atlas control-redocking followed by MM/GBSA on docked poses.
If a local PDBbind/PDBbind-derived tree or zip is available, pass
`--pdbbind-root <dir>` or `--pdbbind-zip <zip>` to copy matching representative
`*_ligand.sdf` and `*_ligand.mol2` files into `<out-dir>/pdbbind_ligands/` and
write `wang2016_pdbbind_ligands.csv`. Representative controls without a matched
PDBbind ligand source, such as `1C5Z` in the PDBbind++ 2020 archive, are skipped
from the batch-ready CSV and recorded separately in
`wang2016_pdbbind_ligands_skipped.csv`.

The Amber manuals page is the canonical current reference for AmberTools
options: https://ambermd.org/Manuals.php. The original MMPBSA.py publication
describes the method and its parallel frame distribution:
https://doi.org/10.1021/ct300418h.

Additional ligand-prep references: Amber's GAFF/antechamber tutorials document
the standard `antechamber` and `parmchk2` workflow
(https://ambermd.org/tutorials/basic/tutorial4b/index.php and
https://ambermd.org/tutorials/advanced/tutorial39/index.php), while tools such
as Schrödinger LigPrep/Epik are commonly used to enumerate and assign
protonation/tautomer states before MD.
