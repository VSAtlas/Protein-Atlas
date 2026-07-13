# Historical SPD90 release backfill

This is an engineering handoff for
`spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503`. It does not approve a
receptor, pose, native redock, score source, or public artifact. Record all
scientific approvals in the [decision register](docking_atlas_decision_register.md)
and explicit annotation files. Those files may remain external to the manifest, but
Atlas hashes their content and stores their source-record provenance during the
build.

## Frozen evidence inventory

Canonical status reports 8,010 of 8,010 distributed chunks complete with no
recorded run failures. The failure-complete release audit found:

- 90 legacy/blank-pH receptor contexts;
- 440,183 master result rows, all with blank `pose_valid_any`;
- 176 unvalidated control rows across 53 contexts;
- 935,313 retained PDBQT files, 413 SDF files, and 40,948 completion JSON files;
- 160 named control-redock PDBQT outputs across 55 targets;
- no `posebusters_all_stages.csv`, `archive_index.json`, or `.tar.zst` archive;
- no structured native-redock RMSD or qualification table.

An engineering checkpoint produced an approximately 439 MB compact public SQLite
projection. That file exceeds Cloudflare Pages' current 25 MiB per-asset limit.
Either serve the complete verified site from a host that accepts the file, or
externalize the download to object storage or another host through a supported link
projection before using Pages.

The database importer also reconstructs 320,695 expected pairs with no master
result as structured missing-result cells. A failed calculation with no pose is not
an invalid pose; preserve that distinction during every backfill.

## Pair pose validation

A full docking rerun is not the first step because pose files survive. Do not run
the quarantined `tools/pose_bust.py` directly: it assumes a
`PDB/variant/pH/stage` layout, while this run uses legacy `PDB/stage` locations and
must join exactly to `(PDB, LEGACY, blank pH, ligand)`.

A supported backfill command does not exist yet. As of 2026-07-13,
`atlas analysis --help` does not list `pose-validate`. Implement and
regression-test a canonical `atlas analysis pose-validate` wrapper before
processing this run. It must:

1. Resolve receptor, ligand, pose, and context paths from frozen run records.
2. Treat `LEGACY` and blank pH as explicit values.
3. Resume without reconverting completed poses.
4. Write one structured record per expected pair, including status, reason, source
   pose, input hashes, validator version, and configuration.
5. Separate absent calculations from present-but-invalid poses.
6. Use no more than 32 workers.

### Design-only pose-validation interface (not implemented)

The intended future interface is shown for design review only. Do not run or
script against it: the current CLI rejects this subcommand.

```bash
# DESIGN ONLY — `atlas analysis pose-validate` is not implemented.
atlas analysis pose-validate \
  spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503 \
  --docked-root outputs/spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503/outputs/docked \
  --post-docked-root outputs/spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503/outputs/post_docked \
  --legacy-context \
  --workers 32 \
  --resume
```

After validation joins have been checked, regenerate the canonical report/table:

```bash
atlas analysis report \
  spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503 \
  --overwrite-master \
  --status-html
```

Only pairs whose required pose artifact is absent or unusable need targeted
docking reruns.

## Native-redock qualification

Historical control-center code computed RMSD transiently but did not persist an
auditable outcome. Do not convert a control score, filename, or apparent pose into
a pass/fail token.

Analysis-only recomputation can avoid redocking when the exact prepared receptor,
native crystal-ligand reference, docked control pose, atom mapping, search box,
configuration, and source hashes all survive. Add a maintained RMSD-evidence
command that stores those inputs and the user-approved calculation method. If an
exact reference or control pose is missing, rerun only that native control.

At least 37 contexts lack a control row in the master export, and named control
artifacts cover only 55 targets. APO contexts or structures without a native ligand
require an explicit user policy. Store reviewed outcomes through the receptor
annotation template; only an approved outcome receives a controlled qualification
token.

## Pair artifact archive and index

Artifact indexing can be recovered without docking. First capture a size baseline:

```bash
atlas artifacts measure \
  --run-id spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503 \
  --root outputs/spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503/outputs/docked/spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503 \
  --root outputs/spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503/outputs/post_docked/spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503 \
  --capture baseline \
  --out outputs/spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503/outputs/data/spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503/artifact_reports
```

Atlas does not yet expose retention through `atlas artifacts`; use the documented
module only after reviewing the dry-run and available capacity:

```bash
python -m post_docking.artifact_retention \
  --run-id spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503 \
  --mode rerun_safe \
  --dry-run \
  --threads 32 \
  --group-workers 1 \
  --docked-root outputs/spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503/outputs/docked \
  --post-docked-root outputs/spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503/outputs/post_docked
```

After explicit review, remove `--dry-run`. `rerun_safe` deletes originals only
after verified archive creation. Verify the result with the same roots:

```bash
python -m post_docking.artifact_retention \
  --run-id spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503 \
  --mode rerun_safe \
  --verify-only \
  --threads 32 \
  --group-workers 1 \
  --docked-root outputs/spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503/outputs/docked \
  --post-docked-root outputs/spd90_fda_fda_dud_spr10_40h_patchdist_20260525_161503/outputs/post_docked
```

The verified `archive_index.json` then becomes the release manifest's
`paths.archive_index`. `atlas reproduce bundle` remains complementary: it indexes
top-level reproducibility/report outputs, not the pair-level pose archive.
