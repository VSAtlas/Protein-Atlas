# Artifact Retention Validation Report

Date: 2026-02-09  
Run ID used: `retention_fixture_20260209`

## Notes

- Existing checked-in run folders in this workspace had almost no `.sdf/.pdbqt/.mol2` artifacts.
- To satisfy objective before/after thresholds, a dedicated fixture run-id was created with representative per-stage ligand files plus CSV outputs.

## Commands and Outputs

### 1) Baseline measurement

Command:

```bash
python measure_artifacts.py \
  --run-id retention_fixture_20260209 \
  --root docked/retention_fixture_20260209 \
  --root post_docked/retention_fixture_20260209 \
  --out size_artifacts/retention_fixture_20260209/baseline_metrics
```

Output (key lines):

```text
Total files: 11
Total bytes: 266
Top extensions by bytes:
  .sdf files=2
  .pdbqt files=3
  .csv files=3
  .mol2 files=1
```

### 2) Dry-run retention

Command:

```bash
python -m post_docking.artifact_retention \
  --run-id retention_fixture_20260209 \
  --mode rerun_safe \
  --dry-run \
  --threads 2
```

Output (key lines):

```text
groups_discovered=4
files_discovered=6
groups_planned=4
files_deleted=0
```

Re-measure command:

```bash
python measure_artifacts.py \
  --run-id retention_fixture_20260209 \
  --root docked/retention_fixture_20260209 \
  --root post_docked/retention_fixture_20260209 \
  --out size_artifacts/retention_fixture_20260209/after_dry_run_metrics
```

Output (key lines):

```text
Total files: 11
Total bytes: 266
```

### 3) Apply `rerun_safe`

Command:

```bash
python -m post_docking.artifact_retention \
  --run-id retention_fixture_20260209 \
  --mode rerun_safe \
  --threads 2
```

Output (key lines):

```text
groups_archived=4
groups_verified=4
files_archived=6
files_deleted=6
groups_failed=0
```

Re-measure command:

```bash
python measure_artifacts.py \
  --run-id retention_fixture_20260209 \
  --root docked/retention_fixture_20260209 \
  --root post_docked/retention_fixture_20260209 \
  --out size_artifacts/retention_fixture_20260209/after_rerun_safe_metrics
```

Output (key lines):

```text
Total files: 14
Top extensions:
  .zst files=4
  .json files=6
  .csv files=3
  .pdbqt/.sdf/.mol2 files=0
```

Archive presence checks:

```text
archive_files=4
manifest_files=4
archive_index.json exists=true
```

Verification command:

```bash
python -m post_docking.artifact_retention \
  --run-id retention_fixture_20260209 \
  --verify-only
```

Output (key lines):

```text
groups_verified=4
groups_failed=0
```

### 4) Idempotency rerun (no overwrite)

Command:

```bash
python -m post_docking.artifact_retention \
  --run-id retention_fixture_20260209 \
  --mode rerun_safe \
  --threads 2
```

Output (key lines):

```text
groups_discovered=0
files_discovered=0
groups_archived=0
files_deleted=0
```

Re-measure command:

```bash
python measure_artifacts.py \
  --run-id retention_fixture_20260209 \
  --root docked/retention_fixture_20260209 \
  --root post_docked/retention_fixture_20260209 \
  --out size_artifacts/retention_fixture_20260209/after_idempotency_metrics
```

Output (key lines):

```text
Total files: 14
Total bytes: 12480
```

### 5) Restore spot-check

Command:

```bash
python -m post_docking.artifact_retention \
  --restore \
  --run-id retention_fixture_20260209 \
  --restore-group "run=retention_fixture_20260209|source=docked|pdb=TEST|variant=HOLO|ph=pH7_0|stage=stage1|mode=unknown"
```

Output:

```text
restored=2
failed=0
```

Manifest consistency check:

```text
manifest_entries=2
restored_exists=2
checksum_mismatch=0
```

### 6) Optional `minimal_disk`

Command:

```bash
python -m post_docking.artifact_retention \
  --run-id retention_fixture_20260209 \
  --mode minimal_disk \
  --threads 2
```

Output (key lines):

```text
extras_deleted=2
groups_failed=0
```

Re-measure command:

```bash
python measure_artifacts.py \
  --run-id retention_fixture_20260209 \
  --root docked/retention_fixture_20260209 \
  --root post_docked/retention_fixture_20260209 \
  --out size_artifacts/retention_fixture_20260209/after_minimal_disk_metrics
```

Output (key lines):

```text
csv files=3 (unchanged)
json files 6 -> 5 (extra non-essential reduced)
```

## Acceptance Thresholds

- Dry-run: no count changes -> PASS
- Rerun-safe: CSV count delta = 0 -> PASS
- Rerun-safe: raw per-ligand artifacts reduced and archive count > 0 -> PASS
- Idempotency: second run without overwrite has no effective changes -> PASS
- Restore: restored files match manifest count + checksum -> PASS
- Minimal disk (optional): additional non-essential files reduced while CSV unchanged -> PASS

## Lint and Types

### Ruff

Command:

```bash
ruff check --fix .
```

Result:

- Fails due large pre-existing repository-wide Ruff debt unrelated to this patch set.
- Focused checks for touched Python files pass:

```bash
ruff check src/post_docking/artifact_retention.py postrun_hooks.py
```

### MyPy

Baseline command:

```bash
mypy src/post_docking/artifact_retention.py postrun_hooks.py main.py
```

Result:

- Fails due environment/module-layout constraints (`types-tqdm` missing and duplicate module discovery).

Focused checks used for touched files:

```bash
mypy --ignore-missing-imports src/post_docking/artifact_retention.py
mypy --ignore-missing-imports postrun_hooks.py
mypy --ignore-missing-imports --follow-imports=skip main.py
```

Result:

- All three succeeded.
