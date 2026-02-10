# Artifact Retention

## Overview

`post_docking.artifact_retention` archives bulky per-ligand artifacts while keeping CSV outputs unarchived.

Key properties:
- Default mode is `rerun_safe`.
- Archive format is per-group `.tar.zst` (not one giant run tarball).
- Verification is fail-closed (no deletion when verification fails).
- Idempotent behavior on repeated runs.

## Control Model

Single config/env key:
- `ARTIFACT_RETENTION`
- allowed values: `rerun_safe`, `minimal_disk`, `off`
- default when missing: `rerun_safe`

CLI flags:
- `--retain` (force enable)
- `--noretain` (force disable)
- `--retainmode <rerun_safe|minimal_disk>` (force mode)

Precedence:
- CLI > ENV > CONFIG > default

ENV uses the same key/value model:
- `ARTIFACT_RETENTION=rerun_safe|minimal_disk|off`

## Timing

Retention now runs per-protein, immediately after successful protein completion in `main.py`.

Behavior:
- Success path: retention invoked for that protein scope (`--pdb-id <PDB>`).
- Failure path: retention is not called.
- Parallel protein mode: invocations are serialized through a lock to avoid archive/index races.

Run-level retention is no longer the primary path.

## What Is Archived

Archive targets (per-ligand bulky artifacts) under `docked/<RUN_ID>` and `post_docked/<RUN_ID>`:
- `*.sdf`, `*.sdf.gz`, `*.sd`
- `*.pdbqt`, `*.pdbqt.gz`
- `*.mol2`, `*.mol2.gz`
- `*.dok`, `*.dlg`
- `*.mae`, `*.maegz`

CSV outputs are always kept unarchived.

## Grouping and Manifest

Archive key:
- `run_id`, `source_root`, `pdb_id`, `variant`, `ph`

Each group has:
- `artifacts.tar.zst`

Run-level manifest/index (single file per run):
- `docked/<RUN_ID>/_artifact_archives/archive_index.json`

The run-level index contains all group records and per-entry metadata
(including `stage_dir` and `mode`), so per-archive manifests are no longer
written for new runs.

Manifest entry fields:
- `run_id`, `pdb_id`, `variant`, `ph`, `stage_dir`, `mode`
- `original_path`, `archive_path`, `member_name`
- `size_bytes`, `mtime`, `sha256`, `file_type`

## Modes

`rerun_safe`:
- archive + verify + delete archived originals after successful verification.
- keep CSV outputs untouched.

`minimal_disk`:
- runs `rerun_safe` behavior first.
- then deletes only extra temp-like files (`.tmp`, `.part`, `.bak`, `.lock`).
- does not delete `_DONE`, `_done`, or `completion_*.json` by default.

`off`:
- skip retention.

## Safety and Edge Cases

- Verified-index mismatch handling: if discovered files do not match verified index members, the group is rebuilt instead of skipped.
- Tool preflight: retention checks `tar`, `zstd`, and `tar --zstd` support before any archive/delete work. On failure, it skips safely.
- Legacy compatibility: verify/restore can still read old per-archive manifests.
- Restore hardening:
  - rejects empty targets
  - rejects directory targets
  - requires target path to resolve under `docked/<run_id>` or `post_docked/<run_id>`
  - verifies checksums after restore

## CLI Usage

Direct retention:

```bash
python -m post_docking.artifact_retention --run-id <RUN_ID> --mode rerun_safe
```

Per-protein scope:

```bash
python -m post_docking.artifact_retention --run-id <RUN_ID> --pdb-id <PDB_ID> --mode rerun_safe
```

Useful flags:
- `--dry-run`
- `--overwrite`
- `--include-glob` / `--exclude-glob`
- `--threads`
- `--verify-only`
- `--restore --archive <path>`
- `--restore --run-id <RUN_ID> --restore-group <group-id>`
- `--restore-stage <stage_dir>`
- `--restore-mode <mode>`
- `--restore-path-glob <glob>` (repeatable)

## Restore

By archive path:

```bash
python -m post_docking.artifact_retention \
  --restore \
  --archive docked/<RUN_ID>/_artifact_archives/<...>/artifacts.tar.zst
```

By group id:

```bash
python -m post_docking.artifact_retention \
  --restore \
  --run-id <RUN_ID> \
  --restore-group "run=<RUN_ID>|source=docked|pdb=...|variant=...|ph=..."
```
