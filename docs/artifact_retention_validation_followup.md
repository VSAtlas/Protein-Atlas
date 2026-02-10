# Artifact Retention Follow-up Validation

Date: 2026-02-09

## Scope

Validated follow-up requirements:
- single-key control model (`ARTIFACT_RETENTION`)
- CLI flags (`--retain`, `--noretain`, `--retainmode`)
- per-protein invocation timing
- failure-path skip behavior
- parallel safety via serialized retention lock
- required bug/edge-case fixes

## Scenario 1: Default behavior (no retention flags)

Fixture run id: `retention_followup_default`

Baseline:

```bash
python measure_artifacts.py \
  --run-id retention_followup_default \
  --root docked/retention_followup_default \
  --root post_docked/retention_followup_default \
  --out size_artifacts/retention_followup_default/baseline
```

Key totals:
- files=9
- `.csv`=4, `.pdbqt`=2, `.sdf`=2

Per-protein retention call (simulates success path call in `main.py`):

```bash
python - <<'PY'
import threading
from postrun_hooks import _maybe_run_artifact_retention_for_pdb
_maybe_run_artifact_retention_for_pdb({}, "retention_followup_default", "PDBA", lock=threading.Lock())
PY
```

After:

```bash
python measure_artifacts.py \
  --run-id retention_followup_default \
  --root docked/retention_followup_default \
  --root post_docked/retention_followup_default \
  --out size_artifacts/retention_followup_default/after_default_per_protein
```

Key totals:
- files=12
- `.csv`=4 (unchanged)
- `.pdbqt`=1 and `.sdf`=1 (only failed/unprocessed protein remains raw)
- `.zst`=2, `.json`=4

Result: PASS

## Scenario 2: `--noretain` skips retention

Fixture run id: `retention_followup_noretain`

Resolved mode check:

```bash
python - <<'PY'
from main import _resolve_artifact_retention_mode
print(_resolve_artifact_retention_mode({}, ["main.py", "--noretain"]))
PY
```

Output:
- `('off', 'CLI(--noretain)')`

No-retain invocation:

```bash
python - <<'PY'
from main import _resolve_artifact_retention_mode
from postrun_hooks import _maybe_run_artifact_retention_for_pdb
mode, _ = _resolve_artifact_retention_mode({}, ["main.py", "--noretain"])
_maybe_run_artifact_retention_for_pdb({"ARTIFACT_RETENTION": mode}, "retention_followup_noretain", "PDBC")
PY
```

Metrics unchanged:
- baseline: files=4, bytes=105
- after: files=4, bytes=105

Result: PASS

## Scenario 3: `--retain --retainmode minimal_disk`

Fixture run id: `retention_followup_minimal`

Resolved mode check:

```bash
python - <<'PY'
from main import _resolve_artifact_retention_mode
print(_resolve_artifact_retention_mode({}, ["main.py", "--retain", "--retainmode", "minimal_disk"]))
PY
```

Output:
- `('minimal_disk', 'CLI(--retainmode)')`

Invocation:

```bash
python - <<'PY'
from main import _resolve_artifact_retention_mode
from postrun_hooks import _maybe_run_artifact_retention_for_pdb
mode, _ = _resolve_artifact_retention_mode({}, ["main.py", "--retain", "--retainmode", "minimal_disk"])
_maybe_run_artifact_retention_for_pdb({"ARTIFACT_RETENTION": mode}, "retention_followup_minimal", "PDBD")
PY
```

After metrics:
- `.csv`: 2 -> 2 (unchanged)
- `.tmp`: 1 -> 0 (extra temp removed)
- `_DONE` count (extension `(none)`): 1 -> 1 (preserved)
- `completion_vina.json` preserved

Result: PASS

## Scenario 4: Dry-run no destructive change

Fixture run id: `retention_followup_dryrun`

Command:

```bash
python -m post_docking.artifact_retention \
  --run-id retention_followup_dryrun \
  --pdb-id PDBE \
  --mode rerun_safe \
  --dry-run \
  --threads 2
```

Metrics:
- baseline: files=4, bytes=105
- after dry-run: files=4, bytes=105

Result: PASS

## Scenario 5: Idempotency (no churn)

Run id: `retention_followup_default`

Second call:

```bash
python - <<'PY'
import threading
from postrun_hooks import _maybe_run_artifact_retention_for_pdb
_maybe_run_artifact_retention_for_pdb({}, "retention_followup_default", "PDBA", lock=threading.Lock())
PY
```

Observed:
- summary reports `files_discovered=0`, `groups_discovered=0`
- metrics unchanged from first retained state

Result: PASS

## Scenario 6: Restore-group checksum match

Run id: `retention_followup_default`

Restore:

```bash
python -m post_docking.artifact_retention \
  --restore \
  --run-id retention_followup_default \
  --restore-group "run=retention_followup_default|source=docked|pdb=PDBA|variant=HOLO|ph=pH7_0|stage=stage1|mode=unknown"
```

Output:
- `restored=1`, `failed=0`

Checksum check:

```bash
python - <<'PY'
import json, hashlib
from pathlib import Path
m = Path("docked/retention_followup_default/_artifact_archives/docked/PDBA/HOLO/pH7_0/stage1/unknown/artifacts.manifest.json")
entries = json.loads(m.read_text())["entries"]
ok = 0
for e in entries:
    p = Path(e["original_path"])
    if p.exists() and hashlib.sha256(p.read_bytes()).hexdigest() == e["sha256"]:
        ok += 1
print({"manifest_entries": len(entries), "checksum_match": ok})
PY
```

Output:
- `{'manifest_entries': 1, 'checksum_match': 1}`

Result: PASS

## Scenario 7: Failure path skips retention call

Run id: `retention_followup_default`

Only `PDBA` retention was invoked. `PDBB` simulates failed protein (no invocation).

Check:

```bash
python - <<'PY'
from pathlib import Path
print({
  "pdbb_docked_pdbqt_exists": any(p.suffix==".pdbqt" for p in Path("docked/retention_followup_default/PDBB").rglob("*")),
  "pdbb_post_sdf_exists": any(p.suffix==".sdf" for p in Path("post_docked/retention_followup_default/PDBB").rglob("*")),
  "pdbb_archives_present": any("PDBB" in str(p) for p in (Path("docked/retention_followup_default")/"_artifact_archives").rglob("*.tar.zst")),
})
PY
```

Output:
- `{'pdbb_docked_pdbqt_exists': True, 'pdbb_post_sdf_exists': True, 'pdbb_archives_present': False}`

Result: PASS

## Scenario 8: Parallel mode safety

Fixture run id: `retention_followup_parallel`

Parallel simulation with shared lock:

```bash
python - <<'PY'
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from postrun_hooks import _maybe_run_artifact_retention_for_pdb
run_id = "retention_followup_parallel"
lock = threading.Lock()
errors = []
with ThreadPoolExecutor(max_workers=3) as pool:
    futs = [pool.submit(_maybe_run_artifact_retention_for_pdb, {}, run_id, pdb, lock=lock) for pdb in ["PX1","PX2","PX3"]]
    for fut in as_completed(futs):
        try:
            fut.result()
        except Exception as exc:
            errors.append(str(exc))
print({"errors": len(errors), "details": errors})
PY
```

Output:
- `{'errors': 0, 'details': []}`

Index integrity:

```bash
python - <<'PY'
import json
from pathlib import Path
idx = Path("docked/retention_followup_parallel/_artifact_archives/archive_index.json")
data = json.loads(idx.read_text())
print({"index_exists": idx.exists(), "group_count": data.get("group_count"), "unique_groups": len({g.get("group_id") for g in data.get("groups", [])})})
PY
```

Output:
- `{'index_exists': True, 'group_count': 6, 'unique_groups': 6}`

Result: PASS

## Edge-case Fix Validation

### Verified-manifest skip mismatch rebuild

Run id: `retention_followup_mismatch`

After first archive, added new files in same group and reran retention:

```bash
python -m post_docking.artifact_retention --run-id retention_followup_mismatch --pdb-id PM1 --mode rerun_safe --threads 2
```

Observed warnings:
- `reason=manifest_mismatch ... action=rebuild`

Result: PASS

### Tool preflight skip-safe

Forced missing tools via `PATH=/tmp`:

```bash
PATH=/tmp /home/michael/atlas/micromamba/envs/docking-env/bin/python -m post_docking.artifact_retention --run-id retention_followup_dryrun --pdb-id PDBE --mode rerun_safe
```

Observed:
- `reason=tooling_preflight_failed detail=tar_missing`
- return code `0`, no file deletion

Result: PASS

### Restore hardening path guard

Tampered manifest target outside allowed roots:

```bash
python -m post_docking.artifact_retention \
  --restore \
  --run-id retention_followup_default \
  --archive docked/retention_followup_default/_artifact_archives/docked/PDBA/HOLO/pH7_0/stage1/unknown/artifacts.tar.zst \
  --manifest docked/retention_followup_default/_artifact_archives/docked/PDBA/HOLO/pH7_0/stage1/unknown/tampered.manifest.json
```

Observed:
- `reason=outside_allowed_roots`
- `failed=1`, `restored=0`, exit code `1`

Result: PASS

### CLI edge case `-test` / `--test`

```bash
python - <<'PY'
from src.cli.cli_utils import _parse_specified_proteins
print(_parse_specified_proteins(["main.py", "--test"], {}))
print(_parse_specified_proteins(["main.py", "-test"], {}))
PY
```

Output:
- `([], '')`
- `([], '')`

Result: PASS

## Quality Gates

### Ruff

Repository-wide:

```bash
ruff check --fix .
```

Result:
- fails with pre-existing unrelated repository errors (DeepCoy/protein_prep).

Focused touched-file lint:

```bash
ruff check src/post_docking/artifact_retention.py postrun_hooks.py main.py src/cli/cli_utils.py --ignore E402
```

Result:
- PASS

### MyPy (touched files)

```bash
mypy --ignore-missing-imports src/post_docking/artifact_retention.py
mypy --ignore-missing-imports postrun_hooks.py
mypy --ignore-missing-imports --follow-imports=skip main.py src/cli/cli_utils.py
```

Result:
- PASS

## Git Diff Summary

Executed:

```bash
git diff --stat
```

Output:
- `main.py`: retention mode precedence/flags, per-protein invocation, no run-end primary call
- `postrun_hooks.py`: new per-pdb retention invoker + lock support + simplified mode handling
- `src/cli/cli_utils.py`: `-test/--test` PDB parsing fix
- new docs + new module `src/post_docking/artifact_retention.py`
