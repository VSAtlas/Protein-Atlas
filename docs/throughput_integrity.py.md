# throughput_integrity.py changes

## What was fixed

1. Prevent false passes when expected ligand set cannot be resolved.
2. Align ligand canonicalization with production `canonical_ligand_base` logic (with fallback).
3. Restrict known-failure attribution to stage-scoped markers and run-window timestamps.
4. Add required protein-prep artifact parity checks per `(pdb, variant, ph)` combo.
5. Remove unused helper code and expand failure-reason reporting in JSON/CSV outputs.

## Why this matters

These changes make throughput benchmarks harder to game accidentally: strict mode now fails when ligand coverage cannot be proven or when receptor prep artifacts are missing, while reducing false attribution from stale failure markers.
