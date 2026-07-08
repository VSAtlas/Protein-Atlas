# tools/architecture_gate.py

Runs the repo-local architecture gates for maintained source code:

```bash
python tools/architecture_gate.py
```

The command runs `tach check` using `tach.toml`, then runs Xenon against `src/`
with `--max-absolute B --max-modules B --max-average A`.

Current pre-existing Xenon violations are listed in `tools/xenon_baseline.txt`.
Remove files from that baseline as complexity is reduced; newly introduced
high-complexity files outside the baseline will fail the gate.

Use `--skip-tach` or `--skip-xenon` when isolating one gate during debugging.
