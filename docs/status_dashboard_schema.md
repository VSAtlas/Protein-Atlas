# Status Dashboard JSON Schema

`atlas status <RUN_ID> --json` emits a stable dashboard document described by:

```text
docs/schemas/status_dashboard.schema.json
```

The top-level object includes:

- `lifecycle`: manifest status and timing.
- `progress`: protein-level scheduled/completed/failed/running counts.
- `stages`: per-stage prep, pocket detection, docking, and postprocessing state.
- `distributed`: chunk plan, result, claim, and phase-marker counts.
- `slurm`: live `squeue` or accounting `sacct` rows when a Slurm job id is known.
- `eta`: current ETA estimate and whether it came from run history, chunks, or current completion rate.
- `eta_history`: compact summary of completed-run history used for calibration.
- `errors`: structured root-cause summary when `--errors` or `--deep` is enabled.
- `paths` and `existing_paths`: expected and currently present artifacts.
- `next` and `explanations`: suggested operator actions.
- `html_path`: present only when `atlas status --html` writes a browser dashboard.

The schema is intentionally permissive with `additionalProperties` so future
dashboards can add fields without breaking consumers that validate the current
contract.
