# Atlas Status Dashboard

`atlas status` is the operational dashboard for active and completed runs.

```bash
atlas status <RUN_ID>
atlas status <RUN_ID> --errors --explain
atlas status <RUN_ID> --watch 30
atlas status <RUN_ID> --json
atlas status <RUN_ID> --html
atlas runs
```

The dashboard reads only run artifacts. It does not mutate docking outputs.

Reported sections:

- run lifecycle from `outputs/manifests/<RUN_ID>/run_manifest.yaml`
- protein completion counts and examples of running/failed entries
- per-stage progress for `prep`, `pocket_detection`, `docking`, and `postprocessing`
- distributed chunk plan/result/claim counts when combo chunking is active
- conservative ETA based on completed-run history, distributed chunks, or current completion rate
- Slurm state from `squeue` while active and `sacct` for accounting when available
- structured root-cause summaries with `--errors`
- key output paths, next commands, and plain-language explanations with `--explain`
- standalone browser output with `--html` or `--html path/to/status.html`

Slurm integration is optional. On machines without Slurm, the dashboard prints
`slurm: unavailable` and continues using manifest and filesystem evidence.

Error summaries are intentionally grouped into broad operational buckets such as
`missing_tool_or_config`, `receptor_prep`, `ligand_prep`, `docking_engine`,
`pose_validation`, `scheduler_or_timeout`, `disk_path_io`, and
`report_generation`.

Root-cause bucket rules live in `analysis/config/root_cause_rules.yaml`, so new
failure signatures can be added without changing Python code.

The JSON output contract is documented in `docs/status_dashboard_schema.md` and
`docs/schemas/status_dashboard.schema.json`. `atlas runs` lists recent run IDs,
progress, failures, age, and whether an HTML report exists.

The HTML status page is generated from the same in-memory status object as
`atlas status --json`; it is a static file and does not require a web server.
