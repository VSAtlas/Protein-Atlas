## Schema: run_manifest

The run manifest is the reproducibility ledger for a pipeline run. It should be enough to identify what was run, which config snapshot was used, and which per-target outputs are expected.

| Field | Type | Required | Meaning | Example | Source |
| --- | --- | --- | --- | --- | --- |
| `run_id` | string | yes | unique run identifier | `20260428_210000` | runtime lifecycle |
| `argv` | array[string] | yes | CLI arguments used | `["--pdb","1M17","--fast"]` | CLI capture |
| `config_snapshot_path` | string | yes | path to frozen run config | `outputs/configs/20260428_210000/config.txt` | run bootstrap |
| `start_time` | string | yes | run start timestamp | `2026-04-28T21:00:00Z` | runtime logging |
| `end_time` | string | no | run end timestamp | `2026-04-28T21:18:22Z` | runtime logging |
| `status` | string | yes | overall status | `completed` | run lifecycle |
| `proteins` | object | yes | per-protein statuses and metadata | `{...}` | manifest runtime writers |
| `input_roots` | object | no | resolved input/runtime roots | `{"input_pdbs":"input_pdbs"}` | config snapshot |
| `tool_versions` | object | no | detected external tool versions | `{"vina":"1.2.5"}` | tool verification |
| `library_map` | object | no | resolved ligand-library mapping | `{"bNJS":"bench_pur2"}` | run profile/config |
| `schema_version` | string | no | manifest schema version | `0.1.0-alpha` | run lifecycle |

Do not store machine-specific secrets in manifests. Absolute local paths may appear in private run ledgers, but public artifact bundles should include enough relative-path context to reproduce without exposing local workstation layout.
