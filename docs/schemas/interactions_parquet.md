## Schema: interactions_parquet

Parquet export for per-pose interaction evidence used by analysis/reporting.

| Field | Type | Required | Meaning | Example | Source |
| --- | --- | --- | --- | --- | --- |
| `run_id` | string | yes | run identifier | `20260428_210000` | runtime context |
| `pdb_id` | string | yes | protein structure ID | `1M17` | docking result |
| `drug_id` | string | yes | ligand identifier | `imatinib` | docking result |
| `pose_rank` | int | yes | pose index/rank | `1` | docking output |
| `interaction_type` | string | yes | interaction class | `hbond` | interaction parser |
| `residue_id` | string | yes | contacted residue | `MET793` | interaction parser |
| `distance_angstrom` | float | no | interaction distance | `2.9` | interaction parser |
| `score_component` | float | no | optional component score | `-0.4` | engine-specific parser |
| `target_id` | string | no | target identifier for report joins | `EGFR` | target mapping |
| `uniprot_id` | string | no | UniProt accession for report joins | `P00533` | target mapping |
| `atlas_score` | float | no | joined plausibility score when exported with report context | `2.41` | report scoring |
| `z_selected` | float | no | joined canonical selected score | `-1.73` | docking summaries |
| `mmgbsa_score` | float | no | joined MM/GBSA support score | `-31.2` | MM/GBSA outputs |
| `ligand_chemotype` | string | no | ligand chemotype family | `anilide` | ligand annotations |

Interaction parquet files are evidence tables, not the canonical heatmap matrix. Report joins may add optional score and annotation fields, but interaction identity fields remain required.
