## Schema: heatmap_input

Heatmap input is the report-facing matrix source. Each row represents one displayed drug-target/structure value plus optional biological context.

| Field | Type | Required | Meaning | Example | Source |
| --- | --- | --- | --- | --- | --- |
| `drug_id` | string | yes | heatmap row key | `imatinib` | reporting join |
| `target_id` | string | yes | heatmap column key | `EGFR` | reporting join |
| `pdb_id` | string | yes | structure backing the score | `1M17` | docking summary |
| `z_selected` | float | yes | canonical numeric value for heatmap fill | `-1.73` | scoring selection |
| `atlas_score` | float | no | optional tooltip/rank metric | `2.41` | report scoring |
| `mmgbsa_score` | float | no | optional tooltip metric | `-31.2` | post-docking |
| `free_cmax` | float | no | ligand exposure context | `0.09` | annotation joins |
| `exposure_plausibility` | string | no | exposure support category | `moderate` | report annotations |
| `tissue_expression` | string | no | target context annotation | `high:liver` | target cache |
| `target_adr_evidence` | string | no | target-level ADR evidence | `literature+label` | evidence join |
| `pathway_link_to_side_effect` | string | no | ADR pathway annotation | `MAPK signaling` | pathway resolver |
| `structure_resolution` | float | no | structure resolution (A) | `2.1` | structure metadata |
| `ligand_chemotype` | string | no | chemotype family | `anilide` | ligand annotations |
| `literature_supported_label` | bool | no | literature support flag | `true` | curated evidence |

Heatmap generation should prefer `z_selected` when multiple score columns are present. Missing optional evidence fields should be rendered as unknown or absent context, not as negative evidence.
