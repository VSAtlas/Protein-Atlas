## Schema: master_rows

Master rows are report-facing merged rows joining docking, target, ligand, and evidence annotations.
They are the preferred long-form table for publication-facing side-effect/ADR analysis.

| Field | Type | Required | Meaning | Example | Source |
| --- | --- | --- | --- | --- | --- |
| `drug_id` | string | yes | ligand identifier | `imatinib` | ligand table |
| `target_id` | string | yes | target identifier | `EGFR` | target table |
| `uniprot_id` | string | no | UniProt accession | `P00533` | target mapping |
| `pdb_id` | string | yes | structure ID | `1M17` | run outputs |
| `atlas_score` | float | no | summary ranking signal | `2.41` | report scoring |
| `z_selected` | float | yes | canonical selected score | `-1.73` | docking summaries |
| `mmgbsa_score` | float | no | MM/GBSA evidence | `-31.2` | MM/GBSA outputs |
| `free_cmax` | float | no | free Cmax (uM) | `0.09` | exposure annotations |
| `exposure_plausibility` | string | no | exposure support category | `moderate` | report annotations |
| `tissue_expression` | string | no | expression summary | `high:liver` | target expression cache |
| `target_adr_evidence` | string | no | target-level ADR evidence | `literature+label` | evidence join |
| `pathway_link_to_side_effect` | string | no | pathway rationale | `MAPK signaling` | pathway resolver |
| `structure_resolution` | float | no | structure resolution (A) | `2.1` | structure metadata |
| `ligand_chemotype` | string | no | chemotype family | `anilide` | ligand annotations |
| `literature_supported_label` | bool | no | literature support flag | `true` | curated evidence |

`z_selected` remains the canonical selected score. `atlas_score` is a relative plausibility/ranking signal and must not be described as an experimental affinity.
