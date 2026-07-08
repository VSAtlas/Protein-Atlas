## Schema: docking_score_summary

Canonical score-summary rows describe one ligand-target/structure score selected for downstream reporting.
`z_selected` is the canonical selected-score column for heatmaps and publication-facing summaries.

| Field | Type | Required | Meaning | Example | Source |
| --- | --- | --- | --- | --- | --- |
| `drug_id` | string | yes | ligand/drug identifier | `imatinib` | ligand metadata/join keys |
| `target_id` | string | yes | target label used by report | `EGFR` | target mapping |
| `uniprot_id` | string | no | UniProt accession | `P00533` | target annotations |
| `pdb_id` | string | yes | protein structure ID | `1M17` | run selection |
| `atlas_score` | float | no | composite plausibility score | `2.41` | report layer |
| `z_selected` | float | yes | canonical selected z-score | `-1.73` | score selection policy |
| `mmgbsa_score` | float | no | MM/GBSA support score | `-31.2` | post-docking mmgbsa |
| `free_cmax` | float | no | free Cmax estimate (uM) | `0.09` | ligand exposure annotations |
| `exposure_plausibility` | string | no | exposure plausibility bucket | `moderate` | reporting annotation |
| `tissue_expression` | string | no | target tissue expression summary | `high:liver` | target expression cache |
| `target_adr_evidence` | string | no | ADR evidence summary at target level | `literature+label` | report evidence join |
| `pathway_link_to_side_effect` | string | no | pathway rationale for ADR link | `MAPK signaling` | pathway resolver |
| `structure_resolution` | float | no | experimental structure resolution (A) | `2.1` | PDB metadata |
| `ligand_chemotype` | string | no | ligand chemotype family | `anilide` | ligand annotations |
| `literature_supported_label` | bool | no | literature support flag | `true` | curation/evidence |

Schema changes to score meaning, score preference, or required keys are scientific-method changes and should be reviewed with reproducibility notes.
