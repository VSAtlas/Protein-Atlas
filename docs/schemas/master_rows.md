## Schema: master_rows

Master rows are report-facing merged rows joining docking, target, ligand, and evidence annotations.
They are the preferred long-form table for publication-facing side-effect/ADR analysis.

| Field | Type | Required | Meaning | Example | Source |
| --- | --- | --- | --- | --- | --- |
| `drug_id` | string | yes | ligand identifier | `imatinib` | ligand table |
| `target_id` | string | yes | target identifier | `EGFR` | target table |
| `uniprot_id` | string | no | UniProt accession | `P00533` | target mapping |
| `pdb_id` | string | yes | structure ID | `1M17` | run outputs |
| `fda_mapping_csv` | string | no | repo-relative mapping path or redacted `external/<basename>` | `chemdb/data/fda_mapping_from_pdbqt.csv` | mapping provenance |
| `fda_mapping_sha256` | string | no | exact SHA-256 of the FDA identity mapping | `4697ca38...` | mapping provenance |
| `atlas_score` | float | no | summary ranking signal | `2.41` | report scoring |
| `z_selected` | float | no | canonical decoy-standardized selected score | `-1.73` | docking summaries |
| `z_selected_source` | string | no | selected-score provenance or missing-null reason | `consensus_z_reconstructed_from_decoys` | report scoring |
| `consensus_score` | float | no | raw within-library consensus rank percentile | `0.84` | docking summaries |
| `mmgbsa_score` | float | no | MM/GBSA evidence | `-31.2` | MM/GBSA outputs |
| `free_cmax` | float | no | free Cmax (uM) | `0.09` | exposure annotations |
| `exposure_plausibility` | string | no | exposure support category | `moderate` | report annotations |
| `tissue_expression` | string | no | expression summary | `high:liver` | target expression cache |
| `target_adr_evidence` | string | no | target-level ADR evidence | `literature+label` | evidence join |
| `pathway_link_to_side_effect` | string | no | pathway rationale | `MAPK signaling` | pathway resolver |
| `structure_resolution` | float | no | structure resolution (A) | `2.1` | structure metadata |
| `ligand_chemotype` | string | no | chemotype family | `anilide` | ligand annotations |
| `literature_supported_label` | bool | no | literature support flag | `true` | curated evidence |

`z_selected` remains the canonical selected score. It must contain only a decoy-standardized score; raw `consensus_score` values must not be substituted when the decoy null is unavailable. In that case `z_selected` remains missing and `z_selected_source` records `missing_consensus_decoy_null`. `atlas_score` is a relative plausibility/ranking signal and must not be described as an experimental affinity.
