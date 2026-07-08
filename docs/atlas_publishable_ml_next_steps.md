Date: 2026-05-05

# Atlas Publishable ML Readiness Plan

This note records the current pilot ML status and the next engineering steps needed before making manuscript-strength ML claims.

## Current Label Sources

The pilot now stages local assay-backed activity labels from:

- SPD: exposure-relevant secondary pharmacology benchmark.
- Papyrus/ChEMBL: broader measured drug-target bioactivity.
- ToxCast/invitrodb v4.3: independent in vitro assay activity benchmark.

ToxCast must be described as an in vitro activity benchmark, not clinical ADR validation. ADR claims should continue to use SIDER, Open Targets Safety, CTD, Reactome, and curated literature labels as mechanism/evaluation layers.

## Implemented ToxCast Filters

The ToxCast source table is now filtered before label generation:

- Assay filter: `export_ready == 1`, `data_usability == 1`, and `cell_viability_assay != 1`.
- Chemical QC filter: analytical QC `pass_or_caution` is `pass` or `caution`.
- Cytotoxicity filter: active hits with `ac50_uM >= cytotox_lower_bound_um` are excluded as likely cytotoxicity-confounded, not converted to inactive.
- DSSTox enrichment: staged ToxCast rows include `toxcast_dsstox_substance_id` where available, with pilot drug-level DSSTox mappings recorded in the extraction summary.

## Exposure Policy

Free-Cmax remains sparse and is explicitly preserved as unknown where missing. Missing exposure must not be treated as negative evidence. Exposure-aware analyses should report the number of rows with nonmissing `free_cmax_um` and should either:

- run an exposure-available subset analysis, or
- include explicit missingness indicators only in sensitivity analyses, not in the primary transparent score.

## Publishable ML Protocol

Before claiming ML performance, run all of the following:

- Deduplicate to one row per drug-target label source decision.
- Report label source provenance for every row.
- Evaluate random split only as a sanity check.
- Treat drug-holdout, target-holdout, and scaffold-holdout as required validation.
- Add source-holdout validation: train on ChEMBL/Papyrus and test on SPD/ToxCast where overlap permits.
- Add time-split validation once evidence/publication or database release years are available.
- Calibrate predicted probabilities and report Brier score plus reliability diagrams.
- Report AUROC, AUPRC, EF@1%, EF@5%, EF@10%, precision@K, and bootstrap confidence intervals.
- Compare against non-ML baselines: atlas score only, consensus/Vina score only, SCORCH subset, priority score, random, and label-shuffled controls.
- Add ablations that remove exposure, assay-source indicators, target evidence, pathway evidence, and mechanism graph score.
- Keep raw `drug_id`, `target_id`, `pdb_id`, side-effect labels, and literature labels out of predictive features.

## Remaining ML Gaps

- Final FDA/RDK mapping should be rerun before freezing the ML table.
- DSSTox matching can still improve if an EPA CTX API key is available for InChIKey or PubChem-CID batch resolution.
- ToxCast assay filtering is conservative but still basic; target-family-specific assay review may be needed for publication panels.
- Free-Cmax needs a dedicated pharmacokinetic source table or curated subset before exposure-aware ML claims.
- Temporal validation needs year columns such as `activity_publication_year`, `database_release_year`, or `label_publication_year`.
- Positive-unlabeled handling is still needed for unlabeled Atlas pairs; no-assay rows should remain unknown.

## Source Transfer And Temporal Validation Update

Source-holdout validation should be interpreted as a domain-transfer stress test, not as the primary headline metric. ChEMBL/Papyrus labels are literature-curated bioactivity records with strong target and publication bias, while ToxCast/invitrodb is a heterogeneous high-throughput in vitro assay resource with a much lower positive rate in the pilot. The training and test source prevalence must therefore be reported alongside AUROC/AUPRC.

The source-holdout trainer now writes:

- `source_transfer_label_balance.csv`
- `source_transfer_prediction_balance.csv`

These files show whether a weak source-transfer result reflects score failure, source label imbalance, or probability miscalibration.

Temporal validation should prefer true evidence dates over database release dates:

1. Use `activity_publication_year` from ChEMBL `document_year` and Papyrus `Year` where available.
2. Exclude rows with missing publication years from publication-year temporal splits rather than assigning them to train or test.
3. Use `database_release_year` only as a coarser release-snapshot sensitivity analysis.
4. Treat ToxCast v4.3 rows as release-year dated unless a comparable per-assay publication/date field is staged.

Current pilot caveat: ToxCast dominates the labeled table but lacks comparable activity-publication years in the local summary extract, so publication-year temporal validation is currently a ChEMBL/Papyrus-backed subset analysis.

## Negative Label Hierarchy

Atlas keeps negative/background concepts separate. The detailed policy is in `docs/atlas_ml_negative_label_strategy.md` and should be treated as the source of truth for mechanism/ADR labels.

Mechanism ML uses a four-state label table:

| value | meaning |
|---:|---|
| `1` | strict positive mechanism triad |
| `0` | measured or reliable negative |
| `NaN` | unknown |
| `-1` | excluded, ambiguous, or conflicting |

Assign `0` only from measured inactive assay pairs, OMOP/OHDSI drug-outcome negative controls, adequately observed FAERS non-signals with no independent positive evidence, or fold-specific reliable negatives mined with external similarity features only. Do not use Atlas score, SCORCH, FDR, heatmap rank, or `mechanism_graph_score` to choose negatives if those quantities are model features.

Decoys remain benchmark-only:

1. Measured assay inactives are the primary supervised negatives for binding/exposure tasks.
2. OMOP/OHDSI and FAERS non-signals are calibration or sensitivity controls for drug-ADR benchmarks.
3. Unlabeled FDA-target pairs may be included only as positive-unlabeled/background rows with low sample weight.
4. Reliable negatives may be mined only inside each training fold using external similarity features.
5. DUD-E/DEKOIS/DeepCoy/diversity-library decoys are technical docking benchmark controls, not FDA/off-target ML labels.

The current main pilot model uses measured actives/inactives plus PU-weighted FDA background:

```bash
python -m analysis.cli.build_ml_dataset \
  --pair-table data/pilotstudy/bioactivity_ml_source.csv \
  --label bioactivity_ml_label \
  --feature-set pilot_nonleaky \
  --out data/pilotstudy/ml_pair_table_bioactivity_pu_weighted.csv \
  --deduplicate-drug-target \
  --include-unlabeled-as-background \
  --unlabeled-weight 0.25

python -m analysis.cli.train_ml_model \
  --dataset data/pilotstudy/ml_pair_table_bioactivity_pu_weighted.csv \
  --label bioactivity_ml_label \
  --feature-set pilot_nonleaky \
  --model logistic_regression \
  --split drug_holdout \
  --pu-mode positive_unlabeled_weighted \
  --out-dir data/pilotstudy/ml_bioactivity_logistic_main_pu_drug_holdout
```

The separated decoy benchmark panel is run independently:

```bash
python -m analysis.cli.run_decoy_benchmark \
  --input data/pilotstudy/master_rows.csv \
  --out-dir data/pilotstudy/decoy_benchmark_panel \
  --scores z_selected final_score SCORCH_score_used z_vs_decoys_blend consensus_score
```

Decoy benchmark outputs:

- `decoy_benchmark_summary.csv`
- `decoy_benchmark_by_target.csv`
- `decoy_benchmark_rows.csv`
- `decoy_benchmark_manifest.json`
- `decoy_benchmark_score_distributions.png`

Reusable negative-evidence ingestion now lives behind:

```bash
python -m analysis.cli.build_negative_evidence
python -m analysis.cli.build_mechanism_label_table
python -m analysis.cli.build_reliable_negatives
```

These commands normalize measured inactive sources, OMOP/OHDSI controls, FAERS observed non-signals, and fold-local reliable negatives into the four-state mechanism policy. Drug-target measured negatives can label pair rows directly. Drug-ADR controls only label rows when an ADR key is present; otherwise they remain calibration/control evidence and are not forced onto target-pair mechanism labels.

Current staged pilot outputs:

- `data/pilotstudy/negative_evidence/negative_evidence.csv`
- `data/pilotstudy/negative_evidence/negative_evidence.summary.json`
- `data/pilotstudy/negative_evidence/mechanism_four_state_labels.csv`
- `data/pilotstudy/negative_evidence/negative_pair_mapping_audit.csv`
- `data/pilotstudy/negative_evidence/negative_pair_mapping_audit_summary.csv`
- `data/pilotstudy/negative_evidence/mechanism_four_state_labels_source_balanced.csv`
- `data/external/faers/ohdsi_omop_faers_disproportionality.csv`
- `data/pilotstudy/moe_experts/moe_expert_performance_summary.csv`
- `data/pilotstudy/moe_experts/moe_expert_label_source_summary.csv`
- `data/pilotstudy/moe_experts/mechanism_label_composition_summary.csv`
- `data/pilotstudy/moe_experts/moe_adr_prioritization_with_mechanism_predictions.csv`

Current negative evidence staging:

| evidence type | rows | role |
|---|---:|---|
| measured inactive assay | 6,210 | pair-level measured/reliable negatives |
| OMOP/OHDSI negative control | 510 | drug-ADR calibration/control evidence |
| FAERS observed non-signal | 117 | weak drug-ADR calibration/control evidence |

The current source-transfer pair table has no ADR key, so OMOP/OHDSI and FAERS controls are not forced onto pair-level mechanism labels. They remain staged for future ADR-keyed calibration/control analysis.

Current four-state mechanism labels:

| mechanism label | rows | interpretation |
|---|---:|---|
| `1` | 16 | strict positive mechanism triad |
| `0` | 6,157 | measured/reliable negative |
| `NaN` | 792 | unknown |
| `-1` | 53 | excluded ambiguous/conflicting |

Negative source mix among the 6,157 pair-level negatives:

| negative source | rows | percent of negatives |
|---|---:|---:|
| ToxCast HTS inactive | 5,836 | 94.786 |
| curated ChEMBL/Papyrus inactive | 321 | 5.214 |

This is intentionally audited because 94.786% ToxCast negatives is too source-dominated for a single headline mechanism classifier. Use the all-negative model as a sensitivity check and the source-balanced dataset when asking whether the signal survives without ToxCast overwhelming the negatives.

Pair-mapping audit:

| check | passed | failed |
|---|---:|---:|
| negative rows joined to source pair | 6,210 | 0 |
| inactive label verified | 6,210 | 0 |
| target metadata present | 6,210 | 0 |
| assay IDs present | 6,210 | 0 |
| unique negative drug-target pairs | 6,210 | 0 |

Current three-expert monitoring table:

| expert | objective | test positives | test negatives | AUROC | AUPRC |
|---|---|---:|---:|---:|---:|
| bioactivity measured activity | measured in-vitro activity/inactivity | 161 | 1,223 | 0.739 | 0.314 |
| SPD exposure relevance | AC50/free-Cmax exposure relevance | 30 | 70 | 0.601 | 0.339 |
| mechanism sensitivity | strict triads vs measured inactive negatives | 1 | 1,247 | 0.709 | 0.003 |
| mechanism source-balanced sensitivity | strict triads vs source-balanced inactive negatives | 6 | 143 | 0.752 | 0.287 |

The mechanism sensitivity model is not publication-ready: the drug-holdout test set has only one strict positive and AUPRC is therefore nearly base-rate limited. Treat it as a pipeline sanity check, not a positive mechanism ML claim.

## Current Pilot ML Sensitivity Outputs

The following pilot analyses are now materialized under `data/pilotstudy/`:

- `ml_bioactivity_logistic_measured_only_drug_holdout/`: measured assay actives/inactives only, no PU background.
- `ml_bioactivity_logistic_main_pu_drug_holdout/`: measured assay labels plus low-weight unlabeled FDA background.
- `ml_bioactivity_logistic_chembl_document_year_temporal/`: ChEMBL-only temporal validation using `document_year` propagated into `activity_publication_year`.
- `ml_bioactivity_source_calibration_toxcast/`: ChEMBL/Papyrus train, ToxCast calibration split, held-out ToxCast test split.
- `exposure_coverage/`: free-Cmax, total Cmax, and fraction-unbound coverage report.

Current pilot metrics:

- Measured-only drug holdout: AUROC `0.733`, AUPRC `0.427`.
- Main PU drug holdout: AUROC `0.682`, AUPRC `0.026`.
- ChEMBL-only document-year temporal holdout: AUROC `0.531`, AUPRC `0.371`.
- Raw ChEMBL/Papyrus to ToxCast transfer: AUROC `0.433`, AUPRC `0.052`, Brier `0.348`, ECE `0.531`.
- ToxCast-split logistic calibration sensitivity: AUROC `0.567`, AUPRC `0.075`, Brier `0.057`, ECE `0.002`.

The ToxCast calibration run is not pure external source transfer because it uses a ToxCast calibration split. It should be described as source-specific calibration/sensitivity analysis, not as independent external validation.

Exposure coverage in the current pilot `bioactivity_ml_source.csv`:

- `free_cmax_um`: 63,833 / 176,682 rows, 36.1%.
- `cmax_um`: 28,908 / 176,682 rows, 16.4%.
- `fraction_unbound_plasma`: 27,360 / 176,682 rows, 15.5%.

Missing exposure remains unknown and is not converted to negative evidence.

## ChEMBL Temporal Drift Diagnostics

Additional diagnostics are staged under `data/pilotstudy/chembl_temporal_diagnostics/`.

The supervised ChEMBL table does not contain DUD/decoy rows. The pilot builder filters `is_decoy` rows before bioactivity labels are joined, so DUD/decoy rows currently support FDR and the separate decoy benchmark only.

ChEMBL-only model checks:

- Random holdout: AUROC `0.533`, AUPRC `0.669`.
- Drug holdout: AUROC `0.533`, AUPRC `0.655`.
- Target holdout: AUROC `0.523`, AUPRC `0.747`.
- Scaffold holdout: AUROC `0.447`, AUPRC `0.596`.
- Document-year temporal holdout: AUROC `0.531`, AUPRC `0.371`.

Interpretation: the temporal split is not uniquely failing by AUROC. ChEMBL-only predictive signal is weak across split types; AUPRC changes largely track prevalence and label distribution.

Temporal label shift:

- Train pre-2016: 304 positives / 68 negatives.
- Test post-2015: 110 positives / 216 negatives.

Balanced resampling of the temporal test predictions gives median AUROC `0.531` and median AUPRC `0.533`, confirming that the low temporal AUPRC is partly prevalence-driven.

Single-score baselines on the same temporal test set:

- `atlas_score`: AUROC `0.489`, AUPRC `0.395`.
- `consensus_score`: AUROC `0.509`, AUPRC `0.373`.
- `free_cmax_um`: AUROC `0.298`, AUPRC `0.257` on rows with exposure.
- `fraction_unbound_plasma`: AUROC `0.601`, AUPRC `0.510` on sparse rows.

Assay metadata drift exists and should be reported as a limitation:

- Train rows are enriched for IC50/Ki and exact `=` relations.
- Test rows are enriched for Kd and censored `>` relations.
- Median ChEMBL assay count shifts from 4 in train to 1 in test.

Assay metadata was not added as a main predictive feature because it is label-generation metadata and is not available for new FDA-target pairs before the pair is assayed. It is retained for drift diagnostics only.

## Reproducible ML Sensitivity Flags

Current ML CLI flags support the publication sensitivity settings directly:

- Exposure-complete training: add `--require-nonmissing-cols free_cmax_um` when building the ML dataset.
- No-PK sensitivity training: add `--exclude-features free_cmax_um cmax_um fraction_unbound_plasma` when building or training.
- Source-specific training: add `--label-source-include chembl papyrus toxcast` when building the ML dataset.
- Temporal evidence filtering: add `--require-year-col activity_publication_year` when building and `--split temporal --temporal-year-col activity_publication_year` when training.

Pilot model panel generated on 2026-05-06:

| model | test rows | test positive rate | AUROC | AUPRC | EF@1% | EF@5% | Brier | ECE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| all measured actives/inactives | 1406 | 0.120 | 0.796 | 0.383 | 4.437 | 4.570 | 0.213 | 0.332 |
| exposure-complete only | 517 | 0.180 | 0.781 | 0.453 | 3.706 | 3.421 | 0.212 | 0.285 |
| SPD exposure relevance | 100 | 0.300 | 0.897 | 0.889 | 3.333 | 3.333 | 0.110 | 0.163 |
| no-PK sensitivity | 1406 | 0.120 | 0.791 | 0.374 | 4.437 | 4.453 | 0.215 | 0.338 |

Interpretation: the exposure-complete subset improves AUPRC and calibration but has fewer drugs/targets and a different prevalence. The no-PK sensitivity model is close to the all-feature model, so current pilot data do not yet prove free-Cmax is the main driver. The SPD result is promising but small and should be rerun after the full SPD protein panel is processed through Atlas.

## Structured PK Source Builder

`python -m analysis.cli.build_pk_table` now builds a local structured PK table with standard columns:

- `drug_id`
- `drug_name`
- `free_cmax_um`
- `cmax_um`
- `fraction_unbound_plasma`
- `exposure_source`
- `pk_source`

Supported local inputs:

- `--existing-pk`
- `--spd-pk`
- `--openfda-pk`
- `--pkdb`
- `--drugbank-cmax`
- `--drugbank-protein-binding`

The builder does not download licensed DrugBank tables. DrugBank structured exports must be supplied locally if available. Missing exposure remains unknown.

## Adversarial ML Strategy Audit

Additional audit outputs are staged under `data/pilotstudy/ml_audit_*`.

Fixed pipeline-level loopholes:

- Numeric missingness is now explicit in model design matrices via `__missing` indicator features instead of being hidden by median imputation alone.
- `build_ml_dataset` can drop rows with flagged source conflicts using `--exclude-flag-cols`.
- SPD leakage checks now forbid direct exposure-label ingredients such as `exposure_margin`, `ac50_nM`, `free_cmax_nM`, and `total_cmax_nM` as predictive features.
- A stricter `pilot_binding_only` feature set is available for primary publication models that avoid PK provenance/missingness effects.

Audit metrics:

| model | test rows | test positive rate | AUROC | AUPRC | EF@1% | EF@5% | Brier | ECE |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| binding-only drug holdout | 1406 | 0.120 | 0.791 | 0.374 | 4.437 | 4.453 | 0.215 | 0.338 |
| binding-only target holdout | 1987 | 0.084 | 0.808 | 0.282 | 5.354 | 3.688 | 0.192 | 0.338 |
| full with missingness indicators, drug holdout | 1406 | 0.120 | 0.823 | 0.379 | 3.882 | 3.984 | 0.197 | 0.314 |
| full target holdout | 1987 | 0.084 | 0.811 | 0.276 | 4.164 | 4.045 | 0.190 | 0.328 |
| no source-disagreement rows | 1397 | 0.115 | 0.816 | 0.357 | 4.366 | 3.867 | 0.198 | 0.319 |
| SPD exposure relevance, no PK features | 100 | 0.300 | 0.601 | 0.339 | 0.000 | 0.000 | 0.245 | 0.195 |
| ChEMBL/Papyrus to ToxCast source holdout | 6056 | 0.025 | 0.507 | 0.026 | 0.649 | 0.131 | 0.343 | 0.555 |

Publication interpretation after audit:

- The most defensible primary ML claim is currently binding-score based pair ranking under drug/target/scaffold holdouts.
- Exposure should be reported as stratification/external validation until the full SPD panel is processed, because SPD labels are mathematically defined using free-Cmax and source coverage is non-random.
- ChEMBL/Papyrus-to-ToxCast transfer is not ready as a positive claim; the source-holdout result is near-random.
- The no-PK/binding-only model being close to the full model means the current pilot does not prove that free-Cmax drives model performance.

Remaining publication gaps:

- Full SPD protein panel has not yet been run through Atlas, so SPD ML validation is small and target-limited.
- Source imbalance remains severe: ToxCast dominates rows and is mostly negative; Papyrus/ChEMBL are much more positive.
- Target prevalence is uneven, with several targets nearly all-positive or all-negative.
- `ligand_chemotype` is only a fallback scaffold proxy; true Murcko scaffold grouping should be added when reliable SMILES are staged.
- Free-Cmax/protein-binding provenance is still incomplete and should be improved through `build_pk_table` plus curated PK-DB/DrugBank/SPD/local exports.
- Clinical ADR prediction should not be claimed from these assay labels alone; SIDER/OpenTargets/CTD/Reactome remain mechanism/evaluation layers, not binding labels.

## Separated Evaluation Targets

The pilot now writes explicit target tables under `data/pilotstudy/evaluation_targets/` with:

```bash
python -m analysis.cli.build_evaluation_targets \
  --out-dir data/pilotstudy/evaluation_targets \
  --pair-table data/pilotstudy/bioactivity_pair_table.csv \
  --bioactivity-table data/pilotstudy/bioactivity_ml_source.csv \
  --spd-table data/pilotstudy/spd_atlas_benchmark_ml_ready.csv
```

Generated pilot target tables:

- `bioactivity_assay_target.csv`: ChEMBL/Papyrus/ToxCast measured assay activity.
- `spd_exposure_relevance_target.csv`: SPD AC50/free-Cmax exposure relevance.
- `mechanism_adr_support_target.csv`: ADR mechanism support label scaffold.
- `evaluation_target_manifest.json`: target-specific feature policy, claims, and gaps.

Current pilot target coverage:

| target | rows | labelable | positives | negatives | primary feature policy |
|---|---:|---:|---:|---:|---|
| bioactivity assay activity | 176682 | 10465 | 1189 | 9276 | binding-only primary; PK sensitivity only |
| SPD exposure relevance | 498 | 498 | 170 | 328 | binding-only primary; PK is label-defining sensitivity |
| ADR mechanism support | 176682 | 0 | 0 | 0 | blocked until mechanism graph evidence is staged |

BigBind note: existing runtime BigBind utilities are already owned by `src/ml/data/bigbind.py` and related `src/ml/*` modules. Atlas analysis now has only a small local pair-table join shim, `analysis.external.bigbind`, for optional Atlas-compatible benchmark joins. BigBind is appropriate for docking/bioactivity evaluation, not exposure or ADR validation.

Revised target strategy:

1. Bioactivity/docking benchmark:
   - Prefer BigBind or source-specific ChEMBL/Papyrus/ToxCast as docking/activity evaluation.
   - Do not use PK as a primary feature for in-vitro assay activity.
   - Required publication checks: source holdout, target holdout, scaffold holdout, and source-prevalence reporting.

2. Exposure relevance:
   - Use SPD as the primary target.
   - Do not claim noncircular prediction from PK/free-Cmax because SPD labels are defined by AC50/free-Cmax.
   - Use PK/free-Cmax as exposure-context stratification and report binding-only SPD performance separately.

3. Mechanism/ADR support:
   - Build labels from mechanism graph evidence only after SIDER/OpenTargets Safety/CTD/Reactome are staged.
   - Missing mechanism evidence remains unknown, not negative.
   - Do not train on `literature_supported_label`, `drug_adr_known`, `target_adr_known`, `target_pathway_adr_link`, or `triad_complete` as predictive features.

Remaining ML publication gaps after target separation:

- Full SPD protein panel still needs to be docked/scored through Atlas.
- BigBind can strengthen the bioactivity/docking benchmark, but it is not yet joined to the pilot.
- Mechanism graph evidence is absent in the current pilot target table, so mechanism/ADR ML is not trainable yet.
- Source transfer remains weak; ChEMBL/Papyrus-to-ToxCast is not a positive claim.
- True Murcko scaffold splits should replace `ligand_chemotype` fallback once reliable SMILES are staged.

## SPD Panel Coverage and MoE Evidence Layer

The target builder now accepts the full SPD supplement:

```bash
python -m analysis.cli.build_evaluation_targets \
  --out-dir data/pilotstudy/evaluation_targets \
  --pair-table data/pilotstudy/bioactivity_pair_table.csv \
  --bioactivity-table data/pilotstudy/bioactivity_ml_source.csv \
  --spd-table data/pilotstudy/spd_atlas_benchmark_ml_ready.csv \
  --spd-full-panel data/external/spd/41467_2023_40064_MOESM4_ESM.xlsx
```

Current pilot SPD coverage after drug-name and target-gene matching:

| item | count |
|---|---:|
| SPD assayed drug-target pairs | 95512 |
| SPD unique drugs | 1948 |
| SPD unique targets | 101 |
| SPD pairs already represented in current Atlas pilot | 4453 |
| SPD pairs missing from current Atlas pilot | 91059 |
| SPD targets absent from current Atlas pilot | 90 |
| SPD drugs absent from current Atlas pilot | 1361 |

The generated run-list artifacts are under `data/pilotstudy/evaluation_targets/spd_full_panel/`. They are coverage/run-planning files, not a substitute for actually docking and rescoring the missing SPD panel.

Mechanism source staging is now wired into `build_external_labels`: SIDER, CTD, OpenTargets Safety, and Reactome paths are reported in `mechanism_source_stage_summary.csv`; if any local source files exist, the same command builds `mechanism_graph/pair_mechanism_scores.csv`.

A conservative mixture-of-evidence table is available with:

```bash
python -m analysis.cli.build_moe_adr_score \
  --pair-table data/pilotstudy/bioactivity_pair_table.csv \
  --bioactivity-predictions data/pilotstudy/ml_panel_main_all_measured_drug_holdout/model_predictions.csv \
  --exposure-predictions data/pilotstudy/ml_panel_spd_exposure_relevance_drug_holdout/model_predictions.csv \
  --mechanism-scores data/pilotstudy/evaluation_targets/mechanism_adr_support_target.csv \
  --out data/pilotstudy/evaluation_targets/moe_adr_prioritization.csv
```

Current pilot MoE coverage:

| expert availability | rows |
|---|---:|
| any MoE score | 2320 |
| bioactivity expert | 2103 |
| exposure expert | 235 |
| mechanism expert | 0 |

Interpretation: an MoE-style architecture is reasonable, but the current output is a conservative evidence-combination layer, not a trained clinical ADR probability. A publishable trained MoE needs calibrated held-out predictions for all three expert tasks: bioactivity, SPD exposure relevance, and mechanism/ADR support.

## Corrected Source-Transfer Audit

The ChEMBL/Papyrus-to-ToxCast transfer weakness is now audited with:

```bash
python -m analysis.cli.audit_source_transfer \
  --dataset data/pilotstudy/ml_pair_table_bioactivity_dedup_leakage_controlled.csv \
  --label bioactivity_ml_label \
  --train-sources chembl papyrus \
  --test-sources toxcast \
  --predictions data/pilotstudy/ml_audit_source_holdout_toxcast_corrected/model_predictions.csv \
  --out-dir data/pilotstudy/source_transfer_audit_toxcast_corrected
```

The label-fusion rule now treats source-level active/inactive disagreement as ambiguous instead of letting a positive source win. After rebuilding, cross-source disagreement rows are no longer labelable for supervised ML.

Corrected transfer diagnostics:

| item | value |
|---|---:|
| train rows | 959 |
| train positive rate | 0.660 |
| ToxCast test rows | 6021 |
| ToxCast test positive rate | 0.023 |
| drug overlap into ToxCast test | 0.082 |
| target overlap into ToxCast test | 1.000 |
| median atlas_score shift, IQR units | -0.700 |
| median consensus_score shift, IQR units | -0.701 |

Corrected pure source-holdout metrics:

| model | AUROC | AUPRC | Brier | ECE |
|---|---:|---:|---:|---:|
| ChEMBL/Papyrus -> ToxCast, no ToxCast calibration | 0.425 | 0.080 | 0.258 | 0.484 |

Corrected ToxCast source-calibration sensitivity:

| model | AUROC | AUPRC | Brier | ECE |
|---|---:|---:|---:|---:|
| ChEMBL/Papyrus train, ToxCast calibration split, ToxCast test | 0.648 | 0.140 | 0.023 | 0.000 |

Interpretation: the weak pure transfer is expected from severe prevalence shift, low drug overlap, lower docking-score distribution in ToxCast rows, and different assay semantics. The calibrated ToxCast sensitivity suggests some within-source ranking signal, but it is not pure external validation because it uses a ToxCast calibration split.

## Source-Specific Benchmark Tables

The fused bioactivity label is no longer the preferred benchmark for source-transfer claims. Build separated objectives with:

```bash
python -m analysis.cli.build_source_benchmark_tables \
  --bioactivity-source data/pilotstudy/bioactivity_ml_source.csv \
  --spd-table data/pilotstudy/spd_atlas_benchmark_dedup_drug_target_ml_annotated.csv \
  --out-dir data/pilotstudy/source_benchmark_tables
```

Outputs:

- `ml_curated_bioactivity_table.csv`: ChEMBL/Papyrus active/inactive objective.
- `ml_toxcast_hts_table.csv`: ToxCast HTS active/inactive objective after local ToxCast staging filters and conflict removal.
- `ml_spd_exposure_table.csv`: SPD AC50/free-Cmax exposure-relevance objective.
- `ml_source_transfer_table.csv`: internal table for curated-to-ToxCast stress tests.
- `ml_matched_source_transfer_table.csv`: source-transfer table restricted to target, chemotype, and score-range overlap.

Run the source-specific comparison with:

```bash
python -m analysis.cli.run_source_benchmark_comparison \
  --bioactivity-source data/pilotstudy/bioactivity_ml_source.csv \
  --spd-table data/pilotstudy/spd_atlas_benchmark_dedup_drug_target_ml_annotated.csv \
  --feature-set pilot_binding_only \
  --model logistic_regression \
  --split drug_holdout \
  --out-dir data/pilotstudy/source_benchmark_comparison
```

Current pilot comparison:

| analysis | AUROC | AUPRC | Brier | ECE |
|---|---:|---:|---:|---:|
| within-source curated bioactivity | 0.439 | 0.575 | 0.253 | 0.127 |
| within-source ToxCast HTS | 0.645 | 0.188 | 0.212 | 0.415 |
| within-source SPD exposure | 0.601 | 0.339 | 0.245 | 0.195 |
| curated -> ToxCast raw transfer | 0.652 | 0.146 | 0.246 | 0.472 |
| curated -> ToxCast source-calibrated | 0.636 | 0.139 | 0.023 | 0.000 |
| curated -> ToxCast matched transfer | 0.674 | 0.071 | 0.248 | 0.471 |

The separated objective improves the raw transfer audit versus the older fused-label runs, but the publication interpretation remains source-aware. ToxCast has about 2.4% positives after filtering, while curated ChEMBL/Papyrus remains about 67% positive. Report ToxCast as an independent HTS benchmark and calibrated sensitivity, not as proof of source-invariant bioactivity prediction.

The source-benchmark comparison now also reports `test_positive_rate`, `test_rows`, `mean_prediction`, `median_prediction`, and `AUPRC_over_prevalence`. These are required for imbalanced sources such as ToxCast, where AUPRC must be interpreted relative to the base positive rate.

The source-specific table builder now writes model-ready feature-restricted tables under:

```text
data/pilotstudy/source_benchmark_comparison*/tables/model_ready/
```

These model-ready tables keep IDs, one label, `source_objective`, `atlas_score`, and `consensus_score`, while the full evidence tables keep assay/status/threshold/AC50 metadata for audit only. Leakage checks now forbid source labels/status fields, assay metadata, activity thresholds, and SPD label-defining PK/AC50 fields as predictive features.

## Normalized Mechanism Evidence

The configured normalized mechanism files are now staged in the current VDS repo:

| source | normalized file | rows |
|---|---|---:|
| SIDER | `data/external/sider/drug_adr.tsv` | 771263 |
| CTD | `data/external/ctd/ctd_edges.tsv` | 773625 |
| OpenTargets Safety | `data/external/opentargets/safety.tsv` | 753 |
| Reactome | `data/external/reactome/pathways.tsv` | 1384 |

Refresh command:

```bash
python -m analysis.cli.stage_mechanism_sources \
  --base-dir data/external \
  --pair-table data/pilotstudy/source_benchmark_comparison_v2/tables/ml_source_transfer_table.csv \
  --fda-mapping chemdb/data/fda_mapping_from_pdbqt.csv
```

Pilot normalized-source graph:

| output | rows |
|---|---:|
| mechanism graph nodes | 13877 |
| mechanism graph edges | 1270332 |
| pair mechanism score rows | 7018 |
| strict mechanism triad positives | 69 |

The normalized mechanism target table is:

```text
data/pilotstudy/evaluation_targets_normalized_mechanism/mechanism_adr_support_target.csv
```

| mechanism label status | rows |
|---|---:|
| labeled mechanism triad | 69 |
| unknown no triad evidence | 6949 |

Interpretation: normalized mechanism staging is no longer blocked. It is still not a full supervised mechanism classifier because the current table has strict positives plus unknown background, not curated negatives. Missing mechanism evidence remains unknown, not negative.

## Four-State Evidence Collapse

The current pilot now has a raw-evidence and collapsed-label layer:

```bash
python -m analysis.cli.build_four_state_evidence \
  --chembl data/external/chembl/bioactivity.tsv \
  --bindingdb data/external/bindingdb/bioactivity.tsv \
  --papyrus data/external/papyrus/bioactivity.tsv \
  --toxcast data/external/toxcast/bioactivity.tsv \
  --tox21 data/external/tox21/bioactivity.tsv \
  --pubchem data/external/pubchem_bioassay/pubchem_bioassay.tsv \
  --iuphar-gtopdb data/external/iuphar_gtopdb/bioactivity.tsv \
  --spd data/pilotstudy/source_benchmark_comparison_v2/tables/ml_spd_exposure_table.csv \
  --omop-ohdsi data/external/ohdsi/omopReferenceSet.csv data/external/ohdsi/ohdsiNegativeControls.csv data/external/ohdsi/ohdsiDevelopmentNegativeControls.csv \
  --faers data/external/faers/ohdsi_omop_faers_disproportionality.csv \
  --feature-table data/pilotstudy/source_benchmark_comparison_v2/tables/ml_source_transfer_table.csv \
  --out-dir data/pilotstudy/four_state_evidence
```

Outputs:

- `data/pilotstudy/four_state_evidence/four_state_raw_evidence.csv`
- `data/pilotstudy/four_state_evidence/four_state_collapsed_labels.csv`
- `data/pilotstudy/four_state_evidence/four_state_joined_feature_table.csv`
- `data/pilotstudy/four_state_evidence/four_state_pair_mapping_audit.csv`

Current pilot collapse:

| item | rows |
|---|---:|
| raw evidence rows | 37293 |
| collapsed canonical pair rows | 11451 |
| joined ML feature rows | 7018 |
| joined trainable rows | 6727 |
| strict positives | 540 |
| measured/reliable negatives | 6187 |
| conflict/excluded joined rows | 154 |

Four-state joined negative mix:

| dataset | negatives | negative rate | ToxCast share | curated share |
|---|---:|---:|---:|---:|
| all four-state joined | 6187 | 92.0% | 94.9% | 5.1% |
| source-balanced sensitivity | 636 | 54.1% | 50.0% | 50.0% |

The 92% negative rate is expected when all ToxCast measured inactive rows are retained. Keep this as the main measured-evidence table only if the claim is explicitly ToxCast-heavy. The source-balanced table is the better sensitivity check for whether the signal survives without one source dominating the negatives.

Current four-state drug-holdout performance with binding-only features:

| model | test rows | positives | negatives | AUROC | AUPRC | Brier | ECE |
|---|---:|---:|---:|---:|---:|---:|---:|
| previous measured bioactivity | 1384 | 161 | 1223 | 0.739 | 0.314 | 0.215 | 0.329 |
| four-state all labels | 1387 | 106 | 1281 | 0.688 | 0.213 | 0.216 | 0.359 |
| four-state source-balanced | 222 | 93 | 129 | 0.663 | 0.643 | 0.233 | 0.072 |

Interpretation: the stricter four-state policy reduces apparent performance on the full imbalanced table because conflicts and gray-zone rows are removed and ToxCast negatives dominate. Source balancing improves calibration and AUPRC but changes prevalence, so it is a sensitivity analysis rather than a replacement for the full measured table.

## External ML Source Staging

The source-staging manifest is:

```bash
python -m analysis.cli.stage_ml_label_sources \
  --out-dir data/external/ml_label_sources \
  --sources chembl bindingdb toxcast tox21 pubchem_bioassay iuphar_gtopdb ncats_inxight \
  --download-small
```

Current staged source status:

| source | status | use policy |
|---|---|---|
| ChEMBL | registered/staged export expected | binding and functional activity labels after assay-type and target-confidence QC |
| BindingDB | full 202605 TSV zip staged | direct binding labels after filtered local export |
| ToxCast/Tox21 | registered/staged export expected | HTS activity benchmark after assay, chemical, cytotoxicity, and target-mapping QC |
| PubChem BioAssay | manifest only | explicit assay active/inactive outcomes by AID after assay-level QC |
| IUPHAR/BPS GtoPdb | registered; small interaction export can be staged | high-precision ligand-target positives; absence/weak affinity are not negatives |
| NCATS Inxight Drugs | registered; FRDB PK/safety context expected | PK, drug-status, toxicity, adverse-event, and DDI context only; not drug-target labels |
| ExCAPE-DB | manifest only | training-only; high overlap risk with ChEMBL/PubChem |
| LIT-PCBA subset | downloaded | benchmark-only; not production truth |
| DUD-E | manifest only | docking benchmark only |
| MUV | manifest only | virtual-screening benchmark only |

The repo is on the VDS work filesystem, not the home-quota-limited path. Large sources may be downloaded when there is a direct source URL. Current direct staged files:

- `data/external/bindingdb/BindingDB_All_202605_tsv.zip`
- `data/external/lit_pcba/LITPCBA_9t_subset.tar.xz`

For query/portal sources such as PubChem BioAssay and ExCAPE, stage filtered exports at the manifest `filtered_path`, then pass them to `build_four_state_evidence`. For NCATS Inxight FRDB, extract `frdb-pk.tsv` from the downloaded archive and build exposure context separately:

```bash
python -m analysis.cli.build_pk_table \
  --out data/external/ncats_inxight/pk_table.tsv \
  --ncats-inxight-pk data/external/ncats_inxight/frdb-pk.tsv
```

Then pass that PK table to `python -m analysis.cli.stage_pilot_bioactivity_sources --pk-table ...`; missing PK remains unknown, not negative.

## ADR-Keyed Control Table

OMOP/OHDSI and FAERS controls are drug-ADR labels, not drug-target labels. The separate ADR-keyed table is:

```bash
python -m analysis.cli.build_adr_mechanism_panel \
  --pair-table data/pilotstudy/source_benchmark_comparison_v2/tables/ml_source_transfer_table.csv \
  --four-state-labels data/pilotstudy/four_state_evidence/four_state_collapsed_labels.csv \
  --mechanism-scores data/pilotstudy/mechanism_graph_normalized_sources/pair_mechanism_scores.csv \
  --out data/pilotstudy/adr_mechanism_panel/adr_mechanism_panel_table.csv
```

Current pilot ADR-control overlap:

| item | rows |
|---|---:|
| drug-ADR control rows | 909 |
| labelable drug-ADR controls | 516 |
| controls with Atlas pilot drug features | 96 |
| positives with Atlas features | 27 |
| negatives with Atlas features | 69 |

ADR-control drug-holdout smoke performance:

| model | test rows | positives | negatives | AUROC | AUPRC | Brier | ECE |
|---|---:|---:|---:|---:|---:|---:|---:|
| ADR panel binding summary | 18 | 4 | 14 | 0.429 | 0.297 | 0.289 | 0.329 |
| ADR panel mechanism summary | 18 | 4 | 14 | 0.339 | 0.272 | 0.295 | 0.340 |

Interpretation: this is a useful calibration/control table, but the pilot overlap is too small for a publishable ADR classifier. The mechanism-summary features do not currently improve the ADR-control smoke model.

## Added ML Diagnostics

Run diagnostics with:

```bash
python -m analysis.cli.run_ml_diagnostics \
  --dataset data/pilotstudy/four_state_evidence/four_state_joined_feature_table.csv \
  --label four_state_ml_label \
  --feature-set pilot_binding_only \
  --model logistic_regression \
  --split drug_holdout \
  --out-dir data/pilotstudy/four_state_evidence/ml_four_state_diagnostics
```

Outputs:

- `learning_curve_metrics.csv`
- `learning_curve_log_loss.png`
- `learning_curve_auc.png`
- `feature_ablation_metrics.csv`
- `score_baseline_metrics.csv`

Current four-state learning curve is stable but shallow:

| train fraction | test log loss | test AUPRC | test AUROC |
|---:|---:|---:|---:|
| 0.10 | 0.636 | 0.209 | 0.683 |
| 0.25 | 0.640 | 0.212 | 0.687 |
| 0.50 | 0.623 | 0.210 | 0.686 |
| 0.75 | 0.620 | 0.210 | 0.685 |
| 1.00 | 0.617 | 0.213 | 0.688 |

Interpretation: more rows improve log loss slightly but do not materially improve AUPRC. That suggests the current pilot features are the bottleneck, not just the amount of training data.

Best-practice checks now covered:

- drug-holdout splits for current smoke models
- source-aware train/test tables
- source-balanced sensitivity analysis
- feature ablation reruns
- score-only baselines
- learning curves with log loss and held-out AUROC/AUPRC
- calibration/reliability outputs from model training

Still needed before publication:

- target-holdout and scaffold-holdout repeats on the final full-sized datasets
- temporal validation where activity/evidence year is available
- external-source holdout claims only when overlap is adequate
- confidence intervals on every final reported AUROC/AUPRC/EF metric
- applicability-domain reporting by drug/scaffold/target class
