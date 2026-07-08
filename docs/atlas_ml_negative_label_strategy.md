# Atlas ML Negative Label Strategy

This note is the durable policy for Atlas mechanism/ADR machine-learning labels. Future agents should use it before adding mechanism labels, negative controls, or decoy-backed training data.

## Core Policy

Atlas has sparse positives and many unknowns. Missing evidence is not negative evidence.

Use a four-state mechanism label:

| value | meaning | allowed use |
|---:|---|---|
| `1` | strict positive mechanism triad | supervised positives, enrichment positives |
| `0` | measured or reliable negative | sensitivity training, calibration/control analyses |
| `NaN` | unknown | PU/background rows, matched-background sampling |
| `-1` | excluded, ambiguous, or conflicting | exclude from supervised training |

The primary publishable mechanism model should remain PU or matched-background. A measured/reliable-negative supervised model is a sensitivity analysis, not the primary claim.

## Negative Evidence Tiers

Only assign `mechanism_label = 0` when at least one strong negative condition is present.

1. **Measured inactive assay pairs**
   - Sources: SPD, ToxCast, ChEMBL, Papyrus, PubChem BioAssay, LIT-PCBA, ExCAPE-DB, BindingDB, Davis/KIBA kinase panels where the pair was actually tested.
   - Use for: bioactivity and exposure-relevance negatives.
   - Caveat: target binding inactivity does not prove the target cannot mediate an ADR.

2. **OMOP/OHDSI drug-outcome negative controls**
   - Sources: OMOP reference set, OHDSI MethodEvaluation negative controls, EU-ADR-style reference sets where available.
   - Use for: drug-ADR calibration and sensitivity controls.
   - Caveat: these are drug-event controls, not pair-level target-mechanism negatives. Published audits report meaningful misclassification, so keep provenance and confidence fields.

3. **FAERS non-signals with adequate reporting opportunity**
   - Required condition: enough drug and event reports to make a non-signal meaningful, no ROR/PRR/EBGM signal under the configured thresholds, and no positive evidence in SIDER/OFFSIDES/CTD/OpenTargets.
   - Use for: weak negatives or calibration background.
   - Caveat: FDA spontaneous reports do not establish causality, and absence of a FAERS signal is not a strict negative.

4. **Reliable-negative mining**
   - Select from unknowns using external drug, target, and ADR similarity only.
   - Recompute inside each training fold; do not preselect reliable negatives globally before splitting.
   - Use for: sensitivity training when measured negatives are too sparse.
   - Caveat: label these `reliable_negative`, not `confirmed_negative`.

5. **Decoys**
   - Sources: DeepCoy, DUD-E, DEKOIS, or the existing diversity-library docked decoys.
   - Use for: docking-recognition benchmark panels only.
   - Do not merge decoys into FDA/off-target ADR mechanism labels.

## Leakage Rules

Never use Atlas-derived ranking/confidence fields to choose negatives when those fields or related scores are model features.

Forbidden for negative selection:

- `atlas_score`
- `SCORCH_score_used`
- `final_score`
- `consensus_score`
- `z_score`
- `z_selected`
- `fdr_q_value`
- heatmap rank or top-percentile status
- `mechanism_graph_score`

Forbidden as predictive features:

- `drug_id`
- `target_id`
- `pdb_id`
- raw UniProt one-hot features
- `literature_supported_label`
- side-effect labels
- source membership labels unless the experiment is explicitly source-stratified

Reliable-negative mining may use external descriptors such as chemical fingerprints, ATC class, target sequence/family similarity, tissue/ontology distance, and ADR ontology distance, but the exact columns must be written to `negative_selection_features`.

## Required Label Columns

Mechanism ML tables should include these audit columns:

- `mechanism_label`
- `mechanism_label_source`
- `mechanism_label_confidence`
- `negative_evidence_type`
- `negative_source`
- `negative_confidence`
- `negative_selection_fold`
- `negative_selection_features`
- `excluded_reason`

Use `NaN` for unknown values and `-1` only when a row is ambiguous, contradictory, or intentionally excluded.

## Recommended Modeling Setup

1. **Primary model**
   - PU or matched-background mechanism model.
   - Unknown rows remain unknown/background with explicit weights.

2. **Sensitivity model**
   - Strict positives versus measured/reliable negatives.
   - Report source composition and class balance.

3. **Calibration/control**
   - OMOP/OHDSI negative controls and FAERS non-signals.
   - Report as drug-ADR controls, not target-mechanism ground truth.

4. **Benchmark only**
   - Decoy panels for docking-recognition.
   - Report separately from FDA/off-target ML.

## Candidate Sources For More Negatives

Prioritize sources that intentionally record tested inactive or negative outcomes:

- PubChem BioAssay confirmatory/dose-response assays with explicit `inactive` outcomes.
- LIT-PCBA, which curates PubChem dose-response actives and confirmed inactives.
- ExCAPE-DB, which integrates ChEMBL and PubChem active/inactive labels for target-prediction benchmarks.
- BindingDB rows with censored weak/non-binding values, such as `Ki`, `Kd`, `IC50`, or `EC50` above a configured threshold.
- Davis and KIBA kinase panels for kinase-specific weak/non-binding controls.
- OMOP/OHDSI reference sets for drug-outcome negative controls.
- SIDER placebo-frequency fields only as weak controls when placebo frequency is comparable to or higher than drug frequency.

Do not treat IUPHAR/GtoPdb absence as negative; it is curated toward active ligand-target knowledge. Do not treat missing NCATS Inxight PK/safety metadata as negative evidence; use it as exposure and safety context only.

## Implemented Ingestion Commands

Normalize local negative-evidence files into one reusable schema:

```bash
python -m analysis.cli.build_negative_evidence \
  --out data/pilotstudy/negative_evidence/negative_evidence.csv \
  --mapping data/mappings/drug_target_mapping.csv \
  --measured data/external/pubchem/pubchem_confirmatory.csv data/external/lit_pcba/lit_pcba.csv \
  --measured-format generic_bioactivity \
  --inactive-threshold-nm 10000 \
  --omop-ohdsi data/external/ohdsi/negative_controls.csv \
  --faers-nonsignal data/external/faers/faers_disproportionality.csv \
  --positive-drug-adr-evidence data/external/sider/drug_adr.tsv
```

`--measured` accepts PubChem BioAssay, LIT-PCBA, ExCAPE-DB, BindingDB, Davis, KIBA, or generic active/inactive bioactivity tables through shared column aliases. Use `--measured-format pubchem_bioassay` when PubChem activity outcomes are encoded as PubChem outcome codes.

Build the four-state mechanism label table:

```bash
python -m analysis.cli.build_mechanism_label_table \
  --pair-table data/pilotstudy/bioactivity_pair_table.csv \
  --mechanism-scores data/pilotstudy/mechanism_graph_normalized_sources/pair_mechanism_scores.csv \
  --negative-evidence data/pilotstudy/negative_evidence/negative_evidence.csv \
  --out data/pilotstudy/mechanism_four_state_labels.csv
```

Mine reliable negatives only inside the selected training fold and only from external descriptors:

```bash
python -m analysis.cli.build_reliable_negatives \
  --dataset data/pilotstudy/mechanism_four_state_labels.csv \
  --label mechanism_label \
  --selection-features chemical_fingerprint_distance target_family_distance adr_ontology_distance \
  --split drug_holdout \
  --out data/pilotstudy/negative_evidence/reliable_negatives_drug_holdout.csv
```

The reliable-negative command refuses Atlas-derived selection fields such as `atlas_score`, `SCORCH_score_used`, `final_score`, `consensus_score`, `z_selected`, `fdr_q_value`, and `mechanism_graph_score`.

## Deep Research Prompt

Use this prompt for a dedicated literature/data-source search:

```text
Find public datasets and papers that intentionally report negative or inactive evidence for drug-target binding, compound-target bioactivity, drug-ADR associations, or target-ADR mechanisms. Prioritize downloadable datasets with explicit inactive/negative outcomes rather than absence of evidence.

For each candidate source, report:
1. dataset/paper name, citation, URL, and license/access constraints;
2. whether negatives are measured inactive assays, curated drug-outcome negative controls, FAERS-style non-signals, reliable-negative mined examples, or docking decoys;
3. exact fields needed to identify actives, inactives, inconclusive rows, and excluded rows;
4. recommended thresholds for active/inactive labels, including how to handle censored values such as ">" or "<";
5. target coverage, drug/compound coverage, and whether FDA-approved drugs can be mapped;
6. whether targets can be mapped to UniProt/PDB/Atlas target IDs;
7. known bias or leakage risks;
8. whether the source is appropriate for core mechanism ML, sensitivity training, calibration/control, or benchmark-only use;
9. implementation steps to ingest the source into a four-state Atlas label table with 1 = strict positive, 0 = measured/reliable negative, NaN = unknown, and -1 = excluded/ambiguous/conflicting.

Do not recommend using Atlas score, SCORCH score, FDR, heatmap rank, or mechanism_graph_score to choose negatives, because those can be model features and would leak the answer.
```
