# ML Mechanistic Feature Workflow

## Ligand feature registry

Feature metadata refreshes materialize the registered RDKit ligand descriptors from
canonical SMILES. The recommended binding ligand branch is
`ligand_primary_functional_branch_no_qed`. It contains the nine individual
physicochemical descriptors without QED, twelve functional-group counts, and six
branch/side-chain counts. The five topology indices remain available through
`ligand_topological_shape` and `ligand_full_interpretable_no_qed` for sensitivity
analysis.

Functional counts:

- `rdkit_aldehyde_count`
- `rdkit_ester_count`
- `rdkit_ketone_count`
- `rdkit_amide_count`
- `rdkit_thiol_count`
- `rdkit_carboxylic_acid_count`
- `rdkit_carboxylate_count`
- `rdkit_amine_count`
- `rdkit_primary_amine_count`
- `rdkit_secondary_amine_count`
- `rdkit_tertiary_amine_count`
- `rdkit_benzene_ring_count`

Branch and side-chain counts:

- `rdkit_heavy_atom_branch_point_count`
- `rdkit_terminal_heavy_atom_count`
- `rdkit_ring_nonring_attachment_bond_count`
- `rdkit_nonring_heavy_atom_count`
- `rdkit_bemis_murcko_outside_heavy_atom_count`
- `rdkit_bemis_murcko_side_chain_fragment_count`

Topology sensitivity descriptors:

- `rdkit_bertz_ct`
- `rdkit_hall_kier_alpha`
- `rdkit_kappa1`
- `rdkit_kappa2`
- `rdkit_kappa3`

## Binding grouped OOF

Named feature sets can be combined with explicit extra columns:

```bash
python -m analysis.cli.run_pair_residual_binding \
  --dataset <model_ready.csv> \
  --out-dir <binding_oof_dir> \
  --ligand-feature-set ligand_primary_functional_branch_no_qed \
  --pair-feature-set spd_binding_pair_final_scores_only \
  --outer-group drug_id \
  --outer-group drug_id+target_id \
  --outer-group drug_id+scaffold_key \
  --variant pair-only --variant ligand-only --variant combined \
  --residual-model logistic --residual-model lightgbm
```

The pair set is exactly `final_score` and `banana_score_normalized`. Ligand-only
features are constant across targets for one drug and must not be described as
pair-specific binding evidence.

## Mechanistic PK and exposure

Materialize copy-only, provenance-preserving PK fields before training:

```bash
python -m analysis.cli.materialize_pk_mechanistic_features \
  --input <model_ready.csv> \
  --out <model_ready_pk_mechanistic.csv>
```

Run cold-drug OOF using registered feature sets:

```bash
python -m analysis.cli.run_spd_exposure_grouped_oof \
  --dataset <model_ready_pk_mechanistic.csv> \
  --out-dir <exposure_oof_dir> \
  --potency-feature-set spd_binding_pair_final_full_structure_no_qed \
  --pk-feature-set ligand_full_interpretable_no_qed \
  --pk-feature pk_mech_fraction_unbound_plasma \
  --group-col drug_id --n-splits 5 --model lightgbm
```

Fraction unbound is a label-related mechanistic sensitivity. Protein-binding
percent is its deterministic complement and is audit-only. Dose, route,
formulation, and regimen are training-eligible only when their endpoint scenario
matches the selected Cmax/free-Cmax record. Missing alignment is not imputed.

## Pose interaction features

Claim-grade pair interactions require immutable selected-pose and receptor paths,
hashes, and model indices. The explicit Stage-1/model-1 mode is sensitivity-only:

```bash
python -m analysis.cli.materialize_pair_interaction_features \
  --dataset <model_ready.csv> \
  --out <model_ready_pair_interactions.csv> \
  --pose-policy stage1-model-1 \
  --frozen-run-id <run_id> \
  --frozen-run-root <run_root>
```

The current provider emits distance-shell, atom-type, hydrophobic-carbon,
polar-proximity, metal-proximity, clash-proxy, and partial-charge Coulomb-proxy
features. It does not claim exact hydrogen bonds, salt bridges, pi interactions,
water mediation, buried surface area, or rigorous electrostatics.
