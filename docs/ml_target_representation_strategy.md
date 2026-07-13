# Target Representation Strategy for Atlas Phase 1

## Problem

The current clean binding model has ligand descriptors and Atlas docking scores but no continuous protein representation. `target_id`, `pdb_id`, UniProt one-hot features, `target_family`, and `protein_class` must remain identifiers, split axes, or audit columns rather than clean predictive inputs. A model with no sequence or pocket representation has little basis for transferring learned interaction rules to an unseen target.

## Recommended implementation order

1. **Deterministic pocket descriptors**
   - Derive residue composition, hydrophobic/aromatic/charged fractions, pocket atom count, residue count, radius of gyration, and simple shape/volume proxies from the existing pocket-PDB map.
   - Cache one row per `pdb_id` and pocket definition. Preserve pocket method, radius, control ligand, chain, and structure quality.
   - This is the fastest interpretable baseline and requires no learned labels.

2. **Frozen ESM-2 sequence embeddings**
   - Map every target to a reviewed human UniProt sequence.
   - Generate one frozen mean-pooled embedding per canonical sequence; do not fine-tune ESM-2 on the small SPD panel.
   - Fit dimensionality reduction inside each training fold only. Start with 16-64 PCA components and retain the frozen full embedding as a cached artifact.

3. **Pocket-residue ESM pooling**
   - Reuse the same residue-level ESM-2 embedding, but pool only residues belonging to the mapped docking pocket.
   - This preserves target sequence context while focusing the representation on the local site that produced the docking score.

4. **Optional structural pocket embedding**
   - Compare a frozen pocket descriptor such as DeeplyTough or a PocketVec-style representation with the deterministic and ESM baselines.
   - Treat this as an ablation, not a required dependency for Phase 1.

## Exact Phase 1 feature contract

The first target-representation ablation uses these cached, label-independent
features from each reviewed pocket PDB:

- composition: `pocket_frac_polar`, `pocket_frac_hydrophobic`,
  `pocket_frac_aromatic`, `pocket_frac_charged_pos`,
  `pocket_frac_charged_neg`, and `pocket_net_charge_proxy`;
- geometry: `pocket_rgyr`, `pocket_bbox_vol`, `pocket_pca_ratio1`, and
  `pocket_pca_ratio2`;
- size/definition: `pocket_atom_count`, `pocket_residue_count`, and
  `pocket_radius_a`;
- bound metal context: `pocket_metal_zn`, `pocket_metal_fe`,
  `pocket_metal_mg`, `pocket_metal_ca`, `pocket_metal_mn`,
  `pocket_metal_cu`, `pocket_metal_co`, and `pocket_metal_ni`.

`pocket_definition_source`, `pocket_control_ligand_id`, and pocket-map status
are provenance/audit columns, not clean predictive features. The descriptor
table is generated with:

```bash
python -m analysis.cli.build_target_pocket_features \
  --pocket-map <primary_pocket_map.csv> [recovery_pocket_map.csv ...] \
  --out-dir data/<run_id>/target_representations
```

ESM-2 is not ESMFold. ESM-2 is a sequence language model; ESMFold adds a
structure-prediction head. Frozen ESM-2 embeddings are useful here because they
provide evolutionary and sequence-context information for unseen targets, while
the pocket descriptors provide the explicit structural signal. Atlas should
compare them separately and together.

PocketVec does not expose a small list of physicochemical variables to copy. Its
descriptor is a fixed vector of rankings produced by docking a reference
lead-like ligand panel against each pocket. Reproducing it would therefore be a
separate inverse-screening baseline. Atlas should not call the deterministic
features above "PocketVec"; a true PocketVec-style vector should be added only
as an explicit, fixed-probe-library ablation.


## Model integration

For the first publishable comparison, concatenate the target representation with the existing nonleaky pair features and fit the existing logistic/elastic-net, EBM, CatBoost, or shallow LightGBM candidates. Run the same immutable splits for:

- ligand/RDKit plus consensus score only;
- deterministic pocket descriptors;
- global frozen ESM-2;
- pocket-residue ESM-2;
- global ESM-2 plus pocket descriptors;
- the best target representation plus BANANA as a separate ablation.

Do not use target identifiers or target-family labels as substitutes for a target representation. BANANA's final probability is a pair score, not a reusable protein representation.

## Required validation

- Keep drug, chemical-cluster/scaffold, and target-ID holdouts.
- Add a sequence-cluster holdout, grouping proteins at no more than 40% sequence identity before assigning train/validation/test.
- Add a pocket-similarity cluster holdout when pocket descriptors are available.
- Keep leave-one-target-family-out as a diagnostic because current family label balance is sparse.
- Fit PCA, scalers, calibration, and feature selection inside the training fold only.
- Report nearest training sequence identity and pocket similarity as applicability-domain fields for every prediction.

## Rationale and references

ESM-2 provides pretrained sequence representations learned without Atlas labels ([Lin et al., Science 2023](https://pubmed.ncbi.nlm.nih.gov/36927031/)). ConPLex demonstrates that pretrained protein language representations can support drug-target prediction and experimental hit recovery ([Singh et al., Nature Machine Intelligence 2023](https://pubmed.ncbi.nlm.nih.gov/37289807/)). Pocket descriptors complement global sequence because ligand recognition occurs locally; PocketVec describes this role and its use in proteome-scale pocket comparison ([La Sala et al., Nature Communications 2024](https://www.nature.com/articles/s41467-024-52146-3)). DeeplyTough is an established alignment-free structural pocket embedding baseline ([Simonovsky and Meyers, JCIM 2020](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00554)). Sequence-identity-aware target splits are necessary for a real unseen-target claim; a recent DTA study uses CD-HIT at 40% identity for its novel-target split ([AdaMBind, Nature Communications 2026](https://www.nature.com/articles/s41467-026-70554-5)).

## Phase 1 recommendation

Implement deterministic pocket descriptors and frozen ESM-2 embeddings first. They are reproducible, cacheable across all ligand rows, and can be evaluated with the current tabular model suite. Defer end-to-end graph or meta-learning models until these representations show improvement under target-ID and sequence-cluster holdouts.
