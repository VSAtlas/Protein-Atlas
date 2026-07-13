from __future__ import annotations


PHYSICHEM_DESCRIPTOR_FEATURES = [
    "rdkit_mol_wt",
    "rdkit_mol_logp",
    "rdkit_tpsa",
    "rdkit_hbd",
    "rdkit_hba",
    "rdkit_rotatable_bonds",
    "rdkit_formal_charge",
    "rdkit_aromatic_rings",
    "rdkit_fraction_csp3",
    "rdkit_qed",
]

PHYSICHEM_DESCRIPTOR_FEATURES_NO_QED = [
    feature for feature in PHYSICHEM_DESCRIPTOR_FEATURES if feature != "rdkit_qed"
]

SHORTCUT_REDUCED_DESCRIPTOR_FEATURES = [
    feature
    for feature in PHYSICHEM_DESCRIPTOR_FEATURES
    if feature not in {"rdkit_mol_wt", "rdkit_mol_logp", "rdkit_qed"}
]

PAIR_FINAL_SCORE_FEATURES = [
    "final_score",
    "banana_score_normalized",
]

PAIR_FINAL_FULL_NO_QED_FEATURES = [
    *PAIR_FINAL_SCORE_FEATURES,
    *PHYSICHEM_DESCRIPTOR_FEATURES_NO_QED,
]

# Default clean feature sets quarantine high-memorization context columns.
# The columns remain in model-ready tables for splitting, audit, and ablation,
# but clean models do not train on them unless a *_with_context feature set is
# explicitly selected.
QUARANTINED_CONTEXT_FEATURES = ["structure_quality", "protein_class", "target_family"]

# Chemistry group identifiers are retained for leakage audits, OOD splits, and
# data-gap reports, but they are not used as predictive inputs. They can encode
# scaffold/chemotype priors directly enough to create memorization under weak
# splits.
QUARANTINED_CHEMISTRY_GROUP_FEATURES = [
    "ligand_chemotype",
    "scaffold_key",
    "chemical_cluster",
]

CLEAN_BINDING_EXPERT_WITH_CONTEXT_FEATURES = [
    "consensus_score",
    "banana_score_normalized",
    "structure_quality",
    "protein_class",
    "target_family",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
]

CLEAN_BINDING_EXPERT_FEATURES = [
    "consensus_score",
    "banana_score_normalized",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
]

ATLAS_PRIOR_BINDING_FEATURES = [
    "atlas_binding_prior",
    "banana_score_normalized",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
]

CONSENSUS_Z_BINDING_FEATURES = [
    "consensus_z_score",
    "banana_score_normalized",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
]

CONSENSUS_Z_BANANA_SCORE_FEATURES = [
    "consensus_z_score",
    "banana_score_normalized",
]

CONSENSUS_Z_POTENCY_FEATURES = [
    "consensus_z_score",
    "banana_score_normalized",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
]

CLEAN_BINDING_EXPERT_NO_TARGET_CLASS_FEATURES = [
    *CLEAN_BINDING_EXPERT_FEATURES,
]

CONSENSUS_CONTEXT_BINDING_FEATURES = [
    "consensus_score",
    "structure_quality",
    "protein_class",
    "target_family",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
]

CONSENSUS_CONTEXT_BINDING_NO_TARGET_CLASS_FEATURES = [
    "consensus_score",
    "structure_quality",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
]

RDKIT_TARGET_METADATA_FEATURES = [
    "structure_quality",
    "protein_class",
    "target_family",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
]

RDKIT_CHEMISTRY_ONLY_FEATURES = [
    "structure_quality",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
]

CONSENSUS_BANANA_SCORE_FEATURES = [
    "consensus_score",
    "banana_score_normalized",
]

FROZEN_BINDING_PRIOR_INPUTS = [
    "banana_binding_probability",
    "banana_score_normalized",
    "atlas_score_normalized",
    "consensus_score_normalized",
    "SCORCH_score_used_normalized",
    "final_score_normalized",
]

FROZEN_BINDING_PRIOR_OUTPUTS = [
    "banana_atlas_blend_score",
    "binding_expert_score",
]

DIRECT_SPD_EXPOSURE_WITH_CONTEXT_FEATURES = [
    "binding_expert_score",
    "structure_quality",
    "protein_class",
    "target_family",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
]

DIRECT_SPD_EXPOSURE_FEATURES = [
    "binding_expert_score",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
]

TISSUE_EXPRESSION_FEATURES = [
    "expression_presence_score",
    "site_specificity_score",
    "expression_concordance_score",
    "hpa_consensus_ntpm",
    "hpa_specificity_score",
    "gtex_median_tpm",
    "bgee_expression_score",
    "ot_expression_value",
    "ot_expression_zscore",
    "distribution_penalty",
    "sensitive_offsite_expression_score",
]

TISSUE_SITE_COMPOSITE_FEATURES = [
    "site_relevance_score",
]

TISSUE_SAFETY_EVIDENCE_FEATURES = [
    "site_safety_risk_score",
    "ot_safety_liability_count",
]

TISSUE_SITE_FEATURES = [
    *TISSUE_SITE_COMPOSITE_FEATURES,
    *TISSUE_EXPRESSION_FEATURES,
    *TISSUE_SAFETY_EVIDENCE_FEATURES,
]

PAIR_TISSUE_SITE_RELEVANCE_WITH_CONTEXT_FEATURES = [
    "binding_expert_score",
    "binding_probability_clean",
    "exposure_relevance_probability",
    "site_relevance_score",
    "expression_presence_score",
    "site_specificity_score",
    "expression_concordance_score",
    "protein_class",
    "target_family",
]

PAIR_TISSUE_SITE_RELEVANCE_FEATURES = [
    "binding_expert_score",
    "binding_probability_clean",
    "exposure_relevance_probability",
    "site_relevance_score",
    "expression_presence_score",
    "site_specificity_score",
    "expression_concordance_score",
]

MECHANISM_EVIDENCE_FEATURES = [
    "mechanism_graph_score",
    "mechanism_path_count",
    "drug_adr_known",
    "target_adr_known",
    "target_pathway_adr_link",
    "triad_complete",
    "target_adr_evidence",
    "pathway_evidence",
]

CLEAN_MECHANISM_WITH_CONTEXT_FEATURES = [
    "binding_expert_score",
    "binding_probability_clean",
    "exposure_relevance_probability",
    "site_relevance_score",
    "tissue_site_relevance_probability",
    "clean_tissue_site_relevance_score",
    "structure_quality",
    "protein_class",
    "target_family",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
]

CLEAN_MECHANISM_FEATURES = [
    "binding_expert_score",
    "binding_probability_clean",
    "exposure_relevance_probability",
    "site_relevance_score",
    "tissue_site_relevance_probability",
    "clean_tissue_site_relevance_score",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
]

FLAT_MECHANISM_BASELINE_WITH_CONTEXT_FEATURES = [
    "consensus_score",
    "banana_score_normalized",
    "structure_quality",
    "protein_class",
    "target_family",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
    *TISSUE_EXPRESSION_FEATURES,
    "site_relevance_score",
]

FLAT_MECHANISM_BASELINE_FEATURES = [
    "consensus_score",
    "banana_score_normalized",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
    *TISSUE_EXPRESSION_FEATURES,
    "site_relevance_score",
]

CARDIAC_QT_MECHANISM_CLEAN_WITH_CONTEXT_FEATURES = [
    "consensus_score",
    "banana_score_normalized",
    "structure_quality",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
    "expression_presence_score",
    "site_specificity_score",
    "expression_concordance_score",
    "hpa_consensus_ntpm",
    "hpa_specificity_score",
    "gtex_median_tpm",
    "bgee_expression_score",
    "ot_expression_value",
    "ot_expression_zscore",
    "distribution_penalty",
    "sensitive_offsite_expression_score",
]

CARDIAC_QT_MECHANISM_CLEAN_FEATURES = [
    "consensus_score",
    "banana_score_normalized",
    *PHYSICHEM_DESCRIPTOR_FEATURES,
    "expression_presence_score",
    "site_specificity_score",
    "expression_concordance_score",
    "hpa_consensus_ntpm",
    "hpa_specificity_score",
    "gtex_median_tpm",
    "bgee_expression_score",
    "ot_expression_value",
    "ot_expression_zscore",
    "distribution_penalty",
    "sensitive_offsite_expression_score",
]

BINDING_EXPERT_FEATURES = CLEAN_BINDING_EXPERT_FEATURES


FEATURE_SETS = {
    "binding_only": ["atlas_score", "mmgbsa_score", "structure_quality", "protein_class", "ligand_chemotype"],
    "pilot_nonleaky": [
        "atlas_score",
        "consensus_score",
    ],
    "pilot_exposure_sensitivity": [
        "atlas_score",
        "consensus_score",
        "free_cmax_um",
        "cmax_um",
        "fraction_unbound_plasma",
    ],
    "pilot_binding_only": [
        "atlas_score",
        "consensus_score",
    ],
    "atlas_only": [
        "atlas_score",
    ],
    "consensus_only": [
        "consensus_score",
    ],
    "scorch_only": [
        "pilot_SCORCH_score_used",
        "SCORCH_score_used",
        "pilot_final_score",
    ],
    "consensus_scorch": [
        "consensus_score",
        "pilot_SCORCH_score_used",
        "SCORCH_score_used",
        "pilot_final_score",
    ],
    "atlas_scorch": [
        "atlas_score",
        "consensus_score",
        "pilot_SCORCH_score_used",
        "SCORCH_score_used",
        "pilot_final_score",
    ],
    "atlas_scorch_banana": [
        "atlas_score",
        "consensus_score",
        "pilot_SCORCH_score_used",
        "SCORCH_score_used",
        "pilot_final_score",
        "banana_score",
        "banana_binding_probability",
        "banana_score_normalized",
        "banana_atlas_blend_score",
        "binding_expert_score",
    ],
    "banana_only": [
        "banana_score",
        "banana_binding_probability",
        "banana_score_normalized",
    ],
    "frozen_binding_prior_inputs": [
        *FROZEN_BINDING_PRIOR_INPUTS,
    ],
    "frozen_binding_prior_baseline": [
        *FROZEN_BINDING_PRIOR_OUTPUTS,
    ],
    "consensus_banana_binding": [
        *CLEAN_BINDING_EXPERT_FEATURES,
    ],
    "atlas_banana": [
        "atlas_score",
        "banana_score",
        "banana_binding_probability",
        "banana_score_normalized",
        "banana_atlas_blend_score",
        "binding_expert_score",
    ],
    "banana_atlas_binding": [
        "atlas_score",
        "consensus_score",
        "banana_score",
        "banana_binding_probability",
        "banana_score_normalized",
        "banana_atlas_blend_score",
        "binding_expert_score",
    ],
    "mechanism_consensus_banana_graph": [
        "consensus_score",
        "banana_score",
        "banana_binding_probability",
        "banana_score_normalized",
        "banana_atlas_blend_score",
        "binding_expert_score",
        "mechanism_graph_score",
        "mechanism_path_count",
        "drug_adr_known",
        "target_adr_known",
        "target_pathway_adr_link",
        "triad_complete",
    ],
    "scorch_binding_only": [
        "atlas_score",
        "consensus_score",
        "pilot_SCORCH_score_used",
        "pilot_final_score",
    ],
    "banana_scorch_binding": [
        "atlas_score",
        "consensus_score",
        "banana_score",
        "banana_binding_probability",
        "banana_score_normalized",
        "binding_expert_score",
        "pilot_SCORCH_score_used",
        "pilot_final_score",
    ],
    "ligand_physchem_descriptors": [
        *PHYSICHEM_DESCRIPTOR_FEATURES,
    ],
    "ligand_physchem_descriptors_no_qed": [
        *PHYSICHEM_DESCRIPTOR_FEATURES_NO_QED,
    ],
    "ligand_physchem_shortcut_reduced": [
        *SHORTCUT_REDUCED_DESCRIPTOR_FEATURES,
    ],
    "spd_binding_pair_final_scores_only": [
        *PAIR_FINAL_SCORE_FEATURES,
    ],
    "spd_binding_pair_final_full_no_qed": [
        *PAIR_FINAL_FULL_NO_QED_FEATURES,
    ],
    "spd_exposure_descriptor_only": [
        *PHYSICHEM_DESCRIPTOR_FEATURES,
            "scaffold_key",
    ],
    "spd_potency_physchem": [
        *CLEAN_BINDING_EXPERT_FEATURES,
    ],
    "spd_potency_atlas_prior": [
        *ATLAS_PRIOR_BINDING_FEATURES,
    ],
    "spd_potency_consensus_z": [
        *CONSENSUS_Z_POTENCY_FEATURES,
    ],
    "spd_binding_nonleaky": [
        *CLEAN_BINDING_EXPERT_FEATURES,
    ],
    "spd_binding_nonleaky_consensus_z": [
        *CONSENSUS_Z_BINDING_FEATURES,
    ],
    "spd_binding_nonleaky_with_context": [
        *CLEAN_BINDING_EXPERT_WITH_CONTEXT_FEATURES,
    ],
    "spd_binding_nonleaky_no_target_class": [
        *CLEAN_BINDING_EXPERT_NO_TARGET_CLASS_FEATURES,
    ],
    "spd_binding_consensus_context_no_banana": [
        *CONSENSUS_CONTEXT_BINDING_FEATURES,
    ],
    "spd_binding_consensus_context_no_banana_no_target_class": [
        *CONSENSUS_CONTEXT_BINDING_NO_TARGET_CLASS_FEATURES,
    ],
    "spd_binding_rdkit_target_metadata": [
        *RDKIT_TARGET_METADATA_FEATURES,
    ],
    "spd_binding_rdkit_chemistry_only": [
        *RDKIT_CHEMISTRY_ONLY_FEATURES,
    ],
    "spd_binding_consensus_banana_scores_only": [
        *CONSENSUS_BANANA_SCORE_FEATURES,
    ],
    "spd_binding_consensus_z_banana_scores_only": [
        *CONSENSUS_Z_BANANA_SCORE_FEATURES,
    ],
    "spd_exposure_nonleaky": [
        *DIRECT_SPD_EXPOSURE_FEATURES,
    ],
    "spd_exposure_nonleaky_with_context": [
        *DIRECT_SPD_EXPOSURE_WITH_CONTEXT_FEATURES,
    ],
    "spd_tissue_site_nonleaky": [
        *PAIR_TISSUE_SITE_RELEVANCE_FEATURES,
    ],
    "spd_tissue_site_nonleaky_with_context": [
        *PAIR_TISSUE_SITE_RELEVANCE_WITH_CONTEXT_FEATURES,
    ],
    "spd_mechanism_nonleaky": [
        *CLEAN_MECHANISM_FEATURES,
    ],
    "spd_mechanism_nonleaky_with_context": [
        *CLEAN_MECHANISM_WITH_CONTEXT_FEATURES,
    ],
    "spd_mechanism_evidence_sensitivity": [
        *CLEAN_MECHANISM_FEATURES,
        *MECHANISM_EVIDENCE_FEATURES,
        *TISSUE_SAFETY_EVIDENCE_FEATURES,
    ],
    "flat_mechanism_baseline": [
        *FLAT_MECHANISM_BASELINE_FEATURES,
    ],
    "flat_mechanism_baseline_with_context": [
        *FLAT_MECHANISM_BASELINE_WITH_CONTEXT_FEATURES,
    ],
    "cardiac_qt_mechanism_clean": [
        *CARDIAC_QT_MECHANISM_CLEAN_FEATURES,
    ],
    "cardiac_qt_mechanism_clean_with_context": [
        *CARDIAC_QT_MECHANISM_CLEAN_WITH_CONTEXT_FEATURES,
    ],
    "tissue_site_score_only": [
        *TISSUE_EXPRESSION_FEATURES,
    ],
    "tissue_site_relevance_score_only": [
        *TISSUE_SITE_COMPOSITE_FEATURES,
    ],
    "mechanism_graph_score_only": [
        *MECHANISM_EVIDENCE_FEATURES,
    ],
    "spd_exposure_pk_sensitivity": [
        "binding_expert_score",
        "structure_quality",
        "protein_class",
        "target_family",
            *PHYSICHEM_DESCRIPTOR_FEATURES,
        "free_cmax_um",
        "cmax_um",
        "fraction_unbound_plasma",
    ],
    "adr_panel_binding_summary": [
        "n_targets_screened",
        "n_pdbs_screened",
        "atlas_score_max",
        "atlas_score_mean",
        "atlas_score_top5_mean",
        "consensus_score_max",
        "consensus_score_mean",
        "consensus_score_top5_mean",
    ],
    "adr_panel_mechanism_summary": [
        "n_targets_screened",
        "n_pdbs_screened",
        "atlas_score_max",
        "atlas_score_mean",
        "atlas_score_top5_mean",
        "consensus_score_max",
        "consensus_score_mean",
        "consensus_score_top5_mean",
        "mechanism_graph_score_max",
        "mechanism_graph_score_mean",
        "mechanism_graph_score_top5_mean",
        "mechanism_path_count_max",
        "mechanism_path_count_mean",
        "target_adr_known_max",
        "target_adr_known_mean",
        "target_pathway_adr_link_max",
        "target_pathway_adr_link_mean",
        "triad_complete_max",
    ],
    "binding_exposure": ["atlas_score", "mmgbsa_score", "free_cmax", "exposure_plausibility", "structure_quality", "protein_class", "ligand_chemotype"],
    "biology_augmented": [
        "atlas_score",
        "mmgbsa_score",
        "free_cmax",
        "exposure_plausibility",
        "tissue_expression",
        "target_adr_evidence",
        "pathway_evidence",
        "structure_quality",
        "protein_class",
        ],
    "full_nonleaky": [
        "atlas_score",
        "mmgbsa_score",
        "tissue_expression",
        "structure_quality",
        "protein_class",
        ],
}
FORBIDDEN_FEATURES = {
    "drug_id",
    "target_id",
    "pdb_id",
    "literature_supported_label",
    *QUARANTINED_CHEMISTRY_GROUP_FEATURES,
}
SPD_EXPOSURE_LABELS = {
    "spd_exposure_label",
    "spd_exposure_relevant",
    "ml_binary_label",
}
SPD_EXPOSURE_DEFINITION_FEATURES = {
    "exposure_margin",
    "spd_exposure_margin",
    "ac50_nM",
    "free_cmax_nM",
    "total_cmax_nM",
    "free_cmax",
    "free_cmax_um",
    "free_cmax_uM",
    "cmax_um",
    "fraction_unbound_plasma",
    "spd_ac50_uM",
    "spd_exposure_relevant",
    "spd_exposure_weak",
    "spd_exposure_unlikely",
    "ml_binary_label",
}
TISSUE_SITE_LABELS = {
    "tissue_site_label",
    "tissue_relevance_label",
    "site_relevance_label",
    "target_site_relevance_label",
}
TISSUE_SITE_DEFINITION_FEATURES = set(TISSUE_SITE_FEATURES) | {
    "adr_site_group",
    "adr_soc",
    "meddra_soc",
    "meddra_pt",
}
MECHANISM_LABELS = {
    "mechanism_ml_label",
    "mechanism_label",
    "four_state_ml_label",
    "drug_target_adr_mechanism_label",
}
MECHANISM_DEFINITION_FEATURES = set(MECHANISM_EVIDENCE_FEATURES) | set(TISSUE_SAFETY_EVIDENCE_FEATURES) | {
    "mechanism_label_status",
    "negative_evidence_type",
    "negative_source",
    "negative_confidence",
}

RETAINED_CONTEXT_AUDIT_COLUMNS = {
    "structure_quality",
    "target_family",
    "protein_class",
    "ligand_chemotype",
    "scaffold_key",
    "chemical_cluster",
    "site_relevance_score",
    "site_relevance_score_by_adr",
    "site_safety_risk_score",
    "ot_safety_liability_count",
    "drug_adr_known",
    "target_adr_known",
    "target_pathway_adr_link",
    "triad_complete",
    "target_adr_evidence",
    "pathway_evidence",
    "mechanism_graph_score",
    "mechanism_path_count",
}


def get_feature_set(name: str) -> list[str]:
    if name not in FEATURE_SETS:
        raise ValueError(f"unknown feature set: {name}")
    return [feature for feature in FEATURE_SETS[name] if feature not in FORBIDDEN_FEATURES]


def normalized_exclude_features(exclude_features: list[str] | None) -> set[str]:
    return {str(f).strip() for f in (exclude_features or []) if str(f).strip()}


def label_definition_exclude_features(label_col: str) -> set[str]:
    label = str(label_col).strip()
    if label in SPD_EXPOSURE_LABELS:
        return set(SPD_EXPOSURE_DEFINITION_FEATURES)
    if label in TISSUE_SITE_LABELS:
        return set(TISSUE_SITE_DEFINITION_FEATURES)
    if label in MECHANISM_LABELS:
        return set(MECHANISM_DEFINITION_FEATURES)
    return set()


def effective_exclude_features(
    label_col: str,
    exclude_features: list[str] | None,
    *,
    allow_label_definition_features: bool = False,
) -> set[str]:
    excluded = normalized_exclude_features(exclude_features)
    if not allow_label_definition_features:
        excluded.update(label_definition_exclude_features(label_col))
    return excluded
