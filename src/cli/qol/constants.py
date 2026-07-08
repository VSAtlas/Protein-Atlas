from __future__ import annotations




DEMO_RUN_ID = "atlas_demo"


SLURM_SUBMIT_DEFAULT_ARRAY = "0-31%4"


SLURM_SUBMIT_DEFAULT_CPUS_PER_TASK = 8


SLURM_SUBMIT_DEFAULT_SCRIPT = "tools/slurm/run_spr.sh"


SLURM_SUBMIT_DEFAULT_JOB_NAME = "atlas"


SLURM_BENCH2_CANARY_ARRAY = "0-2%3"


SLURM_BENCH2_CANARY_CPUS_PER_TASK = "auto"


SLURM_BENCH2_CANARY_SCRIPT = "tools/slurm/run_bench2_canary.sh"


SLURM_BENCH2_CANARY_TIME = "00:15:00"


SLURM_BENCH2_CANARY_JOB_NAME = "atlas_bench2_canary"


SLURM_BENCH2_CANARY_MAIN_ARGS = "-bench2 -fast --run-id {run_id}"


DEMO_ROWS: tuple[dict[str, str], ...] = (
    {
        "ligand_display": "Imatinib",
        "ligand_base": "imatinib",
        "drug_id": "imatinib",
        "target_id": "ABL1|HOLO|pH7_4",
        "pdb_id": "ABL1",
        "variant": "HOLO",
        "ph_label": "pH7_4",
        "target_name": "ABL1 kinase",
        "protein_class": "Kinase",
        "z_selected": "2.45",
        "atlas_score": "0.96",
        "tissue_expression": "high:hematopoietic",
        "pathway_link_to_side_effect": "BCR-ABL signaling",
    },
    {
        "ligand_display": "Imatinib",
        "ligand_base": "imatinib",
        "drug_id": "imatinib",
        "target_id": "KIT|HOLO|pH7_4",
        "pdb_id": "KIT",
        "variant": "HOLO",
        "ph_label": "pH7_4",
        "target_name": "KIT receptor tyrosine kinase",
        "protein_class": "Kinase",
        "z_selected": "2.05",
        "atlas_score": "0.88",
        "tissue_expression": "medium:skin",
        "pathway_link_to_side_effect": "mast-cell signaling",
    },
    {
        "ligand_display": "Imatinib",
        "ligand_base": "imatinib",
        "drug_id": "imatinib",
        "target_id": "EGFR|HOLO|pH7_4",
        "pdb_id": "EGFR",
        "variant": "HOLO",
        "ph_label": "pH7_4",
        "target_name": "EGFR",
        "protein_class": "Kinase",
        "z_selected": "0.72",
        "atlas_score": "0.41",
        "tissue_expression": "medium:epithelium",
        "pathway_link_to_side_effect": "ERBB signaling",
    },
    {
        "ligand_display": "Gefitinib",
        "ligand_base": "gefitinib",
        "drug_id": "gefitinib",
        "target_id": "ABL1|HOLO|pH7_4",
        "pdb_id": "ABL1",
        "variant": "HOLO",
        "ph_label": "pH7_4",
        "target_name": "ABL1 kinase",
        "protein_class": "Kinase",
        "z_selected": "0.31",
        "atlas_score": "0.27",
        "tissue_expression": "high:hematopoietic",
        "pathway_link_to_side_effect": "BCR-ABL signaling",
    },
    {
        "ligand_display": "Gefitinib",
        "ligand_base": "gefitinib",
        "drug_id": "gefitinib",
        "target_id": "KIT|HOLO|pH7_4",
        "pdb_id": "KIT",
        "variant": "HOLO",
        "ph_label": "pH7_4",
        "target_name": "KIT receptor tyrosine kinase",
        "protein_class": "Kinase",
        "z_selected": "0.58",
        "atlas_score": "0.35",
        "tissue_expression": "medium:skin",
        "pathway_link_to_side_effect": "mast-cell signaling",
    },
    {
        "ligand_display": "Gefitinib",
        "ligand_base": "gefitinib",
        "drug_id": "gefitinib",
        "target_id": "EGFR|HOLO|pH7_4",
        "pdb_id": "EGFR",
        "variant": "HOLO",
        "ph_label": "pH7_4",
        "target_name": "EGFR",
        "protein_class": "Kinase",
        "z_selected": "2.39",
        "atlas_score": "0.95",
        "tissue_expression": "medium:epithelium",
        "pathway_link_to_side_effect": "ERBB signaling",
    },
    {
        "ligand_display": "Loratadine",
        "ligand_base": "loratadine",
        "drug_id": "loratadine",
        "target_id": "ABL1|HOLO|pH7_4",
        "pdb_id": "ABL1",
        "variant": "HOLO",
        "ph_label": "pH7_4",
        "target_name": "ABL1 kinase",
        "protein_class": "Kinase",
        "z_selected": "-0.64",
        "atlas_score": "0.08",
        "tissue_expression": "high:hematopoietic",
        "pathway_link_to_side_effect": "BCR-ABL signaling",
    },
    {
        "ligand_display": "Loratadine",
        "ligand_base": "loratadine",
        "drug_id": "loratadine",
        "target_id": "KIT|HOLO|pH7_4",
        "pdb_id": "KIT",
        "variant": "HOLO",
        "ph_label": "pH7_4",
        "target_name": "KIT receptor tyrosine kinase",
        "protein_class": "Kinase",
        "z_selected": "-0.22",
        "atlas_score": "0.13",
        "tissue_expression": "medium:skin",
        "pathway_link_to_side_effect": "mast-cell signaling",
    },
    {
        "ligand_display": "Loratadine",
        "ligand_base": "loratadine",
        "drug_id": "loratadine",
        "target_id": "EGFR|HOLO|pH7_4",
        "pdb_id": "EGFR",
        "variant": "HOLO",
        "ph_label": "pH7_4",
        "target_name": "EGFR",
        "protein_class": "Kinase",
        "z_selected": "-0.48",
        "atlas_score": "0.09",
        "tissue_expression": "medium:epithelium",
        "pathway_link_to_side_effect": "ERBB signaling",
    },
)
