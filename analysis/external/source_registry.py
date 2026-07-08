from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from analysis.external.source_tables import download_to_cache, write_source_manifest


SOURCE_ROLES: dict[str, dict[str, str]] = {
    "spd": {
        "role": "gold_standard_exposure_validation",
        "label_family": "activity_labels,exposure_context",
        "default_claim": "Atlas recovers exposure-relevant secondary pharmacology where standardized SPD assays exist.",
        "reference": "https://www.nature.com/articles/s41467-023-40064-9",
    },
    "papyrus": {
        "role": "large_scale_bioactivity_training",
        "label_family": "activity_labels",
        "default_claim": "Atlas ranks measured drug-target bioactivity in broad public assay corpora.",
        "reference": "https://zenodo.org/records/13987985",
    },
    "chembl": {
        "role": "large_scale_bioactivity_training",
        "label_family": "activity_labels",
        "default_claim": "Atlas ranks curated ChEMBL drug-target bioactivity where assay data exist.",
        "reference": "https://www.ebi.ac.uk/chembl/",
    },
    "bindingdb": {
        "role": "direct_binding_bioactivity_training",
        "label_family": "activity_labels",
        "default_claim": "Atlas ranks experimentally measured BindingDB drug-target affinities where mapped pairs exist.",
        "reference": "https://www.bindingdb.org/",
    },
    "toxcast": {
        "role": "independent_in_vitro_activity_benchmark",
        "label_family": "activity_labels",
        "default_claim": "Atlas recovers experimentally active chemical-target assay pairs beyond SPD.",
        "reference": "https://www.epa.gov/comptox-tools/exploring-toxcast-data",
    },
    "tox21": {
        "role": "independent_in_vitro_toxicology_benchmark",
        "label_family": "activity_labels",
        "default_claim": "Atlas recovers Tox21/ToxCast toxicology-relevant in vitro activity after source-specific QC.",
        "reference": "https://www.epa.gov/comptox-tools/exploring-toxcast-data",
    },
    "pubchem_bioassay": {
        "role": "high_volume_screening_activity_training",
        "label_family": "activity_labels",
        "default_claim": "Atlas ranks PubChem BioAssay active/inactive outcomes after assay-level QC.",
        "reference": "https://pubchem.ncbi.nlm.nih.gov/docs/bioassays",
    },
    "iuphar_gtopdb": {
        "role": "expert_curated_pharmacology_positive_evidence",
        "label_family": "activity_labels",
        "default_claim": "Atlas recovers high-confidence IUPHAR/BPS Guide to Pharmacology ligand-target positives.",
        "reference": "https://www.guidetopharmacology.org/download.jsp",
    },
    "ncats_inxight": {
        "role": "pk_safety_context",
        "label_family": "exposure_context,safety_metadata",
        "default_claim": "Atlas exposure-aware analyses use NCATS Inxight PK/safety metadata as context without turning missing data into negatives.",
        "reference": "https://drugs.ncats.io/downloads-public",
    },
    "sider": {
        "role": "adr_mechanism_evidence",
        "label_family": "mechanism_evidence",
        "default_claim": "Atlas-prioritized pairs are enriched for drug ADR mechanisms.",
        "reference": "http://sideeffects.embl.de/",
    },
    "ctd": {
        "role": "chemical_gene_disease_mechanism_evidence",
        "label_family": "mechanism_evidence",
        "default_claim": "Atlas mechanisms connect chemical-gene-disease evidence.",
        "reference": "https://ctdbase.org/",
    },
    "opentargets_safety": {
        "role": "target_safety_evidence",
        "label_family": "mechanism_evidence",
        "default_claim": "Atlas mechanisms are supported by curated target safety liabilities.",
        "reference": "https://platform-docs.opentargets.org/target/safety",
    },
    "reactome": {
        "role": "pathway_mechanism_evidence",
        "label_family": "mechanism_evidence",
        "default_claim": "Atlas mechanisms are interpretable through target-pathway-ADR links.",
        "reference": "https://reactome.org/dev/content-service/",
    },
}


def source_registry_from_config(config: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    configured = config.get("external_sources")
    registry: dict[str, dict[str, Any]] = {}
    if isinstance(configured, Mapping):
        for name, raw_spec in configured.items():
            spec = dict(raw_spec) if isinstance(raw_spec, Mapping) else {}
            registry[str(name)] = {**SOURCE_ROLES.get(str(name), {}), **spec}
    raw_inputs = config.get("inputs")
    inputs: Mapping[str, Any] = raw_inputs if isinstance(raw_inputs, Mapping) else {}
    legacy_paths = {
        "spd": "spd_file",
        "toxcast": "toxcast_file",
        "papyrus": "papyrus_file",
        "chembl": "chembl_file",
        "bindingdb": "bindingdb_file",
        "sider": "sider_file",
        "ctd": "ctd_file",
        "opentargets_safety": "opentargets_safety_file",
        "reactome": "reactome_file",
        "pubchem_bioassay": "pubchem_bioassay_file",
        "tox21": "tox21_file",
        "iuphar_gtopdb": "iuphar_gtopdb_file",
        "ncats_inxight": "ncats_inxight_file",
    }
    for name, key in legacy_paths.items():
        if name not in registry and inputs.get(key):
            registry[name] = {**SOURCE_ROLES.get(name, {}), "path": inputs[key], "enabled": True}
    return registry


def write_registry_manifest(config: Mapping[str, Any], out_path: str | Path) -> None:
    write_source_manifest(out_path, {"sources": source_registry_from_config(config)})


def refresh_configured_sources(config: Mapping[str, Any], cache_dir: str | Path, *, overwrite: bool = False) -> dict[str, str]:
    downloaded: dict[str, str] = {}
    cache = Path(cache_dir)
    for name, spec in source_registry_from_config(config).items():
        if spec.get("enabled") is False or not spec.get("url"):
            continue
        filename = str(spec.get("cache_name") or Path(str(spec["url"])).name or f"{name}.dat")
        out = download_to_cache(str(spec["url"]), cache / filename, overwrite=overwrite)
        downloaded[name] = str(out)
    return downloaded
