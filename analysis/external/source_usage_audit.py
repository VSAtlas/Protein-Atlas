from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def _classify(path: Path) -> tuple[str, str]:
    text = str(path)
    name = path.name
    if "/raw_extracted/" in text:
        return "raw_source_of_normalized_training_table", "Use normalized cardiac_ion_channel_bioactivity.tsv instead."
    if "/raw_mechanism_sources/" in text:
        return "raw_source_of_normalized_mechanism_table", "Use normalized SIDER/Reactome tables instead."
    if text.endswith("/benchmarks/muv/muv.csv"):
        return "raw_source_of_benchmark_table", "Use normalized MUV benchmark long table."
    if "/benchmarks/muv/muv_benchmark_long.tsv" in text:
        return "benchmark_only", "Included in four-state evidence as benchmark-only; never production truth."
    if "/banana/" in text:
        return "runtime_or_benchmark_support", "BANANA files support binding expert runtime/benchmarking, not labels."
    if "spd_overlap_audit" in text or name.endswith(".cid_properties.csv") or name.endswith(".target_audit.csv"):
        return "audit_support", "Audit/mapping support; not a label source."
    if name.endswith(".chembl_atlas_mapping.csv") or "ml_source_catalog_manifest" in text:
        return "mapping_or_manifest_support", "Consumed as source lineage or mapping support."
    if name.startswith("AID_") and name.endswith(".normalized.tsv"):
        return "production_candidate_drug_target_activity", "Included through PubChem BioAssay four-state input."
    if name.startswith("AID_") and name.endswith((".csv", ".txt")):
        return "raw_source_of_normalized_training_table", "Use the normalized AID table."
    if name == "pubchem.chembl.dataset4publication_inchi_smiles.tsv":
        return "raw_source_of_normalized_training_table", "Staged into data/external/excape/bioactivity.tsv."
    if text.endswith("/excape/bioactivity.tsv"):
        return "normalized_training_label_source", "Included in four-state evidence through --excape."
    if text.endswith("/bindingdb/BindingDB_All_202605_tsv.zip"):
        return "raw_source_of_normalized_training_table", "Staged into data/external/bindingdb/bioactivity.tsv."
    if text.endswith("/bindingdb/bioactivity.tsv"):
        return "normalized_training_label_source", "Included in four-state evidence through --bindingdb."
    if "/tox21/" in text and name == "bioactivity.tsv":
        return "normalized_training_label_source", "Included in four-state evidence through --tox21."
    if "/iuphar_gtopdb/" in text and name in {"bioactivity.tsv", "interactions.tsv", "approved_drug_interactions.tsv"}:
        return "normalized_training_label_source", "Included in four-state evidence through --iuphar-gtopdb; absence is never a negative."
    if "/ncats_inxight/" in text and name in {
        "frdb-drugs.tsv",
        "frdb-pk.tsv",
        "frdb-toxicity.tsv",
        "frdb-adverseevents.tsv",
        "frdb-ddi.tsv",
        "pk_table.tsv",
    }:
        return "pk_safety_context_source", "Use through build_pk_table or ADR/safety-context staging; not direct drug-target truth."
    if text.endswith("/adrecs_target/target_adr_evidence.tsv"):
        return "normalized_mechanism_label_source", "Included as target-ADR positive mechanism evidence."
    if text.endswith("/adrecs_target/drug_target_adr_evidence.tsv"):
        return "normalized_mechanism_label_source", "Included as drug-target-ADR positive mechanism evidence."
    if "/adrecs_target/" in text and name.endswith(".xlsx"):
        return "raw_source_of_normalized_mechanism_table", "Use normalized adrecs_target/target_adr_evidence.tsv."
    production = {
        "bioactivity.tsv",
        "kiba_bioactivity.tsv",
        "cardiac_ion_channel_bioactivity.tsv",
        "pubchem_bioassay.tsv",
        "ohdsi_omop_faers_disproportionality.csv",
        "omopReferenceSet.csv",
        "ohdsiNegativeControls.csv",
        "ohdsiDevelopmentNegativeControls.csv",
        "mitotox_mechanism_evidence.tsv",
        "ctd_edges.tsv",
        "drug_adr.tsv",
        "safety.tsv",
        "pathways.tsv",
    }
    if name in production:
        return "production_or_validation_label_source", "Included in activity or mechanism label/graph pipelines."
    if name in {
        "mitotox_atlas_mapped_activity.tsv",
        "mitotox_records.csv",
        "mitotox_targets.csv",
        "mitotox_compounds.csv",
    }:
        return "raw_or_mapped_mechanism_support", "Consumed by MitoTox mechanism staging; not direct binding truth."
    if name == "P1-09-Target_compound_activity.txt":
        return "production_candidate_drug_target_activity", "Staged into data/external/ttd/ttd_activity.tsv."
    if text.endswith("/ttd/ttd_activity.tsv"):
        return "normalized_training_label_source", "Included in four-state evidence through --ttd-activity."
    if name in {
        "P1-01-TTD_target_download.txt",
        "P1-03-TTD_crossmatching.txt",
        "P2-01-TTD_uniprot_all.txt",
        "P3-07-Approved_smi_inchi.txt",
        "P1-07-Drug-TargetMapping(Sheet1).csv",
    }:
        return "mapping_or_manifest_support", "Consumed by TTD staging or reserved for known-target metadata."
    if name in {"P1-02-TTD_drug_download.txt", "P1-05-Drug_disease.txt", "P1-06-Target_disease.txt"}:
        return "mechanism_or_metadata_source", "Available for mechanism/metadata extension; not direct binding labels."
    if name == "drug.target.interaction.tsv":
        return "production_candidate_drug_target_activity", "Staged into data/external/drugcentral/drugcentral_activity.tsv."
    if text.endswith("/drugcentral/drugcentral_activity.tsv"):
        return "normalized_training_label_source", "Included in four-state evidence through --drugcentral."
    if name == "evebio_computational_result_doc.csv":
        return "documentation_only", "Documentation file only; no measured label rows."
    return "unclassified", "Review required before training use."


def audit_external_source_usage(external_dir: str | Path, out_dir: str | Path) -> pd.DataFrame:
    root = Path(external_dir)
    rows: list[dict[str, str | int]] = []
    for path in sorted(root.rglob("*")):
        if path.suffix.lower() not in {".csv", ".tsv", ".txt"} or not path.is_file():
            continue
        usage, note = _classify(path)
        rows.append(
            {
                "path": str(path),
                "relative_path": str(path.relative_to(root)),
                "bytes": int(path.stat().st_size),
                "usage_class": usage,
                "note": note,
            }
        )
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    table = out / "external_source_usage_manifest.csv"
    frame.to_csv(table, index=False)
    manifest = {
        "external_dir": str(root),
        "n_delimited_files": int(len(frame)),
        "usage_counts": {str(k): int(v) for k, v in frame["usage_class"].value_counts().items()} if not frame.empty else {},
        "unclassified": frame.loc[frame["usage_class"].eq("unclassified"), "relative_path"].tolist() if not frame.empty else [],
        "table": str(table),
    }
    (out / "external_source_usage_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return frame
