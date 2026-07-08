from __future__ import annotations

import json
import zipfile
from dataclasses import asdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd


CANONICAL_EVIDENCE_FIELDS = [
    "relation_domain",
    "canonical_pair_key",
    "compound_id",
    "compound_namespace",
    "inchikey",
    "smiles",
    "drug_name",
    "target_id",
    "target_namespace",
    "uniprot",
    "gene_symbol",
    "target_species",
    "target_variant",
    "adr_id",
    "adr_namespace",
    "adr_term",
    "meddra_level",
    "organ_system",
    "assay_id",
    "assay_type",
    "endpoint_type",
    "assay_mode",
    "tissue_or_cell_context",
    "raw_call",
    "raw_value",
    "raw_relation",
    "raw_unit",
    "standardized_value_nm",
    "exposure_metric",
    "cmax_total",
    "cmax_free",
    "cave_free",
    "css_free",
    "auc",
    "exposure_margin",
    "source_name",
    "source_family",
    "source_version",
    "source_row_id",
    "upstream_source",
    "evidence_type",
    "provenance_priority",
    "quality_flag",
    "license_allows_ml_training",
    "redistribution_allowed",
    "commercial_restriction",
    "benchmark_only",
    "validation_only",
    "training_allowed",
    "label_state",
    "label_confidence",
    "exclude_reason",
]


CONFLICT_RULES = [
    "Collapse by canonical pair plus endpoint family; do not collapse across assay mode, species, tissue, target variant, endpoint type, or exposure metric.",
    "Measured assay positives override decoy, mined-negative, benchmark, and weak non-signal evidence.",
    "Measured assay negatives override benchmark-only assumptions.",
    "Comparable measured active/inactive conflicts become label_state=-1, not a majority vote.",
    "Positive-only curation sources create 1/NaN only, never 0.",
    "Absence of evidence remains NaN.",
    "Benchmark-only decoys stay in benchmark_decoy or excluded state; they never create production negatives.",
    "Lower-priority derived wrappers append provenance but do not create independent training rows when upstream primary evidence is already loaded.",
]


@dataclass(frozen=True)
class SourceSpec:
    source_name: str
    source_family: str
    priority: int
    relation_domains: tuple[str, ...]
    expert_feeds: tuple[str, ...]
    default_use: str
    label_policy: str
    upstream_source: str
    access_path: str
    local_path: str
    license_status: str
    inactive_status: str
    benchmark_only: bool = False
    validation_only: bool = False
    training_allowed_after_verification: bool = False


SOURCE_SPECS = [
    SourceSpec(
        "SPD",
        "panel_complete_secondary_pharmacology",
        1,
        ("drug_target_activity", "drug_target_exposure"),
        ("binding_activity", "exposure_relevance"),
        "training_validation",
        "Use assay-specific thresholds; active measured rows -> 1, explicit inactive/no activity rows -> 0, gray zone -> NaN, conflicts -> -1.",
        "Novartis Secondary Pharmacology Database",
        "local export or publication supplement",
        "data/external/spd/sutherland_2023_spd_supplementary_data_1_15.xlsx",
        "publication supplement; license/use must be confirmed before redistribution",
        "explicit measured inactive where panel row states inactive/no activity",
        training_allowed_after_verification=True,
    ),
    SourceSpec(
        "EvE Bio",
        "panel_complete_safety_pharmacology",
        1,
        ("drug_target_activity",),
        ("binding_activity",),
        "candidate_training_validation",
        "Active rows -> 1, explicit inactive rows -> 0, artifact/cytotoxicity/high-frequency-hitter flags -> -1.",
        "EvE Bio panel reports",
        "requires manual access confirmation",
        "data/external/eve_bio/eve_bio_assays.tsv",
        "unknown",
        "unknown until raw panel headers are inspected",
    ),
    SourceSpec(
        "SafetyScreen44/SAFETYscan47",
        "panel_complete_safety_pharmacology",
        1,
        ("drug_target_activity",),
        ("binding_activity",),
        "candidate_training_validation",
        "Use only raw full-panel reports; active measured rows -> 1, printed inactive rows -> 0, narrative absence -> NaN.",
        "Eurofins/Cerep-style safety panels",
        "requires raw report access",
        "data/external/safetyscan/safetyscan.tsv",
        "unknown",
        "only full-panel printed inactive rows qualify",
    ),
    SourceSpec(
        "BindingDB",
        "direct_measured_bioactivity",
        2,
        ("drug_target_activity",),
        ("binding_activity",),
        "training",
        "Prefer Kd/Ki over IC50; <=1 uM -> 1, >=10 uM measured/censored inactive -> 0, 1-10 uM -> NaN, non-direct/failed target map -> -1.",
        "BindingDB",
        "https://www.bindingdb.org/rwd/bind/downloads/",
        "data/external/bindingdb/BindingDB_All_202605_tsv.zip",
        "requires local file/header/license verification",
        "measured or curated inactive only",
        training_allowed_after_verification=True,
    ),
    SourceSpec(
        "ChEMBL",
        "direct_measured_bioactivity",
        2,
        ("drug_target_activity",),
        ("binding_activity",),
        "training_validation",
        "Single-target high-confidence measured rows only; pChEMBL>=6 or <=1 uM -> 1, pChEMBL<=5 or >=10 uM measured inactive -> 0, gray zone -> NaN.",
        "ChEMBL",
        "https://www.ebi.ac.uk/chembl/",
        "data/external/chembl/bioactivity.tsv",
        "local staged export; ChEMBL licensing must be cited from source release",
        "measured inactive only; missing rows are unknown",
        training_allowed_after_verification=True,
    ),
    SourceSpec(
        "PubChem BioAssay",
        "direct_measured_bioactivity",
        2,
        ("drug_target_activity",),
        ("binding_activity",),
        "training_validation",
        "Confirmatory single-target assays only; Active -> 1, Inactive -> 0, Inconclusive -> NaN, unspecified/probe-only -> -1.",
        "PubChem BioAssay",
        "https://pubchem.ncbi.nlm.nih.gov/docs/bioassays",
        "data/external/pubchem_bioassay/aid_1963048/AID_1963048.normalized.tsv",
        "targeted PUG-REST export; PubChem source provenance belongs to submitters",
        "explicit assay inactive calls",
        training_allowed_after_verification=True,
    ),
    SourceSpec(
        "ToxCast/Tox21",
        "direct_hts_toxicology_activity",
        2,
        ("drug_target_activity",),
        ("binding_activity", "mechanism_adr"),
        "training_validation_source_holdout",
        "Use assay-QC-passed endpoint hit calls only; active/hit -> 1, explicit inactive/non-hit -> 0, inconclusive/cytotoxic/artifact/conflict -> -1 or NaN; missing target coverage is unknown.",
        "EPA invitrodb/ToxCast/Tox21",
        "https://www.epa.gov/comptox-tools/exploring-toxcast-data",
        "data/external/toxcast/bioactivity.tsv",
        "EPA invitrodb local staged export; cite release and assay-QC filters",
        "explicit assay inactive/non-hit calls only after assay and chemical QC",
        training_allowed_after_verification=True,
    ),
    SourceSpec(
        "IUPHAR/BPS Guide to Pharmacology",
        "expert_curated_pharmacology",
        2,
        ("drug_target_activity",),
        ("binding_activity",),
        "high_precision_positive_training",
        "Expert-curated ligand-target interactions create high-confidence positives when quantitative affinity or curated interaction evidence supports activity; absence and weak affinity are never negatives.",
        "IUPHAR/BPS Guide to Pharmacology",
        "https://www.guidetopharmacology.org/download.jsp",
        "data/external/iuphar_gtopdb/bioactivity.tsv",
        "ODbL database with CC BY-SA 4.0 contents; cite GtoPdb release",
        "positive-biased curated source; only explicit inactive/no-activity rows can become negatives",
        training_allowed_after_verification=True,
    ),
    SourceSpec(
        "TTD non-binders",
        "dense_family_matrix",
        3,
        ("drug_target_activity",),
        ("binding_activity",),
        "candidate_training_validation",
        "Non-binders >200 uM -> 0; poor binders 50-200 uM -> weak 0 or NaN depending endpoint; preserve subclass.",
        "Therapeutic Target Database",
        "requires current access confirmation",
        "data/external/ttd/ttd_activity.tsv",
        "TTD 10.1 local files staged from data/external/P1-09 plus mapping support; cite TTD terms",
        "measured target-compound rows with conservative thresholds; not absence-based",
        training_allowed_after_verification=True,
    ),
    SourceSpec(
        "GLASS2/PDSP/Davis/Karaman",
        "dense_family_matrix",
        3,
        ("drug_target_activity",),
        ("binding_activity",),
        "family_specific_training_validation",
        "Human-target measured rows with conservative inactivity thresholds; kinase panels must carry target-family flags.",
        "family-specific measured matrices",
        "requires source-specific local exports",
        "data/external/family_panels/kiba_bioactivity.tsv",
        "Zenodo dataset; cite Tang et al. KIBA source; normalized from data/external/family_panels/kiba.zip",
        "measured matrix inactives only",
    ),
    SourceSpec(
        "Cardiac exposure-liability panels",
        "exposure_linked_measured",
        4,
        ("drug_target_exposure", "organ_toxicity"),
        ("exposure_relevance", "mechanism_adr"),
        "validation_calibration",
        "Preserve endpoint family; do not mix Redfern/Kramer/Gintant/CiPA/HESI/FDA families across train/test.",
        "Jansson-Lofmark/Tsuchitani/Kotani/Goldstein/Liu/Redfern/Kramer/CiPA/Gintant",
        "requires paper/table-specific verification",
        "data/external/cardiac_safety/cardiac_ion_channel_bioactivity.tsv",
        "Zenodo cardiac ion-channel dataset; cite source DOI and upstream composition; normalized from cardiac_ion_channels_raw.rar",
        "measured inactive/control rows only where printed",
        validation_only=True,
    ),
    SourceSpec(
        "ADReCS-Target/MitoTox/DITOP/DART/TPDB",
        "curated_mechanism_positive",
        5,
        ("target_adr_mechanism", "organ_toxicity"),
        ("mechanism_adr",),
        "enrichment_training_positive_only",
        "Target/organ toxicity positives only; no negatives.",
        "curated safety mechanism resources",
        "requires source-specific download verification",
        "data/external/adrecs_target/adr_protein_associations.xlsx",
        "ADReCS-Target page states CC BY-NC-SA 4.0 academic use",
        "positive-only source",
        training_allowed_after_verification=True,
    ),
    SourceSpec(
        "MitoTox",
        "curated_mechanism_positive",
        5,
        ("drug_target_adr_mechanism", "organ_toxicity"),
        ("mechanism_adr",),
        "sensitivity_training_positive_negative",
        "MitoTox result=true creates mitochondrial mechanism positive evidence; result=false is measured context-specific negative evidence, not a global ADR negative.",
        "MitoTox REST API",
        "https://www.mitotox.org/api/",
        "data/external/mitotox/mitotox_mechanism_evidence.tsv",
        "CC BY-NC 4.0 per MitoTox API page",
        "context-specific measured negative where result=false",
        training_allowed_after_verification=True,
    ),
    SourceSpec(
        "AOP-Wiki/AOP-DB",
        "curated_pathway_mechanism",
        5,
        ("pathway_adr",),
        ("mechanism_adr",),
        "enrichment",
        "Pathway-ADR positives; direct drug-target-ADR remains NaN unless independently grounded.",
        "AOP-Wiki/AOP-DB",
        "requires current dump/API verification",
        "data/external/reactome/pathways.tsv",
        "Reactome/curated pathway local export; source license must be cited separately",
        "positive-only source",
    ),
    SourceSpec(
        "OMOP/EU-ADR/OHDSI/Harpaz",
        "clinical_drug_adr_controls",
        6,
        ("drug_adr_clinical",),
        ("mechanism_adr",),
        "calibration_sensitivity",
        "Clinical drug-ADR controls only; do not back-project into target mechanism labels.",
        "OMOP/OHDSI reference sets and related controls",
        "requires local reference-set verification",
        "data/external/ohdsi/omopReferenceSet.csv",
        "local OHDSI/OMOP reference export; source license must be cited separately",
        "clinical negative controls; not biochemical negatives",
        validation_only=True,
    ),
    SourceSpec(
        "Papyrus/ExCAPE/DTC/DrugCentral",
        "derived_wrapper",
        7,
        ("drug_target_activity",),
        ("binding_activity",),
        "secondary_training_enrichment",
        "Use only after upstream lineage is recorded; do not double-count if ChEMBL/PubChem/BindingDB primary evidence is loaded.",
        "aggregated databases",
        "requires source-lineage/header verification",
        "data/external/drugcentral/drugcentral_activity.tsv",
        "local staged export; parent-source lineage must be preserved",
        "depends on upstream source; high overlap risk",
        training_allowed_after_verification=True,
    ),
    SourceSpec(
        "NCATS Inxight Drugs",
        "drug_pk_safety_metadata",
        6,
        ("drug_pk_context", "drug_safety_context", "drug_adr_clinical"),
        ("exposure_relevance", "mechanism_adr"),
        "exposure_context_validation",
        "PK, dose, toxicity, adverse-event, DDI, and drug-status metadata provide exposure/safety context only; do not create drug-target activity labels or negatives from missing metadata.",
        "NCATS Inxight Drugs / Fast Response Database",
        "https://drugs.ncats.io/downloads-public",
        "data/external/ncats_inxight/frdb-v2024-12-30.zip",
        "public-domain NCATS facts; cite Inxight/FRDB release",
        "not a biochemical negative source",
        validation_only=True,
    ),
    SourceSpec(
        "T-ARDIS/ADRtarget/Ietswaart/LAERTES/CEM",
        "weak_hybrid_mined",
        8,
        ("drug_adr_clinical", "target_adr_mechanism", "drug_target_adr_mechanism"),
        ("mechanism_adr",),
        "sensitivity_enrichment",
        "Weak/hybrid/mined evidence only; do not train/evaluate against FAERS/SIDER/OFFSIDES-like labels without source-family holdout.",
        "FAERS/SIDER/literature-mined hybrids",
        "requires source-specific verification",
        "data/external/weak_hybrid/",
        "unknown",
        "weak non-signal only; not hard negatives",
        validation_only=True,
    ),
    SourceSpec(
        "DUD-E/DEKOIS/MUV/decoy BigBind",
        "benchmark_decoy",
        9,
        ("benchmark_decoy",),
        ("binding_activity",),
        "benchmark_only",
        "Actives/decoys remain benchmark namespace or excluded conflict; decoys never become production negatives.",
        "benchmark/decoy sets",
        "requires benchmark-specific local files",
        "data/external/benchmarks/muv/muv_benchmark_long.tsv",
        "benchmark-only local archive; source-specific license required",
        "decoy/inferred benchmark-only",
        benchmark_only=True,
    ),
    SourceSpec(
        "SIDER",
        "clinical_drug_adr_positive",
        5,
        ("drug_adr_clinical",),
        ("mechanism_adr",),
        "enrichment_positive_only",
        "Drug-ADR positives from package-insert labels; no negatives.",
        "SIDER",
        "http://sideeffects.embl.de/download/",
        "data/external/raw_mechanism_sources/sider/meddra_all_se.tsv.gz",
        "local staged export; cite SIDER source terms",
        "positive-only source",
        validation_only=True,
    ),
    SourceSpec(
        "CTD",
        "curated_chemical_gene_disease",
        5,
        ("drug_adr_clinical", "target_adr_mechanism", "drug_target_adr_mechanism"),
        ("mechanism_adr",),
        "enrichment_positive_only",
        "Curated chemical-gene-disease edges; positives only unless CTD explicitly marks a negative/control context.",
        "Comparative Toxicogenomics Database",
        "https://ctdbase.org/downloads/",
        "data/external/ctd/ctd_edges.tsv",
        "local staged export; cite CTD source terms",
        "positive-only source",
        validation_only=True,
    ),
    SourceSpec(
        "OpenTargets Safety",
        "curated_target_safety",
        5,
        ("target_adr_mechanism",),
        ("mechanism_adr",),
        "enrichment_positive_only",
        "Target safety positives/evidence only; missing target safety evidence remains unknown.",
        "Open Targets Platform safety evidence",
        "https://platform.opentargets.org/",
        "data/external/opentargets/safety.tsv",
        "local staged export; cite Open Targets source terms",
        "positive-only source",
        validation_only=True,
    ),
    SourceSpec(
        "Reactome",
        "curated_pathway_mechanism",
        5,
        ("pathway_adr",),
        ("mechanism_adr",),
        "enrichment_positive_only",
        "Target-pathway edges only; pathway association is not a direct ADR negative or positive without ADR evidence.",
        "Reactome",
        "https://reactome.org/download-data",
        "data/external/reactome/pathways.tsv",
        "local staged export; cite Reactome source terms",
        "positive-only source",
        validation_only=True,
    ),
    SourceSpec(
        "FAERS non-signals",
        "weak_clinical_nonsignal",
        6,
        ("drug_adr_clinical",),
        ("mechanism_adr",),
        "calibration_sensitivity",
        "Adequately observed FAERS non-signals can be weak controls only; FAERS absence remains NaN.",
        "openFDA/FAERS disproportionality cache",
        "repo-local openFDA-derived cache",
        "data/external/faers/ohdsi_omop_faers_disproportionality.csv",
        "openFDA-derived local cache; cite FDA/openFDA caveats",
        "weak non-signal only; not a hard causal negative",
        validation_only=True,
    ),
]


def _inspect_header(path: Path) -> tuple[bool, list[str]]:
    if not path.exists():
        return False, []
    if path.is_dir():
        for child in sorted(path.iterdir()):
            exists, header = _inspect_header(child)
            if exists and header:
                return True, header
        return True, []
    try:
        suffixes = "".join(path.suffixes).lower()
        if path.suffix.lower() in {".xlsx", ".xls"}:
            excel = pd.ExcelFile(path)
            if not excel.sheet_names:
                return True, []
            return True, list(pd.read_excel(path, sheet_name=excel.sheet_names[0], nrows=0).columns)
        if path.suffix.lower() == ".zip":
            with zipfile.ZipFile(path) as archive:
                members = [
                    name
                    for name in archive.namelist()
                    if not name.endswith("/")
                    and Path(name).suffix.lower() in {".csv", ".tsv", ".txt", ".xlsx", ".xls"}
                ]
                if not members:
                    return True, []
                if Path(members[0]).suffix.lower() in {".xlsx", ".xls"}:
                    with archive.open(members[0]) as handle:
                        return True, list(pd.read_excel(handle, nrows=0).columns)
                with archive.open(members[0]) as handle:
                    first = handle.readline().decode("utf-8", errors="replace").strip()
                sep = "\t" if Path(members[0]).suffix.lower() in {".tsv", ".txt"} else ","
                return True, [field.strip() for field in first.split(sep) if field.strip()]
        if any(suffixes.endswith(ext) for ext in (".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz", ".csv.xz", ".tsv.xz", ".txt.xz")):
            sep = "\t" if ".tsv" in suffixes or ".txt" in suffixes else ","
            return True, list(pd.read_csv(path, sep=sep, nrows=0, compression="infer").columns)
    except Exception:
        return True, []
    return True, []


def build_source_label_spec(out_dir: str | Path) -> dict[str, pd.DataFrame]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for spec in SOURCE_SPECS:
        row = asdict(spec)
        exists, header = _inspect_header(Path(spec.local_path))
        row["local_exists"] = exists
        row["header_fields_json"] = json.dumps(header)
        row["header_verified"] = bool(header)
        row["verification_status"] = "header_verified" if header else "requires_local_file_header_license_verification"
        row["training_allowed"] = bool(spec.training_allowed_after_verification and header and not spec.benchmark_only and not spec.validation_only)
        row["blocked_reason"] = "" if row["training_allowed"] else "not training-enabled until local file/header/license/upstream lineage are verified"
        rows.append(row)
    sources = pd.DataFrame(rows)
    mapping = sources[
        [
            "source_name",
            "source_family",
            "priority",
            "relation_domains",
            "expert_feeds",
            "default_use",
            "label_policy",
            "training_allowed",
            "validation_only",
            "benchmark_only",
            "blocked_reason",
        ]
    ].copy()
    mapping["relation_domains"] = mapping["relation_domains"].map(";".join)
    mapping["expert_feeds"] = mapping["expert_feeds"].map(";".join)
    rules = pd.DataFrame({"rule_order": range(1, len(CONFLICT_RULES) + 1), "rule": CONFLICT_RULES})
    schema = pd.DataFrame({"field": CANONICAL_EVIDENCE_FIELDS})
    do_not = sources.loc[
        sources["benchmark_only"].eq(True)
        | sources["verification_status"].ne("header_verified")
        | sources["training_allowed"].eq(False),
        ["source_name", "default_use", "verification_status", "blocked_reason", "label_policy"],
    ].copy()
    outputs = {
        "verified_deduplicated_source_table": sources,
        "source_to_label_mapping": mapping,
        "source_priority_conflict_rules": rules,
        "minimum_ingestion_schema": schema,
        "do_not_ingest_yet": do_not,
    }
    for name, frame in outputs.items():
        frame.to_csv(out / f"{name}.csv", index=False)
    manifest = {
        "outputs": {name: str(out / f"{name}.csv") for name in outputs},
        "label_states": {"1": "strict positive", "0": "measured/reliable negative", "NaN": "unknown/gray-zone", "-1": "excluded/ambiguous/conflicting"},
        "absolute_rules": CONFLICT_RULES,
        "warning": "This exporter specifies source policy and verifies local headers when present; it does not certify a source as production-ready without license and upstream-lineage review.",
    }
    (out / "source_label_spec_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return outputs
