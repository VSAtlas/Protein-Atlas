from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.four_state_evidence import (
    activity_source_to_evidence,
    collapse_four_state_evidence,
    drug_adr_control_source_to_evidence,
    drug_adr_positive_source_to_evidence,
    faers_nonsignal_source_to_evidence,
    join_collapsed_labels,
    negative_evidence_to_raw,
    spd_source_to_evidence,
    target_adr_source_to_evidence,
    write_four_state_outputs,
)


def _extend(parts: list[pd.DataFrame], paths: list[Path], **kwargs: Any) -> None:
    for path in paths:
        if path.exists():
            parts.append(activity_source_to_evidence(path, **kwargs))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build raw and collapsed four-state Atlas evidence labels."
    )
    parser.add_argument("--chembl", nargs="*", type=Path, default=[])
    parser.add_argument("--bindingdb", nargs="*", type=Path, default=[])
    parser.add_argument("--pubchem", nargs="*", type=Path, default=[])
    parser.add_argument("--iuphar-gtopdb", "--gtopdb", nargs="*", type=Path, default=[])
    parser.add_argument("--spd", nargs="*", type=Path, default=[])
    parser.add_argument("--papyrus", nargs="*", type=Path, default=[])
    parser.add_argument("--excape", nargs="*", type=Path, default=[])
    parser.add_argument("--ttd-activity", nargs="*", type=Path, default=[])
    parser.add_argument("--drugcentral", nargs="*", type=Path, default=[])
    parser.add_argument("--toxcast", nargs="*", type=Path, default=[])
    parser.add_argument("--tox21", nargs="*", type=Path, default=[])
    parser.add_argument("--family-panel", nargs="*", type=Path, default=[])
    parser.add_argument("--cardiac-safety", nargs="*", type=Path, default=[])
    parser.add_argument("--omop-ohdsi", nargs="*", type=Path, default=[])
    parser.add_argument("--drug-adr-positive", nargs="*", type=Path, default=[])
    parser.add_argument("--faers", nargs="*", type=Path, default=[])
    parser.add_argument("--target-adr", nargs="*", type=Path, default=[])
    parser.add_argument("--lit-pcba", nargs="*", type=Path, default=[])
    parser.add_argument("--dude", nargs="*", type=Path, default=[])
    parser.add_argument("--muv", nargs="*", type=Path, default=[])
    parser.add_argument("--negative-evidence", nargs="*", type=Path, default=[])
    parser.add_argument("--mapping", type=Path, default=None)
    parser.add_argument("--feature-table", type=Path, default=None)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args(argv)

    parts: list[pd.DataFrame] = []
    _extend(
        parts,
        args.chembl,
        source_name="ChEMBL",
        source_priority=1,
        source_role="direct_assay",
        mapping_path=args.mapping,
        confidence=0.9,
    )
    _extend(
        parts,
        args.bindingdb,
        source_name="BindingDB",
        source_priority=1,
        source_role="direct_binding_assay",
        mapping_path=args.mapping,
        confidence=0.92,
    )
    _extend(
        parts,
        args.pubchem,
        source_name="PubChem BioAssay",
        source_priority=1,
        source_role="direct_assay_outcome",
        mapping_path=args.mapping,
        confidence=0.85,
        use_numeric_thresholds=False,
        use_pchembl=False,
    )
    _extend(
        parts,
        args.iuphar_gtopdb,
        source_name="IUPHAR/BPS Guide to Pharmacology",
        source_priority=1,
        source_role="expert_curated_ligand_target_pharmacology",
        mapping_path=args.mapping,
        confidence=0.95,
        positive_threshold_nM=1000.0,
        negative_threshold_nM=10000.0,
        use_pchembl=True,
        allow_numeric_negatives=False,
    )
    for path in args.spd:
        if path.exists():
            parts.append(spd_source_to_evidence(path, mapping_path=args.mapping))
    _extend(
        parts,
        args.papyrus,
        source_name="Papyrus",
        source_priority=4,
        source_role="aggregated_training_secondary",
        mapping_path=args.mapping,
        confidence=0.75,
    )
    _extend(
        parts,
        args.excape,
        source_name="ExCAPE-DB",
        source_priority=4,
        source_role="aggregated_training_only",
        mapping_path=args.mapping,
        confidence=0.7,
    )
    _extend(
        parts,
        args.ttd_activity,
        source_name="TTD",
        source_priority=2,
        source_role="direct_measured_target_compound_activity",
        mapping_path=args.mapping,
        confidence=0.85,
    )
    _extend(
        parts,
        args.drugcentral,
        source_name="DrugCentral",
        source_priority=4,
        source_role="derived_wrapper_target_interaction",
        mapping_path=args.mapping,
        confidence=0.72,
    )
    _extend(
        parts,
        args.toxcast,
        source_name="ToxCast",
        source_priority=1,
        source_role="direct_hts_assay",
        mapping_path=args.mapping,
        confidence=0.85,
    )
    _extend(
        parts,
        args.tox21,
        source_name="Tox21",
        source_priority=1,
        source_role="direct_hts_toxicology_assay",
        mapping_path=args.mapping,
        confidence=0.82,
        use_numeric_thresholds=False,
        use_pchembl=False,
    )
    _extend(
        parts,
        args.family_panel,
        source_name="KIBA/family panel",
        source_priority=2,
        source_role="family_specific_measured_matrix",
        mapping_path=args.mapping,
        confidence=0.8,
    )
    _extend(
        parts,
        args.cardiac_safety,
        source_name="Cardiac ion-channel panel",
        source_priority=2,
        source_role="cardiac_safety_measured_panel",
        mapping_path=args.mapping,
        confidence=0.85,
    )
    for path in args.omop_ohdsi:
        if path.exists():
            parts.append(
                drug_adr_control_source_to_evidence(
                    path,
                    source_name="OMOP/OHDSI",
                    source_priority=2,
                    source_role="curated_clinical_drug_outcome_control",
                    mapping_path=args.mapping,
                    confidence=0.65,
                )
            )
    for path in args.drug_adr_positive:
        if path.exists():
            parts.append(
                drug_adr_positive_source_to_evidence(
                    path,
                    source_name="SIDER/openFDA drug-ADR positive",
                    source_priority=2,
                    source_role="positive_only_clinical_drug_adr_label",
                    mapping_path=args.mapping,
                    confidence=0.7,
                )
            )
    for path in args.faers:
        if path.exists():
            parts.append(faers_nonsignal_source_to_evidence(path, mapping_path=args.mapping))
    for path in args.target_adr:
        if path.exists():
            parts.append(
                target_adr_source_to_evidence(
                    path,
                    source_name="ADReCS-Target",
                    source_priority=2,
                    source_role="curated_target_adr_positive",
                    confidence=0.9,
                )
            )
    _extend(
        parts,
        args.lit_pcba,
        source_name="LIT-PCBA",
        source_priority=5,
        source_role="benchmark_only",
        mapping_path=args.mapping,
        confidence=0.8,
        benchmark_only=True,
    )
    _extend(
        parts,
        args.dude,
        source_name="DUD-E",
        source_priority=5,
        source_role="benchmark_only_decoy",
        mapping_path=args.mapping,
        confidence=0.5,
        benchmark_only=True,
    )
    _extend(
        parts,
        args.muv,
        source_name="MUV",
        source_priority=5,
        source_role="benchmark_only_virtual_screening",
        mapping_path=args.mapping,
        confidence=0.8,
        benchmark_only=True,
    )
    for path in args.negative_evidence:
        if path.exists():
            parts.append(negative_evidence_to_raw(path))

    nonempty_parts = [part for part in parts if not part.empty]
    raw = pd.concat(nonempty_parts, ignore_index=True) if nonempty_parts else pd.DataFrame()
    collapsed = collapse_four_state_evidence(raw)
    joined = None
    if args.feature_table is not None and args.feature_table.exists():
        features = pd.read_csv(args.feature_table, low_memory=False)
        joined = join_collapsed_labels(features, collapsed)
    write_four_state_outputs(raw, collapsed, args.out_dir, joined=joined)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
