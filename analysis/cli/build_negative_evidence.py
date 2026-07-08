from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from analysis.external.negative_evidence import (
    normalize_faers_nonsignal_source,
    normalize_measured_negative_source,
    normalize_omop_ohdsi_negative_source,
    write_negative_evidence,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Normalize measured/reliable negative evidence into the Atlas four-state label schema."
    )
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--mapping", type=Path, default=None)
    parser.add_argument("--measured", nargs="*", type=Path, default=[])
    parser.add_argument("--measured-format", default="generic_bioactivity")
    parser.add_argument("--inactive-threshold-nm", type=float, default=10000.0)
    parser.add_argument("--omop-ohdsi", nargs="*", type=Path, default=[])
    parser.add_argument("--faers-nonsignal", nargs="*", type=Path, default=[])
    parser.add_argument("--positive-drug-adr-evidence", type=Path, default=None)
    parser.add_argument("--min-drug-reports", type=int, default=100)
    parser.add_argument("--min-event-reports", type=int, default=100)
    parser.add_argument("--min-pair-reports", type=int, default=0)
    args = parser.parse_args(argv)
    frames: list[pd.DataFrame] = []
    for path in args.measured:
        frames.append(
            normalize_measured_negative_source(
                path,
                args.mapping,
                source_format=args.measured_format,
                inactive_threshold_nM=args.inactive_threshold_nm,
            )
        )
    for path in args.omop_ohdsi:
        frames.append(normalize_omop_ohdsi_negative_source(path, args.mapping))
    for path in args.faers_nonsignal:
        frames.append(
            normalize_faers_nonsignal_source(
                path,
                args.mapping,
                positive_evidence_path=args.positive_drug_adr_evidence,
                min_drug_reports=args.min_drug_reports,
                min_event_reports=args.min_event_reports,
                min_pair_reports=args.min_pair_reports,
            )
        )
    combined = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()
    write_negative_evidence(
        combined,
        args.out,
        {
            "measured_sources": [str(path) for path in args.measured],
            "omop_ohdsi_sources": [str(path) for path in args.omop_ohdsi],
            "faers_nonsignal_sources": [str(path) for path in args.faers_nonsignal],
            "policy": "0 labels require measured inactive, curated drug-ADR negative control, FAERS observed non-signal, or fold-local reliable negative evidence.",
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

