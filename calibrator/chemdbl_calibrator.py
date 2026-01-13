import argparse
import json
import os
import sys
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from unittest import mock

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))
DEEPCOY_MODULE_ROOT = REPO_ROOT / "DeepCoy_duds"
if str(DEEPCOY_MODULE_ROOT) not in sys.path:
    sys.path.insert(0, str(DEEPCOY_MODULE_ROOT))

from DeepCoy_duds import generate_dud_library  # noqa: E402
from DeepCoy_duds.external_sources import fetch_chembl_labeled_smiles  # noqa: E402
from DeepCoy_duds.generate_dud_library import get_uniprot_and_ec  # noqa: E402

DEFAULT_PDBS = ["1T46", "6LU7", "1Q4X", "2AZR", "1UYG"]
DEFAULT_ACTIVITY_TYPES = ["Ki", "Kd", "IC50", "EC50"]
DEFAULT_LABEL_THRESHOLDS = {
    "pchembl_strong": 7.0,
    "pchembl_weak": 5.0,
    "standard_value_nm_strong": 100.0,
    "standard_value_nm_weak": 10000.0,
}
DEFAULT_TIMEOUT = 20
DEFAULT_RETRIES = 2
DEFAULT_CHEMBL_MAX_PHASE = 4
DEFAULT_OUT_ROOT = REPO_ROOT / "extracted_ligands"
DEFAULT_CACHE_DIR = REPO_ROOT / "calibrator" / ".cache"
DEFAULT_DEEPCOY_ROOT = REPO_ROOT / "extracted_ligands" / "deepcoy"
DEFAULT_LOG_ROOT = REPO_ROOT / "calibrator" / "calibrator_logs"
LABEL_POLICY_DESC = (
    "label_policy=pchembl>=7 strong; pchembl>=5 weak; else non; "
    "fallback requires units=nm relation='=' with <=100 strong, <=10000 weak"
)


def parse_pdbs(raw: str) -> List[str]:
    tokens = []
    for part in raw.replace(";", ",").split(","):
        token = part.strip()
        if token:
            tokens.append(token.upper())
    return tokens


def _default_run_tag() -> str:
    atlas_run_id = os.environ.get("ATLAS_RUN_ID")
    if atlas_run_id:
        return atlas_run_id
    return datetime.utcnow().strftime("%Y%m%d_%H%M%S") + f"_{os.getpid()}"


@contextmanager
def tee_to_log(log_file: Path):
    log_file.parent.mkdir(parents=True, exist_ok=True)
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    log_handle = log_file.open("a")

    class _Tee:
        def __init__(self, stream):
            self.stream = stream

        def write(self, data):
            self.stream.write(data)
            log_handle.write(data)

        def flush(self):
            try:
                self.stream.flush()
            except Exception:
                pass
            log_handle.flush()

    sys.stdout = _Tee(original_stdout)
    sys.stderr = _Tee(original_stderr)
    try:
        yield
    finally:
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        log_handle.flush()
        log_handle.close()


def infer_uniprot_from_deepcoy_dir(
    pdb_id: str, deepcoy_root: Path = DEFAULT_DEEPCOY_ROOT
) -> Optional[str]:
    prefix = f"{pdb_id.upper()}_"
    if not deepcoy_root.is_dir():
        return None
    for child in deepcoy_root.iterdir():
        if not child.is_dir():
            continue
        name = child.name
        if name.upper().startswith(prefix) and "_" in name:
            suffix = name.split("_", 1)[1]
            if suffix:
                return suffix
    return None


def resolve_target(
    pdb_id: str, deepcoy_root: Path = DEFAULT_DEEPCOY_ROOT
) -> Tuple[Optional[str], List[str], Dict]:
    mapping: Dict[str, object] = {"source": "unresolved"}
    uniprot_id: Optional[str] = None
    ec_numbers: List[str] = []
    fallback_used = False

    try:
        original_get = generate_dud_library.requests.get

        def _wrapped_get(url, *args, **kwargs):
            nonlocal fallback_used
            if "uniprot.org" in url:
                fallback_used = True
            return original_get(url, *args, **kwargs)

        with mock.patch(
            "DeepCoy_duds.generate_dud_library.requests.get", side_effect=_wrapped_get
        ):
            uniprot_id, ec_numbers = get_uniprot_and_ec(pdb_id)
        if uniprot_id:
            mapping["source"] = (
                "uniprot_fallback" if fallback_used else "rcsb_polymer_entity"
            )
            mapping["details"] = {"resolver": "get_uniprot_and_ec"}
            mapping["ec_count"] = len(ec_numbers or [])
            return uniprot_id, ec_numbers or [], mapping
        else:
            mapping["details"] = {"reason": "get_uniprot_and_ec returned no match"}
    except Exception as exc:
        mapping["error"] = str(exc)

    inferred_uniprot = infer_uniprot_from_deepcoy_dir(pdb_id, deepcoy_root)
    if inferred_uniprot:
        mapping["source"] = "deepcoy_dir_fallback"
        mapping["details"] = {
            "deepcoy_dir": str(deepcoy_root),
            "inferred_from": inferred_uniprot,
        }
        mapping["ec_count"] = 0
        return inferred_uniprot, [], mapping

    mapping["source"] = "unresolved"
    if mapping.get("error"):
        mapping["details"] = {"reason": mapping["error"]}
    mapping["ec_count"] = 0
    return None, [], mapping


def _write_smiles(path: Path, smiles: List[str]) -> None:
    path.write_text("\n".join(smiles) + ("\n" if smiles else ""))


def _log_source_audit(
    uniprot_id: str, activity_types: List[str], chembl_max_phase: Optional[int]
) -> None:
    types_param = ",".join(activity_types)
    print(
        "[calibrator.source_audit] source=chembl endpoint=/target "
        f"params=target_components__accession={uniprot_id} format=json"
    )
    print(
        "[calibrator.source_audit] source=chembl endpoint=/assay "
        "params=target_chembl_id=<TARGET_ID> assay_type=B relationship_type=D format=json"
    )
    print(
        "[calibrator.source_audit] source=chembl endpoint=/activity "
        f"params=assay_chembl_id=<ASSAY_ID> standard_type__in={types_param} format=json"
    )
    print(f"[calibrator.source_audit] {LABEL_POLICY_DESC}")
    print(f"[calibrator.source_audit] phase_filter chembl_max_phase={chembl_max_phase}")


def _log_chembl_debug(chembl_meta: Dict) -> None:
    debug = chembl_meta.get("debug") or {}
    skips = debug.get("skips") or {}
    if skips:
        for key, val in sorted(skips.items(), key=lambda kv: kv[1], reverse=True):
            print(f"[calibrator.chembl_skip] {key}={val}")
    targets_sample = debug.get("target_ids_sample") or []
    assays_sample = debug.get("assay_ids_sample") or []
    if targets_sample:
        print(
            "[calibrator.source_audit] source=chembl endpoint=/assay "
            f"params=target_chembl_id={','.join(targets_sample)} assay_type=B relationship_type=D format=json"
        )
    if assays_sample:
        print(
            "[calibrator.source_audit] source=chembl endpoint=/activity "
            f"params=assay_chembl_id={','.join(assays_sample)} format=json"
        )
    types_seen = debug.get("activity_types_seen_top") or []
    if types_seen:
        formatted = "; ".join(
            f"{entry.get('type')}={entry.get('count')}" for entry in types_seen
        )
        print(f"[calibrator.chembl_activity_types] {formatted}")


def _format_histogram(hist: Dict[str, int], top_n: int = 3) -> str:
    items = sorted(hist.items(), key=lambda kv: kv[1], reverse=True)
    top_items = items[:top_n]
    return (
        "; ".join(f"{name}={count}" for name, count in top_items) if top_items else ""
    )


def _log_chembl_telemetry(chembl_meta: Dict, max_urls_log: int = 30) -> Dict:
    telemetry = chembl_meta.get("telemetry") or {}
    activity = telemetry.get("activity_sanity") or {}
    labeling = telemetry.get("labeling_sanity") or {}
    molecule = telemetry.get("molecule_sanity") or {}
    relation_top = _format_histogram(activity.get("relation_histogram") or {})
    units_top = _format_histogram(activity.get("units_histogram") or {})
    print(
        "[calibrator.telemetry] activity_sanity "
        f"n_total={activity.get('n_activities_total', 0)} "
        f"with_molecule_id={activity.get('n_with_molecule_chembl_id', 0)} "
        f"with_pchembl={activity.get('n_with_pchembl_value', 0)} "
        f"standard_value_parseable={activity.get('n_standard_value_parseable', 0)} "
        f"relation_top={relation_top} units_top={units_top}"
    )
    print(
        "[calibrator.telemetry] labeling_sanity "
        f"labeled_strong={labeling.get('n_labeled_strong', 0)} "
        f"labeled_weak={labeling.get('n_labeled_weak', 0)} "
        f"labeled_non={labeling.get('n_labeled_non', 0)} "
        f"rejected_relation={labeling.get('n_rejected_by_relation', 0)} "
        f"rejected_units={labeling.get('n_rejected_by_units', 0)} "
        f"rejected_missing_pchembl={labeling.get('n_rejected_by_missing_pchembl', 0)} "
        f"rejected_missing_standard_value={labeling.get('n_rejected_by_missing_standard_value', 0)} "
        f"rejected_value_threshold={labeling.get('n_rejected_by_value_threshold', 0)} "
        f"rejected_type_not_allowed={labeling.get('n_rejected_by_type_not_allowed', 0)} "
        f"rejected_phase_filtered={labeling.get('n_rejected_by_phase_filtered', 0)}"
    )
    print(
        "[calibrator.telemetry] molecule_sanity "
        f"fetch_attempted={molecule.get('n_molecule_fetch_attempted', 0)} "
        f"fetch_failed_http={molecule.get('n_molecule_fetch_failed_http', 0)} "
        f"missing_smiles={molecule.get('n_molecule_missing_smiles', 0)} "
        f"parsed_smiles_ok={molecule.get('n_molecule_parsed_smiles_ok', 0)}"
    )
    request_urls = telemetry.get("request_urls") or []
    truncated = telemetry.get("request_urls_truncated", False)
    for url in request_urls[:max_urls_log]:
        print(f"[calibrator.source] url={url}")
    if request_urls and (len(request_urls) > max_urls_log or truncated):
        print(
            "[calibrator.source] url_sample_truncated=True "
            f"logged={len(request_urls[:max_urls_log])} total_recorded={len(request_urls)}"
        )
    elif truncated:
        print("[calibrator.source] url_sample_truncated=True")
    return telemetry


def _log_rejected_value_samples(telemetry: Dict) -> None:
    labeling = (telemetry or {}).get("labeling_sanity") or {}
    samples = labeling.get("rejected_value_threshold_samples") or []
    if not samples:
        return
    max_samples = labeling.get("rejected_value_threshold_samples_max")
    try:
        max_samples_val = int(max_samples) if max_samples is not None else "n/a"
    except Exception:
        max_samples_val = max_samples
    print(
        "[calibrator.reject_sample] kind=rejected_value_threshold "
        f"count={len(samples)} max={max_samples_val}"
    )
    ordered_keys = [
        "pchembl_value",
        "standard_value",
        "standard_units",
        "standard_relation",
        "standard_type",
        "activity_chembl_id",
        "assay_chembl_id",
        "molecule_chembl_id",
    ]
    for sample in samples:
        record = sample if isinstance(sample, dict) else {}
        parts = [f"{key}={record.get(key)}" for key in ordered_keys]
        print(f"[calibrator.reject_sample] {' '.join(parts)}")


def _diagnose_from_telemetry(telemetry: Dict) -> str:
    activity = telemetry.get("activity_sanity") or {}
    labeling = telemetry.get("labeling_sanity") or {}
    molecule = telemetry.get("molecule_sanity") or {}
    if not (activity or labeling or molecule):
        return "[calibrator.diagnose] telemetry_inconclusive"
    if activity.get("n_with_molecule_chembl_id", 0) == 0:
        return (
            "[calibrator.diagnose] likely_schema_keying_issue "
            "activity_items_missing_molecule_chembl_id"
        )
    if (
        labeling.get("n_labeled_strong", 0) == 0
        and labeling.get("n_labeled_weak", 0) == 0
        and labeling.get("n_labeled_non", 0) == 0
    ):
        return "[calibrator.diagnose] likely_policy_overfiltering no_labels_assigned"
    if (
        molecule.get("n_molecule_fetch_attempted", 0) > 0
        and molecule.get("n_molecule_parsed_smiles_ok", 0) == 0
    ):
        return (
            "[calibrator.diagnose] likely_molecule_parsing_or_endpoint_issue "
            "no_smiles_parsed"
        )
    return "[calibrator.diagnose] telemetry_inconclusive"


def run_calibrator_for_pdb(
    pdb_id: str,
    *,
    out_root: Path = DEFAULT_OUT_ROOT,
    cache_dir: Path = DEFAULT_CACHE_DIR,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = DEFAULT_RETRIES,
    activity_types: Optional[List[str]] = None,
    chembl_max_phase: Optional[int] = DEFAULT_CHEMBL_MAX_PHASE,
    label_thresholds: Dict = DEFAULT_LABEL_THRESHOLDS,
    deepcoy_root: Path = DEFAULT_DEEPCOY_ROOT,
    log_dir: Optional[Path] = None,
    run_tag: Optional[str] = None,
    debug_chembl: bool = False,
    debug_max_ids: int = 25,
    debug_reject_samples: int = 25,
    fetch_fn=None,
) -> Dict:
    pdb_norm = pdb_id.upper()
    out_dir = Path(out_root) / f"{pdb_norm}_calibrator"
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    resolved_activity_types = activity_types or list(DEFAULT_ACTIVITY_TYPES)
    run_tag_value = run_tag or _default_run_tag()
    log_root = Path(log_dir) if log_dir else DEFAULT_LOG_ROOT
    run_log_dir = log_root / run_tag_value
    log_file = run_log_dir / f"{pdb_norm}.log"

    with tee_to_log(log_file):
        print(f"[calibrator.logs] run_tag={run_tag_value} log_file={log_file}")
        uniprot_id, ec_numbers, mapping = resolve_target(pdb_norm, deepcoy_root)
        meta: Dict = {
            "pdb_id": pdb_norm,
            "uniprot_id": uniprot_id,
            "ec_numbers": ec_numbers or [],
            "mapping": mapping,
            "label_thresholds": label_thresholds,
            "activity_types": resolved_activity_types,
            "chembl_max_phase": chembl_max_phase,
            "cache": {},
            "output_dir": str(out_dir),
            "log_file": str(log_file),
            "run_tag": run_tag_value,
        }

        if not uniprot_id:
            meta["status"] = "unresolved"
            meta["counts"] = {}
            meta_path = out_dir / "calibrator_meta.json"
            meta_path.write_text(json.dumps(meta, indent=2))
            audit_path = out_dir / "calibrator_audit.json"
            audit_path.write_text(
                json.dumps(
                    {
                        "run_tag": run_tag_value,
                        "log_file": str(log_file),
                        "pdb_id": pdb_norm,
                        "mapping": mapping,
                        "status": "unresolved",
                    },
                    indent=2,
                )
            )
            return meta

        _log_source_audit(uniprot_id, resolved_activity_types, chembl_max_phase)
        fetch_impl = fetch_fn or fetch_chembl_labeled_smiles
        labels, chembl_meta = fetch_impl(
            uniprot_id,
            pdb_norm,
            cache_dir,
            timeout,
            retries,
            resolved_activity_types,
            chembl_max_phase,
            label_thresholds,
            debug=debug_chembl,
            debug_max_ids=debug_max_ids,
            debug_rejection_samples_max=debug_reject_samples,
        )
        bin_counts = {k: len(v) for k, v in labels.items()}
        chembl_counts = chembl_meta.get("counts", {})
        counts = {**chembl_counts, **bin_counts}
        for name, smiles in labels.items():
            _write_smiles(out_dir / f"{name}_binders.smi", smiles)

        cache_status = chembl_meta.get("cached", False)
        print(
            f"[calibrator.chembl] cached={cache_status} "
            f"targets={chembl_counts.get('targets', 0)} assays={chembl_counts.get('assays', 0)} "
            f"activities={chembl_counts.get('activities', 0)} molecules={chembl_counts.get('molecules', 0)} "
            f"strong={bin_counts.get('strong', 0)} weak={bin_counts.get('weak', 0)} non={bin_counts.get('non', 0)}"
        )
        if chembl_counts.get("targets", 0) == 0:
            print(
                "[calibrator.chembl] ChEMBL target lookup returned zero targets for this UniProt; "
                "likely UniProt mapping too broad or ChEMBL accession mismatch."
            )
        _log_chembl_debug(chembl_meta)
        telemetry = _log_chembl_telemetry(chembl_meta)
        _log_rejected_value_samples(telemetry)
        print(_diagnose_from_telemetry(telemetry))

        meta.update(
            {
                "status": "ok" if not chembl_meta.get("error") else "error",
                "counts": counts,
                "cache": {"chembl": cache_status},
                "provenance": {"fetcher": getattr(fetch_impl, "__name__", "unknown")},
                "timestamp": datetime.utcnow().isoformat() + "Z",
                "telemetry": telemetry,
                "chembl_request_urls_sample": (telemetry.get("request_urls") or []),
                "chembl_request_urls_truncated": telemetry.get(
                    "request_urls_truncated", False
                ),
            }
        )
        if chembl_meta.get("error"):
            meta["error"] = chembl_meta["error"]

        audit_path = out_dir / "calibrator_audit.json"
        audit_payload = {
            "run_tag": run_tag_value,
            "log_file": str(log_file),
            "pdb_id": pdb_norm,
            "uniprot_id": uniprot_id,
            "ec_numbers": ec_numbers or [],
            "mapping": mapping,
            "chembl_meta": chembl_meta,
            "telemetry": telemetry,
        }
        audit_path.write_text(json.dumps(audit_payload, indent=2))

        meta_path = out_dir / "calibrator_meta.json"
        meta_path.write_text(json.dumps(meta, indent=2))
        return meta


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export labeled ChEMBL SMILES per PDB via UniProt mapping."
    )
    parser.add_argument(
        "--pdbs",
        default=",".join(DEFAULT_PDBS),
        help="Comma-separated PDB IDs (whitespace allowed).",
    )
    parser.add_argument(
        "--out-root",
        default=str(DEFAULT_OUT_ROOT),
        help="Root directory for calibrator outputs.",
    )
    parser.add_argument(
        "--cache-dir",
        default=str(DEFAULT_CACHE_DIR),
        help="Cache directory for ChEMBL responses.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=DEFAULT_TIMEOUT,
        help="HTTP timeout (seconds) for ChEMBL queries.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=DEFAULT_RETRIES,
        help="Retries for ChEMBL queries.",
    )
    parser.add_argument(
        "--activity-types",
        default=",".join(DEFAULT_ACTIVITY_TYPES),
        help="Comma-separated allowed standard types (Ki,Kd,IC50,EC50).",
    )
    parser.add_argument(
        "--chembl-max-phase",
        type=int,
        default=DEFAULT_CHEMBL_MAX_PHASE,
        help="Maximum ChEMBL clinical phase to include (None for all).",
    )
    parser.add_argument(
        "--run-tag",
        default=None,
        help="Run tag for logging (default: ATLAS_RUN_ID or timestamp+pid).",
    )
    parser.add_argument(
        "--log-dir",
        default=str(DEFAULT_LOG_ROOT),
        help="Directory to store calibrator logs (default: calibrator/calibrator_logs).",
    )
    parser.add_argument(
        "--debug-chembl",
        action="store_true",
        help="Enable verbose ChEMBL debug metadata (ids samples, skip counters).",
    )
    parser.add_argument(
        "--debug-reject-samples",
        type=int,
        default=25,
        help="Max number of value-threshold rejection samples to log (default: 25).",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    pdbs = parse_pdbs(args.pdbs or "")
    if not pdbs:
        pdbs = list(DEFAULT_PDBS)
    out_root = Path(args.out_root)
    cache_dir = Path(args.cache_dir)
    log_dir = Path(args.log_dir) if args.log_dir else DEFAULT_LOG_ROOT
    run_tag = args.run_tag or _default_run_tag()
    activity_types = [
        t.strip() for t in args.activity_types.split(",") if t.strip()
    ] or list(DEFAULT_ACTIVITY_TYPES)
    for pdb_id in pdbs:
        print(f"[calibrator] Processing {pdb_id.upper()}...")
        meta = run_calibrator_for_pdb(
            pdb_id,
            out_root=out_root,
            cache_dir=cache_dir,
            timeout=args.timeout,
            retries=args.retries,
            activity_types=activity_types,
            chembl_max_phase=args.chembl_max_phase,
            label_thresholds=DEFAULT_LABEL_THRESHOLDS,
            deepcoy_root=DEFAULT_DEEPCOY_ROOT,
            log_dir=log_dir,
            run_tag=run_tag,
            debug_chembl=args.debug_chembl,
            debug_reject_samples=args.debug_reject_samples,
        )
        status = meta.get("status", "unknown")
        mapping_source = (meta.get("mapping") or {}).get("source", "unknown")
        print(
            f"[calibrator] {pdb_id.upper()} status={status} mapping={mapping_source} "
            f"strong={meta.get('counts', {}).get('strong', 0)} "
            f"weak={meta.get('counts', {}).get('weak', 0)} "
            f"non={meta.get('counts', {}).get('non', 0)}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
