import argparse
import hashlib
import inspect
import json
import os
import sys
import shutil
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from unittest import mock

import requests  # type: ignore[import-untyped]

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))
DEEPCOY_MODULE_ROOT = REPO_ROOT / "DeepCoy_duds"
if str(DEEPCOY_MODULE_ROOT) not in sys.path:
    sys.path.insert(0, str(DEEPCOY_MODULE_ROOT))

from calibrator import uniprot_resolver  # noqa: E402
from DeepCoy_duds import generate_dud_library  # noqa: E402
from DeepCoy_duds.external_sources import (  # noqa: E402
    CHEMBL_API_BASE,
    fetch_chembl_labeled_smiles,
    load_cached_smiles,
    save_cached_smiles,
)
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
    chain_segments: Dict[str, Dict] = {}
    selected_chain: Optional[str] = None
    selected_segment: Optional[Dict] = None
    resolver_info: Dict[str, object] = {}

    try:
        result = uniprot_resolver.resolve_chain_uniprot_segments(
            pdb_id, write_file=True, return_info=True
        )
        if isinstance(result, tuple):
            chain_segments, resolver_info = result  # type: ignore[misc]
        else:  # pragma: no cover - defensive
            chain_segments = result
    except FileNotFoundError as exc:
        resolver_info = {"error": str(exc)}
        chain_segments = {}

    selected_chain, selected_segment = uniprot_resolver.select_primary_chain(
        chain_segments
    )
    if selected_segment:
        uniprot_id = selected_segment.get("uniprot")

    mapping.update(
        {
            "source": (selected_segment or {}).get("source", "unresolved"),
            "chain_segments": chain_segments,
            "selected_chain": selected_chain,
            "selected_segment": selected_segment,
            "dbref_count": resolver_info.get("dbref_count"),
            "fallback_used": resolver_info.get("fallback_used"),
        }
    )

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
            lookup_uniprot, ec_numbers = get_uniprot_and_ec(pdb_id)
        if not uniprot_id and lookup_uniprot:
            uniprot_id = lookup_uniprot
            mapping["source"] = (
                "uniprot_fallback" if fallback_used else "rcsb_polymer_entity"
            )
        mapping["details"] = {"resolver": "get_uniprot_and_ec"}
    except Exception as exc:  # pragma: no cover - network fallback safety
        mapping["error"] = str(exc)

    if not uniprot_id:
        inferred_uniprot = infer_uniprot_from_deepcoy_dir(pdb_id, deepcoy_root)
        if inferred_uniprot:
            uniprot_id = inferred_uniprot
            mapping["source"] = "deepcoy_dir_fallback"
            mapping["details"] = {
                "deepcoy_dir": str(deepcoy_root),
                "inferred_from": inferred_uniprot,
            }

    if (
        uniprot_id
        and selected_chain
        and selected_segment
        and not selected_segment.get("uniprot")
    ):
        updated_segment = dict(selected_segment)
        updated_segment["uniprot"] = uniprot_id
        selected_segment = updated_segment
        chain_segments[selected_chain] = updated_segment
        mapping["selected_segment"] = updated_segment
        mapping["chain_segments"] = chain_segments

    if not uniprot_id:
        mapping["details"] = {
            "reason": mapping.get("error") or "no UniProt mapping resolved"
        }
    mapping["ec_count"] = len(ec_numbers or [])
    return uniprot_id, ec_numbers or [], mapping


def _write_smiles(path: Path, smiles: List[str]) -> None:
    path.write_text("\n".join(smiles) + ("\n" if smiles else ""))


def _smiles_hash(smiles: str) -> str:
    return hashlib.sha1(smiles.encode("utf-8")).hexdigest()[:10]


def _ligand_id(label: str, smiles: str) -> str:
    return f"chembl:{label}_{_smiles_hash(smiles)}"


def _write_jsonl(path: Path, rows: List[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=True) + "\n")


def _normalize_supporting_activity(raw: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    payload = {
        "pchembl_value": None,
        "standard_value": None,
        "standard_units": None,
        "standard_type": None,
        "standard_relation": None,
    }
    if isinstance(raw, dict):
        for key in payload:
            payload[key] = raw.get(key)
    return payload


def _label_by_smiles(labels: Dict[str, List[str]]) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for label, smi_list in labels.items():
        for smi in smi_list or []:
            mapping[str(smi)] = str(label)
    return mapping


def _minimal_label_records(labels: Dict[str, List[str]]) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    for label, smi_list in labels.items():
        for smi in smi_list or []:
            records.append(
                {
                    "molecule_chembl_id": None,
                    "canonical_smiles": smi,
                    "inchi_key": None,
                    "label": label,
                    "supporting_activity": _normalize_supporting_activity(None),
                    "assay_chembl_id": None,
                    "activity_chembl_id": None,
                    "document_chembl_id": None,
                }
            )
    return records


def _build_calibrator_ligand_records(
    pdb_id: str,
    selected_chain: Optional[str],
    selected_uniprot: Optional[str],
    labels: Dict[str, List[str]],
    chembl_meta: Optional[Dict[str, Any]],
    activity_types: List[str],
    chembl_max_phase: Optional[int],
) -> List[Dict[str, Any]]:
    label_lookup = _label_by_smiles(labels)
    raw_records: List[Dict[str, Any]] = []
    if isinstance(chembl_meta, dict):
        meta_records = chembl_meta.get("records")
        if isinstance(meta_records, list):
            raw_records = [r for r in meta_records if isinstance(r, dict)]
    if not raw_records:
        raw_records = _minimal_label_records(labels)

    records_out: List[Dict[str, Any]] = []
    for raw in raw_records:
        smi = raw.get("canonical_smiles") or raw.get("smiles")
        if not smi:
            continue
        label = label_lookup.get(str(smi)) or raw.get("label")
        if label not in {"strong", "weak", "non"}:
            continue
        records_out.append(
            {
                "pdb_id": pdb_id,
                "selected_chain": selected_chain,
                "selected_uniprot": selected_uniprot,
                "ligand_id": _ligand_id(label, str(smi)),
                "molecule_chembl_id": raw.get("molecule_chembl_id"),
                "canonical_smiles": smi,
                "inchi_key": raw.get("inchi_key"),
                "label": label,
                "supporting_activity": _normalize_supporting_activity(
                    raw.get("supporting_activity")
                    if isinstance(raw.get("supporting_activity"), dict)
                    else None
                ),
                "assay_chembl_id": raw.get("assay_chembl_id"),
                "activity_chembl_id": raw.get("activity_chembl_id"),
                "document_chembl_id": raw.get("document_chembl_id"),
                "chembl_max_phase": chembl_max_phase,
                "activity_types": list(activity_types),
            }
        )
    return records_out


def _supports_kwarg(func, name: str) -> bool:
    try:
        return name in inspect.signature(func).parameters
    except (TypeError, ValueError):
        return False


def _load_chain_uniprot_map(
    pdb_id: str, fallback: Optional[Dict[str, Dict[str, Any]]] = None
) -> Dict[str, Dict[str, Any]]:
    out_file = REPO_ROOT / "calibrator" / pdb_id / "chain_uniprot.json"
    if out_file.is_file():
        try:
            return json.loads(out_file.read_text())
        except Exception:
            return fallback or {}
    return fallback or {}


def _build_uniprot_chain_index(
    chain_map: Dict[str, Dict[str, Any]],
) -> Dict[str, List[str]]:
    index: Dict[str, List[str]] = {}
    for chain_id, entry in (chain_map or {}).items():
        if not isinstance(entry, dict):
            continue
        uniprot = entry.get("uniprot")
        if uniprot:
            index.setdefault(str(uniprot).upper(), []).append(str(chain_id))
    return index


def _first_present(*values: Any) -> Optional[str]:
    for value in values:
        if value:
            return str(value)
    return None


def _summarize_site_components(
    binding_site: Optional[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    summary: List[Dict[str, Any]] = []
    site_components = (
        binding_site.get("site_components") if isinstance(binding_site, dict) else None
    )
    for comp in site_components or []:
        if not isinstance(comp, dict):
            continue
        component = (
            comp.get("component") if isinstance(comp.get("component"), dict) else {}
        )
        accession = _first_present(
            comp.get("accession"),
            comp.get("component_accession"),
            comp.get("protein_accession"),
            comp.get("uniprot_accession"),
            component.get("accession"),
            component.get("component_accession"),
        )
        summary.append(
            {
                "component_id": comp.get("component_id")
                or comp.get("site_component_id")
                or component.get("component_id"),
                "component_type": comp.get("component_type")
                or comp.get("type")
                or component.get("component_type"),
                "accession": accession,
                "relationship": comp.get("relationship"),
            }
        )
    return summary


def _extract_uniprot_accessions(binding_site: Optional[Dict[str, Any]]) -> List[str]:
    accessions: List[str] = []
    seen: set[str] = set()
    site_components = (
        binding_site.get("site_components") if isinstance(binding_site, dict) else None
    )
    for comp in site_components or []:
        if not isinstance(comp, dict):
            continue
        component = (
            comp.get("component") if isinstance(comp.get("component"), dict) else {}
        )
        accession = _first_present(
            comp.get("accession"),
            comp.get("component_accession"),
            comp.get("protein_accession"),
            comp.get("uniprot_accession"),
            component.get("accession"),
            component.get("component_accession"),
        )
        if accession:
            accession_norm = str(accession).upper()
            if accession_norm not in seen:
                seen.add(accession_norm)
                accessions.append(accession_norm)
    return accessions


def _map_evidence_chains(
    evidence_uniprots: List[str], chain_index: Dict[str, List[str]]
) -> List[str]:
    chains: List[str] = []
    for accession in evidence_uniprots or []:
        for chain_id in chain_index.get(str(accession).upper(), []):
            if chain_id not in chains:
                chains.append(chain_id)
    return chains


def _chembl_get_json(
    url: str,
    params: Dict[str, Any],
    timeout: int,
    retries: int,
    session=None,
) -> Dict[str, Any]:
    last_err = None
    client = session or requests
    for _ in range(retries + 1):
        try:
            resp = client.get(
                url,
                params=params,
                headers={"Accept": "application/json"},
                timeout=timeout,
            )
            if resp.status_code == 200:
                return resp.json()
            last_err = f"status={resp.status_code}"
        except Exception as exc:
            last_err = str(exc)
    raise RuntimeError(last_err or "request failed")


def _load_cached_response(
    cache_dir: Path, source: str, key: Dict[str, Any]
) -> Optional[Any]:
    cached = load_cached_smiles(cache_dir, source, key)
    if cached and isinstance(cached, dict):
        response = cached.get("response")
        if isinstance(response, (dict, list)):
            return response
    return None


def _extract_mechanisms(payload: Any) -> List[Dict[str, Any]]:
    if isinstance(payload, dict):
        mechanisms = (
            payload.get("mechanisms")
            or payload.get("mechanism")
            or payload.get("data")
            or []
        )
    elif isinstance(payload, list):
        mechanisms = payload
    else:
        mechanisms = []
    return [entry for entry in mechanisms if isinstance(entry, dict)]


def _fetch_chembl_mechanisms(
    molecule_chembl_id: str,
    target_chembl_id: str,
    cache_dir: Path,
    timeout: int,
    retries: int,
    *,
    session=None,
    force_refresh: bool = False,
    fetch_fn=None,
) -> List[Dict[str, Any]]:
    if fetch_fn:
        return fetch_fn(molecule_chembl_id, target_chembl_id) or []
    key = {
        "source": "chembl_mechanism",
        "molecule_chembl_id": molecule_chembl_id,
        "target_chembl_id": target_chembl_id,
    }
    if not force_refresh:
        cached = _load_cached_response(cache_dir, "chembl_mechanism", key)
        if cached is not None:
            return _extract_mechanisms(cached)
    response = _chembl_get_json(
        f"{CHEMBL_API_BASE}/mechanism",
        {
            "molecule_chembl_id": molecule_chembl_id,
            "target_chembl_id": target_chembl_id,
            "format": "json",
        },
        timeout,
        retries,
        session=session,
    )
    save_cached_smiles(cache_dir, "chembl_mechanism", key, {"response": response})
    return _extract_mechanisms(response)


def _fetch_binding_site(
    site_id: str,
    cache_dir: Path,
    timeout: int,
    retries: int,
    *,
    session=None,
    force_refresh: bool = False,
    fetch_fn=None,
) -> Optional[Dict[str, Any]]:
    if not site_id:
        return None
    if fetch_fn:
        return fetch_fn(site_id)
    key = {"source": "chembl_binding_site", "site_id": site_id}
    if not force_refresh:
        cached = _load_cached_response(cache_dir, "chembl_binding_site", key)
        if isinstance(cached, dict):
            return cached
    response = _chembl_get_json(
        f"{CHEMBL_API_BASE}/binding_site/{site_id}.json",
        {},
        timeout,
        retries,
        session=session,
    )
    save_cached_smiles(cache_dir, "chembl_binding_site", key, {"response": response})
    return response if isinstance(response, dict) else None


def _build_site_evidence_records(
    pdb_id: str,
    ligand_records: List[Dict[str, Any]],
    target_ids: List[str],
    chain_map: Dict[str, Dict[str, Any]],
    cache_dir: Path,
    timeout: int,
    retries: int,
    force_refresh: bool,
    *,
    session=None,
    mechanism_fetch_fn=None,
    binding_site_fetch_fn=None,
) -> List[Dict[str, Any]]:
    if not ligand_records or not target_ids:
        return []
    chain_index = _build_uniprot_chain_index(chain_map)
    ligand_id_by_mol: Dict[str, Optional[str]] = {}
    for record in ligand_records:
        mol_id = record.get("molecule_chembl_id")
        if mol_id and mol_id not in ligand_id_by_mol:
            ligand_id_by_mol[str(mol_id)] = record.get("ligand_id")

    evidence_records: List[Dict[str, Any]] = []
    seen_pairs: set[tuple[str, str]] = set()
    for record in ligand_records:
        mol_id = record.get("molecule_chembl_id")
        if not mol_id:
            continue
        mol_id_str = str(mol_id)
        for target_id in target_ids:
            target_id_str = str(target_id)
            pair = (mol_id_str, target_id_str)
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            mechanisms = _fetch_chembl_mechanisms(
                mol_id_str,
                target_id_str,
                cache_dir,
                timeout,
                retries,
                session=session,
                force_refresh=force_refresh,
                fetch_fn=mechanism_fetch_fn,
            )
            if not mechanisms:
                evidence_records.append(
                    {
                        "pdb_id": pdb_id,
                        "ligand_id": ligand_id_by_mol.get(mol_id_str),
                        "molecule_chembl_id": mol_id_str,
                        "target_chembl_id": target_id_str,
                        "site_id": None,
                        "binding_site_name": None,
                        "binding_site_comment": None,
                        "site_components": [],
                        "evidence_uniprot_accessions": [],
                        "evidence_chains": [],
                        "evidence_strength": "none",
                    }
                )
                continue
            for mechanism in mechanisms:
                site_id = _first_present(
                    mechanism.get("binding_site_id"),
                    mechanism.get("site_id"),
                    mechanism.get("binding_site"),
                )
                binding_site = (
                    _fetch_binding_site(
                        site_id,
                        cache_dir,
                        timeout,
                        retries,
                        session=session,
                        force_refresh=force_refresh,
                        fetch_fn=binding_site_fetch_fn,
                    )
                    if site_id
                    else None
                )
                site_name = _first_present(
                    binding_site.get("site_name")
                    if isinstance(binding_site, dict)
                    else None,
                    binding_site.get("binding_site_name")
                    if isinstance(binding_site, dict)
                    else None,
                    binding_site.get("name")
                    if isinstance(binding_site, dict)
                    else None,
                    mechanism.get("binding_site_name"),
                    mechanism.get("site_name"),
                )
                site_comment = _first_present(
                    binding_site.get("comment")
                    if isinstance(binding_site, dict)
                    else None,
                    binding_site.get("description")
                    if isinstance(binding_site, dict)
                    else None,
                    mechanism.get("comment"),
                    mechanism.get("mechanism_of_action"),
                )
                site_components = _summarize_site_components(binding_site)
                evidence_uniprots = _extract_uniprot_accessions(binding_site)
                evidence_chains = _map_evidence_chains(evidence_uniprots, chain_index)
                if site_id:
                    evidence_strength = (
                        "site_components_mapped"
                        if evidence_chains
                        else "site_id_present"
                    )
                else:
                    evidence_strength = "none"
                evidence_records.append(
                    {
                        "pdb_id": pdb_id,
                        "ligand_id": ligand_id_by_mol.get(mol_id_str),
                        "molecule_chembl_id": mol_id_str,
                        "target_chembl_id": target_id_str,
                        "site_id": site_id,
                        "binding_site_name": site_name,
                        "binding_site_comment": site_comment,
                        "site_components": site_components,
                        "evidence_uniprot_accessions": evidence_uniprots,
                        "evidence_chains": evidence_chains,
                        "evidence_strength": evidence_strength,
                    }
                )
    return evidence_records


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
    force_refresh: bool = False,
    fetch_fn=None,
    mechanism_fetch_fn=None,
    binding_site_fetch_fn=None,
    enable_site_evidence: bool = True,
) -> Dict:
    pdb_norm = pdb_id.upper()
    out_dir = Path(out_root) / f"{pdb_norm}_calibrator"
    if force_refresh and out_dir.exists():
        shutil.rmtree(out_dir, ignore_errors=True)
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
        chain_segments = (
            mapping.get("chain_segments") if isinstance(mapping, dict) else {}
        )
        selected_chain = (
            mapping.get("selected_chain") if isinstance(mapping, dict) else None
        )
        selected_segment = mapping.get("selected_segment") or {}
        selected_unp_start = (
            selected_segment.get("unp_start")
            if isinstance(selected_segment, dict)
            else None
        )
        selected_unp_end = (
            selected_segment.get("unp_end")
            if isinstance(selected_segment, dict)
            else None
        )
        selected_pdb_start = (
            selected_segment.get("pdb_start")
            if isinstance(selected_segment, dict)
            else None
        )
        selected_pdb_end = (
            selected_segment.get("pdb_end")
            if isinstance(selected_segment, dict)
            else None
        )
        selected_uniprot = uniprot_id or (
            selected_segment.get("uniprot")
            if isinstance(selected_segment, dict)
            else None
        )
        primary_uniprot = selected_uniprot
        meta: Dict = {
            "pdb_id": pdb_norm,
            "uniprot_id": primary_uniprot,
            "ec_numbers": ec_numbers or [],
            "mapping": mapping,
            "selected_chain": selected_chain,
            "selected_segment": selected_segment,
            "label_thresholds": label_thresholds,
            "activity_types": resolved_activity_types,
            "chembl_max_phase": chembl_max_phase,
            "cache": {},
            "output_dir": str(out_dir),
            "log_file": str(log_file),
            "run_tag": run_tag_value,
        }

        if not selected_uniprot:
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
                        "selected_chain": selected_chain,
                        "status": "unresolved",
                    },
                    indent=2,
                )
            )
            return meta

        _log_source_audit(selected_uniprot, resolved_activity_types, chembl_max_phase)
        fetch_impl = fetch_fn or fetch_chembl_labeled_smiles
        fetch_kwargs = {
            "unp_start": selected_unp_start,
            "unp_end": selected_unp_end,
            "debug": debug_chembl,
            "debug_max_ids": debug_max_ids,
            "debug_rejection_samples_max": debug_reject_samples,
        }
        if force_refresh:
            fetch_kwargs["force_refresh"] = True
        if _supports_kwarg(fetch_impl, "return_records"):
            fetch_kwargs["return_records"] = True
        labels, chembl_meta = fetch_impl(
            selected_uniprot,
            pdb_norm,
            cache_dir,
            timeout,
            retries,
            resolved_activity_types,
            chembl_max_phase,
            label_thresholds,
            **fetch_kwargs,
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

        selection_meta = chembl_meta.get("target_selection") or {}
        if not selection_meta:
            telemetry_selection = chembl_meta.get("telemetry", {}).get(
                "target_selection", {}
            )
            if isinstance(telemetry_selection, dict):
                selection_meta = telemetry_selection
        selected_targets = selection_meta.get("selected_ids") or []

        ligand_records = _build_calibrator_ligand_records(
            pdb_norm,
            selected_chain,
            primary_uniprot,
            labels,
            chembl_meta if isinstance(chembl_meta, dict) else None,
            resolved_activity_types,
            chembl_max_phase,
        )
        _write_jsonl(out_dir / "calibrator_ligands.jsonl", ligand_records)

        site_evidence_records: List[Dict[str, Any]] = []
        if enable_site_evidence and selected_targets:
            chain_map = _load_chain_uniprot_map(
                pdb_norm,
                chain_segments if isinstance(chain_segments, dict) else {},
            )
            chembl_session = requests.Session()
            try:
                site_evidence_records = _build_site_evidence_records(
                    pdb_norm,
                    ligand_records,
                    selected_targets,
                    chain_map,
                    cache_dir,
                    timeout,
                    retries,
                    force_refresh,
                    session=chembl_session,
                    mechanism_fetch_fn=mechanism_fetch_fn,
                    binding_site_fetch_fn=binding_site_fetch_fn,
                )
            finally:
                try:
                    chembl_session.close()
                except Exception:
                    pass
        _write_jsonl(out_dir / "calibrator_site_evidence.jsonl", site_evidence_records)

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
                "selected_chain": selected_chain,
                "selected_segment": selected_segment,
                "selected_uniprot": primary_uniprot,
                "selected_targets": selected_targets,
                "target_selection": selection_meta,
            }
        )
        if chembl_meta.get("error"):
            meta["error"] = chembl_meta["error"]

        audit_path = out_dir / "calibrator_audit.json"
        audit_payload = {
            "run_tag": run_tag_value,
            "log_file": str(log_file),
            "pdb_id": pdb_norm,
            "uniprot_id": primary_uniprot,
            "ec_numbers": ec_numbers or [],
            "mapping": mapping,
            "selected_chain": selected_chain,
            "selected_segment": selected_segment,
            "selected_targets": selected_targets,
            "target_selection": selection_meta,
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
