import hashlib
import json
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Set, Tuple, cast

import requests  # type: ignore[import-untyped]


class SourceAudit:
    """
    Lightweight collector for per-request audit details.
    Stores resolved URLs (with params) and optional paging metadata.
    """

    def __init__(self, enabled: bool = False, max_lines: int = 200):
        self.enabled = bool(enabled)
        self.max_lines = max_lines
        self.records: List[Dict[str, Any]] = []
        self._source_order: List[str] = []

    def record(
        self,
        *,
        source: str,
        purpose: str,
        url: str,
        params: Optional[Dict[str, Any]] = None,
        response_meta: Optional[Dict[str, Any]] = None,
    ) -> None:
        if not self.enabled:
            return
        entry: Dict[str, Any] = {
            "source": source,
            "purpose": purpose,
            "url": url,
            "timestamp": time.time(),
        }
        if params:
            entry["params"] = params
        if response_meta:
            entry["response_meta"] = response_meta
        self.records.append(entry)
        if source not in self._source_order:
            self._source_order.append(source)

    def as_json(self) -> List[Dict[str, Any]]:
        return list(self.records)

    def summary_lines(self) -> List[str]:
        if not self.enabled:
            return []
        by_source: Dict[str, List[str]] = defaultdict(list)
        for rec in self.records:
            src = rec.get("source") or "unknown"
            url = rec.get("url")
            if url:
                by_source[src].append(str(url))
        lines: List[str] = []
        for src in self._source_order:
            urls = by_source.get(src, [])
            count = len(urls)
            capped = urls[: self.max_lines]
            suffix = " (truncated)" if len(capped) < count else ""
            urls_part = "; ".join(capped)
            lines.append(
                f"[deepcoy.audit] source={src} requests={count} sample_urls={urls_part}{suffix}"
            )
        return lines


def _resolved_url(url: str, params: Optional[Dict[str, Any]] = None) -> str:
    try:
        req = requests.Request("GET", url, params=params).prepare()
        return req.url
    except Exception:
        if params:
            return f"{url}?{params}"
        return url


def _safe_get(d: Any, key: str, default=None):
    if not isinstance(d, dict):
        return default
    return d.get(key, default)


def _safe_dict(x: Any) -> Dict:
    return x if isinstance(x, dict) else {}


def _safe_list(x: Any) -> List:
    return x if isinstance(x, list) else []


def _increment(skip_counts: Optional[Dict[str, int]], key: str) -> None:
    if skip_counts is None:
        return
    skip_counts[key] = skip_counts.get(key, 0) + 1


def _record_audit(
    audit: Optional[SourceAudit],
    source: str,
    purpose: str,
    url: str,
    params: Optional[Dict[str, Any]] = None,
    response_meta: Optional[Dict[str, Any]] = None,
) -> None:
    if not audit or not audit.enabled:
        return
    audit.record(
        source=source,
        purpose=purpose,
        url=_resolved_url(url, params),
        params=params or {},
        response_meta=response_meta,
    )


def _make_cache_path(cache_dir: Path, source: str, key_dict: Dict) -> Path:
    cache_dir.mkdir(parents=True, exist_ok=True)
    key_str = json.dumps(key_dict, sort_keys=True)
    digest = hashlib.sha256(key_str.encode()).hexdigest()[:16]
    return cache_dir / f"{source}_{digest}.json"


def load_cached_smiles(cache_dir: Path, source: str, key_dict: Dict) -> Optional[Dict]:
    path = _make_cache_path(cache_dir, source, key_dict)
    if path.is_file():
        try:
            with path.open("r") as f:
                return json.load(f)
        except Exception:
            return None
    return None


def save_cached_smiles(
    cache_dir: Path, source: str, key_dict: Dict, payload: Dict
) -> None:
    path = _make_cache_path(cache_dir, source, key_dict)
    tmp = tempfile.NamedTemporaryFile("w", delete=False, dir=cache_dir, suffix=".tmp")
    try:
        json.dump(payload, tmp, indent=2)
        tmp.flush()
        tmp.close()
        Path(tmp.name).replace(path)
    finally:
        try:
            Path(tmp.name).unlink(missing_ok=True)
        except Exception:
            pass


def _get_json(
    url: str,
    params: Dict,
    headers: Dict,
    timeout: int,
    retries: int,
    session: Optional[requests.Session] = None,
    request_log: Optional[Callable[[str], None]] = None,
) -> Dict:
    last_err = None
    sess = session or requests.Session()
    for attempt in range(retries + 1):
        try:
            resolved = _resolved_url(url, params)
            if request_log:
                request_log(resolved)
            resp = sess.get(url, params=params, headers=headers, timeout=timeout)
            if resp.status_code == 200:
                return resp.json()
            last_err = f"status={resp.status_code}"
        except Exception as exc:
            last_err = str(exc)
    raise RuntimeError(last_err or "request failed")


def _get_json_paged(
    url: str,
    params: Dict[str, Any],
    headers: Dict,
    timeout: int,
    retries: int,
    *,
    session: Optional[requests.Session] = None,
    item_key: str,
    max_pages: int = 5,
    max_items: Optional[int] = None,
    audit: Optional[SourceAudit] = None,
    source: str = "",
    purpose: str = "",
    request_logger: Optional[Callable[[str], None]] = None,
) -> Tuple[List[Any], List[Dict[str, Any]]]:
    collected: List[Any] = []
    page_metas: List[Dict[str, Any]] = []
    next_url = url
    next_params: Dict[str, Any] = dict(params or {})
    page_count = 0

    while next_url and page_count < max_pages:
        raw_data = _get_json(
            next_url,
            next_params,
            headers,
            timeout,
            retries,
            session=session,
            request_log=request_logger,
        )
        data = _safe_dict(raw_data)
        page_meta = _safe_dict(_safe_get(data, "page_meta", {}))
        response_meta = {
            "page": page_meta.get("page") or page_count + 1,
            "total_count": page_meta.get("total_count"),
            "next": page_meta.get("next"),
            "offset": page_meta.get("offset"),
            "limit": page_meta.get("limit"),
        }
        _record_audit(
            audit,
            source or "chembl",
            purpose or f"{source}.paged",
            next_url,
            params=next_params if page_count == 0 else {},
            response_meta=response_meta,
        )

        entries_raw = (
            _safe_get(data, item_key) or _safe_get(data, item_key.rstrip("s")) or []
        )
        entries = _safe_list(entries_raw)
        for entry in entries:
            collected.append(entry)
            if max_items is not None and len(collected) >= max_items:
                page_metas.append(page_meta)
                return collected, page_metas
        page_metas.append(page_meta)

        next_link = page_meta.get("next")
        if not next_link:
            break
        next_url = next_link
        next_params = {}
        page_count += 1

    return collected, page_metas


def _limit_smiles(smiles: List[str], limit: int) -> List[str]:
    if limit <= 0:
        return smiles
    return smiles[:limit]


def _normalize_debug(debug_data: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    debug: Dict[str, Any] = {
        "target_ids_sample": [],
        "assay_ids_sample": [],
        "skips": {},
        "activity_types_seen_top": [],
    }
    if isinstance(debug_data, dict):
        for key in debug:
            if key in debug_data:
                debug[key] = debug_data[key]
    return debug


def _serialize_provenance_map(prov: Dict[str, Set[str]]) -> Dict[str, List[str]]:
    return {smi: sorted(labels) for smi, labels in prov.items()}


def _capture_rejection_sample(
    activity: Any,
    samples: Optional[List[Dict[str, Any]]],
    cap: int,
) -> None:
    if samples is None:
        return
    if cap is None:
        return
    try:
        cap_int = int(cap)
    except Exception:
        cap_int = 0
    if cap_int <= 0 or len(samples) >= cap_int:
        return
    safe_activity = _safe_dict(activity)
    samples.append(
        {
            "pchembl_value": safe_activity.get("pchembl_value"),
            "standard_value": safe_activity.get("standard_value"),
            "standard_units": safe_activity.get("standard_units"),
            "standard_relation": safe_activity.get("standard_relation"),
            "standard_type": safe_activity.get("standard_type"),
            "activity_chembl_id": safe_activity.get("activity_chembl_id"),
            "assay_chembl_id": safe_activity.get("assay_chembl_id"),
            "molecule_chembl_id": safe_activity.get("molecule_chembl_id"),
        }
    )


def _chembl_paginated(
    url: str,
    params: Dict,
    headers: Dict,
    timeout: int,
    retries: int,
    result_key: str,
    session: Optional[requests.Session] = None,
    limit: int = 200,
    max_pages: int = 50,
    request_logger: Optional[Callable[[str], None]] = None,
):
    """
    Iterate through a ChEMBL collection endpoint using limit/offset pagination.
    """
    offset = 0
    page_count = 0
    total_count = None
    while page_count < max_pages:
        page_params = dict(params or {})
        page_params.setdefault("limit", limit)
        page_params["offset"] = offset
        raw_data = _get_json(
            url,
            page_params,
            headers,
            timeout,
            retries,
            session=session,
            request_log=request_logger,
        )
        data = _safe_dict(raw_data)
        entries_raw = (
            _safe_get(data, result_key) or _safe_get(data, result_key.rstrip("s")) or []
        )
        entries = _safe_list(entries_raw)
        for entry in entries:
            yield entry
        page_meta = _safe_dict(_safe_get(data, "page_meta") or {})
        total_count = _safe_get(page_meta, "total_count", total_count)
        limit_used = _safe_get(page_meta, "limit") or page_params.get("limit") or limit
        offset = (_safe_get(page_meta, "offset") or offset) + limit_used
        page_count += 1
        if total_count is not None and offset >= total_count:
            break
        if not entries:
            break
        if len(entries) < limit_used:
            break


def parse_chembl_molecule(mol_json: Dict, max_phase: Optional[int]) -> Optional[str]:
    mol = _safe_dict(mol_json)
    structures = _safe_dict(mol.get("molecule_structures"))
    smi = structures.get("canonical_smiles")
    if smi is None:
        return None
    if max_phase is not None:
        max_phase_val = mol.get("max_phase")
        max_phase_num: Optional[float]
        try:
            max_phase_num = (
                float(max_phase_val)
                if isinstance(max_phase_val, (int, float))
                else None
            )
        except Exception:
            max_phase_num = None
        if max_phase_num is not None and max_phase_num > float(max_phase):
            return None
    return smi


def filter_chembl_activity(
    activity: Dict, cutoff_nm: int, allowed_types: List[str]
) -> bool:
    if not isinstance(activity, dict):
        return False
    stype = (activity.get("standard_type") or "").upper()
    if allowed_types and stype and stype not in [t.upper() for t in allowed_types]:
        return False
    units = str(activity.get("standard_units") or "").lower()
    val = activity.get("standard_value")
    if units == "nm":
        try:
            if val is not None and float(val) > float(cutoff_nm):
                return False
        except Exception:
            pass
    return True


def _label_chembl_activity(
    activity: Dict,
    thresholds: Dict[str, float],
    allowed_types: List[str],
    skip_counts: Optional[Dict[str, int]] = None,
    activity_types_seen: Optional[Counter] = None,
    standard_value_override: Optional[float] = None,
    labeling_telemetry: Optional[Dict[str, int]] = None,
    rejected_value_samples: Optional[List[Dict[str, Any]]] = None,
    rejected_value_samples_max: int = 25,
) -> Optional[str]:
    """
    Return one of {"strong","weak","non"} or None if the activity should be ignored.
    Prefers pchembl_value; falls back to standard_value (nM, "=").
    """
    activity = _safe_dict(activity)

    def _inc_labeling(key: str) -> None:
        if labeling_telemetry is not None:
            labeling_telemetry[key] = labeling_telemetry.get(key, 0) + 1

    allowed_upper = [t.upper() for t in allowed_types] if allowed_types else []
    stype = (activity.get("standard_type") or "").upper()
    if activity_types_seen is not None and stype:
        activity_types_seen[stype] += 1
    has_pchembl = activity.get("pchembl_value") is not None
    if not stype:
        if skip_counts is not None:
            skip_counts["skip_missing_standard_type"] = (
                skip_counts.get("skip_missing_standard_type", 0) + 1
            )
            skip_counts_key = "skip_has_pchembl" if has_pchembl else "skip_no_pchembl"
            skip_counts[skip_counts_key] = skip_counts.get(skip_counts_key, 0) + 1
        if not has_pchembl:
            _inc_labeling("n_rejected_by_missing_pchembl")
        return None
    if allowed_upper and stype not in allowed_upper:
        if skip_counts is not None:
            skip_counts["skip_type_not_allowed"] = (
                skip_counts.get("skip_type_not_allowed", 0) + 1
            )
            skip_counts_key = "skip_has_pchembl" if has_pchembl else "skip_no_pchembl"
            skip_counts[skip_counts_key] = skip_counts.get(skip_counts_key, 0) + 1
        _inc_labeling("n_rejected_by_type_not_allowed")
        if not has_pchembl:
            _inc_labeling("n_rejected_by_missing_pchembl")
        return None
    pchembl_val = activity.get("pchembl_value")
    pchembl_float: Optional[float] = None
    try:
        if pchembl_val is not None:
            pchembl_float = float(pchembl_val)
    except Exception:
        pchembl_float = None
        if skip_counts is not None:
            skip_counts["skip_has_pchembl"] = skip_counts.get("skip_has_pchembl", 0) + 1
    if pchembl_float is not None:
        if pchembl_float >= thresholds.get("pchembl_strong", 7.0):
            return "strong"
        if pchembl_float >= thresholds.get("pchembl_weak", 5.0):
            return "weak"
        return "non"
    relation = (activity.get("standard_relation") or "").strip()
    units = str(activity.get("standard_units") or "").lower()
    if relation != "=" or units != "nm":
        if skip_counts is not None:
            if relation != "=":
                skip_counts["skip_relation_not_equal"] = (
                    skip_counts.get("skip_relation_not_equal", 0) + 1
                )
            if units != "nm":
                skip_counts["skip_units_not_nm"] = (
                    skip_counts.get("skip_units_not_nm", 0) + 1
                )
            skip_counts_key = "skip_has_pchembl" if has_pchembl else "skip_no_pchembl"
            skip_counts[skip_counts_key] = skip_counts.get(skip_counts_key, 0) + 1
        if standard_value_override is not None:
            _inc_labeling("n_rejected_by_value_threshold")
            _capture_rejection_sample(
                activity, rejected_value_samples, rejected_value_samples_max
            )
        if relation != "=":
            _inc_labeling("n_rejected_by_relation")
        if units != "nm":
            _inc_labeling("n_rejected_by_units")
        if pchembl_float is None:
            _inc_labeling("n_rejected_by_missing_pchembl")
        return None
    try:
        raw_val = activity.get("standard_value")
        if standard_value_override is not None:
            val = standard_value_override
        else:
            val = float(raw_val) if raw_val is not None else None
    except Exception:
        if skip_counts is not None:
            skip_counts["skip_missing_standard_value"] = (
                skip_counts.get("skip_missing_standard_value", 0) + 1
            )
            skip_counts_key = "skip_has_pchembl" if has_pchembl else "skip_no_pchembl"
            skip_counts[skip_counts_key] = skip_counts.get(skip_counts_key, 0) + 1
        _inc_labeling("n_rejected_by_missing_standard_value")
        if pchembl_float is None:
            _inc_labeling("n_rejected_by_missing_pchembl")
        return None
    if val is None:
        if skip_counts is not None:
            skip_counts["skip_missing_standard_value"] = (
                skip_counts.get("skip_missing_standard_value", 0) + 1
            )
            skip_counts_key = "skip_has_pchembl" if has_pchembl else "skip_no_pchembl"
            skip_counts[skip_counts_key] = skip_counts.get(skip_counts_key, 0) + 1
        _inc_labeling("n_rejected_by_missing_standard_value")
        if pchembl_float is None:
            _inc_labeling("n_rejected_by_missing_pchembl")
        return None
    strong_cutoff = thresholds.get("standard_value_nm_strong", 100.0)
    weak_cutoff = thresholds.get("standard_value_nm_weak", 10000.0)
    if val <= strong_cutoff:
        return "strong"
    if val <= weak_cutoff:
        return "weak"
    return "non"


def parse_iuphar_structure(struct_json: Dict) -> Optional[str]:
    if not isinstance(struct_json, dict):
        return None
    smi = struct_json.get("smiles")
    if smi:
        return smi
    inchi = struct_json.get("inchi")
    if inchi:
        try:
            from rdkit import Chem

            mol = Chem.MolFromInchi(inchi)
            if mol:
                return Chem.MolToSmiles(mol)
        except Exception:
            return None
    return None


def fetch_chembl_smiles(
    uniprot_id: str,
    ec_numbers: List[str],
    pdb_id: str,
    cache_dir: Path,
    timeout: int,
    retries: int,
    max_actives: int,
    affinity_cutoff_nm: int,
    activity_types: List[str],
    max_phase: Optional[int],
    session: Optional[requests.Session] = None,
    audit: Optional[SourceAudit] = None,
    max_pages: int = 5,
) -> Tuple[List[str], Dict]:
    source = "chembl"
    key = {
        "source": source,
        "uniprot": uniprot_id,
        "ecs": sorted(ec_numbers),
        "pdb": pdb_id,
        "types": activity_types,
        "cutoff": affinity_cutoff_nm,
        "max_phase": max_phase,
        "max_pages": max_pages,
    }
    cached = load_cached_smiles(cache_dir, source, key)
    if cached and "smiles" in cached:
        cached_meta = cached.get("meta", {})
        cached_meta = cached_meta if isinstance(cached_meta, dict) else {}
        cached_meta = dict(cached_meta)
        cached_meta["cached"] = True
        return cached["smiles"], cached_meta

    base = "https://www.ebi.ac.uk/chembl/api/data"
    headers = {"Accept": "application/json"}
    smiles: List[str] = []
    seen_smiles: Set[str] = set()
    meta_counts = {"targets": 0, "assays": 0, "activities": 0, "molecules": 0}
    provenance: Dict[str, Set[str]] = {}
    provenance_details: Dict[str, List[Dict[str, Any]]] = {}
    sess = session or requests.Session()
    owns_session = session is None
    max_pages = max(1, int(max_pages or 1))
    max_activity_items = max_actives * 5 if max_actives > 0 else None

    def add_provenance(smiles_val: str, label: str, detail: Optional[Dict[str, Any]]):
        if not smiles_val:
            return
        provenance.setdefault(smiles_val, set()).add(label)
        if detail:
            provenance_details.setdefault(smiles_val, []).append(detail)

    error_meta: Optional[Dict[str, Any]] = None
    try:
        target_params = {
            "target_components__accession": uniprot_id,
            "format": "json",
            "limit": 200,
        }
        targets, _ = _get_json_paged(
            f"{base}/target",
            target_params,
            headers,
            timeout,
            retries,
            session=sess,
            item_key="targets",
            max_pages=max_pages,
            audit=audit,
            source=source,
            purpose="chembl.target_search",
        )
        target_ids = []
        for t in targets:
            t_obj = _safe_dict(t)
            tid = t_obj.get("target_chembl_id")
            if tid and tid not in target_ids:
                target_ids.append(tid)
        meta_counts["targets"] = len(target_ids)
        for tid in target_ids:
            assay_params = {
                "target_chembl_id": tid,
                "assay_type": "B",
                "relationship_type": "D",
                "format": "json",
                "limit": 200,
            }
            assays, _ = _get_json_paged(
                f"{base}/assay",
                assay_params,
                headers,
                timeout,
                retries,
                session=sess,
                item_key="assays",
                max_pages=max_pages,
                audit=audit,
                source=source,
                purpose="chembl.assays_for_target",
            )
            assay_ids = []
            for a in assays:
                assay_obj = _safe_dict(a)
                aid = assay_obj.get("assay_chembl_id")
                if aid and aid not in assay_ids:
                    assay_ids.append(aid)
            meta_counts["assays"] += len(assay_ids)
            for aid in assay_ids:
                if max_actives > 0 and len(smiles) >= max_actives:
                    break
                activity_params = {
                    "assay_chembl_id": aid,
                    "format": "json",
                    "limit": 200,
                }
                if activity_types:
                    activity_params["standard_type__in"] = ",".join(activity_types)
                activities, _ = _get_json_paged(
                    f"{base}/activity",
                    activity_params,
                    headers,
                    timeout,
                    retries,
                    session=sess,
                    item_key="activities",
                    max_pages=max_pages,
                    max_items=max_activity_items,
                    audit=audit,
                    source=source,
                    purpose="chembl.activities_for_assay",
                )
                meta_counts["activities"] += len(activities)
                for act in activities:
                    if max_actives > 0 and len(smiles) >= max_actives:
                        break
                    if not filter_chembl_activity(
                        act, affinity_cutoff_nm, activity_types
                    ):
                        continue
                    mol_id = act.get("molecule_chembl_id")
                    if mol_id is None:
                        continue
                    mol_url = f"{base}/molecule/{mol_id}.json"
                    mol_json = _get_json(
                        mol_url,
                        {},
                        headers,
                        timeout,
                        retries,
                        session=sess,
                    )
                    _record_audit(
                        audit,
                        source,
                        "chembl.molecule_fetch",
                        mol_url,
                        response_meta={"molecule_chembl_id": mol_id},
                    )
                    smi = parse_chembl_molecule(mol_json, max_phase)
                    if not smi or smi in seen_smiles:
                        continue
                    seen_smiles.add(smi)
                    smiles.append(smi)
                    meta_counts["molecules"] += 1
                    label_parts = []
                    if tid:
                        label_parts.append(f"target={tid}")
                    if aid:
                        label_parts.append(f"assay={aid}")
                    act_id = act.get("activity_chembl_id")
                    if act_id:
                        label_parts.append(f"activity={act_id}")
                    label_parts.append(f"molecule={mol_id}")
                    stype = act.get("standard_type")
                    sval = act.get("standard_value")
                    sunits = act.get("standard_units")
                    if stype:
                        type_piece = stype
                        if sval is not None:
                            type_piece += f"={sval}{sunits or ''}"
                        label_parts.append(type_piece)
                    label = "chembl:" + ";".join(label_parts)
                    detail = {
                        "source": source,
                        "target_chembl_id": tid,
                        "assay_chembl_id": aid,
                        "activity_chembl_id": act_id,
                        "molecule_chembl_id": mol_id,
                        "standard_type": act.get("standard_type"),
                        "standard_relation": act.get("standard_relation"),
                        "standard_value": act.get("standard_value"),
                        "standard_units": act.get("standard_units"),
                        "label": label,
                    }
                    add_provenance(smi, label, detail)
                if max_actives > 0 and len(smiles) >= max_actives:
                    break
            if max_actives > 0 and len(smiles) >= max_actives:
                break
    except Exception as exc:
        error_meta = {"error": str(exc), "counts": meta_counts}
    finally:
        if owns_session:
            try:
                sess.close()
            except Exception:
                pass

    if error_meta:
        return [], error_meta

    smiles = _limit_smiles(smiles, max_actives)
    result_meta: Dict[str, Any] = {"cached": False, "counts": meta_counts}
    if provenance:
        result_meta["provenance"] = _serialize_provenance_map(provenance)
    if provenance_details:
        result_meta["provenance_details"] = provenance_details
    save_cached_smiles(cache_dir, source, key, {"smiles": smiles, "meta": result_meta})
    return smiles, result_meta


def fetch_chembl_labeled_smiles(
    uniprot_id: str,
    pdb_id: str,
    cache_dir: Path,
    timeout: int,
    retries: int,
    activity_types: List[str],
    chembl_max_phase: Optional[int],
    label_thresholds: Dict[str, float],
    session: Optional[requests.Session] = None,
    debug: bool = False,
    debug_max_ids: int = 25,
    debug_rejection_samples_max: int = 25,
) -> Tuple[Dict[str, List[str]], Dict]:
    """
    Retrieve ChEMBL activities for a UniProt target and bucket SMILES into
    strong/weak/non binders.
    """
    source = "chembl_labeled"
    cleaned_thresholds = {
        "pchembl_strong": float(label_thresholds.get("pchembl_strong", 7.0)),
        "pchembl_weak": float(label_thresholds.get("pchembl_weak", 5.0)),
        "standard_value_nm_strong": float(
            label_thresholds.get("standard_value_nm_strong", 100.0)
        ),
        "standard_value_nm_weak": float(
            label_thresholds.get("standard_value_nm_weak", 10000.0)
        ),
    }
    activity_types = list(activity_types or [])
    key = {
        "source": source,
        "uniprot": uniprot_id,
        "pdb": pdb_id,
        "activity_types": sorted(activity_types),
        "chembl_max_phase": chembl_max_phase,
        "thresholds": cleaned_thresholds,
        "debug": debug,
        "debug_max_ids": debug_max_ids,
        "debug_rejection_samples_max": debug_rejection_samples_max,
    }
    cached = load_cached_smiles(cache_dir, source, key)
    if cached and "labels" in cached:
        cached_meta = cached.get("meta", {})
        cached_meta["cached"] = True
        cached_meta["debug"] = _normalize_debug(cached_meta.get("debug"))
        cached_meta.setdefault("telemetry", {})
        labeling = cached_meta["telemetry"].get("labeling_sanity")
        if isinstance(labeling, dict):
            labeling.setdefault("rejected_value_threshold_samples", [])
            labeling.setdefault(
                "rejected_value_threshold_samples_max", debug_rejection_samples_max
            )
        return cached["labels"], cached_meta

    base = "https://www.ebi.ac.uk/chembl/api/data"
    headers = {"Accept": "application/json"}
    labels: Dict[str, List[str]] = {"strong": [], "weak": [], "non": []}
    label_map: Dict[str, str] = {}
    label_priority = {"non": 1, "weak": 2, "strong": 3}
    counts = {"targets": 0, "assays": 0, "activities": 0, "molecules": 0}
    molecule_cache: Dict[str, Optional[str]] = {}
    error = None
    debug_info = _normalize_debug({})
    target_ids_sample = cast(List[str], debug_info.get("target_ids_sample", []))
    assay_ids_sample = cast(List[str], debug_info.get("assay_ids_sample", []))
    skip_counts = cast(Dict[str, int], debug_info.get("skips", {}))
    activity_types_seen: Counter = Counter()
    debug_info["target_ids_sample"] = target_ids_sample
    debug_info["assay_ids_sample"] = assay_ids_sample
    debug_info["skips"] = skip_counts
    telemetry_max_urls = 200
    request_urls: List[str] = []
    telemetry: Dict[str, Any] = {
        "activity_sanity": {
            "n_activities_total": 0,
            "n_with_molecule_chembl_id": 0,
            "n_with_pchembl_value": 0,
            "relation_histogram": {},
            "units_histogram": {},
            "n_standard_value_parseable": 0,
        },
        "labeling_sanity": {
            "n_labeled_strong": 0,
            "n_labeled_weak": 0,
            "n_labeled_non": 0,
            "n_rejected_by_relation": 0,
            "n_rejected_by_units": 0,
            "n_rejected_by_missing_pchembl": 0,
            "n_rejected_by_missing_standard_value": 0,
            "n_rejected_by_value_threshold": 0,
            "n_rejected_by_type_not_allowed": 0,
            "n_rejected_by_phase_filtered": 0,
            "rejected_value_threshold_samples": [],
            "rejected_value_threshold_samples_max": debug_rejection_samples_max,
        },
        "molecule_sanity": {
            "n_molecule_fetch_attempted": 0,
            "n_molecule_fetch_failed_http": 0,
            "n_molecule_missing_smiles": 0,
            "n_molecule_parsed_smiles_ok": 0,
        },
        "request_urls": request_urls,
        "request_urls_truncated": False,
    }
    activity_sanity = telemetry["activity_sanity"]
    labeling_sanity = telemetry["labeling_sanity"]
    molecule_sanity = telemetry["molecule_sanity"]
    relation_hist: Counter = Counter()
    units_hist: Counter = Counter()
    rejected_value_samples: List[Dict[str, Any]] = labeling_sanity.get(
        "rejected_value_threshold_samples", []
    )

    def _normalize_relation(rel_val: Any) -> str:
        relation = str(rel_val).strip() if rel_val is not None else ""
        return relation if relation else "<missing>"

    def _normalize_units(units_val: Any) -> str:
        if units_val is None:
            return "<missing>"
        unit_str = str(units_val).strip()
        if not unit_str:
            return "<missing>"
        lower = unit_str.lower()
        if lower in {
            "um",
            "\u00b5m",
            "\u03bcm",
            "micromolar",
            "micromole",
            "micromoles",
            "microm",
        }:
            return "uM"
        if lower == "nm":
            return "nM"
        return unit_str

    def log_request_url(url: Optional[str]) -> None:
        if not url:
            return
        resolved = str(url)
        if len(request_urls) < telemetry_max_urls:
            request_urls.append(resolved)
        else:
            telemetry["request_urls_truncated"] = True

    sess = session or requests.Session()
    owns_session = session is None

    try:
        target_ids = []
        target_params = {"target_components__accession": uniprot_id, "format": "json"}
        for target in _chembl_paginated(
            f"{base}/target",
            target_params,
            headers,
            timeout,
            retries,
            "targets",
            session=sess,
            request_logger=log_request_url,
        ):
            target_obj = _safe_dict(target)
            if not target_obj:
                _increment(skip_counts, "skip_payload_not_dict")
                continue
            tid = target_obj.get("target_chembl_id")
            if tid and tid not in target_ids:
                target_ids.append(tid)
                if len(target_ids_sample) < debug_max_ids:
                    target_ids_sample.append(tid)
        counts["targets"] = len(target_ids)

        assay_ids = []
        for tid in target_ids:
            for assay in _chembl_paginated(
                f"{base}/assay",
                {
                    "target_chembl_id": tid,
                    "assay_type": "B",
                    "relationship_type": "D",
                    "format": "json",
                },
                headers,
                timeout,
                retries,
                "assays",
                session=sess,
                request_logger=log_request_url,
            ):
                assay_obj = _safe_dict(assay)
                if not assay_obj:
                    _increment(skip_counts, "skip_payload_not_dict")
                    continue
                aid = assay_obj.get("assay_chembl_id")
                if aid and aid not in assay_ids:
                    assay_ids.append(aid)
                    if len(assay_ids_sample) < debug_max_ids:
                        assay_ids_sample.append(aid)
        counts["assays"] = len(assay_ids)

        def assign_label(smiles: str, label: str):
            current = label_map.get(smiles)
            if current and label_priority[current] >= label_priority[label]:
                return
            label_map[smiles] = label

        for aid in assay_ids:
            activity_params = {
                "assay_chembl_id": aid,
                "format": "json",
            }
            if activity_types:
                activity_params["standard_type__in"] = ",".join(activity_types)
            activity_seen = False
            for activity in _chembl_paginated(
                f"{base}/activity",
                activity_params,
                headers,
                timeout,
                retries,
                "activities",
                session=sess,
                limit=500,
                request_logger=log_request_url,
            ):
                activity_obj = _safe_dict(activity)
                if not activity_obj:
                    _increment(skip_counts, "skip_payload_not_dict")
                    continue
                activity_seen = True
                counts["activities"] += 1
                relation_hist[
                    _normalize_relation(activity_obj.get("standard_relation"))
                ] += 1
                units_hist[_normalize_units(activity_obj.get("standard_units"))] += 1

                std_val_float: Optional[float] = None
                try:
                    raw_std_val = activity_obj.get("standard_value")
                    if raw_std_val is not None:
                        std_val_float = float(raw_std_val)
                        activity_sanity["n_standard_value_parseable"] += 1
                except Exception:
                    std_val_float = None

                if activity_obj.get("molecule_chembl_id"):
                    activity_sanity["n_with_molecule_chembl_id"] += 1
                try:
                    pchembl_val = activity_obj.get("pchembl_value")
                    if pchembl_val is not None:
                        float(pchembl_val)
                        activity_sanity["n_with_pchembl_value"] += 1
                except Exception:
                    pass
                label = _label_chembl_activity(
                    activity_obj,
                    cleaned_thresholds,
                    activity_types,
                    skip_counts,
                    activity_types_seen,
                    standard_value_override=std_val_float,
                    labeling_telemetry=labeling_sanity,
                    rejected_value_samples=rejected_value_samples,
                    rejected_value_samples_max=debug_rejection_samples_max,
                )
                if not label:
                    continue
                mol_id = activity_obj.get("molecule_chembl_id")
                if not mol_id:
                    _increment(skip_counts, "skip_missing_molecule_chembl_id")
                    continue
                phase_filtered = False
                if mol_id in molecule_cache:
                    smi = molecule_cache[mol_id]
                else:
                    molecule_sanity["n_molecule_fetch_attempted"] += 1
                    try:
                        mol_json = _get_json(
                            f"{base}/molecule/{mol_id}.json",
                            {},
                            headers,
                            timeout,
                            retries,
                            session=sess,
                            request_log=log_request_url,
                        )
                    except Exception:
                        molecule_sanity["n_molecule_fetch_failed_http"] += 1
                        molecule_cache[mol_id] = None
                        continue
                    mol_payload = _safe_dict(mol_json)
                    if not isinstance(mol_json, dict):
                        _increment(skip_counts, "skip_payload_not_dict")
                    mol_structures = _safe_dict(mol_payload.get("molecule_structures"))
                    if not mol_structures:
                        _increment(skip_counts, "skip_molecule_structures_null")
                    max_phase_val = mol_payload.get("max_phase")
                    if chembl_max_phase is not None:
                        phase_value: Optional[float]
                        try:
                            phase_value = (
                                float(max_phase_val)
                                if isinstance(max_phase_val, (int, float))
                                else None
                            )
                        except Exception:
                            phase_value = None
                        if phase_value is None and max_phase_val is not None:
                            _increment(skip_counts, "skip_phase_unknown")
                        if phase_value is not None and phase_value > float(
                            chembl_max_phase
                        ):
                            phase_filtered = True
                    smi = (
                        mol_structures.get("canonical_smiles")
                        if mol_structures
                        else None
                    )
                    if not smi:
                        smi = parse_chembl_molecule(mol_payload, chembl_max_phase)
                    molecule_cache[mol_id] = smi
                if not smi:
                    if phase_filtered:
                        _increment(skip_counts, "skip_phase_filtered")
                        labeling_sanity["n_rejected_by_phase_filtered"] = (
                            labeling_sanity.get("n_rejected_by_phase_filtered", 0) + 1
                        )
                    else:
                        _increment(skip_counts, "skip_molecule_no_smiles")
                        molecule_sanity["n_molecule_missing_smiles"] = (
                            molecule_sanity.get("n_molecule_missing_smiles", 0) + 1
                        )
                    continue
                molecule_sanity["n_molecule_parsed_smiles_ok"] = (
                    molecule_sanity.get("n_molecule_parsed_smiles_ok", 0) + 1
                )
                assign_label(smi, label)
            if not activity_seen:
                _increment(skip_counts, "skip_empty_activity_page")
        for smi, lbl in label_map.items():
            labels[lbl].append(smi)
        counts["molecules"] = len(label_map)
    except Exception as exc:
        error = str(exc)
    finally:
        if owns_session:
            try:
                sess.close()
            except Exception:
                pass

    activity_sanity["n_activities_total"] = counts.get("activities", 0)
    activity_sanity["relation_histogram"] = dict(relation_hist)
    activity_sanity["units_histogram"] = dict(units_hist)

    if activity_types_seen:
        debug_info["activity_types_seen_top"] = [
            {"type": t, "count": c} for t, c in activity_types_seen.most_common(10)
        ]
    bin_counts = {k: len(v) for k, v in labels.items()}
    labeling_sanity["n_labeled_strong"] = bin_counts.get("strong", 0)
    labeling_sanity["n_labeled_weak"] = bin_counts.get("weak", 0)
    labeling_sanity["n_labeled_non"] = bin_counts.get("non", 0)
    result_meta: Dict[str, object] = {
        "cached": False,
        "counts": {**counts, **bin_counts},
        "debug": _normalize_debug(debug_info),
        "telemetry": telemetry,
    }
    if error:
        result_meta["error"] = error
        return labels, result_meta

    save_cached_smiles(cache_dir, source, key, {"labels": labels, "meta": result_meta})
    return labels, result_meta


def _parse_bindingdb_entries(
    entries, affinity_cutoff_nm: int
) -> List[Tuple[str, Dict]]:
    results: List[Tuple[str, Dict]] = []
    seen: Set[str] = set()
    for entry in entries:
        smi = entry.get("smiles") or entry.get("LigandSMILES") or entry.get("SMILES")
        if not smi:
            continue
        affinity = (
            entry.get("Affinity_nM")
            or entry.get("affinity")
            or entry.get("Kd")
            or entry.get("Ki")
        )
        try:
            if affinity is not None and float(affinity) > float(affinity_cutoff_nm):
                continue
        except Exception:
            pass
        if smi in seen:
            continue
        seen.add(smi)
        results.append((smi, entry))
    return results


def fetch_bindingdb_smiles(
    uniprot_id: str,
    ec_numbers: List[str],
    pdb_id: str,
    cache_dir: Path,
    timeout: int,
    retries: int,
    max_actives: int,
    affinity_cutoff_nm: int,
    session: Optional[requests.Session] = None,
    audit: Optional[SourceAudit] = None,
) -> Tuple[List[str], Dict]:
    source = "bindingdb"
    key = {
        "source": source,
        "uniprot": uniprot_id,
        "ecs": sorted(ec_numbers),
        "pdb": pdb_id,
        "cutoff": affinity_cutoff_nm,
    }
    cached = load_cached_smiles(cache_dir, source, key)
    if cached and "smiles" in cached:
        cached_meta = cached.get("meta", {})
        cached_meta = cached_meta if isinstance(cached_meta, dict) else {}
        cached_meta = dict(cached_meta)
        cached_meta["cached"] = True
        return cached["smiles"], cached_meta

    headers = {"Accept": "application/json"}
    base = "https://bindingdb.org"
    meta = {"pdb_hits": 0, "uniprot_hits": 0}
    smiles: List[str] = []
    seen: Set[str] = set()
    provenance: Dict[str, Set[str]] = {}
    provenance_details: Dict[str, List[Dict[str, Any]]] = {}
    sess = session or requests.Session()
    owns_session = session is None

    def add_provenance(smiles_val: str, label: str, detail: Optional[Dict[str, Any]]):
        provenance.setdefault(smiles_val, set()).add(label)
        if detail:
            provenance_details.setdefault(smiles_val, []).append(detail)

    def fetch(url, params, counter_key, query_value):
        try:
            data = _get_json(url, params, headers, timeout, retries, session=sess)
            _record_audit(
                audit,
                source,
                f"bindingdb.{counter_key}",
                url,
                params=params,
                response_meta={"query": query_value},
            )
            if isinstance(data, dict) and "results" in data:
                entries = data.get("results") or []
            else:
                entries = data if isinstance(data, list) else []
            parsed = _parse_bindingdb_entries(entries, affinity_cutoff_nm)
            meta[counter_key] += len(parsed)
            for smi, entry in parsed:
                if len(smiles) >= max_actives:
                    break
                if smi not in seen:
                    smiles.append(smi)
                    seen.add(smi)
                detail = {
                    "source": source,
                    "query": counter_key,
                    "query_value": query_value,
                }
                lig_id = (
                    entry.get("LigandID")
                    or entry.get("ligandId")
                    or entry.get("BindingDBMonomerID")
                )
                if lig_id is not None:
                    detail["ligand_id"] = lig_id
                label = f"bindingdb:{counter_key}"
                detail["label"] = label
                add_provenance(smi, label, detail)
        except Exception:
            return

    fetch(
        f"{base}/rest/getLigandsByPDBs",
        {
            "pdb": pdb_id,
            "cutoff": affinity_cutoff_nm,
            "identity": 100,
            "response": "application/json",
        },
        "pdb_hits",
        pdb_id,
    )
    if len(smiles) < max_actives:
        fetch(
            f"{base}/rest/getLigandsByUniprots",
            {
                "uniprot": uniprot_id,
                "cutoff": affinity_cutoff_nm,
                "response": "application/json",
            },
            "uniprot_hits",
            uniprot_id,
        )

    if owns_session:
        try:
            sess.close()
        except Exception:
            pass

    smiles = _limit_smiles(smiles, max_actives)
    result_meta: Dict[str, Any] = {"cached": False, "counts": meta}
    if provenance:
        result_meta["provenance"] = _serialize_provenance_map(provenance)
    if provenance_details:
        result_meta["provenance_details"] = provenance_details
    save_cached_smiles(cache_dir, source, key, {"smiles": smiles, "meta": result_meta})
    return smiles, result_meta


def fetch_iuphar_smiles(
    uniprot_id: str,
    ec_numbers: List[str],
    pdb_id: str,
    cache_dir: Path,
    timeout: int,
    retries: int,
    max_actives: int,
    approved_only: bool,
    primary_only: bool,
    session: Optional[requests.Session] = None,
    audit: Optional[SourceAudit] = None,
) -> Tuple[List[str], Dict]:
    source = "iuphar"
    key = {
        "source": source,
        "uniprot": uniprot_id,
        "ecs": sorted(ec_numbers),
        "pdb": pdb_id,
        "approved": approved_only,
        "primary": primary_only,
    }
    cached = load_cached_smiles(cache_dir, source, key)
    if cached and "smiles" in cached:
        cached_meta = cached.get("meta", {})
        cached_meta = cached_meta if isinstance(cached_meta, dict) else {}
        cached_meta = dict(cached_meta)
        cached_meta["cached"] = True
        return cached["smiles"], cached_meta

    base = "https://www.guidetopharmacology.org/services"
    headers = {"Accept": "application/json"}
    meta = {"targets": 0, "interactions": 0}
    smiles: List[str] = []
    seen: Set[str] = set()
    provenance: Dict[str, Set[str]] = {}
    provenance_details: Dict[str, List[Dict[str, Any]]] = {}

    def add_provenance(smiles_val: str, label: str, detail: Optional[Dict[str, Any]]):
        provenance.setdefault(smiles_val, set()).add(label)
        if detail:
            provenance_details.setdefault(smiles_val, []).append(detail)

    def get_targets():
        targets = []
        try:
            tgt_json = _get_json(
                f"{base}/targets",
                {"accession": uniprot_id},
                headers,
                timeout,
                retries,
                session=session,
            )
            _record_audit(
                audit,
                source,
                "iuphar.targets_by_accession",
                f"{base}/targets",
                params={"accession": uniprot_id},
            )
            if isinstance(tgt_json, list):
                targets.extend(tgt_json)
        except Exception:
            pass
        for ec in ec_numbers:
            if len(targets) >= max_actives:
                break
            try:
                params = {"ecNumber": ec}
                tgt_json = _get_json(
                    f"{base}/targets",
                    params,
                    headers,
                    timeout,
                    retries,
                    session=session,
                )
                _record_audit(
                    audit,
                    source,
                    "iuphar.targets_by_ec",
                    f"{base}/targets",
                    params=params,
                )
                if isinstance(tgt_json, list):
                    targets.extend(tgt_json)
            except Exception:
                continue
        return targets

    targets = get_targets()
    target_ids = []
    for t in targets:
        tid = t.get("targetId")
        if tid and tid not in target_ids:
            target_ids.append(tid)
    meta["targets"] = len(target_ids)

    for tid in target_ids:
        if len(smiles) >= max_actives:
            break
        params = {
            "approved": str(approved_only).lower(),
            "primaryTarget": str(primary_only).lower(),
        }
        try:
            _record_audit(
                audit,
                source,
                "iuphar.interactions",
                f"{base}/targets/{tid}/interactions",
                params=params,
            )
            interactions = _get_json(
                f"{base}/targets/{tid}/interactions",
                params,
                headers,
                timeout,
                retries,
                session=session,
            )
        except Exception:
            continue
        if not isinstance(interactions, list):
            continue
        meta["interactions"] += len(interactions)
        for inter in interactions:
            if len(smiles) >= max_actives:
                break
            ligand_id = inter.get("ligandId")
            if ligand_id is None:
                continue
            try:
                struct = _get_json(
                    f"{base}/ligands/{ligand_id}/structure",
                    {},
                    headers,
                    timeout,
                    retries,
                    session=session,
                )
                _record_audit(
                    audit,
                    source,
                    "iuphar.ligand_structure",
                    f"{base}/ligands/{ligand_id}/structure",
                )
            except Exception:
                continue
            smi = struct.get("smiles")
            if not smi:
                inchi = struct.get("inchi")
                if inchi:
                    try:
                        from rdkit import Chem

                        mol = Chem.MolFromInchi(inchi)
                        if mol:
                            smi = Chem.MolToSmiles(mol)
                    except Exception:
                        smi = None
            if smi:
                label_parts = [f"target={tid}", f"ligand={ligand_id}"]
                interaction_id = inter.get("interactionId")
                if interaction_id is not None:
                    label_parts.append(f"interaction={interaction_id}")
                label = "iuphar:" + ";".join(label_parts)
                detail = {
                    "source": source,
                    "targetId": tid,
                    "ligandId": ligand_id,
                }
                if interaction_id is not None:
                    detail["interactionId"] = interaction_id
                if smi not in seen:
                    smiles.append(smi)
                    seen.add(smi)
                detail["label"] = label
                add_provenance(smi, label, detail)
    smiles = _limit_smiles(smiles, max_actives)
    result_meta: Dict[str, Any] = {"cached": False, "counts": meta}
    if provenance:
        result_meta["provenance"] = _serialize_provenance_map(provenance)
    if provenance_details:
        result_meta["provenance_details"] = provenance_details
    save_cached_smiles(cache_dir, source, key, {"smiles": smiles, "meta": result_meta})
    return smiles, result_meta


def fetch_drugbank_smiles(
    uniprot_id: str,
    ec_numbers: List[str],
    pdb_id: str,
    cache_dir: Path,
    data_dir: Optional[Path],
) -> Tuple[List[str], Dict]:
    source = "drugbank"
    key = {
        "source": source,
        "uniprot": uniprot_id,
        "ecs": sorted(ec_numbers),
        "pdb": pdb_id,
        "data_dir": str(data_dir) if data_dir else "",
    }
    cached = load_cached_smiles(cache_dir, source, key)
    if cached and "smiles" in cached:
        return cached["smiles"], {"cached": True, "counts": cached.get("meta", {})}
    if not data_dir or not data_dir.exists():
        return [], {"skipped": True, "reason": "no_data_dir"}
    # Placeholder: no local DrugBank parser implemented; skip with cache.
    save_cached_smiles(
        cache_dir, source, key, {"smiles": [], "meta": {"skipped": True}}
    )
    return [], {"skipped": True, "reason": "not_implemented"}
