import hashlib
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import tempfile
import requests


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
) -> Dict:
    last_err = None
    sess = session or requests.Session()
    for attempt in range(retries + 1):
        try:
            resp = sess.get(url, params=params, headers=headers, timeout=timeout)
            if resp.status_code == 200:
                return resp.json()
            last_err = f"status={resp.status_code}"
        except Exception as exc:
            last_err = str(exc)
    raise RuntimeError(last_err or "request failed")


def _limit_smiles(smiles: List[str], limit: int) -> List[str]:
    if limit <= 0:
        return smiles
    return smiles[:limit]


def parse_chembl_molecule(mol_json: Dict, max_phase: Optional[int]) -> Optional[str]:
    structures = (
        mol_json.get("molecule_structures", {}) if isinstance(mol_json, dict) else {}
    )
    smi = structures.get("canonical_smiles")
    if smi is None:
        return None
    if max_phase is not None:
        try:
            if int(mol_json.get("max_phase", 0)) > int(max_phase):
                return None
        except Exception:
            pass
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
            if float(val) > float(cutoff_nm):
                return False
        except Exception:
            pass
    return True


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
    }
    cached = load_cached_smiles(cache_dir, source, key)
    if cached and "smiles" in cached:
        return cached["smiles"], {"cached": True, "counts": cached.get("meta", {})}

    base = "https://www.ebi.ac.uk/chembl/api/data"
    headers = {"Accept": "application/json"}
    smiles: List[str] = []
    meta = {"targets": 0, "assays": 0, "activities": 0, "molecules": 0}

    try:
        targets_json = _get_json(
            f"{base}/target",
            {"target_components__accession": uniprot_id, "format": "json"},
            headers,
            timeout,
            retries,
            session=session,
        )
        targets = targets_json.get("targets") or targets_json.get("target") or []
        target_ids = [
            t.get("target_chembl_id") for t in targets if t.get("target_chembl_id")
        ]
        meta["targets"] = len(target_ids)
        for tid in target_ids:
            assays_json = _get_json(
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
                session=session,
            )
            assays = assays_json.get("assays") or assays_json.get("assay") or []
            assay_ids = [
                a.get("assay_chembl_id") for a in assays if a.get("assay_chembl_id")
            ]
            meta["assays"] += len(assay_ids)
            for aid in assay_ids:
                activities_json = _get_json(
                    f"{base}/activity",
                    {
                        "assay_chembl_id": aid,
                        "standard_type__in": ",".join(activity_types),
                        "format": "json",
                    },
                    headers,
                    timeout,
                    retries,
                    session=session,
                )
                activities = (
                    activities_json.get("activities")
                    or activities_json.get("activity")
                    or []
                )
                meta["activities"] += len(activities)
                for act in activities:
                    if meta["molecules"] >= max_actives:
                        break
                    if not filter_chembl_activity(
                        act, affinity_cutoff_nm, activity_types
                    ):
                        continue
                    mol_id = act.get("molecule_chembl_id")
                    if mol_id is None:
                        continue
                    mol_json = _get_json(
                        f"{base}/molecule/{mol_id}.json",
                        {},
                        headers,
                        timeout,
                        retries,
                        session=session,
                    )
                    smi = parse_chembl_molecule(mol_json, max_phase)
                    if smi and smi not in smiles:
                        smiles.append(smi)
                        meta["molecules"] += 1
                if meta["molecules"] >= max_actives:
                    break
            if meta["molecules"] >= max_actives:
                break
    except Exception as exc:
        return [], {"error": str(exc), "counts": meta}

    smiles = _limit_smiles(smiles, max_actives)
    save_cached_smiles(cache_dir, source, key, {"smiles": smiles, "meta": meta})
    return smiles, {"cached": False, "counts": meta}


def _parse_bindingdb_entries(entries, affinity_cutoff_nm: int) -> List[str]:
    results = []
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
        if smi not in results:
            results.append(smi)
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
        return cached["smiles"], {"cached": True, "counts": cached.get("meta", {})}

    headers = {"Accept": "application/json"}
    base = "https://bindingdb.org"
    meta = {"pdb_hits": 0, "uniprot_hits": 0}
    smiles: List[str] = []

    def fetch(url, params, counter_key):
        try:
            data = _get_json(url, params, headers, timeout, retries, session=session)
            if isinstance(data, dict) and "results" in data:
                entries = data.get("results") or []
            else:
                entries = data if isinstance(data, list) else []
            parsed = _parse_bindingdb_entries(entries, affinity_cutoff_nm)
            meta[counter_key] += len(parsed)
            for smi in parsed:
                if len(smiles) >= max_actives:
                    break
                if smi not in smiles:
                    smiles.append(smi)
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
        )

    smiles = _limit_smiles(smiles, max_actives)
    save_cached_smiles(cache_dir, source, key, {"smiles": smiles, "meta": meta})
    return smiles, {"cached": False, "counts": meta}


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
        return cached["smiles"], {"cached": True, "counts": cached.get("meta", {})}

    base = "https://www.guidetopharmacology.org/services"
    headers = {"Accept": "application/json"}
    meta = {"targets": 0, "interactions": 0}
    smiles: List[str] = []

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
            if isinstance(tgt_json, list):
                targets.extend(tgt_json)
        except Exception:
            pass
        for ec in ec_numbers:
            if len(targets) >= max_actives:
                break
            try:
                tgt_json = _get_json(
                    f"{base}/targets",
                    {"ecNumber": ec},
                    headers,
                    timeout,
                    retries,
                    session=session,
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
            if smi and smi not in smiles:
                smiles.append(smi)
    smiles = _limit_smiles(smiles, max_actives)
    save_cached_smiles(cache_dir, source, key, {"smiles": smiles, "meta": meta})
    return smiles, {"cached": False, "counts": meta}


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
