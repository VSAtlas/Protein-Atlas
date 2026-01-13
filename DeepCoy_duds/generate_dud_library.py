# ruff: noqa: E402
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from external_sources import (
    SourceAudit,
    _record_audit,
    fetch_bindingdb_smiles,
    fetch_chembl_smiles,
    fetch_drugbank_smiles,
    fetch_iuphar_smiles,
)
from input_and_export_functions import load_config

DEFAULT_INPUT_PDB_DIR = "/home/michael/atlas/code/protein_automation/input_pdbs"
DEFAULT_OUT_ROOT = "../extracted_ligands/deepcoy"
DEFAULT_DEEPCOY_PYTHON = "/home/michael/atlas/anaconda3/envs/DeepCoy-env-cpu/bin/python"
PDB_EXTENSIONS = [".pdb", ".ent", ".cif", ".pdb.gz", ".cif.gz"]
CONFIG_FILE_NAME = "config.txt"
EC_PATTERN = re.compile(r"^(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)$")
DEFAULT_SOURCES = ["chembl", "bindingdb", "iuphar", "rcsb", "rhea"]
SOURCE_ORDER = ["chembl", "bindingdb", "iuphar", "drugbank", "rcsb", "rhea"]
VALID_SOURCES = set(SOURCE_ORDER)
SOURCE_ALIASES = {"pdb_ligands": "rcsb", "pdb": "rcsb"}
DEFAULT_DEEPCOY_ALLOWED_ATOMS = ("C", "N", "O", "S", "F", "Cl", "Br", "I", "P")
DEFAULT_DEEPCOY_ALLOWED_ATOMS_STR = ",".join(DEFAULT_DEEPCOY_ALLOWED_ATOMS)
DEFAULT_DEEPCOY_MODEL = Path("models/DeepCoy_DUDE_model_e09.pickle")
DEFAULT_DEEPCOY_PHOS_MODEL = Path("models/DeepCoy_DUDE_phosphorus_model_e10.pickle")

# --- 1. PDB to EC Conversion (PHASE 1) ---


def _add_ec_candidate(target_list, value):
    if not value:
        return
    cleaned = value.strip()
    if cleaned.startswith("EC:"):
        cleaned = cleaned[3:].strip()
    if cleaned.endswith(";"):
        cleaned = cleaned[:-1]
    if EC_PATTERN.match(cleaned) and cleaned not in target_list:
        target_list.append(cleaned)


def _extract_ec_from_string(value, target_list):
    if not isinstance(value, str):
        return
    for token in re.split(r"[;,\\s]+", value):
        _add_ec_candidate(target_list, token)
    # Also capture inline EC= occurrences
    for match in re.finditer(r"EC=([0-9.\\-]+)", value):
        _add_ec_candidate(target_list, match.group(1))


def extract_ec_numbers_from_polymer_entity_json(entity_data):
    ec_numbers = []
    rcsb_entity = (
        entity_data.get("rcsb_polymer_entity", {})
        if isinstance(entity_data, dict)
        else {}
    )

    combined = rcsb_entity.get("rcsb_enzyme_class_combined", [])
    if isinstance(combined, str):
        combined = [combined]
    for entry in combined or []:
        _extract_ec_from_string(entry, ec_numbers)

    enzyme_class = (
        rcsb_entity.get("rcsb_enzyme_class")
        or rcsb_entity.get("rcsb_enzyme_class_list")
        or []
    )
    if isinstance(enzyme_class, str):
        enzyme_class = [enzyme_class]
    for entry in enzyme_class:
        if isinstance(entry, dict):
            for v in entry.values():
                _extract_ec_from_string(v, ec_numbers)
        else:
            _extract_ec_from_string(entry, ec_numbers)

    def recurse(obj):
        if isinstance(obj, dict):
            for v in obj.values():
                recurse(v)
        elif isinstance(obj, list):
            for v in obj:
                recurse(v)
        elif isinstance(obj, str):
            if "EC" in obj:
                _extract_ec_from_string(obj, ec_numbers)

    recurse(entity_data)
    return ec_numbers


def get_uniprot_and_ec(pdb_id):
    """
    (1) Takes a PDB ID and finds the UniProt ID and EC number.
    """
    pdb_id = pdb_id.upper()
    print(f"## 1. Converting PDB {pdb_id} to UniProt and EC number...")
    uniprot_id = None
    ec_numbers = []
    max_entities = 8

    for entity_id in range(1, max_entities + 1):
        annotation_url = (
            f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
        )
        try:
            print(
                f"-> Querying RCSB Annotation API for {pdb_id} (Entity {entity_id})..."
            )
            response = requests.get(annotation_url)
            if response.status_code != 200:
                continue
            entity_data = response.json()

            identifiers = entity_data.get(
                "rcsb_polymer_entity_container_identifiers", {}
            )
            uniprot_ids = identifiers.get("uniprot_ids", [])

            if uniprot_ids and not uniprot_id:
                uniprot_id = uniprot_ids[0]
                print(f"-> Found UniProt ID: {uniprot_id}")

            ecs_here = extract_ec_numbers_from_polymer_entity_json(entity_data)
            if ecs_here:
                for ec in ecs_here:
                    _add_ec_candidate(ec_numbers, ec)
                print(f"-> Found ECs from entity {entity_id}: {ecs_here}")
        except requests.exceptions.RequestException as req_e:
            print(f"Error querying RCSB Annotation API: {req_e}")
        except Exception as e:
            print(f"Error parsing RCSB Annotation API response: {e}")

        if uniprot_id and ec_numbers:
            break

    if not uniprot_id:
        print("Workflow aborted: Failed to retrieve a UniProt ID.")
        return None, None

    # UniProt EC Fallback/Primary Source
    if not ec_numbers:
        print("-> EC not found in PDB data. Querying UniProt API...")
        uniprot_ec_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"

        try:
            uniprot_response = requests.get(uniprot_ec_url)
            uniprot_response.raise_for_status()

            for line in uniprot_response.text.splitlines():
                if line.startswith("DE") and "EC=" in line:
                    for match in re.finditer(r"EC=([0-9.\\-]+)", line):
                        _add_ec_candidate(ec_numbers, match.group(1))

        except requests.exceptions.RequestException as req_e:
            print(f"Error querying UniProt API: {req_e}")

    final_ec_list = list(ec_numbers)
    print(f"-> Found EC Numbers: {final_ec_list if final_ec_list else 'None'}")

    if not final_ec_list:
        print(
            "Warning: Could not find any EC numbers. Phase 2 (Metabolite collection) may be limited."
        )

    return uniprot_id, final_ec_list


# --- 2. Collecting Actives/Metabolites (PHASE 2) ---


def query_external_sources(
    uniprot_id,
    ec_numbers,
    pdb_id,
    sources,
    cache_dir,
    source_opts,
    max_workers=1,
    run_log_path: Optional[Path] = None,
    audit: Optional["SourceAudit"] = None,
):
    """
    (2) Queries external data sources for ligands/metabolites.
    Returns ([(SMILES, source)], smiles_to_sources, provenance_details).
    """
    print("## 2. Querying external data sources for known metabolites...")
    all_actives = []

    smiles_to_sources: Dict[str, Set[str]] = {}
    provenance_details: Dict[str, List[dict]] = {}

    def add_provenance(smiles, label, detail=None):
        if not smiles or not isinstance(smiles, str):
            return
        label = label or "unknown"
        smiles_to_sources.setdefault(smiles, set()).add(label)
        if detail:
            provenance_details.setdefault(smiles, []).append(detail)

    def add_active(smiles, source):
        if not smiles or not isinstance(smiles, str):
            return
        add_provenance(smiles, source)
        if not any(smiles == active[0] for active in all_actives):
            all_actives.append((smiles, source))

    def merge_provenance_meta(meta: Optional[dict], default_label: str):
        if not isinstance(meta, dict):
            return
        prov_map = meta.get("provenance") or {}
        if isinstance(prov_map, dict):
            for smi, labels in prov_map.items():
                if not isinstance(labels, list):
                    continue
                for lbl in labels:
                    if lbl:
                        add_provenance(smi, lbl, None)
        details_map = meta.get("provenance_details") or {}
        if isinstance(details_map, dict):
            for smi, detail_list in details_map.items():
                if not isinstance(detail_list, list):
                    continue
                for detail in detail_list:
                    if not isinstance(detail, dict):
                        continue
                    label = detail.get("label") or default_label
                    add_provenance(smi, label, detail)

    def _print_meta(src, meta):
        if not meta:
            return
        counts = meta.get("counts", {})
        cached = meta.get("cached", False)
        line = None
        if src == "chembl":
            line = (
                "[chembl] targets={targets} assays={assays} activities={activities} "
                "smiles={smiles} cached={cached}"
            ).format(
                targets=counts.get("targets", 0),
                assays=counts.get("assays", 0),
                activities=counts.get("activities", 0),
                smiles=(
                    len(meta.get("smiles", []))
                    if "smiles" in meta
                    else counts.get("molecules", 0) or len(all_actives)
                ),
                cached=cached,
            )
        elif src == "bindingdb":
            line = (
                "[bindingdb] pdb_hits={pdb_hits} uniprot_hits={uniprot_hits} "
                "smiles={smiles} cached={cached}"
            ).format(
                pdb_hits=counts.get("pdb_hits", 0),
                uniprot_hits=counts.get("uniprot_hits", 0),
                smiles=(
                    len(meta.get("smiles", []))
                    if "smiles" in meta
                    else counts.get("pdb_hits", 0) + counts.get("uniprot_hits", 0)
                ),
                cached=cached,
            )
        elif src == "iuphar":
            line = (
                "[iuphar] targets={targets} interactions={interactions} "
                "smiles={smiles} cached={cached}"
            ).format(
                targets=counts.get("targets", 0),
                interactions=counts.get("interactions", 0),
                smiles=(
                    len(meta.get("smiles", []))
                    if "smiles" in meta
                    else counts.get("interactions", 0)
                ),
                cached=cached,
            )
        elif src == "drugbank":
            line = f"[drugbank] smiles={len(meta.get('smiles', [])) if 'smiles' in meta else 0} info={meta}"
        if line:
            print(line)
            _tee_run_log(run_log_path, line)

    cache_dir.mkdir(parents=True, exist_ok=True)

    source_set = set(sources)
    tasks = {}

    def run_rhea():
        local_smiles = []
        print(
            f"-> Querying Rhea/ChEBI for UniProt {uniprot_id} (Metabolic reactions)..."
        )
        rhea_uniprot_url = (
            f"https://www.rhea-db.org/rest/ws/reaction/uniprot/{uniprot_id}"
        )
        try:
            with requests.Session() as sess:
                rhea_response = sess.get(
                    rhea_uniprot_url, timeout=source_opts["timeout"]
                )
                _record_audit(
                    audit,
                    "rhea",
                    "rhea.uniprot_lookup",
                    rhea_uniprot_url,
                    response_meta={"uniprot": uniprot_id},
                )
                if (
                    rhea_response.status_code == 200
                    and rhea_response.text.strip().startswith("[")
                ):
                    rhea_data = rhea_response.json()
                    for reaction in rhea_data:
                        for side in reaction.get("reactionSides", []):
                            for component in side.get("reactionComponents", []):
                                chebi_id = component.get("chebi")
                                if chebi_id:
                                    chebi_smiles_url = f"https://www.ebi.ac.uk/chebi/rest/api/getCompleteEntity?chebiId={chebi_id}"
                                    chebi_xml_response = sess.get(
                                        chebi_smiles_url, timeout=source_opts["timeout"]
                                    )
                                    _record_audit(
                                        audit,
                                        "rhea",
                                        "rhea.chebi_lookup",
                                        chebi_smiles_url,
                                        response_meta={"chebi_id": chebi_id},
                                    )
                                    if (
                                        chebi_xml_response.status_code == 200
                                        and "<smiles>" in chebi_xml_response.text
                                    ):
                                        start_tag = "<smiles>"
                                        end_tag = "</smiles>"
                                        smiles_start = chebi_xml_response.text.find(
                                            start_tag
                                        ) + len(start_tag)
                                        smiles_end = chebi_xml_response.text.find(
                                            end_tag
                                        )
                                        smiles = chebi_xml_response.text[
                                            smiles_start:smiles_end
                                        ].strip()
                                        if smiles:
                                            label = f"Rhea/ChEBI/{chebi_id}"
                                            local_smiles.append((smiles, label))
                                            add_provenance(
                                                smiles,
                                                label,
                                                {
                                                    "source": "rhea",
                                                    "chebi_id": chebi_id,
                                                    "label": label,
                                                },
                                            )
                else:
                    print(
                        "-> Rhea returned an empty or non-JSON response. Skipping ChEBI extraction."
                    )
        except Exception as e:
            print(f"Error processing Rhea/ChEBI data: {e}")
        return local_smiles, {}

    def run_chembl():
        with requests.Session() as sess:
            smiles, meta = fetch_chembl_smiles(
                uniprot_id,
                ec_numbers,
                pdb_id,
                cache_dir,
                source_opts["timeout"],
                source_opts["retries"],
                source_opts["max_actives_per_source"],
                source_opts["affinity_cutoff_nm"],
                source_opts["chembl_activity_types"],
                source_opts["chembl_max_phase"],
                session=sess,
                audit=audit,
                max_pages=source_opts.get("chembl_max_pages", 5),
            )
        return [(s, "ChEMBL") for s in smiles], meta

    def run_bindingdb():
        with requests.Session() as sess:
            smiles, meta = fetch_bindingdb_smiles(
                uniprot_id,
                ec_numbers,
                pdb_id,
                cache_dir,
                source_opts["timeout"],
                source_opts["retries"],
                source_opts["max_actives_per_source"],
                source_opts["affinity_cutoff_nm"],
                session=sess,
                audit=audit,
            )
        return [(s, "BindingDB") for s in smiles], meta

    def run_iuphar():
        with requests.Session() as sess:
            smiles, meta = fetch_iuphar_smiles(
                uniprot_id,
                ec_numbers,
                pdb_id,
                cache_dir,
                source_opts["timeout"],
                source_opts["retries"],
                source_opts["max_actives_per_source"],
                source_opts["iuphar_approved_only"],
                source_opts["iuphar_primary_only"],
                session=sess,
                audit=audit,
            )
        return [(s, "IUPHAR") for s in smiles], meta

    def run_drugbank():
        smiles, meta = fetch_drugbank_smiles(
            uniprot_id,
            ec_numbers,
            pdb_id,
            cache_dir,
            source_opts.get("drugbank_data_dir"),
        )
        return [(s, "DrugBank") for s in smiles], meta

    def run_rcsb():
        local_smiles = []
        print(f"-> Querying RCSB PDB for ligands in PDB {pdb_id}...")
        ligand_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        try:
            with requests.Session() as sess:
                ligand_response = sess.get(
                    ligand_url, timeout=source_opts["timeout"]
                ).json()
                _record_audit(
                    audit,
                    "rcsb",
                    "rcsb.entry",
                    ligand_url,
                    response_meta={"pdb_id": pdb_id},
                )
                for entity in ligand_response.get("nonpolymer_entity", []):
                    comp_id = entity.get("pdbx_nonpolymer_entity", {}).get(
                        "chem_comp_id"
                    )
                    if comp_id:
                        chem_comp_url = (
                            f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}"
                        )
                        comp_response = sess.get(
                            chem_comp_url, timeout=source_opts["timeout"]
                        ).json()
                        _record_audit(
                            audit,
                            "rcsb",
                            "rcsb.chemcomp",
                            chem_comp_url,
                            response_meta={"chem_comp_id": comp_id},
                        )
                        known_smiles = comp_response.get(
                            "pdbx_chem_comp_descriptor", [{}]
                        )
                        for desc in known_smiles:
                            if desc.get("type") == "SMILES":
                                smi = desc.get("descriptor")
                                if smi:
                                    label = f"PDB/ChemComp/{comp_id}"
                                    local_smiles.append((smi, label))
                                    add_provenance(
                                        smi,
                                        label,
                                        {
                                            "source": "rcsb",
                                            "chem_comp_id": comp_id,
                                            "label": label,
                                        },
                                    )
                                break
        except Exception as e:
            print(f"Error querying RCSB ligands: {e}")
        return local_smiles, {}

    runners = {
        "rhea": run_rhea,
        "chembl": run_chembl,
        "bindingdb": run_bindingdb,
        "iuphar": run_iuphar,
        "drugbank": run_drugbank,
        "rcsb": run_rcsb,
    }

    def submit(executor, name):
        return executor.submit(runners[name])

    if max_workers > 1:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            for src in source_set:
                if src in runners:
                    tasks[src] = submit(executor, src)
            for src in SOURCE_ORDER:
                if src not in tasks:
                    continue
                fut = tasks[src]
                try:
                    results, meta = fut.result()
                except Exception as exc:
                    print(f"[{src}] error={exc}")
                    continue
                _print_meta(src, meta)
                merge_provenance_meta(meta, default_label=src)
                for smi, label in results:
                    add_active(smi, label)
    else:
        for src in SOURCE_ORDER:
            if src not in source_set or src not in runners:
                continue
            try:
                results, meta = runners[src]()
            except Exception as exc:
                print(f"[{src}] error={exc}")
                continue
            _print_meta(src, meta)
            merge_provenance_meta(meta, default_label=src)
            for smi, label in results:
                add_active(smi, label)

    print(f"-> Found {len(all_actives)} raw active candidates from external sources.")
    return all_actives, smiles_to_sources, provenance_details


# --- 3. Filtering and Validation (PHASE 3) ---


def validate_and_filter_actives(
    actives_list,
    smiles_to_sources: Optional[Dict[str, Set[str]]] = None,
    provenance_details=None,
    *,
    potency_cutoff_nm: Optional[float] = None,
    potency_keep_unknown: bool = True,
    run_log_path: Optional[Path] = None,
):
    """
    (3) Merges, deduplicates, and filters the list of active metabolites
        based on RDKit validation and size. Returns filtered smiles, source map,
        and provenance details for kept entries.
    """
    print("## 3. Merging, Deduplicating, and Filtering actives...")

    provided_sources = smiles_to_sources or {}
    combined_sources: Dict[str, Set[str]] = {
        smi: set(labels) for smi, labels in provided_sources.items()
    }
    for smiles, source in actives_list:
        if smiles and isinstance(smiles, str):
            combined_sources.setdefault(smiles, set()).add(source or "unknown")

    ordered_smiles: List[str] = []
    for smiles, _ in actives_list:
        if smiles and smiles not in ordered_smiles:
            ordered_smiles.append(smiles)
    for smiles in combined_sources.keys():
        if smiles not in ordered_smiles:
            ordered_smiles.append(smiles)

    unique_smiles = set(ordered_smiles)
    valid_actives = []
    kept_sources: Dict[str, Set[str]] = {}
    kept_details: Dict[str, List[dict]] = {}
    potency_enabled = potency_cutoff_nm is not None
    potency_kept = 0
    potency_drop_weak = 0
    potency_drop_unknown = 0
    potency_samples = []
    MAX_POTENCY_SAMPLES = 10

    print(f"-> Found {len(unique_smiles)} unique SMILES strings.")

    EXCLUDE_SMILES = {"O", "C", "[Na+]", "[Cl-]", "[Mg+2]", "[K+]"}

    prov_details_map = provenance_details or {}

    for smiles in ordered_smiles:
        if potency_enabled:
            details = (
                prov_details_map.get(smiles, [])
                if isinstance(prov_details_map, dict)
                else []
            )
            best_nm = _extract_best_potency_nm(details)
            if best_nm is None:
                if not potency_keep_unknown:
                    potency_drop_unknown += 1
                    if len(potency_samples) < MAX_POTENCY_SAMPLES:
                        potency_samples.append(
                            {
                                "smiles": smiles,
                                "reason": "unknown_potency",
                                "best_nm": None,
                                "sources": sorted(combined_sources.get(smiles, [])),
                            }
                        )
                    continue
            else:
                if potency_cutoff_nm is not None and best_nm > potency_cutoff_nm:
                    potency_drop_weak += 1
                    if len(potency_samples) < MAX_POTENCY_SAMPLES:
                        potency_samples.append(
                            {
                                "smiles": smiles,
                                "reason": "weak_potency",
                                "best_nm": best_nm,
                                "sources": sorted(combined_sources.get(smiles, [])),
                            }
                        )
                    continue
            potency_kept += 1
        if smiles in EXCLUDE_SMILES:
            continue
        if smiles not in unique_smiles:
            continue
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue

            # Heavy atom count (non-hydrogen atoms)
            heavy_atom_count = mol.GetNumHeavyAtoms()

            # Filter: Tiny < 4 heavy atoms
            if heavy_atom_count < 4:
                continue

            valid_actives.append(smiles)
            labels = combined_sources.get(smiles, {"unknown"})
            label_set = {str(lbl) for lbl in labels if lbl} or {"unknown"}
            kept_sources[smiles] = label_set
            details = (
                prov_details_map.get(smiles)
                if isinstance(prov_details_map, dict)
                else None
            )
            if details:
                kept_details[smiles] = list(details)

        except Exception:
            continue

    if potency_enabled:
        summary_line = (
            "[deepcoy.actives.potency_filter] enabled=true "
            f"cutoff_nm={potency_cutoff_nm} keep_unknown={potency_keep_unknown} "
            f"kept={potency_kept} drop_weak={potency_drop_weak} drop_unknown={potency_drop_unknown}"
        )
        print(summary_line)
        _tee_run_log(run_log_path, summary_line)
        if potency_samples:
            sample_parts = []
            for sample in potency_samples:
                part = (
                    f"smiles={sample.get('smiles')} "
                    f"reason={sample.get('reason')} "
                    f"best_nm={sample.get('best_nm')} "
                    f"sources={','.join(sample.get('sources') or [])}"
                )
                sample_parts.append(part)
            sample_line = "[deepcoy.actives.potency_filter.sample] " + " | ".join(
                sample_parts[:MAX_POTENCY_SAMPLES]
            )
            print(sample_line)
            _tee_run_log(run_log_path, sample_line)

    print(
        f"-> Final list contains {len(valid_actives)} filtered, unique active metabolites."
    )
    return valid_actives, kept_sources, kept_details


# --- 4. DeepCoy Execution (PHASE 4) ---


def write_actives_smiles(
    final_actives_smiles: List[str],
    output_dir,
    label,
    names: Optional[List[str]] = None,
):
    """
    (4) Writes actives to a .smi file for DeepCoy.
    """
    print("## 4. Writing SMILES file...")

    output_dir.mkdir(parents=True, exist_ok=True)
    actives_smi_path = output_dir / f"{label}_actives.smi"

    names = names or [f"METAB_{i + 1}" for i in range(len(final_actives_smiles))]
    if len(names) != len(final_actives_smiles):
        names = [f"METAB_{i + 1}" for i in range(len(final_actives_smiles))]

    # Write the Actives SMILES file (SMILES [tab] Name)
    with actives_smi_path.open("w") as f:
        for smiles, name in zip(final_actives_smiles, names):
            f.write(f"{smiles}\t{name}\n")

    print(f"-> Active SMILES written to: {actives_smi_path}")
    return actives_smi_path, names


def _read_smiles_with_labels(smi_path: Path) -> List[Tuple[str, str]]:
    entries: List[Tuple[str, str]] = []
    if not smi_path.is_file():
        return entries
    with smi_path.open() as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split(None, 1)
            smiles = parts[0]
            label = parts[1] if len(parts) > 1 else ""
            entries.append((label, smiles))
    return entries


def build_provenance_urls(details: List[Dict[str, object]]) -> List[str]:
    urls: List[str] = []
    seen: Set[str] = set()
    for detail in details or []:
        if not isinstance(detail, dict):
            continue
        source = str(detail.get("source") or "").lower()
        if source == "chembl":
            for key, tmpl in [
                (
                    "target_chembl_id",
                    "https://www.ebi.ac.uk/chembl/target_report_card/{id}/",
                ),
                (
                    "assay_chembl_id",
                    "https://www.ebi.ac.uk/chembl/assay_report_card/{id}/",
                ),
                ("activity_chembl_id", "https://www.ebi.ac.uk/chembl/activity/{id}"),
                (
                    "molecule_chembl_id",
                    "https://www.ebi.ac.uk/chembl/compound_report_card/{id}/",
                ),
            ]:
                cid = detail.get(key)
                if cid:
                    url = tmpl.format(id=cid)
                    if url not in seen:
                        urls.append(url)
                        seen.add(url)
        elif source == "iuphar":
            ligand_id = detail.get("ligandId") or detail.get("ligand_id")
            if ligand_id:
                url = f"https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId={ligand_id}"
                if url not in seen:
                    urls.append(url)
                    seen.add(url)
        elif source == "bindingdb":
            lig_id = detail.get("ligand_id")
            if lig_id:
                url = f"https://www.bindingdb.org/bind/chemsearch/marvin/MolStructure.jsp?monomerid={lig_id}"
                if url not in seen:
                    urls.append(url)
                    seen.add(url)
        elif source == "rcsb":
            comp_id = detail.get("chem_comp_id")
            if comp_id:
                url = f"https://www.rcsb.org/chemical/{comp_id}"
                if url not in seen:
                    urls.append(url)
                    seen.add(url)
        elif source == "rhea":
            chebi_id = detail.get("chebi_id")
            if chebi_id:
                url = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={chebi_id}"
                if url not in seen:
                    urls.append(url)
                    seen.add(url)
    return urls


def _normalize_units_to_nm(units: str) -> Optional[float]:
    unit = str(units or "").strip().lower()
    unit = unit.replace("µ", "u").replace("μ", "u")
    if not unit:
        return None
    if unit in {"pm", "picomolar"}:
        return 0.001
    if unit in {"nm", "nanomolar"}:
        return 1.0
    if unit in {"um", "micromolar"}:
        return 1000.0
    if unit in {"mm", "millimolar"}:
        return 1_000_000.0
    return None


def _extract_best_potency_nm(details: List[dict]) -> Optional[float]:
    best_nm: Optional[float] = None
    for detail in details or []:
        if not isinstance(detail, dict):
            continue
        val = detail.get("standard_value")
        if val is None:
            val = detail.get("value")
        try:
            num_val = float(val)
        except Exception:
            continue
        units = detail.get("standard_units") or detail.get("units")
        mult = _normalize_units_to_nm(units)
        if mult is None:
            continue
        relation = (
            str(detail.get("standard_relation") or detail.get("relation") or "=")
            .strip()
            .replace(" ", "")
        )
        if relation in {">", ">="}:
            continue
        if relation not in {"", "=", "<", "<="}:
            continue
        candidate = num_val * mult
        if best_nm is None or candidate < best_nm:
            best_nm = candidate
    return best_nm


def write_provenance_files(
    label: str,
    output_dir: Path,
    labeled_smiles: List[Tuple[str, str]],
    smiles_to_sources: Dict[str, Set[str]],
    provenance_details: Dict[str, List[dict]],
) -> Tuple[Optional[Path], Optional[Path]]:
    if not labeled_smiles:
        return None, None
    output_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = output_dir / f"{label}_actives.provenance.tsv"
    json_path = output_dir / f"{label}_actives.provenance.json"
    records: List[Dict[str, object]] = []
    with tsv_path.open("w") as handle:
        handle.write("metab_id\tsmiles\tsources_joined\n")
        for idx, (metab_id, smiles) in enumerate(labeled_smiles, start=1):
            metab = metab_id or f"METAB_{idx}"
            sources = sorted(smiles_to_sources.get(smiles, [])) or ["unknown"]
            handle.write(f"{metab}\t{smiles}\t{';'.join(sources)}\n")
            details = provenance_details.get(smiles, [])
            records.append(
                {
                    "metab_id": metab,
                    "smiles": smiles,
                    "sources": sources,
                    "browse_urls": build_provenance_urls(details),
                    "details": details,
                }
            )
    json_path.write_text(json.dumps(records, indent=2))
    return tsv_path, json_path


def _parse_allowed_atom_list(raw_value: Optional[str]):
    tokens = []
    seen = set()
    raw_value = (
        raw_value if raw_value is not None else DEFAULT_DEEPCOY_ALLOWED_ATOMS_STR
    )
    for token in re.split(r"[;,\s]+", raw_value):
        cleaned = token.strip()
        if not cleaned or cleaned in seen:
            continue
        tokens.append(cleaned)
        seen.add(cleaned)
    if not tokens:
        tokens = list(DEFAULT_DEEPCOY_ALLOWED_ATOMS)
    return tokens, set(tokens)


def _tee_run_log(run_log_path: Optional[Path], line: str) -> None:
    if not run_log_path:
        return
    try:
        run_log_path.parent.mkdir(parents=True, exist_ok=True)
        with run_log_path.open("a") as handle:
            handle.write(line.rstrip("\n") + "\n")
    except Exception:
        pass


def write_source_audit_artifact(
    audit: Optional["SourceAudit"], dest: Path, run_log_path: Optional[Path] = None
) -> Optional[Path]:
    if not audit or not getattr(audit, "enabled", False):
        return None
    try:
        dest.parent.mkdir(parents=True, exist_ok=True)
        payload = {"requests": audit.as_json()}
        dest.write_text(json.dumps(payload, indent=2))
        line = f"[deepcoy.audit.saved] path={dest} requests={len(payload.get('requests', []))}"
        print(line)
        _tee_run_log(run_log_path, line)
        return dest
    except Exception as exc:  # pragma: no cover - defensive
        warn_line = (
            f"[deepcoy.audit.saved] warning=write_failed path={dest} error={exc}"
        )
        print(warn_line)
        _tee_run_log(run_log_path, warn_line)
        return None


def _parse_sources_value(raw_value: Optional[str]) -> List[str]:
    tokens: List[str] = []
    if raw_value is None:
        return tokens
    seen = set()
    for token in re.split(r"[;,\s]+", str(raw_value)):
        cleaned = token.strip().lower()
        mapped = SOURCE_ALIASES.get(cleaned, cleaned)
        if not mapped or mapped in seen:
            continue
        tokens.append(mapped)
        seen.add(mapped)
    return tokens


def _filter_known_sources(sources: List[str]) -> Tuple[List[str], List[str]]:
    resolved: List[str] = []
    unknown: List[str] = []
    seen = set()
    for src in sources:
        if src not in VALID_SOURCES:
            if src not in unknown:
                unknown.append(src)
            continue
        if src in seen:
            continue
        resolved.append(src)
        seen.add(src)
    return resolved, unknown


def resolve_sources(
    args_sources: Optional[str],
    cfg,
) -> Tuple[List[str], str, str]:
    cli_sources = _parse_sources_value(args_sources)
    try:
        cfg_raw_active = cfg.get("DEEPCOY_ACTIVE_SOURCES")
        cfg_raw_legacy = cfg.get("DEEPCOY_SOURCES")
    except Exception:
        cfg_raw_active = None
        cfg_raw_legacy = None
    cfg_active_sources = (
        _parse_sources_value(cfg_raw_active) if cfg_raw_active is not None else []
    )
    cfg_sources = (
        _parse_sources_value(cfg_raw_legacy) if cfg_raw_legacy is not None else []
    )

    if cli_sources:
        return cli_sources, "cli", args_sources if args_sources is not None else ""
    if cfg_active_sources:
        return cfg_active_sources, "config_active", str(cfg_raw_active)
    if cfg_sources:
        return cfg_sources, "config", str(cfg_raw_legacy)
    return list(DEFAULT_SOURCES), "default", ",".join(DEFAULT_SOURCES)


def build_source_audit_lines(
    uniprot_id: str,
    ec_numbers: List[str],
    pdb_id: str,
    source_opts,
    sources: List[str],
) -> List[str]:
    chembl_types_raw = source_opts.get("chembl_activity_types") if source_opts else []
    if isinstance(chembl_types_raw, str):
        chembl_types = [t for t in re.split(r"[;,\s]+", chembl_types_raw) if t.strip()]
    else:
        chembl_types = list(chembl_types_raw or [])
    chembl_type_param = ",".join(chembl_types)
    cutoff = source_opts.get("affinity_cutoff_nm") if source_opts else None
    cutoff_str = str(cutoff) if cutoff is not None else ""
    approved_only = source_opts.get("iuphar_approved_only") if source_opts else None
    primary_only = source_opts.get("iuphar_primary_only") if source_opts else None
    approved_str = str(bool(approved_only)).lower()
    primary_str = str(bool(primary_only)).lower()
    data_dir = source_opts.get("drugbank_data_dir") if source_opts else None

    lines: List[str] = []
    for src in sources:
        if src == "chembl":
            base = "https://www.ebi.ac.uk/chembl/api/data"
            lines.append(
                f"[deepcoy.source_audit] source=chembl url={base}/target?target_components__accession={uniprot_id}&format=json"
            )
            lines.append(
                "[deepcoy.source_audit] source=chembl url_template=https://www.ebi.ac.uk/chembl/api/data/assay?target_chembl_id=<TARGET_CHEMBL_ID>&assay_type=B&relationship_type=D&format=json"
            )
            lines.append(
                f"[deepcoy.source_audit] source=chembl url_template=https://www.ebi.ac.uk/chembl/api/data/activity?assay_chembl_id=<ASSAY_CHEMBL_ID>&standard_type__in={chembl_type_param}&format=json"
            )
            lines.append(
                "[deepcoy.source_audit] source=chembl url_template=https://www.ebi.ac.uk/chembl/api/data/molecule/<MOLECULE_CHEMBL_ID>.json"
            )
        elif src == "bindingdb":
            base = "https://bindingdb.org"
            lines.append(
                f"[deepcoy.source_audit] source=bindingdb url={base}/rest/getLigandsByPDBs?pdb={pdb_id}&cutoff={cutoff_str}&identity=100&response=application/json"
            )
            lines.append(
                f"[deepcoy.source_audit] source=bindingdb url={base}/rest/getLigandsByUniprots?uniprot={uniprot_id}&cutoff={cutoff_str}&response=application/json"
            )
        elif src == "iuphar":
            base = "https://www.guidetopharmacology.org/services"
            lines.append(
                f"[deepcoy.source_audit] source=iuphar url={base}/targets?accession={uniprot_id}"
            )
            for ec in ec_numbers or []:
                lines.append(
                    f"[deepcoy.source_audit] source=iuphar url={base}/targets?ecNumber={ec}"
                )
            lines.append(
                f"[deepcoy.source_audit] source=iuphar url_template={base}/targets/<TARGET_ID>/interactions?approved={approved_str}&primaryTarget={primary_str}"
            )
            lines.append(
                f"[deepcoy.source_audit] source=iuphar url_template={base}/ligands/<LIGAND_ID>/structure"
            )
        elif src == "rcsb":
            base = "https://data.rcsb.org/rest/v1/core"
            lines.append(
                f"[deepcoy.source_audit] source=rcsb url={base}/entry/{pdb_id}"
            )
            lines.append(
                f"[deepcoy.source_audit] source=rcsb url_template={base}/chemcomp/<CHEM_COMP_ID>"
            )
        elif src == "rhea":
            lines.append(
                f"[deepcoy.source_audit] source=rhea url=https://www.rhea-db.org/rest/ws/reaction/uniprot/{uniprot_id}"
            )
            lines.append(
                "[deepcoy.source_audit] source=rhea url_template=https://www.ebi.ac.uk/chebi/rest/api/getCompleteEntity?chebiId=<CHEBI_ID>"
            )
        elif src == "drugbank":
            lines.append(
                f"[deepcoy.source_audit] source=drugbank local_data_dir={data_dir or 'none'}"
            )
    return lines


def filter_actives_for_deepcoy(
    actives_smi_path: Path,
    allowed_atoms: List[str],
    *,
    filter_enabled: bool = True,
    run_log_path: Optional[Path] = None,
):
    allowed_atoms = allowed_atoms or list(DEFAULT_DEEPCOY_ALLOWED_ATOMS)
    allowed_atoms_set = set(allowed_atoms)
    allowed_atoms_str = ",".join(allowed_atoms)
    rejected_path = actives_smi_path.with_name(
        f"{actives_smi_path.stem}.rejected_by_deepcoy{actives_smi_path.suffix}"
    )
    raw_copy_path = actives_smi_path.with_name(
        f"{actives_smi_path.stem}.raw{actives_smi_path.suffix}"
    )

    if not filter_enabled:
        total = 0
        with actives_smi_path.open("r") as handle:
            for raw_line in handle:
                stripped = raw_line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                total += 1
        log_line = (
            f"[deepcoy.filter] disabled mode; input={total} kept={total} rejected=0 "
            f"allowed_atoms={allowed_atoms_str}"
        )
        print(log_line)
        _tee_run_log(run_log_path, log_line)
        return total, total, 0

    shutil.copyfile(actives_smi_path, raw_copy_path)

    total = 0
    kept = 0
    kept_lines = []
    rejected_entries = []

    with raw_copy_path.open("r") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            total += 1
            parts = stripped.split(None, 1)
            smiles = parts[0]
            label = parts[1] if len(parts) > 1 else ""
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                rejected_entries.append((smiles, label, "parse_fail"))
                continue
            atom_symbols = {atom.GetSymbol() for atom in mol.GetAtoms()}
            unsupported = sorted(
                sym for sym in atom_symbols if sym not in allowed_atoms_set
            )
            if unsupported:
                rejected_entries.append(
                    (smiles, label, f"unsupported_atoms:{','.join(unsupported)}")
                )
                continue
            kept += 1
            if label:
                kept_lines.append(f"{smiles}\t{label}\n")
            else:
                kept_lines.append(f"{smiles}\n")

    rejected_path.parent.mkdir(parents=True, exist_ok=True)
    if rejected_entries:
        with rejected_path.open("w") as handle:
            for smiles, label, reason in rejected_entries:
                handle.write(f"{smiles}\t{label}\t{reason}\n")
    else:
        if rejected_path.exists():
            rejected_path.unlink()

    log_line = (
        f"[deepcoy.filter] input={total} kept={kept} rejected={len(rejected_entries)} "
        f"allowed_atoms={allowed_atoms_str} rejected_file={rejected_path}"
    )
    print(log_line)
    _tee_run_log(run_log_path, log_line)

    if not kept_lines:
        raise ValueError(
            f"All actives were rejected by DeepCoy compatibility filter. See {rejected_path} for details."
        )

    actives_smi_path.write_text("".join(kept_lines))
    return total, kept, len(rejected_entries)


def copy_filter_debug_artifacts(
    actives_smi_path: Path,
    artifact_dir: Path,
    run_tag: Optional[str],
    run_log_path: Optional[Path] = None,
) -> None:
    if not run_tag:
        return
    try:
        debug_dir = artifact_dir / "filter_debug" / run_tag
        debug_dir.mkdir(parents=True, exist_ok=True)
        candidates = [
            actives_smi_path,
            actives_smi_path.with_name(
                f"{actives_smi_path.stem}.raw{actives_smi_path.suffix}"
            ),
            actives_smi_path.with_name(
                f"{actives_smi_path.stem}.rejected_by_deepcoy{actives_smi_path.suffix}"
            ),
        ]
        copied = []
        for path in candidates:
            if path.exists():
                shutil.copyfile(path, debug_dir / path.name)
                copied.append(path.name)
        log_line = f"[deepcoy.filter.debug] run_tag={run_tag} dir={debug_dir} copied={','.join(copied) if copied else 'none'}"
        print(log_line)
        _tee_run_log(run_log_path, log_line)
    except Exception:
        pass


def _chunk_contains_atom(actives_lines: List[str], atom_symbol: str) -> bool:
    for raw_line in actives_lines:
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split(None, 1)
        smiles = parts[0]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == atom_symbol:
                return True
    return False


def convert_smi_to_sdf(
    smi_path: Path, sdf_path: Path, *, require_output: bool = False
) -> int:
    if not smi_path.is_file():
        print(f"WARNING: SMILES file not found for SDF conversion: {smi_path}")
        if require_output:
            raise FileNotFoundError(f"Missing SMILES for SDF conversion: {smi_path}")
        return 0
    suppl = []
    with smi_path.open() as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split()
            smi = parts[0]
            name = parts[1] if len(parts) > 1 else "DECOY"
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            try:
                AllChem.Compute2DCoords(mol)
            except Exception:
                pass
            mol.SetProp("_Name", name)
            suppl.append(mol)
    if not suppl:
        print(f"WARNING: No valid molecules to convert for {smi_path}")
        if require_output:
            raise ValueError(f"No valid molecules to convert for {smi_path}")
        return 0
    sdf_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(sdf_path))
    count = 0
    for mol in suppl:
        writer.write(mol)
        count += 1
    writer.close()
    print(
        f"-> Wrote {count} decoys to SDF: {sdf_path} ({sdf_path.stat().st_size} bytes)"
    )
    return count


def resolve_path(path_str, base_dir):
    path = Path(path_str).expanduser()
    if path.is_absolute():
        return path
    return base_dir / path


def find_pdb_file(input_dir, pdb_id):
    if not input_dir.is_dir():
        return None

    name_map = {}
    for entry in input_dir.iterdir():
        if entry.is_file():
            name_map[entry.name.lower()] = entry.name

    pdb_lower = pdb_id.lower()
    for ext in PDB_EXTENSIONS:
        candidate = f"{pdb_lower}{ext}"
        if candidate in name_map:
            return input_dir / name_map[candidate]
    return None


class DeepcoyChunkPlan:
    def __init__(
        self,
        index: int,
        actives_path: Path,
        output_dir: Path,
        artifact_dir: Path,
        line_count: int,
        is_chunk: bool = True,
        model_path: Optional[Path] = None,
        seed: Optional[int] = None,
        dataset: Optional[str] = None,
    ):
        self.index = index
        self.actives_path = actives_path
        self.output_dir = output_dir
        self.artifact_dir = artifact_dir
        self.line_count = line_count
        self.is_chunk = is_chunk
        self.model_path = model_path
        self.seed = seed
        self.dataset = dataset


def build_deepcoy_env(thread_count: int) -> dict:
    env = {}
    if thread_count and thread_count > 0:
        thread_val = str(thread_count)
        for var in [
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS",
        ]:
            env[var] = thread_val
        env.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
        env["DEEPCOY_TF_INTRA_THREADS"] = thread_val
        env["DEEPCOY_TF_INTER_THREADS"] = thread_val
    return env


def _sanitize_run_tag(tag: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(tag))


def build_deepcoy_run_tag(run_tag: Optional[str] = None) -> str:
    if run_tag:
        return _sanitize_run_tag(run_tag)
    parts = []
    atlas_run_id = os.environ.get("ATLAS_RUN_ID")
    if atlas_run_id:
        parts.append(_sanitize_run_tag(atlas_run_id))
    parts.append(time.strftime("%Y%m%d_%H%M%S"))
    parts.append(str(os.getpid()))
    return "_".join([p for p in parts if p])


def resolve_deepcoy_chunk_roots(
    output_dir: Path,
    artifact_dir: Path,
    actives_smi_path: Path,
    keep_chunks: bool,
    run_tag: Optional[str],
) -> Tuple[Path, Path, str]:
    label_name = actives_smi_path.stem
    if label_name.endswith("_actives"):
        label_name = label_name[: -len("_actives")]

    chunk_output_root = output_dir.parent / ".deepcoy_chunks" / label_name
    chunk_artifact_root = artifact_dir / "chunks"

    if keep_chunks:
        tag = build_deepcoy_run_tag(run_tag)
        chunk_output_root = chunk_output_root / tag
        chunk_artifact_root = chunk_artifact_root / tag

    return chunk_output_root, chunk_artifact_root, label_name


def prepare_deepcoy_chunks(
    actives_smi_path: Path,
    output_dir: Path,
    artifact_dir: Path,
    chunk_output_root: Path,
    chunk_artifact_root: Path,
    chunk_size: int,
    use_chunking: bool,
):
    if chunk_size <= 0:
        raise ValueError("deepcoy_chunk_size must be positive.")

    with actives_smi_path.open("r") as f:
        actives_lines = f.readlines()

    if not use_chunking:
        return actives_lines, [
            DeepcoyChunkPlan(
                index=0,
                actives_path=actives_smi_path,
                output_dir=output_dir,
                artifact_dir=artifact_dir,
                line_count=len(actives_lines),
                is_chunk=False,
            )
        ]

    shutil.rmtree(chunk_output_root, ignore_errors=True)
    shutil.rmtree(chunk_artifact_root, ignore_errors=True)
    chunk_output_root.mkdir(parents=True, exist_ok=True)
    chunk_artifact_root.mkdir(parents=True, exist_ok=True)

    chunk_plans = []
    for start in range(0, len(actives_lines), chunk_size):
        chunk_index = start // chunk_size
        chunk_name = f"chunk_{chunk_index:03d}"
        artifact_dir_chunk = chunk_artifact_root / chunk_name
        output_dir_chunk = chunk_output_root / chunk_name
        artifact_dir_chunk.mkdir(parents=True, exist_ok=True)
        output_dir_chunk.mkdir(parents=True, exist_ok=True)

        chunk_actives_path = artifact_dir_chunk / "actives.smi"
        chunk_lines = actives_lines[start : start + chunk_size]
        with chunk_actives_path.open("w") as f:
            f.writelines(chunk_lines)

        chunk_plans.append(
            DeepcoyChunkPlan(
                index=chunk_index,
                actives_path=chunk_actives_path,
                output_dir=output_dir_chunk,
                artifact_dir=artifact_dir_chunk,
                line_count=len(chunk_lines),
            )
        )

    return actives_lines, chunk_plans


def _assign_chunk_metadata(
    chunk_plans,
    phosphorus_model_path: Path,
    default_model_path: Path,
    base_seed: Optional[int],
    seed_per_chunk: bool,
    run_log_path: Optional[Path],
):
    for plan in chunk_plans:
        try:
            with plan.actives_path.open("r") as handle:
                lines = handle.readlines()
        except Exception:
            lines = []
        requires_p = _chunk_contains_atom(lines, "P")
        model = phosphorus_model_path if requires_p else default_model_path
        if not model.is_file():
            raise FileNotFoundError(f"DeepCoy model not found: {model}")
        plan.model_path = model
        dataset = "zinc_phosphorus" if requires_p else "zinc"
        plan.dataset = dataset
        reason = "contains_P" if requires_p else "default"
        log_line = f"[deepcoy.model] chunk={plan.index:03d} model={model} reason={reason} dataset={dataset}"
        print(log_line)
        _tee_run_log(run_log_path, log_line)
        if base_seed is not None:
            plan.seed = base_seed + plan.index if seed_per_chunk else base_seed
            mode = "base+chunk_index" if seed_per_chunk else "base_only"
            seed_line = f"[deepcoy.seed] chunk={plan.index:03d} seed={plan.seed} base={base_seed} mode={mode}"
            print(seed_line)
            _tee_run_log(run_log_path, seed_line)


def assign_chunk_metadata(
    chunk_plans,
    phosphorus_model_path: Path,
    default_model_path: Path,
    base_seed: Optional[int],
    seed_per_chunk: bool,
    run_log_path: Optional[Path],
):
    _assign_chunk_metadata(
        chunk_plans,
        phosphorus_model_path,
        default_model_path,
        base_seed,
        seed_per_chunk,
        run_log_path,
    )


def merge_deepcoy_outputs(
    chunk_plans,
    merged_smi_path: Path,
    decoys_per_active: int,
    total_actives: int,
    keep_chunks: bool,
    chunk_output_root: Optional[Path] = None,
    chunk_artifact_root: Optional[Path] = None,
    run_log_path: Optional[Path] = None,
):
    merged_lines = []
    for plan in sorted(chunk_plans, key=lambda p: p.index):
        smi_path = plan.output_dir / "deepcoy_decoys.smi"
        if not smi_path.is_file():
            raise FileNotFoundError(f"DeepCoy chunk output missing: {smi_path}")
        with smi_path.open("r") as handle:
            for raw in handle:
                stripped = raw.strip()
                if not stripped:
                    continue
                merged_lines.append(raw.rstrip("\n"))

    merged_smi_path.parent.mkdir(parents=True, exist_ok=True)
    merged_smi_path.write_text("\n".join(merged_lines) + ("\n" if merged_lines else ""))

    produced = len(merged_lines)
    expected = (
        total_actives * decoys_per_active if decoys_per_active is not None else None
    )
    if expected is not None and produced < expected:
        warning_line = f"[deepcoy] WARNING: DeepCoy generated {produced} decoys (expected ~{expected})."
        print(warning_line)
        _tee_run_log(run_log_path, warning_line)

    if not keep_chunks:
        if chunk_output_root:
            shutil.rmtree(chunk_output_root, ignore_errors=True)
        if chunk_artifact_root:
            shutil.rmtree(chunk_artifact_root, ignore_errors=True)
    return produced


def run_deepcoy_workflow(
    actives_smi_path,
    output_dir,
    artifact_dir,
    deepcoy_run_sh,
    deepcoy_sdf_sh,
    decoys_per_active,
    deepcoy_python,
    config_path,
    skip_sdf,
    ensure_sdf,
    restrict_data,
    run_log_path: Optional[Path] = None,
    deepcoy_workers=1,
    deepcoy_chunk_size=1,
    deepcoy_threads=1,
    deepcoy_keep_chunks=False,
    default_model_path: Optional[Path] = None,
    phosphorus_model_path: Optional[Path] = None,
    base_seed: Optional[int] = None,
    seed_per_chunk: bool = False,
    use_argmax_generation: Optional[bool] = None,
    try_different_starting: Optional[bool] = None,
    num_different_starting: Optional[int] = None,
    num_samples: Optional[int] = None,
    run_tag: Optional[str] = None,
):
    start_line = "## 5. Executing DeepCoy workflow..."
    print(start_line)
    _tee_run_log(run_log_path, start_line)

    repo_root = Path(__file__).resolve().parent
    default_model_path = default_model_path or repo_root / DEFAULT_DEEPCOY_MODEL
    phosphorus_model_path = (
        phosphorus_model_path or repo_root / DEFAULT_DEEPCOY_PHOS_MODEL
    )
    env_overrides = build_deepcoy_env(deepcoy_threads)
    _ = deepcoy_sdf_sh

    keep_chunks_flag = bool(deepcoy_keep_chunks)
    use_chunking = deepcoy_workers > 1
    run_tag_resolved = (
        build_deepcoy_run_tag(run_tag)
        if keep_chunks_flag
        else (run_tag if run_tag else None)
    )
    chunk_output_root, chunk_artifact_root, label_name = resolve_deepcoy_chunk_roots(
        output_dir, artifact_dir, actives_smi_path, keep_chunks_flag, run_tag_resolved
    )
    run_log_line = (
        f"[deepcoy.run] run_tag={run_tag_resolved or '-'} "
        f"chunk_out={chunk_output_root} chunk_art={chunk_artifact_root} keep_chunks={keep_chunks_flag}"
    )
    print(run_log_line)
    _tee_run_log(run_log_path, run_log_line)
    actives_lines, chunk_plans = prepare_deepcoy_chunks(
        actives_smi_path,
        output_dir,
        artifact_dir,
        chunk_output_root,
        chunk_artifact_root,
        deepcoy_chunk_size,
        use_chunking=use_chunking,
    )
    print(
        f"[deepcoy.chunking] chunk_size={deepcoy_chunk_size} use_chunking={use_chunking} "
        f"chunks={len(chunk_plans)} total_actives={len(actives_lines)}"
    )
    _tee_run_log(
        run_log_path,
        f"[deepcoy.chunking] chunk_size={deepcoy_chunk_size} use_chunking={use_chunking} "
        f"chunks={len(chunk_plans)} total_actives={len(actives_lines)}",
    )
    _assign_chunk_metadata(
        chunk_plans,
        phosphorus_model_path,
        default_model_path,
        base_seed,
        seed_per_chunk,
        run_log_path,
    )

    if not use_chunking:
        plan_single = chunk_plans[0] if chunk_plans else None
        deepcoy_cmd = [
            str(deepcoy_run_sh),
            str(actives_smi_path),
            str(output_dir),
            "--artifact-dir",
            str(artifact_dir),
            "--decoys-per-active",
            str(decoys_per_active),
            "--config",
            str(config_path),
            "--deepcoy-python",
            deepcoy_python,
        ]
        dataset = plan_single.dataset if plan_single and plan_single.dataset else "zinc"
        deepcoy_cmd.extend(["--dataset", dataset])
        if restrict_data and int(restrict_data) > 0:
            deepcoy_cmd.extend(["--restrict-data", str(restrict_data)])
        if plan_single and plan_single.model_path:
            deepcoy_cmd.extend(["--restore-model", str(plan_single.model_path)])
        if plan_single and plan_single.seed is not None:
            deepcoy_cmd.extend(["--random-seed", str(plan_single.seed)])
        if use_argmax_generation is not None:
            deepcoy_cmd.extend(["--use-argmax-generation", str(use_argmax_generation)])
        if try_different_starting is not None:
            deepcoy_cmd.extend(
                ["--try-different-starting", str(try_different_starting)]
            )
        if num_different_starting is not None:
            deepcoy_cmd.extend(
                ["--num-different-starting", str(num_different_starting)]
            )
        if num_samples is not None:
            deepcoy_cmd.extend(["--num-samples", str(num_samples)])
        print(f"Running: {' '.join(deepcoy_cmd)}")
        env = os.environ.copy()
        env.update(env_overrides)
        subprocess.run(deepcoy_cmd, check=True, env=env)

        smi_path = output_dir / "deepcoy_decoys.smi"
        sdf_path = output_dir / "deepcoy_decoys.sdf"
        if ensure_sdf or not skip_sdf:
            try:
                convert_smi_to_sdf(smi_path, sdf_path, require_output=True)
            except Exception as exc:
                raise RuntimeError(f"Decoy SDF conversion failed: {exc}") from exc
        else:
            print("-> Skipping SDF conversion (--skip-sdf).")
        return

    def _caches_ready():
        cache_dir = repo_root / "data" / "cache"
        train_cache = cache_dir / "molecules_train_zinc.preprocessed.pkl"
        valid_cache = cache_dir / "molecules_valid_zinc.preprocessed.pkl"
        return train_cache.exists() and valid_cache.exists()

    def _run_chunk(plan: DeepcoyChunkPlan):
        env = os.environ.copy()
        env.update(env_overrides)
        deepcoy_cmd_local = [
            str(deepcoy_run_sh),
            str(plan.actives_path),
            str(plan.output_dir),
            "--artifact-dir",
            str(plan.artifact_dir),
            "--decoys-per-active",
            str(decoys_per_active),
            "--config",
            str(config_path),
            "--deepcoy-python",
            deepcoy_python,
        ]
        dataset = plan.dataset or "zinc"
        deepcoy_cmd_local.extend(["--dataset", dataset])
        if restrict_data and int(restrict_data) > 0:
            deepcoy_cmd_local.extend(["--restrict-data", str(restrict_data)])
        if plan.model_path:
            deepcoy_cmd_local.extend(["--restore-model", str(plan.model_path)])
        if plan.seed is not None:
            deepcoy_cmd_local.extend(["--random-seed", str(plan.seed)])
        if use_argmax_generation is not None:
            deepcoy_cmd_local.extend(
                ["--use-argmax-generation", str(use_argmax_generation)]
            )
        if try_different_starting is not None:
            deepcoy_cmd_local.extend(
                ["--try-different-starting", str(try_different_starting)]
            )
        if num_different_starting is not None:
            deepcoy_cmd_local.extend(
                ["--num-different-starting", str(num_different_starting)]
            )
        if num_samples is not None:
            deepcoy_cmd_local.extend(["--num-samples", str(num_samples)])
        print(
            f"[deepcoy] Starting chunk {plan.index:03d} ({plan.line_count} actives) -> {plan.output_dir}"
        )
        run_kwargs = {
            "check": True,
            "stdout": subprocess.PIPE,
            "stderr": subprocess.PIPE,
            "universal_newlines": True,
            "env": env,
        }
        try:
            subprocess.run(
                deepcoy_cmd_local,
                **run_kwargs,
            )
        except subprocess.CalledProcessError as exc:
            print(
                f"[deepcoy] Chunk {plan.index:03d} failed with exit code {exc.returncode}",
                file=sys.stderr,
            )
            if exc.stdout:
                print(exc.stdout)
            if exc.stderr:
                print(exc.stderr, file=sys.stderr)
            if exc.returncode == 4:
                skip_line = f"[deepcoy] Chunk {plan.index:03d} rejected by preprocessing (rc=4); skipping chunk."
                print(skip_line)
                _tee_run_log(run_log_path, skip_line)
                return None
            raise
        print(f"[deepcoy] Finished chunk {plan.index:03d}")
        return plan.output_dir / "deepcoy_decoys.smi"

    plans_to_run = list(chunk_plans)
    successful_plans = []
    skipped_plans = []
    if not _caches_ready() and plans_to_run:
        print(
            "[deepcoy] Priming DeepCoy cache with first chunk before parallel execution..."
        )
        while plans_to_run:
            plan = plans_to_run.pop(0)
            result = _run_chunk(plan)
            if result is None:
                skipped_plans.append(plan)
                continue
            successful_plans.append(plan)
            break

    if plans_to_run:
        concurrency = min(deepcoy_workers, len(plans_to_run))
        with ThreadPoolExecutor(max_workers=concurrency) as executor:
            future_map = {
                executor.submit(_run_chunk, plan): plan for plan in plans_to_run
            }
            for future in as_completed(future_map):
                plan = future_map[future]
                try:
                    result = future.result()
                except Exception:
                    for f in future_map:
                        f.cancel()
                    raise
                if result is None:
                    skipped_plans.append(plan)
                else:
                    successful_plans.append(plan)

    skipped_actives = sum(plan.line_count for plan in skipped_plans)
    successful_actives = sum(plan.line_count for plan in successful_plans)
    if skipped_plans:
        warning_line = (
            f"[deepcoy] WARNING: skipped {len(skipped_plans)} chunk(s) "
            f"({skipped_actives} actives) due to preprocessing rc=4."
        )
        print(warning_line)
        _tee_run_log(run_log_path, warning_line)

    if chunk_plans and not successful_plans:
        raise subprocess.CalledProcessError(
            4,
            [str(deepcoy_run_sh), "--chunked-run"],
        )

    final_smi_path = output_dir / "deepcoy_decoys.smi"
    merge_deepcoy_outputs(
        successful_plans,
        final_smi_path,
        decoys_per_active,
        total_actives=successful_actives,
        keep_chunks=keep_chunks_flag,
        chunk_output_root=chunk_output_root,
        chunk_artifact_root=chunk_artifact_root,
        run_log_path=run_log_path,
    )

    if ensure_sdf or not skip_sdf:
        try:
            convert_smi_to_sdf(
                final_smi_path, output_dir / "deepcoy_decoys.sdf", require_output=True
            )
        except Exception as exc:
            raise RuntimeError(f"Decoy SDF conversion failed: {exc}") from exc
    else:
        print("-> Skipping SDF conversion (--skip-sdf).")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a DeepCoy DUD-style decoy library from a target PDB."
    )

    parser.add_argument("--pdb", required=True, help="PDB ID to query (required).")
    parser.add_argument(
        "--input-pdb-dir",
        default=DEFAULT_INPUT_PDB_DIR,
        help="Directory containing input PDBs.",
    )
    parser.add_argument(
        "--out-root", default=DEFAULT_OUT_ROOT, help="Output root for actives/decoys."
    )
    parser.add_argument(
        "--artifact-dir",
        default=None,
        help="Directory for DeepCoy run artifacts (defaults to deepcoy_work/<PDB>).",
    )
    parser.add_argument(
        "--sources",
        default=None,
        help=(
            "Comma-separated sources to query. If omitted, uses config.txt DEEPCOY_ACTIVE_SOURCES (preferred) or DEEPCOY_SOURCES, "
            f"else falls back to built-in default {','.join(DEFAULT_SOURCES)}."
        ),
    )
    parser.add_argument(
        "--http-timeout",
        type=int,
        default=20,
        help="HTTP timeout (seconds) for external source queries.",
    )
    parser.add_argument(
        "--http-retries",
        type=int,
        default=2,
        help="HTTP retries for external source queries.",
    )
    parser.add_argument(
        "--max-actives-per-source",
        type=int,
        default=500,
        help="Maximum SMILES to retain per source.",
    )
    parser.add_argument(
        "--affinity-cutoff-nm",
        type=int,
        default=10000,
        help="Affinity cutoff in nM for potency-filtered sources (e.g., BindingDB, ChEMBL).",
    )
    parser.add_argument(
        "--iuphar-approved-only",
        dest="iuphar_approved_only",
        action="store_true",
        help="Use only approved IUPHAR interactions (default).",
    )
    parser.add_argument(
        "--iuphar-include-nonapproved",
        dest="iuphar_approved_only",
        action="store_false",
        help="Allow non-approved IUPHAR interactions.",
    )
    parser.add_argument(
        "--iuphar-primary-target-only",
        dest="iuphar_primary_only",
        action="store_true",
        help="Use only primary targets in IUPHAR (default).",
    )
    parser.add_argument(
        "--iuphar-include-secondary-targets",
        dest="iuphar_primary_only",
        action="store_false",
        help="Allow secondary targets in IUPHAR.",
    )
    parser.add_argument(
        "--chembl-activity-types",
        default="Ki,Kd,IC50",
        help="Comma-separated ChEMBL activity types (default: Ki,Kd,IC50).",
    )
    parser.add_argument(
        "--chembl-max-phase",
        type=int,
        default=None,
        help="Optional ChEMBL max_phase upper bound to bias toward approved drugs.",
    )
    parser.add_argument(
        "--chembl-max-pages",
        type=int,
        default=5,
        help="Maximum number of ChEMBL pages to follow for paginated endpoints (default: 5).",
    )
    parser.add_argument(
        "--max-source-workers",
        type=int,
        default=1,
        help="Max threads for parallel source queries (default 1 = sequential).",
    )
    parser.add_argument(
        "--drugbank-data-dir",
        default=None,
        help="Path to local DrugBank structured data. If unset, DrugBank is skipped.",
    )

    run_group = parser.add_mutually_exclusive_group()
    run_group.add_argument(
        "--run-deepcoy",
        dest="run_deepcoy",
        action="store_true",
        help="Run DeepCoy after preparing actives.",
    )
    run_group.add_argument(
        "--no-run-deepcoy",
        dest="run_deepcoy",
        action="store_false",
        help="Skip DeepCoy run.",
    )
    run_group.add_argument(
        "--skip-deepcoy",
        dest="run_deepcoy",
        action="store_false",
        help="Alias for --no-run-deepcoy.",
    )
    parser.set_defaults(run_deepcoy=True)
    parser.set_defaults(
        iuphar_approved_only=True, iuphar_primary_only=True, deepcoy_filter_actives=True
    )

    parser.add_argument(
        "--deepcoy-python",
        default=None,
        help=(
            "Path to DeepCoy Python interpreter. Defaults to $DEEPCOY_PYTHON, "
            "config.txt DEEPCOY_PYTHON, or the built-in DeepCoy path."
        ),
    )
    parser.add_argument(
        "--deepcoy-run-sh", default="./deepcoy_run.sh", help="Path to deepcoy_run.sh."
    )
    parser.add_argument(
        "--deepcoy-sdf-sh",
        default="./deepcoy_smiles_to_sdf.sh",
        help="Path to deepcoy_smiles_to_sdf.sh.",
    )
    parser.add_argument(
        "--ensure-sdf",
        action="store_true",
        help="Force SDF generation with RDKit even if --skip-sdf is set.",
    )
    parser.add_argument(
        "--restrict-data",
        type=int,
        default=0,
        help="Restrict DeepCoy data size for faster smoke tests (passed to DeepCoy).",
    )

    parser.add_argument(
        "--skip-sdf",
        action="store_true",
        help="Skip SDF conversion step after DeepCoy.",
    )
    parser.add_argument(
        "--decoys-per-active",
        type=int,
        default=None,
        help="Number of decoys to generate per active (default 50, can be set via config.txt DEEPCOY_DECOYS_PER_ACTIVE).",
    )
    parser.add_argument(
        "--deepcoy-workers",
        type=int,
        default=10,
        help="Parallel DeepCoy workers (default 10).",
    )
    parser.add_argument(
        "--deepcoy-chunk-size",
        type=int,
        default=None,
        help="Number of actives per DeepCoy chunk (default 1).",
    )
    parser.add_argument(
        "--deepcoy-threads",
        type=int,
        default=1,
        help="Threads per DeepCoy worker for BLAS/TF (default 1).",
    )
    parser.add_argument(
        "--deepcoy-keep-chunks",
        dest="deepcoy_keep_chunks",
        action="store_true",
        help="Keep per-chunk DeepCoy outputs instead of cleaning after merge.",
    )
    parser.add_argument(
        "--deepcoy-no-keep-chunks",
        dest="deepcoy_keep_chunks",
        action="store_false",
        help="Do not keep per-chunk DeepCoy outputs after merge.",
    )
    parser.set_defaults(deepcoy_keep_chunks=None)
    parser.add_argument(
        "--deepcoy-base-seed",
        type=int,
        default=None,
        help="Base random seed for DeepCoy (overrides config).",
    )
    parser.add_argument(
        "--deepcoy-seed-per-chunk",
        action="store_true",
        default=None,
        help="Vary DeepCoy random seed per chunk as base+chunk_index.",
    )
    parser.add_argument(
        "--deepcoy-use-argmax-generation",
        dest="deepcoy_use_argmax_generation",
        action="store_true",
        help="Force DeepCoy use_argmax_generation=true (overrides config/base).",
    )
    parser.add_argument(
        "--deepcoy-no-argmax-generation",
        dest="deepcoy_use_argmax_generation",
        action="store_false",
        help="Force DeepCoy use_argmax_generation=false (overrides config/base).",
    )
    parser.add_argument(
        "--deepcoy-try-different-starting",
        dest="deepcoy_try_different_starting",
        action="store_true",
        help="Force DeepCoy try_different_starting=true (overrides config/base).",
    )
    parser.add_argument(
        "--deepcoy-no-try-different-starting",
        dest="deepcoy_try_different_starting",
        action="store_false",
        help="Force DeepCoy try_different_starting=false (overrides config/base).",
    )
    parser.add_argument(
        "--deepcoy-num-different-starting",
        type=int,
        default=None,
        help="Set DeepCoy num_different_starting (overrides config/base).",
    )
    parser.add_argument(
        "--deepcoy-num-samples",
        type=int,
        default=None,
        help="Set DeepCoy num_samples (overrides config/base).",
    )
    parser.set_defaults(
        deepcoy_use_argmax_generation=None,
        deepcoy_try_different_starting=None,
    )
    parser.add_argument(
        "--deepcoy-allowed-atoms",
        default=DEFAULT_DEEPCOY_ALLOWED_ATOMS_STR,
        help=(
            "Comma-separated atom symbols allowed for DeepCoy compatibility filtering "
            f"(default: {DEFAULT_DEEPCOY_ALLOWED_ATOMS_STR})."
        ),
    )
    parser.add_argument(
        "--no-deepcoy-filter-actives",
        dest="deepcoy_filter_actives",
        action="store_false",
        help="Disable DeepCoy compatibility filtering of actives before DeepCoy.",
    )
    parser.add_argument(
        "--fallback-smiles",
        default=None,
        help="Fallback SMILES to use only if no actives are found.",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output.")
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Skip external queries; rely on fallback/cached actives.",
    )
    parser.add_argument(
        "--force",
        dest="force",
        action="store_true",
        help="Force regeneration by removing existing outputs.",
    )
    parser.add_argument(
        "--force-deepcoy",
        dest="force_deepcoy",
        action="store_true",
        help="Force DeepCoy generation even if outputs appear to exist.",
    )

    return parser.parse_args()


def main():
    args = parse_args()
    repo_root = Path(__file__).resolve().parent
    project_root = repo_root.parent
    config_txt_path = project_root / CONFIG_FILE_NAME
    cfg = load_config(config_path=str(config_txt_path), base_dir=project_root)
    atlas_run_id = os.environ.get("ATLAS_RUN_ID")
    run_log_path = (
        project_root / "logs" / f"main_{atlas_run_id}.log" if atlas_run_id else None
    )

    pdb_id = args.pdb.upper()
    offline = bool(args.offline)
    force = bool(args.force or args.force_deepcoy)

    input_pdb_dir = resolve_path(args.input_pdb_dir, repo_root)
    out_root = resolve_path(args.out_root, repo_root)
    artifact_dir = (
        resolve_path(args.artifact_dir, repo_root)
        if args.artifact_dir
        else repo_root / "deepcoy_work" / pdb_id
    )
    cache_dir = artifact_dir / "cache"
    deepcoy_run_sh = resolve_path(args.deepcoy_run_sh, repo_root)
    deepcoy_sdf_sh = resolve_path(args.deepcoy_sdf_sh, repo_root)
    config_path = repo_root / "deepcoy_config.json"

    deepcoy_python = (
        args.deepcoy_python
        or os.environ.get("DEEPCOY_PYTHON")
        or cfg.get("DEEPCOY_PYTHON")
        or DEFAULT_DEEPCOY_PYTHON
    )
    decoys_per_active = (
        args.decoys_per_active
        if args.decoys_per_active is not None
        else cfg.get("DEEPCOY_DECOYS_PER_ACTIVE")
    )
    decoys_source = (
        "cli"
        if args.decoys_per_active is not None
        else "config"
        if cfg.get("DEEPCOY_DECOYS_PER_ACTIVE") is not None
        else "default"
    )
    if decoys_per_active is None or decoys_per_active <= 0:
        decoys_per_active = 50
        if decoys_source != "cli":
            decoys_source = "default"
    deepcoy_workers = max(1, args.deepcoy_workers)
    deepcoy_chunk_size = args.deepcoy_chunk_size
    if deepcoy_chunk_size is None:
        deepcoy_chunk_size = cfg.get("DEEPCOY_CHUNK_SIZE", 1)
    deepcoy_chunk_size = max(1, deepcoy_chunk_size or 1)
    deepcoy_threads = max(1, args.deepcoy_threads)
    deepcoy_keep_chunks = args.deepcoy_keep_chunks
    if deepcoy_keep_chunks is None:
        deepcoy_keep_chunks = cfg.get("DEEPCOY_KEEP_CHUNKS")
        if deepcoy_keep_chunks is None:
            deepcoy_keep_chunks = True
    deepcoy_keep_chunks = bool(deepcoy_keep_chunks)
    base_seed = args.deepcoy_base_seed
    if base_seed is None:
        base_seed = cfg.get("DEEPCOY_BASE_SEED", 0)
    seed_per_chunk = (
        args.deepcoy_seed_per_chunk
        if args.deepcoy_seed_per_chunk is not None
        else cfg.get("DEEPCOY_SEED_PER_CHUNK")
    )
    if seed_per_chunk is None:
        seed_per_chunk = False
    seed_per_chunk = bool(seed_per_chunk)
    use_argmax_generation = (
        args.deepcoy_use_argmax_generation
        if args.deepcoy_use_argmax_generation is not None
        else cfg.get("DEEPCOY_USE_ARGMAX_GENERATION")
    )
    try_different_starting = (
        args.deepcoy_try_different_starting
        if args.deepcoy_try_different_starting is not None
        else cfg.get("DEEPCOY_TRY_DIFFERENT_STARTING")
    )
    num_different_starting = (
        args.deepcoy_num_different_starting
        if args.deepcoy_num_different_starting is not None
        else cfg.get("DEEPCOY_NUM_DIFFERENT_STARTING")
    )
    # Only pass --num-samples when explicitly requested via CLI.
    # Config-driven num_samples has historically caused accidental overrides
    # that collapse output to 1 decoy/active even when DEEPCOY_DECOYS_PER_ACTIVE is larger.
    num_samples = args.deepcoy_num_samples
    allowed_atoms_list, _ = _parse_allowed_atom_list(args.deepcoy_allowed_atoms)
    allowed_atoms_str = ",".join(allowed_atoms_list)
    chembl_max_pages = args.chembl_max_pages
    cfg_chembl_pages = cfg.get("DEEPCOY_CHEMBL_MAX_PAGES")
    if cfg_chembl_pages is not None:
        try:
            chembl_max_pages = int(cfg_chembl_pages)
        except Exception:
            chembl_max_pages = chembl_max_pages
    chembl_max_pages = max(1, chembl_max_pages or 5)
    potency_cutoff_nm: Optional[float] = None
    cfg_potency_raw = cfg.get("DEEPCOY_ACTIVE_POTENCY_CUTOFF_NM")
    if cfg_potency_raw is not None:
        try:
            cutoff_str = str(cfg_potency_raw).strip()
            potency_cutoff_nm = float(cutoff_str) if cutoff_str else None
        except Exception:
            potency_cutoff_nm = None
    potency_keep_unknown = True
    if "DEEPCOY_ACTIVE_POTENCY_KEEP_UNKNOWN" in cfg:
        potency_keep_unknown = bool(cfg.get("DEEPCOY_ACTIVE_POTENCY_KEEP_UNKNOWN"))

    resolved_sources, sources_origin, raw_sources_value = resolve_sources(
        args.sources, cfg
    )
    sources, unknown_sources = _filter_known_sources(resolved_sources)
    source_opts = {
        "timeout": args.http_timeout,
        "retries": args.http_retries,
        "max_actives_per_source": args.max_actives_per_source,
        "affinity_cutoff_nm": args.affinity_cutoff_nm,
        "iuphar_approved_only": args.iuphar_approved_only,
        "iuphar_primary_only": args.iuphar_primary_only,
        "chembl_activity_types": [
            t.strip() for t in args.chembl_activity_types.split(",") if t.strip()
        ],
        "chembl_max_phase": args.chembl_max_phase,
        "chembl_max_pages": chembl_max_pages,
        "drugbank_data_dir": resolve_path(args.drugbank_data_dir, repo_root)
        if args.drugbank_data_dir
        else None,
    }
    if unknown_sources:
        warning_line = f"[deepcoy.sources] warning=unknown_sources ignored={','.join(unknown_sources)}"
        print(warning_line)
        _tee_run_log(run_log_path, warning_line)
    if not sources:
        sources = list(DEFAULT_SOURCES)
        fallback_line = f"[deepcoy.sources] warning=no_valid_sources_fallback resolved={','.join(sources)}"
        print(fallback_line)
        _tee_run_log(run_log_path, fallback_line)
        if sources_origin != "cli":
            sources_origin = "default"
        if not raw_sources_value:
            raw_sources_value = ",".join(DEFAULT_SOURCES)
    sources_log_line = f"[deepcoy.sources] source={sources_origin} value={raw_sources_value} resolved={','.join(sources)}"
    print(sources_log_line)
    _tee_run_log(run_log_path, sources_log_line)
    potency_cfg_line = (
        "[deepcoy.actives.potency_config] "
        f"cutoff_nm={potency_cutoff_nm if potency_cutoff_nm is not None else 'none'} "
        f"keep_unknown={potency_keep_unknown}"
    )
    print(potency_cfg_line)
    _tee_run_log(run_log_path, potency_cfg_line)
    audit_enabled = bool(cfg.get("DEEPCOY_SOURCE_AUDIT", False))
    try:
        audit_max_lines = int(cfg.get("DEEPCOY_SOURCE_AUDIT_MAX_LINES", 200) or 200)
    except Exception:
        audit_max_lines = 200
    audit = SourceAudit(enabled=audit_enabled, max_lines=audit_max_lines)

    if args.verbose:
        print(f"-> Repo root: {repo_root}")
        print(f"-> Project root: {project_root}")
        print(f"-> Config.txt: {config_txt_path}")
        print(f"-> Input PDB dir: {input_pdb_dir}")
        print(f"-> Output root: {out_root}")
        print(f"-> DeepCoy python: {deepcoy_python}")
        print(f"-> Decoys per active: {decoys_per_active}")
        print(f"-> DeepCoy workers: {deepcoy_workers}")
        print(f"-> DeepCoy chunk size: {deepcoy_chunk_size}")
        print(f"-> DeepCoy threads: {deepcoy_threads}")
        print(f"-> DeepCoy keep chunks: {deepcoy_keep_chunks}")
        print(f"-> DeepCoy allowed atoms: {allowed_atoms_str}")
        print(f"-> DeepCoy filter actives: {args.deepcoy_filter_actives}")
        print(f"-> DeepCoy base seed: {base_seed}")
        print(f"-> DeepCoy seed per chunk: {seed_per_chunk}")
        print(f"-> DeepCoy use_argmax_generation: {use_argmax_generation}")
        print(f"-> DeepCoy try_different_starting: {try_different_starting}")
        print(f"-> DeepCoy num_different_starting: {num_different_starting}")
        print(f"-> DeepCoy num_samples (CLI-only): {num_samples}")
        print(f"-> Artifact dir: {artifact_dir}")
        print(f"-> Sources: {sources}")
        print(f"-> Source workers: {args.max_source_workers}")
        print(f"-> Ensure SDF: {args.ensure_sdf}")

    if not input_pdb_dir.is_dir():
        print(
            f"ERROR: input-pdb-dir does not exist or is not a directory: {input_pdb_dir}",
            file=sys.stderr,
        )
        return 2

    pdb_file = find_pdb_file(input_pdb_dir, pdb_id)
    if not pdb_file:
        expected = ", ".join([f"{pdb_id}{ext}" for ext in PDB_EXTENSIONS])
        print(
            f"ERROR: PDB file for {pdb_id} not found in {input_pdb_dir}. "
            f"Looked for: {expected}",
            file=sys.stderr,
        )
        return 2

    if args.verbose:
        print(f"-> Found input PDB file: {pdb_file}")

    fallback_smiles_cli = args.fallback_smiles
    if offline:
        print("## Offline mode enabled: skipping UniProt/EC queries.")
        uniprot_id, ec_numbers = "UNKNOWN", []
    else:
        uniprot_id, ec_numbers = get_uniprot_and_ec(pdb_id)
        if not uniprot_id:
            if fallback_smiles_cli:
                print(
                    "WARNING: UniProt lookup failed; proceeding with fallback SMILES and UNKNOWN tag."
                )
                uniprot_id, ec_numbers = "UNKNOWN", []
            else:
                print(
                    "ERROR: Failed to retrieve UniProt ID. Cannot proceed.",
                    file=sys.stderr,
                )
                return 1
    audit_lines = []
    if not offline:
        audit_lines = build_source_audit_lines(
            uniprot_id or "UNKNOWN",
            ec_numbers or [],
            pdb_id,
            source_opts,
            sources,
        )
        for line in audit_lines:
            print(line)
            _tee_run_log(run_log_path, line)

    all_actives: List[Tuple[str, str]] = []
    smiles_to_sources: Dict[str, Set[str]] = {}
    provenance_details: Dict[str, List[dict]] = {}
    if not offline:
        all_actives, smiles_to_sources, provenance_details = query_external_sources(
            uniprot_id,
            ec_numbers,
            pdb_id,
            sources,
            cache_dir,
            source_opts,
            max_workers=args.max_source_workers,
            run_log_path=run_log_path,
            audit=audit,
        )
    if audit and audit.enabled:
        for line in audit.summary_lines():
            print(line)
            _tee_run_log(run_log_path, line)
    (
        final_actives_smiles,
        smiles_to_sources,
        provenance_details,
    ) = validate_and_filter_actives(
        all_actives,
        smiles_to_sources=smiles_to_sources,
        provenance_details=provenance_details,
        potency_cutoff_nm=potency_cutoff_nm,
        potency_keep_unknown=potency_keep_unknown,
        run_log_path=run_log_path,
    )

    if not final_actives_smiles:
        if args.fallback_smiles:
            print("WARNING: No actives found; using fallback SMILES for smoke testing.")
            (
                final_actives_smiles,
                smiles_to_sources,
                provenance_details,
            ) = validate_and_filter_actives(
                [(args.fallback_smiles, "fallback")],
                smiles_to_sources=smiles_to_sources,
                provenance_details=provenance_details,
                potency_cutoff_nm=potency_cutoff_nm,
                potency_keep_unknown=potency_keep_unknown,
                run_log_path=run_log_path,
            )

        if not final_actives_smiles:
            print(
                "WARNING: No valid active metabolites found. Provide --fallback-smiles to proceed.",
                file=sys.stderr,
            )
            return 3

    label = f"{pdb_id}_{(uniprot_id or 'UNKNOWN').split('-')[0]}"
    run_tag = build_deepcoy_run_tag() if deepcoy_keep_chunks else None
    output_dir = out_root / label
    staging_dir = None
    work_output_dir = output_dir

    def _preserve_staging(src: Path, base_dest: Path, reason: str) -> Optional[Path]:
        if not src.exists():
            return None
        dest = base_dest
        idx = 1
        while dest.exists():
            dest = base_dest.with_name(f"{base_dest.name}_{idx}")
            idx += 1
        src.rename(dest)
        msg = f"[deepcoy.staging.preserved] path={dest} reason={reason}"
        print(msg)
        _tee_run_log(run_log_path, msg)
        return dest

    if args.run_deepcoy:
        staging_dir = output_dir.parent / f".{label}_staging"
        if staging_dir.exists():
            if deepcoy_keep_chunks:
                _preserve_staging(
                    staging_dir,
                    output_dir.parent / f".{label}_staging_prev_{run_tag or 'run'}",
                    "existing_staging",
                )
            else:
                shutil.rmtree(staging_dir, ignore_errors=True)
        work_output_dir = staging_dir
    if force and artifact_dir.exists():
        shutil.rmtree(artifact_dir, ignore_errors=True)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    def _cleanup_staging(reason: str = "cleanup") -> None:
        if staging_dir and staging_dir.exists():
            if deepcoy_keep_chunks:
                _preserve_staging(
                    staging_dir,
                    output_dir.parent / f".{label}_staging_failed_{run_tag or 'run'}",
                    reason,
                )
            else:
                shutil.rmtree(staging_dir, ignore_errors=True)

    if audit and audit.enabled:
        write_source_audit_artifact(
            audit, cache_dir / "source_audit.json", run_log_path=run_log_path
        )

    actives_smi_path, _ = write_actives_smiles(
        final_actives_smiles, work_output_dir, label
    )
    try:
        _, filtered_kept, _ = filter_actives_for_deepcoy(
            actives_smi_path,
            allowed_atoms_list,
            filter_enabled=args.deepcoy_filter_actives,
            run_log_path=run_log_path,
        )
        if deepcoy_keep_chunks:
            copy_filter_debug_artifacts(
                actives_smi_path, artifact_dir, run_tag, run_log_path=run_log_path
            )
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        if deepcoy_keep_chunks:
            copy_filter_debug_artifacts(
                actives_smi_path, artifact_dir, run_tag, run_log_path=run_log_path
            )
        _cleanup_staging()
        return 1
    except FileNotFoundError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        if deepcoy_keep_chunks:
            copy_filter_debug_artifacts(
                actives_smi_path, artifact_dir, run_tag, run_log_path=run_log_path
            )
        _cleanup_staging()
        return 1
    filtered_entries = _read_smiles_with_labels(actives_smi_path)
    filtered_sources = {
        smi: smiles_to_sources.get(smi, {"unknown"}) for _, smi in filtered_entries
    }
    filtered_details = {
        smi: provenance_details.get(smi, []) for _, smi in filtered_entries
    }
    prov_tsv, prov_json = write_provenance_files(
        label, cache_dir, filtered_entries, filtered_sources, filtered_details
    )
    if prov_tsv:
        prov_line = (
            f"[deepcoy.provenance] tsv={prov_tsv} json={prov_json or 'none'} "
            f"n={len(filtered_entries)}"
        )
        print(prov_line)
        _tee_run_log(run_log_path, prov_line)
    actives_sdf_path = work_output_dir / "deepcoy_actives.sdf"
    convert_smi_to_sdf(actives_smi_path, actives_sdf_path)
    actives_log_line = (
        f"[deepcoy.actives.sdf] tag={label} out={actives_sdf_path} n={filtered_kept}"
    )
    print(actives_log_line)
    _tee_run_log(run_log_path, actives_log_line)

    if not args.run_deepcoy:
        if staging_dir and staging_dir.exists():
            shutil.rmtree(staging_dir, ignore_errors=True)
        print("DeepCoy run skipped (--no-run-deepcoy).")
        return 0

    if not deepcoy_run_sh.is_file():
        print(f"ERROR: deepcoy_run.sh not found: {deepcoy_run_sh}", file=sys.stderr)
        _cleanup_staging()
        return 2
    if not config_path.is_file():
        print(f"ERROR: deepcoy_config.json not found: {config_path}", file=sys.stderr)
        _cleanup_staging()
        return 2
    if not args.skip_sdf and not deepcoy_sdf_sh.is_file():
        print(
            f"ERROR: deepcoy_smiles_to_sdf.sh not found: {deepcoy_sdf_sh}",
            file=sys.stderr,
        )
        _cleanup_staging()
        return 2

    if args.verbose:
        print(f"-> DeepCoy python: {deepcoy_python}")
        print(f"-> DeepCoy run script: {deepcoy_run_sh}")
        print(f"-> DeepCoy config: {config_path}")

    config_log = (
        f"[deepcoy.config] decoys_per_active={decoys_per_active} source={decoys_source}"
    )
    print(config_log)
    _tee_run_log(run_log_path, config_log)

    try:
        run_deepcoy_workflow(
            actives_smi_path,
            work_output_dir,
            artifact_dir,
            deepcoy_run_sh,
            deepcoy_sdf_sh,
            decoys_per_active,
            deepcoy_python,
            config_path,
            args.skip_sdf,
            args.ensure_sdf,
            args.restrict_data,
            run_log_path=run_log_path,
            deepcoy_workers=deepcoy_workers,
            deepcoy_chunk_size=deepcoy_chunk_size,
            deepcoy_threads=deepcoy_threads,
            deepcoy_keep_chunks=deepcoy_keep_chunks,
            default_model_path=repo_root / DEFAULT_DEEPCOY_MODEL,
            phosphorus_model_path=repo_root / DEFAULT_DEEPCOY_PHOS_MODEL,
            base_seed=base_seed,
            seed_per_chunk=seed_per_chunk,
            use_argmax_generation=use_argmax_generation,
            try_different_starting=try_different_starting,
            num_different_starting=num_different_starting,
            num_samples=num_samples,
            run_tag=run_tag,
        )
    except FileNotFoundError:
        print(
            "ERROR: DeepCoy wrapper script not found or not executable.",
            file=sys.stderr,
        )
        _cleanup_staging("wrapper_not_found")
        return 1
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        _cleanup_staging("runtime_error")
        return 1
    except subprocess.CalledProcessError as e:
        if e.returncode == 4 and args.fallback_smiles:
            print(
                "WARNING: DeepCoy preprocessing rejected actives; retrying with fallback ligand."
            )
            simple_smiles = "CCN(CC)CCO"
            with actives_smi_path.open("w") as f:
                f.write(f"{simple_smiles}\tFALLBACK\n")
            run_tag_fallback = f"{run_tag}_fallback" if run_tag else None
            try:
                run_deepcoy_workflow(
                    actives_smi_path,
                    work_output_dir,
                    artifact_dir,
                    deepcoy_run_sh,
                    deepcoy_sdf_sh,
                    decoys_per_active,
                    deepcoy_python,
                    config_path,
                    args.skip_sdf,
                    args.ensure_sdf,
                    args.restrict_data,
                    run_log_path=run_log_path,
                    deepcoy_workers=deepcoy_workers,
                    deepcoy_chunk_size=deepcoy_chunk_size,
                    deepcoy_threads=deepcoy_threads,
                    deepcoy_keep_chunks=deepcoy_keep_chunks,
                    default_model_path=repo_root / DEFAULT_DEEPCOY_MODEL,
                    phosphorus_model_path=repo_root / DEFAULT_DEEPCOY_PHOS_MODEL,
                    base_seed=base_seed,
                    seed_per_chunk=seed_per_chunk,
                    use_argmax_generation=use_argmax_generation,
                    try_different_starting=try_different_starting,
                    num_different_starting=num_different_starting,
                    num_samples=num_samples,
                    run_tag=run_tag_fallback,
                )
            except subprocess.CalledProcessError as e2:
                print(
                    f"ERROR: DeepCoy command failed after retry with exit code {e2.returncode}.",
                    file=sys.stderr,
                )
                _cleanup_staging("fallback_retry_failed")
                return e2.returncode or 1
        else:
            print(
                f"ERROR: DeepCoy command failed with exit code {e.returncode}.",
                file=sys.stderr,
            )
            _cleanup_staging("deepcoy_failed")
            return e.returncode or 1

    if staging_dir and staging_dir.exists():
        if output_dir.exists():
            shutil.rmtree(output_dir, ignore_errors=True)
        staging_dir.rename(output_dir)

    final_line = f"Decoys should be generated in: {output_dir}"
    print(final_line)
    _tee_run_log(run_log_path, final_line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
# ruff: noqa: E402
