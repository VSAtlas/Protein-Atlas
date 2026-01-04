import argparse
import os
import sys
from pathlib import Path
import subprocess
import re
import json
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from external_sources import (
    fetch_bindingdb_smiles,
    fetch_chembl_smiles,
    fetch_drugbank_smiles,
    fetch_iuphar_smiles,
)

DEFAULT_INPUT_PDB_DIR = "/home/michael/atlas/code/protein_automation/input_pdbs"
DEFAULT_OUT_ROOT = "../extracted_ligands/deepcoy"
DEFAULT_DEEPCOY_PYTHON = "/home/michael/atlas/anaconda3/envs/DeepCoy-env-run/bin/python"
PDB_EXTENSIONS = [".pdb", ".ent", ".cif", ".pdb.gz", ".cif.gz"]
CONFIG_FILE_NAME = "config.txt"
EC_PATTERN = re.compile(r"^(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)$")
DEFAULT_SOURCES = ["chembl", "bindingdb", "iuphar", "rcsb", "rhea"]
SOURCE_ORDER = ["chembl", "bindingdb", "iuphar", "drugbank", "rcsb", "rhea"]


def read_config_value(config_path, key):
    try:
        with open(config_path, "r") as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                if "=" not in stripped:
                    continue
                k, v = stripped.split("=", 1)
                if k.strip() != key:
                    continue
                value = v.strip().strip('"').strip("'")
                return value
    except FileNotFoundError:
        return None
    except Exception:
        return None
    return None

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
    rcsb_entity = entity_data.get("rcsb_polymer_entity", {}) if isinstance(entity_data, dict) else {}

    combined = rcsb_entity.get("rcsb_enzyme_class_combined", [])
    if isinstance(combined, str):
        combined = [combined]
    for entry in combined or []:
        _extract_ec_from_string(entry, ec_numbers)

    enzyme_class = rcsb_entity.get("rcsb_enzyme_class") or rcsb_entity.get("rcsb_enzyme_class_list") or []
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
        annotation_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
        try:
            print(f"-> Querying RCSB Annotation API for {pdb_id} (Entity {entity_id})...")
            response = requests.get(annotation_url)
            if response.status_code != 200:
                continue
            entity_data = response.json()

            identifiers = entity_data.get("rcsb_polymer_entity_container_identifiers", {})
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
        print("Warning: Could not find any EC numbers. Phase 2 (Metabolite collection) may be limited.")

    return uniprot_id, final_ec_list

# --- 2. Collecting Actives/Metabolites (PHASE 2) ---

def query_external_sources(uniprot_id, ec_numbers, pdb_id, sources, cache_dir, source_opts, max_workers=1):
    """
    (2) Queries external data sources for ligands/metabolites.
    Returns a list of (SMILES, source) tuples.
    """
    print("## 2. Querying external data sources for known metabolites...")
    all_actives = []

    def add_active(smiles, source):
        if smiles and isinstance(smiles, str) and not any(smiles == active[0] for active in all_actives):
            all_actives.append((smiles, source))

    def _print_meta(src, meta):
        if not meta:
            return
        counts = meta.get("counts", {})
        cached = meta.get("cached", False)
        if src == "chembl":
            print(f"[chembl] targets={counts.get('targets', 0)} assays={counts.get('assays', 0)} activities={counts.get('activities', 0)} smiles={len(meta.get('smiles', [])) if 'smiles' in meta else counts.get('molecules', 0) or len(all_actives)} cached={cached}")
        elif src == "bindingdb":
            print(f"[bindingdb] pdb_hits={counts.get('pdb_hits', 0)} uniprot_hits={counts.get('uniprot_hits', 0)} smiles={len(meta.get('smiles', [])) if 'smiles' in meta else counts.get('pdb_hits', 0) + counts.get('uniprot_hits', 0)} cached={cached}")
        elif src == "iuphar":
            print(f"[iuphar] targets={counts.get('targets', 0)} interactions={counts.get('interactions', 0)} smiles={len(meta.get('smiles', [])) if 'smiles' in meta else counts.get('interactions', 0)} cached={cached}")
        elif src == "drugbank":
            print(f"[drugbank] smiles={len(meta.get('smiles', [])) if 'smiles' in meta else 0} info={meta}")

    cache_dir.mkdir(parents=True, exist_ok=True)

    source_set = set(sources)
    tasks = {}

    def run_rhea():
        local_smiles = []
        print(f"-> Querying Rhea/ChEBI for UniProt {uniprot_id} (Metabolic reactions)...")
        rhea_uniprot_url = f"https://www.rhea-db.org/rest/ws/reaction/uniprot/{uniprot_id}"
        try:
            with requests.Session() as sess:
                rhea_response = sess.get(rhea_uniprot_url, timeout=source_opts["timeout"])
                if rhea_response.status_code == 200 and rhea_response.text.strip().startswith("["):
                    rhea_data = rhea_response.json()
                    for reaction in rhea_data:
                        for side in reaction.get("reactionSides", []):
                            for component in side.get("reactionComponents", []):
                                chebi_id = component.get("chebi")
                                if chebi_id:
                                    chebi_smiles_url = (
                                        f"https://www.ebi.ac.uk/chebi/rest/api/getCompleteEntity?chebiId={chebi_id}"
                                    )
                                    chebi_xml_response = sess.get(chebi_smiles_url, timeout=source_opts["timeout"])
                                    if chebi_xml_response.status_code == 200 and "<smiles>" in chebi_xml_response.text:
                                        start_tag = "<smiles>"
                                        end_tag = "</smiles>"
                                        smiles_start = chebi_xml_response.text.find(start_tag) + len(start_tag)
                                        smiles_end = chebi_xml_response.text.find(end_tag)
                                        smiles = chebi_xml_response.text[smiles_start:smiles_end].strip()
                                        if smiles:
                                            local_smiles.append((smiles, "Rhea/ChEBI"))
                else:
                    print("-> Rhea returned an empty or non-JSON response. Skipping ChEBI extraction.")
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
                ligand_response = sess.get(ligand_url, timeout=source_opts["timeout"]).json()
                for entity in ligand_response.get("nonpolymer_entity", []):
                    comp_id = entity.get("pdbx_nonpolymer_entity", {}).get("chem_comp_id")
                    if comp_id:
                        chem_comp_url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}"
                        comp_response = sess.get(chem_comp_url, timeout=source_opts["timeout"]).json()
                        known_smiles = comp_response.get("pdbx_chem_comp_descriptor", [{}])
                        for desc in known_smiles:
                            if desc.get("type") == "SMILES":
                                smi = desc.get("descriptor")
                                if smi:
                                    local_smiles.append((smi, f"PDB/ChemComp/{comp_id}"))
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
            for smi, label in results:
                add_active(smi, label)

    print(f"-> Found {len(all_actives)} raw active candidates from external sources.")
    return all_actives

# --- 3. Filtering and Validation (PHASE 3) ---

def validate_and_filter_actives(actives_list):
    """
    (3) Merges, deduplicates, and filters the list of active metabolites
        based on RDKit validation and size.
    """
    print("## 3. Merging, Deduplicating, and Filtering actives...")

    unique_smiles = {smiles for smiles, source in actives_list if smiles and isinstance(smiles, str)}
    valid_actives = []

    print(f"-> Found {len(unique_smiles)} unique SMILES strings.")

    EXCLUDE_SMILES = {"O", "C", "[Na+]", "[Cl-]", "[Mg+2]", "[K+]"}

    for smiles in unique_smiles:
        if smiles in EXCLUDE_SMILES:
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

        except Exception:
            continue

    print(f"-> Final list contains {len(valid_actives)} filtered, unique active metabolites.")
    return valid_actives

# --- 4. DeepCoy Execution (PHASE 4) ---

def write_actives_smiles(final_actives_smiles, output_dir, label):
    """
    (4) Writes actives to a .smi file for DeepCoy.
    """
    print("## 4. Writing SMILES file...")

    output_dir.mkdir(parents=True, exist_ok=True)
    actives_smi_path = output_dir / f"{label}_actives.smi"

    # Write the Actives SMILES file (SMILES [tab] Name)
    with actives_smi_path.open("w") as f:
        for i, smiles in enumerate(final_actives_smiles):
            name = f"METAB_{i + 1}"
            f.write(f"{smiles}\t{name}\n")

    print(f"-> Active SMILES written to: {actives_smi_path}")
    return actives_smi_path


def convert_smi_to_sdf(smi_path: Path, sdf_path: Path):
    if not smi_path.is_file():
        print(f"WARNING: SMILES file not found for SDF conversion: {smi_path}")
        return
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
        return
    sdf_path.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(sdf_path))
    count = 0
    for mol in suppl:
        writer.write(mol)
        count += 1
    writer.close()
    print(f"-> Wrote {count} decoys to SDF: {sdf_path} ({sdf_path.stat().st_size} bytes)")


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
):
    print("## 5. Executing DeepCoy workflow...")

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
    if restrict_data and int(restrict_data) > 0:
        deepcoy_cmd.extend(["--restrict-data", str(restrict_data)])
    print(f"Running: {' '.join(deepcoy_cmd)}")
    subprocess.run(deepcoy_cmd, check=True)

    smi_path = output_dir / "deepcoy_decoys.smi"
    sdf_path = output_dir / "deepcoy_decoys.sdf"
    if ensure_sdf or not skip_sdf:
        convert_smi_to_sdf(smi_path, sdf_path)
    else:
        print("-> Skipping SDF conversion (--skip-sdf).")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a DeepCoy DUD-style decoy library from a target PDB.")

    parser.add_argument("--pdb", required=True, help="PDB ID to query (required).")
    parser.add_argument(
        "--input-pdb-dir",
        default=DEFAULT_INPUT_PDB_DIR,
        help="Directory containing input PDBs.")
    parser.add_argument(
        "--out-root",
        default=DEFAULT_OUT_ROOT,
        help="Output root for actives/decoys.")
    parser.add_argument(
        "--artifact-dir",
        default=None,
        help="Directory for DeepCoy run artifacts (defaults to deepcoy_work/<PDB>).")
    parser.add_argument(
        "--sources",
        default="chembl,bindingdb,iuphar,rcsb,rhea",
        help="Comma-separated sources to query (default: chembl,bindingdb,iuphar,rcsb,rhea).",
    )
    parser.add_argument(
        "--http-timeout",
        type=int,
        default=20,
        help="HTTP timeout (seconds) for external source queries.")
    parser.add_argument(
        "--http-retries",
        type=int,
        default=2,
        help="HTTP retries for external source queries.")
    parser.add_argument(
        "--max-actives-per-source",
        type=int,
        default=500,
        help="Maximum SMILES to retain per source.")
    parser.add_argument(
        "--affinity-cutoff-nm",
        type=int,
        default=10000,
        help="Affinity cutoff in nM for potency-filtered sources (e.g., BindingDB, ChEMBL).")
    parser.add_argument(
        "--iuphar-approved-only",
        dest="iuphar_approved_only",
        action="store_true",
        help="Use only approved IUPHAR interactions (default).")
    parser.add_argument(
        "--iuphar-include-nonapproved",
        dest="iuphar_approved_only",
        action="store_false",
        help="Allow non-approved IUPHAR interactions.")
    parser.add_argument(
        "--iuphar-primary-target-only",
        dest="iuphar_primary_only",
        action="store_true",
        help="Use only primary targets in IUPHAR (default).")
    parser.add_argument(
        "--iuphar-include-secondary-targets",
        dest="iuphar_primary_only",
        action="store_false",
        help="Allow secondary targets in IUPHAR.")
    parser.add_argument(
        "--chembl-activity-types",
        default="Ki,Kd,IC50",
        help="Comma-separated ChEMBL activity types (default: Ki,Kd,IC50).")
    parser.add_argument(
        "--chembl-max-phase",
        type=int,
        default=None,
        help="Optional ChEMBL max_phase upper bound to bias toward approved drugs.")
    parser.add_argument(
        "--max-source-workers",
        type=int,
        default=1,
        help="Max threads for parallel source queries (default 1 = sequential).")
    parser.add_argument(
        "--drugbank-data-dir",
        default=None,
        help="Path to local DrugBank structured data. If unset, DrugBank is skipped.")

    run_group = parser.add_mutually_exclusive_group()
    run_group.add_argument(
        "--run-deepcoy",
        dest="run_deepcoy",
        action="store_true",
        help="Run DeepCoy after preparing actives.")
    run_group.add_argument(
        "--no-run-deepcoy",
        dest="run_deepcoy",
        action="store_false",
        help="Skip DeepCoy run.")
    run_group.add_argument(
        "--skip-deepcoy",
        dest="run_deepcoy",
        action="store_false",
        help="Alias for --no-run-deepcoy.")
    parser.set_defaults(run_deepcoy=True)
    parser.set_defaults(iuphar_approved_only=True, iuphar_primary_only=True)

    parser.add_argument(
        "--deepcoy-python",
        default=None,
        help=(
            "Path to DeepCoy Python interpreter. Defaults to $DEEPCOY_PYTHON, "
            "config.txt DEEPCOY_PYTHON, or the built-in DeepCoy path."
        ),
    )
    parser.add_argument(
        "--deepcoy-run-sh",
        default="./deepcoy_run.sh",
        help="Path to deepcoy_run.sh.")
    parser.add_argument(
        "--deepcoy-sdf-sh",
        default="./deepcoy_smiles_to_sdf.sh",
        help="Path to deepcoy_smiles_to_sdf.sh.")
    parser.add_argument(
        "--ensure-sdf",
        action="store_true",
        help="Force SDF generation with RDKit even if --skip-sdf is set.")
    parser.add_argument(
        "--restrict-data",
        type=int,
        default=0,
        help="Restrict DeepCoy data size for faster smoke tests (passed to DeepCoy).")

    parser.add_argument(
        "--skip-sdf",
        action="store_true",
        help="Skip SDF conversion step after DeepCoy.")
    parser.add_argument(
        "--decoys-per-active",
        type=int,
        default=50,
        help="Number of decoys to generate per active.")
    parser.add_argument(
        "--fallback-smiles",
        default=None,
        help="Fallback SMILES to use only if no actives are found.")
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output.")
    parser.add_argument(
        "--force-deepcoy",
        dest="force_deepcoy",
        action="store_true",
        help="Force DeepCoy generation even if outputs appear to exist.")

    return parser.parse_args()


def main():
    args = parse_args()
    repo_root = Path(__file__).resolve().parent
    project_root = repo_root.parent
    config_txt_path = project_root / CONFIG_FILE_NAME

    pdb_id = args.pdb.upper()

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
        or read_config_value(config_txt_path, "DEEPCOY_PYTHON")
        or DEFAULT_DEEPCOY_PYTHON
    )

    sources = [s.strip().lower() for s in (args.sources or "").split(",") if s.strip()]
    if not sources:
        sources = DEFAULT_SOURCES
    source_opts = {
        "timeout": args.http_timeout,
        "retries": args.http_retries,
        "max_actives_per_source": args.max_actives_per_source,
        "affinity_cutoff_nm": args.affinity_cutoff_nm,
        "iuphar_approved_only": args.iuphar_approved_only,
        "iuphar_primary_only": args.iuphar_primary_only,
        "chembl_activity_types": [t.strip() for t in args.chembl_activity_types.split(",") if t.strip()],
        "chembl_max_phase": args.chembl_max_phase,
        "drugbank_data_dir": resolve_path(args.drugbank_data_dir, repo_root) if args.drugbank_data_dir else None,
    }

    if args.verbose:
        print(f"-> Repo root: {repo_root}")
        print(f"-> Project root: {project_root}")
        print(f"-> Config.txt: {config_txt_path}")
        print(f"-> Input PDB dir: {input_pdb_dir}")
        print(f"-> Output root: {out_root}")
        print(f"-> DeepCoy python: {deepcoy_python}")
        print(f"-> Artifact dir: {artifact_dir}")
        print(f"-> Sources: {sources}")
        print(f"-> Source workers: {args.max_source_workers}")
        print(f"-> Ensure SDF: {args.ensure_sdf}")

    if not input_pdb_dir.is_dir():
        print(f"ERROR: input-pdb-dir does not exist or is not a directory: {input_pdb_dir}", file=sys.stderr)
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

    uniprot_id, ec_numbers = get_uniprot_and_ec(pdb_id)
    if not uniprot_id:
        print("ERROR: Failed to retrieve UniProt ID. Cannot proceed.", file=sys.stderr)
        return 1

    all_actives = query_external_sources(
        uniprot_id,
        ec_numbers,
        pdb_id,
        sources,
        cache_dir,
        source_opts,
        max_workers=args.max_source_workers,
    )
    final_actives_smiles = validate_and_filter_actives(all_actives)

    if not final_actives_smiles:
        if args.fallback_smiles:
            print("WARNING: No actives found; using fallback SMILES for smoke testing.")
            final_actives_smiles = validate_and_filter_actives([(args.fallback_smiles, "fallback")])

        if not final_actives_smiles:
            print(
                "WARNING: No valid active metabolites found. Provide --fallback-smiles to proceed.",
                file=sys.stderr,
            )
            return 3

    label = f"{pdb_id}_{uniprot_id.split('-')[0]}"
    output_dir = out_root / label
    artifact_dir.mkdir(parents=True, exist_ok=True)

    actives_smi_path = write_actives_smiles(final_actives_smiles, output_dir, label)

    if not args.run_deepcoy:
        print("DeepCoy run skipped (--no-run-deepcoy).")
        return 0

    if not deepcoy_run_sh.is_file():
        print(f"ERROR: deepcoy_run.sh not found: {deepcoy_run_sh}", file=sys.stderr)
        return 2
    if not config_path.is_file():
        print(f"ERROR: deepcoy_config.json not found: {config_path}", file=sys.stderr)
        return 2
    if not args.skip_sdf and not deepcoy_sdf_sh.is_file():
        print(f"ERROR: deepcoy_smiles_to_sdf.sh not found: {deepcoy_sdf_sh}", file=sys.stderr)
        return 2

    if args.verbose:
        print(f"-> DeepCoy python: {deepcoy_python}")
        print(f"-> DeepCoy run script: {deepcoy_run_sh}")
        print(f"-> DeepCoy config: {config_path}")

    try:
        run_deepcoy_workflow(
            actives_smi_path,
            output_dir,
            artifact_dir,
            deepcoy_run_sh,
            deepcoy_sdf_sh,
            args.decoys_per_active,
            deepcoy_python,
            config_path,
            args.skip_sdf,
            args.ensure_sdf,
            args.restrict_data,
        )
    except FileNotFoundError:
        print("ERROR: DeepCoy wrapper script not found or not executable.", file=sys.stderr)
        return 1
    except subprocess.CalledProcessError as e:
        print(f"ERROR: DeepCoy command failed with exit code {e.returncode}.", file=sys.stderr)
        return e.returncode or 1

    print(f"Decoys should be generated in: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
