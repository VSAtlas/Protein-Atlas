import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

# Setup paths
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# Import fallback module
DEEPCOY_PATH = REPO_ROOT / "DeepCoy_duds"
if str(DEEPCOY_PATH) not in sys.path:
    sys.path.insert(0, str(DEEPCOY_PATH))

try:
    from DeepCoy_duds.generate_dud_library import get_uniprot_and_ec
except ImportError:
    get_uniprot_and_ec = None


def _parse_int_field(raw: str) -> Optional[int]:
    raw = raw.strip()
    if not raw:
        return None
    try:
        return int(raw)
    except ValueError:
        return None


def _merge_entry(
    existing: Optional[Dict[str, Any]], new_entry: Dict[str, Any]
) -> Dict[str, Any]:
    if not existing:
        return dict(new_entry)
    merged = dict(existing)
    for key in ("uniprot", "unp_start", "unp_end", "pdb_start", "pdb_end", "source"):
        if merged.get(key) is None and new_entry.get(key) is not None:
            merged[key] = new_entry[key]
    return merged


def get_polymer_chains(pdb_path: Path) -> list[str]:
    chains = set()
    with pdb_path.open("r") as f:
        for line in f:
            if line.startswith("ATOM  "):
                if len(line) > 21:
                    chains.add(line[21])
    return sorted(list(chains))


def _parse_dbref_line(line: str) -> Dict[str, Any]:
    entry: Dict[str, Any] = {
        "chain": None,
        "uniprot": None,
        "unp_start": None,
        "unp_end": None,
        "pdb_start": None,
        "pdb_end": None,
        "source": "dbref",
    }
    if len(line) > 12:
        chain_id = line[12].strip()
        if chain_id:
            entry["chain"] = chain_id
    # PDB residue span (may be missing)
    if len(line) >= 24:
        entry["pdb_start"] = _parse_int_field(line[14:18])
        entry["pdb_end"] = _parse_int_field(line[20:24])

    db_name = line[26:32].strip() if len(line) >= 32 else ""
    if db_name == "UNP":
        accession = line[32:41].strip() if len(line) >= 41 else ""
        if accession:
            entry["uniprot"] = accession
        if len(line) >= 60:
            entry["unp_start"] = _parse_int_field(line[55:60])
        if len(line) >= 67:
            entry["unp_end"] = _parse_int_field(line[62:67])

    # Token fallback for odd spacing
    if not entry["chain"]:
        parts = line.split()
        if len(parts) > 2 and len(parts[2]) == 1:
            entry["chain"] = parts[2]
    if not entry["uniprot"]:
        parts = line.split()
        if "UNP" in parts:
            try:
                idx = parts.index("UNP")
                if idx + 1 < len(parts):
                    acc = parts[idx + 1]
                    if re.match(r"^[A-Z0-9]{6,10}$", acc):
                        entry["uniprot"] = acc
                if entry["unp_start"] is None and idx + 2 < len(parts):
                    maybe_start = _parse_int_field(parts[idx + 2])
                    if maybe_start is not None:
                        entry["unp_start"] = maybe_start
                if entry["unp_end"] is None and idx + 3 < len(parts):
                    maybe_end = _parse_int_field(parts[idx + 3])
                    if maybe_end is not None:
                        entry["unp_end"] = maybe_end
            except ValueError:
                pass
    if entry["pdb_start"] is None or entry["pdb_end"] is None:
        parts = line.split()
        if len(parts) >= 5:
            pdb_start = _parse_int_field(parts[3])
            pdb_end = _parse_int_field(parts[4])
            if pdb_start is not None:
                entry["pdb_start"] = entry["pdb_start"] or pdb_start
            if pdb_end is not None:
                entry["pdb_end"] = entry["pdb_end"] or pdb_end
    return entry


def parse_dbref_segments(pdb_path: Path) -> Tuple[Dict[str, Dict[str, Any]], int]:
    mapping: Dict[str, Dict[str, Any]] = {}
    dbref_count = 0
    pending_dbref1: Dict[str, Dict[str, Any]] = {}

    with pdb_path.open("r") as f:
        for line in f:
            if line.startswith(("DBREF ", "DBREF1", "DBREF2")):
                dbref_count += 1
            if line.startswith("DBREF "):
                entry = _parse_dbref_line(line)
                chain = entry.get("chain")
                if chain and entry.get("uniprot"):
                    mapping[chain] = _merge_entry(mapping.get(chain), entry)
            elif line.startswith("DBREF1"):
                entry = _parse_dbref_line(line)
                chain = entry.get("chain")
                if chain and entry.get("uniprot"):
                    pending_dbref1[chain] = entry
                    mapping[chain] = _merge_entry(mapping.get(chain), entry)
            elif line.startswith("DBREF2"):
                entry = _parse_dbref_line(line)
                chain = entry.get("chain")
                if not chain or not entry.get("uniprot"):
                    continue
                merged = _merge_entry(pending_dbref1.get(chain), entry)
                mapping[chain] = _merge_entry(mapping.get(chain), merged)

    return mapping, dbref_count


def resolve_chain_uniprot_segments(
    pdb_id: str,
    *,
    write_file: bool = True,
    return_info: bool = False,
    pdb_dir: Optional[Path] = None,
) -> Dict[str, Dict[str, Any]]:
    pdb_norm = pdb_id.upper()
    pdb_root = pdb_dir or REPO_ROOT / "input_pdbs"
    pdb_path = pdb_root / f"{pdb_norm}.pdb"
    if not pdb_path.exists():
        raise FileNotFoundError(f"{pdb_path} does not exist")

    mapping, dbref_count = parse_dbref_segments(pdb_path)
    info = {"dbref_count": dbref_count, "fallback_used": False, "pdb_id": pdb_norm}

    if not mapping:
        info["fallback_used"] = True
        fallback_uniprot = None
        if get_uniprot_and_ec:
            try:
                fallback_uniprot, _ = get_uniprot_and_ec(pdb_norm)
            except Exception as exc:
                info["fallback_error"] = str(exc)
        chains = get_polymer_chains(pdb_path)
        for chain_id in chains:
            mapping[chain_id] = {
                "uniprot": fallback_uniprot,
                "unp_start": None,
                "unp_end": None,
                "pdb_start": None,
                "pdb_end": None,
                "source": "rcsb_fallback",
            }
        info["fallback_uniprot"] = fallback_uniprot

    if write_file:
        out_dir = REPO_ROOT / "calibrator" / pdb_norm
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / "chain_uniprot.json"
        with out_file.open("w") as f:
            json.dump(mapping, f, indent=2)
    if return_info:
        return mapping, info  # type: ignore[return-value]
    return mapping


def select_primary_chain(
    chain_map: Dict[str, Dict[str, Any]]
) -> Tuple[Optional[str], Optional[Dict[str, Any]]]:
    if not chain_map:
        return None, None
    dbref_chains = [
        (cid, entry)
        for cid, entry in chain_map.items()
        if entry.get("uniprot") and entry.get("source") == "dbref"
    ]
    if len(dbref_chains) == 1:
        return dbref_chains[0]
    if len(dbref_chains) > 1:
        def _length(entry):
            start = entry.get("pdb_start")
            end = entry.get("pdb_end")
            if start is None or end is None:
                return -1
            return int(end) - int(start) + 1

        dbref_chains.sort(
            key=lambda pair: (-_length(pair[1]), pair[0])
        )
        return dbref_chains[0]

    fallback_chains = [
        (cid, entry) for cid, entry in chain_map.items() if entry.get("uniprot")
    ]
    if fallback_chains:
        fallback_chains.sort(key=lambda pair: (pair[0]))
        return fallback_chains[0]
    # No UniProt accession at all; return first chain deterministically
    sorted_items = sorted(chain_map.items(), key=lambda pair: pair[0])
    return sorted_items[0]


def main():
    parser = argparse.ArgumentParser(description="Resolve UniProt IDs for PDB chains.")
    parser.add_argument("--pdb", required=True, help="PDB ID (e.g. 6LU7)")
    args = parser.parse_args()

    pdb_id = args.pdb.upper()
    try:
        mapping, info = resolve_chain_uniprot_segments(
            pdb_id, write_file=True, return_info=True
        )
    except FileNotFoundError as exc:
        print(f"Error: {exc}")
        sys.exit(1)

    print(f"DBREF-like records found: {info.get('dbref_count')}")
    print(f"UNP mappings found: {len(mapping)}")
    print(f"Fallback used: {info.get('fallback_used')}")

    if not mapping:
        print("Error: Could not resolve any UniProt mappings.")
        sys.exit(1)

    out_file = REPO_ROOT / "calibrator" / pdb_id / "chain_uniprot.json"
    print(f"Wrote mapping to {out_file}")
    sys.exit(0)


if __name__ == "__main__":
    main()
