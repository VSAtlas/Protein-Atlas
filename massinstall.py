import requests
from pathlib import Path
from input_and_export_functions import load_inputs, validate_config

cfg = load_inputs()
validate_config(cfg)
input_dir = Path(cfg["INPUT_DIR"])
input_dir.mkdir(parents=True, exist_ok=True)

# Optional: load PDB IDs from a file if provided
pdb_ids = None
pdb_list_path = str(cfg.get("PDB_ID_LIST", "")).strip()
if pdb_list_path:
    try:
        with open(pdb_list_path, "r", encoding="utf-8") as f:
            pdb_ids = [ln.strip().upper() for ln in f if ln.strip()]
    except Exception:
        pdb_ids = None

# Fallback to the baked-in list if no file provided or read failed
if not pdb_ids:
    pdb_ids = [
        "3EML",
        "2HZI",
        "3BKL",
        "1E66",
        "2E1W",
        "2OI0",
        "2VT4",
        "3NY8",
        "3CQW",
        "3D0E",
        "2HV5",
        "1L2S",
        "2AM9",
        "1S3B",
        "3L5D",
        "3D4Q",
        "1BCD",
        "2CNK",
        "1H00",
        "3BWM",
        "1R9O",
        "3NXU",
        "3KRJ",
        "3ODU",
        "1LRU",
        "3FRJ",
        "2I78",
        "3PBL",
        "3NXO",
        "2RGP",
        "1SJ0",
        "2FSZ",
        "3KL6",
        "1W7X",
        "2NNQ",
        "3BZ3",
        "3C4F",
        "1J4H",
        "3E37",
        "1ZW5",
        "3BQD",
        "2V3F",
        "3KGC",
        "1VSO",
        "3MAX",
        "3F07",
        "3NF7",
        "1XL2",
        "3LAN",
        "3CCW",
        "1UYG",
        "3F9M",
        "2OJ9",
        "2H7L",
        "2ICA",
        "3LPB",
        "3CJO",
        "3G0E",
        "2B8T",
        "2I0E",
        "2OF2",
        "3CHP",
        "3M2W",
        "2AA2",
        "3LQ8",
        "2OJG",
        "2ZDT",
        "2QD9",
        "830C",
        "3EQH",
        "1QW6",
        "1B9V",
        "1KVO",
        "3L3M",
        "1UDT",
        "2OYU",
        "3LN1",
        "2OWB",
        "3BGS",
        "2P54",
        "2ZNP",
        "2GTK",
        "3KBA",
        "2AZR",
        "1NJS",
        "1C8K",
        "1D3G",
        "3G6Z",
        "2ETR",
        "1MV9",
        "1LI4",
        "3EL8",
        "3HMM",
        "1Q4X",
        "1YPE",
        "2AYW",
        "2ZEC",
        "1SYN",
        "1SQT",
        "2P2I",
        "3BIZ",
        "3HL5",
    ]

for pdb_id in pdb_ids:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    r = requests.get(url, timeout=60)
    if r.status_code == 200:
        with open(input_dir / f"{pdb_id}.pdb", "wb") as f:
            f.write(r.content)
        print(f"Downloaded {pdb_id}")
    else:
        print(f"Failed to download {pdb_id}: {r.status_code}")
