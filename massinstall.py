import os, requests
from pathlib import Path
from input_and_export_functions import load_inputs, validate_config

cfg = load_inputs(); validate_config(cfg)
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
    pdb_ids = [ "6JQR", "4RT7", "5I96", "6ADQ", "6U4J", "6O0K", "5L7I", "1T46", "6WTN", "1IEP", "1OPJ", "3OG7", "5L2I", "4U5J", "2XP2", "4XV2", "1M17", "4AG8", "2GQG", "3CS9", "3OXZ", "5L7D", "3LXK", "4XUF", "2HYY", "4I4E", "3QX3", "2RGC", "4ASD", "3DZY", "5VY4", "3WZE", "3WZD", "4AGC", "4TWP", "5MO4", "4R7H", "6GQO", "5C7X", "2ITN", "3ZOS", "2ZGQ", "2WGJ", "3G0E", "4U2P", "6O0L", "5VA0", "3ZBF", "2E2B", "5TQH" ]

for pdb_id in pdb_ids:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    r = requests.get(url, timeout=60)
    if r.status_code == 200:
        with open(input_dir / f"{pdb_id}.pdb", "wb") as f:
            f.write(r.content)
        print(f"Downloaded {pdb_id}")
    else:
        print(f"Failed to download {pdb_id}: {r.status_code}")
