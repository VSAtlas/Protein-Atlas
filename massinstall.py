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
    pdb_ids = [ "5X02", "4RT7", "6JQR", "4XUF", "8XB1", "6IO0", "5DE1", "5LGE", "8HB9", "6VEI", "6VFZ", "6ADI", "8I7Z", "7W7Y", "7W7X", "3CS9", "7RN6", "6BBV", "3UGC", "6WTN", "5TTU", "7Q6H", "6HZV", "1PKG", "6GQJ", "6KLA", "8RRQ", "4YJR", "6HM6", "3KMM", "1LKL", "1QPE", "6O6G", "4OQ5", "5FDR", "5JSB", "4IEH", "4AQ3", "6QGK", "8VXM", "4F3I", "8B98", "7YL2", "6CJ2", "6ZB1", "4UYG", "3S91", "6C7Q", "7USJ", "8ZM8", "6NQM", "4UVB", "5LHG", "6X9K", "7SFC", "6X9J", "6NUQ", "7UC6", "7UC7", "6SM8", "6GGH", "5B7V", "4QQC", "4F64", "6NT2", "9EYU", "4HSG", "2EUF", "1XO2", "6P8G", "5A14", "1W0X", "6XD3", "6Z45", "6W9E", "3MVH", "8Q61", "3ML9", "3B7O", "3A4O", "3APD", "8TSC", "2FQQ", "1RWM", "6CKZ", "3KJQ", "1QTN", "1NW9", "1TUP", "8CRC", "4J52", "4B6L", "4ZS0", "2WTV", "4AF3", "6SLG", "5NHJ", "1IAN", "1ZYJ", "2ZDT", "3V6R", "4G9C", "3TV4", "5WG6", "4MI0", "4LWH", "2FWZ", "5CF0", "5ZTY", "4OO9", "5LPK", "5NU5", "8F1G", "8G3E", "4FOB", "9GBE", "3AOX", "1BZH", "1G7F", "1KAV", "3V8S", "3TV7", "4ZS0", "2WTV", "4AF3", "6SLG", "5NHJ", "1IAN", "1ZYJ", "2ZDT", "3V6R", "4G9C", "3TV4", "5WG6", "4MI0", "4LWH", "2FWZ", "5CF0", "5ZTY", "4OO9", "5LPK", "5NU5", "8F1G", "8G3E", "4FOB", "9GBE", "3AOX", "1BZH", "1G7F", "1KAV", "3V8S", "3TV7", "4Y72", "6GU4"
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
