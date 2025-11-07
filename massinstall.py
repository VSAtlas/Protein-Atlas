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
    pdb_ids = [ "4XUF", "1T46", "3JYX", "3C4F", "3DKF", "3CQW", "3D0E", "1Q5K", "3PY3", "3M2W", "2SRC", "1QCF", "2J0J", "2AEV", "3VO3", "4GRL", "3D4Q", "3KMR", "3ADX", "4LXZ", "3SFF", "1UY6", "4UND", "5IKR", "2AM9", "3ERT", "1A28", "3G6U", "2AA2", "2ITY","2HZI", "2OF2", "1AQ1", "6JON", "5TQY", "4JBV", "6GQ7", "4WAF", "4YHJ", "3DWW","4JVG", "6I83", "6I84", "6G76", "4L00", "3KDP", "4N9Y", "6GQ9", "6GQ6", "3LVP","5I35", "5NTS", "4ZDU", "3O96", "6I82", "6XBX", "3E9X", "4PYP", "3PP0", "5W8L","6EG8", "6EKF", "5N2R", "3SMQ", "5LS6", "6FZH", "4PAX", "6C1X", "4MCH", "6NNU","5ZLU", "6VQN", "5X03", "4FNY", "3TYL", "5W3J", "6KKN", "6IGK", "5I35", "6C6Y","6E6E", "6F6T", "6KRM", "6LFE", "6G2H", "5WRT", "6M99", "5N21", "5VC9", "4UMW"
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
