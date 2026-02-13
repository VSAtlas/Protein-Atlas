import subprocess
import json
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
SCRIPT_PATH = REPO_ROOT / "calibrator" / "uniprot_resolver.py"

def test_uniprot_resolver_6lu7():
    pdb_id = "6LU7"
    out_dir = REPO_ROOT / "calibrator" / pdb_id
    out_file = out_dir / "chain_uniprot.json"
    
    # Cleanup before test
    if out_file.exists():
        out_file.unlink()

    # Run script
    cmd = [sys.executable, str(SCRIPT_PATH), "--pdb", pdb_id]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO_ROOT)
    
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
    
    assert result.returncode == 0, f"Script failed with return code {result.returncode}"
    assert out_file.exists(), "Output JSON file was not created"
    
    with open(out_file, 'r') as f:
        data = json.load(f)
        
    assert data.get("A") == "P0DTD1", "Chain A not mapped to P0DTD1"
    # Chain C is the peptide inhibitor, should not be mapped to P0DTD1 via DBREF
    assert "C" not in data, "Chain C should not be mapped (it is not UNP in DBREF)"
    
    # Check stdout for log messages
    assert "DBREF-like records found" in result.stdout
    assert "UNP mappings found" in result.stdout
    assert "Fallback used: False" in result.stdout

    # Cleanup
    if out_file.exists():
        out_file.unlink()
        try:
            out_dir.rmdir()
        except OSError:
            pass
