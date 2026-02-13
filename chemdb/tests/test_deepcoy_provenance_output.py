import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEEPCOY_DIR = REPO_ROOT / "DeepCoy_duds"
if str(DEEPCOY_DIR) not in sys.path:
    sys.path.insert(0, str(DEEPCOY_DIR))

from DeepCoy_duds.generate_dud_library import write_provenance_files  # noqa: E402


def test_provenance_files_include_sources_and_urls(tmp_path):
    label = "TEST"
    labeled = [("METAB_1", "C"), ("", "CC")]
    sources = {"C": set(["chembl:target=T1"]), "CC": set(["iuphar:ligand=5"])}
    details = {
        "C": [
            {
                "source": "chembl",
                "target_chembl_id": "CHEMBL1",
                "molecule_chembl_id": "CHEMBL2",
                "label": "chembl:target=T1",
            }
        ],
        "CC": [{"source": "iuphar", "ligandId": 5, "label": "iuphar:ligand=5"}],
    }

    tsv_path, json_path = write_provenance_files(
        label, tmp_path, labeled, sources, details
    )

    assert tsv_path and tsv_path.exists()
    lines = tsv_path.read_text().splitlines()
    assert lines[0] == "metab_id\tsmiles\tsources_joined"
    assert "chembl:target=T1" in lines[1]

    payload = json.loads(json_path.read_text())
    assert payload and payload[0]["browse_urls"]
