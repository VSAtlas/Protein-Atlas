import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEEPCOY_DIR = REPO_ROOT / "DeepCoy_duds"
if str(DEEPCOY_DIR) not in sys.path:
    sys.path.insert(0, str(DEEPCOY_DIR))

from generate_dud_library import validate_and_filter_actives  # noqa: E402


def _detail(value, units, relation="="):
    return {
        "standard_value": value,
        "standard_units": units,
        "standard_relation": relation,
    }


def test_potency_filter_respects_cutoff_and_unknown_handling():
    strong = "c1ccccc1"
    weak = "CCCCCCCC"
    unknown = "CCCCC"
    greater = "CCCOC"
    actives_list = [
        (strong, "chembl"),
        (weak, "chembl"),
        (unknown, "chembl"),
        (greater, "chembl"),
    ]
    smiles_to_sources = {s: {src} for s, src in actives_list}
    provenance_details = {
        strong: [_detail(50, "nM")],
        weak: [_detail(10, "\u00b5M")],
        greater: [_detail(5, "nM", relation=">")],
    }

    kept, _, _ = validate_and_filter_actives(
        actives_list,
        smiles_to_sources=smiles_to_sources,
        provenance_details=provenance_details,
        potency_cutoff_nm=1000.0,
        potency_keep_unknown=True,
    )
    assert set(kept) == {strong, unknown, greater}
    assert weak not in kept

    kept_strict, _, _ = validate_and_filter_actives(
        actives_list,
        smiles_to_sources=smiles_to_sources,
        provenance_details=provenance_details,
        potency_cutoff_nm=1000.0,
        potency_keep_unknown=False,
    )
    assert set(kept_strict) == {strong}
