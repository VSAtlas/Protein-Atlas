import sys
import os

# Ensure we can import main from root
sys.path.append(os.getcwd())

from cli import run_profiles


def test_apply_bench2_overrides_vina_scorch_only():
    cfg = {
        "TEST_MODE_ENABLE": "off",
        "USE_GNINA": True,
        "USE_LEDOCK": True,
        "USE_DOCK6": True,
        "USE_SCORCH": False,
        "APO_HOLO_MODE": "apo",
        "PH_ENSEMBLE": False,
        "TEST_LIBRARY_MAP": {"EXIST": "lib_exist"},
        "DISCOVERY_SELECTION_PCTS": "1.0,0.25,0.10,0.05,0.02",
        "POLYPHARM_SELECTION_PCTS": "1.0,0.15,0.10",
        "SCORCH_TOP_FRACTION": 0.15,
        "CNN_TOP_FRACTION": 0.10,
        "PERCENT_SEARCHED": 0.25,
    }

    run_profiles._apply_bench2_overrides(cfg)

    assert cfg["TEST_MODE_ENABLE"] == "dud"
    assert cfg["APO_HOLO_MODE"] == "holo"
    assert cfg["PH_ENSEMBLE"] is True
    assert cfg["USE_GNINA"] is False
    assert cfg["USE_LEDOCK"] is False
    assert cfg["USE_DOCK6"] is False
    assert cfg["USE_SCORCH"] is True

    for pid in ["bNJS", "bOJG", "bNNQ"]:
        assert pid in cfg["SPECIFIED_PROTEINS"]

    assert cfg["TEST_LIBRARY_MAP"]["bNJS"] == "bench_pur2"
    assert cfg["TEST_LIBRARY_MAP"]["bOJG"] == "bench_mk01"
    assert cfg["TEST_LIBRARY_MAP"]["bNNQ"] == "bench_fabp4"
    assert cfg["TEST_LIBRARY_MAP"]["EXIST"] == "lib_exist"
    assert cfg["DISCOVERY_SELECTION_PCTS"] == "1.0,0.25,0.10,0.05,0.02"
    assert cfg["POLYPHARM_SELECTION_PCTS"] == "1.0,0.15,0.10"
    assert cfg["SCORCH_TOP_FRACTION"] == 0.15
    assert cfg["CNN_TOP_FRACTION"] == 0.10
    assert cfg["PERCENT_SEARCHED"] == 0.25
