import sys
import os

# Ensure we can import main from root
sys.path.append(os.getcwd())

import main


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
    }

    main._apply_bench2_overrides(cfg)

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
