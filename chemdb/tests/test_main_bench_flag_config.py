import sys
import os
import pytest
from unittest.mock import MagicMock

# Ensure we can import main from root
sys.path.append(os.getcwd())

import main

def test_bench_enabled():
    assert main._bench_enabled(["main.py", "-bench"]) is True
    assert main._bench_enabled(["main.py", "--bench"]) is True
    assert main._bench_enabled(["main.py", "-b"]) is False
    assert main._bench_enabled(["main.py", "bench"]) is False
    assert main._bench_enabled(["main.py"]) is False

def test_apply_bench_overrides():
    # Setup minimal config
    cfg = {
        "TEST_MODE_ENABLE": "off",
        "USE_GNINA": False,
        "APO_HOLO_MODE": "apo",
        "PH_ENSEMBLE": False,
        "TEST_LIBRARY_MAP": {}
    }
    
    # Apply overrides
    main._apply_bench_overrides(cfg)
    
    # Check overrides
    assert cfg["TEST_MODE_ENABLE"] == "dud"
    assert cfg["USE_GNINA"] is True
    assert cfg["USE_LEDOCK"] is True
    assert cfg["USE_DOCK6"] is True
    assert cfg["USE_SCORCH"] is True
    assert cfg["APO_HOLO_MODE"] == "holo"
    assert cfg["PH_ENSEMBLE"] is True
    
    # Check proteins
    for pid in ["bNJS", "bOJG", "bNNQ"]:
        assert pid in cfg["SPECIFIED_PROTEINS"]
        
    # Check mapping
    assert cfg["TEST_LIBRARY_MAP"]["bNJS"] == "bench_pur2"
    assert cfg["TEST_LIBRARY_MAP"]["bOJG"] == "bench_mk01"
    assert cfg["TEST_LIBRARY_MAP"]["bNNQ"] == "bench_fabp4"

def test_apply_bench_overrides_merges_map():
    # Setup config with existing map
    cfg = {
        "TEST_LIBRARY_MAP": {"EXIST": "lib_exist"}
    }
    
    main._apply_bench_overrides(cfg)
    
    # Check new entries
    assert cfg["TEST_LIBRARY_MAP"]["bNJS"] == "bench_pur2"
    # Check existing entry preserved
    assert cfg["TEST_LIBRARY_MAP"]["EXIST"] == "lib_exist"
