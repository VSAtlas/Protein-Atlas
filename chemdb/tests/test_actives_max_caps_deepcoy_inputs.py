from pathlib import Path
import sys

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
DEEPCOY_ROOT = REPO_ROOT / "DeepCoy_duds"
if str(DEEPCOY_ROOT) not in sys.path:
    sys.path.insert(0, str(DEEPCOY_ROOT))

from DeepCoy_duds import generate_dud_library as deepcoy_lib  # type: ignore  # noqa: E402


@pytest.mark.parametrize("actives_max", [3])
def test_actives_max_limits_actives_smi(tmp_path: Path, monkeypatch, actives_max: int):
    actives = [
        "C",
        "CC",
        "CCC",
        "CCCC",
        "CCCCC",
        "CCCCCC",
        "CCCCCCC",
        "CCCCCCCC",
    ]

    real_load_config = deepcoy_lib.load_config

    def _fake_load_config(config_path="config.txt", base_dir=None):
        cfg = real_load_config(config_path=config_path, base_dir=base_dir)
        cfg["ACTIVES_MAX"] = actives_max
        cfg["TEST_MODE_ENABLE"] = "dud"
        cfg["DEEPCOY_DECOYS_PER_ACTIVE"] = 1
        return cfg

    def _fake_validate_and_filter(
        _all_actives, smiles_to_sources=None, provenance_details=None, **kwargs
    ):
        sources = {smi: {"test"} for smi in actives}
        details = {smi: [{"source": "test"}] for smi in actives}
        return list(actives), sources, details

    def _fake_convert_smi_to_sdf(smi_path, sdf_path, require_output=False):
        sdf_path.parent.mkdir(parents=True, exist_ok=True)
        sdf_path.write_text("")

    monkeypatch.setattr(deepcoy_lib, "load_config", _fake_load_config)
    monkeypatch.setattr(
        deepcoy_lib, "validate_and_filter_actives", _fake_validate_and_filter
    )
    monkeypatch.setattr(deepcoy_lib, "convert_smi_to_sdf", _fake_convert_smi_to_sdf)
    monkeypatch.setattr(
        deepcoy_lib, "run_deepcoy_workflow", lambda *args, **kwargs: None
    )
    monkeypatch.setattr(
        deepcoy_lib, "query_external_sources", lambda *args, **kwargs: ([], {}, {})
    )

    output_root = tmp_path / "out"
    artifact_dir = tmp_path / "artifact"
    argv = [
        "generate_dud_library.py",
        "--pdb",
        "1B9V",
        "--input-pdb-dir",
        str(REPO_ROOT / "input_pdbs"),
        "--out-root",
        str(output_root),
        "--artifact-dir",
        str(artifact_dir),
        "--offline",
        "--skip-sdf",
        "--no-run-deepcoy",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    exit_code = deepcoy_lib.main()
    assert exit_code == 0

    actives_files = list(output_root.glob("*/*_actives.smi"))
    assert len(actives_files) == 1, f"expected one actives file, found: {actives_files}"
    actives_path = actives_files[0]
    lines = [
        line.strip()
        for line in actives_path.read_text().splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]
    assert len(lines) == actives_max
