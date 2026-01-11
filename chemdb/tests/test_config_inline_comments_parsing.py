from pathlib import Path

from input_and_export_functions import load_config


def test_inline_comments_are_stripped_and_typed(tmp_path: Path, monkeypatch) -> None:
    for key in (
        "DEEPCOY_CHUNK_SIZE",
        "DEEPCOY_SEED_PER_CHUNK",
        "DEEPCOY_DECOYS_PER_ACTIVE",
        "EARLY_RECENTER_RATIO",
    ):
        monkeypatch.delenv(key, raising=False)

    cfg_path = tmp_path / "config.txt"
    cfg_path.write_text(
        "\n".join(
            [
                "DEEPCOY_CHUNK_SIZE=1 # chunk size",
                "DEEPCOY_SEED_PER_CHUNK=true # seed",
                "DEEPCOY_DECOYS_PER_ACTIVE=1 # faster",
                "EARLY_RECENTER_RATIO=0.70 # ratio",
                'SOME_PATH="C:/path/with#hash" # comment',
                "",
            ]
        )
    )

    cfg = load_config(config_path=str(cfg_path), base_dir=tmp_path)

    assert cfg["DEEPCOY_CHUNK_SIZE"] == 1
    assert cfg["DEEPCOY_SEED_PER_CHUNK"] is True
    assert cfg["DEEPCOY_DECOYS_PER_ACTIVE"] == 1
    assert isinstance(cfg["EARLY_RECENTER_RATIO"], float)
    assert cfg["EARLY_RECENTER_RATIO"] == 0.70
    assert "with#hash" in cfg["SOME_PATH"]
    assert cfg["SOME_PATH"].strip('"') == "C:/path/with#hash"
