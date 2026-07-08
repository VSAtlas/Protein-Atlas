from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace

from docking import docking_ph_subruns


def test_resolve_ph_tags_string_false_uses_base_receptor(tmp_path: Path, monkeypatch) -> None:
    prepped = tmp_path / "prepped_ligands"
    (prepped / "fda_library").mkdir(parents=True)
    cfg = {
        "OVERALL_DIR": str(tmp_path),
        "PREPPED_LIGANDS_DIR": str(prepped),
        "LIBRARY_SUBDIR_DEFAULT": "fda_library",
        "PH_ENSEMBLE": "false",
        "TEST_MODE_ENABLE": "off",
    }
    paths = SimpleNamespace(pdb_id="TEST")

    monkeypatch.setattr(
        docking_ph_subruns,
        "init_ph_tags_and_manifest",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("should not init PH tags")),
    )

    ph_tags, ph_root, _override, library = docking_ph_subruns.resolve_ph_tags_and_root(
        cfg=cfg,
        paths=paths,  # type: ignore[arg-type]
        logger=logging.getLogger("test"),
        run_mode="fda",
        variant_token=None,
        legacy_mode=True,
    )

    assert ph_tags == [None]
    assert ph_root == prepped / "fda_library"
    assert library == "fda_library"
