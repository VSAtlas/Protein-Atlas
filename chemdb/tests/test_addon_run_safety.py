from __future__ import annotations

import json
from pathlib import Path

import pytest

from analysis.cli.score_reference_run_vina import (
    _persist_final_manifest,
    build_parser,
)
from analysis.ml import addon_run_safety
from analysis.ml.addon_run_safety import ActiveAtlasRun
from cli.qol.pipeline import _is_compare_addon_run


def test_require_idle_atlas_blocks_recent_manifest(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    recent = ActiveAtlasRun(
        run_id="active_run",
        status="running",
        manifest="/repo/manifests/active_run/run_manifest.yaml",
        activity_age_seconds=60.0,
        reason="recent_active_manifest",
    )
    monkeypatch.setattr(
        addon_run_safety,
        "active_atlas_runs",
        lambda _root: [recent],
    )

    with pytest.raises(RuntimeError, match="active_run.*recent_active_manifest"):
        addon_run_safety.require_idle_atlas(Path("/repo"))


def test_require_idle_atlas_accepts_only_empty_active_set(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        addon_run_safety,
        "active_atlas_runs",
        lambda _root: [],
    )

    assert addon_run_safety.require_idle_atlas(Path("/repo")) is None


def test_scoring_cli_rejects_removed_manifest_bypass(tmp_path: Path) -> None:
    pairs = tmp_path / "pairs.csv"
    pairs.write_text("pdb_id,ligand_base,pdbqt_path\n", encoding="utf-8")

    with pytest.raises(SystemExit):
        build_parser().parse_args(
            [
                "--run-id",
                "addon_run",
                "--selected-pairs",
                str(pairs),
                "--compare-run",
                "reference_run",
                "--out-dir",
                str(tmp_path / "out"),
                "--allow-stale-manifest",
            ]
        )


@pytest.mark.parametrize(
    ("schema", "dry_run"),
    [
        ("atlas.reference-addon-plan.v1", True),
        ("atlas.reference-vina-compare.v1", False),
        ("atlas.reference-addon-score.v1", False),
    ],
)
def test_final_manifest_persists_exact_returned_provenance(
    tmp_path: Path,
    schema: str,
    dry_run: bool,
) -> None:
    manifest = {"schema": schema, "outputs": {"scores": "scores.csv"}}
    returned = _persist_final_manifest(
        manifest,
        out_dir=tmp_path,
        lock_path=tmp_path / ".atlas_ml_addon.lock",
        dry_run=dry_run,
    )
    final_path = tmp_path / "reference_addon_manifest.json"
    persisted = json.loads(final_path.read_text(encoding="utf-8"))

    assert persisted == returned
    assert persisted["exclusive_lock"].endswith(".atlas_ml_addon.lock")
    assert persisted["bypassed_stale_manifests"] == []
    assert persisted["active_run_policy"]["manifest_bypass_supported"] is False
    assert persisted["active_run_policy"]["check_result"] == "passed"
    assert persisted["dry_run"] is dry_run
    assert persisted["outputs"]["final_manifest"] == str(final_path)


def test_run_router_recognizes_reference_addon_aliases() -> None:
    assert _is_compare_addon_run(
        ["--selected-pairs=pairs.csv", "--compare-runid=reference"]
    )
    assert not _is_compare_addon_run(["--selected-pairs=pairs.csv"])
    assert not _is_compare_addon_run(["--compare-runid=reference"])
