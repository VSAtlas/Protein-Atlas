from __future__ import annotations

import json
import logging
from pathlib import Path

from cli import distributed_chunk_planner as chunk_planner


class _DistCtx:
    enabled = True
    mode = "slurm_array"
    task_id = 2
    task_count = 4
    task_min_id = 1
    leader_task_id = 1
    is_leader = False


def _write_manifest(root: Path) -> None:
    root.mkdir(parents=True, exist_ok=True)
    payload = {"entries": {}, "filenames": {}}
    (root / "_manifest.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def test_manifest_preflight_rebuilds_only_broken_roots(
    tmp_path: Path, monkeypatch
) -> None:
    run_id = "pref_run"
    good = tmp_path / "prepped_ligands" / "good"
    bad = tmp_path / "prepped_ligands" / "bad"
    _write_manifest(good)
    bad.mkdir(parents=True, exist_ok=True)

    captured: dict[str, list[Path]] = {"roots": []}

    monkeypatch.setattr(
        chunk_planner,
        "_collect_combo_library_roots",
        lambda cfg, combo_items, run_tokens, logger: [good, bad],
    )

    def _fake_rebuild(cfg, *, roots, logger):
        captured["roots"] = list(roots)
        return {
            "roots_total": len(roots),
            "built": len(roots),
            "rebuilt": 0,
            "unchanged": 0,
            "missing_root": 0,
            "failed": 0,
        }

    monkeypatch.setattr(chunk_planner, "_rebuild_library_manifests", _fake_rebuild)

    cfg = {"OVERALL_DIR": str(tmp_path)}
    chunk_planner._ensure_distributed_manifest_preflight(
        cfg,
        run_id=run_id,
        dist_ctx=_DistCtx(),
        variant_label="HOLO",
        combo_items=[("BNJS.pdb", "pH7_0")],
        run_tokens=["dud"],
    )

    assert captured["roots"] == [bad]
    state_path = (
        tmp_path
        / "manifests"
        / run_id
        / "distributed"
        / "manifest_preflight_holo.json"
    )
    payload = json.loads(state_path.read_text(encoding="utf-8"))
    assert payload.get("status") == "ready"
    summary = payload.get("summary") or {}
    assert int(summary.get("roots_total", 0)) == 2
    assert int(summary.get("healthy", 0)) == 1
    assert int(summary.get("broken_total", 0)) == 1


def test_manifest_preflight_wait_logs_on_non_owner(tmp_path: Path, monkeypatch, caplog) -> None:
    run_id = "pref_wait"
    cfg = {"OVERALL_DIR": str(tmp_path)}
    dist_dir = tmp_path / "manifests" / run_id / "distributed"
    state_path = dist_dir / "manifest_preflight_holo.json"
    lock_path = state_path.with_suffix(state_path.suffix + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    lock_path.write_text("", encoding="utf-8")

    monkeypatch.setattr(chunk_planner, "_COMBO_PREP_WAIT_TIMEOUT_SEC", 0.25)
    monkeypatch.setattr(chunk_planner, "_COMBO_PREP_WAIT_POLL_SEC", 0.05)
    monkeypatch.setattr(chunk_planner, "_MANIFEST_PREFLIGHT_WAIT_LOG_SEC", 0.1)
    monkeypatch.setattr(
        chunk_planner,
        "_collect_combo_library_roots",
        lambda cfg, combo_items, run_tokens, logger: [],
    )
    monkeypatch.setattr(
        chunk_planner,
        "_rebuild_library_manifests",
        lambda cfg, roots, logger: {
            "roots_total": 0,
            "built": 0,
            "rebuilt": 0,
            "unchanged": 0,
            "missing_root": 0,
            "failed": 0,
        },
    )

    caplog.set_level(logging.INFO, logger="distributed.chunk")
    chunk_planner._ensure_distributed_manifest_preflight(
        cfg,
        run_id=run_id,
        dist_ctx=_DistCtx(),
        variant_label="HOLO",
        combo_items=[("BNJS.pdb", "pH7_0")],
        run_tokens=["dud"],
    )
    assert any(
        "[distributed.chunk.manifest-preflight.wait]" in rec.getMessage()
        for rec in caplog.records
    )
