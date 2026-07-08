from __future__ import annotations

import os
import time
from pathlib import Path

import main
from cli.distributed_context import (
    DistributedRunContext,
    read_combo_prep_state,
    try_begin_combo_prep,
    write_combo_prep_state,
)


def _ctx(run_id: str, *, task_id: int, task_count: int = 2) -> DistributedRunContext:
    return DistributedRunContext(
        mode="slurm_array",
        enabled=True,
        run_id=run_id,
        task_id=task_id,
        task_count=task_count,
        task_min_id=1,
        leader_task_id=1,
        barrier_timeout_sec=10.0,
        barrier_poll_sec=0.1,
    )


def test_pick_next_combo_index_skips_inflight_pdb() -> None:
    pending = [
        ("BOJG.pdb", "pH7_2"),
        ("BOJG.pdb", "pH7_7"),
        ("BNJS.pdb", "pH7_2"),
    ]
    idx = main._pick_next_combo_index(pending, inflight_pdb_ids={"BOJG"})
    assert idx == 2


def test_pick_next_combo_index_hybrid_unlock_allows_ready_scope() -> None:
    pending = [
        ("BOJG.pdb", "pH7_7"),
        ("BNJS.pdb", "pH7_0"),
    ]
    scope_key = main._local_combo_scope_key(
        pdb_file="BOJG.pdb",
        variant_label="holo",
    )
    idx = main._pick_next_combo_index(
        pending,
        inflight_scope_counts={scope_key: 1},
        prep_ready_scopes={scope_key},
        allow_same_scope_when_ready=True,
        variant_label="holo",
    )
    assert idx == 0


def test_pick_next_combo_index_hybrid_unlock_blocks_unready_scope() -> None:
    pending = [
        ("BOJG.pdb", "pH7_7"),
        ("BNJS.pdb", "pH7_0"),
    ]
    scope_key = main._local_combo_scope_key(
        pdb_file="BOJG.pdb",
        variant_label="holo",
    )
    idx = main._pick_next_combo_index(
        pending,
        inflight_scope_counts={scope_key: 1},
        prep_ready_scopes=set(),
        allow_same_scope_when_ready=True,
        variant_label="holo",
    )
    assert idx == 1


def test_combo_prep_owner_claim_and_ready_state(tmp_path: Path) -> None:
    cfg = {"OVERALL_DIR": str(tmp_path)}
    combo_key = "combo_owner_demo"
    owner = _ctx("run_combo_owner", task_id=1)
    follower = _ctx("run_combo_owner", task_id=2)

    assert try_begin_combo_prep(owner, cfg, combo_key=combo_key, lease_sec=600.0)
    assert not try_begin_combo_prep(follower, cfg, combo_key=combo_key, lease_sec=600.0)

    state_pre = read_combo_prep_state(owner, cfg, combo_key=combo_key)
    assert isinstance(state_pre, dict)
    assert state_pre.get("status") == "preparing"
    assert int(state_pre.get("task_id") or -1) == 1

    write_combo_prep_state(
        owner,
        cfg,
        combo_key=combo_key,
        status="ready",
        payload={"note": "prep complete"},
    )
    state_ready = read_combo_prep_state(follower, cfg, combo_key=combo_key)
    assert isinstance(state_ready, dict)
    assert state_ready.get("status") == "ready"
    assert state_ready.get("note") == "prep complete"
    assert not try_begin_combo_prep(follower, cfg, combo_key=combo_key, lease_sec=600.0)


def test_scope_prep_key_collapses_ph_dimension() -> None:
    combo_key_low = main._combo_prep_key(
        run_id="r1",
        variant_label="holo",
        pdb_file="BOJG.pdb",
        ph_tag="pH7_2",
    )
    combo_key_high = main._combo_prep_key(
        run_id="r1",
        variant_label="holo",
        pdb_file="BOJG.pdb",
        ph_tag="pH7_7",
    )
    scope_key = main._scope_prep_key(
        run_id="r1",
        variant_label="holo",
        pdb_file="BOJG.pdb",
    )
    assert combo_key_low != combo_key_high
    assert len(scope_key) == len(combo_key_low)


def test_combo_receptor_ready_respects_min_mtime(
    monkeypatch, tmp_path: Path
) -> None:
    receptor = tmp_path / "BOJG_receptor.pdbqt"
    receptor.write_text("ATOM\n", encoding="utf-8")
    stale_at = time.time() - 30.0
    os.utime(receptor, (stale_at, stale_at))

    class _DummyPaths:
        pdb_id = "BOJG"

    monkeypatch.setattr(main, "make_paths", lambda *args, **kwargs: _DummyPaths())
    monkeypatch.setattr(main, "receptor_file", lambda *args, **kwargs: receptor)

    assert main._combo_receptor_ready(
        {}, pdb_file="BOJG.pdb", variant_token=None, ph_tag="pH7_2"
    )
    assert not main._combo_receptor_ready(
        {},
        pdb_file="BOJG.pdb",
        variant_token=None,
        ph_tag="pH7_2",
        min_mtime=stale_at + 10.0,
    )

    fresh_at = time.time()
    os.utime(receptor, (fresh_at, fresh_at))
    assert main._combo_receptor_ready(
        {},
        pdb_file="BOJG.pdb",
        variant_token=None,
        ph_tag="pH7_2",
        min_mtime=fresh_at - 1.0,
    )


def test_refresh_local_prep_ready_scopes_uses_ready_floor(
    monkeypatch, tmp_path: Path
) -> None:
    receptor = tmp_path / "BOJG_receptor.pdbqt"
    receptor.write_text("ATOM\n", encoding="utf-8")
    stale_at = time.time() - 20.0
    os.utime(receptor, (stale_at, stale_at))

    class _DummyPaths:
        pdb_id = "BOJG"

    monkeypatch.setattr(main, "make_paths", lambda *args, **kwargs: _DummyPaths())
    monkeypatch.setattr(main, "receptor_file", lambda *args, **kwargs: receptor)

    scope = main._local_combo_scope_key(pdb_file="BOJG.pdb", variant_label="holo")
    ready_scopes: set[str] = set()
    ready_floor = stale_at + 5.0
    main._refresh_local_prep_ready_scopes(
        {},
        variant_token=None,
        variant_label="holo",
        combos=[("BOJG.pdb", "pH7_2")],
        ready_scopes=ready_scopes,
        scope_ready_after={scope: ready_floor},
        default_ready_after=ready_floor,
    )
    assert scope not in ready_scopes

    fresh_at = time.time()
    os.utime(receptor, (fresh_at, fresh_at))
    main._refresh_local_prep_ready_scopes(
        {},
        variant_token=None,
        variant_label="holo",
        combos=[("BOJG.pdb", "pH7_2")],
        ready_scopes=ready_scopes,
        scope_ready_after={scope: ready_floor},
        default_ready_after=ready_floor,
    )
    assert scope in ready_scopes
