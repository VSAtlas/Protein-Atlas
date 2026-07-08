from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from docking.engine_config_emit import _write_bytes_atomic, _write_manifest_payload


def test_write_bytes_atomic_allows_concurrent_same_target(tmp_path: Path) -> None:
    target = tmp_path / "ligand_stage1.txt"
    payloads = [f"payload-{idx}\n".encode("utf-8") for idx in range(32)]

    with ThreadPoolExecutor(max_workers=16) as pool:
        list(pool.map(lambda payload: _write_bytes_atomic(target, payload), payloads))

    assert target.read_bytes() in payloads
    assert not [path for path in tmp_path.iterdir() if path.name.endswith(".part")]


def test_write_manifest_payload_merges_concurrent_entries(tmp_path: Path) -> None:
    manifest = tmp_path / "vina.json"

    def write_entry(idx: int) -> None:
        _write_manifest_payload(
            manifest_path=manifest,
            run_id="run1",
            pdb_id="P123",
            stage_name="stage1",
            variant_token="HOLO",
            ph_label="pH7_0",
            legacy_mode=False,
            lig_base=f"lig_{idx:02d}",
            cfg_path=tmp_path / f"lig_{idx:02d}.txt",
            out_path=tmp_path / f"lig_{idx:02d}.pdbqt",
            receptor_for_config="receptor.pdbqt",
            flex_pdbqt=None,
        )

    with ThreadPoolExecutor(max_workers=16) as pool:
        list(pool.map(write_entry, range(32)))

    data = json.loads(manifest.read_text(encoding="utf-8"))
    ligands = {entry["ligand"] for entry in data["entries"]}
    assert ligands == {f"lig_{idx:02d}" for idx in range(32)}
