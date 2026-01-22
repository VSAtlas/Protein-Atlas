from __future__ import annotations

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import post_docking.rescoring.rescoring_scorch as rescoring_scorch


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def test_decoy_prefix_discovery(tmp_path: Path) -> None:
    previous = rescoring_scorch.DECOY_PREFIX_VALUE
    try:
        rescoring_scorch._set_decoy_prefix("fda_dud")
        run_id = "unitrun"
        pdb_id = "1ABC"
        variant = "HOLO"
        ph = "pH7_0"
        run_combo_root = tmp_path / "docked" / run_id / pdb_id / variant / ph
        _ensure_dir(run_combo_root / "fda_dud_stage3")
        _ensure_dir(run_combo_root / "gnina_fda_dud_stage3")

        post_combo_root = tmp_path / "post_docked" / run_id / pdb_id / variant / ph
        _ensure_dir(post_combo_root / "fda_dud_ledock_pdbqt")

        vina_candidates = rescoring_scorch.stage_dir_candidates("vina", "dud")
        assert vina_candidates[0] == "fda_dud_stage3"

        gnina_candidates = rescoring_scorch.stage_dir_candidates(
            "gnina", "dud", run_combo_root
        )
        assert "gnina_fda_dud_stage3" in gnina_candidates

        combos = rescoring_scorch.discover_combos(
            tmp_path / "docked" / run_id, tmp_path / "post_docked" / run_id
        )
        assert (pdb_id, variant, ph) in combos
    finally:
        rescoring_scorch._set_decoy_prefix(previous)
