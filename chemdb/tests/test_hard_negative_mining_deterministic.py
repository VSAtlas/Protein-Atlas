from __future__ import annotations

import pandas as pd

from ml.hard_negatives import merge_hard_negatives


def _row(smiles: str, active: int, scaffold: str, target: str, pocket: str) -> dict[str, object]:
    return {
        "lig_smiles": smiles,
        "active": active,
        "murcko_scaffold": scaffold,
        "ex_rec_pdb": target,
        "pocket": pocket,
    }


def test_hard_negative_mining_deterministic():
    df = pd.DataFrame(
        [
            _row("c1ccccc1O", 1, "A", "TRN1", "PK1"),
            _row("c1ccncc1", 1, "B", "TRN1", "PK1"),
            _row("CCO", 0, "C", "TRN1", "PK1"),
            _row("CCN", 0, "D", "TRN1", "PK1"),
            _row("CCC", 0, "E", "TRN1", "PK1"),
            _row("c1ccccc1N", 1, "F", "TRN2", "PK2"),
            _row("CCCO", 0, "G", "TRN2", "PK2"),
            _row("CCCC", 0, "H", "TRN2", "PK2"),
        ]
    )

    aug1 = merge_hard_negatives(
        df,
        policy="hybrid",
        ratio=1.0,
        per_active_k=2,
        within_group=True,
        group_key="target",
        max_candidates=100,
        seed=123,
    )
    aug2 = merge_hard_negatives(
        df,
        policy="hybrid",
        ratio=1.0,
        per_active_k=2,
        within_group=True,
        group_key="target",
        max_candidates=100,
        seed=123,
    )

    mined1 = aug1.loc[aug1["_hard_negative_augmented"] == 1, ["lig_smiles", "ex_rec_pdb"]].copy()
    mined2 = aug2.loc[aug2["_hard_negative_augmented"] == 1, ["lig_smiles", "ex_rec_pdb"]].copy()
    assert mined1.reset_index(drop=True).equals(mined2.reset_index(drop=True))
    assert set(mined1["ex_rec_pdb"].unique().tolist()).issubset({"TRN1", "TRN2"})
