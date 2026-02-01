from __future__ import annotations

from typing import List, Dict

from chemdb import calibrator_sampling


def _synthetic_bags() -> List[Dict[str, object]]:
    scaffolds = [
        "C1CCCCC1",
        "C1CCCC1",
        "C1CCC1",
        "C1CCOCC1",
        "C1CCNCC1",
        "C1CCSCC1",
        "c1ccccc1",
        "c1ccncc1",
        "c1ncncc1",
        "c1ccc2ccccc2c1",
    ]
    rows: List[Dict[str, object]] = []
    idx = 0
    for s_idx, base in enumerate(scaffolds):
        for j in range(6):
            label = "strong" if j < 3 else "weak"
            bag_uid = f"bag_{s_idx}_{j}"
            side = "C" * (j + 1)
            smiles = f"{base}{side}"
            inchikey = f"SCF{s_idx:02d}INCHIKEYXX"
            if idx % 7 == 0:
                inchikey = None
            rows.append(
                {
                    "bag_uid": bag_uid,
                    "ligand_uid": bag_uid,
                    "canonical_smiles": smiles,
                    "inchikey": inchikey,
                    "y_binding_class": label,
                }
            )
            idx += 1
    return rows


def test_calibrator_sampling_determinism_and_sets():
    bags = _synthetic_bags()
    seed = 123
    assigned_a = calibrator_sampling.assign_holdout_sets(
        bags, seed=seed, test_fraction=0.20, remainder_eval_fraction=0.0
    )
    assigned_b = calibrator_sampling.assign_holdout_sets(
        bags, seed=seed, test_fraction=0.20, remainder_eval_fraction=0.0
    )

    map_a = {row["bag_uid"]: row["set_name"] for row in assigned_a}
    map_b = {row["bag_uid"]: row["set_name"] for row in assigned_b}
    assert map_a == map_b

    train_pool = [r for r in assigned_a if r["set_name"] == "train_pool"]
    sample_a, _ = calibrator_sampling.select_pocket_eval_sample(
        train_pool, max_n=20, policy="stratified_binding_class", seed=seed
    )
    sample_b, _ = calibrator_sampling.select_pocket_eval_sample(
        train_pool, max_n=20, policy="stratified_binding_class", seed=seed
    )
    sample_ids_a = [r["bag_uid"] for r in sample_a]
    sample_ids_b = [r["bag_uid"] for r in sample_b]
    assert sample_ids_a == sample_ids_b

    strong_pool = sum(1 for r in train_pool if r["y_binding_class"] == "strong")
    strong_sample = sum(1 for r in sample_a if r["y_binding_class"] == "strong")
    pool_ratio = strong_pool / len(train_pool)
    sample_ratio = strong_sample / len(sample_a)
    assert abs(sample_ratio - pool_ratio) <= 0.10

    sample_scaffold, _ = calibrator_sampling.select_pocket_eval_sample(
        train_pool, max_n=20, policy="stratified_scaffold", seed=seed
    )
    scaffolds = {r.get("scaffold_key") for r in sample_scaffold}
    scaffolds.discard(None)
    assert len(scaffolds) >= 8

    pocket_ids = {r["bag_uid"] for r in sample_scaffold}
    train_ids = {r["bag_uid"] for r in train_pool if r["bag_uid"] not in pocket_ids}
    test_ids = {r["bag_uid"] for r in assigned_a if r["set_name"] == "test"}
    remainder_ids = {
        r["bag_uid"] for r in assigned_a if r["set_name"] == "remainder_eval"
    }

    assert pocket_ids.isdisjoint(train_ids)
    assert pocket_ids.isdisjoint(test_ids)
    assert pocket_ids.isdisjoint(remainder_ids)
    assert train_ids.isdisjoint(test_ids)
    assert train_ids.isdisjoint(remainder_ids)
    assert test_ids.isdisjoint(remainder_ids)

    union_ids = pocket_ids | train_ids | test_ids | remainder_ids
    assert len(union_ids) == len(assigned_a)
    assert len(remainder_ids) == 0

    expected = int(round(0.20 * len(assigned_a)))
    max_group = 6
    assert abs(len(test_ids) - expected) <= max_group
