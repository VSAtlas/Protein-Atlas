from __future__ import annotations

import hashlib
import logging
import math
import random
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, TypeVar

try:  # optional dependency
    import pandas as pd  # type: ignore[import-untyped]
except Exception:  # pragma: no cover - optional dependency
    pd = None

try:  # optional dependency
    from rdkit import Chem  # type: ignore[import-untyped]
    from rdkit.Chem.Scaffolds import MurckoScaffold  # type: ignore[import-untyped]
except Exception:  # pragma: no cover - optional dependency
    Chem = None  # type: ignore[assignment]
    MurckoScaffold = None  # type: ignore[assignment]

_ALLOWED_POLICIES = {"stratified_scaffold", "stratified_binding_class", "random"}


def _stable_hash(text: str) -> int:
    return int(hashlib.sha1(text.encode("utf-8")).hexdigest()[:16], 16)


def _seed_for_key(seed: int, key: str) -> int:
    return int((_stable_hash(f"{seed}:{key}") % (2**31 - 1)) + 1)


def _as_records(bags_df: Any) -> Tuple[List[Dict[str, Any]], bool]:
    if pd is not None and isinstance(bags_df, pd.DataFrame):
        return bags_df.to_dict(orient="records"), True
    if isinstance(bags_df, list):
        return [dict(row) for row in bags_df], False
    if isinstance(bags_df, Iterable):
        return [dict(row) for row in bags_df], False
    raise TypeError("bags_df must be a pandas DataFrame or iterable of dicts")


def _return_records(records: List[Dict[str, Any]], as_df: bool) -> Any:
    if as_df and pd is not None:
        return pd.DataFrame(records)
    return records


def _row_uid(row: Dict[str, Any], idx: int) -> str:
    for key in ("bag_uid", "ligand_uid", "ligand_id"):
        val = row.get(key)
        if val:
            return str(val)
    return f"row_{idx}"


def _scaffold_key(smiles: Optional[str], inchikey: Optional[str]) -> str:
    smi = str(smiles).strip() if smiles else ""
    ik = str(inchikey).strip() if inchikey else ""

    if smi and Chem is not None and MurckoScaffold is not None:
        try:
            mol = Chem.MolFromSmiles(smi)
        except Exception:
            mol = None
        if mol is not None:
            try:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            except Exception:
                scaffold = None
            if scaffold is not None:
                try:
                    scaffold_smi = Chem.MolToSmiles(scaffold, isomericSmiles=False)
                except Exception:
                    scaffold_smi = ""
                if scaffold_smi:
                    return scaffold_smi

    if ik:
        return ik[:14] if len(ik) >= 14 else ik
    if smi:
        return f"smiles:{_stable_hash(smi)}"
    return "unknown"


def _normalize_fractions(
    test_fraction: float, remainder_eval_fraction: float
) -> Tuple[float, float]:
    test_fraction = float(test_fraction or 0.0)
    remainder_eval_fraction = float(remainder_eval_fraction or 0.0)
    if test_fraction < 0:
        test_fraction = 0.0
    if remainder_eval_fraction < 0:
        remainder_eval_fraction = 0.0
    if test_fraction + remainder_eval_fraction > 1.0:
        test_fraction = max(0.0, 1.0 - remainder_eval_fraction)
    return test_fraction, remainder_eval_fraction


def assign_holdout_sets(
    bags_df: Any,
    seed: int,
    test_fraction: float,
    remainder_eval_fraction: float,
) -> Any:
    records, as_df = _as_records(bags_df)
    seed = int(seed or 0)
    test_fraction, remainder_eval_fraction = _normalize_fractions(
        test_fraction, remainder_eval_fraction
    )

    groups: Dict[str, List[int]] = {}
    for idx, row in enumerate(records):
        smi = row.get("canonical_smiles") or row.get("smiles")
        ik = row.get("inchikey") or row.get("inchi_key") or row.get("inchiKey")
        if smi and not row.get("canonical_smiles"):
            row["canonical_smiles"] = smi
        row["scaffold_key"] = _scaffold_key(smi, ik)
        groups.setdefault(row["scaffold_key"], []).append(idx)

    assignments: Dict[str, str] = {}
    for scaffold in sorted(groups.keys()):
        hval = _stable_hash(f"{seed}:{scaffold}")
        frac = (hval % 10_000_000) / 10_000_000.0
        if remainder_eval_fraction > 0 and frac < remainder_eval_fraction:
            assignments[scaffold] = "remainder_eval"
        elif frac < remainder_eval_fraction + test_fraction:
            assignments[scaffold] = "test"
        else:
            assignments[scaffold] = "train_pool"

    for row in records:
        row["set_name"] = assignments.get(row["scaffold_key"], "train_pool")

    return _return_records(records, as_df)


def _allocate_targets(class_counts: Dict[str, int], max_n: int) -> Dict[str, int]:
    if not class_counts or max_n <= 0:
        return {k: 0 for k in class_counts}
    total = sum(class_counts.values())
    if total <= 0:
        return {k: 0 for k in class_counts}

    exact = {k: (class_counts[k] / total) * max_n for k in class_counts}
    targets = {k: int(math.floor(v)) for k, v in exact.items()}
    remainder = max_n - sum(targets.values())

    for k, _ in sorted(
        exact.items(), key=lambda kv: (-kv[1] + math.floor(kv[1]), kv[0])
    ):
        if remainder <= 0:
            break
        targets[k] += 1
        remainder -= 1

    classes = [k for k, v in class_counts.items() if v > 0]
    if len(classes) <= max_n:
        for k in classes:
            if targets[k] == 0:
                donor = max(
                    (c for c in targets if targets[c] > 1),
                    key=lambda c: targets[c],
                    default=None,
                )
                if donor is None:
                    continue
                targets[donor] -= 1
                targets[k] = 1

    return targets


T = TypeVar("T")


def _shuffle(rows: List[T], seed: int) -> List[T]:
    rng = random.Random(seed)
    rows = list(rows)
    rng.shuffle(rows)
    return rows


def select_pocket_eval_sample(
    train_pool_df: Any, max_n: Optional[int], policy: str, seed: int
) -> Tuple[Any, Any]:
    records, as_df = _as_records(train_pool_df)
    seed = int(seed or 0)
    if not records:
        return _return_records([], as_df), _return_records([], as_df)

    if max_n is None or int(max_n) <= 0:
        return _return_records(records, as_df), _return_records([], as_df)

    max_n = min(int(max_n), len(records))
    policy_norm = str(policy or "").strip().lower()
    if policy_norm not in _ALLOWED_POLICIES:
        policy_norm = "stratified_binding_class"

    indexed: List[Tuple[int, Dict[str, Any]]] = list(enumerate(records))
    if policy_norm == "random":
        shuffled = _shuffle(
            [row for _, row in indexed], seed=_seed_for_key(seed, "random")
        )
        selected = shuffled[:max_n]
    elif policy_norm == "stratified_binding_class":
        by_class: Dict[str, List[Dict[str, Any]]] = {}
        for idx, row in indexed:
            cls = str(row.get("y_binding_class") or "")
            by_class.setdefault(cls, []).append(row)
        targets = _allocate_targets({k: len(v) for k, v in by_class.items()}, max_n)
        selected = []
        remaining = []
        for cls in sorted(by_class.keys()):
            rows = sorted(by_class[cls], key=lambda r: _row_uid(r, 0))
            rows = _shuffle(rows, seed=_seed_for_key(seed, f"class:{cls}"))
            take = min(targets.get(cls, 0), len(rows))
            selected.extend(rows[:take])
            remaining.extend(rows[take:])
        if len(selected) < max_n and remaining:
            remaining = _shuffle(
                sorted(remaining, key=lambda r: _row_uid(r, 0)),
                seed=_seed_for_key(seed, "class:fill"),
            )
            selected.extend(remaining[: max_n - len(selected)])
    else:  # stratified_scaffold
        by_scaffold: Dict[str, List[Dict[str, Any]]] = {}
        for idx, row in indexed:
            scaffold = row.get("scaffold_key") or _scaffold_key(
                row.get("canonical_smiles") or row.get("smiles"),
                row.get("inchikey") or row.get("inchi_key") or row.get("inchiKey"),
            )
            row["scaffold_key"] = scaffold
            by_scaffold.setdefault(scaffold, []).append(row)

        class_counts: Dict[str, int] = {}
        for row in records:
            cls = str(row.get("y_binding_class") or "")
            class_counts[cls] = class_counts.get(cls, 0) + 1
        targets = _allocate_targets(class_counts, max_n)

        scaffold_order = sorted(by_scaffold.keys())
        scaffold_order = _shuffle(scaffold_order, seed=_seed_for_key(seed, "scaffold"))
        selected = []
        selected_ids: set[str] = set()
        remaining_targets = dict(targets)

        for scaffold in scaffold_order:
            if len(selected) >= max_n:
                break
            rows = by_scaffold.get(scaffold, [])
            if not rows:
                continue
        by_class_scaffold: Dict[str, List[Dict[str, Any]]] = {}
        for row in rows:
            cls = str(row.get("y_binding_class") or "")
            by_class_scaffold.setdefault(cls, []).append(row)
        class_priority = sorted(
            remaining_targets.keys(),
            key=lambda c: (-remaining_targets.get(c, 0), c),
        )
        chosen_row = None
        for cls in class_priority:
            if remaining_targets.get(cls, 0) <= 0:
                continue
            candidates = by_class_scaffold.get(cls, [])
            if not candidates:
                continue
                candidates = sorted(candidates, key=lambda r: _row_uid(r, 0))
                candidates = _shuffle(
                    candidates, seed=_seed_for_key(seed, f"{scaffold}:{cls}")
                )
                chosen_row = candidates[0]
                remaining_targets[cls] -= 1
                break
            if chosen_row is None:
                continue
            row_id = _row_uid(chosen_row, 0)
            if row_id in selected_ids:
                continue
            selected_ids.add(row_id)
            selected.append(chosen_row)

        remaining_rows = [
            row for row in records if _row_uid(row, 0) not in selected_ids
        ]

        for cls in sorted(remaining_targets.keys()):
            need = remaining_targets.get(cls, 0)
            if need <= 0:
                continue
            candidates = [
                r for r in remaining_rows if str(r.get("y_binding_class") or "") == cls
            ]
            if not candidates:
                continue
            candidates = sorted(candidates, key=lambda r: _row_uid(r, 0))
            candidates = _shuffle(candidates, seed=_seed_for_key(seed, f"fill:{cls}"))
            take = min(need, len(candidates))
            selected.extend(candidates[:take])
            selected_ids.update({_row_uid(r, 0) for r in candidates[:take]})
            remaining_rows = [
                r for r in remaining_rows if _row_uid(r, 0) not in selected_ids
            ]

        if len(selected) < max_n and remaining_rows:
            remaining_rows = sorted(remaining_rows, key=lambda r: _row_uid(r, 0))
            remaining_rows = _shuffle(
                remaining_rows, seed=_seed_for_key(seed, "scaffold:fill")
            )
            selected.extend(remaining_rows[: max_n - len(selected)])

    selected_ids = {_row_uid(row, 0) for row in selected}
    train_rows = [row for row in records if _row_uid(row, 0) not in selected_ids]
    return _return_records(selected, as_df), _return_records(train_rows, as_df)


def build_calibrator_sets_table(
    bags_df: Any,
    *,
    pdb_id: str,
    primary_uniprot: Optional[str],
    sampling_policy: str,
    sampling_seed: int,
    split_regime: str = "ligand_scaffold_holdout",
) -> Any:
    records, as_df = _as_records(bags_df)
    seen: Dict[str, str] = {}
    out: List[Dict[str, Any]] = []
    for idx, row in enumerate(records):
        bag_uid = row.get("bag_uid") or row.get("ligand_uid") or row.get("ligand_id")
        if not bag_uid:
            continue
        set_name = row.get("set_name")
        if not set_name:
            continue
        bag_uid = str(bag_uid)
        if bag_uid in seen and seen[bag_uid] != set_name:
            raise ValueError(f"bag_uid {bag_uid} assigned to multiple sets")
        if bag_uid in seen:
            continue
        seen[bag_uid] = str(set_name)
        out.append(
            {
                "bag_uid": bag_uid,
                "pdb_id": pdb_id,
                "primary_uniprot": primary_uniprot,
                "ligand_uid": str(row.get("ligand_uid") or bag_uid),
                "set_name": str(set_name),
                "sampling_policy": str(sampling_policy),
                "sampling_seed": int(sampling_seed),
                "split_regime": str(split_regime),
                "scaffold_key": row.get("scaffold_key"),
                "y_binding_class": row.get("y_binding_class"),
            }
        )
    return _return_records(out, as_df)


def write_calibrator_sets_table(
    out_dir: Path,
    rows: List[Dict[str, Any]],
    *,
    table_format: str,
    logger: Optional[logging.Logger] = None,
) -> Optional[Path]:
    import csv

    log = logger or logging.getLogger(__name__)
    if not rows:
        return None
    out_dir.mkdir(parents=True, exist_ok=True)

    fmt = str(table_format or "parquet").strip().lower()
    if fmt == "parquet":
        if pd is None:
            log.info("[calibrator_sets] parquet.skip reason=pandas_missing")
        else:
            try:
                df = pd.DataFrame(rows)
                out_path = out_dir / "calibrator_sets.parquet"
                df.to_parquet(out_path, index=False)
                return out_path
            except Exception as exc:
                log.info(
                    "[calibrator_sets] parquet.skip reason=write_error err=%s", exc
                )

    out_path = out_dir / "calibrator_sets.csv"
    fieldnames = list(rows[0].keys())
    try:
        with out_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        return out_path
    except Exception as exc:
        log.info("[calibrator_sets] csv.write_failed err=%s", exc)
        return None
