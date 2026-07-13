from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import pandas as pd

from src.ml.pocket_features import POCKET_FEATURE_COLUMNS, compute_pocket_features


POCKET_MAP_FEATURE_COLUMNS = (
    "pocket_atom_count",
    "pocket_residue_count",
    "pocket_radius_a",
)
TARGET_POCKET_FEATURE_COLUMNS = (*POCKET_FEATURE_COLUMNS, *POCKET_MAP_FEATURE_COLUMNS)


def _numeric(value: Any) -> float | None:
    parsed = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    return float(parsed) if pd.notna(parsed) else None


def build_target_pocket_feature_table(
    pocket_map: pd.DataFrame,
    *,
    cache_root: Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build one non-label-derived pocket descriptor row per PDB structure."""

    required = {"pdb_id", "pocket_pdb"}
    missing = required - set(pocket_map.columns)
    if missing:
        raise ValueError(f"pocket map missing required columns: {sorted(missing)}")

    rows: list[dict[str, Any]] = []
    ranked_map = pocket_map.copy()
    ranked_map["_pocket_exists"] = ranked_map["pocket_pdb"].map(
        lambda value: Path(str(value or "")).is_file()
    )
    ranked_map = ranked_map.sort_values(
        ["pdb_id", "_pocket_exists"],
        ascending=[True, False],
    )
    for _, source in ranked_map.drop_duplicates("pdb_id").iterrows():
        pdb_id = str(source.get("pdb_id") or "").strip().upper()
        pocket_path = Path(str(source.get("pocket_pdb") or ""))
        source_status = str(source.get("status") or "")
        ready = bool(pdb_id and pocket_path.is_file())
        row: dict[str, Any] = {
            "pdb_id": pdb_id,
            "pocket_pdb": str(pocket_path) if str(pocket_path) != "." else "",
            "pocket_feature_status": "ready" if ready else "missing_pocket_pdb",
            "pocket_definition_status": source_status,
            "pocket_definition_source": str(source.get("center_source") or ""),
            "pocket_control_ligand_id": str(source.get("center_ligand_id") or ""),
            "pocket_radius_a": _numeric(source.get("radius_a")),
            "pocket_atom_count": _numeric(source.get("pocket_atom_count")),
            "pocket_residue_count": _numeric(source.get("pocket_residue_count")),
        }
        if ready:
            row.update(compute_pocket_features(pocket_path, cache_root=cache_root))
        else:
            row.update({column: pd.NA for column in POCKET_FEATURE_COLUMNS})
        rows.append(row)

    table = pd.DataFrame(rows)
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "rows": int(len(table)),
        "pdb_ids": int(table["pdb_id"].nunique()) if not table.empty else 0,
        "ready_rows": int(table["pocket_feature_status"].eq("ready").sum())
        if not table.empty
        else 0,
        "missing_rows": int(table["pocket_feature_status"].ne("ready").sum())
        if not table.empty
        else 0,
        "feature_columns": list(TARGET_POCKET_FEATURE_COLUMNS),
        "provenance_columns": [
            "pocket_pdb",
            "pocket_feature_status",
            "pocket_definition_status",
            "pocket_definition_source",
            "pocket_control_ligand_id",
        ],
        "policy": (
            "Descriptors are target/PDB features derived without activity labels. "
            "Pocket definition provenance remains audit-only. Missing pockets remain missing."
        ),
    }
    return table, manifest


def write_target_pocket_features(
    pocket_map: pd.DataFrame,
    *,
    out_dir: Path,
    cache_root: Path | None = None,
) -> dict[str, Any]:
    destination = Path(out_dir)
    destination.mkdir(parents=True, exist_ok=True)
    table, manifest = build_target_pocket_feature_table(
        pocket_map,
        cache_root=Path(cache_root or destination / "cache"),
    )
    table_path = destination / "target_pocket_features.csv"
    table.to_csv(table_path, index=False)
    manifest["table"] = str(table_path)
    manifest["cache_root"] = str(Path(cache_root or destination / "cache"))
    (destination / "target_pocket_features_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    return manifest
