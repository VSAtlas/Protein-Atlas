from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path
from typing import Any

import pandas as pd

from analysis.ml.feature_backfill import _finalize_score_priors

_LIGAND_SUFFIX = re.compile(r"\.(pdbqt|mol2)$", flags=re.IGNORECASE)


def _ligand_key(value: Any) -> str:
    return _LIGAND_SUFFIX.sub("", str(value or "").strip())


def _is_true(value: Any) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "t", "yes", "y"}


def _finite_float(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _reference_score_layouts(root: Path) -> list[tuple[str, Path, str]]:
    candidates = sorted(root.rglob("dud_consensus_docking_scores.csv"))
    if not candidates:
        raise ValueError(f"no DUD consensus files found under {root}")

    by_pdb: dict[str, list[tuple[Path, str]]] = {}
    for dud_path in candidates:
        relative = dud_path.relative_to(root)
        if len(relative.parts) < 2:
            raise ValueError(
                "DUD consensus file is not nested under a PDB directory: "
                f"{dud_path}"
            )
        pdb_id = relative.parts[0].strip().upper()
        if not pdb_id:
            raise ValueError(f"empty PDB directory in reference path: {dud_path}")
        layout_parts = relative.parts[1:-1]
        layout = "/".join(layout_parts) if layout_parts else "flat"
        by_pdb.setdefault(pdb_id, []).append((dud_path.parent, layout))

    ambiguous = {
        pdb_id: entries for pdb_id, entries in by_pdb.items() if len(entries) > 1
    }
    if ambiguous:
        details = "; ".join(
            f"{pdb_id}={','.join(layout for _path, layout in entries)}"
            for pdb_id, entries in sorted(ambiguous.items())
        )
        raise ValueError(
            "ambiguous DUD consensus layouts; select a single variant/pH "
            f"reference stratum per PDB: {details}"
        )
    return [
        (pdb_id, entries[0][0], entries[0][1])
        for pdb_id, entries in sorted(by_pdb.items())
    ]


def _read_reference_run(
    run_root: str | Path,
) -> tuple[dict[str, dict[str, Any]], dict[tuple[str, str], float]]:
    root = Path(run_root)
    if not root.is_dir():
        raise FileNotFoundError(f"reference run directory does not exist: {root}")

    nulls: dict[str, dict[str, Any]] = {}
    fda_scores: dict[tuple[str, str], float] = {}
    conflicts: list[tuple[str, str, float, float]] = []

    for pdb_id, pdb_dir, layout in _reference_score_layouts(root):
        dud_path = pdb_dir / "dud_consensus_docking_scores.csv"
        decoy_scores: list[float] = []
        with dud_path.open("r", newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if not _is_true(row.get("is_decoy")):
                    continue
                score = _finite_float(row.get("consensus_score"))
                if score is not None:
                    decoy_scores.append(score)
        if not decoy_scores:
            raise ValueError(f"no explicitly marked decoy consensus scores: {dud_path}")
        mu = sum(decoy_scores) / len(decoy_scores)
        variance = sum((score - mu) ** 2 for score in decoy_scores) / len(decoy_scores)
        sigma = math.sqrt(variance)
        if sigma <= 0 or not math.isfinite(sigma):
            raise ValueError(f"zero or non-finite decoy consensus sigma: {dud_path}")
        nulls[pdb_id] = {
            "pdb_id": pdb_id,
            "n_decoys": len(decoy_scores),
            "n_unique": len(set(decoy_scores)),
            "zero_fraction": sum(score == 0.0 for score in decoy_scores) / len(decoy_scores),
            "mu": mu,
            "sigma": sigma,
            "dud_consensus_file": str(dud_path),
            "reference_layout": layout,
            "reference_score_dir": str(pdb_dir),
        }

        fda_path = pdb_dir / "consensus_docking_scores.csv"
        if not fda_path.is_file():
            raise FileNotFoundError(f"FDA consensus file missing: {fda_path}")
        with fda_path.open("r", newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                score = _finite_float(row.get("consensus_score"))
                ligand = _ligand_key(row.get("ligand"))
                if score is None or not ligand:
                    continue
                key = (pdb_id, ligand)
                prior = fda_scores.get(key)
                if prior is not None and abs(prior - score) > 1e-12:
                    conflicts.append((pdb_id, ligand, prior, score))
                fda_scores[key] = score
    if conflicts:
        first = conflicts[0]
        raise ValueError(
            "conflicting FDA consensus scores for "
            f"{first[0]}/{first[1]}: {first[2]} versus {first[3]}"
        )
    return nulls, fda_scores


def repair_consensus_z_scores(
    dataset_path: str | Path,
    reference_run_root: str | Path,
    out_path: str | Path,
    *,
    allow_reference_fallback: bool = False,
    replace_model_score_columns: bool = False,
) -> dict[str, Any]:
    source = Path(dataset_path)
    out = Path(out_path)
    frame = pd.read_csv(source, low_memory=False)
    required = {"pdb_id", "ligand_base", "consensus_score"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"dataset missing required columns: {', '.join(missing)}")

    nulls, fda_scores = _read_reference_run(reference_run_root)
    raw_consensus_col = (
        "consensus_score_raw"
        if "consensus_score_raw" in frame.columns
        else "consensus_score"
    )
    raw_consensus = pd.to_numeric(frame[raw_consensus_col], errors="coerce")
    pdb_keys = frame["pdb_id"].fillna("").astype(str).str.strip().str.upper()
    ligand_keys = frame["ligand_base"].map(_ligand_key)
    expected_fda = pd.Series(
        [fda_scores.get((pdb_id, ligand)) for pdb_id, ligand in zip(pdb_keys, ligand_keys)],
        index=frame.index,
        dtype="float64",
    )
    same_run_match = (
        expected_fda.notna()
        & raw_consensus.notna()
        & expected_fda.sub(raw_consensus).abs().le(1e-12)
    )

    mu = pdb_keys.map({key: value["mu"] for key, value in nulls.items()})
    sigma = pdb_keys.map({key: value["sigma"] for key, value in nulls.items()})
    z_score = (raw_consensus - mu) / sigma
    reference_fallback = z_score.notna() & ~same_run_match
    if reference_fallback.any() and not allow_reference_fallback:
        raise ValueError(
            f"{int(reference_fallback.sum())} rows require a cross-run reference null; "
            "rerun with allow_reference_fallback=True to mark and retain them"
        )

    source_name = pd.Series("missing_consensus_decoy_null", index=frame.index, dtype="object")
    source_name.loc[same_run_match] = "exact_reference_run_fda_and_dud_consensus_z"
    source_name.loc[reference_fallback] = "borrowed_reference_run_dud_consensus_z"
    missing_z = raw_consensus.notna() & z_score.isna()
    if missing_z.any():
        missing_pdb = sorted(set(pdb_keys.loc[missing_z]))
        raise ValueError(
            f"{int(missing_z.sum())} raw consensus rows lack a reference null for PDBs: "
            + ", ".join(missing_pdb)
        )

    preserved = {
        "consensus_score": "consensus_score_raw",
        "atlas_score": "atlas_score_pre_z_repair",
        "z_selected": "z_selected_pre_z_repair",
        "z_selected_source": "z_selected_source_pre_z_repair",
        "atlas_binding_prior": "atlas_binding_prior_pre_z_repair",
        "atlas_binding_prior_source": "atlas_binding_prior_source_pre_z_repair",
    }
    if replace_model_score_columns:
        for original, backup in preserved.items():
            if original in frame.columns and backup not in frame.columns:
                frame[backup] = frame[original]

    frame["consensus_z_score"] = z_score
    frame["consensus_z_score_source"] = source_name
    frame["consensus_z_decoy_n"] = pdb_keys.map(
        {key: value["n_decoys"] for key, value in nulls.items()}
    )
    frame["consensus_z_decoy_unique"] = pdb_keys.map(
        {key: value["n_unique"] for key, value in nulls.items()}
    )
    frame["consensus_z_decoy_zero_fraction"] = pdb_keys.map(
        {key: value["zero_fraction"] for key, value in nulls.items()}
    )
    frame["consensus_z_decoy_mu"] = mu
    frame["consensus_z_decoy_sigma"] = sigma

    if replace_model_score_columns:
        frame["consensus_score"] = z_score
        frame["atlas_score"] = z_score
        frame["z_selected"] = z_score
        frame["z_selected_source"] = source_name
        frame["consensus_score_feature_source"] = source_name
        frame["z_selected_feature_source"] = source_name
        frame["atlas_score_source_for_ml"] = source_name

        stale = [
            "atlas_score_normalized",
            "consensus_score_normalized",
            "consensus_z_score_normalized",
            "atlas_binding_prior",
            "atlas_binding_prior_source",
            "banana_atlas_blend_score",
            "binding_expert_score",
            "binding_expert_source",
        ]
        frame = frame.drop(columns=[col for col in stale if col in frame.columns])
        frame, prior_summary = _finalize_score_priors(frame)
        frame["atlas_binding_prior_source"] = source_name
    else:
        prior_summary = {
            "status": "not_recomputed",
            "reason": "canonical model score columns were preserved",
        }

    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out, index=False)

    null_quality = pd.DataFrame(nulls.values()).sort_values("pdb_id")
    null_quality_path = out.with_name(f"{out.stem}.consensus_decoy_null_quality.csv")
    null_quality.to_csv(null_quality_path, index=False)
    source_counts = (
        source_name.value_counts(dropna=False)
        .rename_axis("consensus_z_score_source")
        .reset_index(name="n_rows")
    )
    source_counts["fraction"] = source_counts["n_rows"] / len(frame)
    source_counts_path = out.with_name(f"{out.stem}.consensus_z_source_counts.csv")
    source_counts.to_csv(source_counts_path, index=False)

    manifest = {
        "dataset": str(source),
        "reference_run_root": str(reference_run_root),
        "output": str(out),
        "n_rows": int(len(frame)),
        "n_raw_consensus": int(raw_consensus.notna().sum()),
        "n_exact_reference_run_fda_and_dud": int(same_run_match.sum()),
        "n_borrowed_reference_null": int(reference_fallback.sum()),
        "source_by_label_source": (
            pd.DataFrame(
                {
                    "label_source": frame.get(
                        "label_source",
                        pd.Series("missing", index=frame.index),
                    ).fillna("missing"),
                    "consensus_z_score_source": source_name,
                }
            )
            .value_counts()
            .rename("n_rows")
            .reset_index()
            .to_dict("records")
        ),
        "n_missing_z": int(frame["consensus_z_score"].isna().sum()),
        "n_reference_pdbs": int(len(nulls)),
        "allow_reference_fallback": bool(allow_reference_fallback),
        "replace_model_score_columns": bool(replace_model_score_columns),
        "formula": "z = (FDA consensus_score - DUD consensus mean) / DUD consensus population std",
        "decoy_policy": "Only rows with is_decoy == 1 define the null; controls are excluded.",
        "legacy_alias_policy": (
            "canonical aliases replaced by consensus_z_score under explicit opt-in"
            if replace_model_score_columns
            else "canonical score columns preserved; consensus_z_score is additive"
        ),
        "prior_recompute": prior_summary,
        "outputs": {
            "dataset": str(out),
            "null_quality": str(null_quality_path),
            "source_counts": str(source_counts_path),
        },
    }
    manifest_path = out.with_name(f"{out.stem}.consensus_z_repair_manifest.json")
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest
