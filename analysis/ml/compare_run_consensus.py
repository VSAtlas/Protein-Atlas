from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

import pandas as pd


_LIGAND_SUFFIX = re.compile(r"\.(pdbqt|mol2)$", flags=re.IGNORECASE)
_COLUMN_TOKEN = re.compile(r"[^A-Za-z0-9]+")
_CANONICAL_DECOY_PREFIX = "dud_"
_LEGACY_DECOY_PREFIXES = (
    _CANONICAL_DECOY_PREFIX,
    "decoy_",
    "decoys_",
)
_COMPARE_STREAMS: dict[str, dict[str, Any]] = {
    "docked_consensus": {
        "decoy_filename": "dud_consensus_docking_scores.csv",
        "score_filename": "consensus_docking_scores.csv",
        "recursive": False,
        "allow_legacy_prefix": False,
    },
    "post_docked_scorch": {
        "decoy_filename": "dud_consensus_reranked_scorch.csv",
        "score_filename": "consensus_reranked_scorch.csv",
        "recursive": True,
        "allow_legacy_prefix": True,
    },
}


def _ligand_key(value: Any) -> str:
    return _LIGAND_SUFFIX.sub("", str(value or "").strip())


def _finite_float(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _is_decoy(value: Any) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "t", "yes", "y"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_sha256(payload: Any) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _score_output_column(score_column: str) -> str:
    token = _COLUMN_TOKEN.sub("_", str(score_column).strip()).strip("_").lower()
    if not token:
        raise ValueError("comparison score column must contain a letter or digit")
    if token == "consensus_score":
        return "z_vs_compare_run_consensus"
    return f"z_vs_compare_run_{token}"


def _reference_ligand_key(row: dict[str, Any]) -> str:
    return _ligand_key(row.get("ligand") or row.get("ligand_base"))


def _decoy_identification_mode(
    row: dict[str, Any],
    *,
    allow_legacy_prefix: bool,
) -> str | None:
    explicit = str(row.get("is_decoy") or "").strip()
    if explicit:
        return "explicit_is_decoy" if _is_decoy(explicit) else None
    if not allow_legacy_prefix:
        return None
    for key in ("ligand", "ligand_base"):
        identifier = str(row.get(key) or "").strip().lower()
        if identifier.startswith(_LEGACY_DECOY_PREFIXES):
            return "canonical_dud_prefix_fallback"
    return None


def _identification_summary(counts: dict[str, int]) -> str:
    populated = sorted(mode for mode, count in counts.items() if count)
    if not populated:
        return "none"
    if len(populated) == 1:
        return populated[0]
    return "mixed:" + "+".join(populated)


def _reference_paths(root: Path, stream: dict[str, Any]) -> list[Path]:
    filename = str(stream["decoy_filename"])
    if stream["recursive"]:
        return sorted(path for path in root.rglob(filename) if path.is_file())
    return sorted(path for path in root.glob(f"*/{filename}") if path.is_file())


def _pdb_from_reference_path(root: Path, path: Path) -> str:
    relative = path.relative_to(root)
    if len(relative.parts) < 2:
        raise ValueError(f"comparison reference is not under a PDB directory: {path}")
    return relative.parts[0].strip().upper()


def _read_comparison_run(
    run_root: str | Path,
    *,
    run_id: str,
    compare_stream: str = "docked_consensus",
    compare_score_column: str = "consensus_score",
) -> tuple[dict[str, dict[str, Any]], dict[tuple[str, str], float]]:
    root = Path(run_root)
    if not root.is_dir():
        raise FileNotFoundError(f"comparison run directory does not exist: {root}")
    try:
        stream = _COMPARE_STREAMS[compare_stream]
    except KeyError as exc:
        choices = ", ".join(sorted(_COMPARE_STREAMS))
        raise ValueError(
            f"unsupported comparison stream {compare_stream!r}; choose from {choices}"
        ) from exc
    decoy_paths = _reference_paths(root, stream)
    if not decoy_paths:
        raise ValueError(
            f"no {stream['decoy_filename']} reference files found under {root}"
        )
    paths_by_pdb: dict[str, list[Path]] = {}
    for path in decoy_paths:
        paths_by_pdb.setdefault(_pdb_from_reference_path(root, path), []).append(path)

    nulls: dict[str, dict[str, Any]] = {}
    fda_scores: dict[tuple[str, str], float] = {}
    expected_run_id = str(run_id).strip()
    for pdb_id, pdb_decoy_paths in sorted(paths_by_pdb.items()):
        decoy_scores: list[float] = []
        strata: set[tuple[str, str]] = set()
        identification_counts = {
            "explicit_is_decoy": 0,
            "canonical_dud_prefix_fallback": 0,
        }
        decoy_hashes: list[str] = []
        score_hashes: list[str] = []
        score_paths: list[Path] = []
        for dud_path in pdb_decoy_paths:
            decoy_hashes.append(_sha256(dud_path))
            with dud_path.open("r", newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                if compare_score_column not in (reader.fieldnames or []):
                    raise ValueError(
                        "comparison decoy file missing score column "
                        f"{compare_score_column!r}: {dud_path}"
                    )
                for row in reader:
                    mode = _decoy_identification_mode(
                        row,
                        allow_legacy_prefix=bool(stream["allow_legacy_prefix"]),
                    )
                    if mode is None:
                        continue
                    row_run_id = str(row.get("run_id") or "").strip()
                    if row_run_id != expected_run_id:
                        raise ValueError(
                            f"comparison decoy run_id mismatch in {dud_path}: "
                            f"expected {expected_run_id!r}, found {row_run_id!r}"
                        )
                    row_pdb = str(row.get("pdb_id") or "").strip().upper()
                    if row_pdb != pdb_id:
                        raise ValueError(
                            f"comparison decoy PDB mismatch in {dud_path}: "
                            f"directory {pdb_id!r}, row {row_pdb!r}"
                        )
                    strata.add(
                        (
                            str(row.get("variant") or "").strip(),
                            str(row.get("ph_label") or "").strip(),
                        )
                    )
                    score = _finite_float(row.get(compare_score_column))
                    if score is not None:
                        decoy_scores.append(score)
                        identification_counts[mode] += 1
        if len(strata) != 1:
            raise ValueError(
                f"comparison run has {len(strata)} variant/pH strata for {pdb_id}; "
                "the ML table lacks an unambiguous stratum key"
            )
        if not decoy_scores:
            policy = (
                "explicit is_decoy values or blank flags with known legacy decoy prefixes"
                if stream["allow_legacy_prefix"]
                else "explicitly marked decoys"
            )
            raise ValueError(
                f"no finite {compare_score_column} scores from {policy}: "
                + ", ".join(str(path) for path in pdb_decoy_paths)
            )
        mu = sum(decoy_scores) / len(decoy_scores)
        variance = sum((score - mu) ** 2 for score in decoy_scores) / len(decoy_scores)
        sigma = math.sqrt(variance)
        if sigma <= 0.0 or not math.isfinite(sigma):
            raise ValueError(
                f"zero or non-finite decoy {compare_score_column} sigma: "
                + ", ".join(str(path) for path in pdb_decoy_paths)
            )
        variant, ph_label = next(iter(strata))
        decoy_path_values = [str(path) for path in pdb_decoy_paths]
        null_record: dict[str, Any] = {
            "pdb_id": pdb_id,
            "comparison_run_id": expected_run_id,
            "comparison_stream": compare_stream,
            "score_column": compare_score_column,
            "variant": variant,
            "ph_label": ph_label,
            "n_decoys": len(decoy_scores),
            "n_unique": len(set(decoy_scores)),
            "mu": mu,
            "sigma": sigma,
            "decoy_identification_mode": _identification_summary(identification_counts),
            "n_decoys_explicit_is_decoy": identification_counts["explicit_is_decoy"],
            "n_decoys_canonical_dud_prefix_fallback": identification_counts[
                "canonical_dud_prefix_fallback"
            ],
            "canonical_decoy_prefix": _CANONICAL_DECOY_PREFIX,
            "decoy_reference_file": json.dumps(decoy_path_values),
            "decoy_reference_sha256": json.dumps(decoy_hashes),
            # Preserve the original null-source aliases for downstream readers.
            "dud_consensus_file": decoy_path_values[0],
            "dud_consensus_sha256": decoy_hashes[0],
        }

        for dud_path in pdb_decoy_paths:
            fda_path = dud_path.with_name(str(stream["score_filename"]))
            if not fda_path.is_file():
                continue
            score_paths.append(fda_path)
            score_hashes.append(_sha256(fda_path))
            with fda_path.open("r", newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                if compare_score_column not in (reader.fieldnames or []):
                    raise ValueError(
                        "comparison score file missing column "
                        f"{compare_score_column!r}: {fda_path}"
                    )
                for row in reader:
                    row_run_id = str(row.get("run_id") or "").strip()
                    if row_run_id != expected_run_id:
                        raise ValueError(
                            f"comparison score run_id mismatch in {fda_path}: "
                            f"expected {expected_run_id!r}, found {row_run_id!r}"
                        )
                    row_pdb = str(row.get("pdb_id") or "").strip().upper()
                    if row_pdb != pdb_id:
                        raise ValueError(
                            f"comparison score PDB mismatch in {fda_path}: "
                            f"directory {pdb_id!r}, row {row_pdb!r}"
                        )
                    ligand = _reference_ligand_key(row)
                    score = _finite_float(row.get(compare_score_column))
                    if ligand and score is not None:
                        key = (pdb_id, ligand)
                        prior = fda_scores.get(key)
                        if prior is not None and abs(prior - score) > 1e-12:
                            raise ValueError(
                                f"conflicting comparison scores for {pdb_id}/{ligand}"
                            )
                        fda_scores[key] = score
        if score_paths:
            score_path_values = [str(path) for path in score_paths]
            null_record["score_reference_file"] = json.dumps(score_path_values)
            null_record["score_reference_sha256"] = json.dumps(score_hashes)
            null_record["fda_consensus_file"] = score_path_values[0]
            null_record["fda_consensus_sha256"] = score_hashes[0]
        null_record["null_sha256"] = _json_sha256(
            {
                "comparison_run_id": expected_run_id,
                "comparison_stream": compare_stream,
                "score_column": compare_score_column,
                "pdb_id": pdb_id,
                "variant": variant,
                "ph_label": ph_label,
                "n_decoys": len(decoy_scores),
                "mu": mu,
                "sigma": sigma,
                "identification_counts": identification_counts,
                "decoy_reference_sha256": decoy_hashes,
            }
        )
        nulls[pdb_id] = null_record
    return nulls, fda_scores


def compare_consensus_to_run_decoys(
    dataset_path: str | Path,
    comparison_run_root: str | Path,
    out_path: str | Path,
    *,
    comparison_run_id: str,
    compare_stream: str = "docked_consensus",
    compare_score_column: str = "consensus_score",
) -> dict[str, Any]:
    source = Path(dataset_path)
    out = Path(out_path)
    frame = pd.read_csv(source, low_memory=False)
    required = {"pdb_id", "ligand_base", compare_score_column}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"dataset missing required columns: {', '.join(missing)}")

    nulls, fda_scores = _read_comparison_run(
        comparison_run_root,
        run_id=comparison_run_id,
        compare_stream=compare_stream,
        compare_score_column=compare_score_column,
    )
    raw_col = compare_score_column
    if compare_score_column == "consensus_score" and "consensus_score_raw" in frame:
        raw_col = "consensus_score_raw"
    raw = pd.to_numeric(frame[raw_col], errors="coerce")
    pdb_keys = frame["pdb_id"].fillna("").astype(str).str.strip().str.upper()
    ligand_keys = frame["ligand_base"].map(_ligand_key)
    mu = pdb_keys.map({key: value["mu"] for key, value in nulls.items()})
    sigma = pdb_keys.map({key: value["sigma"] for key, value in nulls.items()})
    z_score = (raw - mu) / sigma
    missing_null = raw.notna() & z_score.isna()
    if missing_null.any():
        missing_pdbs = sorted(set(pdb_keys.loc[missing_null]))
        raise ValueError(
            f"{int(missing_null.sum())} score rows lack comparison-run decoy nulls: "
            + ", ".join(missing_pdbs)
        )

    expected_fda = pd.Series(
        [
            fda_scores.get((pdb_id, ligand))
            for pdb_id, ligand in zip(pdb_keys, ligand_keys, strict=False)
        ],
        index=frame.index,
        dtype="float64",
    )
    exact_match = (
        expected_fda.notna() & raw.notna() & expected_fda.sub(raw).abs().le(1e-12)
    )
    score_mismatch = expected_fda.notna() & raw.notna() & ~exact_match
    null_only_source = (
        "comparison_run_dud_null_only"
        if compare_stream == "docked_consensus"
        else "comparison_run_post_docked_scorch_decoy_null_only"
    )
    source_name = pd.Series(
        "missing_comparison_null", index=frame.index, dtype="object"
    )
    source_name.loc[z_score.notna()] = null_only_source
    source_name.loc[exact_match] = "same_ligand_and_score_in_comparison_run"
    source_name.loc[score_mismatch] = "same_ligand_score_mismatch_to_comparison_run"

    prefix = _score_output_column(compare_score_column)
    source_sha256 = _sha256(source)
    run_sha256 = _json_sha256(
        {
            "comparison_run_id": str(comparison_run_id),
            "comparison_stream": compare_stream,
            "score_column": compare_score_column,
            "nulls": [
                {
                    "pdb_id": value["pdb_id"],
                    "null_sha256": value["null_sha256"],
                    "score_reference_sha256": value.get("score_reference_sha256", ""),
                }
                for value in sorted(nulls.values(), key=lambda item: item["pdb_id"])
            ],
        }
    )
    frame[prefix] = z_score
    frame[f"{prefix}_source"] = source_name
    frame[f"{prefix}_run_id"] = str(comparison_run_id)
    frame[f"{prefix}_source_sha256"] = source_sha256
    frame[f"{prefix}_run_sha256"] = run_sha256
    frame[f"{prefix}_null_sha256"] = pdb_keys.map(
        {key: value["null_sha256"] for key, value in nulls.items()}
    )
    frame[f"{prefix}_decoy_n"] = pdb_keys.map(
        {key: value["n_decoys"] for key, value in nulls.items()}
    )
    frame[f"{prefix}_decoy_mu"] = mu
    frame[f"{prefix}_decoy_sigma"] = sigma
    frame[f"{prefix}_decoy_identification_mode"] = pdb_keys.map(
        {key: value["decoy_identification_mode"] for key, value in nulls.items()}
    )
    frame[f"{prefix}_decoy_explicit_n"] = pdb_keys.map(
        {key: value["n_decoys_explicit_is_decoy"] for key, value in nulls.items()}
    )
    frame[f"{prefix}_decoy_prefix_fallback_n"] = pdb_keys.map(
        {
            key: value["n_decoys_canonical_dud_prefix_fallback"]
            for key, value in nulls.items()
        }
    )
    decoy_reference_hashes = pdb_keys.map(
        {key: value["dud_consensus_sha256"] for key, value in nulls.items()}
    )
    frame[f"{prefix}_decoy_sha256"] = decoy_reference_hashes
    frame[f"{prefix}_dud_sha256"] = decoy_reference_hashes

    out.parent.mkdir(parents=True, exist_ok=True)
    null_path = out.with_name(f"{out.stem}.compare_run_nulls.csv")
    pd.DataFrame(nulls.values()).sort_values("pdb_id").to_csv(null_path, index=False)
    nulls_sha256 = _sha256(null_path)
    frame[f"{prefix}_nulls_sha256"] = nulls_sha256
    frame.to_csv(out, index=False)

    explicit_decoy_count = sum(
        int(value["n_decoys_explicit_is_decoy"]) for value in nulls.values()
    )
    fallback_decoy_count = sum(
        int(value["n_decoys_canonical_dud_prefix_fallback"]) for value in nulls.values()
    )
    identification_modes: dict[str, int] = {}
    for value in nulls.values():
        mode = str(value["decoy_identification_mode"])
        identification_modes[mode] = identification_modes.get(mode, 0) + int(
            value["n_decoys"]
        )
    if (
        compare_stream == "docked_consensus"
        and compare_score_column == "consensus_score"
    ):
        formula = (
            "z_vs_compare_run_consensus = (raw consensus_score - comparison-run "
            "explicit-DUD mean) / comparison-run explicit-DUD population std"
        )
    else:
        formula = (
            f"{prefix} = ({raw_col} - same-PDB comparison-run decoy mean) / "
            "same-PDB comparison-run decoy population std"
        )
    manifest: dict[str, Any] = {
        "dataset": str(source),
        "output": str(out),
        "comparison_run_id": str(comparison_run_id),
        "comparison_run_root": str(comparison_run_root),
        "comparison_stream": compare_stream,
        "compare_stream": compare_stream,
        "comparison_score_column": compare_score_column,
        "compare_score_column": compare_score_column,
        "score_input_column": raw_col,
        "score_output_column": prefix,
        "n_rows": int(len(frame)),
        "n_scored": int(z_score.notna().sum()),
        "n_same_ligand_and_score": int(exact_match.sum()),
        "n_same_ligand_score_mismatch": int(score_mismatch.sum()),
        "n_comparison_null_only": int(source_name.eq(null_only_source).sum()),
        "n_reference_pdbs": int(len(nulls)),
        "n_decoys_explicit_is_decoy": explicit_decoy_count,
        "n_decoys_canonical_dud_prefix_fallback": fallback_decoy_count,
        "decoy_identification_modes": identification_modes,
        "decoy_policy": (
            "Populated is_decoy is authoritative; blank legacy post-docked flags "
            "are accepted only for ligand/ligand_base identifiers beginning dud_, decoy_, or decoys_."
            if compare_stream == "post_docked_scorch"
            else "Only rows with explicitly truthy is_decoy define the null."
        ),
        "formula": formula,
        "hashes": {
            "dataset_source_sha256": source_sha256,
            "comparison_run_sha256": run_sha256,
            "nulls_csv_sha256": nulls_sha256,
            "output_dataset_sha256": _sha256(out),
        },
        "canonical_score_columns_overwritten": False,
        "final_score_overwritten": False,
        "outputs": {"dataset": str(out), "nulls": str(null_path)},
    }
    manifest_path = out.with_name(f"{out.stem}.compare_run_manifest.json")
    manifest["outputs"]["manifest"] = str(manifest_path)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest
