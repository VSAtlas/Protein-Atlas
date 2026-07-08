from __future__ import annotations

import csv
import datetime
import hashlib
import json
import logging
import os
import re
import tempfile
from pathlib import Path
from typing import Any, Mapping, Optional

from config.output_paths import runtime_root
from docking.completion_markers import completion_marker_path
from docking.docking_ligands import _chunk_ligand_key


def _sanitize_run_mode_token(token: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_]+", "_", str(token or ""))
    cleaned = cleaned.strip("_")
    return cleaned or "custom"


def _csv_prefix_for_run_mode(run_mode: Optional[str]) -> Optional[str]:
    token = str(run_mode or "").strip().lower()
    if not token or token == "fda":
        return ""
    if token == "dud":
        return "dud_"
    if token == "hmdb":
        return "hmdb_"
    return f"{_sanitize_run_mode_token(token)}_"


def _normalize_ph_tag_token(ph_tag: Optional[str]) -> str:
    token = str(ph_tag or "").strip()
    return token if token else "base"


def _resolve_combo_docking_csv(
    combo_docked_dir: Path,
    library_name: str,
    *,
    run_mode: Optional[str],
    kind: str,
) -> Path:
    suffix = f"_docking_score_{kind}.csv"
    default_path = combo_docked_dir / f"docking_score_{kind}.csv"
    prefix = _csv_prefix_for_run_mode(run_mode)
    if prefix is not None:
        run_mode_path = combo_docked_dir / f"{prefix}docking_score_{kind}.csv"
        if run_mode_path.exists():
            return run_mode_path

    prefixed = sorted(
        p
        for p in combo_docked_dir.glob(f"*{suffix}")
        if p.name != default_path.name
    )

    if not prefixed:
        return default_path

    lib = str(library_name or "").strip()
    if lib:
        wanted = {lib.lower(), _chunk_ligand_key(lib).lower()}
        for candidate in prefixed:
            prefix = candidate.name[: -len(suffix)].strip().lower()
            if prefix in wanted:
                return candidate

    if default_path.exists():
        return default_path
    return prefixed[0]


def _resolve_combo_docking_summary_csv(
    combo_docked_dir: Path,
    library_name: str,
    *,
    run_mode: Optional[str] = None,
) -> Path:
    return _resolve_combo_docking_csv(
        combo_docked_dir,
        library_name,
        run_mode=run_mode,
        kind="summary",
    )


def _resolve_combo_docking_long_csv(
    combo_docked_dir: Path,
    library_name: str,
    *,
    run_mode: Optional[str] = None,
) -> Path:
    return _resolve_combo_docking_csv(
        combo_docked_dir,
        library_name,
        run_mode=run_mode,
        kind="long",
    )


def _load_scored_ligand_keys_from_summary(path: Path) -> set[str]:
    scored: set[str] = set()
    if not path.exists():
        return scored

    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames or []
            ligand_col = next(
                (
                    c
                    for c in ("Ligand", "ligand", "ligand_file", "ligand_id", "Ligand_ID")
                    if c in fieldnames
                ),
                None,
            )
            stage_cols = [c for c in fieldnames if c.lower().startswith("stage")]
            for row in reader:
                ligand_raw = str(row.get(ligand_col or "", "")).strip()
                if not ligand_raw:
                    continue
                key = _chunk_ligand_key(ligand_raw)
                if not key:
                    continue
                if not stage_cols:
                    scored.add(key)
                    continue
                if any(str(row.get(col, "")).strip() for col in stage_cols):
                    scored.add(key)
    except Exception as exc:
        logging.getLogger("distributed.chunk").warning(
            "[distributed.chunk.verify.read] action=score_summary_parse_failed path=%s reason=%s",
            str(path),
            str(exc),
        )
    return scored


def _resolve_combo_post_consensus_csv(combo_post_dir: Path, library_name: str) -> Path:
    suffix = "_consensus_reranked_scorch.csv"
    default_path = combo_post_dir / "consensus_reranked_scorch.csv"
    prefixed = sorted(
        p
        for p in combo_post_dir.glob(f"*{suffix}")
        if p.name != "consensus_reranked_scorch.csv"
    )

    if default_path.exists() and default_path.stat().st_size > 0:
        return default_path

    if not prefixed:
        return default_path

    lib = str(library_name or "").strip()
    if lib:
        wanted = {lib.lower(), _chunk_ligand_key(lib).lower()}
        for candidate in prefixed:
            prefix = candidate.name[: -len(suffix)].strip().lower()
            if prefix in wanted:
                return candidate

    non_dud = [p for p in prefixed if not p.name.lower().startswith("dud_")]
    if non_dud:
        return non_dud[0]
    return prefixed[0]


def _load_post_scored_ligand_keys(path: Path) -> set[str]:
    scored: set[str] = set()
    if not path.exists():
        return scored

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        ligand_col = next(
            (
                c
                for c in ("ligand_base", "ligand", "Ligand", "ligand_file", "Ligand_ID")
                if c in fieldnames
            ),
            None,
        )
        score_col = next(
            (c for c in ("final_score", "consensus_score", "score") if c in fieldnames),
            None,
        )
        for row in reader:
            ligand_raw = str(row.get(ligand_col or "", "")).strip()
            if not ligand_raw:
                continue
            key = _chunk_ligand_key(ligand_raw)
            if not key:
                continue
            if score_col is not None:
                if not str(row.get(score_col, "")).strip():
                    continue
            scored.add(key)
    return scored


def _coverage_snapshot_path(
    cfg: Mapping[str, Any],
    *,
    run_id: str,
    pdb_id: str,
    variant_label: str,
    ph_tag: Optional[str],
) -> Path:
    variant_token = str(variant_label or "BASE").strip().upper() or "BASE"
    ph_token = _normalize_ph_tag_token(ph_tag)
    filename = f"{str(pdb_id).upper()}__{variant_token}__{ph_token}.json"
    return (
        runtime_root(cfg, "MANIFESTS_DIR", "manifests")
        / str(run_id)
        / "coverage"
        / filename
    )


def _load_keyset(payload: Mapping[str, Any], key: str) -> set[str]:
    raw = payload.get(key)
    if not isinstance(raw, list):
        return set()
    out: set[str] = set()
    for item in raw:
        val = _chunk_ligand_key(str(item))
        if val:
            out.add(val)
    return out


def _keys_sha256(keys: set[str]) -> str:
    raw = "\n".join(sorted(keys))
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _resolve_run_scoped_root(root: Path, run_id: str) -> Path:
    token = str(run_id or "").strip()
    if token and root.name != token:
        return root / token
    return root


def _resolve_combo_output_dir(
    cfg: Mapping[str, Any],
    *,
    root_key: str,
    default_name: str,
    run_id: str,
    pdb_id: str,
    variant_label: str,
    ph_tag: Optional[str],
) -> Path:
    repo_root = Path(__file__).resolve().parents[2]
    base_root = Path(
        str(cfg.get(root_key, repo_root / default_name) or repo_root / default_name)
    ).resolve()
    combo_dir = _resolve_run_scoped_root(base_root, str(run_id)) / str(pdb_id).upper()
    variant_token = str(variant_label or "").strip().upper()
    if variant_token and variant_token not in {"LEGACY", "BASE"}:
        combo_dir = combo_dir / variant_token
    if ph_tag:
        combo_dir = combo_dir / str(ph_tag)
    return combo_dir


def _update_combo_coverage_snapshot(
    cfg: Mapping[str, Any],
    *,
    run_id: str,
    pdb_id: str,
    variant_label: str,
    ph_tag: Optional[str],
    library_name: str,
    expected_keys: Optional[set[str]] = None,
    docking_keys: Optional[set[str]] = None,
    post_keys: Optional[set[str]] = None,
) -> None:
    path = _coverage_snapshot_path(
        cfg,
        run_id=run_id,
        pdb_id=pdb_id,
        variant_label=variant_label,
        ph_tag=ph_tag,
    )
    path.parent.mkdir(parents=True, exist_ok=True)

    payload: dict[str, Any] = {}
    if path.exists():
        try:
            loaded = json.loads(path.read_text(encoding="utf-8")) or {}
            if isinstance(loaded, dict):
                payload = loaded
        except Exception:
            payload = {}

    expected = _load_keyset(payload, "expected_ligands")
    docking = _load_keyset(payload, "docking_scored_ligands")
    post = _load_keyset(payload, "post_scored_ligands")

    if expected_keys is not None:
        expected.update({k for k in expected_keys if str(k).strip()})
    if docking_keys is not None:
        docking.update({k for k in docking_keys if str(k).strip()})
    if post_keys is not None:
        post.update({k for k in post_keys if str(k).strip()})

    if expected:
        docking &= expected
        post &= expected

    out = {
        "run_id": str(run_id),
        "pdb_id": str(pdb_id).upper(),
        "variant": str(variant_label or "").strip().upper(),
        "ph": _normalize_ph_tag_token(ph_tag),
        "library": str(library_name or ""),
        "expected_ligands": sorted(expected),
        "docking_scored_ligands": sorted(docking),
        "post_scored_ligands": sorted(post),
        "expected_count": len(expected),
        "docking_count": len(docking),
        "post_count": len(post),
        "expected_sha256": _keys_sha256(expected),
        "docking_sha256": _keys_sha256(docking),
        "post_sha256": _keys_sha256(post),
        "generated_at": datetime.datetime.utcnow().replace(microsecond=0).isoformat()
        + "Z",
    }

    tmp_path: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=path.name + ".",
            suffix=".tmp",
            dir=path.parent,
            delete=False,
        ) as handle:
            tmp_path = Path(handle.name)
            handle.write(json.dumps(out, indent=2, sort_keys=False))
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_path, path)
    finally:
        if tmp_path is not None and tmp_path.exists():
            try:
                tmp_path.unlink()
            except Exception:
                pass


def _verify_chunk_combo_outputs(
    cfg: Mapping[str, Any],
    *,
    run_id: str,
    chunk_id: Optional[str] = None,
    pdb_id: str,
    variant_label: str,
    ph_tag: Optional[str],
    library_name: str,
    chunk_ligand_bases: list[str],
    run_mode: Optional[str] = None,
) -> tuple[bool, str]:
    chunk_keys = {
        _chunk_ligand_key(str(x))
        for x in chunk_ligand_bases
        if str(x).strip()
    }
    if not chunk_keys:
        return False, "empty_chunk_keys"

    combo_dir = _resolve_combo_output_dir(
        cfg,
        root_key="DOCKED_DIR",
        default_name="docked",
        run_id=str(run_id),
        pdb_id=str(pdb_id),
        variant_label=str(variant_label),
        ph_tag=ph_tag,
    )
    marker_missing_reason = ""
    if chunk_id:
        stage_prefix = _csv_prefix_for_run_mode(run_mode) or ""
        stage1_dir = combo_dir / f"{stage_prefix}stage1"
        stage1_marker = completion_marker_path(
            stage1_dir,
            engine="vina",
            chunk_id=str(chunk_id),
        )
        if not stage1_marker.exists():
            marker_missing_reason = f"missing_stage1_chunk_marker:{stage1_marker}"

    summary_csv = _resolve_combo_docking_summary_csv(
        combo_dir,
        str(library_name or ""),
        run_mode=run_mode,
    )
    long_csv = _resolve_combo_docking_long_csv(
        combo_dir,
        str(library_name or ""),
        run_mode=run_mode,
    )
    if not summary_csv.exists():
        if long_csv.exists():
            scored_keys = _load_scored_ligand_keys_from_summary(long_csv)
        elif marker_missing_reason:
            return (
                False,
                f"retryable_{marker_missing_reason}:missing_docking_summary:{summary_csv}",
            )
        else:
            return True, f"degraded_missing_docking_summary:{summary_csv}"
    else:
        scored_keys = _load_scored_ligand_keys_from_summary(summary_csv)
        if long_csv.exists() and long_csv != summary_csv:
            scored_keys |= _load_scored_ligand_keys_from_summary(long_csv)

    missing_scored = chunk_keys - scored_keys
    if missing_scored:
        sample = ",".join(sorted(missing_scored)[:6])
        if chunk_id:
            return (
                False,
                "missing_chunk_scores:"
                f"missing={len(missing_scored)}:"
                f"expected={len(chunk_keys)}:"
                f"scored={len(scored_keys)}:"
                f"sample={sample}",
            )
        return (
            True,
            (
                f"degraded_{marker_missing_reason}:"
                if marker_missing_reason
                else "degraded_missing_chunk_scores:"
            )
            + (
                "missing="
                f"{len(missing_scored)}:"
                f"expected={len(chunk_keys)}:"
                f"scored={len(scored_keys)}:"
                f"sample={sample}"
            ),
        )
    if marker_missing_reason:
        return (
            True,
            f"degraded_{marker_missing_reason}:scored={len(scored_keys)}:expected={len(chunk_keys)}",
        )
    return True, "ok"
