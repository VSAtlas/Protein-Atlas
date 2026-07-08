# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
import csv
import logging
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from chemdb.decoy_control import (
    decoy_role_token,
    row_has_explicit_decoy,
    row_is_control,
)
from config.runtime_config import load_config
from config.output_paths import run_output_dir
from post_docking.rescoring import scorch_prefix_and_stages as _scorch_prefix_mod
from post_docking.rescoring.scorch_coverage import (
    normalize_ph as _normalize_scope_ph,
    normalize_pdb as _normalize_scope_pdb,
    normalize_variant as _normalize_scope_variant,
)
from post_docking.rescoring import rescore_reranker_support as _support
from post_docking.rescoring.rescore_reranker_support import CnnPick, ScorchPick

COMPONENT = "[rescore-reranker]"
SCORCH_WEIGHT_DEFAULT = 0.65
CNN_WEIGHT_DEFAULT = 0.35
SCORCH_WEIGHT_KEY = "SCORCH_WEIGHT"
CNN_WEIGHT_KEY = "CNN_WEIGHT"
DECOY_PREFIX_KEY = "DECOY_PREFIX"
DUD_PREFIX_KEY = "DUD_PREFIX"
DECOY_PREFIX_DEFAULT = "dud"
DECOY_PREFIX_VALUE = DECOY_PREFIX_DEFAULT

NUMERIC_PRIORITY = _support.NUMERIC_PRIORITY
TEXT_PRIORITY = _support.TEXT_PRIORITY
find_consensus_csv = _support.find_consensus_csv
find_dud_consensus_csv = _support.find_dud_consensus_csv
_read_csv = _support._read_csv
_write_csv = _support._write_csv

_DECOY_RE = re.compile(r"\bdecoys?_", re.IGNORECASE)


def _read_consensus_rows_with_chunk_union(
    consensus_csv: Path, logger: logging.Logger
) -> Tuple[List[Dict[str, str]], List[str]]:
    return _support._read_consensus_rows_with_chunk_union(
        consensus_csv,
        logger,
        component=COMPONENT,
        consensus_ligand_base=_consensus_ligand_base,
    )


def _ordered_fields(original_fields: List[str]) -> List[str]:
    return _support._ordered_fields(
        original_fields,
        numeric_priority=NUMERIC_PRIORITY,
        text_priority=TEXT_PRIORITY,
    )


def _sort_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    return _support._sort_rows(rows, as_float=_as_float)


def _decoy_prefix_value() -> str:
    return DECOY_PREFIX_VALUE or DECOY_PREFIX_DEFAULT


def _set_decoy_prefix(
    value: Optional[object], logger: Optional[logging.Logger] = None
) -> str:
    global DECOY_PREFIX_VALUE
    DECOY_PREFIX_VALUE = _scorch_prefix_mod.normalize_decoy_prefix(
        value,
        DECOY_PREFIX_DEFAULT,
    )
    if logger:
        logger.info("%s action=preflight decoy_prefix=%s", COMPONENT, DECOY_PREFIX_VALUE)
    return DECOY_PREFIX_VALUE


def _starts_with_decoy_prefix(text: str) -> bool:
    prefix = _decoy_prefix_value().lower()
    if not prefix:
        return False
    return str(text or "").strip().lower().startswith(f"{prefix}_")


def _is_decoy_scorch_path(path: Path, prefix: str) -> bool:
    name = path.name.lower()
    prefix_lower = prefix.lower()
    if name.startswith(f"{prefix_lower}_"):
        return True
    return prefix_lower != "dud" and name.startswith("dud_")


def _as_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        s = str(x).strip()
        if not s:
            return None
        return float(s)
    except Exception:
        return None


_Z_T_ALIAS_FIELDS: Tuple[Tuple[str, str], ...] = (
    ("z_vs_decoys_consensus", "t_vs_decoys_consensus"),
    ("z_vs_decoys_consensus_pre", "t_vs_decoys_consensus_pre"),
    ("z_vs_decoys_blend", "t_vs_decoys_blend"),
)


def _sync_z_t_aliases(row: Dict[str, Any]) -> None:
    for z_col, t_col in _Z_T_ALIAS_FIELDS:
        z_val = str(row.get(z_col, "")).strip()
        t_val = str(row.get(t_col, "")).strip()
        if z_val:
            row[z_col] = z_val
            row[t_col] = z_val
        elif t_val:
            row[z_col] = t_val
            row[t_col] = t_val
        else:
            row[z_col] = ""
            row[t_col] = ""


def is_decoy_id(s: str) -> bool:
    return bool(_DECOY_RE.search(s or "")) or _starts_with_decoy_prefix(s)


def is_decoy_like(text: str) -> bool:
    return bool(_DECOY_RE.search(text or "")) or _starts_with_decoy_prefix(text)


def _parse_weight(
    cfg: Dict[str, Any], key: str, default: float, logger: logging.Logger
) -> float:
    raw = cfg.get(key)
    if raw is None:
        raw = cfg.get(key.lower())
    if raw is None:
        return default
    try:
        return float(raw)
    except Exception as exc:
        logger.warning(
            "%s action=preflight status=degraded reason=invalid_weight key=%s value=%s error=%s default=%.3f",
            COMPONENT,
            key,
            raw,
            exc,
            default,
        )
        return default


def _load_blend_weights(
    repo_root: Path, logger: logging.Logger
) -> Tuple[float, float, str]:
    cfg: Dict[str, Any] = {}
    try:
        cfg = load_config(config_path=str(repo_root / "config.txt"), base_dir=repo_root)
    except Exception as exc:
        logger.warning(
            "%s action=preflight status=degraded reason=config_load_failed error=%s",
            COMPONENT,
            exc,
        )
        cfg = {}
    scorch_weight = _parse_weight(cfg, SCORCH_WEIGHT_KEY, SCORCH_WEIGHT_DEFAULT, logger)
    cnn_weight = _parse_weight(cfg, CNN_WEIGHT_KEY, CNN_WEIGHT_DEFAULT, logger)
    decoy_prefix = _scorch_prefix_mod.resolve_decoy_prefix_from_config(
        cfg,
        decoy_prefix_key=DECOY_PREFIX_KEY,
        dud_prefix_key=DUD_PREFIX_KEY,
        default_prefix=DECOY_PREFIX_DEFAULT,
    )
    _set_decoy_prefix(decoy_prefix, logger)
    logger.info(
        "%s action=preflight status=ok scorch_weight=%.3f cnn_weight=%.3f",
        COMPONENT,
        scorch_weight,
        cnn_weight,
    )
    return scorch_weight, cnn_weight, decoy_prefix


def _norm_engine(x: Any) -> str:
    return str(x or "").strip().lower()


def _extract_stage_num(s: str) -> int:
    m = re.search(r"stage(\d+)", s)
    if not m:
        return -1
    try:
        return int(m.group(1))
    except Exception:
        return -1


def canonical_ligand_base(lig: str) -> str:
    return _consensus_ligand_base(lig)


def _consensus_ligand_base(lig: str) -> str:
    return _scorch_ligand_base(lig)


def _scorch_ligand_base(lig_id: str) -> str:
    s = str(lig_id or "").strip()
    s = s.replace(".sanitized", "")
    prefixes = {_decoy_prefix_value()}
    if "dud" not in prefixes:
        prefixes.add("dud")
    for prefix in prefixes:
        esc = re.escape(prefix)
        s = re.sub(rf"(_{esc}_gnina_stage\d+)$", "", s)
        s = re.sub(rf"(_gnina_{esc}_stage\d+)$", "", s)
        s = re.sub(rf"(_dock6_{esc}_stage\d+)$", "", s)
        s = re.sub(rf"(_{esc}_dock6_stage\d+)$", "", s)
        s = re.sub(rf"(_{esc}_stage\d+)$", "", s)
        s = re.sub(rf"(__{esc}_dock6_stage\d+)$", "", s)
        s = re.sub(rf"(__dock6_{esc}_stage\d+)$", "", s)
        s = re.sub(rf"(__{esc}_ledock_stage\d+)$", "", s)
        s = re.sub(rf"(__ledock_{esc}_stage\d+)$", "", s)
    s = re.sub(r"(__ledock_stage\d+)$", "", s)
    s = re.sub(r"(__dock6_stage\d+)$", "", s)
    s = re.sub(r"(_gnina_stage\d+)$", "", s)
    s = re.sub(r"(_stage\d+)$", "", s)
    s = re.sub(r"\.(mol2|pdbqt)$", "", s, flags=re.IGNORECASE)
    s = s.replace("__", "_")
    s = re.sub(r"_+$", "", s)
    return s


def _row_identifier(row: Dict[str, Any]) -> str:
    for key in (
        "ligand",
        "ligand_file",
        "Ligand_ID",
        "ligand_id",
        "Ligand",
        "name",
        "molecule",
        "Molecule",
    ):
        val = row.get(key)
        if val:
            return str(val).strip()
    for key, val in row.items():
        if "ligand" in str(key).lower() and val:
            return str(val).strip()
    return ""


def _is_dud_row(row: Dict[str, Any]) -> bool:
    if row_is_control(row):
        return False
    if row_has_explicit_decoy(row):
        return True
    if decoy_role_token(row.get("run_mode")):
        return True
    if decoy_role_token(row.get("library")):
        return True
    return is_decoy_id(_row_identifier(row))


def _row_ligand_base(row: Dict[str, Any]) -> str:
    ident = _row_identifier(row)
    return _scorch_ligand_base(ident) if ident else ""


def _decoy_identifier(row: Dict[str, Any]) -> str:
    for key in ("ligand", "ligand_file", "Ligand_ID", "ligand_base"):
        val = row.get(key)
        if val:
            return str(val).strip()
    return _row_identifier(row)


def _read_consensus_scores(path: Path, logger: logging.Logger) -> List[float]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    try:
        with path.open("r", encoding="utf-8") as handle:
            lines = [
                line
                for line in handle
                if line.strip() and not line.lstrip().startswith("#")
            ]
        if not lines:
            return []
        reader = csv.DictReader(lines)
        fieldnames = reader.fieldnames or []
        if "consensus_score" not in fieldnames:
            logger.warning(
                "%s action=consensus status=skip reason=missing_consensus_score path=%s",
                COMPONENT,
                path,
            )
            return []
        scores: List[float] = []
        for row in reader:
            val = _as_float(row.get("consensus_score"))
            if val is not None and math.isfinite(val):
                scores.append(val)
        return scores
    except Exception as exc:
        logger.warning(
            "%s action=consensus status=skip reason=read_error path=%s error=%s",
            COMPONENT,
            path,
            exc,
        )
        return []


def _cnn_ligand_base(lig_id: str) -> str:
    return _scorch_ligand_base(lig_id)


def _find_repo_root(consensus_csv: Path) -> Optional[Path]:
    for parent in consensus_csv.parents:
        if parent.name == "docked":
            try:
                return parent.parent
            except Exception:
                return None
    return None




def _percentile_map(values: List[float]) -> Dict[float, float]:
    if not values:
        return {}
    if len(values) == 1:
        return {values[0]: 1.0}
    sorted_vals = sorted(values, reverse=True)
    n = len(sorted_vals)
    pct_by_value: Dict[float, float] = {}
    i = 0
    while i < n:
        v = sorted_vals[i]
        j = i
        while j + 1 < n and sorted_vals[j + 1] == v:
            j += 1
        # Average ties: map equal values to the mean rank they span.
        avg_rank = (i + j) / 2.0 + 1.0
        pct_by_value[v] = 1.0 - (avg_rank - 1.0) / (n - 1.0)
        i = j + 1
    return pct_by_value


def rerank_consensus_with_scorch(
    consensus_csv: Path,
    scorch_csvs: List[Path] | Tuple[Path, ...] | Path,
    out_csv: Path,
    logger: logging.Logger,
    *,
    overwrite: bool = False,
    scorch_weight: float = SCORCH_WEIGHT_DEFAULT,
    cnn_weight: float = CNN_WEIGHT_DEFAULT,
    decoy_prefix: Optional[str] = None,
) -> bool:
    if out_csv.exists() and out_csv.stat().st_size > 0 and not overwrite:
        logger.info("%s action=skip reason=exists out=%s", COMPONENT, str(out_csv))
        return True

    if not consensus_csv.exists() or consensus_csv.stat().st_size == 0:
        logger.warning(
            "%s action=skip reason=missing_consensus path=%s",
            COMPONENT,
            str(consensus_csv),
        )
        return False
    prefix = _scorch_prefix_mod.normalize_decoy_prefix(
        decoy_prefix if decoy_prefix is not None else _decoy_prefix_value()
    )
    if decoy_prefix is not None:
        _set_decoy_prefix(prefix)
    if isinstance(scorch_csvs, Path):
        scorch_paths = [scorch_csvs]
    else:
        scorch_paths = list(scorch_csvs)

    if len(scorch_paths) == 1 and scorch_paths[0].name == "scorch_scores_all.csv":
        decoy_candidate = scorch_paths[0].with_name(f"{prefix}_scorch_scores_all.csv")
        if decoy_candidate.exists():
            scorch_paths.append(decoy_candidate)
        if prefix.lower() != "dud":
            legacy_candidate = scorch_paths[0].with_name("dud_scorch_scores_all.csv")
            if legacy_candidate.exists() and legacy_candidate not in scorch_paths:
                scorch_paths.append(legacy_candidate)

    sc_rows_all: List[Dict[str, str]] = []
    counts: Dict[str, int] = {"fda": 0, "dud": 0}
    for path in scorch_paths:
        if not path.exists() or path.stat().st_size == 0:
            logger.warning(
                "%s action=rerank status=degraded reason=missing_or_empty_scorch path=%s",
                COMPONENT,
                str(path),
            )
            continue
        try:
            rows, _ = _read_csv(path)
        except Exception as exc:
            logger.warning(
                "%s action=rerank status=degraded reason=read_error_scorch path=%s error=%s",
                COMPONENT,
                str(path),
                exc,
            )
            continue
        mode_hint = "dud" if _is_decoy_scorch_path(path, prefix) else "fda"
        for row in rows:
            if not row.get("run_mode"):
                row["run_mode"] = mode_hint
            counts[row.get("run_mode", mode_hint)] = (
                counts.get(row.get("run_mode", mode_hint), 0) + 1
            )
            sc_rows_all.append(row)
    logger.info(
        "[rerank.load_scorch] fda_exists=%s dud_exists=%s n_fda=%d n_dud=%d",
        str(counts.get("fda", 0) > 0).lower(),
        str(counts.get("dud", 0) > 0).lower(),
        counts.get("fda", 0),
        counts.get("dud", 0),
    )

    base_dir = scorch_paths[0].parent if scorch_paths else out_csv.parent
    cnn_csv = base_dir / "cnn_rescoring.csv"
    cnn_rows: List[Dict[str, str]] = []
    if not cnn_csv.exists() or cnn_csv.stat().st_size == 0:
        logger.warning(
            "%s action=rerank status=degraded reason=missing_or_empty_cnn path=%s",
            COMPONENT,
            str(cnn_csv),
        )
    else:
        try:
            cnn_rows, _ = _read_csv(cnn_csv)
            if not cnn_rows:
                logger.warning(
                    "%s action=rerank status=degraded reason=missing_or_empty_cnn path=%s",
                    COMPONENT,
                    str(cnn_csv),
                )
        except Exception as exc:
            logger.warning(
                "%s action=rerank status=degraded reason=read_error_cnn path=%s error=%s",
                COMPONENT,
                str(cnn_csv),
                exc,
            )
            cnn_rows = []

    cons_rows, cons_fields = _read_consensus_rows_with_chunk_union(consensus_csv, logger)
    if not cons_rows:
        logger.warning(
            "%s action=skip reason=empty_consensus path=%s",
            COMPONENT,
            str(consensus_csv),
        )
        return False

    first_consensus = cons_rows[0]
    default_pdb_id = _normalize_scope_pdb(first_consensus.get("pdb_id", ""))
    default_variant = _normalize_scope_variant(first_consensus.get("variant", ""))
    default_ph = _normalize_scope_ph(first_consensus.get("ph_label", ""))

    cons_bases: Set[str] = set()
    for row in cons_rows:
        base = _consensus_ligand_base(str(row.get("ligand", "")).strip())
        if base:
            cons_bases.add(base)
    run_id_default = str(first_consensus.get("run_id", "")).strip()

    best_by_source: Dict[Tuple[str, str, str, str, str], ScorchPick] = {}
    best_any: Dict[Tuple[str, str, str, str], ScorchPick] = {}

    for r in sc_rows_all:
        pdb_id = _normalize_scope_pdb(r.get("pdb_id", "")) or default_pdb_id
        variant = _normalize_scope_variant(r.get("variant", ""))
        if not variant and not str(r.get("variant", "")).strip():
            variant = default_variant
        ph = _normalize_scope_ph(r.get("ph", "")) or default_ph
        source = _norm_engine(r.get("source", ""))
        lig_id = str(r.get("Ligand_ID", "")).strip()
        if not lig_id:
            continue

        lig_base = _scorch_ligand_base(lig_id)
        stage_num = _extract_stage_num(lig_id)
        score = _as_float(r.get("SCORCH_score"))
        certainty = _as_float(r.get("SCORCH_certainty"))
        run_mode = str(r.get("run_mode", "")).strip() or "fda"
        selected_stage = str(r.get("selected_stage", "")).strip()
        rescored_stage = str(r.get("rescored_stage", "")).strip()
        selected_docking_score = _as_float(r.get("selected_docking_score"))
        stage_match_flag = str(r.get("stage_match_flag", "")).strip()
        stage_fallback_reason = str(r.get("stage_fallback_reason", "")).strip()

        sc_pick = ScorchPick(
            pdb_id=pdb_id,
            variant=variant,
            ph=ph,
            source=source,
            ligand_id=lig_id,
            ligand_base=lig_base,
            stage_num=stage_num,
            scorch_score=score,
            scorch_certainty=certainty,
            run_mode=run_mode,
            selected_stage=selected_stage,
            rescored_stage=rescored_stage,
            selected_docking_score=selected_docking_score,
            stage_match_flag=stage_match_flag,
            stage_fallback_reason=stage_fallback_reason,
        )

        k1 = (pdb_id, variant, ph, source, lig_base)
        cur = best_by_source.get(k1)
        if cur is None or (
            1 if sc_pick.stage_match_flag == "1" else 0,
            sc_pick.scorch_score or -1e9,
            sc_pick.stage_num,
        ) > (
            1 if cur.stage_match_flag == "1" else 0,
            cur.scorch_score or -1e9,
            cur.stage_num,
        ):
            best_by_source[k1] = sc_pick

        k2 = (pdb_id, variant, ph, lig_base)
        cur2 = best_any.get(k2)
        if cur2 is None or (
            1 if sc_pick.stage_match_flag == "1" else 0,
            sc_pick.scorch_score or -1e9,
            sc_pick.stage_num,
        ) > (
            1 if cur2.stage_match_flag == "1" else 0,
            cur2.scorch_score or -1e9,
            cur2.stage_num,
        ):
            best_any[k2] = sc_pick

    best_cnn: Dict[str, CnnPick] = {}
    for r in cnn_rows:
        lig_id = str(r.get("ligand", "")).strip()
        if not lig_id:
            continue
        lig_base = _cnn_ligand_base(lig_id)
        cnn_score = _as_float(r.get("cnn_score"))
        cnn_affinity = _as_float(r.get("cnn_affinity"))
        if cnn_score is None or cnn_affinity is None:
            continue
        cnn_vs = cnn_score * cnn_affinity
        stage_priority = str(r.get("stage_priority", "")).strip().lower()
        stage_num = _extract_stage_num(stage_priority)

        cnn_pick = CnnPick(
            ligand_base=lig_base,
            stage_num=stage_num,
            cnn_score=cnn_score,
            cnn_affinity=cnn_affinity,
            cnn_vs=cnn_vs,
        )
        current_cnn = best_cnn.get(lig_base)
        if current_cnn is None:
            best_cnn[lig_base] = cnn_pick
        else:
            cur_vs = (
                current_cnn.cnn_vs if current_cnn.cnn_vs is not None else -1e18
            )
            if (
                cnn_vs > cur_vs
                or (cnn_vs == cur_vs and stage_num > current_cnn.stage_num)
            ):
                best_cnn[lig_base] = cnn_pick

    enriched: List[Dict[str, Any]] = []
    for r in cons_rows:
        rr: Dict[str, Any] = dict(r)
        lig = str(r.get("ligand", "")).strip()
        lig_base = _consensus_ligand_base(lig)
        rr["ligand_base"] = lig_base
        library_raw = str(r.get("library", "")).strip()
        rr["library"] = library_raw
        rr["consensus_score_pre"] = str(r.get("consensus_score", "")).strip()
        z_consensus = str(
            r.get("z_vs_decoys_consensus", r.get("t_vs_decoys_consensus", ""))
        ).strip()
        z_consensus_pre = str(
            r.get(
                "z_vs_decoys_consensus_pre",
                r.get("t_vs_decoys_consensus_pre", z_consensus),
            )
        ).strip()
        rr["z_vs_decoys_consensus_pre"] = z_consensus_pre
        rr["z_vs_decoys_consensus"] = str(
            rr.get("z_vs_decoys_consensus", z_consensus)
        ).strip()
        rr["consensus_mu_decoy"] = ""
        rr["consensus_sigma_decoy"] = ""
        rr["consensus_n_decoys"] = ""
        rr["z_vs_decoys_blend"] = ""
        rr["blend_mu_decoy"] = ""
        rr["blend_sigma_decoy"] = ""
        rr["blend_n_decoys"] = ""
        _sync_z_t_aliases(rr)

        pdb_id = _normalize_scope_pdb(r.get("pdb_id", ""))
        variant = _normalize_scope_variant(r.get("variant", ""))
        ph_label = _normalize_scope_ph(r.get("ph_label", ""))
        best_engine = _norm_engine(r.get("best_engine", ""))

        pick: Optional[ScorchPick] = None
        if pdb_id and best_engine:
            pick = best_by_source.get(
                (pdb_id, variant, ph_label, best_engine, lig_base)
            )
        if pick is None and pdb_id:
            pick = best_any.get((pdb_id, variant, ph_label, lig_base))

        run_mode_val = str(r.get("run_mode", "")).strip().lower()
        if not run_mode_val and pick is not None:
            run_mode_val = pick.run_mode
        rr["run_mode"] = run_mode_val

        rr["scorch_source_used"] = pick.source if pick else ""
        rr["SCORCH_score_used"] = (
            ""
            if (not pick or pick.scorch_score is None)
            else f"{pick.scorch_score:.6g}"
        )
        rr["SCORCH_certainty_used"] = (
            ""
            if (not pick or pick.scorch_certainty is None)
            else f"{pick.scorch_certainty:.6g}"
        )
        rr["selected_stage"] = pick.selected_stage if pick else ""
        rr["rescored_stage"] = pick.rescored_stage if pick else ""
        rr["selected_docking_score"] = (
            ""
            if (not pick or pick.selected_docking_score is None)
            else f"{pick.selected_docking_score:.6g}"
        )
        rr["stage_match_flag"] = pick.stage_match_flag if pick else ""
        rr["stage_fallback_reason"] = pick.stage_fallback_reason if pick else ""

        comp = None
        if pick and pick.scorch_score is not None and pick.scorch_certainty is not None:
            comp = pick.scorch_score * pick.scorch_certainty
        rr["scorch_composite"] = "" if comp is None else f"{comp:.6g}"

        cnn_pick_used = best_cnn.get(lig_base)
        rr["cnn_score_used"] = (
            ""
            if (not cnn_pick_used or cnn_pick_used.cnn_score is None)
            else f"{cnn_pick_used.cnn_score:.6g}"
        )
        rr["cnn_affinity_used"] = (
            ""
            if (not cnn_pick_used or cnn_pick_used.cnn_affinity is None)
            else f"{cnn_pick_used.cnn_affinity:.6g}"
        )
        rr["cnn_vs_used"] = (
            ""
            if (not cnn_pick_used or cnn_pick_used.cnn_vs is None)
            else f"{cnn_pick_used.cnn_vs:.6g}"
        )
        rr["cnn_rescored_flag"] = "1" if rr["cnn_vs_used"] != "" else "0"

        rescored = (
            (rr["scorch_composite"] != "")
            or (rr["SCORCH_score_used"] != "")
            or (rr["cnn_vs_used"] != "")
        )
        rr["rescored_flag"] = "1" if rescored else "0"
        if not rescored:
            rr["scorch_source_used"] = ""
            rr["SCORCH_score_used"] = ""
            rr["SCORCH_certainty_used"] = ""
            rr["scorch_composite"] = ""
            rr["selected_stage"] = ""
            rr["rescored_stage"] = ""
            rr["selected_docking_score"] = ""
            rr["stage_match_flag"] = ""
            rr["stage_fallback_reason"] = ""
            rr["scorch_pct"] = ""
            rr["cnn_pct"] = ""
            rr["ml_blend_score"] = ""
            rr["z_vs_decoys_blend"] = ""
        rr["final_score"] = ""
        _sync_z_t_aliases(rr)

        enriched.append(rr)

    for pick in best_any.values():
        if pick.ligand_base in cons_bases:
            continue
        extra_row: Dict[str, Any] = {}
        extra_row["run_id"] = run_id_default
        extra_row["pdb_id"] = pick.pdb_id
        extra_row["variant"] = pick.variant
        extra_row["ph_label"] = pick.ph
        extra_row["ligand"] = pick.ligand_id
        extra_row["ligand_base"] = pick.ligand_base
        extra_row["library"] = ""
        extra_row["consensus_score_pre"] = ""
        extra_row["z_vs_decoys_consensus_pre"] = ""
        extra_row["run_mode"] = pick.run_mode
        extra_row["z_vs_decoys_consensus"] = ""
        extra_row["consensus_mu_decoy"] = ""
        extra_row["consensus_sigma_decoy"] = ""
        extra_row["consensus_n_decoys"] = ""
        extra_row["z_vs_decoys_blend"] = ""
        extra_row["blend_mu_decoy"] = ""
        extra_row["blend_sigma_decoy"] = ""
        extra_row["blend_n_decoys"] = ""
        extra_row["best_engine"] = pick.source
        extra_row["scorch_source_used"] = pick.source
        extra_row["SCORCH_score_used"] = (
            "" if pick.scorch_score is None else f"{pick.scorch_score:.6g}"
        )
        extra_row["SCORCH_certainty_used"] = (
            "" if pick.scorch_certainty is None else f"{pick.scorch_certainty:.6g}"
        )
        extra_row["selected_stage"] = pick.selected_stage
        extra_row["rescored_stage"] = pick.rescored_stage
        extra_row["selected_docking_score"] = (
            ""
            if pick.selected_docking_score is None
            else f"{pick.selected_docking_score:.6g}"
        )
        extra_row["stage_match_flag"] = pick.stage_match_flag
        extra_row["stage_fallback_reason"] = pick.stage_fallback_reason
        comp = None
        if pick.scorch_score is not None and pick.scorch_certainty is not None:
            comp = pick.scorch_score * pick.scorch_certainty
        extra_row["scorch_composite"] = "" if comp is None else f"{comp:.6g}"
        extra_row["cnn_score_used"] = ""
        extra_row["cnn_affinity_used"] = ""
        extra_row["cnn_vs_used"] = ""
        extra_row["cnn_rescored_flag"] = "0"
        extra_row["rescored_flag"] = (
            "1" if extra_row["scorch_composite"] or extra_row["SCORCH_score_used"] else "0"
        )
        extra_row["final_score"] = ""
        _sync_z_t_aliases(extra_row)
        enriched.append(extra_row)

    def group_key(row: Dict[str, Any]) -> Tuple[str, str, str, str]:
        return (
            str(row.get("run_id", "")).strip(),
            str(row.get("pdb_id", "")).strip(),
            str(row.get("variant", "")).strip(),
            str(row.get("ph_label", "")).strip(),
        )

    groups: Dict[Tuple[str, str, str, str], List[int]] = {}
    for i, row in enumerate(enriched):
        groups.setdefault(group_key(row), []).append(i)

    for _, idxs in groups.items():
        scored = []
        for i in idxs:
            cs = _as_float(enriched[i].get("consensus_score"))
            scored.append((cs if cs is not None else -1e18, i))
        scored.sort(key=lambda x: x[0], reverse=True)
        for rank, (_, i) in enumerate(scored, start=1):
            enriched[i]["consensus_rank"] = str(rank)

    dud_rows_for_stats = [row for row in enriched if _is_dud_row(row)]
    decoy_consensus_csv = out_csv.with_name(f"{prefix}_consensus_reranked_scorch.csv")
    legacy_decoy_consensus_csv = out_csv.with_name("dud_consensus_reranked_scorch.csv")
    decoy_scores: List[float] = []
    mu_decoy: Optional[float] = None
    sigma_decoy: Optional[float] = None
    if dud_rows_for_stats:
        consensus_source = "dud_consensus_reranked"
        for row in dud_rows_for_stats:
            val = _as_float(row.get("consensus_score"))
            if val is None or not math.isfinite(val):
                continue
            ident = _decoy_identifier(row)
            if is_decoy_like(ident):
                decoy_scores.append(val)
        if not decoy_scores:
            consensus_source = "dud_role_unlabeled"
            for row in dud_rows_for_stats:
                val = _as_float(row.get("consensus_score"))
                if val is not None and math.isfinite(val):
                    decoy_scores.append(val)
    if not decoy_scores:
        consensus_source = "filename_detection"
        for row in enriched:
            val = _as_float(row.get("consensus_score"))
            if val is None or not math.isfinite(val):
                continue
            ident = _decoy_identifier(row)
            if is_decoy_like(ident):
                decoy_scores.append(val)

    if not decoy_scores:
        logger.info(
            "[z-score.consensus.skip] reason=no_decoy_scores_or_sigma0 n_decoys=0 sigma=nan",
        )
        for row in enriched:
            row["z_vs_decoys_consensus"] = ""
            row["consensus_mu_decoy"] = ""
            row["consensus_sigma_decoy"] = ""
            row["consensus_n_decoys"] = ""
            _sync_z_t_aliases(row)
    else:
        mu_decoy = sum(decoy_scores) / len(decoy_scores)
        variance = sum((v - mu_decoy) ** 2 for v in decoy_scores) / len(decoy_scores)
        sigma_decoy = math.sqrt(variance)
        if (
            sigma_decoy <= 0
            or not math.isfinite(mu_decoy)
            or not math.isfinite(sigma_decoy)
        ):
            logger.info(
                "[z-score.consensus.skip] reason=no_decoy_scores_or_sigma0 n_decoys=%d sigma=%s",
                len(decoy_scores),
                "nan" if not math.isfinite(sigma_decoy) else f"{sigma_decoy:.6g}",
            )
            for row in enriched:
                row["z_vs_decoys_consensus"] = ""
                row["consensus_mu_decoy"] = ""
                row["consensus_sigma_decoy"] = ""
                row["consensus_n_decoys"] = ""
                _sync_z_t_aliases(row)
        else:
            logger.info(
                "[z-score.consensus] source=%s n_decoys=%d mu=%.6g sigma=%.6g",
                consensus_source,
                len(decoy_scores),
                mu_decoy,
                sigma_decoy,
            )
            for row in enriched:
                val = _as_float(row.get("consensus_score"))
                row["consensus_mu_decoy"] = f"{mu_decoy:.6g}"
                row["consensus_sigma_decoy"] = f"{sigma_decoy:.6g}"
                row["consensus_n_decoys"] = str(len(decoy_scores))
                if val is None or not math.isfinite(val):
                    row["z_vs_decoys_consensus"] = ""
                else:
                    z_score = (val - mu_decoy) / sigma_decoy
                    row["z_vs_decoys_consensus"] = f"{z_score:.6g}"
                _sync_z_t_aliases(row)

    for _, idxs in groups.items():
        scorch_values: List[float] = []
        cnn_values: List[float] = []
        for i in idxs:
            comp_val = _as_float(enriched[i].get("scorch_composite"))
            scs_val = _as_float(enriched[i].get("SCORCH_score_used"))
            scorch_val = comp_val if comp_val is not None else scs_val
            if scorch_val is not None:
                scorch_values.append(scorch_val)
            cnn_val = _as_float(enriched[i].get("cnn_vs_used"))
            if cnn_val is not None:
                cnn_values.append(cnn_val)

        scorch_pct_map = _percentile_map(scorch_values)
        cnn_pct_map = _percentile_map(cnn_values)

        for i in idxs:
            rescored = str(enriched[i].get("rescored_flag", "")).strip() == "1"
            if not rescored:
                enriched[i]["scorch_pct"] = ""
                enriched[i]["cnn_pct"] = ""
                enriched[i]["ml_blend_score"] = ""
                continue
            comp_val = _as_float(enriched[i].get("scorch_composite"))
            scs_val = _as_float(enriched[i].get("SCORCH_score_used"))
            scorch_val = comp_val if comp_val is not None else scs_val
            cnn_val = _as_float(enriched[i].get("cnn_vs_used"))
            if scorch_val is None and cnn_val is None:
                enriched[i]["scorch_pct"] = ""
                enriched[i]["cnn_pct"] = ""
                enriched[i]["ml_blend_score"] = ""
                continue
            if scorch_val is None or not scorch_pct_map:
                scorch_pct = 0.5
            else:
                scorch_pct = scorch_pct_map.get(scorch_val, 0.5)

            if cnn_val is None or not cnn_pct_map:
                cnn_pct = 0.5
            else:
                cnn_pct = cnn_pct_map.get(cnn_val, 0.5)

            ml_blend = scorch_weight * scorch_pct + cnn_weight * cnn_pct
            enriched[i]["scorch_pct"] = f"{scorch_pct:.6g}"
            enriched[i]["cnn_pct"] = f"{cnn_pct:.6g}"
            enriched[i]["ml_blend_score"] = f"{ml_blend:.6g}"

    decoy_stats: Dict[str, Dict[str, float]] = {}
    for g_key, idxs in groups.items():
        blend_decoy_scores: List[float] = []
        for i in idxs:
            score_val = _as_float(enriched[i].get("ml_blend_score"))
            if score_val is None or not math.isfinite(score_val):
                continue
            ident = _row_identifier(enriched[i])
            if is_decoy_id(ident):
                blend_decoy_scores.append(score_val)
        if not blend_decoy_scores:
            for i in idxs:
                if not _is_dud_row(enriched[i]):
                    continue
                score_val = _as_float(enriched[i].get("ml_blend_score"))
                if score_val is not None and math.isfinite(score_val):
                    blend_decoy_scores.append(score_val)

        if not blend_decoy_scores:
            logger.info(
                "[z-score.blend.skip] run_id=%s pdb=%s variant=%s ph=%s reason=no_decoy_scores_or_sigma0 n_decoys=0 sigma=nan",
                g_key[0],
                g_key[1],
                g_key[2],
                g_key[3],
            )
            continue

        mu = sum(blend_decoy_scores) / len(blend_decoy_scores)
        variance = sum((v - mu) ** 2 for v in blend_decoy_scores) / len(
            blend_decoy_scores
        )
        sigma = math.sqrt(variance)
        if sigma <= 0 or not math.isfinite(mu) or not math.isfinite(sigma):
            logger.info(
                "[z-score.blend.skip] run_id=%s pdb=%s variant=%s ph=%s reason=no_decoy_scores_or_sigma0 n_decoys=%d sigma=%s",
                g_key[0],
                g_key[1],
                g_key[2],
                g_key[3],
                len(blend_decoy_scores),
                "nan" if not math.isfinite(sigma) else f"{sigma:.6g}",
            )
            continue

        logger.info(
            "[z-score.blend] run_id=%s pdb=%s variant=%s ph=%s n_decoys=%d mu=%.6g sigma=%.6g",
            g_key[0],
            g_key[1],
            g_key[2],
            g_key[3],
            len(blend_decoy_scores),
            mu,
            sigma,
        )
        stats_key = "|".join(g_key)
        decoy_stats[stats_key] = {
            "n_decoys": len(blend_decoy_scores),
            "mu": mu,
            "sigma": sigma,
        }
        for i in idxs:
            val = _as_float(enriched[i].get("ml_blend_score"))
            enriched[i]["blend_mu_decoy"] = f"{mu:.6g}"
            enriched[i]["blend_sigma_decoy"] = f"{sigma:.6g}"
            enriched[i]["blend_n_decoys"] = str(len(blend_decoy_scores))
            if val is None or not math.isfinite(val):
                enriched[i]["z_vs_decoys_blend"] = ""
                continue
            z_score = (val - mu) / sigma
            enriched[i]["z_vs_decoys_blend"] = f"{z_score:.6g}"
            _sync_z_t_aliases(enriched[i])

    for row in enriched:
        blend_val = _as_float(
            row.get("z_vs_decoys_blend", row.get("t_vs_decoys_blend", ""))
        )
        cons_val = _as_float(
            row.get("z_vs_decoys_consensus", row.get("t_vs_decoys_consensus", ""))
        )
        if blend_val is not None and math.isfinite(blend_val):
            row["final_score"] = f"{blend_val:.6g}"
        elif cons_val is not None and math.isfinite(cons_val):
            row["final_score"] = f"{cons_val:.6g}"
        else:
            row["final_score"] = ""

    for _, idxs in groups.items():

        def _final_rank_key(i: int) -> Tuple[float, float, str]:
            fs_val = _as_float(enriched[i].get("final_score"))
            if fs_val is None or not math.isfinite(fs_val):
                fs_val = -1e18
            cs_val = _as_float(enriched[i].get("consensus_score_pre"))
            if cs_val is None or not math.isfinite(cs_val):
                cs_val = _as_float(enriched[i].get("consensus_score"))
            if cs_val is None or not math.isfinite(cs_val):
                cs_val = -1e18
            lig = str(enriched[i].get("ligand", "")).strip().lower()
            return (-fs_val, -cs_val, lig)

        ordered = sorted(idxs, key=_final_rank_key)
        for rank, i in enumerate(ordered, start=1):
            enriched[i]["final_rank"] = str(rank)

    for row in enriched:
        rescored = bool(
            str(row.get("scorch_composite", "")).strip()
            or str(row.get("SCORCH_score_used", "")).strip()
            or str(row.get("cnn_vs_used", "")).strip()
        )
        row["rescored_flag"] = "1" if rescored else "0"
        if not rescored:
            row["scorch_source_used"] = ""
            row["SCORCH_score_used"] = ""
            row["SCORCH_certainty_used"] = ""
            row["scorch_composite"] = ""
            row["selected_stage"] = ""
            row["rescored_stage"] = ""
            row["selected_docking_score"] = ""
            row["stage_match_flag"] = ""
            row["stage_fallback_reason"] = ""
            row["scorch_pct"] = ""
            row["cnn_pct"] = ""
            row["ml_blend_score"] = ""
            row["z_vs_decoys_blend"] = ""
        _sync_z_t_aliases(row)

    out_fields = _ordered_fields(cons_fields)
    pretty_fields = [
        f for f in out_fields if f not in {"run_id", "pdb_id", "variant", "ph_label"}
    ]
    sorted_rows = _sort_rows(enriched)

    _write_csv(out_csv, sorted_rows, out_fields)
    dud_rows = [row for row in sorted_rows if _is_dud_row(row)]
    _write_csv(decoy_consensus_csv, dud_rows, out_fields)
    if prefix.lower() != "dud" and decoy_consensus_csv != legacy_decoy_consensus_csv:
        _write_csv(legacy_decoy_consensus_csv, dud_rows, out_fields)

    meta = [
        f"# run_id={str(cons_rows[0].get('run_id', '')).strip()}\n",
        "# pdb_id={pdb}, variant={variant}, ph={ph}\n".format(
            pdb=str(cons_rows[0].get("pdb_id", "")).strip(),
            variant=str(cons_rows[0].get("variant", "")).strip(),
            ph=str(cons_rows[0].get("ph_label", "")).strip(),
        ),
    ]
    _write_csv(
        out_csv.with_suffix(".pretty.csv"),
        sorted_rows,
        pretty_fields,
        preamble_lines=meta,
    )
    if decoy_stats:
        stats_path = out_csv.with_name("decoy_stats_blend.json")
        try:
            with stats_path.open("w", encoding="utf-8") as handle:
                json.dump(decoy_stats, handle, indent=2)
        except Exception as exc:
            logger.warning(
                "%s action=write status=degraded reason=stats_write_failed path=%s error=%s",
                COMPONENT,
                stats_path,
                exc,
            )

    if (
        consensus_source
        and decoy_scores
        and mu_decoy is not None
        and sigma_decoy is not None
        and sigma_decoy > 0
    ):
        consensus_stats_path = out_csv.with_name("consensus_decoy_stats.json")
        decoy_csv_value = str(decoy_consensus_csv) if decoy_consensus_csv else None
        payload = {
            "n_decoys": len(decoy_scores),
            "mu_decoy": mu_decoy,
            "sigma_decoy": sigma_decoy,
            "source": consensus_source,
            "dud_consensus_csv": decoy_csv_value,
            "dud_consensus_reranked_csv": decoy_csv_value,
            "regular_consensus_csv": str(consensus_csv),
        }
        try:
            with consensus_stats_path.open("w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2)
        except Exception as exc:
            logger.warning(
                "%s action=write status=degraded reason=consensus_stats_write_failed path=%s error=%s",
                COMPONENT,
                consensus_stats_path,
                exc,
            )

    logger.info(
        "%s action=write status=ok consensus=%s scorch=%s out=%s rows=%d",
        COMPONENT,
        str(consensus_csv),
        ",".join(str(p) for p in scorch_paths),
        str(out_csv),
        len(enriched),
    )
    return True


def rerank_run(
    run_id: str,
    repo_root: Path,
    overwrite: bool,
    logger: logging.Logger,
    *,
    scorch_weight: float = SCORCH_WEIGHT_DEFAULT,
    cnn_weight: float = CNN_WEIGHT_DEFAULT,
    decoy_prefix: str = DECOY_PREFIX_DEFAULT,
) -> int:
    prefix = _scorch_prefix_mod.normalize_decoy_prefix(
        decoy_prefix,
        DECOY_PREFIX_DEFAULT,
    )
    _set_decoy_prefix(prefix)
    post_root = run_output_dir(repo_root, "post_docked", run_id)
    dock_root = run_output_dir(repo_root, "docked", run_id)
    if not post_root.exists():
        logger.warning(
            "%s action=discover status=skip reason=missing_post_root path=%s",
            COMPONENT,
            str(post_root),
        )
        return 1
    if not dock_root.exists():
        logger.warning(
            "%s action=discover status=skip reason=missing_dock_root path=%s",
            COMPONENT,
            str(dock_root),
        )
        return 1

    scorch_paths = list(post_root.glob("**/scorch_scores_all.csv"))
    decoy_paths = list(post_root.glob(f"**/{prefix}_scorch_scores_all.csv"))
    legacy_dud_paths = (
        list(post_root.glob("**/dud_scorch_scores_all.csv"))
        if prefix.lower() != "dud"
        else []
    )
    combo_dirs = {p.parent for p in scorch_paths + decoy_paths + legacy_dud_paths}
    if not combo_dirs:
        logger.warning(
            "%s action=discover status=skip reason=no_scorch_scores path=%s",
            COMPONENT,
            str(post_root),
        )
        return 1

    ok = 0
    fail = 0
    for combo_dir in sorted(combo_dirs):
        try:
            ph = combo_dir.name
            variant = combo_dir.parent.name
            pdb_id = combo_dir.parent.parent.name

            dock_combo_dir = dock_root / pdb_id / variant / ph
            consensus_csv = find_consensus_csv(dock_combo_dir)
            if not consensus_csv:
                logger.warning(
                    "%s action=skip reason=missing_consensus pdb=%s variant=%s ph=%s dock_dir=%s",
                    COMPONENT,
                    pdb_id,
                    variant,
                    ph,
                    str(dock_combo_dir),
                )
                fail += 1
                continue

            out_csv = combo_dir / "consensus_reranked_scorch.csv"
            scorch_inputs: List[Path] = []
            fda_csv = combo_dir / "scorch_scores_all.csv"
            decoy_csv = combo_dir / f"{prefix}_scorch_scores_all.csv"
            legacy_dud_csv = combo_dir / "dud_scorch_scores_all.csv"
            if fda_csv.exists():
                scorch_inputs.append(fda_csv)
            if decoy_csv.exists():
                scorch_inputs.append(decoy_csv)
            elif prefix.lower() != "dud" and legacy_dud_csv.exists():
                scorch_inputs.append(legacy_dud_csv)
            if not scorch_inputs:
                logger.warning(
                    "%s action=skip reason=missing_scorch pdb=%s variant=%s ph=%s dir=%s",
                    COMPONENT,
                    pdb_id,
                    variant,
                    ph,
                    combo_dir,
                )
                fail += 1
                continue
            if rerank_consensus_with_scorch(
                consensus_csv,
                scorch_inputs,
                out_csv,
                logger,
                overwrite=overwrite,
                scorch_weight=scorch_weight,
                cnn_weight=cnn_weight,
                decoy_prefix=prefix,
            ):
                ok += 1
            else:
                fail += 1
        except Exception:
            logger.warning(
                "%s action=rerank status=failed reason=exception path=%s",
                COMPONENT,
                str(combo_dir),
                exc_info=True,
            )
            fail += 1

    logger.info(
        "%s action=summary status=ok combos_ok=%d combos_failed=%d", COMPONENT, ok, fail
    )
    return 0 if ok > 0 and fail == 0 else (0 if ok > 0 else 1)


def _configure_logging(verbose: bool) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")
    return logging.getLogger("rescore-reranker")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Rerank consensus docking results using SCORCH rescoring outputs."
    )
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--repo-root", required=True)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    logger = _configure_logging(args.verbose)
    repo_root = Path(args.repo_root).resolve()
    scorch_weight, cnn_weight, decoy_prefix = _load_blend_weights(repo_root, logger)
    return rerank_run(
        args.run_id,
        repo_root,
        args.overwrite,
        logger,
        scorch_weight=scorch_weight,
        cnn_weight=cnn_weight,
        decoy_prefix=decoy_prefix,
    )


if __name__ == "__main__":
    sys.exit(main())
