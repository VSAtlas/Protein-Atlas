from __future__ import annotations

import argparse
import math
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, cast

import pandas as pd  # type: ignore[import-untyped]

from analysis.dud_eval_core.labels import parse_name_and_label
from analysis.dud_eval_core.log import dbg
from analysis.dud_eval_core.schema import guess_ligfile_col
from analysis.dud_eval_core.types import TargetEvaluation, TargetSpec
from src.path_router.path_router import make_paths

FULL_RUN_MIN_LIGANDS = int(os.environ.get("FULL_RUN_MIN_LIGANDS", 10))

CONTROL_CENTERS_RE = re.compile(
    r"\[control-centers\]\s+"
    r"n=(\d+)\s+"
    r"max\?=([0-9]+(?:\.[0-9]+)?)\s*A?\b\s+"
    r"policy=(.+)$",
    re.IGNORECASE,
)

CONTROL_REDOCK_RE = re.compile(
    r"\[control-redock\]\s+"
    r"lig=([^\s]+)\s+"
    r"rmsd=([0-9]+(?:\.[0-9]+)?)\s*A?\b\s+"
    r"score=([0-9.\-]+)",
    re.IGNORECASE,
)

_CONTROL_PATTERNS_LOGGED = False

def _compute_analysis_root(
    args: argparse.Namespace, cfg: Optional[dict], run_id_override: Optional[str] = None
) -> Path:
    """
    Resolve the final analysis root honoring OVERALL_DIR, --out-dir, and --run-id.
    Relative out_dir values are anchored under OVERALL_DIR when available.
    """
    out_root = Path(args.out_dir)
    analysis_root = out_root
    if cfg and "OVERALL_DIR" in cfg:
        base = Path(cfg["OVERALL_DIR"])
        analysis_root = out_root if out_root.is_absolute() else base / out_root
    run_id_val = (
        run_id_override
        if run_id_override is not None
        else getattr(args, "run_id", None)
    )
    if run_id_val:
        analysis_root = analysis_root / str(run_id_val)
    analysis_root.mkdir(parents=True, exist_ok=True)
    return analysis_root


def _format_run_label(run_id: Optional[str]) -> str:
    """
    Sanitize run_id for filenames. Falls back to 'none' when run_id is missing.
    """
    if not run_id:
        return "none"
    label = str(run_id).strip()
    for sep in (os.sep, os.altsep):
        if sep:
            label = label.replace(sep, "_")
    return label.replace(" ", "") or "none"


def select_default_run_id(
    targets: List[TargetSpec],
    csv_paths_by_target: Dict[str, Path],
    lig_col_cli: Optional[str],
    score_col_cli: Optional[str],
) -> Optional[str]:
    run_full_map: Dict[str, Set[str]] = {}
    run_mtimes: Dict[str, float] = {}

    digit_re = re.compile(r"^\d{8,}$")
    iso_re = re.compile(r"^\d{4}-\d{2}-\d{2}([Tt _].*)?$")

    def _parse_isoish(token: str) -> Optional[float]:
        text = token.strip()
        if not text:
            return None
        cleaned = text.rstrip("Z").rstrip("z")
        candidates = [cleaned]
        if "T" in cleaned:
            candidates.append(cleaned.replace("T", " "))
        if "_" in cleaned:
            candidates.append(cleaned.replace("_", " "))
        if cleaned.endswith("T"):
            candidates.append(cleaned.rstrip("T"))
        for cand in candidates:
            cand = cand.strip()
            if not cand:
                continue
            try:
                dt = datetime.fromisoformat(cand)
            except Exception:
                continue
            try:
                return dt.replace(tzinfo=dt.tzinfo or timezone.utc).timestamp()
            except Exception:
                try:
                    return dt.timestamp()
                except Exception:
                    continue
        return None

    for target in targets:
        csv_path = csv_paths_by_target.get(target.target_key)
        if not csv_path or not csv_path.exists():
            continue
        try:
            df = pd.read_csv(csv_path)
        except Exception as exc:
            dbg("WARN", "run", f"pdb={target.pdb_id} auto_select_read_err={exc}")
            continue

        if "run_id" not in df.columns:
            continue

        try:
            lig_col = guess_ligfile_col(df, lig_col_cli)
        except Exception as exc:
            dbg("WARN", "run", f"pdb={target.pdb_id} auto_select_ligcol_err={exc}")
            continue

        lig_series = df[lig_col].astype(str)
        parsed = lig_series.apply(parse_name_and_label)
        lig_ids = parsed.apply(lambda t: t[0])

        run_series = df["run_id"].dropna()
        if run_series.empty:
            continue
        run_ids = run_series.index
        df_subset = df.loc[run_ids].copy()
        df_subset["__run_id"] = run_series.astype(str).str.strip()
        df_subset["__lig_id"] = lig_ids.loc[run_ids]
        df_subset = df_subset[df_subset["__run_id"] != ""]
        if df_subset.empty:
            continue

        lig_counts = df_subset.groupby("__run_id")["__lig_id"].nunique()
        for run_label, count in lig_counts.items():
            if count >= FULL_RUN_MIN_LIGANDS:
                run_full_map.setdefault(run_label, set()).add(target.target_key)
            try:
                stat = csv_path.stat()
            except Exception:
                continue
            run_mtimes[run_label] = max(
                run_mtimes.get(run_label, 0.0), getattr(stat, "st_mtime", 0.0)
            )

    if not run_full_map:
        return None

    scored = {run_id: len(pdbs) for run_id, pdbs in run_full_map.items()}
    best_score = max(scored.values())
    candidates = [run_id for run_id, score in scored.items() if score == best_score]
    if not candidates:
        return None

    parsed_values: Dict[str, float] = {}
    all_parseable = True
    for run_id in candidates:
        token = run_id.strip()
        parsed_val: Optional[float] = None
        if digit_re.match(token):
            try:
                parsed_val = float(int(token))
            except Exception:
                parsed_val = None
        elif iso_re.match(token):
            parsed_val = _parse_isoish(token)
        if parsed_val is None:
            all_parseable = False
            break
        parsed_values[run_id] = parsed_val

    if all_parseable and parsed_values:
        return max(parsed_values, key=lambda k: (parsed_values[k], k))

    return max(candidates, key=lambda k: (run_mtimes.get(k, 0.0), k))


def _load_default_cfg() -> Dict:
    cfg_path = Path("config.txt")
    try:
        from config.runtime_config import load_config, validate_config
    except Exception:
        dbg(
            "WARN",
            "config",
            "input_and_export_functions unavailable; skipping config.txt",
        )
        return {}

    try:
        if not cfg_path.exists():
            dbg("DEBUG", "config", "config.txt not found; using CLI defaults")
            return {}
        cfg = load_config() or {}
        if cfg:
            validate_config(cfg)
            dbg("DEBUG", "config", f"config.txt loaded keys={sorted(cfg.keys())}")
        else:
            dbg("DEBUG", "config", "config.txt empty; using CLI defaults")
        return cfg
    except Exception as exc:
        dbg("WARN", "config", f"failed to load config.txt err={exc}")
        return {}


def _resolve_run_label(
    pdb_id: str,
    meta: Optional[TargetEvaluation],
    cli_run_id: Optional[str],
    active_run_id: Optional[str],
) -> str:
    # Primary: mirror the metrics row's run_id, if it exists.
    if meta is not None and meta.metrics is not None:
        if "run_id" in meta.metrics.index:
            return str(meta.metrics["run_id"])

    # Fallback: previous behaviour, but never join multiple IDs.
    if cli_run_id:
        if meta and meta.has_run_id_column:
            return str(cli_run_id)
        return "(none)"

    if meta and meta.has_run_id_column:
        if meta.run_ids:
            ordered = sorted(str(x) for x in meta.run_ids if x)
            if not ordered:
                return ""
            if active_run_id and active_run_id in ordered:
                return active_run_id
            return ordered[0]
        return ""

    if active_run_id:
        return str(active_run_id)
    return "(none)"


def _candidate_protein_logs(
    pdb_id: str,
    csv_path: Path,
    docked_root: Path,
    log_root_override: Optional[Path],
    cfg: Dict,
) -> Tuple[Optional[Path], List[Path]]:
    candidates: List[Path] = []
    seen: Set[str] = set()

    def _add(path: Path) -> None:
        try:
            key = str(path.resolve())
        except Exception:
            key = str(path)
        if key not in seen:
            seen.add(key)
            candidates.append(path)

    if log_root_override:
        _add(log_root_override / pdb_id / "protein.log")
    else:
        # Prefer the provided docked_root, fall back to legacy docked/ root last.
        _add(docked_root / pdb_id / "protein.log")
        _add(Path("docked") / pdb_id / "protein.log")

    if docked_root.is_dir():
        if docked_root.name.upper() == pdb_id.upper():
            _add(docked_root / "protein.log")
        _add(docked_root / pdb_id / "protein.log")
    else:
        _add(docked_root / pdb_id / "protein.log")

    csv_parent = csv_path.parent
    _add(csv_parent / "protein.log")
    if csv_parent.name.upper() != pdb_id.upper():
        _add(csv_parent.parent / "protein.log")

    if cfg:
        try:
            paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
        except Exception as exc:
            dbg("WARN", "control", f"pdb={pdb_id} router_err={exc}")
        else:
            _add(paths.docked_root / paths.pdb_id / "protein.log")

    selected: Optional[Path] = None
    for cand in candidates:
        if cand.exists():
            selected = cand
            break
    return selected, candidates


def _manifest_library_for_target(
    manifest: dict,
    pdb_id: str,
    variant_label: Optional[str],
    ph_tag: Optional[str],
) -> str:
    """
    Pull the library label for a given target directly from the run manifest.
    Falls back to a pdb_id-only search if the exact variant/pH key is missing.
    """
    if not manifest:
        return ""

    proteins = manifest.get("proteins") or {}
    if not isinstance(proteins, dict):
        return ""

    pdb_norm = str(pdb_id).upper()
    variant_norm = (variant_label or "legacy").strip().upper() or "LEGACY"
    ph_norm = (ph_tag or "base").strip()
    ph_token = ph_norm if ph_norm else "base"
    key = f"{pdb_norm}|{variant_norm}|{ph_token}"

    entry = proteins.get(key)
    if isinstance(entry, dict):
        lib = entry.get("library")
        if isinstance(lib, str) and lib.strip():
            return lib.strip()

    for _k, e in proteins.items():
        if not isinstance(e, dict):
            continue
        if str(e.get("pdb_id", "")).upper() != pdb_norm:
            continue
        lib = e.get("library")
        if isinstance(lib, str) and lib.strip():
            return lib.strip()

    return ""


def _manifest_protein_entries(manifest: dict) -> List[Dict[str, str]]:
    """
    Return the manifest proteins entries in order, without collapsing by pdb_id.
    Each entry is a dict with pdb_id, variant, and ph (strings, may be empty).
    """
    entries: List[Dict[str, str]] = []
    if not manifest:
        return entries
    proteins = manifest.get("proteins") or {}
    if not isinstance(proteins, dict):
        return entries
    for _key, entry in proteins.items():
        if not isinstance(entry, dict):
            continue
        pdb_id = str(entry.get("pdb_id", "")).strip()
        if not pdb_id:
            continue
        variant = str(entry.get("variant", "") or "").strip()
        ph = str(entry.get("ph", "") or "").strip()
        entries.append({"pdb_id": pdb_id, "variant": variant, "ph": ph})
    ph_scoped_keys = {
        (entry["pdb_id"].upper(), entry["variant"].upper())
        for entry in entries
        if entry["ph"]
    }
    filtered: List[Dict[str, str]] = []
    for entry in entries:
        key = (entry["pdb_id"].upper(), entry["variant"].upper())
        if not entry["ph"] and key in ph_scoped_keys:
            continue
        filtered.append(entry)
    return filtered


def _manifest_proteins_by_pdb(manifest: dict) -> Dict[str, Dict[str, str]]:
    """
    Return a mapping:
        pdb_id -> {"variant": <variant or \"\">, "ph": <ph or \"\">}
    using the manifest's proteins section.

    If multiple manifest entries share the same pdb_id, prefer the first.
    We treat the manifest's 'variant' and 'ph' fields as canonical labels
    and do NOT try to parse or split composite pH tags like 'pH7_2+7_7-dup4'.
    """
    out: Dict[str, Dict[str, str]] = {}
    if not manifest:
        return out
    proteins = manifest.get("proteins") or {}
    if not isinstance(proteins, dict):
        return out
    for _key, entry in proteins.items():
        if not isinstance(entry, dict):
            continue
        pdb_id = str(entry.get("pdb_id", "")).strip()
        if not pdb_id:
            continue
        if pdb_id in out:
            continue
        variant = str(entry.get("variant", "") or "").strip()
        ph = str(entry.get("ph", "") or "").strip()
        out[pdb_id] = {"variant": variant, "ph": ph}
    return out


def _build_control_records(
    pdb_id: str,
    run_id_label: str,
    target_name: str,
    library_name: str,
    log_path: Optional[Path],
    variant: Optional[str] = None,
    ph_tag: Optional[str] = None,
) -> Tuple[List[Dict], Dict, Dict[str, int]]:
    long_rows: List[Dict] = []
    controls_found = 0
    max_spread: Optional[float] = None
    policy = ""
    centers_detected = False

    if not log_path or not log_path.exists():
        summary_row = {
            "variant": variant,
            "pH": ph_tag,
            "pdb_id": pdb_id,
            "run_id": run_id_label,
            "target_name": target_name,
            "library_name": library_name,
            "controls_found": 0,
            "control_id": "",
            "rmsd_to_crystal_A": None,
            "redock_best_energy_kcal_mol": None,
            "chosen": 0,
            "max_spread_A": None,
            "policy": "",
        }
        meta = {"centers_found": 0, "redock_lines": 0, "report_rows": 0}
        return long_rows, summary_row, meta

    try:
        with log_path.open("r", encoding="utf-8", errors="ignore") as fh:
            for raw_line in fh:
                line = raw_line.strip()
                if not line:
                    continue
                m_centers = CONTROL_CENTERS_RE.search(line)
                if m_centers:
                    centers_detected = True
                    try:
                        controls_found = int(m_centers.group(1))
                    except Exception:
                        controls_found = 0
                    try:
                        max_spread = float(m_centers.group(2))
                    except Exception:
                        max_spread = None
                    policy = m_centers.group(3)
                    continue
                m_control = CONTROL_REDOCK_RE.search(line)
                if not m_control:
                    continue
                control_id = m_control.group(1)
                try:
                    rmsd = float(m_control.group(2))
                except Exception:
                    rmsd = float("nan")
                try:
                    score = float(m_control.group(3))
                except Exception:
                    score = float("nan")
                long_rows.append(
                    {
                        "variant": variant,
                        "pH": ph_tag,
                        "pdb_id": pdb_id,
                        "run_id": run_id_label,
                        "target_name": target_name,
                        "library_name": library_name,
                        "controls_found": 0,  # fill later
                        "control_id": control_id,
                        "rmsd_to_crystal_A": rmsd,
                        "redock_best_energy_kcal_mol": score,
                        "chosen": 0,
                        "max_spread_A": None,
                        "policy": "",
                    }
                )
    except Exception as exc:
        dbg("WARN", "control", f"pdb={pdb_id} log_read_err={exc}")
        long_rows.clear()
        controls_found = 0
        max_spread = None
        policy = ""

    if controls_found == 0 and long_rows:
        controls_found = len(long_rows)

    chosen_index: Optional[int] = None
    if long_rows:
        # Fill common fields first
        for row in long_rows:
            row["controls_found"] = controls_found
            row["max_spread_A"] = max_spread
            row["policy"] = policy

        def _rmsd_key(val: Optional[float]) -> float:
            if val is None:
                return float("inf")
            try:
                return float(val) if not math.isnan(float(val)) else float("inf")
            except Exception:
                return float("inf")

        ordered = sorted(
            enumerate(long_rows),
            key=lambda pair: (
                _rmsd_key(pair[1]["rmsd_to_crystal_A"]),
                str(pair[1]["control_id"]),
            ),
        )
        if ordered:
            chosen_index = ordered[0][0]
            long_rows[chosen_index]["chosen"] = 1
            chosen_row = long_rows[chosen_index]
            dbg(
                "DEBUG",
                "control",
                (
                    f"pdb={pdb_id} chosen_control={chosen_row['control_id']} "
                    f"best_rmsd_A={chosen_row['rmsd_to_crystal_A']} "
                    f"best_score_kcal={chosen_row['redock_best_energy_kcal_mol']}"
                ),
            )

    summary_row = {
        "variant": variant,
        "pH": ph_tag,
        "pdb_id": pdb_id,
        "run_id": run_id_label,
        "target_name": target_name,
        "library_name": library_name,
        "controls_found": controls_found,
        "control_id": "",
        "rmsd_to_crystal_A": None,
        "redock_best_energy_kcal_mol": None,
        "chosen": 0,
        "max_spread_A": cast(Any, max_spread),
        "policy": policy,
    }

    if chosen_index is not None:
        chosen = long_rows[chosen_index]
        summary_row.update(
            {
                "control_id": chosen["control_id"],
                "rmsd_to_crystal_A": chosen["rmsd_to_crystal_A"],
                "redock_best_energy_kcal_mol": chosen["redock_best_energy_kcal_mol"],
                "chosen": 1,
            }
        )

    meta = {
        "centers_found": 1 if centers_detected else 0,
        "redock_lines": len(long_rows),
        "report_rows": len(long_rows),
    }

    return long_rows, summary_row, meta
