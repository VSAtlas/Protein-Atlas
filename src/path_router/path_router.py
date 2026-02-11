# -*- coding: utf-8 -*-
# path_router.py
from __future__ import annotations

import json
import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional


# ---------------------------
# Helpers
# ---------------------------
def _canon(val: Optional[str]) -> str:
    """
    Canonicalize a string for loose comparisons:
    - lowercased
    - remove spaces, dashes, underscores
    - keep only alphanumerics
    """
    if val is None:
        return ""
    s = str(val).lower()
    s = s.replace(" ", "").replace("-", "").replace("_", "")
    return "".join(ch for ch in s if ch.isalnum())


def expand_variants(mode: Optional[str]) -> list[Optional[str]]:
    """
    Expand a user/config 'mode' into concrete variants to iterate.

    Returns:
      - [None]          for legacy/no-variant (mode None/"", or strings like "none"/"legacy"/"null")
      - ["APO"]         for "apo"
      - ["HOLO"]        for "holo"
      - ["APO","HOLO"]  for "apo_vs_holo" (any loose spelling like "apo vs holo", "apo-vs-holo", "both")

    Notes:
      * Single-variant builders (receptor_dir, receptor_cleaned_pdb, receptor_pdbqt, docked_variant_root,
        docked_stage_dir, configs_stage_dir) should be called with concrete variants only (APO/HOLO/None).
        If you have a multi-variant mode, call this function and iterate.
    """
    key = _canon(mode)
    expanded: list[Optional[str]]
    if key in {"", "none", "legacy", "null"}:
        expanded = [None]
    elif key == "apo":
        expanded = ["APO"]
    elif key == "holo":
        expanded = ["HOLO"]
    elif key in {"apovsholo", "both", "apoandholo"}:
        expanded = ["APO", "HOLO"]
    else:
        # Fallback: treat odd casing like "Apo" or "Holo"
        v = _norm_variant(mode)
        expanded = [v] if v is not None else [None]

    logger = logging.getLogger("path_router")
    logger.info("[router.debug] expand_variants.in=%r -> %r", mode, expanded)
    return expanded


def _norm_variant(variant: Optional[str]) -> Optional[str]:
    """
    Normalize a *single* variant:
      - returns "APO" | "HOLO" | None
      - strings like "none"/"legacy"/"null" are treated as None
      - raises on multi-variant modes (e.g., "apo_vs_holo") to avoid ambiguity in single-variant builders
    """
    if variant is None:
        return None
    key = _canon(variant)
    if key in {"", "none", "legacy", "null"}:
        return None
    if key == "apo":
        return "APO"
    if key == "holo":
        return "HOLO"
    if key in {"apovsholo", "both", "apoandholo"}:
        raise ValueError(
            "variant denotes a multi-variant mode; expand via expand_variants(mode) and iterate APO/HOLO."
        )
    raise ValueError(f"variant must be APO|HOLO or None, got {variant!r}")


# ---------------------------
# pH  Helpers
# ---------------------------


def _ph_manifest_candidates(pdb_id: str, variant: Optional[str]) -> list[Path]:
    """
    Return the only two allowed search locations for ensemble.json, in order:
      1) processed_pdbs/<PDB>/<VARIANT>/receptor/ph_ensemble/ensemble.json (when variant given)
      2) processed_pdbs/<PDB>/receptor/ph_ensemble/ensemble.json
    """
    roots = _ensure_router_roots()
    token = _norm_pdb_id(pdb_id)
    v = _norm_variant(variant)  # returns "APO"|"HOLO"|None

    base = roots.processed / token
    out: list[Path] = []
    if v:
        out.append(base / v / "receptor" / "ph_ensemble" / "ensemble.json")
    out.append(base / "receptor" / "ph_ensemble" / "ensemble.json")
    return out


# ---------------------------
# pH ensemble public helper
# ---------------------------


def ph_ensemble_dir(
    pdb_id: str,
    variant: Optional[str] = None,
    legacy: bool = False,
) -> Path:
    """
    Return the canonical directory for pH ensemble artifacts.

    If variant in {APO, HOLO}: processed_pdbs/<PDB>/<VARIANT>/receptor/ph_ensemble/
    Else (legacy):             processed_pdbs/<PDB>/receptor/ph_ensemble/
    """
    roots = _ensure_router_roots()
    token = _norm_pdb_id(pdb_id)
    base = roots.processed / token
    v = _norm_variant(variant)
    logger = logging.getLogger("path_router")
    logger.info(
        "[router.debug] ph_ensemble_dir variant_in=%r norm=%r legacy_flag=%s",
        variant,
        v,
        legacy,
    )
    if v and not legacy:
        base = base / v
    return base / "receptor" / "ph_ensemble"


# ---------------------------
# Router roots & stateless helpers
# ---------------------------


@dataclass(frozen=True)
class RouterRoots:
    overall: Path
    processed: Path
    docked: Path
    configs: Path


_ROUTER_ROOTS: Optional[RouterRoots] = None


def _norm_pdb_id(pdb_id: str) -> str:
    return str(pdb_id).strip().upper()


def _default_config_path() -> Path:
    cfg_env = os.environ.get("ATLAS_CONFIG")
    if cfg_env:
        return Path(cfg_env).expanduser()
    # path_router.py is in src/path_router/, so repo root is two levels up
    return Path(__file__).resolve().parents[2] / "config.txt"


def _load_router_roots() -> RouterRoots:
    cfg_path = _default_config_path()
    expanded: Dict[str, str]
    try:
        from input_and_export_functions import load_config

        expanded = load_config(config_path=str(cfg_path), base_dir=cfg_path.parent)
    except Exception:
        # Keep router functional even if shared config import fails.
        expanded = {}
    for key in ("OVERALL_DIR", "OUTPUT_DIR", "DOCKED_DIR"):
        env_val = os.environ.get(key)
        if env_val:
            expanded[key] = env_val
    run_id = (os.environ.get("ATLAS_RUN_ID") or "").strip()
    cfg_env = os.environ.get("CONFIGS_DIR")
    if cfg_env:
        expanded["CONFIGS_DIR"] = cfg_env

    over_raw = expanded.get("OVERALL_DIR")
    overall = (
        Path(over_raw).expanduser()
        if over_raw
        else Path(__file__).resolve().parents[2]
    )

    out_raw = expanded.get("OUTPUT_DIR")
    processed = Path(out_raw).expanduser() if out_raw else overall / "processed_pdbs"

    dock_raw = expanded.get("DOCKED_DIR")
    docked = Path(dock_raw).expanduser() if dock_raw else overall / "docked"
    if run_id:
        try:
            if docked.name != run_id:
                docked = docked / run_id
        except Exception:
            docked = docked / run_id

    cfg_raw = expanded.get("CONFIGS_DIR")
    configs = Path(cfg_raw).expanduser() if cfg_raw else overall / "configs"

    if os.environ.get("ROUTER_DEBUG"):
        print(
            "[router.debug.env]"
            f" OVERALL_DIR={over_raw!r}"
            f" OUTPUT_DIR={out_raw!r}"
            f" DOCKED_DIR={dock_raw!r}"
            f" CONFIGS_DIR={cfg_raw!r}"
        )

    return RouterRoots(
        overall=overall, processed=processed, docked=docked, configs=configs
    )


def _set_router_roots(
    overall: Path, processed: Path, docked: Path, configs: Optional[Path] = None
) -> None:
    global _ROUTER_ROOTS
    _ROUTER_ROOTS = RouterRoots(
        overall=Path(overall).expanduser(),
        processed=Path(processed).expanduser(),
        docked=Path(docked).expanduser(),
        configs=Path(configs).expanduser()
        if configs
        else Path(overall).expanduser() / "configs",
    )


def _ensure_router_roots() -> RouterRoots:
    global _ROUTER_ROOTS
    if _ROUTER_ROOTS is not None:
        return _ROUTER_ROOTS

    try:
        _ROUTER_ROOTS = _load_router_roots()
    except Exception:
        cwd = Path.cwd().resolve()
        _ROUTER_ROOTS = RouterRoots(
            overall=cwd,
            processed=cwd / "processed_pdbs",
            docked=cwd / "docked",
            configs=cwd / "configs",
        )

    return _ROUTER_ROOTS


# ---------------------------
# Public stateless path helpers
# ---------------------------


def run_logs_dir(cfg: Optional[Dict[str, object]] = None) -> Path:
    roots = _ensure_router_roots()
    docked_root = roots.docked
    run_id = (os.environ.get("ATLAS_RUN_ID") or "").strip() or None
    if cfg:
        try:
            cfg_run_id = str(cfg.get("RUN_ID") or "").strip()
            if cfg_run_id:
                run_id = cfg_run_id
        except Exception:
            pass
        try:
            cfg_docked = cfg.get("DOCKED_DIR")
            if cfg_docked:
                docked_root = Path(str(cfg_docked)).expanduser()
        except Exception:
            pass
    if run_id and docked_root.name != run_id:
        docked_root = docked_root / run_id
    log_dir = docked_root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("path_router")
    logger.info("[router.debug] kind=run_logs_dir path=%s", log_dir)
    return log_dir


def receptor_dir(
    pdb_id: str,
    variant: Optional[str] = None,
    ph_tag: Optional[str] = None,
    legacy: bool = False,
) -> Path:
    roots = _ensure_router_roots()
    token = _norm_pdb_id(pdb_id)
    base = roots.processed / token
    v = _norm_variant(variant)
    logger = logging.getLogger("path_router")
    logger.info(
        "[router.debug] variant_in=%r norm=%r legacy_flag=%s", variant, v, legacy
    )
    if v:
        base = base / v
    dir_path = base / "receptor"
    if ph_tag:
        dir_path = dir_path / "ph_ensemble"
    return dir_path


def receptor_file(
    pdb_id: str,
    variant: Optional[str] = None,
    ph_tag: Optional[str] = None,
    legacy: bool = False,
) -> Path:
    dir_path = receptor_dir(pdb_id, variant=variant, ph_tag=ph_tag, legacy=legacy)
    stem = _norm_pdb_id(pdb_id)
    suffix = f"_{ph_tag}" if ph_tag else ""
    return dir_path / f"{stem}{suffix}.pdbqt"


def variant_root(pdb_id: str, variant: Optional[str] = None) -> Path:
    """Return the processed/ tree root for a given PDB/variant."""
    roots = _ensure_router_roots()
    token = _norm_pdb_id(pdb_id)
    base = roots.processed / token
    v = _norm_variant(variant)
    if v:
        base = base / v
    logger = logging.getLogger("path_router")
    logger.info(
        "[router.debug] kind=variant_root variant=%s path=%s", v or "None", base
    )
    return base


def raw_dir_v(pdb_id: str, variant: Optional[str] = None) -> Path:
    """Variant-aware raw/ directory."""
    base = variant_root(pdb_id, variant=variant)
    v = _norm_variant(variant)
    path = base / "raw"
    logger = logging.getLogger("path_router")
    logger.info("[router.debug] kind=raw variant=%s path=%s", v or "None", path)
    return path


def nolig_dir_v(pdb_id: str, variant: Optional[str] = None) -> Path:
    """Variant-aware nolig/ directory."""
    base = variant_root(pdb_id, variant=variant)
    v = _norm_variant(variant)
    path = base / "nolig"
    logger = logging.getLogger("path_router")
    logger.info("[router.debug] kind=nolig variant=%s path=%s", v or "None", path)
    return path


def work_dir_v(pdb_id: str, variant: Optional[str] = None) -> Path:
    """Variant-aware work/ directory."""
    base = variant_root(pdb_id, variant=variant)
    v = _norm_variant(variant)
    path = base / "work"
    logger = logging.getLogger("path_router")
    logger.info("[router.debug] kind=work variant=%s path=%s", v or "None", path)
    return path


_raw_dir_v_fn = raw_dir_v
_nolig_dir_v_fn = nolig_dir_v
_work_dir_v_fn = work_dir_v


def docked_dir(
    pdb_id: str,
    variant: Optional[str] = None,
    ph_tag: Optional[str] = None,
    legacy: bool = False,
) -> Path:
    """
    Run-scoped docked path:
      - with ATLAS_RUN_ID or run-scoped roots: DOCKED_DIR/<RUN_ID>/<PDB>/[variant]/[ph_tag]
      - without run id: legacy DOCKED_DIR/<PDB>/[variant]/[ph_tag]
    """
    roots = _ensure_router_roots()
    token = _norm_pdb_id(pdb_id)
    base = roots.docked / token
    logger = logging.getLogger("path_router")
    if legacy and not ph_tag:
        logger.info(
            "[router.debug] variant_in=%r norm=%r legacy_flag=%s", variant, None, legacy
        )
        return base
    v = _norm_variant(variant)
    logger.info(
        "[router.debug] variant_in=%r norm=%r legacy_flag=%s", variant, v, legacy
    )
    if v:
        base = base / v
    if ph_tag:
        base = base / str(ph_tag)
    return base


# ---------------------------
# Config helpers (stateless)
# ---------------------------
def _norm_stage(stage: str) -> str:
    return str(stage).strip() or "stage"


def config_dir(
    run_id: str,
    pdb_id: str,
    stage: str,
    variant: Optional[str] = None,
    ph_tag: Optional[str] = None,
    legacy: bool = False,
) -> Path:
    roots = _ensure_router_roots()
    token = _norm_pdb_id(pdb_id)
    base = roots.configs / str(run_id) / token
    v = _norm_variant(variant)
    logger = logging.getLogger("path_router")
    logger.info(
        "[router.debug] variant_in=%r norm=%r legacy_flag=%s", variant, v, legacy
    )
    if v:
        base = base / v
    if ph_tag:
        base = base / str(ph_tag)
    return base / _norm_stage(stage)


def config_file(
    run_id: str,
    pdb_id: str,
    stage: str,
    variant: Optional[str] = None,
    ph_tag: Optional[str] = None,
    name: str = "vina.json",
    legacy: bool = False,
) -> Path:
    return (
        config_dir(
            run_id,
            pdb_id,
            stage,
            variant=variant,
            ph_tag=ph_tag,
            legacy=legacy,
        )
        / name
    )


# --- full replacement for load_ph_tags() ---
def load_ph_tags(pdb_id: str, variant: Optional[str] = None) -> list[str]:
    """
    Read ensemble.json and return ordered pH tags for the given protein/variant.

    Search order:
      1) processed_pdbs/<PDB>/<VARIANT>/receptor/ph_ensemble/ensemble.json  (if variant set)
      2) processed_pdbs/<PDB>/receptor/ph_ensemble/ensemble.json
    """
    logger = logging.getLogger("path_router")

    tags: list[str] = []
    seen: set[str] = set()
    prefix = f"{_norm_pdb_id(pdb_id)}_"

    for manifest in _ph_manifest_candidates(pdb_id, variant):
        try:
            if not manifest.exists():
                continue
            payload = json.loads(manifest.read_text(encoding="utf-8", errors="ignore"))
        except Exception:
            continue

        members = payload.get("members") if isinstance(payload, dict) else None
        if not isinstance(members, list):
            continue

        for entry in members:
            if not isinstance(entry, dict):
                continue
            raw = entry.get("label") or entry.get("ph_label")
            if raw is not None:
                label = str(raw)
            else:
                receptor_path = (
                    entry.get("pdbqt")
                    or entry.get("receptor_pdbqt")
                    or entry.get("output_pdbqt")
                    or entry.get("path")
                    or entry.get("receptor")
                )
                if not receptor_path:
                    continue
                stem = Path(str(receptor_path)).stem
                label = stem[len(prefix) :] if stem.startswith(prefix) else stem

            if label in seen:
                continue
            seen.add(label)
            tags.append(label)

        if tags:
            logger.info(
                "[router.ph] using manifest=%s tags=%s", manifest, ",".join(tags)
            )
            break

    return tags


def print_pathmap(
    *,
    pdb_id: str,
    variant: Optional[str] = None,
    ph_tag: Optional[str] = None,
    legacy: bool = False,
) -> None:
    logger = logging.getLogger("path_router")
    logger.info(
        "[router.debug] variant_in=%r norm=%r legacy_flag=%s",
        variant,
        _norm_variant(variant),
        legacy,
    )
    rec_dir = receptor_dir(pdb_id, variant=variant, ph_tag=ph_tag, legacy=legacy)
    print(f"receptor_dir={rec_dir}")
    rec_file = receptor_file(pdb_id, variant=variant, ph_tag=ph_tag, legacy=legacy)
    print(f"receptor_file={rec_file}")
    ph_dir = ph_ensemble_dir(pdb_id, variant=variant, legacy=legacy)
    print(f"ph_ensemble_dir={ph_dir}")
    dock_dir = docked_dir(pdb_id, variant=variant, ph_tag=ph_tag, legacy=legacy)
    print(f"docked_dir={dock_dir}")


# ---------------------------
# Core dataclass
# ---------------------------
@dataclass(frozen=True)
class Paths:
    """
    Centralized per-PDB path router.

    Back-compat notes:
      * If variant=None, receptor/config/docked paths use the legacy "no-variant" layout.
      * pH is expressed as an optional *filename token* (no water tokens by design).
    """

    # Identity
    pdb_id: str  # e.g., "1T46" (uppercased)
    pdb_file: str  # e.g., "1T46.pdb" (as provided)

    # Roots (from cfg)
    over_root: Path  # OVERALL_DIR
    input_root: Path  # INPUT_DIR  (input_pdbs/)
    processed_root: Path  # OUTPUT_DIR (processed_pdbs/)
    docked_root: Path  # DOCKED_DIR (docked/[RUN_ID]/)
    prepped_root: Path  # PREPPED_LIGANDS_DIR
    ligands_mol2_root: Path  # LIGANDS_MOL2_DIR
    configs_root: Path  # default OVERALL_DIR/configs (unless CONFIGS_DIR supplied)

    # Canonical per-PDB roots (variant-agnostic)
    root_pdb_dir: Path  # processed_pdbs/<PDB>/

    # Shared subfolders under processed_pdbs/<PDB>/
    raw_dir: Path  # processed_pdbs/<PDB>/raw/          # copies/ingest
    ligands_raw_dir: (
        Path  # processed_pdbs/<PDB>/ligands_raw/  # extracted crystal ligands (*.pdb)
    )
    nolig_dir: Path  # processed_pdbs/<PDB>/nolig/        # <PDB>_nolig.pdb intermediate
    work_dir: Path  # processed_pdbs/<PDB>/work/         # scratch & intermediates

    # Inputs
    input_pdb_path: Path  # input_pdbs/<pdb_file>              # source PDB to process

    # Shared expected filenames (variant-agnostic)
    nolig_pdb_path: Path  # processed_pdbs/<PDB>/nolig/<PDB>_nolig.pdb

    # ---------- Variant-aware receptor paths (back-compat with variant=None) ----------
    def receptor_dir(self, variant: Optional[str]) -> Path:
        """
        If variant in {APO,HOLO}:
          Dir: processed_pdbs/<PDB>/<VARIANT>/receptor/
          Files: <PDB>_cleaned.pdb, <PDB>.pdbqt, optional pH-tokenized PDBQTs.
        If variant is None (legacy):
          Dir: processed_pdbs/<PDB>/receptor/
          Files: <PDB>_cleaned.pdb, <PDB>.pdbqt, optional pH-tokenized PDBQTs.
        """
        p = receptor_dir(self.pdb_id, variant=variant)
        p.mkdir(parents=True, exist_ok=True)
        return p

    def receptor_cleaned_pdb(self, variant: Optional[str]) -> Path:
        """
        File:
          * With variant: processed_pdbs/<PDB>/<VARIANT>/receptor/<PDB>_cleaned.pdb
          * No variant : processed_pdbs/<PDB>/receptor/<PDB>_cleaned.pdb
        """
        return self.receptor_dir(variant) / f"{self.pdb_id}_cleaned.pdb"

    def receptor_pdbqt(
        self, variant: Optional[str], ph_token: Optional[str] = None
    ) -> Path:
        """
        File (no pH):            <receptor_dir>/<PDB>.pdbqt
        File (with pH token):    <receptor_dir>/ph_ensemble/<PDB>_<ph_token>.pdbqt
        Examples of ph_token you may pass:
          "pH8_0+8_1+8_2"   -> 6LYZ_pH8_0+8_1+8_2.pdbqt
          "pH8_0"           -> 6LYZ_pH8_0.pdbqt
        """
        path = receptor_file(self.pdb_id, variant=variant, ph_tag=ph_token)
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def raw_dir_v(self, variant: Optional[str]) -> Path:
        path = _raw_dir_v_fn(self.pdb_id, variant=variant)
        path.mkdir(parents=True, exist_ok=True)
        return path

    def nolig_dir_v(self, variant: Optional[str]) -> Path:
        path = _nolig_dir_v_fn(self.pdb_id, variant=variant)
        path.mkdir(parents=True, exist_ok=True)
        return path

    def work_dir_v(self, variant: Optional[str]) -> Path:
        path = _work_dir_v_fn(self.pdb_id, variant=variant)
        path.mkdir(parents=True, exist_ok=True)
        return path

    # ---------------------------
    # Ligand prep / libraries
    # ---------------------------
    @property
    def ligand_output_dir(self) -> Path:
        """
        Dir: processed_pdbs/<PDB>/ligands_raw/
        Expected: extracted crystal ligands (*.pdb), possibly '*.sanitized.pdb'
        """
        self.ligands_raw_dir.mkdir(parents=True, exist_ok=True)
        return self.ligands_raw_dir

    @property
    def ligands_mol2_dir(self) -> Path:
        """
        Dir: ligands_mol2/<PDB>/
        Expected: MOL2 intermediates (*.mol2) during ligand prep
        """
        d = self.ligands_mol2_root / self.pdb_id
        d.mkdir(parents=True, exist_ok=True)
        return d

    @property
    def prepped_ligands_dir(self) -> Path:
        """
        Dir: prepped_ligands/<PDB>/
        Expected: prepared ligand library (*.pdbqt), incl. controls; optional quarantine/
        """
        d = self.prepped_root / self.pdb_id
        d.mkdir(parents=True, exist_ok=True)
        return d

    # Optional legacy global extracted ligands root (if present in cfg).
    # Historically some code used EXTRACTED_LIGANDS_DIR/<PDB>/.
    def legacy_extracted_ligands_dir(self) -> Optional[Path]:
        """
        Legacy (optional): <EXTRACTED_LIGANDS_DIR>/<PDB>/
        Note: canonical location is processed_pdbs/<PDB>/ligands_raw/.
        """
        return None  # kept as a stub to avoid reintroducing legacy in new code

    # ---------------------------
    # Docking outputs
    # ---------------------------
    def docked_pdb_root(self) -> Path:
        """
        Dir: docked/<RUN_ID>/<PDB>/ (or docked/<PDB>/ if no run_id)
        Expected (top-level analysis):
          - docking_score_long.csv
          - docking_score_summary.csv
          - bench_pocketX_single/ (optional benchmark folders; X is an index)
          - <VARIANT>/          (subdirs when variant-aware outputs are used)
        """
        d = self.docked_root / self.pdb_id
        d.mkdir(parents=True, exist_ok=True)
        return d

    def docking_score_long_csv(self) -> Path:
        """File: docked/<RUN_ID>/<PDB>/docking_score_long.csv"""
        return self.docked_pdb_root() / "docking_score_long.csv"

    def docking_score_summary_csv(self) -> Path:
        """File: docked/<RUN_ID>/<PDB>/docking_score_summary.csv"""
        return self.docked_pdb_root() / "docking_score_summary.csv"

    def bench_pocket_dir(self, pocket_index: int) -> Path:
        """
        Dir: docked/<RUN_ID>/<PDB>/bench_pocketX_single/
        Example: docked/<RUN_ID>/1T46/bench_pocket1_single/
        """
        d = self.docked_pdb_root() / f"bench_pocket{int(pocket_index)}_single"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def docked_variant_root(
        self, variant: Optional[str], ph_label: Optional[str] = None
    ) -> Path:
        """
        If variant in {APO,HOLO}:
          Dir: docked/<RUN_ID>/<PDB>/<VARIANT>/
          Expected: stage1/, stage2/, stage3/ ... (poses, logs, scores)
        If variant is None (legacy):
          Dir: docked/<RUN_ID>/<PDB>/         (stages under top-level; back-compat)
        """
        base = docked_dir(self.pdb_id, variant=variant, ph_tag=ph_label)
        base.mkdir(parents=True, exist_ok=True)
        return base

    def docked_stage_dir(
        self,
        variant: Optional[str],
        stage: str,
        ph_label: Optional[str] = None,
    ) -> Path:
        """
        If variant in {APO,HOLO}:
          Dir: docked/<RUN_ID>/<PDB>/<VARIANT>/<stage>/
        If variant is None:
          Dir: docked/<RUN_ID>/<PDB>/<stage>/
        Expected files: <ligand>_<stage>.pdbqt (poses), vina logs, per-stage CSVs
        """
        d = self.docked_variant_root(variant, ph_label) / stage
        d.mkdir(parents=True, exist_ok=True)
        return d

    # ---------------------------
    # Config emission
    # ---------------------------
    def configs_stage_dir(
        self,
        run_id: str,
        variant: Optional[str],
        stage: str,
        ph_label: Optional[str] = None,
    ) -> Path:
        """
        If variant in {APO,HOLO}:
          Dir: configs/<RUN>/<PDB>/<VARIANT>/<stage>/
        If variant is None:
          Dir: configs/<RUN>/<PDB>/<stage>/
        Expected: <ligand>_<stage>.txt  (Vina config files per ligand)
        """
        v = _norm_variant(variant)
        base = self.configs_root / run_id / self.pdb_id
        d = (base / v) if v else base
        if ph_label:
            d = d / str(ph_label)
        d = d / stage
        d.mkdir(parents=True, exist_ok=True)
        return d

    def expand_variants(self, mode: Optional[str]) -> list[Optional[str]]:
        # Delegate to the module-level function to keep a single source of truth
        return expand_variants(mode)


# Back-compat alias: some modules still expect RouterPaths
RouterPaths = Paths


# ---------------------------
# Factory (required signature)
# ---------------------------
def make_paths(cfg: Dict, base_id: str, pdb_file: str) -> Paths:
    """
    Build per-PDB router from config (no global state).

    Ensures canonical per-PDB directories exist immediately:
      processed_pdbs/<PDB>/{raw, ligands_raw, nolig, work}/
      (variant receptor dirs are created lazily upon first access)
    """
    pdb_id = str(base_id).upper()
    # Normalize here too (idempotent if already clean)
    raw_id = str(base_id)
    pdb_id = re.sub(r"(?i)(_nolig(_cleaned)?|_cleaned)$", "", raw_id).upper()

    # Required roots
    over_root = Path(cfg["OVERALL_DIR"])
    input_root = Path(cfg["INPUT_DIR"])
    processed_root = Path(cfg["OUTPUT_DIR"])  # processed_pdbs/
    docked_root = Path(cfg["DOCKED_DIR"])  # docked/
    run_id = str(cfg.get("RUN_ID") or "").strip()
    if run_id:
        os.environ["ATLAS_RUN_ID"] = run_id
        if docked_root.name != run_id:
            docked_root = docked_root / run_id

    prepped_root = Path(cfg.get("PREPPED_LIGANDS_DIR", over_root / "prepped_ligands"))

    ligands_mol2_root = Path(cfg.get("LIGANDS_MOL2_DIR", over_root / "ligands_mol2"))

    # Configs root: default under OVERALL_DIR/configs (no CONFIGS_DIR in your cfg)
    configs_root = Path(cfg.get("CONFIGS_DIR", over_root / "configs"))

    # Canonical per-PDB root and shared subfolders
    root_pdb_dir = processed_root / pdb_id  # processed_pdbs/<PDB>/
    raw_dir = root_pdb_dir / "raw"  # processed_pdbs/<PDB>/raw/
    lig_raw_dir = root_pdb_dir / "ligands_raw"  # processed_pdbs/<PDB>/ligands_raw/
    nolig_dir = root_pdb_dir / "nolig"  # processed_pdbs/<PDB>/nolig/
    work_dir = root_pdb_dir / "work"  # processed_pdbs/<PDB>/work/

    # Ensure shared subfolders exist
    for d in (root_pdb_dir, raw_dir, lig_raw_dir, nolig_dir, work_dir):
        d.mkdir(parents=True, exist_ok=True)

    # Inputs / shared outputs
    input_pdb_path = input_root / pdb_file  # input_pdbs/<pdb_file>
    nolig_pdb_path = (
        nolig_dir / f"{pdb_id}_nolig.pdb"
    )  # processed_pdbs/<PDB>/nolig/<PDB>_nolig.pdb

    _set_router_roots(over_root, processed_root, docked_root, configs_root)

    return Paths(
        pdb_id=pdb_id,
        pdb_file=pdb_file,
        over_root=over_root,
        input_root=input_root,
        processed_root=processed_root,
        docked_root=docked_root,
        prepped_root=prepped_root,
        ligands_mol2_root=ligands_mol2_root,
        configs_root=configs_root,
        root_pdb_dir=root_pdb_dir,
        raw_dir=raw_dir,
        ligands_raw_dir=lig_raw_dir,
        nolig_dir=nolig_dir,
        work_dir=work_dir,
        input_pdb_path=input_pdb_path,
        nolig_pdb_path=nolig_pdb_path,
    )


if __name__ == "__main__":
    scenarios = [
        (
            "pH ON, apo/HOLO OFF",
            dict(pdb_id="3CS9", variant=None, ph_tag="pH6_0", legacy=False),
        ),
        (
            "apo/HOLO ON, pH ON",
            dict(pdb_id="3CS9", variant="APO", ph_tag="pH6_0+8_0", legacy=False),
        ),
        (
            "apo/HOLO ON, pH OFF",
            dict(pdb_id="3CS9", variant="HOLO", ph_tag=None, legacy=False),
        ),
        ("legacy", dict(pdb_id="3CS9", variant=None, ph_tag=None, legacy=True)),
    ]
    for label, params in scenarios:
        print(f"scenario={label}")
        print_pathmap(**params)  # type: ignore[arg-type]
        print()
