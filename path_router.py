# -*- coding: utf-8 -*-
# path_router.py
from __future__ import annotations

import os, re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Union

# ---------------------------
# Helpers
# ---------------------------
from typing import Iterable  # already imported Optional above; Iterable is harmless if unused elsewhere

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
    if key in {"", "none", "legacy", "null"}:
        return [None]
    if key == "apo":
        return ["APO"]
    if key == "holo":
        return ["HOLO"]
    if key in {"apovsholo", "both", "apoandholo"}:
        return ["APO", "HOLO"]

    # Fallback: treat odd casing like "Apo" or "Holo"
    v = _norm_variant(mode)
    return [v] if v is not None else [None]


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
    pdb_id: str                 # e.g., "1T46" (uppercased)
    pdb_file: str               # e.g., "1T46.pdb" (as provided)

    # Roots (from cfg)
    over_root: Path             # OVERALL_DIR
    input_root: Path            # INPUT_DIR  (input_pdbs/)
    processed_root: Path        # OUTPUT_DIR (processed_pdbs/)
    docked_root: Path           # DOCKED_DIR (docked/)
    prepped_root: Path          # PREPPED_LIGANDS_DIR / OUTPUT_LIGANDS_DIR / PREPPED_LIGANDS_ROOT
    ligands_mol2_root: Path     # LIGANDS_MOL2_DIR
    configs_root: Path          # default OVERALL_DIR/configs (unless CONFIGS_DIR supplied)

    # Canonical per-PDB roots (variant-agnostic)
    root_pdb_dir: Path          # processed_pdbs/<PDB>/

    # Shared subfolders under processed_pdbs/<PDB>/
    raw_dir: Path               # processed_pdbs/<PDB>/raw/          # copies/ingest
    ligands_raw_dir: Path       # processed_pdbs/<PDB>/ligands_raw/  # extracted crystal ligands (*.pdb)
    nolig_dir: Path             # processed_pdbs/<PDB>/nolig/        # <PDB>_nolig.pdb intermediate
    work_dir: Path              # processed_pdbs/<PDB>/work/         # scratch & intermediates

    # Inputs
    input_pdb_path: Path        # input_pdbs/<pdb_file>              # source PDB to process

    # Shared expected filenames (variant-agnostic)
    nolig_pdb_path: Path        # processed_pdbs/<PDB>/nolig/<PDB>_nolig.pdb

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
        v = _norm_variant(variant)
        p = (self.root_pdb_dir / v / "receptor") if v else (self.root_pdb_dir / "receptor")
        p.mkdir(parents=True, exist_ok=True)
        return p

    def receptor_cleaned_pdb(self, variant: Optional[str]) -> Path:
        """
        File:
          * With variant: processed_pdbs/<PDB>/<VARIANT>/receptor/<PDB>_cleaned.pdb
          * No variant : processed_pdbs/<PDB>/receptor/<PDB>_cleaned.pdb
        """
        return self.receptor_dir(variant) / f"{self.pdb_id}_cleaned.pdb"

    def receptor_pdbqt(self, variant: Optional[str], ph_token: Optional[str] = None) -> Path:
        """
        File (no pH):            <receptor_dir>/<PDB>.pdbqt
        File (with pH token):    <receptor_dir>/<PDB>_<ph_token>.pdbqt
        Examples of ph_token you may pass:
          "pH8_0+8_1+8_2"   -> 6LYZ_pH8_0+8_1+8_2.pdbqt
          "pH8_0"           -> 6LYZ_pH8_0.pdbqt
        """
        base = self.pdb_id if not ph_token else f"{self.pdb_id}_{ph_token}"
        return self.receptor_dir(variant) / f"{base}.pdbqt"

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
    # Historically some code used EXTRACTED_LIGANDS_DIR/<PDB>/ or LIGAND_EXTRACTED_DIR/<PDB>/.
    def legacy_extracted_ligands_dir(self) -> Optional[Path]:
        """
        Legacy (optional): <EXTRACTED_LIGANDS_DIR or LIGAND_EXTRACTED_DIR>/<PDB>/
        Note: canonical location is processed_pdbs/<PDB>/ligands_raw/.
        """
        return None  # kept as a stub to avoid reintroducing legacy in new code

    # ---------------------------
    # Docking outputs
    # ---------------------------
    def docked_pdb_root(self, ph_label: Optional[str] = None) -> Path:
        """
        Dir (legacy):        docked/<PDB>/
        Dir (with pH):      docked/<PDB>/<PH>/

        ``ph_label`` should be a filesystem-safe token (e.g., "pH6_0"). When
        omitted, the historical directory layout is preserved.
        """
        base = self.docked_root / self.pdb_id
        base.mkdir(parents=True, exist_ok=True)
        if ph_label is None:
            return base

        ph_dir = base / str(ph_label)
        ph_dir.mkdir(parents=True, exist_ok=True)
        return ph_dir

    def docking_score_long_csv(self, ph_label: Optional[str] = None) -> Path:
        """File: docked/<PDB>/<PH?>/docking_score_long.csv"""
        return self.docked_pdb_root(ph_label=ph_label) / "docking_score_long.csv"

    def docking_score_summary_csv(self, ph_label: Optional[str] = None) -> Path:
        """File: docked/<PDB>/<PH?>/docking_score_summary.csv"""
        return self.docked_pdb_root(ph_label=ph_label) / "docking_score_summary.csv"

    def bench_pocket_dir(self, pocket_index: int) -> Path:
        """
        Dir: docked/<PDB>/bench_pocketX_single/
        Example: docked/1T46/bench_pocket1_single/
        """
        d = self.docked_pdb_root() / f"bench_pocket{int(pocket_index)}_single"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def docked_variant_root(self, variant: Optional[str], ph_label: Optional[str] = None) -> Path:
        """
        If variant in {APO,HOLO}:
          Dir: docked/<PDB>/<VARIANT>/              (legacy)
               docked/<PDB>/<PH>/<VARIANT>/        (with pH)
        If variant is None (legacy):
          Dir: docked/<PDB>/                     (legacy)
               docked/<PDB>/<PH>/                (with pH)
        """
        v = _norm_variant(variant)
        base = self.docked_pdb_root(ph_label=ph_label)
        d = base / v if v else base
        d.mkdir(parents=True, exist_ok=True)
        return d

    def docked_stage_dir(
        self,
        variant: Optional[str],
        stage: Union[str, int],
        ph_label: Optional[str] = None,
    ) -> Path:
        """
        If variant in {APO,HOLO}:
          Dir: docked/<PDB>/<VARIANT>/<stage>/                 (legacy)
               docked/<PDB>/<PH>/<VARIANT>/<stage>/           (with pH)
        If variant is None:
          Dir: docked/<PDB>/<stage>/                          (legacy)
               docked/<PDB>/<PH>/<stage>/                    (with pH)
        Expected files: <ligand>_<stage>.pdbqt (poses), vina logs, per-stage CSVs
        """
        stage_token = stage if isinstance(stage, str) else f"stage{int(stage)}"
        d = self.docked_variant_root(variant, ph_label=ph_label) / str(stage_token)
        d.mkdir(parents=True, exist_ok=True)
        return d

    # ---------------------------
    # Config emission
    # ---------------------------
    def configs_stage_dir(self, run_id: str, variant: Optional[str], stage: str) -> Path:
        """
        If variant in {APO,HOLO}:
          Dir: configs/<RUN>/<PDB>/<VARIANT>/<stage>/
        If variant is None:
          Dir: configs/<RUN>/<PDB>/<stage>/
        Expected: <ligand>_<stage>.txt  (Vina config files per ligand)
        """
        v = _norm_variant(variant)
        base = self.configs_root / run_id / self.pdb_id
        d = (base / v / stage) if v else (base / stage)
        d.mkdir(parents=True, exist_ok=True)
        return d
    def expand_variants(self, mode: Optional[str]) -> list[Optional[str]]:
        # Delegate to the module-level function to keep a single source of truth
        return expand_variants(mode)


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
    pdb_id = re.sub(r'(?i)(_nolig(_cleaned)?|_cleaned)$', '', raw_id).upper()

    # Required roots
    over_root = Path(cfg["OVERALL_DIR"])
    input_root = Path(cfg["INPUT_DIR"])
    processed_root = Path(cfg["OUTPUT_DIR"])      # processed_pdbs/
    docked_root = Path(cfg["DOCKED_DIR"])         # docked/

    # Prepped ligands: pick the first present among synonymous keys
    prepped_root = Path(
        cfg.get("PREPPED_LIGANDS_DIR")
        or cfg.get("OUTPUT_LIGANDS_DIR")
        or cfg.get("PREPPED_LIGANDS_ROOT")
        or (over_root / "prepped_ligands")
    )

    ligands_mol2_root = Path(cfg.get("LIGANDS_MOL2_DIR", over_root / "ligands_mol2"))

    # Configs root: default under OVERALL_DIR/configs (no CONFIGS_DIR in your cfg)
    configs_root = Path(cfg.get("CONFIGS_DIR", over_root / "configs"))

    # Canonical per-PDB root and shared subfolders
    root_pdb_dir = processed_root / pdb_id  # processed_pdbs/<PDB>/
    raw_dir      = root_pdb_dir / "raw"            # processed_pdbs/<PDB>/raw/
    lig_raw_dir  = root_pdb_dir / "ligands_raw"    # processed_pdbs/<PDB>/ligands_raw/
    nolig_dir    = root_pdb_dir / "nolig"          # processed_pdbs/<PDB>/nolig/
    work_dir     = root_pdb_dir / "work"           # processed_pdbs/<PDB>/work/

    # Ensure shared subfolders exist
    for d in (root_pdb_dir, raw_dir, lig_raw_dir, nolig_dir, work_dir):
        d.mkdir(parents=True, exist_ok=True)

    # Inputs / shared outputs
    input_pdb_path = input_root / pdb_file                     # input_pdbs/<pdb_file>
    nolig_pdb_path = nolig_dir / f"{pdb_id}_nolig.pdb"         # processed_pdbs/<PDB>/nolig/<PDB>_nolig.pdb

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
