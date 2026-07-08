from __future__ import annotations

import logging
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Any, Dict

from config.tool_resolver import repo_ambertools_roots
from config.value_access import to_bool, to_float, to_int

_LOG = logging.getLogger(__name__)
_EMITTED_ALIAS_WARNINGS: set[tuple[str, str]] = set()
_EMITTED_PORTABILITY_WARNINGS: set[str] = set()
_VAR_PATTERN = re.compile(r"\$\{([A-Za-z0-9_]+)\}")

# Canonical key -> legacy aliases accepted from config/env.
_ALIASES: dict[str, tuple[str, ...]] = {
    "EXTRACTED_LIGANDS_DIR": ("LIGAND_EXTRACTED_DIR",),
    "PREPPED_LIGANDS_DIR": ("OUTPUT_LIGANDS_DIR", "PREPPED_LIGANDS_ROOT"),
    "PREPARE_LIGAND_SCRIPT": ("PREPARE_LIGAND4",),
    "MGLTOOLS_PATH": ("MGLTOOLS_DIR",),
    "PYMOL_EXE": ("PYMOL_PATH",),
    "REDUCE_EXE": ("REDUCE_LOCAL_CANDIDATE",),
    "SCORCH": ("SCORCH_SCRIPT",),
    "CPU": ("MAX_PARALLEL_JOBS",),
}

_OVERALL_DERIVED: dict[str, str] = {
    "INPUT_DIR": "input_pdbs",
    "LIGAND_DIR": "input_ligands",
    "EXTRACTED_LIGANDS_DIR": "extracted_ligands",
    "LIGANDS_MOL2_DIR": "ligands_mol2",
    "PREPPED_LIGANDS_DIR": "prepped_ligands",
    "OUTPUT_DIR": "outputs/processed_pdbs",
    "PDBQT_DIR": "pdbqts",
    "DOCKED_DIR": "outputs/docked",
    "POST_DOCKED_DIR": "outputs/post_docked",
    "P2RANK_OUTPUT_DIR": "outputs/p2rank_out",
    "CONFIGS_DIR": "outputs/configs",
    "DATA_DIR": "outputs/data",
    "LOGS_DIR": "outputs/logs",
    "MANIFESTS_DIR": "outputs/manifests",
    "PROTEIN_DIR": "pdbqts",
}

_BOOL_KEYS = {
    "CPU_ONLY",
    "FORCE_REPROCESS",
    "ALLOW_BOX_EXPAND",
    "QUIET_CONSOLE",
    "FILTER_INVALID",
    "USE_MEEKO",
    "USE_DOCK6",
    "use_dock6",
    "USE_GNINA",
    "USE_LEDOCK",
    "USE_SCORCH",
    "REPAIR_SCORCH_SELECTION",
    "ATLAS_REPAIR_SCORCH_SELECTION",
    "PH_ENSEMBLE",
    "RESET_CONFIGS",
    "FAST_MODE",
    "NO_LIBRARY_DOCKING",
    "CONTROL_CONSENSUS",
    "CHECKPOINT_ENABLE",
    "MMGBSA_STRIP_METALS",
    "MMGBSA_WATER_USE_ALIASES",
    "MMGBSA_METAL_USE_ALIASES",
    "MMGBSA_TOPOLOGY_PREP_ENABLED",
    "MMGBSA_TLEAP_RUN",
    "MMGBSA_CPPTRAJ_ENABLED",
    "MMGBSA_CPPTRAJ_RUN",
    "MMGBSA_MMPBSA_ENABLED",
    "MMGBSA_MMPBSA_RUN",
    "MMGBSA_MMPBSA_MPI_ENABLED",
    "MMGBSA_ENABLED",
    "MMGBSA_AUTO_RUN",
    "MMGBSA_MD_ENABLED",
    "MMGBSA_MMPBSA_USE_TRAJ_FRAMES",
    "MMGBSA_MD_RESTRAIN_PROTEIN_HEAVY",
    "MMGBSA_MD_COPY_BEST_REPLICATE",
    "MMGBSA_QC_ENABLED",
    "MMGBSA_QC_RUN_CPPTRAJ",
    "MMGBSA_CHEMISTRY_REVIEW_REQUIRED",
    "MMGBSA_PRODUCTION_STRICT",
    "MMGBSA_EXPLICIT_SOLVENT_STRIP_TRAJ",
    "MMGBSA_KEEP_WATERS",
    "MMGBSA_KEEP_METALS",
    "MMGBSA_RECEPTOR_FORCE",
    "MMGBSA_FORCE",
    "MMGBSA_STRICT",
    "MMGBSA_TLEAP_ENABLED",
    "MMGBSA_TLEAP_FORCE",
    "MMGBSA_LIGAND_BCC_SWEEP_INCLUDE_PLUSMINUS2",
    "MMGBSA_LIGAND_FORCE",
    "MMGBSA_RDKit_VALIDATE",
    "MMGBSA_LIGAND_CHEMISTRY_STRICT",
    "MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE",
    "MMGBSA_CHEMISTRY_SUBGATES_REQUIRED",
    "TOOL_VERIFY_ON_START",
    "ENABLE_GLOBAL_SCHEDULER",
}

_INT_KEYS = {
    "CPU",
    "GLOBAL_SCHEDULER_CPUS",
    "GLOBAL_MAX_PROTEINS",
    "GLOBAL_MIN_PROTEINS",
    "MAX_RECENTER_ATTEMPTS",
    "EARLY_RECENTER_MIN_EVAL",
    "CTRL_REDOCK_EXHAUSTIVENESS",
    "CTRL_REDOCK_NMODES",
    "CTRL_REDOCK_MAX_CONTROLS",
    "THREADS_PER_VINA_CTRL",
    "MMGBSA_TRAJ_STARTFRAME",
    "MMGBSA_TRAJ_ENDFRAME",
    "MMGBSA_TRAJ_INTERVAL",
    "MMGBSA_MMPBSA_STARTFRAME",
    "MMGBSA_MMPBSA_ENDFRAME",
    "MMGBSA_MMPBSA_INTERVAL",
    "MMGBSA_MMPBSA_VERBOSE",
    "MMGBSA_GB_IGB",
    "MMGBSA_LIGAND_SQM_LEVEL",
    "MMGBSA_RDKit_RADICAL_LOWCONF_THRESHOLD",
    "MMGBSA_MAX_LIGANDS",
    "MMGBSA_RERANKED_TOP_PCT",
    "MMGBSA_GENERAL_STARTFRAME",
    "MMGBSA_GENERAL_ENDFRAME",
    "MMGBSA_GENERAL_INTERVAL",
    "MMGBSA_GENERAL_VERBOSE",
    "MMGBSA_MD_NREPLICATES",
    "MMGBSA_MD_MPI_RANKS",
    "MMGBSA_MMPBSA_MPI_RANKS",
    "MMGBSA_OPENMM_CPU_THREADS",
    "MMGBSA_ANALYSIS_STARTFRAME",
    "MMGBSA_ANALYSIS_ENDFRAME",
    "MMGBSA_ANALYSIS_INTERVAL",
    "MMGBSA_QC_MIN_FRAMES",
    "MMGBSA_QC_MIN_BLOCKS",
    "MMGBSA_QC_BLOCK_SIZE",
}

_FLOAT_KEYS = {
    "EARLY_RECENTER_RATIO",
    "EARLY_RECENTER_FAR_A",
    "EARLY_RECENTER_MEDIAN_A",
    "CONTROL_CENTER_CLOSE_MAX_A",
    "MMGBSA_WATER_KEEP_RADIUS_A",
    "MMGBSA_WATER_KEEP_RADIUS",
    "MMGBSA_ACTIVE_SITE_RADIUS_FALLBACK",
    "MMGBSA_RERANKED_TOP_PCT",
    "MMGBSA_GB_SALTCON",
    "MMGBSA_MD_DT_PS",
    "MMGBSA_MD_HEAT_PS",
    "MMGBSA_MD_EQUIL_PS",
    "MMGBSA_MD_PROD_PS",
    "MMGBSA_MD_FRAME_STRIDE_PS",
    "MMGBSA_MD_TEMP0",
    "MMGBSA_MD_GAMMA_LN",
    "MMGBSA_MD_RESTRAINT_WT",
    "MMGBSA_ANALYSIS_START_PS",
    "MMGBSA_ANALYSIS_END_PS",
    "MMGBSA_QC_MAX_SEM_KCAL",
    "MMGBSA_QC_MAX_BLOCK_RANGE_KCAL",
}


def which_or_exists(candidates) -> str | None:
    """
    Return the first path that exists (if absolute),
    or the first candidate found on PATH.
    """
    for cand in candidates:
        path = Path(str(cand))
        if path.is_absolute() and path.exists():
            return str(path)
        env_candidate = Path(sys.executable).resolve().parent / path.name
        if env_candidate.exists():
            return str(env_candidate)
        found = shutil.which(path.name)
        if found:
            return found
    return None


def default_base_dir(reference_file: str | os.PathLike[str] | None = None) -> Path:
    """
    Base directory for default config resolution.

    Precedence:
    1) PROTEIN_AUTOMATION_DIR env var
    2) parent directory of `reference_file` if provided
    3) repository root inferred from this module location
    """
    env_dir = os.environ.get("PROTEIN_AUTOMATION_DIR")
    if env_dir:
        return Path(env_dir)
    if reference_file is not None:
        return Path(reference_file).resolve().parent
    return Path(__file__).resolve().parents[2]


def _mgltools_root_from_pythonsh(pythonsh: str | None) -> Path | None:
    if _is_blank(pythonsh):
        return None
    py_path = Path(str(pythonsh)).expanduser()
    if not py_path.exists():
        return None
    for parent in py_path.resolve().parents:
        if (parent / "MGLToolsPckgs").exists() or (parent / "AutoDockTools").exists():
            return parent
    return None


def _mgltools_utility(root: Path | None, script_name: str) -> str | None:
    if root is None:
        return None
    candidate = root / "MGLToolsPckgs" / "AutoDockTools" / "Utilities24" / script_name
    return str(candidate) if candidate.exists() else None


def default_dirs(base: Path) -> Dict[str, str]:
    """
    Minimal directory defaults; derived paths are computed in normalize_config().
    """
    return {
        "OVERALL_DIR": str(base),
        "PHENIX_CLEAN_SCRIPT": str(base / "src" / "cli" / "phenix_clean_cli.py"),
    }


def default_tools(base: Path) -> Dict[str, Any]:
    """
    Tool defaults prefer environment and PATH discovery.
    """
    scorch_root = base.parent.parent / "tools" / "SCORCH"
    scorch_script_default: str | None
    if (scorch_root / "scorch.py").exists():
        scorch_script_default = str(scorch_root / "scorch.py")
    else:
        scorch_script_default = which_or_exists(["scorch.py"])
    mgltools_python = (
        os.environ.get("MGL_PYTHON")
        or os.environ.get("MGLTOOLS_PYTHON")
        or which_or_exists(["pythonsh"])
    )
    mgltools_root = (
        Path(os.environ["MGLTOOLS_PATH"])
        if os.environ.get("MGLTOOLS_PATH")
        else Path(os.environ["MGLTOOLS_DIR"])
        if os.environ.get("MGLTOOLS_DIR")
        else _mgltools_root_from_pythonsh(mgltools_python)
    )

    return {
        "VINA_PATH": which_or_exists(["vina"]),
        "VINA_EXE": which_or_exists(["vina"]),
        "GNINA_EXE": which_or_exists(["gnina"]),
        "OPENBABEL_PATH": which_or_exists(["obabel"]),
        "PROPKA_EXE": os.environ.get("PROPKA_EXE")
        or which_or_exists(["propka3"]),
        "PYMOL_EXE": which_or_exists(["pymol"]),
        "REDUCE_EXE": which_or_exists(["reduce"]),
        "PDB2PQR_EXE": which_or_exists(["pdb2pqr"]),
        "DOCK6_EXE": which_or_exists(["dock6"]),
        "MGLTOOLS_PATH": str(mgltools_root) if mgltools_root else None,
        "MGLTOOLS_PYTHON": mgltools_python,
        "PREPARE_LIGAND_SCRIPT": os.environ.get("PREPARE_LIGAND_SCRIPT")
        or which_or_exists(["prepare_ligand4.py"])
        or _mgltools_utility(mgltools_root, "prepare_ligand4.py"),
        "PREPARE_RECEPTOR_SCRIPT": os.environ.get("PREPARE_RECEPTOR_SCRIPT")
        or which_or_exists(["prepare_receptor4.py"])
        or _mgltools_utility(mgltools_root, "prepare_receptor4.py"),
        "P2RANK_PATH": shutil.which("prank") or "prank",
        "SCORCH": os.environ.get("SCORCH")
        or os.environ.get("SCORCH_SCRIPT")
        or scorch_script_default,
        "SCORCH_ENV_PREFIX": os.environ.get("SCORCH_ENV_PREFIX"),
        "SCORCH_ENV": os.environ.get("SCORCH_ENV") or "scorch-env",
    }


def default_runtime() -> Dict[str, Any]:
    return {
        "CPU_ONLY": True,
        "CPU": os.cpu_count() or 8,
        "ENABLE_GLOBAL_SCHEDULER": True,
        "GLOBAL_SCHEDULER_CPUS": os.cpu_count() or 8,
        "GLOBAL_MAX_PROTEINS": 0,
        "GLOBAL_MIN_PROTEINS": 1,
        "GLOBAL_ADMISSION_POLICY": "adaptive",
        "FORCE_REPROCESS": False,
        # CLI-only gate in main.py; default remains false in config.
        "TOOL_VERIFY_ON_START": False,
        "DOCKING_MODE": "discovery",
        "QUIET_CONSOLE": False,
        "FILTER_INVALID": False,
        # docking/recenter knobs
        "EARLY_RECENTER_RATIO": 0.70,
        "EARLY_RECENTER_MIN_EVAL": 10,
        "EARLY_RECENTER_FAR_A": 15.0,
        "EARLY_RECENTER_MEDIAN_A": 10.0,
        "ALLOW_BOX_EXPAND": True,
        "MAX_RECENTER_ATTEMPTS": 3,
        # control-centering knobs
        "CONTROL_CENTER_POLICY": "best_redock",
        "CONTROL_CENTER_CLOSE_MAX_A": 8.0,
        "CTRL_REDOCK_EXHAUSTIVENESS": 32,
        "CTRL_REDOCK_NMODES": 9,
        "CTRL_REDOCK_MAX_CONTROLS": 12,
        "THREADS_PER_VINA_CTRL": 4,
        # --- benchmark policy knobs ---
        "BENCH_ENFORCE_HARDCODED_CONTROLS_ONLY": True,
        "BENCH_ALLOW_WHITELIST_FALLBACK_IF_CONTROLS_MISSING": False,
        # --- MMGBSA receptor prep ---
        "MMGBSA_STRIP_METALS": True,
        "MMGBSA_WATER_POLICY": "ACTIVE_SITE",
        "MMGBSA_WATER_KEEP_RADIUS_A": 6.0,
        "MMGBSA_WATER_USE_ALIASES": True,
        "MMGBSA_METAL_USE_ALIASES": True,
        "MMGBSA_TOPOLOGY_PREP_ENABLED": True,
        "MMGBSA_TLEAP_RUN": True,
        "MMGBSA_TOPOLOGY_DIRNAME": "mmgbsa_topologies",
        "MMGBSA_AMBERTOOLS_PREFIX": "",
        "MMGBSA_CPPTRAJ_ENABLED": True,
        "MMGBSA_CPPTRAJ_RUN": True,
        "MMGBSA_TRAJOUT_NAME": "mdcrd",
        "MMGBSA_TRAJOUT_FORMAT": "mdcrd",
        "MMGBSA_TRAJIN_SOURCE": "INPCRD",
        "MMGBSA_TRAJIN_PATH": "",
        "MMGBSA_TRAJIN_FORMAT": "inpcrd",
        "MMGBSA_TRAJ_STARTFRAME": 1,
        "MMGBSA_TRAJ_ENDFRAME": 1,
        "MMGBSA_TRAJ_INTERVAL": 1,
        "MMGBSA_MMPBSA_ENABLED": True,
        "MMGBSA_MMPBSA_RUN": True,
        "MMGBSA_MMPBSA_STARTFRAME": 1,
        "MMGBSA_MMPBSA_ENDFRAME": 1,
        "MMGBSA_MMPBSA_INTERVAL": 1,
        "MMGBSA_MMPBSA_VERBOSE": 2,
        "MMGBSA_GB_IGB": 5,
        "MMGBSA_GB_SALTCON": 0.150,
        "MMGBSA_MMPBSA_INPUT_NAME": "mmpbsa.in",
        "MMGBSA_MMPBSA_LOG_NAME": "mmpbsa.log",
        "MMGBSA_MMPBSA_OUT_DAT": "FINAL_RESULTS_MMPBSA.dat",
        "MMGBSA_MMPBSA_OUT_CSV": "FINAL_RESULTS_MMPBSA.csv",
        "MMGBSA_DEFAULT_TRAJ_NAME": "mdcrd",
        "MMGBSA_ENABLED": True,
        "MMGBSA_AUTO_RUN": False,
        "MMGBSA_PROTOCOL": "screening",
        "MMGBSA_PRODUCTION_STRICT": True,
        "MMGBSA_MD_ENABLED": False,
        "MMGBSA_MD_ENGINE": "sander",
        "MMGBSA_OPENMM_PLATFORM": "auto",
        "MMGBSA_OPENMM_DEVICE_INDEX": "",
        "MMGBSA_OPENMM_DEFAULT_DEVICE_INDEX": "0",
        "MMGBSA_OPENMM_PRECISION": "mixed",
        "MMGBSA_OPENMM_CPU_THREADS": 1,
        "MMGBSA_OPENMM_START_STAGE": "prod",
        "MMGBSA_MD_SOLVENT_MODEL": "implicit",
        "MMGBSA_MD_DT_PS": 0.002,
        "MMGBSA_MD_HEAT_PS": 10.0,
        "MMGBSA_MD_EQUIL_PS": 20.0,
        "MMGBSA_MD_PROD_PS": 100.0,
        "MMGBSA_MD_FRAME_STRIDE_PS": 2.0,
        "MMGBSA_MD_NREPLICATES": 1,
        "MMGBSA_MD_REP_AGG": "mean",
        "MMGBSA_MD_COPY_BEST_REPLICATE": False,
        "MMGBSA_MD_TEMP0": 300.0,
        "MMGBSA_MD_NTT": 3,
        "MMGBSA_MD_GAMMA_LN": 2.0,
        "MMGBSA_MD_NTC": 2,
        "MMGBSA_MD_NTF": 2,
        "MMGBSA_MD_IGB": 5,
        "MMGBSA_MD_SALTCON": 0.150,
        "MMGBSA_MD_RESTRAIN_PROTEIN_HEAVY": True,
        "MMGBSA_MD_RESTRAINT_WT": 5.0,
        "MMGBSA_MD_RESTRAINT_MASK": ":1-999999 & !@H=",
        "MMGBSA_MMPBSA_USE_TRAJ_FRAMES": False,
        "MMGBSA_ANALYSIS_START_PS": 0.0,
        "MMGBSA_ANALYSIS_END_PS": 0.0,
        "MMGBSA_ANALYSIS_STARTFRAME": 0,
        "MMGBSA_ANALYSIS_ENDFRAME": 0,
        "MMGBSA_ANALYSIS_INTERVAL": 1,
        "MMGBSA_QC_ENABLED": True,
        "MMGBSA_QC_RUN_CPPTRAJ": True,
        "MMGBSA_QC_MIN_FRAMES": 1,
        "MMGBSA_QC_MIN_BLOCKS": 1,
        "MMGBSA_QC_BLOCK_SIZE": 0,
        "MMGBSA_QC_MAX_SEM_KCAL": 999999.0,
        "MMGBSA_QC_MAX_BLOCK_RANGE_KCAL": 999999.0,
        "MMGBSA_CHEMISTRY_REVIEW_REQUIRED": False,
        "MMGBSA_CHEMISTRY_REVIEW_STATUS": "unreviewed",
        "MMGBSA_CHEMISTRY_REVIEW_NOTES": "",
        "MMGBSA_CHEMISTRY_SUBGATES_REQUIRED": False,
        "MMGBSA_PROTONATION_REVIEW_STATUS": "unreviewed",
        "MMGBSA_TAUTOMER_REVIEW_STATUS": "unreviewed",
        "MMGBSA_STEREOCHEMISTRY_REVIEW_STATUS": "unreviewed",
        "MMGBSA_NET_CHARGE_REVIEW_STATUS": "unreviewed",
        "MMGBSA_PARAMETER_REVIEW_STATUS": "unreviewed",
        "MMGBSA_EXPLICIT_TRAJ_PATH": "",
        "MMGBSA_EXPLICIT_SOLVENT_STRIP_TRAJ": True,
        "MMGBSA_KEEP_WATERS": True,
        "MMGBSA_WATER_KEEP_RADIUS": 6.0,
        "MMGBSA_KEEP_METALS": False,
        "MMGBSA_METAL_RETAIN_TOKENS": "",
        "MMGBSA_WATER_RETAIN_TOKENS": "",
        "MMGBSA_RECEPTOR_FORCE": False,
        "MMGBSA_FORCE": False,
        "MMGBSA_STRICT": False,
        "MMGBSA_LIGAND_AT": "gaff2",
        "MMGBSA_LIGAND_CHARGE_METHOD": "bcc",
        "MMGBSA_LIGAND_PRIMARY_CHARGE_METHOD": "bcc",
        "MMGBSA_LIGAND_FALLBACK_CHARGE_METHOD": "gas",
        "MMGBSA_LIGAND_NOMINAL_NET_CHARGE": "auto",
        "MMGBSA_LIGAND_BCC_CHARGE_SWEEP": "-1,1",
        "MMGBSA_LIGAND_BCC_SWEEP_INCLUDE_PLUSMINUS2": True,
        "MMGBSA_LIGAND_FORCE": False,
        "MMGBSA_RDKit_VALIDATE": True,
        "MMGBSA_RDKit_RADICAL_LOWCONF_THRESHOLD": 1,
        "MMGBSA_LIGAND_CHEMISTRY_STRICT": False,
        "MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE": False,
        "MMGBSA_INPUT_CHEMISTRY_SOURCE": "unreviewed_sdf",
        "MMGBSA_INPUT_CHEMISTRY_NOTES": "",
        "MMGBSA_LIGAND_NET_CHARGE": 0,
        "MMGBSA_LIGAND_SQM_LEVEL": 2,
        "MMGBSA_INPUT_STAGE_DIR": "stage1",
        "MMGBSA_MAX_LIGANDS": 1,
        "MMGBSA_RERANKED_TOP_PCT": 0.0,
        "MMGBSA_ACTIVE_SITE_RADIUS_FALLBACK": 6.0,
        "MMGBSA_TLEAP_ENABLED": True,
        "MMGBSA_TLEAP_FORCE": False,
        "MMGBSA_GENERAL_STARTFRAME": 1,
        "MMGBSA_GENERAL_ENDFRAME": 1,
        "MMGBSA_GENERAL_INTERVAL": 1,
        "MMGBSA_GENERAL_VERBOSE": 2,
    }


def _is_blank(value: Any) -> bool:
    return value is None or (isinstance(value, str) and value.strip() == "")


def _warn_alias(alias: str, canonical: str) -> None:
    pair = (alias, canonical)
    if pair in _EMITTED_ALIAS_WARNINGS:
        return
    _EMITTED_ALIAS_WARNINGS.add(pair)
    _LOG.warning(
        "Config key '%s' is deprecated; use '%s' instead.", alias, canonical
    )


def _warn_portability_once(key: str, message: str) -> None:
    if key in _EMITTED_PORTABILITY_WARNINGS:
        return
    _EMITTED_PORTABILITY_WARNINGS.add(key)
    _LOG.warning(message)


def _to_bool(value: Any, default: Any = None) -> Any:
    return to_bool(value, default=default)


def _to_int(value: Any, default: Any = None) -> Any:
    return to_int(value, default=default)


def _to_float(value: Any, default: Any = None) -> Any:
    return to_float(value, default=default)


def _explicit_keys(cfg: Dict[str, Any]) -> set[str]:
    explicit: set[str] = set()
    file_cfg = cfg.get("_FILE_CFG")
    has_file_cfg = False
    if isinstance(file_cfg, dict):
        has_file_cfg = True
        for key, value in file_cfg.items():
            if not _is_blank(value):
                explicit.add(str(key).upper())
    tracked = set(_ALIASES.keys())
    for aliases in _ALIASES.values():
        tracked.update(aliases)
    tracked.update(_OVERALL_DERIVED.keys())
    tracked.update(
        {
            "MGLTOOLS_PYTHON",
            "PREPARE_RECEPTOR_SCRIPT",
            "PREPARE_LIGAND_SCRIPT",
            "MMGBSA_AMBERTOOLS_PREFIX",
            "AMBERTOOLS_PREFIX",
            "MMGBSA_PROTOCOL",
            "MMGBSA_MD_ENABLED",
            "MMGBSA_MD_SOLVENT_MODEL",
            "MMGBSA_MD_NREPLICATES",
            "MMGBSA_MD_HEAT_PS",
            "MMGBSA_MD_EQUIL_PS",
            "MMGBSA_MD_PROD_PS",
            "MMGBSA_MD_FRAME_STRIDE_PS",
            "MMGBSA_MMPBSA_USE_TRAJ_FRAMES",
            "MMGBSA_ANALYSIS_START_PS",
            "MMGBSA_ANALYSIS_END_PS",
            "MMGBSA_ANALYSIS_STARTFRAME",
            "MMGBSA_ANALYSIS_ENDFRAME",
            "MMGBSA_ANALYSIS_INTERVAL",
            "MMGBSA_QC_MIN_FRAMES",
            "MMGBSA_QC_MIN_BLOCKS",
            "MMGBSA_QC_MAX_SEM_KCAL",
            "MMGBSA_QC_MAX_BLOCK_RANGE_KCAL",
            "MMGBSA_CHEMISTRY_REVIEW_REQUIRED",
            "MMGBSA_CHEMISTRY_REVIEW_STATUS",
            "MMGBSA_PRODUCTION_STRICT",
            "MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE",
            "MMGBSA_INPUT_CHEMISTRY_SOURCE",
            "MMGBSA_INPUT_CHEMISTRY_NOTES",
        }
    )
    if not has_file_cfg:
        for key in tracked:
            if key in cfg and not _is_blank(cfg.get(key)):
                explicit.add(key)
    for key in tracked:
        env_val = os.environ.get(key)
        if env_val is not None and str(env_val).strip() != "":
            explicit.add(key)
    return explicit


def _adopt_aliases(cfg: Dict[str, Any]) -> None:
    for canonical, aliases in _ALIASES.items():
        canonical_value = cfg.get(canonical)
        if _is_blank(canonical_value):
            for alias in aliases:
                alias_value = cfg.get(alias)
                if _is_blank(alias_value):
                    continue
                cfg[canonical] = alias_value
                _warn_alias(alias, canonical)
                break
        else:
            for alias in aliases:
                alias_value = cfg.get(alias)
                if not _is_blank(alias_value):
                    _warn_alias(alias, canonical)


def _derive_from_overall(cfg: Dict[str, Any]) -> None:
    overall = cfg.get("OVERALL_DIR")
    if _is_blank(overall):
        return
    base = Path(str(overall)).expanduser()
    cfg["OVERALL_DIR"] = str(base)
    for key, subdir in _OVERALL_DERIVED.items():
        if _is_blank(cfg.get(key)):
            cfg[key] = str(base / subdir)
    if _is_blank(cfg.get("PHENIX_CLEAN_SCRIPT")):
        cfg["PHENIX_CLEAN_SCRIPT"] = str(base / "src" / "cli" / "phenix_clean_cli.py")


def _repair_missing_overall_dir(cfg: Dict[str, Any], explicit: set[str]) -> None:
    overall = cfg.get("OVERALL_DIR")
    if _is_blank(overall):
        cfg["OVERALL_DIR"] = str(default_base_dir())
        return

    overall_path = Path(str(overall)).expanduser()
    if overall_path.exists():
        cfg["OVERALL_DIR"] = str(overall_path)
        return

    fallback = default_base_dir()
    _LOG.warning(
        "Configured OVERALL_DIR does not exist (%s); falling back to %s",
        overall_path,
        fallback,
    )
    cfg["OVERALL_DIR"] = str(fallback)

    # Treat repaired OVERALL_DIR as non-explicit for derived path keys so they can
    # follow the corrected repo root when stale machine-specific values are present.
    explicit.discard("OVERALL_DIR")
    for key in _OVERALL_DERIVED:
        explicit.discard(key)
    explicit.discard("PHENIX_CLEAN_SCRIPT")


def _derive_mgltools(cfg: Dict[str, Any], explicit: set[str]) -> None:
    mgl_root = cfg.get("MGLTOOLS_PATH")
    mgl_explicit = "MGLTOOLS_PATH" in explicit or "MGLTOOLS_DIR" in explicit
    if not _is_blank(mgl_root):
        mgl_path = Path(str(mgl_root)).expanduser()
        cfg["MGLTOOLS_PATH"] = str(mgl_path)
        if "MGLTOOLS_PYTHON" in explicit or "PREPARE_RECEPTOR_SCRIPT" in explicit:
            _warn_portability_once(
                "mgltools_derived",
                "MGLTOOLS_PATH is set; deriving MGLTOOLS_PYTHON and PREPARE_RECEPTOR_SCRIPT from it for portability.",
            )
        # Canonical behavior: derive ADT paths from MGLTOOLS_PATH for portability.
        cfg["MGLTOOLS_PYTHON"] = str(mgl_path / "bin" / "pythonsh")
        cfg["PREPARE_RECEPTOR_SCRIPT"] = str(
            mgl_path
            / "MGLToolsPckgs"
            / "AutoDockTools"
            / "Utilities24"
            / "prepare_receptor4.py"
        )
        if (
            "PREPARE_LIGAND_SCRIPT" not in explicit
            and "PREPARE_LIGAND4" not in explicit
        ):
            cfg["PREPARE_LIGAND_SCRIPT"] = str(
                mgl_path
                / "MGLToolsPckgs"
                / "AutoDockTools"
                / "Utilities24"
                / "prepare_ligand4.py"
            )
    else:
        if _is_blank(cfg.get("MGLTOOLS_PATH")) and mgl_explicit:
            cfg["MGLTOOLS_PATH"] = str(mgl_root).strip()

    # Runtime env variables intentionally remain highest-priority overrides.
    env_mgl_python = os.environ.get("MGLTOOLS_PYTHON") or os.environ.get("MGL_PYTHON")
    if not _is_blank(env_mgl_python):
        cfg["MGLTOOLS_PYTHON"] = env_mgl_python
    env_prepare_receptor = os.environ.get("PREPARE_RECEPTOR_SCRIPT")
    if not _is_blank(env_prepare_receptor):
        cfg["PREPARE_RECEPTOR_SCRIPT"] = env_prepare_receptor

    mgl_python = cfg.get("MGLTOOLS_PYTHON")
    if _is_blank(cfg.get("MGLTOOLS_PATH")) and not _is_blank(mgl_python):
        inferred_root = _mgltools_root_from_pythonsh(str(mgl_python))
        if inferred_root is not None:
            cfg["MGLTOOLS_PATH"] = str(inferred_root)
            receptor_script = _mgltools_utility(inferred_root, "prepare_receptor4.py")
            ligand_script = _mgltools_utility(inferred_root, "prepare_ligand4.py")
            if (
                "PREPARE_RECEPTOR_SCRIPT" not in explicit
                and _is_blank(cfg.get("PREPARE_RECEPTOR_SCRIPT"))
                and receptor_script
            ):
                cfg["PREPARE_RECEPTOR_SCRIPT"] = receptor_script
            if (
                "PREPARE_LIGAND_SCRIPT" not in explicit
                and "PREPARE_LIGAND4" not in explicit
                and _is_blank(cfg.get("PREPARE_LIGAND_SCRIPT"))
                and ligand_script
            ):
                cfg["PREPARE_LIGAND_SCRIPT"] = ligand_script

    if _is_blank(cfg.get("MGLTOOLS_PYTHON")):
        cfg["MGLTOOLS_PYTHON"] = (
            env_mgl_python
            or shutil.which("pythonsh")
        )
    if _is_blank(cfg.get("PREPARE_RECEPTOR_SCRIPT")):
        cfg["PREPARE_RECEPTOR_SCRIPT"] = (
            env_prepare_receptor
            or shutil.which("prepare_receptor4.py")
        )
    if _is_blank(cfg.get("PREPARE_LIGAND_SCRIPT")):
        cfg["PREPARE_LIGAND_SCRIPT"] = (
            os.environ.get("PREPARE_LIGAND_SCRIPT")
            or shutil.which("prepare_ligand4.py")
        )


def _derive_reduce(cfg: Dict[str, Any]) -> None:
    if _is_blank(cfg.get("REDUCE_EXE")):
        cfg["REDUCE_EXE"] = shutil.which("reduce")

    reduce_exe = cfg.get("REDUCE_EXE")
    if _is_blank(cfg.get("REDUCE_LOCAL_CANDIDATE")) and not _is_blank(reduce_exe):
        cfg["REDUCE_LOCAL_CANDIDATE"] = reduce_exe

    if not _is_blank(cfg.get("REDUCE_HET_DICT")):
        return

    if _is_blank(reduce_exe):
        return

    exe_path = Path(str(reduce_exe)).expanduser()
    parent = exe_path.parent
    candidates = [
        parent / "reduce_wwPDB_het_dict.txt",
        parent.parent / "reduce_wwPDB_het_dict.txt",
        parent.parent.parent / "reduce_wwPDB_het_dict.txt",
    ]
    for cand in candidates:
        if cand.exists():
            cfg["REDUCE_HET_DICT"] = str(cand)
            return
    # Deterministic fallback for common tree: <root>/reduce_src/reduce
    cfg["REDUCE_HET_DICT"] = str(parent.parent / "reduce_wwPDB_het_dict.txt")


def _derive_dock6(cfg: Dict[str, Any]) -> None:
    if _is_blank(cfg.get("DOCK6_EXE")):
        cfg["DOCK6_EXE"] = shutil.which("dock6")

    dock6_exe = cfg.get("DOCK6_EXE")
    if _is_blank(dock6_exe):
        for key, tool in [
            ("DMS_EXE", "dms"),
            ("SPHGEN_EXE", "sphgen"),
            ("SPHERE_SELECTOR_EXE", "sphere_selector"),
            ("SHOWBOX_EXE", "showbox"),
            ("GRID_EXE", "grid"),
        ]:
            if _is_blank(cfg.get(key)):
                cfg[key] = shutil.which(tool)
        return

    exe_path = Path(str(dock6_exe)).expanduser()
    bin_dir = exe_path.parent
    dock6_root = bin_dir.parent

    derived = {
        "SPHGEN_EXE": bin_dir / "sphgen",
        "SPHERE_SELECTOR_EXE": bin_dir / "sphere_selector",
        "SHOWBOX_EXE": bin_dir / "showbox",
        "GRID_EXE": bin_dir / "grid",
        "DOCK6_VDW_DEFN_FILE": dock6_root / "parameters" / "vdw_AMBER_parm99.defn",
    }
    for key, path in derived.items():
        if _is_blank(cfg.get(key)):
            cfg[key] = str(path)

    if _is_blank(cfg.get("DMS_EXE")):
        sibling = dock6_root.parent / "dms" / "bin" / "dms"
        cfg["DMS_EXE"] = str(sibling) if sibling.exists() else (shutil.which("dms") or str(sibling))


def _candidate_ambertools_prefixes(cfg: Dict[str, Any]) -> list[Path]:
    candidates: list[Path] = list(repo_ambertools_roots(cfg))
    candidates.append(Path.home() / "micromamba" / "envs" / "AmberTools25")
    seen: set[str] = set()
    out: list[Path] = []
    for path in candidates:
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        out.append(path)
    return out


def _derive_ambertools(cfg: Dict[str, Any], explicit: set[str]) -> None:
    if "MMGBSA_AMBERTOOLS_PREFIX" in explicit or "AMBERTOOLS_PREFIX" in explicit:
        return
    if not _is_blank(cfg.get("MMGBSA_AMBERTOOLS_PREFIX")) or not _is_blank(
        cfg.get("AMBERTOOLS_PREFIX")
    ):
        return
    required = ("MMPBSA.py", "cpptraj", "tleap", "antechamber", "parmchk2")
    for prefix in _candidate_ambertools_prefixes(cfg):
        if all((prefix / "bin" / tool).is_file() for tool in required):
            cfg["MMGBSA_AMBERTOOLS_PREFIX"] = str(prefix)
            cfg["AMBERTOOLS_PREFIX"] = str(prefix)
            return


def _set_if_not_explicit(
    cfg: Dict[str, Any], explicit: set[str], key: str, value: Any
) -> None:
    if key.upper() not in explicit:
        cfg[key] = value


def _derive_mmgbsa_protocol(cfg: Dict[str, Any], explicit: set[str]) -> None:
    protocol = str(cfg.get("MMGBSA_PROTOCOL", "screening") or "screening").strip().lower()
    if protocol in {"production", "publication", "publish"}:
        cfg["MMGBSA_PROTOCOL"] = "production"
        defaults = {
            "MMGBSA_MD_ENABLED": True,
            "MMGBSA_MD_SOLVENT_MODEL": "explicit",
            "MMGBSA_MD_NREPLICATES": 3,
            "MMGBSA_MD_MIN_STEPS": 10000,
            "MMGBSA_MD_RESCUE_ENABLED": True,
            "MMGBSA_MD_RESCUE_STEPS": 50000,
            "MMGBSA_MD_HEAT_PS": 100.0,
            "MMGBSA_MD_EQUIL_PS": 1000.0,
            "MMGBSA_MD_PROD_PS": 10000.0,
            "MMGBSA_MD_FRAME_STRIDE_PS": 10.0,
            "MMGBSA_MMPBSA_USE_TRAJ_FRAMES": True,
            "MMGBSA_KEEP_WATERS": True,
            "MMGBSA_WATER_POLICY": "ACTIVE_SITE",
            "MMGBSA_WATER_KEEP_RADIUS_A": 4.0,
            "MMGBSA_WATER_KEEP_RADIUS": 4.0,
            "MMGBSA_TLEAP_WATER_MODEL": "tip3p",
            "MMGBSA_EXPLICIT_WATER_BUFFER_A": 8.0,
            "MMGBSA_ANALYSIS_START_PS": 2000.0,
            "MMGBSA_ANALYSIS_END_PS": 0.0,
            "MMGBSA_ANALYSIS_INTERVAL": 1,
            "MMGBSA_FRAME_AGG": "mean",
            "MMGBSA_FRAME_LIMIT": 0,
            "MMGBSA_MD_REP_AGG": "mean",
            "MMGBSA_QC_MIN_FRAMES": 200,
            "MMGBSA_QC_MIN_BLOCKS": 5,
            "MMGBSA_QC_MAX_SEM_KCAL": 2.0,
            "MMGBSA_QC_MAX_BLOCK_RANGE_KCAL": 5.0,
            "MMGBSA_CHEMISTRY_REVIEW_REQUIRED": True,
            "MMGBSA_CHEMISTRY_SUBGATES_REQUIRED": True,
            "MMGBSA_LIGAND_CHEMISTRY_STRICT": True,
            "MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE": False,
            "MMGBSA_INPUT_CHEMISTRY_SOURCE": "unreviewed_sdf",
            "MMGBSA_INPUT_CHEMISTRY_NOTES": "",
            "MMGBSA_PRODUCTION_STRICT": True,
            "MMGBSA_LIGAND_ANTECHAMBER_TIMEOUT_SEC": 1800,
        }
        for key, value in defaults.items():
            _set_if_not_explicit(cfg, explicit, key, value)
    elif protocol in {"screen", "screening", "triage", "fast"}:
        cfg["MMGBSA_PROTOCOL"] = "screening"
    else:
        raise ValueError(
            f"MMGBSA_PROTOCOL must be screening or production (got {protocol})"
        )


def _repair_missing_tool_paths(cfg: Dict[str, Any]) -> None:
    tool_keys = {
        "VINA_EXE": "vina",
        "GNINA_EXE": "gnina",
        "OPENBABEL_PATH": "obabel",
        "PYMOL_EXE": "pymol",
        "REDUCE_EXE": "reduce",
        "DOCK6_EXE": "dock6",
        "DMS_EXE": "dms",
        "SPHGEN_EXE": "sphgen",
        "SPHERE_SELECTOR_EXE": "sphere_selector",
        "SHOWBOX_EXE": "showbox",
        "GRID_EXE": "grid",
        "LEPRO_EXE": "lepro",
        "PROPKA_EXE": "propka3",
        "P2RANK_PATH": "prank",
        "MGLTOOLS_PYTHON": "pythonsh",
        "PREPARE_RECEPTOR_SCRIPT": "prepare_receptor4.py",
        "PREPARE_LIGAND_SCRIPT": "prepare_ligand4.py",
        "SCORCH": "scorch.py",
    }
    which_cache: dict[str, str | None] = {}

    for key, binary in tool_keys.items():
        raw = cfg.get(key)
        if _is_blank(raw):
            continue
        text = str(raw).strip()
        path = Path(text).expanduser()
        if not path.is_absolute():
            continue
        if path.exists():
            cfg[key] = str(path)
            continue
        if key == "VINA_EXE" and path.name == "vina1.sh":
            vina_resolved = shutil.which(binary) or binary
            _LOG.warning(
                "Configured VINA_EXE uses retired wrapper (%s); falling back to %s",
                path,
                vina_resolved,
            )
            cfg[key] = vina_resolved
            continue
        if binary not in which_cache:
            which_cache[binary] = shutil.which(binary)
        resolved = which_cache[binary]
        if resolved:
            _LOG.warning(
                "Configured %s path does not exist (%s); falling back to PATH (%s)",
                key,
                path,
                resolved,
            )
            cfg[key] = resolved


def _expand_value(value: str, cfg: Dict[str, Any]) -> str:
    def _replace(match: re.Match[str]) -> str:
        key = (match.group(1) or "").upper()
        repl = cfg.get(key)
        if repl is None:
            return match.group(0)
        return str(repl)

    out = _VAR_PATTERN.sub(_replace, value)
    # Back-compat for older templates.
    over = cfg.get("OVERALL_DIR")
    if over is not None:
        over_s = str(over)
        out = out.replace("{OVERALL_DIR}", over_s)
        out = out.replace("$OVERALL_DIR", over_s)
    return out


def _expand_all(cfg: Dict[str, Any]) -> None:
    for _ in range(5):
        changed = False
        for key, value in list(cfg.items()):
            if not isinstance(value, str):
                continue
            new_val = _expand_value(value, cfg)
            if new_val != value:
                cfg[key] = new_val
                changed = True
        if not changed:
            break


def _normalize_types(cfg: Dict[str, Any]) -> None:
    for key in _BOOL_KEYS:
        if key in cfg:
            cfg[key] = _to_bool(cfg[key], cfg[key])

    if ("USE_DOCK6" in cfg) or ("use_dock6" in cfg):
        dock6_flag = _to_bool(cfg.get("USE_DOCK6", cfg.get("use_dock6")), False)
        cfg["USE_DOCK6"] = dock6_flag
        cfg["use_dock6"] = dock6_flag

    for key in _INT_KEYS:
        if key in cfg:
            cfg[key] = _to_int(cfg[key], cfg[key])

    for key in _FLOAT_KEYS:
        if key in cfg:
            cfg[key] = _to_float(cfg[key], cfg[key])

    cfg["DOCKING_MODE"] = str(cfg.get("DOCKING_MODE", "discovery")).lower()


def _propagate_compat_aliases(cfg: Dict[str, Any]) -> None:
    for canonical, aliases in _ALIASES.items():
        if canonical == "CPU":
            # Keep MAX_PARALLEL_JOBS as input-only compatibility.
            continue
        value = cfg.get(canonical)
        if _is_blank(value):
            continue
        for alias in aliases:
            cfg[alias] = value


def normalize_config(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """
    Return a normalized config dict with canonical keys and deterministic derivations.

    This function is idempotent.
    """
    normalized: Dict[str, Any] = dict(cfg or {})
    explicit = _explicit_keys(normalized)
    _adopt_aliases(normalized)
    _repair_missing_overall_dir(normalized, explicit)
    _derive_from_overall(normalized)
    _derive_mgltools(normalized, explicit)
    _derive_reduce(normalized)
    _derive_dock6(normalized)
    _derive_ambertools(normalized, explicit)
    _derive_mmgbsa_protocol(normalized, explicit)
    _repair_missing_tool_paths(normalized)

    if _is_blank(normalized.get("PYMOL_EXE")):
        normalized["PYMOL_EXE"] = shutil.which("pymol")

    _expand_all(normalized)
    _repair_missing_tool_paths(normalized)
    _normalize_types(normalized)
    _propagate_compat_aliases(normalized)
    return normalized
