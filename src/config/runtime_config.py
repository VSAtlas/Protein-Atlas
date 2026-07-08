from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Any, Dict

import sitecustomize  # noqa: F401  # ensure HOME is writable for micromamba/pytest sandboxes
from . import normalize as _cfg_norm
from .normalize import (
    default_base_dir as _default_base_dir,
    default_dirs as _default_dirs,
    default_runtime as _default_runtime,
    default_tools as _default_tools,
    normalize_config,
)

_to_bool = _cfg_norm._to_bool
_to_int = _cfg_norm._to_int
_to_float = _cfg_norm._to_float
which_or_exists = _cfg_norm.which_or_exists

# Config loading & validation
# -------------------------
_ALLOWED_ENV_OVERRIDES = {
    "VINA_EXE",
    "VINA_PATH",
    "GNINA_EXE",
    "OPENBABEL_PATH",
    "MGLTOOLS_PATH",
    "MGLTOOLS_DIR",
    "MGLTOOLS_PYTHON",
    "PREPARE_LIGAND_SCRIPT",
    "PREPARE_RECEPTOR_SCRIPT",
    "PREPARE_LIGAND4",
    "PYMOL_EXE",
    "PYMOL_PATH",
    "P2RANK_PATH",
    "PHENIX_DIR",
    "PHENIX_LIB_PATH",
    "PHENIX_CLEAN_SCRIPT",
    "INPUT_DIR",
    "OUTPUT_DIR",
    "PDBQT_DIR",
    "DOCKED_DIR",
    "POST_DOCKED_DIR",
    "DATA_DIR",
    "LOGS_DIR",
    "MANIFESTS_DIR",
    "LIGAND_DIR",
    "EXTRACTED_LIGANDS_DIR",
    "LIGAND_EXTRACTED_DIR",
    "LIGANDS_MOL2_DIR",
    "PREPPED_LIGANDS_DIR",
    "OUTPUT_LIGANDS_DIR",
    "PREPPED_LIGANDS_ROOT",
    "REDUCE_LOCAL_CANDIDATE",
    "REDUCE_HET_DICT",
    "DOCK6_EXE",
    "DMS_EXE",
    "SPHGEN_EXE",
    "SPHERE_SELECTOR_EXE",
    "SHOWBOX_EXE",
    "GRID_EXE",
    "DOCK6_VDW_DEFN_FILE",
    "P2RANK_OUTPUT_DIR",
    "CONFIGS_DIR",
    "CPU",
    "CPU_ONLY",
    "DOCKING_MODE",
    "TEST_MODE_ENABLE",
    "TEST_FDA_LIBRARY_SUBDIR",
    "REDUCE_EXE",
    "USE_MEEKO",
    "OVERALL_DIR",
    "PROTEIN_DIR",
    # ---  allow ENV override for the knobs ---
    "BENCH_ENFORCE_HARDCODED_CONTROLS_ONLY",
    "BENCH_ALLOW_WHITELIST_FALLBACK_IF_CONTROLS_MISSING",
    "CONTROL_CENTER_POLICY",
    "CONTROL_CENTER_CLOSE_MAX_A",
    "CTRL_REDOCK_EXHAUSTIVENESS",
    "CTRL_REDOCK_NMODES",
    "CTRL_REDOCK_MAX_CONTROLS",
    "THREADS_PER_VINA_CTRL",
    # logging/topic gates (opt-in; safe to ignore if unset)
    "LOG_TOPICS",
    "LOG_LEVEL_FILE",
    "LOG_LEVEL_CONSOLE",
    # GNINA follow-up toggle
    "USE_GNINA",
    # LeDock follow-up toggle
    "USE_LEDOCK",
    # DOCK6 follow-up toggle
    "USE_DOCK6",
    # SCORCH rescoring knobs
    "USE_SCORCH",
    "SCORCH",
    "SCORCH_SCRIPT",
    "SCORCH_ENV_PREFIX",
    "SCORCH_ENV",
    "SCORCH_TOP_FRACTION",
    "SCORCH_PARALLEL_PROFILE",
    "SCORCH_PROVISIONAL_ENABLE",
    "SCORCH_PROVISIONAL_REUSE",
    "SCORCH_PROVISIONAL_CACHE_WRITE",
    "SCORCH_PROVISIONAL_TOP_FRACTION",
    "SCORCH_PROVISIONAL_MIN_PROGRESS",
    "SCORCH_PROVISIONAL_RESERVED_WORKERS",
    "SCORCH_PROVISIONAL_CPU_BUDGET",
    "SCORCH_PROVISIONAL_MAX_INFLIGHT",
    "SCORCH_PROVISIONAL_MAX_PER_COMBO",
    "SCORCH_ENV_PREFIX_GPU",
    "SCORCH_ENV_PREFIX_AMD",
    "SCORCH_ENV_PREFIX_NVIDIA",
    "SCORCH_ENV_GPU",
    "SCORCH_ENV_AMD",
    "SCORCH_ENV_NVIDIA",
    "REPAIR_SCORCH_SELECTION",
    "ATLAS_REPAIR_SCORCH_SELECTION",
    # MMGBSA receptor prep knobs
    "MMGBSA_STRIP_METALS",
    "MMGBSA_WATER_POLICY",
    "MMGBSA_WATER_KEEP_RADIUS_A",
    "MMGBSA_WATER_USE_ALIASES",
    "MMGBSA_METAL_USE_ALIASES",
    "MMGBSA_TOPOLOGY_PREP_ENABLED",
    "MMGBSA_TLEAP_RUN",
    "MMGBSA_TOPOLOGY_DIRNAME",
    "MMGBSA_AMBERTOOLS_PREFIX",
    "MMGBSA_CPPTRAJ_ENABLED",
    "MMGBSA_CPPTRAJ_RUN",
    "MMGBSA_TRAJOUT_NAME",
    "MMGBSA_TRAJOUT_FORMAT",
    "MMGBSA_TRAJIN_SOURCE",
    "MMGBSA_TRAJIN_PATH",
    "MMGBSA_TRAJIN_FORMAT",
    "MMGBSA_TRAJ_STARTFRAME",
    "MMGBSA_TRAJ_ENDFRAME",
    "MMGBSA_TRAJ_INTERVAL",
    "MMGBSA_MMPBSA_ENABLED",
    "MMGBSA_MMPBSA_RUN",
    "MMGBSA_MMPBSA_STARTFRAME",
    "MMGBSA_MMPBSA_ENDFRAME",
    "MMGBSA_MMPBSA_INTERVAL",
    "MMGBSA_MMPBSA_VERBOSE",
    "MMGBSA_GB_IGB",
    "MMGBSA_GB_SALTCON",
    "MMGBSA_MMPBSA_INPUT_NAME",
    "MMGBSA_MMPBSA_LOG_NAME",
    "MMGBSA_MMPBSA_OUT_DAT",
    "MMGBSA_MMPBSA_OUT_CSV",
    "MMGBSA_MMPBSA_MPI_ENABLED",
    "MMGBSA_MMPBSA_MPI_RANKS",
    "MMGBSA_DEFAULT_TRAJ_NAME",
    "MMGBSA_ENABLED",
    "MMGBSA_AUTO_RUN",
    "MMGBSA_PROTOCOL",
    "MMGBSA_PRODUCTION_STRICT",
    "MMGBSA_SELECTED_LIGANDS",
    "MMGBSA_INPUT_PH_LABELS",
    "MMGBSA_MD_ENABLED",
    "MMGBSA_MD_ENGINE",
    "MMGBSA_OPENMM_PLATFORM",
    "MMGBSA_OPENMM_DEVICE_INDEX",
    "MMGBSA_OPENMM_DEFAULT_DEVICE_INDEX",
    "MMGBSA_OPENMM_PRECISION",
    "MMGBSA_OPENMM_CPU_THREADS",
    "MMGBSA_OPENMM_START_STAGE",
    "MMGBSA_MD_SOLVENT_MODEL",
    "MMGBSA_MD_NREPLICATES",
    "MMGBSA_MD_MIN_STEPS",
    "MMGBSA_MD_RESCUE_ENABLED",
    "MMGBSA_MD_RESCUE_STEPS",
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
    "MMGBSA_QC_ENABLED",
    "MMGBSA_QC_RUN_CPPTRAJ",
    "MMGBSA_QC_MIN_FRAMES",
    "MMGBSA_QC_MIN_BLOCKS",
    "MMGBSA_QC_BLOCK_SIZE",
    "MMGBSA_QC_MAX_SEM_KCAL",
    "MMGBSA_QC_MAX_BLOCK_RANGE_KCAL",
    "MMGBSA_PREHEAT_QC_MAX_GMAX",
    "MMGBSA_PREHEAT_QC_MAX_RMS",
    "MMGBSA_PREHEAT_QC_MAX_ABS_ENERGY",
    "MMGBSA_PREHEAT_QC_MAX_ABS_VDW",
    "MMGBSA_PREHEAT_QC_MAX_NONFINITE_MARKERS",
    "MMGBSA_CHEMISTRY_REVIEW_REQUIRED",
    "MMGBSA_CHEMISTRY_REVIEW_STATUS",
    "MMGBSA_CHEMISTRY_REVIEW_NOTES",
    "MMGBSA_CHEMISTRY_SUBGATES_REQUIRED",
    "MMGBSA_PROTONATION_REVIEW_STATUS",
    "MMGBSA_TAUTOMER_REVIEW_STATUS",
    "MMGBSA_STEREOCHEMISTRY_REVIEW_STATUS",
    "MMGBSA_NET_CHARGE_REVIEW_STATUS",
    "MMGBSA_PARAMETER_REVIEW_STATUS",
    "MMGBSA_EXPLICIT_TRAJ_PATH",
    "MMGBSA_EXPLICIT_SOLVENT_STRIP_TRAJ",
    "MMGBSA_KEEP_WATERS",
    "MMGBSA_WATER_KEEP_RADIUS",
    "MMGBSA_KEEP_METALS",
    "MMGBSA_METAL_RETAIN_TOKENS",
    "MMGBSA_WATER_RETAIN_TOKENS",
    "MMGBSA_RECEPTOR_FORCE",
    "MMGBSA_FORCE",
    "MMGBSA_STRICT",
    "MMGBSA_LIGAND_AT",
    "MMGBSA_LIGAND_CHARGE_METHOD",
    "MMGBSA_LIGAND_PRIMARY_CHARGE_METHOD",
    "MMGBSA_LIGAND_FALLBACK_CHARGE_METHOD",
    "MMGBSA_LIGAND_NOMINAL_NET_CHARGE",
    "MMGBSA_LIGAND_BCC_CHARGE_SWEEP",
    "MMGBSA_LIGAND_BCC_SWEEP_INCLUDE_PLUSMINUS2",
    "MMGBSA_LIGAND_FORCE",
    "MMGBSA_RDKit_VALIDATE",
    "MMGBSA_RDKit_RADICAL_LOWCONF_THRESHOLD",
    "MMGBSA_LIGAND_CHEMISTRY_STRICT",
    "MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE",
    "MMGBSA_INPUT_CHEMISTRY_SOURCE",
    "MMGBSA_INPUT_CHEMISTRY_NOTES",
    "MMGBSA_LIGAND_NET_CHARGE",
    "MMGBSA_LIGAND_SQM_LEVEL",
    "MMGBSA_INPUT_STAGE_DIR",
    "MMGBSA_MAX_LIGANDS",
    "MMGBSA_RERANKED_TOP_PCT",
    "MMGBSA_ACTIVE_SITE_RADIUS_FALLBACK",
    "MMGBSA_TLEAP_ENABLED",
    "MMGBSA_TLEAP_FORCE",
    "MMGBSA_GENERAL_STARTFRAME",
    "MMGBSA_GENERAL_ENDFRAME",
    "MMGBSA_GENERAL_INTERVAL",
    "MMGBSA_GENERAL_VERBOSE",
}

def _extract_brace_block(text: str, start_idx: int, open_char="{", close_char="}"):
    """Return the brace-balanced substring starting at the first open_char after start_idx."""
    i = text.find(open_char, start_idx)
    if i == -1:
        return None
    level = 0
    j = i
    while j < len(text):
        c = text[j]
        if c == open_char:
            level += 1
        elif c == close_char:
            level -= 1
            if level == 0:
                return text[i : j + 1]
        j += 1
    return None  # unbalanced


def _strip_inline_comment(value: str) -> str:
    """Strip inline comments starting with #, except when inside quotes."""
    in_single = False
    in_double = False
    escaped = False
    for idx, ch in enumerate(value):
        if escaped:
            escaped = False
            continue
        if ch == "\\":
            escaped = True
            continue
        if ch == "'" and not in_double:
            in_single = not in_single
            continue
        if ch == '"' and not in_single:
            in_double = not in_double
            continue
        if ch == "#" and not in_single and not in_double:
            return value[:idx]
    return value


def _parse_kv_config(path: Path) -> Dict[str, str]:
    """
    Lightweight key=value parser; allows comments starting with '#'.
    """
    out: Dict[str, str] = {}
    if not path.exists():
        return out
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                k, v = line.split("=", 1)
                val = _strip_inline_comment(v).strip()
                out[k.strip().upper()] = val
    return out


def _resolve_config_path(config_path: str, base_dir: Path) -> Path:
    requested = Path(config_path)
    candidates: list[Path] = []
    if requested.is_absolute():
        candidates.append(requested)
    else:
        candidates.append(requested)
        candidates.append(base_dir / requested)

    for cand in candidates:
        if cand.exists():
            return cand

    # Backward-compatible fallback for new clones that only have config.example.txt.
    if requested.name == "config.txt":
        alt_name = "config.example.txt"
        alt_candidates = [cand.with_name(alt_name) for cand in candidates]
        for cand in alt_candidates:
            if cand.exists():
                return cand

    # Keep old behavior: return requested path even if missing.
    return candidates[0]


def load_config(
    config_path: str = "config.txt", base_dir: Path | None = None
) -> Dict[str, Any]:
    """
    Load configuration with precedence:
      1) defaults (BRCF-aware)
      2) config.txt (key=value, optional)
      3) environment variables (allowed set only)
    Then normalize derived paths and variable expansions.
    If config.txt is missing, fall back to config.example.txt when present.
    """
    base = base_dir or _default_base_dir(reference_file=__file__)
    resolved_config_path = _resolve_config_path(config_path, base)

    cfg: Dict[str, Any] = {}
    cfg.update(_default_dirs(base))
    cfg.update(_default_tools(base))
    cfg.update(_default_runtime())

    # overlay config file
    file_cfg = _parse_kv_config(resolved_config_path)

    # --- preserve multi-line TEST_LIBRARY_MAP block ---
    try:
        raw = file_cfg.get("TEST_LIBRARY_MAP")
        if isinstance(raw, str):
            s = raw.strip()
            # If it looks like a dict but was truncated to one line, re-extract the full brace block
            if s.startswith("{") and not s.endswith("}"):
                txt = resolved_config_path.read_text(encoding="utf-8", errors="ignore")
                key_idx = txt.find("TEST_LIBRARY_MAP")
                if key_idx != -1:
                    block = _extract_brace_block(txt, key_idx, "{", "}")
                    if block:
                        file_cfg["TEST_LIBRARY_MAP"] = block

    except OSError:
        # Non-fatal: if anything goes wrong, keep the original single-line value
        pass

    cfg.update(file_cfg)
    # Preserve raw file-sourced config (before env overrides) for precedence checks.
    cfg["_FILE_CFG"] = dict(file_cfg)

    # overlay env vars (uppercased keys only)
    for k, v in os.environ.items():
        K = k.upper()
        if (K in cfg or K in _ALLOWED_ENV_OVERRIDES) and v:
            cfg[K] = v

    cfg = normalize_config(cfg)

    return cfg


def _selection_pct_error(config_key: str, details: str) -> ValueError:
    return ValueError(
        f"{config_key}: {details}. "
        "Use comma-separated values as either fractions (e.g. 1.0,0.15,0.10) "
        "or percents (e.g. 100,15,10 or 100%,15%,10%). "
        "Do not mix unsuffixed fraction and percent scales in one list."
    )


def _parse_selection_pct_list(raw_val: Any, config_key: str) -> list[float]:
    if isinstance(raw_val, str):
        text = raw_val.strip()
        if not text:
            return []
        tokens = [tok.strip() for tok in text.split(",")]
    elif isinstance(raw_val, (list, tuple)):
        if not raw_val:
            return []
        tokens = [str(tok).strip() for tok in raw_val]
    else:
        raise _selection_pct_error(
            config_key,
            f"expected comma-separated string/list, got {type(raw_val).__name__}",
        )

    if any(tok == "" for tok in tokens):
        raise _selection_pct_error(config_key, "contains empty token")

    parsed: list[float] = []
    unsuffixed_scales: set[str] = set()
    saw_explicit_percent = False
    saw_unsuffixed_fraction = False

    for idx, token in enumerate(tokens, start=1):
        explicit_percent = token.endswith("%")
        num_text = token[:-1].strip() if explicit_percent else token
        try:
            value = float(num_text)
        except ValueError as exc:
            raise _selection_pct_error(
                config_key, f"token #{idx} '{token}' is not numeric ({exc})"
            ) from exc

        if not math.isfinite(value):
            raise _selection_pct_error(config_key, f"token #{idx} '{token}' is NaN/inf")
        if value <= 0.0:
            raise _selection_pct_error(
                config_key, f"token #{idx} '{token}' must be > 0"
            )

        if explicit_percent:
            saw_explicit_percent = True
            if value > 100.0:
                raise _selection_pct_error(
                    config_key, f"token #{idx} '{token}' percent exceeds 100"
                )
            parsed.append(value / 100.0)
            continue

        if value <= 1.0:
            unsuffixed_scales.add("fraction")
            saw_unsuffixed_fraction = True
            parsed.append(value)
        elif value <= 100.0:
            unsuffixed_scales.add("percent")
            parsed.append(value / 100.0)
        else:
            raise _selection_pct_error(
                config_key, f"token #{idx} '{token}' exceeds maximum of 100"
            )

    if len(unsuffixed_scales) > 1:
        raise _selection_pct_error(
            config_key, f"mixed unsuffixed scales are ambiguous (raw={raw_val!r})"
        )
    if saw_explicit_percent and saw_unsuffixed_fraction:
        raise _selection_pct_error(
            config_key,
            "cannot mix explicit % tokens with unsuffixed fraction tokens (e.g. 100%,0.1 is ambiguous; use 0.1% or 10%)",
        )
    return parsed


def _validate_selection_schedule_config(cfg: Dict[str, Any]) -> None:
    for key in ("DISCOVERY_SELECTION_PCTS", "POLYPHARM_SELECTION_PCTS"):
        raw = cfg.get(key)
        if raw in (None, "", [], ()):
            continue
        parsed = _parse_selection_pct_list(raw, key)
        if not parsed:
            raise _selection_pct_error(key, "no values provided")
        if not math.isclose(parsed[0], 1.0, rel_tol=0.0, abs_tol=1e-12):
            raise _selection_pct_error(
                key,
                f"first value must resolve to 1.0 (100%) because stage1 always runs on the full ligand pool (got {parsed[0]:.12g})",
            )


def validate_config(cfg: Dict[str, Any]):
    """
    Require OVERALL_DIR and INPUT_DIR to exist.
    Auto-create typical output directories if missing.
    """
    required_keys = [
        "OUTPUT_DIR",
        "INPUT_DIR",
        "PDBQT_DIR",
        "DOCKED_DIR",
        "OVERALL_DIR",
        "VINA_EXE",
        "CPU",
        "PYMOL_EXE",
    ]
    missing = [k for k in required_keys if not cfg.get(k)]
    if missing:
        raise ValueError(
            "Missing required config keys or values: "
            f"{missing}. Check config.txt (or config.example.txt fallback) "
            "and ensure these keys are set. Run `atlas --doctor` to inspect runtime config "
            "and `atlas --verify-tools` to validate tool paths."
        )

    # Must exist
    for path_key in ["OVERALL_DIR", "INPUT_DIR"]:
        p = Path(cfg[path_key])
        if not p.exists():
            raise FileNotFoundError(
                f"Required path check failed for {path_key}: {p}. "
                f"Create/fix this directory (for INPUT_DIR this is typically your input_pdbs path), "
                "or correct it in config.txt. Run `atlas --doctor` to inspect resolved paths."
            )

    # Create-if-missing for outputs/caches
    create_keys = [
        "OUTPUT_DIR",
        "PDBQT_DIR",
        "DOCKED_DIR",
        "PREPPED_LIGANDS_DIR",
        "EXTRACTED_LIGANDS_DIR",
        "LIGANDS_MOL2_DIR",
        "P2RANK_OUTPUT_DIR",
        "CONFIGS_DIR",
    ]
    for path_key in create_keys:
        p = Path(cfg[path_key])
        try:
            p.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            raise FileNotFoundError(
                f"Could not create directory for {path_key}: {p} ({e})"
            )

    _validate_selection_schedule_config(cfg)




def load_inputs():
    """Legacy name used by existing scripts. Now just calls load_config()."""
    return load_config()


def get_default_config(prompt: bool = False):
    """Compatibility helper that returns the merged runtime configuration."""
    return load_config()
