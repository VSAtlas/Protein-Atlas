from __future__ import annotations

import os
import re
import sys
from typing import Any, Iterable, Mapping

from cli.cli_utils import _cli_has, _cli_val
from cli.run_context import ConfigDict
from docking.library_mode import _coerce_test_map, test_libraries_value_with_token

BENCH_PDB_IDS = ["bNJS", "bOJG", "bNNQ"]
BENCH_TEST_LIBRARY_MAP = {
    "bNJS": "bench_pur2",
    "bOJG": "bench_mk01",
    "bNNQ": "bench_fabp4",
}
BENCH_SMALL_LIBRARY = "test_library_20"


def _env_int(name: str, default: int, *, min_value: int = 1) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        return max(min_value, int(str(raw).strip()))
    except Exception:
        return default


def _env_bool(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    token = str(raw).strip().lower()
    if token in {"1", "true", "t", "yes", "y", "on"}:
        return True
    if token in {"0", "false", "f", "no", "n", "off"}:
        return False
    return default


def _throughput_bench_seed() -> int | None:
    if not _env_bool("ATLAS_THROUGHPUT_BENCH", False):
        return None
    return _env_int("ATLAS_THROUGHPUT_BENCH_SEED", 1337)


def _env_bench_small_pdb_ids() -> list[str]:
    raw = os.environ.get("ATLAS_BENCH_SMALL_PDBS", "")
    if not raw.strip():
        return list(BENCH_PDB_IDS)
    allowed = {pid.upper(): pid for pid in BENCH_PDB_IDS}
    selected: list[str] = []
    seen: set[str] = set()
    for token in re.split(r"[,\s]+", raw):
        key = token.strip().upper()
        if not key or key in seen or key not in allowed:
            continue
        selected.append(allowed[key])
        seen.add(key)
    return selected or list(BENCH_PDB_IDS)


BENCH_SMALL_PDB_IDS = _env_bench_small_pdb_ids()
BENCH_VINA_SEED = _throughput_bench_seed()
BENCH_SMALL_TARGET_LIGANDS = _env_int("ATLAS_BENCH_SMALL_TARGET_LIGANDS", 96)
BENCH_SMALL_TEST_LIBRARY_MAP = {pid: BENCH_SMALL_LIBRARY for pid in BENCH_SMALL_PDB_IDS}
WATER_BENCH_PDB_IDS = ["3EML", "1UYG", "1XL2"]
WATER_BENCH_TEST_LIBRARY_MAP = {
    "3EML": "water_aa2ar",
    "1UYG": "water_hs90a",
    "1XL2": "water_hivpr",
}
WATER_BENCH_TARGET_LIGANDS = int(os.environ.get("ATLAS_WATER_BENCH_MAX_LIGANDS", "999"))
_DUD_LIBRARY_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def _set_throughput_benchmark_seed(cfg: ConfigDict) -> None:
    if BENCH_VINA_SEED is None:
        return
    cfg["_BENCH_VINA_SEED"] = int(BENCH_VINA_SEED)
    cfg["_THROUGHPUT_BENCH_SEEDED"] = True


def _bench_enabled(argv: list[str]) -> bool:
    """Check if -bench or --bench is present in arguments."""
    return _cli_has(argv, "-bench") or _cli_has(argv, "--bench")


def _bench2_enabled(argv: list[str]) -> bool:
    """Check if -bench2 or --bench2 is present in arguments."""
    return _cli_has(argv, "-bench2") or _cli_has(argv, "--bench2")


def _bench_small_enabled(argv: list[str]) -> bool:
    """Check if -bench-small or --bench-small is present in arguments."""
    return _cli_has(argv, "-bench-small") or _cli_has(argv, "--bench-small")


def _bench_micro_enabled(argv: list[str]) -> bool:
    """Check if -bench-micro or --bench-micro is present in arguments."""
    return _cli_has(argv, "-bench-micro") or _cli_has(argv, "--bench-micro")


def _water_bench_enabled(argv: list[str]) -> bool:
    """Check if a water benchmark flag is present in arguments."""
    return (
        _cli_has(argv, "-water-bench")
        or _cli_has(argv, "--water-bench")
        or _cli_has(argv, "-waterbench")
        or _cli_has(argv, "--waterbench")
    )


def _dude_enabled(argv: list[str]) -> bool:
    """Check if -dude or --dude is present in arguments."""
    return _cli_has(argv, "-dude") or _cli_has(argv, "--dude")


def _dud_library_value(argv: list[str]) -> str:
    raw = _cli_val(argv, "--dud-library")
    if raw is None:
        if _cli_has(argv, "--dud-library"):
            raise ValueError("--dud-library requires a value")
        return ""
    value = str(raw).strip()
    if not value:
        raise ValueError("--dud-library requires a non-empty library subdirectory")
    if not _DUD_LIBRARY_RE.match(value):
        raise ValueError(
            "--dud-library expects a prepped_ligands subdirectory name "
            "using letters, numbers, '.', '_', or '-'"
        )
    return value


def _apply_dud_library_override(
    cfg: ConfigDict, pdb_ids: Iterable[str], dud_library: str
) -> int:
    """Use one configured DUD/decoy library for every selected PDB in this run."""
    library = str(dud_library).strip()
    if not library:
        return 0

    current_map = _coerce_test_map(cfg.get("TEST_LIBRARY_MAP", {}))
    updated = dict(current_map)
    count = 0
    for raw_pid in pdb_ids:
        pid = str(raw_pid or "").strip().upper()
        if not pid:
            continue
        updated[pid] = library
        count += 1

    if count:
        cfg["TEST_LIBRARY_MAP"] = updated
        effective_mode = test_libraries_value_with_token(
            cfg.get("TEST_MODE_ENABLE", "off"), "dud"
        )
        cfg["TEST_MODE_ENABLE"] = effective_mode
        cfg["_TEST_MODE_ENABLE_OVERRIDE"] = effective_mode
        cfg["_DUD_LIBRARY_OVERRIDE"] = library
    return count


def _apply_bench_overrides(cfg: ConfigDict) -> None:
    """Force benchmark configuration settings."""
    cfg["TEST_MODE_ENABLE"] = "dud"
    cfg["_TEST_MODE_ENABLE_OVERRIDE"] = "dud"
    _set_throughput_benchmark_seed(cfg)

    cfg["USE_GNINA"] = True
    cfg["USE_LEDOCK"] = True
    cfg["USE_DOCK6"] = True
    cfg["USE_SCORCH"] = True

    cfg["APO_HOLO_MODE"] = "holo"
    cfg["PH_ENSEMBLE"] = True

    cfg["SPECIFIED_PROTEINS"] = ",".join(BENCH_PDB_IDS)

    current_map = _coerce_test_map(cfg.get("TEST_LIBRARY_MAP", {}))
    if hasattr(current_map, "copy"):
        new_map = dict(current_map)
    else:
        new_map = {}

    for k, v in BENCH_TEST_LIBRARY_MAP.items():
        new_map[k] = v

    cfg["TEST_LIBRARY_MAP"] = new_map
    # Intentionally do not override stage-selection controls for -bench/-bench2.
    # These runs should honor configured discovery/polypharm percentages.


def _apply_bench2_overrides(cfg: ConfigDict) -> None:
    """Force benchmark2 configuration settings (Vina + SCORCH only)."""
    _apply_bench_overrides(cfg)
    cfg["USE_GNINA"] = False
    cfg["USE_LEDOCK"] = False
    cfg["USE_DOCK6"] = False
    cfg["USE_SCORCH"] = True


def _apply_bench_small_overrides(cfg: ConfigDict) -> None:
    """Force small-throughput benchmark profile (separate from bench2)."""
    pdb_ids = list(BENCH_SMALL_PDB_IDS)
    cfg["TEST_MODE_ENABLE"] = BENCH_SMALL_LIBRARY
    cfg["_TEST_MODE_ENABLE_OVERRIDE"] = BENCH_SMALL_LIBRARY
    _set_throughput_benchmark_seed(cfg)
    cfg["USE_GNINA"] = False
    cfg["USE_LEDOCK"] = False
    cfg["USE_DOCK6"] = False
    cfg["USE_SCORCH"] = True
    cfg["APO_HOLO_MODE"] = "holo"
    cfg["PH_ENSEMBLE"] = True
    cfg["SPECIFIED_PROTEINS"] = ",".join(pdb_ids)
    cfg["TEST_LIBRARY_MAP"] = {pid: BENCH_SMALL_LIBRARY for pid in pdb_ids}
    cfg["ENABLE_DISTRIBUTED_COMBO_CHUNKS"] = True
    cfg["_BENCH_SMALL_TARGET_LIGANDS"] = int(BENCH_SMALL_TARGET_LIGANDS)
    cfg["_BENCH_SMALL_SEED"] = 1337
    # Bench parity mode: keep full ligand pass-through across stages and rescoring.
    cfg["DISCOVERY_SELECTION_PCTS"] = "1,1,1,1,1"
    cfg["POLYPHARM_SELECTION_PCTS"] = "1,1,1"
    cfg["SCORCH_TOP_FRACTION"] = 1.0
    cfg["CNN_TOP_FRACTION"] = 1.0
    cfg["PERCENT_SEARCHED"] = 1.0


def _apply_bench_micro_overrides(cfg: ConfigDict) -> None:
    """Force a tiny throughput profile for rapid scheduler iteration."""
    _apply_bench_small_overrides(cfg)
    cfg["SPECIFIED_PROTEINS"] = os.environ.get("ATLAS_BENCH_MICRO_PDBS", "bOJG")
    cfg["TEST_LIBRARY_MAP"] = {
        pid: BENCH_SMALL_LIBRARY
        for pid in [
            token.strip()
            for token in re.split(r"[,\s]+", str(cfg["SPECIFIED_PROTEINS"]))
            if token.strip()
        ]
    }
    cfg["_BENCH_SMALL_TARGET_LIGANDS"] = _env_int(
        "ATLAS_BENCH_MICRO_TARGET_LIGANDS", 6
    )
    cfg["_BENCH_CACHED_PREP_REQUIRED"] = _env_bool(
        "ATLAS_BENCH_CACHED_PREP_REQUIRED", True
    )
    if bool(cfg["_BENCH_CACHED_PREP_REQUIRED"]):
        cfg["FORCE_REPROCESS"] = False


def _apply_water_bench_overrides(cfg: ConfigDict) -> None:
    """Force a compact DUD-E water-policy enrichment benchmark profile."""
    cfg["TEST_MODE_ENABLE"] = "dud"
    cfg["_TEST_MODE_ENABLE_OVERRIDE"] = "dud"
    _set_throughput_benchmark_seed(cfg)
    cfg["USE_GNINA"] = False
    cfg["USE_LEDOCK"] = False
    cfg["USE_DOCK6"] = False
    cfg["USE_SCORCH"] = True
    cfg["APO_HOLO_MODE"] = "holo"
    cfg["PH_ENSEMBLE"] = True
    cfg["SPECIFIED_PROTEINS"] = ",".join(WATER_BENCH_PDB_IDS)
    cfg["TEST_LIBRARY_MAP"] = dict(WATER_BENCH_TEST_LIBRARY_MAP)
    cfg["WATER_POLICY"] = os.environ.get("ATLAS_WATER_BENCH_POLICY", "site_only")
    cfg["WATER_BENCH_STAGE1_ONLY"] = True
    cfg["ENABLE_DISTRIBUTED_COMBO_CHUNKS"] = True
    cfg["_BENCH_SMALL_TARGET_LIGANDS"] = int(WATER_BENCH_TARGET_LIGANDS)
    cfg["_WATER_BENCH_TARGET_LIGANDS"] = int(WATER_BENCH_TARGET_LIGANDS)
    cfg["DISCOVERY_SELECTION_PCTS"] = "1"
    cfg["POLYPHARM_SELECTION_PCTS"] = "1"
    cfg["SCORCH_TOP_FRACTION"] = 1.0
    cfg["CNN_TOP_FRACTION"] = 1.0
    cfg["PERCENT_SEARCHED"] = 1.0


def apply_benchmark_vina_seeds(
    cfg: Mapping[str, Any], stages: Iterable[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    """Return stage configs with fixed Vina seeds for benchmark profiles only."""
    if not bool(cfg.get("_THROUGHPUT_BENCH_SEEDED")) or "_BENCH_VINA_SEED" not in cfg:
        return [dict(stage) for stage in stages]
    try:
        base_seed = int(cfg.get("_BENCH_VINA_SEED", 1337))
    except Exception:
        base_seed = 1337
    seeded: list[dict[str, Any]] = []
    for idx, stage in enumerate(stages):
        item = dict(stage)
        item["seed"] = int(base_seed + idx * 1009)
        seeded.append(item)
    return seeded


def _apply_dude_overrides(cfg: ConfigDict) -> None:
    """Force DUD-only runtime mode without touching on-disk config."""
    cfg["TEST_MODE_ENABLE"] = "dud"
    cfg["_TEST_MODE_ENABLE_OVERRIDE"] = "dud"


def _normalize_artifact_retention_mode(raw: Any) -> str:
    token = str(raw or "").strip().lower().replace("-", "_")
    if token in {"", "rerunsafe", "rerun_safe"}:
        return "rerun_safe"
    if token in {"minimaldisk", "minimal_disk"}:
        return "minimal_disk"
    if token in {"off", "none", "false", "0", "disabled", "disable"}:
        return "off"
    return "rerun_safe"


def _resolve_artifact_retention_mode(
    cfg: Mapping[str, Any], argv: list[str]
) -> tuple[str, str]:
    cli_retain = _cli_has(argv, "--retain")
    cli_noretain = _cli_has(argv, "--noretain")
    cli_has_mode = _cli_has(argv, "--retainmode")
    cli_mode_raw = _cli_val(argv, "--retainmode")
    if cli_has_mode and not cli_mode_raw:
        print("ERROR: --retainmode requires a value", file=sys.stderr)
        sys.exit(2)

    if cli_noretain and (cli_retain or cli_mode_raw):
        print(
            "ERROR: --noretain cannot be combined with --retain or --retainmode",
            file=sys.stderr,
        )
        sys.exit(2)

    if cli_mode_raw:
        mode = _normalize_artifact_retention_mode(cli_mode_raw)
        if mode == "off":
            print(
                "ERROR: --retainmode must be rerun_safe or minimal_disk",
                file=sys.stderr,
            )
            sys.exit(2)
        return mode, "CLI(--retainmode)"

    if cli_noretain:
        return "off", "CLI(--noretain)"

    if cli_retain:
        return "rerun_safe", "CLI(--retain)"

    env_raw = os.environ.get("ARTIFACT_RETENTION")
    if env_raw is not None and str(env_raw).strip():
        return _normalize_artifact_retention_mode(env_raw), "ENV(ARTIFACT_RETENTION)"

    cfg_raw = cfg.get("ARTIFACT_RETENTION")
    if cfg_raw is not None and str(cfg_raw).strip():
        return _normalize_artifact_retention_mode(cfg_raw), "CFG(ARTIFACT_RETENTION)"

    return "rerun_safe", "DEFAULT"
