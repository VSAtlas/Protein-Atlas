from __future__ import annotations

import logging
import os
import shlex
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

from docking.apo_holo_runtime import _debug_normalize_mode_token, resolve_apo_holo_mode
from cli import run_profiles as _run_profiles
from cli.cli_utils import (
    _cli_has,
    _cli_val,
    _norm_pdb_id,
    _parse_fast_flag,
    _parse_single_from_cli,
    _parse_specified_proteins,
)
from cli.distributed_context import (
    remove_own_marker_if_present,
    resolve_distributed_context,
    shard_pdb_files_for_context,
)
from cli.distributed_mode_policy import resolve_execution_mode
from cli.run_context import ConfigDict
from cli.run_bootstrap_helpers import define_docking_stages, init_config_run_dir
from config.normalize import _to_bool
from config.output_paths import output_root, run_scoped_root
from cli.runtime_logging import ensure_file_handler
from path_router.path_router import run_logs_dir
from protein_prep.pdb_fixer_runtime import get_atom_rules
from cli.postrun_hooks_runtime import _start_scorch_mixed_queue_service
from cli.run_manifest_runtime import (
    init_run_manifest,
    load_run_manifest,
    update_manifest_for_config_hash,
    update_manifest_for_run_config,
    update_manifest_for_scheduled_proteins,
)

from cli.distributed_chunk_planner import (
    _chunk_planner_manifest_mode,
    _chunk_planner_manifest_sample_size,
    _chunk_planner_scan_workers,
    _log_distributed_cpu_node,
    _manifest_prebuild_enabled,
    _normalize_ph_tag_token,
    _resolve_runtime_cpu_settings,
)
from docking.library_mode import _coerce_test_map, parse_test_libraries

_DEFAULT_TEST_FDA_LIBRARY_SUBDIR = "fda_test_library_10"


def _clear_test_library_resolution_caches(cfg: ConfigDict) -> None:
    for key in (
        "_TEST_LIBRARY_MAP_CACHE_KEY",
        "_TEST_LIBRARY_MAP_CACHE_VAL",
        "_ALLOWED_LIBRARY_ROOTS_CACHE",
        "_LIB_ROOTS_MAP_LOG_FINGERPRINT",
    ):
        cfg.pop(key, None)


def _retarget_default_test_library_map_for_test_fda(
    cfg: ConfigDict,
    test_fda_library: str,
) -> bool:
    """Retarget only the bundled TEST placeholder map when -test-fda changes library."""
    test_fda_library = str(test_fda_library or "").strip()
    if not test_fda_library or test_fda_library == _DEFAULT_TEST_FDA_LIBRARY_SUBDIR:
        return False
    test_map = _coerce_test_map(cfg.get("TEST_LIBRARY_MAP", {}))
    if test_map.get("TEST") != _DEFAULT_TEST_FDA_LIBRARY_SUBDIR:
        return False
    updated = dict(test_map)
    updated["TEST"] = test_fda_library
    cfg["TEST_LIBRARY_MAP"] = updated
    _clear_test_library_resolution_caches(cfg)
    logging.getLogger("run").info(
        "[config] TEST_LIBRARY_MAP.TEST retargeted from %s to %s for test FDA mode",
        _DEFAULT_TEST_FDA_LIBRARY_SUBDIR,
        test_fda_library,
    )
    return True


def select_pdb_files_for_run(
    cfg: ConfigDict,
    argv_for_parsing: list[str],
    *,
    is_resume: bool = False,
    resume_protein_ids: list[str] | None = None,
) -> list[str]:
    """Resolve the input PDB filenames for a run and apply test-library policy."""
    cfg.setdefault("SPECIFIED_PROTEINS", "")
    requested_ids, _ = _parse_specified_proteins(argv_for_parsing, cfg)
    cfg["_EFFECTIVE_SPECIFIED_PROTEINS"] = requested_ids
    print(
        f"[config] SPECIFIED_PROTEINS effective={requested_ids} (precedence: CLI>ENV>CFG)"
    )

    pdb_files = [
        f
        for f in os.listdir(cfg["INPUT_DIR"])
        if f.lower().endswith(".pdb") and "_nolig" not in f.lower()
    ]

    id_index: dict[str, str] = {}
    for f in pdb_files:
        base = os.path.splitext(f)[0].replace("_cleaned", "")
        nid = _norm_pdb_id(base)
        if nid:
            id_index.setdefault(nid, f)

    req = list(cfg.get("_EFFECTIVE_SPECIFIED_PROTEINS", []) or [])
    if is_resume and not req and resume_protein_ids:
        req = list(resume_protein_ids)
    if req:
        normalized_req: list[str] = []
        seen_req: set[str] = set()
        for raw in req:
            nid = _norm_pdb_id(str(raw))
            if not nid or nid in seen_req:
                continue
            normalized_req.append(nid)
            seen_req.add(nid)
        req = normalized_req
    if req:
        hits = [nid for nid in req if nid in id_index]
        miss = [nid for nid in req if nid not in id_index]

        print(
            f"[filter.proteins] mode=on requested={len(req)} present={len(hits)} missing={len(miss)} ? {hits}"
        )
        for m in miss:
            print(
                f"WARNING: requested PDB '{m}' not found under INPUT_DIR={cfg['INPUT_DIR']} or was excluded (_nolig)."
            )

        if not hits:
            print(
                "ERROR: No requested proteins found. Exiting with status 2 to avoid a no-op run."
            )
            sys.exit(2)

        pdb_files = [id_index[nid] for nid in hits]
        print("Selected proteins (Specified Proteins Mode): " + ", ".join(hits))
    else:
        print(
            f"[filter.proteins] mode=off requested=0 present={len(pdb_files)} missing=0 ? []"
        )

    selected_ids = sorted(
        {
            str(nid).upper()
            for nid in (_norm_pdb_id(f) for f in pdb_files)
            if nid is not None
        }
    )
    dud_library_override = str(cfg.get("_DUD_LIBRARY_OVERRIDE") or "").strip()
    if dud_library_override:
        mapped_count = _run_profiles._apply_dud_library_override(
            cfg, selected_ids, dud_library_override
        )
        print(
            "[config] DUD_LIBRARY_OVERRIDE "
            f"library={dud_library_override} mapped_pdbs={mapped_count}"
        )
        logging.info(
            "[dud-library] override library=%s mapped_pdbs=%d",
            dud_library_override,
            mapped_count,
        )

    tokens = parse_test_libraries(cfg)
    raw_map = cfg.get("TEST_LIBRARY_MAP", {})
    test_map = _coerce_test_map(raw_map)
    try:
        cfg["_TEST_LIBRARY_CANONICAL"] = {
            str(k).upper(): str(v) for k, v in getattr(test_map, "items", lambda: [])()
        }
    except Exception:
        cfg["_TEST_LIBRARY_CANONICAL"] = {}

    if "dud" in tokens:
        test_keys = set()
        for k in getattr(test_map, "keys", lambda: [])():
            nid = _norm_pdb_id(str(k))
            if nid:
                test_keys.add(nid)

        token_set = {
            str(tok).strip().lower()
            for tok in tokens
            if str(tok).strip()
        }
        strict_dud_only = token_set == {"dud"}
        missing_map_ids = [pid for pid in selected_ids if pid not in test_keys]

        if strict_dud_only:
            if not test_keys:
                print(
                    "ERROR: TEST_MODE_ENABLE includes 'dud' but TEST_LIBRARY_MAP is empty/invalid. "
                    "Populate TEST_LIBRARY_MAP for every selected PDB before running distributed DUD mode.",
                    file=sys.stderr,
                )
                sys.exit(2)

            if missing_map_ids:
                print(
                    "ERROR: DUD mode requires TEST_LIBRARY_MAP coverage for all selected PDBs. "
                    f"Missing mappings for: {', '.join(missing_map_ids)}",
                    file=sys.stderr,
                )
                sys.exit(2)

            kept = [f for f in pdb_files if (_norm_pdb_id(f) or "") in test_keys]
            if len(kept) != len(pdb_files):
                print(
                    f"[test-mode] Enabled tokens={'+'.join(tokens)}; restricting to "
                    f"{len(kept)} PDBs from TEST_LIBRARY_MAP keys."
                )
            pdb_files = kept
        elif missing_map_ids:
            print(
                "[test-mode] WARNING: mixed mode includes 'dud' but some PDBs are not "
                "mapped in TEST_LIBRARY_MAP. Unmapped proteins will run non-DUD "
                f"libraries only: {', '.join(missing_map_ids)}",
                file=sys.stderr,
            )

    return pdb_files


def _resolve_test_mode(cfg: Mapping[str, Any]) -> str:
    raw = (
        os.environ.get("TEST_MODE_ENABLE")
        if "TEST_MODE_ENABLE" in os.environ
        else cfg.get("TEST_MODE_ENABLE", "off")
    )
    if isinstance(raw, bool) and not raw:
        return "off"
    token = str(raw).strip()
    if not token or token.lower() in {"0", "false", "no", "off", "none", "null"}:
        return "off"

    tokens = parse_test_libraries(cfg)
    reserved = {"dud", "fda", "hmdb"}
    if tokens and all(tok in reserved for tok in tokens):
        token_set = set(tokens)
        if token_set == {"dud"}:
            return "dud"
        if token_set == {"fda"}:
            return "fda"
        if token_set == {"hmdb"}:
            return "hmdb"
        if token_set == {"dud", "fda"}:
            return "fda+dud"
        if token_set == {"hmdb", "dud"}:
            return "hmdb+dud"
        if token_set == {"hmdb", "fda"}:
            return "hmdb+fda"
        if token_set == {"dud", "fda", "hmdb"}:
            return "fda+dud+hmdb"
    return "+".join(tokens) if tokens else "off"


@dataclass
class BootstrapState:
    cfg: ConfigDict
    argv_for_parsing: list[str]
    dist_ctx: Any
    completed_combo_lookup: set[tuple[str, str, str]]
    params: Any
    stages: Any
    pdb_files: list[str]
    tokens: list[str]
    test_mode: str
    dist_combo_chunk_mode: bool
    dist_combo_hybrid_mode: bool
    execution_mode_decision: Any
    plan_only: bool
    mode: str
    variants: list[Any]
    start: float
    global_start: float
    failed_root: str
    failed_entries: list[tuple[str, str, str, str, str]]
    run_scope_completion_times: dict[tuple[str, str, str], float]
    run_chunk_rebalance_count: int
    recorded_terminal_chunk_ids: set[str]
    distributed_assigned_pdb_ids: set[str]
    retention_lock: threading.Lock
    pending_retention_combos: set[tuple[str, str, str]]
    pending_coverage_refresh: dict[tuple[str, str, str], dict[str, Any]]
    pending_retention_lock: threading.Lock
    scorch_queue_service: Any


@dataclass
class WorkloadSelectionState:
    params: Any
    stages: Any
    pdb_files: list[str]
    all_pdb_files: list[str]
    tokens: list[str]
    test_mode: str
    dist_combo_chunk_mode: bool
    dist_combo_hybrid_mode: bool
    execution_mode_decision: Any


@dataclass
class VariantBootstrapState:
    start: float
    plan_only: bool
    mode: str
    variants: list[Any]
    failed_root: str


def _write_profile_snapshot(
    cfg: Mapping[str, Any],
    *,
    filename: str,
    header: str,
    use_gnina: bool,
    use_ledock: bool,
    use_dock6: bool,
    use_scorch: bool,
    test_library_map: Mapping[str, str],
    extra_lines: Sequence[str] = (),
) -> None:
    if "CONFIG_RUN_DIR" not in cfg:
        return
    try:
        snap_path = Path(str(cfg["CONFIG_RUN_DIR"])) / str(filename)
        lines = [
            str(header),
            f"TEST_MODE_ENABLE={cfg.get('TEST_MODE_ENABLE')}",
            f"USE_GNINA={str(bool(use_gnina)).lower()}",
            f"USE_LEDOCK={str(bool(use_ledock)).lower()}",
            f"USE_DOCK6={str(bool(use_dock6)).lower()}",
            f"USE_SCORCH={str(bool(use_scorch)).lower()}",
            f"SCORCH_TOP_FRACTION={cfg.get('SCORCH_TOP_FRACTION')}",
            "APO_HOLO_MODE=holo",
            "PH_ENSEMBLE=true",
            f"BENCH_PDBS={cfg.get('SPECIFIED_PROTEINS')}",
        ]
        lines.extend(str(x) for x in extra_lines)
        for k, v in test_library_map.items():
            lines.append(f"TEST_LIBRARY_MAP.{k}={v}")
        snap_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        logging.info("[bench] wrote snapshot config to %s", snap_path)
    except Exception:
        logging.warning("[bench] failed to write snapshot config", exc_info=True)


def _setup_distributed_runtime(cfg: ConfigDict, *, run_id: str) -> Any:
    cfg_cpu_raw = 0
    alloc_cpu = 0
    effective_cpu = 0
    effective_scheduler_cpus = 0
    dist_ctx = resolve_distributed_context(cfg, run_id)
    cfg_cpu_raw, alloc_cpu, effective_cpu, effective_scheduler_cpus = (
        _resolve_runtime_cpu_settings(cfg, distributed_enabled=dist_ctx.enabled)
    )
    cfg["CPU"] = effective_cpu
    cfg["GLOBAL_SCHEDULER_CPUS"] = effective_scheduler_cpus

    if dist_ctx.enabled:
        logging.info(
            "[cpu.alloc.distributed] cfg_cpu=%d alloc_cpu=%d effective_cpu=%d effective_scheduler_cpus=%d",
            cfg_cpu_raw,
            alloc_cpu,
            effective_cpu,
            effective_scheduler_cpus,
        )
    elif effective_cpu != cfg_cpu_raw:
        logging.info(
            "[cpu.alloc.cap] cfg_cpu=%d alloc_cpu=%d effective_cpu=%d",
            cfg_cpu_raw,
            alloc_cpu,
            effective_cpu,
        )

    if dist_ctx.enabled:
        logging.info(
            "[distributed.mode] enabled=true mode=%s task_id=%d task_count=%d task_min_id=%d leader_task_id=%d is_leader=%s",
            dist_ctx.mode,
            dist_ctx.task_id,
            dist_ctx.task_count,
            dist_ctx.task_min_id,
            dist_ctx.leader_task_id,
            str(dist_ctx.is_leader).lower(),
        )
        _log_distributed_cpu_node(
            run_id=run_id,
            dist_ctx=dist_ctx,
            cfg_cpu_raw=cfg_cpu_raw,
            alloc_cpu=alloc_cpu,
            effective_cpu=effective_cpu,
            effective_scheduler_cpus=effective_scheduler_cpus,
        )
        remove_own_marker_if_present(dist_ctx, cfg)
    else:
        logging.info("[distributed.mode] enabled=false")
    return dist_ctx


def _log_atom_alias_summary() -> None:
    rules = get_atom_rules()
    alias_sets = getattr(rules, "alias_sets", None)
    waters_set = set(getattr(rules, "waters", set()))
    if not waters_set and alias_sets is not None:
        waters_set = {
            str(tok).strip().upper()
            for tok in getattr(alias_sets, "waters", set())
            if str(tok).strip()
        }
    cofactors_set = {
        str(tok).strip().upper()
        for tok in getattr(rules, "cofactors", set())
        if str(tok).strip()
    }
    elements_set = {
        str(tok).strip().upper()
        for tok in getattr(rules, "elem_tokens_canonical", set())
        if str(tok).strip()
    }
    logging.info(
        "[aliases.summary] mode=%s keep={waters:%d, cofactors:%d, elements:%d}",
        getattr(rules, "policy_mode", "LEGACY"),
        len(waters_set),
        len(cofactors_set),
        len(elements_set),
    )


def _apply_ph_ensemble_settings(cfg: ConfigDict) -> None:
    cfg.setdefault("PH_ENSEMBLE", False)
    cfg.setdefault("PH_RADIUS", 1000000.0)

    env_ph_flag = os.environ.get("PH_ENSEMBLE")
    if env_ph_flag is not None:
        cfg.PH_ENSEMBLE = _to_bool(env_ph_flag)
    else:
        cfg.PH_ENSEMBLE = _to_bool(cfg.PH_ENSEMBLE)

    scope_env = os.environ.get("PH_SCOPE")
    if scope_env is not None:
        cfg.PH_SCOPE = scope_env.strip()
    elif "PH_SCOPE" in cfg:
        cfg.PH_SCOPE = cfg["PH_SCOPE"]

    radius_env = os.environ.get("PH_RADIUS")
    if radius_env is not None:
        try:
            cfg.PH_RADIUS = float(radius_env)
        except Exception:
            cfg.PH_RADIUS = 1000000.0
    else:
        try:
            cfg.PH_RADIUS = float(cfg.PH_RADIUS)
        except Exception:
            cfg.PH_RADIUS = 100000000.0
    if cfg.PH_RADIUS <= 0:
        cfg.PH_RADIUS = 10000000.0

    log = logging.getLogger("ph_ensemble")
    log.info(
        "[ph_ensemble.mode] enabled=%s scope=%s radius=%s",
        cfg.PH_ENSEMBLE,
        getattr(cfg, "PH_SCOPE", "auto"),
        getattr(cfg, "PH_RADIUS", 1000000.0),
    )


def _hydrate_resume_state(
    *,
    cfg: ConfigDict,
    run_id: str,
    is_resume: bool,
    argv_for_parsing: list[str],
    resume_protein_ids: list[str],
    completed_combo_lookup: set[tuple[str, str, str]],
) -> list[str]:
    if not is_resume:
        return argv_for_parsing

    resume_manifest = load_run_manifest(cfg, run_id)
    if resume_manifest is None:
        print(
            f"ERROR: -resume requested but run_manifest.yaml for run_id={run_id} could not be loaded.",
            file=sys.stderr,
        )
        sys.exit(2)

    command_section = resume_manifest.get("command") or {}
    stored_argv = command_section.get("argv")
    if stored_argv:
        try:
            parsed_args = shlex.split(str(stored_argv))
            if parsed_args:
                argv_for_parsing = parsed_args
                logging.info("[resume] using stored argv=%s", parsed_args)
        except Exception:
            logging.warning(
                "[resume] Failed to parse stored argv=%r",
                stored_argv,
                exc_info=True,
            )

    listed = command_section.get("pdb_list") or []
    if isinstance(listed, (list, tuple)):
        resume_protein_ids.extend(str(x).strip().upper() for x in listed if str(x).strip())

    if not resume_protein_ids:
        summary_section = resume_manifest.get("summary") or {}
        summary_list = summary_section.get("total_protein_list") or []
        if isinstance(summary_list, (list, tuple)):
            resume_protein_ids.extend(
                str(x).strip().upper() for x in summary_list if str(x).strip()
            )

    proteins = resume_manifest.get("proteins") or {}
    completed_entries = 0
    if isinstance(proteins, Mapping):
        cfg["_RESUME_MANIFEST_PROTEINS"] = {
            str(key): dict(entry)
            for key, entry in proteins.items()
            if isinstance(entry, Mapping)
        }
        for key, entry in proteins.items():
            try:
                raw_key = str(key)
                parts = raw_key.split("|", 2)
                if len(parts) != 3:
                    continue
                pdb_part, variant_part, _ph = parts
                pdb_token = (pdb_part or "").strip().upper()
                if not pdb_token:
                    continue
                resume_variant_token = (variant_part or "legacy").strip().lower() or "legacy"
                ph_token = _normalize_ph_tag_token(parts[2] if len(parts) > 2 else None)
                if not resume_protein_ids:
                    resume_protein_ids.append(pdb_token)
                status = (
                    ((entry or {}).get("status") if isinstance(entry, Mapping) else "") or ""
                ).strip().lower()
                if status == "completed":
                    completed_combo_lookup.add((pdb_token, resume_variant_token, ph_token))
                    completed_entries += 1
            except Exception:
                logging.warning(
                    "[resume.lookup.skip] key=%r entry=%r",
                    key,
                    entry,
                    exc_info=True,
                )

    if resume_protein_ids:
        seen_ids = set()
        deduped = []
        for pid in resume_protein_ids:
            token = pid.strip().upper()
            if not token or token in seen_ids:
                continue
            seen_ids.add(token)
            deduped.append(token)
        resume_protein_ids[:] = deduped
    logging.info(
        "[resume] loaded manifest for run_id=%s; completed protein entries=%d completed_combos=%d",
        run_id,
        completed_entries,
        len(completed_combo_lookup),
    )
    return argv_for_parsing


def _apply_runtime_config_defaults(
    *,
    cfg: ConfigDict,
    run_id: str,
    argv_for_parsing: list[str],
    is_resume: bool,
    dist_ctx: Any,
    is_bench: bool,
    is_bench2: bool,
    is_bench_small: bool,
    bench_test_library_map: Mapping[str, str],
    bench_small_test_library_map: Mapping[str, str],
    bench_small_target_ligands: int,
    resolve_artifact_retention_mode: Callable[[Mapping[str, Any], list[str]], tuple[str, str]],
    log_path: str,
) -> float:
    log = logging.getLogger("ph_ensemble")
    cfg.setdefault("PH_LIGAND_MODE", "off")
    cli_ph_mode = _cli_val(argv_for_parsing, "--ph-ligand-mode")
    env_ph_mode = os.environ.get("PH_LIGAND_MODE")
    cfg_ph_mode = str(cfg.get("PH_LIGAND_MODE", "off"))
    effective_ph_mode = next(
        (
            m
            for m in (cli_ph_mode, env_ph_mode, cfg_ph_mode)
            if m is not None and str(m).strip() != ""
        ),
        "off",
    )
    cfg["PH_LIGAND_MODE"] = effective_ph_mode
    log.info(
        "[ph_ligand.mode] effective=%r (source=%s)",
        effective_ph_mode,
        "CLI" if cli_ph_mode else "ENV" if env_ph_mode else "CFG",
    )

    cfg.setdefault("CONFIGS_DIR", str(output_root(Path(cfg["OVERALL_DIR"]), "configs")))
    cfg.setdefault("RESET_CONFIGS", True)
    cli_cfg_dir = _cli_val(argv_for_parsing, "--configs-dir")
    cli_no_reset = _cli_has(argv_for_parsing, "--no-reset-configs")
    if cli_cfg_dir:
        cfg["CONFIGS_DIR"] = cli_cfg_dir
    cfg["RUN_ID"] = run_id
    default_processed_root = output_root(Path(cfg["OVERALL_DIR"]), "processed_pdbs")
    cfg["OUTPUT_DIR"] = str(
        run_scoped_root(
            Path(str(cfg.get("OUTPUT_DIR") or default_processed_root)),
            run_id,
        )
    )
    if is_resume:
        cfg["RESET_CONFIGS"] = False
    elif dist_ctx.enabled:
        cfg["RESET_CONFIGS"] = False
    elif cli_no_reset:
        cfg["RESET_CONFIGS"] = False

    init_config_run_dir(
        cfg,
        run_id=cfg.get("RUN_ID"),
        reset=cfg.get("RESET_CONFIGS"),
        logger=logging.getLogger("run"),
    )
    run_log_dir = run_logs_dir(cfg)
    pipeline_log_path = run_log_dir / "pipeline.log"
    ensure_file_handler(logging.getLogger(), pipeline_log_path, cfg, level=logging.INFO)
    logging.info("[run.log] pipeline_log=%s", pipeline_log_path)

    if is_bench:
        _write_profile_snapshot(
            cfg,
            filename="bench_config.txt",
            header="# Benchmark Configuration Snapshot",
            use_gnina=True,
            use_ledock=True,
            use_dock6=True,
            use_scorch=True,
            test_library_map=bench_test_library_map,
        )
    if is_bench2:
        _write_profile_snapshot(
            cfg,
            filename="bench2_config.txt",
            header="# Benchmark2 Configuration Snapshot",
            use_gnina=False,
            use_ledock=False,
            use_dock6=False,
            use_scorch=True,
            test_library_map=bench_test_library_map,
        )
    if is_bench_small:
        _write_profile_snapshot(
            cfg,
            filename="bench_small_config.txt",
            header="# Bench Small Configuration Snapshot",
            use_gnina=False,
            use_ledock=False,
            use_dock6=False,
            use_scorch=True,
            test_library_map=bench_small_test_library_map,
            extra_lines=(
                "BENCH_SMALL_TARGET_LIGANDS="
                f"{int(cfg.get('_BENCH_SMALL_TARGET_LIGANDS', bench_small_target_ligands) or bench_small_target_ligands)}",
                f"ENABLE_DISTRIBUTED_COMBO_CHUNKS={cfg.get('ENABLE_DISTRIBUTED_COMBO_CHUNKS')}",
            ),
        )

    print(f"[cfg.run] run_id={cfg['RUN_ID']} run_dir={cfg['CONFIG_RUN_DIR']}")
    print(f"[ph.mode] PH_ENSEMBLE={cfg.get('PH_ENSEMBLE', False)}")

    global_start = time.time()
    if not is_resume and (not dist_ctx.enabled or dist_ctx.is_leader):
        try:
            init_run_manifest(cfg, run_id, sys.argv, log_path)
        except Exception:
            logging.warning("Failed to initialize run_manifest.yaml", exc_info=True)
        else:
            try:
                update_manifest_for_config_hash(cfg, run_id)
            except Exception:
                logging.warning(
                    "[run-manifest] Failed to record config hash run_id=%s",
                    run_id,
                    exc_info=True,
                )
    elif not is_resume and dist_ctx.enabled:
        logging.info(
            "[run-manifest] init skipped for non-leader distributed worker task_id=%d",
            dist_ctx.task_id,
        )
    else:
        logging.info("[resume] Skipping manifest initialization for run_id=%s", run_id)

    cfg.setdefault("SINGLE_LIGAND", "")
    cfg.setdefault("SINGLE_LIGAND_SUGGESTIONS", 5)
    cfg.setdefault("ALLOW_FDA_FALLBACK", False)
    cfg.setdefault("LIBRARY_MANIFEST_FILENAME", "_manifest.json")
    cfg.setdefault("LIBRARY_MANIFEST_BUILD_ON_SCAN", True)
    env_build_flag = os.environ.get("LIBRARY_MANIFEST_BUILD_ON_SCAN")
    if env_build_flag is not None:
        try:
            cfg["LIBRARY_MANIFEST_BUILD_ON_SCAN"] = _to_bool(env_build_flag)
        except Exception:
            pass

    logging.getLogger("lib-manifest").info(
        "[lib-manifest.scan] build_on_scan=%s",
        str(bool(cfg.get("LIBRARY_MANIFEST_BUILD_ON_SCAN", True))).lower(),
    )
    logging.getLogger("distributed.chunk").info(
        "[distributed.chunk.planner] manifest_mode=%s sample_n=%d scan_workers=%d prebuild=%s",
        _chunk_planner_manifest_mode(cfg),
        _chunk_planner_manifest_sample_size(cfg),
        _chunk_planner_scan_workers(cfg),
        str(_manifest_prebuild_enabled(cfg, distributed_enabled=bool(dist_ctx.enabled))).lower(),
    )
    cfg.setdefault(
        "FDA_MAPPING_CSV",
        str(Path(__file__).resolve().parents[2] / "chemdb" / "data" / "fda_mapping_from_pdbqt.csv"),
    )

    cli_single = _parse_single_from_cli(argv_for_parsing)
    env_single = os.environ.get("SINGLE_LIGAND", "").strip()
    cfg_single = str(cfg.get("SINGLE_LIGAND", "")).strip()
    effective_single = next((x for x in (cli_single, env_single, cfg_single) if x), "")
    cfg["_EFFECTIVE_SINGLE_LIGAND"] = effective_single
    if effective_single:
        print(
            f"[config] SINGLE_LIGAND effective='{effective_single}' "
            f"(order=CLI>{'ENV' if env_single else ''}>{'CFG' if cfg_single else ''})"
        )

    cfg.setdefault("LIBRARY_SUBDIR_DEFAULT", "fda_library")
    cfg.setdefault("HMDB_LIBRARY_SUBDIR", "hmdb")
    cfg.setdefault("TEST_MODE_ENABLE", "off")
    if "TEST_LIBRARY_MAP" not in cfg:
        cfg["TEST_LIBRARY_MAP"] = {}

    cfg.setdefault("TEST_FDA_LIBRARY_SUBDIR", _DEFAULT_TEST_FDA_LIBRARY_SUBDIR)
    if _cli_has(argv_for_parsing, "-test") or _cli_has(argv_for_parsing, "--test"):
        fixture_input_root = (
            Path(__file__).resolve().parents[2]
            / "chemdb"
            / "tests"
            / "fixtures"
            / "input_pdbs"
        )
        if fixture_input_root.joinpath("TEST.pdb").exists() and "INPUT_DIR" not in os.environ:
            cfg["INPUT_DIR"] = str(fixture_input_root)
        cfg["LIBRARY_SUBDIR_DEFAULT"] = cfg.get(
            "TEST_FDA_LIBRARY_SUBDIR",
            _DEFAULT_TEST_FDA_LIBRARY_SUBDIR,
        )
        _retarget_default_test_library_map_for_test_fda(
            cfg,
            str(cfg["LIBRARY_SUBDIR_DEFAULT"]),
        )
    if _cli_has(argv_for_parsing, "-test-fda") or _cli_has(argv_for_parsing, "--test-fda"):
        cfg["LIBRARY_SUBDIR_DEFAULT"] = cfg.get(
            "TEST_FDA_LIBRARY_SUBDIR",
            _DEFAULT_TEST_FDA_LIBRARY_SUBDIR,
        )
        _retarget_default_test_library_map_for_test_fda(
            cfg,
            str(cfg["LIBRARY_SUBDIR_DEFAULT"]),
        )
        print(
            f"[config] TEST_FDA_LIBRARY enabled: LIBRARY_SUBDIR_DEFAULT={cfg['LIBRARY_SUBDIR_DEFAULT']}"
        )

    cfg["FAST_MODE"] = _parse_fast_flag(argv_for_parsing) or bool(cfg.get("FAST_MODE", False))
    if cfg["FAST_MODE"]:
        print("[config] FAST_MODE effective=True (exhaustiveness=1)")

    cfg.setdefault("CONTROL_CONSENSUS", False)
    control_consensus = False
    try:
        control_consensus = _to_bool(cfg.get("CONTROL_CONSENSUS", False))
    except Exception:
        control_consensus = bool(cfg.get("CONTROL_CONSENSUS", False))
    if _cli_has(argv_for_parsing, "-control-consensus") or _cli_has(
        argv_for_parsing, "--control-consensus"
    ):
        control_consensus = True
    cfg["CONTROL_CONSENSUS"] = bool(control_consensus)
    if cfg["CONTROL_CONSENSUS"]:
        print("[config] CONTROL_CONSENSUS effective=True (multi-engine control docking enabled)")

    cfg.setdefault("NO_LIBRARY_DOCKING", False)
    cli_no_dock = _cli_has(argv_for_parsing, "--no-docking")
    env_no_dock = os.environ.get("NO_LIBRARY_DOCKING")
    effective_no_dock = False
    if cli_no_dock:
        effective_no_dock = True
    elif env_no_dock is not None:
        try:
            effective_no_dock = _to_bool(env_no_dock)
        except Exception:
            effective_no_dock = False
    else:
        try:
            effective_no_dock = bool(cfg.get("NO_LIBRARY_DOCKING", False))
        except Exception:
            effective_no_dock = False
    cfg["NO_LIBRARY_DOCKING"] = effective_no_dock
    if effective_no_dock:
        print(
            "[config] NO_LIBRARY_DOCKING=True (controls-only; skip DUD/FDA docking, but still enumerate ligands)"
        )

    retention_mode, retention_source = resolve_artifact_retention_mode(cfg, argv_for_parsing)
    cfg["ARTIFACT_RETENTION"] = retention_mode
    print(f"[config] ARTIFACT_RETENTION effective={retention_mode} source={retention_source}")

    cfg.setdefault("CENTER_MODE", "control-first")
    cfg.setdefault("CONTROL_BLACKLIST", "GOL,EDO,PG4,MPD,ACT,SO4,PO4,CL,NA,CA")
    cfg.setdefault("CONTROL_MIN_HEAVY_ATOMS", 10)
    cfg.setdefault("CONTROL_ANCHOR_MIN_VALID_RATE", 0.10)
    cfg.setdefault("ALLOW_SWITCH_FROM_CONTROL", True)
    cfg.setdefault("REQUIRE_CONTROL_FAILURE_TO_SWITCH", False)
    cfg.setdefault("SWITCH_AWAY_FROM_CONTROL_MIN_BOOST", 2.5)
    cfg.setdefault("CONTROL_LOCK_SCORE_MAX", -6.0)
    cfg.setdefault("CONTROL_LOCK_MIN_HITS", 1)
    cfg.setdefault("CONTROL_LOCK_CENTER_MAX_DIST", 4.0)
    cfg.setdefault("CLUSTER_EPS_ANG", 3.5)
    cfg.setdefault("SWITCH_VALID_RATE_MIN", 0.40)
    cfg.setdefault("SWITCH_SCORE_IMPROVE_MIN", 1.5)
    cfg.setdefault("SWITCH_LE_GAIN_MIN", 0.02)
    cfg.setdefault("SWITCH_SCORE_THRESHOLD", 0.70)
    cfg.setdefault("SWITCH_SCORE_HYSTERESIS", 0.50)
    cfg.setdefault("MAX_GLOBAL_CENTER_SWITCHES", 2)
    cfg.setdefault("THREADS_PER_VINA", 1)
    cfg.setdefault("RMSD_FILTER_ANG", 2.0)
    cfg.setdefault("RMSD_MAX_MODELS", 3)
    cfg.setdefault("ENABLE_GLOBAL_SCHEDULER", True)
    cfg.setdefault("GLOBAL_SCHEDULER_CPUS", int(cfg.get("CPU", os.cpu_count() or 1)))
    cfg.setdefault("GLOBAL_MAX_PROTEINS", 0)
    cfg.setdefault("GLOBAL_MIN_PROTEINS", 1)
    cfg.setdefault("GLOBAL_ADMISSION_POLICY", "adaptive")
    cfg.setdefault("RETRY_NEAR_MISS", True)
    cfg.setdefault("RETRY_EXHAUST_MULT", 2)
    cfg.setdefault("RETRY_DIST_THRESH", 7.0)
    cfg.setdefault("RETRY_SCORE_THRESH", -7.5)
    cfg.setdefault("ADAPTIVE_SHRINK_ENABLE", True)
    cfg.setdefault("ADAPTIVE_SHRINK_MEDIAN_MAX", 4.0)
    cfg.setdefault("ADAPTIVE_SHRINK_DEC", 4.0)
    cfg.setdefault("ADAPTIVE_SHRINK_MIN_BOX", 14.0)
    cfg.setdefault("CHECKPOINT_ENABLE", True)
    cfg.setdefault("PAINS_NAMES", "")
    cfg.setdefault("RECEPTOR_SANITY_CHECK", True)
    cfg.setdefault("AUDIT_JSON", True)
    cfg.setdefault("QUIET_CONSOLE", False)
    cfg.setdefault("VINA_VERBOSITY", 2)
    cfg.setdefault("FILTER_VINA_STDOUT", False)
    cfg.setdefault("SELF_RMSD_MAX_ANG", 2.0)
    cfg.setdefault("SELF_RMSD_REQUIRE_FOR_CONTROLS", True)
    cfg.setdefault("EARLY_EXIT_MAX_MODELS", 3)
    cfg.setdefault("MAX_RETRY_SECONDS_PER_LIGAND", 300)
    return global_start


def _select_workload(
    *,
    cfg: ConfigDict,
    argv_for_parsing: list[str],
    is_resume: bool,
    resume_protein_ids: list[str],
    dist_ctx: Any,
    select_pdb_files_for_run: Callable[..., list[str]],
    get_recenter_params: Callable[[dict], Any],
) -> WorkloadSelectionState:
    params = get_recenter_params(cfg)
    stages = define_docking_stages(cfg.get("DOCKING_MODE", "discovery").lower())
    stages = _run_profiles.apply_benchmark_vina_seeds(cfg, stages)
    if _to_bool(cfg.get("WATER_BENCH_STAGE1_ONLY", False)):
        stages = stages[:1]
        logging.info("[water-bench] stage_policy=stage1_only")
    logging.info("[config] docking_mode=%s", cfg.get("DOCKING_MODE"))

    pdb_files = select_pdb_files_for_run(
        cfg,
        argv_for_parsing,
        is_resume=is_resume,
        resume_protein_ids=resume_protein_ids,
    )
    all_pdb_files = list(pdb_files)
    dist_combo_chunks_disabled = False
    try:
        dist_combo_chunks_disabled = _to_bool(os.environ.get("ATLAS_DISABLE_COMBO_CHUNKS", "0"))
    except Exception:
        dist_combo_chunks_disabled = False

    requested_execution_mode = os.environ.get("ATLAS_EXECUTION_MODE")
    execution_mode_decision = resolve_execution_mode(
        cfg=cfg,
        distributed_enabled=bool(dist_ctx.enabled),
        task_count=int(getattr(dist_ctx, "task_count", 1) or 1),
        pdb_count=int(len(all_pdb_files)),
        combo_disabled=bool(dist_combo_chunks_disabled),
        requested_mode_raw=requested_execution_mode,
    )
    dist_combo_chunk_mode = bool(dist_ctx.enabled and execution_mode_decision.combo_enabled)
    dist_combo_hybrid_mode = bool(dist_ctx.enabled and execution_mode_decision.hybrid_enabled)
    dist_per_protein_mode = bool(
        dist_ctx.enabled and execution_mode_decision.effective_mode == "distributed_per_protein"
    )
    logging.info(
        "[execution.mode] requested=%s effective=%s reason=%s pdb_count=%d distributed=%s task_count=%d imbalance=%.3f",
        execution_mode_decision.requested_mode,
        execution_mode_decision.effective_mode,
        execution_mode_decision.reason,
        int(execution_mode_decision.pdb_count),
        str(bool(dist_ctx.enabled)).lower(),
        int(getattr(dist_ctx, "task_count", 1) or 1),
        float(execution_mode_decision.imbalance_ratio),
    )
    if dist_per_protein_mode:
        pdb_files = shard_pdb_files_for_context(dist_ctx, pdb_files)
    elif dist_combo_chunk_mode:
        pdb_files = list(all_pdb_files)

    tokens = parse_test_libraries(cfg)
    test_mode = _resolve_test_mode(cfg)

    if dist_per_protein_mode:
        if dist_combo_chunks_disabled:
            logging.info(
                "[distributed.combo-chunk] status=disabled reason=ATLAS_DISABLE_COMBO_CHUNKS"
            )
        print(
            f"Proteins queued (global): {len(all_pdb_files)} | assigned to task {dist_ctx.task_id}: {len(pdb_files)}"
        )
    elif dist_combo_chunk_mode:
        print(
            f"Proteins queued (global): {len(all_pdb_files)} | task {dist_ctx.task_id} uses {'combo-hybrid' if dist_combo_hybrid_mode else 'combo'} assignment"
        )
    else:
        print(f"Proteins queued: {len(pdb_files)}")

    return WorkloadSelectionState(
        params=params,
        stages=stages,
        pdb_files=list(pdb_files),
        all_pdb_files=list(all_pdb_files),
        tokens=list(tokens),
        test_mode=str(test_mode),
        dist_combo_chunk_mode=bool(dist_combo_chunk_mode),
        dist_combo_hybrid_mode=bool(dist_combo_hybrid_mode),
        execution_mode_decision=execution_mode_decision,
    )


def _record_scheduled_proteins_for_manifest(
    *,
    cfg: ConfigDict,
    run_id: str,
    is_resume: bool,
    dist_ctx: Any,
    all_pdb_files: list[str],
) -> None:
    if is_resume or (dist_ctx.enabled and not dist_ctx.is_leader):
        return
    try:
        scheduled_ids: list[str] = []
        for f in all_pdb_files:
            nid = _norm_pdb_id(f)
            if nid:
                scheduled_ids.append(nid.upper())
        update_manifest_for_scheduled_proteins(cfg, run_id, scheduled_ids)
    except Exception:
        print(
            "[run-manifest] WARNING: failed to record scheduled proteins in manifest",
            file=sys.stderr,
        )


def _prepare_variant_bootstrap(
    *,
    cfg: ConfigDict,
    run_id: str,
    dist_ctx: Any,
) -> VariantBootstrapState:
    start = time.time()
    plan_only = os.environ.get("A2_PLAN_ONLY") == "1"
    mode, variants = resolve_apo_holo_mode(cfg)
    router_legacy = mode == "legacy"
    ph_enabled = bool(cfg.get("PH_ENSEMBLE"))
    cfg_raw_mode = cfg.get("APO_HOLO_MODE")
    logging.info(
        "[apo-holo] cfg_token_raw=%r resolved_mode=%s variants=%s",
        cfg_raw_mode,
        mode,
        variants,
    )
    logging.info("[apo-holo] router_legacy=%s ph_enabled=%s", router_legacy, ph_enabled)
    logging.info("[apo-holo] normalized_mode_token=%s", _debug_normalize_mode_token(cfg_raw_mode))
    cfg["_ROUTER_LEGACY"] = router_legacy
    cfg["_RESOLVED_APO_HOLO_MODE"] = mode
    logging.info(
        "[apo-holo] resolved mode=%s variants=%s cfg_token=%r",
        mode,
        variants,
        cfg_raw_mode,
    )
    if not dist_ctx.enabled or dist_ctx.is_leader:
        try:
            update_manifest_for_run_config(cfg, run_id)
        except Exception:
            logging.warning(
                "[run-manifest] Failed to record run-level apo/holo + water policy",
                exc_info=True,
            )

    overall_dir = cfg.get("OVERALL_DIR", ".")
    failed_root = os.path.join(overall_dir, "failed")
    os.makedirs(failed_root, exist_ok=True)
    logging.info("[apo-holo] failed log directory: %s", failed_root)
    return VariantBootstrapState(
        start=float(start),
        plan_only=bool(plan_only),
        mode=str(mode),
        variants=list(variants),
        failed_root=str(failed_root),
    )


def bootstrap_pre_variant_state(
    *,
    cfg: ConfigDict,
    run_id: str,
    argv_for_parsing: list[str],
    is_resume: bool,
    is_bench: bool,
    is_bench2: bool,
    is_bench_small: bool,
    resume_protein_ids: list[str],
    completed_combo_lookup: set[tuple[str, str, str]],
    log_path: str,
    bench_test_library_map: Mapping[str, str],
    bench_small_test_library_map: Mapping[str, str],
    bench_small_target_ligands: int,
    resolve_artifact_retention_mode: Callable[[Mapping[str, Any], list[str]], tuple[str, str]],
    select_pdb_files_for_run: Callable[..., list[str]],
    get_recenter_params: Callable[[dict], Any],
) -> BootstrapState:
    dist_ctx = _setup_distributed_runtime(cfg, run_id=run_id)
    _log_atom_alias_summary()
    _apply_ph_ensemble_settings(cfg)
    argv_for_parsing = _hydrate_resume_state(
        cfg=cfg,
        run_id=run_id,
        is_resume=is_resume,
        argv_for_parsing=argv_for_parsing,
        resume_protein_ids=resume_protein_ids,
        completed_combo_lookup=completed_combo_lookup,
    )
    if is_resume and not str(cfg.get("_DUD_LIBRARY_OVERRIDE") or "").strip():
        try:
            restored_dud_library = _run_profiles._dud_library_value(argv_for_parsing)
        except ValueError:
            raise
        except Exception:
            restored_dud_library = ""
        if restored_dud_library:
            cfg["_DUD_LIBRARY_OVERRIDE"] = restored_dud_library
            logging.info(
                "[resume] restored --dud-library override from stored argv: %s",
                restored_dud_library,
            )
    global_start = _apply_runtime_config_defaults(
        cfg=cfg,
        run_id=run_id,
        argv_for_parsing=argv_for_parsing,
        is_resume=is_resume,
        dist_ctx=dist_ctx,
        is_bench=is_bench,
        is_bench2=is_bench2,
        is_bench_small=is_bench_small,
        bench_test_library_map=bench_test_library_map,
        bench_small_test_library_map=bench_small_test_library_map,
        bench_small_target_ligands=bench_small_target_ligands,
        resolve_artifact_retention_mode=resolve_artifact_retention_mode,
        log_path=log_path,
    )
    workload = _select_workload(
        cfg=cfg,
        argv_for_parsing=argv_for_parsing,
        is_resume=is_resume,
        resume_protein_ids=resume_protein_ids,
        dist_ctx=dist_ctx,
        select_pdb_files_for_run=select_pdb_files_for_run,
        get_recenter_params=get_recenter_params,
    )
    _record_scheduled_proteins_for_manifest(
        cfg=cfg,
        run_id=run_id,
        is_resume=is_resume,
        dist_ctx=dist_ctx,
        all_pdb_files=workload.all_pdb_files,
    )
    variant_state = _prepare_variant_bootstrap(cfg=cfg, run_id=run_id, dist_ctx=dist_ctx)

    failed_entries: list[tuple[str, str, str, str, str]] = []
    run_scope_completion_times: dict[tuple[str, str, str], float] = {}
    run_chunk_rebalance_count = 0
    recorded_terminal_chunk_ids: set[str] = set()
    distributed_assigned_pdb_ids: set[str] = set()
    retention_lock = threading.Lock()
    pending_retention_combos: set[tuple[str, str, str]] = set()
    pending_coverage_refresh: dict[tuple[str, str, str], dict[str, Any]] = {}
    pending_retention_lock = threading.Lock()
    scorch_queue_service = _start_scorch_mixed_queue_service(cfg, run_id, verbose=False)

    return BootstrapState(
        cfg=cfg,
        argv_for_parsing=list(argv_for_parsing),
        dist_ctx=dist_ctx,
        completed_combo_lookup=set(completed_combo_lookup),
        params=workload.params,
        stages=workload.stages,
        pdb_files=list(workload.pdb_files),
        tokens=list(workload.tokens),
        test_mode=str(workload.test_mode),
        dist_combo_chunk_mode=bool(workload.dist_combo_chunk_mode),
        dist_combo_hybrid_mode=bool(workload.dist_combo_hybrid_mode),
        execution_mode_decision=workload.execution_mode_decision,
        plan_only=bool(variant_state.plan_only),
        mode=str(variant_state.mode),
        variants=list(variant_state.variants),
        start=float(variant_state.start),
        global_start=float(global_start),
        failed_root=str(variant_state.failed_root),
        failed_entries=failed_entries,
        run_scope_completion_times=run_scope_completion_times,
        run_chunk_rebalance_count=int(run_chunk_rebalance_count),
        recorded_terminal_chunk_ids=recorded_terminal_chunk_ids,
        distributed_assigned_pdb_ids=distributed_assigned_pdb_ids,
        retention_lock=retention_lock,
        pending_retention_combos=pending_retention_combos,
        pending_coverage_refresh=pending_coverage_refresh,
        pending_retention_lock=pending_retention_lock,
        scorch_queue_service=scorch_queue_service,
    )
