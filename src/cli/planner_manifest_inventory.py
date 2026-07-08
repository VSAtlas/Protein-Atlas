from __future__ import annotations

import json
import logging
import os
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from docking.active_preserving_cap import stable_active_preserving_limit
from docking.docking_ligands import _chunk_ligand_key
from docking.library_mode import compute_allowed_library_roots

_RESERVED_TEST_MODE_TOKENS = {"dud", "fda", "hmdb"}
_PH_DIR_NUM_PATTERN = re.compile(r"^[+-]?\d+(?:\.\d+)?$")
_CHUNK_MANIFEST_SAMPLE_SIZE = 64
_CHUNK_PLANNER_SCAN_WORKERS_DEFAULT = 8


def _ph_ligand_mode_enabled(cfg: Mapping[str, Any]) -> bool:
    raw = str(cfg.get("PH_LIGAND_MODE", "off")).strip().lower()
    return raw not in ("", "off", "none", "false", "0")


def _parse_ph_dir_value(dir_name: str) -> Optional[float]:
    lowered = str(dir_name or "").strip().lower()
    if not lowered:
        return None

    def _coerce(candidate: str) -> Optional[float]:
        try:
            value = float(candidate)
        except Exception:
            return None
        return value if 0.0 <= value <= 14.0 else None

    if lowered.startswith("ph"):
        remainder = lowered[2:].lstrip("_- ")
        if remainder and remainder[0] in "+-.0123456789":
            for candidate in (remainder, remainder.replace("_", ".")):
                if _PH_DIR_NUM_PATTERN.fullmatch(candidate):
                    parsed = _coerce(candidate)
                    if parsed is not None:
                        return parsed
            match = re.search(r"([+-]?\d+(?:\.\d+)?)", remainder)
            if match:
                parsed = _coerce(match.group(1))
                if parsed is not None:
                    return parsed

    cleaned = lowered.strip("_- ")
    cleaned_alt = cleaned.replace("_", ".")
    for candidate in (cleaned, cleaned_alt):
        if candidate and _PH_DIR_NUM_PATTERN.fullmatch(candidate):
            parsed = _coerce(candidate)
            if parsed is not None:
                return parsed
    return None


def _is_under_ph_subdir(path: Path, root: Path) -> bool:
    try:
        rel = path.resolve().relative_to(root.resolve())
    except Exception:
        return False
    if not rel.parts:
        return False
    return _parse_ph_dir_value(str(rel.parts[0])) is not None


def _keep_plannable_ligand(path: Path, *, root: Path, ph_ligand_mode_on: bool) -> bool:
    if ph_ligand_mode_on:
        return True
    if any(part.lower() == "microstates" for part in path.parts):
        return False
    if _is_under_ph_subdir(path, root):
        return False
    if "pH" in path.stem:
        return False
    return True


def _keep_plannable_ligand_rel(
    rel_norm: str, *, ph_ligand_mode_on: bool
) -> bool:
    if ph_ligand_mode_on:
        return True
    rel_path = Path(rel_norm)
    parts_lower = [str(p).lower() for p in rel_path.parts]
    if any(part == "microstates" for part in parts_lower):
        return False
    if rel_path.parts:
        if _parse_ph_dir_value(str(rel_path.parts[0])) is not None:
            return False
    if "pH" in rel_path.stem:
        return False
    return True


def _child_dir_mtimes(root: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    try:
        for child in root.iterdir():
            if not child.is_dir():
                continue
            try:
                out[child.name] = float(child.stat().st_mtime)
            except Exception:
                continue
    except Exception:
        return {}
    return out


def _manifest_inventory_cache_key(
    root: Path, *, ph_ligand_mode_on: bool, manifest_mode: str
) -> str:
    manifest_path = root / "_manifest.json"
    try:
        root_key = str(root.resolve())
    except Exception:
        root_key = str(root)
    try:
        st = manifest_path.stat()
        manifest_sig = f"{int(st.st_mtime_ns)}:{int(st.st_size)}"
    except Exception:
        manifest_sig = "missing"
    return (
        f"{root_key}|sig={manifest_sig}|ph_mode={int(ph_ligand_mode_on)}|"
        f"mode={manifest_mode}"
    )


def _read_library_manifest_inventory(
    cfg: Mapping[str, Any],
    root: Path,
    *,
    ph_ligand_mode_on: bool,
    manifest_mode: str,
) -> tuple[list[str], list[Path]]:
    cache: Optional[dict[str, tuple[list[str], list[Path]]]] = None
    if isinstance(cfg, dict):
        cached = cfg.get("_PLANNER_MANIFEST_INVENTORY_CACHE")
        if isinstance(cached, dict):
            cache = cached
        else:
            cache = {}
            cfg["_PLANNER_MANIFEST_INVENTORY_CACHE"] = cache
    cache_key = _manifest_inventory_cache_key(
        root,
        ph_ligand_mode_on=ph_ligand_mode_on,
        manifest_mode=manifest_mode,
    )
    if isinstance(cache, dict):
        cached_payload = cache.get(cache_key)
        if isinstance(cached_payload, tuple) and len(cached_payload) == 2:
            if isinstance(cfg, dict):
                cfg["_PLANNER_MANIFEST_CACHE_HITS"] = int(
                    cfg.get("_PLANNER_MANIFEST_CACHE_HITS", 0) or 0
                ) + 1
            return list(cached_payload[0]), list(cached_payload[1])

    if isinstance(cfg, dict):
        cfg["_PLANNER_MANIFEST_CACHE_MISSES"] = int(
            cfg.get("_PLANNER_MANIFEST_CACHE_MISSES", 0) or 0
        ) + 1

    manifest_path = root / "_manifest.json"
    if not manifest_path.exists():
        return [], []
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8")) or {}
    except Exception:
        return [], []
    rel_values: list[str] = []
    raw_files = payload.get("filenames")
    if isinstance(raw_files, dict):
        rel_values.extend(str(v) for v in raw_files.values())
    raw_entries = payload.get("entries")
    if isinstance(raw_entries, dict):
        rel_values.extend(str(v) for v in raw_entries.values())
    if not rel_values:
        return [], []

    trust_manifest_rel_only = False
    if manifest_mode == "trust" and isinstance(payload, dict):
        recorded_mtime = payload.get("mtime")
        recorded_child_dir_mtimes = payload.get("child_dir_mtimes")
        try:
            current_mtime = float(root.stat().st_mtime)
        except Exception:
            current_mtime = None
        child_dirs_match = (
            isinstance(recorded_child_dir_mtimes, dict)
            and {
                str(k): float(v)
                for k, v in recorded_child_dir_mtimes.items()
            }
            == {
                str(k): float(v)
                for k, v in _child_dir_mtimes(root).items()
            }
        )
        if (
            recorded_mtime is not None
            and current_mtime is not None
            and float(recorded_mtime) == float(current_mtime)
            and child_dirs_match
        ):
            trust_manifest_rel_only = True

    candidate_paths: list[Path] = []
    bases: list[str] = []
    seen_rel: set[str] = set()
    for rel in rel_values:
        rel_norm = str(rel).replace("\\", "/").strip()
        if not rel_norm or rel_norm in seen_rel:
            continue
        seen_rel.add(rel_norm)
        rel_path = root / rel_norm
        candidate_paths.append(rel_path)
        if trust_manifest_rel_only:
            if not _keep_plannable_ligand_rel(
                rel_norm,
                ph_ligand_mode_on=ph_ligand_mode_on,
            ):
                continue
            base = _chunk_ligand_key(Path(rel_norm).name)
        else:
            try:
                if not rel_path.exists() or rel_path.stat().st_size <= 100:
                    continue
            except Exception:
                continue
            if not _keep_plannable_ligand(
                rel_path, root=root, ph_ligand_mode_on=ph_ligand_mode_on
            ):
                continue
            base = _chunk_ligand_key(rel_path.name)
        if base:
            bases.append(base)
    result = (sorted(set(bases)), candidate_paths)
    if isinstance(cache, dict):
        cache[cache_key] = (list(result[0]), list(result[1]))
    return result


def _normalize_chunk_manifest_mode(raw: Any) -> str:
    token = str(raw or "").strip().lower().replace("-", "_")
    if token in {"", "trust", "manifest", "manifest_trust"}:
        return "trust"
    if token in {"sample", "manifest_sample"}:
        return "sample"
    if token in {"strict", "manifest_strict", "verify"}:
        return "strict"
    return "trust"


def _chunk_planner_manifest_mode(cfg: Mapping[str, Any]) -> str:
    env_val = os.environ.get("ATLAS_CHUNK_PLANNER_MANIFEST_MODE")
    if env_val is not None and str(env_val).strip():
        return _normalize_chunk_manifest_mode(env_val)
    return _normalize_chunk_manifest_mode(cfg.get("CHUNK_PLANNER_MANIFEST_MODE", "trust"))


def _chunk_planner_manifest_sample_size(cfg: Mapping[str, Any]) -> int:
    env_val = os.environ.get("ATLAS_CHUNK_PLANNER_MANIFEST_SAMPLE_SIZE")
    raw = env_val if env_val is not None else cfg.get(
        "CHUNK_PLANNER_MANIFEST_SAMPLE_SIZE", _CHUNK_MANIFEST_SAMPLE_SIZE
    )
    try:
        return max(1, int(raw))
    except Exception:
        return int(_CHUNK_MANIFEST_SAMPLE_SIZE)


def _chunk_planner_scan_workers(cfg: Mapping[str, Any]) -> int:
    env_val = os.environ.get("ATLAS_CHUNK_PLANNER_SCAN_WORKERS")
    raw = env_val if env_val is not None else cfg.get(
        "CHUNK_PLANNER_SCAN_WORKERS", _CHUNK_PLANNER_SCAN_WORKERS_DEFAULT
    )
    try:
        parsed = int(raw)
    except Exception:
        parsed = int(_CHUNK_PLANNER_SCAN_WORKERS_DEFAULT)
    parsed = max(1, parsed)
    return max(1, min(parsed, int(os.cpu_count() or 1)))


def _manifest_sample_looks_valid(
    *,
    candidate_paths: list[Path],
    sample_size: int,
) -> bool:
    if not candidate_paths:
        return False
    probes = candidate_paths[: max(1, min(int(sample_size), len(candidate_paths)))]
    ok = 0
    for path in probes:
        try:
            if path.exists() and path.stat().st_size > 100:
                ok += 1
        except Exception:
            continue
    return ok > 0 and ok == len(probes)


def _scan_library_file_bases(
    cfg: Mapping[str, Any], root: Path, *, ph_ligand_mode_on: bool
) -> list[str]:
    cache: Optional[dict[str, list[str]]] = None
    if isinstance(cfg, dict):
        cached = cfg.get("_PLANNER_FILE_BASES_CACHE")
        if isinstance(cached, dict):
            cache = cached
        else:
            cache = {}
            cfg["_PLANNER_FILE_BASES_CACHE"] = cache
    try:
        cache_key = f"{str(root.resolve())}|ph_mode={int(ph_ligand_mode_on)}"
    except Exception:
        cache_key = f"{str(root)}|ph_mode={int(ph_ligand_mode_on)}"
    if isinstance(cache, dict) and cache_key in cache:
        return list(cache[cache_key])

    out: list[str] = []
    for path in root.rglob("*.pdbqt"):
        try:
            if not path.exists() or path.stat().st_size <= 100:
                continue
        except Exception:
            continue
        if not _keep_plannable_ligand(path, root=root, ph_ligand_mode_on=ph_ligand_mode_on):
            continue
        base = _chunk_ligand_key(path.name)
        if base:
            out.append(base)
    result = sorted(set(out))
    if isinstance(cache, dict):
        cache[cache_key] = list(result)
    return result


def _scan_library_file_bases_parallel(
    cfg: Mapping[str, Any],
    *,
    roots: Sequence[Path],
    ph_ligand_mode_on: bool,
) -> dict[Path, list[str]]:
    unique_roots: list[Path] = []
    seen: set[str] = set()
    for root in roots:
        try:
            key = str(root.resolve())
        except Exception:
            key = str(root)
        if key in seen:
            continue
        seen.add(key)
        unique_roots.append(Path(root))
    if not unique_roots:
        return {}
    workers = min(_chunk_planner_scan_workers(cfg), len(unique_roots))
    if workers <= 1:
        return {
            root: _scan_library_file_bases(
                cfg, root, ph_ligand_mode_on=ph_ligand_mode_on
            )
            for root in unique_roots
        }
    out: dict[Path, list[str]] = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        fut_to_root = {
            pool.submit(
                _scan_library_file_bases,
                cfg,
                root,
                ph_ligand_mode_on=ph_ligand_mode_on,
            ): root
            for root in unique_roots
        }
        for fut, root in fut_to_root.items():
            try:
                out[root] = fut.result()
            except Exception:
                out[root] = []
    return out


def _pdbqt_complexity(path: Path) -> tuple[int, int] | None:
    heavy_atoms = 0
    rotatable = 0
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith("BRANCH"):
                    rotatable += 1
                    continue
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                element = line[76:78].strip() if len(line) >= 78 else ""
                if not element:
                    atom_name = line[12:16].strip() if len(line) >= 16 else ""
                    element = "".join(ch for ch in atom_name if ch.isalpha())[:1]
                if element.upper() != "H":
                    heavy_atoms += 1
    except Exception:
        return None
    if heavy_atoms <= 0:
        return None
    return heavy_atoms, rotatable


def _rdkit_complexity(path: Path) -> tuple[int, int] | None:
    try:
        from rdkit import Chem  # type: ignore[import-untyped]
        from rdkit.Chem import Lipinski  # type: ignore[import-untyped]
    except Exception:
        return None
    mol = None
    suffix = path.suffix.lower()
    try:
        if suffix in {".sdf", ".mol"}:
            mol = Chem.MolFromMolFile(str(path), sanitize=False, removeHs=False)
        elif suffix == ".mol2":
            mol = Chem.MolFromMol2File(str(path), sanitize=False, removeHs=False)
        elif suffix == ".pdb":
            mol = Chem.MolFromPDBFile(str(path), sanitize=False, removeHs=False)
    except Exception:
        mol = None
    if mol is None:
        return None
    try:
        heavy_atoms = int(mol.GetNumHeavyAtoms())
        rotatable = int(Lipinski.NumRotatableBonds(mol))
    except Exception:
        return None
    if heavy_atoms <= 0:
        return None
    return heavy_atoms, rotatable


def _ligand_complexity_weight(path: Path) -> float:
    parsed = _pdbqt_complexity(path) if path.suffix.lower() == ".pdbqt" else None
    if parsed is None:
        parsed = _rdkit_complexity(path)
    if parsed is None:
        return 1.0
    heavy_atoms, rotatable = parsed
    return max(0.25, 1.0 + (float(heavy_atoms) / 25.0) + (float(rotatable) / 5.0))


def _cache_ligand_complexity_weights(
    cfg: Mapping[str, Any],
    *,
    selected_bases: Sequence[str],
    candidate_paths: Sequence[Path],
) -> None:
    if not isinstance(cfg, dict) or not selected_bases or not candidate_paths:
        return
    selected = {str(base) for base in selected_bases if str(base).strip()}
    if not selected:
        return
    weights = cfg.get("_PLANNER_LIGAND_COMPLEXITY_WEIGHTS")
    if not isinstance(weights, dict):
        weights = {}
        cfg["_PLANNER_LIGAND_COMPLEXITY_WEIGHTS"] = weights
    path_cache = cfg.get("_PLANNER_LIGAND_COMPLEXITY_PATH_CACHE")
    if not isinstance(path_cache, dict):
        path_cache = {}
        cfg["_PLANNER_LIGAND_COMPLEXITY_PATH_CACHE"] = path_cache
    for path in candidate_paths:
        if all(base in weights for base in selected):
            break
        base = _chunk_ligand_key(path.name)
        if not base or base not in selected or base in weights:
            continue
        try:
            key = str(path.resolve())
        except Exception:
            key = str(path)
        cached = path_cache.get(key)
        if cached is None:
            cached = _ligand_complexity_weight(path)
            path_cache[key] = float(cached)
        weights[base] = float(cached)


def _path_identity(path: Path) -> str:
    try:
        return str(path.resolve())
    except Exception:
        return str(path)


def _stream_tokens(run_tokens: Sequence[str]) -> list[str]:
    tokens = [str(token).strip().lower() for token in run_tokens if str(token).strip()]
    if not tokens:
        return ["fda"]
    out: list[str] = []
    seen: set[str] = set()
    for token in tokens:
        if token in {"default", "off", "none", "null"}:
            token = "fda"
        if token in seen:
            continue
        seen.add(token)
        out.append(token)
    return out or ["fda"]


def _stream_root_specs(
    cfg: Mapping[str, Any],
    *,
    pdb_id: str,
    run_tokens: Sequence[str],
    logger: logging.Logger,
) -> list[dict[str, str]]:
    cfg_for_roots: dict[str, Any]
    if isinstance(cfg, dict):
        cfg_for_roots = cfg
    else:
        cfg_for_roots = dict(cfg)

    specs: list[dict[str, str]] = []
    seen_roots: set[str] = set()
    for run_mode in _stream_tokens(run_tokens):
        roots = compute_allowed_library_roots(
            cfg_for_roots,
            pdb_id,
            logger,
            tokens_override=[run_mode],
        )
        for root in roots:
            if not root.exists():
                continue
            root_key = _path_identity(root)
            if root_key in seen_roots:
                logger.info(
                    "[distributed.chunk.stream-root] pdb=%s run_mode=%s root=%s action=skip_duplicate",
                    pdb_id,
                    run_mode,
                    root,
                )
                continue
            seen_roots.add(root_key)
            specs.append(
                {
                    "run_mode": str(run_mode),
                    "library_name": str(root.name),
                    "library_root": root_key,
                }
            )
    return specs


def _select_ligand_bases_for_root(
    cfg: Mapping[str, Any],
    *,
    pdb_id: str,
    root: Path,
    manifest_mode: str,
    sample_size: int,
    ph_ligand_mode_on: bool,
    logger: logging.Logger,
) -> tuple[list[str], list[Path]]:
    manifest_bases, manifest_paths = _read_library_manifest_inventory(
        cfg,
        root,
        ph_ligand_mode_on=ph_ligand_mode_on,
        manifest_mode=manifest_mode,
    )
    needs_scan = False
    if manifest_mode == "strict":
        needs_scan = True
    elif manifest_mode == "sample":
        needs_scan = not _manifest_sample_looks_valid(
            candidate_paths=manifest_paths,
            sample_size=sample_size,
        )
    else:
        needs_scan = not bool(manifest_bases)

    file_bases: list[str] = []
    if needs_scan:
        file_bases = _scan_library_file_bases(
            cfg,
            root,
            ph_ligand_mode_on=ph_ligand_mode_on,
        )

    manifest_set = set(manifest_bases)
    file_set = set(file_bases)
    selected_source = "manifest_trust"
    reason = f"mode={manifest_mode}"
    if needs_scan:
        selected = file_bases
        selected_source = "filesystem_scan"
        if not manifest_set:
            reason = "manifest_missing_or_empty"
        elif not file_set:
            reason = "scan_empty"
        else:
            reason = "scan_used"
    else:
        selected = manifest_bases
        if not manifest_set:
            selected_source = "manifest_unusable"
            reason = "manifest_empty_after_filter"

    logger.info(
        "[distributed.chunk.ligand-source] pdb=%s root=%s manifest_n=%d files_n=%d selected=%s reason=%s",
        pdb_id,
        root,
        len(manifest_set),
        len(file_set),
        selected_source,
        reason,
    )
    return sorted(set(selected)), list(manifest_paths)


def _resolve_combo_ligand_streams(
    cfg: Mapping[str, Any],
    *,
    pdb_id: str,
    run_tokens: list[str],
    target_ligands: int,
    sample_seed: int,
) -> list[dict[str, Any]]:
    logger = logging.getLogger("distributed.chunk")
    specs = _stream_root_specs(
        cfg,
        pdb_id=pdb_id,
        run_tokens=run_tokens,
        logger=logger,
    )
    manifest_mode = _chunk_planner_manifest_mode(cfg)
    sample_size = _chunk_planner_manifest_sample_size(cfg)
    ph_ligand_mode_on = _ph_ligand_mode_enabled(cfg)
    streams: list[dict[str, Any]] = []
    for spec in specs:
        root = Path(spec["library_root"])
        selected, candidate_paths = _select_ligand_bases_for_root(
            cfg,
            pdb_id=pdb_id,
            root=root,
            manifest_mode=manifest_mode,
            sample_size=sample_size,
            ph_ligand_mode_on=ph_ligand_mode_on,
            logger=logger,
        )
        available_n = len(selected)
        effective_target = (
            available_n if target_ligands <= 0 else min(target_ligands, available_n)
        )
        if target_ligands > 0 and available_n < target_ligands:
            logger.warning(
                "[distributed.chunk.target-clamp] pdb=%s run_mode=%s library=%s configured=%d available=%d effective=%d policy=auto_clamp_continue",
                pdb_id,
                spec["run_mode"],
                spec["library_name"],
                int(target_ligands),
                int(available_n),
                int(effective_target),
            )
        if effective_target > 0 and len(selected) > effective_target:
            selected = stable_active_preserving_limit(
                selected,
                limit=effective_target,
                seed=sample_seed,
                scope=f"{pdb_id}|{spec['run_mode']}|{spec['library_name']}|{spec['library_root']}",
                key=str,
            )
        _cache_ligand_complexity_weights(
            cfg,
            selected_bases=selected,
            candidate_paths=candidate_paths,
        )
        if not selected:
            continue
        streams.append(
            {
                "run_mode": spec["run_mode"],
                "library_name": spec["library_name"],
                "library_root": spec["library_root"],
                "ligand_bases": list(selected),
            }
        )

    if streams:
        return streams

    token_custom = [
        t for t in _stream_tokens(run_tokens) if t not in _RESERVED_TEST_MODE_TOKENS
    ]
    fallback_library = (
        token_custom[0]
        if token_custom
        else str(cfg.get("LIBRARY_SUBDIR_DEFAULT", "fda_library"))
    )
    return [
        {
            "run_mode": token_custom[0] if token_custom else "fda",
            "library_name": fallback_library,
            "library_root": "",
            "ligand_bases": [],
        }
    ]


def _validate_chunk_stream_resolution(
    cfg: Mapping[str, Any],
    *,
    chunks: Sequence[Mapping[str, Any]],
    sample_chunks_per_stream: int = 32,
    max_loss_fraction: float = 0.05,
) -> None:
    if not chunks:
        return
    logger = logging.getLogger("distributed.chunk")
    manifest_mode = _chunk_planner_manifest_mode(cfg)
    sample_size = _chunk_planner_manifest_sample_size(cfg)
    ph_ligand_mode_on = _ph_ligand_mode_enabled(cfg)
    grouped: dict[tuple[str, str, str], list[Mapping[str, Any]]] = {}
    for chunk in chunks:
        run_mode = str(chunk.get("run_mode") or "").strip()
        library_name = str(chunk.get("library_name") or "").strip()
        library_root = str(chunk.get("library_root") or "").strip()
        if not run_mode or not library_root:
            raise RuntimeError(
                "chunk_resolution_preflight_failed:"
                f" reason=missing_stream_identity chunk_id={chunk.get('chunk_id')}"
            )
        grouped.setdefault((run_mode, library_name, library_root), []).append(chunk)

    sample_n = max(1, int(sample_chunks_per_stream))
    allowed_loss = max(0.0, min(1.0, float(max_loss_fraction)))
    for (run_mode, library_name, library_root), stream_chunks in grouped.items():
        root = Path(library_root)
        available, _paths = _select_ligand_bases_for_root(
            cfg,
            pdb_id="chunk_preflight",
            root=root,
            manifest_mode=manifest_mode,
            sample_size=sample_size,
            ph_ligand_mode_on=ph_ligand_mode_on,
            logger=logger,
        )
        available_set = {str(base) for base in available}
        sampled = stream_chunks[: min(sample_n, len(stream_chunks))]
        empty_chunks = 0
        planned_total = 0
        resolved_total = 0
        missing_sample: list[str] = []
        for chunk in sampled:
            ligand_bases = [
                str(x)
                for x in (chunk.get("ligand_bases") or [])
                if str(x).strip()
            ]
            planned = len(ligand_bases)
            resolved = sum(1 for base in ligand_bases if base in available_set)
            planned_total += planned
            resolved_total += resolved
            if planned > 0 and resolved == 0:
                empty_chunks += 1
            if resolved < planned and len(missing_sample) < 8:
                for base in ligand_bases:
                    if base not in available_set:
                        missing_sample.append(base)
                        if len(missing_sample) >= 8:
                            break
        loss_fraction = (
            0.0
            if planned_total <= 0
            else float(planned_total - resolved_total) / float(planned_total)
        )
        if empty_chunks > 0 or loss_fraction > allowed_loss:
            raise RuntimeError(
                "chunk_resolution_preflight_failed:"
                f" run_mode={run_mode} library={library_name or 'unknown'}"
                f" root={library_root} sampled_chunks={len(sampled)}"
                f" planned={planned_total} resolved={resolved_total}"
                f" empty_chunks={empty_chunks} loss_fraction={loss_fraction:.3f}"
                f" sample_missing={','.join(missing_sample[:8]) or 'none'}"
            )
        logger.info(
            "[distributed.chunk.preflight] run_mode=%s library=%s root=%s sampled_chunks=%d planned=%d resolved=%d loss_fraction=%.3f",
            run_mode,
            library_name or "unknown",
            library_root,
            len(sampled),
            planned_total,
            resolved_total,
            loss_fraction,
        )


def _resolve_combo_ligand_bases(
    cfg: Mapping[str, Any],
    *,
    pdb_id: str,
    run_tokens: list[str],
    target_ligands: int,
    sample_seed: int,
) -> tuple[str, list[str]]:
    streams = _resolve_combo_ligand_streams(
        cfg,
        pdb_id=pdb_id,
        run_tokens=run_tokens,
        target_ligands=target_ligands,
        sample_seed=sample_seed,
    )
    if not streams:
        return str(cfg.get("LIBRARY_SUBDIR_DEFAULT", "fda_library")), []
    if len(streams) == 1:
        stream = streams[0]
        return str(stream.get("library_name") or ""), list(stream.get("ligand_bases") or [])
    library_name = "+".join(str(stream.get("library_name") or "") for stream in streams)
    ligands: list[str] = []
    for stream in streams:
        ligands.extend(str(x) for x in (stream.get("ligand_bases") or []) if str(x).strip())
    return library_name, sorted(set(ligands))


def _collect_combo_library_roots(
    cfg: Mapping[str, Any],
    *,
    combo_items: Sequence[tuple[str, Optional[str]]],
    run_tokens: list[str],
    logger: logging.Logger,
) -> list[Path]:
    seen: set[str] = set()
    out: list[Path] = []
    for pdb_file, _ph in combo_items:
        pdb_id = os.path.splitext(os.path.basename(str(pdb_file)))[0].upper()
        roots = compute_allowed_library_roots(
            dict(cfg), pdb_id, logger, tokens_override=run_tokens
        )
        for root in roots:
            try:
                key = str(root.resolve())
            except Exception:
                key = str(root)
            if key in seen:
                continue
            seen.add(key)
            out.append(Path(root))
    return out
