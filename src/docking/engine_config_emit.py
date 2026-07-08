from __future__ import annotations

import fcntl
import json
import os
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional

from path_router.path_router import (
    config_dir as router_config_dir,
    config_file as router_config_file,
    docked_dir as router_docked_dir,
    make_paths,
    receptor_file as router_receptor_file,
)


def _emit_log(
    logger: Optional[Any],
    level: str,
    message: str,
    *args: object,
) -> None:
    if logger:
        getattr(logger, level)(message, *args)
        return
    if args:
        print(message % args)
    else:
        print(message)


def _write_bytes_atomic(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path: Optional[str] = None
    try:
        with tempfile.NamedTemporaryFile(
            "wb",
            delete=False,
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".part",
        ) as handle:
            tmp_path = handle.name
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_path, path)
        tmp_path = None
    finally:
        if tmp_path:
            try:
                os.unlink(tmp_path)
            except FileNotFoundError:
                pass


def _build_vina_config_lines(
    *,
    receptor_for_config: str,
    ligand_path: str,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    cpu_per_job: int,
    stage_info: Dict[str, Any],
    out_path: Path,
    fast_mode: bool,
    flex_pdbqt: object,
) -> list[str]:
    exhaustiveness = int(stage_info.get("exhaustiveness", 8))
    num_modes = int(stage_info.get("num_modes", 4))
    if fast_mode:
        exhaustiveness = 1
        num_modes = 1

    lines = [
        f"receptor = {receptor_for_config}",
        f"ligand   = {ligand_path}",
        f"center_x = {center[0]:.3f}",
        f"center_y = {center[1]:.3f}",
        f"center_z = {center[2]:.3f}",
        f"size_x   = {box_size[0]:.3f}",
        f"size_y   = {box_size[1]:.3f}",
        f"size_z   = {box_size[2]:.3f}",
        f"cpu      = {int(cpu_per_job)}",
        f"exhaustiveness = {exhaustiveness}",
        f"energy_range   = {int(stage_info.get('energy_range', 4))}",
        f"num_modes      = {num_modes}",
        f"verbosity      = {int(stage_info.get('verbosity', 0))}",
        f"out = {out_path}",
    ]
    if flex_pdbqt:
        lines.append(f"flex = {flex_pdbqt}")
    if "seed" in stage_info:
        lines.append(f"seed = {int(stage_info['seed'])}")
    return lines


def _load_manifest_entries_map(
    manifest_path: Path,
) -> Optional[Dict[str, Dict[str, Any]]]:
    if not manifest_path.exists():
        return {}
    try:
        manifest_data = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return None

    entries = manifest_data.get("entries") if isinstance(manifest_data, dict) else None
    if not isinstance(entries, list):
        return {}

    entries_map: Dict[str, Dict[str, Any]] = {}
    for item in entries:
        if not isinstance(item, dict):
            continue
        lig = str(item.get("ligand", ""))
        if lig:
            entries_map[lig] = item
    return entries_map


def _write_manifest_payload(
    *,
    manifest_path: Path,
    run_id: str,
    pdb_id: str,
    stage_name: str,
    variant_token: Optional[str],
    ph_label: Optional[str],
    legacy_mode: bool,
    lig_base: str,
    cfg_path: Path,
    out_path: Path,
    receptor_for_config: str,
    flex_pdbqt: object,
    logger: Optional[Any] = None,
) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = manifest_path.with_name(f".{manifest_path.name}.lock")
    with open(lock_path, "a", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            entries_map = _load_manifest_entries_map(manifest_path)
            if entries_map is None:
                _emit_log(
                    logger,
                    "warning",
                    "[cfg.emit] manifest unreadable (corrupt JSON); skipping update path=%s",
                    str(manifest_path),
                )
                return
            entry: Dict[str, Any] = {
                "ligand": lig_base,
                "config": str(cfg_path),
                "out": str(out_path),
                "receptor": receptor_for_config,
            }
            if flex_pdbqt:
                entry["flex"] = str(flex_pdbqt)
            entries_map[lig_base] = entry
            manifest_data = {
                "run_id": run_id,
                "pdb_id": pdb_id,
                "stage": stage_name,
                "variant": variant_token,
                "ph": ph_label,
                "legacy": legacy_mode,
                "entries": [entries_map[key] for key in sorted(entries_map)],
            }
            payload = (
                json.dumps(manifest_data, indent=2, sort_keys=True) + "\n"
            ).encode("utf-8")
            _write_bytes_atomic(manifest_path, payload)
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def _optional_upper(value: object) -> Optional[str]:
    if value is None:
        return None
    return str(value).strip().upper() or None


def _optional_text(value: object) -> Optional[str]:
    if value is None:
        return None
    return str(value).strip() or None


def _manifest_path_for(
    *,
    run_id: str,
    pdb_id: str,
    stage_name: str,
    variant_token: Optional[str],
    ph_label: Optional[str],
    manifest_name: Optional[str],
    legacy_mode: bool,
) -> Path | None:
    if not manifest_name:
        return None
    return router_config_file(
        run_id,
        pdb_id,
        stage_name,
        variant=variant_token,
        ph_tag=ph_label,
        name=manifest_name,
        legacy=legacy_mode,
    )


def _receptor_for_stage(
    stage_info: Dict[str, Any], receptor_pdbqt: str, expected_receptor: Path
) -> str:
    override = (
        stage_info.get("receptor_pdbqt_override")
        or stage_info.get("receptor_override")
        or receptor_pdbqt
    )
    return str(override or expected_receptor)


def _log_missing_receptor(
    *,
    logger: Optional[Any],
    receptor_exists: bool,
    pdb_id: str,
    variant_display: str,
    ph_display: str,
    expected_receptor: Path,
) -> None:
    if receptor_exists:
        return
    _emit_log(
        logger,
        "error",
        "[router.error] missing receptor for pdb=%s variant=%s ph=%s -> %s",
        pdb_id,
        variant_display,
        ph_display,
        str(expected_receptor),
    )


def _log_vina_config_summary(
    *,
    logger: Optional[Any],
    engine_log_prefix: str,
    ligand_path: str,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
) -> None:
    if not logger:
        return
    cx, cy, cz = center
    sx, sy, sz = box_size
    logger.info(
        "[%s.cfg] lig=%s center=(%.3f,%.3f,%.3f) size=(%.1f,%.1f,%.1f)",
        engine_log_prefix,
        os.path.basename(str(ligand_path)),
        cx,
        cy,
        cz,
        sx,
        sy,
        sz,
    )


def _write_manifest_safely(
    *,
    manifest_path: Path | None,
    logger: Optional[Any],
    run_id: str,
    pdb_id: str,
    stage_name: str,
    variant_token: Optional[str],
    ph_label: Optional[str],
    legacy_mode: bool,
    lig_base: str,
    cfg_path: Path,
    out_path: Path,
    receptor_for_config: str,
    flex_pdbqt: object,
) -> None:
    if manifest_path is None:
        return
    try:
        _write_manifest_payload(
            manifest_path=manifest_path,
            run_id=run_id,
            pdb_id=pdb_id,
            stage_name=stage_name,
            variant_token=variant_token,
            ph_label=ph_label,
            legacy_mode=legacy_mode,
            lig_base=lig_base,
            cfg_path=cfg_path,
            out_path=out_path,
            receptor_for_config=receptor_for_config,
            flex_pdbqt=flex_pdbqt,
            logger=logger,
        )
    except FileNotFoundError:
        _emit_log(
            logger,
            "warning",
            "[cfg.emit] manifest_tmp missing; skipping manifest update tmp=%s dest=%s stage=%s",
            str(manifest_path.with_suffix(".part")),
            str(manifest_path),
            stage_name,
        )
    except Exception as exc:
        _emit_log(
            logger,
            "warning",
            "[cfg.emit] manifest update failed; continuing without manifest stage=%s path=%s err=%s",
            stage_name,
            str(manifest_path),
            exc,
        )


def emit_engine_config(
    cfg: Dict[str, Any],
    pdb_id: str,
    receptor_pdbqt: str,
    center: tuple[float, float, float],
    box_size: tuple[float, float, float],
    ligand_path: str,
    stage_name: str,
    stage_info: Dict[str, Any],
    cpu_per_job: int,
    logger: Optional[Any] = None,
    *,
    variant: Optional[str] = None,
    ph_token: Optional[str] = None,
    legacy: bool = False,
    manifest_name: Optional[str],
    engine_log_prefix: str,
    fast_mode: bool = False,
) -> tuple[str, str]:
    make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")

    variant_token = _optional_upper(variant)
    ph_label = _optional_text(ph_token)
    legacy_mode = bool(legacy)

    lig_base = Path(ligand_path).stem
    run_id = cfg["RUN_ID"]

    cfg_dir = router_config_dir(
        run_id,
        pdb_id,
        stage_name,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    cfg_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = _manifest_path_for(
        run_id=run_id,
        pdb_id=pdb_id,
        stage_name=stage_name,
        variant_token=variant_token,
        ph_label=ph_label,
        manifest_name=manifest_name,
        legacy_mode=legacy_mode,
    )

    cfg_path = cfg_dir / f"{lig_base}_{stage_name}.txt"
    stage_root = router_docked_dir(
        pdb_id,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    out_dir = stage_root / stage_name
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{lig_base}_{stage_name}.pdbqt"

    expected_receptor = router_receptor_file(
        pdb_id,
        variant=variant_token,
        ph_tag=ph_label,
        legacy=legacy_mode,
    )
    receptor_for_config = _receptor_for_stage(
        stage_info, receptor_pdbqt, expected_receptor
    )
    receptor_exists = Path(receptor_for_config).exists()
    flex_pdbqt = stage_info.get("flex_pdbqt") or stage_info.get("flex")

    variant_display = variant_token or "None"
    ph_display = ph_label or "None"
    breadcrumb = (
        "[cfg.emit] run=%s pdb=%s stage=%s variant=%s ph=%s\n"
        "           cfg_dir=%s receptor=%s out_root=%s"
    )
    _emit_log(
        logger,
        "info",
        breadcrumb,
        run_id,
        pdb_id,
        stage_name,
        variant_display,
        ph_display,
        str(cfg_dir),
        receptor_for_config,
        str(stage_root),
    )

    _log_missing_receptor(
        logger=logger,
        receptor_exists=receptor_exists,
        pdb_id=pdb_id,
        variant_display=variant_display,
        ph_display=ph_display,
        expected_receptor=expected_receptor,
    )

    lines = _build_vina_config_lines(
        receptor_for_config=receptor_for_config,
        ligand_path=ligand_path,
        center=center,
        box_size=box_size,
        cpu_per_job=cpu_per_job,
        stage_info=stage_info,
        out_path=out_path,
        fast_mode=fast_mode,
        flex_pdbqt=flex_pdbqt,
    )

    _log_vina_config_summary(
        logger=logger,
        engine_log_prefix=engine_log_prefix,
        ligand_path=ligand_path,
        center=center,
        box_size=box_size,
    )

    payload = ("\n".join(lines)).encode("utf-8")
    overwrite = cfg_path.exists()
    _write_bytes_atomic(cfg_path, payload)

    _write_manifest_safely(
        manifest_path=manifest_path,
        logger=logger,
        run_id=run_id,
        pdb_id=pdb_id,
        stage_name=stage_name,
        variant_token=variant_token,
        ph_label=ph_label,
        legacy_mode=legacy_mode,
        lig_base=lig_base,
        cfg_path=cfg_path,
        out_path=out_path,
        receptor_for_config=receptor_for_config,
        flex_pdbqt=flex_pdbqt,
    )

    _emit_log(
        logger,
        "info",
        "[cfg.emit] run=%s pdb=%s variant=%s ph=%s stage=%s ligand=%s cfg_dir=%s docked_root=%s path=%s overwrite=%s bytes=%d",
        run_id,
        pdb_id,
        variant_display,
        ph_display,
        stage_name,
        lig_base,
        str(cfg_dir),
        str(stage_root),
        str(cfg_path),
        str(overwrite).lower(),
        len(payload),
    )
    return str(cfg_path), str(out_path)
