from __future__ import annotations

import os
import logging
from pathlib import Path
from typing import Any, Dict, Optional

from docking.engine_config_emit import emit_engine_config
from path_router.path_router import make_paths
from docking.score_io import write_score_summary_to_csv
from docking.pose_validation import compute_self_rmsd

# Engine-specific home for Vina config emission and CSV writing.


def _ph_variant_for_paths(
    ph_label: Optional[str], variant: Optional[str]
) -> tuple[Optional[str], Optional[str]]:
    """``(ph_token, variant_token)`` for routed docked paths (APO_HOLO from env when needed)."""
    ph_token = (ph_label or "").strip() or None
    variant_env = (
        (variant or os.environ.get("APO_HOLO_VARIANT", "") or "").strip().upper()
    )
    variant_token = variant_env or None
    return ph_token, variant_token


def emit_vina_config(
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
    skip_manifest_if_exists: bool = False,
):
    result = emit_engine_config(
        cfg,
        pdb_id,
        receptor_pdbqt,
        center,
        box_size,
        ligand_path,
        stage_name,
        stage_info,
        cpu_per_job,
        logger,
        variant=variant,
        ph_token=ph_token,
        legacy=legacy,
        manifest_name=None if skip_manifest_if_exists else "vina.json",
        engine_log_prefix="vina",
        fast_mode=bool(cfg.get("FAST_MODE")),
    )
    _log_cfg_emit_path_check(pdb_id, receptor_pdbqt, variant, legacy)
    return result


def _log_cfg_emit_path_check(
    pdb_id: str,
    receptor_path: str,
    variant: Optional[str],
    legacy: bool,
) -> None:
    variant_token = (
        (str(variant).strip().upper() or None) if variant is not None else None
    )
    variant_label = variant_token or "legacy"
    contains_variant = bool(variant_token and variant_token in str(receptor_path))
    logging.info(
        "[cfg.emit.check] pdb=%s variant=%s receptor_path=%s path_contains_variant=%s",
        pdb_id,
        variant_label,
        receptor_path,
        contains_variant,
    )
    if variant_token and not contains_variant and not legacy:
        logging.warning(
            "[variant.mismatch] expected_variant=%s wrote_legacy_path=%s action=fail_ci",
            variant_token,
            receptor_path,
        )


def _pose_path_for(
    csv_cfg: Dict,
    pdb_id: str,
    stage_name: str,
    lig_path: str,
    ph_label: Optional[str] = None,
    variant: Optional[str] = None,
) -> str:
    """Build the expected pose path for a ligand at a given stage."""
    paths = make_paths(csv_cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_token, variant_token = _ph_variant_for_paths(ph_label, variant)
    stage_dir = paths.docked_stage_dir(variant_token, stage_name, ph_token)
    return str(stage_dir / f"{Path(lig_path).stem}_{stage_name}.pdbqt")


def write_scores_csv(
    cfg: Dict,
    pdb_id: str,
    score_history: Dict[str, Dict[str, Dict]],
    ph_label: Optional[str] = None,
    variant: Optional[str] = None,
    *,
    csv_prefix: str = "",
    summary_basename: str = "docking_score_summary.csv",
    long_basename: str = "docking_score_long.csv",
) -> str:
    import csv as _csv
    import errno as _errno
    import math as _math
    import time as _time

    paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    ph_token, variant_token = _ph_variant_for_paths(ph_label, variant)
    variant_env = variant_token or ""
    dock_dir = paths.docked_variant_root(variant_token, ph_token)
    dock_dir.mkdir(parents=True, exist_ok=True)

    run_id_value = str(cfg.get("RUN_ID") or "")
    variant_value = variant_env
    include_variant = bool(variant_value)

    long_name = f"{csv_prefix}{long_basename}"
    summary_name = f"{csv_prefix}{summary_basename}"
    csv_out_wide_path = dock_dir / summary_name
    csv_out_wide = str(csv_out_wide_path)
    csv_out_long_path = dock_dir / long_name

    flat: Dict[str, Dict[str, Any]] = {}
    for stage_name, stage_map in score_history.items():
        flat[stage_name] = {}
        for lig, rec in stage_map.items():
            lig_key = os.path.basename(lig)
            s = rec.get("score", None)
            if rec.get("valid", False):
                flat[stage_name][lig_key] = s if s is not None else ""
            else:
                flat[stage_name][lig_key] = (
                    f"{s:.2f} (invalid)" if isinstance(s, (int, float)) else "(invalid)"
                )

    def _long_header() -> list[str]:
        header = ["run_id"]
        if include_variant:
            header.append("variant")
        header.extend(
            [
                "stage",
                "ligand",
                "score",
                "valid",
                "reason",
                "heavy_atoms",
                "le",
                "self_rmsd",
                "pains_flag",
            ]
        )
        return header

    def _build_long_row_map() -> dict[tuple[str, str], dict[str, str]]:
        out: dict[tuple[str, str], dict[str, str]] = {}
        for stage_name, stage_map in score_history.items():
            for lig, rec in stage_map.items():
                lig_key = os.path.basename(lig)
                score = rec.get("score", None)
                valid = bool(rec.get("valid", False))
                reason = rec.get("reason", "")

                ha = rec.get("heavy_atoms", None)
                le = rec.get("le", None)
                pains_hit = rec.get("pains_flag", False)

                pose_path = _pose_path_for(
                    cfg,
                    pdb_id,
                    stage_name,
                    lig,
                    ph_label=ph_token,
                    variant=variant_token,
                )
                if os.path.exists(pose_path):
                    try:
                        sr = compute_self_rmsd(pose_path)
                        sr_str = (
                            f"{sr:.2f}"
                            if isinstance(sr, (int, float)) and _math.isfinite(sr)
                            else ""
                        )
                    except Exception:
                        sr_str = ""
                else:
                    sr_str = ""

                score_str = f"{score:.2f}" if isinstance(score, (int, float)) else ""
                ha_str = str(int(ha)) if isinstance(ha, (int, float)) else ""
                le_str = f"{le:.4f}" if isinstance(le, (int, float)) else ""
                reason_str = str(reason) if reason is not None else ""

                row = {
                    "run_id": run_id_value,
                    "stage": str(stage_name),
                    "ligand": str(lig_key),
                    "score": score_str,
                    "valid": str(int(valid)),
                    "reason": reason_str,
                    "heavy_atoms": ha_str,
                    "le": le_str,
                    "self_rmsd": sr_str,
                    "pains_flag": str(int(bool(pains_hit))),
                }
                if include_variant:
                    row["variant"] = variant_value
                out[(str(stage_name), str(lig_key))] = row
        return out

    def _load_summary_flat(path: Path) -> dict[str, dict[str, str]]:
        loaded: dict[str, dict[str, str]] = {}
        if not path.exists():
            return loaded
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = _csv.DictReader(handle)
            fieldnames = reader.fieldnames or []
            ligand_col = "Ligand" if "Ligand" in fieldnames else "ligand"
            stage_cols = [c for c in fieldnames if c.lower().startswith("stage")]
            for row in reader:
                ligand = str(row.get(ligand_col, "")).strip()
                if not ligand:
                    continue
                for stage_col in stage_cols:
                    cell = row.get(stage_col, "")
                    if cell is None or str(cell).strip() == "":
                        continue
                    loaded.setdefault(stage_col, {})[ligand] = str(cell)
        return loaded

    def _load_long_rows(
        path: Path, header: list[str]
    ) -> dict[tuple[str, str], dict[str, str]]:
        rows: dict[tuple[str, str], dict[str, str]] = {}
        if not path.exists():
            return rows
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = _csv.DictReader(handle)
            for raw in reader:
                stage = str(raw.get("stage", "")).strip()
                ligand = str(raw.get("ligand", "")).strip()
                if not stage or not ligand:
                    continue
                normalized = {col: str(raw.get(col, "")) for col in header}
                rows[(stage, ligand)] = normalized
        return rows

    def _write_long_rows_atomic(
        path: Path,
        header: list[str],
        rows: dict[tuple[str, str], dict[str, str]],
    ) -> None:
        tmp_path = path.with_suffix(path.suffix + ".tmp")
        with tmp_path.open("w", newline="", encoding="utf-8") as handle:
            writer = _csv.DictWriter(handle, fieldnames=header)
            writer.writeheader()
            for key in sorted(rows.keys()):
                writer.writerow(rows[key])
        os.replace(tmp_path, path)

    def _chunk_merge_enabled() -> bool:
        raw = cfg.get("_CHUNK_LIGAND_KEYS")
        return isinstance(raw, (list, tuple, set)) and any(str(x).strip() for x in raw)

    def _acquire_lock(lock_path: Path) -> int:
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        start = _time.time()
        timeout_sec = 300.0
        stale_sec = 1800.0
        while True:
            try:
                fd = os.open(
                    str(lock_path),
                    os.O_CREAT | os.O_EXCL | os.O_WRONLY,
                    0o644,
                )
                try:
                    os.write(
                        fd,
                        f"{os.getpid()} {_time.time():.3f}\n".encode("utf-8"),
                    )
                except Exception:
                    pass
                return fd
            except OSError as exc:
                if exc.errno != _errno.EEXIST:
                    raise
                try:
                    age = _time.time() - lock_path.stat().st_mtime
                    if age > stale_sec:
                        lock_path.unlink()
                        continue
                except FileNotFoundError:
                    continue
                if (_time.time() - start) >= timeout_sec:
                    raise TimeoutError(
                        f"Timed out waiting for CSV lock: {lock_path}"
                    ) from exc
                _time.sleep(0.05)

    def _release_lock(fd: Optional[int], lock_path: Path) -> None:
        if fd is None:
            return
        try:
            os.close(fd)
        finally:
            try:
                lock_path.unlink()
            except FileNotFoundError:
                pass

    def _write_non_chunk() -> None:
        write_score_summary_to_csv(
            flat,
            output_path=csv_out_wide,
            run_id=run_id_value,
            variant=variant_value if include_variant else None,
        )
        header = _long_header()
        row_map = _build_long_row_map()
        _write_long_rows_atomic(csv_out_long_path, header, row_map)

    if not _chunk_merge_enabled():
        _write_non_chunk()
        return csv_out_wide

    lock_path = Path(f"{csv_out_wide}.lock")
    lock_fd: Optional[int] = None
    try:
        lock_fd = _acquire_lock(lock_path)

        merged_flat = _load_summary_flat(csv_out_wide_path)
        for stage_name, stage_map in flat.items():
            target = merged_flat.setdefault(str(stage_name), {})
            for ligand_name, value in stage_map.items():
                lig_key = str(ligand_name).strip()
                if lig_key:
                    target[lig_key] = "" if value is None else str(value)
        write_score_summary_to_csv(
            merged_flat,
            output_path=csv_out_wide,
            run_id=run_id_value,
            variant=variant_value if include_variant else None,
        )

        header = _long_header()
        merged_long_rows = _load_long_rows(csv_out_long_path, header)
        merged_long_rows.update(_build_long_row_map())
        _write_long_rows_atomic(csv_out_long_path, header, merged_long_rows)
    finally:
        _release_lock(lock_fd, lock_path)

    return csv_out_wide
