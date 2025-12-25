from __future__ import annotations

import logging
import math
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from path_router import ph_ensemble_dir
import prep_for_ledock


def _resolve_dms_exe(cfg: dict, logger: logging.Logger) -> str:
    """
    Resolve the `dms` executable path.

    Precedence:
      1) cfg["DMS_EXE"] if set
      2) cfg["_FILE_CFG"]["DMS_EXE"] if present
      3) hard-coded /home/michael/atlas/tools/dms/bin/dms
      4) default to "dms" on PATH
    """
    candidate = None
    if isinstance(cfg, dict):
        candidate = cfg.get("DMS_EXE") or None
        file_cfg = cfg.get("_FILE_CFG")
        if not candidate and isinstance(file_cfg, dict):
            candidate = file_cfg.get("DMS_EXE") or None
    if not candidate:
        candidate = "/home/michael/atlas/tools/dms/bin/dms"

    # If the preferred path is missing, fall back to PATH to avoid hard failures in dev sandboxes.
    try:
        if not Path(candidate).exists():
            candidate = "dms"
    except Exception:
        candidate = "dms"

    logger.info("[dock6.dms.exe] path=%s", candidate)
    return str(candidate)


def _resolve_sphgen_exe(cfg: dict, logger: logging.Logger) -> str:
    """
    Resolve the `sphgen` executable path.

    Precedence:
      1) cfg["SPHGEN_EXE"] if set
      2) cfg["_FILE_CFG"]["SPHGEN_EXE"] if present
      3) hard-coded /home/michael/atlas/tools/dock6/bin/sphgen
      4) default to "sphgen" on PATH
    """
    candidate = None
    if isinstance(cfg, dict):
        candidate = cfg.get("SPHGEN_EXE") or None
        file_cfg = cfg.get("_FILE_CFG")
        if not candidate and isinstance(file_cfg, dict):
            candidate = file_cfg.get("SPHGEN_EXE") or None
    if not candidate:
        candidate = "/home/michael/atlas/tools/dock6/bin/sphgen"

    # If the preferred path is missing, fall back to PATH to avoid hard failures in dev sandboxes.
    try:
        if not Path(candidate).exists():
            candidate = "sphgen"
    except Exception:
        candidate = "sphgen"

    logger.info("[dock6.sphgen.exe] path=%s", candidate)
    return str(candidate)


def _resolve_sphere_selector_exe(cfg: dict, logger: logging.Logger) -> str:
    candidate = None
    if isinstance(cfg, dict):
        candidate = cfg.get("SPHERE_SELECTOR_EXE") or None
        file_cfg = cfg.get("_FILE_CFG")
        if not candidate and isinstance(file_cfg, dict):
            candidate = file_cfg.get("SPHERE_SELECTOR_EXE") or None
    if not candidate:
        candidate = "/home/michael/atlas/tools/dock6/bin/sphere_selector"
    logger.info("[dock6.sphere_selector.exe] path=%s", candidate)
    return str(candidate)


def _resolve_showbox_exe(cfg: dict, logger: logging.Logger) -> str:
    candidate = None
    if isinstance(cfg, dict):
        candidate = cfg.get("SHOWBOX_EXE") or None
        file_cfg = cfg.get("_FILE_CFG")
        if not candidate and isinstance(file_cfg, dict):
            candidate = file_cfg.get("SHOWBOX_EXE") or None
    if not candidate:
        candidate = "/home/michael/atlas/tools/dock6/bin/showbox"
    logger.info("[dock6.showbox.exe] path=%s", candidate)
    return str(candidate)


def _resolve_grid_exe(cfg: dict, logger: logging.Logger) -> str:
    candidate = None
    if isinstance(cfg, dict):
        candidate = cfg.get("GRID_EXE") or None
        file_cfg = cfg.get("_FILE_CFG")
        if not candidate and isinstance(file_cfg, dict):
            candidate = file_cfg.get("GRID_EXE") or None
    if not candidate:
        candidate = "/home/michael/atlas/tools/dock6/bin/grid"

    try:
        if not Path(candidate).exists():
            candidate = "grid"
    except Exception:
        candidate = "grid"

    logger.info("[dock6.grid.exe] path=%s", candidate)
    return str(candidate)


def _resolve_vdw_defn_path(cfg: dict, logger: logging.Logger) -> str:
    candidate = None
    if isinstance(cfg, dict):
        candidate = cfg.get("DOCK6_VDW_DEFN_FILE") or None
        file_cfg = cfg.get("_FILE_CFG")
        if not candidate and isinstance(file_cfg, dict):
            candidate = file_cfg.get("DOCK6_VDW_DEFN_FILE") or None
    if not candidate:
        candidate = "/home/michael/atlas/tools/dock6/parameters/vdw_AMBER_parm99.defn"

    try:
        if not Path(candidate).exists():
            logger.warning("[dock6.vdw_defn.missing] path=%s", candidate)
    except Exception:
        logger.warning("[dock6.vdw_defn.missing] path=%s", candidate)

    logger.info("[dock6.vdw_defn] path=%s", candidate)
    return str(candidate)


def _radius_from_box_size(size: tuple[float, float, float]) -> float:
    """
    Convert full box lengths (size_x/size_y/size_z) into a sphere radius.
    Vina-style sizes are full edge lengths, so radius is half the diagonal.
    """
    sx, sy, sz = size
    return 0.5 * math.sqrt(sx * sx + sy * sy + sz * sz)


def strip_hydrogens(in_pdb: Path, out_pdb: Path, logger: logging.Logger) -> None:
    """
    Strip hydrogens from ATOM/HETATM records, preserving other lines.
    Uses both element column and atom name as fallbacks.
    """
    lines_out: list[str] = []
    with in_pdb.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                lines_out.append(line)
                continue

            atomname = line[12:16].strip()
            elem = line[76:78].strip()

            if elem == "H":
                continue
            if not elem and atomname.startswith("H"):
                continue

            lines_out.append(line)

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    with out_pdb.open("w", encoding="utf-8") as out:
        out.writelines(lines_out)

    logger.info("[dock6.noH] in=%s out=%s", in_pdb, out_pdb)


def _ensure_sphgen_spheres(
    cfg: dict,
    dock6_root: Path,
    logger: logging.Logger,
) -> Optional[Path]:
    """
    Ensure DOCK6 sphgen spheres are generated in `dock6_root`.

    Requires:
    - rec.ms in dock6_root (created earlier by dms)

    Produces:
    - INSPH       (control file for sphgen)
    - rec.sph     (sphere file)
    - OUTSPH      (sphgen log)

    Returns:
    - Path to rec.sph on success
    - None on failure
    """
    ms_path = dock6_root / "rec.ms"
    if not ms_path.exists() or ms_path.stat().st_size == 0:
        logger.warning("[dock6.sphgen.skip] reason=missing_rec_ms dir=%s", dock6_root)
        return None

    insph_path = dock6_root / "INSPH"
    rec_sph_path = dock6_root / "rec.sph"
    outsph_path = dock6_root / "OUTSPH"

    # Idempotent: if rec.sph already exists and is non-empty, reuse it.
    if rec_sph_path.exists() and rec_sph_path.stat().st_size > 0:
        if outsph_path.exists() and outsph_path.stat().st_size > 0:
            logger.info(
                "[dock6.sphgen.reuse] rec_sph=%s outsph=%s",
                rec_sph_path,
                outsph_path,
            )
            return rec_sph_path

    insph_contents = "\n".join(
        [
            "rec.ms",
            "R",
            "X",
            "0.0",
            "4.0",
            "1.4",
            "rec.sph",
        ]
    ) + "\n"

    dock6_root.mkdir(parents=True, exist_ok=True)
    insph_path.write_text(insph_contents, encoding="utf-8")
    logger.info("[dock6.sphgen.insp] path=%s", insph_path)

    sphgen_exe = _resolve_sphgen_exe(cfg, logger)

    try:
        # sphgen reads INSPH from the current directory; we only need to set cwd.
        subprocess.run(
            [sphgen_exe],
            check=True,
            cwd=str(dock6_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except Exception as exc:
        logger.error(
            "[dock6.sphgen.error] dir=%s exe=%s reason=%s",
            dock6_root,
            sphgen_exe,
            exc,
        )
        return None

    if not rec_sph_path.exists() or rec_sph_path.stat().st_size == 0:
        logger.error("[dock6.sphgen.empty] rec_sph=%s", rec_sph_path)
        return None

    if not outsph_path.exists() or outsph_path.stat().st_size == 0:
        logger.warning(
            "[dock6.sphgen.no_outsph] OUTSPH missing or empty at %s",
            outsph_path,
        )

    logger.info("[dock6.sphgen.ready] rec_sph=%s outsph=%s", rec_sph_path, outsph_path)
    return rec_sph_path


def _ensure_selected_spheres(
    cfg: dict,
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    dock6_root: Path,
    logger: logging.Logger,
) -> Optional[Path]:
    """
    Ensure active-site spheres are selected around the pocket center.

    Requires:
      - rec.sph present in dock6_root
      - Active-site center/box available from docking logic
    Produces:
      - site_center.mol2 (one atom at the pocket center)
      - selected_spheres.sph (subset around the center)
    """
    rec_sph_path = dock6_root / "rec.sph"
    if not rec_sph_path.exists() or rec_sph_path.stat().st_size == 0:
        logger.warning("[dock6.selector.skip] reason=missing_rec_sph dir=%s", dock6_root)
        return None

    selected_path = dock6_root / "selected_spheres.sph"
    site_center_path = dock6_root / "site_center.mol2"

    if (
        selected_path.exists()
        and selected_path.stat().st_size > 0
        and site_center_path.exists()
        and site_center_path.stat().st_size > 0
    ):
        logger.info("[dock6.selector.rebuild] selected_existing=%s", selected_path)

    try:
        from docking import get_active_site_center_and_size
    except Exception as exc:
        logger.error("[dock6.selector.error] reason=import_failed err=%s", exc)
        return None

    cs = get_active_site_center_and_size(cfg, pdb_id, variant, ph_label, logger)
    if cs is None:
        logger.warning(
            "[dock6.selector.skip] reason=no_active_site pdb=%s variant=%s ph=%s",
            pdb_id,
            variant,
            ph_label,
        )
        return None

    center, size = cs
    cx, cy, cz = center
    radius = _radius_from_box_size(size)

    site_center_contents = (
        "@<TRIPOS>MOLECULE\n"
        "site_center\n"
        " 1 0 0 0 0\n"
        "SMALL\n"
        "NO_CHARGES\n"
        "\n"
        "@<TRIPOS>ATOM\n"
        f"{1:7d} C1{cx:11.4f}{cy:11.4f}{cz:11.4f} C.3{1:8d} SITE{0.0:11.4f}\n"
    )

    dock6_root.mkdir(parents=True, exist_ok=True)
    site_center_path.write_text(site_center_contents, encoding="utf-8")
    logger.info(
        "[dock6.selector.site_center] path=%s center=(%.3f, %.3f, %.3f) radius=%.3f",
        site_center_path,
        cx,
        cy,
        cz,
        radius,
    )

    sphere_selector_exe = _resolve_sphere_selector_exe(cfg, logger)

    try:
        subprocess.run(
            [sphere_selector_exe, str(rec_sph_path), str(site_center_path), f"{radius:.3f}"],
            check=True,
            cwd=str(dock6_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except Exception as exc:
        logger.error(
            "[dock6.selector.error] pdb=%s variant=%s ph=%s exe=%s reason=%s",
            pdb_id,
            variant,
            ph_label,
            sphere_selector_exe,
            exc,
        )
        return None

    if not selected_path.exists() or selected_path.stat().st_size == 0:
        logger.error("[dock6.selector.empty] path=%s", selected_path)
        return None

    logger.info("[dock6.selector.ready] selected=%s", selected_path)
    return selected_path


def _ensure_site_box(
    cfg: dict,
    dock6_root: Path,
    logger: logging.Logger,
) -> Optional[Path]:
    """
    Ensure DOCK6 site box (site_box.pdb) exists for the selected spheres.
    """
    selected_path = dock6_root / "selected_spheres.sph"
    box_path = dock6_root / "site_box.pdb"

    if not selected_path.exists() or selected_path.stat().st_size == 0:
        logger.warning("[dock6.showbox.skip] reason=missing_selected_spheres dir=%s", dock6_root)
        return None

    showbox_exe = _resolve_showbox_exe(cfg, logger)
    showbox_input = "\n".join(
        [
            "Y",  # automatically construct box
            "5.0",  # margin in Å
            "selected_spheres.sph",
            "1",  # cluster number
            "site_box.pdb",
        ]
    ) + "\n"

    try:
        subprocess.run(
            [showbox_exe],
            input=showbox_input.encode("ascii"),
            check=True,
            cwd=str(dock6_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except Exception as exc:
        logger.error(
            "[dock6.showbox.error] dir=%s exe=%s reason=%s",
            dock6_root,
            showbox_exe,
            exc,
        )
        return None

    if not box_path.exists() or box_path.stat().st_size == 0:
        logger.error("[dock6.showbox.empty] box=%s", box_path)
        return None

    logger.info("[dock6.showbox.ready] box=%s", box_path)
    return box_path


def _center_from_pdb_atoms(pdb_path: Path) -> Optional[tuple[float, float, float]]:
    coords: list[tuple[float, float, float]] = []
    try:
        with pdb_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                coords.append((x, y, z))
    except FileNotFoundError:
        return None
    except Exception:
        return None

    if not coords:
        return None

    xs, ys, zs = zip(*coords)
    return (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))


def _log_site_alignment(
    cfg: dict,
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    dock6_root: Path,
    logger: logging.Logger,
) -> None:
    try:
        from docking import get_active_site_center_and_size
    except Exception:
        return

    cs = get_active_site_center_and_size(cfg, pdb_id, variant, ph_label, logger)
    if cs is None:
        return
    center, _ = cs
    box_center = _center_from_pdb_atoms(dock6_root / "site_box.pdb")
    if not box_center:
        return
    delta = math.dist(center, box_center)
    logger.info(
        "[dock6.site.check] pdb=%s variant=%s ph=%s delta_center=%.3f",
        pdb_id,
        variant or "HOLO",
        ph_label or "base",
        delta,
    )


def _write_grid_in(
    grid_in_path: Path,
    receptor_pdb: Path,
    vdw_defn: str,
    grid_prefix: str,
    grid_spacing: float,
) -> None:
    rec_name = receptor_pdb.name if receptor_pdb.parent == grid_in_path.parent else str(receptor_pdb)
    grid_contents = "\n".join(
        [
            "compute_grids yes",
            f"grid_spacing {grid_spacing}",
            "output_molecule yes",
            "",
            "contact_score no",
            "chemical_score no",
            "energy_score yes",
            "energy_cutoff_distance 999",
            "",
            "atom_model a",
            "attractive_exponent 6",
            "repulsive_exponent 9",
            "distance_dielectric yes",
            "dielectric_factor 4",
            "bump_filter yes",
            "bump_overlap 0.75",
            "",
            "allow_non_integral_charges yes",
            "",
            f"receptor_file {rec_name}",
            "box_file site_box.pdb",
            f"vdw_definition_file {vdw_defn}",
            "",
            f"score_grid_prefix {grid_prefix}",
            f"receptor_out_file {receptor_pdb.stem}_{grid_prefix}.grid.pdb",
        ]
    ) + "\n"
    grid_in_path.write_text(grid_contents, encoding="utf-8")


def _run_grid(grid_exe: str, grid_in_path: Path, grid_out_path: Path, dock6_root: Path, logger: logging.Logger) -> bool:
    try:
        subprocess.run(
            [grid_exe, "-i", grid_in_path.name, "-o", grid_out_path.name],
            check=True,
            cwd=str(dock6_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return True
    except Exception as exc:
        logger.error(
            "[dock6.grid.error] dir=%s exe=%s reason=%s",
            dock6_root,
            grid_exe,
            exc,
        )
        return False


def _ensure_dock6_grids(
    cfg: dict,
    dock6_root: Path,
    receptor_pdb: Path,
    logger: logging.Logger,
    *,
    grid_prefix: str,
    grid_spacing: float,
) -> Optional[Path]:
    """
    Ensure DOCK6 grid files exist in dock6_root for the given receptor and site_box.
    """
    box_path = dock6_root / "site_box.pdb"
    if not box_path.exists() or box_path.stat().st_size == 0:
        logger.warning("[dock6.grid.skip] reason=missing_site_box dir=%s", dock6_root)
        return None

    if not receptor_pdb.exists() or receptor_pdb.stat().st_size == 0:
        logger.warning(
            "[dock6.grid.skip] reason=missing_receptor dir=%s receptor=%s",
            dock6_root,
            receptor_pdb,
        )
        return None

    grid_nrg = dock6_root / f"{grid_prefix}.nrg"
    grid_bmp = dock6_root / f"{grid_prefix}.bmp"
    grid_in_path = dock6_root / f"{grid_prefix}.in"
    grid_out_path = dock6_root / f"{grid_prefix}.out"

    if grid_nrg.exists() and grid_nrg.stat().st_size > 0 and grid_bmp.exists() and grid_bmp.stat().st_size > 0:
        logger.info("[dock6.grid.reuse] dir=%s prefix=%s nrg=%s bmp=%s", dock6_root, grid_prefix, grid_nrg, grid_bmp)
        return grid_nrg

    vdw_defn = _resolve_vdw_defn_path(cfg, logger)
    dock6_root.mkdir(parents=True, exist_ok=True)
    _write_grid_in(grid_in_path, receptor_pdb, vdw_defn, grid_prefix, grid_spacing)
    logger.info("[dock6.grid.in] path=%s", grid_in_path)

    grid_exe = _resolve_grid_exe(cfg, logger)
    logger.info("[DOCK6_GRID] building %s spacing=%s", grid_prefix, grid_spacing)
    ok = _run_grid(grid_exe, grid_in_path, grid_out_path, dock6_root, logger)
    if not ok:
        return None

    if not grid_nrg.exists() or grid_nrg.stat().st_size == 0:
        logger.error("[dock6.grid.empty] nrg=%s", grid_nrg)
        return None
    if not grid_bmp.exists() or grid_bmp.stat().st_size == 0:
        logger.warning("[dock6.grid.no_bmp] bmp missing or empty at %s", grid_bmp)

    logger.info("[DOCK6_GRID] ready %s", grid_prefix)
    return grid_nrg



def ensure_dock6_surface(
    cfg: dict,
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    logger: logging.Logger,
) -> Optional[Path]:
    """
    Ensure that a DOCK6-ready surface exists for the given target.

    Steps:
      - Resolve variant/pH tokens using the same logic as LeDock prep.
      - Locate the withH receptor in the ph_ensemble folder.
      - Strip hydrogens to <base>_noH.pdb under dock6/.
      - Run dms to produce rec.ms under dock6/.
    """
    legacy_mode = bool(cfg.get("_ROUTER_LEGACY", False))
    variant_token = (str(variant).strip().upper() or None) if variant is not None else None
    variant_for_ph = prep_for_ledock._variant_for_ph(variant_token, legacy_mode)
    ph_token = prep_for_ledock._normalize_ph_label(ph_label)
    if not ph_token:
        logger.warning("[dock6.surface.skip] reason=missing_ph_label pdb=%s", pdb_id)
        return None

    if bool(cfg.get("FAST_MODE")):
        grid_plan = [("grid_fast", 0.60)]
    else:
        grid_plan = [("grid_s1", 0.50), ("grid_s2", 0.35), ("grid_s3", 0.30)]

    ensemble_dir = ph_ensemble_dir(pdb_id, variant=variant_for_ph, legacy=legacy_mode)
    withh_path = prep_for_ledock._resolve_withh_from_manifest(
        ensemble_dir, pdb_id, ph_token, logger
    )
    if withh_path is None:
        prefix = f"{str(pdb_id).upper()}_"
        fallback_label = ph_token
        if not fallback_label.startswith(prefix):
            fallback_label = f"{prefix}{fallback_label}"
        candidate = ensemble_dir / f"{fallback_label}.withH.pdb"
        if candidate.exists():
            withh_path = candidate
        else:
            logger.warning(
                "[dock6.surface.skip] reason=withH_missing pdb=%s ph=%s dir=%s",
                pdb_id,
                ph_token,
                ensemble_dir,
            )
            return None

    dock6_root = ensemble_dir / "dock6"
    dock6_root.mkdir(parents=True, exist_ok=True)

    base = withh_path.name
    if base.endswith(".withH.pdb"):
        base = base[: -len(".withH.pdb")]
    elif base.endswith(".pdb"):
        base = base[: -len(".pdb")]

    noH_path = dock6_root / f"{base}_noH.pdb"
    ms_path = dock6_root / "rec.ms"

    def _build_grids_once():
        if not grid_plan:
            return
        with ThreadPoolExecutor(max_workers=len(grid_plan)) as pool:
            futures = {}
            for prefix, spacing in grid_plan:
                futures[
                    pool.submit(
                        _ensure_dock6_grids,
                        cfg=cfg,
                        dock6_root=dock6_root,
                        receptor_pdb=noH_path,
                        logger=logger,
                        grid_prefix=prefix,
                        grid_spacing=spacing,
                    )
                ] = (prefix, spacing)
            for fut in as_completed(futures):
                prefix, _spacing = futures[fut]
                res = None
                try:
                    res = fut.result()
                except Exception as exc:
                    logger.error("[dock6.grid.error] prefix=%s dir=%s reason=%s", prefix, dock6_root, exc)
                    raise
                if res is None:
                    raise RuntimeError(f"dock6 grid failed for prefix {prefix}")

    if (
        noH_path.exists()
        and noH_path.stat().st_size > 0
        and ms_path.exists()
        and ms_path.stat().st_size > 0
    ):
        logger.info(
            "[dock6.surface.reuse] pdb=%s variant=%s ph=%s noH=%s ms=%s",
            pdb_id,
            variant_for_ph or "HOLO",
            ph_token,
            noH_path,
            ms_path,
        )
        try:
            _ensure_sphgen_spheres(cfg=cfg, dock6_root=dock6_root, logger=logger)
        except Exception as exc:
            logger.error(
                "[dock6.surface.sphgen.error] pdb=%s variant=%s ph=%s dir=%s reason=%s",
                pdb_id,
                variant_for_ph or "HOLO",
                ph_token,
                dock6_root,
                exc,
            )
        try:
            _ensure_selected_spheres(
                cfg=cfg,
                pdb_id=pdb_id,
                variant=variant_for_ph,
                ph_label=ph_label,
                dock6_root=dock6_root,
                logger=logger,
            )
            _ensure_site_box(cfg=cfg, dock6_root=dock6_root, logger=logger)
            _log_site_alignment(
                cfg=cfg,
                pdb_id=pdb_id,
                variant=variant_for_ph,
                ph_label=ph_label,
                dock6_root=dock6_root,
                logger=logger,
            )
            try:
                _build_grids_once()
            except Exception as exc:
                logger.error(
                    "[dock6.surface.grid.error] pdb=%s variant=%s ph=%s dir=%s reason=%s",
                    pdb_id,
                    variant_for_ph or "HOLO",
                    ph_token,
                    dock6_root,
                    exc,
                )
                return None
        except Exception as exc:
            logger.error(
                "[dock6.surface.site.error] pdb=%s variant=%s ph=%s dir=%s reason=%s",
                pdb_id,
                variant_for_ph or "HOLO",
                ph_token,
                dock6_root,
                exc,
            )
        return ms_path

    strip_hydrogens(withh_path, noH_path, logger=logger)

    from pathlib import Path

    dms_exe = _resolve_dms_exe(cfg, logger)

    # Prefer running dms from its install root so it can find lib/dms/radii and lib/dms/dmsd
    dms_cwd = str(dock6_root)
    try:
        dms_path = Path(dms_exe)
        if dms_path.is_absolute() and dms_path.exists():
            dms_root = dms_path.resolve().parent.parent  # .../dms
            radii = dms_root / "lib" / "dms" / "radii"
            dmsd = dms_root / "lib" / "dms" / "dmsd"
            if radii.exists() and dmsd.exists():
                dms_cwd = str(dms_root)
    except Exception:
        pass

    try:
        subprocess.run(
            [dms_exe, str(noH_path), "-n", "-w", "1.4", "-v", "-o", str(ms_path)],
            check=True,
            cwd=dms_cwd,
        )
    except Exception as exc:
        logger.error(
            "[dock6.surface.dms.error] pdb=%s variant=%s ph=%s exe=%s cwd=%s reason=%s",
            pdb_id,
            variant_for_ph or "HOLO",
            ph_token,
            dms_exe,
            dms_cwd,
            exc,
        )
        return None

    if not ms_path.exists() or ms_path.stat().st_size == 0:
        logger.error(
            "[dock6.surface.dms.empty] pdb=%s variant=%s ph=%s path=%s",
            pdb_id,
            variant_for_ph or "HOLO",
            ph_token,
            ms_path,
        )
        return None

    try:
        _ensure_sphgen_spheres(cfg=cfg, dock6_root=dock6_root, logger=logger)
    except Exception as exc:
        logger.error(
            "[dock6.surface.sphgen.error] pdb=%s variant=%s ph=%s dir=%s reason=%s",
            pdb_id,
            variant_for_ph or "HOLO",
            ph_token,
            dock6_root,
            exc,
        )

    try:
        _ensure_selected_spheres(
            cfg=cfg,
            pdb_id=pdb_id,
            variant=variant_for_ph,
            ph_label=ph_label,
            dock6_root=dock6_root,
            logger=logger,
        )
        _ensure_site_box(cfg=cfg, dock6_root=dock6_root, logger=logger)
        _log_site_alignment(
            cfg=cfg,
            pdb_id=pdb_id,
            variant=variant_for_ph,
            ph_label=ph_label,
            dock6_root=dock6_root,
            logger=logger,
        )
        try:
            _build_grids_once()
        except Exception as exc:
            logger.error(
                "[dock6.surface.grid.error] pdb=%s variant=%s ph=%s dir=%s reason=%s",
                pdb_id,
                variant_for_ph or "HOLO",
                ph_token,
                dock6_root,
                exc,
            )
            return None
    except Exception as exc:
        logger.error(
            "[dock6.surface.site.error] pdb=%s variant=%s ph=%s dir=%s reason=%s",
            pdb_id,
            variant_for_ph or "HOLO",
            ph_token,
            dock6_root,
            exc,
        )

    logger.info(
        "[dock6.surface.ready] pdb=%s variant=%s ph=%s noH=%s ms=%s",
        pdb_id,
        variant_for_ph or "HOLO",
        ph_token,
        noH_path,
        ms_path,
    )
    return ms_path


def ensure_dock6_site(
    cfg: dict,
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    logger: logging.Logger,
) -> Optional[Path]:
    """
    Ensure full DOCK6 site prep (surface + spheres + box) exists.
    """
    return ensure_dock6_surface(cfg, pdb_id, variant, ph_label, logger)
