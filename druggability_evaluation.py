from __future__ import annotations

import logging
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from installation import load_config
from path_router import make_paths

# Load config once at import to mirror other helpers.
config = load_config()

# Defaults
_DEFAULT_FPOCKET_EXE = Path("/home/michael/atlas/tools/fpocket/bin/fpocket")
_DEFAULT_FPOCKET_OUTPUT_ROOT = Path(__file__).resolve().parent / "fpocket"


def _cfg_float(key: str, default: float) -> float:
    raw = config.get(key)
    try:
        return float(raw)
    except Exception:
        return default


@dataclass
class DruggabilityMetrics:
    pdb_id: str
    variant: Optional[str]
    ph_label: Optional[str]
    druggability: float
    volume: float
    openness: float
    polar_fraction: float
    has_metal: bool
    tier: str
    triggers: List[str]
    info_path: Path


def get_fpocket_exe(cfg: Dict[str, Any]) -> Optional[Path]:
    """Return a validated fpocket executable path if present; otherwise None."""
    raw = cfg.get("FPOCKET_EXE") or _DEFAULT_FPOCKET_EXE
    try:
        exe_path = Path(raw)
        if not exe_path.exists():
            logging.warning(
                "[druggability.fpocket.missing] exe=%s reason=missing", exe_path
            )
            return None
        if not os.access(exe_path, os.X_OK):
            logging.warning(
                "[druggability.fpocket.missing] exe=%s reason=not_executable", exe_path
            )
            return None
        return exe_path
    except Exception as exc:
        logging.warning("[druggability.fpocket.error] stage=exe_lookup err=%s", exc)
        return None


def get_fpocket_output_root(cfg: Dict[str, Any]) -> Path:
    """Resolve and ensure the consolidated fpocket output root."""
    raw = cfg.get("FPOCKET_OUTPUT_ROOT") or _DEFAULT_FPOCKET_OUTPUT_ROOT
    try:
        root = Path(raw)
        root.mkdir(parents=True, exist_ok=True)
        return root
    except Exception:
        logging.warning(
            "[druggability.fpocket.output_root.warn] root=%s reason=mkdir_failed",
            raw,
            exc_info=True,
        )
        root = Path(raw)
        return root


def _detect_metal_near_center(
    receptor_pdb: Union[str, Path],
    center: Tuple[float, float, float],
    radius: float = 8.0,
) -> bool:
    """
    Return True if a metal/cofactor atom is found within `radius` Å of the active-site center.
    """
    metals = {
        "ZN",
        "FE",
        "MG",
        "MN",
        "CO",
        "NI",
        "CU",
        "CA",
        "NA",
        "K",
    }
    try:
        path = Path(receptor_pdb)
        if not path.exists():
            return False
        r2 = float(radius) * float(radius)
        cx, cy, cz = (float(center[0]), float(center[1]), float(center[2]))
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except Exception:
                    continue
                dx = x - cx
                dy = y - cy
                dz = z - cz
                if (dx * dx + dy * dy + dz * dz) > r2:
                    continue
                elem = (line[76:78] or "").strip().upper()
                resname = (line[17:20] or "").strip().upper()
                if elem in metals or resname in metals:
                    return True
    except Exception as exc:
        logging.debug("[druggability.fpocket.metal.detect.error] path=%s err=%s", receptor_pdb, exc)
    return False


def _classify_druggability(
    druggability: float,
    volume: float,
    openness: float,
    polar_frac: float,
    has_metal: bool,
) -> Tuple[str, List[str]]:
    """
    Return (tier_label, triggers).
    """
    c_druggability_max = _cfg_float("DRUGGABILITY_TIER_C_DRUGGABILITY_MAX", 0.30)
    c_openness_min = _cfg_float("DRUGGABILITY_TIER_C_OPENNESS_MIN", 0.65)
    c_polar_frac_min = _cfg_float("DRUGGABILITY_TIER_C_POLAR_FRAC_MIN", 0.60)

    a_druggability_min = _cfg_float("DRUGGABILITY_TIER_A_DRUGGABILITY_MIN", 0.50)
    a_openness_max = _cfg_float("DRUGGABILITY_TIER_A_OPENNESS_MAX", 0.55)
    a_polar_frac_max = _cfg_float("DRUGGABILITY_TIER_A_POLAR_FRAC_MAX", 0.45)

    b_druggability_min = _cfg_float("DRUGGABILITY_TIER_B_DRUGGABILITY_MIN", 0.30)
    b_druggability_max = _cfg_float("DRUGGABILITY_TIER_B_DRUGGABILITY_MAX", 0.50)
    b_openness_min = _cfg_float("DRUGGABILITY_TIER_B_OPENNESS_MIN", 0.55)
    b_openness_max = _cfg_float("DRUGGABILITY_TIER_B_OPENNESS_MAX", 0.65)
    b_polar_frac_min = _cfg_float("DRUGGABILITY_TIER_B_POLAR_FRAC_MIN", 0.45)
    b_polar_frac_max = _cfg_float("DRUGGABILITY_TIER_B_POLAR_FRAC_MAX", 0.60)

    triggers_c: List[str] = []
    if druggability < c_druggability_max:
        triggers_c.append(f"druggability<{c_druggability_max}")
    if polar_frac > c_polar_frac_min:
        triggers_c.append(f"polar_frac>{c_polar_frac_min}")
    if openness >= c_openness_min:
        triggers_c.append(f"open>={c_openness_min}")
    if has_metal:
        triggers_c.append("has_metal")
    if triggers_c:
        return "C", triggers_c

    if (
        druggability >= a_druggability_min
        and openness <= a_openness_max
        and polar_frac <= a_polar_frac_max
    ):
        return "A", ["tierA_all_criteria"]

    triggers_b: List[str] = []
    if b_druggability_min <= druggability < b_druggability_max:
        triggers_b.append(f"{b_druggability_min}<=druggability<{b_druggability_max}")
    if b_openness_min < openness < b_openness_max:
        triggers_b.append(f"{b_openness_min}<open<{b_openness_max}")
    if b_polar_frac_min < polar_frac <= b_polar_frac_max:
        triggers_b.append(f"{b_polar_frac_min}<polar_frac<={b_polar_frac_max}")
    if triggers_b:
        return "B", triggers_b

    return "B", ["fallback_uncertain"]


def _summarize_fpocket_info(
    info_path: Path,
    receptor_pdb: Path,
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    center: Tuple[float, float, float],
    logger: logging.Logger,
) -> Optional[DruggabilityMetrics]:
    """
    Parse the first pocket block from <stem>_info.txt, compute metrics, classify into Tier A/B/C.
    """
    if not info_path.exists() or info_path.stat().st_size == 0:
        logger.warning(
            "[druggability.fpocket.info.missing] pdb_id=%s variant=%s ph=%s path=%s",
            pdb_id,
            variant,
            ph_label,
            info_path,
        )
        return None

    druggability: Optional[float] = None
    volume: Optional[float] = None
    total_sasa: Optional[float] = None
    polar_sasa: Optional[float] = None
    openness: Optional[float] = None
    in_pocket = False

    try:
        with info_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith("Pocket "):
                    if in_pocket:
                        break
                    in_pocket = True
                    continue
                if not in_pocket:
                    continue
                if line.startswith("Druggability Score"):
                    try:
                        druggability = float(line.split()[-1])
                    except Exception:
                        pass
                    continue
                if line.startswith("Volume score"):
                    continue
                if line.startswith("Volume :"):
                    try:
                        volume = float(line.split()[-1])
                    except Exception:
                        pass
                    continue
                if line.startswith("Total SASA"):
                    try:
                        total_sasa = float(line.split()[-1])
                    except Exception:
                        pass
                    continue
                if line.startswith("Polar SASA"):
                    try:
                        polar_sasa = float(line.split()[-1])
                    except Exception:
                        pass
                    continue
                if line.startswith("Mean alp. sph. solvent access"):
                    try:
                        openness = float(line.split()[-1])
                    except Exception:
                        pass
                    continue
    except Exception as exc:
        logger.warning(
            "[druggability.fpocket.info.error] pdb_id=%s variant=%s ph=%s path=%s err=%s",
            pdb_id,
            variant,
            ph_label,
            info_path,
            exc,
        )
        return None

    if (
        druggability is None
        or volume is None
        or openness is None
        or total_sasa is None
        or polar_sasa is None
        or total_sasa <= 0.0
    ):
        logger.warning(
            "[druggability.fpocket.info.incomplete] pdb_id=%s variant=%s ph=%s path=%s",
            pdb_id,
            variant,
            ph_label,
            info_path,
        )
        return None

    polar_frac = polar_sasa / total_sasa if total_sasa > 0 else 0.0
    has_metal = _detect_metal_near_center(receptor_pdb, center)
    tier, triggers = _classify_druggability(druggability, volume, openness, polar_frac, has_metal)

    metrics = DruggabilityMetrics(
        pdb_id=pdb_id,
        variant=variant,
        ph_label=ph_label,
        druggability=druggability,
        volume=volume,
        openness=openness,
        polar_fraction=polar_frac,
        has_metal=has_metal,
        tier=tier,
        triggers=triggers,
        info_path=info_path,
    )

    logger.info(
        "[druggability.fpocket.tier] pdb_id=%s variant=%s ph=%s tier=%s druggability=%.3f volume=%.1f open=%.3f polar_frac=%.3f has_metal=%d triggers=%s",
        pdb_id,
        variant or "HOLO",
        ph_label or "base",
        metrics.tier,
        metrics.druggability,
        metrics.volume,
        metrics.openness,
        metrics.polar_fraction,
        1 if metrics.has_metal else 0,
        ",".join(metrics.triggers),
    )
    return metrics


def collect_pocket_residues_from_center(
    receptor_pdb: Union[str, Path],
    center: Tuple[float, float, float],
    box_size: Tuple[float, float, float],
    radius_margin: float = 1.5,
    max_radius: float = 10.0,
) -> List[Tuple[int, str, str]]:
    """
    Returns a sorted list of residues (res_seq, i_code, chain_id) that define the active-site pocket.
    """
    try:
        receptor_path = Path(receptor_pdb)
        half_max = max(float(s) for s in box_size) / 2.0
        radius = min(half_max + float(radius_margin), float(max_radius))
        radius_sq = radius * radius

        residues: set[Tuple[str, int, str]] = set()
        if not receptor_path.exists():
            logging.warning(
                "[druggability.fpocket.residues.skip] reason=missing_receptor path=%s",
                receptor_path,
            )
            return []

        with receptor_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue

                element = (line[76:78] or "").strip().upper()
                atom_name = (line[12:16] or "").strip()
                if element:
                    if element.startswith("H"):
                        continue
                else:
                    if atom_name.upper().startswith("H") or atom_name.startswith(" D"):
                        continue

                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue

                dx = x - float(center[0])
                dy = y - float(center[1])
                dz = z - float(center[2])
                if (dx * dx + dy * dy + dz * dz) > radius_sq:
                    continue

                chain_id = (line[21:22] or " ").strip() or " "
                try:
                    res_seq = int((line[22:26] or "0").strip() or "0")
                except ValueError:
                    res_seq = 0
                i_code = (line[26:27] or "").strip() or "-"
                residues.add((chain_id, res_seq, i_code))

        ordered = sorted(residues, key=lambda r: (r[0], r[1], r[2]))
        return [(res_seq, i_code, chain_id) for chain_id, res_seq, i_code in ordered]
    except Exception as exc:
        logging.warning(
            "[druggability.fpocket.residues.error] pdb=%s err=%s", receptor_pdb, exc
        )
        return []


def format_fpocket_pocket_spec(residues: List[Tuple[int, str, str]]) -> str:
    """
    Given a list of (res_seq, i_code, chain_id), build the -P argument.
    """
    parts: List[str] = []
    for res_seq, i_code, chain_id in sorted(residues, key=lambda r: (r[2], r[0], r[1])):
        icode_token = i_code if i_code else "-"
        chain_token = chain_id if chain_id else " "
        parts.append(f"{int(res_seq)}:{icode_token}:{chain_token}")
    return ".".join(parts)


def run_fpocket_for_explicit_pocket(
    cfg: Dict[str, Any],
    receptor_pdb: Union[str, Path],
    pocket_residues: List[Tuple[int, str, str]],
    logger: logging.Logger,
) -> Optional[Path]:
    """
    Run fpocket with -P on receptor_pdb and pocket_residues.
    """
    receptor_path = Path(receptor_pdb)
    stem = receptor_path.stem
    fpocket_exe = get_fpocket_exe(cfg)
    if fpocket_exe is None:
        logger.warning("[druggability.fpocket.skip] pdb=%s reason=missing_exe", receptor_path)
        return None
    if not pocket_residues:
        logger.warning(
            "[druggability.fpocket.skip] pdb=%s reason=no_residues", receptor_path
        )
        return None

    pocket_spec = format_fpocket_pocket_spec(pocket_residues)
    cmd = [
        str(fpocket_exe),
        "-f",
        str(receptor_path),
        "-P",
        pocket_spec,
    ]
    logger.info(
        "[druggability.fpocket.run] pdb=%s pocket_residues=%d cmd=%s",
        receptor_path,
        len(pocket_residues),
        cmd,
    )

    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except Exception as exc:
        logger.warning(
            "[druggability.fpocket.error] pdb=%s reason=fpocket_failed err=%s",
            receptor_path,
            exc,
        )
        return None

    src_out_dir = receptor_path.parent / f"{stem}_out"
    if not src_out_dir.exists():
        logger.warning(
            "[druggability.fpocket.skip] pdb=%s reason=missing_output path=%s",
            receptor_path,
            src_out_dir,
        )
        return None

    dst_root = get_fpocket_output_root(cfg)
    dst_out_dir = dst_root / f"{stem}_out"
    try:
        shutil.rmtree(dst_out_dir, ignore_errors=True)
        shutil.copytree(src_out_dir, dst_out_dir)
        shutil.rmtree(src_out_dir, ignore_errors=True)
    except Exception as exc:
        logger.warning(
            "[druggability.fpocket.move.error] src=%s dst=%s err=%s",
            src_out_dir,
            dst_out_dir,
            exc,
        )
        return None

    logger.info(
        "[druggability.fpocket.move] src=%s dst=%s removed_src=1",
        src_out_dir,
        dst_out_dir,
    )
    return dst_out_dir


def evaluate_druggability_for_active_site(
    cfg: Dict[str, Any],
    pdb_id: str,
    variant: Optional[str],
    ph_label: Optional[str],
    center: Tuple[float, float, float],
    box_size: Tuple[float, float, float],
    logger: logging.Logger,
) -> Optional[Path]:
    """
    Evaluate fpocket-based druggability for the already-identified active site.
    """
    try:
        paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    except Exception as exc:
        logger.warning(
            "[druggability.fpocket.skip] pdb=%s variant=%s ph=%s reason=paths err=%s",
            pdb_id,
            variant,
            ph_label,
            exc,
        )
        return None

    preferred = getattr(paths, "input_pdb_path", None)
    receptor_pdb = Path(preferred) if preferred and Path(preferred).exists() else Path(
        paths.receptor_cleaned_pdb(variant)
    )
    logger.info(
        "[druggability.fpocket.receptor] pdb_id=%s variant=%s ph=%s receptor_pdb=%s",
        pdb_id,
        variant,
        ph_label,
        receptor_pdb,
    )

    residues = collect_pocket_residues_from_center(receptor_pdb, center, box_size)
    if not residues:
        logger.warning(
            "[druggability.fpocket.residues.empty] pdb=%s variant=%s ph=%s",
            pdb_id,
            variant,
            ph_label,
        )
        return None

    out_dir = run_fpocket_for_explicit_pocket(cfg, receptor_pdb, residues, logger)
    if not out_dir:
        return None

    stem = Path(receptor_pdb).stem
    info_path = out_dir / f"{stem}_info.txt"
    _summarize_fpocket_info(
        info_path=info_path,
        receptor_pdb=receptor_pdb,
        pdb_id=pdb_id,
        variant=variant,
        ph_label=ph_label,
        center=center,
        logger=logger,
    )
    return out_dir


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Evaluate fpocket druggability for a PDB.")
    parser.add_argument("pdb", help="PDB ID or path to a receptor PDB file.")
    args = parser.parse_args()

    pdb_arg = Path(args.pdb)
    if pdb_arg.suffix.lower() == ".pdb":
        pdb_id = pdb_arg.stem
    else:
        pdb_id = str(pdb_arg).upper()

    try:
        from docking import get_active_site_center_and_size  # type: ignore
    except Exception as exc:  # pragma: no cover - defensive; avoids circular import issues
        print(f"Failed to import docking helpers: {exc}")
        raise SystemExit(1)

    cfg = config
    logger = logging.getLogger("druggability.cli")
    logging.basicConfig(level=logging.INFO)

    try:
        paths = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    except Exception as exc:
        print(f"Failed to build paths for {pdb_id}: {exc}")
        raise SystemExit(1)

    preferred = getattr(paths, "input_pdb_path", None)
    receptor_path = Path(preferred) if preferred and Path(preferred).exists() else Path(
        paths.receptor_cleaned_pdb("HOLO")
    )

    center_box = get_active_site_center_and_size(cfg, pdb_id, "HOLO", None, logger)
    if not center_box:
        print(f"No active site center/box found for {pdb_id}")
        raise SystemExit(1)

    center, box = center_box
    out_dir = evaluate_druggability_for_active_site(
        cfg,
        pdb_id=pdb_id,
        variant="HOLO",
        ph_label=None,
        center=center,
        box_size=box,
        logger=logger,
    )
    if out_dir:
        print(f"fpocket output: {out_dir}")
        metrics = _summarize_fpocket_info(
            info_path=out_dir / f"{receptor_path.stem}_info.txt",
            receptor_pdb=receptor_path,
            pdb_id=pdb_id,
            variant="HOLO",
            ph_label=None,
            center=center,
            logger=logger,
        )
        if metrics:
            print(f"tier={metrics.tier} druggability={metrics.druggability:.3f}")
    else:
        print("fpocket did not produce output")
