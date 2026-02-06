from __future__ import annotations

import logging
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

try:
    import requests  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    requests = None  # type: ignore

try:
    import yaml  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    yaml = None  # type: ignore

from installation import load_config
try:
    from chemdb.path_router import make_paths
except Exception:
    try:
        from path_router import make_paths
    except Exception:
        from src.path_router.path_router import make_paths

# Load config once at import to mirror other helpers.
config = load_config()

# Defaults
_DEFAULT_FPOCKET_EXE = Path("/home/michael/atlas/tools/fpocket/bin/fpocket")
_DEFAULT_FPOCKET_OUTPUT_ROOT = Path(__file__).resolve().parent / "fpocket"
_PROTEIN_CLASS_RULES: Optional[Dict[str, Any]] = None


def _cfg_float(key: str, default: float) -> float:
    raw = config.get(key)
    try:
        return float(raw)
    except Exception:
        return default


def _cfg_str(key: str, default: str = "") -> str:
    raw = config.get(key)
    if raw is None:
        return default
    return str(raw)


def _cfg_str_list(key: str, default: str = "") -> List[str]:
    """
    Split a semi-colon or comma separated string into a cleaned list of lowercased tokens.
    """
    value = _cfg_str(key, default)
    if not value:
        return []
    tokens: List[str] = []
    for chunk in value.replace(",", ";").split(";"):
        t = chunk.strip().lower()
        if t:
            tokens.append(t)
    return tokens


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

    # New metadata fields (best-effort; may be None/empty)
    uniprot_id: Optional[str] = None
    ec_numbers: List[str] = field(default_factory=list)
    protein_name: Optional[str] = None
    protein_family: Optional[str] = None
    protein_class: Optional[str] = None
    # Family-based prior difficulty (optional)
    family_prior_tier: Optional[str] = None
    family_prior_triggers: List[str] = field(default_factory=list)


def _fetch_uniprot_and_ec_for_pdb(
    pdb_id: str,
    logger: logging.Logger,
    entity_id: int = 1,
    timeout: float = 10.0,
) -> Tuple[Optional[str], List[str], Optional[str]]:
    """
    Best-effort lookup of UniProt ID and EC numbers for a PDB entry via RCSB + UniProt.

    Returns (uniprot_id, ec_numbers, pdb_description).
    """
    pdb_id = pdb_id.upper()
    if requests is None:
        logger.debug(
            "[druggability.family.skip] pdb_id=%s reason=missing_requests", pdb_id
        )
        return None, [], None

    uniprot_id: Optional[str] = None
    ec_numbers: List[str] = []
    pdb_description: Optional[str] = None

    annotation_url = (
        f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
    )
    try:
        logger.info(
            "[druggability.family.rcsb.query] pdb_id=%s entity=%d url=%s",
            pdb_id,
            entity_id,
            annotation_url,
        )
        resp = requests.get(annotation_url, timeout=timeout)
        resp.raise_for_status()
        entity_data = resp.json()

        identifiers = entity_data.get("rcsb_polymer_entity_container_identifiers", {})
        uniprot_ids = identifiers.get("uniprot_ids", []) or []
        if uniprot_ids:
            uniprot_id = str(uniprot_ids[0]).strip() or None

        ec_list = entity_data.get("entity", {}).get("rcsb_enzyme_class_list", [])
        for ec_entry in ec_list:
            ec = ec_entry.get("ec")
            if ec:
                ec_clean = str(ec).strip().rstrip(".")
                if ec_clean and ec_clean not in ec_numbers:
                    ec_numbers.append(ec_clean)

        pdb_description = entity_data.get("rcsb_polymer_entity", {}).get(
            "pdbx_description", None
        )
    except Exception as exc:
        logger.warning(
            "[druggability.family.rcsb.error] pdb_id=%s url=%s err=%s",
            pdb_id,
            annotation_url,
            exc,
        )

    if not uniprot_id:
        logger.warning(
            "[druggability.family.uniprot.missing] pdb_id=%s reason=no_uniprot_from_rcsb",
            pdb_id,
        )
        return None, ec_numbers, pdb_description

    if not ec_numbers:
        uniprot_txt_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"
        try:
            logger.info(
                "[druggability.family.uniprot.query] pdb_id=%s uniprot=%s url=%s",
                pdb_id,
                uniprot_id,
                uniprot_txt_url,
            )
            u_resp = requests.get(uniprot_txt_url, timeout=timeout)
            u_resp.raise_for_status()
            for line in u_resp.text.splitlines():
                if line.startswith("DE   EC="):
                    chunk = line.split("EC=", 1)[-1]
                    for token in chunk.split(";"):
                        ec_candidate = token.strip()
                        if ec_candidate and "." in ec_candidate and ec_candidate != "-":
                            ec_clean = ec_candidate.rstrip(".")
                            if ec_clean and ec_clean not in ec_numbers:
                                ec_numbers.append(ec_clean)
        except Exception as exc:
            logger.warning(
                "[druggability.family.uniprot.error] pdb_id=%s uniprot=%s url=%s err=%s",
                pdb_id,
                uniprot_id,
                uniprot_txt_url,
                exc,
            )

    return uniprot_id, ec_numbers, pdb_description


def _fetch_protein_name_and_family_from_uniprot(
    uniprot_id: str,
    logger: logging.Logger,
    timeout: float = 10.0,
    pdb_description: Optional[str] = None,
) -> Tuple[Optional[str], Optional[str]]:
    """
    Best-effort extraction of (protein_name, protein_family) from UniProt text entry.

    For now, protein_family is a simple textual label you can later postprocess;
    we prefer DE RecName and fall back to RCSB description or the UniProt ID.
    """
    if requests is None:
        logger.debug(
            "[druggability.family.uniprot.skip] uniprot=%s reason=missing_requests",
            uniprot_id,
        )
        return None, pdb_description or None

    protein_name: Optional[str] = None
    family_label: Optional[str] = None

    txt_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"
    try:
        logger.info(
            "[druggability.family.uniprot.name.query] uniprot=%s url=%s",
            uniprot_id,
            txt_url,
        )
        resp = requests.get(txt_url, timeout=timeout)
        resp.raise_for_status()
        lines = resp.text.splitlines()

        for line in lines:
            if line.startswith("DE   RecName: Full=") and not protein_name:
                try:
                    chunk = line.split("Full=", 1)[1]
                    protein_name = chunk.strip().rstrip(";.")
                except Exception:
                    continue
            if line.startswith("ID   ") and protein_name is None:
                try:
                    token = line.split()[1]
                    protein_name = token.strip()
                except Exception:
                    continue

        if protein_name:
            family_label = protein_name
        elif pdb_description:
            family_label = pdb_description
        else:
            family_label = uniprot_id

    except Exception as exc:
        logger.warning(
            "[druggability.family.uniprot.name.error] uniprot=%s url=%s err=%s",
            uniprot_id,
            txt_url,
            exc,
        )
        if pdb_description:
            family_label = pdb_description
        else:
            family_label = uniprot_id

    return protein_name, family_label


def _load_protein_class_rules(logger: logging.Logger) -> Dict[str, Any]:
    global _PROTEIN_CLASS_RULES
    if _PROTEIN_CLASS_RULES is not None:
        return _PROTEIN_CLASS_RULES

    rules: Dict[str, Any] = {}
    if yaml is None:
        logger.debug("[druggability.family.yaml.skip] reason=missing_yaml_module")
        _PROTEIN_CLASS_RULES = rules
        return rules

    try:
        root = Path(__file__).resolve().parent
        yaml_path = root / "chemdb" / "protein_names.yaml"
        if not yaml_path.is_file():
            logger.warning("[druggability.family.yaml.missing] path=%s", yaml_path)
            _PROTEIN_CLASS_RULES = rules
            return rules

        with yaml_path.open("r", encoding="utf-8") as fh:
            loaded = yaml.safe_load(fh) or {}
        if not isinstance(loaded, dict):
            logger.warning(
                "[druggability.family.yaml.invalid] path=%s type=%s",
                yaml_path,
                type(loaded),
            )
            _PROTEIN_CLASS_RULES = rules
            return rules

        rules = loaded
        classes = (loaded.get("classes") or {}).keys()
        logger.info(
            "[druggability.family.yaml.loaded] path=%s classes=%s",
            yaml_path,
            ",".join(classes),
        )
    except Exception as exc:
        logger.warning("[druggability.family.yaml.error] err=%s", exc)
        rules = {}

    _PROTEIN_CLASS_RULES = rules
    return rules


def _extract_pdb_header_text(receptor_pdb: Path, logger: logging.Logger) -> str:
    """
    Extract informative header/keyword lines from a PDB file.
    """
    if not receptor_pdb.is_file():
        return ""

    lines: List[str] = []
    try:
        with receptor_pdb.open("r", encoding="utf-8", errors="ignore") as fh:
            for i, line in enumerate(fh):
                if i > 500:
                    break
                tag = line[:6].strip().upper()
                if tag in {"HEADER", "TITLE", "KEYWDS", "COMPND"}:
                    lines.append(line[6:].strip())
    except Exception as exc:
        logger.warning(
            "[druggability.family.pdb_header.error] path=%s err=%s",
            receptor_pdb,
            exc,
        )
        return ""

    return " ".join(lines)


def _infer_family_prior_from_metadata(
    uniprot_id: Optional[str],
    protein_name: Optional[str],
    protein_family: Optional[str],
    ec_numbers: List[str],
    pdb_header_text: str,
    logger: logging.Logger,
) -> Tuple[Optional[str], List[str], Optional[str]]:
    """
    Use YAML-driven canonical classes to infer a family-based difficulty prior.

    Returns (family_prior_tier, triggers, protein_class_name).
    """
    rules = _load_protein_class_rules(logger)
    classes = rules.get("classes") or {}
    if not isinstance(classes, dict) or not classes:
        return None, [], None

    parts: List[str] = []
    for val in (pdb_header_text, protein_name, protein_family, uniprot_id):
        if val:
            parts.append(str(val))
    if ec_numbers:
        parts.extend(ec_numbers)
    text = " ".join(parts).lower()

    best_class: Optional[str] = None
    best_tier: Optional[str] = None
    triggers: List[str] = []
    best_val: Optional[int] = None

    def tier_to_val(t: str) -> int:
        return {"A": 0, "B": 1, "C": 2}.get(t.upper(), 1)

    for class_name, info in classes.items():
        if not isinstance(info, dict):
            continue
        cname = str(class_name).strip()
        tier_prior = str(info.get("tier_prior", "B")).upper()
        if tier_prior not in {"A", "B", "C"}:
            tier_prior = "B"
        class_val = tier_to_val(tier_prior)

        text_keywords = info.get("text_keywords") or []
        ec_prefixes = info.get("ec_prefixes") or []

        class_triggers: List[str] = []

        for raw_kw in text_keywords:
            kw = str(raw_kw).strip().lower()
            if kw and kw in text:
                class_triggers.append(f"{cname}:kw:{kw}")

        for raw_prefix in ec_prefixes:
            prefix = str(raw_prefix).strip().lower()
            if not prefix:
                continue
            for ec in ec_numbers:
                if ec.lower().startswith(prefix):
                    class_triggers.append(f"{cname}:ec:{prefix}")
                    break

        if not class_triggers:
            continue

        if best_val is None or class_val > best_val:
            best_val = class_val
            best_class = cname
            best_tier = tier_prior
            triggers = class_triggers

    if best_class and best_tier:
        logger.info(
            "[druggability.family.prior] uniprot=%s class=%s tier=%s triggers=%s",
            uniprot_id or "",
            best_class,
            best_tier,
            ",".join(triggers),
        )

    return best_tier, triggers, best_class


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
        logging.debug(
            "[druggability.fpocket.metal.detect.error] path=%s err=%s",
            receptor_pdb,
            exc,
        )
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
    tier, triggers = _classify_druggability(
        druggability, volume, openness, polar_frac, has_metal
    )

    pdb_header_text = _extract_pdb_header_text(receptor_pdb, logger)

    # --- New: protein family metadata ---
    uniprot_id: Optional[str] = None
    ec_numbers: List[str] = []
    protein_name: Optional[str] = None
    protein_family: Optional[str] = None
    try:
        uniprot_id, ec_numbers, pdb_desc = _fetch_uniprot_and_ec_for_pdb(
            pdb_id=pdb_id,
            logger=logger,
        )
        if uniprot_id:
            protein_name, protein_family = _fetch_protein_name_and_family_from_uniprot(
                uniprot_id=uniprot_id,
                logger=logger,
                pdb_description=pdb_desc,
            )
    except Exception as exc:
        logger.warning(
            "[druggability.fpocket.family.error] pdb_id=%s variant=%s ph=%s err=%s",
            pdb_id,
            variant,
            ph_label,
            exc,
        )
    family_prior_tier: Optional[str] = None
    family_prior_triggers: List[str] = []
    protein_class: Optional[str] = None
    try:
        (
            family_prior_tier,
            family_prior_triggers,
            protein_class,
        ) = _infer_family_prior_from_metadata(
            uniprot_id=uniprot_id,
            protein_name=protein_name,
            protein_family=protein_family,
            ec_numbers=ec_numbers or [],
            pdb_header_text=pdb_header_text,
            logger=logger,
        )
    except Exception as exc:
        logger.warning(
            "[druggability.fpocket.family_prior.error] pdb_id=%s variant=%s ph=%s err=%s",
            pdb_id,
            variant,
            ph_label,
            exc,
        )

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
        uniprot_id=uniprot_id,
        ec_numbers=ec_numbers or [],
        protein_name=protein_name,
        protein_family=protein_family,
        protein_class=protein_class,
        family_prior_tier=family_prior_tier,
        family_prior_triggers=family_prior_triggers,
    )

    logger.info(
        "[druggability.fpocket.tier] pdb_id=%s variant=%s ph=%s tier=%s druggability=%.3f volume=%.1f open=%.3f polar_frac=%.3f has_metal=%d uniprot=%s class=%s family=%s family_prior=%s family_triggers=%s ec=%s triggers=%s",
        pdb_id,
        variant or "HOLO",
        ph_label or "base",
        metrics.tier,
        metrics.druggability,
        metrics.volume,
        metrics.openness,
        metrics.polar_fraction,
        1 if metrics.has_metal else 0,
        metrics.uniprot_id or "",
        metrics.protein_class or "",
        metrics.protein_family or "",
        metrics.family_prior_tier or "",
        ";".join(metrics.family_prior_triggers or []),
        ";".join(metrics.ec_numbers or []),
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
        logger.warning(
            "[druggability.fpocket.skip] pdb=%s reason=missing_exe", receptor_path
        )
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
    receptor_pdb = (
        Path(preferred)
        if preferred and Path(preferred).exists()
        else Path(paths.receptor_cleaned_pdb(variant))
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

    parser = argparse.ArgumentParser(
        description="Evaluate fpocket druggability for a PDB."
    )
    parser.add_argument("pdb", help="PDB ID or path to a receptor PDB file.")
    args = parser.parse_args()

    pdb_arg = Path(args.pdb)
    if pdb_arg.suffix.lower() == ".pdb":
        pdb_id = pdb_arg.stem
    else:
        pdb_id = str(pdb_arg).upper()

    try:
        from docking.docking import get_active_site_center_and_size  # type: ignore
    except (
        Exception
    ) as exc:  # pragma: no cover - defensive; avoids circular import issues
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
    receptor_path = (
        Path(preferred)
        if preferred and Path(preferred).exists()
        else Path(paths.receptor_cleaned_pdb("HOLO"))
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
