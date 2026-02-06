"""Non-standard residue stripping and ion policy helpers."""

from __future__ import annotations

import logging
from collections import Counter
from pathlib import Path
from typing import Iterable, Optional, Tuple, Union

from protein_prep.aliases_policy import (
    _ALIASES_BIND_LOGGED,
    _COFACTOR_CANONICAL,
    _COFACTOR_DROP_LOGGED,
    _COFACTOR_NAMES,
    _COFACTOR_RAW_ALL,
    _ELEM_CANON,
    _ION_AUDIT_METALS,
    _ION_KEEP_LOGGED,
    _ION_AUDIT_SIMPLE_IONS,
    _METAL_RESNAMES,
    _POLICY_MODE,
    _RETAIN_VARIANT,
    _RETAIN_VARIANT_CANONICAL,
    _SALT_RESNAMES,
    _WATER_NAMES,
    _normalize_resname,
)
from protein_prep.receptor_prep import _normalize_ion_policy


def _hydrate_legacy_globals() -> None:
    import automate_protein_prep as legacy

    g = globals()
    for name, value in legacy.__dict__.items():
        g.setdefault(name, value)

def _bucket_counts(counter: Counter) -> dict[str, int]:
    _hydrate_legacy_globals()
    keys = ["ZN", "MG", "NA", "K", "CA", "MN", "FE", "CL"]
    out = {k: 0 for k in keys}
    other = 0
    for resn, count in counter.items():
        token = resn.upper()
        if token in out:
            out[token] += count
        else:
            other += count
    out["OTHER"] = other
    return out

def _format_counts(counter: Counter) -> str:
    _hydrate_legacy_globals()
    bucketed = _bucket_counts(counter)
    keys = ["ZN", "MG", "NA", "K", "CA", "MN", "FE", "CL", "OTHER"]
    parts = [f"{k}={bucketed.get(k, 0)}" for k in keys]
    return " ".join(parts)

def _maybe_strip_ions(
    pdb_path: Union[str, Path],
    cfg: Optional[dict] = None,
    *,
    variant: Optional[str] = None,
    pocket_center: Optional[tuple[float, float, float]] = None,
    extra_keep: Optional[Iterable[str]] = None,
) -> int:
    _hydrate_legacy_globals()
    path = Path(pdb_path)
    if not path.exists():
        logging.warning("[ions] skip_missing file=%s", path)
        return 0

    allow_tokens, _ = _load_retain_allowlist(cfg)
    allow_canonical: set[str] = set()
    for token in allow_tokens:
        if not token:
            continue
        text = str(token).strip().upper()
        if not text:
            continue
        canonical = _normalize_resname(text)
        allow_canonical.add(canonical or text)

    extra_canonical: set[str] = set()
    if extra_keep:
        for token in extra_keep:
            if not token:
                continue
            text = str(token).strip().upper()
            if not text:
                continue
            canonical = _normalize_resname(text)
            extra_canonical.add(canonical or text)

    holo_keep_tokens = allow_canonical | extra_canonical
    apo_keep_tokens = set(extra_canonical)

    policy = _normalize_ion_policy(cfg)
    variant_token = _resolve_variant_token(cfg, variant)
    variant_label = variant_token or "legacy"

    log_pre_variant = getattr(
        _activesite_mod, "log_pre_variant_policy_breadcrumb", None
    )
    if log_pre_variant is not None:
        try:
            log_pre_variant(path)
        except Exception:
            pass

    emit_ion_audit_probe("pre_variant_policy", path, variant=variant_token)

    radius_cfg = 6.0
    if cfg is not None:
        try:
            radius_cfg = float(cfg.get("HOLO_SALT_STRIP_RADIUS", 6.0) or 0.0)
        except Exception:
            radius_cfg = 6.0
    radius_cfg = max(0.0, radius_cfg)
    radius = radius_cfg if (policy == "by_variant" and variant_token == "HOLO") else 0.0
    radius_term = f"{radius:.2f}" if radius > 0.0 else "none"

    # [ions] stage=clean instrumentation
    drop_free_ions = policy != "never_strip"
    if variant_token == "HOLO":
        allow_display_set = holo_keep_tokens
        block_candidates = {
            tok
            for tok in _SALT_RESNAMES
            if _normalize_resname(tok) not in holo_keep_tokens
        }
    elif variant_token == "APO":
        allow_display_set = apo_keep_tokens
        block_candidates = {
            tok
            for tok in (_METAL_RESNAMES | _SALT_RESNAMES)
            if _normalize_resname(tok) not in apo_keep_tokens
        }
    else:
        allow_display_set = allow_canonical
        block_candidates = set()

    logging.info(
        "[ions.policy] stage=central variant=%s drop_free_ions=%s allow_list=%s block_list=%s radius=%s policy=%s file=%s",
        variant_label,
        str(drop_free_ions).lower(),
        _format_token_list(allow_display_set),
        _format_token_list(block_candidates),
        radius_term,
        policy,
        path,
    )

    try:
        text = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception as exc:
        logging.warning("[ions] read_failed file=%s err=%s", path, exc)
        return 0

    residues: dict[tuple[str, str, str, str], list[tuple[int, str]]] = {}
    for idx, line in enumerate(text):
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip().upper()
        key = (line[21], line[22:26], line[26], resname)
        residues.setdefault(key, []).append((idx, line))

    before = Counter()
    stripped = Counter()
    kept_counter = Counter()
    category_totals: dict[str, dict[str, int]] = {
        "metal": {"kept": 0, "stripped": 0},
        "salt": {"kept": 0, "stripped": 0},
        "other": {"kept": 0, "stripped": 0},
    }

    warn_missing_center = False
    remove_indices: set[int] = set()

    for key, atoms in residues.items():
        if len(atoms) != 1:
            continue
        chain, resi, icode, resname = key
        res_token = resname.upper()
        canonical_res = _normalize_resname(res_token)
        line_idx, line = atoms[0]
        before[res_token] += 1

        remove = False
        reason = ""
        elem = (line[76:78].strip() or res_token).upper()
        canonical_elem = _normalize_resname(elem)

        is_metal = (canonical_res in _METAL_RESNAMES) or (
            canonical_elem in _METAL_RESNAMES
        )
        is_salt = (canonical_res in _SALT_RESNAMES) or (
            canonical_elem in _SALT_RESNAMES
        )
        category = "metal" if is_metal else "salt" if is_salt else "other"

        if policy == "never_strip":
            remove = False
        elif policy == "always_strip":
            if canonical_res in holo_keep_tokens:
                remove = False
            else:
                remove = True
                reason = "policy_always_strip"
        else:  # policy == by_variant
            if variant_token == "HOLO":
                if (canonical_res in holo_keep_tokens) or is_metal:
                    remove = False
                elif is_salt:
                    if radius > 0.0:
                        dist = _distance_from_center(line, pocket_center)
                        if dist is None:
                            warn_missing_center = True
                            remove = False
                        elif dist >= radius:
                            remove = True
                            reason = "salt_far"
                        else:
                            remove = False
                    else:
                        remove = False
                else:
                    remove = False
            else:
                if canonical_res in apo_keep_tokens:
                    remove = False
                    reason = "apo_keep_override"
                else:
                    if is_metal:
                        remove = True
                        reason = "apo_strip_metal"
                    elif is_salt:
                        remove = True
                        reason = "apo_strip_salt"
                    else:
                        remove = True
                        reason = "apo_strip_other"

        if remove:
            stripped[res_token] += 1
            remove_indices.add(line_idx)
            logging.debug(
                "[ions.remove] elem=%s resname=%s serial=%s chain=%s resi=%s reason=%s",
                elem,
                res_token,
                line[6:11].strip(),
                chain.strip() or "-",
                (resi or "0").strip() or "0",
                reason,
            )
            category_totals[category]["stripped"] += 1
        else:
            kept_counter[res_token] += 1
            category_totals[category]["kept"] += 1

    logging.info(
        "[ions.counts.before] stage=clean variant=%s file=%s detail=%s",
        variant_label,
        path,
        _format_counts(before),
    )
    kept = before - stripped
    logging.info(
        "[ions.counts.after] stage=clean variant=%s file=%s detail=%s",
        variant_label,
        path,
        _format_counts(kept),
    )
    total_kept = sum(kept_counter.values())
    total_stripped = sum(stripped.values())
    logging.info(
        "[ions.summary] variant=%s kept=%d stripped=%d metals_kept=%d metals_stripped=%d salts_kept=%d salts_stripped=%d",
        variant_label,
        total_kept,
        total_stripped,
        category_totals["metal"]["kept"],
        category_totals["metal"]["stripped"],
        category_totals["salt"]["kept"],
        category_totals["salt"]["stripped"],
    )

    if warn_missing_center and radius > 0.0 and variant_token == "HOLO":
        logging.warning(
            "[ions] holo_salt_radius_set_but_no_center action=keep_salts radius=%.2f",
            radius,
        )

    if not remove_indices:
        logging.debug("[ions] no_monoatomic_hits remove=0")
        return 0

    for idx in sorted(remove_indices):
        text[idx] = None  # type: ignore

    rewritten = [ln for ln in text if ln is not None]
    try:
        path.write_text("\n".join(rewritten) + "\n", encoding="utf-8")
    except Exception as exc:
        logging.warning("[ions] write_failed file=%s err=%s", path, exc)
        return 0

    return len(remove_indices)

def _is_element_token(sym):
    _hydrate_legacy_globals()
    canonical = _normalize_resname(sym)
    return bool(canonical) and (canonical in _ELEM_CANON)

def _is_retained_ion(resname: str) -> bool:
    _hydrate_legacy_globals()
    canonical = _normalize_resname(resname)
    return bool(canonical) and (canonical in _ELEM_CANON)

def _cofactor_policy_keep(resname: str) -> bool:
    _hydrate_legacy_globals()
    """
    YAML-only policy:
      • Keep anything listed under 'retain_in_receptor_resnames'
      • Keep standard residues (handled elsewhere)
      • Water handling is done in strip_nonstandard_residues(), so return False here for waters.
      • Everything else → remove
    """
    rn_upper = (resname or "").strip().upper()
    if not rn_upper:
        return False
    if rn_upper in _WATER_NAMES:
        return False
    canonical = _normalize_resname(rn_upper)
    if canonical and canonical in _COFACTOR_CANONICAL:
        return True
    return rn_upper in _COFACTOR_NAMES

def strip_nonstandard_residues(
    input_pdb: Union[str, Path],
    output_pdb: Union[str, Path],
    *,
    variant: Optional[str] = None,
) -> Tuple[int, str]:
    _hydrate_legacy_globals()
    """
    Remove nonstandard residues while keeping what the YAML says to keep
    (retain_in_receptor_resnames) and standard amino acids.
    Water handling still honors your numeric policy (radius/B-factor) but
    the water names themselves come from the YAML.
    """
    variant_token = _resolve_variant_token(config, variant)
    variant_label = variant_token or "legacy"
    global _ALIASES_BIND_LOGGED
    if not _ALIASES_BIND_LOGGED:
        logging.info(
            "[aliases.bind] water_set=%d cofactor_set=%d elem_tokens=%d",
            len(_WATER_NAMES),
            len(_COFACTOR_NAMES),
            len(_ELEM_CANON),
        )
        _ALIASES_BIND_LOGGED = True
    before_map = _scan_metal_map(input_pdb, ALIASES)
    before_summary = _summarize_ions_file(input_pdb)
    logging.info(
        "[stripnsr.before] variant=%s ions_res=%s ions_elem=%s",
        variant_label,
        before_summary.get("res_hist", "none"),
        before_summary.get("elem_hist", "none"),
    )
    standard_residues = {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "HID",
        "HIE",
        "HIP",
        "SEC",
        "PYL",
        "MSE",
    }

    removed: Set[str] = set()
    removed_hits: List[Tuple[str, str, str, str]] = []
    kept_lines: List[str] = []

    # Estimate pocket center from any YAML-retained cofactors/metals present
    cofm_xyz: List[Tuple[float, float, float]] = []
    with open(input_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip().upper()
            canonical_res = _normalize_resname(resname)
            if (resname in _RETAIN_VARIANT) or (
                canonical_res in _RETAIN_VARIANT_CANONICAL
            ):
                xyz = _parse_xyz(line)
                if xyz:
                    cofm_xyz.append(xyz)

    pocket_center: Optional[Tuple[float, float, float]] = None
    if cofm_xyz:
        from statistics import fmean

        xs, ys, zs = zip(*cofm_xyz)
        pocket_center = (fmean(xs), fmean(ys), fmean(zs))

    water_policy = (config.get("WATER_POLICY", "site_only") or "site_only").lower()
    water_radius = float(config.get("WATER_SITE_RADIUS_ANG", 6.0))
    water_bmax = float(config.get("WATER_MAX_BFACTOR", 60.0))
    logging.info(
        "[water.policy] mode=%s policy=%s waters_kept_rule_applied=true",
        _POLICY_MODE,
        water_policy,
    )

    with open(input_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ATOM  "):
                resname = line[17:20].strip().upper()
                if resname in standard_residues:
                    kept_lines.append(line)
                else:
                    removed.add(resname)
                continue

            if line.startswith("HETATM"):
                resname = line[17:20].strip().upper()
                chain = (line[21] or "-").strip() or "-"
                resseq = (line[22:26] or "0").strip() or "0"
                elem_token = (line[76:78].strip() or resname).upper()
                reason = None
                canonical = _normalize_resname(resname)
                canonical_token = canonical or resname

                # Water names come from YAML (no hardcoded list here)
                if resname in _WATER_NAMES:
                    keep_line = False
                    if water_policy == "keep_all":
                        keep_line = True
                    elif water_policy == "remove_all":
                        removed.add(resname)
                        reason = "water_remove_all"
                    elif pocket_center is not None:
                        xyz = _parse_xyz(line)
                        keep = False
                        if xyz:
                            dx = xyz[0] - pocket_center[0]
                            dy = xyz[1] - pocket_center[1]
                            dz = xyz[2] - pocket_center[2]
                            dist2 = dx * dx + dy * dy + dz * dz
                            if dist2 <= water_radius * water_radius:
                                try:
                                    b = float(line[60:66])
                                    keep = b <= water_bmax
                                except Exception:
                                    keep = True
                        if keep:
                            keep_line = True
                        else:
                            removed.add(resname)
                            reason = "water_far"
                    else:
                        removed.add(resname)
                        reason = "water_no_center"

                    if keep_line:
                        kept_lines.append(line)
                        continue

                    if reason and (
                        resname in _ION_AUDIT_METALS
                        or resname in _ION_AUDIT_SIMPLE_IONS
                        or elem_token in _ION_AUDIT_METALS
                        or elem_token in _ION_AUDIT_SIMPLE_IONS
                    ):
                        removed_hits.append((resname, chain, resseq, reason))
                        logging.info(
                            "[stripnsr.hit] resname=%s chain=%s resSeq=%s reason=%s variant=%s",
                            resname,
                            chain,
                            resseq,
                            reason,
                            variant_label,
                        )
                    continue

                keep_line = False
                if canonical_token and canonical_token in _ELEM_CANON:
                    keep_line = True
                    if canonical_token not in _ION_KEEP_LOGGED:
                        logging.info(
                            "[ion.keep] token=%s source=elem_tokens_canonical",
                            canonical_token,
                        )
                        _ION_KEEP_LOGGED.add(canonical_token)
                elif _cofactor_policy_keep(resname):
                    keep_line = True
                elif resname in _RETAIN_VARIANT or (
                    canonical and canonical in _RETAIN_VARIANT_CANONICAL
                ):
                    keep_line = True

                if keep_line:
                    kept_lines.append(line)
                else:
                    removed.add(resname)
                    reason = "not_in_retain"
                    if _POLICY_MODE == "APO" and (
                        resname in _COFACTOR_RAW_ALL
                        or (canonical and canonical in _COFACTOR_RAW_ALL)
                    ):
                        reason = "apo_policy"
                        drop_key = canonical or resname
                        if drop_key and drop_key not in _COFACTOR_DROP_LOGGED:
                            logging.info(
                                "[cofactor.drop] mode=APO resname=%s reason=apo_policy",
                                drop_key,
                            )
                            _COFACTOR_DROP_LOGGED.add(drop_key)
                if not keep_line and (
                    resname in _ION_AUDIT_METALS
                    or resname in _ION_AUDIT_SIMPLE_IONS
                    or elem_token in _ION_AUDIT_METALS
                    or elem_token in _ION_AUDIT_SIMPLE_IONS
                ):
                    removed_hits.append(
                        (resname, chain, resseq, reason or "not_in_retain")
                    )
                    logging.info(
                        "[stripnsr.hit] resname=%s chain=%s resSeq=%s reason=%s variant=%s",
                        resname,
                        chain,
                        resseq,
                        reason or "not_in_retain",
                        variant_label,
                    )
                continue

            # non-coordinate records pass through
            kept_lines.append(line)

    with open(output_pdb, "w", encoding="utf-8") as f:
        f.writelines(kept_lines)
    _post_write_element_guard("strip_nonstandard", output_pdb)

    after_map = _scan_metal_map(output_pdb, ALIASES)
    after_summary = _summarize_ions_file(output_pdb)
    logging.info(
        "[stripnsr.after] variant=%s ions_res=%s ions_elem=%s removed_total=%d",
        variant_label,
        after_summary.get("res_hist", "none"),
        after_summary.get("elem_hist", "none"),
        len(removed_hits),
    )
    kept_res = sorted([res for res, count in after_map.items() if count > 0])
    stripped_res = sorted(
        {res for res, count in before_map.items() if after_map.get(res, 0) < count}
    )
    logging.info(
        "[stripnsr.diff] kept=%s stripped=%s",
        ",".join(kept_res) if kept_res else "none",
        ",".join(stripped_res) if stripped_res else "none",
    )
    retain_targets = {tok.upper() for tok in getattr(ALIASES, "retain_resnames", [])}
    retain_targets_canonical = {
        _normalize_resname(tok)
        for tok in getattr(ALIASES, "retain_resnames", [])
        if _normalize_resname(tok)
    }
    for res in stripped_res:
        if (res in retain_targets) or (
            _normalize_resname(res) in retain_targets_canonical
        ):
            logging.warning("[stripnsr.violation] resname=%s", res)

    logging.info("Removed nonstandard residues (YAML-driven): %s", sorted(removed))
    if pocket_center:
        logging.info("Estimated pocket center: (%.2f, %.2f, %.2f)", *pocket_center)
    return len(removed), str(output_pdb)
