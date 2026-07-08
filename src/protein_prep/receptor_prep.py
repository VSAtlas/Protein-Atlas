"""Receptor PDBQT preparation with Meeko and related cleanup helpers."""

from __future__ import annotations

import logging
import json
import os
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Iterable, Optional, Sequence, Union

import protein_prep.prep_utils as prep_utils
from protein_prep.aliases_policy import (
    ALIASES,
    RULES,
    _ELEM_CANON,
    _METAL_RESNAMES,
    _RETAIN_VARIANT,
    _SALT_RESNAMES,
    _flatten_semicolons,
    _normalize_resname,
)

from protein_prep.element_guard import _meeko_preflight_or_fail
from protein_prep.pdb_records import line_xyz
from protein_prep.prep_utils import (
    _cfg_bool,
    _load_retain_allowlist,
    _persist_subproc,
    _resolve_variant_token,
)
from protein_prep.tool_runners import (
    build_meeko_base_cmd,
    build_meeko_legacy_cmd,
    build_meeko_modern_cmd,
    run_subprocess_capture,
)
from protein_prep.geometry.repair import repair_terminal_oxt_geometry_in_pdb
from protein_prep.meeko.normalization import (
    meeko_pruned_stage_label,
    write_meeko_clash_pruned_copy,
    write_meeko_normalized_copy,
)
from protein_prep.meeko.safe_export import (
    meeko_safe_stage_label,
    parse_unmatched_residue_keys,
    write_meeko_safe_receptor_copy,
)
from protein_prep.meeko.stage_review import (
    meeko_stage_requires_review as _meeko_stage_requires_review,
)
_MEEKO_DROP_IONS = set(_flatten_semicolons(RULES.get("meeko_drop_free_ions", []))) or {
    "NA",
    "K",
    "LI",
}


def _is_element_token(sym):
    canonical = _normalize_resname(sym)
    return bool(canonical) and (canonical in _ELEM_CANON)


def _ion_pairs_from_records(
    file_path: Union[str, Path], prefixes: Sequence[str]
) -> set[tuple[str, str]]:
    path = Path(file_path)
    pairs: set[tuple[str, str]] = set()
    if not path.exists():
        return pairs
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not any(line.startswith(prefix) for prefix in prefixes):
                continue
            res = line[17:20].strip().upper()
            if not res or not _is_element_token(res):
                continue
            chain = (line[21:22] or "-").strip() or "-"
            resseq = (line[22:26] or "0").strip() or "0"
            pairs.add((res, f"{chain}:{resseq}"))
    return pairs


def _ion_pairs_from_pdb(file_path: Union[str, Path]) -> set[tuple[str, str]]:
    return _ion_pairs_from_records(file_path, ("HETATM",))


def _ion_pairs_from_pdbqt(file_path: Union[str, Path]) -> set[tuple[str, str]]:
    return _ion_pairs_from_records(file_path, ("ATOM  ", "HETATM"))


def _format_ion_pairs(pairs: Iterable[tuple[str, str]]) -> str:
    seq = sorted(pairs)
    if not seq:
        return "none"
    return ",".join(f"{res}:{loc}" for res, loc in seq)


def _log_ion_diff(
    tag: str, pdb_ions: set[tuple[str, str]], pdbqt_ions: set[tuple[str, str]]
) -> None:
    missing = pdb_ions - pdbqt_ions
    kept = pdb_ions & pdbqt_ions
    logging.info(
        "[ions.diff.pdb↔pdbqt.%s] kept=%d stripped=%d detail=%s",
        tag,
        len(kept),
        len(missing),
        _format_ion_pairs(pdb_ions),
    )
    if pdb_ions or missing:
        logging.warning(
            "[ion diff] present_pdb=%s missing_in_pdbqt=%s",
            sorted(pdb_ions),
            sorted(missing),
        )
    retained = {
        str(tok).strip().upper()
        for tok in getattr(ALIASES, "elem_tokens_canonical", set())
        if str(tok).strip()
    }
    if not retained:
        retained = set(_ELEM_CANON)
    lost_retained = sorted(
        [item for item in missing if _normalize_resname(item[0]) in retained]
    )
    if lost_retained:
        logging.warning("[ion lost] %s", lost_retained)


def _ion_candidate_tokens(rules_obj=ALIASES) -> set[str]:
    tokens: set[str] = set()
    canonical = getattr(rules_obj, "elem_tokens_canonical", None) or set()
    for tok in canonical:
        if tok is None:
            continue
        text = str(tok).strip()
        if text:
            tokens.add(text.upper())
    alias_sets = getattr(rules_obj, "alias_sets", None)
    if alias_sets and getattr(alias_sets, "elem_tokens_canonical", None):
        for tok in getattr(alias_sets, "elem_tokens_canonical", set()):
            if tok is None:
                continue
            text = str(tok).strip()
            if text:
                tokens.add(text.upper())
    return tokens


def _scan_metal_map(path: Union[str, Path], rules=ALIASES) -> dict[str, int]:
    file_path = Path(path)
    if not file_path.exists():
        return {}
    tokens = _ion_candidate_tokens(rules)
    counts: Counter[str] = Counter()
    try:
        with file_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.startswith(("HETATM", "ATOM  ")):
                    continue
                resname = line[17:20].strip().upper()
                if not resname:
                    continue
                if resname in tokens:
                    counts[resname] += 1
    except Exception as exc:
        logging.warning("[ions.probe] stage=scan_fail file=%s err=%s", file_path, exc)
        return dict(counts)
    return dict(counts)


def _format_diff_map(data: dict[str, int]) -> str:
    if not data:
        return "none"
    filtered = [(res, count) for res, count in data.items() if count > 0]
    if not filtered:
        return "none"
    items = sorted(filtered, key=lambda kv: (-kv[1], kv[0]))
    return ",".join(f"{res}:{cnt}" for res, cnt in items)


def _log_pdb_pdbqt_counts_diff(
    pdb_path: Union[str, Path],
    pdbqt_path: Union[str, Path],
    *,
    tool: str | None = None,
) -> None:
    pdb_file = Path(pdb_path)
    pdbqt_file = Path(pdbqt_path)
    if not pdb_file.exists() or not pdbqt_file.exists():
        logging.warning(
            "[iondiff.pdb_pdbqt] action=skip reason=missing tool=%s file_pdb=%s exists_pdb=%s file_pdbqt=%s exists_pdbqt=%s",
            tool or "unknown",
            pdb_file,
            pdb_file.exists(),
            pdbqt_file,
            pdbqt_file.exists(),
        )
        return

    pdb_counts = Counter(_scan_metal_map(pdb_file, ALIASES))
    pdbqt_counts = Counter(_scan_metal_map(pdbqt_file, ALIASES))

    kept: dict[str, int] = {}
    lost: dict[str, int] = {}
    gained: dict[str, int] = {}

    for resname, count in pdb_counts.items():
        matched = min(count, pdbqt_counts.get(resname, 0))
        if matched > 0:
            kept[resname] = matched
        delta = count - pdbqt_counts.get(resname, 0)
        if delta > 0:
            lost[resname] = delta

    for resname, count in pdbqt_counts.items():
        delta = count - pdb_counts.get(resname, 0)
        if delta > 0:
            gained[resname] = delta

    logging.info(
        "[iondiff.pdb_pdbqt] tool=%s kept=%s lost=%s added=%s",
        tool or "unknown",
        _format_diff_map(kept),
        _format_diff_map(lost),
        _format_diff_map(gained),
    )


def _normalize_ion_policy(cfg: Optional[dict]) -> str:
    raw = "by_variant"
    if cfg is not None:
        raw = (
            str(cfg.get("ION_STRIP_POLICY", "by_variant")).strip().lower()
            or "by_variant"
        )
    if raw not in {"by_variant", "always_strip", "never_strip"}:
        logging.warning("[ions.cfg] unsupported_policy=%s fallback=by_variant", raw)
        return "by_variant"
    return raw


def _distance_from_center(
    line: str, center: Optional[tuple[float, float, float]]
) -> float | None:
    if center is None:
        return None
    xyz = line_xyz(line)
    if xyz is None:
        return None
    x, y, z = xyz
    dx = x - center[0]
    dy = y - center[1]
    dz = z - center[2]
    return (dx * dx + dy * dy + dz * dz) ** 0.5


def _first_last_lines(p: Path, n: int = 50) -> tuple[list, list]:
    try:
        lines = (p.read_text(encoding="utf-8", errors="ignore")).splitlines()
        return lines[:n], lines[-n:]
    except Exception:
        return [], []


def _ok_receptor_file(p: Path) -> bool:
    return p.exists() and p.stat().st_size > 0


def _classify_and_rename_histidines(
    pdb_in: Union[str, Path], pdb_out: Union[str, Path]
) -> None:
    """
    Inspect each HIS residue's side-chain hydrogens and rename:
      - HD1 only  -> HID
      - HE2 only  -> HIE
      - both      -> HIP
      - neither   -> keep HIS (caller may rewrite later)
    """
    pdb_in, pdb_out = str(pdb_in), str(pdb_out)

    residues = defaultdict(list)
    other_lines = []
    with open(pdb_in, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM  ", "HETATM")):
                chain = ln[21]
                resseq = ln[22:26]
                icode = ln[26]
                resn = ln[17:20]
                key = (chain, resseq, icode, resn)
                residues[key].append(ln)
            else:
                other_lines.append(ln)

    out_lines = []
    for key, atoms in residues.items():
        chain, resseq, icode, resn = key
        resname = resn.strip()
        if resname != "HIS":
            out_lines.extend(atoms)
            continue

        names = {ln[12:16].strip().upper() for ln in atoms}
        has_hd1 = "HD1" in names  # proton on ND1
        has_he2 = "HE2" in names  # proton on NE2

        if has_hd1 and has_he2:
            new = "HIP"
        elif has_hd1 and not has_he2:
            new = "HID"
        elif has_he2 and not has_hd1:
            new = "HIE"
        else:
            new = None

        if new:
            for ln in atoms:
                if ln.startswith(("ATOM  ", "HETATM")):
                    out_lines.append(ln[:17] + f"{new:>3}" + ln[20:])
                else:
                    out_lines.append(ln)
        else:
            out_lines.extend(atoms)

    with open(pdb_out, "w", encoding="utf-8") as w:
        for ln in other_lines:
            w.write(ln)
        for ln in out_lines:
            w.write(ln)


def _rewrite_his_default(
    pdb_in: Union[str, Path], pdb_out: Union[str, Path], default: str = "HIE"
) -> None:
    pdb_in, pdb_out = str(pdb_in), str(pdb_out)
    default = default.upper()
    assert default in {"HIE", "HID", "HIP"}
    with (
        open(pdb_in, "r", encoding="utf-8", errors="ignore") as f,
        open(pdb_out, "w", encoding="utf-8") as w,
    ):
        for ln in f:
            if ln.startswith(("ATOM  ", "HETATM")) and ln[17:20] == "HIS":
                w.write(ln[:17] + default + ln[20:])
            else:
                w.write(ln)


def _with_his_mapping(cmd: list[str], mapping: str) -> list[str]:
    if not mapping:
        return cmd
    return cmd + ["-n", mapping]


def _with_default_altloc(cmd: list[str], default_altloc: str) -> list[str]:
    token = default_altloc.strip()
    if not token or "--default_altloc" in cmd:
        return cmd
    return cmd + ["--default_altloc", token]


def _drop_free_ions_for_meeko(
    pdb_in: str | Path,
    pdb_out: str | Path,
    banlist: set[str] | None = None,
    *,
    cfg: Optional[dict] = None,
    variant: Optional[str] = None,
    pocket_center: Optional[tuple[float, float, float]] = None,
) -> int:
    """
    Remove HETATM entries for simple ions that Meeko chokes on (e.g., Na, K, Li).
    Writes to pdb_out. Returns #atoms dropped.
    Never drops ions explicitly retained by YAML (retain_in_receptor_resnames) that are elemental tokens.
    Set MEEKO_DROP_VERBOSE=1 to log each dropped residue position.
    """
    cfg_obj: Optional[dict] = None
    if isinstance(cfg, dict):
        cfg_obj = cfg
    elif isinstance(prep_utils.config, dict):
        cfg_obj = prep_utils.config

    allow_tokens, _ = _load_retain_allowlist(cfg_obj)
    allow_set = {str(tok).strip().upper() for tok in allow_tokens if str(tok).strip()}
    retain_ions = {r for r in _RETAIN_VARIANT if _is_element_token(r)}
    base_ban = set(banlist) if banlist else _MEEKO_DROP_IONS

    policy = _normalize_ion_policy(cfg_obj)
    variant_token = _resolve_variant_token(cfg_obj, variant)
    variant_label = variant_token or "legacy"

    radius_cfg = 0.0
    if cfg_obj is not None:
        try:
            radius_cfg = float(cfg_obj.get("HOLO_SALT_STRIP_RADIUS", 0.0) or 0.0)
        except Exception:
            radius_cfg = 0.0
    radius_cfg = max(0.0, radius_cfg)
    radius = radius_cfg if (policy == "by_variant" and variant_token == "HOLO") else 0.0
    radius_term = f"{radius:.2f}" if radius > 0.0 else "none"

    drop_metals = False
    drop_salts = False
    if policy == "always_strip":
        drop_metals = True
        drop_salts = True
    elif policy == "by_variant":
        if variant_token == "HOLO":
            drop_metals = False
            drop_salts = radius > 0.0
        else:
            drop_metals = True
            drop_salts = True
    elif policy == "never_strip":
        drop_metals = False
        drop_salts = False

    if (
        drop_salts
        and radius > 0.0
        and pocket_center is None
        and variant_token == "HOLO"
    ):
        logging.warning(
            "[ions] holo_salt_radius_set_but_no_center stage=meeko action=keep_salts radius=%.2f",
            radius,
        )
        drop_salts = False

    logging.info(
        "[ions.policy] stage=meeko variant=%s drop_metals=%s drop_salts=%s radius=%s allowlist=%d",
        variant_label,
        drop_metals,
        drop_salts,
        radius_term,
        len(allow_set),
    )

    verbose = os.environ.get("MEEKO_DROP_VERBOSE", "0") not in ("0", "false", "False")

    dropped = 0
    kept_lines: list[str] = []
    warn_missing_center = False
    with open(pdb_in, "r", encoding="utf-8", errors="ignore") as fh:
        for ln in fh:
            if not ln.startswith("HETATM"):
                kept_lines.append(ln)
                continue
            res = ln[17:20].strip().upper()
            elem = (ln[76:78].strip() or res).upper()
            token = res or elem
            if token in allow_set or elem in allow_set:
                kept_lines.append(ln)
                continue

            is_metal = res in _METAL_RESNAMES or elem in _METAL_RESNAMES
            is_salt = (
                res in _SALT_RESNAMES
                or elem in _SALT_RESNAMES
                or res in base_ban
                or elem in base_ban
            )

            remove = False
            reason = ""
            if policy == "never_strip":
                remove = False
            elif is_metal and drop_metals and (token not in retain_ions):
                remove = True
                reason = "meeko_strip_metal"
            elif is_salt and drop_salts:
                if radius > 0.0:
                    dist = _distance_from_center(ln, pocket_center)
                    if dist is None:
                        warn_missing_center = True
                        remove = False
                    elif dist >= radius:
                        remove = True
                        reason = "meeko_salt_far"
                elif (res in base_ban or elem in base_ban) and policy != "never_strip":
                    remove = True
                    reason = "meeko_salt_policy"
            elif (res in base_ban or elem in base_ban) and policy == "always_strip":
                remove = True
                reason = "meeko_policy"

            if remove:
                dropped += 1
                if verbose:
                    resi = ln[22:26].strip()
                    chain = ln[21]
                    logging.info(
                        "[meeko drop] res=%s chain=%s resi=%s elem=%s reason=%s",
                        res,
                        chain,
                        resi,
                        elem or res,
                        reason or "policy",
                    )
                continue

            kept_lines.append(ln)

    if warn_missing_center and radius > 0.0 and variant_token == "HOLO":
        logging.warning(
            "[ions] holo_salt_radius_set_but_no_center stage=meeko action=keep_salts radius=%.2f",
            radius,
        )

    Path(pdb_out).write_text("".join(kept_lines), encoding="utf-8")
    logging.info(
        "[ions.touch] stage=meeko variant=%s dropped=%d kept=%d file=%s",
        variant_label,
        dropped,
        len(kept_lines),
        pdb_out,
    )
    if dropped:
        logging.warning(
            "Meeko pre-sanitize: dropped %d monoatomics (ban=%s).",
            dropped,
            ",".join(sorted(base_ban - retain_ions)),
        )
    return dropped


def _clean_receptor_pdbqt(pdbqt_path: Union[str, Path]) -> None:
    """
    Drop Reduce USER MOD lines and strip bogus AutoDock 'std' suffixes on
    receptor PDBQT atom records. Safe no-op on missing files.
    """
    p = Path(pdbqt_path)
    if not p.exists():
        return
    try:
        lines = p.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception as exc:
        logging.warning("[receptor.clean] read failed for %s: %s", p, exc)
        return

    cleaned: list[str] = []
    user_removed = 0
    std_fixed = 0
    std_pattern = re.compile(r"\s(?:std|new)\s*$", re.IGNORECASE)

    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith("USER  MOD"):
            user_removed += 1
            continue

        if line.startswith(("ATOM", "HETATM")):
            new_line, n = std_pattern.subn("", line)
            if n:
                std_fixed += n
            cleaned.append(new_line)
        else:
            cleaned.append(line)

    try:
        p.write_text("\n".join(cleaned) + "\n", encoding="utf-8", errors="ignore")
    except Exception as exc:
        logging.warning("[receptor.clean] write failed for %s: %s", p, exc)
        return

    if user_removed or std_fixed:
        logging.info(
            "[receptor.clean] cleaned receptor PDBQT %s: user_mod_removed=%d std_suffix_fixed=%d",
            str(p),
            user_removed,
            std_fixed,
        )


def _open_meeko_fallback_candidates(output_pdbqt: Union[str, Path]) -> list[Path]:
    protein_root = Path(output_pdbqt).resolve().parent.parent
    work_dir = protein_root / "work"
    pdb_id = protein_root.name
    names = (
        f"{pdb_id}_openmm_repaired.pdb",
        f"{pdb_id}_validated.pdb",
        f"{pdb_id}_elemfix.pdb",
        f"{pdb_id}_stripped.pdb",
    )
    return [work_dir / name for name in names if (work_dir / name).exists()]


def _prepare_meeko_fallback_input(
    candidate: Path,
    work_dir: Path,
    cfg: dict,
    variant_token: Optional[str],
    drop_policy: bool,
) -> Path | None:
    prepared = work_dir / f"meeko_fallback_{candidate.stem}.pdb"
    try:
        _classify_and_rename_histidines(str(candidate), str(prepared))
    except Exception as exc:
        logging.warning(
            "[meeko.fallback] HIS classify failed candidate=%s err=%s",
            candidate,
            exc,
        )
        try:
            shutil.copy2(candidate, prepared)
        except Exception:
            return None
    if drop_policy:
        _drop_free_ions_for_meeko(
            prepared,
            prepared,
            cfg=cfg,
            variant=variant_token,
            pocket_center=None,
        )
    try:
        return Path(_meeko_preflight_or_fail(prepared, work_dir))
    except Exception as exc:
        logging.warning(
            "[meeko.fallback] preflight failed candidate=%s err=%s",
            candidate,
            exc,
        )
        return None


def _try_meeko_open_fallbacks(
    output_pdbqt: str,
    cfg: dict,
    variant_token: Optional[str],
    drop_policy: bool,
    base_meeko: list[str] | None,
) -> tuple[bool, str]:
    if not base_meeko:
        return False, ""
    work_dir = Path(output_pdbqt).resolve().parent.parent / "work"
    for candidate in _open_meeko_fallback_candidates(output_pdbqt):
        prepared = _prepare_meeko_fallback_input(
            candidate,
            work_dir,
            cfg,
            variant_token,
            drop_policy,
        )
        if prepared is None:
            continue
        label_base = f"meeko_fallback_{candidate.stem}"
        attempts = (
            ("modern", build_meeko_modern_cmd(base_meeko, prepared, output_pdbqt)),
            ("legacy", build_meeko_legacy_cmd(base_meeko, prepared, output_pdbqt)),
        )
        for mode, cmd in attempts:
            cp = run_subprocess_capture(cmd)
            _persist_subproc(
                f"{label_base}_{mode}",
                cmd,
                cp,
                work_dir,
                Path(output_pdbqt),
            )
            if cp.returncode == 0 and _ok_receptor_file(Path(output_pdbqt)):
                _clean_receptor_pdbqt(output_pdbqt)
                return True, candidate.name
    return False, ""


def _run_meeko_output_attempts(
    attempts: Sequence[tuple[str, list[str]]],
    *,
    log_prefix: str,
    work_dir: Path,
    output_path: Path,
) -> bool:
    for mode, cmd in attempts:
        cp = run_subprocess_capture(cmd)
        _persist_subproc(f"{log_prefix}_{mode}", cmd, cp, work_dir, output_path)
        if cp.returncode == 0 and _ok_receptor_file(output_path):
            _clean_receptor_pdbqt(str(output_path))
            return True
    return False


def _try_meeko_normalized_fallback(
    input_pdb: str,
    output_pdbqt: str,
    base_meeko: list[str] | None,
    his_mapping: str = "",
    *,
    drop_waters: bool,
) -> tuple[bool, str]:
    if not base_meeko:
        return False, ""
    output_path = Path(output_pdbqt)
    work_dir = output_path.resolve().parent.parent / "work"
    label = "normalized_dewatered_noh_retry" if drop_waters else "normalized_noh_retry"
    normalized = work_dir / f"meeko_fallback_{output_path.stem}_{label}.pdb"
    try:
        write_meeko_normalized_copy(
            Path(input_pdb),
            normalized,
            drop_waters=drop_waters,
            drop_hydrogens=True,
        )
    except Exception as exc:
        logging.warning("[meeko.%s] prepare_failed err=%s", label, exc)
        return False, ""
    attempts = (
        (
            "modern",
            _with_his_mapping(
                build_meeko_modern_cmd(base_meeko, normalized, output_pdbqt),
                his_mapping,
            ),
        ),
        (
            "legacy",
            _with_his_mapping(
                build_meeko_legacy_cmd(base_meeko, normalized, output_pdbqt),
                his_mapping,
            ),
        ),
    )
    if _run_meeko_output_attempts(
        attempts,
        log_prefix=f"meeko_{label}",
        work_dir=work_dir,
        output_path=output_path,
    ):
        return True, label
    return False, label


def _try_meeko_clash_pruned_fallback(
    input_pdb: str,
    output_pdbqt: str,
    base_meeko: list[str] | None,
    his_mapping: str = "",
    *,
    allow_bad_res: bool = False,
) -> tuple[bool, str]:
    if not base_meeko:
        return False, ""
    output_path = Path(output_pdbqt)
    work_dir = output_path.resolve().parent.parent / "work"
    pruned = work_dir / f"meeko_fallback_{output_path.stem}_clash_pruned.pdb"
    try:
        dropped = write_meeko_clash_pruned_copy(Path(input_pdb), pruned)
    except Exception as exc:
        logging.warning("[meeko.clash_pruned] prepare_failed err=%s", exc)
        return False, ""
    if dropped <= 0 or not pruned.exists():
        return False, ""
    clean_label = meeko_pruned_stage_label(pruned, dropped=dropped)
    clean_attempts = (
        (
            "modern",
            _with_his_mapping(
                build_meeko_modern_cmd(base_meeko, pruned, output_pdbqt),
                his_mapping,
            ),
        ),
        (
            "legacy",
            _with_his_mapping(
                build_meeko_legacy_cmd(base_meeko, pruned, output_pdbqt),
                his_mapping,
            ),
        ),
    )
    if _run_meeko_output_attempts(
        clean_attempts,
        log_prefix="meeko_clash_pruned",
        work_dir=work_dir,
        output_path=output_path,
    ):
        return True, clean_label
    if not allow_bad_res:
        return False, clean_label
    label = meeko_pruned_stage_label(pruned, dropped=dropped, allow_bad=True)
    allow_bad_attempts = (
        ("modern_allow_bad", clean_attempts[0][1] + ["-a"]),
        ("legacy_allow_bad", clean_attempts[1][1] + ["-a"]),
    )
    if _run_meeko_output_attempts(
        allow_bad_attempts,
        log_prefix="meeko_clash_pruned",
        work_dir=work_dir,
        output_path=output_path,
    ):
        return True, label
    return False, label


def _run_meeko_output_attempts_with_error(
    attempts: Sequence[tuple[str, list[str]]],
    *,
    log_prefix: str,
    work_dir: Path,
    output_path: Path,
) -> tuple[bool, str]:
    last_error = ""
    for mode, cmd in attempts:
        cp = run_subprocess_capture(cmd)
        _persist_subproc(f"{log_prefix}_{mode}", cmd, cp, work_dir, output_path)
        if cp.returncode == 0 and _ok_receptor_file(output_path):
            _clean_receptor_pdbqt(str(output_path))
            return True, ""
        last_error = str(getattr(cp, "stderr", "") or getattr(cp, "stdout", "") or last_error)
    return False, last_error


def _try_meeko_safe_receptor_fallback(
    input_pdb: str,
    output_pdbqt: str,
    base_meeko: list[str] | None,
    his_mapping: str = "",
    *,
    error_context: str = "",
) -> tuple[bool, str]:
    if not base_meeko:
        return False, ""
    output_path = Path(output_pdbqt)
    work_dir = output_path.resolve().parent.parent / "work"
    residue_keys = parse_unmatched_residue_keys(error_context)
    last_error = error_context
    for attempt in range(33):
        suffix = "protein_only" if attempt == 0 else f"template_pruned_{attempt}"
        safe_pdb = work_dir / f"meeko_fallback_{output_path.stem}_{suffix}.pdb"
        try:
            summary = write_meeko_safe_receptor_copy(
                Path(input_pdb),
                safe_pdb,
                drop_residue_keys=residue_keys,
            )
        except Exception as exc:
            logging.warning("[meeko.safe_export] prepare_failed err=%s", exc)
            return False, ""
        attempts = (
            (
                "modern",
                _with_his_mapping(
                    build_meeko_modern_cmd(base_meeko, safe_pdb, output_pdbqt),
                    his_mapping,
                ),
            ),
            (
                "legacy",
                _with_his_mapping(
                    build_meeko_legacy_cmd(base_meeko, safe_pdb, output_pdbqt),
                    his_mapping,
                ),
            ),
        )
        ok, last_error = _run_meeko_output_attempts_with_error(
            attempts,
            log_prefix=f"meeko_safe_export_{suffix}",
            work_dir=work_dir,
            output_path=output_path,
        )
        if ok:
            return True, meeko_safe_stage_label(summary)
        pruned_pdb = work_dir / f"meeko_fallback_{output_path.stem}_{suffix}_clash_pruned.pdb"
        dropped = write_meeko_clash_pruned_copy(safe_pdb, pruned_pdb)
        if dropped:
            pruned_attempts = (
                (
                    "modern",
                    _with_his_mapping(
                        build_meeko_modern_cmd(base_meeko, pruned_pdb, output_pdbqt),
                        his_mapping,
                    ),
                ),
                (
                    "legacy",
                    _with_his_mapping(
                        build_meeko_legacy_cmd(base_meeko, pruned_pdb, output_pdbqt),
                        his_mapping,
                    ),
                ),
            )
            ok, last_error = _run_meeko_output_attempts_with_error(
                pruned_attempts,
                log_prefix=f"meeko_safe_export_{suffix}_clash_pruned",
                work_dir=work_dir,
                output_path=output_path,
            )
            if ok:
                label = (
                    f"{meeko_safe_stage_label(summary)}:"
                    f"clash_pruned_retry:dropped_atoms={dropped}"
                )
                return True, label
        next_keys = parse_unmatched_residue_keys(last_error)
        if not next_keys or next_keys.issubset(residue_keys):
            break
        residue_keys |= next_keys
    return False, "protein_only_sanitized_retry"


def _write_meeko_forced_sanitized_copy(
    input_pdb: Path,
    output_pdb: Path,
    *,
    drop_oxt: bool,
) -> int:
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    dropped = 0
    dropped_records: list[dict[str, str]] = []
    kept: list[str] = []
    with input_pdb.open("r", encoding="utf-8", errors="ignore") as src:
        for line in src:
            if line.startswith("END"):
                continue
            if line.startswith("HETATM"):
                res = line[17:20].strip().upper()
                elem = (line[76:78].strip() or res).upper()
                if res in _MEEKO_DROP_IONS or elem in _MEEKO_DROP_IONS:
                    dropped += 1
                    dropped_records.append(
                        {
                            "record": line[:6].strip(),
                            "chain": line[21:22].strip() or "-",
                            "resseq": line[22:26].strip(),
                            "resname": res,
                            "atom_name": line[12:16].strip(),
                            "reason": "meeko_salt_ion_strip",
                        }
                    )
                    continue
            if (
                drop_oxt
                and line.startswith("ATOM  ")
                and line[12:16].strip().upper() == "OXT"
            ):
                dropped += 1
                dropped_records.append(
                    {
                        "record": line[:6].strip(),
                        "chain": line[21:22].strip() or "-",
                        "resseq": line[22:26].strip(),
                        "resname": line[17:20].strip(),
                        "atom_name": line[12:16].strip(),
                        "reason": "meeko_terminal_oxt_strip",
                    }
                )
                continue
            kept.append(line)
    output_pdb.write_text("".join(kept).rstrip() + "\nEND\n", encoding="utf-8")
    if not drop_oxt:
        repair_terminal_oxt_geometry_in_pdb(
            output_pdb,
            audit_path=output_pdb.with_suffix(".oxt_repair_audit.json"),
        )
    output_pdb.with_suffix(".drop_audit.json").write_text(
        json.dumps(
            {
                "input_pdb": str(input_pdb),
                "output_pdb": str(output_pdb),
                "dropped_atom_count": dropped,
                "dropped_atoms": dropped_records,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return dropped


def _try_meeko_forced_sanitized_fallback(
    input_pdb: str,
    output_pdbqt: str,
    base_meeko: list[str] | None,
    *,
    label_token: str,
    drop_oxt: bool,
    allow_bad_res: bool = False,
) -> tuple[bool, str]:
    if not base_meeko:
        return False, ""
    output_path = Path(output_pdbqt)
    work_dir = output_path.resolve().parent.parent / "work"
    stripped = work_dir / f"meeko_fallback_{output_path.stem}_{label_token}.pdb"
    try:
        dropped = _write_meeko_forced_sanitized_copy(
            Path(input_pdb),
            stripped,
            drop_oxt=drop_oxt,
        )
    except Exception as exc:
        logging.warning("[meeko.%s] prepare_failed err=%s", label_token, exc)
        return False, ""
    if dropped <= 0:
        return False, ""
    his_tmp = stripped.with_suffix(".his.pdb")
    try:
        _classify_and_rename_histidines(str(stripped), str(his_tmp))
        _rewrite_his_default(his_tmp, stripped, default="HIE")
    except Exception as exc:
        logging.warning("[meeko.%s] his_classify_failed err=%s", label_token, exc)
        try:
            _rewrite_his_default(stripped, his_tmp, default="HIE")
            shutil.move(str(his_tmp), str(stripped))
        except Exception as rewrite_exc:
            logging.warning("[meeko.%s] his_default_failed err=%s", label_token, rewrite_exc)
    try:
        prepared = Path(_meeko_preflight_or_fail(stripped, work_dir))
    except Exception as exc:
        logging.warning("[meeko.%s] preflight_failed err=%s", label_token, exc)
        return False, f"{label_token}:dropped_atoms={dropped}"

    clean_label = f"{label_token}_retry:dropped_atoms={dropped}"
    clean_attempts = (
        ("modern", build_meeko_modern_cmd(base_meeko, prepared, output_pdbqt)),
        ("legacy", build_meeko_legacy_cmd(base_meeko, prepared, output_pdbqt)),
    )
    for mode, cmd in clean_attempts:
        cp = run_subprocess_capture(cmd)
        _persist_subproc(
            f"meeko_{label_token}_{mode}",
            cmd,
            cp,
            work_dir,
            output_path,
        )
        if cp.returncode == 0 and _ok_receptor_file(output_path):
            _clean_receptor_pdbqt(output_pdbqt)
            return True, clean_label
    if not allow_bad_res:
        return False, clean_label
    label = f"{label_token}_allow_bad_retry:dropped_atoms={dropped}"
    for mode, cmd in (
        ("modern_allow_bad", clean_attempts[0][1] + ["-a"]),
        ("legacy_allow_bad", clean_attempts[1][1] + ["-a"]),
    ):
        cp = run_subprocess_capture(cmd)
        _persist_subproc(
            f"meeko_{label_token}_{mode}",
            cmd,
            cp,
            work_dir,
            output_path,
        )
        if cp.returncode == 0 and _ok_receptor_file(output_path):
            _clean_receptor_pdbqt(output_pdbqt)
            return True, label
    return False, label


def _write_meeko_stage_audit(
    *,
    input_pdb: str,
    output_pdbqt: str,
    stage: str,
    ok: bool,
    error: str = "",
) -> None:
    output_path = Path(output_pdbqt).resolve()
    work_dir = output_path.parent.parent / "work"
    review_required = _meeko_stage_requires_review(stage)
    publication_review_required = _meeko_stage_requires_review(
        stage,
        profile="publication",
    )
    payload = {
        "input_pdb": str(input_pdb),
        "output_pdbqt": str(output_pdbqt),
        "stage": stage,
        "ok": ok,
        "review_required": review_required,
        "publication_review_required": publication_review_required,
        "error": error[:1000],
    }
    paths = (work_dir / "meeko_input_stage.json", output_path.with_suffix(".meeko_input_stage.json"))
    for path in paths:
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        except Exception as exc:
            logging.warning("[meeko.stage] sidecar_write_failed path=%s err=%s", path, exc)
    logging.info(
        "[meeko.stage] stage=%s ok=%s review_required=%s",
        stage,
        ok,
        review_required,
    )


def run_prepare_receptor(
    input_pdb: Union[str, Path], output_pdbqt: Union[str, Path], cfg: dict
) -> bool:
    """
    Robust Meeko wrapper with:
      • Early check for HIS tautomer ambiguity on the modern (--read_pdb) call
      • Auto -n mapping for all truly ambiguous HIS (no HD1/HE2)
      • Optional -a (allow_bad_res) retry on template-mismatch
      • Heavy-handed HIS?<default> fallback
    """
    input_pdb = str(input_pdb)
    output_pdbqt = str(output_pdbqt)
    variant_token = _resolve_variant_token(cfg, cfg.get("_CURRENT_VARIANT"))
    variant_label = variant_token or "legacy"
    # Persist all attempt logs here (processed_pdbs/<RUN_ID>/<PDB>/work).
    work_dir = Path(output_pdbqt).resolve().parent.parent / "work"
    # --- Helium preflight on the EXACT file we’re about to feed into Meeko pipeline
    _work_dir = Path(output_pdbqt).resolve().parent.parent / "work"
    skip_meeko = False
    try:
        input_pdb = str(_meeko_preflight_or_fail(input_pdb, _work_dir))
    except RuntimeError as e:
        if str(e) == "helium_preflight_failed":
            logging.error("Helium persists pre-Meeko; receptor prep will fail.")
            skip_meeko = True
        else:
            raise

    drop_policy = _cfg_bool("MEEKO_DROP_FREE_IONS", True)
    drop_list = ",".join(sorted(_MEEKO_DROP_IONS)) if _MEEKO_DROP_IONS else "none"
    logging.info(
        "[meeko.policy] drop_free_ions=%s variant=%s list=%s",
        str(bool(drop_policy)).lower(),
        variant_label,
        drop_list,
    )

    his_default = str(cfg.get("HIS_DEFAULT", "HIE")).upper()
    if his_default not in {"HIE", "HID", "HIP"}:
        his_default = "HIE"

    # Meeko --allow_bad_res ignores unresolved residues; keep it opt-in.
    allow_bad_res = str(cfg.get("MEEKO_ALLOW_BAD_RES", "false")).lower() in (
        "1",
        "true",
        "yes",
    )
    if allow_bad_res:
        logging.warning(
            "[meeko.policy] allow_bad_res=true; receptor export may ignore unresolved residues."
        )
    default_altloc = (cfg.get("MEEKO_DEFAULT_ALTLOC") or "").strip()  # e.g. "A" or ""

    def _meeko_cmd() -> list[str]:
        return build_meeko_base_cmd()

    def _build_his_override_mapping(pdb_path: str) -> str:
        """Return a single -n mapping like 'A:27,A:94=HIE' for HIS with no HD1/HE2."""
        by_res = defaultdict(set)
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if not ln.startswith(("ATOM  ", "HETATM")):
                    continue
                if ln[17:20] != "HIS":
                    continue
                chain = ln[21]
                resi = ln[22:26].strip()
                aname = ln[12:16].strip().upper()
                by_res[(chain, resi)].add(aname)
        targets = []
        for (chain, resi), names in by_res.items():
            if ("HD1" not in names) and ("HE2" not in names):
                try:
                    rnum = int(resi)
                except Exception:
                    rnum = int(resi.strip() or "0")
                targets.append(f"{chain}:{rnum}")
        return (",".join(sorted(targets)) + f"={his_default}") if targets else ""

    # --- Step 0: classify HIS by existing hydrogens (best-case, no coord change)
    work_dir.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile("w", suffix=".pdb", dir=work_dir, delete=False) as tmp1:
        tmp1_path = tmp1.name
    try:
        _classify_and_rename_histidines(
            input_pdb, tmp1_path
        )  # leaves non-diagnostic HIS as HIS
    except Exception as e:
        logging.warning(
            "HIS classify/rename step failed (continuing with original): %s", e
        )
        tmp1_path = input_pdb

    if drop_policy:
        logging.info(
            "[ions.touch] stage=meeko-pre variant=%s action=pre_sanitize file=%s",
            variant_label,
            tmp1_path,
        )
        # Never drop YAML-retained elemental ions (banlist computed inside)
        _drop_free_ions_for_meeko(
            tmp1_path,
            tmp1_path,
            cfg=cfg,
            variant=variant_token,
            pocket_center=None,
        )

    # Preflight again on the HIS-classified file we’re actually handing to Meeko now
    if not skip_meeko:
        try:
            _meeko_preflight_or_fail(
                tmp1_path, Path(output_pdbqt).resolve().parent.parent / "work"
            )
        except RuntimeError as e:
            if str(e) == "helium_preflight_failed":
                logging.error("Helium persists after HIS/ion tweaks; receptor prep will fail.")
                skip_meeko = True
            else:
                raise

    work_dir = Path(output_pdbqt).resolve().parent.parent / "work"

    # Build explicit commands for each attempt
    try:
        base_meeko = _meeko_cmd()
    except FileNotFoundError as e:
        base_meeko = None
        logging.warning("Meeko CLI not found at prebuild stage: %s", e)

    his_mapping = _build_his_override_mapping(tmp1_path)

    # Meeko modern (--read_pdb) uses the HIS-classified, ion-sanitized tmp1_path
    if base_meeko:
        meeko_cmd = _with_default_altloc(
            _with_his_mapping(
                build_meeko_modern_cmd(base_meeko, tmp1_path, output_pdbqt),
                his_mapping,
            ),
            default_altloc,
        )
        # Legacy single-shot form (we still capture rc/file size for logging)
        meeko_cmd_legacy = _with_default_altloc(
            _with_his_mapping(
                build_meeko_legacy_cmd(base_meeko, tmp1_path, output_pdbqt),
                his_mapping,
            ),
            default_altloc,
        )
    else:
        meeko_cmd = None
        meeko_cmd_legacy = None

    # --- Step 1: Modern Meeko attempt
    cp: Any
    if not skip_meeko and meeko_cmd:
        cp = run_subprocess_capture(meeko_cmd)
        _persist_subproc("meeko_modern", meeko_cmd, cp, work_dir, Path(output_pdbqt))
    else:
        # If Meeko is skipped or unavailable, synthesize a "failed" result object
        from types import SimpleNamespace

        cp = SimpleNamespace(
            returncode=127,
            stdout="",
            stderr=("meeko_skipped" if skip_meeko else "meeko_not_found"),
        )

    if cp.returncode == 0 and _ok_receptor_file(Path(output_pdbqt)):
        _clean_receptor_pdbqt(output_pdbqt)
        logging.info("Prepared receptor PDBQT with modern Meeko.")
        _write_meeko_stage_audit(
            input_pdb=input_pdb,
            output_pdbqt=output_pdbqt,
            stage=("modern_his_template_input" if his_mapping else "modern_clean_input"),
            ok=True,
        )
        # --- ION DIFF BLOCK: summarize ions kept vs lost in PDBQT
        try:
            _pdb_ions = _ion_pairs_from_pdb(input_pdb)
            _pdbqt_ions = _ion_pairs_from_pdbqt(output_pdbqt)
            _log_ion_diff("meeko", _pdb_ions, _pdbqt_ions)
        except Exception as _e:
            logging.warning("[ion diff] skipped note=%s", _e)
        try:
            _log_pdb_pdbqt_counts_diff(input_pdb, output_pdbqt, tool="Meeko")
        except Exception as diff_exc:
            logging.warning(
                "[iondiff.pdb_pdbqt] action=skip tool=Meeko reason=%s", diff_exc
            )
        return True

    if Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size == 0:
        head, tail = _first_last_lines(work_dir / "meeko_modern.stderr.txt")
        logging.warning(
            "[meeko_modern] receptor PDBQT is empty. head=%r tail=%r", head, tail
        )

    # Decide if we need an -n retry (histidine tie) based on se0
    se0 = cp.stderr or ""
    his_tie = ("tied for fewest missing h" in se0.lower()) and (
        "hie" in se0.lower() and "hid" in se0.lower()
    )
    mapping = ""
    if his_tie and not his_mapping:
        mapping = _build_his_override_mapping(tmp1_path)
        if not mapping:
            m = re.search(r"residue_key='([A-Za-z]):(\d+)'", se0)
            if m:
                mapping = f"{m.group(1)}:{int(m.group(2))}={his_default}"

    # --- Step 1b: Retry with -n if needed
    if mapping and meeko_cmd:
        meeko_cmd_with_n = _with_his_mapping(meeko_cmd, mapping)
        cp = run_subprocess_capture(meeko_cmd_with_n)
        _persist_subproc(
            "meeko_retry_nmap", meeko_cmd_with_n, cp, work_dir, Path(output_pdbqt)
        )
        if cp.returncode == 0 and _ok_receptor_file(Path(output_pdbqt)):
            _clean_receptor_pdbqt(output_pdbqt)
            logging.info("Prepared receptor after HIS -n mapping.")
            _write_meeko_stage_audit(
                input_pdb=input_pdb,
                output_pdbqt=output_pdbqt,
                stage="his_nmap_retry",
                ok=True,
            )
            # --- ION DIFF BLOCK: summarize ions kept vs lost in PDBQT
            try:
                _pdb_ions = _ion_pairs_from_pdb(input_pdb)
                _pdbqt_ions = _ion_pairs_from_pdbqt(output_pdbqt)
                _log_ion_diff("meeko_retry", _pdb_ions, _pdbqt_ions)
            except Exception as _e:
                logging.warning("[ion diff] skipped note=%s", _e)
            try:
                _log_pdb_pdbqt_counts_diff(input_pdb, output_pdbqt, tool="Meeko-nmap")
            except Exception as diff_exc:
                logging.warning(
                    "[iondiff.pdb_pdbqt] action=skip tool=Meeko-nmap reason=%s",
                    diff_exc,
            )

            return True

        if Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size == 0:
            head, tail = _first_last_lines(work_dir / "meeko_retry_nmap.stderr.txt")
            logging.warning(
                "[meeko_retry_nmap] receptor PDBQT is empty. head=%r tail=%r",
                head,
                tail,
            )

    # --- Step 1c: Retry a Meeko-only normalized receptor before any allow_bad path
    if not skip_meeko and base_meeko:
        for drop_waters in (False, True):
            normalized_ok, normalized_label = _try_meeko_normalized_fallback(
                tmp1_path,
                output_pdbqt,
                base_meeko,
                his_mapping,
                drop_waters=drop_waters,
            )
            if not normalized_ok:
                continue
            logging.info("Prepared receptor with Meeko normalized fallback: %s", normalized_label)
            _write_meeko_stage_audit(
                input_pdb=input_pdb,
                output_pdbqt=output_pdbqt,
                stage=normalized_label,
                ok=True,
            )
            try:
                _pdb_ions = _ion_pairs_from_pdb(input_pdb)
                _pdbqt_ions = _ion_pairs_from_pdbqt(output_pdbqt)
                _log_ion_diff(normalized_label, _pdb_ions, _pdbqt_ions)
            except Exception as _e:
                logging.warning("[ion diff] skipped note=%s", _e)
            try:
                _log_pdb_pdbqt_counts_diff(
                    input_pdb,
                    output_pdbqt,
                    tool=f"Meeko-normalized:{normalized_label}",
                )
            except Exception as diff_exc:
                logging.warning(
                    "[iondiff.pdb_pdbqt] action=skip tool=Meeko-normalized reason=%s",
                    diff_exc,
                )
            return True

    if not skip_meeko and base_meeko:
        safe_ok, safe_label = _try_meeko_safe_receptor_fallback(
            tmp1_path,
            output_pdbqt,
            base_meeko,
            his_mapping,
            error_context=se0,
        )
        if safe_ok:
            logging.info("Prepared receptor with Meeko safe-export fallback: %s", safe_label)
            _write_meeko_stage_audit(
                input_pdb=input_pdb,
                output_pdbqt=output_pdbqt,
                stage=safe_label,
                ok=True,
            )
            try:
                _pdb_ions = _ion_pairs_from_pdb(input_pdb)
                _pdbqt_ions = _ion_pairs_from_pdbqt(output_pdbqt)
                _log_ion_diff("meeko_safe_export", _pdb_ions, _pdbqt_ions)
            except Exception as _e:
                logging.warning("[ion diff] skipped note=%s", _e)
            try:
                _log_pdb_pdbqt_counts_diff(
                    input_pdb,
                    output_pdbqt,
                    tool=f"Meeko-safe-export:{safe_label}",
                )
            except Exception as diff_exc:
                logging.warning(
                    "[iondiff.pdb_pdbqt] action=skip tool=Meeko-safe-export reason=%s",
                    diff_exc,
                )
            return True

    # --- Step 1d: Retry with -a allow_bad_res (and default_altloc if hinted)
    templ_fail = (
        ("no template matched for residue_key" in se0.lower())
        or ("template matching failed" in se0.lower())
        or ("template matched failed" in se0.lower())
    )
    extra = []
    if (
        meeko_cmd
        and allow_bad_res
        and (
            templ_fail
            or "allow_bad_res" in se0.lower()
            or "recommendations" in se0.lower()
        )
    ):
        extra = ["-a"]
        if (not default_altloc) and (
            "default_altloc" in se0.lower() or "altloc" in se0.lower()
        ):
            default_altloc = "A"
            logging.warning("Assuming --default_altloc A based on Meeko hint.")
        if default_altloc and "--default_altloc" not in meeko_cmd:
            extra += ["--default_altloc", default_altloc]
        meeko_cmd_allow_bad = meeko_cmd + extra
        cp = run_subprocess_capture(meeko_cmd_allow_bad)
        _persist_subproc(
            "meeko_retry_allow_bad_res",
            meeko_cmd_allow_bad,
            cp,
            work_dir,
            Path(output_pdbqt),
        )
        if cp.returncode == 0 and _ok_receptor_file(Path(output_pdbqt)):
            _clean_receptor_pdbqt(output_pdbqt)
            logging.info("Prepared receptor after -a allow_bad_res.")
            _write_meeko_stage_audit(
                input_pdb=input_pdb,
                output_pdbqt=output_pdbqt,
                stage="allow_bad_res_retry",
                ok=True,
            )
            # --- ION DIFF BLOCK: summarize ions kept vs lost in PDBQT
            try:
                _pdb_ions = _ion_pairs_from_pdb(input_pdb)
                _pdbqt_ions = _ion_pairs_from_pdbqt(output_pdbqt)
                _log_ion_diff("meeko_allow_bad_res", _pdb_ions, _pdbqt_ions)
            except Exception as _e:
                logging.warning("[ion diff] skipped note=%s", _e)
            try:
                _log_pdb_pdbqt_counts_diff(
                    input_pdb, output_pdbqt, tool="Meeko-allow_bad_res"
                )
            except Exception as diff_exc:
                logging.warning(
                    "[iondiff.pdb_pdbqt] action=skip tool=Meeko-allow_bad_res reason=%s",
                    diff_exc,
                )

            return True

        if Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size == 0:
            head, tail = _first_last_lines(
                work_dir / "meeko_retry_allow_bad_res.stderr.txt"
            )
            logging.warning(
                "[meeko_retry_allow_bad_res] receptor PDBQT is empty. head=%r tail=%r",
                head,
                tail,
            )
    # --- Step 2: Legacy Meeko
    if not skip_meeko and meeko_cmd_legacy:
        cp = run_subprocess_capture(meeko_cmd_legacy)
        _persist_subproc(
            "meeko_legacy", meeko_cmd_legacy, cp, work_dir, Path(output_pdbqt)
        )
        if cp.returncode == 0 and _ok_receptor_file(Path(output_pdbqt)):
            _clean_receptor_pdbqt(output_pdbqt)
            logging.info("Prepared receptor with legacy Meeko.")
            _write_meeko_stage_audit(
                input_pdb=input_pdb,
                output_pdbqt=output_pdbqt,
                stage="legacy_cli_retry",
                ok=True,
            )
            # --- ION DIFF BLOCK: summarize ions kept vs lost in PDBQT
            try:
                _pdb_ions = _ion_pairs_from_pdb(input_pdb)
                _pdbqt_ions = _ion_pairs_from_pdbqt(output_pdbqt)
                _log_ion_diff("meeko_legacy", _pdb_ions, _pdbqt_ions)
            except Exception as _e:
                logging.warning("[ion diff] skipped note=%s", _e)
            try:
                _log_pdb_pdbqt_counts_diff(input_pdb, output_pdbqt, tool="Meeko-legacy")
            except Exception as diff_exc:
                logging.warning(
                    "[iondiff.pdb_pdbqt] action=skip tool=Meeko-legacy reason=%s",
                    diff_exc,
                )

            return True

        if Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size == 0:
            head, tail = _first_last_lines(work_dir / "meeko_legacy.stderr.txt")
            logging.warning(
                "[meeko_legacy] receptor PDBQT is empty. head=%r tail=%r", head, tail
            )

    fallback_ok, fallback_name = _try_meeko_open_fallbacks(
        output_pdbqt,
        cfg,
        variant_token,
        drop_policy,
        base_meeko,
    )
    if fallback_ok:
        logging.info("Prepared receptor with Meeko open fallback: %s", fallback_name)
        _write_meeko_stage_audit(
            input_pdb=input_pdb,
            output_pdbqt=output_pdbqt,
            stage=f"open_fallback:{fallback_name}",
            ok=True,
        )
        try:
            _pdb_ions = _ion_pairs_from_pdb(input_pdb)
            _pdbqt_ions = _ion_pairs_from_pdbqt(output_pdbqt)
            _log_ion_diff("meeko_open_fallback", _pdb_ions, _pdbqt_ions)
        except Exception as _e:
            logging.warning("[ion diff] skipped note=%s", _e)
        try:
            _log_pdb_pdbqt_counts_diff(
                input_pdb,
                output_pdbqt,
                tool=f"Meeko-open-fallback:{fallback_name}",
            )
        except Exception as diff_exc:
            logging.warning(
                "[iondiff.pdb_pdbqt] action=skip tool=Meeko-open-fallback reason=%s",
                diff_exc,
            )
        return True

    for label_token, drop_oxt in (("forced_salt_strip", False),):
        salt_ok, salt_label = _try_meeko_forced_sanitized_fallback(
            input_pdb,
            output_pdbqt,
            base_meeko,
            label_token=label_token,
            drop_oxt=drop_oxt,
            allow_bad_res=allow_bad_res,
        )
        if not salt_ok:
            continue
        logging.info(
            "Prepared receptor with Meeko forced sanitize fallback: %s",
            salt_label,
        )
        _write_meeko_stage_audit(
            input_pdb=input_pdb,
            output_pdbqt=output_pdbqt,
            stage=salt_label,
            ok=True,
        )
        try:
            _pdb_ions = _ion_pairs_from_pdb(input_pdb)
            _pdbqt_ions = _ion_pairs_from_pdbqt(output_pdbqt)
            _log_ion_diff(label_token, _pdb_ions, _pdbqt_ions)
        except Exception as _e:
            logging.warning("[ion diff] skipped note=%s", _e)
        try:
            _log_pdb_pdbqt_counts_diff(
                input_pdb,
                output_pdbqt,
                tool=f"Meeko-forced-sanitize:{salt_label}",
            )
        except Exception as diff_exc:
            logging.warning(
                "[iondiff.pdb_pdbqt] action=skip tool=Meeko-forced-sanitize reason=%s",
                diff_exc,
            )
        return True

    clash_ok, clash_label = _try_meeko_clash_pruned_fallback(
        input_pdb,
        output_pdbqt,
        base_meeko,
        his_mapping,
        allow_bad_res=allow_bad_res,
    )
    if clash_ok:
        logging.info("Prepared receptor with Meeko clash-pruned fallback: %s", clash_label)
        _write_meeko_stage_audit(
            input_pdb=input_pdb,
            output_pdbqt=output_pdbqt,
            stage=clash_label,
            ok=True,
        )
        try:
            _pdb_ions = _ion_pairs_from_pdb(input_pdb)
            _pdbqt_ions = _ion_pairs_from_pdbqt(output_pdbqt)
            _log_ion_diff("meeko_clash_pruned", _pdb_ions, _pdbqt_ions)
        except Exception as _e:
            logging.warning("[ion diff] skipped note=%s", _e)
        try:
            _log_pdb_pdbqt_counts_diff(
                input_pdb,
                output_pdbqt,
                tool=f"Meeko-clash-pruned:{clash_label}",
            )
        except Exception as diff_exc:
            logging.warning(
                "[iondiff.pdb_pdbqt] action=skip tool=Meeko-clash-pruned reason=%s",
                diff_exc,
            )
        return True

    if _ok_receptor_file(Path(output_pdbqt)):
        _clean_receptor_pdbqt(output_pdbqt)
        logging.warning(
            "[receptor] proceeding with existing receptor PDBQT despite prep errors."
        )
        _write_meeko_stage_audit(
            input_pdb=input_pdb,
            output_pdbqt=output_pdbqt,
            stage="existing_receptor_after_errors",
            ok=True,
            error=str(getattr(cp, "stderr", "") or getattr(cp, "stdout", "")),
        )
        return True
    _write_meeko_stage_audit(
        input_pdb=input_pdb,
        output_pdbqt=output_pdbqt,
        stage="failed",
        ok=False,
        error=str(getattr(cp, "stderr", "") or getattr(cp, "stdout", "")),
    )
    return False
