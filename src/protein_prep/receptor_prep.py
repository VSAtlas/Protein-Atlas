"""Receptor PDBQT preparation (Meeko/ADT) and related cleanup helpers."""

from __future__ import annotations

import logging
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Iterable, List, Optional, Sequence, Union

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
from protein_prep.prep_utils import (
    _cfg,
    _cfg_bool,
    _load_retain_allowlist,
    _persist_subproc,
    _resolve_variant_token,
)
from protein_prep.tool_runners import (
    build_adt_prepare_receptor_cmd,
    build_meeko_base_cmd,
    build_meeko_legacy_cmd,
    build_meeko_modern_cmd,
    run_subprocess_capture,
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
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except Exception:
        return None
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


def _strip_reduce_user_lines(pdbqt_path: Union[str, Path]) -> None:
    """Back-compat shim to the unified receptor cleaner."""
    _clean_receptor_pdbqt(pdbqt_path)


MGLTOOLS_PYTHON = _cfg("MGLTOOLS_PYTHON", "")
PREPARE_RECEPTOR_SCRIPT = _cfg("PREPARE_RECEPTOR_SCRIPT", "")


def run_prepare_receptor(
    input_pdb: Union[str, Path], output_pdbqt: Union[str, Path], cfg: dict
) -> bool:
    """
    Robust Meeko/ADT wrapper with:
      • Early check for HIS tautomer ambiguity on the modern (--read_pdb) call
      • Auto -n mapping for all truly ambiguous HIS (no HD1/HE2)
      • Optional -a (allow_bad_res) retry on template-mismatch
      • Heavy-handed HIS?<default> fallback
      • ADT fallback
    """
    input_pdb = str(input_pdb)
    output_pdbqt = str(output_pdbqt)
    variant_token = _resolve_variant_token(cfg, cfg.get("_CURRENT_VARIANT"))
    variant_label = variant_token or "legacy"
    # Persist all attempt logs here (processed_pdbs/<PDB>/work) NOTE DIFF FROM _WORK_DIR
    work_dir = Path(output_pdbqt).resolve().parent.parent / "work"
    # --- Helium preflight on the EXACT file we’re about to feed into Meeko pipeline
    _work_dir = Path(output_pdbqt).resolve().parent.parent / "work"
    skip_meeko = False
    try:
        input_pdb = str(_meeko_preflight_or_fail(input_pdb, _work_dir))
    except RuntimeError as e:
        if str(e) == "helium_preflight_failed":
            logging.error(
                "Helium persists pre-Meeko; skipping Meeko and trying ADT fallback."
            )
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

    # improved truthiness parsing
    allow_bad_res = str(cfg.get("MEEKO_ALLOW_BAD_RES", "true")).lower() in (
        "1",
        "true",
        "yes",
    )
    default_altloc = (cfg.get("MEEKO_DEFAULT_ALTLOC") or "").strip()  # e.g. "A" or ""

    def _meeko_cmd() -> list[str]:
        return build_meeko_base_cmd()

    def _run(cmd: List[str]):
        logging.info("Meeko: %s", " ".join(cmd))
        r = run_subprocess_capture(cmd)
        (logging.info if r.returncode == 0 else logging.warning)(
            "rc=%s\nSTDOUT:\n%s\nSTDERR:\n%s",
            r.returncode,
            (r.stdout or ""),
            (r.stderr or ""),
        )
        return r

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
    with NamedTemporaryFile("w", suffix=".pdb", delete=False) as tmp1:
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
                logging.error(
                    "Helium persists after HIS/ion tweaks; skipping Meeko and trying ADT fallback."
                )
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

    # Meeko modern (--read_pdb) uses the HIS-classified, ion-sanitized tmp1_path
    if base_meeko:
        meeko_cmd = build_meeko_modern_cmd(base_meeko, tmp1_path, output_pdbqt)
        # Legacy single-shot form (we still capture rc/file size for logging)
        meeko_cmd_legacy = build_meeko_legacy_cmd(base_meeko, tmp1_path, output_pdbqt)
    else:
        meeko_cmd = None
        meeko_cmd_legacy = None
    adt_cmd = None
    # Only attempt ADT if BOTH are explicitly configured and present
    if not (
        MGLTOOLS_PYTHON
        and Path(MGLTOOLS_PYTHON).is_file()
        and os.access(MGLTOOLS_PYTHON, os.X_OK)
        and PREPARE_RECEPTOR_SCRIPT
        and Path(PREPARE_RECEPTOR_SCRIPT).is_file()
    ):
        logging.info(
            "[receptor] ADT fallback disabled (missing MGLTOOLS_PYTHON or PREPARE_RECEPTOR_SCRIPT)"
        )
    else:
        # --- Step 3: ADT prepare_receptor4 fallback (only if both keys are valid)
        mgltools_python = MGLTOOLS_PYTHON
        prepare_script = PREPARE_RECEPTOR_SCRIPT
        if (
            mgltools_python
            and Path(mgltools_python).is_file()
            and os.access(mgltools_python, os.X_OK)
            and prepare_script
            and Path(prepare_script).is_file()
        ):
            allow_tokens, _ = _load_retain_allowlist(cfg)
            logging.info(
                "[ions.policy] stage=adt variant=%s drop_metals=%s drop_salts=%s radius=none allowlist=%d",
                variant_label,
                False,
                False,
                len(
                    {
                        str(tok).strip().upper()
                        for tok in allow_tokens
                        if str(tok).strip()
                    }
                ),
            )
            adt_cmd = build_adt_prepare_receptor_cmd(
                mgltools_python,
                prepare_script,
                tmp1_path,
                output_pdbqt,
            )
            cp = run_subprocess_capture(adt_cmd)
            _persist_subproc(
                "adt_prepare_receptor4", adt_cmd, cp, work_dir, Path(output_pdbqt)
            )
            if cp.returncode == 0 and _ok_receptor_file(Path(output_pdbqt)):
                _clean_receptor_pdbqt(output_pdbqt)
                logging.info("Prepared receptor with ADT prepare_receptor4.py.")
                # --- ION DIFF BLOCK: summarize ions kept vs lost in PDBQT
                try:
                    _pdb_ions = _ion_pairs_from_pdb(input_pdb)
                    _pdbqt_ions = _ion_pairs_from_pdbqt(output_pdbqt)
                    _log_ion_diff("adt", _pdb_ions, _pdbqt_ions)
                except Exception as _e:
                    logging.warning("[ion diff] skipped note=%s", _e)
                try:
                    _log_pdb_pdbqt_counts_diff(input_pdb, output_pdbqt, tool="ADT")
                except Exception as diff_exc:
                    logging.warning(
                        "[iondiff.pdb_pdbqt] action=skip tool=ADT reason=%s",
                        diff_exc,
                    )

                return True
            if Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size == 0:
                head, tail = _first_last_lines(
                    work_dir / "adt_prepare_receptor4.stderr.txt"
                )
                logging.warning(
                    "[adt_prepare_receptor4] receptor PDBQT is empty. head=%r tail=%r",
                    head,
                    tail,
                )
            elif Path(output_pdbqt).exists():
                _clean_receptor_pdbqt(output_pdbqt)
        else:
            logging.info(
                "[receptor] ADT fallback disabled (MGLTOOLS_PYTHON/PREPARE_RECEPTOR_SCRIPT not set)"
            )

    # --- Step 1: Modern Meeko attempt
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
    if his_tie:
        mapping = _build_his_override_mapping(tmp1_path)
        if not mapping:
            m = re.search(r"residue_key='([A-Za-z]):(\d+)'", se0)
            if m:
                mapping = f"{m.group(1)}:{int(m.group(2))}={his_default}"

    # --- Step 1b: Retry with -n if needed
    if mapping and meeko_cmd:
        meeko_cmd_with_n = meeko_cmd + ["-n", mapping]
        cp = run_subprocess_capture(meeko_cmd_with_n)
        _persist_subproc(
            "meeko_retry_nmap", meeko_cmd_with_n, cp, work_dir, Path(output_pdbqt)
        )
        if cp.returncode == 0 and _ok_receptor_file(Path(output_pdbqt)):
            _clean_receptor_pdbqt(output_pdbqt)
            logging.info("Prepared receptor after HIS -n mapping.")
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

    # --- Step 1c: Retry with -a allow_bad_res (and default_altloc if hinted)
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
        if default_altloc:
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

    # --- Step 3: ADT prepare_receptor4 fallback
    if adt_cmd:
        cp = run_subprocess_capture(adt_cmd)
        _persist_subproc(
            "adt_prepare_receptor4", adt_cmd, cp, work_dir, Path(output_pdbqt)
        )
        if cp.returncode == 0 and _ok_receptor_file(Path(output_pdbqt)):
            _clean_receptor_pdbqt(output_pdbqt)
            logging.info("Prepared receptor with ADT prepare_receptor4.py.")
            # --- ION DIFF BLOCK: summarize ions kept vs lost in PDBQT
            try:
                _pdb_ions = _ion_pairs_from_pdb(input_pdb)
                _pdbqt_ions = _ion_pairs_from_pdbqt(output_pdbqt)
                _log_ion_diff("adt", _pdb_ions, _pdbqt_ions)
            except Exception as _e:
                logging.warning("[ion diff] skipped note=%s", _e)
            try:
                _log_pdb_pdbqt_counts_diff(input_pdb, output_pdbqt, tool="ADT")
            except Exception as diff_exc:
                logging.warning(
                    "[iondiff.pdb_pdbqt] action=skip tool=ADT reason=%s",
                    diff_exc,
                )

            return True

        if Path(output_pdbqt).exists() and Path(output_pdbqt).stat().st_size == 0:
            head, tail = _first_last_lines(
                work_dir / "adt_prepare_receptor4.stderr.txt"
            )
            logging.warning(
                "[adt_prepare_receptor4] receptor PDBQT is empty. head=%r tail=%r",
                head,
                tail,
            )
        elif Path(output_pdbqt).exists():
            _clean_receptor_pdbqt(output_pdbqt)

    if _ok_receptor_file(Path(output_pdbqt)):
        _clean_receptor_pdbqt(output_pdbqt)
        logging.warning(
            "[receptor] proceeding with existing receptor PDBQT despite prep errors."
        )
        return True
    return False
