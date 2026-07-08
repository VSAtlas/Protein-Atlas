"""Ion auditing and ion-diff diagnostics helpers."""

from __future__ import annotations

import json
import logging
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Optional, Sequence, Union

from protein_prep import pdb_fixer_runtime as _activesite_mod
from protein_prep.aliases_policy import (
    ALIASES,
    _ION_AUDIT_ALIAS_MAP,
    _ION_AUDIT_DISABLED_VALUES,
    _ION_AUDIT_ENABLED_VALUES,
    _ION_AUDIT_METALS,
    _ION_AUDIT_SIMPLE_IONS,
    _ION_AUDIT_WATERS,
    _ION_BREADCRUMB_METAL_ORDER,
    _ION_BREADCRUMB_SIMPLE_ORDER,
    _METAL_RESNAMES,
    _SALT_RESNAMES,
    _normalize_resname,
)
from protein_prep.prep_utils import _short_path_for_log

_ION_PIPE_AUDIT: dict = {}
_ION_PIPE_WARNED = False

_ION_AUDIT_STACK: list["_IonAuditManager"] = []


def _safe_short_path_for_log(path: Path) -> str:
    try:
        return str(_short_path_for_log(path))
    except Exception:
        return str(path)


def _get_activesite_module():
    return _activesite_mod

def _format_histogram(counter: Counter[str]) -> str:
    if not counter:
        return "none"
    items = sorted(counter.items(), key=lambda kv: (-kv[1], kv[0]))
    return ",".join(f"{name}:{count}" for name, count in items)

def _format_token_list(tokens: Iterable[str], *, limit: int = 8) -> str:
    unique = []
    seen: set[str] = set()
    for token in sorted(str(tok).strip().upper() for tok in tokens if str(tok).strip()):
        if token in seen:
            continue
        seen.add(token)
        unique.append(token)
    if not unique:
        return "none"
    head = unique[:limit]
    remaining = len(unique) - len(head)
    if remaining > 0:
        head.append(f"+{remaining}")
    return ",".join(head)

def _summarize_ions_file(file_path: Union[str, Path]) -> dict[str, object]:
    path = Path(file_path)
    if not path.exists():
        return {
            "res_hist": "missing",
            "elem_hist": "missing",
            "res_counts": {},
            "elem_counts": {},
            "metals_present": False,
            "salts_present": False,
            "error": "missing",
        }

    res_counts: Counter[str] = Counter()
    elem_counts: Counter[str] = Counter()
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if not line.startswith("HETATM"):
                    continue
                res = line[17:20].strip().upper()
                elem = (line[76:78].strip() or res).upper()
                if res:
                    res_counts[res] += 1
                if elem:
                    elem_counts[elem] += 1
    except Exception as exc:
        return {
            "res_hist": "error",
            "elem_hist": "error",
            "res_counts": {},
            "elem_counts": {},
            "metals_present": False,
            "salts_present": False,
            "error": str(exc),
        }

    metals_present = any(
        token in _ION_AUDIT_METALS and res_counts[token] > 0 for token in res_counts
    ) or any(
        token in _ION_AUDIT_METALS and elem_counts[token] > 0 for token in elem_counts
    )
    salts_present = any(
        token in _ION_AUDIT_SIMPLE_IONS and res_counts[token] > 0
        for token in res_counts
    ) or any(
        token in _ION_AUDIT_SIMPLE_IONS and elem_counts[token] > 0
        for token in elem_counts
    )

    return {
        "res_hist": _format_histogram(res_counts),
        "elem_hist": _format_histogram(elem_counts),
        "res_counts": dict(res_counts),
        "elem_counts": dict(elem_counts),
        "metals_present": metals_present,
        "salts_present": salts_present,
        "error": None,
    }

def _format_ion_pairs(pairs: Iterable[tuple[str, str]]) -> str:
    sorted_pairs = sorted(pairs)
    if not sorted_pairs:
        return "none"
    return ",".join(f"{res}:{loc}" for res, loc in sorted_pairs)

def _ion_audit_enabled() -> bool:
    raw = os.environ.get("ION_AUDIT", "")
    return raw.strip().lower() in _ION_AUDIT_ENABLED_VALUES

def _canon_ion_resname(resname: str) -> str:
    key = resname.strip().upper()
    return _ION_AUDIT_ALIAS_MAP.get(key, key)

def _format_breadcrumb_counts(counts: dict[str, int], order: Sequence[str]) -> str:
    parts: list[str] = []
    for key in order:
        parts.append(f"{key}:{int(counts.get(key, 0))}")
    extras = [key for key in sorted(counts) if key not in order]
    for key in extras:
        parts.append(f"{key}:{int(counts.get(key, 0))}")
    return "{" + ",".join(parts) + "}"

def _breadcrumbs_enabled() -> bool:
    raw = os.environ.get("ION_AUDIT")
    if raw is None:
        return True
    text = raw.strip()
    if not text:
        return True
    lowered = text.lower()
    if lowered in _ION_AUDIT_DISABLED_VALUES:
        return False
    return lowered in _ION_AUDIT_ENABLED_VALUES

def _emit_ion_breadcrumb(stage: str, file_path: Union[str, Path]) -> None:
    if not _breadcrumbs_enabled():
        return
    summarizer = getattr(_get_activesite_module(), "summarize_ions", None)
    if summarizer is None:
        return
    try:
        summary = summarizer(file_path)
    except Exception:
        summary = None
    short_path = _safe_short_path_for_log(Path(file_path))
    if not summary:
        logging.info(
            "[ions.breadcrumb] stage=%s file=%s missing=true", stage, short_path
        )
        return
    if summary.get("missing"):
        logging.info(
            "[ions.breadcrumb] stage=%s file=%s missing=true", stage, short_path
        )
        return
    metals = summary.get("metals", {}) or {}
    simple = summary.get("simple_ions", {}) or {}
    waters = int(summary.get("waters", 0) or 0)
    other = int(summary.get("other_het", 0) or 0)
    logging.info(
        "[ions.breadcrumb] stage=%s file=%s metals=%s waters=%d simple_ions=%s other_het=%d",
        stage,
        short_path,
        _format_breadcrumb_counts(metals, _ION_BREADCRUMB_METAL_ORDER),
        waters,
        _format_breadcrumb_counts(simple, _ION_BREADCRUMB_SIMPLE_ORDER),
        other,
    )

def _serialize_counts(counts: dict[str, int]) -> str:
    if not counts:
        return "{}"
    ordered = {key: counts[key] for key in sorted(counts)}
    return json.dumps(ordered, sort_keys=True)

def _legacy_ion_global_missing() -> None:
    global _ION_PIPE_WARNED
    if not _ION_PIPE_WARNED:
        logging.warning("[ion.audit.warn] disabled=legacy_global_missing")
        _ION_PIPE_WARNED = True

def _gather_ion_counts(path: Path) -> tuple[dict[str, int], dict[str, int], int, int]:
    metals: dict[str, set[tuple[str, str, str]]] = defaultdict(set)
    simple: dict[str, set[tuple[str, str, str]]] = defaultdict(set)
    waters: set[tuple[str, str, str]] = set()
    other: set[tuple[str, str, str, str]] = set()

    if not path.exists():
        return {}, {}, 0, 0

    try:
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.startswith("HETATM"):
                    continue
                resname_raw = line[17:20]
                resname = _canon_ion_resname(resname_raw)
                if not resname:
                    continue
                chain = (line[21:22] or "-").strip() or "-"
                resseq = (line[22:26] or "0").strip() or "0"
                icode = (line[26:27] or "").strip()
                resid = (chain, resseq, icode)
                if resname in _ION_AUDIT_METALS:
                    metals[resname].add(resid)
                elif resname in _ION_AUDIT_SIMPLE_IONS:
                    simple[resname].add(resid)
                elif resname in _ION_AUDIT_WATERS:
                    waters.add(resid)
                else:
                    other.add((resname, *resid))
    except Exception as exc:
        logging.warning("[ions.probe] stage=read_error file=%s err=%s", path, exc)
        return {}, {}, 0, 0

    metals_counts = {key: len(val) for key, val in metals.items() if val}
    simple_counts = {key: len(val) for key, val in simple.items() if val}
    waters_count = len(waters)
    other_count = len(other)
    return metals_counts, simple_counts, waters_count, other_count

def audit_ions(
    pdb_path: Union[str, Path],
    pdb_id: str,
    stage: str,
    variant: Optional[str],
    logger: logging.Logger | None = None,
) -> Optional[dict[str, object]]:
    if not _ion_audit_enabled():
        return None

    log = logger or logging.getLogger(__name__)
    path = Path(pdb_path)
    metals, simple, waters, other = _gather_ion_counts(path)
    stage_tag = "[ions.probe]"
    if stage == "input":
        stage_tag = "[ions.input.counts]"
    elif stage == "final_cleaned":
        stage_tag = "[ions.clean.counts]"

    log.info(
        "%s stage=%s file=%s present_pdb.metals=%s present_pdb.simple_ions=%s present_pdb.waters=%d present_pdb.other_het=%d",
        stage_tag,
        stage,
        _safe_short_path_for_log(path),
        _serialize_counts(metals),
        _serialize_counts(simple),
        waters,
        other,
    )

    try:
        resolved = str(path.resolve(strict=False))
    except Exception:
        resolved = str(path)

    record: dict[str, object] = {
        "stage": stage,
        "pdb": pdb_id,
        "variant": (variant or "NONE").upper(),
        "file": resolved,
        "counts": {
            "metals": {key: metals[key] for key in sorted(metals)},
            "simple_ions": {key: simple[key] for key in sorted(simple)},
            "waters": waters,
            "other_het": other,
        },
    }
    return record

def diff_ions(
    prev: dict[str, object],
    curr: dict[str, object],
    logger: logging.Logger | None = None,
) -> dict[str, object]:
    log = logger or logging.getLogger(__name__)
    prev_counts = prev.get("counts", {}) if isinstance(prev, dict) else {}
    curr_counts = curr.get("counts", {}) if isinstance(curr, dict) else {}

    def _as_dict(section: str) -> dict[str, int]:
        counts: dict[str, int] = {}
        for source in (prev_counts, curr_counts):
            if not isinstance(source, dict):
                continue
            mapping = source.get(section)
            if isinstance(mapping, dict):
                for key, value in mapping.items():
                    if isinstance(key, str) and isinstance(value, int):
                        counts.setdefault(key, 0)
        delta: dict[str, int] = {}
        for key in sorted(counts):
            prev_val = 0
            curr_val = 0
            prev_section = (
                prev_counts.get(section) if isinstance(prev_counts, dict) else {}
            )
            curr_section = (
                curr_counts.get(section) if isinstance(curr_counts, dict) else {}
            )
            if isinstance(prev_section, dict):
                prev_val = int(prev_section.get(key, 0))
            if isinstance(curr_section, dict):
                curr_val = int(curr_section.get(key, 0))
            delta_val = curr_val - prev_val
            if delta_val:
                delta[key] = delta_val
        return delta

    delta_metals = _as_dict("metals")
    delta_simple = _as_dict("simple_ions")
    prev_waters = (
        int(prev_counts.get("waters", 0)) if isinstance(prev_counts, dict) else 0
    )
    curr_waters = (
        int(curr_counts.get("waters", 0)) if isinstance(curr_counts, dict) else 0
    )
    prev_other = (
        int(prev_counts.get("other_het", 0)) if isinstance(prev_counts, dict) else 0
    )
    curr_other = (
        int(curr_counts.get("other_het", 0)) if isinstance(curr_counts, dict) else 0
    )

    delta_waters = curr_waters - prev_waters
    delta_other = curr_other - prev_other

    log.info(
        "[ions.diff] stage=%s delta.metals=%s delta.simple_ions=%s delta.waters=%+d delta.other_het=%+d",
        curr.get("stage", "unknown"),
        json.dumps(delta_metals, sort_keys=True),
        json.dumps(delta_simple, sort_keys=True),
        delta_waters,
        delta_other,
    )

    return {
        "metals": delta_metals,
        "simple_ions": delta_simple,
        "waters": delta_waters,
        "other_het": delta_other,
    }

class _IonAuditManager:
    def __init__(
        self, pdb_id: str, work_dir: Union[str, Path], variant: Optional[str]
    ) -> None:
        self.pdb_id = pdb_id
        self.work_dir = Path(work_dir)
        self.variant = (variant or "NONE").upper()
        self.enabled = _ion_audit_enabled()
        self._prev_record: dict[str, object] | None = None
        self._audit_dir = self.work_dir / "audits"
        if self.enabled:
            try:
                self._audit_dir.mkdir(parents=True, exist_ok=True)
            except Exception as exc:
                logging.warning(
                    "[ions.probe] stage=mkdir_failed dir=%s err=%s",
                    self._audit_dir,
                    exc,
                )
                self.enabled = False

    def probe(
        self,
        stage: str,
        file_path: Union[str, Path],
        *,
        variant_override: Optional[str] = None,
    ) -> None:
        if not self.enabled:
            return
        logger = logging.getLogger(__name__)
        record = audit_ions(
            file_path, self.pdb_id, stage, variant_override or self.variant, logger
        )
        if record is None:
            return
        if self._prev_record is not None:
            deltas = diff_ions(self._prev_record, record, logger)
            record["delta_vs_prev"] = deltas
        elif stage != "input":
            record["delta_vs_prev"] = {
                "metals": {},
                "simple_ions": {},
                "waters": 0,
                "other_het": 0,
            }

        if self.enabled:
            stage_token = stage.replace("/", "_")
            out_path = self._audit_dir / f"{stage_token}_ions.json"
            try:
                with out_path.open("w", encoding="utf-8") as fh:
                    json.dump(record, fh, indent=2, sort_keys=True)
                    fh.write("\n")
            except Exception as exc:
                logging.warning(
                    "[ions.probe] stage=json_write_failed file=%s err=%s", out_path, exc
                )

        self._prev_record = record

def _push_ion_audit_manager(manager: _IonAuditManager) -> None:
    if not manager.enabled:
        return
    _ION_AUDIT_STACK.append(manager)

def _pop_ion_audit_manager(manager: _IonAuditManager) -> None:
    if not manager.enabled:
        return
    if _ION_AUDIT_STACK and _ION_AUDIT_STACK[-1] is manager:
        _ION_AUDIT_STACK.pop()
        return
    try:
        _ION_AUDIT_STACK.remove(manager)
    except ValueError:
        pass

def _current_ion_audit_manager() -> Optional[_IonAuditManager]:
    if not _ION_AUDIT_STACK:
        return None
    return _ION_AUDIT_STACK[-1]

def emit_ion_audit_probe(
    stage: str, file_path: Union[str, Path], *, variant: Optional[str] = None
) -> None:
    manager = _current_ion_audit_manager()
    if manager is not None:
        manager.probe(stage, file_path, variant_override=variant)

def _reset_ion_probe(pdb_id: str) -> None:
    """Reset ion probe cache for a PDB identifier (case-normalized)."""
    if not pdb_id:
        return
    try:
        _ION_PIPE_AUDIT[pdb_id.upper()] = {}
    except NameError:
        _legacy_ion_global_missing()

def _ion_candidate_tokens(rules_obj=None) -> set[str]:
    if rules_obj is None:
        rules_obj = ALIASES
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

def _scan_metal_map(path: Union[str, Path], rules=None) -> dict[str, int]:
    if rules is None:
        rules = ALIASES
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

def _format_probe_counts(counts: dict[str, int], *, limit: int | None = None) -> str:
    if not counts:
        return "none"
    items = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    if limit is not None:
        items = items[:limit]
    return ",".join(f"{res}:{cnt}" for res, cnt in items)

def _log_ions_probe(
    pdb_id: str,
    stage: str,
    file_path: Union[str, Path],
    *,
    phase: str | None = None,
) -> dict[str, int]:
    counts = _scan_metal_map(file_path, ALIASES)
    total = sum(counts.values())
    sample = _format_probe_counts(counts, limit=5)
    if phase:
        logging.info(
            "[ions.probe] stage=%s phase=%s file=%s totals=%d sample=%s",
            stage,
            phase,
            file_path,
            total,
            sample,
        )
        key = f"{stage}:{phase}"
    else:
        logging.info(
            "[ions.probe] stage=%s file=%s totals=%d sample=%s",
            stage,
            file_path,
            total,
            sample,
        )
        key = stage
    try:
        bucket = _ION_PIPE_AUDIT.setdefault(pdb_id.upper() if pdb_id else "UNKNOWN", {})
    except NameError:
        _legacy_ion_global_missing()
    else:
        bucket[key] = counts
    return counts

def _format_diff_map(data: dict[str, int]) -> str:
    if not data:
        return "none"
    filtered = [(res, count) for res, count in data.items() if count > 0]
    if not filtered:
        return "none"
    items = sorted(filtered, key=lambda kv: (-kv[1], kv[0]))
    return ",".join(f"{res}:{cnt}" for res, cnt in items)

def get_ion_probe_map(pdb_id: str) -> dict[str, dict[str, int]]:
    try:
        bucket = _ION_PIPE_AUDIT.get((pdb_id or "").upper(), {})
    except NameError:
        _legacy_ion_global_missing()
        return {}
    return {k: dict(v) for k, v in bucket.items()}

def _bucket_counts(counter: Counter) -> dict[str, int]:
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
    bucketed = _bucket_counts(counter)
    keys = ["ZN", "MG", "NA", "K", "CA", "MN", "FE", "CL", "OTHER"]
    parts = [f"{k}={bucketed.get(k, 0)}" for k in keys]
    return " ".join(parts)

def _collect_monoatomic_records(pdb_path: Union[str, Path]) -> tuple[Counter, Counter]:
    path = Path(pdb_path)
    counts: Counter[str] = Counter()
    detail: Counter[tuple[str, str, str]] = Counter()
    if not path.exists():
        return counts, detail
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.startswith("HETATM"):
                    continue
                resname = line[17:20].strip().upper()
                elem = (line[76:78].strip() or resname).upper()
                canonical_res = _normalize_resname(resname)
                canonical_elem = _normalize_resname(elem)
                if not resname and not elem:
                    continue
                token = canonical_res or canonical_elem or resname or elem
                if (
                    canonical_res in _METAL_RESNAMES
                    or canonical_elem in _METAL_RESNAMES
                    or canonical_res in _SALT_RESNAMES
                    or canonical_elem in _SALT_RESNAMES
                ):
                    chain = (line[21] or "-").strip() or "-"
                    resseq = (line[22:26] or "0").strip() or "0"
                    counts[token] += 1
                    detail[(token, chain, resseq)] += 1
    except Exception as exc:  # logging-only helper
        logging.warning("[ions.probe] file=%s err=%s", path, exc)
    return counts, detail

def _format_ion_hist(counter: Counter) -> str:
    if not counter:
        return "none"
    parts = [f"{token}:{counter[token]}" for token in sorted(counter)]
    return ",".join(parts)

def _diff_detail_records(before: Counter, after: Counter) -> list[str]:
    missing = before - after
    out: list[str] = []
    for (token, chain, resseq), count in sorted(missing.items()):
        for _ in range(count):
            out.append(f"{token}:{chain}:{resseq}")
    return out
