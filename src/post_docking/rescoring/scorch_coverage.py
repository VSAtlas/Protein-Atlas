from __future__ import annotations

import csv
import json
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence, Tuple


ComboKey = Tuple[str, str, str]


def normalize_variant(value: Any) -> str:
    token = str(value or "").strip()
    if token.upper() in {"", "*", "BASE", "LEGACY", "NONE", "NULL"}:
        return ""
    return token.upper()


def normalize_ph(value: Any) -> str:
    token = str(value or "").strip()
    if token.lower() in {"", "*", "base", "none", "null"}:
        return ""
    return token


def normalize_pdb(value: Any) -> str:
    return str(value or "").strip().upper()


def normalize_combo(pdb_id: Any, variant: Any = "", ph: Any = "") -> ComboKey:
    return normalize_pdb(pdb_id), normalize_variant(variant), normalize_ph(ph)


@dataclass(frozen=True)
class ScorchStreamCoverage:
    stream: str
    consensus_csv: Optional[str]
    expected_rows: int
    raw_rows: int
    final_rescored_rows: int
    complete: bool
    invalid_rows: int = 0
    quarantined_rows: int = 0

    def to_dict(self) -> dict[str, Any]:
        return {
            "stream": self.stream,
            "consensus_csv": self.consensus_csv,
            "expected_rows": int(self.expected_rows),
            "raw_rows": int(self.raw_rows),
            "final_rescored_rows": int(self.final_rescored_rows),
            "invalid_rows": int(self.invalid_rows),
            "quarantined_rows": int(self.quarantined_rows),
            "complete": bool(self.complete),
        }


@dataclass(frozen=True)
class ScorchCoverageSummary:
    pdb_id: str
    variant: str
    ph: str
    top_fraction: float
    dock_dir: Optional[str]
    post_dir: Optional[str]
    fda: ScorchStreamCoverage
    dud: ScorchStreamCoverage
    all_complete: bool
    reasons: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "pdb_id": self.pdb_id,
            "variant": self.variant,
            "ph": self.ph,
            "top_fraction": float(self.top_fraction),
            "dock_dir": self.dock_dir,
            "post_dir": self.post_dir,
            "streams": {
                "fda": self.fda.to_dict(),
                "dud": self.dud.to_dict(),
            },
            "all_complete": bool(self.all_complete),
            "reasons": list(self.reasons),
        }


def _safe_fraction(value: Any, default: float = 0.10) -> float:
    try:
        fraction = float(value)
    except Exception:
        return default
    if not (0.0 < fraction <= 1.0):
        return default
    return fraction


def resolve_top_fraction(cfg: Optional[Mapping[str, Any]] = None, default: float = 0.10) -> float:
    cfg = cfg or {}
    for key in ("SCORCH_TOP_FRACTION", "ATLAS_SCORCH_TOP_FRACTION"):
        if key in cfg:
            return _safe_fraction(cfg.get(key), default=default)
    return default


def _candidate_combo_dirs(root: Path, combo: ComboKey) -> list[Path]:
    pdb_id, variant, ph = combo
    candidates: list[Path] = []
    pdb_root = root / pdb_id
    if variant and ph:
        candidates.append(pdb_root / variant / ph)
    if variant:
        candidates.append(pdb_root / variant)
    if ph:
        candidates.append(pdb_root / ph)
    candidates.append(pdb_root)
    if not variant:
        if ph:
            candidates.append(pdb_root / "LEGACY" / ph)
        candidates.append(pdb_root / "LEGACY" / "base")
        candidates.append(pdb_root / "LEGACY")

    deduped: list[Path] = []
    seen: set[str] = set()
    for path in candidates:
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(path)
    return deduped


def find_combo_dir(root: Path, combo: ComboKey, filenames: Sequence[str]) -> Optional[Path]:
    for candidate in _candidate_combo_dirs(root, combo):
        for filename in filenames:
            path = candidate / filename
            if path.is_file() and path.stat().st_size > 0:
                return candidate
    return None


def combo_done_sentinel(post_run_root: Path, combo: ComboKey) -> Path:
    pdb_id, variant, ph = combo
    scope = Path(post_run_root) / pdb_id
    if variant:
        scope = scope / variant
    if ph:
        scope = scope / ph
    return scope / "scorch" / "_DONE"


def combo_coverage_summary_path(post_run_root: Path, combo: ComboKey) -> Path:
    return combo_done_sentinel(post_run_root, combo).parent / "scorch_coverage_summary.json"


def _count_csv_rows(path: Optional[Path]) -> int:
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return 0
    with path.open(newline="", encoding="utf-8") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def _read_csv_rows(path: Optional[Path]) -> Iterable[dict[str, str]]:
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _count_matching(path: Optional[Path], predicate) -> int:
    return sum(1 for row in _read_csv_rows(path) if predicate(row))


def _truthy_flag(value: Any) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y"}


def accept_quarantined_coverage(cfg: Optional[Mapping[str, Any]] = None) -> bool:
    cfg = cfg or {}
    for key in (
        "ATLAS_SCORCH_ACCEPT_QUARANTINED_COVERAGE",
        "SCORCH_ACCEPT_QUARANTINED_COVERAGE",
        "ATLAS_SCORCH_ALLOW_QUARANTINED_COVERAGE",
        "SCORCH_ALLOW_QUARANTINED_COVERAGE",
    ):
        raw = os.environ.get(key)
        if raw is None:
            raw = cfg.get(key)
        if raw is not None:
            return _truthy_flag(raw)
    return True


_INVALID_SCORCH_STATUSES = {
    "error",
    "failed",
    "failure",
    "missing",
    "quarantine",
    "quarantined",
    "skipped",
    "timeout",
    "timed_out",
}


def _valid_scorch_row(row: Mapping[str, Any]) -> bool:
    """Return False for placeholder rows that represent failed SCORCH work."""
    status = str(row.get("scorch_status") or row.get("status") or "").strip().lower()
    if status in _INVALID_SCORCH_STATUSES:
        return False
    if _truthy_flag(row.get("scorch_unscorable_flag")):
        return False
    if str(row.get("scorch_failure_reason") or "").strip():
        return False
    returncode = str(row.get("scorch_returncode") or "").strip()
    if returncode and returncode not in {"0", "0.0"}:
        return False
    return True


_ACCEPTED_QUARANTINE_REASONS = {
    "explicit_ligand_quarantine",
    "global_ligand_timeout",
}


def _accepted_quarantined_scorch_row(row: Mapping[str, Any]) -> bool:
    """Return True for intentionally quarantined rows accepted by policy."""
    status = str(row.get("scorch_status") or row.get("status") or "").strip().lower()
    if status not in {"quarantine", "quarantined"}:
        return False
    if not _truthy_flag(row.get("scorch_unscorable_flag")):
        return False
    reason = str(row.get("scorch_failure_reason") or "").strip().lower()
    return reason in _ACCEPTED_QUARANTINE_REASONS


_RUN_MODE_STAGE_SUFFIX_RE = re.compile(
    r"_(?:dud|decoy|decoys|fda|prod|production|fda_dud|dud_fda|test_library_\d+|dry_bench_\d+)_(?:stage|pose)\d+$",
    re.IGNORECASE,
)
_STAGE_SUFFIX_RE = re.compile(r"_(?:stage|pose)\d+$", re.IGNORECASE)


def _ligand_key(value: Any) -> str:
    token = Path(str(value or "").strip()).name
    if not token:
        return ""
    for suffix in (".pdbqt", ".pdb", ".sdf", ".mol2", ".csv"):
        if token.lower().endswith(suffix):
            token = token[: -len(suffix)]
            break
    token = token.replace(".sanitized", "")
    token = _RUN_MODE_STAGE_SUFFIX_RE.sub("", token)
    return _STAGE_SUFFIX_RE.sub("", token).strip().lower()


def _row_ligand_key(row: Mapping[str, Any]) -> str:
    for field in ("ligand", "ligand_file", "Ligand_ID", "ligand_id"):
        key = _ligand_key(row.get(field))
        if key:
            return key
    return ""


def _row_mode(row: Mapping[str, Any], default: str = "") -> str:
    return str(row.get("run_mode") or default or "").strip().lower()


def _is_decoy_row(row: Mapping[str, Any]) -> bool:
    return _truthy_flag(row.get("is_decoy")) or _row_mode(row) in {"dud", "decoy", "decoys"}


def _stream_hint(row: Mapping[str, Any], default: str = "") -> str:
    fields = (
        str(row.get("run_mode") or ""),
        str(row.get("library") or ""),
        str(row.get("stream") or ""),
    )
    tokens = {field.strip().lower() for field in fields if field.strip()}
    if _truthy_flag(row.get("is_decoy")) or tokens.intersection({"dud", "dud_fda", "decoy", "decoys"}):
        return "dud"
    if tokens.intersection({"fda", "production", "prod"}):
        return "fda"
    return default if default in {"fda", "dud"} else ""


def _expected_rows(consensus_csv: Optional[Path], top_fraction: float) -> int:
    rows = _count_csv_rows(consensus_csv)
    return int(math.ceil(float(rows) * float(top_fraction))) if rows > 0 else 0


def _expected_ligands(
    consensus_csv: Optional[Path],
    top_fraction: float,
    *,
    default_stream: str,
    exclude_bases: Optional[set[str]] = None,
) -> dict[str, set[str]]:
    rows = list(_read_csv_rows(consensus_csv))
    if not rows:
        return {"fda": set(), "dud": set()}
    limit = int(math.ceil(float(len(rows)) * float(top_fraction)))
    expected = {"fda": set(), "dud": set()}
    excluded = exclude_bases or set()
    for row in rows[:limit]:
        key = _row_ligand_key(row)
        if not key:
            continue
        if key in excluded:
            continue
        stream = _stream_hint(row, default=default_stream)
        if stream not in expected:
            stream = default_stream
        expected[stream].add(key)
    return expected


def _first_existing(parent: Optional[Path], filenames: Sequence[str]) -> Optional[Path]:
    if parent is None:
        return None
    for filename in filenames:
        candidate = parent / filename
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    return None


def _first_matching(parent: Optional[Path], patterns: Sequence[str], exclude: set[str]) -> Optional[Path]:
    if parent is None:
        return None
    for pattern in patterns:
        for candidate in sorted(parent.glob(pattern)):
            if candidate.name in exclude:
                continue
            if candidate.is_file() and candidate.stat().st_size > 0:
                return candidate
    return None


def _max_glob_rows(parent: Optional[Path], patterns: Sequence[str], exclude: set[str]) -> int:
    if parent is None:
        return 0
    best = 0
    for pattern in patterns:
        for candidate in parent.glob(pattern):
            if candidate.name in exclude:
                continue
            best = max(best, _count_csv_rows(candidate))
    return best


def _max_glob_rescored_rows(
    parent: Optional[Path],
    patterns: Sequence[str],
    exclude: set[str],
) -> int:
    if parent is None:
        return 0
    best = 0
    for pattern in patterns:
        for candidate in parent.glob(pattern):
            if candidate.name in exclude:
                continue
            best = max(
                best,
                _count_matching(
                    candidate,
                    lambda row: _truthy_flag(row.get("rescored_flag")),
                ),
            )
    return best


def _max_rows(parent: Optional[Path], filenames: Sequence[str]) -> int:
    if parent is None:
        return 0
    return max((_count_csv_rows(parent / filename) for filename in filenames), default=0)


def _max_rescored_rows(parent: Optional[Path], filenames: Sequence[str]) -> int:
    if parent is None:
        return 0
    return max(
        (
            _count_matching(
                parent / filename,
                lambda row: _truthy_flag(row.get("rescored_flag")),
            )
            for filename in filenames
        ),
        default=0,
    )


def _csv_candidates(parent: Optional[Path], filenames: Sequence[str]) -> list[Path]:
    if parent is None:
        return []
    paths: list[Path] = []
    seen: set[str] = set()
    for filename in filenames:
        path = parent / filename
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        if path.is_file() and path.stat().st_size > 0:
            paths.append(path)
    return paths


def _glob_csv_candidates(parent: Optional[Path], patterns: Sequence[str], exclude: set[str]) -> list[Path]:
    if parent is None:
        return []
    paths: list[Path] = []
    seen: set[str] = set()
    for pattern in patterns:
        for path in sorted(parent.glob(pattern)):
            if path.name in exclude:
                continue
            key = str(path)
            if key in seen:
                continue
            seen.add(key)
            if path.is_file() and path.stat().st_size > 0:
                paths.append(path)
    return paths


def _count_expected_hits(
    paths: Sequence[Path],
    *,
    expected_fda: set[str],
    expected_dud: set[str],
    default_stream: str = "",
    require_rescored: bool = False,
    accept_quarantined: bool = False,
) -> dict[str, int]:
    hits = {"fda": set(), "dud": set()}
    expected_present = bool(expected_fda or expected_dud)
    for path in paths:
        for row in _read_csv_rows(path):
            accepted_quarantine = bool(
                accept_quarantined and _accepted_quarantined_scorch_row(row)
            )
            if not (_valid_scorch_row(row) or accepted_quarantine):
                continue
            if require_rescored and not _truthy_flag(row.get("rescored_flag")):
                continue
            key = _row_ligand_key(row)
            if not key:
                continue
            if key in expected_fda:
                hits["fda"].add(key)
                continue
            if key in expected_dud:
                hits["dud"].add(key)
                continue
            if expected_present:
                continue
            stream = _stream_hint(row, default=default_stream)
            if stream in hits:
                hits[stream].add(key)
    return {stream: len(values) for stream, values in hits.items()}


def _count_invalid_expected_hits(
    paths: Sequence[Path],
    *,
    expected_fda: set[str],
    expected_dud: set[str],
    default_stream: str = "",
    accept_quarantined: bool = False,
) -> dict[str, int]:
    hits = {"fda": set(), "dud": set()}
    expected_present = bool(expected_fda or expected_dud)
    for path in paths:
        for row in _read_csv_rows(path):
            if _valid_scorch_row(row):
                continue
            if accept_quarantined and _accepted_quarantined_scorch_row(row):
                continue
            key = _row_ligand_key(row)
            if not key:
                continue
            if key in expected_fda:
                hits["fda"].add(key)
                continue
            if key in expected_dud:
                hits["dud"].add(key)
                continue
            if expected_present:
                continue
            stream = _stream_hint(row, default=default_stream)
            if stream in hits:
                hits[stream].add(key)
    return {stream: len(values) for stream, values in hits.items()}


def _count_quarantined_expected_hits(
    paths: Sequence[Path],
    *,
    expected_fda: set[str],
    expected_dud: set[str],
    default_stream: str = "",
) -> dict[str, int]:
    hits = {"fda": set(), "dud": set()}
    expected_present = bool(expected_fda or expected_dud)
    for path in paths:
        for row in _read_csv_rows(path):
            if not _accepted_quarantined_scorch_row(row):
                continue
            key = _row_ligand_key(row)
            if not key:
                continue
            if key in expected_fda:
                hits["fda"].add(key)
                continue
            if key in expected_dud:
                hits["dud"].add(key)
                continue
            if expected_present:
                continue
            stream = _stream_hint(row, default=default_stream)
            if stream in hits:
                hits[stream].add(key)
    return {stream: len(values) for stream, values in hits.items()}


def _derive_manifest_run_root(docked_run_root: Path) -> Optional[Path]:
    run_root = Path(docked_run_root)
    run_id = run_root.name
    if not run_id:
        return None
    docked_parent = run_root.parent
    if docked_parent.name != "docked":
        return None
    outputs = docked_parent.parent
    return outputs / "manifests" / run_id


def _load_shard_expected_ligands(
    manifest_run_root: Optional[Path],
    combo: ComboKey,
    *,
    decoy_prefix: str,
    control_bases: set[str],
) -> tuple[dict[str, set[str]], set[str]]:
    expected = {"fda": set(), "dud": set()}
    streams_seen: set[str] = set()
    if manifest_run_root is None:
        return expected, streams_seen
    plans_dir = Path(manifest_run_root) / "distributed" / "scorch_shard_plans"
    if not plans_dir.is_dir():
        return expected, streams_seen
    prefix = str(decoy_prefix or "dud").strip().lower() or "dud"
    matching_payloads: list[dict[str, Any]] = []
    fallback_payloads: list[dict[str, Any]] = []
    for plan_path in sorted(plans_dir.glob("*.json")):
        try:
            payload = json.loads(plan_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        if not isinstance(payload, dict):
            continue
        plan_prefix = str(payload.get("decoy_prefix") or prefix).strip().lower() or prefix
        raw_combo = payload.get("combo")
        if not isinstance(raw_combo, list):
            continue
        plan_combo = normalize_combo(
            raw_combo[0] if len(raw_combo) > 0 else "",
            raw_combo[1] if len(raw_combo) > 1 else "",
            raw_combo[2] if len(raw_combo) > 2 else "",
        )
        if plan_combo != combo:
            continue
        if plan_prefix == prefix:
            matching_payloads.append(payload)
        else:
            fallback_payloads.append(payload)

    # Older mixed FDA/DUD-FDA runs can record SCORCH shard plans with the
    # physical file prefix ("dud") while the active config reports the logical
    # test-library prefix ("dud_fda"). Prefer exact-prefix plans when present,
    # but do not drop otherwise matching combo plans and fall back to consensus
    # selection; that fallback can misclassify control ligands as FDA top hits.
    selected_payloads = matching_payloads or fallback_payloads
    for payload in selected_payloads:
        shards = payload.get("shards")
        if not isinstance(shards, list):
            continue
        for raw_record in shards:
            if not isinstance(raw_record, dict):
                continue
            run_mode = str(raw_record.get("run_mode") or "").strip().lower()
            if run_mode in {"dud", "decoy", "decoys", prefix}:
                stream = "dud"
            elif run_mode in {"fda", "production", "prod"}:
                stream = "fda"
            else:
                continue
            streams_seen.add(stream)
            record_controls = {
                _ligand_key(base)
                for base in (raw_record.get("control_bases") or [])
                if _ligand_key(base)
            }
            excluded = set(control_bases) | record_controls
            for base in raw_record.get("allowed_bases") or []:
                key = _ligand_key(base)
                if key and key not in excluded:
                    expected[stream].add(key)
    return expected, streams_seen


def evaluate_scorch_coverage(
    *,
    docked_run_root: Path,
    post_run_root: Path,
    pdb_id: Any,
    variant: Any = "",
    ph: Any = "",
    top_fraction: float = 0.10,
    decoy_prefix: str = "dud",
    control_bases: Optional[Iterable[str]] = None,
    manifest_run_root: Optional[Path] = None,
    accept_quarantined: Optional[bool] = None,
) -> ScorchCoverageSummary:
    combo = normalize_combo(pdb_id, variant, ph)
    fraction = _safe_fraction(top_fraction)
    prefix = str(decoy_prefix or "dud").strip() or "dud"
    accept_quarantine = (
        accept_quarantined_coverage() if accept_quarantined is None else bool(accept_quarantined)
    )
    excluded_controls = {
        _ligand_key(base) for base in (control_bases or ()) if _ligand_key(base)
    }
    dock_dir = find_combo_dir(
        Path(docked_run_root),
        combo,
        (
            "consensus_docking_scores.csv",
            "dud_consensus_docking_scores.csv",
            f"{prefix}_consensus_docking_scores.csv",
            "decoys_consensus_docking_scores.csv",
        ),
    )
    post_dir = find_combo_dir(
        Path(post_run_root),
        combo,
        (
            "scorch_scores_all.csv",
            f"{prefix}_scorch_scores_all.csv",
            "dud_scorch_scores_all.csv",
            "decoys_scorch_scores_all.csv",
            "consensus_reranked_scorch.csv",
            f"{prefix}_consensus_reranked_scorch.csv",
            "dud_consensus_reranked_scorch.csv",
            "decoys_consensus_reranked_scorch.csv",
        ),
    )

    fda_consensus = _first_existing(dock_dir, ("consensus_docking_scores.csv",))
    dud_consensus = _first_existing(
        dock_dir,
        (
            f"{prefix}_consensus_docking_scores.csv",
            "dud_consensus_docking_scores.csv",
            "decoys_consensus_docking_scores.csv",
        ),
    )
    if dud_consensus is None:
        dud_consensus = _first_matching(
            dock_dir,
            ("*_consensus_docking_scores.csv",),
            {"consensus_docking_scores.csv"},
        )
    fda_expected_ligands = _expected_ligands(
        fda_consensus,
        fraction,
        default_stream="fda",
        exclude_bases=excluded_controls,
    )
    dud_expected_ligands = _expected_ligands(
        dud_consensus,
        fraction,
        default_stream="dud",
        exclude_bases=excluded_controls,
    )
    fallback_fda = set(fda_expected_ligands["fda"]) | set(dud_expected_ligands["fda"])
    fallback_dud = set(fda_expected_ligands["dud"]) | set(dud_expected_ligands["dud"])
    shard_expected, shard_streams = _load_shard_expected_ligands(
        manifest_run_root or _derive_manifest_run_root(Path(docked_run_root)),
        combo,
        decoy_prefix=prefix,
        control_bases=excluded_controls,
    )
    expected_fda = (
        set(shard_expected["fda"]) if "fda" in shard_streams else fallback_fda
    )
    expected_dud = (
        set(shard_expected["dud"]) if "dud" in shard_streams else fallback_dud
    )
    fda_expected = len(expected_fda)
    dud_expected = len(expected_dud)
    has_consensus_rows = (
        _count_csv_rows(fda_consensus) + _count_csv_rows(dud_consensus)
    ) > 0

    raw_paths = _csv_candidates(
        post_dir,
        (
            "scorch_scores_all.csv",
            f"{prefix}_scorch_scores_all.csv",
            "dud_scorch_scores_all.csv",
            "decoys_scorch_scores_all.csv",
        ),
    )
    raw_paths.extend(
        _glob_csv_candidates(
            post_dir,
            ("*_scorch_scores_all.csv",),
            {
                "scorch_scores_all.csv",
                f"{prefix}_scorch_scores_all.csv",
                "dud_scorch_scores_all.csv",
                "decoys_scorch_scores_all.csv",
            },
        )
    )
    raw_hits = _count_expected_hits(
        raw_paths,
        expected_fda=expected_fda,
        expected_dud=expected_dud,
        accept_quarantined=accept_quarantine,
    )
    raw_fda = raw_hits["fda"]
    raw_dud = raw_hits["dud"]
    raw_invalid_hits = _count_invalid_expected_hits(
        raw_paths,
        expected_fda=expected_fda,
        expected_dud=expected_dud,
        accept_quarantined=accept_quarantine,
    )
    raw_quarantine_hits = _count_quarantined_expected_hits(
        raw_paths,
        expected_fda=expected_fda,
        expected_dud=expected_dud,
    )

    final_paths = _csv_candidates(
        post_dir,
        (
            "consensus_reranked_scorch.csv",
            f"{prefix}_consensus_reranked_scorch.csv",
            "dud_consensus_reranked_scorch.csv",
            "decoys_consensus_reranked_scorch.csv",
        ),
    )
    final_paths.extend(
        _glob_csv_candidates(
            post_dir,
            ("*_consensus_reranked_scorch.csv",),
            {
                "consensus_reranked_scorch.csv",
                f"{prefix}_consensus_reranked_scorch.csv",
                "dud_consensus_reranked_scorch.csv",
                "decoys_consensus_reranked_scorch.csv",
            },
        )
    )
    final_hits = _count_expected_hits(
        final_paths,
        expected_fda=expected_fda,
        expected_dud=expected_dud,
        require_rescored=True,
        accept_quarantined=accept_quarantine,
    )
    final_fda = final_hits["fda"]
    final_dud = final_hits["dud"]
    final_invalid_hits = _count_invalid_expected_hits(
        final_paths,
        expected_fda=expected_fda,
        expected_dud=expected_dud,
        accept_quarantined=accept_quarantine,
    )
    final_quarantine_hits = _count_quarantined_expected_hits(
        final_paths,
        expected_fda=expected_fda,
        expected_dud=expected_dud,
    )
    invalid_fda = raw_invalid_hits["fda"] + final_invalid_hits["fda"]
    invalid_dud = raw_invalid_hits["dud"] + final_invalid_hits["dud"]
    quarantine_fda = raw_quarantine_hits["fda"] + final_quarantine_hits["fda"]
    quarantine_dud = raw_quarantine_hits["dud"] + final_quarantine_hits["dud"]

    fda_complete = fda_expected == 0 or (
        raw_fda >= fda_expected and final_fda >= fda_expected
    )
    dud_complete = dud_expected == 0 or (
        raw_dud >= dud_expected and final_dud >= dud_expected
    )
    if invalid_fda > 0:
        fda_complete = False
    if invalid_dud > 0:
        dud_complete = False
    reasons: list[str] = []
    if not has_consensus_rows:
        reasons.append("missing_docking_consensus")
        fda_complete = False
        dud_complete = False
    if fda_expected > 0 and raw_fda < fda_expected:
        reasons.append(f"fda_raw_under_coverage:{raw_fda}/{fda_expected}")
    if fda_expected > 0 and final_fda < fda_expected:
        reasons.append(f"fda_final_under_coverage:{final_fda}/{fda_expected}")
    if invalid_fda > 0:
        reasons.append(f"fda_invalid_scorch_rows:{invalid_fda}")
    if dud_expected > 0 and raw_dud < dud_expected:
        reasons.append(f"dud_raw_under_coverage:{raw_dud}/{dud_expected}")
    if dud_expected > 0 and final_dud < dud_expected:
        reasons.append(f"dud_final_under_coverage:{final_dud}/{dud_expected}")
    if invalid_dud > 0:
        reasons.append(f"dud_invalid_scorch_rows:{invalid_dud}")

    fda = ScorchStreamCoverage(
        stream="fda",
        consensus_csv=str(fda_consensus) if fda_consensus else None,
        expected_rows=fda_expected,
        raw_rows=raw_fda,
        final_rescored_rows=final_fda,
        complete=fda_complete,
        invalid_rows=invalid_fda,
        quarantined_rows=quarantine_fda,
    )
    dud = ScorchStreamCoverage(
        stream="dud",
        consensus_csv=str(dud_consensus) if dud_consensus else None,
        expected_rows=dud_expected,
        raw_rows=raw_dud,
        final_rescored_rows=final_dud,
        complete=dud_complete,
        invalid_rows=invalid_dud,
        quarantined_rows=quarantine_dud,
    )
    return ScorchCoverageSummary(
        pdb_id=combo[0],
        variant=combo[1],
        ph=combo[2],
        top_fraction=fraction,
        dock_dir=str(dock_dir) if dock_dir else None,
        post_dir=str(post_dir) if post_dir else None,
        fda=fda,
        dud=dud,
        all_complete=bool(fda_complete and dud_complete),
        reasons=tuple(reasons),
    )


def write_coverage_summary(post_run_root: Path, summary: ScorchCoverageSummary) -> Path:
    path = combo_coverage_summary_path(
        Path(post_run_root),
        (summary.pdb_id, summary.variant, summary.ph),
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(summary.to_dict(), sort_keys=True, indent=2) + "\n",
        encoding="utf-8",
    )
    return path


def remove_done_sentinel(post_run_root: Path, summary: ScorchCoverageSummary) -> bool:
    path = combo_done_sentinel(
        Path(post_run_root),
        (summary.pdb_id, summary.variant, summary.ph),
    )
    if not path.exists():
        return False
    path.unlink()
    return True
