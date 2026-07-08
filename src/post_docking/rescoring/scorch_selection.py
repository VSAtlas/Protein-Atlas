from __future__ import annotations

import csv
import logging
import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

from post_docking.rescoring.scorch_prefix_and_stages import stage_dir_candidates
from post_docking.rescoring.scorch_types import (
    PoseSelectionResult,
    SelectionResult,
    StageSpec,
)

CONSENSUS_INPUT_ORDER_SENTINEL = "__atlas_consensus_input_order__"


def _decoy_prefix_aliases(decoy_prefix: str) -> Tuple[str, ...]:
    prefix = str(decoy_prefix or "dud").strip() or "dud"
    aliases = [prefix]
    for alias in ("fda_dud", "dud"):
        if alias not in aliases:
            aliases.append(alias)
    return tuple(sorted(aliases, key=len, reverse=True))


def _row_stream_hint(row: Dict[str, object], *, decoy_prefix: str) -> str:
    tokens = {
        str(row.get(field) or "").strip().lower()
        for field in ("run_mode", "library", "stream")
        if str(row.get(field) or "").strip()
    }
    decoy_aliases = set(_decoy_prefix_aliases(decoy_prefix))
    if str(row.get("is_decoy") or "").strip().lower() in {"1", "true", "yes", "y"}:
        return "dud"
    if tokens.intersection({"dud", "decoy", "decoys", *decoy_aliases}):
        return "dud"
    if tokens.intersection({"fda", "production", "prod"}):
        return "fda"
    return ""


def _generic_consensus_stream(combo_root: Path, *, decoy_prefix: str) -> str:
    consensus = combo_root / "consensus_docking_scores.csv"
    if not consensus.is_file() or consensus.stat().st_size <= 0:
        return ""
    seen: set[str] = set()
    try:
        with consensus.open(newline="", encoding="utf-8") as handle:
            for idx, row in enumerate(csv.DictReader(handle)):
                stream = _row_stream_hint(row, decoy_prefix=decoy_prefix)
                if stream:
                    seen.add(stream)
                if idx >= 49 or len(seen) > 1:
                    break
    except OSError:
        return ""
    return next(iter(seen)) if len(seen) == 1 else ""


def collect_stage_pdbqts(ph_root: Path, stage_dir: str) -> List[Path]:
    candidates = [
        ph_root / stage_dir,
        ph_root / "retries" / stage_dir,
        ph_root / stage_dir / "retries",
        ph_root / f"{stage_dir}_retries",
        ph_root / f"{stage_dir}_retry",
    ]

    for child in ph_root.iterdir() if ph_root.exists() else []:
        if (
            child.is_dir()
            and child.name.startswith(stage_dir)
            and "retry" in child.name
        ):
            candidates.append(child)

    collected: Set[Path] = set()
    for candidate in candidates:
        if not candidate.exists():
            continue
        if candidate.is_file() and candidate.suffix == ".pdbqt":
            collected.add(candidate.resolve())
            continue
        if candidate.is_dir():
            for pdbqt in candidate.rglob("*.pdbqt"):
                collected.add(pdbqt.resolve())

    return sorted(collected)


def pose_base_from_path(p: Path, *, decoy_prefix: str) -> str:
    stem = p.stem
    stem = stem.replace(".sanitized", "")
    for alias in _decoy_prefix_aliases(decoy_prefix):
        prefix = re.escape(alias)
        stem = re.sub(rf"(_gnina_{prefix}_stage\d+)$", "", stem)
        stem = re.sub(rf"(_{prefix}_gnina_stage\d+)$", "", stem)
        stem = re.sub(rf"(_dock6_{prefix}_stage\d+)$", "", stem)
        stem = re.sub(rf"(_{prefix}_dock6_stage\d+)$", "", stem)
        stem = re.sub(rf"(_{prefix}_stage\d+)$", "", stem)
        stem = re.sub(rf"(__{prefix}_ledock_stage\d+)$", "", stem)
        stem = re.sub(rf"(__ledock_{prefix}_stage\d+)$", "", stem)
        stem = re.sub(rf"(__{prefix}_dock6_stage\d+)$", "", stem)
        stem = re.sub(rf"(__dock6_{prefix}_stage\d+)$", "", stem)
    stem = re.sub(r"(__ledock_stage\d+)$", "", stem)
    stem = re.sub(r"(__dock6_stage\d+)$", "", stem)
    stem = re.sub(r"(_gnina_stage\d+)$", "", stem)
    stem = re.sub(r"(_stage\d+)$", "", stem)
    stem = re.sub(r"\.(mol2|pdbqt)$", "", stem, flags=re.IGNORECASE)
    stem = stem.replace("__", "_")
    stem = re.sub(r"_+$", "", stem)
    return stem


def control_base_from_path(p: Path) -> str:
    stem = p.stem.split("_stage")[0]
    return stem.split(".sanitized")[0]


def load_control_bases(
    processed_root: Path, pdb_id: str, logger: logging.Logger, *, component: str
) -> Set[str]:
    roots = [
        processed_root / pdb_id / "ligands_raw",
        processed_root / f"{pdb_id}_NOLIG" / "ligands_raw",
    ]
    bases: Set[str] = set()
    for root in roots:
        if not root.exists():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            base = control_base_from_path(path)
            if base:
                bases.add(base)
    logger.debug(
        "%s action=controls status=ok pdb_id=%s bases=%d roots=%s",
        component,
        pdb_id,
        len(bases),
        ",".join(str(r) for r in roots),
    )
    return bases


def stage_priority(stage_dir: str, pdbqt: Optional[Path] = None) -> int:
    name_lower = stage_dir.lower()
    if name_lower.endswith("stage3") or "stage3" in name_lower:
        return 3
    if name_lower.endswith("stage2") or "stage2" in name_lower:
        return 2
    if name_lower.endswith("stage1") or "stage1" in name_lower:
        return 1
    if pdbqt is not None:
        stem = pdbqt.stem
        match = re.search(r"__dock6_stage(\d+)", stem, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
        match = re.search(r"__ledock_stage(\d+)", stem, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
    return 0


def _priority_to_stage_label(priority: int, fallback: str = "") -> str:
    if priority > 0:
        return f"stage{priority}"
    return fallback


def _resolved_stage_label(stage_dir: str, pdbqt: Path, fallback: str = "") -> str:
    stem = pdbqt.stem.lower()
    for pattern in (
        r"(gnina_stage\d+)",
        r"(dock6_stage\d+)",
        r"(ledock_stage\d+)",
        r"(stage\d+)",
    ):
        match = re.search(pattern, stem)
        if match:
            return match.group(1)
    lower = stage_dir.lower()
    if "gnina" in lower:
        priority = stage_priority(stage_dir, pdbqt)
        return f"gnina_stage{priority}" if priority > 0 else fallback
    if "dock6" in lower:
        priority = stage_priority(stage_dir, pdbqt)
        return f"dock6_stage{priority}" if priority > 0 else fallback
    if "ledock" in lower:
        priority = stage_priority(stage_dir, pdbqt)
        return f"ledock_stage{priority}" if priority > 0 else fallback
    return _priority_to_stage_label(stage_priority(stage_dir, pdbqt), fallback)


def _normalize_stage_label(raw: object, ligand_val: str = "") -> str:
    text = str(raw or "").strip()
    if text:
        lower = text.lower()
        if lower.isdigit():
            return f"stage{lower}"
        match = re.search(
            r"(gnina_stage\d+|dock6_stage\d+|ledock_stage\d+|stage\d+)",
            lower,
        )
        if match:
            return match.group(1)
    if ligand_val:
        lower_lig = str(ligand_val).strip().lower()
        match = re.search(
            r"(gnina_stage\d+|dock6_stage\d+|ledock_stage\d+|stage\d+)",
            lower_lig,
        )
        if match:
            return match.group(1)
    return ""


def collect_best_pose_per_base(
    ph_root: Path,
    stage_dirs: Sequence[str],
    allowed_bases: Optional[Set[str]],
    logger: logging.Logger,
    *,
    component: str,
    decoy_prefix: str,
    preferred_stage_by_base: Optional[Dict[str, str]] = None,
) -> PoseSelectionResult:
    preferred_stage_by_base = preferred_stage_by_base or {}
    candidates_by_base: Dict[str, List[Tuple[int, Path, str]]] = {}
    total_candidates = 0
    for stage_dir in stage_dirs:
        stage_pdbqts = collect_stage_pdbqts(ph_root, stage_dir)
        total_candidates += len(stage_pdbqts)
        for pdbqt in stage_pdbqts:
            base = pose_base_from_path(pdbqt, decoy_prefix=decoy_prefix)
            priority = stage_priority(stage_dir, pdbqt)
            candidates_by_base.setdefault(base, []).append((priority, pdbqt, stage_dir))

    available_bases = set(candidates_by_base.keys())
    if allowed_bases is not None:
        filtered_candidates = {
            base: entries
            for base, entries in candidates_by_base.items()
            if base in allowed_bases
        }
    else:
        filtered_candidates = candidates_by_base

    stage_counts: Dict[int, int] = {}
    rescored_stage_by_base: Dict[str, str] = {}
    stage_fallback_reason_by_base: Dict[str, str] = {}
    chosen: Dict[str, Tuple[int, Path, str]] = {}
    for base, entries in filtered_candidates.items():
        ordered_entries = sorted(entries, key=lambda t: (-t[0], str(t[1])))
        preferred_stage = preferred_stage_by_base.get(base, "")
        preferred_priority = stage_priority(preferred_stage) if preferred_stage else 0
        selected_entry: Optional[Tuple[int, Path, str]] = None
        if preferred_priority > 0:
            matching = [entry for entry in ordered_entries if entry[0] == preferred_priority]
            if matching:
                selected_entry = sorted(matching, key=lambda t: str(t[1]))[0]
            else:
                selected_entry = ordered_entries[0]
                stage_fallback_reason_by_base[base] = f"preferred_missing:{preferred_stage}"
        else:
            selected_entry = ordered_entries[0]
        chosen[base] = selected_entry
        rescored_stage_by_base[base] = _resolved_stage_label(
            selected_entry[2],
            selected_entry[1],
            fallback=preferred_stage or "",
        )
    for priority, _path, _stage_dir in chosen.values():
        stage_counts[priority] = stage_counts.get(priority, 0) + 1

    ligands = [
        entry[1]
        for entry in sorted(chosen.values(), key=lambda t: (-t[0], str(t[1])))
    ]
    logger.debug(
        "%s action=collect_best status=ok stage_dirs=%s candidates=%d selected=%d stage3=%d stage2=%d stage1=%d stage0=%d",
        component,
        ",".join(stage_dirs),
        total_candidates,
        len(ligands),
        stage_counts.get(3, 0),
        stage_counts.get(2, 0),
        stage_counts.get(1, 0),
        stage_counts.get(0, 0),
    )
    return PoseSelectionResult(
        ligands=ligands,
        stage_counts=stage_counts,
        total_candidates=total_candidates,
        available_bases=available_bases,
        rescored_stage_by_base=rescored_stage_by_base,
        stage_fallback_reason_by_base=stage_fallback_reason_by_base,
    )


def load_consensus_top_bases(
    consensus_csv: Path,
    frac: float,
    logger: logging.Logger,
    *,
    component: str,
    control_bases: Optional[Set[str]] = None,
    decoy_prefix: str,
) -> Tuple[Set[str], Set[str], int, int, int, int, int]:
    rows: List[Dict[str, str]] = []
    try:
        with consensus_csv.open() as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames or []
            missing = {"ligand", "consensus_score"} - set(fieldnames)
            if missing:
                logger.error(
                    "%s action=select status=failed reason=missing_fields path=%s missing=%s",
                    component,
                    consensus_csv,
                    ",".join(sorted(missing)),
                )
                return set(), set(), 0, 0, 0, 0, 0
            for row in reader:
                rows.append(row)
    except Exception as exc:
        logger.error(
            "%s action=select status=failed reason=read_error path=%s error=%s",
            component,
            consensus_csv,
            exc,
        )
        return set(), set(), 0, 0, 0, 0, 0

    if not rows:
        return set(), set(), 0, 0, 0, 0, 0

    def _score(row: Dict[str, str]) -> float:
        try:
            return float(row.get("consensus_score", ""))
        except Exception:
            return float("-inf")

    control_set = control_bases or set()
    scored_rows: List[Tuple[float, str, bool, str]] = []
    controls_in_consensus: Set[str] = set()
    has_library = "library" in (rows[0].keys() if rows else [])
    unknown_count = 0
    for row in rows:
        lig = str(row.get("ligand", "")).strip()
        if not lig:
            continue
        base = pose_base_from_path(Path(lig), decoy_prefix=decoy_prefix)
        if not base:
            continue
        is_ctrl = base in control_set
        if is_ctrl:
            controls_in_consensus.add(base)
        lib = str(row.get("library", "")).strip().upper() if has_library else ""
        if lib == "":
            lib = "UNKNOWN"
        scored_rows.append((_score(row), base, is_ctrl, lib))
        if lib == "UNKNOWN" and not is_ctrl:
            unknown_count += 1

    if not scored_rows:
        empty = set(control_set)
        return (
            empty,
            empty,
            len(rows),
            0,
            0,
            len(control_set),
            len(controls_in_consensus),
        )

    if not has_library:
        non_controls = [
            (score, base) for score, base, is_ctrl, _ in scored_rows if not is_ctrl
        ]
        non_controls_sorted = sorted(
            non_controls, key=lambda item: item[0], reverse=True
        )
        k = (
            max(1, math.ceil(len(non_controls_sorted) * frac))
            if non_controls_sorted
            else 0
        )

        allowed: Set[str] = set(control_set)
        for _, base in non_controls_sorted[:k]:
            allowed.add(base)

        logger.info(
            "[select.top_frac] library=absent total_rows=%d non_controls=%d k=%d controls=%d controls_in_consensus=%d frac=%.3f",
            len(rows),
            len(non_controls_sorted),
            k,
            len(control_set),
            len(controls_in_consensus),
            frac,
        )
        return (
            allowed,
            allowed,
            len(rows),
            k,
            len(non_controls_sorted),
            len(control_set),
            len(controls_in_consensus),
        )

    buckets: Dict[str, List[Tuple[float, str]]] = {"FDA": [], "DECOY": []}
    non_controls_total = 0
    for score, base, is_ctrl, lib in scored_rows:
        if is_ctrl:
            continue
        non_controls_total += 1
        if lib in buckets:
            buckets[lib].append((score, base))
    allowed_fda: Set[str] = set(control_set)
    allowed_decoy: Set[str] = set(control_set)
    k_fda = 0
    k_decoy = 0
    for lib, bucket in buckets.items():
        sorted_bucket = sorted(bucket, key=lambda item: item[0], reverse=True)
        k = max(1, math.ceil(len(sorted_bucket) * frac)) if sorted_bucket else 0
        if lib == "FDA":
            k_fda = k
        elif lib == "DECOY":
            k_decoy = k
        for _, base in sorted_bucket[:k]:
            if lib == "FDA":
                allowed_fda.add(base)
            elif lib == "DECOY":
                allowed_decoy.add(base)

    logger.info(
        "[select.top_frac] library=stratified total_rows=%d non_controls=%d non_controls_fda=%d k_fda=%d non_controls_decoy=%d k_decoy=%d controls=%d controls_in_consensus=%d unknown=%d frac=%.3f",
        len(rows),
        non_controls_total,
        len(buckets["FDA"]),
        k_fda,
        len(buckets["DECOY"]),
        k_decoy,
        len(control_set),
        len(controls_in_consensus),
        unknown_count,
        frac,
    )
    return (
        allowed_fda,
        allowed_decoy,
        len(rows),
        k_fda + k_decoy,
        non_controls_total,
        len(control_set),
        len(controls_in_consensus),
    )


def select_top_bases_from_score_csv(
    score_csv: Path,
    frac: float,
    control_bases: Optional[Set[str]],
    *,
    higher_is_better: bool,
    logger: logging.Logger,
    component: str,
    score_cols: Optional[Sequence[str]] = None,
    decoy_prefix: str,
) -> SelectionResult:
    control_set = control_bases or set()
    score_cols_list = list(score_cols or [])
    preserve_input_order = CONSENSUS_INPUT_ORDER_SENTINEL in score_cols_list
    score_cols_effective = [
        col for col in score_cols_list if col != CONSENSUS_INPUT_ORDER_SENTINEL
    ]
    if not score_csv.exists():
        logger.info(
            "%s action=select status=skip reason=missing_score_csv path=%s",
            component,
            score_csv,
        )
        return SelectionResult(
            allowed_bases=set(control_set),
            n_pool=0,
            k=0,
            controls_total=len(control_set),
            selected_stage_by_base={},
            selected_score_by_base={},
        )

    ligand_keys = (
        "ligand",
        "Ligand",
        "ligand_file",
        "Ligand_file",
        "Ligand_ID",
        "ligand_id",
        "ligand_name",
        "name",
        "molecule",
        "Molecule",
    )
    best_by_base: Dict[str, float] = {}
    best_stage_by_base: Dict[str, str] = {}
    noncontrol_seen: Set[str] = set()
    ordered_seen: List[str] = []
    n_rows = 0
    try:
        with score_csv.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames or []
            for row in reader:
                n_rows += 1
                lig_val = ""
                for key in ligand_keys:
                    if row.get(key):
                        lig_val = str(row.get(key, "")).strip()
                        break
                if not lig_val:
                    continue
                base = pose_base_from_path(Path(lig_val), decoy_prefix=decoy_prefix)
                if not base:
                    continue
                if base in control_set:
                    continue
                if base not in noncontrol_seen:
                    ordered_seen.append(base)
                    noncontrol_seen.add(base)

                row_candidates: List[Tuple[float, str]] = []
                if score_cols_effective:
                    for col in score_cols_effective:
                        raw = row.get(col, "")
                        try:
                            val = float(str(raw).strip())
                        except Exception:
                            continue
                        if math.isfinite(val):
                            row_candidates.append(
                                (
                                    val,
                                    _normalize_stage_label(col, lig_val),
                                )
                            )
                else:
                    for col in fields:
                        if col in ligand_keys or col in {
                            "run_id",
                            "variant",
                            "stage",
                            "valid",
                            "reason",
                        }:
                            continue
                        raw = row.get(col, "")
                        try:
                            val = float(str(raw).strip())
                        except Exception:
                            continue
                        if math.isfinite(val):
                            row_candidates.append(
                                (
                                    val,
                                    _normalize_stage_label(row.get("stage", ""), lig_val),
                                )
                            )
                if not row_candidates and not preserve_input_order:
                    continue
                if row_candidates:
                    if higher_is_better:
                        score_val, stage_label = max(
                            row_candidates, key=lambda item: item[0]
                        )
                    else:
                        score_val, stage_label = min(
                            row_candidates, key=lambda item: item[0]
                        )
                    current = best_by_base.get(base)
                    if current is None:
                        best_by_base[base] = score_val
                        best_stage_by_base[base] = stage_label
                    else:
                        if higher_is_better and score_val > current:
                            best_by_base[base] = score_val
                            best_stage_by_base[base] = stage_label
                        elif (not higher_is_better) and score_val < current:
                            best_by_base[base] = score_val
                            best_stage_by_base[base] = stage_label
    except Exception as exc:
        logger.warning(
            "%s action=select status=skip reason=read_error path=%s error=%s",
            component,
            score_csv,
            exc,
        )
        return SelectionResult(
            allowed_bases=set(control_set),
            n_pool=0,
            k=0,
            controls_total=len(control_set),
            selected_stage_by_base={},
            selected_score_by_base={},
        )

    pool_bases = sorted(noncontrol_seen)
    k = max(1, math.ceil(len(pool_bases) * frac)) if pool_bases else 0

    if preserve_input_order:
        ordered_bases = [base for base in ordered_seen if base in noncontrol_seen]
    else:
        ordered_scored = sorted(
            best_by_base.items(),
            key=lambda item: item[1],
            reverse=higher_is_better,
        )
        ordered_bases = [base for base, _ in ordered_scored]
        if len(ordered_bases) < len(pool_bases):
            remaining = [base for base in pool_bases if base not in best_by_base]
            ordered_bases.extend(remaining)

    allowed: Set[str] = set(control_set)
    selected_stage_subset: Dict[str, str] = {}
    selected_score_subset: Dict[str, float] = {}
    for base in ordered_bases[:k]:
        allowed.add(base)
        if base in best_stage_by_base and best_stage_by_base[base]:
            selected_stage_subset[base] = best_stage_by_base[base]
        if base in best_by_base:
            selected_score_subset[base] = best_by_base[base]

    logger.info(
        "%s action=select status=ok path=%s rows=%d n_pool=%d k=%d controls=%d allowed_total=%d scored=%d",
        component,
        score_csv,
        n_rows,
        len(pool_bases),
        k,
        len(control_set),
        len(allowed),
        len(best_by_base),
    )
    return SelectionResult(
        allowed_bases=allowed,
        n_pool=len(pool_bases),
        k=k,
        controls_total=len(control_set),
        selected_stage_by_base=selected_stage_subset,
        selected_score_by_base=selected_score_subset,
    )


def score_csv_for_spec(
    run_root: Path,
    combo: Tuple[str, str, str],
    spec: StageSpec,
    run_mode: str,
    *,
    decoy_prefix: str,
) -> Tuple[Optional[Path], List[str], bool]:
    pdb_id, variant, ph = combo
    combo_root = run_root / pdb_id / variant / ph
    dud = run_mode == "dud"
    prefix = decoy_prefix
    higher_is_better = False

    if dud:
        for candidate in (
            combo_root / f"{prefix}_consensus_docking_scores.csv",
            combo_root / "dud_consensus_docking_scores.csv",
            combo_root / "decoys_consensus_docking_scores.csv",
        ):
            if candidate.exists():
                return (
                    candidate,
                    ["consensus_score"],
                    True,
                )
        generic_consensus = combo_root / "consensus_docking_scores.csv"
        if (
            generic_consensus.exists()
            and _generic_consensus_stream(combo_root, decoy_prefix=decoy_prefix)
            == "dud"
        ):
            return (
                generic_consensus,
                ["consensus_score"],
                True,
            )
    else:
        consensus = combo_root / "consensus_docking_scores.csv"
        if consensus.exists():
            return (
                consensus,
                ["consensus_score"],
                True,
            )

    if spec.source == "vina":
        aliases = _decoy_prefix_aliases(prefix) if dud else (prefix,)
        summary = combo_root / "docking_score_summary.csv"
        long_csv = combo_root / "docking_score_long.csv"
        if dud:
            summary = next(
                (
                    combo_root / f"{alias}_docking_score_summary.csv"
                    for alias in aliases
                    if (combo_root / f"{alias}_docking_score_summary.csv").exists()
                ),
                combo_root / f"{prefix}_docking_score_summary.csv",
            )
            long_csv = next(
                (
                    combo_root / f"{alias}_docking_score_long.csv"
                    for alias in aliases
                    if (combo_root / f"{alias}_docking_score_long.csv").exists()
                ),
                combo_root / f"{prefix}_docking_score_long.csv",
            )
        if summary.exists():
            selected_prefix = summary.name[: -len("_docking_score_summary.csv")]
            cols = (
                [
                    f"{selected_prefix}_stage3",
                    f"{selected_prefix}_stage2",
                    f"{selected_prefix}_stage1",
                ]
                if dud
                else ["stage3", "stage2", "stage1"]
            )
            return summary, cols, higher_is_better
        return long_csv, ["score"], higher_is_better

    if spec.source == "gnina":
        aliases = _decoy_prefix_aliases(prefix) if dud else (prefix,)
        summary = combo_root / "gnina_docking_score_summary.csv"
        long_csv = combo_root / "gnina_docking_score_long.csv"
        if dud:
            summary = next(
                (
                    combo_root / f"{alias}_gnina_docking_score_summary.csv"
                    for alias in aliases
                    if (combo_root / f"{alias}_gnina_docking_score_summary.csv").exists()
                ),
                combo_root / f"{prefix}_gnina_docking_score_summary.csv",
            )
            long_csv = next(
                (
                    combo_root / f"{alias}_gnina_docking_score_long.csv"
                    for alias in aliases
                    if (combo_root / f"{alias}_gnina_docking_score_long.csv").exists()
                ),
                combo_root / f"{prefix}_gnina_docking_score_long.csv",
            )
        if summary.exists():
            selected_prefix = summary.name[: -len("_gnina_docking_score_summary.csv")]
            stage_prefix = f"gnina_{selected_prefix}_stage" if dud else "gnina_stage"
            cols = [f"{stage_prefix}3", f"{stage_prefix}2", f"{stage_prefix}1"]
            return summary, cols, higher_is_better
        return (
            long_csv,
            ["gnina_primary_score", "gnina_minimized_affinity_kcal"],
            higher_is_better,
        )

    if spec.source == "ledock":
        aliases = _decoy_prefix_aliases(prefix) if dud else (prefix,)
        summary = combo_root / "ledock_docking_score_summary.csv"
        long_csv = combo_root / "ledock_docking_score_long.csv"
        if dud:
            summary = next(
                (
                    combo_root / f"{alias}_ledock_docking_score_summary.csv"
                    for alias in aliases
                    if (combo_root / f"{alias}_ledock_docking_score_summary.csv").exists()
                ),
                combo_root / f"{prefix}_ledock_docking_score_summary.csv",
            )
            long_csv = next(
                (
                    combo_root / f"{alias}_ledock_docking_score_long.csv"
                    for alias in aliases
                    if (combo_root / f"{alias}_ledock_docking_score_long.csv").exists()
                ),
                combo_root / f"{prefix}_ledock_docking_score_long.csv",
            )
        if summary.exists():
            selected_prefix = summary.name[: -len("_ledock_docking_score_summary.csv")]
            cols = (
                [
                    f"{selected_prefix}_stage3",
                    f"{selected_prefix}_stage2",
                    f"{selected_prefix}_stage1",
                ]
                if dud
                else ["stage3", "stage2", "stage1"]
            )
            return summary, cols, higher_is_better
        return long_csv, ["ledock_best_score_kcal"], higher_is_better

    if spec.source == "dock6":
        aliases = _decoy_prefix_aliases(prefix) if dud else (prefix,)
        summary = combo_root / "dock6_docking_score_summary.csv"
        long_csv = combo_root / "dock6_docking_score_long.csv"
        if dud:
            summary = next(
                (
                    combo_root / f"{alias}_dock6_docking_score_summary.csv"
                    for alias in aliases
                    if (combo_root / f"{alias}_dock6_docking_score_summary.csv").exists()
                ),
                combo_root / f"{prefix}_dock6_docking_score_summary.csv",
            )
            long_csv = next(
                (
                    combo_root / f"{alias}_dock6_docking_score_long.csv"
                    for alias in aliases
                    if (combo_root / f"{alias}_dock6_docking_score_long.csv").exists()
                ),
                combo_root / f"{prefix}_dock6_docking_score_long.csv",
            )
        if summary.exists():
            selected_prefix = summary.name[: -len("_dock6_docking_score_summary.csv")]
            cols = (
                [
                    f"{selected_prefix}_stage3",
                    f"{selected_prefix}_stage2",
                    f"{selected_prefix}_stage1",
                ]
                if dud
                else ["stage3", "stage2", "stage1"]
            )
            return summary, cols, higher_is_better
        return long_csv, ["dock6_grid_score"], higher_is_better

    return None, [], higher_is_better


def discover_mode_dirs(
    combo: Tuple[str, str, str],
    run_root: Path,
    post_root: Path,
    specs: Sequence[StageSpec],
    logger: logging.Logger,
    *,
    decoy_prefix: str,
) -> Dict[str, List[Path]]:
    pdb_id, variant, ph = combo
    mode_dirs: Dict[str, Set[Path]] = {"fda": set(), "dud": set()}
    for spec in specs:
        root = (
            run_root / pdb_id / variant / ph
            if spec.source in {"vina", "gnina"}
            else post_root / pdb_id / variant / ph
        )
        for mode in ("fda", "dud"):
            for stage_dir in stage_dir_candidates(
                spec.source,
                mode,
                root,
                decoy_prefix=decoy_prefix,
            ):
                candidate = root / stage_dir
                if candidate.exists():
                    mode_dirs[mode].add(candidate)
    generic_stream = _generic_consensus_stream(
        run_root / pdb_id / variant / ph,
        decoy_prefix=decoy_prefix,
    )
    if generic_stream == "dud":
        generic_dirs: Set[Path] = set()
        for spec in specs:
            root = (
                run_root / pdb_id / variant / ph
                if spec.source in {"vina", "gnina"}
                else post_root / pdb_id / variant / ph
            )
            for stage_dir in stage_dir_candidates(
                spec.source,
                "fda",
                root,
                decoy_prefix=decoy_prefix,
            ):
                candidate = root / stage_dir
                if candidate.exists():
                    generic_dirs.add(candidate)
        if generic_dirs:
            mode_dirs["dud"].update(generic_dirs)
            mode_dirs["fda"].difference_update(generic_dirs)
    fda_list = sorted(mode_dirs["fda"])
    dud_list = sorted(mode_dirs["dud"])
    logger.info(
        "[scorch.discover] pdb_id=%s variant=%s ph=%s fda=%s dud=%s",
        pdb_id,
        variant,
        ph,
        ",".join(str(p) for p in fda_list) if fda_list else "None",
        ",".join(str(p) for p in dud_list) if dud_list else "None",
    )
    return {"fda": fda_list, "dud": dud_list}


def log_score_csv_coverage(
    source: str,
    combo: Tuple[str, str, str],
    score_csv: Path,
    allowed_bases: Set[str],
    logger: logging.Logger,
    *,
    component: str,
    decoy_prefix: str,
) -> None:
    if not score_csv.exists():
        logger.debug(
            "%s action=select source=%s pdb_id=%s variant=%s ph=%s status=skip reason=missing_score_csv path=%s",
            component,
            source,
            combo[0],
            combo[1],
            combo[2],
            score_csv,
        )
        return
    try:
        with score_csv.open() as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or "ligand" not in reader.fieldnames:
                logger.debug(
                    "%s action=select source=%s pdb_id=%s variant=%s ph=%s status=skip reason=missing_ligand_column path=%s",
                    component,
                    source,
                    combo[0],
                    combo[1],
                    combo[2],
                    score_csv,
                )
                return
            bases_in_csv = {
                base
                for row in reader
                for base in [
                    pose_base_from_path(
                        Path(str(row.get("ligand", "")).strip()),
                        decoy_prefix=decoy_prefix,
                    )
                ]
                if row.get("ligand") and base
            }
    except Exception as exc:
        logger.debug(
            "%s action=select source=%s pdb_id=%s variant=%s ph=%s status=skip reason=read_error path=%s error=%s",
            component,
            source,
            combo[0],
            combo[1],
            combo[2],
            score_csv,
            exc,
        )
        return

    covered = len(allowed_bases & bases_in_csv)
    logger.info(
        "%s action=select source=%s pdb_id=%s variant=%s ph=%s allowed=%d score_rows=%d covered_by_scores=%d path=%s | score CSV coverage only",
        component,
        source,
        combo[0],
        combo[1],
        combo[2],
        len(allowed_bases),
        len(bases_in_csv),
        covered,
        score_csv,
    )
