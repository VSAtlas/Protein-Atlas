from __future__ import annotations

import csv
import json
import logging
import math
import re
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple


def _mmgbsa_normalize_ligand_base(name: str) -> str:
    s = str(name or "").strip()
    if not s:
        return ""
    s = re.sub(r"\.(pdbqt|mol2|sdf)$", "", s, flags=re.IGNORECASE)
    s = s.replace(".sanitized", "")
    return s


def _mmgbsa_sdf_base_candidates(stem: str, stage_dir_label: str) -> List[str]:
    candidates: set[str] = set()
    base = stem.replace(".sanitized", "")
    candidates.add(base)

    if stage_dir_label:
        for sep in ("_", "__"):
            suffix = f"{sep}{stage_dir_label}"
            if base.endswith(suffix):
                candidates.add(base[: -len(suffix)])

    for pattern in (
        r"(_+gnina_stage\d+)$",
        r"(_+ledock_stage\d+)$",
        r"(_+dock6_stage\d+)$",
        r"(_+stage\d+)$",
    ):
        stripped = re.sub(pattern, "", base, flags=re.IGNORECASE)
        if stripped != base:
            candidates.add(stripped)

    cleaned = set()
    for cand in candidates:
        cleaned.add(re.sub(r"_+$", "", cand))
    return sorted(cleaned)


def _mmgbsa_default_pose_group_regexes() -> List[str]:
    return [
        r"(_pose\d+)$",
        r"(_rank\d+)$",
        r"(_conf\d+)$",
        r"(_model\d+)$",
        r"(_p\d+)$",
        r"(_stage\d+)(_pose\d+)$",
        r"(_gnina_stage\d+)(_pose\d+)$",
        r"(_ledock_stage\d+)(_pose\d+)$",
        r"(_dock6_stage\d+)(_pose\d+)$",
    ]


def _mmgbsa_pose_group_regexes(
    cfg: Mapping[str, Any], logger: logging.Logger
) -> List[re.Pattern]:
    raw = str(cfg.get("MMGBSA_POSE_GROUP_REGEXES", "") or "").strip()
    patterns: List[str] = []
    if raw:
        for part in raw.split(","):
            part = part.strip()
            if part:
                patterns.append(part)
    if not patterns:
        patterns = _mmgbsa_default_pose_group_regexes()

    compiled: List[re.Pattern] = []
    for pat in patterns:
        try:
            compiled.append(re.compile(pat, flags=re.IGNORECASE))
        except re.error as exc:
            logger.warning(
                "[mmgbsa.pipeline] selection=pose_regex_invalid regex=%s err=%s",
                pat,
                exc,
            )
    if not compiled:
        compiled = [
            re.compile(pat, flags=re.IGNORECASE)
            for pat in _mmgbsa_default_pose_group_regexes()
        ]
    return compiled


def _mmgbsa_pose_candidates(
    stem: str, stage_dir_label: str, regexes: List[re.Pattern]
) -> List[str]:
    base = stem.replace(".sanitized", "")
    seeds: set[str] = {base}
    for regex in regexes:
        stripped = regex.sub("", base)
        if stripped != base:
            seeds.add(stripped)

    candidates: set[str] = set()
    for seed in seeds:
        for cand in _mmgbsa_sdf_base_candidates(seed, stage_dir_label):
            candidates.add(cand)
        for regex in regexes:
            stripped = regex.sub("", seed)
            if stripped != seed:
                for cand in _mmgbsa_sdf_base_candidates(stripped, stage_dir_label):
                    candidates.add(cand)

    cleaned: set[str] = set()
    for cand in candidates:
        cleaned.add(re.sub(r"_+$", "", cand))
    return sorted([c for c in cleaned if c])


def _mmgbsa_canonical_ligand_id(
    stem: str,
    stage_dir_label: str,
    regexes: List[re.Pattern],
    logger: logging.Logger,
) -> str:
    candidates = _mmgbsa_pose_candidates(stem, stage_dir_label, regexes)
    if not candidates:
        return _mmgbsa_normalize_ligand_base(stem)

    min_len = min(len(cand) for cand in candidates)
    shortest = sorted([cand for cand in candidates if len(cand) == min_len])
    chosen = shortest[0]
    if len(shortest) > 1:
        logger.warning(
            "[mmgbsa.pipeline] selection=ligand_id_collision stem=%s candidates=%s chosen=%s",
            stem,
            shortest,
            chosen,
        )
    return chosen


def _mmgbsa_group_pose_sdfs(
    sdfs: List[Path],
    stage_dir_label: str,
    regexes: List[re.Pattern],
    logger: logging.Logger,
) -> Dict[str, List[Path]]:
    grouped: Dict[str, List[Path]] = {}
    for sdf in sdfs:
        ligand_id = _mmgbsa_canonical_ligand_id(
            sdf.stem, stage_dir_label, regexes, logger
        )
        if not ligand_id:
            ligand_id = sdf.stem
        grouped.setdefault(ligand_id, []).append(sdf)
    return grouped


def _mmgbsa_pose_numeric_rank(stem: str) -> Optional[int]:
    primary = re.search(r"(?:pose|rank|conf|model|p)(\d+)$", stem, flags=re.IGNORECASE)
    if primary:
        try:
            return int(primary.group(1))
        except Exception:
            return None

    trailing = re.search(r"(\d+)$", stem)
    if trailing:
        try:
            return int(trailing.group(1))
        except Exception:
            return None
    return None


def _mmgbsa_authoritative_pose_rank(path: Path) -> int:
    meta_path = path.with_suffix(".mmgbsa_pose.json")
    if not meta_path.is_file():
        return 1
    try:
        payload = json.loads(meta_path.read_text(encoding="utf-8"))
    except Exception:
        return 1
    if not isinstance(payload, dict):
        return 1
    authoritative = payload.get("chemistry_authoritative") is True
    reconstructed = payload.get("pdbqt_to_sdf_reconstruction") is True
    return 0 if authoritative and not reconstructed else 1


def _mmgbsa_pose_sort_key(path: Path, mode: str) -> Tuple[int, int, object, str]:
    stem = path.stem.lower()
    authority_rank = _mmgbsa_authoritative_pose_rank(path)
    if mode == "name_numeric":
        rank = _mmgbsa_pose_numeric_rank(path.stem)
        if rank is not None:
            return (authority_rank, 0, rank, stem)
        return (authority_rank, 1, stem, stem)
    return (authority_rank, 0, stem, stem)


def _mmgbsa_load_reranked_bases(csv_path: Path, logger: logging.Logger) -> List[str]:
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return []

    rows: List[Tuple[Optional[int], int, str]] = []
    try:
        with csv_path.open(
            "r", encoding="utf-8", errors="ignore", newline=""
        ) as handle:
            reader = csv.DictReader(handle)
            for idx, row in enumerate(reader):
                lig_base = row.get("ligand_base") or _mmgbsa_normalize_ligand_base(
                    row.get("ligand", "")
                )
                if not lig_base:
                    continue
                rank_val = None
                rank_raw = str(row.get("final_rank", "")).strip()
                if rank_raw:
                    try:
                        rank_val = int(float(rank_raw))
                    except Exception:
                        rank_val = None
                rows.append((rank_val, idx, lig_base))
    except Exception as exc:
        logger.warning(
            "[mmgbsa.pipeline] selection=reranked_top_pct failed reason=csv_read_error path=%s err=%s",
            csv_path,
            exc,
        )
        return []

    if not rows:
        return []

    rows.sort(key=lambda x: (x[0] if x[0] is not None else 1_000_000_000, x[1]))
    seen: set[str] = set()
    ordered: List[str] = []
    for _, _, base in rows:
        if base in seen:
            continue
        seen.add(base)
        ordered.append(base)
    return ordered


def _mmgbsa_select_sdfs(
    stage_dir: Path,
    max_ligands: Optional[int],
    poses_per_ligand: int,
    pose_sort_mode: str,
    pose_group_regexes: List[re.Pattern],
    reranked_csv: Optional[Path] = None,
    top_pct: Optional[float] = None,
    selected_ligand_ids: Optional[Sequence[str]] = None,
    logger: Optional[logging.Logger] = None,
) -> Tuple[List[Path], bool, bool, List[str]]:
    sdfs = sorted(stage_dir.glob("*.sdf"))
    if not sdfs:
        return [], False, False, []

    if logger is None:
        logger = logging.getLogger("mmgbsa.pipeline")

    used_reranked = False
    used_actives = False
    stage_label = stage_dir.name
    pose_groups = _mmgbsa_group_pose_sdfs(sdfs, stage_label, pose_group_regexes, logger)
    selected_group_map = pose_groups
    selected_ligands: List[str] = []
    requested_ligands = [
        str(item).strip()
        for item in (selected_ligand_ids or [])
        if str(item).strip()
    ]
    explicit_selection = bool(requested_ligands)

    if explicit_selection:
        requested_ids: set[str] = set()
        for token in requested_ligands:
            requested_ids.add(_mmgbsa_normalize_ligand_base(token))
            requested_ids.add(
                _mmgbsa_canonical_ligand_id(
                    Path(token).stem, stage_label, pose_group_regexes, logger
                )
            )
        selected_ligands = [
            ligand_id for ligand_id in sorted(pose_groups) if ligand_id in requested_ids
        ]
        missing = sorted(requested_ids - set(selected_ligands))
        if missing:
            logger.warning(
                "[mmgbsa.pipeline] selection=explicit_ligands_missing stage_dir=%s missing=%s",
                stage_label,
                ",".join(missing),
            )
        logger.info(
            "[mmgbsa.pipeline] selection=explicit_ligands stage_dir=%s requested=%d matched=%d",
            stage_label,
            len(requested_ligands),
            len(selected_ligands),
        )

    if (
        not selected_ligands
        and not explicit_selection
        and reranked_csv is not None
        and top_pct is not None
        and top_pct > 0
    ):
        if not reranked_csv.exists():
            logger.warning(
                "[mmgbsa.pipeline] selection=reranked_top_pct failed reason=missing_reranked_csv path=%s",
                reranked_csv,
            )
        else:
            ordered_bases = _mmgbsa_load_reranked_bases(reranked_csv, logger)
            if ordered_bases:
                pct_val = float(top_pct)
                if pct_val > 100.0:
                    pct_val = 100.0
                if pct_val <= 0.0:
                    pct_val = 0.0
                ordered_ligands: List[str] = []
                seen: set[str] = set()
                for base in ordered_bases:
                    ligand_id = _mmgbsa_canonical_ligand_id(
                        base, stage_label, pose_group_regexes, logger
                    )
                    if ligand_id in seen or ligand_id not in pose_groups:
                        continue
                    seen.add(ligand_id)
                    ordered_ligands.append(ligand_id)

                total = len(ordered_ligands)
                target = max(1, int(math.ceil(total * pct_val / 100.0))) if total else 0
                if ordered_ligands and target:
                    selected_ligands = ordered_ligands[:target]
                    used_reranked = True
                    logger.info(
                        "[mmgbsa.pipeline] selection=reranked_top_pct pct=%.3g total=%d target=%d selected=%d path=%s",
                        pct_val,
                        total,
                        target,
                        len(selected_ligands),
                        reranked_csv,
                    )
                else:
                    logger.warning(
                        "[mmgbsa.pipeline] selection=reranked_top_pct failed reason=no_matching_sdfs path=%s stage_dir=%s",
                        reranked_csv,
                        stage_dir,
                    )
            else:
                logger.warning(
                    "[mmgbsa.pipeline] selection=reranked_top_pct failed reason=empty_csv path=%s",
                    reranked_csv,
                )

    if not selected_ligands and not explicit_selection:
        actives = [p for p in sdfs if "actives_final" in p.name]
        if actives:
            used_actives = True
            selected_group_map = _mmgbsa_group_pose_sdfs(
                actives, stage_label, pose_group_regexes, logger
            )
            selected_ligands = sorted(selected_group_map.keys())
        else:
            selected_group_map = pose_groups
            selected_ligands = sorted(pose_groups.keys())

    if max_ligands is not None and max_ligands > 0:
        selected_ligands = selected_ligands[:max_ligands]

    selected_sdfs: List[Path] = []
    for ligand_id in selected_ligands:
        group = selected_group_map.get(ligand_id, [])
        if not group:
            continue
        sorted_group = sorted(
            group, key=lambda p: _mmgbsa_pose_sort_key(p, pose_sort_mode)
        )
        if poses_per_ligand > 0:
            sorted_group = sorted_group[:poses_per_ligand]
        selected_sdfs.extend(sorted_group)

    logger.info(
        "[mmgbsa.pipeline] selection=poses stage_dir=%s ligand_ids=%d poses_per_ligand=%d poses_selected=%d reranked=%s actives_only=%s",
        stage_label,
        len(selected_ligands),
        poses_per_ligand,
        len(selected_sdfs),
        used_reranked,
        used_actives,
    )
    return selected_sdfs, used_reranked, used_actives, selected_ligands
