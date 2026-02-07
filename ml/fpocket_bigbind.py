from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
from pathlib import Path
import threading
from datetime import datetime, timezone
from typing import Any

import pandas as pd

from druggability_evaluation import run_fpocket_for_explicit_pocket
from druggability_orchestrator import load_fpocket_metrics_for_receptor_pdb
from ml.data.bigbind import resolve_pocket_path
from ml.pocket_features import POCKET_FEATURE_COLUMNS, compute_pocket_features


FPOCKET_METRIC_COLUMNS = (
    "fpocket_druggability",
    "fpocket_volume",
    "fpocket_openness",
    "fpocket_polar_fraction",
    "fpocket_has_metal",
    "fpocket_tier",
)
FPOCKET_CENTER_COLUMNS = ("pocket_center_x", "pocket_center_y", "pocket_center_z")


def _coerce_positive_int(value: Any) -> int | None:
    if value in (None, "", "none", "null"):
        return None
    try:
        parsed = int(value)
    except Exception:
        return None
    return parsed if parsed > 0 else None


def _resolve_worker_count(
    atlas_cfg: dict[str, Any],
    *,
    max_workers: int | None = None,
) -> int:
    explicit = _coerce_positive_int(max_workers)
    if explicit is not None:
        return explicit
    for key in ("ML_FPOCKET_WORKERS", "ML_WORKERS", "CPU"):
        parsed = _coerce_positive_int(atlas_cfg.get(key))
        if parsed is not None:
            return parsed
    return 1


def _coerce_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, (int, float)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "no", "n", "off", "", "none", "null"}:
        return False
    return bool(text)


def _write_no_pocket_sentinel(
    sentinel_path: Path,
    *,
    receptor_pdb: Path,
    reason: str,
) -> None:
    sentinel_path.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).isoformat()
    sentinel_path.write_text(
        f"created_utc={timestamp}\nreceptor={receptor_pdb}\nreason={reason}\n",
        encoding="utf-8",
    )


def parse_pocket_residues(pocket_pdb: Path) -> list[tuple[int, str, str]]:
    if not pocket_pdb.exists():
        return []

    residues: set[tuple[int, str, str]] = set()
    with pocket_pdb.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            chain = (line[21:22] or " ").strip() or " "
            try:
                resseq = int((line[22:26] or "0").strip() or "0")
            except ValueError:
                continue
            icode = (line[26:27] or "").strip() or "-"
            residues.add((resseq, icode, chain))
    return sorted(residues, key=lambda entry: (entry[2], entry[0], entry[1]))


def _infer_chain_by_residue_from_receptor(receptor_pdb: Path) -> dict[tuple[int, str], str]:
    if not receptor_pdb.exists():
        return {}
    by_key: dict[tuple[int, str], set[str]] = {}
    with receptor_pdb.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                resseq = int((line[22:26] or "0").strip() or "0")
            except ValueError:
                continue
            icode = (line[26:27] or "").strip() or "-"
            chain = (line[21:22] or " ").strip() or " "
            key = (resseq, icode)
            if key not in by_key:
                by_key[key] = set()
            by_key[key].add(chain)
    inferred: dict[tuple[int, str], str] = {}
    for key, chains in by_key.items():
        ordered = sorted(chains)
        if ordered:
            inferred[key] = ordered[0]
    return inferred


def _normalize_pocket_residues_for_receptor(
    residues: list[tuple[int, str, str]],
    receptor_pdb: Path,
) -> list[tuple[int, str, str]]:
    if not residues:
        return []
    chain_by_residue = _infer_chain_by_residue_from_receptor(receptor_pdb)
    normalized: set[tuple[int, str, str]] = set()
    for resseq, icode, chain in residues:
        icode_norm = (icode or "").strip() or "-"
        chain_norm = (chain or "").strip()
        if not chain_norm or chain_norm == "-":
            chain_norm = (
                chain_by_residue.get((resseq, icode_norm))
                or chain_by_residue.get((resseq, "-"))
                or "A"
            )
        normalized.add((resseq, icode_norm, chain_norm))
    return sorted(normalized, key=lambda entry: (entry[2], entry[0], entry[1]))


def pocket_center_from_pdb(pocket_pdb: Path) -> tuple[float, float, float] | None:
    if not pocket_pdb.exists():
        return None

    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    with pocket_pdb.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            try:
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
            except Exception:
                continue
    if not xs:
        return None
    n = float(len(xs))
    return (sum(xs) / n, sum(ys) / n, sum(zs) / n)


def derive_receptor_pdb(pocket_pdb: Path) -> Path:
    pocket_name = pocket_pdb.name
    candidates: list[Path] = []

    if pocket_name.endswith("_rec_pocket.pdb"):
        candidates.append(pocket_pdb.with_name(pocket_name.replace("_rec_pocket.pdb", "_rec.pdb")))
        candidates.append(
            pocket_pdb.with_name(pocket_name.replace("_rec_pocket.pdb", "_rec_nofix.pdb"))
        )
    elif pocket_name.endswith("_rec_pocket.mmtf"):
        candidates.append(
            pocket_pdb.with_name(pocket_name.replace("_rec_pocket.mmtf", "_rec.pdb"))
        )
        candidates.append(
            pocket_pdb.with_name(pocket_name.replace("_rec_pocket.mmtf", "_rec_nofix.pdb"))
        )
    elif pocket_name.endswith("_pocket.pdb"):
        candidates.append(pocket_pdb.with_name(pocket_name.replace("_pocket.pdb", ".pdb")))

    if not candidates:
        return pocket_pdb
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def _normalize_pocket_values(df: pd.DataFrame) -> list[str]:
    if "pocket" not in df.columns:
        return []
    values: list[str] = []
    for raw in df["pocket"].tolist():
        text = str(raw).strip()
        if not text:
            continue
        if text.lower() in {"none", "null", "nan"}:
            continue
        values.append(text)
    return sorted(set(values))


def _resolve_pocket_path_for_value(
    df: pd.DataFrame,
    bigbind_root: Path,
    raw_pocket: str,
) -> Path:
    candidate = resolve_pocket_path(bigbind_root, raw_pocket)
    if candidate.exists() and candidate.is_file():
        return candidate.resolve()

    if "ex_rec_pocket_file" in df.columns:
        matches = df.loc[df["pocket"].astype(str) == raw_pocket, "ex_rec_pocket_file"]
        for raw_alt in matches.tolist():
            alt = str(raw_alt).strip()
            if not alt or alt.lower() in {"none", "null", "nan"}:
                continue
            alt_path = resolve_pocket_path(bigbind_root, alt)
            if alt_path.exists():
                return alt_path.resolve()

    if "ex_rec_pdb" in df.columns:
        pdb_matches = df.loc[df["pocket"].astype(str) == raw_pocket, "ex_rec_pdb"]
        for pdb_value in pdb_matches.tolist():
            pdb_id = str(pdb_value).strip()
            if not pdb_id or pdb_id.lower() in {"none", "null", "nan"}:
                continue
            implied = resolve_pocket_path(
                bigbind_root,
                f"{raw_pocket}/{pdb_id}_rec_pocket.pdb",
            )
            if implied.exists():
                return implied.resolve()

    if candidate.exists() and candidate.is_dir():
        for pattern in ("*_rec_pocket.pdb", "*_pocket.pdb", "*.pdb"):
            matches = sorted(p for p in candidate.glob(pattern) if p.is_file())
            if matches:
                return matches[0].resolve()

    return candidate.resolve()


def _ensure_fpocket_output_root(atlas_cfg: dict[str, Any]) -> Path:
    if "FPOCKET_OUTPUT_ROOT" not in atlas_cfg or not str(atlas_cfg.get("FPOCKET_OUTPUT_ROOT")):
        repo_root = Path(__file__).resolve().parents[1]
        atlas_cfg["FPOCKET_OUTPUT_ROOT"] = str(repo_root / "fpocket")
    return Path(str(atlas_cfg["FPOCKET_OUTPUT_ROOT"]))


def precompute_fpocket_for_bigbind_df(
    df: pd.DataFrame,
    bigbind_root: Path,
    atlas_cfg: dict[str, Any],
    run_dir: Path,
    max_unique_pockets: int | None = None,
    max_workers: int | None = None,
) -> pd.DataFrame:
    logger = logging.getLogger("ml.fpocket_bigbind")
    cfg = atlas_cfg
    output_root = _ensure_fpocket_output_root(cfg)
    output_root.mkdir(parents=True, exist_ok=True)

    unique_pockets = _normalize_pocket_values(df)
    if max_unique_pockets is not None and max_unique_pockets > 0:
        unique_pockets = unique_pockets[: int(max_unique_pockets)]

    receptor_locks: dict[str, threading.Lock] = {}
    receptor_locks_guard = threading.Lock()

    def _lock_for_receptor(receptor_pdb: Path) -> threading.Lock:
        key = str(receptor_pdb.resolve())
        with receptor_locks_guard:
            existing = receptor_locks.get(key)
            if existing is not None:
                return existing
            created = threading.Lock()
            receptor_locks[key] = created
            return created

    def _row_for_pocket(raw_pocket: str) -> dict[str, Any]:
        pocket_abs_path = _resolve_pocket_path_for_value(df, bigbind_root, raw_pocket)
        if pocket_abs_path.is_dir():
            logger.info(
                "[ml.fpocket_bigbind.skip] pocket=%s path=%s reason=resolved_to_directory",
                raw_pocket,
                pocket_abs_path,
            )
            dir_row = {
                "pocket": raw_pocket,
                "pocket_abs_path": str(pocket_abs_path),
                "receptor_pdb": str(pocket_abs_path),
                "pocket_center_x": float("nan"),
                "pocket_center_y": float("nan"),
                "pocket_center_z": float("nan"),
            }
            for col in FPOCKET_METRIC_COLUMNS:
                dir_row[col] = float("nan")
            for col in POCKET_FEATURE_COLUMNS:
                dir_row[col] = 0.0
            return dir_row

        pocket_feature_values = compute_pocket_features(pocket_abs_path)
        receptor_pdb = derive_receptor_pdb(pocket_abs_path)
        center = pocket_center_from_pdb(pocket_abs_path)
        pocket_residues = _normalize_pocket_residues_for_receptor(
            parse_pocket_residues(pocket_abs_path),
            receptor_pdb,
        )

        stem = receptor_pdb.stem
        out_dir = output_root / f"{stem}_out"
        info_path = out_dir / f"{stem}_info.txt"
        sentinel_path = out_dir / ".no_pocket"
        force_retry_empty = _coerce_bool(cfg.get("FPOCKET_FORCE_RETRY_EMPTY", False))

        def _has_nonempty_info() -> bool:
            return info_path.exists() and info_path.stat().st_size > 0

        def _retry_fpocket_once(reason: str) -> bool:
            if not receptor_pdb.exists() or not pocket_residues:
                logger.info(
                    "[ml.fpocket_bigbind.skip] pocket=%s receptor=%s reason=missing_receptor_or_residues",
                    raw_pocket,
                    receptor_pdb,
                )
                return False
            logger.info(
                "[ml.fpocket_bigbind.retry] pocket=%s receptor=%s reason=%s",
                raw_pocket,
                receptor_pdb,
                reason,
            )
            run_fpocket_for_explicit_pocket(
                cfg=cfg,
                receptor_pdb=receptor_pdb,
                pocket_residues=pocket_residues,
                logger=logger,
            )
            return _has_nonempty_info()

        if not _has_nonempty_info():
            with _lock_for_receptor(receptor_pdb):
                has_info = _has_nonempty_info()
                has_empty_info = info_path.exists() and info_path.stat().st_size <= 0
                has_sentinel = sentinel_path.exists()
                attempted_retry = False
                retry_succeeded = False

                if has_info:
                    pass
                elif has_sentinel and not force_retry_empty:
                    logger.debug(
                        "[ml.fpocket_bigbind.cache] pocket=%s receptor=%s action=skip_sentinel path=%s",
                        raw_pocket,
                        receptor_pdb,
                        sentinel_path,
                    )
                else:
                    if has_sentinel and force_retry_empty:
                        sentinel_path.unlink(missing_ok=True)
                    if has_empty_info:
                        logger.info(
                            "[ml.fpocket_bigbind.cache] pocket=%s receptor=%s action=retry_empty_info path=%s",
                            raw_pocket,
                            receptor_pdb,
                            info_path,
                        )
                        info_path.unlink(missing_ok=True)
                        attempted_retry = True
                        retry_succeeded = _retry_fpocket_once("empty_info")
                    elif not has_info:
                        attempted_retry = True
                        retry_succeeded = _retry_fpocket_once("missing_info")

                if attempted_retry and not retry_succeeded:
                    _write_no_pocket_sentinel(
                        sentinel_path,
                        receptor_pdb=receptor_pdb,
                        reason="empty_or_no_output_after_single_retry",
                    )
                if retry_succeeded and sentinel_path.exists():
                    sentinel_path.unlink(missing_ok=True)

        metrics = None
        if center is not None:
            metrics = load_fpocket_metrics_for_receptor_pdb(
                cfg=cfg,
                receptor_pdb=receptor_pdb,
                center=center,
                logger=logger,
            )

        row: dict[str, Any] = {
            "pocket": raw_pocket,
            "pocket_abs_path": str(pocket_abs_path),
            "receptor_pdb": str(receptor_pdb),
            "pocket_center_x": center[0] if center is not None else float("nan"),
            "pocket_center_y": center[1] if center is not None else float("nan"),
            "pocket_center_z": center[2] if center is not None else float("nan"),
        }
        for col in FPOCKET_METRIC_COLUMNS:
            row[col] = (
                float(metrics[col])
                if metrics is not None and col in metrics
                else float("nan")
            )
        for col in POCKET_FEATURE_COLUMNS:
            row[col] = float(pocket_feature_values.get(col, 0.0))
        return row

    workers = max(1, _resolve_worker_count(cfg, max_workers=max_workers))
    workers = min(workers, max(1, len(unique_pockets)))
    rows: list[dict[str, Any]] = []
    if workers <= 1 or len(unique_pockets) <= 1:
        for raw_pocket in unique_pockets:
            rows.append(_row_for_pocket(raw_pocket))
    else:
        logger.info(
            "[ml.fpocket_bigbind.parallel] workers=%d unique_pockets=%d",
            workers,
            len(unique_pockets),
        )
        ordered_rows: list[dict[str, Any] | None] = [None] * len(unique_pockets)
        with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="ml-fpocket") as pool:
            future_to_index = {
                pool.submit(_row_for_pocket, raw_pocket): idx
                for idx, raw_pocket in enumerate(unique_pockets)
            }
            for future in as_completed(future_to_index):
                idx = future_to_index[future]
                ordered_rows[idx] = future.result()
        rows = [row for row in ordered_rows if row is not None]

    columns = [
        "pocket",
        "pocket_abs_path",
        "receptor_pdb",
        *FPOCKET_CENTER_COLUMNS,
        *FPOCKET_METRIC_COLUMNS,
        *POCKET_FEATURE_COLUMNS,
    ]
    metrics_df = pd.DataFrame(rows, columns=columns)
    run_dir.mkdir(parents=True, exist_ok=True)
    metrics_df.to_csv(run_dir / "fpocket_metrics.csv", index=False)
    return metrics_df


def merge_fpocket_metrics_on_pocket(df: pd.DataFrame, metrics_df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    if "pocket" not in df.columns or metrics_df.empty or "pocket" not in metrics_df.columns:
        return df.copy()

    result = df.copy()
    dedup = metrics_df.drop_duplicates(subset=["pocket"]).set_index("pocket")
    for col in (
        "pocket_abs_path",
        "receptor_pdb",
        *FPOCKET_CENTER_COLUMNS,
        *FPOCKET_METRIC_COLUMNS,
        *POCKET_FEATURE_COLUMNS,
    ):
        if col not in dedup.columns:
            continue
        mapped = result["pocket"].map(dedup[col])
        if col in result.columns:
            result[col] = mapped.where(mapped.notna(), result[col])
        else:
            result[col] = mapped
    return result
