from __future__ import annotations

import csv
import math
import re
from pathlib import Path
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Set, Tuple

from post_docking.rescoring.scorch_types import StageSpec

_RUN_MODE_STAGE_SUFFIX_RE = re.compile(
    r"_(?:dud|decoy|decoys|fda|prod|production)_(?:stage|pose)\d+$",
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
    for field in ("Ligand_ID", "ligand", "ligand_file", "ligand_id"):
        key = _ligand_key(row.get(field))
        if key:
            return key
    return ""


def aggregate_combo(
    post_root: Path,
    specs: Sequence[StageSpec],
    combo: Tuple[str, str, str],
    logger: Any,
    *,
    component: str,
    run_mode: str = "fda",
    output_name: str = "scorch_scores_all.csv",
    decoy_prefix: str = "dud",
    input_csvs: Optional[Sequence[Path | str]] = None,
    allowed_ligand_bases: Optional[Set[str] | Sequence[str]] = None,
) -> Optional[Path]:
    def _safe_float(value: object) -> Optional[float]:
        try:
            parsed = float(str(value).strip())
        except Exception:
            return None
        if not math.isfinite(parsed):
            return None
        return parsed

    def _dedup_key(row: Mapping[str, str]) -> Optional[Tuple[str, str]]:
        lig_raw = str(row.get("Ligand_ID", "") or row.get("ligand", "")).strip()
        if not lig_raw:
            return None
        source = str(row.get("source", "")).strip().lower()
        lig_norm = Path(lig_raw).name
        if not lig_norm:
            return None
        return source, lig_norm

    def _dedup_rank(row: Mapping[str, str]) -> float:
        score = _safe_float(row.get("SCORCH_score") or row.get("SCORCH_score_used"))
        cert = _safe_float(
            row.get("SCORCH_certainty") or row.get("SCORCH_certainty_used")
        )
        if score is None:
            return float("-inf")
        if cert is None:
            return score
        return score * cert

    pdb_id, variant, ph = combo
    combo_dir = post_root / pdb_id / variant / ph
    rows: List[Dict[str, str]] = []
    fields: List[str] = []

    paths: list[Path] = []
    seen_paths: set[Path] = set()
    if input_csvs is not None:
        for raw_path in input_csvs:
            p = Path(raw_path)
            if not p.exists() or p.stat().st_size <= 0:
                continue
            try:
                resolved = p.resolve()
            except Exception:
                resolved = p
            if resolved in seen_paths:
                continue
            seen_paths.add(resolved)
            paths.append(p)
        logger.info(
            "%s action=aggregate_input status=plan_bound pdb_id=%s variant=%s ph=%s run_mode=%s input_csvs=%d",
            component,
            pdb_id,
            variant,
            ph,
            run_mode,
            len(paths),
        )
    else:
        for spec in specs:
            base_name = (
                spec.output_name if run_mode == "fda" else f"{decoy_prefix}_{spec.output_name}"
            )
            base_path = combo_dir / base_name
            if base_path.exists() and base_path.stat().st_size > 0:
                try:
                    resolved = base_path.resolve()
                except Exception:
                    resolved = base_path
                seen_paths.add(resolved)
                paths.append(base_path)
            shard_globs = [
                base_name.replace(".csv", ".part*.csv"),
                base_name.replace(".csv", ".tailt*.csv"),
                base_name.replace(".csv", ".rt*.csv"),
            ]
            for shard_glob in shard_globs:
                for p in sorted(combo_dir.glob(shard_glob)):
                    if not p.exists() or p.stat().st_size <= 0:
                        continue
                    try:
                        resolved = p.resolve()
                    except Exception:
                        resolved = p
                    if resolved in seen_paths:
                        continue
                    seen_paths.add(resolved)
                    paths.append(p)

    for csv_path in paths:
        try:
            with csv_path.open() as handle:
                reader = csv.DictReader(handle)
                if not reader.fieldnames:
                    continue
                for row in reader:
                    rows.append(row)
                for fn in reader.fieldnames:
                    if fn not in fields:
                        fields.append(fn)
        except Exception as exc:
            logger.error(
                "%s action=aggregate status=failed path=%s reason=read_error error=%s",
                component,
                csv_path,
                exc,
            )
            continue

    filter_to_allowed = allowed_ligand_bases is not None
    allowed_keys = {
        _ligand_key(base) for base in (allowed_ligand_bases or set()) if _ligand_key(base)
    }
    filtered_rows = 0
    unkeyed_filtered_rows = 0
    if filter_to_allowed:
        kept_rows: List[Dict[str, str]] = []
        for row in rows:
            row_key = _row_ligand_key(row)
            if row_key and row_key in allowed_keys:
                kept_rows.append(row)
                continue
            filtered_rows += 1
            if not row_key:
                unkeyed_filtered_rows += 1
        if filtered_rows:
            logger.info(
                "%s action=aggregate_filter status=ok pdb_id=%s variant=%s ph=%s run_mode=%s allowed=%d filtered=%d unkeyed_filtered=%d",
                component,
                pdb_id,
                variant,
                ph,
                run_mode,
                len(allowed_keys),
                filtered_rows,
                unkeyed_filtered_rows,
            )
        rows = kept_rows

    if not rows:
        logger.warning(
            "%s action=aggregate status=skip reason=no_rows pdb_id=%s variant=%s ph=%s",
            component,
            pdb_id,
            variant,
            ph,
        )
        return None

    deduped_rows: List[Dict[str, str]] = []
    deduped_rank: List[float] = []
    deduped_index: Dict[Tuple[str, str], int] = {}
    duplicates_removed = 0
    for row in rows:
        key = _dedup_key(row)
        if key is None:
            deduped_rows.append(row)
            deduped_rank.append(float("-inf"))
            continue
        new_rank = _dedup_rank(row)
        existing_idx = deduped_index.get(key)
        if existing_idx is None:
            deduped_index[key] = len(deduped_rows)
            deduped_rows.append(row)
            deduped_rank.append(new_rank)
            continue
        duplicates_removed += 1
        if new_rank > deduped_rank[existing_idx]:
            deduped_rows[existing_idx] = row
            deduped_rank[existing_idx] = new_rank
    rows = deduped_rows

    out_path = combo_dir / output_name
    try:
        with out_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
    except Exception as exc:
        logger.error(
            "%s action=aggregate status=failed path=%s reason=write_error error=%s",
            component,
            out_path,
            exc,
        )
        return None

    logger.info(
        "%s action=aggregate status=ok rows=%d output=%s pdb_id=%s variant=%s ph=%s run_mode=%s",
        component,
        len(rows),
        out_path,
        pdb_id,
        variant,
        ph,
        run_mode,
    )
    if duplicates_removed > 0:
        logger.info(
            "%s action=aggregate_dedup status=ok pdb_id=%s variant=%s ph=%s run_mode=%s duplicates_removed=%d rows_out=%d",
            component,
            pdb_id,
            variant,
            ph,
            run_mode,
            int(duplicates_removed),
            len(rows),
        )
    return out_path


def compute_best_composites(
    rows: List[Dict[str, str]],
    *,
    pose_base_from_path: Callable[[Path], str],
) -> Dict[str, float]:
    def _as_float_local(val: Any) -> Optional[float]:
        try:
            if val is None:
                return None
            s = str(val).strip()
            if not s:
                return None
            out = float(s)
            return out
        except Exception:
            return None

    best: Dict[str, float] = {}
    for r in rows:
        lig_id = str(r.get("Ligand_ID", "") or r.get("ligand", "")).strip()
        if not lig_id:
            continue
        base = pose_base_from_path(Path(lig_id))
        if not base:
            continue
        score_val = _as_float_local(r.get("SCORCH_score") or r.get("SCORCH_score_used"))
        cert_val = _as_float_local(
            r.get("SCORCH_certainty") or r.get("SCORCH_certainty_used")
        )
        if score_val is None or not math.isfinite(score_val):
            continue
        comp = (
            score_val
            if cert_val is None or not math.isfinite(cert_val)
            else score_val * cert_val
        )
        if not math.isfinite(comp):
            continue
        if base not in best or comp > best[base]:
            best[base] = comp
    return best


def annotate_scorch_file(
    csv_path: Path,
    best_map: Dict[str, float],
    mu_decoy: Optional[float],
    sigma_decoy: Optional[float],
    n_decoys: int,
    *,
    pose_base_from_path: Callable[[Path], str],
) -> None:
    rows: List[Dict[str, str]] = []
    with csv_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        rows = [dict(r) for r in reader]

    extra_fields = [
        "scorch_best_composite",
        "scorch_mu_decoy",
        "scorch_sigma_decoy",
        "scorch_n_decoys",
        "z_vs_decoys_scorch",
        "t_vs_decoys_scorch",
    ]
    for ef in extra_fields:
        if ef not in fields:
            fields.append(ef)

    for r in rows:
        lig_id = str(r.get("Ligand_ID", "") or r.get("ligand", "")).strip()
        base = pose_base_from_path(Path(lig_id)) if lig_id else ""
        comp = best_map.get(base) if base else None
        r["scorch_best_composite"] = "" if comp is None else f"{comp:.6g}"
        r["scorch_mu_decoy"] = "" if mu_decoy is None else f"{mu_decoy:.6g}"
        r["scorch_sigma_decoy"] = "" if sigma_decoy is None else f"{sigma_decoy:.6g}"
        r["scorch_n_decoys"] = "" if n_decoys <= 0 else str(n_decoys)
        if comp is None or mu_decoy is None or sigma_decoy is None or sigma_decoy <= 0:
            r["z_vs_decoys_scorch"] = ""
            r["t_vs_decoys_scorch"] = ""
        else:
            z_val = f"{(comp - mu_decoy) / sigma_decoy:.6g}"
            r["z_vs_decoys_scorch"] = z_val
            r["t_vs_decoys_scorch"] = z_val

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def annotate_scorch_z_scores(
    fda_csv: Path,
    dud_csv: Path,
    logger: Any,
    *,
    pose_base_from_path: Callable[[Path], str],
) -> None:
    if not fda_csv.exists() or not dud_csv.exists():
        logger.info(
            "[z-score.schorch.skip] reason=missing_inputs fda_exists=%s dud_exists=%s",
            fda_csv.exists(),
            dud_csv.exists(),
        )
        return

    def _read(path: Path) -> List[Dict[str, str]]:
        with path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            return [dict(r) for r in reader]

    fda_rows = _read(fda_csv)
    dud_rows = _read(dud_csv)
    best_fda = compute_best_composites(
        fda_rows, pose_base_from_path=pose_base_from_path
    )
    best_dud = compute_best_composites(
        dud_rows, pose_base_from_path=pose_base_from_path
    )
    decoy_values = [v for v in best_dud.values() if v is not None and math.isfinite(v)]
    if not decoy_values:
        logger.info(
            "[z-score.schorch.skip] reason=no_decoys n_decoys=0 fda=%s dud=%s",
            fda_csv,
            dud_csv,
        )
        mu = sigma = None
        n_decoys = 0
    else:
        mu = sum(decoy_values) / len(decoy_values)
        variance = sum((v - mu) ** 2 for v in decoy_values) / len(decoy_values)
        sigma = math.sqrt(variance)
        n_decoys = len(decoy_values)
        if sigma <= 0 or not math.isfinite(mu) or not math.isfinite(sigma):
            logger.info(
                "[z-score.schorch.skip] reason=no_decoys_or_sigma0 n_decoys=%d mu=%s sigma=%s",
                n_decoys,
                mu,
                sigma,
            )
            mu = sigma = None
        else:
            logger.info(
                "[z-score.schorch] n_decoys=%d mu=%.6g sigma=%.6g fda=%s dud=%s",
                n_decoys,
                mu,
                sigma,
                fda_csv,
                dud_csv,
            )

    annotate_mu = mu if mu is not None else None
    annotate_sigma = sigma if sigma is not None else None
    annotate_n = n_decoys if mu is not None and sigma is not None else 0
    annotate_scorch_file(
        fda_csv,
        best_fda,
        annotate_mu,
        annotate_sigma,
        annotate_n,
        pose_base_from_path=pose_base_from_path,
    )
    annotate_scorch_file(
        dud_csv,
        best_dud,
        annotate_mu,
        annotate_sigma,
        annotate_n,
        pose_base_from_path=pose_base_from_path,
    )


def annotate_scorch_t_scores(
    fda_csv: Path,
    dud_csv: Path,
    logger: Any,
    *,
    pose_base_from_path: Callable[[Path], str],
) -> None:
    annotate_scorch_z_scores(
        fda_csv,
        dud_csv,
        logger,
        pose_base_from_path=pose_base_from_path,
    )
