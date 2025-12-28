# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

COMPONENT = "[rescore-reranker]"
NUMERIC_PRIORITY = [
    "consensus_score",
    "final_score",
    "scorch_composite",
    "SCORCH_score_used",
    "SCORCH_certainty_used",
    "rescored_flag",
    "consensus_rank",
    "final_rank",
    "p_energy",
    "p_vina",
    "p_gnina_energy",
    "p_ledock",
    "p_dock6",
    "p_cnn",
    "n_engines_with_data",
]
TEXT_PRIORITY = [
    "best_engine",
    "best_signal",
    "scorch_source_used",
]


def _as_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        s = str(x).strip()
        if not s:
            return None
        return float(s)
    except Exception:
        return None


def _norm_engine(x: Any) -> str:
    return str(x or "").strip().lower()


def _extract_stage_num(s: str) -> int:
    m = re.search(r"stage(\d+)", s)
    if not m:
        return -1
    try:
        return int(m.group(1))
    except Exception:
        return -1


def _consensus_ligand_base(lig: str) -> str:
    s = str(lig or "").strip()
    s = re.sub(r"\.(pdbqt|mol2|sdf)$", "", s, flags=re.IGNORECASE)
    s = s.replace(".sanitized", "")
    return s


def _scorch_ligand_base(lig_id: str) -> str:
    s = str(lig_id or "").strip()
    s = re.sub(r"(__ledock_stage\d+)$", "", s)
    s = re.sub(r"(__dock6_stage\d+)$", "", s)
    s = re.sub(r"(_gnina_stage\d+)$", "", s)
    s = re.sub(r"(_stage\d+)$", "", s)
    s = re.sub(r"\.(mol2|pdbqt)$", "", s, flags=re.IGNORECASE)
    s = s.replace(".sanitized", "")
    s = s.replace("__", "_")
    s = re.sub(r"_+$", "", s)
    return s


def find_consensus_csv(docked_combo_dir: Path) -> Optional[Path]:
    candidates = [
        docked_combo_dir / "consensus_docking_scores.csv",
        docked_combo_dir / "consensus.csv",
    ]
    for p in candidates:
        if p.exists() and p.stat().st_size > 0:
            return p
    for p in sorted(docked_combo_dir.glob("consensus*.csv")):
        if p.exists() and p.stat().st_size > 0:
            return p
    return None


@dataclass(frozen=True)
class ScorchPick:
    pdb_id: str
    variant: str
    ph: str
    source: str
    ligand_base: str
    stage_num: int
    scorch_score: Optional[float]
    scorch_certainty: Optional[float]


def _read_csv(path: Path) -> Tuple[List[Dict[str, str]], List[str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        rows = [dict(r) for r in reader]
        fields = list(reader.fieldnames or [])
    return rows, fields


def _write_csv(
    path: Path, rows: List[Dict[str, Any]], fieldnames: List[str], preamble_lines: Optional[List[str]] = None
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        if preamble_lines:
            for line in preamble_lines:
                f.write(line)
                if not line.endswith("\n"):
                    f.write("\n")
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in fieldnames})


def _is_numeric_field(name: str, rows: List[Dict[str, Any]]) -> bool:
    lowered = name.lower()
    if lowered in {n.lower() for n in NUMERIC_PRIORITY}:
        return True
    values = [r.get(name, "") for r in rows]
    non_empty = [v for v in values if str(v).strip() != ""]
    if not non_empty:
        return False
    success = 0
    for v in non_empty:
        try:
            float(str(v).strip())
            success += 1
        except Exception:
            continue
    return success >= 0.8 * len(non_empty)


def _ordered_fields(original_fields: List[str]) -> List[str]:
    # Exclude internal fields
    base_fields = [f for f in original_fields if f not in {"ligand_base", "scorch_stage_used"}]
    # Always include added fields if not present
    for extra in [
        "scorch_source_used",
        "SCORCH_score_used",
        "SCORCH_certainty_used",
        "scorch_composite",
        "final_score",
        "rescored_flag",
        "consensus_rank",
        "final_rank",
    ]:
        if extra not in base_fields:
            base_fields.append(extra)

    # Build numeric/qualitative sets
    dummy_rows: List[Dict[str, Any]] = []
    numeric_fields = []
    qualitative_fields = []
    for f in base_fields:
        if _is_numeric_field(f, dummy_rows):
            numeric_fields.append(f)
        else:
            qualitative_fields.append(f)

    # Respect ligand first
    ordered: List[str] = []
    if "ligand" in base_fields:
        ordered.append("ligand")
    # Numeric priority
    for f in NUMERIC_PRIORITY:
        if f in base_fields and f not in ordered:
            ordered.append(f)
    for f in base_fields:
        if f in numeric_fields and f not in ordered:
            ordered.append(f)
    # Qualitative priority
    for f in TEXT_PRIORITY:
        if f in base_fields and f not in ordered:
            ordered.append(f)
    for f in base_fields:
        if f not in ordered:
            ordered.append(f)
    return ordered


def _sort_rows(rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    def _rank_key(r: Dict[str, Any]) -> Tuple[int, float, float, float]:
        try:
            fr = int(str(r.get("final_rank", "")).strip())
            return (fr, 0.0, 0.0, 0.0)
        except Exception:
            comp = _as_float(r.get("scorch_composite")) or -1e18
            scs = _as_float(r.get("SCORCH_score_used")) or -1e18
            cs = _as_float(r.get("consensus_score")) or -1e18
            # large negative final_rank signals fallback path; ensure fr high to deprioritize
            return (int(1e9), -comp, -scs, -cs)

    return sorted(rows, key=_rank_key)


def rerank_consensus_with_scorch(
    consensus_csv: Path,
    scorch_csv: Path,
    out_csv: Path,
    logger: logging.Logger,
    *,
    overwrite: bool = False,
    ) -> bool:
    if out_csv.exists() and out_csv.stat().st_size > 0 and not overwrite:
        logger.info("%s action=skip reason=exists out=%s", COMPONENT, str(out_csv))
        return True

    if not consensus_csv.exists() or consensus_csv.stat().st_size == 0:
        logger.warning("%s action=skip reason=missing_consensus path=%s", COMPONENT, str(consensus_csv))
        return False
    degraded = False
    sc_rows: List[Dict[str, str]] = []
    if not scorch_csv.exists() or scorch_csv.stat().st_size == 0:
        degraded = True
        logger.warning("%s action=rerank status=degraded reason=missing_or_empty_scorch path=%s", COMPONENT, str(scorch_csv))
    else:
        try:
            sc_rows, _ = _read_csv(scorch_csv)
            if not sc_rows:
                degraded = True
                logger.warning(
                    "%s action=rerank status=degraded reason=missing_or_empty_scorch path=%s",
                    COMPONENT,
                    str(scorch_csv),
                )
        except Exception as exc:
            degraded = True
            logger.warning(
                "%s action=rerank status=degraded reason=read_error_scorch path=%s error=%s",
                COMPONENT,
                str(scorch_csv),
                exc,
            )
            sc_rows = []

    cons_rows, cons_fields = _read_csv(consensus_csv)

    if not cons_rows:
        logger.warning("%s action=skip reason=empty_consensus path=%s", COMPONENT, str(consensus_csv))
        return False

    best_by_source: Dict[Tuple[str, str, str, str, str], ScorchPick] = {}
    best_any: Dict[Tuple[str, str, str, str], ScorchPick] = {}

    for r in sc_rows:
        pdb_id = str(r.get("pdb_id", "")).strip()
        variant = str(r.get("variant", "")).strip()
        ph = str(r.get("ph", "")).strip()
        source = _norm_engine(r.get("source", ""))
        lig_id = str(r.get("Ligand_ID", "")).strip()
        if not lig_id:
            continue

        lig_base = _scorch_ligand_base(lig_id)
        stage_num = _extract_stage_num(lig_id)
        score = _as_float(r.get("SCORCH_score"))
        certainty = _as_float(r.get("SCORCH_certainty"))

        pick = ScorchPick(
            pdb_id=pdb_id,
            variant=variant,
            ph=ph,
            source=source,
            ligand_base=lig_base,
            stage_num=stage_num,
            scorch_score=score,
            scorch_certainty=certainty,
        )

        k1 = (pdb_id, variant, ph, source, lig_base)
        cur = best_by_source.get(k1)
        if cur is None or (pick.stage_num, pick.scorch_score or -1e9) > (cur.stage_num, cur.scorch_score or -1e9):
            best_by_source[k1] = pick

        k2 = (pdb_id, variant, ph, lig_base)
        cur2 = best_any.get(k2)
        if cur2 is None or (pick.stage_num, pick.scorch_score or -1e9) > (cur2.stage_num, cur2.scorch_score or -1e9):
            best_any[k2] = pick

    enriched: List[Dict[str, Any]] = []
    for r in cons_rows:
        rr: Dict[str, Any] = dict(r)
        lig = str(r.get("ligand", "")).strip()
        lig_base = _consensus_ligand_base(lig)
        rr["ligand_base"] = lig_base

        pdb_id = str(r.get("pdb_id", "")).strip()
        variant = str(r.get("variant", "")).strip()
        ph_label = str(r.get("ph_label", "")).strip()
        best_engine = _norm_engine(r.get("best_engine", ""))

        pick = None
        if pdb_id and variant and ph_label and best_engine:
            pick = best_by_source.get((pdb_id, variant, ph_label, best_engine, lig_base))
        if pick is None and pdb_id and variant and ph_label:
            pick = best_any.get((pdb_id, variant, ph_label, lig_base))

        rr["scorch_source_used"] = pick.source if pick else ""
        rr["SCORCH_score_used"] = "" if (not pick or pick.scorch_score is None) else f"{pick.scorch_score:.6g}"
        rr["SCORCH_certainty_used"] = "" if (not pick or pick.scorch_certainty is None) else f"{pick.scorch_certainty:.6g}"

        comp = None
        if pick and pick.scorch_score is not None and pick.scorch_certainty is not None:
            comp = pick.scorch_score * pick.scorch_certainty
        rr["scorch_composite"] = "" if comp is None else f"{comp:.6g}"

        rescored = (rr["scorch_composite"] != "") or (rr["SCORCH_score_used"] != "")
        rr["rescored_flag"] = "1" if rescored else "0"
        if rescored and rr["scorch_composite"] != "":
            rr["final_score"] = rr["scorch_composite"]
        elif rescored and rr["SCORCH_score_used"] != "":
            rr["final_score"] = rr["SCORCH_score_used"]
        else:
            rr["final_score"] = str(rr.get("consensus_score", "")).strip()

        enriched.append(rr)

    def group_key(row: Dict[str, Any]) -> Tuple[str, str, str, str]:
        return (
            str(row.get("run_id", "")).strip(),
            str(row.get("pdb_id", "")).strip(),
            str(row.get("variant", "")).strip(),
            str(row.get("ph_label", "")).strip(),
        )

    groups: Dict[Tuple[str, str, str, str], List[int]] = {}
    for i, row in enumerate(enriched):
        groups.setdefault(group_key(row), []).append(i)

    for _, idxs in groups.items():
        scored = []
        for i in idxs:
            cs = _as_float(enriched[i].get("consensus_score"))
            scored.append((cs if cs is not None else -1e18, i))
        scored.sort(key=lambda x: x[0], reverse=True)
        for rank, (_, i) in enumerate(scored, start=1):
            enriched[i]["consensus_rank"] = str(rank)

    for _, idxs in groups.items():
        rescored_idxs: List[int] = []
        nonrescored_idxs: List[int] = []
        for i in idxs:
            comp = enriched[i].get("scorch_composite")
            scs = enriched[i].get("SCORCH_score_used")
            if str(comp).strip() != "" or str(scs).strip() != "":
                rescored_idxs.append(i)
            else:
                nonrescored_idxs.append(i)

        def _rescored_key(i: int) -> Tuple[float, float, float]:
            comp = _as_float(enriched[i].get("scorch_composite")) or -1e18
            scs = _as_float(enriched[i].get("SCORCH_score_used")) or -1e18
            cs = _as_float(enriched[i].get("consensus_score")) or -1e18
            return (comp, scs, cs)

        rescored_sorted = sorted(rescored_idxs, key=lambda idx: _rescored_key(idx), reverse=True)
        nonrescored_sorted = sorted(
            nonrescored_idxs,
            key=lambda idx: _as_float(enriched[idx].get("consensus_score")) or -1e18,
            reverse=True,
        )

        final_order = rescored_sorted + nonrescored_sorted
        for rank, i in enumerate(final_order, start=1):
            enriched[i]["final_rank"] = str(rank)

    for row in enriched:
        rescored = bool(str(row.get("scorch_composite", "")).strip() or str(row.get("SCORCH_score_used", "")).strip())
        row["rescored_flag"] = "1" if rescored else "0"
        fr = _as_float(row.get("final_rank"))
        row["final_score"] = "" if fr is None else f"{(-fr):.6g}"

    out_fields = _ordered_fields(cons_fields)
    pretty_fields = [f for f in out_fields if f not in {"run_id", "pdb_id", "variant", "ph_label"}]
    sorted_rows = _sort_rows(enriched)

    _write_csv(out_csv, sorted_rows, out_fields)

    meta = [
        f"# run_id={str(cons_rows[0].get('run_id', '')).strip()}\n",
        "# pdb_id={pdb}, variant={variant}, ph={ph}\n".format(
            pdb=str(cons_rows[0].get("pdb_id", "")).strip(),
            variant=str(cons_rows[0].get("variant", "")).strip(),
            ph=str(cons_rows[0].get("ph_label", "")).strip(),
        ),
    ]
    _write_csv(out_csv.with_suffix(".pretty.csv"), sorted_rows, pretty_fields, preamble_lines=meta)

    logger.info(
        "%s action=write status=ok consensus=%s scorch=%s out=%s rows=%d",
        COMPONENT,
        str(consensus_csv),
        str(scorch_csv),
        str(out_csv),
        len(enriched),
    )
    return True


def rerank_run(run_id: str, repo_root: Path, overwrite: bool, logger: logging.Logger) -> int:
    post_root = repo_root / "post_docked" / run_id
    dock_root = repo_root / "docked" / run_id
    if not post_root.exists():
        logger.warning("%s action=discover status=skip reason=missing_post_root path=%s", COMPONENT, str(post_root))
        return 1
    if not dock_root.exists():
        logger.warning("%s action=discover status=skip reason=missing_dock_root path=%s", COMPONENT, str(dock_root))
        return 1

    scorch_paths = list(post_root.glob("**/scorch_scores_all.csv"))
    if not scorch_paths:
        logger.warning("%s action=discover status=skip reason=no_scorch_scores path=%s", COMPONENT, str(post_root))
        return 1

    ok = 0
    fail = 0
    for scorch_csv in sorted(scorch_paths):
        try:
            combo_dir = scorch_csv.parent
            ph = combo_dir.name
            variant = combo_dir.parent.name
            pdb_id = combo_dir.parent.parent.name

            dock_combo_dir = dock_root / pdb_id / variant / ph
            consensus_csv = find_consensus_csv(dock_combo_dir)
            if not consensus_csv:
                logger.warning(
                    "%s action=skip reason=missing_consensus pdb=%s variant=%s ph=%s dock_dir=%s",
                    COMPONENT,
                    pdb_id,
                    variant,
                    ph,
                    str(dock_combo_dir),
                )
                fail += 1
                continue

            out_csv = combo_dir / "consensus_reranked_scorch.csv"
            if rerank_consensus_with_scorch(consensus_csv, scorch_csv, out_csv, logger, overwrite=overwrite):
                ok += 1
            else:
                fail += 1
        except Exception:
            logger.warning(
                "%s action=rerank status=failed reason=exception path=%s",
                COMPONENT,
                str(scorch_csv),
                exc_info=True,
            )
            fail += 1

    logger.info("%s action=summary status=ok combos_ok=%d combos_failed=%d", COMPONENT, ok, fail)
    return 0 if ok > 0 and fail == 0 else (0 if ok > 0 else 1)


def _configure_logging(verbose: bool) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(message)s")
    return logging.getLogger("rescore-reranker")


def main() -> int:
    ap = argparse.ArgumentParser(description="Rerank consensus docking results using SCORCH rescoring outputs.")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--repo-root", required=True)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    logger = _configure_logging(args.verbose)
    repo_root = Path(args.repo_root).resolve()
    return rerank_run(args.run_id, repo_root, args.overwrite, logger)


if __name__ == "__main__":
    sys.exit(main())
