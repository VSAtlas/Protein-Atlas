from __future__ import annotations

import csv
import hashlib
import logging
import os
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Set

from post_docking.rescoring.scorch_selection import pose_base_from_path
from post_docking.rescoring.scorch_types import AnnotateResult, MaterializeResult


COMPONENT = "[scorch-rescore]"
_DECOY_PREFIX_VALUE = "dud"


def set_decoy_prefix(value: str) -> None:
    global _DECOY_PREFIX_VALUE
    _DECOY_PREFIX_VALUE = value


def _planned_materialized_names(ph_root: Path, ligands: List[Path]) -> Dict[Path, str]:
    existing_names: Set[str] = set()
    planned: Dict[Path, str] = {}
    for src in ligands:
        name = src.name
        if name in existing_names:
            try:
                rel_str = str(src.relative_to(ph_root))
            except ValueError:
                rel_str = str(src)
            digest = hashlib.sha1(rel_str.encode("utf-8")).hexdigest()[:8]
            name = f"{digest}_{name}"
        existing_names.add(name)
        planned[src] = name
    return planned


def materialize_inputs(
    ph_root: Path,
    combo_post_root: Path,
    stage_dir: str,
    ligands: List[Path],
    overwrite: bool,
    logger: logging.Logger,
    *,
    run_mode: str,
    chunk_tag: Optional[str] = None,
) -> MaterializeResult:
    if not ligands:
        return MaterializeResult(
            input_dir=None,
            materialized_count=0,
            failed_count=0,
            failed_examples=(),
        )

    input_dir = combo_post_root / ".scorch_inputs" / run_mode / stage_dir
    if chunk_tag:
        input_dir = input_dir / chunk_tag
    if input_dir.exists() and not overwrite:
        existing = list(input_dir.glob("*.pdbqt"))
        expected_names = set(_planned_materialized_names(ph_root, ligands).values())
        existing_names = {path.name for path in existing}
        if existing and existing_names == expected_names:
            return MaterializeResult(
                input_dir=input_dir,
                materialized_count=len(existing),
                failed_count=0,
                failed_examples=(),
            )
        logger.warning(
            "%s action=inputs status=rebuild stage=%s run_mode=%s chunk=%s reason=stale_materialized_inputs existing=%d expected=%d",
            COMPONENT,
            stage_dir,
            run_mode,
            chunk_tag or "all",
            len(existing_names),
            len(expected_names),
        )
        shutil.rmtree(input_dir)
    elif input_dir.exists():
        shutil.rmtree(input_dir)

    input_dir.mkdir(parents=True, exist_ok=True)

    materialized_count = 0
    failed_sources: List[str] = []
    planned_names = _planned_materialized_names(ph_root, ligands)
    for src, name in planned_names.items():
        dst = input_dir / name
        if not src.exists():
            logger.warning(
                "%s action=link status=degraded stage=%s source=%s reason=missing_source",
                COMPONENT,
                stage_dir,
                src,
            )
            failed_sources.append(str(src))
            continue
        try:
            os.link(src, dst)
            materialized_count += 1
        except OSError:
            try:
                shutil.copy2(src, dst)
                materialized_count += 1
            except Exception as exc:
                try:
                    dst.unlink(missing_ok=True)
                except Exception:
                    pass
                logger.warning(
                    "%s action=link status=degraded stage=%s source=%s reason=copy_error error=%s",
                    COMPONENT,
                    stage_dir,
                    src,
                    exc,
                )
                failed_sources.append(str(src))
                continue

    failed_examples = tuple(failed_sources[:5])
    if materialized_count <= 0:
        return MaterializeResult(
            input_dir=None,
            materialized_count=0,
            failed_count=len(failed_sources),
            failed_examples=failed_examples,
        )
    return MaterializeResult(
        input_dir=input_dir,
        materialized_count=materialized_count,
        failed_count=len(failed_sources),
        failed_examples=failed_examples,
    )





def annotate_csv(
    csv_path: Path,
    metadata: Dict[str, str],
    row_metadata_by_base: Optional[Dict[str, Dict[str, str]]],
    logger: logging.Logger,
) -> AnnotateResult:
    try:
        with csv_path.open() as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames or []
            rows = list(reader)
    except Exception as exc:
        logger.error(
            "%s action=annotate status=failed path=%s reason=read_error error=%s",
            COMPONENT,
            csv_path,
            exc,
        )
        return AnnotateResult(
            ok=False,
            degraded=True,
            reason="annotate_read_error",
        )

    if not rows:
        logger.warning(
            "%s action=annotate status=warning reason=no_rows path=%s",
            COMPONENT,
            csv_path,
        )

    ligand_keys = (
        "ligand",
        "Ligand",
        "Ligand_ID",
        "ligand_id",
        "ligand_file",
        "Ligand_file",
        "file",
        "filename",
        "name",
        "molecule",
        "Molecule",
    )
    out_rows: List[Dict[str, str]] = []
    for row in rows:
        ligand_val = ""
        for key in ligand_keys:
            if row.get(key):
                ligand_val = row[key]
                break
        ligand_val = Path(ligand_val).name if ligand_val else ""
        row.update(metadata)
        lig_base = pose_base_from_path(Path(ligand_val), decoy_prefix=_DECOY_PREFIX_VALUE) if ligand_val else ""
        if lig_base and row_metadata_by_base:
            row.update(row_metadata_by_base.get(lig_base, {}))
        row["ligand_file"] = ligand_val or metadata.get("ligand_file", "")
        out_rows.append(row)

    row_meta_fields: List[str] = []
    if row_metadata_by_base:
        seen_row_meta: Set[str] = set()
        for values in row_metadata_by_base.values():
            for key in values.keys():
                if key not in seen_row_meta:
                    seen_row_meta.add(key)
                    row_meta_fields.append(key)

    meta_fields = list(metadata.keys()) + row_meta_fields + ["ligand_file"]
    out_fields = meta_fields + [fn for fn in fieldnames if fn not in meta_fields]

    try:
        with csv_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=out_fields)
            writer.writeheader()
            for row in out_rows:
                writer.writerow(row)
    except Exception as exc:
        logger.error(
            "%s action=annotate status=failed path=%s reason=write_error error=%s",
            COMPONENT,
            csv_path,
            exc,
        )
        return AnnotateResult(
            ok=False,
            degraded=True,
            reason="annotate_write_error",
        )
    return AnnotateResult(ok=True, degraded=False, reason="ok")

