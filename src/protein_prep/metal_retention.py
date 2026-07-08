"""Metal and cofactor retention helpers for protein preparation."""

from __future__ import annotations

import logging
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from protein_prep.pdb_records import atom_identity, line_element, line_serial, line_xyz
from protein_prep.pdb_fixer_runtime import (
    fix_pdb_elements,
    load_canonical_cofactors,
    load_canonical_metals,
)

TAIL_RECORD_PREFIXES = ("CONECT", "MASTER", "END", "ENDMDL")
WATER_NAMES = {"HOH", "WAT", "H2O", "DOD"}
METAL_BOUND_HET_MAX_A = 3.0


@dataclass(frozen=True)
class RetainedHetRecord:
    kind: str
    line: str
    identity: tuple[str, str, str, str, str, str]
    token: str
    serial: int | None


def _canonical_sets(
    cfg: Mapping[str, object] | None = None,
) -> tuple[set[str], set[str]]:
    try:
        return set(load_canonical_metals(cfg)), set(load_canonical_cofactors(cfg))
    except Exception:
        return set(load_canonical_metals(None)), set(load_canonical_cofactors(None))


_line_element = line_element
_line_identity = atom_identity
_line_serial = line_serial


def _is_retained_het(
    line: str, metals: set[str], cofactors: set[str]
) -> tuple[str, str] | None:
    if not line.startswith("HETATM"):
        return None
    resname = line[17:20].strip().upper()
    element = _line_element(line)
    if resname in metals or element in metals:
        return "metal", element or resname
    if resname in cofactors:
        return "cofactor", resname
    return None


def _collect_retained_records(
    path: Path, metals: set[str], cofactors: set[str]
) -> list[RetainedHetRecord]:
    records: list[RetainedHetRecord] = []
    if not path.exists():
        return records
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            classification = _is_retained_het(line, metals, cofactors)
            if classification is None:
                continue
            kind, token = classification
            records.append(
                RetainedHetRecord(
                    kind=kind,
                    line=line.rstrip("\n"),
                    identity=_line_identity(line),
                    token=token,
                    serial=_line_serial(line),
                )
            )
    return records


def _count_records(records: Iterable[RetainedHetRecord]) -> tuple[int, int]:
    metals = 0
    cofactors = 0
    for record in records:
        if record.kind == "metal":
            metals += 1
        elif record.kind == "cofactor":
            cofactors += 1
    return metals, cofactors


def _existing_serials(lines: Sequence[str]) -> set[int]:
    serials: set[int] = set()
    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        serial = _line_serial(line)
        if serial is not None:
            serials.add(serial)
    return serials


def _rewrite_serial(line: str, serial: int) -> str:
    padded = line.rstrip("\n").ljust(80)
    return f"{padded[:6]}{serial:5d}{padded[11:]}".rstrip()


def _tail_insert_index(lines: Sequence[str]) -> int:
    idx = len(lines)
    while idx > 0:
        prefix = lines[idx - 1][:6].strip().upper()
        if prefix not in TAIL_RECORD_PREFIXES:
            break
        idx -= 1
    return idx


def _insert_records_before_tail(
    target_path: Path, records: Sequence[RetainedHetRecord]
) -> int:
    lines = target_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    existing = {
        _line_identity(line)
        for line in lines
        if line.startswith(("ATOM  ", "HETATM"))
    }
    serials = _existing_serials(lines)
    next_serial = max(serials, default=0) + 1
    inserted_lines: list[str] = []
    seen: set[tuple[str, str, str, str, str, str]] = set()
    for record in records:
        if record.identity in existing or record.identity in seen:
            continue
        serial = record.serial
        if serial is None or serial in serials:
            serial = next_serial
            next_serial += 1
        serials.add(serial)
        inserted_lines.append(_rewrite_serial(record.line, serial))
        seen.add(record.identity)
        existing.add(record.identity)

    if not inserted_lines:
        return 0

    insert_at = _tail_insert_index(lines)
    updated = [*lines[:insert_at], *inserted_lines, *lines[insert_at:]]
    target_path.write_text("\n".join(updated).rstrip() + "\n", encoding="utf-8")
    return len(inserted_lines)


def rescue_retained_hets(
    *,
    source_pdb: str | Path,
    target_pdb: str | Path,
    stage: str,
    cfg: Mapping[str, object] | None = None,
    logger: logging.Logger | None = None,
    audit_path: str | Path | None = None,
) -> int:
    """Reinsert source metal/cofactor HET records missing from a prepared PDB."""

    log = logger or logging.getLogger(__name__)
    source = Path(source_pdb)
    target = Path(target_pdb)
    metals, cofactors = _canonical_sets(cfg)
    source_records = _collect_retained_records(source, metals, cofactors)
    target_before = _collect_retained_records(target, metals, cofactors)
    before_metals, before_cofactors = _count_records(target_before)
    inserted = 0
    if source_records and target.exists():
        inserted = _insert_records_before_tail(target, source_records)
        if inserted:
            try:
                fix_pdb_elements(target)
            except Exception as exc:
                log.warning(
                    "[metal.retention.elemfix] stage=%s file=%s err=%s",
                    stage,
                    target,
                    exc,
                )
    target_after = _collect_retained_records(target, metals, cofactors)
    after_metals, after_cofactors = _count_records(target_after)
    _write_rescue_audit(
        audit_path=Path(audit_path) if audit_path else _default_audit_path(target, stage),
        source=source,
        target=target,
        stage=stage,
        source_records=source_records,
        target_before=target_before,
        target_after=target_after,
        inserted=inserted,
    )
    log.info(
        "[metal.retention] stage=%s source=%s target=%s candidates=%d inserted=%d metals=%d->%d cofactors=%d->%d",
        stage,
        source,
        target,
        len(source_records),
        inserted,
        before_metals,
        after_metals,
        before_cofactors,
        after_cofactors,
    )
    return inserted


def rescue_metal_bound_hets(
    *,
    source_pdb: str | Path,
    target_pdb: str | Path,
    stage: str,
    cfg: Mapping[str, object] | None = None,
    logger: logging.Logger | None = None,
    audit_path: str | Path | None = None,
    radius_a: float = METAL_BOUND_HET_MAX_A,
) -> int:
    """Reinsert non-water HET residues that coordinate source metals.

    This is intended for repair stages that temporarily produce protein-only
    PDBs, such as restrained minimization feasibility checks.  It operates on
    the already-prepared source PDB, so competitive ligands previously removed
    by APO/HOLO or docking-box policy are not reintroduced.
    """

    log = logger or logging.getLogger(__name__)
    source = Path(source_pdb)
    target = Path(target_pdb)
    metals, _ = _canonical_sets(cfg)
    source_records = _collect_metal_bound_het_records(
        source,
        metals=metals,
        radius_a=radius_a,
    )
    inserted = 0
    if source_records and target.exists():
        inserted = _insert_records_before_tail(target, source_records)
        if inserted:
            try:
                fix_pdb_elements(target)
            except Exception as exc:
                log.warning(
                    "[metal.bound_het.elemfix] stage=%s file=%s err=%s",
                    stage,
                    target,
                    exc,
                )
    target_after = _records_present_in_target(target, source_records)
    _write_rescue_audit(
        audit_path=Path(audit_path) if audit_path else _default_audit_path(target, stage),
        source=source,
        target=target,
        stage=stage,
        source_records=source_records,
        target_before=[],
        target_after=target_after,
        inserted=inserted,
    )
    log.info(
        "[metal.bound_het.retention] stage=%s source=%s target=%s candidates=%d inserted=%d",
        stage,
        source,
        target,
        len(source_records),
        inserted,
    )
    return inserted


def _collect_metal_bound_het_records(
    path: Path,
    *,
    metals: set[str],
    radius_a: float,
) -> list[RetainedHetRecord]:
    if not path.exists():
        return []
    metal_points, het_by_residue, het_points = _scan_metal_bound_het_candidates(
        path,
        metals=metals,
    )
    if not metal_points or not het_points:
        return []
    radius2 = radius_a * radius_a
    bound_keys = {
        key
        for key, xyz in het_points
        if any(_distance2(xyz, metal_xyz) <= radius2 for metal_xyz in metal_points)
    }
    records: list[RetainedHetRecord] = []
    for key in sorted(bound_keys):
        records.extend(het_by_residue.get(key, ()))
    return records


def _scan_metal_bound_het_candidates(
    path: Path,
    *,
    metals: set[str],
) -> tuple[
    list[tuple[float, float, float]],
    dict[tuple[str, str, str, str], list[RetainedHetRecord]],
    list[tuple[tuple[str, str, str, str], tuple[float, float, float]]],
]:
    metal_points: list[tuple[float, float, float]] = []
    het_by_residue: dict[tuple[str, str, str, str], list[RetainedHetRecord]] = {}
    het_points: list[tuple[tuple[str, str, str, str], tuple[float, float, float]]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            _record_metal_bound_het_candidate(
                line,
                metals=metals,
                metal_points=metal_points,
                het_by_residue=het_by_residue,
                het_points=het_points,
            )
    return metal_points, het_by_residue, het_points


def _record_metal_bound_het_candidate(
    line: str,
    *,
    metals: set[str],
    metal_points: list[tuple[float, float, float]],
    het_by_residue: dict[tuple[str, str, str, str], list[RetainedHetRecord]],
    het_points: list[tuple[tuple[str, str, str, str], tuple[float, float, float]]],
) -> None:
    if not line.startswith("HETATM"):
        return
    resname = line[17:20].strip().upper()
    xyz = line_xyz(line)
    if xyz is None:
        return
    if resname in metals or _line_element(line) in metals:
        metal_points.append(xyz)
        return
    if resname in WATER_NAMES:
        return
    key = _residue_key(line)
    record = RetainedHetRecord(
        kind="metal_bound_het",
        line=line.rstrip("\n"),
        identity=_line_identity(line),
        token=resname,
        serial=_line_serial(line),
    )
    het_by_residue.setdefault(key, []).append(record)
    het_points.append((key, xyz))


def _records_present_in_target(
    target: Path,
    source_records: Sequence[RetainedHetRecord],
) -> list[RetainedHetRecord]:
    if not target.exists() or not source_records:
        return []
    wanted = {record.identity: record for record in source_records}
    present: list[RetainedHetRecord] = []
    with target.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("HETATM"):
                continue
            identity = _line_identity(line)
            if identity in wanted:
                present.append(wanted[identity])
    return present


def _residue_key(line: str) -> tuple[str, str, str, str]:
    return (
        line[17:20].strip().upper(),
        line[21:22].strip() or "-",
        line[22:26].strip(),
        line[26:27].strip(),
    )


def _distance2(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def _default_audit_path(target: Path, stage: str) -> Path:
    safe_stage = "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in stage)
    return target.with_suffix(f".{safe_stage}.retained_hets_audit.json")


def _write_rescue_audit(
    *,
    audit_path: Path,
    source: Path,
    target: Path,
    stage: str,
    source_records: Sequence[RetainedHetRecord],
    target_before: Sequence[RetainedHetRecord],
    target_after: Sequence[RetainedHetRecord],
    inserted: int,
) -> None:
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    before_metals, before_cofactors = _count_records(target_before)
    after_metals, after_cofactors = _count_records(target_after)
    missing_after = _missing_records(source_records, target_after)
    payload = {
        "stage": stage,
        "source_pdb": str(source),
        "target_pdb": str(target),
        "candidate_count": len(source_records),
        "inserted_count": inserted,
        "before": {"metals": before_metals, "cofactors": before_cofactors},
        "after": {"metals": after_metals, "cofactors": after_cofactors},
        "missing_after_count": len(missing_after),
        "missing_after": [_record_audit_row(record) for record in missing_after],
        "candidates": [_record_audit_row(record) for record in source_records],
    }
    audit_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _missing_records(
    source_records: Sequence[RetainedHetRecord],
    target_records: Sequence[RetainedHetRecord],
) -> list[RetainedHetRecord]:
    target_ids = {record.identity for record in target_records}
    return [record for record in source_records if record.identity not in target_ids]


def _record_audit_row(record: RetainedHetRecord) -> dict[str, object]:
    record_type, chain, resseq, icode, resname, atom_name = record.identity
    return {
        "kind": record.kind,
        "token": record.token,
        "record": record_type,
        "chain": chain,
        "resseq": resseq,
        "icode": icode,
        "resname": resname,
        "atom_name": atom_name,
        "serial": record.serial,
    }
