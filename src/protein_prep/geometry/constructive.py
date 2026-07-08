"""Targeted constructive geometry repair for binding-site prep."""

from __future__ import annotations

import json
import logging
import math
import os
import shutil
from pathlib import Path
from typing import Any, Mapping, Sequence

from protein_prep.geometry.audit import (
    apply_conservative_geometry_fixes,
    audit_geometry,
)
from protein_prep.openmm_repair import repair_with_pdbfixer
from protein_prep.pdb_records import line_xyz

RESIDUAL_BINDING_SITE_RADIUS_A = 8.0
RESIDUAL_SEVERE_CLASH_A = 1.2


def _heavy_atom_coord_from_line(line: str) -> tuple[float, float, float] | None:
    if not line.startswith(("ATOM  ", "HETATM")):
        return None
    if _line_element(line) == "H":
        return None
    return line_xyz(line)


def _heavy_atom_coords_from_pdb(path: Path) -> list[tuple[float, float, float]]:
    coords: list[tuple[float, float, float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            xyz = _heavy_atom_coord_from_line(line)
            if xyz is not None:
                coords.append(xyz)
    return coords


def ligand_center_from_pdb(path: Path | None) -> tuple[float, float, float] | None:
    """Return the heavy-atom centroid of an extracted ligand PDB."""

    if path is None or not path.exists():
        return None
    coords = _heavy_atom_coords_from_pdb(path)
    if not coords:
        return None
    return (
        sum(x for x, _, _ in coords) / len(coords),
        sum(y for _, y, _ in coords) / len(coords),
        sum(z for _, _, z in coords) / len(coords),
    )


def apply_targeted_binding_site_geometry_policy(
    input_pdb: Path,
    output_pdb: Path,
    *,
    reference_pdb: Path | None,
    ligand_pdb: Path | None,
    work_dir: Path,
    target_ph: float,
    audit_prefix: str,
    logger: logging.Logger | None = None,
) -> dict[str, object]:
    """
    Apply publication-oriented geometry cleanup without global sidechain repacking.

    Policy:
    - remove clashing waters conservatively;
    - preserve resolved sidechains;
    - if newly added sidechain atoms clash near the ligand, try restrained repair;
    - if repair cannot improve the site, leave a review sidecar instead of
      silently deleting binding-site sidechain atoms.
    """

    log = logger or logging.getLogger(__name__)
    work_dir.mkdir(parents=True, exist_ok=True)
    ligand_center = ligand_center_from_pdb(ligand_pdb)
    pre_audit = audit_geometry(
        input_pdb,
        sidecar_path=work_dir / f"{audit_prefix}_geometry_clash_audit_pre.json",
    )
    water_only = _apply_water_only_cleanup(
        input_pdb,
        work_dir / f"{audit_prefix}_water_geometry_fixed.pdb",
        work_dir=work_dir,
        audit_prefix=audit_prefix,
    )
    full_fix = apply_conservative_geometry_fixes(
        input_pdb,
        work_dir / f"{audit_prefix}_full_geometry_fixed.pdb",
        reference_pdb=reference_pdb,
        sidecar_path=work_dir / f"{audit_prefix}_geometry_cleanup_full.json",
        ligand_center=ligand_center,
    )
    chosen = _choose_initial_candidate(input_pdb, water_only, full_fix)
    constructive = _maybe_constructive_site_repair(
        input_pdb,
        full_fix,
        work_dir=work_dir,
        audit_prefix=audit_prefix,
        target_ph=target_ph,
        ligand_center=ligand_center,
        log=log,
    )
    if constructive.get("geometry_fix_constructive_repair_accepted"):
        chosen = Path(str(constructive["geometry_fix_constructive_repair_pdb"]))

    if chosen != input_pdb and chosen.exists():
        output_pdb.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(chosen, output_pdb)
        final_input = output_pdb
    else:
        final_input = input_pdb
    final_audit_path = work_dir / f"{audit_prefix}_geometry_clash_audit_final.json"
    final_audit = audit_geometry(
        final_input,
        sidecar_path=final_audit_path,
    )
    residual = _maybe_escalate_residual_clashes(
        final_input,
        output_pdb,
        final_audit=final_audit,
        final_audit_path=final_audit_path,
        reference_pdb=reference_pdb,
        ligand_pdb=ligand_pdb,
        work_dir=work_dir,
        audit_prefix=audit_prefix,
        target_ph=target_ph,
        log=log,
    )
    if residual.get("geometry_residual_repair_accepted"):
        accepted_status = str(residual.get("geometry_residual_repair_status", ""))
        final_input = Path(str(residual["geometry_residual_repair_pdb"]))
        if final_input != output_pdb and final_input.exists():
            output_pdb.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(final_input, output_pdb)
            final_input = output_pdb
        final_audit = audit_geometry(final_input, sidecar_path=final_audit_path)
        residual = _classify_residual_clashes(
            final_audit_path,
            ligand_pdb=ligand_pdb,
            attempted_repair=True,
            repair_status=accepted_status,
            repair_pdb=str(final_input),
        )
        residual["geometry_residual_repair_accepted"] = True
        residual_sidecar = work_dir / f"{audit_prefix}_geometry_residual_escalation.json"
        residual["geometry_residual_escalation_audit_path"] = str(residual_sidecar)
        residual_sidecar.write_text(
            json.dumps(residual, indent=2, sort_keys=True),
            encoding="utf-8",
        )
    summary: dict[str, object] = {
        "binding_site_geometry_policy_status": _policy_status(
            input_pdb,
            final_input,
            full_fix,
            constructive,
            residual,
        ),
        "binding_site_geometry_policy_output_pdb": str(output_pdb)
        if final_input == output_pdb
        else "",
        "binding_site_geometry_ligand_pdb": str(ligand_pdb) if ligand_pdb else "",
        "binding_site_geometry_ligand_center": _round_center(ligand_center),
        "binding_site_geometry_pre_clash_count": pre_audit.get(
            "geometry_clash_count", 0
        ),
        "binding_site_geometry_final_clash_count": final_audit.get(
            "geometry_clash_count", 0
        ),
        "binding_site_geometry_residual_gate_status": residual.get(
            "geometry_residual_gate_status", ""
        ),
        "binding_site_geometry_residual_near_ligand_clash_count": residual.get(
            "geometry_residual_near_ligand_clash_count", 0
        ),
        "binding_site_geometry_residual_severe_clash_count": residual.get(
            "geometry_residual_severe_clash_count", 0
        ),
        "binding_site_geometry_residual_offsite_clash_count": residual.get(
            "geometry_residual_offsite_clash_count", 0
        ),
        "binding_site_geometry_residual_repair_status": residual.get(
            "geometry_residual_repair_status", ""
        ),
        "binding_site_geometry_water_cleanup_status": water_only.get(
            "geometry_fix_status", ""
        ),
        "binding_site_geometry_full_cleanup_status": full_fix.get(
            "geometry_fix_status", ""
        ),
        "binding_site_geometry_removed_water_count": water_only.get(
            "geometry_fix_removed_water_count", 0
        ),
        "binding_site_geometry_candidate_dropped_sidechain_atom_count": full_fix.get(
            "geometry_fix_dropped_sidechain_atom_count", 0
        ),
        "binding_site_geometry_site_sidechain_atom_count": full_fix.get(
            "geometry_fix_binding_site_sidechain_count", 0
        ),
        **constructive,
        **residual,
    }
    sidecar = work_dir / f"{audit_prefix}_binding_site_geometry_policy.json"
    sidecar.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    summary["binding_site_geometry_policy_audit_path"] = str(sidecar)
    return summary


def _maybe_escalate_residual_clashes(
    input_pdb: Path,
    output_pdb: Path,
    *,
    final_audit: Mapping[str, object],
    final_audit_path: Path,
    reference_pdb: Path | None,
    ligand_pdb: Path | None,
    work_dir: Path,
    audit_prefix: str,
    target_ph: float,
    log: logging.Logger,
) -> dict[str, object]:
    residual = _classify_residual_clashes(final_audit_path, ligand_pdb=ligand_pdb)
    residual_sidecar = work_dir / f"{audit_prefix}_geometry_residual_escalation.json"
    if str(final_audit.get("geometry_status", "")) != "clashes":
        residual["geometry_residual_escalation_audit_path"] = str(residual_sidecar)
        residual_sidecar.write_text(
            json.dumps(residual, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        return residual
    best_residual: dict[str, object] | None = None
    attempted_statuses: list[str] = []
    for label, restraint_k, max_iterations in _residual_minimization_protocols():
        repaired = work_dir / f"{audit_prefix}_residual_geometry_minimized_{label}.pdb"
        status = run_restrained_minimization_repair(
            input_pdb,
            repaired,
            target_ph=target_ph,
            log=log,
            restraint_k=restraint_k,
            max_iterations=max_iterations,
        )
        attempted_statuses.append(f"{label}:{status}")
        if status != "ok" or not repaired.exists():
            continue
        candidate_residual = _evaluate_residual_repair_candidate(
            repaired,
            original_pdb=input_pdb,
            reference_pdb=reference_pdb,
            ligand_pdb=ligand_pdb,
            work_dir=work_dir,
            audit_prefix=f"{audit_prefix}_{label}",
            repair_status=";".join(attempted_statuses),
        )
        if not _residual_repair_improved(candidate_residual, residual):
            continue
        if best_residual is None or _residual_repair_improved(
            candidate_residual,
            best_residual,
        ):
            best_residual = candidate_residual
        if candidate_residual.get("geometry_residual_gate_status") == "ok":
            break
    if best_residual is not None:
        return best_residual
    residual["geometry_residual_repair_status"] = "rejected_no_geometry_improvement"
    if attempted_statuses:
        residual["geometry_residual_repair_status"] = (
            str(residual["geometry_residual_repair_status"])
            + ":"
            + ";".join(attempted_statuses)
        )
    residual["geometry_residual_escalation_audit_path"] = str(residual_sidecar)
    residual_sidecar.write_text(
        json.dumps(residual, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return residual


def _residual_minimization_protocols() -> list[tuple[str, float, int]]:
    protocols = [
        ("standard", 10.0, 50),
        ("polish", 5.0, 125),
    ]
    if _relaxed_geometry_minimization_enabled():
        protocols.append(("relaxed", 2.0, 250))
    return protocols


def _relaxed_geometry_minimization_enabled() -> bool:
    override = os.environ.get("ATLAS_GEOMETRY_RELAXED_MINIMIZATION")
    if override is not None:
        return override.strip().lower() not in {"0", "false", "no", "off", "none"}
    if _geometry_test_context_enabled():
        return False
    return True


def _geometry_test_context_enabled() -> bool:
    test_mode = os.environ.get("TEST_MODE_ENABLE", "").strip().lower()
    if test_mode and test_mode not in {"0", "false", "no", "off", "none"}:
        return True
    if os.environ.get("ATLAS_GEOMETRY_TEST_PDB", "").strip():
        return True
    requested = {
        token.strip().upper()
        for raw in (
            os.environ.get("PDB_ID", ""),
            os.environ.get("PDB", ""),
            os.environ.get("PDB_CODE", ""),
            os.environ.get("SPECIFIED_PROTEINS", ""),
            os.environ.get("PDBS", ""),
        )
        for token in raw.replace(";", ",").split(",")
        if token.strip()
    }
    return bool(requested & {"TEST", "TEST_PDB"})


def _evaluate_residual_repair_candidate(
    repaired: Path,
    *,
    original_pdb: Path,
    reference_pdb: Path | None,
    ligand_pdb: Path | None,
    work_dir: Path,
    audit_prefix: str,
    repair_status: str,
) -> dict[str, object]:
    cleaned = work_dir / f"{audit_prefix}_residual_geometry_cleaned.pdb"
    cleanup = apply_conservative_geometry_fixes(
        repaired,
        cleaned,
        reference_pdb=reference_pdb or original_pdb,
        sidecar_path=work_dir / f"{audit_prefix}_geometry_cleanup_residual.json",
        ligand_center=ligand_center_from_pdb(ligand_pdb),
    )
    if _summary_int(cleanup, "geometry_fix_dropped_sidechain_atom_count") > 0:
        return {
            "geometry_residual_gate_status": "rejected_sidechain_atom_drop",
            "geometry_residual_clash_count": 999999,
            "geometry_residual_near_ligand_clash_count": 999999,
            "geometry_residual_severe_clash_count": 999999,
            "geometry_residual_offsite_clash_count": 999999,
            "geometry_residual_repair_attempted": True,
            "geometry_residual_repair_status": "rejected_sidechain_atom_drop:"
            + str(cleanup.get("geometry_fix_status", "unknown")),
            "geometry_residual_repair_accepted": False,
            "geometry_residual_repair_pdb": str(repaired),
            "geometry_residual_cleanup_status": str(
                cleanup.get("geometry_fix_status", "")
            ),
        }
    candidate = _cleanup_candidate(repaired, cleaned, cleanup)
    candidate_audit_path = work_dir / f"{audit_prefix}_geometry_clash_audit_residual.json"
    audit_geometry(candidate, sidecar_path=candidate_audit_path)
    candidate_residual = _classify_residual_clashes(
        candidate_audit_path,
        ligand_pdb=ligand_pdb,
        attempted_repair=True,
        repair_status=repair_status,
        repair_pdb=str(candidate),
    )
    candidate_residual["geometry_residual_repair_accepted"] = True
    candidate_residual["geometry_residual_repair_pdb"] = str(candidate)
    candidate_residual["geometry_residual_cleanup_status"] = str(
        cleanup.get("geometry_fix_status", "")
    )
    return candidate_residual


def _classify_residual_clashes(
    audit_sidecar: Path,
    *,
    ligand_pdb: Path | None,
    attempted_repair: bool = False,
    repair_status: str = "",
    repair_pdb: str = "",
) -> dict[str, object]:
    contacts = _load_actionable_clashes(audit_sidecar)
    ligand_coords = _heavy_atom_coords_from_ligand_context(ligand_pdb)
    near_ligand = 0
    severe = 0
    offsite = 0
    min_distance: float | None = None
    min_ligand_distance: float | None = None
    for contact in contacts:
        distance = _contact_distance(contact)
        if distance is not None:
            min_distance = (
                round(distance, 3)
                if min_distance is None
                else min(min_distance, round(distance, 3))
            )
            if distance < RESIDUAL_SEVERE_CLASH_A:
                severe += 1
        ligand_distance = _contact_min_ligand_distance(contact, ligand_coords)
        if ligand_distance is not None:
            min_ligand_distance = (
                round(ligand_distance, 3)
                if min_ligand_distance is None
                else min(min_ligand_distance, round(ligand_distance, 3))
            )
            if ligand_distance <= RESIDUAL_BINDING_SITE_RADIUS_A:
                near_ligand += 1
                continue
        offsite += 1
    if not contacts:
        gate_status = "ok"
    elif severe:
        gate_status = "severe_geometry_review"
    elif near_ligand:
        gate_status = "binding_site_geometry_review"
    else:
        gate_status = "offsite_geometry_review"
    return {
        "geometry_residual_gate_status": gate_status,
        "geometry_residual_clash_count": len(contacts),
        "geometry_residual_near_ligand_clash_count": near_ligand,
        "geometry_residual_severe_clash_count": severe,
        "geometry_residual_offsite_clash_count": offsite,
        "geometry_residual_binding_site_radius_a": RESIDUAL_BINDING_SITE_RADIUS_A,
        "geometry_residual_severe_cutoff_a": RESIDUAL_SEVERE_CLASH_A,
        "geometry_residual_min_clash_distance_a": (
            "" if min_distance is None else min_distance
        ),
        "geometry_residual_min_ligand_distance_a": (
            "" if min_ligand_distance is None else min_ligand_distance
        ),
        "geometry_residual_repair_attempted": attempted_repair,
        "geometry_residual_repair_status": repair_status,
        "geometry_residual_repair_accepted": False,
        "geometry_residual_repair_pdb": repair_pdb,
    }


def _load_actionable_clashes(audit_sidecar: Path) -> list[dict[str, object]]:
    try:
        payload = json.loads(audit_sidecar.read_text(encoding="utf-8"))
    except Exception:
        return []
    rows = payload.get("clashes", [])
    if not isinstance(rows, list):
        return []
    return [
        row
        for row in rows
        if isinstance(row, dict) and not bool(row.get("ignored", False))
    ]


def _heavy_atom_coords_from_ligand_context(
    ligand_pdb: Path | None,
) -> list[tuple[float, float, float]]:
    if ligand_pdb is None or not ligand_pdb.exists():
        return []
    ligand_paths = [ligand_pdb]
    try:
        siblings = sorted(
            path
            for path in ligand_pdb.parent.glob("*.pdb")
            if path.is_file() and path != ligand_pdb
        )
        ligand_paths.extend(siblings)
    except Exception:
        pass
    coords: list[tuple[float, float, float]] = []
    for path in ligand_paths:
        coords.extend(_heavy_atom_coords_from_pdb(path))
    return coords


def _contact_distance(contact: Mapping[str, object]) -> float | None:
    try:
        return float(str(contact.get("distance_a", "")))
    except (TypeError, ValueError):
        return None


def _contact_min_ligand_distance(
    contact: Mapping[str, object],
    ligand_coords: Sequence[tuple[float, float, float]],
) -> float | None:
    if not ligand_coords:
        return None
    distances: list[float] = []
    for atom_key in ("atom_a", "atom_b"):
        atom = contact.get(atom_key)
        if not isinstance(atom, Mapping):
            continue
        xyz = _xyz_from_atom_record(atom)
        if xyz is None:
            continue
        distances.append(min(math.dist(xyz, ligand_xyz) for ligand_xyz in ligand_coords))
    return min(distances) if distances else None


def _xyz_from_atom_record(atom: Mapping[str, Any]) -> tuple[float, float, float] | None:
    raw = atom.get("xyz")
    if not isinstance(raw, Sequence) or len(raw) != 3:
        return None
    try:
        return float(raw[0]), float(raw[1]), float(raw[2])
    except (TypeError, ValueError):
        return None


def _residual_blocker_count(summary: Mapping[str, object]) -> int:
    return _summary_int(summary, "geometry_residual_near_ligand_clash_count") + _summary_int(
        summary,
        "geometry_residual_severe_clash_count",
    )


def _residual_repair_improved(
    candidate: Mapping[str, object],
    original: Mapping[str, object],
) -> bool:
    candidate_rank = _residual_gate_rank(candidate)
    original_rank = _residual_gate_rank(original)
    if candidate_rank > original_rank:
        return False
    candidate_blockers = _residual_blocker_count(candidate)
    original_blockers = _residual_blocker_count(original)
    if candidate_blockers > original_blockers:
        return False
    if candidate_rank < original_rank or candidate_blockers < original_blockers:
        return True
    return _summary_int(candidate, "geometry_residual_clash_count") < _summary_int(
        original,
        "geometry_residual_clash_count",
    )


def _residual_gate_rank(summary: Mapping[str, object]) -> int:
    gate = str(summary.get("geometry_residual_gate_status", ""))
    return {
        "ok": 0,
        "offsite_geometry_review": 1,
        "binding_site_geometry_review": 2,
        "severe_geometry_review": 3,
    }.get(gate, 4)


def _apply_water_only_cleanup(
    input_pdb: Path,
    output_pdb: Path,
    *,
    work_dir: Path,
    audit_prefix: str,
) -> dict[str, object]:
    return apply_conservative_geometry_fixes(
        input_pdb,
        output_pdb,
        reference_pdb=None,
        sidecar_path=work_dir / f"{audit_prefix}_geometry_cleanup_water_only.json",
        ligand_center=None,
    )


def _choose_initial_candidate(
    input_pdb: Path,
    water_only: Mapping[str, object],
    full_fix: Mapping[str, object],
) -> Path:
    sidechain_drops = _summary_int(
        full_fix,
        "geometry_fix_dropped_sidechain_atom_count",
    )
    if sidechain_drops == 0:
        candidate = str(full_fix.get("geometry_fix_output_pdb", "") or "")
        return Path(candidate) if candidate else input_pdb
    candidate = str(water_only.get("geometry_fix_output_pdb", "") or "")
    return Path(candidate) if candidate else input_pdb


def _maybe_constructive_site_repair(
    input_pdb: Path,
    full_fix: Mapping[str, object],
    *,
    work_dir: Path,
    audit_prefix: str,
    target_ph: float,
    ligand_center: tuple[float, float, float] | None,
    log: logging.Logger,
) -> dict[str, object]:
    if _summary_int(full_fix, "geometry_fix_binding_site_sidechain_count") <= 0:
        return {
            "geometry_fix_constructive_repair_status": "not_needed",
            "geometry_fix_constructive_repair_accepted": False,
            "geometry_fix_constructive_repair_pdb": "",
        }
    repaired = work_dir / f"{audit_prefix}_constructive_geometry_repaired.pdb"
    status = run_restrained_minimization_repair(
        input_pdb,
        repaired,
        target_ph=target_ph,
        log=log,
    )
    summary: dict[str, object] = {
        "geometry_fix_constructive_repair_status": status,
        "geometry_fix_constructive_repair_accepted": False,
        "geometry_fix_constructive_repair_pdb": str(repaired) if repaired.exists() else "",
    }
    if status != "ok" or not repaired.exists():
        return summary
    cleaned = work_dir / f"{audit_prefix}_constructive_geometry_cleaned.pdb"
    cleanup = apply_conservative_geometry_fixes(
        repaired,
        cleaned,
        reference_pdb=input_pdb,
        sidecar_path=work_dir / f"{audit_prefix}_geometry_cleanup_constructive.json",
        ligand_center=ligand_center,
    )
    if _summary_int(cleanup, "geometry_fix_dropped_sidechain_atom_count") > 0:
        summary["geometry_fix_constructive_repair_status"] = (
            "failed_constructive_cleanup:"
            + str(cleanup.get("geometry_fix_status", "unknown"))
        )
        return summary
    candidate = _cleanup_candidate(repaired, cleaned, cleanup)
    final_audit = audit_geometry(
        candidate,
        sidecar_path=work_dir / f"{audit_prefix}_geometry_clash_audit_constructive.json",
    )
    if str(final_audit.get("geometry_status", "")) == "clashes":
        summary["geometry_fix_constructive_repair_status"] = "failed_geometry_clashes"
        return summary
    summary["geometry_fix_constructive_repair_accepted"] = True
    summary["geometry_fix_constructive_repair_pdb"] = str(candidate)
    summary["geometry_fix_constructive_cleanup_status"] = str(
        cleanup.get("geometry_fix_status", "")
    )
    return summary


def run_restrained_minimization_repair(
    input_pdb: Path,
    output_pdb: Path,
    *,
    target_ph: float,
    log: logging.Logger | None = None,
    restraint_k: float = 10.0,
    max_iterations: int = 50,
) -> str:
    """Run a short restrained OpenMM cleanup and restore current receptor HETATM."""

    logger = log or logging.getLogger(__name__)
    protein_only = output_pdb.with_suffix(".protein_only.pdb")
    _write_protein_minimization_input(input_pdb, protein_only)
    repaired_input = output_pdb.with_suffix(".protein_repaired.pdb")
    if repair_with_pdbfixer(
        protein_only,
        repaired_input,
        target_ph=target_ph,
        cfg={"PDBFIXER_ADD_MISSING_RESIDUES": False},
        add_hydrogens=False,
        logger=logger,
    ):
        protein_only = repaired_input
    try:
        from openmm import CustomExternalForce, LangevinIntegrator, Platform, unit
        from openmm.app import ForceField, Modeller, PDBFile, Simulation

        pdb = PDBFile(str(protein_only))
        forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(forcefield, pH=float(target_ph))
        system = forcefield.createSystem(modeller.topology, constraints=None)
        restraint = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        restraint.addGlobalParameter("k", float(restraint_k))
        restraint.addPerParticleParameter("x0")
        restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")
        positions = modeller.positions
        for index, atom in enumerate(modeller.topology.atoms()):
            if atom.element is None or atom.element.symbol == "H":
                continue
            pos = positions[index].value_in_unit(unit.nanometer)
            restraint.addParticle(index, [pos.x, pos.y, pos.z])
        system.addForce(restraint)
        integrator = LangevinIntegrator(
            300 * unit.kelvin,
            1 / unit.picosecond,
            0.002 * unit.picoseconds,
        )
        platform, properties = _openmm_minimization_platform(
            Platform,
            atom_count=sum(1 for _ in modeller.topology.atoms()),
        )
        simulation = Simulation(
            modeller.topology,
            system,
            integrator,
            platform,
            properties,
        )
        simulation.context.setPositions(positions)
        simulation.minimizeEnergy(maxIterations=int(max_iterations))
        state = simulation.context.getState(getPositions=True)
        output_pdb.parent.mkdir(parents=True, exist_ok=True)
        with output_pdb.open("w", encoding="utf-8") as handle:
            PDBFile.writeFile(
                simulation.topology,
                state.getPositions(asNumpy=True),
                handle,
                keepIds=True,
            )
        _restore_hetatm_records(input_pdb, output_pdb)
        return "ok" if output_pdb.exists() and output_pdb.stat().st_size > 0 else "empty_output"
    except Exception as exc:
        return f"failed:{str(exc)[:180]}"


def _openmm_minimization_platform(
    platform_cls: Any,
    *,
    atom_count: int,
) -> tuple[Any, dict[str, str]]:
    """Prefer OpenMM CPU; Reference is too slow for large stress-test receptors."""

    max_threads = _bounded_openmm_threads()
    try:
        platform = platform_cls.getPlatformByName("CPU")
        return platform, {"Threads": str(max_threads)}
    except Exception:
        pass
    reference_max_atoms = _env_int("ATLAS_GEOMETRY_REFERENCE_MAX_ATOMS", 12000)
    if atom_count > reference_max_atoms:
        raise RuntimeError(
            "OpenMM CPU platform unavailable and receptor too large for Reference "
            f"minimization: atoms={atom_count} limit={reference_max_atoms}"
        )
    return platform_cls.getPlatformByName("Reference"), {}


def _bounded_openmm_threads() -> int:
    requested = _env_int("ATLAS_GEOMETRY_OPENMM_THREADS", 0)
    if requested <= 0:
        allocation = _env_int("SLURM_CPUS_PER_TASK", 0)
        if allocation <= 0:
            allocation = _env_int("CPU", 1)
        requested = min(max(1, allocation), 2)
    return max(1, min(requested, 8))


def _env_int(name: str, default: int) -> int:
    try:
        return int(float(str(os.environ.get(name, default)).strip()))
    except (TypeError, ValueError):
        return default


def _write_protein_minimization_input(input_pdb: Path, output_pdb: Path) -> None:
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with input_pdb.open("r", encoding="utf-8", errors="ignore") as src:
        with output_pdb.open("w", encoding="utf-8") as dst:
            for line in src:
                if line.startswith("ATOM  ") and _line_element(line) != "H":
                    dst.write(line)
            dst.write("END\n")


def _restore_hetatm_records(source_pdb: Path, target_pdb: Path) -> None:
    records: list[str] = []
    with source_pdb.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("HETATM"):
                records.append(line)
    if records:
        _append_records_before_end(target_pdb, records)


def _append_records_before_end(path: Path, records: Sequence[str]) -> None:
    original = path.read_text(encoding="utf-8", errors="ignore").splitlines(
        keepends=True
    )
    output: list[str] = []
    inserted = False
    for line in original:
        if line.startswith("END") and not inserted:
            output.extend(records)
            inserted = True
        output.append(line)
    if not inserted:
        output.extend(records)
        output.append("END\n")
    path.write_text("".join(output), encoding="utf-8")


def _cleanup_candidate(
    repaired: Path,
    cleaned: Path,
    cleanup: Mapping[str, object],
) -> Path:
    if str(cleanup.get("geometry_fix_status", "")) == "unchanged":
        return repaired
    return cleaned if cleaned.exists() else repaired


def _policy_status(
    input_pdb: Path,
    final_pdb: Path,
    full_fix: Mapping[str, object],
    constructive: Mapping[str, object],
    residual: Mapping[str, object],
) -> str:
    residual_gate = str(residual.get("geometry_residual_gate_status", ""))
    if bool(residual.get("geometry_residual_repair_accepted")):
        if residual_gate == "ok":
            return "residual_geometry_repair_applied"
        if residual_gate == "offsite_geometry_review":
            return "residual_geometry_repair_applied_offsite_review"
    if bool(constructive.get("geometry_fix_constructive_repair_accepted")):
        return "constructive_binding_site_repair_applied"
    if _summary_int(full_fix, "geometry_fix_binding_site_sidechain_count") > 0:
        return "binding_site_sidechain_review_required"
    if residual_gate in {"severe_geometry_review", "binding_site_geometry_review"}:
        return residual_gate
    if residual_gate == "offsite_geometry_review":
        return "offsite_geometry_review"
    return "geometry_cleanup_applied" if final_pdb != input_pdb else "unchanged"


def _summary_int(summary: Mapping[str, object], key: str) -> int:
    try:
        return int(str(summary.get(key, 0) or 0))
    except (TypeError, ValueError):
        return 0


def _line_element(line: str) -> str:
    element = line[76:78].strip() if len(line) >= 78 else ""
    if element:
        return element.upper()
    return line[12:16].strip()[:1].upper()


def _round_center(
    center: tuple[float, float, float] | None,
) -> list[float] | str:
    if center is None:
        return ""
    return [round(value, 3) for value in center]
