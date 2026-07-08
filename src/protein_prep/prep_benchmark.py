"""Open receptor-preparation benchmark harness."""

from __future__ import annotations

import json
import logging
import math
import os
import shutil
import sys
import urllib.request
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

from protein_prep.pdb_records import (
    line_resseq_int as _line_resseq_int,
    line_xyz as _line_xyz,
    residue_key as _residue_key,
)
from protein_prep.autodock4zn import (
    build_autodock4zn_plan,
    run_autodock4zn_plan,
    summarize_autodock4zn_plan,
    write_autodock4zn_plan,
)
from protein_prep.het_atom.ccd import write_ccd_instance_sdf
from protein_prep.binding_site_network_apply import (
    apply_binding_site_network_recommendations,
)
from protein_prep.geometry.audit import (
    apply_conservative_geometry_fixes,
    audit_geometry,
)
from protein_prep.het_state_selection import select_het_states
from protein_prep.metal_chemistry import (
    audit_metal_chemistry,
    summarize_metal_chemistry,
)
from protein_prep.metal_publication_policy import classify_metal_publication_policy
from protein_prep.metal_parameterization import (
    summarize_parameterization_plan,
    write_parameterization_plan,
)
from protein_prep.openmm_repair import get_target_ph_for_prep, repair_with_pdbfixer
from protein_prep.metal_retention import rescue_metal_bound_hets, rescue_retained_hets
from protein_prep.protonation import (
    REDUCE_EXE,
    _protonate_with_pdb2pqr_if_available,
    assign_protonation_states,
)
from protein_prep.tool_runners import (
    build_meeko_base_cmd,
    build_meeko_legacy_cmd,
    build_meeko_modern_cmd,
    run_subprocess_capture,
)
from protein_prep.meeko.normalization import (
    meeko_pruned_stage_label,
    write_meeko_clash_pruned_copy,
    write_meeko_normalized_copy,
    write_meeko_ordered_copy,
)
from protein_prep.meeko.safe_export import (
    meeko_safe_stage_label,
    parse_unmatched_residue_keys,
    write_meeko_safe_receptor_copy,
)
from protein_prep.water_policy import audit_water_policy, write_supported_water_receptor
from protein_prep.water_evidence import augment_water_policy_evidence
from protein_prep.water_redock_sensitivity import (
    audit_water_redock_sensitivity,
    empty_water_redock_sensitivity_summary,
)
from protein_prep.meeko.stage_review import meeko_stage_requires_review


def _active_env_tool(name: str) -> Path:
    return Path(sys.executable).resolve().parent / name


def _tool_available(name: str, configured: str | None = None) -> bool:
    if configured:
        configured_path = Path(str(configured))
        if configured_path.is_absolute() and configured_path.exists():
            return True
        if shutil.which(str(configured)):
            return True
    env_candidate = _active_env_tool(name)
    return bool(
        env_candidate.exists()
        or shutil.which(name)
    )


DEFAULT_PDB_IDS: tuple[str, ...] = (
    "1IEP",
    "1HPX",
    "1HVR",
    "3PTB",
    "4DFR",
    "1M17",
    "3ERK",
    "2JDU",
    "2SRC",
    "3K5V",
    "4AG8",
    "4J8M",
    "5EW8",
    "3PBL",
    "2ZV2",
    "1OYT",
    "1BJU",
    "2RGP",
    "1K1I",
)

WATER_NAMES = {"HOH", "WAT", "H2O", "DOD"}
COMMON_IONS = {
    "AG",
    "AL",
    "BA",
    "BR",
    "CA",
    "CD",
    "CL",
    "CO",
    "CS",
    "CU",
    "FE",
    "HG",
    "IOD",
    "K",
    "LI",
    "MG",
    "MN",
    "NA",
    "NI",
    "RB",
    "SR",
    "ZN",
}
COORDINATION_METALS = {
    "CA",
    "CD",
    "CO",
    "CU",
    "FE",
    "HG",
    "MG",
    "MN",
    "NI",
    "ZN",
}
COFACTOR_NAMES = {"ADP", "ATP", "FAD", "FMN", "HEM", "NAD", "NAP", "SAM"}
PH_STATE_NAMES = {"ASH", "GLH", "HID", "HIE", "HIP", "LYN", "CYM", "TYM"}
PROTEIN_DONOR_ELEMENTS = {"N", "O", "S"}
METAL_DONOR_ELEMENTS = {"N", "O", "S"}
CAP_RESNAMES = {"ACE", "NME"}
PROTEIN_RECORD_RESNAMES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}


@dataclass
class StructureStats:
    heavy_atoms: int = 0
    hydrogens: int = 0
    waters: int = 0
    metals: int = 0
    cofactors: int = 0
    ligands: int = 0
    ph_state_residues: str = ""
    het_counts: str = ""


BenchmarkRow = dict[str, object]


@dataclass(frozen=True)
class HetResidue:
    resname: str
    chain: str
    resseq: str
    heavy_atoms: int


@dataclass(frozen=True)
class AtomPoint:
    serial: str
    resname: str
    chain: str
    resseq: str
    record: str
    atom_name: str
    element: str
    xyz: tuple[float, float, float]


@dataclass(frozen=True)
class MetalBoundHet:
    residue_key: tuple[str, str, str, str]
    classification: str
    retained_in_baseline: bool
    expected_removed: bool
    donor_count: int
    metal_ids: tuple[str, ...]
    reason: str


def _is_heavy_atom(line: str) -> bool:
    return _line_element(line) != "H"


def _line_element(line: str) -> str:
    return (line[76:78].strip() if len(line) >= 78 else line[12:16].strip()[:1]).upper()


def _format_counter(counter: Counter[str]) -> str:
    if not counter:
        return ""
    return ";".join(f"{key}:{counter[key]}" for key in sorted(counter))


def collect_structure_stats(path: Path) -> StructureStats:
    hydrogens = 0
    heavy_atoms = 0
    waters: set[tuple[str, str, str, str]] = set()
    ligands: set[tuple[str, str, str, str]] = set()
    metals = 0
    cofactors: set[tuple[str, str, str, str]] = set()
    ph_states: Counter[str] = Counter()
    het_counts: Counter[str] = Counter()

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            resname = line[17:20].strip().upper()
            is_heavy = _is_heavy_atom(line)
            if is_heavy:
                heavy_atoms += 1
            else:
                hydrogens += 1
            if resname in PH_STATE_NAMES:
                ph_states[resname] += 1
            if not line.startswith("HETATM"):
                continue
            key = _residue_key(line)
            het_counts[resname] += 1
            if resname in WATER_NAMES:
                waters.add(key)
            elif resname in COORDINATION_METALS:
                metals += 1 if is_heavy else 0
            elif resname in COFACTOR_NAMES:
                cofactors.add(key)
            else:
                ligands.add(key)

    return StructureStats(
        heavy_atoms=heavy_atoms,
        hydrogens=hydrogens,
        waters=len(waters),
        metals=metals,
        cofactors=len(cofactors),
        ligands=len(ligands),
        ph_state_residues=_format_counter(ph_states),
        het_counts=_format_counter(het_counts),
    )


def download_pdb(pdb_id: str, out_path: Path, *, timeout: int = 60) -> bool:
    if out_path.exists() and out_path.stat().st_size > 0:
        return True
    out_path.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            out_path.write_bytes(response.read())
        return out_path.stat().st_size > 0
    except Exception as exc:
        logging.warning("[benchmark] download failed pdb=%s err=%s", pdb_id, exc)
        return False


def write_water_variant(input_pdb: Path, output_pdb: Path, *, keep_water: bool) -> None:
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with input_pdb.open("r", encoding="utf-8", errors="ignore") as src:
        with output_pdb.open("w", encoding="utf-8") as dst:
            for line in src:
                if line.startswith("HETATM") and line[17:20].strip().upper() in WATER_NAMES:
                    if keep_water:
                        dst.write(line)
                    continue
                dst.write(line)


def extract_largest_ligand(input_pdb: Path, output_pdb: Path) -> tuple[str, int]:
    groups: dict[tuple[str, str, str, str], list[str]] = defaultdict(list)
    with input_pdb.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("HETATM") or not _is_heavy_atom(line):
                continue
            resname = line[17:20].strip().upper()
            if resname in WATER_NAMES or resname in COMMON_IONS:
                continue
            groups[_residue_key(line)].append(line)

    if not groups:
        return "", 0

    key, lines = max(groups.items(), key=lambda item: len(item[1]))
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with output_pdb.open("w", encoding="utf-8") as handle:
        handle.writelines(lines)
        handle.write("END\n")
    return key[0], len(lines)


def _het_status(status: str, charge: object = "", tool: str = "") -> dict[str, object]:
    return {
        "het_chemistry_status": status,
        "het_formal_charge": charge,
        "het_tool": tool,
    }


def _obabel_ligand_to_sdf(ligand_pdb: Path, ligand_sdf: Path) -> str:
    if not shutil.which("obabel"):
        return "missing_obabel"
    result = run_subprocess_capture(
        ["obabel", str(ligand_pdb), "-O", str(ligand_sdf), "-h"],
        timeout=120,
    )
    if result.returncode != 0 or not ligand_sdf.exists():
        return "obabel_failed"
    return "ok"


def _rdkit_sdf_charge(ligand_sdf: Path) -> tuple[str, object]:
    try:
        from rdkit import Chem

        supplier = Chem.SDMolSupplier(str(ligand_sdf), sanitize=True, removeHs=False)
        mol = next((item for item in supplier if item is not None), None)
        if mol is None:
            return "rdkit_sanitize_failed", ""
        charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        return "ok", charge
    except Exception as exc:
        return f"rdkit_error:{str(exc)[:160]}", ""


def download_ccd_sdf(resname: str, out_path: Path) -> bool:
    token = resname.strip().upper()
    if not token:
        return False
    if out_path.exists() and out_path.stat().st_size > 0:
        return True
    url = f"https://files.rcsb.org/ligands/download/{token}_ideal.sdf"
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            out_path.write_bytes(response.read())
        return out_path.stat().st_size > 0
    except Exception as exc:
        logging.info("[benchmark] CCD ligand SDF unavailable res=%s err=%s", token, exc)
        return False


def validate_ligand_chemistry(
    ligand_pdb: Path,
    work_dir: Path,
    ligand_resname: str,
) -> dict[str, object]:
    if not ligand_pdb.exists() or ligand_pdb.stat().st_size == 0:
        return _het_status("no_ligand")
    ccd_sdf = work_dir / f"{ligand_resname.upper()}_ideal.sdf"
    if download_ccd_sdf(ligand_resname, ccd_sdf):
        ccd_status, charge = _rdkit_sdf_charge(ccd_sdf)
        if ccd_status == "ok":
            return _het_status("ok", charge=charge, tool="rcsb_ccd+rdkit")

    ligand_sdf = work_dir / "crystal_ligand_chemistry.sdf"
    obabel_status = _obabel_ligand_to_sdf(ligand_pdb, ligand_sdf)
    if obabel_status != "ok":
        tool = "obabel" if obabel_status != "missing_obabel" else ""
        return _het_status(obabel_status, tool=tool)
    rdkit_status, charge = _rdkit_sdf_charge(ligand_sdf)
    return _het_status(rdkit_status, charge=charge, tool="obabel+rdkit")


def _rdkit_mol_from_sdf(sdf_path: Path) -> Any | None:
    try:
        from rdkit import Chem

        supplier = Chem.SDMolSupplier(str(sdf_path), sanitize=True, removeHs=False)
        return next((item for item in supplier if item is not None), None)
    except Exception:
        return None


def _mol_formal_charge(mol: Any) -> int:
    return int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))


def _canonical_smiles(mol: Any) -> str:
    from rdkit import Chem

    return str(Chem.MolToSmiles(mol, isomericSmiles=True))


def _enumerate_rdkit_tautomers(smiles: str, *, limit: int = 64) -> set[str]:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return set()
    enumerator = rdMolStandardize.TautomerEnumerator()
    enumerator.SetMaxTautomers(int(limit))
    return {_canonical_smiles(tautomer) for tautomer in enumerator.Enumerate(mol)}


def _enumerate_dimorphite_states(smiles: str, target_ph: float, *, limit: int = 32) -> tuple[set[str], str]:
    try:
        from dimorphite_dl import protonate_smiles

        variants = protonate_smiles(
            smiles,
            ph_min=max(0.0, float(target_ph) - 0.5),
            ph_max=min(14.0, float(target_ph) + 0.5),
            max_variants=int(limit),
            validate_output=True,
        )
        return {str(item).strip() for item in variants if str(item).strip()}, "dimorphite-dl"
    except Exception as exc:
        return set(), f"dimorphite_unavailable:{str(exc)[:120]}"


def _collect_het_residues(path: Path) -> list[HetResidue]:
    grouped: dict[tuple[str, str, str, str], int] = defaultdict(int)
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("HETATM") or not _is_heavy_atom(line):
                continue
            resname = line[17:20].strip().upper()
            if resname in WATER_NAMES or resname in COMMON_IONS:
                continue
            grouped[_residue_key(line)] += 1
    return [
        HetResidue(
            resname=key[0],
            chain=key[1],
            resseq=key[2],
            heavy_atoms=count,
        )
        for key, count in sorted(grouped.items(), key=lambda item: item[0])
    ]


def _audit_one_het_state(
    resname: str,
    work_dir: Path,
    target_ph: float,
) -> dict[str, object]:
    ccd_sdf = work_dir / f"{resname.upper()}_ideal.sdf"
    if not download_ccd_sdf(resname, ccd_sdf):
        return _het_state_failure(resname, "missing_ccd_sdf", "missing", "")
    mol = _rdkit_mol_from_sdf(ccd_sdf)
    if mol is None:
        return _het_state_failure(resname, "rdkit_sanitize_failed", "rcsb_ccd", "rdkit")

    base_smiles = _canonical_smiles(mol)
    candidate_smiles, dimorphite_status = _enumerate_dimorphite_states(base_smiles, target_ph)
    if not candidate_smiles:
        candidate_smiles = {base_smiles}
    states = _enumerate_state_smiles(candidate_smiles)
    charges = _state_charge_set(states)
    tool = "rcsb_ccd+rdkit"
    if dimorphite_status == "dimorphite-dl":
        tool += "+dimorphite-dl"
    return {
        "resname": resname,
        "status": "ok" if states else "no_states",
        "bond_order_source": "rcsb_ccd",
        "state_count": len(states),
        "charge_set": ",".join(str(charge) for charge in sorted(charges)),
        "tool": tool,
        "dimorphite_status": dimorphite_status,
        "representative_smiles": sorted(states)[:5],
    }


def _het_state_failure(
    resname: str,
    status: str,
    bond_order_source: str,
    tool: str,
) -> dict[str, object]:
    return {
        "resname": resname,
        "status": status,
        "bond_order_source": bond_order_source,
        "state_count": 0,
        "charge_set": "",
        "tool": tool,
    }


def _enumerate_state_smiles(candidate_smiles: set[str]) -> set[str]:
    states: set[str] = set()
    for smiles in candidate_smiles:
        tautomers = _enumerate_rdkit_tautomers(smiles)
        states.update(tautomers or {smiles})
    return states


def _state_charge_set(states: set[str]) -> set[int]:
    from rdkit import Chem

    charges: set[int] = set()
    for smiles in states:
        state_mol = Chem.MolFromSmiles(smiles)
        if state_mol is not None:
            charges.add(_mol_formal_charge(state_mol))
    return charges


def audit_het_states(path: Path, work_dir: Path, target_ph: float) -> tuple[dict[str, object], list[dict[str, object]]]:
    residues = _collect_het_residues(path)
    if not residues:
        return _empty_het_state_summary(), []

    by_resname = {residue.resname for residue in residues}
    audits = [
        _audit_one_het_state(resname, work_dir, target_ph)
        for resname in sorted(by_resname)
    ]
    detailed = [
        {
            **item,
            "residue_instances": sum(1 for residue in residues if residue.resname == item.get("resname")),
        }
        for item in audits
    ]
    status, ok_count = _het_audit_status(audits)
    return {
        "het_state_status": status,
        "het_state_residue_count": len(residues),
        "het_state_resname_count": len(by_resname),
        "het_state_variant_count": _het_total_states(audits),
        "het_state_charge_set": ",".join(_het_charge_set(audits)),
        "het_state_failures": ";".join(_het_failures(audits)),
        "het_state_tool": "rcsb_ccd+rdkit+dimorphite-dl",
        "het_bond_order_source": "rcsb_ccd" if ok_count else "missing_or_failed",
    }, detailed


def _empty_het_state_summary() -> dict[str, object]:
    return {
        "het_state_status": "no_het",
        "het_state_residue_count": 0,
        "het_state_resname_count": 0,
        "het_state_variant_count": 0,
        "het_state_charge_set": "",
        "het_state_failures": "",
        "het_state_tool": "",
        "het_bond_order_source": "",
    }


def _het_audit_status(audits: list[dict[str, object]]) -> tuple[str, int]:
    ok_count = Counter(str(item.get("status", "")) for item in audits).get("ok", 0)
    if ok_count == len(audits):
        return "ok", ok_count
    if ok_count:
        return "partial", ok_count
    return "failed", ok_count


def _het_total_states(audits: list[dict[str, object]]) -> int:
    return sum(int(str(item.get("state_count", 0) or 0)) for item in audits)


def _het_charge_set(audits: list[dict[str, object]]) -> list[str]:
    return sorted(
        {
            charge
            for item in audits
            for charge in str(item.get("charge_set", "")).split(",")
            if charge
        }
    )


def _het_failures(audits: list[dict[str, object]]) -> list[str]:
    return [
        f"{item.get('resname')}:{item.get('status')}"
        for item in audits
        if item.get("status") != "ok"
    ]


def _distance2(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def _collect_coordination_atoms(
    path: Path,
) -> tuple[
    list[tuple[str, tuple[str, str, str, str], tuple[float, float, float]]],
    list[tuple[str, str, tuple[float, float, float]]],
]:
    metals: list[tuple[str, tuple[str, str, str, str], tuple[float, float, float]]] = []
    donors: list[tuple[str, str, tuple[float, float, float]]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")) or not _is_heavy_atom(line):
                continue
            xyz = _line_xyz(line)
            if xyz is None:
                continue
            elem = _line_element(line)
            resname = line[17:20].strip().upper()
            if line.startswith("HETATM") and resname in COORDINATION_METALS:
                metals.append((elem or resname, _residue_key(line), xyz))
            elif elem in PROTEIN_DONOR_ELEMENTS:
                donors.append((line[12:16].strip().upper(), resname, xyz))
    return metals, donors


def _coordination_part(
    elem: str,
    key: tuple[str, str, str, str],
    xyz: tuple[float, float, float],
    donors: Sequence[tuple[str, str, tuple[float, float, float]]],
    radius2: float,
) -> str:
    nearby = Counter(
        resname for _, resname, donor_xyz in donors if _distance2(xyz, donor_xyz) <= radius2
    )
    donor_text = ",".join(f"{res}:{nearby[res]}" for res in sorted(nearby)) or "none"
    return f"{elem}:{key[1]}:{key[2]}:{sum(nearby.values())}[{donor_text}]"


def metal_coordination_signature(path: Path, *, radius: float = 3.0) -> str:
    metals, donors = _collect_coordination_atoms(path)
    if not metals:
        return ""
    radius2 = radius * radius
    return ";".join(_coordination_part(elem, key, xyz, donors, radius2) for elem, key, xyz in metals)


def metal_coordination_donor_count(path: Path, *, radius: float = 3.0) -> int:
    metals, donors = _collect_coordination_atoms(path)
    if not metals:
        return 0
    radius2 = radius * radius
    return sum(
        1
        for _, _, metal_xyz in metals
        for _, _, donor_xyz in donors
        if _distance2(metal_xyz, donor_xyz) <= radius2
    )


def audit_benchmark_metal_chemistry(
    path: Path,
    *,
    radius: float = 3.2,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    atoms = _collect_atom_points(path)
    metals = [
        atom
        for atom in atoms
        if atom.record == "HETATM" and atom.resname in COORDINATION_METALS
    ]
    rows = [_benchmark_metal_chemistry_row(metal, atoms, radius) for metal in metals]
    return summarize_metal_chemistry(rows), rows


def _benchmark_metal_chemistry_row(
    metal: AtomPoint,
    atoms: Sequence[AtomPoint],
    radius: float,
) -> dict[str, object]:
    donors = _benchmark_metal_donors(metal, atoms, radius)
    chemistry = audit_metal_chemistry(
        element=metal.element or metal.resname,
        metal_xyz=metal.xyz,
        donors=donors,
    )
    return {
        "id": f"{metal.element or metal.resname}:{metal.chain}:{metal.resseq}",
        "atom_serial": metal.serial,
        "resname": metal.resname,
        "chain": metal.chain,
        "resseq": metal.resseq,
        "atom_name": metal.atom_name,
        "coords": list(metal.xyz),
        "donors": donors,
        **chemistry,
    }


def _benchmark_metal_donors(
    metal: AtomPoint,
    atoms: Sequence[AtomPoint],
    radius: float,
) -> list[dict[str, object]]:
    radius2 = radius * radius
    donors: list[dict[str, object]] = []
    for atom in atoms:
        if atom.element not in METAL_DONOR_ELEMENTS:
            continue
        if atom is metal or _distance2(metal.xyz, atom.xyz) > radius2:
            continue
        donors.append(
            {
                "category": _benchmark_donor_category(atom),
                "resname": atom.resname,
                "chain": atom.chain,
                "resseq": atom.resseq,
                "atom_name": atom.atom_name,
                "element": atom.element,
                "coords": list(atom.xyz),
            }
        )
    return donors


def _benchmark_donor_category(atom: AtomPoint) -> str:
    if atom.record == "ATOM":
        return "protein"
    if atom.resname in WATER_NAMES:
        return "water"
    return "ligand"


def _point_inside_box(
    xyz: tuple[float, float, float],
    center: tuple[float, float, float] | None,
    box_size: float,
) -> bool:
    if center is None:
        return False
    half = box_size / 2.0
    return all(abs(xyz[idx] - center[idx]) <= half for idx in range(3))


def _collect_metal_bound_het_atoms(
    path: Path,
) -> tuple[
    list[tuple[str, tuple[str, str, str, str], tuple[float, float, float]]],
    dict[tuple[str, str, str, str], list[AtomPoint]],
]:
    metals: list[tuple[str, tuple[str, str, str, str], tuple[float, float, float]]] = []
    het_atoms: dict[tuple[str, str, str, str], list[AtomPoint]] = defaultdict(list)
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("HETATM") or not _is_heavy_atom(line):
                continue
            xyz = _line_xyz(line)
            if xyz is None:
                continue
            key = _residue_key(line)
            resname = key[0]
            element = _line_element(line)
            atom = AtomPoint(
                serial=line[6:11].strip(),
                resname=resname,
                chain=key[1],
                resseq=key[2],
                record="HETATM",
                atom_name=line[12:16].strip(),
                element=element,
                xyz=xyz,
            )
            if resname in COORDINATION_METALS:
                metals.append((element or resname, key, xyz))
            else:
                het_atoms[key].append(atom)
    return metals, het_atoms


def _classify_metal_bound_het(
    key: tuple[str, str, str, str],
    atoms: Sequence[AtomPoint],
    *,
    ligand_resname: str,
    ligand_center: tuple[float, float, float] | None,
    box_size: float,
) -> tuple[str, bool, str]:
    resname = key[0]
    in_box = any(_point_inside_box(atom.xyz, ligand_center, box_size) for atom in atoms)
    if resname in WATER_NAMES:
        return "bridging_water", False, "metal-bound water"
    if resname in COMMON_IONS:
        return "buffer/salt", True, "simple ion/salt"
    if resname in COFACTOR_NAMES:
        return "structural_cofactor", False, "known cofactor"
    if resname == ligand_resname:
        return "competitive_ligand", True, "matches extracted co-crystal ligand"
    if in_box:
        return "competitive_ligand", True, "overlaps docking box"
    return "unknown_structural_ligand", False, "metal-bound HET outside docking box"


def audit_metal_bound_hets(
    input_pdb: Path,
    baseline_pdb: Path,
    ligand_pdb: Path,
    ligand_resname: str,
    *,
    box_size: float,
    radius: float = 3.0,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    metals, het_atoms = _collect_metal_bound_het_atoms(input_pdb)
    _, baseline_hets = _collect_metal_bound_het_atoms(baseline_pdb)
    baseline_keys = set(baseline_hets)
    ligand_center = _centroid(ligand_pdb) if ligand_pdb.exists() else None
    bound_by_residue = _find_metal_bound_het_entries(
        metals,
        het_atoms,
        radius * radius,
    )
    rows = _metal_bound_het_rows(
        bound_by_residue,
        het_atoms,
        baseline_keys,
        ligand_resname=ligand_resname,
        ligand_center=ligand_center,
        box_size=box_size,
    )
    return _metal_bound_het_summary(rows), _metal_bound_het_details(rows)


def _find_metal_bound_het_entries(
    metals: Sequence[tuple[str, tuple[str, str, str, str], tuple[float, float, float]]],
    het_atoms: dict[tuple[str, str, str, str], list[AtomPoint]],
    radius2: float,
) -> dict[tuple[str, str, str, str], dict[str, object]]:
    bound_by_residue: dict[tuple[str, str, str, str], dict[str, object]] = {}
    for metal_elem, metal_key, metal_xyz in metals:
        metal_id = f"{metal_elem}:{metal_key[1]}:{metal_key[2]}:{metal_key[3] or '-'}"
        for key, atoms in het_atoms.items():
            donor_atoms = [
                atom
                for atom in atoms
                if atom.element in METAL_DONOR_ELEMENTS
                and _distance2(metal_xyz, atom.xyz) <= radius2
            ]
            if not donor_atoms:
                continue
            entry = bound_by_residue.setdefault(
                key,
                {
                    "residue_key": key,
                    "donor_count": 0,
                    "metal_ids": set(),
                },
            )
            entry["donor_count"] = int(str(entry["donor_count"])) + len(donor_atoms)
            metal_ids = entry["metal_ids"]
            if isinstance(metal_ids, set):
                metal_ids.add(metal_id)
    return bound_by_residue


def _metal_bound_het_rows(
    bound_by_residue: dict[tuple[str, str, str, str], dict[str, object]],
    het_atoms: dict[tuple[str, str, str, str], list[AtomPoint]],
    baseline_keys: set[tuple[str, str, str, str]],
    *,
    ligand_resname: str,
    ligand_center: tuple[float, float, float] | None,
    box_size: float,
) -> list[MetalBoundHet]:
    rows: list[MetalBoundHet] = []
    for key, entry in sorted(bound_by_residue.items()):
        classification, expected_removed, reason = _classify_metal_bound_het(
            key,
            het_atoms.get(key, []),
            ligand_resname=ligand_resname,
            ligand_center=ligand_center,
            box_size=box_size,
        )
        metal_ids_obj = entry.get("metal_ids", set())
        metal_ids = tuple(sorted(metal_ids_obj)) if isinstance(metal_ids_obj, set) else ()
        rows.append(
            MetalBoundHet(
                residue_key=key,
                classification=classification,
                retained_in_baseline=key in baseline_keys,
                expected_removed=expected_removed,
                donor_count=int(str(entry.get("donor_count", 0))),
                metal_ids=metal_ids,
                reason=reason,
            )
        )
    return rows


def _metal_bound_het_summary(rows: Sequence[MetalBoundHet]) -> dict[str, object]:
    class_counts = Counter(row.classification for row in rows)
    expected_removed_count = _metal_bound_het_count(rows, _is_expected_removed_het)
    unexpected_removed_count = _metal_bound_het_count(rows, _is_unexpected_removed_het)
    retained_count = _metal_bound_het_count(rows, _is_retained_het)
    retained_competitive_count = _metal_bound_het_count(
        rows,
        _is_retained_competitive_ligand,
    )
    summary = {
        "metal_bound_het_status": _metal_bound_het_status(
            rows,
            unexpected_removed_count,
        ),
        "metal_bound_het_counts": _format_counter(class_counts),
        "metal_bound_het_retained_count": retained_count,
        "metal_bound_het_expected_removed_count": expected_removed_count,
        "metal_bound_het_unexpected_removed_count": unexpected_removed_count,
        "metal_bound_het_retained_donor_count": _metal_bound_het_donor_count(
            rows,
            _is_retained_het,
        ),
        "metal_bound_het_expected_removed_donor_count": _metal_bound_het_donor_count(
            rows,
            _is_expected_removed_het,
        ),
        "metal_bound_het_unexpected_removed_donor_count": _metal_bound_het_donor_count(
            rows,
            _is_unexpected_removed_het,
        ),
        "metal_bound_competitive_ligand_retained_count": retained_competitive_count,
        "metal_bound_competitive_ligand_retained_donor_count": _metal_bound_het_donor_count(
            rows,
            _is_retained_competitive_ligand,
        ),
    }
    return summary


def _is_expected_removed_het(row: MetalBoundHet) -> bool:
    return row.expected_removed and not row.retained_in_baseline


def _is_unexpected_removed_het(row: MetalBoundHet) -> bool:
    return not row.expected_removed and not row.retained_in_baseline


def _is_retained_het(row: MetalBoundHet) -> bool:
    return row.retained_in_baseline


def _is_retained_competitive_ligand(row: MetalBoundHet) -> bool:
    return row.classification == "competitive_ligand" and row.retained_in_baseline


def _metal_bound_het_count(
    rows: Sequence[MetalBoundHet],
    predicate: Callable[[MetalBoundHet], bool],
) -> int:
    return sum(1 for row in rows if predicate(row))


def _metal_bound_het_donor_count(
    rows: Sequence[MetalBoundHet],
    predicate: Callable[[MetalBoundHet], bool],
) -> int:
    return sum(row.donor_count for row in rows if predicate(row))


def _metal_bound_het_details(
    rows: Sequence[MetalBoundHet],
) -> list[dict[str, object]]:
    return [
        {
            "resname": row.residue_key[0],
            "chain": row.residue_key[1],
            "resseq": row.residue_key[2],
            "icode": row.residue_key[3],
            "classification": row.classification,
            "retained_in_baseline": row.retained_in_baseline,
            "expected_removed": row.expected_removed,
            "donor_count": row.donor_count,
            "metal_ids": list(row.metal_ids),
            "reason": row.reason,
        }
        for row in rows
    ]


def _metal_bound_het_status(
    rows: Sequence[MetalBoundHet],
    unexpected_removed_count: int,
) -> str:
    if unexpected_removed_count:
        return "unexpected_loss"
    if not rows:
        return "none"
    if any(row.expected_removed and not row.retained_in_baseline for row in rows):
        return "ok_expected_ligand_or_salt_loss"
    return "ok"


def _collect_chain_residues_and_caps(
    path: Path,
) -> tuple[dict[str, list[tuple[int, str]]], set[tuple[str, int, str]]]:
    chains: dict[str, list[tuple[int, str]]] = defaultdict(list)
    caps: set[tuple[str, int, str]] = set()
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            resname = line[17:20].strip().upper()
            chain = line[21:22].strip() or "-"
            resseq = _line_resseq_int(line)
            if resname in CAP_RESNAMES:
                caps.add((chain, resseq, resname))
                continue
            if not line.startswith("ATOM  "):
                continue
            if resname not in PROTEIN_RECORD_RESNAMES:
                continue
            if resseq == 0:
                continue
            pair = (resseq, resname)
            if pair not in chains[chain]:
                chains[chain].append(pair)
    return chains, caps


def audit_termini(path: Path) -> dict[str, object]:
    chains, caps = _collect_chain_residues_and_caps(path)
    termini = sum(2 for residues in chains.values() if residues)
    details, uncapped, capped_terms = _termini_details(chains, caps)
    status = _termini_status(termini, capped_terms)
    return {
        "termini_count": termini,
        "termini_capped_count": capped_terms,
        "termini_uncapped_count": max(0, termini - capped_terms),
        "termini_status": status,
        "termini_uncapped_sites": ";".join(uncapped[:20]),
        "termini_cap_plan": "no_action"
        if not uncapped
        else "audit_only_consider_explicit_caps_for_truncated_or_binding_site_termini",
        "termini_details": details,
    }


def _termini_details(
    chains: dict[str, list[tuple[int, str]]],
    caps: set[tuple[str, int, str]],
) -> tuple[list[dict[str, object]], list[str], int]:
    capped_terms = 0
    uncapped: list[str] = []
    details: list[dict[str, object]] = []
    for chain, residues in sorted(chains.items()):
        if not residues:
            continue
        detail = _terminus_detail(chain, residues, caps)
        capped_terms += int(bool(detail["n_capped"])) + int(bool(detail["c_capped"]))
        if not detail["n_capped"]:
            uncapped.append(f"{chain}:N:{detail['n_term']}")
        if not detail["c_capped"]:
            uncapped.append(f"{chain}:C:{detail['c_term']}")
        details.append(detail)
    return details, uncapped, capped_terms


def _terminus_detail(
    chain: str,
    residues: list[tuple[int, str]],
    caps: set[tuple[str, int, str]],
) -> dict[str, object]:
    ordered = sorted(residues)
    n_resseq, n_resname = ordered[0]
    c_resseq, c_resname = ordered[-1]
    n_capped = _has_cap(caps, chain, "ACE", n_resseq, before=True)
    c_capped = _has_cap(caps, chain, "NME", c_resseq, before=False)
    return {
        "chain": chain,
        "n_term": f"{n_resname}{n_resseq}",
        "c_term": f"{c_resname}{c_resseq}",
        "n_capped": n_capped,
        "c_capped": c_capped,
    }


def _has_cap(
    caps: set[tuple[str, int, str]],
    chain: str,
    cap_name: str,
    resseq: int,
    *,
    before: bool,
) -> bool:
    return any(
        cap_chain == chain
        and name == cap_name
        and (cap_resseq <= resseq if before else cap_resseq >= resseq)
        for cap_chain, cap_resseq, name in caps
    )


def _termini_status(termini: int, capped_terms: int) -> str:
    if termini == 0:
        return "no_protein_chains"
    if capped_terms == termini:
        return "fully_capped"
    if capped_terms:
        return "partially_capped"
    return "uncapped"


def run_probe_validation(path: Path) -> dict[str, object]:
    probe = shutil.which("probe")
    if not probe:
        return {"probe_status": "missing_probe", "probe_bad_contacts": ""}
    result = run_subprocess_capture(
        [probe, "-Quiet", "-Countdots", "-ONELINE", str(path)],
        timeout=120,
    )
    status = "ok" if result.returncode == 0 else "failed"
    output = (result.stdout or result.stderr or "").strip().splitlines()
    return {
        "probe_status": status,
        "probe_bad_contacts": output[-1][:300] if output else "",
    }


def _collect_atom_points(path: Path) -> list[AtomPoint]:
    atoms: list[AtomPoint] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")) or not _is_heavy_atom(line):
                continue
            xyz = _line_xyz(line)
            if xyz is None:
                continue
            atoms.append(
                AtomPoint(
                    serial=line[6:11].strip(),
                    resname=line[17:20].strip().upper(),
                    chain=line[21:22].strip() or "-",
                    resseq=line[22:26].strip(),
                    record=line[:6].strip(),
                    atom_name=line[12:16].strip(),
                    element=_line_element(line),
                    xyz=xyz,
                )
            )
    return atoms


def _grid_cell(atom: AtomPoint, cell_size: float) -> tuple[int, int, int]:
    return (
        math.floor(atom.xyz[0] / cell_size),
        math.floor(atom.xyz[1] / cell_size),
        math.floor(atom.xyz[2] / cell_size),
    )


def _neighbor_cells(cell: tuple[int, int, int]) -> list[tuple[int, int, int]]:
    return [
        (cell[0] + dx, cell[1] + dy, cell[2] + dz)
        for dx in (-1, 0, 1)
        for dy in (-1, 0, 1)
        for dz in (-1, 0, 1)
    ]


def write_protein_minimization_input(input_pdb: Path, output_pdb: Path) -> None:
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with input_pdb.open("r", encoding="utf-8", errors="ignore") as src:
        with output_pdb.open("w", encoding="utf-8") as dst:
            for line in src:
                if line.startswith("ATOM  ") and _line_element(line) != "H":
                    dst.write(line)
            dst.write("END\n")


def run_restrained_minimization_check(
    input_pdb: Path,
    output_pdb: Path,
    *,
    target_ph: float,
) -> str:
    work_input = output_pdb.with_suffix(".protein_only.pdb")
    write_protein_minimization_input(input_pdb, work_input)
    repaired_input = output_pdb.with_suffix(".protein_repaired.pdb")
    if repair_with_pdbfixer(
        work_input,
        repaired_input,
        target_ph=target_ph,
        cfg={"PDBFIXER_ADD_MISSING_RESIDUES": False},
        add_hydrogens=False,
    ):
        work_input = repaired_input
    try:
        from openmm import CustomExternalForce, LangevinIntegrator, Platform, unit
        from openmm.app import ForceField, Modeller, PDBFile, Simulation

        pdb = PDBFile(str(work_input))
        forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(forcefield, pH=float(target_ph))
        system = forcefield.createSystem(modeller.topology, constraints=None)
        restraint = CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        restraint.addGlobalParameter("k", 10.0)
        restraint.addPerParticleParameter("x0")
        restraint.addPerParticleParameter("y0")
        restraint.addPerParticleParameter("z0")
        positions = modeller.positions
        for index, atom in enumerate(modeller.topology.atoms()):
            if atom.element is not None and atom.element.symbol != "H":
                pos = positions[index].value_in_unit(unit.nanometer)
                restraint.addParticle(index, [pos.x, pos.y, pos.z])
        system.addForce(restraint)
        integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        platform = Platform.getPlatformByName("Reference")
        simulation = Simulation(modeller.topology, system, integrator, platform)
        simulation.context.setPositions(positions)
        simulation.minimizeEnergy(maxIterations=50)
        state = simulation.context.getState(getPositions=True)
        with output_pdb.open("w", encoding="utf-8") as handle:
            PDBFile.writeFile(
                simulation.topology,
                state.getPositions(asNumpy=True),
                handle,
                keepIds=True,
            )
        return "ok" if output_pdb.exists() and output_pdb.stat().st_size > 0 else "empty_output"
    except Exception as exc:
        return f"failed:{str(exc)[:180]}"


def try_constructive_geometry_repair(
    input_pdb: Path,
    target_dir: Path,
    *,
    target_ph: float,
    ligand_center: tuple[float, float, float] | None,
) -> tuple[Path | None, dict[str, object]]:
    """Use restrained minimization to repair clashing added sidechains."""

    minimized = target_dir / "baseline_constructive_geometry_repaired.pdb"
    status = run_restrained_minimization_check(
        input_pdb,
        minimized,
        target_ph=target_ph,
    )
    summary: dict[str, object] = {
        "geometry_fix_constructive_repair_status": status,
        "geometry_fix_constructive_repair_pdb": str(minimized) if minimized.exists() else "",
    }
    if status != "ok" or not minimized.exists():
        return None, summary
    rescue_retained_hets(
        source_pdb=input_pdb,
        target_pdb=minimized,
        stage="benchmark_constructive_geometry_repair",
        logger=logging.getLogger("protein_prep_benchmark"),
    )
    rescue_metal_bound_hets(
        source_pdb=input_pdb,
        target_pdb=minimized,
        stage="benchmark_constructive_geometry_repair_metal_bound",
        logger=logging.getLogger("protein_prep_benchmark"),
    )
    restore_metal_bound_waters(input_pdb, minimized)
    cleaned = target_dir / "baseline_constructive_geometry_repaired_clean.pdb"
    cleanup = apply_conservative_geometry_fixes(
        minimized,
        cleaned,
        reference_pdb=input_pdb,
        sidecar_path=target_dir / "geometry_cleanup_audit_constructive_repair.json",
        ligand_center=ligand_center,
    )
    candidate = _constructive_candidate(minimized, cleaned, cleanup)
    if candidate is None:
        summary["geometry_fix_constructive_repair_status"] = (
            "failed_constructive_cleanup:"
            + str(cleanup.get("geometry_fix_status", "unknown"))
        )
        return None, summary | _constructive_cleanup_summary(cleanup)
    audit = audit_geometry(
        candidate,
        sidecar_path=target_dir / "geometry_clash_audit_constructive_repair.json",
    )
    summary |= _constructive_cleanup_summary(cleanup)
    summary["geometry_fix_constructive_repair_pdb"] = str(candidate)
    if str(audit.get("geometry_status", "")) == "clashes":
        summary["geometry_fix_constructive_repair_status"] = "failed_geometry_clashes"
        return None, summary
    return candidate, summary


def _constructive_candidate(
    minimized: Path,
    cleaned: Path,
    cleanup: Mapping[str, object],
) -> Path | None:
    if _row_int(cleanup, "geometry_fix_dropped_sidechain_atom_count") > 0:
        return None
    if str(cleanup.get("geometry_fix_status", "")) == "unchanged":
        return minimized
    return cleaned if cleaned.exists() else None


def _constructive_cleanup_summary(cleanup: Mapping[str, object]) -> dict[str, object]:
    return {
        "geometry_fix_constructive_cleanup_status": str(
            cleanup.get("geometry_fix_status", "")
        ),
        "geometry_fix_constructive_removed_water_count": _row_int(
            cleanup,
            "geometry_fix_removed_water_count",
        ),
    }


def restore_metal_bound_waters(
    source_pdb: Path,
    target_pdb: Path,
    *,
    radius: float = 3.0,
) -> int:
    """Restore source waters that directly coordinate retained metals."""

    if not source_pdb.exists() or not target_pdb.exists():
        return 0
    metals, het_atoms = _collect_metal_bound_het_atoms(source_pdb)
    bound = _find_metal_bound_het_entries(metals, het_atoms, radius * radius)
    water_keys = {key for key in bound if key[0] in WATER_NAMES}
    if not water_keys:
        return 0
    existing = _water_residue_keys(target_pdb)
    missing = water_keys - existing
    if not missing:
        return 0
    lines = _source_residue_lines(source_pdb, missing)
    if not lines:
        return 0
    _append_records_before_end(target_pdb, lines)
    return len(missing)


def _water_residue_keys(path: Path) -> set[tuple[str, str, str, str]]:
    keys: set[tuple[str, str, str, str]] = set()
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("HETATM") and line[17:20].strip().upper() in WATER_NAMES:
                keys.add(_residue_key(line))
    return keys


def _source_residue_lines(
    source_pdb: Path,
    residue_keys: set[tuple[str, str, str, str]],
) -> list[str]:
    lines: list[str] = []
    with source_pdb.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("HETATM") and _residue_key(line) in residue_keys:
                lines.append(line)
    return lines


def _append_records_before_end(path: Path, records: Sequence[str]) -> None:
    original = path.read_text(encoding="utf-8", errors="ignore").splitlines(keepends=True)
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


def _short_process_error(result: object) -> str:
    text = ""
    if hasattr(result, "stderr") or hasattr(result, "stdout"):
        text = (getattr(result, "stderr", "") or getattr(result, "stdout", "") or "").strip()
    lines = [line for line in text.splitlines() if line.strip()]
    key_lines = [line for line in lines if "residue_key=" in line]
    selected = [*key_lines[-4:], *lines[-8:]]
    deduped = list(dict.fromkeys(selected))
    return "\n".join(deduped)[-1500:]


def _centroid(path: Path) -> tuple[float, float, float] | None:
    if not path.exists():
        return None
    coords = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")) and _is_heavy_atom(line):
                xyz = _line_xyz(line)
                if xyz:
                    coords.append(xyz)
    if not coords:
        return None
    n = float(len(coords))
    return (
        sum(point[0] for point in coords) / n,
        sum(point[1] for point in coords) / n,
        sum(point[2] for point in coords) / n,
    )


def _binding_site_clash_counts(
    sidecar_path: Path,
    ligand_center: tuple[float, float, float] | None,
    *,
    radius: float = 8.0,
) -> tuple[int, int]:
    try:
        audit = json.loads(sidecar_path.read_text(encoding="utf-8"))
    except Exception:
        return (0, 0)
    clashes = [
        row
        for row in audit.get("clashes", [])
        if isinstance(row, dict) and not row.get("ignored")
    ]
    if ligand_center is None:
        return (len(clashes), 0)
    radius2 = radius * radius
    binding_site = 0
    for row in clashes:
        atom_rows = [row.get("atom_a"), row.get("atom_b")]
        if any(
            isinstance(atom, dict)
            and _point_distance2(atom.get("xyz"), ligand_center) <= radius2
            for atom in atom_rows
        ):
            binding_site += 1
    return (binding_site, max(0, len(clashes) - binding_site))


def _point_distance2(
    raw_xyz: object,
    point: tuple[float, float, float],
) -> float:
    if not isinstance(raw_xyz, (list, tuple)) or len(raw_xyz) != 3:
        return math.inf
    try:
        xyz = (float(raw_xyz[0]), float(raw_xyz[1]), float(raw_xyz[2]))
    except Exception:
        return math.inf
    return (
        (xyz[0] - point[0]) ** 2
        + (xyz[1] - point[1]) ** 2
        + (xyz[2] - point[2]) ** 2
    )


def run_open_baseline(
    input_pdb: Path,
    work_dir: Path,
    target_ph: float,
    *,
    add_missing_residues: bool,
) -> tuple[Path, Path, bool, bool]:
    repaired = work_dir / "baseline_pdbfixer.pdb"
    protonated = work_dir / "baseline_protonated.pdb"
    repair_with_pdbfixer(
        input_pdb,
        repaired,
        target_ph=target_ph,
        cfg={
            "PDBFIXER_REPAIR": True,
            "PDBFIXER_ADD_MISSING_RESIDUES": add_missing_residues,
        },
        logger=logging.getLogger("protein_prep_benchmark"),
        add_hydrogens=False,
    )
    if not repaired.exists():
        shutil.copy2(input_pdb, repaired)
    rescue_retained_hets(
        source_pdb=input_pdb,
        target_pdb=repaired,
        stage="benchmark_pdbfixer",
        logger=logging.getLogger("protein_prep_benchmark"),
    )

    old_target_ph = os.environ.get("TARGET_PH")
    os.environ["TARGET_PH"] = f"{float(target_ph):.2f}"
    try:
        p2p_path, pdb2pqr_used, _ = _protonate_with_pdb2pqr_if_available(
            str(repaired),
            work_dir,
            variant="HOLO",
            logger=logging.getLogger("protein_prep_benchmark"),
        )
    finally:
        if old_target_ph is None:
            os.environ.pop("TARGET_PH", None)
        else:
            os.environ["TARGET_PH"] = old_target_ph
    reduce_exe = REDUCE_EXE if Path(str(REDUCE_EXE)).exists() else None
    try:
        assign_protonation_states(p2p_path, protonated, reduce_exe=reduce_exe)
        rescue_retained_hets(
            source_pdb=repaired,
            target_pdb=protonated,
            stage="benchmark_post_reduce",
            logger=logging.getLogger("protein_prep_benchmark"),
        )
        reduce_used = True
    except Exception as exc:
        logging.warning("[benchmark] reduce/fallback failed; copying baseline: %s", exc)
        shutil.copy2(p2p_path, protonated)
        reduce_used = False
    return protonated if protonated.exists() else Path(p2p_path), repaired, pdb2pqr_used, reduce_used


def run_meeko_with_repair_fallback(
    protonated_pdb: Path,
    repaired_pdb: Path,
    output_pdbqt: Path,
) -> tuple[bool, str, str, str]:
    ok, error = run_meeko_receptor(protonated_pdb, output_pdbqt)
    if ok:
        stage = "pdb2pqr_reduce" if not error else f"pdb2pqr_reduce:{error}"
        return True, error, stage, ""
    primary_error = error
    if repaired_pdb == protonated_pdb or not repaired_pdb.exists():
        return False, error, "pdb2pqr_reduce", primary_error
    ok, error = run_meeko_receptor(repaired_pdb, output_pdbqt)
    stage = "pdbfixer_repair_fallback" if ok else "pdb2pqr_reduce"
    return ok, error, stage, primary_error


def _run_meeko_attempts(
    attempts: Sequence[Sequence[str]], output_pdbqt: Path
) -> tuple[bool, str]:
    last_error = ""
    for cmd in attempts:
        try:
            result = run_subprocess_capture(cmd, timeout=180)
        except Exception as exc:
            last_error = last_error or str(exc)
            continue
        if result.returncode == 0 and output_pdbqt.exists() and output_pdbqt.stat().st_size > 0:
            return True, ""
        last_error = last_error or _short_process_error(result)
    return False, last_error


def _meeko_cmds_for(
    base_cmd: Sequence[str],
    path: Path,
    output_pdbqt: Path,
) -> tuple[list[str], list[str]]:
    return (
        build_meeko_modern_cmd(base_cmd, path, output_pdbqt),
        build_meeko_legacy_cmd(base_cmd, path, output_pdbqt),
    )


def _try_meeko_prepared_input(
    input_pdb: Path,
    prepared_pdb: Path,
    output_pdbqt: Path,
    base_cmd: Sequence[str],
    label: str,
    writer: Callable[[Path, Path], object],
) -> tuple[bool, str, str]:
    writer(input_pdb, prepared_pdb)
    ok, error = _run_meeko_attempts(
        _meeko_cmds_for(base_cmd, prepared_pdb, output_pdbqt),
        output_pdbqt,
    )
    return ok, label if ok else "", error


def _try_meeko_clash_pruned_input(
    input_pdb: Path,
    output_pdbqt: Path,
    base_cmd: Sequence[str],
) -> tuple[bool, str, str]:
    pruned = output_pdbqt.with_suffix(".meeko_clash_pruned_input.pdb")
    dropped = write_meeko_clash_pruned_copy(input_pdb, pruned)
    if not dropped:
        return False, "", ""
    clean_attempts = _meeko_cmds_for(base_cmd, pruned, output_pdbqt)
    ok, error = _run_meeko_attempts(clean_attempts, output_pdbqt)
    label = meeko_pruned_stage_label(pruned, dropped=dropped) if ok else ""
    if ok or os.environ.get("MEEKO_ALLOW_BAD_RES", "").lower() not in {"1", "true", "yes"}:
        return ok, label, error
    ok, error = _run_meeko_attempts(tuple(cmd + ["-a"] for cmd in clean_attempts), output_pdbqt)
    label = (
        meeko_pruned_stage_label(pruned, dropped=dropped, allow_bad=True) if ok else ""
    )
    return ok, label, error


def _try_meeko_safe_stage(
    input_pdb: Path,
    output_pdbqt: Path,
    base_cmd: Sequence[str],
    *,
    suffix: str,
    residue_keys: Iterable[tuple[str, str]] = (),
) -> tuple[Path, dict[str, object], bool, str]:
    safe = output_pdbqt.with_suffix(suffix)
    summary = write_meeko_safe_receptor_copy(
        input_pdb,
        safe,
        drop_residue_keys=residue_keys,
    )
    ok, error = _run_meeko_attempts(
        _meeko_cmds_for(base_cmd, safe, output_pdbqt),
        output_pdbqt,
    )
    return safe, summary, ok, error


def _try_meeko_safe_clash_stage(
    safe_input: Path,
    output_pdbqt: Path,
    base_cmd: Sequence[str],
    summary: dict[str, object],
    *,
    suffix: str,
) -> tuple[bool, str, str]:
    pruned = output_pdbqt.with_suffix(suffix)
    dropped = write_meeko_clash_pruned_copy(safe_input, pruned)
    if not dropped:
        return False, "", ""
    ok, error = _run_meeko_attempts(
        _meeko_cmds_for(base_cmd, pruned, output_pdbqt),
        output_pdbqt,
    )
    if not ok:
        return False, "", error
    label = f"{meeko_safe_stage_label(summary)}:clash_pruned_retry:dropped_atoms={dropped}"
    return True, label, ""


def _try_meeko_template_pruned_inputs(
    input_pdb: Path,
    output_pdbqt: Path,
    base_cmd: Sequence[str],
    residue_keys: set[tuple[str, str]],
    last_error: str,
) -> tuple[bool, str, str]:
    for attempt in range(1, 33):
        if not residue_keys:
            break
        safe, summary, ok, last_error = _try_meeko_safe_stage(
            input_pdb,
            output_pdbqt,
            base_cmd,
            suffix=f".meeko_template_pruned_input{attempt}.pdb",
            residue_keys=residue_keys,
        )
        if ok:
            return True, meeko_safe_stage_label(summary), ""
        ok, label, clash_error = _try_meeko_safe_clash_stage(
            safe,
            output_pdbqt,
            base_cmd,
            summary,
            suffix=f".meeko_template_pruned_clash_input{attempt}.pdb",
        )
        if ok:
            return True, label, ""
        last_error = clash_error or last_error
        next_keys = parse_unmatched_residue_keys(last_error)
        if next_keys.issubset(residue_keys):
            break
        residue_keys |= next_keys
    return False, "", last_error


def _try_meeko_safe_protein_input(
    input_pdb: Path,
    output_pdbqt: Path,
    base_cmd: Sequence[str],
    *,
    error_context: str,
) -> tuple[bool, str, str]:
    safe, summary, ok, error = _try_meeko_safe_stage(
        input_pdb,
        output_pdbqt,
        base_cmd,
        suffix=".meeko_protein_only_input.pdb",
    )
    if ok:
        return True, meeko_safe_stage_label(summary), ""
    ok, label, pruned_error = _try_meeko_safe_clash_stage(
        safe,
        output_pdbqt,
        base_cmd,
        summary,
        suffix=".meeko_protein_only_clash_pruned_input.pdb",
    )
    if ok:
        return True, label, ""
    error = pruned_error or error
    residue_keys = parse_unmatched_residue_keys(f"{error_context}\n{error}")
    return _try_meeko_template_pruned_inputs(
        input_pdb,
        output_pdbqt,
        base_cmd,
        residue_keys,
        error,
    )


def run_meeko_receptor(input_pdb: Path, output_pdbqt: Path) -> tuple[bool, str]:
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    try:
        base_cmd = build_meeko_base_cmd()
    except Exception as exc:
        return False, f"meeko_unavailable:{exc}"

    ok, last_error = _run_meeko_attempts(
        _meeko_cmds_for(base_cmd, input_pdb, output_pdbqt),
        output_pdbqt,
    )
    if ok:
        return True, ""

    retries: tuple[tuple[Path, str, Callable[[Path, Path], object]], ...] = (
        (
            output_pdbqt.with_suffix(".meeko_ordered_input.pdb"),
            "ordered_copy_retry",
            write_meeko_ordered_copy,
        ),
        (
            output_pdbqt.with_suffix(".meeko_normalized_input.pdb"),
            "normalized_dewatered_noh_retry",
            write_meeko_normalized_copy,
        ),
    )
    for prepared, label, writer in retries:
        ok, retry_label, retry_error = _try_meeko_prepared_input(
            input_pdb,
            prepared,
            output_pdbqt,
            base_cmd,
            label,
            writer,
        )
        if ok:
            return True, retry_label
        last_error = retry_error or last_error

    ok, label, pruned_error = _try_meeko_clash_pruned_input(
        input_pdb,
        output_pdbqt,
        base_cmd,
    )
    if ok:
        return True, label
    last_error = pruned_error or last_error

    ok, label, safe_error = _try_meeko_safe_protein_input(
        input_pdb,
        output_pdbqt,
        base_cmd,
        error_context=last_error,
    )
    if ok:
        return True, label
    last_error = safe_error or last_error
    return False, last_error


def _write_meeko_stage_sidecar(
    path: Path,
    *,
    ok: bool,
    stage: str,
    error: str,
) -> None:
    stage_label = stage or "direct"
    path.write_text(
        json.dumps(
            {
                "ok": ok,
                "stage": stage_label,
                "error": error,
                "review_required": meeko_stage_requires_review(stage_label),
                "publication_review_required": meeko_stage_requires_review(
                    stage_label,
                    profile="publication",
                ),
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )


def _meeko_ligand_cmd() -> list[str] | None:
    try:
        import meeko  # type: ignore[import-untyped]  # noqa: F401

        return [sys.executable, "-m", "meeko.cli.mk_prepare_ligand"]
    except Exception:
        exe = shutil.which("mk_prepare_ligand.py") or shutil.which("mk_prepare_ligand")
        return [exe] if exe else None


def _redock_prereq_status(ligand_pdb: Path, receptor_pdbqt: Path) -> tuple[str, str | None]:
    if not receptor_pdbqt.exists():
        return "no_receptor_pdbqt", None
    if not ligand_pdb.exists():
        return "no_crystal_ligand", None
    if not shutil.which("obabel"):
        return "missing_obabel", None
    vina = shutil.which("vina")
    if not vina:
        return "missing_vina", None
    if _meeko_ligand_cmd() is None:
        return "missing_meeko_ligand_cli", None
    return "ok", vina


def _prepare_redock_ligand(
    ligand_pdb: Path,
    out_dir: Path,
    ligand_resname: str = "",
) -> tuple[Path | None, str]:
    ligand_sdf = out_dir / "crystal_ligand.sdf"
    ligand_pdbqt = out_dir / "crystal_ligand.pdbqt"
    sdf_status = _write_redock_ligand_sdf(ligand_pdb, ligand_sdf, out_dir, ligand_resname)
    if not sdf_status.startswith("ok"):
        return None, sdf_status

    ligand_cmd = _meeko_ligand_cmd()
    if ligand_cmd is None:
        return None, "missing_meeko_ligand_cli"
    ligprep = run_subprocess_capture(
        list(ligand_cmd) + ["-i", str(ligand_sdf), "-o", str(ligand_pdbqt)],
        timeout=180,
    )
    if ligprep.returncode != 0 or not ligand_pdbqt.exists():
        return None, f"meeko_ligand_failed:{sdf_status}"
    return ligand_pdbqt, sdf_status


def _write_redock_ligand_sdf(
    ligand_pdb: Path,
    ligand_sdf: Path,
    out_dir: Path,
    ligand_resname: str,
) -> str:
    ccd_sdf = out_dir / "crystal_ligand_ccd_instance.sdf"
    selected_smiles = _selected_het_smiles(out_dir, ligand_resname)
    ccd_status = write_ccd_instance_sdf(
        ligand_pdb,
        ligand_resname,
        out_dir,
        ccd_sdf,
        selected_smiles=selected_smiles,
    )
    if ccd_status.get("ccd_instance_status") == "ok":
        shutil.copyfile(ccd_sdf, ligand_sdf)
        _write_redock_ligand_audit(
            out_dir,
            ccd_status | {"selected_het_smiles_used": bool(selected_smiles)},
        )
        return "ok_ccd_instance"

    obabel = run_subprocess_capture(
        ["obabel", str(ligand_pdb), "-O", str(ligand_sdf), "-h"],
        timeout=120,
    )
    if obabel.returncode != 0 or not ligand_sdf.exists():
        _write_redock_ligand_audit(out_dir, ccd_status | {"fallback_status": "obabel_failed"})
        return "obabel_failed"
    _write_redock_ligand_audit(out_dir, ccd_status | {"fallback_status": "ok_obabel"})
    return "ok_obabel"


def _selected_het_smiles(out_dir: Path, ligand_resname: str) -> str:
    token = ligand_resname.strip().upper()
    if not token:
        return ""
    for base in (out_dir, *out_dir.parents):
        sidecar = base / "het_state_selection.json"
        if not sidecar.exists():
            continue
        try:
            data = json.loads(sidecar.read_text(encoding="utf-8"))
        except Exception:
            continue
        for row in data.get("selections", []):
            if not isinstance(row, Mapping):
                continue
            if str(row.get("resname", "")).strip().upper() == token:
                if str(row.get("confidence", "")).strip().lower() == "low":
                    return ""
                return str(row.get("selected_smiles", "") or "").strip()
    return ""


def _write_redock_ligand_audit(out_dir: Path, payload: Mapping[str, object]) -> None:
    (out_dir / "redock_ligand_prep.json").write_text(
        json.dumps(dict(payload), indent=2, sort_keys=True),
        encoding="utf-8",
    )


def _run_vina_redock(
    vina: str,
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    docked: Path,
    center: tuple[float, float, float],
    *,
    box_size: float,
    exhaustiveness: int,
    num_modes: int = 9,
    energy_range: float = 6.0,
    seed: int = 1,
) -> str:
    cmd = [
        vina,
        "--receptor",
        str(receptor_pdbqt),
        "--ligand",
        str(ligand_pdbqt),
        "--center_x",
        f"{center[0]:.3f}",
        "--center_y",
        f"{center[1]:.3f}",
        "--center_z",
        f"{center[2]:.3f}",
        "--size_x",
        str(box_size),
        "--size_y",
        str(box_size),
        "--size_z",
        str(box_size),
        "--exhaustiveness",
        str(exhaustiveness),
        "--num_modes",
        str(num_modes),
        "--energy_range",
        str(energy_range),
        "--cpu",
        str(_benchmark_cpu_count()),
        "--seed",
        str(seed),
        "--out",
        str(docked),
    ]
    dock = run_subprocess_capture(cmd, timeout=900)
    return "ok" if dock.returncode == 0 and docked.exists() else "vina_failed"


def _benchmark_cpu_count() -> int:
    return max(1, min(32, os.cpu_count() or 1))


def run_optional_redock(
    ligand_pdb: Path,
    receptor_pdbqt: Path,
    out_dir: Path,
    *,
    ligand_resname: str,
    box_size: float,
    exhaustiveness: int,
) -> tuple[float | None, str, str]:
    status, vina = _redock_prereq_status(ligand_pdb, receptor_pdbqt)
    if status != "ok" or vina is None:
        return None, status, "not_run"
    docked = out_dir / "redock_out.pdbqt"
    center = _centroid(ligand_pdb)
    if center is None:
        return None, "no_ligand_centroid", "not_run"

    ligand_pdbqt, ligand_prep_status = _prepare_redock_ligand(
        ligand_pdb,
        out_dir,
        ligand_resname,
    )
    if ligand_pdbqt is None:
        return None, ligand_prep_status, ligand_prep_status
    status = _run_vina_redock(
        vina,
        receptor_pdbqt,
        ligand_pdbqt,
        docked,
        center,
        box_size=box_size,
        exhaustiveness=exhaustiveness,
    )
    if status != "ok":
        return None, status, ligand_prep_status

    try:
        from docking.pose_validation import compute_redock_rmsd

        rmsd = _compute_prepared_ligand_redock_rmsd(
            ligand_pdb,
            ligand_pdbqt,
            docked,
            ligand_prep_status,
        )
        if rmsd is None:
            rmsd = compute_redock_rmsd(str(ligand_pdb), str(docked))
        if rmsd is None or rmsd > 2.5:
            high_rmsd, high_status = _run_high_search_redock(
                vina,
                receptor_pdbqt,
                ligand_pdbqt,
                ligand_pdb,
                out_dir,
                center,
                ligand_prep_status,
                box_size=box_size,
                compute_fallback_rmsd=compute_redock_rmsd,
            )
            if high_rmsd is not None and (rmsd is None or high_rmsd < rmsd):
                return high_rmsd, high_status, ligand_prep_status
        return rmsd, "ok", ligand_prep_status
    except Exception as exc:
        return None, f"rmsd_failed:{exc}", ligand_prep_status


def _run_high_search_redock(
    vina: str,
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    ligand_pdb: Path,
    out_dir: Path,
    center: tuple[float, float, float],
    ligand_prep_status: str,
    *,
    box_size: float,
    compute_fallback_rmsd: object,
) -> tuple[float | None, str]:
    high_docked = out_dir / "redock_out_high_search.pdbqt"
    status = _run_vina_redock(
        vina,
        receptor_pdbqt,
        ligand_pdbqt,
        high_docked,
        center,
        box_size=box_size,
        exhaustiveness=64,
        num_modes=20,
        energy_range=9.0,
        seed=1337,
    )
    if status != "ok":
        return None, status
    rmsd = _compute_best_prepared_ligand_redock_rmsd(
        ligand_pdbqt,
        high_docked,
        ligand_prep_status,
    )
    if rmsd is None and callable(compute_fallback_rmsd):
        rmsd = compute_fallback_rmsd(str(ligand_pdb), str(high_docked))
    return rmsd, "ok_high_search_best_pose"


def _water_redock_variants(
    *,
    receptor_pdbqt: Path,
    selected_water_pdbqt: Path,
    selected_water_meeko_ok: bool,
    remove_all_pdbqt: Path,
    remove_all_meeko_ok: bool,
) -> tuple[tuple[str, Path], ...]:
    variants: list[tuple[str, Path]] = [("baseline", receptor_pdbqt)]
    if selected_water_meeko_ok:
        variants.append(("selected", selected_water_pdbqt))
    if remove_all_meeko_ok:
        variants.append(("remove_all", remove_all_pdbqt))
    return tuple(variants)


def _compute_prepared_ligand_redock_rmsd(
    ligand_pdb: Path,
    ligand_pdbqt: Path,
    docked_pdbqt: Path,
    ligand_prep_status: str,
) -> float | None:
    if ligand_prep_status != "ok_ccd_instance":
        return None
    del ligand_pdb
    ref = _pdbqt_heavy_coords(ligand_pdbqt)
    pose = _pdbqt_heavy_coords(docked_pdbqt)
    if len(ref) < 5 or len(ref) != len(pose):
        return None
    ordered_rmsd = _kabsch_rmsd(ref, pose)
    symmetry_rmsd = _compute_rdkit_pose_rmsd(
        ligand_pdbqt.with_suffix(".sdf"),
        docked_pdbqt,
    )
    if symmetry_rmsd is None:
        return ordered_rmsd
    return min(ordered_rmsd, symmetry_rmsd)


def _compute_best_prepared_ligand_redock_rmsd(
    ligand_pdbqt: Path,
    docked_pdbqt: Path,
    ligand_prep_status: str,
) -> float | None:
    if ligand_prep_status != "ok_ccd_instance":
        return None
    reference_sdf = ligand_pdbqt.with_suffix(".sdf")
    ref = _pdbqt_heavy_coords(ligand_pdbqt)
    if len(ref) < 5:
        return None
    ref_mol = _load_redock_reference_mol(reference_sdf)
    if ref_mol is not None and ref_mol.GetNumAtoms() < 5:
        ref_mol = None
    scores: list[float] = []
    for pose in _pdbqt_heavy_coords_by_model(docked_pdbqt):
        if len(pose) != len(ref):
            continue
        ordered = _kabsch_rmsd(ref, pose)
        symmetry = (
            _rdkit_best_rms_for_pose_coords(ref_mol, pose) if ref_mol is not None else None
        )
        scores.append(min(ordered, symmetry) if symmetry is not None else ordered)
    return min(scores) if scores else None


def _compute_rdkit_pose_rmsd(reference_sdf: Path, docked_pdbqt: Path) -> float | None:
    pose_coords = _pdbqt_heavy_coords(docked_pdbqt)
    return _compute_rdkit_pose_rmsd_from_coords(reference_sdf, pose_coords)


def _load_redock_reference_mol(reference_sdf: Path) -> Any:
    try:
        from rdkit import Chem

        return next(
            (
                mol
                for mol in Chem.SDMolSupplier(
                    str(reference_sdf),
                    sanitize=True,
                    removeHs=True,
                )
                if mol is not None
            ),
            None,
        )
    except Exception:
        return None


def _rdkit_best_rms_for_pose_coords(
    reference_mol: Any,
    pose_coords: Sequence[tuple[float, float, float]],
) -> float | None:
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolAlign
        from rdkit.Geometry import Point3D

        ref = reference_mol
        if ref is None or ref.GetNumAtoms() < 5:
            return None
        if len(pose_coords) != ref.GetNumAtoms():
            return None
        pose = Chem.Mol(ref)
        pose.RemoveAllConformers()
        conf = Chem.Conformer(pose.GetNumAtoms())
        for idx, xyz in enumerate(pose_coords):
            conf.SetAtomPosition(idx, Point3D(float(xyz[0]), float(xyz[1]), float(xyz[2])))
        pose.AddConformer(conf, assignId=True)
        return float(rdMolAlign.GetBestRMS(ref, pose))
    except Exception:
        return None


def _compute_rdkit_pose_rmsd_from_coords(
    reference_sdf: Path,
    pose_coords: Sequence[tuple[float, float, float]],
) -> float | None:
    ref_mol = _load_redock_reference_mol(reference_sdf)
    if ref_mol is None or ref_mol.GetNumAtoms() < 5:
        return None
    return _rdkit_best_rms_for_pose_coords(ref_mol, pose_coords)


def _pdbqt_heavy_coords(path: Path) -> list[tuple[float, float, float]]:
    coords: list[tuple[float, float, float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("ENDMDL") and coords:
                break
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            element = _pdbqt_element(line)
            if element == "H":
                continue
            xyz = _line_xyz(line)
            if xyz is not None:
                coords.append(xyz)
    return coords


def _pdbqt_heavy_coords_by_model(path: Path) -> list[list[tuple[float, float, float]]]:
    models: list[list[tuple[float, float, float]]] = []
    current: list[tuple[float, float, float]] = []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith("MODEL") and current:
                models.append(current)
                current = []
                continue
            if line.startswith("ENDMDL"):
                if current:
                    models.append(current)
                    current = []
                continue
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if _pdbqt_element(line) == "H":
                continue
            xyz = _line_xyz(line)
            if xyz is not None:
                current.append(xyz)
    if current:
        models.append(current)
    return models


def _pdbqt_element(line: str) -> str:
    raw = line[77:79].strip() if len(line) >= 79 else ""
    if not raw:
        raw = line[12:16].strip()
    return raw[:1].upper()


def _kabsch_rmsd(
    ref: Sequence[tuple[float, float, float]],
    pose: Sequence[tuple[float, float, float]],
) -> float:
    import numpy as np

    ref_arr = np.asarray(ref, dtype=float)
    pose_arr = np.asarray(pose, dtype=float)
    ref_centered = ref_arr - ref_arr.mean(axis=0)
    pose_centered = pose_arr - pose_arr.mean(axis=0)
    covariance = ref_centered.T @ pose_centered
    left, _, right_t = np.linalg.svd(covariance)
    sign = np.sign(np.linalg.det(right_t.T @ left.T))
    correction = np.diag([1.0, 1.0, sign])
    rotated = ref_centered @ (left @ correction @ right_t)
    diff = rotated - pose_centered
    return float(np.sqrt((diff * diff).sum() / ref_arr.shape[0]))


def _finite_or_none(value: float | None) -> float | None:
    if value is None or not math.isfinite(float(value)):
        return None
    return float(value)


def _row_int(row: Mapping[str, object], key: str) -> int:
    raw = row.get(key, 0)
    try:
        return int(float(str(raw or 0)))
    except Exception:
        return 0


def _row_float(row: BenchmarkRow, key: str) -> float | None:
    raw = row.get(key, "")
    try:
        value = float(str(raw))
    except Exception:
        return None
    return value if math.isfinite(value) else None


def classify_validity(row: BenchmarkRow) -> tuple[bool, str]:
    failures: list[str] = []
    if _row_int(row, "baseline_heavy_atoms") <= 0:
        failures.append("no_baseline_heavy_atoms")
    if _row_int(row, "baseline_hydrogens") <= 0:
        failures.append("no_baseline_hydrogens")
    if not bool(row.get("meeko_ok")):
        failures.append("meeko_receptor_failed")
    if str(row.get("het_chemistry_status", "")) not in {"", "no_ligand", "ok"}:
        failures.append("het_chemistry_failed")
    if str(row.get("het_state_status", "")) in {"failed"}:
        failures.append("het_state_enumeration_failed")
    if _row_int(row, "metals_baseline") < _row_int(row, "metals_input"):
        failures.append("metal_count_decreased")
    if _row_int(row, "metal_bound_het_unexpected_removed_count") > 0:
        failures.append("unexpected_metal_bound_het_removed")
    if str(row.get("probe_status", "")) == "failed":
        failures.append("probe_failed")
    if (
        str(row.get("baseline_geometry_status", "")) == "clashes"
        and _row_int(row, "baseline_geometry_binding_site_clash_count") > 0
    ):
        failures.append("geometry_clashes")
    return not failures, ";".join(failures)


def classify_publication_readiness(
    row: BenchmarkRow,
    *,
    require_redock: bool,
    redock_rmsd_max: float,
    require_minimization: bool,
) -> dict[str, object]:
    """Strict publication gate plus disclosure-only limitations."""

    failures: list[str] = []
    disclosures: list[str] = []
    _add_publication_structure_gates(row, failures, disclosures)
    _add_publication_geometry_gates(row, failures, disclosures)
    _add_publication_water_gates(row, failures, disclosures)
    _add_publication_het_gates(row, failures, disclosures)
    _add_publication_hbond_disclosure(row, failures, disclosures)
    _add_publication_metal_gates(row, failures, disclosures)
    _add_publication_termini_disclosure(row, disclosures)
    _add_publication_redock_gate(
        row,
        failures,
        disclosures,
        require_redock=require_redock,
        redock_rmsd_max=redock_rmsd_max,
    )
    _add_publication_minimization_gate(
        row,
        failures,
        disclosures,
        require_minimization=require_minimization,
    )
    return {
        "publication_ready": not failures,
        "publication_readiness_failures": ";".join(failures),
        "publication_readiness_gate_count": len(failures),
        "publication_readiness_disclosures": ";".join(disclosures),
    }


def _add_publication_structure_gates(
    row: BenchmarkRow,
    failures: list[str],
    disclosures: list[str],
) -> None:
    if not bool(row.get("structure_valid")):
        failures.append(
            "structure_validity_failed:"
            + str(row.get("validity_failures", "") or "unknown")
        )
    stage = str(row.get("meeko_input_stage", "") or "")
    if meeko_stage_requires_review(stage, profile="publication"):
        failures.append(f"meeko_parser_recovery_review:{stage}")
    elif "normalized" in stage:
        disclosures.append(f"meeko_parser_normalized:{stage}")


def _add_publication_geometry_gates(
    row: BenchmarkRow,
    failures: list[str],
    disclosures: list[str],
) -> None:
    if _row_int(row, "baseline_geometry_binding_site_clash_count") > 0:
        failures.append("binding_site_geometry_clashes")
    elif _row_int(row, "baseline_geometry_nonbinding_clash_count") > 0:
        disclosures.append("remote_geometry_clashes_disclosed")
    if bool(row.get("geometry_fix_requires_constructive_repair")):
        failures.append("binding_site_sidechain_constructive_repair_needed")
    elif str(row.get("geometry_fix_status", "")) != "unchanged":
        disclosures.append("geometry_cleanup_applied")


def _add_publication_water_gates(
    row: BenchmarkRow,
    failures: list[str],
    disclosures: list[str],
) -> None:
    if str(row.get("water_policy_status", "")) == "review":
        if (
            str(row.get("water_policy_selection_status", "")) == "selected"
            and bool(row.get("water_policy_selected_meeko_ok"))
        ):
            if bool(row.get("water_redock_selected_supported")):
                if bool(row.get("water_evidence_supported")):
                    disclosures.append("water_policy_selected_by_density_conservation_and_redocking")
                else:
                    failures.append("water_density_conservation_evidence_missing")
            else:
                failures.append(
                    "water_redock_sensitivity_review:"
                    + str(row.get("water_redock_sensitivity_status", "unknown"))
                )
        else:
            failures.append("water_policy_review")
    if bool(row.get("water_keep_all_meeko_ok")) != bool(
        row.get("water_remove_all_meeko_ok")
    ):
        failures.append("water_meeko_sensitivity_review")


def _add_publication_het_gates(
    row: BenchmarkRow,
    failures: list[str],
    disclosures: list[str],
) -> None:
    het_status = str(row.get("het_state_status", ""))
    selection_status = str(row.get("het_state_selection_status", ""))
    if het_status in {"partial", "failed"}:
        failures.append("het_state_selection_review")
    elif selection_status in {"partial", "failed"}:
        failures.append(f"het_state_context_selection_{selection_status}")
    elif _row_int(row, "het_state_variant_count") > _row_int(
        row,
        "het_state_resname_count",
    ) and selection_status != "selected":
        failures.append("het_multiple_states_not_contextually_selected")
    elif selection_status == "selected":
        disclosures.append("het_state_selected_by_binding_site_network_score")
        labels = str(row.get("het_state_selection_review_labels", ""))
        if labels:
            disclosures.append(f"het_state_selection_review_labels:{labels}")
        if _row_int(row, "het_state_selection_low_confidence_count") > 0:
            disclosures.append("het_state_selection_low_confidence")


def _add_publication_hbond_disclosure(
    row: BenchmarkRow,
    failures: list[str],
    disclosures: list[str],
) -> None:
    if str(row.get("hbond_network_status", "")) != "reduce_probe_available":
        failures.append("hbond_network_tooling_partial")
    else:
        disclosures.append("hbond_network_local_not_global")
    if _row_int(row, "network_application_applied_count") > 0:
        disclosures.append("binding_site_network_local_modes_applied")
    if _row_int(row, "water_network_application_water_hydrogen_added_count") > 0:
        disclosures.append("selected_water_hydrogens_oriented_by_local_network")


def _add_publication_metal_gates(
    row: BenchmarkRow,
    failures: list[str],
    disclosures: list[str],
) -> None:
    failure = str(row.get("metal_publication_failure", ""))
    if failure:
        failures.append(failure)
    disclosure = str(row.get("metal_publication_disclosure", ""))
    if disclosure:
        disclosures.append(disclosure)


def _add_publication_termini_disclosure(
    row: BenchmarkRow,
    disclosures: list[str],
) -> None:
    if _row_int(row, "termini_uncapped_count") > 0:
        disclosures.append("termini_uncapped_default")


def _add_publication_redock_gate(
    row: BenchmarkRow,
    failures: list[str],
    disclosures: list[str],
    *,
    require_redock: bool,
    redock_rmsd_max: float,
) -> None:
    if require_redock:
        rmsd = _row_float(row, "redock_rmsd_a")
        status = str(row.get("redock_status", ""))
        if status != "ok" or rmsd is None:
            failures.append(f"redock_missing_or_failed:{status or 'unknown'}")
        elif rmsd > redock_rmsd_max:
            failures.append(f"redock_rmsd_gt_{redock_rmsd_max:.2f}A")
    else:
        disclosures.append("redocking_not_required_for_default_gate")


def _add_publication_minimization_gate(
    row: BenchmarkRow,
    failures: list[str],
    disclosures: list[str],
    *,
    require_minimization: bool,
) -> None:
    if require_minimization and str(row.get("geometry_fix_status", "")) != "unchanged":
        status = str(row.get("restrained_minimization_status", ""))
        if status != "ok":
            failures.append(f"restrained_minimization_missing_after_cleanup:{status}")
    elif not require_minimization:
        disclosures.append("restrained_minimization_optional")


def _out_of_scope_parameterization_plan() -> dict[str, object]:
    return {
        "status": "out_of_scope_default_docking",
        "site_count": 0,
        "missing_required_tools": [],
        "optional_missing_tools": [],
        "reason": (
            "Amber/MCPB metal-center parameterization is reserved for the "
            "MMGBSA/MD workflow and is not a default docking-prep gate."
        ),
        "mcpb_preclean": {"status": "not_run_default_docking"},
        "mcpb_cofactor_parameterization": {
            "status": "not_run_default_docking",
            "cofactor_residue_count": 0,
            "parameterized_count": 0,
            "unresolved_count": 0,
            "entries": [],
        },
    }


def run_one(
    pdb_id: str,
    out_dir: Path,
    *,
    add_missing_residues: bool,
    run_minimize: bool,
    run_redock: bool,
    run_autodock4zn: bool,
    run_water_analysis: bool,
    publication_mode: bool,
    publication_redock_rmsd_max: float,
    box_size: float,
    exhaustiveness: int,
) -> BenchmarkRow:
    target_dir = out_dir / pdb_id.upper()
    raw_pdb = target_dir / f"{pdb_id.upper()}.pdb"
    target_dir.mkdir(parents=True, exist_ok=True)
    if not download_pdb(pdb_id, raw_pdb):
        raise RuntimeError(f"download_failed:{pdb_id}")

    ligand_pdb = target_dir / "crystal_ligand.pdb"
    ligand_resname, ligand_heavy_atoms = extract_largest_ligand(raw_pdb, ligand_pdb)
    ligand_center = _centroid(ligand_pdb)
    target_ph = get_target_ph_for_prep(raw_pdb)
    baseline_pdb, repaired_pdb, pdb2pqr_used, reduce_used = run_open_baseline(
        raw_pdb,
        target_dir,
        target_ph,
        add_missing_residues=add_missing_residues,
    )
    pre_fix_geometry = audit_geometry(
        baseline_pdb,
        sidecar_path=target_dir / "geometry_clash_audit_pre_fix.json",
    )
    pre_geometry_cleanup_pdb = baseline_pdb
    geometry_fix = apply_conservative_geometry_fixes(
        baseline_pdb,
        target_dir / "baseline_geometry_fixed.pdb",
        reference_pdb=raw_pdb,
        sidecar_path=target_dir / "geometry_cleanup_audit.json",
        ligand_center=ligand_center,
    )
    if _row_int(geometry_fix, "geometry_fix_dropped_sidechain_atom_count") > 0:
        constructive_pdb, constructive_summary = try_constructive_geometry_repair(
            pre_geometry_cleanup_pdb,
            target_dir,
            target_ph=target_ph,
            ligand_center=ligand_center,
        )
        geometry_fix |= constructive_summary
        if constructive_pdb is not None:
            geometry_fix |= {
                "geometry_fix_status": (
                    "restrained_constructive_repair:"
                    + str(constructive_summary.get("geometry_fix_constructive_cleanup_status", ""))
                ),
                "geometry_fix_original_dropped_sidechain_atom_count": _row_int(
                    geometry_fix,
                    "geometry_fix_dropped_sidechain_atom_count",
                ),
                "geometry_fix_dropped_sidechain_atom_count": 0,
                "geometry_fix_requires_constructive_repair": False,
                "geometry_fix_output_pdb": str(constructive_pdb),
            }
    if str(geometry_fix.get("geometry_fix_status", "")) != "unchanged":
        baseline_pdb = Path(str(geometry_fix["geometry_fix_output_pdb"]))

    input_stats = collect_structure_stats(raw_pdb)
    ligand_chemistry = validate_ligand_chemistry(ligand_pdb, target_dir, ligand_resname)
    het_state_summary, het_state_details = audit_het_states(raw_pdb, target_dir, target_ph)
    (target_dir / "het_state_audit.json").write_text(
        json.dumps(het_state_details, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    het_selection_summary, _ = select_het_states(
        raw_pdb,
        target_dir,
        target_ph,
        sidecar_path=target_dir / "het_state_selection.json",
    )
    network_application_summary = apply_binding_site_network_recommendations(
        baseline_pdb,
        target_dir / "baseline_binding_site_network_applied.pdb",
        selection_sidecar=target_dir / "het_state_selection.json",
        ligand_pdb=ligand_pdb,
        apply_protein_modes=True,
        apply_water_hydrogens=False,
        audit_path=target_dir / "binding_site_network_application.json",
    )
    if _row_int(network_application_summary, "network_application_applied_count") > 0:
        baseline_pdb = Path(
            str(network_application_summary["network_application_output_pdb"])
        )
    baseline_stats = collect_structure_stats(baseline_pdb)
    input_metal_coord = metal_coordination_signature(raw_pdb)
    baseline_metal_coord = metal_coordination_signature(baseline_pdb)
    input_metal_donor_count = metal_coordination_donor_count(raw_pdb)
    baseline_metal_donor_count = metal_coordination_donor_count(baseline_pdb)
    metal_chemistry_summary, metal_chemistry_details = audit_benchmark_metal_chemistry(
        baseline_pdb
    )
    (target_dir / "metal_chemistry_audit.json").write_text(
        json.dumps(metal_chemistry_details, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    metal_bound_het_summary, metal_bound_het_details = audit_metal_bound_hets(
        raw_pdb,
        baseline_pdb,
        ligand_pdb,
        ligand_resname,
        box_size=box_size,
    )
    (target_dir / "metal_bound_het_audit.json").write_text(
        json.dumps(metal_bound_het_details, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    termini_audit = audit_termini(baseline_pdb)
    termini_details = termini_audit.pop("termini_details", [])
    (target_dir / "termini_audit.json").write_text(
        json.dumps(termini_details, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    input_geometry = audit_geometry(
        raw_pdb,
        sidecar_path=target_dir / "input_geometry_clash_audit.json",
    )
    baseline_geometry = audit_geometry(
        baseline_pdb,
        sidecar_path=target_dir / "geometry_clash_audit.json",
    )
    (
        baseline_geometry_binding_site_clash_count,
        baseline_geometry_nonbinding_clash_count,
    ) = _binding_site_clash_counts(
        target_dir / "geometry_clash_audit.json",
        ligand_center,
    )
    probe_audit = run_probe_validation(baseline_pdb)

    receptor_pdbqt = target_dir / "baseline_receptor.pdbqt"
    meeko_ok, meeko_error, meeko_input_stage, meeko_primary_error = (
        run_meeko_with_repair_fallback(baseline_pdb, repaired_pdb, receptor_pdbqt)
    )
    _write_meeko_stage_sidecar(
        receptor_pdbqt.with_suffix(".meeko_input_stage.json"),
        ok=meeko_ok,
        stage=meeko_input_stage if meeko_ok else "failed",
        error=meeko_error if not meeko_ok else meeko_primary_error,
    )
    ad4zn_ligand_pdbqt: Path | None = None
    ad4zn_ligand_status = "not_requested"
    if run_autodock4zn:
        ad4zn_ligand_pdbqt, ad4zn_ligand_status = _prepare_redock_ligand(
            ligand_pdb,
            target_dir,
            ligand_resname,
        )
    autodock4zn_plan = build_autodock4zn_plan(
        receptor_pdbqt=receptor_pdbqt,
        metal_rows=metal_chemistry_details,
        work_dir=target_dir / "autodock4zn",
        ligand_pdbqt=ad4zn_ligand_pdbqt,
        center=ligand_center,
        box_size=box_size,
        exhaustiveness=exhaustiveness,
    )
    autodock4zn_plan_path = target_dir / "autodock4zn_plan.json"
    autodock4zn_plan["plan_path"] = str(autodock4zn_plan_path)
    autodock4zn_plan["ligand_prep_status"] = ad4zn_ligand_status
    autodock4zn_result = (
        run_autodock4zn_plan(autodock4zn_plan)
        if run_autodock4zn
        else {"status": "disabled"}
    )
    autodock4zn_plan["execution_result"] = autodock4zn_result
    write_autodock4zn_plan(autodock4zn_plan_path, autodock4zn_plan)

    parameterization_plan_path = target_dir / "metal_parameterization_plan.json"
    parameterization_plan = _out_of_scope_parameterization_plan()
    parameterization_plan["plan_path"] = str(parameterization_plan_path)
    write_parameterization_plan(parameterization_plan_path, parameterization_plan)

    keep_all = target_dir / "water_keep_all.pdb"
    remove_all = target_dir / "water_remove_all.pdb"
    write_water_variant(raw_pdb, keep_all, keep_water=True)
    write_water_variant(raw_pdb, remove_all, keep_water=False)
    keep_stats = collect_structure_stats(keep_all)
    remove_stats = collect_structure_stats(remove_all)
    if run_water_analysis:
        keep_ok, _ = run_meeko_receptor(keep_all, target_dir / "water_keep_all.pdbqt")
        remove_ok, _ = run_meeko_receptor(
            remove_all,
            target_dir / "water_remove_all.pdbqt",
        )
    else:
        keep_ok = False
        remove_ok = False
    baseline_remove_all = target_dir / "water_baseline_remove_all.pdb"
    baseline_remove_all_pdbqt = target_dir / "water_baseline_remove_all.pdbqt"
    write_water_variant(baseline_pdb, baseline_remove_all, keep_water=False)
    if run_water_analysis:
        baseline_remove_all_meeko_ok, baseline_remove_all_meeko_error = (
            run_meeko_receptor(
                baseline_remove_all,
                baseline_remove_all_pdbqt,
            )
        )
    else:
        baseline_remove_all_meeko_ok = False
        baseline_remove_all_meeko_error = "water_analysis_skipped"
    water_policy_summary, water_policy_rows = audit_water_policy(
        raw_pdb,
        baseline_pdb,
        ligand_pdb,
        sidecar_path=target_dir / "water_policy_audit.json",
    )
    water_selection_summary = write_supported_water_receptor(
        baseline_pdb,
        target_dir / "water_policy_selected.pdb",
        water_policy_rows,
        source_water_pdb=raw_pdb,
    )
    water_evidence_summary, water_policy_rows = augment_water_policy_evidence(
        raw_pdb=raw_pdb,
        prepared_pdb=baseline_pdb,
        ligand_pdb=ligand_pdb,
        rows=water_policy_rows,
        target_dir=target_dir,
        pdb_id=pdb_id,
    )
    selected_water_receptor = Path(str(water_selection_summary["water_policy_selected_receptor"]))
    water_network_summary = apply_binding_site_network_recommendations(
        selected_water_receptor,
        target_dir / "water_policy_selected_network_applied.pdb",
        selection_sidecar=target_dir / "het_state_selection.json",
        ligand_pdb=ligand_pdb,
        water_policy_rows=water_policy_rows,
        apply_protein_modes=False,
        apply_water_hydrogens=True,
        audit_path=target_dir / "water_binding_site_network_application.json",
    )
    water_network_row = {
        f"water_{key}": value for key, value in water_network_summary.items()
    }
    if _row_int(water_network_summary, "network_application_applied_count") > 0:
        selected_water_receptor = Path(
            str(water_network_summary["network_application_output_pdb"])
        )
        water_selection_summary["water_policy_selected_receptor"] = str(
            selected_water_receptor
        )
    selected_water_pdbqt = target_dir / "water_policy_selected.pdbqt"
    selected_water_meeko_ok = False
    if run_water_analysis and selected_water_receptor.exists():
        selected_water_meeko_ok, selected_water_meeko_error = run_meeko_receptor(
            selected_water_receptor,
            selected_water_pdbqt,
        )
    elif not run_water_analysis:
        selected_water_meeko_error = "water_analysis_skipped"
    else:
        selected_water_meeko_error = "missing_selected_water_receptor"
    water_selection_summary["water_policy_selected_meeko_ok"] = selected_water_meeko_ok
    water_selection_summary["water_policy_selected_meeko_error"] = selected_water_meeko_error
    _write_meeko_stage_sidecar(
        target_dir / "water_policy_selected.meeko_input_stage.json",
        ok=selected_water_meeko_ok,
        stage=selected_water_meeko_error if selected_water_meeko_ok else "failed",
        error="" if selected_water_meeko_ok else selected_water_meeko_error,
    )

    redock_rmsd: float | None = None
    redock_status = "disabled"
    redock_ligand_prep_status = "disabled"
    water_redock_summary = empty_water_redock_sensitivity_summary()
    if run_redock:
        redock_receptor_pdbqt = (
            selected_water_pdbqt
            if publication_mode and selected_water_meeko_ok
            else receptor_pdbqt
        )
        redock_rmsd, redock_status, redock_ligand_prep_status = run_optional_redock(
            ligand_pdb,
            redock_receptor_pdbqt,
            target_dir,
            ligand_resname=ligand_resname,
            box_size=box_size,
            exhaustiveness=exhaustiveness,
        )
        water_redock_summary = audit_water_redock_sensitivity(
            ligand_pdb=ligand_pdb,
            ligand_resname=ligand_resname,
            variants=_water_redock_variants(
                receptor_pdbqt=receptor_pdbqt,
                selected_water_pdbqt=selected_water_pdbqt,
                selected_water_meeko_ok=selected_water_meeko_ok,
                remove_all_pdbqt=baseline_remove_all_pdbqt,
                remove_all_meeko_ok=baseline_remove_all_meeko_ok,
            ),
            out_dir=target_dir / "water_redock_sensitivity",
            sidecar_path=target_dir / "water_redock_sensitivity.json",
            redock_runner=run_optional_redock,
            box_size=box_size,
            exhaustiveness=exhaustiveness,
            redock_rmsd_max=publication_redock_rmsd_max,
        )

    minimization_status = "disabled"
    needs_cleanup_minimization = str(geometry_fix.get("geometry_fix_status", "")) != "unchanged"
    if run_minimize or (publication_mode and needs_cleanup_minimization):
        minimization_status = run_restrained_minimization_check(
            baseline_pdb,
            target_dir / "baseline_restrained_minimized.pdb",
            target_ph=target_ph,
        )

    unexpected_metal_bound_het_loss = _row_int(
        metal_bound_het_summary,
        "metal_bound_het_unexpected_removed_count",
    )
    row: BenchmarkRow = {
        "pdb_id": pdb_id.upper(),
        "target_ph": target_ph,
        "input_heavy_atoms": input_stats.heavy_atoms,
        "baseline_heavy_atoms": baseline_stats.heavy_atoms,
        "heavy_atom_delta": baseline_stats.heavy_atoms - input_stats.heavy_atoms,
        "input_hydrogens": input_stats.hydrogens,
        "baseline_hydrogens": baseline_stats.hydrogens,
        "input_waters": input_stats.waters,
        "baseline_waters": baseline_stats.waters,
        "metals_input": input_stats.metals,
        "metals_baseline": baseline_stats.metals,
        "cofactors_input": input_stats.cofactors,
        "cofactors_baseline": baseline_stats.cofactors,
        "ligand_resname": ligand_resname,
        "ligand_heavy_atoms": ligand_heavy_atoms,
        **ligand_chemistry,
        **het_state_summary,
        **het_selection_summary,
        **network_application_summary,
        "metal_coord_input": input_metal_coord,
        "metal_coord_baseline": baseline_metal_coord,
        "metal_coord_input_donor_count": input_metal_donor_count,
        "metal_coord_baseline_donor_count": baseline_metal_donor_count,
        "metal_coord_donor_delta": baseline_metal_donor_count
        - input_metal_donor_count,
        "metal_coord_raw_first_shell_retained": (
            baseline_metal_donor_count >= input_metal_donor_count
        ),
        "metal_coord_retained": (
            baseline_stats.metals >= input_stats.metals
            and (not input_metal_coord or bool(baseline_metal_coord))
            and unexpected_metal_bound_het_loss == 0
        ),
        **metal_chemistry_summary,
        **summarize_autodock4zn_plan(autodock4zn_plan),
        "autodock4zn_plan_path": str(autodock4zn_plan_path),
        "autodock4zn_execution_status": str(autodock4zn_result.get("status", "")),
        "autodock4zn_ligand_prep_status": ad4zn_ligand_status,
        **summarize_parameterization_plan(parameterization_plan),
        "metal_parameterization_plan_path": str(parameterization_plan_path),
        **metal_bound_het_summary,
        **termini_audit,
        "input_geometry_status": input_geometry["geometry_status"],
        "input_geometry_clash_count": input_geometry["geometry_clash_count"],
        "input_geometry_min_nonbonded_distance_a": input_geometry[
            "geometry_min_nonbonded_distance_a"
        ],
        "input_geometry_total_close_contact_count": input_geometry[
            "geometry_total_close_contact_count"
        ],
        "input_geometry_ignored_close_contact_count": input_geometry[
            "geometry_ignored_close_contact_count"
        ],
        "baseline_geometry_status": baseline_geometry["geometry_status"],
        "baseline_geometry_clash_count": baseline_geometry["geometry_clash_count"],
        "baseline_geometry_min_nonbonded_distance_a": baseline_geometry[
            "geometry_min_nonbonded_distance_a"
        ],
        "baseline_geometry_pre_fix_clash_count": pre_fix_geometry[
            "geometry_clash_count"
        ],
        "baseline_geometry_pre_fix_min_nonbonded_distance_a": pre_fix_geometry[
            "geometry_min_nonbonded_distance_a"
        ],
        "baseline_geometry_pre_fix_clash_audit_path": pre_fix_geometry[
            "geometry_clash_audit_path"
        ],
        "baseline_geometry_total_close_contact_count": baseline_geometry[
            "geometry_total_close_contact_count"
        ],
        "baseline_geometry_ignored_close_contact_count": baseline_geometry[
            "geometry_ignored_close_contact_count"
        ],
        "baseline_geometry_binding_site_clash_count": (
            baseline_geometry_binding_site_clash_count
        ),
        "baseline_geometry_nonbinding_clash_count": (
            baseline_geometry_nonbinding_clash_count
        ),
        "baseline_geometry_expected_metal_coordination_count": baseline_geometry[
            "geometry_expected_metal_coordination_count"
        ],
        "baseline_geometry_peptide_covalent_contact_count": baseline_geometry[
            "geometry_peptide_covalent_contact_count"
        ],
        "baseline_geometry_metal_ligand_adjacent_contact_count": baseline_geometry[
            "geometry_metal_ligand_adjacent_contact_count"
        ],
        "baseline_geometry_clash_audit_path": baseline_geometry[
            "geometry_clash_audit_path"
        ],
        **geometry_fix,
        "ph_state_residues": baseline_stats.ph_state_residues,
        "pdb2pqr_used": pdb2pqr_used,
        "reduce_or_fallback_used": reduce_used,
        "reduce_available": _tool_available("reduce", str(REDUCE_EXE)),
        "probe_available": _tool_available("probe"),
        "hbond_network_status": "reduce_probe_available"
        if _tool_available("reduce", str(REDUCE_EXE)) and _tool_available("probe")
        else "partial_no_reduce_probe",
        **probe_audit,
        "restrained_minimization_status": minimization_status,
        "meeko_ok": meeko_ok,
        "meeko_error": meeko_error,
        "meeko_input_stage": meeko_input_stage,
        "meeko_primary_error": meeko_primary_error,
        "water_keep_all_heavy_atoms": keep_stats.heavy_atoms,
        "water_keep_all_meeko_ok": keep_ok,
        "water_remove_all_heavy_atoms": remove_stats.heavy_atoms,
        "water_remove_all_meeko_ok": remove_ok,
        "water_baseline_remove_all_meeko_ok": baseline_remove_all_meeko_ok,
        "water_baseline_remove_all_meeko_error": baseline_remove_all_meeko_error,
        **water_policy_summary,
        **water_evidence_summary,
        **water_selection_summary,
        **water_network_row,
        **water_redock_summary,
        "redock_rmsd_a": _finite_or_none(redock_rmsd),
        "redock_status": redock_status,
        "redock_ligand_prep_status": redock_ligand_prep_status,
    }
    row.update(classify_metal_publication_policy(row))
    valid, failures = classify_validity(row)
    row["structure_valid"] = valid
    row["validity_failures"] = failures
    publication = classify_publication_readiness(
        row,
        require_redock=publication_mode,
        redock_rmsd_max=publication_redock_rmsd_max,
        require_minimization=publication_mode,
    )
    readiness_path = target_dir / "publication_readiness.json"
    row.update(publication)
    row["publication_readiness_path"] = str(readiness_path)
    readiness_path.write_text(
        json.dumps(
            {
                "pdb_id": pdb_id.upper(),
                **publication,
                "publication_mode": publication_mode,
                "redock_rmsd_max_a": publication_redock_rmsd_max,
                "structure_valid": valid,
                "validity_failures": failures,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return row


def main(argv: Sequence[str] | None = None) -> int:
    from protein_prep.benchmark.cli import main as cli_main

    return cli_main(argv)


if __name__ == "__main__":
    raise SystemExit(main())
