"""Context-scored HET tautomer/protomer selection for prep benchmarks."""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from protein_prep.network.optimizer import optimize_binding_site_network
from protein_prep.het_atom.ccd import download_ccd_sdf
from protein_prep.pdb_records import line_xyz, residue_key

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
ACIDIC_RESIDUES = {"ASP", "GLU"}
BASIC_RESIDUES = {"ARG", "LYS", "HIP"}
FLIPPABLE_RESIDUES = {"ASN", "GLN", "HIS", "HID", "HIE", "HIP"}
BORDERLINE_TITRATABLE_RESIDUES = {"ASP", "GLU", "HIS", "HID", "HIE", "HIP"}
DONOR_ELEMENTS = {"N", "S"}
ACCEPTOR_ELEMENTS = {"O", "S"}
POLAR_ELEMENTS = {"N", "O", "S"}
METAL_ELEMENTS = {"CA", "CD", "CO", "CU", "FE", "MG", "MN", "NI", "ZN"}


@dataclass(frozen=True)
class HetInstance:
    key: tuple[str, str, str, str]
    resname: str
    atoms: tuple[tuple[str, tuple[float, float, float]], ...]


@dataclass(frozen=True)
class ContactAtom:
    resname: str
    element: str
    xyz: tuple[float, float, float]
    record: str


def select_het_states(
    pdb_path: Path,
    work_dir: Path,
    target_ph: float,
    *,
    sidecar_path: Path | None = None,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Select one HET state per CCD residue name using local receptor contacts."""

    instances = _collect_het_instances(pdb_path)
    if not instances:
        summary = _empty_summary()
        _write_sidecar(sidecar_path, pdb_path, target_ph, summary, [])
        return summary, []

    atoms = _collect_contact_atoms(pdb_path)
    details = [
        _select_one_resname(
            pdb_path,
            resname,
            group,
            atoms,
            work_dir,
            target_ph,
        )
        for resname, group in sorted(_group_instances(instances).items())
    ]
    summary = _summarize(details, sidecar_path)
    _write_sidecar(sidecar_path, pdb_path, target_ph, summary, details)
    return summary, details


def _select_one_resname(
    pdb_path: Path,
    resname: str,
    instances: list[HetInstance],
    atoms: list[ContactAtom],
    work_dir: Path,
    target_ph: float,
) -> dict[str, object]:
    mol = _read_ccd_mol(resname, work_dir)
    if mol is None:
        return _failed_row(resname, instances, "missing_or_invalid_ccd")
    candidates = _candidate_mols(mol, target_ph)
    if not candidates:
        return _failed_row(resname, instances, "no_candidate_states")
    context = _context_summary(instances, atoms)
    candidate_features = {
        smiles: _candidate_feature_row(smiles, candidate)
        for smiles, candidate in candidates.items()
    }
    network = optimize_binding_site_network(
        pdb_path=pdb_path,
        het_resname=resname,
        het_instances=instances,
        candidates=list(candidate_features.values()),
        target_ph=target_ph,
    )
    scored = sorted(
        (
            _score_candidate(feature, context, network)
            for feature in candidate_features.values()
        ),
        key=_score_value,
        reverse=True,
    )
    best = scored[0]
    margin = (
        _score_value(best) - _score_value(scored[1])
        if len(scored) > 1
        else math.inf
    )
    confidence = _confidence(len(scored), margin)
    review = _review_labels(
        candidate_count=len(scored),
        confidence=confidence,
        margin=margin,
        context=context,
    )
    return {
        "resname": resname,
        "status": "selected",
        "confidence": confidence,
        "instance_count": len(instances),
        "candidate_count": len(scored),
        "selected_smiles": best["smiles"],
        "selected_formal_charge": best["formal_charge"],
        "selected_score": round(_score_value(best), 3),
        "score_margin": "inf" if math.isinf(margin) else round(margin, 3),
        "selection_basis": "ccd+rdkit+dimorphite-dl+binding_site_network_score",
        "network_scope": "bounded_binding_site_water_metal_hnq_titratable_component",
        "network_optimizer": _network_public_summary(network),
        **review,
        "context": context,
        "top_candidates": scored[:8],
    }


def _candidate_mols(mol: Any, target_ph: float) -> dict[str, Any]:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    base_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    seeds, _ = _dimorphite_smiles(base_smiles, target_ph)
    seeds.add(base_smiles)
    enumerator = rdMolStandardize.TautomerEnumerator()
    enumerator.SetMaxTautomers(64)
    candidates: dict[str, Any] = {}
    for smiles in seeds:
        seed_mol = Chem.MolFromSmiles(smiles)
        if seed_mol is None:
            continue
        for tautomer in enumerator.Enumerate(seed_mol):
            key = Chem.MolToSmiles(tautomer, isomericSmiles=True)
            candidates.setdefault(key, tautomer)
            if len(candidates) >= 128:
                return candidates
    return candidates


def _dimorphite_smiles(smiles: str, target_ph: float) -> tuple[set[str], str]:
    try:
        from dimorphite_dl import protonate_smiles

        variants = protonate_smiles(
            smiles,
            ph_min=max(0.0, float(target_ph) - 0.5),
            ph_max=min(14.0, float(target_ph) + 0.5),
            max_variants=64,
            validate_output=True,
        )
        return {str(item).strip() for item in variants if str(item).strip()}, "ok"
    except Exception as exc:
        return set(), f"dimorphite_unavailable:{str(exc)[:120]}"


def _candidate_feature_row(smiles: str, mol: Any) -> dict[str, object]:
    from rdkit.Chem import rdMolDescriptors

    charge = int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms()))
    hbd = int(rdMolDescriptors.CalcNumHBD(mol))
    hba = int(rdMolDescriptors.CalcNumHBA(mol))
    return {
        "smiles": smiles,
        "formal_charge": charge,
        "hbd": hbd,
        "hba": hba,
    }


def _score_candidate(
    candidate: Mapping[str, object],
    context: dict[str, int],
    network: Mapping[str, object],
) -> dict[str, object]:
    smiles = str(candidate.get("smiles", ""))
    charge = _int_field(candidate, "formal_charge")
    hbd = _int_field(candidate, "hbd")
    hba = _int_field(candidate, "hba")
    acidic = context["acidic_contacts"]
    basic = context["basic_contacts"]
    metals = context["metal_contacts"]
    donor_contacts = context["protein_donor_contacts"]
    acceptor_contacts = context["protein_acceptor_contacts"]
    polar_contacts = context["protein_polar_contacts"]
    water_contacts = context["water_contacts"]
    flippable_contacts = context["flippable_sidechain_contacts"]
    titration_contacts = context["borderline_titration_contacts"]
    counter_contacts = acidic + basic + metals
    score = 0.0
    score += charge * acidic * 0.75
    score -= charge * basic * 0.75
    score -= charge * metals * 1.10
    score += min(hbd, acceptor_contacts) * 0.35
    score += min(hba, donor_contacts + metals) * 0.35
    score += min(hbd + hba, polar_contacts) * 0.08
    score += min(hbd + hba, water_contacts) * 0.12
    score += min(hbd + hba, flippable_contacts) * 0.06
    score += min(abs(charge), titration_contacts) * 0.10
    if abs(charge) and not counter_contacts:
        score -= abs(charge) * 0.85
    if abs(charge) > 2:
        score -= (abs(charge) - 2) * 0.50
    network_row = _network_candidate_row(network, smiles)
    score += _float_field(network_row, "network_score")
    return {
        "smiles": smiles,
        "formal_charge": charge,
        "hbd": hbd,
        "hba": hba,
        "score": round(score, 3),
        "network_score": network_row.get("network_score", 0.0),
        "network_mode_explanation": network_row.get("mode_explanation", ""),
    }


def _int_field(row: Mapping[str, object], key: str) -> int:
    try:
        return int(str(row.get(key, 0)))
    except Exception:
        return 0


def _float_field(row: Mapping[str, object], key: str) -> float:
    try:
        return float(str(row.get(key, 0.0)))
    except Exception:
        return 0.0


def _network_candidate_row(
    network: Mapping[str, object],
    smiles: str,
) -> dict[str, object]:
    raw = network.get("candidate_scores", {})
    if not isinstance(raw, Mapping):
        return {}
    row = raw.get(smiles, {})
    return dict(row) if isinstance(row, Mapping) else {}


def _network_public_summary(network: Mapping[str, object]) -> dict[str, object]:
    keys = (
        "status",
        "method",
        "not_protoss_equivalent",
        "radius_a",
        "hbond_cutoff_a",
        "variable_site_count",
        "site_type_counts",
        "component_count",
        "largest_component_site_count",
        "mode_count_estimate",
        "mode_search",
        "recommended_site_modes",
    )
    return {key: network.get(key) for key in keys}


def _score_value(row: dict[str, object]) -> float:
    raw = row.get("score", 0.0)
    try:
        return float(str(raw))
    except Exception:
        return 0.0


def _review_labels(
    *,
    candidate_count: int,
    confidence: str,
    margin: float,
    context: Mapping[str, int],
) -> dict[str, object]:
    labels: list[str] = ["open_binding_site_network_optimizer_not_protoss_equivalent"]
    if candidate_count > 1:
        labels.append("multiple_state_candidates")
    if confidence == "low":
        labels.append("low_margin_state_selection")
    if _context_total(context) <= 0:
        labels.append("no_local_receptor_context")
    if _context_metal_contacts(context) > 0:
        labels.append("metal_adjacent_het_state")
    return {
        "review_required": False,
        "review_label": ";".join(labels),
        "review_reasons": labels,
        "score_margin_numeric": "" if math.isinf(margin) else round(margin, 3),
    }


def _context_total(context: Mapping[str, int]) -> int:
    return sum(int(value) for value in context.values())


def _context_metal_contacts(context: Mapping[str, int]) -> int:
    return int(context.get("metal_contacts", 0))


def _context_summary(
    instances: Iterable[HetInstance],
    atoms: list[ContactAtom],
    *,
    radius: float = 4.0,
) -> dict[str, int]:
    radius2 = radius * radius
    counts: Counter[str] = Counter()
    het_points = [xyz for instance in instances for _, xyz in instance.atoms]
    for atom in atoms:
        if _atom_near_het(atom, het_points, radius2):
            _update_context_counts(counts, atom)
    return _context_counts_dict(counts)


def _atom_near_het(
    atom: ContactAtom,
    het_points: Iterable[tuple[float, float, float]],
    radius2: float,
) -> bool:
    return any(_distance2(atom.xyz, xyz) <= radius2 for xyz in het_points)


def _update_context_counts(counts: Counter[str], atom: ContactAtom) -> None:
    if atom.record == "HETATM" and atom.element in METAL_ELEMENTS:
        counts["metal_contacts"] += 1
    if atom.record == "HETATM" and atom.resname in WATER_NAMES:
        counts["water_contacts"] += 1
    if atom.record != "ATOM":
        return
    if atom.resname in ACIDIC_RESIDUES:
        counts["acidic_contacts"] += 1
    if atom.resname in BASIC_RESIDUES:
        counts["basic_contacts"] += 1
    if atom.resname in FLIPPABLE_RESIDUES:
        counts["flippable_sidechain_contacts"] += 1
    if atom.resname in BORDERLINE_TITRATABLE_RESIDUES:
        counts["borderline_titration_contacts"] += 1
    _update_polar_context_counts(counts, atom)


def _update_polar_context_counts(counts: Counter[str], atom: ContactAtom) -> None:
    if atom.element in POLAR_ELEMENTS:
        counts["protein_polar_contacts"] += 1
    if atom.element in DONOR_ELEMENTS:
        counts["protein_donor_contacts"] += 1
    if atom.element in ACCEPTOR_ELEMENTS:
        counts["protein_acceptor_contacts"] += 1


def _context_counts_dict(counts: Counter[str]) -> dict[str, int]:
    return {
        "acidic_contacts": counts["acidic_contacts"],
        "basic_contacts": counts["basic_contacts"],
        "metal_contacts": counts["metal_contacts"],
        "protein_polar_contacts": counts["protein_polar_contacts"],
        "protein_donor_contacts": counts["protein_donor_contacts"],
        "protein_acceptor_contacts": counts["protein_acceptor_contacts"],
        "water_contacts": counts["water_contacts"],
        "flippable_sidechain_contacts": counts["flippable_sidechain_contacts"],
        "borderline_titration_contacts": counts["borderline_titration_contacts"],
    }


def _collect_het_instances(path: Path) -> list[HetInstance]:
    grouped: dict[tuple[str, str, str, str], list[tuple[str, tuple[float, float, float]]]] = defaultdict(list)
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("HETATM") or _line_element(line) == "H":
                continue
            resname = line[17:20].strip().upper()
            if resname in WATER_NAMES or resname in COMMON_IONS:
                continue
            xyz = line_xyz(line)
            if xyz is None:
                continue
            grouped[residue_key(line)].append((line[12:16].strip().upper(), xyz))
    return [
        HetInstance(key=key, resname=key[0], atoms=tuple(atoms))
        for key, atoms in sorted(grouped.items())
    ]


def _collect_contact_atoms(path: Path) -> list[ContactAtom]:
    atoms: list[ContactAtom] = []
    if not path.exists():
        return atoms
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")) or _line_element(line) == "H":
                continue
            xyz = line_xyz(line)
            if xyz is None:
                continue
            atoms.append(
                ContactAtom(
                    resname=line[17:20].strip().upper(),
                    element=_line_element(line),
                    xyz=xyz,
                    record=line[:6].strip(),
                )
            )
    return atoms


def _group_instances(instances: Iterable[HetInstance]) -> dict[str, list[HetInstance]]:
    grouped: dict[str, list[HetInstance]] = defaultdict(list)
    for instance in instances:
        grouped[instance.resname].append(instance)
    return grouped


def _read_ccd_mol(resname: str, work_dir: Path) -> Any | None:
    from rdkit import Chem

    sdf_path = download_ccd_sdf(resname, work_dir)
    if sdf_path is None:
        return None
    supplier = Chem.SDMolSupplier(str(sdf_path), sanitize=True, removeHs=False)
    return next((item for item in supplier if item is not None), None)


def _summarize(
    details: list[dict[str, object]],
    sidecar_path: Path | None,
) -> dict[str, object]:
    failed = [row for row in details if row.get("status") != "selected"]
    low_confidence = [
        row for row in details if str(row.get("confidence", "")) == "low"
    ]
    review_labels = _selection_review_labels(details)
    status = "selected"
    if failed and len(failed) == len(details):
        status = "failed"
    elif failed:
        status = "partial"
    return {
        "het_state_selection_status": status,
        "het_state_selected_resname_count": len(details) - len(failed),
        "het_state_selection_failed_count": len(failed),
        "het_state_selection_low_confidence_count": len(low_confidence),
        "het_state_selection_review_label_count": len(review_labels),
        "het_state_selection_review_labels": ";".join(review_labels),
        "het_state_context_variant_count": sum(
            _candidate_count(row) for row in details
        ),
        "het_state_selection_tool": (
            "ccd+rdkit+dimorphite-dl+binding_site_network_optimizer"
        ),
        "het_state_selection_path": str(sidecar_path) if sidecar_path else "",
    }


def _candidate_count(row: dict[str, object]) -> int:
    try:
        return int(str(row.get("candidate_count") or 0))
    except Exception:
        return 0


def _selection_review_labels(details: Sequence[Mapping[str, object]]) -> list[str]:
    labels: set[str] = set()
    for row in details:
        raw = row.get("review_reasons", [])
        if isinstance(raw, list):
            labels.update(str(item) for item in raw if str(item))
    return sorted(labels)


def _empty_summary() -> dict[str, object]:
    return {
        "het_state_selection_status": "no_het",
        "het_state_selected_resname_count": 0,
        "het_state_selection_failed_count": 0,
        "het_state_selection_low_confidence_count": 0,
        "het_state_context_variant_count": 0,
        "het_state_selection_tool": "",
        "het_state_selection_path": "",
    }


def _failed_row(
    resname: str,
    instances: list[HetInstance],
    status: str,
) -> dict[str, object]:
    return {
        "resname": resname,
        "status": status,
        "confidence": "none",
        "instance_count": len(instances),
        "candidate_count": 0,
    }


def _confidence(candidate_count: int, margin: float) -> str:
    if candidate_count <= 1:
        return "single"
    if margin >= 1.0:
        return "high"
    if margin >= 0.35:
        return "medium"
    return "low"


def _distance2(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
) -> float:
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def _line_element(line: str) -> str:
    return (line[76:78].strip() if len(line) >= 78 else line[12:16].strip()[:1]).upper()


def _write_sidecar(
    sidecar_path: Path | None,
    pdb_path: Path,
    target_ph: float,
    summary: dict[str, object],
    details: list[dict[str, object]],
) -> None:
    if sidecar_path is None:
        return
    sidecar_path.parent.mkdir(parents=True, exist_ok=True)
    sidecar_path.write_text(
        json.dumps(
            {
                "pdb_path": str(pdb_path),
                "target_ph": target_ph,
                "summary": summary,
                "selections": details,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
