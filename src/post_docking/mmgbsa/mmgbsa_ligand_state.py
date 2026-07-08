from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable, Mapping


COMMON_ORGANIC_ELEMENTS = {
    "B",
    "BR",
    "C",
    "CL",
    "F",
    "H",
    "I",
    "N",
    "O",
    "P",
    "S",
    "SE",
    "SI",
}
METAL_ELEMENTS = {
    "AL",
    "CA",
    "CD",
    "CO",
    "CU",
    "FE",
    "HG",
    "K",
    "LI",
    "MG",
    "MN",
    "NA",
    "NI",
    "ZN",
}
SEVERITY_RANK = {"none": 0, "info": 1, "warning": 2, "blocker": 3, "critical": 4}


def analyze_ligand_state(sdf_path: str | Path, cfg: Mapping[str, Any] | None) -> dict[str, Any]:
    """Return auditable ligand-state provenance and suspicious-chemistry flags."""
    path = Path(sdf_path)
    mol, parse_status, sanitize_error = _read_mol(path)
    ph_context = _resolve_ph_context(path, cfg)
    ph_label = str(ph_context.get("ph_label", "") or "")
    ph_value = ph_context.get("target_ph")
    provenance: dict[str, Any] = {
        "schema": "mmgbsa_ligand_state_v1",
        "input_sdf": str(path),
        "input_format": "sdf",
        "ph_label": ph_label,
        "ph_value": ph_value,
        "ph_context": ph_context,
        "parse_status": parse_status,
        "sanitize_error": sanitize_error,
        "protomer_source": "input_structure_preserved",
        "tautomer_source": "input_structure_preserved",
    }
    suspicious = {
        "schema": "mmgbsa_suspicious_chemistry_v1",
        "severity": "none",
        "reasons": [],
        "details": {},
    }
    if mol is None:
        _add_reason(
            suspicious,
            "blocker",
            "rdkit_parse_failed",
            {"parse_status": parse_status, "sanitize_error": sanitize_error},
        )
        return {
            "ligand_state_provenance": provenance,
            "suspicious_chemistry": suspicious,
        }

    try:
        _populate_rdkit_provenance(provenance, suspicious, mol)
    except Exception as exc:
        message = str(exc)[:240]
        provenance["state_analysis_error"] = message
        _add_reason(suspicious, "blocker", "rdkit_state_analysis_failed", message)
        return {
            "ligand_state_provenance": provenance,
            "suspicious_chemistry": suspicious,
        }
    _classify_ph_context(ph_context, suspicious)
    smiles = str(provenance.get("canonical_isomeric_smiles", "") or "")
    if smiles:
        provenance["rdkit_tautomer"] = _enumerate_rdkit_tautomers(smiles)
        provenance["dimorphite"] = _enumerate_dimorphite_states(smiles, ph_value, cfg)
        provenance["state_inference"] = _infer_state_selection(
            path,
            smiles,
            provenance,
            suspicious,
            cfg,
        )
    else:
        provenance["rdkit_tautomer"] = {"enabled": False, "status": "missing_smiles"}
        provenance["dimorphite"] = {"enabled": False, "status": "missing_smiles"}
        provenance["state_inference"] = {"status": "missing_smiles"}
    return {
        "ligand_state_provenance": provenance,
        "suspicious_chemistry": suspicious,
    }


def suspicious_blocks_publication(payload: Mapping[str, Any] | None) -> bool:
    if not isinstance(payload, Mapping):
        return False
    severity = str(payload.get("severity", "none") or "none").lower()
    return SEVERITY_RANK.get(severity, 0) >= SEVERITY_RANK["blocker"]


def _cfg_get(cfg: Mapping[str, Any] | None, key: str, default: Any = None) -> Any:
    if isinstance(cfg, Mapping):
        return cfg.get(key, default)
    return default


def _cfg_bool(cfg: Mapping[str, Any] | None, key: str, default: bool) -> bool:
    raw = _cfg_get(cfg, key, default)
    text = str(raw).strip().lower()
    if text in {"1", "true", "yes", "on"}:
        return True
    if text in {"0", "false", "no", "off"}:
        return False
    return bool(default)


def _cfg_int(cfg: Mapping[str, Any] | None, key: str, default: int) -> int:
    try:
        return int(_cfg_get(cfg, key, default))
    except Exception:
        return default


def _production_mode(cfg: Mapping[str, Any] | None) -> bool:
    return str(_cfg_get(cfg, "MMGBSA_PROTOCOL", "") or "").strip().lower() in {
        "production",
        "publication",
    }


def _resolve_ph_context(sdf_path: Path, cfg: Mapping[str, Any] | None) -> dict[str, Any]:
    label = str(_cfg_get(cfg, "MMGBSA_INPUT_PH_LABEL", "") or "")
    label_values = _parse_existing_ph_label(label)
    receptor = _existing_path(
        _cfg_get(cfg, "MMGBSA_RECEPTOR_PDB_CONTEXT")
        or _cfg_get(cfg, "MMGBSA_RECEPTOR_PDB")
        or _cfg_get(cfg, "MMGBSA_INPUT_RECEPTOR_PDB")
    )
    receptor_context = _protein_prep_ph_context(receptor)
    target = label_values[0] if label_values else receptor_context.get("target_ph", 7.0)
    source = "mmgbsa_ph_label" if label_values else receptor_context.get("source", "default")
    return {
        "schema": "mmgbsa_ligand_ph_context_v1",
        "source": source,
        "ph_label": label,
        "label_values": label_values,
        "target_ph": float(target),
        "receptor_pdb": str(receptor) if receptor else "",
        "protein_prep": receptor_context,
        "ligand_sdf": str(sdf_path),
    }


def _parse_existing_ph_label(label: str) -> list[float]:
    try:
        from docking.ph_ensemble_docking import _parse_ph_values_from_label

        return [float(item) for item in _parse_ph_values_from_label(label)]
    except Exception:
        return []


def _existing_path(raw: Any) -> Path | None:
    if raw is None:
        return None
    path = Path(str(raw)).expanduser()
    return path if path.is_file() else None


def _protein_prep_ph_context(receptor: Path | None) -> dict[str, Any]:
    if receptor is None:
        return {"source": "default", "target_ph": 7.0, "status": "missing_receptor"}
    payload: dict[str, Any] = {"source": "protein_prep_context", "status": "ok"}
    try:
        from path_router.context_ph import select_ph_from_pdb

        selected = select_ph_from_pdb(str(receptor))
        payload["context_ph"] = selected
    except Exception as exc:
        payload["context_ph_error"] = str(exc)[:160]
    try:
        from protein_prep.openmm_repair import get_target_ph_for_prep

        payload["target_ph"] = float(get_target_ph_for_prep(receptor))
    except Exception as exc:
        payload["target_ph"] = 7.0
        payload["status"] = "fallback_default"
        payload["target_ph_error"] = str(exc)[:160]
    return payload


def _classify_ph_context(
    ph_context: Mapping[str, Any], suspicious: dict[str, Any]
) -> None:
    label_values = ph_context.get("label_values")
    protein_target = _nested_float(ph_context, ("protein_prep", "target_ph"))
    target = _as_float(ph_context.get("target_ph"))
    if not isinstance(label_values, list) or not label_values or protein_target is None:
        return
    if target is not None and abs(float(target) - float(protein_target)) > 0.35:
        _add_reason(
            suspicious,
            "warning",
            "ligand_ph_context_mismatch",
            {
                "mmgbsa_target_ph": target,
                "protein_prep_target_ph": protein_target,
                "policy": "using_mmgbsa_ph_label_for_ligand_state",
            },
        )


def _nested_float(payload: Mapping[str, Any], keys: tuple[str, ...]) -> float | None:
    value: Any = payload
    for key in keys:
        if not isinstance(value, Mapping):
            return None
        value = value.get(key)
    return _as_float(value)


def _as_float(value: Any) -> float | None:
    try:
        return float(value)
    except Exception:
        return None


def _read_mol(path: Path) -> tuple[Any, str, str]:
    try:
        from rdkit import Chem

        supplier = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False)
        mol = supplier[0] if supplier and len(supplier) else None
        if mol is not None:
            return mol, "sanitized", ""
        supplier = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)
        mol = supplier[0] if supplier and len(supplier) else None
        if mol is None:
            return None, "parse_failed", ""
        try:
            Chem.SanitizeMol(mol)
            return mol, "sanitize_recovered", ""
        except Exception as exc:
            return mol, "unsanitized", str(exc)
    except Exception as exc:
        return None, "rdkit_unavailable_or_failed", str(exc)


def _populate_rdkit_provenance(
    provenance: dict[str, Any], suspicious: dict[str, Any], mol: Any
) -> None:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    atoms = list(mol.GetAtoms())
    elements, metals, unusual = _element_groups(atoms)
    formal_charge = int(Chem.GetFormalCharge(mol))
    radicals = int(sum(atom.GetNumRadicalElectrons() for atom in atoms))
    fragments = int(len(Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)))
    provenance.update(
        {
            "canonical_isomeric_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "atom_count": int(mol.GetNumAtoms()),
            "heavy_atom_count": int(mol.GetNumHeavyAtoms()),
            "formal_charge": formal_charge,
            "radical_electrons": radicals,
            "fragment_count": fragments,
            "elements": elements,
            "metal_elements": metals,
            "unusual_elements": unusual,
        }
    )
    _classify_rdkit_findings(
        provenance,
        suspicious,
        formal_charge=formal_charge,
        radicals=radicals,
        fragments=fragments,
        metals=metals,
        unusual=unusual,
    )


def _element_groups(atoms: list[Any]) -> tuple[list[str], list[str], list[str]]:
    elements = sorted({atom.GetSymbol().upper() for atom in atoms})
    metals = sorted(element for element in elements if element in METAL_ELEMENTS)
    unusual = sorted(
        element
        for element in elements
        if element not in COMMON_ORGANIC_ELEMENTS and element not in METAL_ELEMENTS
    )
    return elements, metals, unusual


def _classify_rdkit_findings(
    provenance: Mapping[str, Any],
    suspicious: dict[str, Any],
    *,
    formal_charge: int,
    radicals: int,
    fragments: int,
    metals: list[str],
    unusual: list[str],
) -> None:
    if provenance.get("sanitize_error"):
        _add_reason(suspicious, "blocker", "rdkit_sanitize_failed", provenance["sanitize_error"])
    if radicals:
        _add_reason(suspicious, "blocker", "radical_electrons_present", radicals)
    if abs(formal_charge) > 2:
        _add_reason(suspicious, "warning", "large_absolute_formal_charge", formal_charge)
    if fragments > 1:
        _add_reason(suspicious, "warning", "multi_fragment_ligand", fragments)
    if metals:
        _add_reason(suspicious, "warning", "metal_in_ligand_graph", metals)
    if unusual:
        _add_reason(suspicious, "warning", "unusual_elements", unusual)


def _add_reason(payload: dict[str, Any], severity: str, reason: str, detail: Any) -> None:
    current = str(payload.get("severity", "none") or "none").lower()
    if SEVERITY_RANK.get(severity, 0) > SEVERITY_RANK.get(current, 0):
        payload["severity"] = severity
    payload.setdefault("reasons", []).append(reason)
    payload.setdefault("details", {})[reason] = detail


def _enumerate_rdkit_tautomers(smiles: str) -> dict[str, Any]:
    try:
        from rdkit import Chem
        from rdkit.Chem.MolStandardize import rdMolStandardize

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"enabled": True, "status": "mol_from_smiles_failed", "state_count": 0}
        enumerator = rdMolStandardize.TautomerEnumerator()
        enumerator.SetMaxTautomers(64)
        states = {
            Chem.MolToSmiles(tautomer, isomericSmiles=True)
            for tautomer in enumerator.Enumerate(mol)
        }
        return {
            "enabled": True,
            "status": "ok",
            "state_count": len(states),
            "candidate_smiles": sorted(states),
            "representative_smiles": sorted(states)[:8],
        }
    except Exception as exc:
        return {"enabled": False, "status": f"rdkit_tautomer_failed:{str(exc)[:120]}"}


def _enumerate_dimorphite_states(
    smiles: str, ph_value: float | None, cfg: Mapping[str, Any] | None
) -> dict[str, Any]:
    if not _cfg_bool(cfg, "MMGBSA_DIMORPHITE_ENABLED", True):
        return {"enabled": False, "status": "disabled"}
    target_ph = 7.0 if ph_value is None else float(ph_value)
    window = float(_cfg_get(cfg, "MMGBSA_DIMORPHITE_PH_WINDOW", 0.5) or 0.5)
    max_variants = _cfg_int(cfg, "MMGBSA_DIMORPHITE_MAX_VARIANTS", 32)
    try:
        from dimorphite_dl import protonate_smiles

        variants = protonate_smiles(
            smiles,
            ph_min=max(0.0, target_ph - window),
            ph_max=min(14.0, target_ph + window),
            max_variants=max(1, max_variants),
            validate_output=True,
        )
        states = sorted({str(item).strip() for item in variants if str(item).strip()})
        return {
            "enabled": True,
            "status": "ok",
            "tool": "dimorphite-dl",
            "target_ph": target_ph,
            "ph_min": max(0.0, target_ph - window),
            "ph_max": min(14.0, target_ph + window),
            "state_count": len(states),
            "candidate_smiles": states,
            "representative_smiles": states[:8],
        }
    except Exception as exc:
        return {
            "enabled": False,
            "status": f"dimorphite_unavailable_or_failed:{str(exc)[:120]}",
        }


def _infer_state_selection(
    sdf_path: Path,
    input_smiles: str,
    provenance: Mapping[str, Any],
    suspicious: dict[str, Any],
    cfg: Mapping[str, Any] | None,
) -> dict[str, Any]:
    pose_meta = _read_pose_metadata(sdf_path)
    authoritative = _authoritative_input(cfg, pose_meta)
    candidates = _rank_state_candidates(input_smiles, provenance, authoritative)
    selected, margin = _select_state_candidate(candidates, authoritative)
    status = "selected" if margin > 10.0 or len(candidates) <= 1 else "ambiguous"
    if status == "ambiguous":
        severity = "blocker" if _production_mode(cfg) else "warning"
        _add_reason(
            suspicious,
            severity,
            "protomer_tautomer_state_ambiguous",
            {"top_margin": margin, "candidate_count": len(candidates)},
        )
    return {
        "schema": "mmgbsa_ligand_state_inference_v1",
        "status": status,
        "selection_policy": (
            "preserve_authoritative_input" if authoritative else "rank_heuristic_candidates"
        ),
        "selected_state": selected,
        "top_margin": margin,
        "candidate_count": len(candidates),
        "candidates": candidates[:12],
        "evidence": {
            "solution_population_proxy": "dimorphite_pH_window_membership",
            "input_geometry_consistency": "input_smiles_identity",
            "pose_metadata": _pose_evidence(pose_meta),
            "docking_redock_evidence": _pose_score_evidence(pose_meta),
            "short_relaxation_or_mmgbsa": {"status": "not_run"},
        },
    }


def _read_pose_metadata(sdf_path: Path) -> dict[str, Any]:
    meta_path = sdf_path.with_suffix(".mmgbsa_pose.json")
    if not meta_path.is_file():
        return {}
    try:
        data = json.loads(meta_path.read_text(encoding="utf-8"))
        return dict(data) if isinstance(data, Mapping) else {}
    except Exception:
        return {}


def _authoritative_input(
    cfg: Mapping[str, Any] | None, pose_meta: Mapping[str, Any]
) -> bool:
    return _cfg_bool(cfg, "MMGBSA_INPUT_CHEMISTRY_AUTHORITATIVE", False) or bool(
        pose_meta.get("chemistry_authoritative")
        or pose_meta.get("input_chemistry_authoritative")
    )


def _rank_state_candidates(
    input_smiles: str, provenance: Mapping[str, Any], authoritative: bool
) -> list[dict[str, Any]]:
    dimorphite_states = _candidate_set(provenance.get("dimorphite"))
    tautomer_states = _candidate_set(provenance.get("rdkit_tautomer"))
    candidate_smiles: list[str] = []
    seen: set[str] = set()
    for smiles in [input_smiles, *sorted(dimorphite_states | tautomer_states)]:
        if smiles and smiles not in seen:
            seen.add(smiles)
            candidate_smiles.append(smiles)
    ranked = [
        _score_candidate(
            smiles,
            input_smiles=input_smiles,
            dimorphite_states=dimorphite_states,
            tautomer_states=tautomer_states,
            authoritative=authoritative,
        )
        for smiles in candidate_smiles
    ]
    return sorted(ranked, key=lambda item: float(item.get("score", 0.0)), reverse=True)


def _candidate_set(payload: Any) -> set[str]:
    if not isinstance(payload, Mapping):
        return set()
    values = payload.get("candidate_smiles") or payload.get("representative_smiles") or []
    return {str(item).strip() for item in values if str(item).strip()}


def _score_candidate(
    smiles: str,
    *,
    input_smiles: str,
    dimorphite_states: set[str],
    tautomer_states: set[str],
    authoritative: bool,
) -> dict[str, Any]:
    quality = _screen_state_smiles(smiles)
    terms = {
        "preserve_input": 70.0 if authoritative and smiles == input_smiles else 0.0,
        "input_identity": 25.0 if smiles == input_smiles else 0.0,
        "dimorphite_pH_window": 20.0 if smiles in dimorphite_states else 0.0,
        "rdkit_tautomer_candidate": 8.0 if smiles in tautomer_states else 0.0,
        "quality_penalty": -float(quality.get("penalty", 0.0)),
    }
    if not quality.get("ok", False):
        terms["quality_penalty"] -= 100.0
    return {
        "smiles": smiles,
        "score": round(sum(terms.values()), 3),
        "terms": terms,
        "quality": quality,
        "is_input_state": smiles == input_smiles,
    }


def _screen_state_smiles(smiles: str) -> dict[str, Any]:
    try:
        from rdkit import Chem

        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            return {"ok": False, "reasons": ["mol_from_smiles_failed"], "penalty": 100.0}
        charge = int(Chem.GetFormalCharge(mol))
        radicals = int(sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()))
        fragments = int(len(Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)))
        reasons = _state_quality_reasons(charge, radicals, fragments)
        return {
            "ok": not any(reason.startswith("blocker:") for reason in reasons),
            "reasons": reasons,
            "formal_charge": charge,
            "radical_electrons": radicals,
            "fragment_count": fragments,
            "penalty": _state_quality_penalty(charge, radicals, fragments),
        }
    except Exception as exc:
        return {"ok": False, "reasons": [f"screen_failed:{str(exc)[:120]}"], "penalty": 100.0}


def _state_quality_reasons(charge: int, radicals: int, fragments: int) -> list[str]:
    reasons: list[str] = []
    if radicals:
        reasons.append("blocker:radical_electrons")
    if fragments > 1:
        reasons.append("warning:multi_fragment")
    if abs(charge) > 2:
        reasons.append("warning:large_absolute_charge")
    return reasons


def _state_quality_penalty(charge: int, radicals: int, fragments: int) -> float:
    return (100.0 if radicals else 0.0) + max(0, fragments - 1) * 15.0 + max(
        0, abs(charge) - 2
    ) * 10.0


def _select_state_candidate(
    candidates: list[dict[str, Any]], authoritative: bool
) -> tuple[dict[str, Any], float]:
    if not candidates:
        return {}, 0.0
    selected = _input_candidate(candidates) if authoritative else candidates[0]
    second_score = _best_other_score(candidates, selected)
    margin = float(selected.get("score", 0.0)) - second_score
    return selected, round(margin, 3)


def _input_candidate(candidates: Iterable[dict[str, Any]]) -> dict[str, Any]:
    for candidate in candidates:
        if candidate.get("is_input_state"):
            return candidate
    return {}


def _best_other_score(candidates: list[dict[str, Any]], selected: Mapping[str, Any]) -> float:
    selected_smiles = selected.get("smiles")
    other_scores = [
        float(item.get("score", 0.0))
        for item in candidates
        if item.get("smiles") != selected_smiles
    ]
    return max(other_scores) if other_scores else float(selected.get("score", 0.0)) - 999.0


def _pose_evidence(pose_meta: Mapping[str, Any]) -> dict[str, Any]:
    if not pose_meta:
        return {"status": "not_available"}
    return {
        "status": "available",
        "chemistry_authoritative": bool(
            pose_meta.get("chemistry_authoritative")
            or pose_meta.get("input_chemistry_authoritative")
        ),
        "source_kind": pose_meta.get("source_kind", ""),
        "source_used": pose_meta.get("source_used", ""),
        "pdbqt_to_sdf_reconstruction": bool(pose_meta.get("pdbqt_to_sdf_reconstruction")),
    }


def _pose_score_evidence(pose_meta: Mapping[str, Any]) -> dict[str, Any]:
    keys = ("docking_score", "score", "redock_rmsd", "rmsd")
    values = {key: pose_meta.get(key) for key in keys if key in pose_meta}
    return {"status": "available", **values} if values else {"status": "not_available"}
