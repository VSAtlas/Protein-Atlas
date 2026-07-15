from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import pytest
from rdkit import Chem
from rdkit.Geometry import Point3D

from docking import docking_ligands
from docking.native_redock_rmsd import (
    NATIVE_REDOCK_RMSD_LITERATURE_DOIS,
    NATIVE_REDOCK_RMSD_METHOD_ID,
    RECEPTOR_TRANSFORM_SERIALIZATION_TOLERANCE_A,
    canonical_json_sha256,
    evaluate_native_redock_rmsd,
    sha256_file,
)
from docking.pose_validation import compute_redock_rmsd


_SMILES = "c1ccccc1"
_NATIVE_COORDS = (
    (1.400, 0.000, 0.000),
    (0.700, 1.212, 0.000),
    (-0.700, 1.212, 0.000),
    (-1.400, 0.000, 0.000),
    (-0.700, -1.212, 0.000),
    (0.700, -1.212, 0.000),
)
_TRANSLATION = [
    [1.0, 0.0, 0.0, 3.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
]
_MAPPING = "REMARK SMILES IDX 1 1 2 2 3 3 4 4 5 5 6 6\n"


@dataclass
class _Case:
    reference: Path
    graph: Path
    prepared_pdbqt: Path
    poses: Path
    source_receptor: Path
    prepared_receptor: Path
    extraction: Path
    frame: Path
    selection: Path
    frame_sha256: str
    selection_sha256: str
    model_blocks: tuple[bytes, ...]
    model_scores: tuple[float, ...]


def _translated(
    coords: Sequence[tuple[float, float, float]],
    dx: float,
) -> tuple[tuple[float, float, float], ...]:
    return tuple((x + dx, y, z) for x, y, z in coords)


def _write_sdf(
    path: Path,
    *,
    smiles: str = _SMILES,
    coords: Sequence[tuple[float, float, float]] = _NATIVE_COORDS,
    reverse_atoms: bool = False,
) -> None:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    assert mol.GetNumAtoms() == len(coords)
    conformer = Chem.Conformer(mol.GetNumAtoms())
    conformer.Set3D(True)
    for atom_index, coord in enumerate(coords):
        conformer.SetAtomPosition(atom_index, Point3D(*coord))
    mol.AddConformer(conformer, assignId=True)
    if reverse_atoms:
        mol = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()


def _atom_line(
    index: int,
    coord: tuple[float, float, float],
    *,
    atom_type: str = "A",
) -> str:
    x, y, z = coord
    return (
        f"ATOM  {index:5d}  {f'C{index}':<3} LIG     1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  0.00  0.00     0.000 {atom_type}\n"
    )


def _receptor_atom_line(
    index: int,
    atom_name: str,
    residue_number: int,
    coord: tuple[float, float, float],
) -> str:
    x, y, z = coord
    return (
        f"ATOM  {index:5d} {atom_name:>4} ALA A{residue_number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n"
    )


def _pdbqt_body(
    coords: Sequence[tuple[float, float, float]],
    *,
    mapping: str = _MAPPING,
    atom_type: str = "A",
) -> bytes:
    text = f"REMARK SMILES {_SMILES}\n{mapping}ROOT\n"
    text += "".join(
        _atom_line(index, coord, atom_type=atom_type)
        for index, coord in enumerate(coords, start=1)
    )
    text += "ENDROOT\nTORSDOF 0\n"
    return text.encode("ascii")


def _model_block(
    label: int,
    coords: Sequence[tuple[float, float, float]],
    *,
    score: float,
    mapping: str = _MAPPING,
    atom_type: str = "A",
) -> bytes:
    return (
        f"MODEL {label}\n".encode("ascii")
        + f"REMARK VINA RESULT: {score:.3f} 0.000 0.000\n".encode("ascii")
        + _pdbqt_body(coords, mapping=mapping, atom_type=atom_type)
        + b"ENDMDL\n"
    )


def _write_json(path: Path, payload: dict[str, Any]) -> str:
    path.write_text(
        json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    digest = sha256_file(path)
    assert digest is not None
    return digest


def _write_extraction_and_frame(
    case: _Case,
    *,
    transform: list[list[float]] = _TRANSLATION,
    transform_sha256: str | None = None,
) -> None:
    source_sha = sha256_file(case.source_receptor)
    prepared_sha = sha256_file(case.prepared_receptor)
    reference_sha = sha256_file(case.reference)
    assert source_sha and prepared_sha and reference_sha
    extraction_sha = _write_json(
        case.extraction,
        {
            "schema_version": 1,
            "source_receptor_sha256": source_sha,
            "reference_ligand_sha256": reference_sha,
            "extraction_method": "pdb_component_to_sdf_with_bond_orders",
            "ligand_id": "BEN:A:100",
        },
    )
    case.frame_sha256 = _write_json(
        case.frame,
        {
            "schema_version": 1,
            "source_receptor_sha256": source_sha,
            "prepared_receptor_sha256": prepared_sha,
            "reference_extraction_manifest_sha256": extraction_sha,
            "source_to_prepared_transform": transform,
            "source_to_prepared_transform_sha256": (
                transform_sha256 or canonical_json_sha256(transform)
            ),
        },
    )


def _write_selection(
    case: _Case,
    *,
    selected_ordinal: Any = 1,
    selected_label: Any = 1,
    selected_sha256: str | None = None,
    score_value: float | None = None,
) -> None:
    poses_sha = sha256_file(case.poses)
    assert poses_sha is not None
    block_index = selected_ordinal - 1 if type(selected_ordinal) is int else 0
    default_model_sha = (
        hashlib.sha256(case.model_blocks[block_index]).hexdigest()
        if 0 <= block_index < len(case.model_blocks)
        else "0" * 64
    )
    default_score = (
        case.model_scores[block_index]
        if 0 <= block_index < len(case.model_scores)
        else -8.0
    )
    case.selection_sha256 = _write_json(
        case.selection,
        {
            "schema_version": 1,
            "docked_poses_sha256": poses_sha,
            "engine": "vina",
            "stage": "native_redock",
            "score_source": "vina_affinity_kcal_mol",
            "score_value": default_score if score_value is None else score_value,
            "score_direction": "lower_is_better",
            "tie_break_rule": "model_ordinal_ascending",
            "selected_model_ordinal": selected_ordinal,
            "selected_model_label": selected_label,
            "selected_model_sha256": selected_sha256 or default_model_sha,
        },
    )


def _make_case(
    tmp_path: Path,
    *,
    pose_coords: Sequence[Sequence[tuple[float, float, float]]] | None = None,
    pose_scores: Sequence[float] | None = None,
) -> _Case:
    coords = tuple(pose_coords or (_translated(_NATIVE_COORDS, 3.0),))
    scores = tuple(
        pose_scores
        if pose_scores is not None
        else (-8.0 + index for index in range(len(coords)))
    )
    assert len(scores) == len(coords)
    reference = tmp_path / "reference.sdf"
    graph = tmp_path / "prepared_graph.sdf"
    prepared_pdbqt = tmp_path / "prepared.pdbqt"
    poses = tmp_path / "poses.pdbqt"
    source_receptor = tmp_path / "source_receptor.pdb"
    prepared_receptor = tmp_path / "prepared_receptor.pdbqt"
    extraction = tmp_path / "reference_extraction.json"
    frame = tmp_path / "coordinate_frame.json"
    selection = tmp_path / "pose_selection.json"
    _write_sdf(reference)
    # Deliberately reverse this graph's atom order. RMSD reconstruction must use
    # Meeko's explicit SMILES-IDX mapping rather than SDF/PDBQT element order.
    _write_sdf(graph, reverse_atoms=True)
    prepared_pdbqt.write_bytes(_pdbqt_body(_translated(_NATIVE_COORDS, 3.0)))
    blocks = tuple(
        _model_block(label, model_coords, score=scores[label - 1])
        for label, model_coords in enumerate(coords, start=1)
    )
    poses.write_bytes(b"".join(blocks))
    receptor_coords = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    receptor_names = ("N", "CA", "C")
    source_receptor.write_text(
        "".join(
            _receptor_atom_line(index, name, 1, coord)
            for index, (name, coord) in enumerate(
                zip(receptor_names, receptor_coords), start=1
            )
        ),
        encoding="ascii",
    )
    prepared_receptor.write_text(
        "".join(
            _receptor_atom_line(index, name, 1, coord)
            for index, (name, coord) in enumerate(
                zip(receptor_names, _translated(receptor_coords, 3.0)), start=1
            )
        ),
        encoding="ascii",
    )
    case = _Case(
        reference=reference,
        graph=graph,
        prepared_pdbqt=prepared_pdbqt,
        poses=poses,
        source_receptor=source_receptor,
        prepared_receptor=prepared_receptor,
        extraction=extraction,
        frame=frame,
        selection=selection,
        frame_sha256="",
        selection_sha256="",
        model_blocks=blocks,
        model_scores=scores,
    )
    _write_extraction_and_frame(case)
    _write_selection(case)
    return case


def _evaluate(case: _Case, **kwargs: Any):
    return evaluate_native_redock_rmsd(
        case.reference,
        case.graph,
        case.prepared_pdbqt,
        case.poses,
        source_receptor_path=case.source_receptor,
        prepared_receptor_path=case.prepared_receptor,
        coordinate_frame_manifest_path=case.frame,
        coordinate_frame_manifest_sha256=case.frame_sha256,
        reference_extraction_manifest_path=case.extraction,
        pose_selection_manifest_path=case.selection,
        pose_selection_manifest_sha256=case.selection_sha256,
        **kwargs,
    )


def test_receptor_frame_transform_score_order_and_mapping_are_verified(
    tmp_path: Path,
) -> None:
    case = _make_case(tmp_path)

    result = _evaluate(case)

    assert result.status == "qualified"
    assert result.qualified is True
    assert result.failure_reason is None
    assert result.top_ranked_rmsd_a == pytest.approx(0.0, abs=1e-3)
    assert result.best_generated_rmsd_a == pytest.approx(0.0, abs=1e-3)
    assert result.coordinate_transform_is_identity is False
    assert result.coordinate_transform_verification_status == "verified"
    assert result.coordinate_frame_relation == (
        "verified_by_exact_common_receptor_atom_correspondence"
    )
    assert result.coordinate_transform_verification_tolerance_a == (
        RECEPTOR_TRANSFORM_SERIALIZATION_TOLERANCE_A
    )
    assert result.common_receptor_atom_count == 3
    assert result.coordinate_transform_max_residual_a == pytest.approx(0.0)
    assert result.selection_score_verification_status == "verified"
    assert result.pose_rmsds[0].selection_score == pytest.approx(-8.0)
    assert result.coordinate_frame_manifest_sha256 == case.frame_sha256
    assert result.pose_selection_manifest_sha256 == case.selection_sha256
    assert result.meeko_mapping_sha256 is not None
    assert result.heavy_atom_count == 6
    assert result.method_id == NATIVE_REDOCK_RMSD_METHOD_ID
    assert result.literature_dois == NATIVE_REDOCK_RMSD_LITERATURE_DOIS
    assert "10.1186/s13321-019-0362-7" in result.literature_dois


def test_qualification_uses_selected_top_pose_not_best_generated_pose(
    tmp_path: Path,
) -> None:
    case = _make_case(
        tmp_path,
        pose_coords=(
            _translated(_NATIVE_COORDS, 6.0),
            _translated(_NATIVE_COORDS, 3.0),
        ),
    )

    result = _evaluate(case)

    assert result.status == "not_qualified"
    assert result.qualified is False
    assert result.top_ranked_pose_ordinal == 1
    assert result.top_ranked_rmsd_a == pytest.approx(3.0, abs=1e-3)
    assert result.best_generated_pose_ordinal == 2
    assert result.best_generated_rmsd_a == pytest.approx(0.0, abs=1e-3)


@pytest.mark.parametrize(
    ("mapping", "reason"),
    [
        ("", "meeko_smiles_idx_missing"),
        (
            "REMARK SMILES IDX 1 1 1 2 3 3 4 4 5 5 6 6\n",
            "meeko_smiles_idx_duplicate_smiles_atom",
        ),
        (
            "REMARK SMILES IDX 1 1 2 2 3 3 4 4 5 5\n",
            "meeko_mapping_heavy_smiles_coverage_invalid",
        ),
    ],
)
def test_prepared_meeko_mapping_must_be_one_to_one_and_complete(
    tmp_path: Path,
    mapping: str,
    reason: str,
) -> None:
    case = _make_case(tmp_path)
    case.prepared_pdbqt.write_bytes(
        _pdbqt_body(_translated(_NATIVE_COORDS, 3.0), mapping=mapping)
    )

    result = _evaluate(case)

    assert result.status == "invalid"
    assert result.failure_reason == reason


@pytest.mark.parametrize(
    ("payload", "reason"),
    [
        (
            b"MODEL 1\nMODEL 2\nENDMDL\nENDMDL\n",
            "pdbqt_model_nested",
        ),
        (
            b"MODEL 1\n" + _pdbqt_body(_translated(_NATIVE_COORDS, 3.0)),
            "pdbqt_model_missing_endmdl",
        ),
        (
            _model_block(2, _translated(_NATIVE_COORDS, 3.0), score=-8.0),
            "pdbqt_model_label_nonsequential:2:expected_1",
        ),
        (b"MODEL 1\nENDMDL\n", "pdbqt_model_empty:1"),
    ],
)
def test_explicit_model_parser_rejects_malformed_or_ambiguous_models(
    tmp_path: Path,
    payload: bytes,
    reason: str,
) -> None:
    case = _make_case(tmp_path)
    case.poses.write_bytes(payload)
    case.model_blocks = (payload,)
    _write_selection(case, selected_sha256="0" * 64)

    result = _evaluate(case)

    assert result.status == "invalid"
    assert result.failure_reason == reason


def test_frame_manifest_expected_hash_is_mandatory(tmp_path: Path) -> None:
    case = _make_case(tmp_path)

    result = evaluate_native_redock_rmsd(
        case.reference,
        case.graph,
        case.prepared_pdbqt,
        case.poses,
        source_receptor_path=case.source_receptor,
        prepared_receptor_path=case.prepared_receptor,
        coordinate_frame_manifest_path=case.frame,
        coordinate_frame_manifest_sha256="0" * 64,
        reference_extraction_manifest_path=case.extraction,
        pose_selection_manifest_path=case.selection,
        pose_selection_manifest_sha256=case.selection_sha256,
    )

    assert result.status == "invalid"
    assert result.failure_reason == "coordinate_frame_manifest_sha256_mismatch"


def test_frame_manifest_binds_receptor_content(tmp_path: Path) -> None:
    case = _make_case(tmp_path)
    case.source_receptor.write_text("TAMPERED SOURCE RECEPTOR\n", encoding="ascii")

    result = _evaluate(case)

    assert result.status == "invalid"
    assert result.failure_reason == ("coordinate_frame_source_receptor_sha256_mismatch")


def test_frame_manifest_binds_transform_content(tmp_path: Path) -> None:
    case = _make_case(tmp_path)
    _write_extraction_and_frame(case, transform_sha256="0" * 64)

    result = _evaluate(case)

    assert result.status == "invalid"
    assert result.failure_reason == "coordinate_transform_sha256_mismatch"


def test_pose_selection_manifest_binds_exact_selected_model(tmp_path: Path) -> None:
    case = _make_case(tmp_path)
    _write_selection(case, selected_sha256="0" * 64)

    result = _evaluate(case)

    assert result.status == "invalid"
    assert result.failure_reason == "pose_selection_model_sha256_mismatch"


def test_unverified_receptor_transform_is_evidence_pending(tmp_path: Path) -> None:
    case = _make_case(tmp_path)
    moved = ((3.0, 0.0, 0.0), (4.0, 0.0, 0.0), (3.0, 1.01, 0.0))
    case.prepared_receptor.write_text(
        "".join(
            _receptor_atom_line(index, name, 1, coord)
            for index, (name, coord) in enumerate(
                zip(("N", "CA", "C"), moved), start=1
            )
        ),
        encoding="ascii",
    )
    _write_extraction_and_frame(case)

    result = _evaluate(case)

    assert result.status == "evidence_pending"
    assert result.qualified is None
    assert result.coordinate_transform_verification_status == "evidence_pending"
    assert result.coordinate_transform_max_residual_a == pytest.approx(0.01)
    assert result.failure_reason is not None
    assert "coordinate_residual_exceeds" in result.failure_reason


def test_manifest_selected_pose_must_be_score_verified_top_ranked(
    tmp_path: Path,
) -> None:
    case = _make_case(
        tmp_path,
        pose_coords=(
            _translated(_NATIVE_COORDS, 3.0),
            _translated(_NATIVE_COORDS, 4.0),
        ),
        pose_scores=(-8.0, -7.0),
    )
    _write_selection(case, selected_ordinal=2, selected_label=2)

    result = _evaluate(case)

    assert result.status == "evidence_pending"
    assert result.qualified is None
    assert result.manifest_selected_pose_ordinal == 2
    assert result.top_ranked_pose_ordinal is None
    assert result.failure_reason == (
        "evidence_pending:pose_selection_not_verified_top_ranked"
    )


def test_every_competing_vina_model_requires_exact_score_evidence(
    tmp_path: Path,
) -> None:
    case = _make_case(
        tmp_path,
        pose_coords=(
            _translated(_NATIVE_COORDS, 3.0),
            _translated(_NATIVE_COORDS, 4.0),
        ),
    )
    second = case.model_blocks[1].replace(
        b"REMARK VINA RESULT: -7.000 0.000 0.000\n",
        b"",
    )
    case.model_blocks = (case.model_blocks[0], second)
    case.poses.write_bytes(b"".join(case.model_blocks))
    _write_selection(case)

    result = _evaluate(case)

    assert result.status == "evidence_pending"
    assert result.qualified is None
    assert result.top_ranked_pose_ordinal is None
    assert result.failure_reason == (
        "evidence_pending:vina_result_record_count_not_one:model_2:0"
    )


def test_raw_pdb_is_not_accepted_as_authoritative_ligand_topology(
    tmp_path: Path,
) -> None:
    case = _make_case(tmp_path)
    raw_reference = tmp_path / "reference.pdb"
    raw_reference.write_text("HETATM    1  C1  BEN A 100       0.000   0.000   0.000\n")
    case.reference = raw_reference
    _write_extraction_and_frame(case)

    result = _evaluate(case)

    assert result.status == "invalid"
    assert result.failure_reason == (
        "reference_ligand_authoritative_topology_format_unsupported"
    )


@pytest.mark.parametrize(
    ("kwargs", "reason"),
    [
        ({"threshold_a": True}, "qualification_threshold_invalid"),
        ({"threshold_a": 0.0}, "qualification_threshold_invalid"),
        ({"max_graph_mappings": True}, "graph_mapping_limit_invalid"),
        ({"max_graph_mappings": 1}, "graph_mapping_limit_exceeded"),
        (
            {"max_graph_mappings": 100_001},
            "graph_mapping_limit_exceeds_supported_max",
        ),
    ],
)
def test_numeric_policy_inputs_fail_closed_before_coercion(
    tmp_path: Path,
    kwargs: dict[str, Any],
    reason: str,
) -> None:
    case = _make_case(tmp_path)

    result = _evaluate(case, **kwargs)

    assert result.status == "invalid"
    assert result.failure_reason == reason


def test_pose_selection_boolean_model_index_is_rejected(tmp_path: Path) -> None:
    case = _make_case(tmp_path)
    _write_selection(case, selected_ordinal=True, selected_label=1)

    result = _evaluate(case)

    assert result.status == "invalid"
    assert result.failure_reason == "pose_selection_model_ordinal_invalid"


@pytest.mark.parametrize(
    ("atom_type", "reason"),
    [
        ("a", "pdbqt_atom_type_case_invalid:a"),
        ("G0", "pdbqt_pseudoatom_type:G0"),
        ("XX", "pdbqt_atom_type_unknown:XX"),
    ],
)
def test_unknown_case_invalid_and_pseudo_atom_types_are_rejected(
    tmp_path: Path,
    atom_type: str,
    reason: str,
) -> None:
    case = _make_case(tmp_path)
    case.prepared_pdbqt.write_bytes(
        _pdbqt_body(_translated(_NATIVE_COORDS, 3.0), atom_type=atom_type)
    )

    result = _evaluate(case)

    assert result.status == "invalid"
    assert result.failure_reason == reason


def test_unspecified_potential_stereochemistry_is_policy_pending(tmp_path: Path) -> None:
    case = _make_case(tmp_path)
    case.reference = tmp_path / "unspecified_stereo.sdf"
    _write_sdf(
        case.reference,
        smiles="CC(F)Cl",
        coords=((0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (2.0, 1.0, 0.0), (2.0, -1.0, 0.0)),
    )
    _write_extraction_and_frame(case)

    result = _evaluate(case)

    assert result.status == "evidence_pending"
    assert result.qualified is None
    assert result.failure_reason == (
        "policy_pending:unspecified_potential_stereochemistry:reference_ligand"
    )


def test_threshold_is_inclusive_at_two_point_five_angstrom(tmp_path: Path) -> None:
    case = _make_case(
        tmp_path,
        pose_coords=(_translated(_NATIVE_COORDS, 5.5),),
    )

    result = _evaluate(case)

    assert result.top_ranked_rmsd_a == pytest.approx(2.5, abs=1e-3)
    assert result.status == "qualified"
    assert result.qualified is True


def test_legacy_wrapper_refuses_unprovenanced_two_path_call() -> None:
    assert compute_redock_rmsd("reference.sdf", "poses.pdbqt") is None


def test_native_control_never_falls_back_to_ligand_fitted_rmsd(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference = tmp_path / "reference.sdf"
    reference.write_text("reference", encoding="utf-8")

    monkeypatch.setattr(
        docking_ligands, "compute_redock_rmsd", lambda *_args, **_kwargs: None
    )

    def reject_legacy(*_args: Any, **_kwargs: Any) -> float:
        raise AssertionError("legacy ligand-fitted RMSD must not qualify a control")

    monkeypatch.setattr(docking_ligands, "compute_rmsd", reject_legacy)

    assert not docking_ligands.validate_ligand(
        "native_control", "poses.pdbqt", str(reference), rmsd_thresh=2.5
    )


def test_native_control_default_cutoff_is_inclusive_two_point_five(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reference = tmp_path / "reference.sdf"
    reference.write_text("reference", encoding="utf-8")
    monkeypatch.setattr(
        docking_ligands,
        "compute_redock_rmsd",
        lambda *_args, **_kwargs: 2.5,
    )

    assert docking_ligands.validate_ligand(
        "native_control",
        "poses.pdbqt",
        str(reference),
    )
