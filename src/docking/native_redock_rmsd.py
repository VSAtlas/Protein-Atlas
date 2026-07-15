"""Content-bound, fail-closed native-ligand redocking RMSD.

The calculation is performed in the prepared-receptor coordinate frame without
independent ligand superposition. Meeko SMILES/SMILES-IDX remarks provide the
durable topology-to-coordinate map for every pose.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Sequence

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolAlign
from rdkit.Geometry import Point3D

from protein_prep.pdb_records import atom_name, pdbqt_line_xyz

NATIVE_REDOCK_RMSD_SCHEMA_VERSION = 2
NATIVE_REDOCK_RMSD_METHOD_ID = (
    "atlas.native_redock.meeko_mapped_in_place_rdkit_calcrms.v2"
)
NATIVE_REDOCK_RMSD_FORMULA = (
    "min_{pi in Iso(G_pose,G_ref)} "
    "sqrt((1/N) * sum_i ||T_source_to_prepared(x_ref_pi(i))-x_pose_i||^2)"
)
NATIVE_REDOCK_RMSD_THRESHOLD_A = 2.5
NATIVE_REDOCK_RMSD_LITERATURE_DOIS = ("10.1186/s13321-019-0362-7",)
DEFAULT_MAX_GRAPH_MAPPINGS = 10_000
MAX_GRAPH_MAPPINGS = 100_000
MANIFEST_SCHEMA_VERSION = 1
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_POSITIVE_INT_RE = re.compile(r"^[1-9][0-9]*$")
_MODEL_RE = re.compile(r"^MODEL[ \t]+([1-9][0-9]*)[ \t]*$")
_ENDMDL_RE = re.compile(r"^ENDMDL[ \t]*$")
_VINA_RESULT_RE = re.compile(
    r"^REMARK[ \t]+VINA[ \t]+RESULT:[ \t]+"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
    r"[ \t]+([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
    r"[ \t]+([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
    r"[ \t]*$"
)
_TIE_RULE = "model_ordinal_ascending"
RECEPTOR_TRANSFORM_SERIALIZATION_TOLERANCE_A = 0.002
RECEPTOR_TRANSFORM_VERIFICATION_METHOD = (
    "exact_common_pdb_atom_keys_all_coordinate_residuals_v1"
)
VINA_POSE_SCORE_VERIFICATION_METHOD = "vina_remark_vina_result_all_models_v1"
_ATOM_TYPE_ELEMENTS = {
    "H": "H",
    "HD": "H",
    "B": "B",
    "C": "C",
    "A": "C",
    "N": "N",
    "NA": "N",
    "OA": "O",
    "F": "F",
    "Mg": "MG",
    "Si": "SI",
    "P": "P",
    "S": "S",
    "SA": "S",
    "Cl": "CL",
    "Ca": "CA",
    "Mn": "MN",
    "Fe": "FE",
    "Zn": "ZN",
    "Br": "BR",
    "I": "I",
}


@dataclass(frozen=True)
class PoseRmsd:
    pose_ordinal: int
    model_label: int | None
    model_sha256: str
    selection_score: float
    rmsd_a: float


@dataclass(frozen=True)
class NativeRedockRmsdResult:
    schema_version: int
    method_id: str
    formula: str
    literature_dois: tuple[str, ...]
    rdkit_version: str
    atom_policy: str
    alignment_policy: str
    symmetry_policy: str
    bond_identity_policy: str
    stereo_policy: str
    symmetrize_conjugated_terminal_groups: bool
    qualification_threshold_a: float | None
    qualification_rule: str
    graph_mapping_limit: int | None
    status: str
    qualified: bool | None
    failure_reason: str | None
    reference_ligand_path: str
    prepared_ligand_graph_path: str
    prepared_ligand_pdbqt_path: str
    docked_poses_path: str
    reference_ligand_sha256: str | None
    prepared_ligand_graph_sha256: str | None
    prepared_ligand_pdbqt_sha256: str | None
    docked_poses_sha256: str | None
    source_receptor_path: str
    prepared_receptor_path: str
    source_receptor_sha256: str | None
    prepared_receptor_sha256: str | None
    coordinate_frame_manifest_path: str
    coordinate_frame_manifest_sha256: str | None
    reference_extraction_manifest_path: str
    reference_extraction_manifest_sha256: str | None
    pose_selection_manifest_path: str
    pose_selection_manifest_sha256: str | None
    source_to_prepared_transform: tuple[tuple[float, ...], ...] | None
    source_to_prepared_transform_sha256: str | None
    coordinate_transform_is_identity: bool | None
    coordinate_frame_relation: str
    coordinate_transform_verification_method: str
    coordinate_transform_verification_tolerance_a: float
    coordinate_transform_verification_status: str
    source_receptor_atom_count: int
    prepared_receptor_atom_count: int
    common_receptor_atom_count: int
    coordinate_transform_max_residual_a: float | None
    reference_extraction_method: str | None
    reference_ligand_id: str | None
    identity_method: str
    reference_identity: str | None
    prepared_identity: str | None
    meeko_smiles: str | None
    meeko_mapping_sha256: str | None
    heavy_atom_count: int | None
    graph_mapping_count: int
    prepared_atom_signature_sha256: str | None
    generated_pose_count: int
    evaluated_pose_count: int
    selection_engine: str | None
    selection_stage: str | None
    selection_score_source: str | None
    selection_score_value: float | None
    selection_score_direction: str | None
    selection_tie_break_rule: str | None
    selection_score_verification_method: str | None
    selection_score_verification_status: str
    manifest_selected_pose_ordinal: int | None
    manifest_selected_model_label: int | None
    manifest_selected_model_sha256: str | None
    top_ranked_pose_ordinal: int | None
    top_ranked_model_label: int | None
    top_ranked_model_sha256: str | None
    top_ranked_rmsd_a: float | None
    best_generated_pose_ordinal: int | None
    best_generated_rmsd_a: float | None
    pose_rmsds: tuple[PoseRmsd, ...]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


class _NativeRedockInputError(ValueError):
    pass


class _NativeRedockEvidencePending(ValueError):
    pass


@dataclass(frozen=True)
class _PdbqtAtom:
    name: str
    atom_type: str
    element: str
    coords: tuple[float, float, float]


@dataclass(frozen=True)
class _PdbqtModel:
    ordinal: int
    label: int | None
    raw_sha256: str
    lines: tuple[str, ...]
    atoms: tuple[_PdbqtAtom, ...]


@dataclass(frozen=True)
class _ReceptorAtom:
    key: tuple[str, str, str, str, str]
    coords: tuple[float, float, float]


@dataclass(frozen=True)
class _MeekoMapping:
    smiles: str
    pairs: tuple[tuple[int, int], ...]
    mapping_sha256: str
    smiles_mol: Chem.Mol
    heavy_mol: Chem.Mol
    original_to_heavy: dict[int, int]


@dataclass(frozen=True)
class _PoseSelection:
    engine: str
    stage: str
    score_source: str
    score_value: float
    score_direction: str
    tie_break_rule: str
    selected_model_ordinal: int
    selected_model_label: int | None
    selected_model_sha256: str


def sha256_file(path: str | Path) -> str | None:
    try:
        digest = hashlib.sha256()
        with Path(path).open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()
    except OSError:
        return None


def canonical_json_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _validate_sha256(value: Any, label: str) -> str:
    if not isinstance(value, str) or _SHA256_RE.fullmatch(value) is None:
        raise _NativeRedockInputError(f"{label}_sha256_invalid")
    return value


def _reject_duplicate_json_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise _NativeRedockInputError(f"manifest_duplicate_key:{key}")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> None:
    raise _NativeRedockInputError(f"manifest_nonfinite_number:{value}")


def _load_hashed_manifest(
    path: Path,
    expected_sha256: Any,
    label: str,
) -> tuple[dict[str, Any], str]:
    expected = _validate_sha256(expected_sha256, label)
    if not path.is_file():
        raise _NativeRedockInputError(f"{label}_manifest_missing")
    actual = sha256_file(path)
    if actual != expected:
        raise _NativeRedockInputError(f"{label}_manifest_sha256_mismatch")
    try:
        payload = json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_keys,
            parse_constant=_reject_json_constant,
        )
    except _NativeRedockInputError:
        raise
    except Exception as exc:
        raise _NativeRedockInputError(
            f"{label}_manifest_json_invalid:{type(exc).__name__}"
        ) from exc
    if not isinstance(payload, dict):
        raise _NativeRedockInputError(f"{label}_manifest_root_not_object")
    if type(payload.get("schema_version")) is not int:
        raise _NativeRedockInputError(f"{label}_manifest_schema_invalid")
    if payload["schema_version"] != MANIFEST_SCHEMA_VERSION:
        raise _NativeRedockInputError(f"{label}_manifest_schema_unsupported")
    return payload, actual


def _require_nonempty_string(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise _NativeRedockInputError(f"{label}_missing")
    return value.strip()


def _load_one_heavy_molecule(path: Path, label: str) -> Chem.Mol:
    suffix = path.suffix.lower()
    if suffix not in {".sdf", ".sd", ".mol", ".mol2"}:
        raise _NativeRedockInputError(
            f"{label}_authoritative_topology_format_unsupported"
        )
    try:
        if suffix in {".sdf", ".sd"}:
            records = [
                mol
                for mol in Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False)
                if mol is not None
            ]
            if len(records) != 1:
                raise _NativeRedockInputError(
                    f"{label}_record_count_not_one:{len(records)}"
                )
            mol = records[0]
        elif suffix == ".mol2":
            mol = Chem.MolFromMol2File(str(path), sanitize=True, removeHs=False)
        else:
            mol = Chem.MolFromMolFile(str(path), sanitize=True, removeHs=False)
    except _NativeRedockInputError:
        raise
    except Exception as exc:
        raise _NativeRedockInputError(
            f"{label}_load_failed:{type(exc).__name__}"
        ) from exc
    if mol is None:
        raise _NativeRedockInputError(f"{label}_load_failed")
    try:
        heavy = Chem.RemoveHs(mol, sanitize=True)
        Chem.AssignStereochemistry(heavy, cleanIt=True, force=True)
    except Exception as exc:
        raise _NativeRedockInputError(
            f"{label}_sanitization_failed:{type(exc).__name__}"
        ) from exc
    if heavy.GetNumAtoms() < 1:
        raise _NativeRedockInputError(f"{label}_has_no_heavy_atoms")
    _require_fully_specified_stereo(heavy, label)
    return heavy


def _require_fully_specified_stereo(mol: Chem.Mol, label: str) -> None:
    try:
        stereo = Chem.FindPotentialStereo(mol, cleanIt=True, flagPossible=True)
    except Exception as exc:
        raise _NativeRedockInputError(
            f"{label}_stereochemistry_audit_failed:{type(exc).__name__}"
        ) from exc
    if any(str(item.specified).endswith("Unspecified") for item in stereo):
        raise _NativeRedockEvidencePending(
            f"policy_pending:unspecified_potential_stereochemistry:{label}"
        )


def _canonical_identity(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(
        mol,
        canonical=True,
        isomericSmiles=True,
        allHsExplicit=False,
        allBondsExplicit=True,
    )


def _atom_type_element(atom_type: str) -> str:
    if atom_type == "W" or re.fullmatch(r"G[0-9]+", atom_type):
        raise _NativeRedockInputError(f"pdbqt_pseudoatom_type:{atom_type}")
    if re.fullmatch(r"CG[0-9]+", atom_type):
        return "C"
    element = _ATOM_TYPE_ELEMENTS.get(atom_type)
    if element is not None:
        return element
    if any(atom_type.lower() == known.lower() for known in _ATOM_TYPE_ELEMENTS):
        raise _NativeRedockInputError(f"pdbqt_atom_type_case_invalid:{atom_type}")
    if re.fullmatch(r"(?:cg|g)[0-9]+", atom_type, flags=re.IGNORECASE):
        raise _NativeRedockInputError(f"pdbqt_atom_type_case_invalid:{atom_type}")
    raise _NativeRedockInputError(f"pdbqt_atom_type_unknown:{atom_type}")


def _build_pdbqt_model(
    raw_lines: Sequence[bytes],
    ordinal: int,
    label: int | None,
) -> _PdbqtModel:
    raw = b"".join(raw_lines)
    text_lines: list[str] = []
    atoms: list[_PdbqtAtom] = []
    try:
        decoded = [line.decode("ascii") for line in raw_lines]
    except UnicodeDecodeError as exc:
        raise _NativeRedockInputError("pdbqt_non_ascii_content") from exc
    for line in decoded:
        text = line.rstrip("\r\n")
        text_lines.append(text)
        if not text.startswith(("ATOM", "HETATM")):
            continue
        xyz = pdbqt_line_xyz(text)
        if xyz is None:
            raise _NativeRedockInputError(
                f"pdbqt_atom_coordinates_malformed:model_{ordinal}"
            )
        if not all(math.isfinite(value) for value in xyz):
            raise _NativeRedockInputError(
                f"pdbqt_atom_coordinates_nonfinite:model_{ordinal}"
            )
        fields = text.split()
        if not fields:
            raise _NativeRedockInputError(f"pdbqt_atom_type_missing:model_{ordinal}")
        atom_type = fields[-1]
        element = _atom_type_element(atom_type)
        name = atom_name(text).strip()
        if not name:
            raise _NativeRedockInputError(f"pdbqt_atom_name_missing:model_{ordinal}")
        atoms.append(
            _PdbqtAtom(
                name=name,
                atom_type=atom_type,
                element=element,
                coords=(float(xyz[0]), float(xyz[1]), float(xyz[2])),
            )
        )
    if not atoms:
        raise _NativeRedockInputError(f"pdbqt_model_empty:{ordinal}")
    return _PdbqtModel(
        ordinal=ordinal,
        label=label,
        raw_sha256=hashlib.sha256(raw).hexdigest(),
        lines=tuple(text_lines),
        atoms=tuple(atoms),
    )


def _parse_pdbqt_models(path: Path) -> tuple[_PdbqtModel, ...]:
    try:
        raw_lines = path.read_bytes().splitlines(keepends=True)
    except OSError as exc:
        raise _NativeRedockInputError("pdbqt_read_failed") from exc
    if not raw_lines:
        raise _NativeRedockInputError("pdbqt_file_empty")
    try:
        texts = [line.decode("ascii").rstrip("\r\n") for line in raw_lines]
    except UnicodeDecodeError as exc:
        raise _NativeRedockInputError("pdbqt_non_ascii_content") from exc
    explicit = any(text.startswith("MODEL") for text in texts)
    if not explicit:
        if any(text.startswith("ENDMDL") for text in texts):
            raise _NativeRedockInputError("pdbqt_unexpected_endmdl")
        return (_build_pdbqt_model(raw_lines, 1, None),)

    models: list[_PdbqtModel] = []
    current: list[bytes] | None = None
    current_label: int | None = None
    for raw_line, text in zip(raw_lines, texts):
        if text.startswith("MODEL"):
            match = _MODEL_RE.fullmatch(text)
            if match is None:
                raise _NativeRedockInputError("pdbqt_model_record_malformed")
            if current is not None:
                raise _NativeRedockInputError("pdbqt_model_nested")
            label = int(match.group(1))
            expected = len(models) + 1
            if label != expected:
                raise _NativeRedockInputError(
                    f"pdbqt_model_label_nonsequential:{label}:expected_{expected}"
                )
            current = [raw_line]
            current_label = label
            continue
        if text.startswith("ENDMDL"):
            if _ENDMDL_RE.fullmatch(text) is None:
                raise _NativeRedockInputError("pdbqt_endmdl_record_malformed")
            if current is None:
                raise _NativeRedockInputError("pdbqt_unexpected_endmdl")
            current.append(raw_line)
            models.append(_build_pdbqt_model(current, len(models) + 1, current_label))
            current = None
            current_label = None
            continue
        if current is None:
            if text.strip():
                raise _NativeRedockInputError("pdbqt_content_outside_model")
            continue
        current.append(raw_line)
    if current is not None:
        raise _NativeRedockInputError("pdbqt_model_missing_endmdl")
    if not models:
        raise _NativeRedockInputError("pdbqt_models_missing")
    return tuple(models)


def _parse_positive_index(token: str, label: str) -> int:
    if _POSITIVE_INT_RE.fullmatch(token) is None:
        raise _NativeRedockInputError(f"{label}_invalid:{token}")
    return int(token)


def _parse_meeko_mapping(model: _PdbqtModel) -> _MeekoMapping:
    smiles_values: list[str] = []
    index_tokens: list[str] = []
    for line in model.lines:
        if line.startswith("REMARK SMILES IDX"):
            fields = line.split()
            if fields[:3] != ["REMARK", "SMILES", "IDX"]:
                raise _NativeRedockInputError("meeko_smiles_idx_record_malformed")
            index_tokens.extend(fields[3:])
        elif line.startswith("REMARK SMILES"):
            fields = line.split()
            if len(fields) != 3 or fields[:2] != ["REMARK", "SMILES"]:
                raise _NativeRedockInputError("meeko_smiles_record_malformed")
            smiles_values.append(fields[2])
    if len(smiles_values) != 1:
        raise _NativeRedockInputError(
            f"meeko_smiles_record_count_not_one:{len(smiles_values)}"
        )
    if not index_tokens:
        raise _NativeRedockInputError("meeko_smiles_idx_missing")
    if len(index_tokens) % 2:
        raise _NativeRedockInputError("meeko_smiles_idx_odd_length")
    pairs = tuple(
        (
            _parse_positive_index(index_tokens[offset], "meeko_smiles_atom_index"),
            _parse_positive_index(
                index_tokens[offset + 1], "meeko_pdbqt_coordinate_index"
            ),
        )
        for offset in range(0, len(index_tokens), 2)
    )
    smiles_indices = [pair[0] for pair in pairs]
    coordinate_indices = [pair[1] for pair in pairs]
    if len(smiles_indices) != len(set(smiles_indices)):
        raise _NativeRedockInputError("meeko_smiles_idx_duplicate_smiles_atom")
    if len(coordinate_indices) != len(set(coordinate_indices)):
        raise _NativeRedockInputError("meeko_smiles_idx_duplicate_coordinate")

    smiles = smiles_values[0]
    try:
        smiles_mol = Chem.MolFromSmiles(smiles, sanitize=True)
    except Exception as exc:
        raise _NativeRedockInputError(
            f"meeko_smiles_parse_failed:{type(exc).__name__}"
        ) from exc
    if smiles_mol is None:
        raise _NativeRedockInputError("meeko_smiles_parse_failed")
    Chem.AssignStereochemistry(smiles_mol, cleanIt=True, force=True)
    _require_fully_specified_stereo(smiles_mol, "meeko_smiles")
    heavy_original = [
        atom.GetIdx() for atom in smiles_mol.GetAtoms() if atom.GetAtomicNum() != 1
    ]
    original_to_heavy = {
        original_index: heavy_index
        for heavy_index, original_index in enumerate(heavy_original)
    }
    heavy_mol = Chem.RemoveHs(smiles_mol, sanitize=True)

    mapped_heavy_smiles: set[int] = set()
    mapped_heavy_coords: set[int] = set()
    for smiles_index, coordinate_index in pairs:
        if smiles_index > smiles_mol.GetNumAtoms():
            raise _NativeRedockInputError("meeko_smiles_idx_out_of_range")
        if coordinate_index > len(model.atoms):
            raise _NativeRedockInputError("meeko_pdbqt_coordinate_idx_out_of_range")
        graph_atom = smiles_mol.GetAtomWithIdx(smiles_index - 1)
        pdbqt_atom = model.atoms[coordinate_index - 1]
        graph_element = graph_atom.GetSymbol().upper()
        if graph_element != pdbqt_atom.element:
            raise _NativeRedockInputError(
                f"meeko_mapping_element_mismatch:{smiles_index}:{coordinate_index}"
            )
        if graph_atom.GetAtomicNum() != 1:
            mapped_heavy_smiles.add(smiles_index)
            mapped_heavy_coords.add(coordinate_index)
    expected_heavy_smiles = {index + 1 for index in heavy_original}
    expected_heavy_coords = {
        index + 1 for index, atom in enumerate(model.atoms) if atom.element != "H"
    }
    if mapped_heavy_smiles != expected_heavy_smiles:
        raise _NativeRedockInputError("meeko_mapping_heavy_smiles_coverage_invalid")
    if mapped_heavy_coords != expected_heavy_coords:
        raise _NativeRedockInputError("meeko_mapping_heavy_coordinate_coverage_invalid")
    mapping_payload = {"smiles": smiles, "pairs": [list(pair) for pair in pairs]}
    return _MeekoMapping(
        smiles=smiles,
        pairs=pairs,
        mapping_sha256=canonical_json_sha256(mapping_payload),
        smiles_mol=smiles_mol,
        heavy_mol=heavy_mol,
        original_to_heavy=original_to_heavy,
    )


def _atom_signature_sha256(model: _PdbqtModel) -> str:
    signature = [[atom.name, atom.atom_type, atom.element] for atom in model.atoms]
    return canonical_json_sha256(signature)


def _pose_from_meeko_mapping(
    model: _PdbqtModel,
    mapping: _MeekoMapping,
) -> Chem.Mol:
    coords: list[tuple[float, float, float] | None] = [
        None
    ] * mapping.heavy_mol.GetNumAtoms()
    for smiles_index, coordinate_index in mapping.pairs:
        original_index = smiles_index - 1
        heavy_index = mapping.original_to_heavy.get(original_index)
        if heavy_index is None:
            continue
        coords[heavy_index] = model.atoms[coordinate_index - 1].coords
    if any(coord is None for coord in coords):
        raise _NativeRedockInputError("meeko_pose_reconstruction_incomplete")
    pose = Chem.Mol(mapping.heavy_mol)
    pose.RemoveAllConformers()
    conf = Chem.Conformer(pose.GetNumAtoms())
    conf.Set3D(True)
    for atom_index, coord in enumerate(coords):
        if coord is None:
            raise _NativeRedockInputError("meeko_pose_reconstruction_incomplete")
        conf.SetAtomPosition(atom_index, Point3D(*coord))
    pose.AddConformer(conf, assignId=True)
    return pose


def _validate_reference_coordinates(mol: Chem.Mol) -> None:
    if mol.GetNumConformers() != 1:
        raise _NativeRedockInputError(
            f"reference_conformer_count_not_one:{mol.GetNumConformers()}"
        )
    conf = mol.GetConformer()
    if not conf.Is3D():
        raise _NativeRedockInputError("reference_conformer_not_3d")
    for atom_index in range(mol.GetNumAtoms()):
        point = conf.GetAtomPosition(atom_index)
        if not all(math.isfinite(value) for value in (point.x, point.y, point.z)):
            raise _NativeRedockInputError("reference_coordinates_nonfinite")


def _parse_rigid_transform(
    value: Any,
) -> tuple[tuple[float, ...], ...]:
    if not isinstance(value, list) or len(value) != 4:
        raise _NativeRedockInputError("coordinate_transform_shape_invalid")
    rows: list[tuple[float, ...]] = []
    for row in value:
        if not isinstance(row, list) or len(row) != 4:
            raise _NativeRedockInputError("coordinate_transform_shape_invalid")
        parsed: list[float] = []
        for item in row:
            if isinstance(item, bool) or not isinstance(item, (int, float)):
                raise _NativeRedockInputError("coordinate_transform_value_invalid")
            number = float(item)
            if not math.isfinite(number):
                raise _NativeRedockInputError("coordinate_transform_value_nonfinite")
            parsed.append(number)
        rows.append(tuple(parsed))
    matrix = tuple(rows)
    tolerance = 1e-6
    if any(
        abs(matrix[3][index] - expected) > tolerance
        for index, expected in enumerate((0.0, 0.0, 0.0, 1.0))
    ):
        raise _NativeRedockInputError("coordinate_transform_not_affine")
    rotation = [row[:3] for row in matrix[:3]]
    for left in range(3):
        for right in range(3):
            dot = sum(rotation[left][axis] * rotation[right][axis] for axis in range(3))
            expected = 1.0 if left == right else 0.0
            if abs(dot - expected) > tolerance:
                raise _NativeRedockInputError("coordinate_transform_not_rigid")
    determinant = (
        rotation[0][0]
        * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
        - rotation[0][1]
        * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
        + rotation[0][2]
        * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0])
    )
    if abs(determinant - 1.0) > tolerance:
        raise _NativeRedockInputError("coordinate_transform_not_proper_rotation")
    return matrix


def _transform_is_identity(matrix: Sequence[Sequence[float]]) -> bool:
    return all(
        abs(matrix[row][column] - (1.0 if row == column else 0.0)) <= 1e-12
        for row in range(4)
        for column in range(4)
    )


def _transform_reference(
    reference: Chem.Mol,
    matrix: Sequence[Sequence[float]],
) -> Chem.Mol:
    transformed = Chem.Mol(reference)
    conf = transformed.GetConformer()
    for atom_index in range(transformed.GetNumAtoms()):
        point = conf.GetAtomPosition(atom_index)
        vector = (point.x, point.y, point.z, 1.0)
        output = tuple(
            sum(matrix[row][column] * vector[column] for column in range(4))
            for row in range(4)
        )
        if abs(output[3] - 1.0) > 1e-6:
            raise _NativeRedockInputError("coordinate_transform_homogeneous_w_invalid")
        if not all(math.isfinite(value) for value in output):
            raise _NativeRedockInputError("coordinate_transform_output_nonfinite")
        conf.SetAtomPosition(atom_index, Point3D(output[0], output[1], output[2]))
    return transformed


def _load_receptor_atoms(
    path: Path,
    label: str,
) -> dict[tuple[str, str, str, str, str], _ReceptorAtom]:
    atoms: dict[tuple[str, str, str, str, str], _ReceptorAtom] = {}
    try:
        with path.open("r", encoding="ascii") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.startswith("ATOM"):
                    continue
                if len(line) < 54:
                    raise _NativeRedockEvidencePending(
                        f"evidence_pending:{label}_atom_record_short:{line_number}"
                    )
                xyz = pdbqt_line_xyz(line)
                if xyz is None or not all(math.isfinite(value) for value in xyz):
                    raise _NativeRedockEvidencePending(
                        f"evidence_pending:{label}_atom_coordinates_invalid:{line_number}"
                    )
                key = (
                    line[12:16].strip(),
                    line[16:17].strip(),
                    line[17:20].strip(),
                    line[21:22].strip(),
                    line[22:27].strip(),
                )
                if not key[0] or not key[2] or not key[4]:
                    raise _NativeRedockEvidencePending(
                        f"evidence_pending:{label}_atom_key_incomplete:{line_number}"
                    )
                if key in atoms:
                    raise _NativeRedockEvidencePending(
                        f"evidence_pending:{label}_atom_key_duplicate:{line_number}"
                    )
                atoms[key] = _ReceptorAtom(
                    key=key,
                    coords=(float(xyz[0]), float(xyz[1]), float(xyz[2])),
                )
    except (OSError, UnicodeError) as exc:
        raise _NativeRedockEvidencePending(
            f"evidence_pending:{label}_receptor_read_failed:{type(exc).__name__}"
        ) from exc
    if not atoms:
        raise _NativeRedockEvidencePending(
            f"evidence_pending:{label}_protein_atom_records_missing"
        )
    return atoms


def _has_noncollinear_points(
    points: Sequence[tuple[float, float, float]],
    tolerance_a: float,
) -> bool:
    if len(points) < 3:
        return False
    origin = points[0]
    for second in points[1:]:
        first_vector = tuple(second[axis] - origin[axis] for axis in range(3))
        first_norm = math.sqrt(sum(value * value for value in first_vector))
        if first_norm <= tolerance_a:
            continue
        for third in points[1:]:
            second_vector = tuple(third[axis] - origin[axis] for axis in range(3))
            cross = (
                first_vector[1] * second_vector[2]
                - first_vector[2] * second_vector[1],
                first_vector[2] * second_vector[0]
                - first_vector[0] * second_vector[2],
                first_vector[0] * second_vector[1]
                - first_vector[1] * second_vector[0],
            )
            if math.sqrt(sum(value * value for value in cross)) > tolerance_a**2:
                return True
    return False


def _receptor_transform_verification(
    *,
    source_path: Path,
    prepared_path: Path,
    matrix: Sequence[Sequence[float]],
    tolerance_a: float,
) -> tuple[dict[str, Any], str | None]:
    source_atoms = _load_receptor_atoms(source_path, "source_receptor")
    prepared_atoms = _load_receptor_atoms(prepared_path, "prepared_receptor")
    common_keys = sorted(set(source_atoms) & set(prepared_atoms))
    evidence: dict[str, Any] = {
        "source_receptor_atom_count": len(source_atoms),
        "prepared_receptor_atom_count": len(prepared_atoms),
        "common_receptor_atom_count": len(common_keys),
        "coordinate_transform_max_residual_a": None,
    }
    if len(common_keys) < 3:
        return evidence, "evidence_pending:receptor_transform_common_atoms_fewer_than_three"
    source_points = [source_atoms[key].coords for key in common_keys]
    if not _has_noncollinear_points(source_points, tolerance_a):
        return evidence, "evidence_pending:receptor_transform_common_atoms_collinear"
    residuals: list[float] = []
    for key in common_keys:
        source = source_atoms[key].coords
        prepared = prepared_atoms[key].coords
        transformed = tuple(
            sum(matrix[row][column] * (*source, 1.0)[column] for column in range(4))
            for row in range(3)
        )
        residuals.append(
            math.sqrt(
                sum(
                    (transformed[axis] - prepared[axis]) ** 2
                    for axis in range(3)
                )
            )
        )
    maximum = max(residuals)
    evidence["coordinate_transform_max_residual_a"] = maximum
    if maximum > tolerance_a:
        return (
            evidence,
            "evidence_pending:receptor_transform_coordinate_residual_exceeds_"
            f"serialization_tolerance:{maximum:.6g}>{tolerance_a:.6g}",
        )
    return evidence, None


def _verify_frame_and_extraction_manifests(
    *,
    frame_path: Path,
    frame_expected_sha256: Any,
    extraction_path: Path,
    source_receptor_sha256: str,
    prepared_receptor_sha256: str,
    reference_ligand_sha256: str,
) -> tuple[
    tuple[tuple[float, ...], ...],
    str,
    bool,
    str,
    str,
    str,
]:
    frame, _ = _load_hashed_manifest(
        frame_path, frame_expected_sha256, "coordinate_frame"
    )
    for key, actual in (
        ("source_receptor_sha256", source_receptor_sha256),
        ("prepared_receptor_sha256", prepared_receptor_sha256),
    ):
        declared = _validate_sha256(frame.get(key), f"coordinate_frame_{key}")
        if declared != actual:
            raise _NativeRedockInputError(f"coordinate_frame_{key}_mismatch")
    extraction_expected = _validate_sha256(
        frame.get("reference_extraction_manifest_sha256"),
        "coordinate_frame_reference_extraction_manifest",
    )
    extraction, extraction_actual = _load_hashed_manifest(
        extraction_path,
        extraction_expected,
        "reference_extraction",
    )
    extraction_source = _validate_sha256(
        extraction.get("source_receptor_sha256"),
        "reference_extraction_source_receptor",
    )
    if extraction_source != source_receptor_sha256:
        raise _NativeRedockInputError(
            "reference_extraction_source_receptor_sha256_mismatch"
        )
    extraction_ligand = _validate_sha256(
        extraction.get("reference_ligand_sha256"),
        "reference_extraction_reference_ligand",
    )
    if extraction_ligand != reference_ligand_sha256:
        raise _NativeRedockInputError(
            "reference_extraction_reference_ligand_sha256_mismatch"
        )
    extraction_method = _require_nonempty_string(
        extraction.get("extraction_method"),
        "reference_extraction_method",
    )
    ligand_id = _require_nonempty_string(
        extraction.get("ligand_id"),
        "reference_extraction_ligand_id",
    )
    raw_transform = frame.get("source_to_prepared_transform")
    matrix = _parse_rigid_transform(raw_transform)
    transform_expected = _validate_sha256(
        frame.get("source_to_prepared_transform_sha256"),
        "coordinate_frame_transform",
    )
    try:
        transform_actual = canonical_json_sha256(raw_transform)
    except (TypeError, ValueError) as exc:
        raise _NativeRedockInputError("coordinate_transform_hash_failed") from exc
    if transform_actual != transform_expected:
        raise _NativeRedockInputError("coordinate_transform_sha256_mismatch")
    return (
        matrix,
        transform_actual,
        _transform_is_identity(matrix),
        extraction_actual,
        extraction_method,
        ligand_id,
    )


def _parse_pose_selection_manifest(
    *,
    path: Path,
    expected_sha256: Any,
    docked_poses_sha256: str,
) -> tuple[_PoseSelection, str]:
    payload, actual_sha256 = _load_hashed_manifest(
        path, expected_sha256, "pose_selection"
    )
    declared_poses = _validate_sha256(
        payload.get("docked_poses_sha256"), "pose_selection_docked_poses"
    )
    if declared_poses != docked_poses_sha256:
        raise _NativeRedockInputError("pose_selection_docked_poses_sha256_mismatch")
    engine = _require_nonempty_string(payload.get("engine"), "pose_selection_engine")
    stage = _require_nonempty_string(payload.get("stage"), "pose_selection_stage")
    score_source = _require_nonempty_string(
        payload.get("score_source"), "pose_selection_score_source"
    )
    score_raw = payload.get("score_value")
    if isinstance(score_raw, bool) or not isinstance(score_raw, (int, float)):
        raise _NativeRedockInputError("pose_selection_score_value_invalid")
    score_value = float(score_raw)
    if not math.isfinite(score_value):
        raise _NativeRedockInputError("pose_selection_score_value_nonfinite")
    score_direction = payload.get("score_direction")
    if score_direction not in {"lower_is_better", "higher_is_better"}:
        raise _NativeRedockInputError("pose_selection_score_direction_invalid")
    tie_break_rule = payload.get("tie_break_rule")
    if tie_break_rule != _TIE_RULE:
        raise _NativeRedockInputError("pose_selection_tie_break_rule_unsupported")
    ordinal = payload.get("selected_model_ordinal")
    if type(ordinal) is not int or ordinal < 1:
        raise _NativeRedockInputError("pose_selection_model_ordinal_invalid")
    label = payload.get("selected_model_label")
    if label is not None and (type(label) is not int or label < 1):
        raise _NativeRedockInputError("pose_selection_model_label_invalid")
    model_sha256 = _validate_sha256(
        payload.get("selected_model_sha256"), "pose_selection_model"
    )
    return (
        _PoseSelection(
            engine=engine,
            stage=stage,
            score_source=score_source,
            score_value=score_value,
            score_direction=score_direction,
            tie_break_rule=tie_break_rule,
            selected_model_ordinal=ordinal,
            selected_model_label=label,
            selected_model_sha256=model_sha256,
        ),
        actual_sha256,
    )


def _vina_score_from_model(model: _PdbqtModel) -> float:
    records = [
        line for line in model.lines if line.startswith("REMARK VINA RESULT")
    ]
    if len(records) != 1:
        raise _NativeRedockEvidencePending(
            "evidence_pending:vina_result_record_count_not_one:"
            f"model_{model.ordinal}:{len(records)}"
        )
    match = _VINA_RESULT_RE.fullmatch(records[0])
    if match is None:
        raise _NativeRedockEvidencePending(
            f"evidence_pending:vina_result_record_malformed:model_{model.ordinal}"
        )
    values = tuple(float(token) for token in match.groups())
    if not all(math.isfinite(value) for value in values):
        raise _NativeRedockEvidencePending(
            f"evidence_pending:vina_result_nonfinite:model_{model.ordinal}"
        )
    return values[0]


def _verify_vina_pose_score_order(
    selection: _PoseSelection,
    models: Sequence[_PdbqtModel],
) -> tuple[float, ...]:
    if (
        selection.engine.strip().lower() != "vina"
        or selection.score_source != "vina_affinity_kcal_mol"
        or selection.score_direction != "lower_is_better"
    ):
        raise _NativeRedockEvidencePending(
            "evidence_pending:pose_score_source_not_independently_supported"
        )
    scores = tuple(_vina_score_from_model(model) for model in models)
    selected_index = selection.selected_model_ordinal - 1
    if selected_index >= len(scores):
        raise _NativeRedockInputError("pose_selection_model_ordinal_out_of_range")
    if selection.score_value != scores[selected_index]:
        raise _NativeRedockEvidencePending(
            "evidence_pending:pose_selection_score_value_mismatch"
        )
    expected = min(
        range(len(scores)),
        key=lambda index: (scores[index], models[index].ordinal),
    )
    if selected_index != expected:
        raise _NativeRedockEvidencePending(
            "evidence_pending:pose_selection_not_verified_top_ranked"
        )
    return scores


def _graph_maps(
    prepared_graph: Chem.Mol,
    reference: Chem.Mol,
    limit: int,
) -> tuple[tuple[tuple[int, int], ...], ...]:
    matches = reference.GetSubstructMatches(
        prepared_graph,
        uniquify=False,
        useChirality=True,
        maxMatches=limit + 1,
    )
    if len(matches) > limit:
        raise _NativeRedockInputError("graph_mapping_limit_exceeded")
    if not matches:
        raise _NativeRedockInputError("no_exact_graph_isomorphism")
    return tuple(
        tuple((pose_index, ref_index) for pose_index, ref_index in enumerate(match))
        for match in matches
    )


def evaluate_native_redock_rmsd(
    reference_ligand_path: str | Path,
    prepared_ligand_graph_path: str | Path,
    prepared_ligand_pdbqt_path: str | Path,
    docked_poses_path: str | Path,
    *,
    source_receptor_path: str | Path,
    prepared_receptor_path: str | Path,
    coordinate_frame_manifest_path: str | Path,
    coordinate_frame_manifest_sha256: str,
    reference_extraction_manifest_path: str | Path,
    pose_selection_manifest_path: str | Path,
    pose_selection_manifest_sha256: str,
    threshold_a: float = NATIVE_REDOCK_RMSD_THRESHOLD_A,
    max_graph_mappings: int = DEFAULT_MAX_GRAPH_MAPPINGS,
) -> NativeRedockRmsdResult:
    """Evaluate a native redock only when every content-bound contract verifies."""

    reference_path = Path(reference_ligand_path)
    graph_path = Path(prepared_ligand_graph_path)
    prepared_pdbqt_path = Path(prepared_ligand_pdbqt_path)
    poses_path = Path(docked_poses_path)
    source_receptor = Path(source_receptor_path)
    prepared_receptor = Path(prepared_receptor_path)
    frame_manifest = Path(coordinate_frame_manifest_path)
    extraction_manifest = Path(reference_extraction_manifest_path)
    selection_manifest = Path(pose_selection_manifest_path)
    threshold_value = (
        float(threshold_a)
        if not isinstance(threshold_a, bool)
        and isinstance(threshold_a, (int, float))
        and math.isfinite(float(threshold_a))
        else None
    )
    graph_limit_value = max_graph_mappings if type(max_graph_mappings) is int else None
    base: dict[str, Any] = {
        "schema_version": NATIVE_REDOCK_RMSD_SCHEMA_VERSION,
        "method_id": NATIVE_REDOCK_RMSD_METHOD_ID,
        "formula": NATIVE_REDOCK_RMSD_FORMULA,
        "literature_dois": NATIVE_REDOCK_RMSD_LITERATURE_DOIS,
        "rdkit_version": str(rdBase.rdkitVersion),
        "atom_policy": "meeko_mapped_heavy_atoms_only",
        "alignment_policy": (
            "common_receptor_atom_verified_transform_no_ligand_fit"
        ),
        "symmetry_policy": "all_exact_graph_isomorphisms_with_chirality",
        "bond_identity_policy": "exact_bond_orders_and_aromaticity_required",
        "stereo_policy": "unspecified_potential_stereochemistry_is_policy_pending",
        "symmetrize_conjugated_terminal_groups": False,
        "qualification_threshold_a": threshold_value,
        "qualification_rule": (
            "score_verified_top_ranked_rmsd_a <= qualification_threshold_a"
        ),
        "graph_mapping_limit": graph_limit_value,
        "reference_ligand_path": str(reference_path),
        "prepared_ligand_graph_path": str(graph_path),
        "prepared_ligand_pdbqt_path": str(prepared_pdbqt_path),
        "docked_poses_path": str(poses_path),
        "reference_ligand_sha256": sha256_file(reference_path),
        "prepared_ligand_graph_sha256": sha256_file(graph_path),
        "prepared_ligand_pdbqt_sha256": sha256_file(prepared_pdbqt_path),
        "docked_poses_sha256": sha256_file(poses_path),
        "source_receptor_path": str(source_receptor),
        "prepared_receptor_path": str(prepared_receptor),
        "source_receptor_sha256": sha256_file(source_receptor),
        "prepared_receptor_sha256": sha256_file(prepared_receptor),
        "coordinate_frame_manifest_path": str(frame_manifest),
        "coordinate_frame_manifest_sha256": sha256_file(frame_manifest),
        "reference_extraction_manifest_path": str(extraction_manifest),
        "reference_extraction_manifest_sha256": sha256_file(extraction_manifest),
        "pose_selection_manifest_path": str(selection_manifest),
        "pose_selection_manifest_sha256": sha256_file(selection_manifest),
        "source_to_prepared_transform": None,
        "source_to_prepared_transform_sha256": None,
        "coordinate_transform_is_identity": None,
        "coordinate_frame_relation": "unverified",
        "coordinate_transform_verification_method": (
            RECEPTOR_TRANSFORM_VERIFICATION_METHOD
        ),
        "coordinate_transform_verification_tolerance_a": (
            RECEPTOR_TRANSFORM_SERIALIZATION_TOLERANCE_A
        ),
        "coordinate_transform_verification_status": "evidence_pending",
        "source_receptor_atom_count": 0,
        "prepared_receptor_atom_count": 0,
        "common_receptor_atom_count": 0,
        "coordinate_transform_max_residual_a": None,
        "reference_extraction_method": None,
        "reference_ligand_id": None,
        "identity_method": "exact_graph_plus_meeko_smiles_idx_one_to_one_mapping",
        "reference_identity": None,
        "prepared_identity": None,
        "meeko_smiles": None,
        "meeko_mapping_sha256": None,
        "heavy_atom_count": None,
        "graph_mapping_count": 0,
        "prepared_atom_signature_sha256": None,
        "generated_pose_count": 0,
        "evaluated_pose_count": 0,
        "selection_engine": None,
        "selection_stage": None,
        "selection_score_source": None,
        "selection_score_value": None,
        "selection_score_direction": None,
        "selection_tie_break_rule": None,
        "selection_score_verification_method": None,
        "selection_score_verification_status": "evidence_pending",
        "manifest_selected_pose_ordinal": None,
        "manifest_selected_model_label": None,
        "manifest_selected_model_sha256": None,
        "top_ranked_pose_ordinal": None,
        "top_ranked_model_label": None,
        "top_ranked_model_sha256": None,
        "top_ranked_rmsd_a": None,
        "best_generated_pose_ordinal": None,
        "best_generated_rmsd_a": None,
        "pose_rmsds": (),
    }

    def finish(
        status: str,
        *,
        qualified: bool | None,
        failure_reason: str | None,
        **updates: Any,
    ) -> NativeRedockRmsdResult:
        return NativeRedockRmsdResult(
            **(base | updates),
            status=status,
            qualified=qualified,
            failure_reason=failure_reason,
        )

    try:
        if (
            isinstance(threshold_a, bool)
            or not isinstance(threshold_a, (int, float))
            or not math.isfinite(float(threshold_a))
            or float(threshold_a) <= 0
        ):
            raise _NativeRedockInputError("qualification_threshold_invalid")
        if type(max_graph_mappings) is not int or max_graph_mappings < 1:
            raise _NativeRedockInputError("graph_mapping_limit_invalid")
        if max_graph_mappings > MAX_GRAPH_MAPPINGS:
            raise _NativeRedockInputError("graph_mapping_limit_exceeds_supported_max")
        required_paths = (
            ("reference_ligand", reference_path),
            ("prepared_ligand_graph", graph_path),
            ("prepared_ligand_pdbqt", prepared_pdbqt_path),
            ("docked_poses", poses_path),
            ("source_receptor", source_receptor),
            ("prepared_receptor", prepared_receptor),
        )
        for label, path in required_paths:
            if not path.is_file():
                raise _NativeRedockInputError(f"{label}_file_missing")
        source_sha = base["source_receptor_sha256"]
        prepared_receptor_sha = base["prepared_receptor_sha256"]
        reference_sha = base["reference_ligand_sha256"]
        poses_sha = base["docked_poses_sha256"]
        if not all(
            isinstance(value, str)
            for value in (source_sha, prepared_receptor_sha, reference_sha, poses_sha)
        ):
            raise _NativeRedockInputError("required_content_sha256_unavailable")

        (
            transform,
            transform_sha,
            transform_is_identity,
            extraction_sha,
            extraction_method,
            ligand_id,
        ) = _verify_frame_and_extraction_manifests(
            frame_path=frame_manifest,
            frame_expected_sha256=coordinate_frame_manifest_sha256,
            extraction_path=extraction_manifest,
            source_receptor_sha256=source_sha,
            prepared_receptor_sha256=prepared_receptor_sha,
            reference_ligand_sha256=reference_sha,
        )
        base.update(
            {
                "source_to_prepared_transform": transform,
                "source_to_prepared_transform_sha256": transform_sha,
                "coordinate_transform_is_identity": transform_is_identity,
                "coordinate_frame_relation": (
                    "content_bound_declared_transform_pending_independent_verification"
                ),
                "reference_extraction_manifest_sha256": extraction_sha,
                "reference_extraction_method": extraction_method,
                "reference_ligand_id": ligand_id,
            }
        )
        transform_evidence, transform_pending_reason = (
            _receptor_transform_verification(
                source_path=source_receptor,
                prepared_path=prepared_receptor,
                matrix=transform,
                tolerance_a=RECEPTOR_TRANSFORM_SERIALIZATION_TOLERANCE_A,
            )
        )
        base.update(transform_evidence)
        if transform_pending_reason is not None:
            raise _NativeRedockEvidencePending(transform_pending_reason)
        base.update(
            {
                "coordinate_transform_verification_status": "verified",
                "coordinate_frame_relation": (
                    "verified_by_exact_common_receptor_atom_correspondence"
                ),
            }
        )
        selection, selection_sha = _parse_pose_selection_manifest(
            path=selection_manifest,
            expected_sha256=pose_selection_manifest_sha256,
            docked_poses_sha256=poses_sha,
        )
        base.update(
            {
                "pose_selection_manifest_sha256": selection_sha,
                "selection_engine": selection.engine,
                "selection_stage": selection.stage,
                "selection_score_source": selection.score_source,
                "selection_score_value": selection.score_value,
                "selection_score_direction": selection.score_direction,
                "selection_tie_break_rule": selection.tie_break_rule,
                "manifest_selected_pose_ordinal": selection.selected_model_ordinal,
                "manifest_selected_model_label": selection.selected_model_label,
                "manifest_selected_model_sha256": selection.selected_model_sha256,
            }
        )

        reference = _load_one_heavy_molecule(reference_path, "reference_ligand")
        prepared_graph = _load_one_heavy_molecule(graph_path, "prepared_ligand_graph")
        _validate_reference_coordinates(reference)
        transformed_reference = _transform_reference(reference, transform)
        reference_identity = _canonical_identity(reference)
        prepared_identity = _canonical_identity(prepared_graph)
        base["reference_identity"] = reference_identity
        base["prepared_identity"] = prepared_identity
        if reference.GetNumAtoms() != prepared_graph.GetNumAtoms():
            raise _NativeRedockInputError("heavy_atom_count_mismatch")
        if reference_identity != prepared_identity:
            raise _NativeRedockInputError("ligand_identity_mismatch")

        prepared_models = _parse_pdbqt_models(prepared_pdbqt_path)
        if len(prepared_models) != 1:
            raise _NativeRedockInputError(
                f"prepared_ligand_model_count_not_one:{len(prepared_models)}"
            )
        prepared_model = prepared_models[0]
        prepared_mapping = _parse_meeko_mapping(prepared_model)
        meeko_identity = _canonical_identity(prepared_mapping.heavy_mol)
        if meeko_identity != reference_identity:
            raise _NativeRedockInputError("meeko_smiles_ligand_identity_mismatch")
        base.update(
            {
                "meeko_smiles": prepared_mapping.smiles,
                "meeko_mapping_sha256": prepared_mapping.mapping_sha256,
                "heavy_atom_count": prepared_mapping.heavy_mol.GetNumAtoms(),
                "prepared_atom_signature_sha256": _atom_signature_sha256(
                    prepared_model
                ),
            }
        )

        pose_models = _parse_pdbqt_models(poses_path)
        base["generated_pose_count"] = len(pose_models)
        prepared_signature = _atom_signature_sha256(prepared_model)
        pose_molecules: list[Chem.Mol] = []
        for pose_model in pose_models:
            pose_mapping = _parse_meeko_mapping(pose_model)
            if (
                pose_mapping.smiles != prepared_mapping.smiles
                or pose_mapping.pairs != prepared_mapping.pairs
                or pose_mapping.mapping_sha256 != prepared_mapping.mapping_sha256
            ):
                raise _NativeRedockInputError(
                    f"docked_pose_meeko_mapping_mismatch:{pose_model.ordinal}"
                )
            if _atom_signature_sha256(pose_model) != prepared_signature:
                raise _NativeRedockInputError(
                    f"docked_pose_atom_signature_mismatch:{pose_model.ordinal}"
                )
            pose_molecules.append(_pose_from_meeko_mapping(pose_model, pose_mapping))

        if selection.selected_model_ordinal > len(pose_models):
            raise _NativeRedockInputError("pose_selection_model_ordinal_out_of_range")
        selected_model = pose_models[selection.selected_model_ordinal - 1]
        if selected_model.label != selection.selected_model_label:
            raise _NativeRedockInputError("pose_selection_model_label_mismatch")
        if selected_model.raw_sha256 != selection.selected_model_sha256:
            raise _NativeRedockInputError("pose_selection_model_sha256_mismatch")

        if (
            selection.engine.strip().lower() == "vina"
            and selection.score_source == "vina_affinity_kcal_mol"
            and selection.score_direction == "lower_is_better"
        ):
            base["selection_score_verification_method"] = (
                VINA_POSE_SCORE_VERIFICATION_METHOD
            )
        pose_scores = _verify_vina_pose_score_order(selection, pose_models)
        base.update(
            {
                "selection_score_verification_status": "verified",
                "top_ranked_pose_ordinal": selection.selected_model_ordinal,
                "top_ranked_model_label": selection.selected_model_label,
                "top_ranked_model_sha256": selection.selected_model_sha256,
            }
        )

        graph_maps = _graph_maps(
            prepared_mapping.heavy_mol,
            transformed_reference,
            max_graph_mappings,
        )
        base["graph_mapping_count"] = len(graph_maps)
        rdkit_maps = [list(mapping) for mapping in graph_maps]
        rmsds: list[PoseRmsd] = []
        for model, pose_mol, selection_score in zip(
            pose_models, pose_molecules, pose_scores
        ):
            value = float(
                rdMolAlign.CalcRMS(
                    pose_mol,
                    transformed_reference,
                    map=rdkit_maps,
                    maxMatches=max_graph_mappings,
                    symmetrizeConjugatedTerminalGroups=False,
                )
            )
            if not math.isfinite(value) or value < 0:
                raise _NativeRedockInputError("rmsd_result_invalid")
            rmsds.append(
                PoseRmsd(
                    pose_ordinal=model.ordinal,
                    model_label=model.label,
                    model_sha256=model.raw_sha256,
                    selection_score=selection_score,
                    rmsd_a=value,
                )
            )
        top = rmsds[selection.selected_model_ordinal - 1]
        best = min(rmsds, key=lambda item: (item.rmsd_a, item.pose_ordinal))
        qualified = top.rmsd_a <= float(threshold_a)
        return finish(
            "qualified" if qualified else "not_qualified",
            qualified=qualified,
            failure_reason=None,
            evaluated_pose_count=len(rmsds),
            top_ranked_rmsd_a=top.rmsd_a,
            best_generated_pose_ordinal=best.pose_ordinal,
            best_generated_rmsd_a=best.rmsd_a,
            pose_rmsds=tuple(rmsds),
        )
    except _NativeRedockEvidencePending as exc:
        return finish("evidence_pending", qualified=None, failure_reason=str(exc))
    except _NativeRedockInputError as exc:
        return finish("invalid", qualified=None, failure_reason=str(exc))
    except Exception as exc:
        return finish(
            "invalid",
            qualified=None,
            failure_reason=f"rmsd_calculation_failed:{type(exc).__name__}",
        )


__all__ = [
    "DEFAULT_MAX_GRAPH_MAPPINGS",
    "MANIFEST_SCHEMA_VERSION",
    "MAX_GRAPH_MAPPINGS",
    "NATIVE_REDOCK_RMSD_FORMULA",
    "NATIVE_REDOCK_RMSD_LITERATURE_DOIS",
    "NATIVE_REDOCK_RMSD_METHOD_ID",
    "NATIVE_REDOCK_RMSD_SCHEMA_VERSION",
    "NATIVE_REDOCK_RMSD_THRESHOLD_A",
    "NativeRedockRmsdResult",
    "PoseRmsd",
    "canonical_json_sha256",
    "evaluate_native_redock_rmsd",
    "sha256_file",
]
