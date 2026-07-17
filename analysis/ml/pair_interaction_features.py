"""Provenance-first, PDBQT-defensible pair-interaction features.

This module deliberately does not infer which pose contributed to a docking or
rescoring result.  Claim-grade callers must supply an exact, content-bound pose
and receptor.  The only fallback policy is an explicit Stage-1/model-1
sensitivity analysis tied to a caller-supplied frozen Atlas run.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Callable, Mapping, Protocol, Sequence

from protein_prep.pdb_records import pdbqt_line_xyz
from src.config.output_paths import run_output_dir


PAIR_INTERACTION_SCHEMA_VERSION = "atlas_pair_interactions_v1"
GEOMETRY_PROVIDER_NAME = "pdbqt_geometry"
GEOMETRY_PROVIDER_VERSION = "pdbqt_geometry_v1"
DEFAULT_PAIR_KEY_COLUMNS = (
    "run_id",
    "pdb_id",
    "variant",
    "ph_label",
    "ligand_base",
)
EXACT_POSE_COLUMNS = (
    "pose_path",
    "pose_sha256",
    "pose_model_index",
    "receptor_path",
    "receptor_sha256",
)
POSE_POLICY_EXACT = "exact-supplied"
POSE_POLICY_STAGE1_MODEL1 = "stage1-model-1"
_MISSING_TEXT = {"", "na", "n/a", "nan", "none", "null"}
_SHA256_RE = re.compile(r"^(?:sha256:)?([0-9a-f]{64})$")
_SAFE_TOKEN_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.+-]*$")
_MODEL_RE = re.compile(r"^MODEL[ \t]+([1-9][0-9]*)[ \t]*$")
_ENDMDL_RE = re.compile(r"^ENDMDL[ \t]*$")

# AutoDock4 atom types accepted by Atlas' strict native-redock parser.
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
_METAL_ELEMENTS = {"CA", "FE", "MG", "MN", "ZN"}
_POLAR_ELEMENTS = {"N", "O", "S"}
_HYDROPHOBIC_CARBON_TYPES = {"A", "C"}
_VDW_RADII_A = {
    "B": 1.92,
    "BR": 1.85,
    "C": 1.70,
    "CA": 2.31,
    "CL": 1.75,
    "F": 1.47,
    "FE": 2.00,
    "I": 1.98,
    "MG": 1.73,
    "MN": 2.00,
    "N": 1.55,
    "O": 1.52,
    "P": 1.80,
    "S": 1.80,
    "SI": 2.10,
    "ZN": 1.39,
}


class PairInteractionError(ValueError):
    """Base error for fail-closed pair-interaction materialization."""


class PairProvenanceError(PairInteractionError):
    """Raised when pose, receptor, or pair identity is not exact."""


class PairFeatureProviderError(PairInteractionError):
    """Raised when a provider cannot defensibly parse its input artifacts."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json_sha256(payload: Any) -> str:
    encoded = json.dumps(
        payload,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _clean(value: Any) -> str:
    return "" if value is None else str(value).strip()


def _is_missing(value: Any) -> bool:
    return _clean(value).lower() in _MISSING_TEXT


def _canonical_sha256(value: Any, label: str) -> str:
    match = _SHA256_RE.fullmatch(_clean(value))
    if match is None:
        raise PairProvenanceError(f"{label} must be a lowercase SHA-256 digest")
    return match.group(1)


def _positive_int(value: Any, label: str) -> int:
    text = _clean(value)
    if not re.fullmatch(r"[1-9][0-9]*", text):
        raise PairProvenanceError(f"{label} must be a positive integer")
    return int(text)


def _safe_token(value: Any, label: str) -> str:
    token = _clean(value)
    if not token or _SAFE_TOKEN_RE.fullmatch(token) is None:
        raise PairProvenanceError(f"{label} is missing or unsafe: {token!r}")
    if Path(token).name != token or token in {".", ".."}:
        raise PairProvenanceError(f"{label} must be one path-safe token")
    return token


def _resolve_file(path_value: Any, *, base_dir: Path, label: str) -> Path:
    text = _clean(path_value)
    if not text:
        raise PairProvenanceError(f"{label} is missing")
    path = Path(text).expanduser()
    if not path.is_absolute():
        path = base_dir / path
    try:
        resolved = path.resolve(strict=True)
    except OSError as exc:
        raise PairProvenanceError(
            f"{label} does not resolve to a file: {path}"
        ) from exc
    if not resolved.is_file():
        raise PairProvenanceError(f"{label} is not a file: {resolved}")
    if resolved.suffix.lower() != ".pdbqt":
        raise PairProvenanceError(f"{label} must be a PDBQT file: {resolved}")
    return resolved


@dataclass(frozen=True)
class PairKey:
    columns: tuple[str, ...]
    values: tuple[str, ...]

    def as_dict(self) -> dict[str, str]:
        return dict(zip(self.columns, self.values))

    @property
    def canonical_json(self) -> str:
        return json.dumps(
            self.as_dict(),
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.canonical_json.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class ResolvedPairPose:
    pair_key: PairKey
    pose_path: Path
    pose_sha256: str
    pose_model_index: int
    receptor_path: Path
    receptor_sha256: str
    pose_policy: str
    claim_status: str
    frozen_run_id: str | None = None


@dataclass(frozen=True)
class ProviderResult:
    features: Mapping[str, int | float | None]
    provenance: Mapping[str, str | int | float | None]
    unavailable_reasons: tuple[str, ...] = ()


class PairFeatureProvider(Protocol):
    name: str
    version: str
    config_sha256: str

    @property
    def config_payload(self) -> Mapping[str, Any]: ...

    def compute(self, pair: ResolvedPairPose) -> ProviderResult: ...


ProviderFactory = Callable[[Mapping[str, Any] | None], PairFeatureProvider]


class ProviderRegistry:
    """Small collision-safe registry for independently versioned feature blocks."""

    def __init__(self) -> None:
        self._factories: dict[str, ProviderFactory] = {}

    def register(self, name: str, factory: ProviderFactory) -> None:
        token = _clean(name)
        if not token or token in self._factories:
            raise PairInteractionError(
                f"provider name is empty or registered: {name!r}"
            )
        self._factories[token] = factory

    def names(self) -> tuple[str, ...]:
        return tuple(sorted(self._factories))

    def create(
        self,
        name: str,
        config: Mapping[str, Any] | None = None,
    ) -> PairFeatureProvider:
        try:
            factory = self._factories[name]
        except KeyError as exc:
            raise PairInteractionError(
                f"unknown pair feature provider {name!r}; available={self.names()}"
            ) from exc
        return factory(config)


@dataclass(frozen=True)
class PdbqtGeometryConfig:
    shell_edges_a: tuple[float, ...] = (2.5, 3.5, 5.0)
    vdw_overlap_scale: float = 0.75
    hydrophobic_carbon_cutoff_a: float = 4.0
    polar_nos_cutoff_a: float = 3.5
    metal_polar_cutoff_a: float = 3.0
    coulomb_cutoff_a: float = 8.0
    max_contact_feature_columns: int = 2048

    def __post_init__(self) -> None:
        edges = self.shell_edges_a
        if not edges or any(not math.isfinite(edge) or edge <= 0 for edge in edges):
            raise PairInteractionError("shell edges must be positive finite values")
        if tuple(sorted(set(edges))) != edges:
            raise PairInteractionError("shell edges must be strictly increasing")
        positive = (
            self.vdw_overlap_scale,
            self.hydrophobic_carbon_cutoff_a,
            self.polar_nos_cutoff_a,
            self.metal_polar_cutoff_a,
            self.coulomb_cutoff_a,
        )
        if any(not math.isfinite(value) or value <= 0 for value in positive):
            raise PairInteractionError("geometry cutoffs and scales must be positive")
        if self.coulomb_cutoff_a < max(edges):
            raise PairInteractionError("Coulomb cutoff must cover every contact shell")
        if self.max_contact_feature_columns < 1:
            raise PairInteractionError("max_contact_feature_columns must be positive")

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any] | None) -> PdbqtGeometryConfig:
        if not value:
            return cls()
        allowed = set(cls.__dataclass_fields__)
        unknown = sorted(set(value) - allowed)
        if unknown:
            raise PairInteractionError(f"unknown geometry config fields: {unknown}")
        kwargs = dict(value)
        if "shell_edges_a" in kwargs:
            kwargs["shell_edges_a"] = tuple(
                float(item) for item in kwargs["shell_edges_a"]
            )
        return cls(**kwargs)


@dataclass(frozen=True)
class _PdbqtAtom:
    atom_type: str
    element: str
    x: float
    y: float
    z: float
    charge: float | None


@dataclass(frozen=True)
class _ParsedModel:
    ordinal: int
    atoms: tuple[_PdbqtAtom, ...]


def _parse_charge(line: str) -> float | None:
    if len(line) < 76:
        return None
    token = line[70:76].strip()
    if not token:
        return None
    try:
        value = float(token)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def _build_model(lines: Sequence[str], ordinal: int) -> _ParsedModel:
    atoms: list[_PdbqtAtom] = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        xyz = pdbqt_line_xyz(line)
        if xyz is None or not all(math.isfinite(value) for value in xyz):
            raise PairFeatureProviderError(
                f"malformed PDBQT coordinate in model {ordinal}"
            )
        fields = line.split()
        atom_type = fields[-1] if fields else ""
        element = _ATOM_TYPE_ELEMENTS.get(atom_type)
        if element is None:
            raise PairFeatureProviderError(
                f"unsupported AutoDock atom type {atom_type!r} in model {ordinal}"
            )
        if element == "H":
            continue
        atoms.append(
            _PdbqtAtom(
                atom_type=atom_type,
                element=element,
                x=float(xyz[0]),
                y=float(xyz[1]),
                z=float(xyz[2]),
                charge=_parse_charge(line),
            )
        )
    if not atoms:
        raise PairFeatureProviderError(f"PDBQT model {ordinal} has no heavy atoms")
    return _ParsedModel(ordinal=ordinal, atoms=tuple(atoms))


def _parse_pdbqt_models(path: Path) -> tuple[_ParsedModel, ...]:
    try:
        lines = path.read_text(encoding="ascii").splitlines()
    except (OSError, UnicodeDecodeError) as exc:
        raise PairFeatureProviderError(f"cannot read ASCII PDBQT: {path}") from exc
    if not lines:
        raise PairFeatureProviderError(f"empty PDBQT file: {path}")
    explicit = any(line.startswith("MODEL") for line in lines)
    if not explicit:
        if any(line.startswith("ENDMDL") for line in lines):
            raise PairFeatureProviderError(f"unexpected ENDMDL in {path}")
        return (_build_model(lines, 1),)

    models: list[_ParsedModel] = []
    current: list[str] | None = None
    for line in lines:
        if line.startswith("MODEL"):
            match = _MODEL_RE.fullmatch(line)
            if match is None or current is not None:
                raise PairFeatureProviderError(f"malformed or nested MODEL in {path}")
            expected = len(models) + 1
            if int(match.group(1)) != expected:
                raise PairFeatureProviderError(
                    f"nonsequential MODEL label in {path}: expected {expected}"
                )
            current = [line]
        elif line.startswith("ENDMDL"):
            if _ENDMDL_RE.fullmatch(line) is None or current is None:
                raise PairFeatureProviderError(
                    f"malformed or unexpected ENDMDL in {path}"
                )
            current.append(line)
            models.append(_build_model(current, len(models) + 1))
            current = None
        elif current is not None:
            current.append(line)
        elif line.strip():
            raise PairFeatureProviderError(f"content outside MODEL blocks in {path}")
    if current is not None or not models:
        raise PairFeatureProviderError(f"unterminated or absent MODEL block in {path}")
    return tuple(models)


def _distance(left: _PdbqtAtom, right: _PdbqtAtom) -> float:
    return math.sqrt(
        (left.x - right.x) ** 2 + (left.y - right.y) ** 2 + (left.z - right.z) ** 2
    )


def _slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def _angstrom_token(value: float) -> str:
    return f"{value:.3f}".rstrip("0").rstrip(".").replace(".", "p")


def _shell_labels(edges: Sequence[float]) -> tuple[str, ...]:
    lower = 0.0
    labels: list[str] = []
    for upper in edges:
        labels.append(f"{_angstrom_token(lower)}_{_angstrom_token(upper)}a")
        lower = upper
    return tuple(labels)


class _SpatialHash:
    def __init__(self, atoms: Sequence[_PdbqtAtom], cell_size: float) -> None:
        self.atoms = tuple(atoms)
        self.cell_size = cell_size
        cells: dict[tuple[int, int, int], list[int]] = {}
        for index, atom in enumerate(self.atoms):
            cells.setdefault(self._cell(atom.x, atom.y, atom.z), []).append(index)
        self.cells = cells

    def _cell(self, x: float, y: float, z: float) -> tuple[int, int, int]:
        return (
            math.floor(x / self.cell_size),
            math.floor(y / self.cell_size),
            math.floor(z / self.cell_size),
        )

    def neighbors(self, atom: _PdbqtAtom) -> Sequence[_PdbqtAtom]:
        cx, cy, cz = self._cell(atom.x, atom.y, atom.z)
        indices: list[int] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    indices.extend(self.cells.get((cx + dx, cy + dy, cz + dz), ()))
        return tuple(self.atoms[index] for index in indices)


class PdbqtGeometryProvider:
    name = GEOMETRY_PROVIDER_NAME
    version = GEOMETRY_PROVIDER_VERSION

    def __init__(self, config: PdbqtGeometryConfig | None = None) -> None:
        self.config = config or PdbqtGeometryConfig()
        self.config_sha256 = canonical_json_sha256(self.config_payload)
        self._receptor_cache: dict[str, tuple[_ParsedModel, _SpatialHash]] = {}

    @property
    def config_payload(self) -> Mapping[str, Any]:
        return {
            "provider": self.name,
            "provider_version": self.version,
            **asdict(self.config),
            "atom_type_vocabulary": sorted(_ATOM_TYPE_ELEMENTS),
            "distance_shell_semantics": "lower_exclusive_upper_inclusive",
            "polar_proxy_semantics": "N/O/S distance only; not a hydrogen bond",
            "coulomb_proxy_semantics": "sum(q_ligand*q_receptor/r_angstrom)",
        }

    def _receptor(self, pair: ResolvedPairPose) -> tuple[_ParsedModel, _SpatialHash]:
        cached = self._receptor_cache.get(pair.receptor_sha256)
        if cached is not None:
            return cached
        models = _parse_pdbqt_models(pair.receptor_path)
        if len(models) != 1:
            raise PairFeatureProviderError(
                "receptor PDBQT must contain exactly one coordinate model"
            )
        model = models[0]
        index = _SpatialHash(model.atoms, self.config.coulomb_cutoff_a)
        self._receptor_cache[pair.receptor_sha256] = (model, index)
        return model, index

    def compute(self, pair: ResolvedPairPose) -> ProviderResult:
        pose_models = _parse_pdbqt_models(pair.pose_path)
        if pair.pose_model_index > len(pose_models):
            raise PairFeatureProviderError(
                f"pose model {pair.pose_model_index} absent from {pair.pose_path}"
            )
        ligand = pose_models[pair.pose_model_index - 1]
        receptor, receptor_index = self._receptor(pair)
        return self._compute_geometry(ligand, receptor, receptor_index)

    def _compute_geometry(
        self,
        ligand: _ParsedModel,
        receptor: _ParsedModel,
        receptor_index: _SpatialHash,
    ) -> ProviderResult:
        config = self.config
        shell_labels = _shell_labels(config.shell_edges_a)
        counts: dict[str, int] = {
            f"pi_geom_total_contact_pairs_{shell}": 0 for shell in shell_labels
        }
        minimum_distance = math.inf
        minimum_metal_polar = math.inf
        vdw_overlap_count = 0
        hydrophobic_count = 0
        polar_count = 0
        metal_count = 0
        coulomb_sum = 0.0
        coulomb_pairs = 0
        max_shell = config.shell_edges_a[-1]
        all_charges = all(atom.charge is not None for atom in ligand.atoms) and all(
            atom.charge is not None for atom in receptor.atoms
        )

        for ligand_atom in ligand.atoms:
            for receptor_atom in receptor_index.neighbors(ligand_atom):
                distance = _distance(ligand_atom, receptor_atom)
                if distance > config.coulomb_cutoff_a:
                    continue
                minimum_distance = min(minimum_distance, distance)
                if distance <= max_shell:
                    shell_index = next(
                        index
                        for index, edge in enumerate(config.shell_edges_a)
                        if distance <= edge
                    )
                    shell = shell_labels[shell_index]
                    element_key = (
                        "pi_geom_element_"
                        f"{_slug(ligand_atom.element)}__{_slug(receptor_atom.element)}_"
                        f"contacts_{shell}"
                    )
                    adtype_key = (
                        "pi_geom_adtype_"
                        f"{_slug(ligand_atom.atom_type)}__{_slug(receptor_atom.atom_type)}_"
                        f"contacts_{shell}"
                    )
                    total_key = f"pi_geom_total_contact_pairs_{shell}"
                    for key in (element_key, adtype_key, total_key):
                        counts[key] = counts.get(key, 0) + 1
                ligand_radius = _VDW_RADII_A.get(ligand_atom.element)
                receptor_radius = _VDW_RADII_A.get(receptor_atom.element)
                if (
                    ligand_radius is not None
                    and receptor_radius is not None
                    and distance
                    < config.vdw_overlap_scale * (ligand_radius + receptor_radius)
                ):
                    vdw_overlap_count += 1
                if (
                    distance <= config.hydrophobic_carbon_cutoff_a
                    and ligand_atom.atom_type in _HYDROPHOBIC_CARBON_TYPES
                    and receptor_atom.atom_type in _HYDROPHOBIC_CARBON_TYPES
                ):
                    hydrophobic_count += 1
                if (
                    distance <= config.polar_nos_cutoff_a
                    and ligand_atom.element in _POLAR_ELEMENTS
                    and receptor_atom.element in _POLAR_ELEMENTS
                ):
                    polar_count += 1
                metal_polar = (
                    ligand_atom.element in _METAL_ELEMENTS
                    and receptor_atom.element in _POLAR_ELEMENTS
                ) or (
                    receptor_atom.element in _METAL_ELEMENTS
                    and ligand_atom.element in _POLAR_ELEMENTS
                )
                if metal_polar:
                    minimum_metal_polar = min(minimum_metal_polar, distance)
                    if distance <= config.metal_polar_cutoff_a:
                        metal_count += 1
                if all_charges:
                    assert ligand_atom.charge is not None
                    assert receptor_atom.charge is not None
                    if distance <= 0:
                        raise PairFeatureProviderError(
                            "coincident atoms make Coulomb proxy undefined"
                        )
                    coulomb_sum += ligand_atom.charge * receptor_atom.charge / distance
                    coulomb_pairs += 1

        if math.isinf(minimum_distance):
            minimum_distance = min(
                _distance(ligand_atom, receptor_atom)
                for ligand_atom in ligand.atoms
                for receptor_atom in receptor.atoms
            )
        if math.isinf(minimum_metal_polar):
            metal_polar_distances = [
                _distance(ligand_atom, receptor_atom)
                for ligand_atom in ligand.atoms
                for receptor_atom in receptor.atoms
                if (
                    ligand_atom.element in _METAL_ELEMENTS
                    and receptor_atom.element in _POLAR_ELEMENTS
                )
                or (
                    receptor_atom.element in _METAL_ELEMENTS
                    and ligand_atom.element in _POLAR_ELEMENTS
                )
            ]
            if metal_polar_distances:
                minimum_metal_polar = min(metal_polar_distances)

        ligand_charge_count = sum(atom.charge is not None for atom in ligand.atoms)
        receptor_charge_count = sum(atom.charge is not None for atom in receptor.atoms)
        hydrophobic_token = _angstrom_token(config.hydrophobic_carbon_cutoff_a)
        polar_token = _angstrom_token(config.polar_nos_cutoff_a)
        metal_token = _angstrom_token(config.metal_polar_cutoff_a)
        coulomb_token = _angstrom_token(config.coulomb_cutoff_a)
        features: dict[str, int | float | None] = {
            **counts,
            "pi_geom_ligand_heavy_atom_count": len(ligand.atoms),
            "pi_geom_receptor_heavy_atom_count": len(receptor.atoms),
            "pi_geom_ligand_partial_charge_coverage_fraction": (
                ligand_charge_count / len(ligand.atoms)
            ),
            "pi_geom_receptor_partial_charge_coverage_fraction": (
                receptor_charge_count / len(receptor.atoms)
            ),
            "pi_geom_min_heavy_atom_distance_a": minimum_distance,
            "pi_geom_vdw_overlap_proxy_count": vdw_overlap_count,
            f"pi_geom_hydrophobic_carbon_contact_proxy_count_le_{hydrophobic_token}a": hydrophobic_count,
            f"pi_geom_polar_nos_proximity_proxy_count_le_{polar_token}a": polar_count,
            f"pi_geom_metal_polar_distance_proxy_count_le_{metal_token}a": metal_count,
            "pi_geom_min_metal_polar_distance_a": (
                None if math.isinf(minimum_metal_polar) else minimum_metal_polar
            ),
            f"pi_geom_partial_charge_coulomb_proxy_sum_qiqj_over_r_le_{coulomb_token}a": (
                coulomb_sum if all_charges else None
            ),
            f"pi_geom_partial_charge_coulomb_proxy_pair_count_le_{coulomb_token}a": (
                coulomb_pairs if all_charges else None
            ),
        }
        dynamic_count = sum(
            key.startswith(("pi_geom_element_", "pi_geom_adtype_")) for key in features
        )
        if dynamic_count > config.max_contact_feature_columns:
            raise PairFeatureProviderError(
                "contact feature vocabulary exceeds configured safety bound"
            )
        unavailable: list[str] = []
        if not all_charges:
            unavailable.append(
                "partial_charge_coulomb_proxy_unavailable:incomplete_heavy_atom_charges"
            )
        if math.isinf(minimum_metal_polar):
            unavailable.append("metal_distance_proxy_unavailable:no_metal_polar_pair")
        return ProviderResult(
            features=features,
            provenance={
                "pi_geom_provider_version": self.version,
                "pi_geom_provider_config_sha256": self.config_sha256,
                "pi_geom_pose_model_index": ligand.ordinal,
                "pi_geom_distance_units": "angstrom",
                "pi_geom_coulomb_proxy_units": "pdbqt_charge_squared_per_angstrom",
            },
            unavailable_reasons=tuple(unavailable),
        )


def _geometry_factory(config: Mapping[str, Any] | None) -> PairFeatureProvider:
    return PdbqtGeometryProvider(PdbqtGeometryConfig.from_mapping(config))


PAIR_FEATURE_PROVIDER_REGISTRY = ProviderRegistry()
PAIR_FEATURE_PROVIDER_REGISTRY.register(GEOMETRY_PROVIDER_NAME, _geometry_factory)


def _read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        if not fields:
            raise PairInteractionError(f"CSV has no header: {path}")
        duplicate_fields = sorted(
            {field for field in fields if fields.count(field) > 1}
        )
        if duplicate_fields:
            raise PairInteractionError(f"CSV has duplicate columns: {duplicate_fields}")
        rows = [dict(row) for row in reader]
        if not rows:
            raise PairInteractionError(f"CSV has no data rows: {path}")
        malformed = [index for index, row in enumerate(rows, start=2) if None in row]
        if malformed:
            raise PairInteractionError(
                f"CSV rows have more values than header columns: {malformed[:10]}"
            )
        return fields, rows


def _canonical_pair_key(
    row: Mapping[str, Any],
    columns: Sequence[str],
    row_number: int,
) -> PairKey:
    values: list[str] = []
    for column in columns:
        value = _clean(row.get(column))
        if column == "pdb_id":
            value = value.upper()
        elif column == "variant":
            value = value.upper() or "LEGACY"
        if not value and column != "ph_label":
            raise PairProvenanceError(
                f"row {row_number}: pair key column {column!r} is missing"
            )
        values.append(value)
    return PairKey(tuple(columns), tuple(values))


def _require_unique_keys(keys: Sequence[PairKey]) -> None:
    seen: dict[tuple[str, ...], int] = {}
    for row_number, key in enumerate(keys, start=2):
        prior = seen.get(key.values)
        if prior is not None:
            raise PairProvenanceError(
                "duplicate canonical pair key; first-row deduplication is forbidden: "
                f"rows {prior} and {row_number}, key={key.canonical_json}"
            )
        seen[key.values] = row_number


def _verified_hash(
    path: Path,
    expected: str,
    cache: dict[Path, str],
    label: str,
) -> str:
    actual = cache.get(path)
    if actual is None:
        actual = sha256_file(path)
        cache[path] = actual
    if actual != expected:
        raise PairProvenanceError(
            f"{label} SHA-256 mismatch for {path}: expected {expected}, observed {actual}"
        )
    return actual


def _exact_pair_pose(
    row: Mapping[str, Any],
    key: PairKey,
    *,
    base_dir: Path,
    hash_cache: dict[Path, str],
) -> ResolvedPairPose:
    pose_path = _resolve_file(
        row.get("pose_path"), base_dir=base_dir, label="pose_path"
    )
    receptor_path = _resolve_file(
        row.get("receptor_path"), base_dir=base_dir, label="receptor_path"
    )
    pose_sha256 = _canonical_sha256(row.get("pose_sha256"), "pose_sha256")
    receptor_sha256 = _canonical_sha256(row.get("receptor_sha256"), "receptor_sha256")
    _verified_hash(pose_path, pose_sha256, hash_cache, "pose")
    _verified_hash(receptor_path, receptor_sha256, hash_cache, "receptor")
    return ResolvedPairPose(
        pair_key=key,
        pose_path=pose_path,
        pose_sha256=pose_sha256,
        pose_model_index=_positive_int(row.get("pose_model_index"), "pose_model_index"),
        receptor_path=receptor_path,
        receptor_sha256=receptor_sha256,
        pose_policy=POSE_POLICY_EXACT,
        claim_status="exact_pose_identity_supplied",
    )


def _context_variant(value: Any) -> str:
    token = _clean(value).upper()
    return (
        "LEGACY" if token in {"", "LEGACY", "NONE"} else _safe_token(token, "variant")
    )


def _exploratory_paths(
    row: Mapping[str, Any],
    *,
    frozen_run_root: Path,
    frozen_run_id: str,
) -> tuple[Path, Path, str, str]:
    pdb_id = _safe_token(_clean(row.get("pdb_id")).upper(), "pdb_id")
    ligand_base = _safe_token(row.get("ligand_base"), "ligand_base")
    variant = _context_variant(row.get("variant"))
    ph_label = _clean(row.get("ph_label"))
    if ph_label:
        ph_label = _safe_token(ph_label, "ph_label")
    docked_root = run_output_dir(frozen_run_root, "docked", frozen_run_id)
    processed_root = run_output_dir(frozen_run_root, "processed_pdbs", frozen_run_id)
    context_parts = [] if variant == "LEGACY" else [variant]
    if ph_label:
        context_parts.append(ph_label)
    pose_path = docked_root.joinpath(
        pdb_id, *context_parts, "stage1", f"{ligand_base}_stage1.pdbqt"
    )
    receptor_dir = processed_root.joinpath(
        pdb_id, *([] if variant == "LEGACY" else [variant]), "receptor"
    )
    receptor_path = (
        receptor_dir / "ph_ensemble" / f"{pdb_id}_{ph_label}.pdbqt"
        if ph_label
        else receptor_dir / f"{pdb_id}.pdbqt"
    )
    return pose_path, receptor_path, variant, ph_label


def _exploratory_pair_pose(
    row: dict[str, str],
    key_columns: Sequence[str],
    *,
    frozen_run_root: Path,
    frozen_run_id: str,
    hash_cache: dict[Path, str],
    row_number: int,
) -> tuple[PairKey, ResolvedPairPose]:
    existing_run = _clean(row.get("run_id"))
    if existing_run and existing_run != frozen_run_id:
        raise PairProvenanceError(
            f"row {row_number}: run_id {existing_run!r} does not match frozen run"
        )
    pose_candidate, receptor_candidate, variant, ph_label = _exploratory_paths(
        row,
        frozen_run_root=frozen_run_root,
        frozen_run_id=frozen_run_id,
    )
    pose_path = _resolve_file(
        pose_candidate, base_dir=frozen_run_root, label="stage1 pose_path"
    )
    receptor_path = _resolve_file(
        receptor_candidate, base_dir=frozen_run_root, label="frozen receptor_path"
    )
    pose_sha256 = hash_cache.get(pose_path)
    if pose_sha256 is None:
        pose_sha256 = sha256_file(pose_path)
        hash_cache[pose_path] = pose_sha256
    receptor_sha256 = hash_cache.get(receptor_path)
    if receptor_sha256 is None:
        receptor_sha256 = sha256_file(receptor_path)
        hash_cache[receptor_path] = receptor_sha256
    generated = {
        "run_id": frozen_run_id,
        "variant": variant,
        "ph_label": ph_label,
        "pose_path": str(pose_path),
        "pose_sha256": pose_sha256,
        "pose_model_index": "1",
        "receptor_path": str(receptor_path),
        "receptor_sha256": receptor_sha256,
    }
    for column, value in generated.items():
        existing = _clean(row.get(column))
        if existing and existing != value:
            raise PairProvenanceError(
                f"row {row_number}: exploratory policy conflicts with {column}={existing!r}"
            )
        row[column] = value
    key = _canonical_pair_key(row, key_columns, row_number)
    return key, ResolvedPairPose(
        pair_key=key,
        pose_path=pose_path,
        pose_sha256=pose_sha256,
        pose_model_index=1,
        receptor_path=receptor_path,
        receptor_sha256=receptor_sha256,
        pose_policy=POSE_POLICY_STAGE1_MODEL1,
        claim_status="sensitivity_stage1_model1",
        frozen_run_id=frozen_run_id,
    )


def _serialize_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, float) and not math.isfinite(value):
        raise PairInteractionError("provider produced a nonfinite feature")
    return value


def _write_csv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def materialize_pair_interaction_features(
    dataset: Path,
    output: Path,
    manifest_output: Path,
    *,
    pose_policy: str = POSE_POLICY_EXACT,
    pair_key_columns: Sequence[str] = DEFAULT_PAIR_KEY_COLUMNS,
    provider_names: Sequence[str] = (GEOMETRY_PROVIDER_NAME,),
    provider_configs: Mapping[str, Mapping[str, Any]] | None = None,
    frozen_run_id: str | None = None,
    frozen_run_root: Path | None = None,
    overwrite: bool = False,
    registry: ProviderRegistry = PAIR_FEATURE_PROVIDER_REGISTRY,
) -> dict[str, Any]:
    """Append pair features while preserving exact row order and cardinality."""

    dataset = Path(dataset).expanduser().resolve(strict=True)
    if not dataset.is_file():
        raise PairInteractionError(f"dataset is not a file: {dataset}")
    output = Path(output).expanduser().resolve()
    manifest_output = Path(manifest_output).expanduser().resolve()
    if output == dataset or manifest_output == dataset or output == manifest_output:
        raise PairInteractionError("dataset, output, and manifest paths must differ")
    if not overwrite and (output.exists() or manifest_output.exists()):
        raise FileExistsError("output exists; pass overwrite=True only after review")
    if pose_policy not in {POSE_POLICY_EXACT, POSE_POLICY_STAGE1_MODEL1}:
        raise PairInteractionError(f"unsupported pose policy: {pose_policy}")
    if not pair_key_columns or len(set(pair_key_columns)) != len(pair_key_columns):
        raise PairInteractionError("pair key columns must be nonempty and unique")
    if not provider_names or len(set(provider_names)) != len(provider_names):
        raise PairInteractionError("provider names must be nonempty and unique")
    unused_configs = sorted(set(provider_configs or {}) - set(provider_names))
    if unused_configs:
        raise PairInteractionError(
            f"configuration supplied for unselected providers: {unused_configs}"
        )
    if pose_policy == POSE_POLICY_STAGE1_MODEL1:
        if not frozen_run_id or frozen_run_root is None:
            raise PairProvenanceError(
                "stage1-model-1 requires frozen_run_id and frozen_run_root"
            )
        frozen_run_id = _safe_token(frozen_run_id, "frozen_run_id")
        frozen_run_root = Path(frozen_run_root).expanduser().resolve(strict=True)
    elif frozen_run_id is not None or frozen_run_root is not None:
        raise PairProvenanceError(
            "frozen run arguments are valid only with stage1-model-1"
        )

    input_fields, input_rows = _read_csv(dataset)
    missing_keys = sorted(set(pair_key_columns) - set(input_fields))
    if pose_policy == POSE_POLICY_STAGE1_MODEL1:
        missing_keys = [
            column
            for column in missing_keys
            if column not in {"run_id", "variant", "ph_label"}
        ]
    if missing_keys:
        raise PairProvenanceError(
            f"dataset is missing pair key columns: {missing_keys}"
        )
    if pose_policy == POSE_POLICY_EXACT:
        missing_exact = sorted(set(EXACT_POSE_COLUMNS) - set(input_fields))
        if missing_exact:
            raise PairProvenanceError(
                f"exact-supplied policy requires columns: {missing_exact}"
            )

    rows = [dict(row) for row in input_rows]
    hash_cache: dict[Path, str] = {}
    keys: list[PairKey] = []
    resolved_pairs: list[ResolvedPairPose] = []
    if pose_policy == POSE_POLICY_EXACT:
        for row_number, row in enumerate(rows, start=2):
            key = _canonical_pair_key(row, pair_key_columns, row_number)
            keys.append(key)
            resolved_pairs.append(
                _exact_pair_pose(
                    row,
                    key,
                    base_dir=dataset.parent,
                    hash_cache=hash_cache,
                )
            )
    else:
        assert frozen_run_root is not None
        assert frozen_run_id is not None
        for row_number, row in enumerate(rows, start=2):
            key, pair = _exploratory_pair_pose(
                row,
                pair_key_columns,
                frozen_run_root=frozen_run_root,
                frozen_run_id=frozen_run_id,
                hash_cache=hash_cache,
                row_number=row_number,
            )
            keys.append(key)
            resolved_pairs.append(pair)
    _require_unique_keys(keys)

    providers = [
        registry.create(name, (provider_configs or {}).get(name))
        for name in provider_names
    ]
    added_rows: list[dict[str, Any]] = []
    all_feature_columns: set[str] = set()
    all_provenance_columns: set[str] = set()
    provider_feature_columns: dict[str, set[str]] = {
        provider.name: set() for provider in providers
    }
    count_columns: set[str] = set()
    for row, pair in zip(rows, resolved_pairs):
        additions: dict[str, Any] = {
            "pair_interaction_pair_key_json": pair.pair_key.canonical_json,
            "pair_interaction_pair_key_sha256": pair.pair_key.sha256,
            "pair_interaction_pose_policy": pair.pose_policy,
            "pair_interaction_claim_status": pair.claim_status,
            "pair_interaction_pose_path": str(pair.pose_path),
            "pair_interaction_pose_sha256": pair.pose_sha256,
            "pair_interaction_pose_model_index": pair.pose_model_index,
            "pair_interaction_receptor_path": str(pair.receptor_path),
            "pair_interaction_receptor_sha256": pair.receptor_sha256,
            "pair_interaction_schema_version": PAIR_INTERACTION_SCHEMA_VERSION,
        }
        unavailable: list[str] = []
        for provider in providers:
            result = provider.compute(pair)
            for source in (result.features, result.provenance):
                collisions = sorted(set(additions) & set(source))
                if collisions:
                    raise PairInteractionError(
                        f"provider output collision for {provider.name}: {collisions}"
                    )
                additions.update(source)
            all_feature_columns.update(result.features)
            provider_feature_columns[provider.name].update(result.features)
            all_provenance_columns.update(result.provenance)
            count_columns.update(
                column
                for column in result.features
                if "_contacts_" in column or "_contact_pairs_" in column
            )
            unavailable.extend(
                f"{provider.name}:{reason}" for reason in result.unavailable_reasons
            )
        additions["pair_interaction_provider_names"] = "|".join(
            provider.name for provider in providers
        )
        additions["pair_interaction_provider_versions"] = "|".join(
            f"{provider.name}:{provider.version}" for provider in providers
        )
        additions["pair_interaction_provider_config_sha256"] = canonical_json_sha256(
            {provider.name: provider.config_sha256 for provider in providers}
        )
        additions["pair_interaction_status"] = "ready"
        additions["pair_interaction_missing_reasons"] = ""
        additions["pair_interaction_unavailable_reasons"] = "|".join(
            sorted(unavailable)
        )
        collisions = sorted(set(row) & set(additions))
        if collisions:
            raise PairInteractionError(
                f"refusing to overwrite existing materialized columns: {collisions}"
            )
        added_rows.append(additions)

    for provider in providers:
        if provider.name != GEOMETRY_PROVIDER_NAME:
            continue
        contact_columns = {
            column
            for column in provider_feature_columns[provider.name]
            if column.startswith(("pi_geom_element_", "pi_geom_adtype_"))
        }
        limit = int(provider.config_payload["max_contact_feature_columns"])
        if len(contact_columns) > limit:
            raise PairFeatureProviderError(
                "table-wide contact feature vocabulary exceeds configured safety "
                f"bound: observed={len(contact_columns)}, limit={limit}"
            )

    materialized_columns = (
        sorted(set().union(*(set(additions) for additions in added_rows)))
        if added_rows
        else []
    )
    output_fields = list(input_fields)
    for column in (*pair_key_columns, *EXACT_POSE_COLUMNS):
        if column not in output_fields and any(column in row for row in rows):
            output_fields.append(column)
    output_fields.extend(materialized_columns)
    output_rows: list[dict[str, Any]] = []
    for row, additions in zip(rows, added_rows):
        merged: dict[str, Any] = {**row, **additions}
        for column in all_feature_columns:
            if column not in merged:
                merged[column] = 0 if column in count_columns else ""
        output_rows.append(
            {
                column: _serialize_value(merged.get(column, ""))
                for column in output_fields
            }
        )

    _write_csv(output, output_fields, output_rows)
    manifest: dict[str, Any] = {
        "status": "complete",
        "schema_version": PAIR_INTERACTION_SCHEMA_VERSION,
        "dataset": str(dataset),
        "dataset_sha256": sha256_file(dataset),
        "output": str(output),
        "output_sha256": sha256_file(output),
        "rows": len(output_rows),
        "input_columns": len(input_fields),
        "output_columns": len(output_fields),
        "pair_key_columns": list(pair_key_columns),
        "pair_key_set_sha256": canonical_json_sha256(
            [key.canonical_json for key in keys]
        ),
        "pair_keys_unique": True,
        "row_order_preserved": True,
        "first_row_deduplication_allowed": False,
        "pose_policy": pose_policy,
        "claim_status": (
            "sensitivity_stage1_model1"
            if pose_policy == POSE_POLICY_STAGE1_MODEL1
            else "exact_pose_identity_supplied"
        ),
        "selected_final_score_pose_inferred": False,
        "frozen_run_id": frozen_run_id,
        "frozen_run_root": str(frozen_run_root) if frozen_run_root else None,
        "providers": [
            {
                "name": provider.name,
                "version": provider.version,
                "config_sha256": provider.config_sha256,
                "config": provider.config_payload,
            }
            for provider in providers
        ],
        "feature_columns": sorted(all_feature_columns),
        "feature_schema_sha256": canonical_json_sha256(sorted(all_feature_columns)),
        "provider_provenance_columns": sorted(all_provenance_columns),
        "materialized_columns": materialized_columns,
        "limitations": [
            "PDBQT geometry does not establish hydrogen-bond donor/acceptor angles.",
            "Hydrophobic, polar, metal, VDW, and Coulomb values are explicitly proxies.",
            "No selected-final-score pose is inferred from stage or filename.",
            "Stage-1/model-1 output is sensitivity-only even when hashes verify.",
        ],
    }
    manifest_output.parent.mkdir(parents=True, exist_ok=True)
    manifest_output.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=True, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


__all__ = [
    "DEFAULT_PAIR_KEY_COLUMNS",
    "EXACT_POSE_COLUMNS",
    "GEOMETRY_PROVIDER_NAME",
    "PAIR_FEATURE_PROVIDER_REGISTRY",
    "PAIR_INTERACTION_SCHEMA_VERSION",
    "POSE_POLICY_EXACT",
    "POSE_POLICY_STAGE1_MODEL1",
    "PairFeatureProvider",
    "PairFeatureProviderError",
    "PairInteractionError",
    "PairProvenanceError",
    "PdbqtGeometryConfig",
    "PdbqtGeometryProvider",
    "ProviderRegistry",
    "ResolvedPairPose",
    "materialize_pair_interaction_features",
]
