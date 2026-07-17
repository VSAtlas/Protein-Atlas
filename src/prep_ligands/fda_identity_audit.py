"""Offline, structure-first audit of the Atlas FDA ligand mapping.

The audit never rewrites its mapping, source SDF, or prepared ligand library.
It records exact and standardized-parent chemistry separately so expected salts
and solvates are not confused with true structure-association errors.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from rdkit import Chem, rdBase
from rdkit.Chem import inchi, rdFMCS, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from prep_ligands.ligand_library_manager import validate_fda_filter_manifest


STANDARDIZATION_VERSION = "atlas-fda-parent-v2-stereo"

VERIFIED_EXACT = "verified_exact"
VERIFIED_PARENT = "verified_parent"
CONNECTIVITY_EXACT = "connectivity_consistent_exact"
CONNECTIVITY_PARENT = "connectivity_consistent_parent"
CONNECTIVITY_DISAGREES = "connectivity_disagrees"
CONNECTIVITY_UNAVAILABLE = "connectivity_unavailable"
PROBABLE_MISMATCH = "probable_mismatch"
AMBIGUOUS = "ambiguous"
UNVERIFIABLE = "unverifiable"

ACTIVE = "active_ingredient"
SALT_OR_SOLVATE = "salt_or_solvate"
ADDITIVE = "additive_or_nonmedication"
COMBINATION = "multi_active_combination"
METAL_COMPLEX = "metal_complex"
AMBIGUOUS_MIXTURE = "ambiguous_mixture"

CONFIRMED_ACTIVE_CATEGORY = "confirmed FDA active ingredient"
CONFIRMED_SALT_PARENT_CATEGORY = "confirmed FDA salt whose parent was docked"
RETAINED_SALT_CATEGORY = "retained salt/counterion"
COLLAPSED_COMBINATION_CATEGORY = "collapsed combination product"
ADDITIVE_CATEGORY = "additive/excipient"
INCORRECT_MAPPING_CATEGORY = "incorrect name–structure mapping"
UNVERIFIABLE_CATEGORY = "unverifiable"
FDA_CROSS_REFERENCE_CATEGORIES = (
    CONFIRMED_ACTIVE_CATEGORY,
    CONFIRMED_SALT_PARENT_CATEGORY,
    RETAINED_SALT_CATEGORY,
    COLLAPSED_COMBINATION_CATEGORY,
    ADDITIVE_CATEGORY,
    INCORRECT_MAPPING_CATEGORY,
    UNVERIFIABLE_CATEGORY,
)

_RDK_ID_RE = re.compile(r"rdk[_-]?(\d+)", re.IGNORECASE)
_VENDOR_CATALOG_ID_RE = re.compile(
    r"^(?:T(?=[A-Za-z0-9_-]*\d)[A-Za-z0-9_-]+|Fr\d+)$",
    re.IGNORECASE,
)
_PDBQT_SMILES_RE = re.compile(
    r"^REMARK\s+(?:SMILES|SMI)\s*(?::|=)?\s*(\S+)", re.IGNORECASE
)
_PDBQT_INCHIKEY_RE = re.compile(r"^REMARK\s+INCHIKEY\s*(?::|=)?\s*(\S+)", re.IGNORECASE)
_NONMED_NAME_RE = re.compile(
    r"\b(?:"
    r"fd&c|d&c|d and c|lake|tartrazine|carminic acid|carmine|"
    r"chocolate brown|brilliant blue|allura red|sunset yellow|erythrosine|"
    r"indigo carmine|fast green|quinoline yellow|ponceau|amaranth|annatto|"
    r"magnesium stearate|microcrystalline cellulose|cellulose|lactose|"
    r"starch|sucrose|sorbitol|xanthan gum|povidone|crospovidone|"
    r"polyethylene glycol|propylene glycol|polysorbate|hypromellose|talc|"
    r"silicon dioxide|titanium dioxide|iron oxide|ferric oxide|zinc oxide"
    r")\b",
    re.IGNORECASE,
)
_COMPLEX_METAL_ATOMIC_NUMBERS = frozenset(
    {
        13,
        21,
        22,
        23,
        24,
        25,
        26,
        27,
        28,
        29,
        30,
        31,
        39,
        40,
        44,
        45,
        46,
        47,
        48,
        50,
        57,
        78,
        79,
        80,
        81,
        82,
        83,
    }
)
_SALT_OR_SOLVATE_NAME_RE = re.compile(
    r"\b(?:hydrochloride|hydrobromide|sodium|potassium|calcium|magnesium|"
    r"besylate|besilate|mesylate|hydrate|hemihydrate|monohydrate|dihydrate|"
    r"solvate|sulfate|sulphate|hemisulfate|hemisulphate|oxalate|malate|"
    r"tartrate|tosylate|tosilate|tosilas|esylate|isethionate|isethiolate|"
    r"edisylate|epolamine|pamoate|embonate|xinafoate|salicylate|gluconate|"
    r"digluconate|hyclate|diaceturate|benzathine|benzathin|pidolate|etiprate|"
    r"hippurate|erbumine|tromethamine|meglumine|olamine|propanediol|propandiol|"
    r"ethanolate|hemiethanolate|acetone|camphorsulfonate|nicotinate|"
    r"dipropionate|methyl sulfate|metilsulfate|metilsolfato|napadisylate|"
    r"dimethyl sulfoxide|sulfosalicylate|diolamine|diolamin|hydrogen oxalate|"
    r"fendizoate)\b",
    re.IGNORECASE,
)
_NAME_SUPPORT_FIELDS = (
    "pubchem_record_title",
    "pubchem_name",
    "generic_name",
    "drugcentral_generic_name",
)
_NAME_AUXILIARY_FIELDS = ("rxnorm_generic_name", "openfda_generic_name")
_TRUSTED_NAME_QUALIFIER_TOKENS = frozenset(
    {
        "hydrochloride",
        "hydrobromide",
        "chloride",
        "bromide",
        "sodium",
        "potassium",
        "calcium",
        "magnesium",
        "besylate",
        "besilate",
        "mesylate",
        "maleate",
        "fumarate",
        "citrate",
        "acetate",
        "hydrate",
        "hemihydrate",
        "monohydrate",
        "dihydrate",
        "trihydrate",
        "anhydrous",
        "solvate",
        "salt",
    }
)
_NAME_COMPONENT_TOKENS = frozenset({"and", "with", "plus"})
_COMBINATION_NAME_RE = re.compile(
    r"(?:\b(?:and|with|plus|combination)\b|[/+])",
    re.IGNORECASE,
)
_KNOWN_ORGANIC_COUNTERION_SMILES = (
    "CC(=O)O",
    "O=C(O)c1ccccc1",
    "O=C(O)C=CC(=O)O",
    "O=C(O)CCC(=O)O",
    "O=C(O)CC(O)(CC(=O)O)C(=O)O",
    "O=C(O)C(O)C(O)C(=O)O",
    "CC(O)C(=O)O",
    "CS(=O)(=O)O",
    "O=S(=O)(O)c1ccccc1",
    "NCCCC[C@H](N)C(=O)O",
    "N=C(N)NCCC[C@H](N)C(=O)O",
    "CN[C@@H](CO)[C@@H](O)[C@H](O)[C@H](O)CO",
    "NC(CO)(CO)CO",
    "C[N+](C)(C)CCO",
    "O=C(O)CC(O)C(=O)O",
    "Cc1ccc(S(=O)(=O)O)cc1",
    "O=S(=O)(O)CCO",
    "O=C(O)C(=O)O",
    "NCCN",
    "CC(C)(N)CO",
    "CCCCCCCCCCCCOS(=O)(=O)O",
)
_SOURCE_PROPERTY_KIND = {
    "id": "id",
    "structid": "id",
    "structureid": "id",
    "drugcentralid": "id",
    "cas": "cas",
    "casrn": "cas",
    "casnumber": "cas",
    "inchikey": "inchikey",
    "name": "name",
    "preferredname": "name",
    "genericname": "name",
    "drugname": "name",
    "inn": "name",
    "synonyms": "name",
}
_ELIGIBILITY_CONFLICT_REASONS = frozenset(
    {
        "pdbqt_stated_inchikey_conflict",
        "source_sidecar_invalid",
        "source_sidecar_sdf_mismatch",
        "source_sidecar_record_name_mismatch",
        "prepared_candidate_structure_disagreement",
        "pdbqt_sdf_heavy_atom_mismatch",
        "duplicate_ligand_id",
        "possible_additive_or_excipient_name",
    }
)


@dataclass(frozen=True)
class FDAIdentityAuditConfig:
    mapping_csv: Path
    output_dir: Path
    source_sdf: Path | None = None
    approval_manifest: Path | None = None
    library_dir: Path | None = None
    classified_nonmed_csv: Path | None = None
    score_csvs: tuple[Path, ...] = ()


@dataclass(frozen=True)
class FDAIdentityAuditResult:
    row_audit_csv: Path
    evidence_jsonl: Path
    summary_json: Path
    review_csv: Path
    quarantine_manifest_csv: Path
    eligible_manifest_csv: Path
    source_catalog_csv: Path
    approval_cross_reference_csv: Path
    total_rows: int
    blocker_count: int
    review_count: int
    repaired_mapping_csv: Path | None = None
    score_reuse_manifest_csv: Path | None = None
    redock_delta_csv: Path | None = None
    repair_quarantine_csv: Path | None = None

    @property
    def has_blockers(self) -> bool:
        return self.blocker_count > 0


@dataclass(frozen=True)
class _ResolvedConfig:
    mapping_csv: Path
    output_dir: Path
    source_sdf: Path | None
    approval_manifest: Path | None
    library_dir: Path | None
    classified_nonmed_csv: Path | None
    score_csvs: tuple[Path, ...]


@dataclass(frozen=True)
class _AuditOutputPaths:
    row_audit_csv: Path
    evidence_jsonl: Path
    summary_json: Path
    review_csv: Path
    quarantine_csv: Path
    eligible_csv: Path
    source_catalog_csv: Path
    cross_reference_csv: Path


@dataclass
class FragmentEvidence:
    smiles: str
    normalized_smiles: str
    formula: str
    heavy_atoms: int
    formal_charge: int
    organic: bool
    drug_sized: bool
    contains_metal: bool


@dataclass
class StructureEvidence:
    source: str
    raw_smiles: str = ""
    canonical_smiles: str = ""
    exact_inchikey: str = ""
    parent_smiles: str = ""
    parent_inchikey: str = ""
    parent_connectivity_key: str = ""
    stereo_specified: bool = False
    formula: str = ""
    heavy_atoms: int = 0
    formal_charge: int = 0
    fragment_count: int = 0
    fragments: list[FragmentEvidence] = field(default_factory=list)
    error: str = ""


@dataclass
class PDBQTConnectivityEvidence:
    source: str
    method: str
    verdict: str
    heavy_atoms: int = 0
    bonds: int = 0
    error: str = ""


@dataclass
class _NameAssessment:
    support_codes: list[str]
    reason_codes: list[str]
    corroborated: bool
    blocker: bool


@dataclass
class _SourceSDFRecord:
    index: int
    name: str
    names: tuple[str, ...]
    structure: StructureEvidence
    identifiers: dict[str, str]


@dataclass
class _SourceSDFIndex:
    path: Path | None
    by_index: dict[int, _SourceSDFRecord]
    by_identifier: dict[tuple[str, str], list[_SourceSDFRecord]]
    approval_provenance_verified: bool
    approval_provenance_reason: str
    approval_manifest: Path | None


@dataclass
class _AuditedRow:
    row_number: int
    rdk_id: str
    display_name: str
    mapping_path: str
    match_method: str
    source_match_method: str
    identity_verdict: str
    substance_class: str
    name_corroborated: bool
    approval_verified: bool
    fda_drug_eligible: bool
    reason_codes: list[str]
    support_codes: list[str]
    mapping_structure: StructureEvidence
    source_structure: StructureEvidence | None
    source_record: _SourceSDFRecord | None
    display_name_source_record: _SourceSDFRecord | None
    actual_structures: list[StructureEvidence]
    actual_connectivity: list[PDBQTConnectivityEvidence]
    actual_paths: list[str]
    prepared_form_state: str
    curated_nonmed: bool
    mapping_row: dict[str, str]
    cross_reference_category: str = UNVERIFIABLE_CATEGORY
    cross_reference_basis: str = "insufficient_verified_evidence"

    @property
    def blocker(self) -> bool:
        return self.identity_verdict == PROBABLE_MISMATCH

    @property
    def review(self) -> bool:
        return bool(
            self.blocker
            or self.identity_verdict
            in {AMBIGUOUS, UNVERIFIABLE, CONNECTIVITY_EXACT, CONNECTIVITY_PARENT}
            or self.substance_class
            in {
                SALT_OR_SOLVATE,
                ADDITIVE,
                COMBINATION,
                METAL_COMPLEX,
                AMBIGUOUS_MIXTURE,
            }
            or not self.name_corroborated
            or not self.approval_verified
            or any(code.startswith("possible_additive") for code in self.reason_codes)
            or bool(set(self.reason_codes) & _ELIGIBILITY_CONFLICT_REASONS)
        )


def _clean(value: object) -> str:
    text = str(value or "").strip()
    return "" if text.lower() == "nan" else text


def _norm_name(value: object) -> str:
    return re.sub(
        r"\s+", " ", re.sub(r"[^a-z0-9]+", " ", _clean(value).casefold())
    ).strip()


def _name_agrees(left: str, right: str) -> bool:
    left_key = _norm_name(left)
    right_key = _norm_name(right)
    padded_left = f" {left_key} "
    padded_right = f" {right_key} "
    return bool(
        left_key
        and right_key
        and (
            left_key == right_key
            or padded_left in padded_right
            or padded_right in padded_left
        )
    )


def _trusted_name_agrees(left: str, right: str) -> bool:
    left_key = _norm_name(left)
    right_key = _norm_name(right)
    if not left_key or not right_key:
        return False
    if _has_component_conjunction(left_key) or _has_component_conjunction(right_key):
        return False
    left_core = _trusted_name_core(left_key)
    right_core = _trusted_name_core(right_key)
    if not left_core or not right_core:
        return False
    if left_key == right_key:
        return True
    return left_core == right_core


def _has_component_conjunction(normalized_name: str) -> bool:
    return bool(set(normalized_name.split()) & _NAME_COMPONENT_TOKENS)


def _trusted_name_core(normalized_name: str) -> str:
    tokens = normalized_name.split()
    while tokens and tokens[-1] in _TRUSTED_NAME_QUALIFIER_TOKENS:
        tokens.pop()
    return " ".join(tokens)


def _rdk_id(row: dict[str, str], row_number: int) -> str:
    for value in (row.get("path"), row.get("rdk_id"), row.get("id")):
        match = _RDK_ID_RE.search(_clean(value))
        if match:
            return f"rdk_{int(match.group(1)):07d}"
    if _clean(row.get("scheme")).casefold() == "rdk" and _clean(row.get("file_num")):
        try:
            return f"rdk_{int(float(_clean(row.get('file_num')))):07d}"
        except ValueError:
            pass
    return f"mapping_row_{row_number:07d}"


def _safe_mol(mol: Chem.Mol | None) -> Chem.Mol | None:
    if mol is None:
        return None
    copied = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(copied)
    except Exception:
        try:
            copied.UpdatePropertyCache(strict=False)
        except Exception:
            return None
    return copied


def _parent_mol(mol: Chem.Mol) -> Chem.Mol:
    parent = Chem.Mol(mol)
    for transform in (
        rdMolStandardize.Cleanup,
        rdMolStandardize.FragmentParent,
        rdMolStandardize.ChargeParent,
    ):
        try:
            parent = transform(parent)
        except Exception:
            continue
    return parent


def _normalized_fragment_smiles(mol: Chem.Mol) -> str:
    normalized = Chem.Mol(mol)
    for transform in (rdMolStandardize.Cleanup, rdMolStandardize.ChargeParent):
        try:
            normalized = transform(normalized)
        except Exception:
            continue
    return Chem.MolToSmiles(normalized, canonical=True, isomericSmiles=False)


@lru_cache(maxsize=1)
def _known_counterion_smiles() -> frozenset[str]:
    values: set[str] = set()
    for smiles in _KNOWN_ORGANIC_COUNTERION_SMILES:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            values.add(_normalized_fragment_smiles(mol))
    return frozenset(values)


def _fragment_evidence(mol: Chem.Mol) -> FragmentEvidence:
    heavy_atoms = int(mol.GetNumHeavyAtoms())
    atomic_numbers = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    organic = 6 in atomic_numbers
    contains_metal = bool(atomic_numbers & _COMPLEX_METAL_ATOMIC_NUMBERS)
    return FragmentEvidence(
        smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
        normalized_smiles=_normalized_fragment_smiles(mol),
        formula=rdMolDescriptors.CalcMolFormula(mol),
        heavy_atoms=heavy_atoms,
        formal_charge=sum(atom.GetFormalCharge() for atom in mol.GetAtoms()),
        organic=organic,
        drug_sized=organic and heavy_atoms >= 14,
        contains_metal=contains_metal,
    )


def _structure_from_mol(mol: Chem.Mol | None, source: str) -> StructureEvidence:
    sanitized = _safe_mol(mol)
    if sanitized is None:
        return StructureEvidence(source=source, error="rdkit_parse_failed")
    try:
        fragments = [
            _fragment_evidence(fragment)
            for fragment in Chem.GetMolFrags(sanitized, asMols=True, sanitizeFrags=True)
        ]
        parent = _parent_mol(sanitized)
        parent_inchikey = inchi.MolToInchiKey(parent)
        return StructureEvidence(
            source=source,
            raw_smiles=Chem.MolToSmiles(
                sanitized, canonical=False, isomericSmiles=True
            ),
            canonical_smiles=Chem.MolToSmiles(
                sanitized, canonical=True, isomericSmiles=True
            ),
            exact_inchikey=inchi.MolToInchiKey(sanitized),
            parent_smiles=Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True),
            parent_inchikey=parent_inchikey,
            parent_connectivity_key=parent_inchikey.split("-", 1)[0]
            if parent_inchikey
            else "",
            stereo_specified=_has_specified_stereo(parent),
            formula=rdMolDescriptors.CalcMolFormula(sanitized),
            heavy_atoms=int(sanitized.GetNumHeavyAtoms()),
            formal_charge=sum(atom.GetFormalCharge() for atom in sanitized.GetAtoms()),
            fragment_count=len(fragments),
            fragments=fragments,
        )
    except Exception as exc:
        return StructureEvidence(
            source=source, error=f"structure_evidence_failed:{exc}"
        )


def _has_specified_stereo(mol: Chem.Mol) -> bool:
    atom_stereo = any(
        atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        for atom in mol.GetAtoms()
    )
    bond_stereo = any(
        bond.GetStereo() != Chem.BondStereo.STEREONONE for bond in mol.GetBonds()
    )
    return atom_stereo or bond_stereo


def _structure_from_smiles(smiles: str, source: str) -> StructureEvidence:
    text = _clean(smiles)
    if not text:
        return StructureEvidence(source=source, error="missing_smiles")
    try:
        mol = Chem.MolFromSmiles(text)
    except Exception:
        mol = None
    evidence = _structure_from_mol(mol, source)
    evidence.raw_smiles = text
    return evidence


def _mapping_structure(row: dict[str, str]) -> StructureEvidence:
    for field_name in ("smiles", "remark_smiles", "smiles_neutral"):
        value = _clean(row.get(field_name))
        if value:
            return _structure_from_smiles(value, f"mapping:{field_name}")
    return StructureEvidence(source="mapping", error="mapping_smiles_missing")


def _source_value(kind: str, value: object) -> str:
    text = _clean(value)
    if kind in {"inchikey", "parent_inchikey"}:
        return text.upper()
    if kind == "name":
        return _norm_name(text)
    return text.casefold()


def _source_identifiers(
    mol: Chem.Mol | None,
    name: str,
    structure: StructureEvidence,
) -> dict[str, str]:
    identifiers: dict[str, str] = {}
    if name:
        identifiers["name"] = _source_value("name", name)
    if structure.exact_inchikey:
        identifiers["inchikey"] = structure.exact_inchikey.upper()
    if structure.parent_inchikey:
        identifiers["parent_inchikey"] = structure.parent_inchikey.upper()
    if mol is None:
        return identifiers
    for property_name in mol.GetPropNames():
        normalized = re.sub(r"[^a-z0-9]+", "", property_name.casefold())
        kind = _SOURCE_PROPERTY_KIND.get(normalized)
        if kind and kind not in identifiers:
            value = _source_value(kind, mol.GetProp(property_name))
            if value:
                identifiers[kind] = value
    return identifiers


def _source_names(mol: Chem.Mol | None, primary_name: str) -> tuple[str, ...]:
    names = {_clean(primary_name)} if _clean(primary_name) else set()
    if mol is not None:
        for property_name in mol.GetPropNames():
            normalized = re.sub(r"[^a-z0-9]+", "", property_name.casefold())
            if _SOURCE_PROPERTY_KIND.get(normalized) == "name":
                names.update(_split_names(mol.GetProp(property_name)))
    return tuple(sorted(name for name in names if name))


def _read_source_sdf(
    path: Path | None,
    approval_manifest: Path | None,
) -> _SourceSDFIndex:
    if path is None:
        return _SourceSDFIndex(
            path=None,
            by_index={},
            by_identifier={},
            approval_provenance_verified=False,
            approval_provenance_reason="approval_source_not_supplied",
            approval_manifest=None,
        )
    if not path.is_file():
        raise ValueError(f"source SDF not found: {path}")
    provenance_verified, provenance_reason, resolved_manifest = (
        validate_fda_filter_manifest(path, approval_manifest)
    )
    records: dict[int, _SourceSDFRecord] = {}
    by_identifier: dict[tuple[str, str], list[_SourceSDFRecord]] = defaultdict(list)
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
    for index, mol in enumerate(supplier, start=1):
        name = ""
        if mol is not None and mol.HasProp("_Name"):
            name = _clean(mol.GetProp("_Name"))
        structure = _structure_from_mol(mol, f"source_sdf:{index}")
        names = _source_names(mol, name)
        record = _SourceSDFRecord(
            index=index,
            name=name,
            names=names,
            structure=structure,
            identifiers=_source_identifiers(mol, name, structure),
        )
        records[index] = record
        for kind, value in record.identifiers.items():
            by_identifier[(kind, value)].append(record)
        for source_name in record.names:
            key = ("name", _source_value("name", source_name))
            if record not in by_identifier[key]:
                by_identifier[key].append(record)
    return _SourceSDFIndex(
        path=path.resolve(),
        by_index=records,
        by_identifier=dict(by_identifier),
        approval_provenance_verified=provenance_verified,
        approval_provenance_reason=provenance_reason,
        approval_manifest=resolved_manifest,
    )


def _mapping_source_identifiers(
    row: dict[str, str],
    mapping: StructureEvidence,
) -> list[tuple[str, str]]:
    candidates = (
        ("inchikey", mapping.exact_inchikey or row.get("inchikey")),
        ("parent_inchikey", mapping.parent_inchikey),
        ("id", row.get("drugcentral_id")),
        ("id", row.get("id")),
        ("cas", row.get("cas")),
    )
    return list(
        dict.fromkeys(
            (kind, normalized)
            for kind, value in candidates
            if (normalized := _source_value(kind, value))
        )
    )


def _source_record_for_row(
    row: dict[str, str],
    mapping: StructureEvidence,
    source_index: _SourceSDFIndex,
) -> tuple[_SourceSDFRecord | None, str]:
    matched, conflict = _stable_source_matches(row, mapping, source_index)
    if conflict:
        return None, conflict
    record_indices = {record.index for _, record in matched}
    if len(record_indices) > 1:
        return None, "conflicting_stable_identifiers"
    if matched:
        matched_kinds = {kind for kind, _ in matched}
        if "inchikey" in matched_kinds:
            matched_kinds.discard("parent_inchikey")
        kinds = "_".join(sorted(matched_kinds))
        unmatched_suffix = (
            "_unmatched_id" if _has_unmatched_claimed_id(row, source_index) else ""
        )
        return matched[0][1], f"stable_{kinds}{unmatched_suffix}"
    return _fallback_source_record(row, source_index)


def _has_unmatched_claimed_id(
    row: dict[str, str],
    source_index: _SourceSDFIndex,
) -> bool:
    claimed_ids = {
        value
        for field_name in ("drugcentral_id", "id")
        if (value := _source_value("id", row.get(field_name)))
        and not (
            field_name == "id"
            and _VENDOR_CATALOG_ID_RE.fullmatch(_clean(row.get(field_name)))
        )
    }
    return bool(
        claimed_ids
        and any(
            ("id", value) not in source_index.by_identifier for value in claimed_ids
        )
    )


def _stable_source_matches(
    row: dict[str, str],
    mapping: StructureEvidence,
    source_index: _SourceSDFIndex,
) -> tuple[list[tuple[str, _SourceSDFRecord]], str]:
    matched: list[tuple[str, _SourceSDFRecord]] = []
    for kind, value in _mapping_source_identifiers(row, mapping):
        matches = source_index.by_identifier.get((kind, value), [])
        if len(matches) > 1:
            return [], f"ambiguous_{kind}"
        if matches:
            matched.append((kind, matches[0]))
    return matched, ""


def _fallback_source_record(
    row: dict[str, str],
    source_index: _SourceSDFIndex,
) -> tuple[_SourceSDFRecord | None, str]:
    name = _source_value("name", row.get("sdf_title"))
    name_matches = source_index.by_identifier.get(("name", name), []) if name else []
    if len(name_matches) == 1:
        return name_matches[0], "name_only"
    if len(name_matches) > 1:
        return None, "ambiguous_name"
    raw_index = _clean(row.get("sdf_index"))
    if not raw_index:
        return None, "missing"
    try:
        index = int(float(raw_index))
    except ValueError:
        return None, "invalid_index"
    record = source_index.by_index.get(index)
    return (record, "ordinal") if record else (None, "index_not_found")


def _library_index(library_dir: Path | None) -> dict[str, list[Path]]:
    if library_dir is None:
        return {}
    if not library_dir.is_dir():
        raise ValueError(f"prepared library directory not found: {library_dir}")
    index: dict[str, list[Path]] = defaultdict(list)
    for path in sorted(library_dir.rglob("*.pdbqt")):
        match = _RDK_ID_RE.search(path.stem)
        key = f"rdk_{int(match.group(1)):07d}" if match else path.stem
        index[key].append(path)
    return {
        key: _preferred_pdbqt_candidates(key, paths) for key, paths in index.items()
    }


def _preferred_pdbqt_candidates(key: str, paths: Sequence[Path]) -> list[Path]:
    canonical = [path for path in paths if path.stem.casefold() == key.casefold()]
    return canonical or list(paths)


def _read_pdbqt_lines(path: Path) -> list[str] | None:
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            return list(handle)
    except OSError:
        return None


def _first_remark(lines: Sequence[str], pattern: re.Pattern[str]) -> str:
    for line in lines:
        match = pattern.match(line)
        if match:
            return match.group(1).strip()
    return ""


def _pdbqt_structure(path: Path) -> tuple[StructureEvidence | None, str, str]:
    lines = _read_pdbqt_lines(path)
    if lines is None:
        return None, "pdbqt_unreadable", ""
    smiles = _first_remark(lines, _PDBQT_SMILES_RE)
    if not smiles:
        return None, "pdbqt_structure_remark_missing", ""
    evidence = _structure_from_smiles(smiles, f"pdbqt:{path.name}")
    stated_inchikey = _first_remark(lines, _PDBQT_INCHIKEY_RE).upper()
    if _pdbqt_key_conflicts(stated_inchikey, evidence):
        return evidence, "pdbqt_stated_inchikey_conflict", stated_inchikey
    return evidence, "", stated_inchikey


def _pdbqt_connectivity_evidence(
    path: Path,
    mapping: StructureEvidence,
) -> PDBQTConnectivityEvidence:
    lines = _read_pdbqt_lines(path)
    if lines is None:
        return _unavailable_connectivity(path, "pdbqt_unreadable")
    mol = _pdbqt_coordinate_mol(lines)
    method = "rdkit_pdb_proximity_graph"
    if mol is None:
        mol = _obabel_pdbqt_mol(path)
        method = "openbabel_pdbqt_graph"
    if mol is None or mol.GetNumAtoms() == 0:
        return _unavailable_connectivity(path, "pdbqt_connectivity_parse_failed")
    verdict = _connectivity_verdict(mol, mapping)
    return PDBQTConnectivityEvidence(
        source=str(path),
        method=method,
        verdict=verdict,
        heavy_atoms=mol.GetNumAtoms(),
        bonds=mol.GetNumBonds(),
    )


def _unavailable_connectivity(
    path: Path,
    error: str,
) -> PDBQTConnectivityEvidence:
    return PDBQTConnectivityEvidence(
        source=str(path),
        method="",
        verdict=CONNECTIVITY_UNAVAILABLE,
        error=error,
    )


def _pdbqt_coordinate_mol(lines: Sequence[str]) -> Chem.Mol | None:
    atom_lines = [
        f"{line[:66]}\n" for line in lines if line.startswith(("ATOM", "HETATM"))
    ]
    if not atom_lines:
        return None
    try:
        mol = Chem.MolFromPDBBlock(
            "".join(atom_lines) + "END\n",
            sanitize=False,
            removeHs=False,
            proximityBonding=True,
        )
    except Exception:
        mol = None
    return _heavy_only_mol(mol)


def _heavy_only_mol(mol: Chem.Mol | None) -> Chem.Mol | None:
    if mol is None:
        return None
    editable = Chem.RWMol(mol)
    hydrogen_indices = [
        atom.GetIdx() for atom in editable.GetAtoms() if atom.GetAtomicNum() == 1
    ]
    for atom_index in reversed(hydrogen_indices):
        editable.RemoveAtom(atom_index)
    return editable.GetMol()


@lru_cache(maxsize=1)
def _obabel_executable() -> str:
    return shutil.which("obabel") or ""


def _obabel_pdbqt_mol(path: Path) -> Chem.Mol | None:
    executable = _obabel_executable()
    if not executable or not path.is_file():
        return None
    try:
        if path.stat().st_size == 0:
            return None
        completed = subprocess.run(
            [executable, "-ipdbqt", str(path), "-osmi"],
            check=False,
            capture_output=True,
            text=True,
            timeout=15,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    smiles = completed.stdout.split(maxsplit=1)[0] if completed.stdout.strip() else ""
    return _heavy_only_mol(Chem.MolFromSmiles(smiles)) if smiles else None


def _connectivity_verdict(
    actual: Chem.Mol,
    mapping: StructureEvidence,
) -> str:
    exact = Chem.MolFromSmiles(mapping.canonical_smiles)
    if _same_element_connectivity(actual, exact):
        return CONNECTIVITY_EXACT
    parent = Chem.MolFromSmiles(mapping.parent_smiles)
    if _same_element_connectivity(actual, parent):
        return CONNECTIVITY_PARENT
    return CONNECTIVITY_DISAGREES


def _same_element_connectivity(
    actual: Chem.Mol | None,
    expected: Chem.Mol | None,
) -> bool:
    if actual is None or expected is None:
        return False
    expected = _heavy_only_mol(expected)
    if expected is None:
        return False
    if actual.GetNumAtoms() != expected.GetNumAtoms():
        return False
    if actual.GetNumBonds() != expected.GetNumBonds():
        return False
    try:
        match = rdFMCS.FindMCS(
            [actual, expected],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            ringMatchesRingOnly=False,
            completeRingsOnly=False,
            timeout=2,
        )
    except Exception:
        return False
    return bool(
        not match.canceled
        and match.numAtoms == actual.GetNumAtoms()
        and match.numBonds == actual.GetNumBonds()
    )


def _pdbqt_key_conflicts(stated_inchikey: str, evidence: StructureEvidence) -> bool:
    return bool(
        stated_inchikey
        and evidence.exact_inchikey
        and stated_inchikey != evidence.exact_inchikey.upper()
    )


def _stated_key_matches_source(
    stated_inchikey: str,
    source_structure: StructureEvidence | None,
) -> bool:
    if not stated_inchikey or source_structure is None:
        return False
    return stated_inchikey in {
        source_structure.exact_inchikey.upper(),
        source_structure.parent_inchikey.upper(),
    }


def _read_sidecar_payload(sidecar: Path) -> tuple[dict[str, object] | None, str]:
    if not sidecar.is_file():
        return None, ""
    try:
        payload = json.loads(sidecar.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None, "source_sidecar_invalid"
    return (
        (payload, "") if isinstance(payload, dict) else (None, "source_sidecar_invalid")
    )


def _sidecar_record(
    payload: dict[str, object], source_index: _SourceSDFIndex
) -> tuple[_SourceSDFRecord | None, str]:
    if payload.get("chemistry_authoritative") is not True:
        return None, "source_sidecar_not_authoritative"
    payload_source = Path(_clean(payload.get("source_sdf"))).expanduser()
    if source_index.path is None:
        return None, "source_sidecar_record_unavailable"
    if payload_source.name and payload_source.name != source_index.path.name:
        return None, "source_sidecar_sdf_mismatch"
    try:
        record_index = int(str(payload.get("source_record_index")))
    except (TypeError, ValueError):
        return None, "source_sidecar_record_index_invalid"
    record = source_index.by_index.get(record_index)
    if record is None:
        return None, "source_sidecar_record_unavailable"
    payload_name = _clean(payload.get("source_record_name"))
    if (
        payload_name
        and record.name
        and _norm_name(payload_name) != _norm_name(record.name)
    ):
        return None, "source_sidecar_record_name_mismatch"
    return record, ""


def _sidecar_source_structure(
    pdbqt_path: Path,
    source_index: _SourceSDFIndex,
) -> tuple[StructureEvidence | None, str]:
    payload, reason = _read_sidecar_payload(
        pdbqt_path.with_suffix(".ligprep_source.json")
    )
    if payload is None:
        return None, reason
    record, reason = _sidecar_record(payload, source_index)
    if record is None:
        return None, reason
    structure = copy.deepcopy(record.structure)
    structure.source = f"source_sidecar:{pdbqt_path.name}:{record.index}"
    return structure, ""


def _classification_name(
    row: dict[str, str],
    display_name: str,
) -> str:
    values = [
        display_name,
        _clean(row.get("pubchem_record_title")),
        _clean(row.get("generic_name")),
    ]
    deduplicated: dict[str, str] = {}
    for value in values:
        if value:
            deduplicated.setdefault(value.casefold(), value)
    return " | ".join(deduplicated.values())


def _classify_substance(
    structure: StructureEvidence,
    display_name: str,
    curated_nonmed: bool,
) -> str:
    if curated_nonmed:
        return ADDITIVE
    fragments = structure.fragments
    if any(fragment.contains_metal for fragment in fragments):
        return METAL_COMPLEX
    if len(fragments) <= 1:
        return _single_fragment_substance_class(display_name)
    active_fragments = [
        fragment for fragment in fragments if not _is_counterion(fragment)
    ]
    return _multifragment_substance_class(
        fragments,
        active_fragments,
        display_name,
    )


def _single_fragment_substance_class(display_name: str) -> str:
    if _SALT_OR_SOLVATE_NAME_RE.search(display_name):
        return SALT_OR_SOLVATE
    return ACTIVE


def _multifragment_substance_class(
    fragments: Sequence[FragmentEvidence],
    active_fragments: Sequence[FragmentEvidence],
    display_name: str,
) -> str:
    if len(active_fragments) >= 2:
        return _multiple_active_fragment_class(
            fragments,
            active_fragments,
            display_name,
        )
    if len(active_fragments) == 1:
        return SALT_OR_SOLVATE
    return AMBIGUOUS_MIXTURE


def _multiple_active_fragment_class(
    fragments: Sequence[FragmentEvidence],
    active_fragments: Sequence[FragmentEvidence],
    display_name: str,
) -> str:
    unique_active_structures = {
        fragment.normalized_smiles for fragment in active_fragments
    }
    has_counterion_or_solvate = len(fragments) > len(active_fragments)
    if len(unique_active_structures) == 1 and has_counterion_or_solvate:
        return SALT_OR_SOLVATE
    salt_form_name = bool(_SALT_OR_SOLVATE_NAME_RE.search(display_name))
    explicit_combination_name = bool(_COMBINATION_NAME_RE.search(display_name))
    if salt_form_name and not explicit_combination_name:
        return SALT_OR_SOLVATE
    return COMBINATION


def _is_counterion(fragment: FragmentEvidence) -> bool:
    if fragment.contains_metal:
        return False
    if not fragment.organic:
        return True
    return fragment.normalized_smiles in _known_counterion_smiles()


def _compare_structures(
    mapping: StructureEvidence,
    candidates: Iterable[StructureEvidence],
) -> str:
    usable = [candidate for candidate in candidates if candidate.exact_inchikey]
    if not mapping.exact_inchikey or not usable:
        return ""
    if _has_exact_match(mapping, usable):
        return VERIFIED_EXACT
    if _has_parent_match(mapping, usable):
        return VERIFIED_PARENT
    return _connectivity_match_verdict(mapping, usable)


def _has_exact_match(
    mapping: StructureEvidence, candidates: Sequence[StructureEvidence]
) -> bool:
    return any(
        candidate.exact_inchikey == mapping.exact_inchikey for candidate in candidates
    )


def _has_parent_match(
    mapping: StructureEvidence, candidates: Sequence[StructureEvidence]
) -> bool:
    return bool(
        mapping.parent_inchikey
        and any(
            candidate.parent_inchikey == mapping.parent_inchikey
            for candidate in candidates
        )
    )


def _connectivity_match_verdict(
    mapping: StructureEvidence, candidates: Sequence[StructureEvidence]
) -> str:
    matches = [
        candidate
        for candidate in candidates
        if mapping.parent_connectivity_key
        and candidate.parent_connectivity_key == mapping.parent_connectivity_key
    ]
    if not matches:
        return PROBABLE_MISMATCH
    both_specified = mapping.stereo_specified and all(
        candidate.stereo_specified for candidate in matches
    )
    return PROBABLE_MISMATCH if both_specified else AMBIGUOUS


def _int_or_none(value: object) -> int | None:
    text = _clean(value)
    if not text:
        return None
    try:
        return int(float(text))
    except ValueError:
        return None


def _split_names(value: object) -> list[str]:
    return [
        part.strip()
        for part in re.split(r"[;|\r\n]+", _clean(value))
        if part.strip()
    ]


def _field_values(row: dict[str, str], fields: Sequence[str]) -> list[tuple[str, str]]:
    return [
        (field_name, value)
        for field_name in fields
        for value in _split_names(row.get(field_name))
    ]


def _support_codes(display_name: str, values: Sequence[tuple[str, str]]) -> list[str]:
    return sorted(
        {
            f"{field_name}_supports_display"
            for field_name, value in values
            if _name_agrees(display_name, value)
        }
    )


def _trusted_support_codes(
    display_name: str,
    values: Sequence[tuple[str, str]],
) -> list[str]:
    return sorted(
        {
            f"{field_name}_supports_display"
            for field_name, value in values
            if _trusted_name_agrees(display_name, value)
        }
    )


def _name_assessment(
    row: dict[str, str],
    display_name: str,
    mapping: StructureEvidence,
    source_record: _SourceSDFRecord | None,
    source_match_method: str,
    approval_verified: bool,
    source_index: _SourceSDFIndex,
    display_name_source_record: _SourceSDFRecord | None,
) -> _NameAssessment:
    general_values = _field_values(
        row,
        (*_NAME_SUPPORT_FIELDS, "brand_names", "pubchem_synonyms"),
    )
    support = _support_codes(display_name, general_values)
    source_values = [
        ("approved_source_name", name)
        for name in (source_record.names if source_record else ())
    ]
    source_support = _trusted_support_codes(display_name, source_values)
    source_is_linked = bool(
        approval_verified and source_match_method.startswith("stable_")
    )
    trusted_support = bool(source_is_linked and source_support)
    blocker = bool(
        source_index.approval_provenance_verified
        and _display_name_resolves_to_different_structure(
            mapping,
            source_record,
            display_name_source_record,
        )
    )
    reasons = _primary_name_reasons(display_name, trusted_support, blocker)
    reasons.extend(_unverified_mapping_name_reasons(row, display_name))
    reasons.extend(_auxiliary_name_reasons(row, display_name))
    return _NameAssessment(
        support_codes=sorted(set([*support, *source_support])),
        reason_codes=sorted(set(reasons)),
        corroborated=trusted_support,
        blocker=blocker,
    )


def _unique_display_name_source_record(
    display_name: str,
    source_index: _SourceSDFIndex,
) -> _SourceSDFRecord | None:
    name_key = _source_value("name", display_name)
    matches = source_index.by_identifier.get(("name", name_key), []) if name_key else []
    return matches[0] if len(matches) == 1 else None


def _display_name_resolves_to_different_structure(
    mapping: StructureEvidence,
    matched_record: _SourceSDFRecord | None,
    display_name_record: _SourceSDFRecord | None,
) -> bool:
    if display_name_record is None:
        return False
    if matched_record is not None and display_name_record.index == matched_record.index:
        return False
    return (
        _compare_structures(mapping, [display_name_record.structure])
        == PROBABLE_MISMATCH
    )


def _unverified_mapping_name_reasons(
    row: dict[str, str],
    display_name: str,
) -> list[str]:
    if not display_name:
        return []
    source_fields = (
        (
            "pubchem",
            "pubchem_cid_resolved",
            ("pubchem_record_title", "pubchem_name", "pubchem_synonyms"),
        ),
        (
            "drugcentral",
            "drugcentral_id",
            ("drugcentral_generic_name", "drugcentral_brand_names"),
        ),
    )
    reasons: list[str] = []
    for label, identifier_field, name_fields in source_fields:
        values = _field_values(row, name_fields)
        if (
            _clean(row.get(identifier_field))
            and values
            and not _support_codes(display_name, values)
        ):
            reasons.append(f"unverified_mapping_{label}_name_conflict")
    return reasons


def _primary_name_reasons(
    display_name: str, trusted_support: bool, blocker: bool
) -> list[str]:
    if not display_name:
        return ["display_name_missing"]
    if blocker:
        return ["probable_name_structure_mismatch"]
    if not trusted_support:
        return ["display_name_not_structure_corroborated"]
    return []


def _auxiliary_name_reasons(row: dict[str, str], display_name: str) -> list[str]:
    reasons: list[str] = []
    for field_name in _NAME_AUXILIARY_FIELDS:
        value = _clean(row.get(field_name))
        if value and not _name_agrees(display_name, value):
            reasons.append(f"aux_name_conflict:{field_name}")
    return reasons


def _load_curated_nonmed(path: Path | None) -> dict[str, tuple[str, ...]]:
    if path is None:
        return {}
    if not path.is_file():
        raise ValueError(f"classified non-medication CSV not found: {path}")
    accepted_statuses = {"valid", "valid fda non medication substance"}
    names_by_id: dict[str, set[str]] = defaultdict(set)
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            rdk_id = _clean(row.get("rdk_id"))
            display_name = _clean(row.get("display_name"))
            status = _norm_name(row.get("status"))
            if rdk_id and display_name and status in accepted_statuses:
                names_by_id[rdk_id].add(display_name)
    return {
        rdk_id: tuple(sorted(display_names))
        for rdk_id, display_names in names_by_id.items()
    }


def _curated_nonmed_agrees(
    curated_nonmed: dict[str, tuple[str, ...]],
    rdk_id: str,
    display_name: str,
) -> bool:
    return any(
        _name_agrees(display_name, curated_name)
        for curated_name in curated_nonmed.get(rdk_id, ())
    )


@dataclass
class _ActualEvidence:
    structures: list[StructureEvidence]
    connectivity: list[PDBQTConnectivityEvidence]
    paths: list[str]
    reason_codes: list[str]


def _recorded_pdbqt_path(row: dict[str, str], mapping_dir: Path) -> Path | None:
    raw_path = _clean(row.get("path"))
    if not raw_path:
        return None
    path = Path(raw_path).expanduser()
    if not path.is_absolute():
        path = mapping_dir / path
    return path if path.is_file() and path.suffix.casefold() == ".pdbqt" else None


def _pdbqt_evidence(
    path: Path,
    source_index: _SourceSDFIndex,
    mapping: StructureEvidence,
) -> tuple[
    list[StructureEvidence],
    list[PDBQTConnectivityEvidence],
    list[str],
]:
    structures: list[StructureEvidence] = []
    connectivity: list[PDBQTConnectivityEvidence] = []
    reasons: list[str] = []
    sidecar_structure, sidecar_reason = _sidecar_source_structure(path, source_index)
    if sidecar_reason:
        reasons.append(sidecar_reason)
    if sidecar_structure and sidecar_structure.exact_inchikey:
        structures.append(sidecar_structure)
    pdbqt_structure, pdbqt_reason, stated_key = _pdbqt_structure(path)
    if pdbqt_reason == "pdbqt_stated_inchikey_conflict" and _stated_key_matches_source(
        stated_key, sidecar_structure
    ):
        pdbqt_reason = ""
    if pdbqt_reason:
        reasons.append(pdbqt_reason)
    if pdbqt_structure and pdbqt_structure.exact_inchikey and not pdbqt_reason:
        structures.append(pdbqt_structure)
    else:
        connectivity.append(_pdbqt_connectivity_evidence(path, mapping))
    return structures, connectivity, reasons


def _collect_actual_evidence(
    row: dict[str, str],
    rdk_id: str,
    prepared_index: dict[str, list[Path]],
    source_index: _SourceSDFIndex,
    mapping_dir: Path,
    mapping: StructureEvidence,
) -> _ActualEvidence:
    candidates = list(prepared_index.get(rdk_id, []))
    recorded = _recorded_pdbqt_path(row, mapping_dir)
    if recorded and recorded not in candidates:
        candidates.append(recorded)
    reasons = ["multiple_prepared_pdbqt_candidates"] if len(candidates) > 1 else []
    structures: list[StructureEvidence] = []
    connectivity: list[PDBQTConnectivityEvidence] = []
    for path in candidates:
        path_structures, path_connectivity, path_reasons = _pdbqt_evidence(
            path,
            source_index,
            mapping,
        )
        structures.extend(path_structures)
        connectivity.extend(path_connectivity)
        reasons.extend(path_reasons)
    if _clean(row.get("remark_smiles")) or _clean(row.get("remark_inchikey")):
        reasons.append("mapping_embedded_remark_not_independent")
    return _ActualEvidence(
        structures=structures,
        connectivity=connectivity,
        paths=[str(path) for path in candidates],
        reason_codes=reasons,
    )


def _actual_structure_verdict(
    mapping: StructureEvidence, actual: _ActualEvidence
) -> tuple[str, list[str]]:
    verdicts = {
        _compare_structures(mapping, [candidate]) for candidate in actual.structures
    } - {""}
    reasons: list[str] = []
    if verdicts & {VERIFIED_EXACT, VERIFIED_PARENT} and _connectivity_has_disagreement(
        actual.connectivity
    ):
        return AMBIGUOUS, ["prepared_high_fidelity_connectivity_conflict"]
    parent_keys = {
        candidate.parent_inchikey
        for candidate in actual.structures
        if candidate.parent_inchikey
    }
    if len(parent_keys) > 1:
        reasons.append("prepared_candidate_structure_disagreement")
        return AMBIGUOUS, reasons
    if verdicts == {VERIFIED_EXACT, VERIFIED_PARENT}:
        return VERIFIED_PARENT, reasons
    if len(verdicts) > 1:
        reasons.append("prepared_candidate_structure_disagreement")
        return AMBIGUOUS, reasons
    if not verdicts:
        return _connectivity_identity_verdict(actual.connectivity)
    return next(iter(verdicts)), reasons


def _connectivity_has_disagreement(
    evidence: Sequence[PDBQTConnectivityEvidence],
) -> bool:
    return any(item.verdict == CONNECTIVITY_DISAGREES for item in evidence)


def _connectivity_identity_verdict(
    evidence: Sequence[PDBQTConnectivityEvidence],
) -> tuple[str, list[str]]:
    verdicts = {
        item.verdict for item in evidence if item.verdict != CONNECTIVITY_UNAVAILABLE
    }
    if CONNECTIVITY_DISAGREES in verdicts:
        return AMBIGUOUS, ["pdbqt_connectivity_topology_disagreement"]
    if verdicts == {CONNECTIVITY_EXACT, CONNECTIVITY_PARENT}:
        return CONNECTIVITY_PARENT, ["pdbqt_connectivity_only"]
    if verdicts == {CONNECTIVITY_EXACT}:
        return CONNECTIVITY_EXACT, ["pdbqt_connectivity_only"]
    if verdicts == {CONNECTIVITY_PARENT}:
        return CONNECTIVITY_PARENT, ["pdbqt_connectivity_only"]
    if len(verdicts) > 1:
        return AMBIGUOUS, ["prepared_candidate_connectivity_disagreement"]
    return UNVERIFIABLE, ["pdbqt_connectivity_unavailable"] if evidence else []


def _source_approval_assessment(
    mapping: StructureEvidence,
    source_structure: StructureEvidence | None,
    source_match_method: str,
    source_index: _SourceSDFIndex,
) -> tuple[bool, str, list[str]]:
    if source_structure is None:
        return _missing_source_approval(source_match_method)
    source_verdict = _compare_structures(mapping, [source_structure])
    if not source_match_method.startswith("stable_"):
        return _nonstable_source_approval(source_verdict)
    if "unmatched_id" in source_match_method:
        return False, AMBIGUOUS, ["source_identifier_unmatched"]
    if source_verdict == PROBABLE_MISMATCH:
        return False, PROBABLE_MISMATCH, ["mapping_source_sdf_structure_mismatch"]
    if source_verdict == AMBIGUOUS:
        return False, AMBIGUOUS, ["approved_source_stereo_ambiguous"]
    if _source_record_has_unsafe_parent_role(source_structure, source_verdict):
        return False, AMBIGUOUS, ["approved_source_parent_is_mixture"]
    if not source_index.approval_provenance_verified:
        return (
            False,
            "",
            [source_index.approval_provenance_reason or "approval_manifest_unverified"],
        )
    if source_verdict in {VERIFIED_EXACT, VERIFIED_PARENT}:
        return True, "", ["approved_source_structure_verified"]
    return False, "", ["approved_source_structure_unverifiable"]


def _source_record_has_unsafe_parent_role(
    source_structure: StructureEvidence,
    source_verdict: str,
) -> bool:
    if source_verdict != VERIFIED_PARENT:
        return False
    source_class = _classify_substance(source_structure, "", False)
    return source_class in {COMBINATION, METAL_COMPLEX, AMBIGUOUS_MIXTURE}


def _missing_source_approval(
    source_match_method: str,
) -> tuple[bool, str, list[str]]:
    if source_match_method == "missing_source":
        return False, "", ["approval_source_not_supplied"]
    if source_match_method.startswith(("ambiguous_", "conflicting_")):
        return False, AMBIGUOUS, ["source_identifier_conflict"]
    return False, "", ["approved_source_record_not_resolved"]


def _nonstable_source_approval(
    source_verdict: str,
) -> tuple[bool, str, list[str]]:
    reasons = ["source_sdf_ordinal_only"]
    if source_verdict == PROBABLE_MISMATCH:
        reasons.append("ordinal_source_structure_disagrees")
    return False, "", reasons


def _match_method_reason(method: str) -> str:
    normalized = method.casefold()
    if normalized == "index-fallback":
        return "ordinal_index_fallback_unverified"
    if normalized.startswith("composition"):
        return "composition_match_not_identity_proof"
    return "historical_match_method_not_independent"


def _mapping_key_integrity(
    row: dict[str, str], mapping: StructureEvidence
) -> tuple[str, list[str]]:
    stored_key = _clean(row.get("inchikey")).upper()
    recalculated_key = mapping.exact_inchikey.upper()
    if not stored_key or not recalculated_key or stored_key == recalculated_key:
        return "", []
    if stored_key.split("-", 1)[0] == recalculated_key.split("-", 1)[0]:
        return AMBIGUOUS, ["mapping_stored_inchikey_stereo_or_protonation_mismatch"]
    return PROBABLE_MISMATCH, ["mapping_stored_inchikey_connectivity_mismatch"]


def _heavy_atom_assessment(
    row: dict[str, str], mapping: StructureEvidence, verdict: str
) -> tuple[str, list[str]]:
    reasons: list[str] = []
    stored_sdf_heavy = _int_or_none(row.get("sdf_heavy_atoms"))
    pdbqt_heavy = _int_or_none(row.get("pdbqt_heavy_atoms"))
    if (
        stored_sdf_heavy
        and mapping.heavy_atoms
        and stored_sdf_heavy != mapping.heavy_atoms
    ):
        reasons.append("source_form_parent_heavy_atom_difference")
    if stored_sdf_heavy and pdbqt_heavy and stored_sdf_heavy != pdbqt_heavy:
        reasons.append("pdbqt_sdf_heavy_atom_mismatch")
        if verdict == UNVERIFIABLE:
            verdict = AMBIGUOUS
    return verdict, reasons


def _substance_reasons(
    substance_class: str, display_name: str, curated_nonmed: bool
) -> list[str]:
    reasons = {
        SALT_OR_SOLVATE: "salt_or_solvate_parent_required",
        ADDITIVE: "non_drug_substance",
        COMBINATION: "multiple_active_fragments",
        METAL_COMPLEX: "metal_complex_requires_review",
        AMBIGUOUS_MIXTURE: "mixture_role_ambiguous",
    }
    output = [reasons[substance_class]] if substance_class in reasons else []
    if not curated_nonmed and _NONMED_NAME_RE.search(display_name):
        output.append("possible_additive_or_excipient_name")
    return output


def _actual_verdict_with_method(
    mapping: StructureEvidence, actual: _ActualEvidence, match_method: str
) -> tuple[str, list[str]]:
    verdict, reasons = _actual_structure_verdict(mapping, actual)
    if verdict == UNVERIFIABLE:
        reasons.append(_match_method_reason(match_method))
    return verdict, reasons


def _prepared_form_state(verdict: str) -> str:
    return {
        VERIFIED_EXACT: "high_fidelity_exact",
        VERIFIED_PARENT: "high_fidelity_parent",
        CONNECTIVITY_EXACT: "coordinate_connectivity_exact",
        CONNECTIVITY_PARENT: "coordinate_connectivity_parent",
        PROBABLE_MISMATCH: "high_fidelity_disagrees",
        AMBIGUOUS: "ambiguous",
        UNVERIFIABLE: "unverifiable",
    }.get(verdict, "unverifiable")


def _reconcile_source_and_key_conflicts(
    verdict: str, source_conflict: str, key_conflict: str
) -> str:
    if key_conflict == PROBABLE_MISMATCH or source_conflict == PROBABLE_MISMATCH:
        return PROBABLE_MISMATCH
    if (
        key_conflict == AMBIGUOUS or source_conflict == AMBIGUOUS
    ) and verdict != PROBABLE_MISMATCH:
        return AMBIGUOUS
    return verdict


def _combination_adjusted_verdict(verdict: str, substance_class: str) -> str:
    if substance_class == COMBINATION and verdict in {
        VERIFIED_PARENT,
        CONNECTIVITY_PARENT,
    }:
        return AMBIGUOUS
    return verdict


def _row_is_eligible(
    verdict: str,
    substance_class: str,
    mapping: StructureEvidence,
    name: _NameAssessment,
    approval_verified: bool,
    reason_codes: Sequence[str],
) -> bool:
    return bool(
        verdict in {VERIFIED_EXACT, VERIFIED_PARENT}
        and substance_class in {ACTIVE, SALT_OR_SOLVATE}
        and mapping.parent_inchikey
        and name.corroborated
        and approval_verified
        and not (set(reason_codes) & _ELIGIBILITY_CONFLICT_REASONS)
    )


def _audit_row(
    row_number: int,
    row: dict[str, str],
    *,
    source_index: _SourceSDFIndex,
    prepared_index: dict[str, list[Path]],
    curated_nonmed_entries: dict[str, tuple[str, ...]],
    mapping_dir: Path,
) -> _AuditedRow:
    rdk_id = _rdk_id(row, row_number)
    display_name = _clean(row.get("display_name")) or _clean(row.get("generic_name"))
    mapping = _mapping_structure(row)
    source_record, source_match_method = _source_record_for_row(
        row, mapping, source_index
    )
    if source_index.path is None:
        source_match_method = "missing_source"
    source_structure = source_record.structure if source_record else None
    actual = _collect_actual_evidence(
        row,
        rdk_id,
        prepared_index,
        source_index,
        mapping_dir,
        mapping,
    )
    verdict, reason_codes = _actual_verdict_with_method(
        mapping, actual, _clean(row.get("match_method"))
    )
    prepared_form_state = _prepared_form_state(verdict)
    reason_codes.extend(actual.reason_codes)
    approval_verified, source_conflict, source_reasons = _source_approval_assessment(
        mapping, source_structure, source_match_method, source_index
    )
    reason_codes.extend(source_reasons)
    key_conflict, key_reasons = _mapping_key_integrity(row, mapping)
    reason_codes.extend(key_reasons)
    verdict = _reconcile_source_and_key_conflicts(
        verdict,
        source_conflict,
        key_conflict,
    )
    verdict, heavy_reasons = _heavy_atom_assessment(row, mapping, verdict)
    reason_codes.extend(heavy_reasons)

    curated_nonmed = _curated_nonmed_agrees(
        curated_nonmed_entries,
        rdk_id,
        display_name,
    )
    classification_name = display_name
    if mapping.fragment_count > 1:
        classification_name = _classification_name(
            row,
            display_name,
        )
    substance_class = _classify_substance(
        mapping,
        classification_name,
        curated_nonmed,
    )
    adjusted_verdict = _combination_adjusted_verdict(verdict, substance_class)
    if adjusted_verdict != verdict:
        reason_codes.append("prepared_structure_lost_combination_component")
    verdict = adjusted_verdict
    reason_codes.extend(
        _substance_reasons(substance_class, display_name, curated_nonmed)
    )
    display_name_source_record = _unique_display_name_source_record(
        display_name,
        source_index,
    )
    name = _name_assessment(
        row,
        display_name,
        mapping,
        source_record,
        source_match_method,
        approval_verified,
        source_index,
        display_name_source_record,
    )
    reason_codes.extend(name.reason_codes)
    if name.blocker:
        verdict = PROBABLE_MISMATCH
    eligible = _row_is_eligible(
        verdict,
        substance_class,
        mapping,
        name,
        approval_verified,
        reason_codes,
    )
    if not eligible:
        reason_codes.append("not_docking_eligible")
    return _AuditedRow(
        row_number=row_number,
        rdk_id=rdk_id,
        display_name=display_name,
        mapping_path=_clean(row.get("path")),
        match_method=_clean(row.get("match_method")),
        source_match_method=source_match_method,
        identity_verdict=verdict,
        substance_class=substance_class,
        name_corroborated=name.corroborated,
        approval_verified=approval_verified,
        fda_drug_eligible=eligible,
        reason_codes=sorted(set(reason_codes)),
        support_codes=name.support_codes,
        mapping_structure=mapping,
        source_structure=source_structure,
        source_record=source_record,
        display_name_source_record=display_name_source_record,
        actual_structures=actual.structures,
        actual_connectivity=actual.connectivity,
        actual_paths=actual.paths,
        prepared_form_state=prepared_form_state,
        curated_nonmed=curated_nonmed,
        mapping_row=row,
    )


_ROW_FIELDS = [
    "row_number",
    "rdk_id",
    "display_name",
    "category",
    "category_basis",
    "identity_verdict",
    "substance_class",
    "name_corroborated",
    "approval_verified",
    "fda_drug_eligible",
    "match_method",
    "source_match_method",
    "reason_codes",
    "support_codes",
    "mapping_exact_inchikey",
    "mapping_stored_inchikey",
    "mapping_parent_inchikey",
    "mapping_parent_connectivity_key",
    "mapping_stereo_specified",
    "mapping_parent_smiles",
    "mapping_formula",
    "mapping_heavy_atoms",
    "mapping_formal_charge",
    "mapping_fragment_count",
    "source_exact_inchikey",
    "source_parent_inchikey",
    "actual_exact_inchikeys",
    "actual_parent_inchikeys",
    "prepared_connectivity_verdicts",
    "prepared_connectivity_methods",
    "prepared_connectivity_heavy_atoms",
    "prepared_connectivity_bonds",
    "prepared_paths",
    "mapping_path",
]


def _flat_row(row: _AuditedRow) -> dict[str, object]:
    source = row.source_structure
    return {
        "row_number": row.row_number,
        "rdk_id": row.rdk_id,
        "display_name": row.display_name,
        "category": row.cross_reference_category,
        "category_basis": row.cross_reference_basis,
        "identity_verdict": row.identity_verdict,
        "substance_class": row.substance_class,
        "name_corroborated": str(row.name_corroborated).lower(),
        "approval_verified": str(row.approval_verified).lower(),
        "fda_drug_eligible": str(row.fda_drug_eligible).lower(),
        "match_method": row.match_method,
        "source_match_method": row.source_match_method,
        "reason_codes": ";".join(row.reason_codes),
        "support_codes": ";".join(row.support_codes),
        "mapping_exact_inchikey": row.mapping_structure.exact_inchikey,
        "mapping_stored_inchikey": _clean(row.mapping_row.get("inchikey")),
        "mapping_parent_inchikey": row.mapping_structure.parent_inchikey,
        "mapping_parent_connectivity_key": row.mapping_structure.parent_connectivity_key,
        "mapping_stereo_specified": str(row.mapping_structure.stereo_specified).lower(),
        "mapping_parent_smiles": row.mapping_structure.parent_smiles,
        "mapping_formula": row.mapping_structure.formula,
        "mapping_heavy_atoms": row.mapping_structure.heavy_atoms,
        "mapping_formal_charge": row.mapping_structure.formal_charge,
        "mapping_fragment_count": row.mapping_structure.fragment_count,
        "source_exact_inchikey": _structure_value(source, "exact_inchikey"),
        "source_parent_inchikey": _structure_value(source, "parent_inchikey"),
        "actual_exact_inchikeys": _structure_values(
            row.actual_structures,
            "exact_inchikey",
        ),
        "actual_parent_inchikeys": _structure_values(
            row.actual_structures,
            "parent_inchikey",
        ),
        "prepared_connectivity_verdicts": _connectivity_values(
            row.actual_connectivity,
            "verdict",
        ),
        "prepared_connectivity_methods": _connectivity_values(
            row.actual_connectivity,
            "method",
        ),
        "prepared_connectivity_heavy_atoms": _connectivity_values(
            row.actual_connectivity,
            "heavy_atoms",
        ),
        "prepared_connectivity_bonds": _connectivity_values(
            row.actual_connectivity,
            "bonds",
        ),
        "prepared_paths": ";".join(row.actual_paths),
        "mapping_path": row.mapping_path,
    }


def _structure_value(
    structure: StructureEvidence | None,
    field_name: str,
) -> str:
    return str(getattr(structure, field_name, "")) if structure else ""


def _structure_values(
    structures: Sequence[StructureEvidence],
    field_name: str,
) -> str:
    values = {str(getattr(item, field_name, "")) for item in structures}
    return ";".join(sorted(value for value in values if value))


def _connectivity_values(
    evidence: Sequence[PDBQTConnectivityEvidence],
    field_name: str,
) -> str:
    return ";".join(str(getattr(item, field_name, "")) for item in evidence)


_CROSS_REFERENCE_FIELDS = [
    "category",
    "category_basis",
    "row_number",
    "rdk_id",
    "scheme",
    "file_num",
    "sdf_index",
    "sdf_title",
    "display_name",
    "mapping_vendor_id",
    "mapping_drugcentral_id",
    "mapping_cas",
    "match_method",
    "source_match_method",
    "linked_source_record_index",
    "linked_drugcentral_id",
    "linked_drugcentral_preferred_name",
    "display_name_source_record_index",
    "display_name_drugcentral_id",
    "display_name_drugcentral_preferred_name",
    "mapping_smiles",
    "mapping_canonical_smiles",
    "mapping_exact_inchikey",
    "mapping_parent_smiles",
    "mapping_parent_inchikey",
    "substance_class",
    "mapping_fragment_count",
    "mapping_formula",
    "mapping_formal_charge",
    "prepared_form_state",
    "prepared_connectivity_state",
    "prepared_exact_inchikeys",
    "prepared_parent_inchikeys",
    "prepared_connectivity_verdicts",
    "prepared_connectivity_methods",
    "prepared_paths",
    "identity_verdict",
    "approval_verified",
    "name_corroborated",
    "fda_drug_eligible",
    "reason_codes",
    "support_codes",
    "mapping_path",
]


def _record_identifier(record: _SourceSDFRecord | None, kind: str) -> str:
    return record.identifiers.get(kind, "") if record is not None else ""


def _record_preferred_name(record: _SourceSDFRecord | None) -> str:
    return record.name if record is not None else ""


def _cross_reference_row(row: _AuditedRow) -> dict[str, object]:
    method_is_linked = bool(
        row.source_match_method.startswith("stable_")
        and "unmatched" not in row.source_match_method
        and "conflict" not in row.source_match_method
    )
    linked_record = (
        row.source_record if method_is_linked and row.approval_verified else None
    )
    display_record = row.display_name_source_record
    return {
        "category": row.cross_reference_category,
        "category_basis": row.cross_reference_basis,
        "row_number": row.row_number,
        "rdk_id": row.rdk_id,
        "scheme": _clean(row.mapping_row.get("scheme")),
        "file_num": _clean(row.mapping_row.get("file_num")),
        "sdf_index": _clean(row.mapping_row.get("sdf_index")),
        "sdf_title": _clean(row.mapping_row.get("sdf_title")),
        "display_name": row.display_name,
        "mapping_vendor_id": _clean(row.mapping_row.get("id")),
        "mapping_drugcentral_id": _clean(row.mapping_row.get("drugcentral_id")),
        "mapping_cas": _clean(row.mapping_row.get("cas")),
        "match_method": row.match_method,
        "source_match_method": row.source_match_method,
        "linked_source_record_index": linked_record.index if linked_record else "",
        "linked_drugcentral_id": _record_identifier(linked_record, "id"),
        "linked_drugcentral_preferred_name": _record_preferred_name(linked_record),
        "display_name_source_record_index": display_record.index
        if display_record
        else "",
        "display_name_drugcentral_id": _record_identifier(display_record, "id"),
        "display_name_drugcentral_preferred_name": _record_preferred_name(
            display_record
        ),
        "mapping_smiles": _clean(row.mapping_row.get("smiles")),
        "mapping_canonical_smiles": row.mapping_structure.canonical_smiles,
        "mapping_exact_inchikey": row.mapping_structure.exact_inchikey,
        "mapping_parent_smiles": row.mapping_structure.parent_smiles,
        "mapping_parent_inchikey": row.mapping_structure.parent_inchikey,
        "substance_class": row.substance_class,
        "mapping_fragment_count": row.mapping_structure.fragment_count,
        "mapping_formula": row.mapping_structure.formula,
        "mapping_formal_charge": row.mapping_structure.formal_charge,
        "prepared_form_state": row.prepared_form_state,
        "prepared_connectivity_state": _prepared_connectivity_state(row),
        "prepared_exact_inchikeys": _structure_values(
            row.actual_structures, "exact_inchikey"
        ),
        "prepared_parent_inchikeys": _structure_values(
            row.actual_structures, "parent_inchikey"
        ),
        "prepared_connectivity_verdicts": _connectivity_values(
            row.actual_connectivity, "verdict"
        ),
        "prepared_connectivity_methods": _connectivity_values(
            row.actual_connectivity, "method"
        ),
        "prepared_paths": ";".join(row.actual_paths),
        "identity_verdict": row.identity_verdict,
        "approval_verified": str(row.approval_verified).lower(),
        "name_corroborated": str(row.name_corroborated).lower(),
        "fda_drug_eligible": str(row.fda_drug_eligible).lower(),
        "reason_codes": ";".join(row.reason_codes),
        "support_codes": ";".join(row.support_codes),
        "mapping_path": row.mapping_path,
    }


def _write_csv(
    path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _eligible_parent_rows(rows: Sequence[_AuditedRow]) -> list[dict[str, object]]:
    grouped: dict[str, list[_AuditedRow]] = defaultdict(list)
    for row in rows:
        if row.fda_drug_eligible:
            grouped[row.mapping_structure.parent_inchikey].append(row)
    output: list[dict[str, object]] = []
    for parent_key in sorted(grouped):
        group = grouped[parent_key]
        output.append(
            {
                "parent_inchikey": parent_key,
                "parent_smiles": group[0].mapping_structure.parent_smiles,
                "display_names": ";".join(
                    sorted({row.display_name for row in group if row.display_name})
                ),
                "source_rdk_ids": ";".join(sorted({row.rdk_id for row in group})),
                "source_substance_classes": ";".join(
                    sorted({row.substance_class for row in group})
                ),
                "standardization_version": STANDARDIZATION_VERSION,
            }
        )
    return output


def _source_catalog_rows(rows: Sequence[_AuditedRow]) -> list[dict[str, object]]:
    return [
        {
            "rdk_id": row.rdk_id,
            "display_name": row.display_name,
            "substance_class": row.substance_class,
            "identity_verdict": row.identity_verdict,
            "prepared_connectivity_state": _prepared_connectivity_state(row),
            "name_corroborated": str(row.name_corroborated).lower(),
            "approval_verified": str(row.approval_verified).lower(),
            "fda_drug_eligible": str(row.fda_drug_eligible).lower(),
            "exact_inchikey": row.mapping_structure.exact_inchikey,
            "parent_inchikey": row.mapping_structure.parent_inchikey,
            "parent_smiles": row.mapping_structure.parent_smiles,
            "reason_codes": ";".join(row.reason_codes),
        }
        for row in rows
    ]


def _evidence_payload(row: _AuditedRow) -> dict[str, object]:
    mapping_fields = {
        key: _clean(row.mapping_row.get(key))
        for key in (
            "scheme",
            "file_num",
            "sdf_index",
            "sdf_title",
            "smiles",
            "inchikey",
            "display_name",
            "generic_name",
            "pubchem_cid_resolved",
            "pubchem_record_title",
            "rxnorm_rxcui",
            "rxnorm_generic_name",
            "drugcentral_id",
            "drugcentral_generic_name",
        )
        if _clean(row.mapping_row.get(key))
    }
    return {
        "row_number": row.row_number,
        "rdk_id": row.rdk_id,
        "identity_verdict": row.identity_verdict,
        "substance_class": row.substance_class,
        "name_corroborated": row.name_corroborated,
        "approval_verified": row.approval_verified,
        "source_match_method": row.source_match_method,
        "fda_drug_eligible": row.fda_drug_eligible,
        "reason_codes": row.reason_codes,
        "support_codes": row.support_codes,
        "mapping_fields": mapping_fields,
        "mapping_structure": asdict(row.mapping_structure),
        "source_structure": asdict(row.source_structure)
        if row.source_structure
        else None,
        "prepared_structures": [asdict(item) for item in row.actual_structures],
        "prepared_connectivity": [asdict(item) for item in row.actual_connectivity],
        "prepared_paths": row.actual_paths,
        "standardization_version": STANDARDIZATION_VERSION,
    }


def _resolve_optional_path(path: Path | None) -> Path | None:
    if path is None:
        return None
    return path.expanduser().resolve()


def _resolve_approval_manifest(
    source_sdf: Path | None,
    approval_manifest: Path | None,
) -> Path | None:
    resolved = _resolve_optional_path(approval_manifest)
    if source_sdf is None or resolved is not None:
        return resolved
    adjacent_manifest = source_sdf.parent / "fda_filter_manifest.json"
    if not adjacent_manifest.is_file():
        return None
    return adjacent_manifest.resolve()


def _resolve_score_csvs(score_csvs: Sequence[Path]) -> tuple[Path, ...]:
    resolved = tuple(
        dict.fromkeys(path.expanduser().resolve() for path in score_csvs)
    )
    missing_score_csvs = [path for path in resolved if not path.is_file()]
    if not missing_score_csvs:
        return resolved
    missing = ", ".join(str(path) for path in missing_score_csvs)
    raise ValueError(f"FDA score CSV not found: {missing}")


def _resolve_config(config: FDAIdentityAuditConfig) -> _ResolvedConfig:
    source_sdf = _resolve_optional_path(config.source_sdf)
    approval_manifest = _resolve_approval_manifest(
        source_sdf, config.approval_manifest
    )
    return _ResolvedConfig(
        mapping_csv=Path(config.mapping_csv).expanduser().resolve(),
        output_dir=Path(config.output_dir).expanduser().resolve(),
        source_sdf=source_sdf,
        approval_manifest=approval_manifest,
        library_dir=_resolve_optional_path(config.library_dir),
        classified_nonmed_csv=_resolve_optional_path(config.classified_nonmed_csv),
        score_csvs=_resolve_score_csvs(config.score_csvs),
    )


def _read_mapping_rows(mapping_csv: Path) -> list[dict[str, str]]:
    if not mapping_csv.is_file():
        raise ValueError(f"FDA mapping CSV not found: {mapping_csv}")
    with mapping_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            raise ValueError(f"FDA mapping CSV has no header: {mapping_csv}")
        if not any(name in fieldnames for name in ("path", "rdk_id", "file_num")):
            raise ValueError("mapping has no ligand identifier column")
        return [{key: _clean(value) for key, value in row.items()} for row in reader]


def _audit_mapping_rows(
    rows: Sequence[dict[str, str]], config: _ResolvedConfig
) -> list[_AuditedRow]:
    with rdBase.BlockLogs():
        source_index = _read_source_sdf(
            config.source_sdf,
            config.approval_manifest,
        )
        prepared_index = _library_index(config.library_dir)
        curated_nonmed_entries = _load_curated_nonmed(config.classified_nonmed_csv)
        return [
            _audit_row(
                row_number,
                row,
                source_index=source_index,
                prepared_index=prepared_index,
                curated_nonmed_entries=curated_nonmed_entries,
                mapping_dir=config.mapping_csv.parent,
            )
            for row_number, row in enumerate(rows, start=1)
        ]


def _mark_duplicate_ids(rows: Sequence[_AuditedRow]) -> None:
    duplicate_ids = {
        key for key, count in Counter(row.rdk_id for row in rows).items() if count > 1
    }
    for row in rows:
        if row.rdk_id not in duplicate_ids:
            continue
        row.reason_codes = sorted(set([*row.reason_codes, "duplicate_ligand_id"]))
        row.fda_drug_eligible = False
        if row.identity_verdict != PROBABLE_MISMATCH:
            row.identity_verdict = AMBIGUOUS


_PREPARED_EXACT_STATES = frozenset(
    {"high_fidelity_exact", "coordinate_connectivity_exact"}
)
_PREPARED_PARENT_STATES = frozenset(
    {"high_fidelity_parent", "coordinate_connectivity_parent"}
)
_PREPARED_CONSISTENT_STATES = _PREPARED_EXACT_STATES | _PREPARED_PARENT_STATES


def _cross_reference_category(row: _AuditedRow) -> tuple[str, str]:
    classifiers = (
        _incorrect_mapping_category,
        _curated_additive_category,
        _collapsed_combination_category,
        _retained_salt_category,
        _confirmed_salt_parent_category,
        _confirmed_active_category,
    )
    for classifier in classifiers:
        assignment = classifier(row)
        if assignment is not None:
            return assignment
    return UNVERIFIABLE_CATEGORY, _unverifiable_category_basis(row)


def _incorrect_mapping_category(row: _AuditedRow) -> tuple[str, str] | None:
    if "probable_name_structure_mismatch" in row.reason_codes:
        return (
            INCORRECT_MAPPING_CATEGORY,
            "display_name_resolves_to_different_approved_structure",
        )
    if "mapping_source_sdf_structure_mismatch" in row.reason_codes:
        return (
            INCORRECT_MAPPING_CATEGORY,
            "stable_identifier_resolves_to_different_approved_structure",
        )
    return None


def _curated_additive_category(row: _AuditedRow) -> tuple[str, str] | None:
    if row.curated_nonmed:
        return ADDITIVE_CATEGORY, "curated_status_and_name_match"
    return None


def _collapsed_combination_category(row: _AuditedRow) -> tuple[str, str] | None:
    if (
        row.substance_class == COMBINATION
        and row.prepared_form_state in _PREPARED_PARENT_STATES
    ):
        return (
            COLLAPSED_COMBINATION_CATEGORY,
            f"multi_active_{row.prepared_form_state}",
        )
    return None


def _is_explicit_salt(row: _AuditedRow) -> bool:
    return bool(
        row.substance_class == SALT_OR_SOLVATE
        and row.mapping_structure.fragment_count > 1
    )


def _retained_salt_category(row: _AuditedRow) -> tuple[str, str] | None:
    if _is_explicit_salt(row) and row.prepared_form_state in _PREPARED_EXACT_STATES:
        return (
            RETAINED_SALT_CATEGORY,
            f"explicit_multifragment_salt_{row.prepared_form_state}",
        )
    return None


def _is_manifest_name_confirmed(row: _AuditedRow) -> bool:
    return row.approval_verified and row.name_corroborated


def _confirmed_salt_parent_category(row: _AuditedRow) -> tuple[str, str] | None:
    if (
        _is_manifest_name_confirmed(row)
        and _is_explicit_salt(row)
        and row.prepared_form_state in _PREPARED_PARENT_STATES
    ):
        return (
            CONFIRMED_SALT_PARENT_CATEGORY,
            f"manifest_name_confirmed_explicit_salt_{row.prepared_form_state}",
        )
    return None


def _confirmed_active_category(row: _AuditedRow) -> tuple[str, str] | None:
    if (
        _is_manifest_name_confirmed(row)
        and row.substance_class == ACTIVE
        and row.prepared_form_state in _PREPARED_CONSISTENT_STATES
    ):
        return (
            CONFIRMED_ACTIVE_CATEGORY,
            f"manifest_name_confirmed_active_{row.prepared_form_state}",
        )
    return None


def _unverifiable_category_basis(row: _AuditedRow) -> str:
    tokens: list[str] = []
    if not row.approval_verified:
        tokens.append("manifest_fda_structure_not_verified")
    if not row.name_corroborated:
        tokens.append("drugcentral_name_not_corroborated")
    tokens.append(f"prepared_{row.prepared_form_state}")
    if (
        row.substance_class == SALT_OR_SOLVATE
        and row.mapping_structure.fragment_count <= 1
    ):
        tokens.append("salt_name_without_explicit_counterion_fragment")
    if row.substance_class in {COMBINATION, METAL_COMPLEX, AMBIGUOUS_MIXTURE}:
        tokens.append("metal_or_mixture_requires_review")
    source_method = re.sub(
        r"[^a-z0-9]+",
        "_",
        row.source_match_method.casefold(),
    ).strip("_")
    tokens.append(f"source_method_{source_method or 'missing'}")
    return ";".join(dict.fromkeys(tokens))


def _assign_cross_reference_categories(rows: Sequence[_AuditedRow]) -> None:
    for row in rows:
        category, basis = _cross_reference_category(row)
        if category not in FDA_CROSS_REFERENCE_CATEGORIES:
            raise RuntimeError(f"unknown FDA cross-reference category: {category}")
        row.cross_reference_category = category
        row.cross_reference_basis = basis


def _audit_output_paths(output_dir: Path) -> _AuditOutputPaths:
    return _AuditOutputPaths(
        row_audit_csv=output_dir / "fda_identity_audit.csv",
        evidence_jsonl=output_dir / "fda_identity_evidence.jsonl",
        summary_json=output_dir / "fda_identity_summary.json",
        review_csv=output_dir / "fda_identity_review.csv",
        quarantine_csv=output_dir / "fda_identity_quarantine.csv",
        eligible_csv=output_dir / "fda_docking_eligible_active_parents.csv",
        source_catalog_csv=output_dir / "fda_source_substance_catalog.csv",
        cross_reference_csv=output_dir / "fda_approval_cross_reference.csv",
    )


def _validate_output_scope(config: _ResolvedConfig, outputs: _AuditOutputPaths) -> None:
    input_files = {
        path
        for path in (
            config.mapping_csv,
            config.source_sdf,
            config.approval_manifest,
            config.classified_nonmed_csv,
        )
        if path is not None
    }
    if any(path.resolve() in input_files for path in asdict(outputs).values()):
        raise ValueError("audit output path collides with an input file")
    if config.library_dir and config.output_dir.is_relative_to(config.library_dir):
        raise ValueError(
            "audit output directory must not be inside the prepared library"
        )


def _write_evidence_jsonl(path: Path, rows: Sequence[_AuditedRow]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(_evidence_payload(row), sort_keys=True) + "\n")


def _write_audit_artifacts(
    rows: Sequence[_AuditedRow], outputs: _AuditOutputPaths
) -> None:
    flat_rows = [_flat_row(row) for row in rows]
    _write_csv(outputs.row_audit_csv, flat_rows, _ROW_FIELDS)
    _write_csv(
        outputs.review_csv,
        [_flat_row(row) for row in rows if row.review],
        _ROW_FIELDS,
    )
    _write_csv(
        outputs.quarantine_csv,
        [_flat_row(row) for row in rows if row.blocker],
        _ROW_FIELDS,
    )
    eligible_fields = [
        "parent_inchikey",
        "parent_smiles",
        "display_names",
        "source_rdk_ids",
        "source_substance_classes",
        "standardization_version",
    ]
    _write_csv(outputs.eligible_csv, _eligible_parent_rows(rows), eligible_fields)
    catalog_fields = [
        "rdk_id",
        "display_name",
        "substance_class",
        "identity_verdict",
        "prepared_connectivity_state",
        "name_corroborated",
        "approval_verified",
        "fda_drug_eligible",
        "exact_inchikey",
        "parent_inchikey",
        "parent_smiles",
        "reason_codes",
    ]
    _write_csv(outputs.source_catalog_csv, _source_catalog_rows(rows), catalog_fields)
    _write_csv(
        outputs.cross_reference_csv,
        [_cross_reference_row(row) for row in rows],
        _CROSS_REFERENCE_FIELDS,
    )
    _write_evidence_jsonl(outputs.evidence_jsonl, rows)


def _sha256(path: Path | None) -> str:
    if path is None or not path.is_file():
        return ""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _summary_counts(rows: Sequence[_AuditedRow]) -> dict[str, object]:
    verdict_counts: Counter[str] = Counter()
    connectivity_counts: Counter[str] = Counter()
    class_counts: Counter[str] = Counter()
    category_counts: Counter[str] = Counter(
        {category: 0 for category in FDA_CROSS_REFERENCE_CATEGORIES}
    )
    reason_counts: Counter[str] = Counter()
    scalar_counts: Counter[str] = Counter()
    for row in rows:
        verdict_counts[row.identity_verdict] += 1
        connectivity_state = _prepared_connectivity_state(row)
        connectivity_counts[connectivity_state] += 1
        class_counts[row.substance_class] += 1
        if row.cross_reference_category not in FDA_CROSS_REFERENCE_CATEGORIES:
            raise RuntimeError(
                f"unknown FDA cross-reference category: {row.cross_reference_category}"
            )
        category_counts[row.cross_reference_category] += 1
        reason_counts.update(row.reason_codes)
        scalar_counts["blocker"] += int(row.blocker)
        scalar_counts["review"] += int(row.review)
        scalar_counts["name_corroborated"] += int(row.name_corroborated)
        scalar_counts["approval_verified"] += int(row.approval_verified)
        scalar_counts["eligible"] += int(row.fda_drug_eligible)
        scalar_counts["mapping_structure_parsed"] += int(
            bool(row.mapping_structure.canonical_smiles)
        )
        scalar_counts["prepared_pdbqt"] += int(bool(row.actual_paths))
        scalar_counts["high_fidelity_prepared_structure"] += int(
            bool(row.actual_structures)
        )
        scalar_counts["connectivity_consistent"] += int(
            connectivity_state in {CONNECTIVITY_EXACT, CONNECTIVITY_PARENT}
        )
        scalar_counts["connectivity_evaluable"] += int(
            connectivity_state
            in {
                CONNECTIVITY_EXACT,
                CONNECTIVITY_PARENT,
                CONNECTIVITY_DISAGREES,
            }
        )
    serialized_category_counts = {
        category: int(category_counts[category])
        for category in FDA_CROSS_REFERENCE_CATEGORIES
    }
    if (
        len(serialized_category_counts) != 7
        or sum(serialized_category_counts.values()) != len(rows)
    ):
        raise RuntimeError("FDA cross-reference categories are not mutually exhaustive")
    return {
        "verdict_counts": dict(sorted(verdict_counts.items())),
        "connectivity_counts": dict(sorted(connectivity_counts.items())),
        "class_counts": dict(sorted(class_counts.items())),
        "category_counts": serialized_category_counts,
        "reason_counts": dict(sorted(reason_counts.items())),
        "scalar_counts": dict(scalar_counts),
    }


def _prepared_connectivity_state(row: _AuditedRow) -> str:
    verdicts = {item.verdict for item in row.actual_connectivity}
    if CONNECTIVITY_DISAGREES in verdicts:
        return CONNECTIVITY_DISAGREES
    if verdicts == {CONNECTIVITY_EXACT, CONNECTIVITY_PARENT}:
        return CONNECTIVITY_PARENT
    if CONNECTIVITY_EXACT in verdicts:
        return CONNECTIVITY_EXACT
    if CONNECTIVITY_PARENT in verdicts:
        return CONNECTIVITY_PARENT
    if CONNECTIVITY_UNAVAILABLE in verdicts:
        return CONNECTIVITY_UNAVAILABLE
    if row.actual_structures:
        return "high_fidelity_structure_available"
    return "not_supplied"


def _summary_output_paths(outputs: _AuditOutputPaths) -> dict[str, str]:
    return {
        "row_audit_csv": str(outputs.row_audit_csv),
        "evidence_jsonl": str(outputs.evidence_jsonl),
        "review_csv": str(outputs.review_csv),
        "quarantine_manifest_csv": str(outputs.quarantine_csv),
        "eligible_manifest_csv": str(outputs.eligible_csv),
        "source_catalog_csv": str(outputs.source_catalog_csv),
        "approval_cross_reference_csv": str(outputs.cross_reference_csv),
    }


def _prepared_pdbqt_inventory(library_dir: Path | None) -> dict[str, object]:
    if library_dir is None:
        return {"count": 0, "sha256": ""}
    paths = sorted(path for path in library_dir.rglob("*.pdbqt") if path.is_file())
    digest = hashlib.sha256()
    for path in paths:
        relative = path.relative_to(library_dir).as_posix()
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        digest.update(b"\0")
    return {"count": len(paths), "sha256": digest.hexdigest()}


def _summary_payload(
    rows: Sequence[_AuditedRow],
    config: _ResolvedConfig,
    outputs: _AuditOutputPaths,
    repair_summary: Mapping[str, object] | None = None,
) -> dict[str, object]:
    counts = _summary_counts(rows)
    scalar = counts["scalar_counts"]
    assert isinstance(scalar, dict)
    return {
        "schema_version": 4,
        "fda_repair_schema_version": 1,
        "fda_approval_cross_reference_schema_version": 1,
        "standardization_version": STANDARDIZATION_VERSION,
        "rdkit_version": rdBase.rdkitVersion,
        "mapping_csv": str(config.mapping_csv),
        "source_sdf": str(config.source_sdf or ""),
        "approval_manifest": str(config.approval_manifest or ""),
        "library_dir": str(config.library_dir or ""),
        "classified_nonmed_csv": str(config.classified_nonmed_csv or ""),
        "input_sha256": {
            "mapping_csv": _sha256(config.mapping_csv),
            "source_sdf": _sha256(config.source_sdf),
            "approval_manifest": _sha256(config.approval_manifest),
            "classified_nonmed_csv": _sha256(config.classified_nonmed_csv),
        },
        "total_rows": len(rows),
        "blocker_count": scalar.get("blocker", 0),
        "review_count": scalar.get("review", 0),
        "name_corroborated_rows": scalar.get("name_corroborated", 0),
        "approval_verified_rows": scalar.get("approval_verified", 0),
        "docking_eligible_source_rows": scalar.get("eligible", 0),
        "docking_eligible_unique_parents": len(_eligible_parent_rows(rows)),
        "identity_verdict_counts": counts["verdict_counts"],
        "prepared_pdbqt_connectivity_counts": counts["connectivity_counts"],
        "prepared_pdbqt_inventory": _prepared_pdbqt_inventory(config.library_dir),
        "mapping_structure_parsed_rows": scalar.get("mapping_structure_parsed", 0),
        "prepared_pdbqt_rows": scalar.get("prepared_pdbqt", 0),
        "high_fidelity_prepared_structure_rows": scalar.get(
            "high_fidelity_prepared_structure", 0
        ),
        "prepared_connectivity_consistent_rows": scalar.get(
            "connectivity_consistent", 0
        ),
        "prepared_connectivity_evaluable_rows": scalar.get("connectivity_evaluable", 0),
        "substance_class_counts": counts["class_counts"],
        "fda_approval_cross_reference_category_counts": counts["category_counts"],
        "reason_counts": counts["reason_counts"],
        "outputs": _summary_output_paths(outputs),
        "repair": dict(repair_summary or {}),
        "policy": {
            "canonical_mapping_mutated": False,
            "quarantine_verdicts": [PROBABLE_MISMATCH],
            "eligible_substance_classes": [ACTIVE, SALT_OR_SOLVATE],
            "eligibility_requires_name_corroboration": True,
            "eligibility_requires_stable_approved_source_record": True,
            "eligibility_requires_valid_fda_filter_manifest": True,
            "pdbqt_connectivity_only_cannot_establish_eligibility": True,
            "index_fallback_is_verified": False,
        },
    }


def run_fda_identity_audit(config: FDAIdentityAuditConfig) -> FDAIdentityAuditResult:
    """Audit a mapping and write deterministic evidence without mutating inputs."""

    resolved = _resolve_config(config)
    mapping_rows = _read_mapping_rows(resolved.mapping_csv)
    audited = _audit_mapping_rows(mapping_rows, resolved)
    _mark_duplicate_ids(audited)
    _assign_cross_reference_categories(audited)
    outputs = _audit_output_paths(resolved.output_dir)
    _validate_output_scope(resolved, outputs)
    resolved.output_dir.mkdir(parents=True, exist_ok=True)
    _write_audit_artifacts(audited, outputs)
    with rdBase.BlockLogs():
        source_index = _read_source_sdf(
            resolved.source_sdf,
            resolved.approval_manifest,
        )
    source_class_by_index = {
        int(record.index): _classify_substance(record.structure, record.name, False)
        for record in source_index.by_index.values()
    }
    from prep_ligands.fda_repair_manifest import write_fda_repair_artifacts

    repair_outputs, repair_summary = write_fda_repair_artifacts(
        audited_rows=audited,
        source_index=source_index,
        source_class_by_index=source_class_by_index,
        output_dir=resolved.output_dir,
        score_csvs=resolved.score_csvs,
        protected_input_paths=tuple(
            path
            for path in (
                resolved.mapping_csv,
                resolved.source_sdf,
                resolved.approval_manifest,
                resolved.classified_nonmed_csv,
            )
            if path is not None
        ),
    )
    summary = _summary_payload(audited, resolved, outputs, repair_summary)
    outputs.summary_json.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return FDAIdentityAuditResult(
        row_audit_csv=outputs.row_audit_csv,
        evidence_jsonl=outputs.evidence_jsonl,
        summary_json=outputs.summary_json,
        review_csv=outputs.review_csv,
        quarantine_manifest_csv=outputs.quarantine_csv,
        eligible_manifest_csv=outputs.eligible_csv,
        source_catalog_csv=outputs.source_catalog_csv,
        approval_cross_reference_csv=outputs.cross_reference_csv,
        total_rows=len(audited),
        blocker_count=int(str(summary["blocker_count"])),
        review_count=int(str(summary["review_count"])),
        repaired_mapping_csv=repair_outputs.repaired_mapping_csv,
        score_reuse_manifest_csv=repair_outputs.score_reuse_manifest_csv,
        redock_delta_csv=repair_outputs.redock_delta_csv,
        repair_quarantine_csv=repair_outputs.repair_quarantine_csv,
    )


__all__ = [
    "FDAIdentityAuditConfig",
    "FDAIdentityAuditResult",
    "run_fda_identity_audit",
]
