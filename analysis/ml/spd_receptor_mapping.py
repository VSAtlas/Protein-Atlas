"""Versioned, fail-closed SPD receptor-target mapping contract."""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import pandas as pd


POLICY_VERSION = "spd_receptor_target_mapping_v1.0.0"
DEFAULT_CONTRACT_PATH = (
    Path(__file__).resolve().parents[1]
    / "gene_list"
    / "spd_receptor_target_mapping_v1.json"
)
OUTPUT_PREFIX = "spd_receptor_mapping_"
OUTPUT_FIELDS = (
    "policy_version",
    "status",
    "strict_eligible",
    "exclusion_reason",
    "contract_source",
    "contract_sha256",
    "changed",
    "requires_target_derived_rebuild",
    "legacy_target_id",
    "legacy_target_gene",
    "legacy_target_uniprot",
    "revised_target_id",
    "revised_target_gene",
    "revised_target_uniprot",
    "selected_chain",
    "site_chain",
)
OUTPUT_COLUMNS = tuple(f"{OUTPUT_PREFIX}{field}" for field in OUTPUT_FIELDS)

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_REQUIRED_EVIDENCE_KINDS = frozenset(
    {
        "audit_mapping_review",
        "pdbe_sifts_snapshot",
        "uniprot_snapshot",
        "local_structure",
        "prepared_receptor",
    }
)
_LOCKED_MAPPINGS: dict[str, dict[str, object]] = {
    "9I52": {
        "status": "validated_revision",
        "strict_eligible": True,
        "changed": True,
        "requires_target_derived_rebuild": True,
        "legacy_target_id": "ADRB2",
        "legacy_target_gene": "ADRB2",
        "legacy_target_uniprot": "P07550",
        "revised_target_id": "DRD1",
        "revised_target_gene": "DRD1",
        "revised_target_uniprot": "P21728",
        "selected_chain": "R",
        "site_chain": "R",
        "prepared_chains": ("A", "B", "G", "R"),
        "site_relevant_chains": ("R",),
    },
    "6HUJ": {
        "status": "unresolved",
        "strict_eligible": False,
        "changed": False,
        "requires_target_derived_rebuild": False,
        "legacy_target_id": "GABRA1",
        "legacy_target_gene": "GABRA1",
        "legacy_target_uniprot": "P14867",
        "revised_target_id": None,
        "revised_target_gene": None,
        "revised_target_uniprot": None,
        "selected_chain": None,
        "site_chain": None,
        "prepared_chains": ("A", "B", "C", "D", "E", "G"),
        "site_relevant_chains": ("A", "B", "C", "D", "E"),
    },
    "8HCQ": {
        "status": "unresolved",
        "strict_eligible": False,
        "changed": False,
        "requires_target_derived_rebuild": False,
        "legacy_target_id": "EDNRA",
        "legacy_target_gene": "EDNRA",
        "legacy_target_uniprot": "P25101",
        "revised_target_id": None,
        "revised_target_gene": None,
        "revised_target_uniprot": None,
        "selected_chain": None,
        "site_chain": None,
        "prepared_chains": ("A", "B", "E", "G", "L", "R"),
        "site_relevant_chains": ("L", "R"),
    },
}
_LOCKED_9I52_SITE_EVIDENCE: dict[str, object] = {
    "ligand": "A1I",
    "ligand_atom_count": 26,
    "contact_cutoff_angstrom": 5.0,
    "target_gene": "DRD1",
    "target_uniprot": "P21728",
    "chain": "R",
    "ligand_atoms_within_cutoff": 26,
    "atom_pairs_within_cutoff": 242,
    "minimum_distance_angstrom": 2.9937418,
}


class ReceptorMappingContractError(ValueError):
    """Raised when the receptor-target mapping contract is not safe to use."""


@dataclass(frozen=True)
class EvidenceReference:
    kind: str
    path: str
    sha256: str


@dataclass(frozen=True)
class SPDReceptorMapping:
    pdb_id: str
    status: str
    strict_eligible: bool
    exclusion_reason: str
    changed: bool
    requires_target_derived_rebuild: bool
    legacy_target_id: str
    legacy_target_gene: str
    legacy_target_uniprot: str
    revised_target_id: str | None
    revised_target_gene: str | None
    revised_target_uniprot: str | None
    selected_chain: str | None
    site_chain: str | None
    prepared_chains: tuple[str, ...]
    site_relevant_chains: tuple[str, ...]
    ligand_site_evidence: Mapping[str, object] | None
    evidence: tuple[EvidenceReference, ...]


@dataclass(frozen=True)
class SPDReceptorMappingContract:
    schema_version: int
    policy_version: str
    evidence_retrieved_at_utc: str
    mappings: tuple[SPDReceptorMapping, ...]
    source: str
    sha256: str

    def by_pdb(self) -> dict[str, SPDReceptorMapping]:
        return {mapping.pdb_id: mapping for mapping in self.mappings}


def _canonical_payload_sha256(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _required_text(row: Mapping[str, Any], field: str, *, context: str) -> str:
    value = row.get(field)
    if not isinstance(value, str) or not value.strip():
        raise ReceptorMappingContractError(
            f"{context}: {field} must be a nonempty string"
        )
    return value.strip()


def _optional_text(row: Mapping[str, Any], field: str, *, context: str) -> str | None:
    value = row.get(field)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ReceptorMappingContractError(
            f"{context}: {field} must be null or a nonempty string"
        )
    return value.strip()


def _required_bool(row: Mapping[str, Any], field: str, *, context: str) -> bool:
    value = row.get(field)
    if not isinstance(value, bool):
        raise ReceptorMappingContractError(f"{context}: {field} must be boolean")
    return value


def _required_chain_list(
    row: Mapping[str, Any],
    field: str,
    *,
    context: str,
) -> tuple[str, ...]:
    value = row.get(field)
    if not isinstance(value, list) or not value:
        raise ReceptorMappingContractError(
            f"{context}: {field} must be a nonempty list"
        )
    chains: list[str] = []
    for chain in value:
        if not isinstance(chain, str) or not chain.strip():
            raise ReceptorMappingContractError(
                f"{context}: {field} entries must be nonempty strings"
            )
        chains.append(chain.strip())
    if len(chains) != len(set(chains)):
        raise ReceptorMappingContractError(
            f"{context}: {field} cannot contain duplicate chains"
        )
    return tuple(chains)


def _validate_evidence(
    values: Any,
    *,
    context: str,
) -> tuple[EvidenceReference, ...]:
    if not isinstance(values, list) or not values:
        raise ReceptorMappingContractError(
            f"{context}: evidence must be a nonempty list"
        )
    evidence: list[EvidenceReference] = []
    seen_kinds: set[str] = set()
    for index, raw in enumerate(values):
        evidence_context = f"{context}.evidence[{index}]"
        if not isinstance(raw, Mapping):
            raise ReceptorMappingContractError(
                f"{evidence_context}: evidence entry must be an object"
            )
        kind = _required_text(raw, "kind", context=evidence_context)
        path = _required_text(raw, "path", context=evidence_context)
        sha256 = _required_text(raw, "sha256", context=evidence_context).lower()
        if kind in seen_kinds:
            raise ReceptorMappingContractError(
                f"{context}: duplicate evidence kind {kind!r}"
            )
        if not _SHA256_RE.fullmatch(sha256):
            raise ReceptorMappingContractError(
                f"{evidence_context}: sha256 must be 64 lowercase hexadecimal characters"
            )
        seen_kinds.add(kind)
        evidence.append(EvidenceReference(kind=kind, path=path, sha256=sha256))

    missing = sorted(_REQUIRED_EVIDENCE_KINDS - seen_kinds)
    unexpected = sorted(seen_kinds - _REQUIRED_EVIDENCE_KINDS)
    if missing or unexpected:
        raise ReceptorMappingContractError(
            f"{context}: invalid evidence kinds; missing={missing}, unexpected={unexpected}"
        )
    return tuple(evidence)


def _validate_locked_identity(
    pdb_id: str,
    values: Mapping[str, object],
    *,
    context: str,
) -> None:
    expected = _LOCKED_MAPPINGS.get(pdb_id)
    if expected is None:
        raise ReceptorMappingContractError(
            f"{context}: unexpected PDB identity {pdb_id!r}"
        )
    for field, expected_value in expected.items():
        if values[field] != expected_value:
            raise ReceptorMappingContractError(
                f"{context}: unexpected {field} for {pdb_id}; "
                f"expected {expected_value!r}, observed {values[field]!r}"
            )


def validate_spd_receptor_mapping_contract(
    payload: Mapping[str, Any],
    *,
    source: str = "<memory>",
    contract_sha256: str | None = None,
) -> SPDReceptorMappingContract:
    """Validate and parse a mapping payload, failing closed on any policy drift."""

    if not isinstance(payload, Mapping):
        raise ReceptorMappingContractError("contract payload must be an object")
    schema_version = payload.get("schema_version")
    if schema_version != 1:
        raise ReceptorMappingContractError(
            f"unsupported schema_version {schema_version!r}; expected 1"
        )
    policy_version = _required_text(payload, "policy_version", context="contract")
    if policy_version != POLICY_VERSION:
        raise ReceptorMappingContractError(
            f"unexpected policy_version {policy_version!r}; expected {POLICY_VERSION!r}"
        )
    retrieved_at = _required_text(
        payload,
        "evidence_retrieved_at_utc",
        context="contract",
    )
    raw_mappings = payload.get("mappings")
    if not isinstance(raw_mappings, list) or not raw_mappings:
        raise ReceptorMappingContractError("contract.mappings must be a nonempty list")

    mappings: list[SPDReceptorMapping] = []
    seen_pdbs: set[str] = set()
    for index, raw in enumerate(raw_mappings):
        context = f"contract.mappings[{index}]"
        if not isinstance(raw, Mapping):
            raise ReceptorMappingContractError(f"{context} must be an object")
        pdb_id = _required_text(raw, "pdb_id", context=context).upper()
        if pdb_id in seen_pdbs:
            raise ReceptorMappingContractError(f"duplicate PDB mapping {pdb_id}")
        seen_pdbs.add(pdb_id)

        status = _required_text(raw, "status", context=context)
        strict_eligible = _required_bool(raw, "strict_eligible", context=context)
        exclusion_reason = raw.get("exclusion_reason")
        if not isinstance(exclusion_reason, str):
            raise ReceptorMappingContractError(
                f"{context}: exclusion_reason must be a string"
            )
        changed = _required_bool(raw, "changed", context=context)
        rebuild = _required_bool(
            raw,
            "requires_target_derived_rebuild",
            context=context,
        )
        values: dict[str, object] = {
            "status": status,
            "strict_eligible": strict_eligible,
            "changed": changed,
            "requires_target_derived_rebuild": rebuild,
            "legacy_target_id": _required_text(
                raw, "legacy_target_id", context=context
            ),
            "legacy_target_gene": _required_text(
                raw, "legacy_target_gene", context=context
            ),
            "legacy_target_uniprot": _required_text(
                raw, "legacy_target_uniprot", context=context
            ),
            "revised_target_id": _optional_text(
                raw, "revised_target_id", context=context
            ),
            "revised_target_gene": _optional_text(
                raw, "revised_target_gene", context=context
            ),
            "revised_target_uniprot": _optional_text(
                raw, "revised_target_uniprot", context=context
            ),
            "selected_chain": _optional_text(
                raw, "selected_chain", context=context
            ),
            "site_chain": _optional_text(raw, "site_chain", context=context),
            "prepared_chains": _required_chain_list(
                raw, "prepared_chains", context=context
            ),
            "site_relevant_chains": _required_chain_list(
                raw, "site_relevant_chains", context=context
            ),
        }
        site_evidence = raw.get("ligand_site_evidence")
        if site_evidence is not None and not isinstance(site_evidence, Mapping):
            raise ReceptorMappingContractError(
                f"{context}: ligand_site_evidence must be null or an object"
            )

        revised_identity = (
            values["revised_target_id"],
            values["revised_target_gene"],
            values["revised_target_uniprot"],
        )
        if status == "unresolved":
            if any(value is not None for value in revised_identity):
                raise ReceptorMappingContractError(
                    f"{context}: unresolved mapping cannot declare a revised identity"
                )
            if strict_eligible:
                raise ReceptorMappingContractError(
                    f"{context}: unresolved mapping cannot be strict eligible"
                )
            if not exclusion_reason.strip():
                raise ReceptorMappingContractError(
                    f"{context}: unresolved mapping requires an exclusion_reason"
                )
        elif status == "validated_revision":
            if any(value is None for value in revised_identity):
                raise ReceptorMappingContractError(
                    f"{context}: validated revision requires a complete revised identity"
                )
            if not strict_eligible:
                raise ReceptorMappingContractError(
                    f"{context}: validated revision must be strict eligible"
                )
            if not values["selected_chain"] or not values["site_chain"]:
                raise ReceptorMappingContractError(
                    f"{context}: validated revision requires selected and site chains"
                )
            if dict(site_evidence or {}) != _LOCKED_9I52_SITE_EVIDENCE:
                raise ReceptorMappingContractError(
                    f"{context}: validated 9I52 revision requires the locked ligand-site evidence"
                )
        else:
            raise ReceptorMappingContractError(
                f"{context}: unexpected status {status!r}"
            )

        _validate_locked_identity(pdb_id, values, context=context)
        if status == "unresolved" and site_evidence is not None:
            raise ReceptorMappingContractError(
                f"{context}: unresolved mapping cannot claim resolved ligand-site evidence"
            )
        evidence = _validate_evidence(raw.get("evidence"), context=context)
        mappings.append(
            SPDReceptorMapping(
                pdb_id=pdb_id,
                status=status,
                strict_eligible=strict_eligible,
                exclusion_reason=exclusion_reason.strip(),
                changed=changed,
                requires_target_derived_rebuild=rebuild,
                legacy_target_id=str(values["legacy_target_id"]),
                legacy_target_gene=str(values["legacy_target_gene"]),
                legacy_target_uniprot=str(values["legacy_target_uniprot"]),
                revised_target_id=values["revised_target_id"],  # type: ignore[arg-type]
                revised_target_gene=values["revised_target_gene"],  # type: ignore[arg-type]
                revised_target_uniprot=values["revised_target_uniprot"],  # type: ignore[arg-type]
                selected_chain=values["selected_chain"],  # type: ignore[arg-type]
                site_chain=values["site_chain"],  # type: ignore[arg-type]
                prepared_chains=values["prepared_chains"],  # type: ignore[arg-type]
                site_relevant_chains=values["site_relevant_chains"],  # type: ignore[arg-type]
                ligand_site_evidence=dict(site_evidence) if site_evidence else None,
                evidence=evidence,
            )
        )

    expected_pdbs = set(_LOCKED_MAPPINGS)
    if seen_pdbs != expected_pdbs:
        raise ReceptorMappingContractError(
            "contract PDB scope mismatch; "
            f"missing={sorted(expected_pdbs - seen_pdbs)}, "
            f"unexpected={sorted(seen_pdbs - expected_pdbs)}"
        )

    digest = contract_sha256 or _canonical_payload_sha256(payload)
    if not _SHA256_RE.fullmatch(digest):
        raise ReceptorMappingContractError("contract_sha256 is not a valid SHA-256")
    return SPDReceptorMappingContract(
        schema_version=1,
        policy_version=policy_version,
        evidence_retrieved_at_utc=retrieved_at,
        mappings=tuple(mappings),
        source=source,
        sha256=digest,
    )


def load_spd_receptor_mapping_contract(
    path: str | Path = DEFAULT_CONTRACT_PATH,
) -> SPDReceptorMappingContract:
    """Load the versioned contract and validate all identities and provenance."""

    contract_path = Path(path)
    try:
        payload = json.loads(contract_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ReceptorMappingContractError(
            f"unable to load receptor mapping contract {contract_path}: {exc}"
        ) from exc
    return validate_spd_receptor_mapping_contract(
        payload,
        source=str(contract_path),
        contract_sha256=_file_sha256(contract_path),
    )


def get_spd_receptor_mapping(
    pdb_id: object,
    *,
    contract: SPDReceptorMappingContract | None = None,
) -> SPDReceptorMapping | None:
    """Return a validated mapping for a PDB, or ``None`` outside contract scope."""

    active_contract = contract or load_spd_receptor_mapping_contract()
    normalized = str(pdb_id).strip().upper()
    return active_contract.by_pdb().get(normalized)


def _validate_scoped_input_identities(
    source: pd.DataFrame,
    *,
    pdb_col: str,
    contract: SPDReceptorMappingContract,
) -> None:
    """Reject a strict rebuild when its input no longer matches the contract."""

    normalized_pdb = source[pdb_col].astype("string").str.strip().str.upper()
    for mapping in contract.mappings:
        rows = normalized_pdb.eq(mapping.pdb_id).fillna(False)
        if not bool(rows.any()):
            continue
        for field in ("target_id", "target_gene", "target_uniprot"):
            if field not in source:
                continue
            legacy = str(getattr(mapping, f"legacy_{field}")).strip().upper()
            revised_value = getattr(mapping, f"revised_{field}")
            allowed = {legacy}
            if revised_value is not None:
                allowed.add(str(revised_value).strip().upper())
            observed = (
                source.loc[rows, field]
                .dropna()
                .astype(str)
                .str.strip()
                .str.upper()
            )
            unexpected = sorted(
                value for value in observed.unique().tolist() if value and value not in allowed
            )
            if unexpected:
                raise ReceptorMappingContractError(
                    f"{mapping.pdb_id}: input {field} conflicts with the receptor "
                    f"mapping contract; allowed={sorted(allowed)}, "
                    f"observed={unexpected[:5]}"
                )


def apply_spd_receptor_mapping_contract(
    frame: pd.DataFrame,
    *,
    pdb_col: str = "pdb_id",
    contract: SPDReceptorMappingContract | None = None,
) -> pd.DataFrame:
    """Annotate rows without transferring or clearing target-derived labels.

    Revised identity columns are authoritative only when ``strict_eligible`` is
    true. Callers remain responsible for rebuilding fields identified by
    ``requires_target_derived_rebuild`` before using strict rows.
    """

    if pdb_col not in frame.columns:
        raise KeyError(f"missing PDB column {pdb_col!r}")
    active_contract = contract or load_spd_receptor_mapping_contract()
    mapping_by_pdb = active_contract.by_pdb()
    normalized_pdb = frame[pdb_col].astype("string").str.strip().str.upper()

    output = frame.copy()
    records: Sequence[SPDReceptorMapping | None] = tuple(
        mapping_by_pdb.get(pdb_id) if pd.notna(pdb_id) else None
        for pdb_id in normalized_pdb
    )
    common_values: dict[str, object] = {
        "policy_version": active_contract.policy_version,
        "contract_source": active_contract.source,
        "contract_sha256": active_contract.sha256,
    }
    mapping_fields = (
        "status",
        "strict_eligible",
        "exclusion_reason",
        "changed",
        "requires_target_derived_rebuild",
        "legacy_target_id",
        "legacy_target_gene",
        "legacy_target_uniprot",
        "revised_target_id",
        "revised_target_gene",
        "revised_target_uniprot",
        "selected_chain",
        "site_chain",
    )
    for field in OUTPUT_FIELDS:
        column = f"{OUTPUT_PREFIX}{field}"
        if field in common_values:
            output[column] = common_values[field]
        elif field in mapping_fields:
            unaffected_default: object = pd.NA
            if field == "status":
                unaffected_default = "unaffected"
            elif field == "strict_eligible":
                unaffected_default = True
            elif field in ("changed", "requires_target_derived_rebuild"):
                unaffected_default = False
            elif field == "exclusion_reason":
                unaffected_default = ""
            output[column] = [
                getattr(record, field) if record is not None else unaffected_default
                for record in records
            ]
        else:  # pragma: no cover - OUTPUT_FIELDS and mapping_fields are constants
            raise AssertionError(f"unhandled receptor mapping output field {field}")

    for field in (
        "strict_eligible",
        "changed",
        "requires_target_derived_rebuild",
    ):
        column = f"{OUTPUT_PREFIX}{field}"
        output[column] = output[column].astype("boolean")
    return output


def apply_receptor_target_contract(
    source: pd.DataFrame,
    *,
    mode: str = "legacy",
    contract_path: str | Path = DEFAULT_CONTRACT_PATH,
    pdb_col: str = "pdb_id",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Apply the reconciled identity policy and return the frame plus summary.

    ``legacy`` is an exact schema and identity pass-through. ``strict`` annotates
    all rows and applies the validated 9I52 revision. Unresolved 6HUJ/8HCQ rows
    retain their legacy identity but are marked strict-ineligible. This function
    does not transfer, clear, or rebuild target-derived labels.
    """

    normalized_mode = str(mode).strip().lower()
    if normalized_mode not in {"legacy", "strict"}:
        raise ValueError("mode must be 'legacy' or 'strict'")
    contract = load_spd_receptor_mapping_contract(contract_path)
    if normalized_mode == "strict":
        _validate_scoped_input_identities(
            source,
            pdb_col=pdb_col,
            contract=contract,
        )
    annotated = apply_spd_receptor_mapping_contract(
        source,
        pdb_col=pdb_col,
        contract=contract,
    )
    status_col = f"{OUTPUT_PREFIX}status"
    strict_col = f"{OUTPUT_PREFIX}strict_eligible"
    changed_col = f"{OUTPUT_PREFIX}changed"
    rebuild_col = f"{OUTPUT_PREFIX}requires_target_derived_rebuild"
    in_contract = annotated[status_col].ne("unaffected")

    if normalized_mode == "strict":
        output = annotated
        eligible_revision = in_contract & output[strict_col].fillna(False)
        identity_fields = ("target_id", "target_gene", "target_uniprot")
        for field in identity_fields:
            revised_column = f"{OUTPUT_PREFIX}revised_{field}"
            if field not in output.columns:
                output[field] = pd.NA
            output.loc[eligible_revision, field] = output.loc[
                eligible_revision, revised_column
            ]
    else:
        output = source.copy()

    pdb_counts = (
        annotated.loc[in_contract, pdb_col]
        .astype("string")
        .str.strip()
        .str.upper()
        .value_counts(dropna=False)
        .sort_index()
    )
    summary: dict[str, Any] = {
        "mode": normalized_mode,
        "policy_version": contract.policy_version,
        "contract_source": contract.source,
        "contract_sha256": contract.sha256,
        "rows_total": int(len(annotated)),
        "rows_in_contract": int(in_contract.sum()),
        "rows_unaffected": int((~in_contract).sum()),
        "rows_strict_eligible": int(annotated[strict_col].fillna(False).sum()),
        "rows_strict_ineligible": int(
            (~annotated[strict_col].fillna(False)).sum()
        ),
        "rows_changed": int(annotated[changed_col].fillna(False).sum()),
        "rows_requiring_target_derived_rebuild": int(
            annotated[rebuild_col].fillna(False).sum()
        ),
        "contract_pdb_row_counts": {
            str(pdb_id): int(count) for pdb_id, count in pdb_counts.items()
        },
    }
    return output, summary


__all__ = [
    "DEFAULT_CONTRACT_PATH",
    "OUTPUT_COLUMNS",
    "POLICY_VERSION",
    "EvidenceReference",
    "ReceptorMappingContractError",
    "SPDReceptorMapping",
    "SPDReceptorMappingContract",
    "apply_receptor_target_contract",
    "apply_spd_receptor_mapping_contract",
    "get_spd_receptor_mapping",
    "load_spd_receptor_mapping_contract",
    "validate_spd_receptor_mapping_contract",
]
