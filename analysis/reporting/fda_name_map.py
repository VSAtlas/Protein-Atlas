# -*- coding: utf-8 -*-
import csv
import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

from analysis._common import clean_text
from analysis.reporting.report_config import (
    read_config_keys,
    report_config_candidates,
)
from analysis.reporting.value_utils import strip_quotes
from prep_ligands.fda_mapping_identity import (
    has_resolved_identity,
    resolved_preferred_name,
    suppresses_legacy_aliases,
)

_PRIMARY_NAME_FIELDS = (
    "display_name",
    "pubchem_record_title",
    "drugcentral_generic_name",
    "rxnorm_generic_name",
    "generic_name",
    "pubchem_name",
    "pubchem_iupac_name",
    "remark_name",
    "sdf_title",
)
_FALLBACK_NAME_FIELDS = (
    "fda_name",
    "name",
    "drug_name",
    "preferred_name",
    "international_nonproprietary_name",
    "inn",
    "brand_name",
    "brand_names",
    "trade_name",
    "label_name",
)
_SYNONYM_FIELDS = (
    "synonyms",
    "alias",
    "alts",
    "pubchem_synonyms",
    "brand_names",
    "rxnorm_brand_names",
    "drugcentral_brand_names",
)

_RDK_RE = re.compile(r"rdk[_-]?(\d+)", re.IGNORECASE)
@dataclass
class LibraryIndex:
    id_to_name: Dict[str, str] = field(default_factory=dict)
    name_index: Dict[str, List[str]] = field(default_factory=dict)
    resolved_identity_ids: set[str] = field(default_factory=set)

    def add(
        self,
        rdk_id: str,
        canonical_name: str,
        aliases: List[str],
        *,
        resolved_identity: bool = False,
    ) -> None:
        rid = _extract_rdk_id(rdk_id) or rdk_id
        rid = str(rid).strip()
        if not rid:
            return
        name = str(canonical_name or rid).strip() or rid
        self.id_to_name[rid] = name
        if resolved_identity:
            self.resolved_identity_ids.add(rid)
        tokens = [name, *aliases, rid]
        for token in tokens:
            key = _norm(token)
            if not key:
                continue
            self.name_index.setdefault(key, [])
            if rid not in self.name_index[key]:
                self.name_index[key].append(rid)


def _norm(value: Any) -> str:
    text = str(value or "").strip().lower()
    text = re.sub(r"[-_/\\,;:\|\[\]\(\)\{\}\.\+\*'\"]+", " ", text)
    text = re.sub(r"\s+", " ", text)
    return text


def _tokenize(value: Any) -> List[str]:
    return [tok for tok in re.split(r"[^a-z0-9\+]+", _norm(value)) if tok]


def _extract_rdk_id(value: Any) -> Optional[str]:
    match = _RDK_RE.search(str(value or ""))
    if not match:
        return None
    return f"rdk_{match.group(1).zfill(7)}"


def _clean_cell(row: Dict[str, Any], key: str) -> str:
    return clean_text(row.get(key))


def _pick_first(row: Dict[str, Any], keys: List[str]) -> str:
    for key in keys:
        value = _clean_cell(row, key)
        if value:
            return value
    return ""


def load_library_index(mapping_csv: str) -> LibraryIndex:
    idx = LibraryIndex()
    with open(mapping_csv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            lowered = {str(k).strip().lower(): (v or "") for k, v in row.items()}

            raw_id = _pick_first(
                lowered,
                ["rdk_id", "rdkid", "ligand_id", "id", "path", "sdf_title", "remark_name"],
            )
            rdk_id = _extract_rdk_id(raw_id)

            if not rdk_id and str(lowered.get("scheme", "")).strip().lower() == "rdk":
                file_num = str(lowered.get("file_num", "")).strip()
                if file_num:
                    try:
                        rdk_id = f"rdk_{int(float(file_num)):07d}"
                    except ValueError:
                        rdk_id = None

            if not rdk_id:
                continue

            resolved_name = resolved_preferred_name(lowered)
            preferred_name = resolved_name or _pick_first(
                lowered, list(_PRIMARY_NAME_FIELDS + _FALLBACK_NAME_FIELDS)
            )
            if not preferred_name:
                preferred_name = rdk_id

            synonyms_raw: List[str] = [resolved_name] if resolved_name else []
            if not suppresses_legacy_aliases(lowered):
                for key in _SYNONYM_FIELDS:
                    value = _clean_cell(lowered, key)
                    if value:
                        synonyms_raw.append(value)

            path_val = _clean_cell(lowered, "path")
            if path_val:
                stem = Path(path_val).stem
                if stem:
                    synonyms_raw.append(stem)

            synonyms: List[str] = []
            for raw in synonyms_raw:
                synonyms.extend(
                    token.strip()
                    for token in re.split(r"[|,;]", raw)
                    if token.strip()
                )

            idx.add(
                rdk_id,
                preferred_name,
                synonyms,
                resolved_identity=has_resolved_identity(lowered),
            )
    return idx


def resolve_corresponding_name_for_rdk(
    rdk_id: str, fda_index: LibraryIndex
) -> Optional[str]:
    rid = _extract_rdk_id(rdk_id) or str(rdk_id or "").strip()
    if not rid:
        return None
    return fda_index.id_to_name.get(rid)


def resolve_corresponding_name_from_text(
    text: str, fda_index: LibraryIndex
) -> Optional[str]:
    key = _norm(text)
    exact_ids = fda_index.name_index.get(key, [])
    if len(exact_ids) == 1:
        rid = exact_ids[0]
        return fda_index.id_to_name.get(rid)
    if len(exact_ids) > 1:
        return None

    rdk_id = _extract_rdk_id(text)
    if rdk_id:
        return fda_index.id_to_name.get(rdk_id)

    for token in _tokenize(text):
        token_ids = fda_index.name_index.get(token, [])
        if len(token_ids) == 1:
            rid = token_ids[0]
            name = fda_index.id_to_name.get(rid)
            if name:
                return name
        elif len(token_ids) > 1:
            return None
    return None


def mapping_file_sha256(mapping_csv: Optional[Path]) -> str:
    if mapping_csv is None or not mapping_csv.is_file():
        return ""
    digest = hashlib.sha256()
    with mapping_csv.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def mapping_file_provenance_path(
    mapping_csv: Optional[Path],
    repo_root: Path,
) -> str:
    """Return source-safe mapping provenance; the file hash remains authoritative."""

    if mapping_csv is None:
        return ""
    mapping = Path(mapping_csv).expanduser().resolve()
    root = Path(repo_root).expanduser().resolve()
    try:
        return mapping.relative_to(root).as_posix()
    except ValueError:
        return f"external/{mapping.name}"


def resolve_authoritative_ligand_display_name(
    ligand_base: str,
    ligand_file: str,
    fda_index: Optional[LibraryIndex],
) -> Optional[str]:
    if fda_index is None:
        return None
    exact_ids = {
        f"rdk_{match.group(1).zfill(7)}"
        for text in (ligand_base, ligand_file)
        if (match := _RDK_RE.search(text or "")) is not None
    }
    if len(exact_ids) != 1:
        return None
    rdk_id = next(iter(exact_ids))
    if rdk_id not in fda_index.resolved_identity_ids:
        return None
    return clean_text(resolve_corresponding_name_for_rdk(rdk_id, fda_index)) or None


def _resolve_path(value: str, base_dir: Path) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return path


def resolve_mapping_csv_path(
    repo_root: Path, run_id: str, cli_value: Optional[str] = None
) -> Optional[Path]:
    if cli_value is not None:
        candidate = strip_quotes(str(cli_value))
        if candidate:
            path = _resolve_path(candidate, repo_root)
            if path.exists():
                return path

    config_candidates = report_config_candidates(repo_root, run_id)

    # Priority for implicit defaults:
    # 1) FDA_MAPPING_CSV (current config key)
    # 2) backward-compatible aliases used by older flows
    key_candidates = [
        "FDA_MAPPING_CSV",
        "fda_mapping_csv",
        "fda_mapping_filled.csv",
        "IDENTIFY_MAPPING_CSV",
    ]
    for path in config_candidates:
        values = read_config_keys(path, key_candidates)
        for key_name in key_candidates:
            value = values.get(key_name)
            if not value:
                continue
            resolved = _resolve_path(value, path.parent)
            if resolved.exists():
                return resolved

    fallbacks = [
        repo_root / "fda_mapping_filled.csv",
        repo_root / "chemdb" / "data" / "fda_mapping_from_pdbqt.csv",
        repo_root / "data" / "fda_mapping_from_pdbqt.csv",
    ]
    for path in fallbacks:
        if path.exists():
            return path
    return None


def try_load_fda_index(mapping_csv: Optional[Path]) -> Optional[LibraryIndex]:
    if mapping_csv is None:
        return None
    try:
        return load_library_index(str(mapping_csv))
    except (OSError, UnicodeDecodeError, csv.Error, ValueError):
        return None


def resolve_ligand_display_name(
    ligand_base: str, ligand_file: str, fda_index: Optional[LibraryIndex]
) -> str:
    base = clean_text(ligand_base)
    lig_file = clean_text(ligand_file)
    if fda_index is None:
        return base

    exact_ids = {
        f"rdk_{match.group(1).zfill(7)}"
        for text in (base, lig_file)
        if (match := _RDK_RE.search(text or "")) is not None
    }
    if len(exact_ids) > 1:
        return base
    if exact_ids:
        rdk_id = next(iter(exact_ids))
        name = clean_text(resolve_corresponding_name_for_rdk(rdk_id, fda_index))
        if name:
            return name

    for text in (base, Path(lig_file).stem if lig_file else ""):
        name = clean_text(resolve_corresponding_name_from_text(text, fda_index))
        if name:
            return name

    return base
