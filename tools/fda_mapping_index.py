"""Shared FDA mapping CSV index for benchmark tooling."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

try:
    from chemdb.tools.fill_fda_mapping_names import display_name_quality_reasons
except Exception:

    def display_name_quality_reasons(
        name: str,
        pubchem_iupac_name: str | None = None,
        *,
        strict: bool = False,
    ) -> list[str]:
        _ = pubchem_iupac_name, strict
        return ["empty"] if not str(name or "").strip() else []


def _clean(value: object) -> str:
    text = str(value or "").strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def _norm(value: object) -> str:
    text = _clean(value).lower()
    text = re.sub(r"[-_/\\,;:\|\[\]\(\)\{\}\.\+\*'\"]+", " ", text)
    text = re.sub(r"\s+", " ", text)
    return text


def _is_stereo_prefix_extension(candidate_key: str, display_key: str) -> bool:
    if not display_key or not candidate_key.endswith(display_key):
        return False
    prefix = candidate_key[: -len(display_key)]
    if not prefix:
        return False
    return re.fullmatch(r"(?:\d+)?[rsdl]+", prefix) is not None


def _tokenize(value: object) -> list[str]:
    return [tok for tok in re.split(r"[^a-z0-9\+]+", _norm(value)) if tok]


def _split_values(value: object) -> list[str]:
    return [part.strip() for part in re.split(r"[|,;]", _clean(value)) if part.strip()]


_RDK_RE = re.compile(r"rdk[_-]?(\d+)", re.IGNORECASE)


def extract_rdk_id(value: object) -> Optional[str]:
    match = _RDK_RE.search(str(value or ""))
    if not match:
        return None
    return f"rdk_{match.group(1).zfill(7)}"


@dataclass
class MappingRow:
    path: str
    display_name: str = ""
    generic_name: str = ""
    brand_names: str = ""
    pubchem_name: str = ""
    pubchem_record_title: str = ""
    pubchem_iupac_name: str = ""
    pubchem_synonyms: str = ""
    rxnorm_generic_name: str = ""
    rxnorm_brand_names: str = ""
    drugcentral_generic_name: str = ""
    drugcentral_brand_names: str = ""
    remark_name: str = ""
    sdf_title: str = ""
    inchikey: str = ""

    @classmethod
    def from_csv_row(cls, row: dict[str, object]) -> "MappingRow":
        return cls(
            path=_clean(row.get("path")),
            display_name=_clean(row.get("display_name")),
            generic_name=_clean(row.get("generic_name")),
            brand_names=_clean(row.get("brand_names")),
            pubchem_name=_clean(row.get("pubchem_name")),
            pubchem_record_title=_clean(row.get("pubchem_record_title")),
            pubchem_iupac_name=_clean(row.get("pubchem_iupac_name")),
            pubchem_synonyms=_clean(row.get("pubchem_synonyms")),
            rxnorm_generic_name=_clean(row.get("rxnorm_generic_name")),
            rxnorm_brand_names=_clean(row.get("rxnorm_brand_names")),
            drugcentral_generic_name=_clean(row.get("drugcentral_generic_name")),
            drugcentral_brand_names=_clean(row.get("drugcentral_brand_names")),
            remark_name=_clean(row.get("remark_name")),
            sdf_title=_clean(row.get("sdf_title")),
            inchikey=_clean(row.get("inchikey")),
        )

    def all_name_fields(self) -> list[tuple[str, str]]:
        return [
            ("rxnorm_generic_name", self.rxnorm_generic_name),
            ("drugcentral_generic_name", self.drugcentral_generic_name),
            ("pubchem_record_title", self.pubchem_record_title),
            ("generic_name", self.generic_name),
            ("display_name", self.display_name),
            ("brand_names", self.brand_names),
            ("pubchem_name", self.pubchem_name),
            ("pubchem_synonyms", self.pubchem_synonyms),
            ("rxnorm_brand_names", self.rxnorm_brand_names),
            ("drugcentral_brand_names", self.drugcentral_brand_names),
            ("remark_name", self.remark_name),
            ("sdf_title", self.sdf_title),
        ]

    def preferred_name(self) -> str:
        display = _clean(self.display_name)
        if display and not self.display_quality_reasons() and not self._semantic_reasons():
            return display
        for candidate in self._better_name_candidates():
            if candidate and not display_name_quality_reasons(
                candidate, self.pubchem_iupac_name, strict=True
            ):
                return candidate
        return display

    def _better_name_candidates(self) -> list[str]:
        candidates: list[str] = []
        for field in [
            self.pubchem_record_title,
            self.generic_name,
            self.pubchem_name,
            self.drugcentral_generic_name,
            self.rxnorm_generic_name,
            self.remark_name,
            self.sdf_title,
        ]:
            candidate = _clean(field)
            if candidate:
                candidates.append(candidate)
        for field in [
            self.brand_names,
            self.drugcentral_brand_names,
            self.rxnorm_brand_names,
        ]:
            candidates.extend(_split_values(field))
        return list(dict.fromkeys(candidate for candidate in candidates if candidate))

    def display_quality_reasons(self) -> list[str]:
        reasons = display_name_quality_reasons(
            self.display_name,
            self.pubchem_iupac_name,
            strict=True,
        )
        display_key = _norm(self.display_name).replace(" ", "")
        trusted_keys = {
            _norm(self.pubchem_record_title).replace(" ", ""),
            _norm(self.drugcentral_generic_name).replace(" ", ""),
            _norm(self.rxnorm_generic_name).replace(" ", ""),
        }
        if display_key and display_key in {key for key in trusted_keys if key}:
            return [reason for reason in reasons if reason != "iupac_like"]
        return reasons

    def audit_reasons(self) -> list[str]:
        reasons = list(self.display_quality_reasons())
        reasons.extend(self._semantic_reasons())
        return list(dict.fromkeys(reason for reason in reasons if reason))

    def _semantic_reasons(self) -> list[str]:
        reasons: list[str] = []
        display = _clean(self.display_name)
        if not display:
            return reasons
        display_key = _norm(display).replace(" ", "")
        if not display_key:
            return reasons
        trusted_keys = {
            _norm(self.pubchem_record_title).replace(" ", ""),
            _norm(self.drugcentral_generic_name).replace(" ", ""),
            _norm(self.rxnorm_generic_name).replace(" ", ""),
        }
        if display_key in {key for key in trusted_keys if key}:
            return reasons
        candidates = [
            candidate
            for candidate in self._better_name_candidates()
            if not display_name_quality_reasons(candidate, self.pubchem_iupac_name, strict=True)
        ]
        candidate_keys = [_norm(candidate).replace(" ", "") for candidate in candidates]
        if len(display_key) <= 4 and any(len(key) > len(display_key) for key in candidate_keys):
            reasons.append("short_alias_display_name")
        if any(
            display_key != key
            and len(display_key) >= 5
            and key.endswith(display_key)
            and len(key) - len(display_key) <= 3
            and not _is_stereo_prefix_extension(key, display_key)
            for key in candidate_keys
        ):
            reasons.append("truncated_display_name")
        return reasons


class MappingIndex:
    def __init__(self, csv_path: Optional[Path]):
        self.csv_path = Path(csv_path) if csv_path else None
        self.rows: list[MappingRow] = []
        if not self.csv_path or not self.csv_path.is_file():
            return
        with self.csv_path.open("r", encoding="utf-8", errors="ignore", newline="") as f:
            reader = csv.DictReader(f)
            for raw in reader:
                row = MappingRow.from_csv_row(raw)
                if row.path:
                    self.rows.append(row)

    def rows_by_rdk_id(self, rdk_id: str) -> list[MappingRow]:
        rid = (extract_rdk_id(rdk_id) or _clean(rdk_id)).lower()
        if not rid:
            return []
        return [row for row in self.rows if rid in Path(row.path).stem.lower()]

    def resolve_name_for_rdk(self, rdk_id: str) -> Optional[str]:
        for row in self.rows_by_rdk_id(rdk_id):
            name = row.preferred_name()
            if name:
                return name
        return None

    def audit_reasons_for_rdk(self, rdk_id: str) -> list[str]:
        reasons: list[str] = []
        for row in self.rows_by_rdk_id(rdk_id):
            reasons.extend(row.audit_reasons())
        return list(dict.fromkeys(reasons))

    def search(
        self,
        hints: Sequence[str],
        inchikey: Optional[str] = None,
        max_results: int = 6,
    ) -> list[tuple[MappingRow, int, str]]:
        hints_norm = [_norm(h) for h in hints if h]
        hint_tokens = {token for hint in hints_norm for token in _tokenize(hint)}
        raw_het_codes = {h.strip().upper() for h in hints if h and 2 <= len(h.strip()) <= 5}
        out: list[tuple[MappingRow, int, str]] = []

        for row in self.rows:
            best = 0
            why = ""
            path = Path(row.path)
            stem_upper = path.stem.upper()
            base_upper = path.name.upper()
            parents_upper = " ".join(parent.name.upper() for parent in path.parents)

            if inchikey and row.inchikey and row.inchikey.strip().upper() == inchikey.strip().upper():
                best, why = 100, "inchikey_exact"

            if best < 100 and raw_het_codes:
                for code in raw_het_codes:
                    if re.search(rf"\b{re.escape(code)}\b", stem_upper) or re.search(
                        rf"\b{re.escape(code)}\b", parents_upper
                    ):
                        score = 93
                        reason = f"path_token:{code}"
                    elif code in base_upper:
                        score = 88
                        reason = f"path_substr:{code}"
                    else:
                        score = 0
                        reason = ""
                    if score > best:
                        best, why = score, reason

            if best < 100 and hints_norm:
                for field, value in row.all_name_fields():
                    value_norm = _norm(value)
                    if not value_norm:
                        continue
                    if value_norm in hints_norm:
                        score = 95
                        reason = f"{field}_exact"
                    elif any(hint in value_norm for hint in hints_norm if len(hint) >= 3):
                        score = 85
                        reason = f"{field}_substr"
                    else:
                        tokens = set(_tokenize(value_norm))
                        overlap = len(tokens & hint_tokens)
                        score = 70 if overlap >= 2 else (65 if overlap == 1 else 0)
                        reason = f"{field}_tokens:{overlap}"
                    if score > best:
                        best, why = score, reason

            if best > 0:
                out.append((row, best, why))

        out.sort(key=lambda item: item[1], reverse=True)
        return out[:max_results]


def iter_audited_rows(rows: Iterable[MappingRow]) -> list[tuple[MappingRow, list[str]]]:
    audited: list[tuple[MappingRow, list[str]]] = []
    for row in rows:
        reasons = row.audit_reasons()
        if reasons:
            audited.append((row, reasons))
    return audited
