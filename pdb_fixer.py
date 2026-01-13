import os
import re
import logging
import yaml
from pathlib import Path
from typing import (
    Any,
    Dict,
    Iterable,
    List,
    Mapping,
    NamedTuple,
    Optional,
    Set,
    Tuple,
    Union,
)
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select
from types import SimpleNamespace
from installation import load_config

logger = logging.getLogger(__name__)

# Load config
config = load_config()

BASE_DIR = Path(__file__).resolve().parent

ALIASES_PATH = (
    config.get("ALIASES_PATH")
    or os.environ.get("ALIASES_YAML")
    or str(BASE_DIR / "aliases.yaml")
)

# if not found, try canonical chemdb/aliases.yaml automatically
if not os.path.exists(ALIASES_PATH):
    for candidate in (
        BASE_DIR / "chemdb" / "aliases.yaml",
        BASE_DIR / "activesite" / "aliases.yaml",
        BASE_DIR / "activesite" / "chemdb" / "aliases.yaml",
        BASE_DIR.parent / "chemdb" / "aliases.yaml",
    ):
        try_path = candidate.resolve() if hasattr(candidate, "resolve") else candidate
        if Path(try_path).exists():
            ALIASES_PATH = str(try_path)
            break

_aliases_cache = None
_rules_cache = None


class AliasSets(NamedTuple):
    waters: Set[str]
    cofactors: Set[str]
    element_tokens: Set[str]
    element_alias: Dict[str, str]
    elem_tokens_canonical: Set[str]


_ALIAS_POLICY_LOGGED = False
_ALIAS_WARNED_KEYS: set[str] = set()


_ION_BREADCRUMB_METAL_ORDER = (
    "ZN",
    "HG",
    "MG",
    "FE",
    "MN",
    "CO",
    "NI",
    "CU",
    "CD",
    "CA",
)
_ION_BREADCRUMB_SIMPLE_ORDER = ("NA", "K", "CL", "BR", "I")
_ION_BREADCRUMB_WATERS = {"HOH", "WAT"}
_ION_BREADCRUMB_ALIAS_MAP = {
    "ZN1": "ZN",
    "ZN2": "ZN",
    "ZN3": "ZN",
    "ZN+": "ZN",
    "ZN+2": "ZN",
    "ZN2+": "ZN",
    "MG1": "MG",
    "MG2": "MG",
    "MG+": "MG",
    "MN2": "MN",
    "MN3": "MN",
    "FE2": "FE",
    "FE3": "FE",
    "CO2": "CO",
    "NI2": "NI",
    "CU1": "CU",
    "CU2": "CU",
    "CD2": "CD",
    "HG2": "HG",
    "CA1": "CA",
    "CA2": "CA",
    "NA1": "NA",
    "K1": "K",
    "K+": "K",
    "CL-": "CL",
    "BR-": "BR",
    "I-": "I",
}
_ION_BREADCRUMB_ENABLED_VALUES = {"1", "true", "yes"}
_ION_BREADCRUMB_DISABLED_VALUES = {"0", "false", "no"}


def _ion_audit_env_enabled() -> bool:
    raw = os.environ.get("ION_AUDIT")
    if raw is None:
        return True
    text = raw.strip()
    if not text:
        return True
    lowered = text.lower()
    if lowered in _ION_BREADCRUMB_DISABLED_VALUES:
        return False
    return lowered in _ION_BREADCRUMB_ENABLED_VALUES


def _short_path_for_log(path: Union[str, Path]) -> str:
    try:
        return str(Path(path).resolve(strict=False).relative_to(Path.cwd()))
    except Exception:
        try:
            return str(Path(path).resolve(strict=False))
        except Exception:
            return str(Path(path))


def _format_breadcrumb_counts(counts: dict[str, int], order: Iterable[str]) -> str:
    parts: list[str] = []
    for token in order:
        parts.append(f"{token}:{int(counts.get(token, 0))}")
    extras = [tok for tok in sorted(counts) if tok not in order]
    for token in extras:
        parts.append(f"{token}:{int(counts.get(token, 0))}")
    return "{" + ",".join(parts) + "}"


def summarize_ions(pdb_path: Union[str, Path]) -> dict[str, object]:
    path = Path(pdb_path)
    metals_counts = {token: 0 for token in _ION_BREADCRUMB_METAL_ORDER}
    simple_counts = {token: 0 for token in _ION_BREADCRUMB_SIMPLE_ORDER}
    if not path.exists():
        return {
            "metals": metals_counts,
            "simple_ions": simple_counts,
            "waters": 0,
            "other_het": 0,
            "missing": True,
        }

    metal_hits = {token: set() for token in _ION_BREADCRUMB_METAL_ORDER}
    simple_hits = {token: set() for token in _ION_BREADCRUMB_SIMPLE_ORDER}
    water_hits: set[str] = set()
    other_hits: set[tuple[str, str]] = set()

    try:
        with path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.startswith("HETATM"):
                    continue
                resname_raw = line[17:20].strip().upper()
                if not resname_raw:
                    continue
                chain = (line[21:22] or "-").strip() or "-"
                resseq = (line[22:26] or "0").strip() or "0"
                icode = (line[26:27] or "").strip() or "-"
                canonical = _ION_BREADCRUMB_ALIAS_MAP.get(resname_raw, resname_raw)
                token = canonical.upper()
                residue_key = f"{chain}:{resseq}:{icode}"
                if token in _ION_BREADCRUMB_METAL_ORDER:
                    metal_hits.setdefault(token, set()).add(residue_key)
                    continue
                if token in _ION_BREADCRUMB_SIMPLE_ORDER:
                    simple_hits.setdefault(token, set()).add(residue_key)
                    continue
                if token in _ION_BREADCRUMB_WATERS:
                    water_hits.add(residue_key)
                    continue
                other_hits.add((token, residue_key))
    except FileNotFoundError:
        return {
            "metals": metals_counts,
            "simple_ions": simple_counts,
            "waters": 0,
            "other_het": 0,
            "missing": True,
        }
    except Exception as exc:
        logging.debug(
            "[ions.breadcrumb] stage=summarize error=%s file=%s",
            exc,
            path,
        )

    for token in metals_counts:
        metals_counts[token] = len(metal_hits.get(token, set()))
    for token in simple_counts:
        simple_counts[token] = len(simple_hits.get(token, set()))

    return {
        "metals": metals_counts,
        "simple_ions": simple_counts,
        "waters": len(water_hits),
        "other_het": len(other_hits),
        "missing": False,
    }


def log_pre_variant_policy_breadcrumb(pdb_path: Union[str, Path]) -> None:
    if not _ion_audit_env_enabled():
        return
    summary = summarize_ions(pdb_path)
    short_path = _short_path_for_log(pdb_path)
    if summary.get("missing"):
        logging.info(
            "[ions.breadcrumb] stage=pre_variant_policy file=%s missing=true",
            short_path,
        )
        return
    logging.info(
        "[ions.breadcrumb] stage=pre_variant_policy file=%s metals=%s waters=%d simple_ions=%s other_het=%d",
        short_path,
        _format_breadcrumb_counts(
            summary.get("metals", {}), _ION_BREADCRUMB_METAL_ORDER
        ),
        int(summary.get("waters", 0) or 0),
        _format_breadcrumb_counts(
            summary.get("simple_ions", {}), _ION_BREADCRUMB_SIMPLE_ORDER
        ),
        int(summary.get("other_het", 0) or 0),
    )


def _log_alias_tokens(key: str, tokens: set[str]) -> None:
    sample = ",".join(sorted(tokens)[:10]) if tokens else "none"
    logging.info(
        "[aliases.tokens] key=%s count=%d sample=%s source=%s",
        key,
        len(tokens),
        sample,
        ALIASES_PATH,
    )


def _normalize_alias_token(
    token: str,
    alias_map: Dict[str, str],
    canonical_targets: Optional[Set[str]] = None,
) -> str:
    raw = (token or "").strip().upper()
    if not raw:
        return ""
    mapped = alias_map.get(raw)
    if mapped:
        return mapped
    if canonical_targets is not None and raw not in canonical_targets:
        if raw not in _ALIAS_WARNED_KEYS and any(
            ch.isdigit() or ch in "+-" for ch in raw
        ):
            logging.debug(
                "[aliases.alias.warn] key=%s had no mapping in retain_element_alias_map",
                raw,
            )
            _ALIAS_WARNED_KEYS.add(raw)
    return raw


def _format_alias_sample(tokens: Iterable[str], limit: int = 6) -> str:
    deduped: list[str] = []
    seen: set[str] = set()
    for tok in sorted(str(t).strip().upper() for t in tokens if str(t).strip()):
        if tok in seen:
            continue
        deduped.append(tok)
        seen.add(tok)
    if not deduped:
        return "none"
    if len(deduped) > limit:
        trimmed = deduped[:limit]
        trimmed.append(f"+{len(deduped) - limit}")
        return ",".join(trimmed)
    return ",".join(deduped)


def _derive_alias_sets(
    aliases_root: dict,
    *,
    as_set,
    element_symbol_hints: Set[str],
) -> tuple[AliasSets, set[str], bool]:
    alias_cfg = aliases_root or {}

    alias_raw = alias_cfg.get("retain_element_alias_map", {})
    element_alias: Dict[str, str] = {}
    if isinstance(alias_raw, dict):
        for k, v in alias_raw.items():
            key = str(k or "").strip().upper()
            val = str(v or "").strip().upper()
            if not key or not val:
                continue
            element_alias[key] = val
    if not element_alias:
        element_alias = {
            k.upper(): v.upper() for k, v in _ION_BREADCRUMB_ALIAS_MAP.items()
        }

    waters = as_set(
        alias_cfg.get("retain_water_resnames"), section="retain_water_resnames"
    )
    cofactors = as_set(
        alias_cfg.get("retain_cofactor_resnames"), section="retain_cofactor_resnames"
    )
    element_tokens = as_set(
        alias_cfg.get("retain_element_tokens"), section="retain_element_tokens"
    )

    # TODO(aliases-migration): uses legacy retain_in_receptor_resnames.
    legacy = as_set(alias_cfg.get("retain_in_receptor_resnames", []))
    used_backcompat = False
    if not (waters or cofactors or element_tokens):
        if legacy:
            used_backcompat = True
            waters = {tok for tok in legacy if tok in _ION_BREADCRUMB_WATERS}
            canonical_elements: set[str] = set()
            for tok in legacy:
                canonical = _normalize_alias_token(tok, element_alias)
                if not canonical:
                    continue
                if (len(canonical) <= 2) or (canonical in element_symbol_hints):
                    canonical_elements.add(canonical)
            element_tokens = set(canonical_elements)
            cofactors = {
                tok
                for tok in legacy
                if tok not in waters
                and _normalize_alias_token(tok, element_alias) not in canonical_elements
            }
        else:
            waters = set()
            cofactors = set()
            element_tokens = set()

    alias_values = {
        str(v or "").strip().upper()
        for v in element_alias.values()
        if str(v or "").strip()
    }
    canonical_hints = set(element_tokens) | alias_values
    canonical_elem_tokens: Set[str] = set()
    for tok in set(element_tokens) | alias_values:
        canonical = _normalize_alias_token(
            tok,
            element_alias,
            canonical_targets=canonical_hints,
        )
        if canonical:
            canonical_elem_tokens.add(canonical)

    alias_sets = AliasSets(
        waters=set(waters),
        cofactors=set(cofactors),
        element_tokens=set(element_tokens),
        element_alias=dict(element_alias),
        elem_tokens_canonical=set(canonical_elem_tokens),
    )
    return alias_sets, legacy, used_backcompat


# ---  YAML loader with encoding fallbacks & punctuation cleanup ---
def _load_aliases_yaml():
    import yaml

    global _aliases_cache
    if _aliases_cache is not None:
        return _aliases_cache

    path = ALIASES_PATH
    try:
        # try UTF-8 first (fast path)
        with open(path, "r", encoding="utf-8") as fh:
            _aliases_cache = yaml.safe_load(fh) or {}
            return _aliases_cache
    except Exception as e_utf8:
        logging.warning("[aliases] UTF-8 load failed for %s: %s", path, e_utf8)

    try:
        # fallback: cp1252 decode + punctuation normalization → UTF-8
        raw = Path(path).read_bytes()
        text = raw.decode("cp1252", errors="strict")

        # normalize common Windows punctuation to ASCII
        repl = {
            "\u2018": "'",
            "\u2019": "'",
            "\u201c": '"',
            "\u201d": '"',
            "\u2013": "-",
            "\u2014": "-",
            "\u2026": "...",
            "\u00a0": " ",
            "\u200b": "",
            "\ufeff": "",
        }
        for k, v in repl.items():
            text = text.replace(k, v)

        # fix common mojibake 'â†’' if present
        text = text.replace("â†’", "->")

        text = text.replace("\r\n", "\n").replace("\r", "\n")
        _aliases_cache = yaml.safe_load(text) or {}
        logging.info("[aliases] loaded with cp1252 fallback and sanitized punctuation")
        return _aliases_cache
    except Exception as e_cp:
        logging.warning("[aliases] cp1252 fallback failed for %s: %s", path, e_cp)
        _aliases_cache = {}
        return _aliases_cache


def get_atom_rules():
    """
    Return a SimpleNamespace of normalized sets/maps used by element fixing,
    plus compatibility views so legacy code can still do RULES["element_sets"].
    """
    global _rules_cache
    if _rules_cache is not None:
        return _rules_cache

    a = _load_aliases_yaml() or {}
    logging.info(
        "[aliases.load] source=%s keys=%d",
        ALIASES_PATH,
        len(a),
    )
    es = a.get("element_sets", {}) or {}
    ligand_sets = a.get("ligand_sets", {}) or {}
    meeko_cfg = a.get("meeko", {}) or {}

    token_splitter = re.compile(r"[,\s]+")

    def _flatten(items):
        for x in items or []:
            if isinstance(x, (list, tuple, set)):
                yield from _flatten(x)
            else:
                text = "" if x is None else str(x)
                if not text:
                    continue
                for semi in text.split(";"):
                    segment = semi.strip()
                    if not segment:
                        continue
                    for token in token_splitter.split(segment):
                        tok = token.strip()
                        if tok:
                            yield tok

    def _as_set(items, up=True, section=None):
        tokens = list(_flatten(items))
        normalized = []
        ignored = []
        for tok in tokens:
            norm = tok.strip()
            if not norm:
                ignored.append(tok)
                continue
            normalized.append(norm.upper() if up else norm)
        result = set(normalized)
        if section and ignored:
            filtered = sorted({t.strip() for t in ignored if t and t.strip()})[:5]
            if filtered:
                logging.debug(
                    "[aliases.tokens.ignored] section=%s tokens=%s",
                    section,
                    ",".join(filtered),
                )
        return result

    # element/name logic
    peptide_like = _as_set(es.get("peptide_like_names"))
    one_letter = _as_set(es.get("one_letter_elements"))
    two_letter = _as_set(es.get("two_letter_elements"))
    halide_resnames = _as_set(es.get("halide_resnames"))
    default_element = (es.get("default_element") or "C").upper()
    treat_backbone_ca = bool(es.get("treat_backbone_CA_as_C", True))

    # AD4 types (support both top-level and nested)
    ad4_types_yaml = a.get("ad4_types") or (a.get("pdbqt_types") or {}).get("ad4_types")
    ad4_types = _as_set(ad4_types_yaml) if ad4_types_yaml else set()

    #  ligand & Meeko lists from YAML
    nucleotide_like_resnames = _as_set(ligand_sets.get("nucleotide_like_resnames"))
    meeko_drop_free_ions = _as_set(
        meeko_cfg.get("drop_free_ions"), section="meeko.drop_free_ions"
    )
    strip_raw = a.get("strip_in_receptor_resnames", []) or []
    strip_tokens = _as_set(strip_raw, section="strip_in_receptor_resnames")
    strip_tokens |= meeko_drop_free_ions
    _log_alias_tokens("strip_in_receptor_resnames", strip_tokens)

    alias_sets, legacy_tokens, used_backcompat = _derive_alias_sets(
        a,
        as_set=_as_set,
        element_symbol_hints=one_letter | two_letter,
    )

    _log_alias_tokens("retain_water_resnames", set(alias_sets.waters))
    _log_alias_tokens("retain_cofactor_resnames", set(alias_sets.cofactors))
    _log_alias_tokens("retain_element_tokens", set(alias_sets.element_tokens))

    allow_raw = a.get("allow_in_receptor_resnames", []) or []
    allow_tokens = _as_set(allow_raw, section="allow_in_receptor_resnames")
    _log_alias_tokens("allow_in_receptor_resnames", allow_tokens)

    waters_set = set(alias_sets.waters)
    cofactors_set = set(alias_sets.cofactors)
    element_tokens_raw = set(alias_sets.element_tokens)
    canonical_elements = set(alias_sets.elem_tokens_canonical)
    element_alias = dict(alias_sets.element_alias)

    logging.info(
        "[aliases.loaded] waters=%d cofactors=%d elements_raw=%d aliases=%d",
        len(waters_set),
        len(cofactors_set),
        len(element_tokens_raw),
        len(element_alias),
    )
    logging.info(
        "[aliases.canonical] elem_tokens=%d samples=[%s]",
        len(canonical_elements),
        _format_alias_sample(canonical_elements),
    )
    if used_backcompat:
        logging.info(
            "[aliases.backcompat] using retain_in_receptor_resnames; split not provided"
        )

    required_metals = ["ZN", "MG", "CA", "FE", "MN", "CU", "CO", "NI", "NA", "K"]
    coverage = {tok: (tok in canonical_elements) for tok in required_metals}
    logging.info("[aliases.audit] metal_core_coverage=%s", coverage)

    mode_raw = (
        (os.environ.get("APO_HOLO_MODE") or a.get("APO_HOLO_MODE", "") or "")
        .strip()
        .upper()
    )
    if mode_raw == "APO":
        policy_mode = "APO"
    elif mode_raw == "HOLO":
        policy_mode = "HOLO"
    else:
        policy_mode = "LEGACY"

    if policy_mode == "APO":
        cofactors_policy = set()
    else:
        cofactors_policy = set(cofactors_set)

    retain_res = set(waters_set) | set(canonical_elements) | set(cofactors_policy)

    logging.info(
        "[aliases.policy] mode=%s keep_sets=waters{n=%d} cofactors{n=%d} elements{n=%d} final_retained=%d",
        policy_mode,
        len(waters_set),
        len(cofactors_policy),
        len(canonical_elements),
        len(retain_res),
    )
    logging.info(
        "[aliases.samples] waters=%s",
        _format_alias_sample(waters_set),
    )
    logging.info(
        "[aliases.samples] cofactors=%s",
        _format_alias_sample(cofactors_policy),
    )
    logging.info(
        "[aliases.samples] elements=%s",
        _format_alias_sample(canonical_elements),
    )
    global _ALIAS_POLICY_LOGGED
    if not _ALIAS_POLICY_LOGGED:
        logging.info(
            "[aliases.summary] mode=%s keep={waters:%d, cofactors:%d, elements:%d}",
            policy_mode,
            len(waters_set),
            len(cofactors_policy),
            len(canonical_elements),
        )
        _ALIAS_POLICY_LOGGED = True

    metal_tokens = sorted(tok for tok in canonical_elements if len(tok) <= 2)
    logging.info(
        "[aliases.section] name=element_resnames.metals size=%d sample=[%s]",
        len(metal_tokens),
        _format_alias_sample(metal_tokens),
    )

    compat_element_sets = dict(es)
    compat_retain_list = sorted(retain_res)

    def _normalize_resname_for_rules(token: str) -> str:
        return _normalize_alias_token(
            token,
            element_alias,
            canonical_targets=canonical_elements,
        )

    _rules_cache = SimpleNamespace(
        # normalized sets / maps
        peptide_like=peptide_like,
        one_letter=one_letter,
        two_letter=two_letter,
        halide_resnames=halide_resnames,
        default_element=default_element,
        treat_backbone_ca=treat_backbone_ca,
        retain_resnames=set(retain_res),
        # alias policy exposure
        alias_sets=alias_sets,
        waters=waters_set,
        cofactors=set(cofactors_policy),
        cofactors_all=cofactors_set,
        element_tokens=element_tokens_raw,
        element_alias=element_alias,
        elem_tokens_canonical=canonical_elements,
        policy_mode=policy_mode,
        normalize_resname=_normalize_resname_for_rules,
        # NEW exports used elsewhere
        nucleotide_like_resnames=sorted(nucleotide_like_resnames),
        meeko_drop_free_ions=sorted(meeko_drop_free_ions),
        # name/alias maps (uppercased keys/values)
        prefix_map={
            (k or "").upper(): (v or "").upper()
            for k, v in (es.get("derive_prefix_map") or {}).items()
        },
        special_names={
            (k or "").upper(): (v or "").upper()
            for k, v in (es.get("special_atom_names") or {}).items()
        },
        halide_aliases={
            (k or "").upper(): (v or "").upper()
            for k, v in (es.get("halide_resname_aliases") or {}).items()
        },
        cation_aliases={
            (k or "").upper(): (v or "").upper()
            for k, v in (es.get("cation_resname_aliases") or {}).items()
        },
        # compatibility views for older call sites
        element_sets=compat_element_sets,
        # TODO(aliases-migration): uses legacy retain_in_receptor_resnames.
        retain_in_receptor_resnames=compat_retain_list,
        ad4_types=ad4_types,
        legacy_retain_tokens=legacy_tokens,
    )

    metals_probe = [
        "ZN",
        "HG",
        "MG",
        "FE",
        "MN",
        "CA",
        "CU",
        "CO",
        "NI",
        "NA",
        "K",
        "CL",
    ]
    includes = {tok: (tok in retain_res) for tok in metals_probe}
    detected_metals = sorted([tok for tok in retain_res if tok in metals_probe])
    logging.info(
        "[aliases.audit] retain_in_receptor_resnames size=%d includes=%s",
        len(retain_res),
        includes,
    )
    logging.info(
        "[aliases.audit] metal_tokens_detected=%s",
        ",".join(detected_metals) if detected_metals else "none",
    )
    variant_env = (os.environ.get("APO_HOLO_VARIANT") or "").strip().upper() or "LEGACY"
    logging.info(
        "[activesite.retention] variant=%s retain_resnames_size=%d contains=%s",
        variant_env,
        len(retain_res),
        len(includes),
    )
    return _rules_cache


def scan_helium_counts(text_or_lines) -> int:
    """
    Count occurrences that look like element/ADT 'He' in PDB or PDBQT.
    - PDB: fixed-width columns 77–78
    - PDBQT: trailing ADT token
    """
    n = 0
    if isinstance(text_or_lines, str):
        lines = text_or_lines.splitlines()
    else:
        lines = list(text_or_lines)
    for ln in lines:
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        # PDB fixed-width first
        if len(ln) >= 78 and ln[76:78].strip().upper() == "HE":
            n += 1
            continue
        # PDBQT trailing token
        parts = ln.split()
        if parts and parts[-1].strip().upper() == "HE":
            n += 1
    return n


def rules_version() -> str:
    """
    Return a short identifier for the current YAML/alias ruleset if available,
    else a static placeholder. Keep this stable for log grepping.
    """
    try:
        r = get_atom_rules()
        for k in ("version", "rules_version", "hash", "source_hash"):
            v = (
                getattr(r, k, None)
                if hasattr(r, k)
                else r.get(k)
                if isinstance(r, dict)
                else None
            )
            if v:
                return str(v)
    except Exception:
        pass
    return "rules:v1"


def format_pdb_atom_debug(line: str, line_no: int | None = None) -> str:
    """
    Render a fixed-width debug string exposing key PDB fields + element cols 77–78.
    """
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return line.strip()
    resn = line[17:20].strip()
    chain = line[21].strip() or "_"
    resi = line[22:26].strip() or "0"
    icode = line[26].strip() or " "
    name = line[12:16].strip()
    alt = line[16] if len(line) > 16 else " "
    elem = (line[76:78] if len(line) >= 78 else "  ").strip() or "?"
    ln = f"line{line_no}" if line_no is not None else "line?"
    return f"{ln} {chain}:{resi}{icode}:{resn} {name} alt={alt} elem={elem}  cols77-78='{(line[76:78] if len(line) >= 78 else '  ')}'"


def scan_helium_counts_with_hits(text: str, max_hits: int = 5) -> tuple[int, list[str]]:
    """
    Wrapper around scan_helium_counts that also returns up to `max_hits` formatted offenders.
    """
    n = scan_helium_counts(text)
    if n == 0:
        return 0, []
    hits = []
    for i, ln in enumerate(text.splitlines(), start=1):
        if ln.startswith(("ATOM", "HETATM")) and (ln[76:78].strip().upper() == "HE"):
            hits.append(format_pdb_atom_debug(ln, i))
            if len(hits) >= max_hits:
                break
    return n, hits


def assert_no_helium_in_hydrogen_names(pdb_text: str) -> Tuple[str, int]:
    """
    If an ATOM/HETATM line has an atom *name* beginning with 'H' but element/ADT is 'He',
    rewrite the element to 'H'. Returns (new_text, n_fixes).
    """
    out: List[str] = []
    fixes = 0
    for ln in pdb_text.splitlines(True):
        if ln.startswith(("ATOM  ", "HETATM")) and len(ln) >= 78:
            aname = ln[12:16].strip().upper()
            # ADT usually in the last whitespace token for your PDBQT; for PDB use columns 77-78
            adt_or_elem = ln.split()[-1].strip().upper()
            if aname.startswith("H") and adt_or_elem == "HE":
                # fix element columns if fixed-width PDB, else re-map the trailing ADT token
                # Prefer fixed-width column 77-78 if present:
                if ln[76:78].strip().upper() in {"HE", "H", ""}:
                    ln = ln[:76] + f"{ 'H':>2}" + ln[78:]
                else:
                    # trailing token path: reassemble line by replacing the last token
                    parts = ln.rstrip("\n").split()
                    parts[-1] = "H"
                    ln = " ".join(parts) + ("\n" if ln.endswith("\n") else "")
                fixes += 1
        out.append(ln)
    return "".join(out), fixes


def assert_no_helium_in_pdbqt(
    lines: Iterable[str], ligand_name: str
) -> Tuple[List[str], int, str]:
    """
    Post-write guardrail for PDBQT:
      - If ADT='He' appears where atom NAME begins with 'H' → auto-correct to 'H', log fixes.
      - If ADT='He' appears on a non-H-named atom → quarantine is required; return reason.
    Returns (new_lines, fixes_count, quarantine_reason_or_empty).
    """
    new_lines: List[str] = []
    fixes = 0
    quarantine_reason = ""
    for ln in lines:
        if ln.startswith(("ATOM", "HETATM")):
            parts = ln.split()
            if parts:
                adt = parts[-1].strip().upper()
                if adt == "HE":
                    aname = ln[12:16].strip().upper() if len(ln) >= 16 else ""
                    if aname.startswith("H"):
                        parts[-1] = "H"
                        ln = " ".join(parts) + ("\n" if not ln.endswith("\n") else "")
                        fixes += 1
                    else:
                        quarantine_reason = f"adt_helium_inconsistent: non-H name '{aname}' has ADT=He (ligand={ligand_name})"
        new_lines.append(ln)
    return new_lines, fixes, quarantine_reason


def _normalize_element_token(elem_raw: str, rules=None) -> str | None:
    """Return canonical element token if the input looks valid; otherwise None."""
    if rules is None:
        rules = get_atom_rules()
    elem = (elem_raw or "").strip()
    if not elem:
        return None
    elem_up = elem.upper()
    if elem_up == "H":
        return "H"
    if elem_up in rules.one_letter:
        return elem_up
    if elem_up in rules.two_letter:
        return elem_up
    return None


def _format_element_token(elem_token: str, upper: bool = True) -> str:
    token = (elem_token or "").strip()
    if not token:
        return ""
    token_up = token.upper()
    if upper:
        return token_up
    if len(token_up) == 1:
        return token_up
    return token_up[0] + token_up[1].lower()


def derive_element(aname: str, resname: str, is_het: bool, rules=None) -> str:
    """Infer element symbol from atom/residue context using YAML-driven rules."""
    if rules is None:
        rules = get_atom_rules()

    an = (aname or "").strip().upper()
    rn = (resname or "").strip().upper()

    # special atom names (e.g., OXT)
    if an in rules.special_names:
        logging.debug(
            "[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
            an,
            rn,
            is_het,
            "special_names",
            rules.special_names.get(an, "?"),
        )
        return rules.special_names[an]

    # normalize residue-name aliases (e.g., IOD->I, CL- -> CL)
    if rn in rules.halide_aliases:
        logging.debug(
            "[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
            an,
            rn,
            is_het,
            "halide_alias",
            rules.halide_aliases[rn],
        )
        rn = rules.halide_aliases[rn]
    if rn in rules.cation_aliases:
        logging.debug(
            "[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
            an,
            rn,
            is_het,
            "cation_alias",
            rules.cation_aliases[rn],
        )
        rn = rules.cation_aliases[rn]

    # derive by leading functional prefix (OE1, NE2, OD1, ND2, SD, ...)
    pref = an[:2]
    if pref in rules.prefix_map:
        logging.debug(
            "[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
            an,
            rn,
            is_het,
            "prefix_map",
            rules.prefix_map.get(pref, "?"),
        )
        return rules.prefix_map[pref]
    # halide ions by residue name for single-atom HETATMs
    if is_het and rn in rules.halide_resnames:
        logging.debug(
            "[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
            an,
            rn,
            is_het,
            "halide_resname",
            rn,
        )
        if rn in {"CL", "BR"}:
            return rn[0] + rn[1].lower()
        return rn  # I, F

    # CA special-case (avoid backbone CA -> Calcium)
    if len(an) >= 2 and an[:2] == "CA":
        # Protein / peptide backbone CA: always treat as carbon
        if not is_het and rules.treat_backbone_ca:
            logging.debug(
                "[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                an,
                rn,
                is_het,
                "backbone_CA_guard",
                "C",
            )
            return "C"
        # Simple Ca2+ ions: residue names CA/CAL (pure ion residues)
        if is_het and rn in {"CA", "CAL"}:
            logging.debug(
                "[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                an,
                rn,
                is_het,
                "Ca_ion",
                "Ca",
            )
            return "Ca"

    # two-letter elements at name start (Cl, Br, Na, Mg, ...) — HETs only, and only when the
    # atom name itself looks like a stand-alone element token (length==2 or 3rd char not alpha).
    if is_het:
        two = an[:2].upper()
        looks_like_standalone = (len(an) == 2) or (len(an) >= 3 and not an[2].isalpha())
        if two in rules.two_letter and looks_like_standalone:
            # avoid mislabeling organic HET alpha carbons "CA" as Calcium.
            # If the residue is *not* a simple Ca ion (CA/CAL), and the atom name
            # is exactly "CA" (or equivalent), prefer carbon over calcium.
            if two == "CA" and rn not in {"CA", "CAL"}:
                logging.debug(
                    "[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                    an,
                    rn,
                    is_het,
                    "het_CA_non_ion_guard",
                    "C",
                )
                return "C"

            t = two
            logging.debug(
                "[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                an,
                rn,
                is_het,
                "two_letter",
                t[0] + t[1].lower(),
            )
            return t[0] + t[1].lower()

        # --- Guard: avoid mislabeling organic HET atom names like "CAK","NAA" as Ca/Na ---
        # If this is a HET but NOT a known simple ion residue, and the atom name
        # continues with another alpha character (e.g., "CAK", "NAA"), prefer a
        # one-letter element guess (C/N/...) rather than a two-letter metal.
    if is_het:
        _ion_res = {
            "LI",
            "NA",
            "K",
            "RB",
            "CS",
            "MG",
            "CA",
            "SR",
            "BA",
            "ZN",
            "CU",
            "NI",
            "CO",
            "FE",
            "MN",
            "CD",
            "AL",
            "HG",
            "AG",
            "PB",
            "PT",
            "PD",
            "AU",
            "RU",
            "IR",
            "OS",
        }
        if rn not in _ion_res:
            if len(an) >= 3 and an[0].isalpha() and an[1].isalpha() and an[2].isalpha():
                c = an[0].upper()
                if c in rules.one_letter:
                    logging.debug(
                        "[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                        an,
                        rn,
                        is_het,
                        "het_three_letters_guard",
                        c,
                    )
                    return c
    # hydrogens (H, 1H, 2H...)
    if an.startswith("H") or (an[:1].isdigit() and len(an) >= 2 and an[1] == "H"):
        logging.debug(
            "[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
            an,
            rn,
            is_het,
            "hydrogen_name",
            "H",
        )
        return "H"

    # one-letter defaults by first alpha
    if an and an[0].isalpha():
        c = an[0].upper()
        if c in rules.one_letter:
            logging.debug(
                "[element] aname=%s resn=%s is_het=%s -> via=%s => %s",
                an,
                rn,
                is_het,
                "one_letter",
                c,
            )
            return c

    for ch in an:
        if ch.isalpha():
            c = ch.upper()
            return c if c in rules.one_letter else rules.default_element

    return rules.default_element


def ensure_model_records(pdb_input_path: str, pdb_output_path: str):
    with open(pdb_input_path, "r") as f:
        lines = f.readlines()

    has_model = any(line.startswith("MODEL") for line in lines)

    if has_model:
        with open(pdb_output_path, "w") as f:
            f.writelines(lines)
        logging.info(f"MODEL record found in {pdb_input_path}. File copied unchanged.")
        return

    atom_start_idx = None
    atom_end_idx = None

    for i, line in enumerate(lines):
        if line.startswith(("ATOM  ", "HETATM")):
            atom_start_idx = i
            break

    for i in range(len(lines) - 1, -1, -1):
        if lines[i].startswith(("ATOM  ", "HETATM")):
            atom_end_idx = i
            break

    if atom_start_idx is None or atom_end_idx is None:
        raise ValueError(f"No ATOM or HETATM lines found in {pdb_input_path}.")

    lines.insert(atom_start_idx, "MODEL        1\n")
    lines.insert(atom_end_idx + 2, "ENDMDL\n")

    with open(pdb_output_path, "w") as f:
        f.writelines(lines)

    logging.info(
        f"MODEL/ENDMDL added in {pdb_output_path} between lines {atom_start_idx + 1} and {atom_end_idx + 3}."
    )


def remove_unparsable_hetatms(pdb_path):
    cleaned_lines = []
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                if "UNK" in atom_name or "UNX" in line:
                    continue
            cleaned_lines.append(line)
    with open(pdb_path, "w") as f:
        f.writelines(cleaned_lines)
    logging.info(f"Unparsable HETATM entries removed from {pdb_path}.")


class ElementFixer(Select):
    """Rewrite element symbols safely, with peptide-like HET awareness (YAML-driven)."""

    def __init__(self):
        super().__init__()
        self.rules = get_atom_rules()

    def _is_peptidic_like(self, residue) -> bool:
        try:
            names = [a.get_name().strip().upper() for a in residue.get_atoms()]
        except Exception:
            return False
        if not names:
            return False
        hits = sum(
            (n in self.rules.peptide_like) or (n[:2] in {"OE", "NE", "OD", "ND", "SD"})
            for n in names
        )
        return hits >= max(4, int(0.6 * len(names)))

    def get_atom_element(self, atom):
        aname = atom.get_name()
        res = atom.get_parent()
        rname = getattr(res, "get_resname", lambda: "")()

        try:
            hetflag = res.get_id()[0] if res is not None else " "
        except Exception:
            hetflag = " "
        is_het = hetflag != " "

        elem_raw = getattr(atom, "element", "")
        elem_token = _normalize_element_token(elem_raw, self.rules)
        if elem_token is not None:
            return _format_element_token(elem_token, upper=False)

        # if residue looks peptide-like, still derive by name (prevents ion mislabels)
        if self._is_peptidic_like(res):
            el = derive_element(aname, rname, is_het, self.rules)
        else:
            el = derive_element(aname, rname, is_het, self.rules)

        if aname.strip().upper().startswith("H") and str(el).strip() in {
            "He",
            "HE",
            "he",
        }:
            el = "H"
        return el

    def accept_atom(self, atom):
        atom.element = self.get_atom_element(atom)
        return True


def fix_pdb_elements(input_path, output_path=None):
    parser = PDBParser(QUIET=True)
    io = PDBIO()
    structure = parser.get_structure("structure", input_path)
    io.set_structure(structure)

    fixer = ElementFixer()

    if output_path is None:
        output_path = input_path

    with open(output_path, "w") as fh:
        io.save(fh, fixer)
    print(f"Saved fixed PDB to {output_path}")


def fix_ligand_element_columns(lines: List[str]) -> List[str]:
    """Rewrite element cols (77–78) for HET ligands using YAML rules."""
    rules = get_atom_rules()
    out = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 78:
            aname = line[12:16]
            resn = line[17:20]
            is_het = line.startswith("HETATM")
            elem_token = _normalize_element_token(line[76:78], rules)
            if elem_token is None:
                # Fall back only when the element column is missing or invalid.
                el = derive_element(aname, resn, is_het, rules)
                # Correct PTR-style hydrogens mislabeled as Helium (HE1/HE2 → H) without touching real metals.
                if aname.strip().upper().startswith("H") and str(el).strip() in {
                    "He",
                    "HE",
                    "he",
                }:
                    el = "H"
            else:
                el = elem_token
            # [elem-normalize] Force uppercase 2-char element slot before write.
            el_clean = _format_element_token(el, upper=True)
            el_clean = el_clean[:2]
            line = line[:76] + f"{el_clean:>2}" + line[78:]

        out.append(line)
    return out


def fix_element_columns_in_file(src_path, dst_path=None, rewrite_atoms=False):
    rules = get_atom_rules()
    src_path = str(src_path)
    dst_path = src_path if dst_path is None else str(dst_path)
    logging.info("[elemfix] in=%s rewrite_atoms=%s", src_path, rewrite_atoms)
    out_lines = []
    with open(src_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            is_atom = line.startswith("ATOM  ")
            is_het = line.startswith("HETATM")
            if len(line) >= 78 and (is_het or (rewrite_atoms and is_atom)):
                aname = line[12:16]
                resn = line[17:20]
                elem_token = _normalize_element_token(line[76:78], rules)
                if elem_token is None:
                    # Fall back only when the element column is missing or invalid.
                    el = derive_element(aname, resn, is_het, rules)
                    # Correct PTR-style hydrogens mislabeled as Helium (HE1/HE2 → H) without touching real metals.
                    if str(aname).strip().upper().startswith("H") and str(
                        el
                    ).strip() in {"He", "HE", "he"}:
                        el = "H"
                else:
                    el = elem_token
                # [elem-normalize] Align on uppercase 2-char element column for files.
                el_clean = _format_element_token(el, upper=True)
                el_clean = el_clean[:2]
                line = line[:76] + f"{el_clean:>2}" + line[78:]

            out_lines.append(line)
    with open(dst_path, "w", encoding="utf-8") as out:
        out.writelines(out_lines)
    logging.info("[elemfix] done=%s", dst_path or src_path)
    return dst_path


# --- Canonical residue loader helpers (aliases-migration) ---
_TOKEN_SPLIT = re.compile(r"[;\s,]+")


def _iter_values(payload: Any) -> Iterable[Any]:
    if payload is None:
        return []
    if isinstance(payload, Mapping):
        items: list[Any] = []
        for key, value in payload.items():
            items.append(key)
            items.extend(_iter_values(value))
        return items
    if isinstance(payload, (str, bytes)):
        return [payload]
    if isinstance(payload, Iterable):
        items: list[Any] = []
        for entry in payload:
            items.extend(_iter_values(entry))
        return items
    return [payload]


def _load_default_alias_cfg() -> Mapping[str, Any]:
    global _aliases_cache
    if _aliases_cache is not None:
        cached = _aliases_cache
        if isinstance(cached, Mapping):
            return cached
    try:
        with open(ALIASES_PATH, "r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}
    except Exception:
        data = {}
    if not isinstance(data, Mapping):
        data = {}
    _aliases_cache = data
    return data


def _normalize_tokens(raw_tokens: Iterable[Any]) -> Set[str]:
    tokens: Set[str] = set()
    for raw in raw_tokens:
        if raw is None:
            continue
        text = str(raw).strip()
        if not text:
            continue
        for piece in _TOKEN_SPLIT.split(text):
            token = piece.strip()
            if token:
                tokens.add(token.upper())
    return tokens


def _collect_alias_synonyms(alias_map: Any, canonical: Set[str]) -> Set[str]:
    if not isinstance(alias_map, Mapping):
        return set()
    synonyms: Set[str] = set()
    for alias, target in alias_map.items():
        alias_tokens = _normalize_tokens([alias])
        target_tokens = _normalize_tokens([target])
        if canonical & target_tokens:
            synonyms |= alias_tokens
    return synonyms


def _resolve_alias_cfg(cfg: Mapping[str, Any] | None) -> Mapping[str, Any]:
    """
    Decide which mapping to treat as the alias config for canonical metals/cofactors/waters.

    - If the caller passes a mapping that already looks like an aliases.yaml payload
      (has canonical_* keys or element_sets), use it.
    - Otherwise, fall back to the canonical chemdb/aliases.yaml on disk.
    """
    if isinstance(cfg, Mapping):
        if any(
            key in cfg
            for key in (
                "canonical_metals",
                "canonical_cofactors",
                "canonical_waters",
                "element_sets",
            )
        ):
            return cfg
    return _load_default_alias_cfg()


# Public loaders --------------------------------------------------------------
def load_canonical_metals(cfg: Mapping[str, Any] | None) -> Set[str]:
    cfg_map = _resolve_alias_cfg(cfg)
    canonical = _normalize_tokens(_iter_values(cfg_map.get("canonical_metals")))

    element_sets = cfg_map.get("element_sets")
    if isinstance(element_sets, Mapping):
        canonical |= _collect_alias_synonyms(
            element_sets.get("cation_resname_aliases"), canonical
        )
        canonical |= _collect_alias_synonyms(
            element_sets.get("halide_resname_aliases"), canonical
        )

    canonical |= _collect_alias_synonyms(
        cfg_map.get("retain_element_alias_map"), canonical
    )
    return {token.strip().upper() for token in canonical if token.strip()}


def load_canonical_cofactors(cfg: Mapping[str, Any] | None) -> Set[str]:
    cfg_map = _resolve_alias_cfg(cfg)
    return _normalize_tokens(_iter_values(cfg_map.get("canonical_cofactors")))


def load_canonical_waters(cfg: Mapping[str, Any] | None) -> Set[str]:
    cfg_map = _resolve_alias_cfg(cfg)
    return _normalize_tokens(_iter_values(cfg_map.get("canonical_waters")))
