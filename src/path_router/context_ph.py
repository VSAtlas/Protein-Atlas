"""
Context-driven pH selector for receptor preparation.

Usage:
  python -m protein_automation.context_ph /path/to/structure.pdb [--buffer-ph 7.4] [--compartment lysosome] [--assay-ph 6.8] [--print-provenance]

Programmatic:
  from protein_automation.context_ph import select_ph_from_pdb
  result = select_ph_from_pdb(pdb_path, buffer_ph=7.4, project_compartment="lysosome")
  print(result["target_pH"], result["provenance"], result["confidence"], result["ensemble"])
"""

import os
import re
import json
import argparse
from typing import Dict, List, Tuple, Optional
import logging

# >>> PATHS IMPORT START
from .path_router import make_paths
# >>> PATHS IMPORT END

logger = logging.getLogger("context_ph")
if not logger.handlers:
    _h = logging.StreamHandler()
    _h.setFormatter(logging.Formatter("[context_ph] %(levelname)s - %(message)s"))
    logger.addHandler(_h)
logger.setLevel(logging.INFO)
# -------------------------
# Presets loading
# -------------------------


def _project_root() -> str:
    here = os.path.abspath(os.path.dirname(__file__))
    return os.path.abspath(os.path.join(here))


def _presets_path() -> str:
    # chemdb/context_ph_presets.json relative to this module
    return os.path.join(_project_root(), "chemdb", "context_ph_presets.json")


def load_presets(path: Optional[str] = None) -> Dict:
    p = path or _presets_path()
    try:
        with open(p, "r") as f:
            data = json.load(f)
    except Exception as e:
        logger.error("Failed to load presets from %s: %s", p, e)
        # Minimal safe fallback
        return {
            "compartment_pH": {"extracellular": [7.35, 7.45]},
            "aliases": {},
            "weights": {"default": 0.5},
            "default_pH": 7.0,
        }

    # Minimal schema validation + normalization
    if "compartment_pH" not in data or not isinstance(data["compartment_pH"], dict):
        logger.warning("Presets missing 'compartment_pH'; adding fallback.")
        data["compartment_pH"] = {"extracellular": [7.35, 7.45]}

    # Normalize ranges and sanity check values
    clean_map = {}
    for k, v in data["compartment_pH"].items():
        try:
            lo, hi = float(v[0]), float(v[1])
            if hi < lo:
                lo, hi = hi, lo
            # clip and warn if extreme
            lo_c = max(0.0, min(14.0, lo))
            hi_c = max(0.0, min(14.0, hi))
            if (lo_c, hi_c) != (lo, hi):
                logger.debug(
                    "Clipped pH range for %s from (%.2f, %.2f) to (%.2f, %.2f)",
                    k,
                    lo,
                    hi,
                    lo_c,
                    hi_c,
                )
            clean_map[_norm(k)] = [lo_c, hi_c]
        except Exception:
            logger.warning("Ignoring invalid pH range for %r: %r", k, v)
    data["compartment_pH"] = clean_map

    if "aliases" in data and isinstance(data["aliases"], dict):
        data["aliases"] = {_norm(k): _norm(v) for k, v in data["aliases"].items()}
    else:
        data["aliases"] = {}

    if "weights" not in data or not isinstance(data["weights"], dict):
        data["weights"] = {"default": 0.5}
    if "default_pH" not in data:
        data["default_pH"] = 7.0
    # validate  location_tokens(phrases + singles)  from JSON
    lt = data.get("location_tokens", {})
    phrases = lt.get("phrases", []) if isinstance(lt, dict) else []
    singles = lt.get("singles", []) if isinstance(lt, dict) else []
    if not isinstance(phrases, list):
        phrases = []
    if not isinstance(singles, list):
        singles = []

    # normalize to lower, dedupe, drop empties
    def _norm_list(xs):
        out = []
        seen = set()
        for x in xs:
            s = (x or "").strip().lower()
            if not s or s in seen:
                continue
            seen.add(s)
            out.append(s)
        return out

    data["location_tokens"] = {
        "phrases": _norm_list(phrases),
        "singles": _norm_list(singles),
    }

    logger.debug(
        "Loaded presets: %d compartments, %d aliases, %d phrases, %d singles",
        len(data["compartment_pH"]),
        len(data["aliases"]),
        len(data["location_tokens"]["phrases"]),
        len(data["location_tokens"]["singles"]),
    )
    return data


# -------------------------
# Small utilities
# -------------------------


def _mid(a: float, b: float) -> float:
    return (a + b) / 2.0


def _clip(x: float, lo: float = 3.0, hi: float = 10.5) -> float:
    return max(lo, min(hi, x))


def _round01(x: float) -> float:
    return round(x, 1)


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")


def _cand(name: str, value: float, weight: float, note: str) -> Dict:
    return {
        "source": name,
        "pH": float(round(_clip(value), 3)),
        "w": float(weight),
        "note": str(note),
    }


# -------------------------
# Header parsing (PDB / AlphaFold)
# -------------------------

_PH_RE = re.compile(r"\bPH\s*([0-9]+(?:\.[0-9]+)?)", re.IGNORECASE)


def parse_pdb_header_text(pdb_path: str, max_lines: int = 2000) -> str:
    """
    Read header-only (stop at first ATOM/HETATM) with robust error handling.
    """
    lines: List[str] = []
    try:
        with open(pdb_path, "r", errors="ignore") as fh:
            for i, line in enumerate(fh):
                if line.startswith(("ATOM", "HETATM")) or i >= max_lines:
                    break
                lines.append(line.rstrip("\n"))
    except Exception as e:
        logger.warning("Failed reading header from %s: %s", pdb_path, e)
        return ""
    hdr = "\n".join(lines)
    logger.debug("Header lines read: %d", len(lines))
    return hdr


def header_find_explicit_ph(header_text: str) -> Optional[float]:
    """
    Look for crystallization/buffer pH in REMARK 200/280 etc. Pattern: "PH 7.5"
    """
    m = _PH_RE.search(header_text)
    return float(m.group(1)) if m else None


def header_detect_alphafold(header_text: str) -> bool:
    """
    AlphaFold PDBs typically contain 'AlphaFold' in TITLE/HEADER/COMPND, and DBREF/DBSOURCE UNP.
    """
    ht = header_text.lower()
    # Common indicators
    if "alphafold" in ht or "model generated by alphafold" in ht:
        return True
    if "dbref" in ht and "uniprot" in ht:
        # Many AF PDBs include DBREF UNP mapping
        return True
    return False


def header_extract_uniprot_accessions(header_text: str) -> List[str]:
    """
    Extract UniProt accessions from DBREF/COMPND/DBSOURCE-like lines.
    Very simple heuristic; returns list of candidate accessions.
    """
    accs = set()
    # DBREF ... UNP <ACCN> <NAME>
    for line in header_text.splitlines():
        if line.startswith(
            ("DBREF", "DBREF1", "DBREF2", "DBREF3", "COMPND", "DBSOURCE", "REMARK")
        ):
            # Look for UniProt-style accessions: 1-2 letters + 4-5 alphanum, possibly with a dash for isoforms
            for m in re.finditer(
                r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])(-\d+)?\b",
                line,
            ):
                accs.add(m.group(0))
            # Also detect tokens after 'UNP ' or 'UNIPROT'
            m2 = re.search(
                r"\bUNP(?:ROTK?)?\s*[:=]?\s*([A-Z0-9\-]+)", line, re.IGNORECASE
            )
            if m2:
                token = m2.group(1).strip()
                if re.match(r"^[A-Z0-9\-]{5,10}$", token):
                    accs.add(token)
    return sorted(accs)


def header_extract_location_keywords(header_text: str, presets: Dict) -> List[str]:
    """
    Scrape subcellular/location words from COMPND/KEYWDS/TITLE/SOURCE blocks,
    using word-boundary regex and the token lists provided in presets["location_tokens"].
    """
    loc_blocks: List[str] = []
    for line in header_text.splitlines():
        prefix = line[:6].strip()
        if prefix in ("COMPND", "KEYWDS", "TITLE", "SOURCE"):
            loc_blocks.append(line[10:].strip().lower())

    blob = " ".join(loc_blocks)
    # Normalize: non-alnum -> spaces so \b works reliably
    norm = re.sub(r"[^a-z0-9]+", " ", blob)

    lt = presets.get("location_tokens", {})
    phrases: List[str] = lt.get("phrases", []) or []
    singles: List[str] = lt.get("singles", []) or []

    hits: List[str] = []

    # Match phrases as whole phrases via regex word boundaries
    for ph in phrases:
        # e.g., "cell surface" -> r"\bcell\s+surface\b"
        pat = r"\b" + re.sub(r"\s+", r"\\s+", ph) + r"\b"
        if re.search(pat, norm):
            hits.append(ph)

    # Match single tokens with strict word boundaries
    for tok in singles:
        if re.search(rf"\b{re.escape(tok)}\b", norm):
            hits.append(tok)

    dedup = sorted(set(hits))
    if dedup:
        logger.debug("Location keyword hits: %s", dedup)
    return dedup


# -------------------------
# Mapping locations ? pH (presets)
# -------------------------


def ph_from_compartment(name: str, presets: Dict) -> Optional[float]:
    comp_map = presets.get("compartment_pH", {})
    key = _norm(name)
    if key in comp_map:
        lo, hi = comp_map[key]
        return round(_mid(lo, hi), 2)
    # Fuzzy aliases
    aliases = presets.get("aliases", {})
    if key in aliases:
        alias_key = aliases[key]
        if alias_key in comp_map:
            lo, hi = comp_map[alias_key]
            return round(_mid(lo, hi), 2)
    return None


# -------------------------
# Main selection logic
# -------------------------


def select_ph(signals: Dict, presets: Dict) -> Dict:
    """
    Combine heterogeneous signals to produce target pH, provenance, confidence, and optional ensemble.

    signals:
      TARGET_PH: Optional[float]
      ASSAY_PH / BUFFER_PH: Optional[float]
      PDB_HEADER_TEXT: Optional[str]
      PROJECT_COMPARTMENT: Optional[str]
      UNIPROT_LOC_HINTS: Optional[List[str]]  # extra location terms you found elsewhere
      ENZYME_OPT_RANGE: Optional[Tuple[float,float]]
      DEFAULT_PH: Optional[float]

    returns:
      { "target_pH": float, "provenance": List[dict], "confidence": float, "ensemble": Optional[List[float]] }
    """
    W = presets["weights"]
    cands: List[Dict] = []

    # 1) Explicit override
    if (v := signals.get("TARGET_PH")) is not None:
        return {
            "target_pH": float(v),
            "provenance": [_cand("override", float(v), W["override"], "explicit")],
            "confidence": 1.0,
            "ensemble": None,
        }

    # 2) Assay/buffer
    for key in ("ASSAY_PH", "BUFFER_PH"):
        if (v := signals.get(key)) is not None:
            cands.append(_cand("assay", float(v), W["assay"], key))

    # 3) PDB header explicit PH
    header = signals.get("PDB_HEADER_TEXT") or ""
    if header:
        if (v := header_find_explicit_ph(header)) is not None:
            cands.append(_cand("pdb_header", float(v), W["pdb_header"], "REMARK PH"))

    # 4) Compartment from project hint
    if pc := signals.get("PROJECT_COMPARTMENT"):
        v = ph_from_compartment(pc, presets)
        if v is not None:
            cands.append(_cand("project_compartment", v, W["project_compartment"], pc))

    # 5) Location keywords → compartment mapping (JSON-driven)
    loc_terms = list(signals.get("UNIPROT_LOC_HINTS") or [])
    # Also scrape from header automatically using presets tokens
    loc_terms += header_extract_location_keywords(header, presets) if header else []
    for term in sorted(set(loc_terms)):
        v = ph_from_compartment(term, presets)
        if v is not None:
            cands.append(_cand("uniprot_go", v, W["uniprot_go"], term))

    # 6) Enzyme optimum range (if known)
    if rng := signals.get("ENZYME_OPT_RANGE"):
        lo, hi = rng
        cands.append(
            _cand(
                "enzyme_optimum",
                _mid(float(lo), float(hi)),
                W["enzyme_optimum"],
                f"{lo}-{hi}",
            )
        )

    # 7) Default if still nothing
    if not cands:
        v = float(signals.get("DEFAULT_PH", presets.get("default_pH", 7.0)))
        return {
            "target_pH": v,
            "provenance": [_cand("default", v, W["default"], "no signals")],
            "confidence": W["default"],
            "ensemble": None,
        }

    # Choose winner / blend
    cands_sorted = sorted(cands, key=lambda d: d["w"], reverse=True)
    if logger.isEnabledFor(logging.DEBUG):
        logger.debug(
            "Candidates (n=%d): %s",
            len(cands_sorted),
            "; ".join(f"{d['source']}@{d['pH']} (w={d['w']})" for d in cands_sorted),
        )

    top = cands_sorted[0]
    second = cands_sorted[1] if len(cands_sorted) > 1 else None

    # If a clearly dominant source exists, pick it
    if second and (top["w"] - second["w"] >= 0.2):
        target = _round01(top["pH"])
        logger.info(
            "select_ph: dominant=%s pH=%.1f (w=%.2f)", top["source"], target, top["w"]
        )
        return {
            "target_pH": target,
            "provenance": cands_sorted,
            "confidence": min(1.0, top["w"]),
            "ensemble": None,
        }

    # Otherwise weighted blend
    wsum = sum(d["w"] for d in cands_sorted)
    target = _round01(
        sum(d["pH"] * d["w"] for d in cands_sorted) / (wsum if wsum else 1.0)
    )
    span = max(d["pH"] for d in cands_sorted) - min(d["pH"] for d in cands_sorted)

    if span <= 0.3:
        logger.info("select_ph: blended pH=%.1f span=%.2f (no ensemble)", target, span)
        return {
            "target_pH": target,
            "provenance": cands_sorted,
            "confidence": min(0.9, wsum / max(1.0, len(cands_sorted))),
            "ensemble": None,
        }
    elif span <= 1.0:
        alt = target + 0.5 if target < top["pH"] else target - 0.5
        ens = sorted({_clip(target), _clip(_round01(alt))})
        logger.info(
            "select_ph: blended pH=%.1f span=%.2f → ensemble=%s",
            target,
            span,
            [float(x) for x in ens],
        )
        return {
            "target_pH": target,
            "provenance": cands_sorted,
            "confidence": 0.7,
            "ensemble": [float(x) for x in ens],
        }
    else:
        ens = sorted(
            {
                _clip(_round01(target - 0.5)),
                _clip(target),
                _clip(_round01(target + 0.5)),
            }
        )
        logger.info(
            "select_ph: blended pH=%.1f span=%.2f → ensemble=%s",
            target,
            span,
            [float(x) for x in ens],
        )
        return {
            "target_pH": target,
            "provenance": cands_sorted,
            "confidence": 0.6,
            "ensemble": [float(x) for x in ens],
        }


# -------------------------
# High-level convenience
# -------------------------


def select_ph_from_pdb(
    pdb_path: str,
    buffer_ph: Optional[float] = None,
    assay_ph: Optional[float] = None,
    project_compartment: Optional[str] = None,
    enzyme_opt_range: Optional[Tuple[float, float]] = None,
    default_ph: float = 7.0,
    extra_location_hints: Optional[List[str]] = None,
    presets_path: Optional[str] = None,
) -> Dict:
    """
    One-shot: parse the file, harvest signals, and select pH.
    """
    presets = load_presets(presets_path)
    header = parse_pdb_header_text(pdb_path)
    is_af = header_detect_alphafold(header)
    uniprot_accs = header_extract_uniprot_accessions(header) if is_af else []

    signals = {
        "ASSAY_PH": assay_ph if assay_ph is not None else buffer_ph,
        "BUFFER_PH": buffer_ph,
        "PDB_HEADER_TEXT": header,
        "PROJECT_COMPARTMENT": project_compartment,
        "ENZYME_OPT_RANGE": enzyme_opt_range,
        "DEFAULT_PH": default_ph,
        "UNIPROT_LOC_HINTS": extra_location_hints or [],
    }

    # NOTE: We don't web-fetch UniProt; if you have a sidecar annotator,
    # pass its terms via extra_location_hints.
    # Still, AF header scraping may contribute location tokens already.

    result = select_ph(signals, presets)
    # Attach a couple of introspection fields
    result["is_alphafold"] = bool(is_af)
    result["uniprot_accessions"] = uniprot_accs
    return result


# >>> CONTEXTPH PATHS PATCH START
def select_ph_for_pdbid(cfg, pdb_id: str, **kwargs) -> Dict:
    """
    Resolve input PDB via path router and call select_ph_from_pdb.
    kwargs are forwarded to select_ph_from_pdb (buffer_ph, assay_ph, project_compartment, ...).
    """
    p = make_paths(cfg, base_id=pdb_id, pdb_file=f"{pdb_id}.pdb")
    return select_ph_from_pdb(str(p.input_pdb_path), **kwargs)


# >>> CONTEXTPH PATHS PATCH END

# -------------------------
# CLI
# -------------------------


def _build_cli():
    ap = argparse.ArgumentParser(
        description="Infer environmental pH from PDB/AlphaFold headers + context hints."
    )
    ap.add_argument("pdb", help="Path to PDB file (AlphaFold or experimental).")
    ap.add_argument(
        "--buffer-ph", type=float, default=None, help="Assay/buffer pH hint."
    )
    ap.add_argument(
        "--assay-ph",
        type=float,
        default=None,
        help="Assay pH hint (overrides buffer if both).",
    )
    ap.add_argument(
        "--compartment",
        type=str,
        default=None,
        help="Project compartment hint (e.g., lysosome, cytosol).",
    )
    ap.add_argument(
        "--enzyme-opt",
        type=str,
        default=None,
        help="Enzyme optimum range 'lo,hi' (e.g., '5.5,6.5').",
    )
    ap.add_argument(
        "--default-ph",
        type=float,
        default=7.0,
        help="Fallback pH when no signals found.",
    )
    ap.add_argument(
        "--presets",
        type=str,
        default=None,
        help="Path to context_ph_presets.json (optional).",
    )
    ap.add_argument(
        "--print-provenance", action="store_true", help="Print detailed provenance."
    )
    ap.add_argument(
        "--verbose", "-v", action="store_true", help="Enable DEBUG logging."
    )
    ap.add_argument("--json", action="store_true", help="Emit JSON result to stdout.")

    return ap


def _parse_opt_range(s: Optional[str]) -> Optional[Tuple[float, float]]:
    if not s:
        return None
    m = re.match(r"\s*([0-9.]+)\s*,\s*([0-9.]+)\s*$", s)
    if not m:
        raise SystemExit(f"--enzyme-opt must be 'lo,hi', got {s!r}")
    lo, hi = float(m.group(1)), float(m.group(2))
    if hi < lo:
        lo, hi = hi, lo
    return (lo, hi)


def select_ph_values_for_protonation(pdb_path: str, **kwargs) -> List[float]:
    """
    Returns a list of pH values to use (ensemble if present, else [target]).
    kwargs are passed to select_ph_from_pdb (buffer_ph, assay_ph, project_compartment, ...).
    """
    res = select_ph_from_pdb(pdb_path, **kwargs)
    return list(res.get("ensemble") or [res["target_pH"]])


def main():
    ap = _build_cli()
    args = ap.parse_args()
    rng = _parse_opt_range(args.enzyme_opt)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        logger.debug("Verbose logging enabled")

    res = select_ph_from_pdb(
        pdb_path=args.pdb,
        buffer_ph=args.buffer_ph,
        assay_ph=args.assay_ph,
        project_compartment=args.compartment,
        enzyme_opt_range=rng,
        default_ph=args.default_ph,
        presets_path=args.presets,
    )
    if args.json:
        import json as _json

        print(_json.dumps(res, indent=2, sort_keys=True))
        return
    # Minimal, stable text output
    print(f"target_pH: {res['target_pH']}")
    print(f"confidence: {res['confidence']:.2f}")
    if res.get("ensemble"):
        print("ensemble:", ", ".join(str(x) for x in res["ensemble"]))
    print(f"is_alphafold: {res.get('is_alphafold', False)}")
    if res.get("uniprot_accessions"):
        print("uniprot_accessions:", ", ".join(res["uniprot_accessions"]))
    if args.print_provenance:
        print("provenance:")
        for d in res["provenance"]:
            print(f"  - {d['source']}: pH={d['pH']} w={d['w']} note={d['note']}")


if __name__ == "__main__":
    main()
