#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Atlas Log Resurfacer - single entry (v2, ASCII-safe)

Default flow (no flags):
  1) Ensure ./logs exists
  2) Run benchmark_mode.py, tee stdout/stderr to ./logs/bench_<YYYYmmdd_HHMMSS>.log
  3) Resurface with:
     --proteins ./docked
     --ligands  ./prepped_ligands
     --since    <bench ts as YYYY-MM-DD HH:MM:SS>
     (write JSON + TSV + MD to log_resurfacer/out/<ts>/)
  4) Print concise console summary and [paths] footer
  5) Exit policy: --fail-on hard (non-zero if any hard/fatal)

Examples
--------
python log_resurfacer/resurface.py
python log_resurfacer/resurface.py --no-benchmark --bench-log logs/bench_20251007_153012.log
python log_resurfacer/resurface.py -- --quick --pdbs 1T46,5MO4
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import subprocess, os
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

# -----------------------
# Default rules (embedded)
# -----------------------
DEFAULT_RULES = {
    "rules": [
        # A. Docking/scoring & render pipeline (protein.log)
        {
            "id": "dock.no_score",
            "severity": "hard",
            "regex": r"^.*WARNING\s*-\s*No score for (?P<lig>[^\s]+?\.pdbqt)",
            "hint": "Check ligand_prep_status.tsv for this ligand; verify docking ran; inspect pocket center/grid logs.",
            "priority": 100,
        },
        {
            "id": "render.no_ctrl_rdk_paths",
            "severity": "soft",
            "regex": r"\[render-skip].*no ctrl/rdk pose paths",
            "hint": "Why no valid poses? Correlate with 'No score'/validation; ensure control selection produced a pose.",
            "priority": 10,
        },
        {
            "id": "render.enqueue_none",
            "severity": "soft",
            "regex": r"\[render-enqueue].*ctrl=None\s*\|\s*rdk=None",
            "hint": "Control selection failed or disabled; check 'hint_count' and control whitelist.",
            "priority": 10,
        },
        {
            "id": "recentering.skipped_zero",
            "severity": "soft",
            "regex": r"Early recenter skipped: evaluated=0\s*<\s*threshold",
            "hint": "No evaluated poses to guide recentering; widen first-pass or ensure at least one pose completes.",
            "priority": 8,
        },
        {
            "id": "cap.limiting",
            "severity": "info",
            "regex": r"^\[CAP]\s+.*limiting non-controls.*",
            "hint": "If controls fail, few whitelist ligands might still succeed; not an error by itself.",
            "priority": 1,
        },
        # B. Ligand preparation (benchmark & TSV corroboration)
        {
            "id": "ligprep.rdkit_valence",
            "severity": "hard",
            "regex": r"Explicit valence for atom|Can't kekulize mol|sanitize .*\.pdb",
            "hint": "Quarantine; consider pre-sanitize with RDKit, fix aromaticity, add hydrogens earlier.",
            "priority": 90,
        },
        {
            "id": "ligprep.obabel_h_charge",
            "severity": "soft",
            "regex": r"no explicit hydrogens.*needed for formal charge estimation",
            "hint": "Run Reduce/OpenBabel addH; if Reduce already ran, try -nohyd branch then addH fallback.",
            "priority": 5,
        },
        {
            "id": "ligprep.mgl_missing_mol2",
            "severity": "fatal",
            "regex": r"doesn'?t exist .*\.mol2",
            "hint": "OBabel write failed or source missing; retry fallback path; if still missing, quarantine.",
            "priority": 110,
        },
        {
            "id": "ligprep.mgl_partial_write",
            "severity": "hard",
            "regex": r"were not written",
            "hint": "Conversion partial; capture stderr, retry fallback, then quarantine if persistent.",
            "priority": 80,
        },
        {
            "id": "ligprep.pipeline_bug_scope",
            "severity": "hard",
            "regex": r"Could not create intermediates symlink: local variable 'pdb_file' referenced before assignment",
            "hint": "Initialize pdb_file before symlink block; this masks real prep outcomes.",
            "priority": 70,
        },
        # TSV-reported failures are handled programmatically and labeled ligprep.tsv_failure
    ]
}

SEVERITY_ORDER = {"fatal": 3, "hard": 2, "soft": 1, "info": 0}
TIMESTAMP_RE = re.compile(r"(?P<ts>(?:20\d{2}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2}))")


# Map rule IDs to pipeline steps
RULE_TO_STEP = {
    "ligprep.rdkit_valence":     "ligprep.sanitize",
    "ligprep.mgl_partial_write": "ligprep.write_mol2",
    "ligprep.obabel_h_charge":   "ligprep.add_h",
    "ligprep.mgl_missing_mol2":  "ligprep.write_mol2",
    "ligprep.tsv_failure":       "ligprep.parse",
    "dock.no_score":             "dock.score",
    "render.no_ctrl_rdk_paths":  "render",
    "render.enqueue_none":       "render",
    "recentering.skipped_zero":  "dock.run",
    "cap.limiting":              "dock.run",
}






# -----------------------
# Data structures
# -----------------------
@dataclass
class Rule:
    id: str
    severity: str
    regex: re.Pattern
    hint: str
    priority: int = 0


@dataclass
class FileSample:
    path: str
    samples: List[str] = field(default_factory=list)


@dataclass
class GroupItem:
    rule_id: str
    severity: str
    hint: str
    pdb: str
    ligand: Optional[str]
    file_basename: str
    count: int = 0
    first_seen: Optional[datetime] = None
    last_seen: Optional[datetime] = None
    files: Dict[str, FileSample] = field(default_factory=dict)
    joins: List[Dict[str, str]] = field(default_factory=list)
    priority: int = 0
    had_timed: bool = False
    step: str = ""
    paths: Dict[str, List[str]] = field(default_factory=dict)

    def add_hit(self, path: str, line: str, ts: Optional[datetime]):
        self.count += 1
        if path not in self.files:
            self.files[path] = FileSample(path=path, samples=[])
        if len(self.files[path].samples) < 3:
            self.files[path].samples.append(line.strip())
        if ts:
            self.had_timed = True
            if not self.first_seen or ts < self.first_seen:
                self.first_seen = ts
            if not self.last_seen or ts > self.last_seen:
                self.last_seen = ts

    @property
    def recency(self) -> str:
        if self.last_seen:
            return self.last_seen.strftime("%Y-%m-%d %H:%M:%S")
        return ""


# -----------------------
# Helpers
# -----------------------
def _parse_rules(path: Optional[str]) -> List[Rule]:
    data = DEFAULT_RULES if not path else json.loads(Path(path).read_text(encoding="utf-8"))
    return [
        Rule(
            id=r["id"],
            severity=r["severity"],
            regex=re.compile(r["regex"]),
            hint=r.get("hint", ""),
            priority=int(r.get("priority", 0)),
        )
        for r in data.get("rules", [])
    ]


def _iter_lines(path: Path) -> Iterable[Tuple[Optional[datetime], str]]:
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                ts = None
                m = TIMESTAMP_RE.search(line)
                if m:
                    try:
                        ts = datetime.fromisoformat(m.group("ts").replace("T", " "))
                    except Exception:
                        ts = None
                yield ts, line.rstrip("\n")
    except FileNotFoundError:
        return


def _load_ligand_tsv(tsv_path: Path) -> Dict[str, Dict[str, str]]:
    results: Dict[str, Dict[str, str]] = {}
    if not tsv_path.exists():
        return results
    with open(tsv_path, "r", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = {h.lower(): h for h in (reader.fieldnames or [])}
        lid = headers.get("ligand_id") or headers.get("ligand") or "ligand_id"
        status_col = headers.get("status") or "status"
        reason_col = headers.get("reason") or "reason"
        for row in reader:
            ligand = (row.get(lid) or "").strip()
            if not ligand:
                continue
            base = Path(ligand).name
            base_noext = Path(base).stem
            status = (row.get(status_col) or "").strip().lower()
            reason = (row.get(reason_col) or "").strip()
            payload = {"prep_status": status, "reason": reason}
            results[base] = payload
            results[base_noext] = payload
    return results


def _severity_at_or_above(sev: str, threshold: str) -> bool:
    return SEVERITY_ORDER.get(sev, -1) >= SEVERITY_ORDER.get(threshold, 99)

def _shorten(s: str, width: int = 64) -> str:
    if not s:
        return ""
    if len(s) <= width:
        return s
    # center-ellipsis to keep head/tail informative
    keep = width - 3
    head = keep // 2
    tail = keep - head
    return f"{s[:head]}...{s[-tail:]}"

def _pick_representative_path(g) -> str:
    """
    Prefer post-prep artifacts, then raw, then dock, then bench log.
    """
    if not getattr(g, "paths", None):
        return ""
    order = [
        "post_pdbqt",
        "post_mol2",
        "pre_raw_pdb_sanitized",
        "pre_raw_pdb",
        "protein_log",
        "dock_inputs",
        "bench_logs",
    ]
    for key in order:
        hits = g.paths.get(key, [])
        if hits:
            return hits[0]
    return ""


import re

_PDB_RX = re.compile(r"/(?:processed_pdbs|docked)/([^/]+)/")
_LIG_RX = re.compile(r"([A-Za-z0-9_.+-]+)\.(?:pdbqt|mol2|pdb)(?:\s|$)")
# also catch sanitized variants
_LIG_SAN_RX = re.compile(r"([A-Za-z0-9_.+-]+)\.sanitized\.pdb(?:\s|$)")

def infer_missing_ids(items: list):
    """Fill in g.pdb / g.ligand from sample lines and file paths if missing."""
    for g in items:
        if g.pdb and g.ligand:
            continue
        # Walk all samples from all files for this group
        samples = []
        for fs in g.files.values():
            samples.extend(fs.samples)
            samples.append(fs.path)  # also try the file path itself
        # Infer PDB
        if not g.pdb:
            for s in samples:
                m = _PDB_RX.search(s)
                if m:
                    g.pdb = m.group(1)
                    break
        # Infer ligand
        if not g.ligand:
            for s in samples:
                m = _LIG_RX.search(s) or _LIG_SAN_RX.search(s)
                if m:
                    g.ligand = m.group(1)
                    break
    # Final pass: if still missing PDB but have ligand, reverse-lookup by ligand
    project_root = Path(__file__).resolve().parent.parent
    for g in items:
        if not g.pdb and g.ligand:
            lig_base = Path(g.ligand).stem
            pdb_guess = _reverse_lookup_pdb_by_ligand(project_root, lig_base)
            if pdb_guess:
                g.pdb = pdb_guess



# -----------------------
# Scanner
# -----------------------
def scan(
    proteins_root: Path,
    ligands_root: Path,
    rules: List[Rule],
    pdb_filter: Optional[str] = None,
    since: Optional[datetime] = None,
    bench_root: Optional[Path] = None,
    bench_logs: Optional[List[Path]] = None,
) -> Tuple[Dict[str, GroupItem], Dict[str, int], Dict[str, List[str]], bool]:
    groups: Dict[str, GroupItem] = {}
    sev_counts = defaultdict(int)
    scanned = {"proteins": [], "ligand_tsvs": [], "bench_logs": []}
    saw_timed_since = False

    def make_key(rule: Rule, pdb: str, ligand: Optional[str], file_base: str) -> str:
        return f"{rule.id}|{pdb}|{ligand or ''}|{file_base}"

    # Walk PDBs under proteins_root
    for pdb_dir in sorted(proteins_root.iterdir()):
        if not pdb_dir.is_dir():
            continue
        pdb = pdb_dir.name
        if pdb_filter and pdb.upper() != pdb_filter.upper():
            continue

        protein_log = pdb_dir / "protein.log"
        lig_tsv = ligands_root / pdb / "ligand_prep_status.tsv"

        scanned["proteins"].append(str(protein_log))
        scanned["ligand_tsvs"].append(str(lig_tsv))

        lig_index = _load_ligand_tsv(lig_tsv)

        # 1) Apply regex rules to protein.log
        for ts, line in _iter_lines(protein_log):
            if since and ts and ts < since:
                continue
            if ts and (not since or ts >= since):
                saw_timed_since = True
            for rule in rules:
                m = rule.regex.search(line)
                if not m:
                    continue
                ligand = Path(m.group("lig")).name if "lig" in m.groupdict() else None
                key = make_key(rule, pdb, ligand, protein_log.name)
                if key not in groups:
                    groups[key] = GroupItem(
                        rule_id=rule.id,
                        severity=rule.severity,
                        hint=rule.hint,
                        pdb=pdb,
                        ligand=ligand,
                        file_basename=protein_log.name,
                        priority=rule.priority,
                        step=RULE_TO_STEP.get(rule.id, ""),
                    )

                groups[key].add_hit(str(protein_log), line, ts)
                sev_counts[rule.severity] += 1

        # 2) TSV hard failures
        for lig_key, payload in lig_index.items():
            status = payload.get("prep_status", "")
            if status and status not in ("ok", "success", "passed", "done"):
                rule_id = "ligprep.tsv_failure"
                key = f"{rule_id}|{pdb}|{lig_key}|ligand_prep_status.tsv"
                if key not in groups:
                    groups[key] = GroupItem(
                        rule_id=rule_id,
                        severity="hard",
                        hint="Treat TSV as ground truth; review the 'reason' column and quarantine failing ligand.",
                        pdb=pdb,
                        ligand=lig_key,
                        file_basename="ligand_prep_status.tsv",
                        priority=95,
                        step=RULE_TO_STEP.get(rule_id, ""),
                    )
                line = f"TSV failure for {lig_key}: status={payload.get('prep_status')} reason={payload.get('reason')}"
                groups[key].add_hit(str(lig_tsv), line, None)
                sev_counts["hard"] += 1

        # 3) Correlate TSV with no_score
        for g in list(groups.values()):
            if g.pdb != pdb or g.rule_id != "dock.no_score" or not g.ligand:
                continue
            base = Path(g.ligand).name
            base_noext = Path(base).stem
            joined = lig_index.get(base) or lig_index.get(base_noext)
            if joined:
                g.joins.append({"ligand": base, **joined})

    # Optional bench directory
    if bench_root and bench_root.exists():
        for logf in sorted(bench_root.glob("*.log")):
            scanned["bench_logs"].append(str(logf))
            for ts, line in _iter_lines(logf):
                if since and ts and ts < since:
                    continue
                if ts and (not since or ts >= since):
                    saw_timed_since = True
                for rule in rules:
                    m = rule.regex.search(line)
                    if not m:
                        continue
                    pdb = _extract_pdb_from_line(line) or _extract_pdb_from_path(logf)
                    ligand = Path(m.group("lig")).name if "lig" in m.groupdict() else None
                    key = f"{rule.id}|{pdb or ''}|{ligand or ''}|{logf.name}"
                    if key not in groups:
                        groups[key] = GroupItem(
                            rule_id=rule.id,
                            severity=rule.severity,
                            hint=rule.hint,
                            pdb=pdb or "",
                            ligand=ligand,
                            file_basename=logf.name,  # or Path(bl).name
                            priority=rule.priority,
                            step=RULE_TO_STEP.get(rule.id, ""),
                        )
                    groups[key].add_hit(str(logf), line, ts)
                    sev_counts[rule.severity] += 1

    # Explicit bench logs
    for bl in (bench_logs or []):
        if not bl or not Path(bl).exists():
            continue
        scanned["bench_logs"].append(str(bl))
        for ts, line in _iter_lines(Path(bl)):
            if since and ts and ts < since:
                continue
            if ts and (not since or ts >= since):
                saw_timed_since = True
            for rule in rules:
                m = rule.regex.search(line)
                if not m:
                    continue
                pdb = _extract_pdb_from_line(line) or _extract_pdb_from_path(Path(bl))
                ligand = Path(m.group("lig")).name if "lig" in m.groupdict() else None
                key = f"{rule.id}|{pdb or ''}|{ligand or ''}|{Path(bl).name}"
                if key not in groups:
                    groups[key] = GroupItem(
                        rule_id=rule.id,
                        severity=rule.severity,
                        hint=rule.hint,
                        pdb=pdb or "",
                        ligand=ligand,
                        file_basename=Path(bl).name,
                        priority=rule.priority,
                        step=RULE_TO_STEP.get(rule.id, ""),
                    )
                groups[key].add_hit(str(bl), line, ts)
                sev_counts[rule.severity] += 1

    return groups, sev_counts, scanned, saw_timed_since

def _norm(p: Path) -> str:
    return str(p.resolve()).replace("\\", "/")

def enrich_paths(items: List[GroupItem], project_root: Path,
                 proteins_root: Path, ligands_root: Path, logs_dir: Path) -> List[GroupItem]:
    """
    Populate g.paths with existing files for each item, using your tree:
      processed_pdbs/<PDB>/ligands_raw/*[.sanitized].pdb
      prepped_ligands/<PDB>/*.pdbqt  (+ ligand_prep_status.tsv)
      docked/<PDB>/bench_pocket1_single/*{ligand_base}*, docked/<PDB>/protein.log
      logs/bench_*.log (that matched)
    Also: optionally probe for <ligand_base>.mol2 within those same three roots;
    include only if found (we don't assume a canonical mol2 location).
    """
    processed_root = project_root / "processed_pdbs"
    prepped_root   = ligands_root
    docked_root    = proteins_root

    for g in items:
        paths: Dict[str, List[str]] = {}
        pdb = (g.pdb or "").strip()
        lig = (g.ligand or "").strip()
        lig_base = Path(lig).stem if lig else ""

        # processed_pdbs/<PDB>/ligands_raw/*
        if pdb:
            raw_dir = processed_root / pdb / "ligands_raw"
            if raw_dir.exists():
                pre = [p for p in raw_dir.glob(f"{lig_base}*.pdb")] if lig_base else list(raw_dir.glob("*.pdb"))
                san = [p for p in raw_dir.glob(f"{lig_base}*.sanitized.pdb")] if lig_base else [p for p in pre if str(p).endswith(".sanitized.pdb")]
                pre_raw = [_norm(p) for p in pre if not str(p).endswith(".sanitized.pdb")]
                pre_san = [_norm(p) for p in san]
                if pre_raw: paths["pre_raw_pdb"] = pre_raw
                if pre_san: paths["pre_raw_pdb_sanitized"] = pre_san

        # prepped_ligands/<PDB>/*.pdbqt (+ TSV)
        if pdb:
            prepped_dir = prepped_root / pdb
            if prepped_dir.exists():
                pdbqt_hits = list(prepped_dir.glob(f"{lig_base}*.pdbqt")) if lig_base else list(prepped_dir.glob("*.pdbqt"))
                if pdbqt_hits:
                    paths["post_pdbqt"] = [_norm(p) for p in pdbqt_hits]
                tsv = prepped_dir / "ligand_prep_status.tsv"
                if tsv.exists():
                    paths["ligand_prep_status_tsv"] = [_norm(tsv)]

                # Optional: try to find a .mol2 sibling if present
                mol2_hits = list(prepped_dir.glob(f"{lig_base}*.mol2")) if lig_base else []
                if mol2_hits:
                    paths["post_mol2"] = [_norm(p) for p in mol2_hits]

        # docked/<PDB>/bench_pocket1_single/*
        if pdb:
            dock_dir = docked_root / pdb / "bench_pocket1_single"
            if dock_dir.exists():
                dock_hits = list(dock_dir.glob(f"*{lig_base}*")) if lig_base else []
                if dock_hits:
                    paths["dock_inputs"] = [_norm(p) for p in dock_hits]
            prot_log = docked_root / pdb / "protein.log"
            if prot_log.exists():
                paths["protein_log"] = [_norm(prot_log)]

            # Optional: .mol2 might also be under docked (rare but cheap to probe)
            if lig_base and dock_dir.exists():
                mol2_hits_docked = list(dock_dir.glob(f"*{lig_base}*.mol2"))
                if mol2_hits_docked:
                    paths.setdefault("post_mol2", [])
                    paths["post_mol2"].extend(_norm(p) for p in mol2_hits_docked)

        # Bench logs that matched this group (from g.files)
            bench_matches = []
            for fp in g.files.keys():
                try:
                    p = Path(fp).resolve()
                except Exception:
                    p = Path(fp)
                if p.suffix == ".log" and p.parent.resolve() == logs_dir.resolve():
                    bench_matches.append(_norm(p))
            if bench_matches:
                paths["bench_logs"] = sorted(set(bench_matches))

            g.paths = paths
    return items


def _extract_pdb_from_line(line: str) -> Optional[str]:
    m = re.search(r"([0-9](?:[A-Za-z][A-Za-z0-9]{2}|[A-Za-z0-9][A-Za-z][A-Za-z0-9]|[A-Za-z0-9]{2}[A-Za-z]))", line)
    return m.group(1).upper() if m else None


def _extract_pdb_from_path(p: Path) -> Optional[str]:
    m = re.search(r"([0-9](?:[A-Za-z][A-Za-z0-9]{2}|[A-Za-z0-9][A-Za-z][A-Za-z0-9]|[A-Za-z0-9]{2}[A-Za-z]))",p.stem)
    return m.group(1).upper() if m else None

def _reverse_lookup_pdb_by_ligand(project_root: Path, lig_base: str) -> Optional[str]:
    # Try prepped_ligands/<PDB>/<lig_base>*
    prepped = project_root / "prepped_ligands"
    for pdir in prepped.iterdir() if prepped.exists() else []:
        if not pdir.is_dir():
            continue
        if list(pdir.glob(f"{lig_base}*")):
            return pdir.name.upper()
    # Try processed_pdbs/<PDB>/ligands_raw/<lig_base>*
    processed = project_root / "processed_pdbs"
    for pdir in processed.iterdir() if processed.exists() else []:
        rawdir = pdir / "ligands_raw"
        if rawdir.exists() and list(rawdir.glob(f"{lig_base}*")):
            return pdir.name.upper()
    return None

# -----------------------
# Ranking and reporting
# -----------------------
def rank_groups(groups: Dict[str, GroupItem]) -> List[GroupItem]:
    items = list(groups.values())
    def sort_key(g: GroupItem):
        sev = SEVERITY_ORDER.get(g.severity, -1)
        last = g.last_seen.timestamp() if g.last_seen else 0
        return (-sev, -g.count, -last, -g.priority)
    items.sort(key=sort_key)
    for i, g in enumerate(items, 1):
        setattr(g, "rank", i)
    return items


def compute_health(stats: Dict[str, int]) -> int:
    score = 100
    score -= stats.get("fatal", 0) * 30
    score -= stats.get("hard", 0) * 10
    score -= stats.get("soft", 0) * 2
    return max(0, min(100, score))


# -----------------------
# Writers (render + write)
# -----------------------
def compose_json(items: List[GroupItem], stats: Dict[str, int], scanned: Dict[str, List[str]], config: Dict[str, str]):
    return {
        "schema_version": 1,
        "generator": "atlas.log_resurfacer.v2",
        "config": config,
        "scanned": scanned,
        "stats": stats,
        "items": [
            {
                "rank": getattr(g, "rank", None),
                "severity": g.severity,
                "rule_id": g.rule_id,
                "pdb": g.pdb,
                "ligand": g.ligand,
                "message": _group_message(g),
                "hint": g.hint,
                "count": g.count,
                "first_seen": g.first_seen.strftime("%Y-%m-%d %H:%M:%S") if g.first_seen else None,
                "last_seen": g.last_seen.strftime("%Y-%m-%d %H:%M:%S") if g.last_seen else None,
                "files": [{"path": fs.path, "samples": fs.samples} for fs in g.files.values()],
                "joins": g.joins,
                "step": getattr(g, "step", ""),
                "paths": getattr(g, "paths", {}),
            }
            for g in items
        ],
    }


def write_json(outdir: Path, payload: Dict):
    (outdir / "resurfaced.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")


def render_tsv(items: List[GroupItem]) -> str:
    lines = ["\t".join([
        "rank", "severity", "rule_id", "pdb", "ligand", "message", "count", "recency", "file", "sample", "hint", "step",
        "paths"
    ])]
    for g in items:
        for fs in g.files.values():
            sample = fs.samples[0] if fs.samples else ""
            paths_str = ";".join(sum([v for v in getattr(g, "paths", {}).values()], []))
            row = [
                str(getattr(g, "rank", "")), g.severity, g.rule_id, g.pdb, g.ligand or "",
                _group_message(g), str(g.count), g.recency, fs.path, sample, g.hint,
                getattr(g, "step", ""), paths_str
            ]
            lines.append("\t".join(row))
    return "\n".join(lines) + "\n"


def write_tsv(outdir: Path, tsv_text: str):
    (outdir / "resurfaced.tsv").write_text(tsv_text, encoding="utf-8")


def render_markdown(items: List[GroupItem], stats: Dict[str, int], stale_note: Optional[str] = None) -> str:
    md: List[str] = []
    md.append("# Atlas Log Resurfacer - Summary\n\n")
    if stale_note:
        md.append(f"> Note: {stale_note}\n\n")
    md.append(f"*Fatal:* {stats.get('fatal',0)}  *Hard:* {stats.get('hard',0)}  *Soft:* {stats.get('soft',0)}  *Info:* {stats.get('info',0)}\n\n")
    md.append(f"**Health score:** {compute_health(stats)} / 100\n\n")
    for g in items:
        md.append(f"## {getattr(g,'rank', '?')}. [{g.severity.upper()}] {g.rule_id} - {g.pdb}{' - ' + (g.ligand or '') if g.ligand else ''}\n")
        md.append(f"- **Message:** {_group_message(g)}\n")
        md.append(f"- **Count:** {g.count}  |  **Recency:** {g.recency}\n")
        md.append(f"- **Hint:** {g.hint}\n")
        for fs in g.files.values():
            md.append(f"  - **File:** `{fs.path}`\n")
            for s in fs.samples:
                md.append(f"    - `{s}`\n")
        if g.joins:
            md.append("- **Joins:**\n")
            for j in g.joins:
                md.append(f"  - ligand={j.get('ligand')} prep_status={j.get('prep_status')} reason={j.get('reason')}\n")
        if getattr(g, "step", ""):
            md.append(f"- **Step:** {g.step}\n")
        if getattr(g, "paths", {}):
            md.append("- **Paths:**\n")
            for key, plist in g.paths.items():
                for p in plist:
                    md.append(f"  - {key}: `{p}`\n")
        md.append("\n")
    return "".join(md)


def write_markdown(outdir: Path, md_text: str):
    (outdir / "resurfaced.md").write_text(md_text, encoding="utf-8")


def _group_message(g: GroupItem) -> str:
    mapping = {
        "dock.no_score": "No score for {n} ligand(s)",
        "render.no_ctrl_rdk_paths": "Screenshots skipped: no ctrl/rdk pose paths",
        "render.enqueue_none": "Controls/RDK not enqueued",
        "recentering.skipped_zero": "Early recenter skipped: evaluated=0",
        "cap.limiting": "CAP limiting non-controls",
        "ligprep.rdkit_valence": "RDKit sanitize/valence/kekulize error(s)",
        "ligprep.obabel_h_charge": "OpenBabel H/charge warning(s)",
        "ligprep.mgl_missing_mol2": "Missing .mol2 during conversion (fatal)",
        "ligprep.mgl_partial_write": "Partial write during conversion",
        "ligprep.pipeline_bug_scope": "Pipeline bug: uninitialized variable masks prep outcome",
        "ligprep.tsv_failure": "Ligand prep failed per TSV",
    }
    templ = mapping.get(g.rule_id, g.rule_id)
    return templ.replace("{n}", str(g.count))


# -----------------------
# Console summary
# -----------------------
def print_console(items: List[GroupItem], stats: Dict[str, int], top: int):
    print("\n== Per-severity counts ==")
    for sev in ("fatal", "hard", "soft", "info"):
        print(f"{sev:>5}: {stats.get(sev,0)}")
    print(f"Health score: {compute_health(stats)} / 100\n")

    headers = ["#", "sev", "rule", "pdb", "lig", "count", "recency", "file", "sample", "step", "rep_path"]
    rows = []
    for g in items[:top]:
        fs = next(iter(g.files.values())) if g.files else FileSample(path="", samples=[""])
        sample = fs.samples[0] if fs.samples else ""
        # representative path (post_pdbqt > post_mol2 > pre_raw_* > protein_log > dock_inputs > bench_logs)
        rep = _pick_representative_path(g)
        rows.append([
            getattr(g, "rank", ""),
            g.severity,
            g.rule_id,
            g.pdb,
            g.ligand or "",
            g.count,
            g.recency,
            Path(fs.path).name if fs.path else "",
            sample,  
            getattr(g, "step", ""),
            rep,
        ])
    _print_table(headers, rows)
def print_paths_footer(items: List[GroupItem], top: int):
    print("\n[paths] top (rank • step → rep_path)")
    for g in items[:top]:
        rep = _pick_representative_path(g)
        step = getattr(g, "step", "")
        rank = getattr(g, "rank", "?")
        if rep or step:
            print(f"  {rank}. {step or '-'} → {rep or '-'}")


def _truncate(s: str, n: int) -> str:
    return s if len(s) <= n else s[: n - 1] + "..."


def _print_table(headers: List[str], rows: List[List[object]]):
    """
    Pretty, fixed-width console table:
      - clamps each column to a sensible max,
      - adapts to terminal width,
      - ellipsizes overly long cells so rows never wrap.
    """
    import shutil

    # Column order for headers:
    # 0:# 1:sev 2:rule 3:pdb 4:lig 5:count 6:recency 7:file 8:sample 9:step 10:rep_path
    # Set *target caps* for fixed columns; we'll size 'rep_path' dynamically
    caps = {
        0: 3,   # #
        1: 4,   # sev
        2: 25,  # rule
        3: 4,   # pdb
        4: 24,  # lig
        5: 5,   # count
        6: 19,  # recency
        7: 22,  # file
        8: 48,  # sample (we'll shrink first if needed)
        9: 22,  # step
        # 10 rep_path -> dynamic
    }
    min_sample = 18
    min_rep    = 24

    # Terminal width (fallback if not a TTY)
    term = shutil.get_terminal_size((140, 24)).columns

    # Compute static sum (all except rep_path) + spaces
    static_sum = sum(caps.values())
    gaps = len(headers) - 1
    # Desired width for rep_path with a decent minimum
    desired_rep = max(min_rep, term - static_sum - gaps)

    # If negative space, shrink sample first, then rep_path down to minima
    sample_w = caps[8]
    rep_w = desired_rep
    overflow = (static_sum + gaps + desired_rep) - term
    if overflow > 0:
        take = min(overflow, sample_w - min_sample)
        sample_w -= take
        overflow -= take
    if overflow > 0:
        rep_w = max(min_rep, desired_rep - overflow)
        overflow = 0

    # Final widths array (with clamps) and alignments
    widths = [0] * len(headers)
    aligns = ["<"] * len(headers)  # left by default
    numeric_right = {0, 5}         # #, count right align
    for i in numeric_right:
        aligns[i] = ">"
    # recency right-ish looks nicer in narrow view
    aligns[6] = ">"
    # file, sample, rep_path left
    widths[10] = rep_w
    for i in range(len(headers)):
        if i == 8:
            widths[i] = sample_w
        elif i in caps:
            widths[i] = caps[i]
        # else: already set for 10

    def _ellipsize(s: str, w: int) -> str:
        s = "" if s is None else str(s).replace("\n", " ").replace("\t", " ")
        if len(s) <= w:
            return s
        if w <= 3:
            return s[:w]
        # center-ellipsis for long paths and names
        keep = w - 3
        head = keep // 2
        tail = keep - head
        return f"{s[:head]}...{s[-tail:]}"

    # Print header
    fmt = " ".join([f"{{:{aligns[i]}{widths[i]}}}" for i in range(len(headers))])
    print(fmt.format(*[ _ellipsize(h, widths[i]) for i, h in enumerate(headers) ]))
    print(" ".join("-" * widths[i] for i in range(len(headers))))

    # Print rows, clamped
    for r in rows:
        cells = [ _ellipsize(r[i], widths[i]) for i in range(len(headers)) ]
        print(fmt.format(*cells))



# -----------------------
# Runner (benchmark + resurface)
# -----------------------
# --- at top of file if not present ---
from pathlib import Path
import os
import subprocess

def tee_run(cmd, log_path, cwd=None, env=None):
    """
    Run a subprocess, stream stdout to console AND to log_path.
    Robust to non-UTF-8 bytes.
    """
    # Normalize and create parent dir
    log_path = str(log_path)
    os.makedirs(os.path.dirname(log_path), exist_ok=True)

    print(f"[tee] opening log file: {log_path}")
    with open(log_path, "a", encoding="utf-8", errors="replace") as flog:
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd) if cwd is not None else None,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
        for line in proc.stdout:
            print(line, end="")
            flog.write(line)
        rc = proc.wait()
    return rc

# -----------------------
# CLI
# -----------------------
def main(argv: Optional[List[str]] = None) -> int:
    # Split on passthrough marker '--'
    argv = sys.argv[1:] if argv is None else argv
    if "--" in argv:
        idx = argv.index("--")
        own_args = argv[:idx]
        bench_passthrough = argv[idx + 1 :]
    else:
        own_args, bench_passthrough = argv, []

    ap = argparse.ArgumentParser(description="Atlas Log Resurfacer (single entry)", add_help=True)
    ap.add_argument("--no-benchmark", action="store_true", help="skip launching benchmark_mode.py; only resurface")
    ap.add_argument("--bench-log", action="append", help="path to a bench log to include; can repeat")
    ap.add_argument("--proteins", type=Path, default=Path("./docked"), help="root of docked/<PDB>/protein.log")
    ap.add_argument("--ligands", type=Path, default=Path("./prepped_ligands"), help="root of prepped_ligands/<PDB>/ligand_prep_status.tsv")
    ap.add_argument("--include-bench", type=Path, help="optional root logs directory to include (*.log)")
    ap.add_argument("--rules", type=Path, help="optional JSON rule file to override defaults")
    ap.add_argument("--out-root", type=Path, default=Path(__file__).parent / "out", help="output root directory")
    ap.add_argument("--pdb", help="limit to a single PDB (e.g., 1T46)")
    ap.add_argument("--since", help="timestamp floor 'YYYY-MM-DD HH:MM:SS' (defaults to benchmark start)")
    ap.add_argument("--top", type=int, default=30, help="top N issues to print to console")
    ap.add_argument("--markdown", action="store_true", help="write Markdown summary")
    ap.add_argument("--tsv", action="store_true", help="write TSV summary")
    ap.add_argument("--fail-on", default="hard", choices=list(SEVERITY_ORDER.keys()), help="exit 1 if >= severity present")
    ap.add_argument("--force-resurface", action="store_true", help="if benchmark fails, still exit based on resurfacer findings")

    args = ap.parse_args(own_args)

    # Resolve paths
    project_root = Path(__file__).resolve().parent.parent  # one up from log_resurfacer/
    logs_dir = project_root / "logs"
    bench_script = project_root / "benchmark_mode.py"

    # Establish timestamp for this session
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    bench_log_created: Optional[Path] = None

    # Launch benchmark unless requested not to
    bench_rc = 0
    if not args.no_benchmark:
        logs_dir.mkdir(parents=True, exist_ok=True)
        bench_log_created = logs_dir / f"bench_{ts}.log"
        cmd = [sys.executable, str(bench_script)] + bench_passthrough
        print(f"[launcher] running: {' '.join(cmd)}\n[launcher] tee -> {bench_log_created}")
        bench_rc = tee_run(cmd, bench_log_created, cwd=project_root)
        print(f"[launcher] benchmark exited with rc={bench_rc}")

    # Determine since floor
    since_str = args.since or None
    if since_str is None:
        # Default to benchmark-start ts if we ran it; else now
        since_str = datetime.strptime(ts, "%Y%m%d_%H%M%S").strftime("%Y-%m-%d %H:%M:%S")
    try:
        since_dt = datetime.fromisoformat(since_str)
    except Exception:
        since_dt = None

    # Bench logs list (explicit + created)
    bench_logs: List[Path] = []
    if bench_log_created:
        bench_logs.append(bench_log_created)
    for bl in (args.bench_log or []):
        bench_logs.append(Path(bl))

    # Parse rules
    rules = _parse_rules(str(args.rules) if args.rules else None)

    # Scan & rank
    groups, sev_counts, scanned, saw_timed_since = scan(
        proteins_root=args.proteins,
        ligands_root=args.ligands,
        rules=rules,
        pdb_filter=args.pdb,
        since=since_dt,
        bench_root=args.include_bench,
        bench_logs=bench_logs,
    )
    items = rank_groups(groups)
    # Enrich each item with discovered filesystem paths
    infer_missing_ids(items)

    items = enrich_paths(
        items=items,
        project_root=project_root,
        proteins_root=args.proteins,
        ligands_root=args.ligands,
        logs_dir=logs_dir,
    )
    # Output directory (use same ts so bench & reports align)
    outdir = (Path(args.out_root) if args.out_root else Path(__file__).parent / "out") / ts
    outdir.mkdir(parents=True, exist_ok=True)

    # Config echo for JSON
    config_echo = {
        "fail_on": args.fail_on,
        "top": args.top,
        "since": since_str,
        "rules_path": str(args.rules) if args.rules else "<embedded>",
        "pdb_filter": args.pdb or "",
    }

    # Reports
    payload = compose_json(items, sev_counts, scanned, config_echo)
    tsv_text = render_tsv(items)
    stale_note = None if (since_dt and saw_timed_since) else "No timestamped lines >= --since were found; results may include older untimed lines."
    md_text = render_markdown(items, sev_counts, stale_note)

    write_json(outdir, payload)
    if args.tsv or (not args.no_benchmark):
        write_tsv(outdir, tsv_text)
    if args.markdown or (not args.no_benchmark):
        write_markdown(outdir, md_text)

    # Console
    print_console(items, sev_counts, args.top)
    print_paths_footer(items, args.top) 
    if bench_log_created:
        print(f"\n[paths] bench log: {bench_log_created}")
    print(f"[paths] reports: {outdir}")

    # Exit code policy
    fail_present = any(_severity_at_or_above(sev, args.fail_on) and sev_counts.get(sev, 0) > 0 for sev in SEVERITY_ORDER)

    if not args.no_benchmark and bench_rc not in (0, None):
        if args.force_resurface:
            return 1 if fail_present else 2
        return 2

    return 1 if fail_present else 0


if __name__ == "__main__":
    sys.exit(main())
