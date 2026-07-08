"""Small fixed-column PDB/PDBQT record helpers."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Iterable

PDB_ATOM_PREFIXES = ("ATOM", "HETATM")
_UNIPROT_ACCESSION_RE = re.compile(
    r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?\b"
)


def is_atom_record(line: str) -> bool:
    return line.startswith(PDB_ATOM_PREFIXES)


def record_type(line: str) -> str:
    return line[:6].strip().upper()


def padded_record(line: str) -> str:
    return line.rstrip("\n").ljust(80)


def atom_name(line: str) -> str:
    return padded_record(line)[12:16].strip()


def residue_name(line: str) -> str:
    return padded_record(line)[17:20].strip().upper()


def chain_id(line: str, default: str = "-") -> str:
    return padded_record(line)[21:22].strip() or default


def resseq(line: str, default: str = "0") -> str:
    return padded_record(line)[22:26].strip() or default


def insertion_code(line: str, default: str = "-") -> str:
    return padded_record(line)[26:27].strip() or default


def line_element(line: str) -> str:
    padded = padded_record(line)
    element = padded[76:78].strip().upper()
    if element:
        return element
    atom_name = padded[12:16].strip().upper()
    return "".join(ch for ch in atom_name[:2] if ch.isalpha()).upper()


def line_serial(line: str) -> int | None:
    text = line[6:11].strip()
    return int(text) if text.isdigit() else None


def line_xyz(line: str) -> tuple[float, float, float] | None:
    if not is_atom_record(line):
        return None
    try:
        return float(line[30:38]), float(line[38:46]), float(line[46:54])
    except (TypeError, ValueError):
        return None


def pdbqt_line_xyz(line: str) -> tuple[float, float, float] | None:
    fixed = line_xyz(line)
    if fixed is not None:
        return fixed
    if not is_atom_record(line):
        return None
    parts = line.split()
    if len(parts) < 8:
        return None
    try:
        return float(parts[5]), float(parts[6]), float(parts[7])
    except (TypeError, ValueError):
        return None


def pdbqt_line_element(line: str) -> str:
    fields = line.split()
    if fields:
        token = re.sub(r"[^A-Za-z]", "", fields[-1]).upper()
        if token:
            if token.startswith("CL"):
                return "Cl"
            if token.startswith("BR"):
                return "Br"
            if token.startswith("H"):
                return "H"
            return token[0].upper()
    name = re.sub(r"[^A-Za-z]", "", atom_name(line)).upper()
    return name[0].upper() if name else ""


def residue_key(line: str) -> tuple[str, str, str, str]:
    return (
        residue_name(line),
        chain_id(line),
        resseq(line, default=""),
        insertion_code(line, default=""),
    )


def residue_key_text(line: str) -> str:
    _, chain, resseq, icode = residue_key(line)
    return f"{chain}:{resseq or '0'}:{icode or '-'}"


def atom_identity(line: str) -> tuple[str, str, str, str, str, str]:
    return (
        atom_name(line),
        padded_record(line)[16:17].strip(),
        residue_name(line),
        chain_id(line, default=""),
        resseq(line, default=""),
        insertion_code(line, default=""),
    )


def line_resseq_int(line: str) -> int:
    try:
        return int(padded_record(line)[22:26].strip())
    except (TypeError, ValueError):
        return 0


def centroid(coords: Iterable[tuple[float, float, float]]) -> tuple[float, float, float] | None:
    sx = sy = sz = 0.0
    count = 0
    for x, y, z in coords:
        sx += float(x)
        sy += float(y)
        sz += float(z)
        count += 1
    if count == 0:
        return None
    return sx / count, sy / count, sz / count


def pdbqt_model_blocks(lines: Iterable[str]) -> list[list[str]]:
    blocks: list[list[str]] = []
    cur: list[str] = []
    in_model = False
    for line in lines:
        if line.startswith("MODEL"):
            if cur:
                blocks.append(cur)
                cur = []
            in_model = True
            cur.append(line)
        elif line.startswith("ENDMDL"):
            cur.append(line)
            blocks.append(cur)
            cur = []
            in_model = False
        else:
            cur.append(line)
    if in_model:
        blocks.append(cur)
    elif not blocks:
        blocks = [cur]
    return blocks


def pdbqt_model_blocks_from_path(path: str | Path) -> list[list[str]]:
    with Path(path).open("r", encoding="utf-8", errors="ignore") as handle:
        return pdbqt_model_blocks(handle.readlines())


def local_uniprots_from_pdb(repo_root: Path, pdb_id: str) -> list[str]:
    pdb_path = Path(repo_root) / "input_pdbs" / f"{str(pdb_id).strip().upper()}.pdb"
    if not pdb_path.exists():
        return []
    hits: set[str] = set()
    try:
        with pdb_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                rec = record_type(line)
                if rec in {"ATOM", "HETATM", "MODEL"}:
                    break
                if rec not in {"HEADER", "TITLE", "COMPND", "SOURCE", "DBREF", "DBSOURCE", "REMARK"}:
                    continue
                for match in _UNIPROT_ACCESSION_RE.finditer(line):
                    hits.add(match.group(0).split("-", 1)[0].upper())
                token_match = re.search(
                    r"\bUNP(?:ROTK?)?\s*[:=]?\s*([A-Z0-9\-]+)", line, re.IGNORECASE
                )
                if token_match:
                    token = token_match.group(1).strip().upper()
                    if _UNIPROT_ACCESSION_RE.fullmatch(token):
                        hits.add(token.split("-", 1)[0].upper())
    except OSError:
        return []
    return sorted(hits)
