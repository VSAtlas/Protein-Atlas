from __future__ import annotations

import csv
import hashlib
import math
import re
import shutil
import urllib.request
import zipfile
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
from xml.etree import ElementTree

from post_docking.mmgbsa._atomic_io import (
    tmp_path as _tmp_path,
    write_csv_rows_atomic,
    write_json_atomic,
)

WANG2016_SUPPLEMENT_URL = (
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/5018451/bin/"
    "NIHMS806653-supplement-Supp_info.docx"
)
WANG2016_SUPPLEMENT_URLS = (
    WANG2016_SUPPLEMENT_URL,
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC5018451/bin/"
    "NIHMS806653-supplement-Supp_info.docx",
)

WANG2016_TABLE1 = (
    ("CDK+PKA", "INP=1", 5.81, 0.81, 0.18),
    ("CDK+PKA", "INP=2", 3.42, 0.71, 0.29),
    ("glucosidase", "INP=1", 4.43, 0.57, 0.16),
    ("glucosidase", "INP=2", 1.94, 0.30, 0.23),
    ("thrombin", "INP=1", 11.02, 0.87, 0.23),
    ("thrombin", "INP=2", 4.26, 0.86, 0.42),
    ("trypsin", "INP=1", 9.89, 0.74, 0.18),
    ("trypsin", "INP=2", 3.30, 0.53, 0.32),
    ("urokinase", "INP=1", 7.19, 0.73, 0.22),
    ("urokinase", "INP=2", 2.61, 0.53, 0.33),
    ("Factor Xa", "INP=1", 5.99, 0.75, 0.23),
    ("Factor Xa", "INP=2", 3.26, 0.67, 0.27),
)

WANG2016_TABLE8 = (
    ("glucosidase", "INP=1", 9.73),
    ("glucosidase", "INP=2", 8.55),
    ("thrombin", "INP=1", 23.90),
    ("thrombin", "INP=2", 4.85),
    ("trypsin", "INP=1", 18.50),
    ("trypsin", "INP=2", 6.16),
    ("urokinase", "INP=1", 20.70),
    ("urokinase", "INP=2", 4.70),
    ("Factor Xa", "INP=1", 29.48),
    ("Factor Xa", "INP=2", 6.69),
)

WANG2016_TABLE3_REPRESENTATIVE_PB = (
    ("1Q8W", "CDK+PKA", -8.22, 38.22, 38.05, 37.66),
    ("2J75", "glucosidase", -10.01, 18.45, 18.16, 17.15),
    ("1VZQ", "thrombin", -17.52, 3.75, 3.36, 2.26),
    ("1O2Q", "trypsin", -18.87, 19.87, 19.67, 18.94),
    ("1C5Z", "urokinase", -13.23, 38.41, 38.32, 37.97),
    ("1LPZ", "Factor Xa", -20.19, 12.93, 12.77, 12.31),
)

WANG2016_CONVERGENCE_10NS_VS_1NS = (
    ("CDK+PKA", "10ns", 3.65, 0.72, 0.30),
    ("CDK+PKA", "1ns", 3.51, 0.73, 0.29),
    ("glucosidase", "10ns", 1.77, 0.47, 0.38),
    ("glucosidase", "1ns", 1.70, 0.43, 0.34),
    ("thrombin", "10ns", 4.58, 0.86, 0.40),
    ("thrombin", "1ns", 4.46, 0.88, 0.41),
    ("trypsin", "10ns", 3.10, 0.67, 0.37),
    ("trypsin", "1ns", 3.25, 0.66, 0.36),
    ("urokinase", "10ns", 2.23, 0.67, 0.47),
    ("urokinase", "1ns", 2.17, 0.60, 0.43),
    ("Factor Xa", "10ns", 3.00, 0.70, 0.30),
    ("Factor Xa", "1ns", 2.85, 0.71, 0.31),
)

WANG2016_REPRESENTATIVE_COMPLEXES = (
    ("1Q8W", "CDK+PKA", "representative PB-energy/convergence complex"),
    ("2J75", "glucosidase", "representative PB-energy/convergence complex"),
    ("1VZQ", "thrombin", "representative PB-energy/convergence complex"),
    ("1O2Q", "trypsin", "representative PB-energy/convergence complex"),
    ("1C5Z", "urokinase", "representative PB-energy/convergence complex"),
    ("1LPZ", "Factor Xa", "representative PB-energy/convergence complex"),
)

_PMC_POW_COOKIE = "cloudpmc-viewer-pow"


def _pmc_pow_cookie(data: bytes) -> str:
    text = data.decode("utf-8", errors="ignore")
    challenge = re.search(r'POW_CHALLENGE\s*=\s*"([^"]+)"', text)
    difficulty = re.search(r'POW_DIFFICULTY\s*=\s*"([^"]+)"', text)
    cookie_name = re.search(r'POW_COOKIE_NAME\s*=\s*"([^"]+)"', text)
    if not challenge or not difficulty:
        return ""
    try:
        leading_zeros = max(1, int(difficulty.group(1)))
    except ValueError:
        leading_zeros = 4
    target = "0" * leading_zeros
    challenge_text = challenge.group(1)
    for nonce in range(10_000_000):
        digest = hashlib.sha256(f"{challenge_text}{nonce}".encode()).hexdigest()
        if digest.startswith(target):
            name = cookie_name.group(1) if cookie_name else _PMC_POW_COOKIE
            return f"{name}={challenge_text},{nonce}"
    return ""


def _download_url(url: str, *, cookie: str = "") -> bytes:
    headers = {"User-Agent": "Atlas-MMGBSA-benchmark/1.0"}
    if cookie:
        headers["Cookie"] = cookie
    request = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    write_json_atomic(path, payload)


def _write_rows(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    write_csv_rows_atomic(path, rows)


def download_supplement(out_dir: str | Path, url: str = WANG2016_SUPPLEMENT_URL) -> Path:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    dest = out_path / "wang2016_supplement.docx"
    data = _download_url(url)
    if not data.startswith(b"PK"):
        cookie = _pmc_pow_cookie(data)
        if cookie:
            data = _download_url(url, cookie=cookie)
    if len(data) < 1024:
        raise RuntimeError(f"downloaded supplement is unexpectedly small: {len(data)} bytes")
    if not data.startswith(b"PK"):
        blocked_path = out_path / "wang2016_supplement_download_response.html"
        blocked_path.write_bytes(data)
        raise RuntimeError(
            "downloaded supplement is not a DOCX/ZIP file; "
            f"saved response to {blocked_path}"
        )
    tmp_path = _tmp_path(dest)
    tmp_path.write_bytes(data)
    tmp_path.replace(dest)
    return dest


def download_pdb(pdb_id: str, out_dir: str | Path, timeout: int = 60) -> dict[str, Any]:
    pdb = str(pdb_id).strip().upper()
    if not re.fullmatch(r"[0-9][A-Z0-9]{3}", pdb):
        raise ValueError(f"invalid PDB ID: {pdb_id}")
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    dest = out_path / f"{pdb}.pdb"
    if dest.exists() and dest.stat().st_size > 0:
        return {"pdb_id": pdb, "path": str(dest), "ok": True, "skipped": True}
    url = f"https://files.rcsb.org/download/{pdb}.pdb"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "Atlas-MMGBSA-benchmark/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            data = response.read()
        if len(data) < 1024 or b"ATOM" not in data:
            raise RuntimeError(f"unexpected PDB payload size={len(data)}")
        tmp_path = _tmp_path(dest)
        tmp_path.write_bytes(data)
        tmp_path.replace(dest)
        return {"pdb_id": pdb, "path": str(dest), "ok": True, "url": url}
    except Exception as exc:
        return {"pdb_id": pdb, "path": str(dest), "ok": False, "url": url, "error": str(exc)}


def _download_first_supplement(out_dir: Path, urls: Sequence[str]) -> tuple[Path | None, str]:
    errors: list[str] = []
    for url in urls:
        try:
            return download_supplement(out_dir, url), ""
        except Exception as exc:
            errors.append(f"{url}: {exc}")
    return None, " | ".join(errors)


def _cell_text(cell: ElementTree.Element, namespace: dict[str, str]) -> str:
    values = [
        node.text or ""
        for node in cell.findall(".//w:t", namespace)
        if node.text is not None
    ]
    return " ".join(part.strip() for part in values if part.strip()).strip()


def parse_docx_tables(docx_path: str | Path) -> list[list[list[str]]]:
    path = Path(docx_path)
    namespace = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    with zipfile.ZipFile(path) as archive:
        xml_payload = archive.read("word/document.xml")
    root = ElementTree.fromstring(xml_payload)
    tables: list[list[list[str]]] = []
    for table in root.findall(".//w:tbl", namespace):
        rows: list[list[str]] = []
        for tr in table.findall("./w:tr", namespace):
            cells = [_cell_text(tc, namespace) for tc in tr.findall("./w:tc", namespace)]
            if any(cells):
                rows.append(cells)
        if rows:
            tables.append(rows)
    return tables


def _slug(text: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(text).strip().lower()).strip("_")


def _table_rows(table: Sequence[Sequence[str]], table_index: int) -> list[dict[str, str]]:
    if not table:
        return []
    header = [_slug(cell) or f"col_{idx + 1}" for idx, cell in enumerate(table[0])]
    rows: list[dict[str, str]] = []
    for raw in table[1:]:
        row = {"table_index": str(table_index)}
        for idx, value in enumerate(raw):
            key = header[idx] if idx < len(header) else f"col_{idx + 1}"
            row[key] = value
        rows.append(row)
    return rows


def _all_table_rows(tables: Sequence[Sequence[Sequence[str]]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for idx, table in enumerate(tables, start=1):
        rows.extend(_table_rows(table, idx))
    return rows


def _looks_like_pdb(value: object) -> bool:
    return bool(re.fullmatch(r"[0-9][A-Za-z0-9]{3}", str(value).strip()))


def _candidate_reference_rows(rows: Sequence[Mapping[str, str]]) -> list[dict[str, str]]:
    candidates: list[dict[str, str]] = []
    for row in rows:
        values = list(row.values())
        if not any(_looks_like_pdb(value) for value in values):
            continue
        text = " ".join(str(value).lower() for value in values)
        if "dg" in text or "gbind" in text or "mmpbsa" in text or "kcal" in text:
            candidates.append(dict(row))
    return candidates


def write_static_targets(out_dir: str | Path) -> dict[str, str]:
    out_path = Path(out_dir)
    rel_rows = [
        {
            "benchmark_kind": "relative_binding_affinity",
            "family": family,
            "setting": setting,
            "rmsd_kcal_mol": rmsd,
            "pearson_r": r_value,
            "slope": slope,
        }
        for family, setting, rmsd, r_value, slope in WANG2016_TABLE1
    ]
    abs_rows = [
        {
            "benchmark_kind": "absolute_binding_affinity",
            "family": family,
            "setting": setting,
            "rmsd_kcal_mol": rmsd,
        }
        for family, setting, rmsd in WANG2016_TABLE8
    ]
    representative_rows = [
        {
            "pdb_id": pdb_id,
            "family": family,
            "experimental_delta_g_kcal_mol": experimental,
            "pb_delta_g_space_0_25": pb_025,
            "pb_delta_g_space_0_50": pb_050,
            "pb_delta_g_space_1_00": pb_100,
            "note": "Representative PB electrostatic table row, not full per-complex MMGBSA.",
        }
        for pdb_id, family, experimental, pb_025, pb_050, pb_100 in WANG2016_TABLE3_REPRESENTATIVE_PB
    ]
    convergence_rows = [
        {
            "family": family,
            "trajectory_window": window,
            "rmsd_kcal_mol": rmsd,
            "pearson_r": r_value,
            "slope": slope,
        }
        for family, window, rmsd, r_value, slope in WANG2016_CONVERGENCE_10NS_VS_1NS
    ]
    relative_path = out_path / "wang2016_table1_relative_targets.csv"
    absolute_path = out_path / "wang2016_table8_absolute_targets.csv"
    representative_path = out_path / "wang2016_table3_representative_pb.csv"
    convergence_path = out_path / "wang2016_10ns_vs_1ns_targets.csv"
    _write_rows(relative_path, rel_rows)
    _write_rows(absolute_path, abs_rows)
    _write_rows(representative_path, representative_rows)
    _write_rows(convergence_path, convergence_rows)
    return {
        "relative_targets": str(relative_path),
        "absolute_targets": str(absolute_path),
        "representative_pb_targets": str(representative_path),
        "convergence_targets": str(convergence_path),
    }


def write_representative_complexes(out_dir: str | Path) -> str:
    rows = [
        {
            "pdb_id": pdb_id,
            "family": family,
            "role": role,
            "benchmark_mode": "crystal_pose_and_control_redock",
        }
        for pdb_id, family, role in WANG2016_REPRESENTATIVE_COMPLEXES
    ]
    path = Path(out_dir) / "wang2016_representative_complexes.csv"
    _write_rows(path, rows)
    return str(path)


def download_representative_pdbs(out_dir: str | Path) -> dict[str, Any]:
    root = Path(out_dir)
    input_dir = root / "input_pdbs"
    rows = []
    downloads = []
    for pdb_id, family, role in WANG2016_REPRESENTATIVE_COMPLEXES:
        result = download_pdb(pdb_id, input_dir)
        downloads.append(result)
        rows.append(
            {
                "pdb_id": pdb_id,
                "family": family,
                "role": role,
                "input_pdb": result.get("path", ""),
                "download_ok": result.get("ok", False),
                "download_error": result.get("error", ""),
            }
        )
    manifest_csv = root / "wang2016_downloaded_pdbs.csv"
    _write_rows(manifest_csv, rows)
    return {
        "input_dir": str(input_dir),
        "manifest_csv": str(manifest_csv),
        "downloads": downloads,
        "ok": all(bool(item.get("ok")) for item in downloads),
    }


def _pdbbind_ligand_names(pdb_id: str) -> tuple[str, str]:
    pdb = str(pdb_id).strip().lower()
    return f"{pdb}_ligand.sdf", f"{pdb}_ligand.mol2"


def _pdbbind_root_ligand_paths(
    root: Path,
    pdb_id: str,
) -> tuple[Path | None, Path | None]:
    pdb = str(pdb_id).strip().lower()
    sdf_name, mol2_name = _pdbbind_ligand_names(pdb)
    candidate_dirs = (
        root / pdb,
        root / "pbpp-2020" / pdb,
        root / pdb.upper(),
        root / "pbpp-2020" / pdb.upper(),
    )
    sdf_path = None
    mol2_path = None
    for directory in candidate_dirs:
        if sdf_path is None and (directory / sdf_name).is_file():
            sdf_path = directory / sdf_name
        if mol2_path is None and (directory / mol2_name).is_file():
            mol2_path = directory / mol2_name
    return sdf_path, mol2_path


def _copy_pdbbind_ligands_from_root(
    *,
    pdb_id: str,
    root: Path,
    dest_dir: Path,
) -> dict[str, Any]:
    sdf_path, mol2_path = _pdbbind_root_ligand_paths(root, pdb_id)
    copied_sdf = ""
    copied_mol2 = ""
    dest_dir.mkdir(parents=True, exist_ok=True)
    if sdf_path is not None:
        copied_sdf = str(dest_dir / f"{pdb_id.upper()}_ligand.sdf")
        shutil.copy2(sdf_path, copied_sdf)
    if mol2_path is not None:
        copied_mol2 = str(dest_dir / f"{pdb_id.upper()}_ligand.mol2")
        shutil.copy2(mol2_path, copied_mol2)
    ok = bool(copied_sdf or copied_mol2)
    return {
        "pdb_id": pdb_id.upper(),
        "source": str(root),
        "source_type": "root",
        "sdf_path": copied_sdf,
        "mol2_path": copied_mol2,
        "ok": ok,
        "error": "" if ok else "pdbbind_ligand_not_found",
    }


def _zip_member_by_suffix(archive: zipfile.ZipFile, suffix: str) -> str | None:
    suffix_lc = suffix.lower()
    matches = [name for name in archive.namelist() if name.lower().endswith(suffix_lc)]
    return sorted(matches, key=lambda name: (name.count("/"), name))[0] if matches else None


def _copy_zip_member(archive: zipfile.ZipFile, member: str, dest: Path) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = _tmp_path(dest)
    tmp_path.write_bytes(archive.read(member))
    tmp_path.replace(dest)
    return str(dest)


def _copy_pdbbind_ligands_from_zip(
    *,
    pdb_id: str,
    zip_path: Path,
    dest_dir: Path,
) -> dict[str, Any]:
    pdb = pdb_id.lower()
    sdf_name, mol2_name = _pdbbind_ligand_names(pdb)
    with zipfile.ZipFile(zip_path) as archive:
        sdf_member = _zip_member_by_suffix(archive, f"/{pdb}/{sdf_name}")
        mol2_member = _zip_member_by_suffix(archive, f"/{pdb}/{mol2_name}")
        copied_sdf = (
            _copy_zip_member(archive, sdf_member, dest_dir / f"{pdb_id.upper()}_ligand.sdf")
            if sdf_member
            else ""
        )
        copied_mol2 = (
            _copy_zip_member(archive, mol2_member, dest_dir / f"{pdb_id.upper()}_ligand.mol2")
            if mol2_member
            else ""
        )
    ok = bool(copied_sdf or copied_mol2)
    return {
        "pdb_id": pdb_id.upper(),
        "source": str(zip_path),
        "source_type": "zip",
        "sdf_path": copied_sdf,
        "mol2_path": copied_mol2,
        "ok": ok,
        "error": "" if ok else "pdbbind_ligand_not_found",
    }


def _missing_pdbbind_ligand_result(
    pdb_id: str,
    source_root: Path | None,
    source_zip: Path | None,
) -> dict[str, Any]:
    return {
        "pdb_id": pdb_id,
        "source": str(source_root or source_zip or ""),
        "source_type": "",
        "sdf_path": "",
        "mol2_path": "",
        "ok": False,
        "error": "missing_pdbbind_source",
    }


def _resolve_one_representative_pdbbind_ligand(
    *,
    pdb_id: str,
    ligand_root: Path,
    source_root: Path | None,
    source_zip: Path | None,
) -> dict[str, Any]:
    target_dir = ligand_root / pdb_id.upper()
    if source_root and source_root.exists():
        return _copy_pdbbind_ligands_from_root(
            pdb_id=pdb_id,
            root=source_root,
            dest_dir=target_dir,
        )
    if source_zip and source_zip.exists():
        return _copy_pdbbind_ligands_from_zip(
            pdb_id=pdb_id,
            zip_path=source_zip,
            dest_dir=target_dir,
        )
    return _missing_pdbbind_ligand_result(pdb_id, source_root, source_zip)


def _pdbbind_ligand_manifest_row(
    pdb_id: str,
    family: str,
    role: str,
    result: Mapping[str, Any],
) -> dict[str, Any]:
    pdb_code = pdb_id.upper()
    return {
        "pdb_id": pdb_code,
        "family": family,
        "role": role,
        "ligand_id": f"{pdb_code}_ligand",
        "stage_dir": "ctrl_redock",
        "ph_label": "pH7_0",
        "source": result.get("source", ""),
        "source_type": result.get("source_type", ""),
        "sdf_path": result.get("sdf_path", ""),
        "mol2_path": result.get("mol2_path", ""),
        "ok": result.get("ok", False),
        "error": result.get("error", ""),
    }


def resolve_representative_pdbbind_ligands(
    *,
    out_dir: str | Path,
    pdbbind_root: str | Path | None = None,
    pdbbind_zip: str | Path | None = None,
) -> dict[str, Any]:
    root = Path(out_dir)
    ligand_root = root / "pdbbind_ligands"
    rows: list[dict[str, Any]] = []
    skipped_rows: list[dict[str, Any]] = []
    source_root = Path(pdbbind_root).expanduser() if pdbbind_root else None
    source_zip = Path(pdbbind_zip).expanduser() if pdbbind_zip else None
    for pdb_id, family, role in WANG2016_REPRESENTATIVE_COMPLEXES:
        result = _resolve_one_representative_pdbbind_ligand(
            pdb_id=pdb_id,
            ligand_root=ligand_root,
            source_root=source_root,
            source_zip=source_zip,
        )
        row = _pdbbind_ligand_manifest_row(pdb_id, family, role, result)
        if row.get("ok"):
            rows.append(row)
        else:
            skipped_rows.append(row)
    manifest_csv = root / "wang2016_pdbbind_ligands.csv"
    skipped_csv = root / "wang2016_pdbbind_ligands_skipped.csv"
    _write_rows(manifest_csv, rows)
    _write_rows(skipped_csv, skipped_rows)
    return {
        "ligand_dir": str(ligand_root),
        "manifest_csv": str(manifest_csv),
        "skipped_csv": str(skipped_csv),
        "ok": bool(rows),
        "rows": rows,
        "skipped_rows": skipped_rows,
    }


def _read_csv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _first_column(row: Mapping[str, str], names: Iterable[str]) -> str:
    lowered = {_slug(key): value for key, value in row.items()}
    for name in names:
        value = lowered.get(_slug(name))
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def _float_or_none(value: object) -> float | None:
    try:
        text = str(value).strip().replace("−", "-")
        if not text:
            return None
        return float(text)
    except Exception:
        return None


def _rank(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda idx: values[idx])
    ranks = [0.0] * len(values)
    idx = 0
    while idx < len(order):
        end = idx + 1
        while end < len(order) and values[order[end]] == values[order[idx]]:
            end += 1
        rank = (idx + end + 1) / 2.0
        for pos in order[idx:end]:
            ranks[pos] = rank
        idx = end
    return ranks


def _pearson(xs: Sequence[float], ys: Sequence[float]) -> float | None:
    if len(xs) < 2 or len(xs) != len(ys):
        return None
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    var_x = sum((x - mean_x) ** 2 for x in xs)
    var_y = sum((y - mean_y) ** 2 for y in ys)
    denom = math.sqrt(var_x * var_y)
    if denom == 0:
        return None
    return cov / denom


def _metrics(pairs: Sequence[tuple[float, float]]) -> dict[str, float | int | None]:
    if not pairs:
        return {"n": 0, "pearson_r": None, "spearman_r": None, "rmse": None, "mae": None}
    atlas = [pair[0] for pair in pairs]
    ref = [pair[1] for pair in pairs]
    diffs = [a - r for a, r in pairs]
    return {
        "n": len(pairs),
        "pearson_r": _pearson(atlas, ref),
        "spearman_r": _pearson(_rank(atlas), _rank(ref)),
        "rmse": math.sqrt(sum(diff * diff for diff in diffs) / len(diffs)),
        "mae": sum(abs(diff) for diff in diffs) / len(diffs),
        "bias": sum(diffs) / len(diffs),
    }


def _row_key(row: Mapping[str, str]) -> str:
    pdb = _first_column(row, ("pdb_id", "pdb", "pdbid", "structure_id"))
    ligand = _first_column(row, ("ligand_id", "ligand", "ligand_base", "compound_id"))
    return f"{pdb.upper()}::{ligand.lower()}" if ligand else pdb.upper()


def compare_against_reference(
    *,
    atlas_results: str | Path,
    reference_rows: Sequence[Mapping[str, str]],
    out_dir: str | Path,
    atlas_energy_column: str,
    reference_energy_column: str,
) -> dict[str, Any]:
    atlas_rows = _read_csv(atlas_results)
    reference_by_key = {_row_key(row): row for row in reference_rows if _row_key(row)}
    matched: list[dict[str, Any]] = []
    pairs: list[tuple[float, float]] = []
    for atlas_row in atlas_rows:
        key = _row_key(atlas_row)
        ref_row = reference_by_key.get(key)
        if not ref_row:
            continue
        atlas_value = _float_or_none(atlas_row.get(atlas_energy_column))
        ref_value = _float_or_none(ref_row.get(reference_energy_column))
        if atlas_value is None or ref_value is None:
            continue
        pairs.append((atlas_value, ref_value))
        matched.append(
            {
                "match_key": key,
                "atlas_delta_g": atlas_value,
                "reference_delta_g": ref_value,
                "difference": atlas_value - ref_value,
            }
        )
    out_path = Path(out_dir)
    matched_path = out_path / "atlas_vs_wang2016_matched.csv"
    metrics_path = out_path / "atlas_vs_wang2016_metrics.json"
    _write_rows(matched_path, matched)
    payload = {
        "atlas_results": str(atlas_results),
        "atlas_energy_column": atlas_energy_column,
        "reference_energy_column": reference_energy_column,
        "matched_csv": str(matched_path),
        "metrics": _metrics(pairs),
    }
    _write_json(metrics_path, payload)
    return {**payload, "metrics_json": str(metrics_path)}


def build_wang2016_benchmark(
    *,
    out_dir: str | Path,
    atlas_results: str | Path | None = None,
    reference_csv: str | Path | None = None,
    atlas_energy_column: str = "delta_total",
    reference_energy_column: str = "delta_g_binding",
    download: bool = True,
    download_pdbs: bool = False,
    pdbbind_root: str | Path | None = None,
    pdbbind_zip: str | Path | None = None,
    supplement_url: str | None = WANG2016_SUPPLEMENT_URL,
) -> dict[str, Any]:
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    urls = (supplement_url,) if supplement_url else WANG2016_SUPPLEMENT_URLS
    outputs: dict[str, Any] = {
        "schema": "atlas_mmgbsa_wang2016_benchmark_v1",
        "source_urls": list(urls),
        "static_targets": write_static_targets(out_path),
        "representative_complexes_csv": write_representative_complexes(out_path),
    }
    table_rows: list[dict[str, str]] = []
    if download:
        supplement, error = _download_first_supplement(out_path, urls)
        if supplement is None:
            outputs["download_error"] = error
        else:
            tables = parse_docx_tables(supplement)
            table_rows = _all_table_rows(tables)
            raw_tables_path = out_path / "wang2016_supplement_tables.json"
            raw_rows_path = out_path / "wang2016_supplement_rows.csv"
            candidates_path = out_path / "wang2016_reference_candidates.csv"
            _write_json(raw_tables_path, {"tables": tables})
            _write_rows(raw_rows_path, table_rows)
            _write_rows(candidates_path, _candidate_reference_rows(table_rows))
            outputs.update(
                {
                    "supplement_docx": str(supplement),
                    "raw_tables_json": str(raw_tables_path),
                    "raw_rows_csv": str(raw_rows_path),
                    "reference_candidates_csv": str(candidates_path),
                    "table_count": len(tables),
                    "raw_row_count": len(table_rows),
                }
            )
    ref_rows = _read_csv(reference_csv) if reference_csv else _candidate_reference_rows(table_rows)
    if atlas_results:
        outputs["comparison"] = compare_against_reference(
            atlas_results=atlas_results,
            reference_rows=ref_rows,
            out_dir=out_path,
            atlas_energy_column=atlas_energy_column,
            reference_energy_column=reference_energy_column,
        )
    if download_pdbs:
        outputs["representative_pdb_downloads"] = download_representative_pdbs(out_path)
    if pdbbind_root or pdbbind_zip:
        outputs["representative_pdbbind_ligands"] = resolve_representative_pdbbind_ligands(
            out_dir=out_path,
            pdbbind_root=pdbbind_root,
            pdbbind_zip=pdbbind_zip,
        )
    outputs["protocol_gap_report"] = _protocol_gap_report()
    _write_json(out_path / "wang2016_benchmark_manifest.json", outputs)
    return outputs


def _protocol_gap_report() -> dict[str, Any]:
    return {
        "reference_protocol": {
            "source": "Wang et al. 2016 J Comput Chem 37:2436-2446",
            "families": [
                "trypsin beta",
                "thrombin alpha",
                "CDK+PKA",
                "urokinase-type plasminogen activator",
                "beta-glucosidase A",
                "factor Xa",
            ],
            "protein_force_field": "Amber ff14SB",
            "ligand_force_field": "GAFF",
            "ligand_charges": "AM1-BCC",
            "water_model": "TIP3P",
            "water_retention": "structural waters within 4 A retained",
            "production_md": "10 ns, snapshots every 10 ps",
            "mmpbsa": "single-trajectory MMPBSA; PB/GB variants surveyed",
        },
        "atlas_comparison_notes": [
            "Use per-family rank correlation and RMSE/MAE as the primary acceptance metrics.",
            "Do not expect exact energy equality unless receptor construction, retained waters/ions, charges, PB/GB settings, and frame windows are matched.",
            "The static target CSVs capture the paper-level expected error bands when per-complex reference values are unavailable.",
        ],
    }
