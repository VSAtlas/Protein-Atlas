#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
from typing import Iterable


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a small DUD-style dry-bench ligand library and SPD PDB canary list."
    )
    parser.add_argument(
        "--spd-selected-csv",
        default="analysis/gene_list/spd_targets_final_93_v2_selected.csv",
    )
    parser.add_argument("--pdb-limit", type=int, default=90)
    parser.add_argument("--out-prefix", default="analysis/gene_list/spd_targets_dry90")
    parser.add_argument("--source-library", default="prepped_ligands/bench_mk01")
    parser.add_argument("--library-name", default="dry_bench_64")
    parser.add_argument("--actives", type=int, default=16)
    parser.add_argument("--decoys", type=int, default=48)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _read_spd_rows(path: Path, limit: int) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            pdb_id = str(row.get("pdb_id", "")).strip().upper()
            if not pdb_id or pdb_id in seen:
                continue
            seen.add(pdb_id)
            row["pdb_id"] = pdb_id
            rows.append(dict(row))
            if len(rows) >= int(limit):
                break
    if not rows:
        raise RuntimeError(f"No PDB IDs found in {path}")
    return rows


def _write_spd_outputs(rows: list[dict[str, str]], out_prefix: Path) -> tuple[Path, Path]:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    selected_csv = out_prefix.with_name(f"{out_prefix.name}_selected.csv")
    pdb_list = out_prefix.with_name(f"{out_prefix.name}_pdbs.txt")
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with selected_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    pdb_list.write_text("\n".join(row["pdb_id"] for row in rows) + "\n", encoding="utf-8")
    return selected_csv, pdb_list


def _pick_ligands(source_root: Path, subdir: str, count: int) -> list[Path]:
    candidates = sorted((source_root / subdir).glob("*.pdbqt"))
    if len(candidates) < int(count):
        raise RuntimeError(
            f"Source library {source_root / subdir} has {len(candidates)} ligands; need {count}"
        )
    return candidates[: int(count)]


def _safe_symlink(source: Path, dest: Path, *, force: bool) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        if not force:
            return
        dest.unlink()
    rel_source = os.path.relpath(source.resolve(), dest.parent.resolve())
    dest.symlink_to(rel_source)


def _manifest_for(root: Path, rel_paths: Iterable[Path]) -> dict[str, object]:
    entries: dict[str, str] = {}
    filenames: dict[str, str] = {}
    for rel_path in rel_paths:
        rel = str(rel_path).replace("\\", "/")
        name = Path(rel).name
        stem = Path(name).stem.lower()
        filenames.setdefault(name.lower(), rel)
        entries.setdefault(stem, rel)
        base = stem.split("_stage")[0]
        if base:
            entries.setdefault(base, rel)
    child_mtimes: dict[str, float] = {}
    for child in root.iterdir():
        if child.is_dir():
            child_mtimes[child.name] = float(child.stat().st_mtime)
    return {
        "mtime": float(root.stat().st_mtime),
        "child_dir_mtimes": child_mtimes,
        "entries": entries,
        "filenames": filenames,
    }


def _create_library(
    *,
    source_root: Path,
    dest_root: Path,
    active_count: int,
    decoy_count: int,
    force: bool,
) -> Path:
    active_sources = _pick_ligands(source_root, "actives", active_count)
    decoy_sources = _pick_ligands(source_root, "decoys", decoy_count)
    rel_paths: list[Path] = []
    for source in active_sources:
        rel = Path("actives") / source.name
        _safe_symlink(source, dest_root / rel, force=force)
        rel_paths.append(rel)
    for source in decoy_sources:
        rel = Path("decoys") / source.name
        _safe_symlink(source, dest_root / rel, force=force)
        rel_paths.append(rel)
    manifest = _manifest_for(dest_root, rel_paths)
    manifest_path = dest_root / "_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest_path


def main() -> int:
    args = _parse_args()
    root = _repo_root()
    rows = _read_spd_rows(root / args.spd_selected_csv, args.pdb_limit)
    selected_csv, pdb_list = _write_spd_outputs(rows, root / args.out_prefix)
    source_root = root / args.source_library
    dest_root = root / "prepped_ligands" / str(args.library_name)
    manifest = _create_library(
        source_root=source_root,
        dest_root=dest_root,
        active_count=int(args.actives),
        decoy_count=int(args.decoys),
        force=bool(args.force),
    )
    summary = {
        "pdb_count": len(rows),
        "pdb_list": str(pdb_list.relative_to(root)),
        "selected_csv": str(selected_csv.relative_to(root)),
        "library": str(dest_root.relative_to(root)),
        "manifest": str(manifest.relative_to(root)),
        "actives": int(args.actives),
        "decoys": int(args.decoys),
    }
    summary_path = root / "outputs" / "data" / "dry_bench_canary_assets.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
