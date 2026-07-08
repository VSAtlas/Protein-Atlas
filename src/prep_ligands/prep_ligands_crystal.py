import logging
import os
import re
import shutil
import contextlib
import json
from pathlib import Path
from typing import Dict, List, Set, cast

from rdkit import Chem, rdBase
from protein_prep import pdb_fixer_runtime as _pdb_fixer

from prep_ligands.prep_ligands_common_conversion import (
    _buffer_like_by_counts_from_mol,
    _buffer_like_by_counts_from_pdbfile,
    _is_polyacidic_buffer_like,
    _looks_like_buffer_salt,
    _matches_counterion,
    _polyacidic_by_counts_from_pdbfile,
)
from prep_ligands.prep_ligands_meeko import (
    prepare_sdf_to_pdbqt_with_meeko,
    restore_extracted_pdb_to_sdf,
)
from prep_ligands.prep_ligands_reporting import (
    _append_prep_status,
    _log_elem_fix_summary,
    _log_malformed,
)
from prep_ligands.prep_ligands_common_validation import (
    EXCLUDE_CRYSTAL_ADDITIVES,
    MIN_ATOMS_FOR_DOCKING,
    QUARANTINE_DIRNAME,
    STANDARD_AMINO_ACIDS,
    collapse_sanitized_path,
    is_valid_ligand,
)
from prep_ligands.prep_ligands_runtime import init_prep_workspace


_SAFE_NAME_RE = re.compile(r"[^A-Za-z0-9_.+-]+")
_EXTRACTED_LIGAND_NAME_EXCLUDES = ("nolig", "phenix_clean", "phenix-clean")


def assert_no_helium_in_hydrogen_names(text: str) -> tuple[str, int]:
    return cast(
        tuple[str, int],
        _pdb_fixer.assert_no_helium_in_hydrogen_names(text),
    )


def fix_pdb_elements(path: str) -> None:
    _pdb_fixer.fix_pdb_elements(path)


def fix_element_columns_in_file(path: str, rewrite_atoms: bool = False) -> None:
    _pdb_fixer.fix_element_columns_in_file(path, rewrite_atoms=rewrite_atoms)


def _rdkit_quiet_logs():
    try:
        return rdBase.BlockLogs()
    except Exception:
        return contextlib.nullcontext()


def _base_stem_from_pdb_path(pdb_path: Path) -> str:
    return pdb_path.stem.split(".sanitized")[0]


def _ligprep_source_sidecar(pdbqt_path: Path) -> Path:
    return pdbqt_path.with_suffix(".ligprep_source.json")


def _write_extracted_source_sidecar(
    *,
    pdbqt_path: Path,
    source_sdf: Path,
    ligand_stem: str,
    ligand_id: str,
    source_kind: str,
    restored_bond_orders: bool,
    detail: str = "",
) -> None:
    if not pdbqt_path.exists() or not source_sdf.exists():
        return
    payload = {
        "schema_version": 1,
        "writer": "meeko",
        "pdbqt_path": str(pdbqt_path),
        "ligand_stem": ligand_stem,
        "source_format": "sdf",
        "source_sdf": str(source_sdf),
        "source_record_index": 1,
        "source_record_name": ligand_id,
        "source_kind": source_kind,
        "bond_orders_restored": bool(restored_bond_orders),
        "chemistry_authoritative": bool(restored_bond_orders),
        "notes": (
            "Extracted crystal-control PDBQT was generated from this SDF; "
            "downstream docked pose SDFs should reuse this chemistry and only replace coordinates."
        ),
        "detail": detail,
    }
    sidecar = _ligprep_source_sidecar(pdbqt_path)
    sidecar.parent.mkdir(parents=True, exist_ok=True)
    tmp = sidecar.with_suffix(sidecar.suffix + ".part")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, sidecar)


def _first_hetatm_resname(pdb_path: Path) -> str:
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith("HETATM"):
                    return (line[17:20] if len(line) >= 20 else "").strip().upper()
    except Exception:
        pass
    return ""


def _looks_like_extracted_ligand_pdb(pdb_path: Path) -> tuple[bool, str]:
    name_lc = pdb_path.name.lower()
    for token in _EXTRACTED_LIGAND_NAME_EXCLUDES:
        if token in name_lc:
            return False, f"name_excluded:{token}"

    has_atom = False
    has_hetatm = False
    try:
        with open(pdb_path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("ATOM"):
                    has_atom = True
                    break
                if line.startswith("HETATM"):
                    has_hetatm = True
    except Exception as exc:
        return False, f"read_error:{exc.__class__.__name__}"

    if has_atom:
        return False, "has_atom_records"
    if not has_hetatm:
        return False, "no_hetatm_records"
    return True, "ok"


def _filter_extracted_ligand_pdbs(pdb_files: List[Path]) -> List[Path]:
    excluded_counts: Dict[str, int] = {}
    kept: List[Path] = []
    for pdb_path in pdb_files:
        ok, reason = _looks_like_extracted_ligand_pdb(pdb_path)
        logging.info(
            "[ligprep-extracted] filter %s ok=%s reason=%s",
            pdb_path.name,
            ok,
            reason,
        )
        if ok:
            kept.append(pdb_path)
        else:
            excluded_counts[reason] = excluded_counts.get(reason, 0) + 1

    logging.info(
        "[ligprep-extracted] filter_summary kept=%d excluded=%d",
        len(kept),
        len(pdb_files) - len(kept),
    )
    for reason, count in sorted(excluded_counts.items()):
        logging.info("[ligprep-extracted] filtered reason=%s count=%d", reason, count)
    return kept


def _dedupe_prefer_sanitized(pdb_files: List[Path]) -> List[Path]:
    grouped: Dict[str, List[Path]] = {}
    for pdb_path in pdb_files:
        grouped.setdefault(_base_stem_from_pdb_path(pdb_path), []).append(pdb_path)

    deduped: List[Path] = []
    for base_stem, paths in sorted(grouped.items()):
        if len(paths) == 1:
            deduped.append(paths[0])
            continue

        sanitized = [p for p in paths if p.stem.endswith(".sanitized")]
        keep = sorted(sanitized or paths)[0]
        dropped = [p for p in sorted(paths) if p != keep]
        logging.info(
            "[ligprep-extracted] dedupe base=%s keep=%s drop=%s",
            base_stem,
            keep.name,
            ",".join(p.name for p in dropped),
        )
        deduped.append(keep)

    return sorted(deduped)


def prep_ligands_from_pdb(
    ligand_output_dir: Path,
    ligands_mol2_dir: Path,
    prepped_ligands_dir: Path,
    dry_run: bool = False,
):
    """
    Crystal-safe path to prepare ligands that were extracted from PDBs (processed_pdbs/*/ligands_raw/*.pdb).
    Publication-oriented flow:
      - treat extracted PDB as coordinates, not chemistry
      - restore coordinate-preserving SDF chemistry from CCD/template when possible
      - write PDBQT with the shared Meeko backend
      - concise one-line summary per ligand
      - test-mode subset selection via EXTRACT_ONLY/EXTRACT_TEST (and mirrored CLI)
      - dry_run plans outputs without invoking external toolchains
    """
    logging.info(
        "Starting ligand preparation from PDB files "
        "(crystal-safe, restored-SDF Meeko path)"
    )

    # --- discovery root diagnostics ---
    root = ligand_output_dir.resolve()
    print(
        f"[ligprep-extracted] discovery_root={root} exists={root.exists()} is_dir={root.is_dir()}"
    )

    # --- collect EXTRACT_ONLY tokens (basename / filename / full path) ---
    raw_only = os.environ.get("EXTRACT_ONLY", "").strip()
    requested_only: Set[str] = set()
    if raw_only:
        for tok in raw_only.replace(",", " ").split():
            t = tok.strip()
            if t:
                requested_only.add(t)

    # If EXTRACT_ONLY contains any *existing* .pdb paths, use them directly
    direct_paths: List[Path] = []
    for tok in list(requested_only):
        p = Path(tok)
        if p.suffix.lower() == ".pdb" and p.exists():
            direct_paths.append(p.resolve())

    missing: Set[str] = set()
    if direct_paths:
        pdb_files = sorted(set(direct_paths))
        print(
            f"[ligprep-extracted] direct-path mode count={len(pdb_files)} keep={','.join(p.stem for p in pdb_files)}"
        )
    else:
        # Fall back to recursive discovery under the configured root
        pdb_files = sorted(root.rglob("*.pdb"))
        logging.info(f"Found {len(pdb_files)} PDB ligand file(s) under {root}")

        # If subset requested by name/stem, filter against discovered set
        if requested_only:
            by_stem: Dict[str, List[Path]] = {}
            by_name: Dict[str, List[Path]] = {}
            by_abs: Dict[str, Path] = {}
            for p in pdb_files:
                by_stem.setdefault(p.stem, []).append(p)
                by_name.setdefault(p.name, []).append(p)
                by_abs[str(p.resolve())] = p

            keep: List[Path] = []
            for r in requested_only:
                if r in by_abs:
                    keep.append(by_abs[r])
                    continue
                if r in by_name and by_name[r]:
                    keep.extend(by_name[r])
                    continue
                stem = Path(r).stem if Path(r).suffix else r
                if stem in by_stem and by_stem[stem]:
                    keep.extend(by_stem[stem])
                    continue
                missing.add(r)

            pdb_files = sorted(set(keep))

        print(
            f"[ligprep-extracted] test-mode enabled count={len(pdb_files)} missing={len(missing)} keep={','.join(Path(p).stem for p in pdb_files)}"
        )
        if missing:
            print(
                f"[ligprep-extracted] warn: requested_not_found={','.join(sorted(missing))}"
            )

    logging.info("[ligprep-extracted] discovery_count=%d", len(pdb_files))
    pdb_files = _filter_extracted_ligand_pdbs(pdb_files)
    pdb_files = _dedupe_prefer_sanitized(pdb_files)
    logging.info(
        "[ligprep-extracted] selected_after_filter count=%d keep=%s",
        len(pdb_files),
        ",".join(p.stem for p in pdb_files),
    )

    if dry_run:
        planned_outputs = [
            (
                prepped_ligands_dir / f"{_base_stem_from_pdb_path(p)}.sanitized.pdbqt"
            ).name
            for p in pdb_files
        ]
        return {
            "selected_inputs": [p.name for p in pdb_files],
            "planned_outputs": planned_outputs,
        }

    # Verbose toggle like bulk
    if str(os.environ.get("EXTRACT_TEST", "0")).lower() not in {"0", "false", "no"}:
        os.environ["LIGPREP_DEBUG"] = "1"

    status_log = init_prep_workspace(
        prepped_ligands_dir,
        source_link_target=ligand_output_dir,
        create_reference_dirs=True,
    )

    for pdb_file in pdb_files:
        normalized_src = collapse_sanitized_path(pdb_file)
        if normalized_src != pdb_file:
            try:
                pdb_file.rename(normalized_src)
                logging.info(
                    "[sanitize] normalized ligand filename %s -> %s",
                    pdb_file.name,
                    normalized_src.name,
                )
                pdb_file = normalized_src
            except Exception as e:
                logging.warning(
                    "[sanitize] unable to normalize ligand filename %s -> %s: %s",
                    pdb_file.name,
                    normalized_src.name,
                    e,
                )

        logging.info(f"Processing: {pdb_file.name}")

        # Ensure the ligand ID is based on the canonical sanitized stem to keep output names stable.
        base_stem = pdb_file.stem.split(".sanitized")[0]
        lig_id = f"{base_stem}.sanitized"

        residue_name = _first_hetatm_resname(pdb_file) or lig_id.split("_")[0].upper()

        # Skip standard amino acids and known crystallization additives
        if residue_name in STANDARD_AMINO_ACIDS:
            logging.info(f"Skipping standard amino acid residue: {pdb_file.name}")
            continue
        if residue_name in EXCLUDE_CRYSTAL_ADDITIVES:
            _log_malformed(
                pdb_file,
                f"excluded_crystal_additive:{residue_name}",
                log_dir=prepped_ligands_dir,
            )
            logging.info(f"Skipping crystallization additive: {pdb_file.name}")
            continue

        # Quick size gate
        with open(pdb_file, "r", encoding="utf-8", errors="ignore") as f:
            atom_lines = [line for line in f if line.startswith(("HETATM", "ATOM"))]
        if len(atom_lines) < MIN_ATOMS_FOR_DOCKING:
            _log_malformed(
                pdb_file,
                f"tiny_ligand_fewer_than_{MIN_ATOMS_FOR_DOCKING}_atoms",
                log_dir=prepped_ligands_dir,
            )
            logging.info(f"Skipping tiny ligand: {pdb_file.name}")
            continue

        # Sanitize PDB: keep only the ligand's HETATM triplet; drop all ATOM lines (protein)
        sanitized = collapse_sanitized_path(pdb_file.with_suffix(".sanitized.pdb"))

        def _sanitize_pdb(in_pdb: Path, out_pdb: Path) -> bool:
            target = None  # (resn, chain, resi)
            wrote_any = False

            def _triplet(ln: str):
                resn = (ln[17:20] if len(ln) >= 20 else "").strip()
                chain = (ln[21] if len(ln) >= 22 else " ").strip()
                resi = (ln[22:26] if len(ln) >= 26 else "").strip()
                return resn, chain, resi

            with open(in_pdb, "r", encoding="utf-8", errors="ignore") as fin:
                lines = fin.readlines()

            for ln in lines:
                if ln.startswith("HETATM"):
                    target = _triplet(ln)
                    break

            if target:
                logging.info(
                    f"[sanitize] {in_pdb.name}: keeping ligand {target[0]} chain {target[1]} resi {target[2]}"
                )
            else:
                logging.warning(
                    f"[sanitize] {in_pdb.name}: no HETATM found; nothing to keep"
                )

            with open(out_pdb, "w", encoding="utf-8") as fout:
                for ln in lines:
                    if not ln.startswith(("ATOM", "HETATM")):
                        fout.write(ln)
                        continue
                    if ln.startswith("ATOM"):
                        continue
                    if target is None:
                        continue
                    resn, chain, resi = _triplet(ln)
                    if (resn, chain, resi) != target:
                        continue
                    altloc = ln[16] if len(ln) > 16 else " "
                    keep_alt = {"", " ", "A"} | set("0123456789")
                    a = altloc.strip()
                    if a and (a not in keep_alt and a.upper() not in keep_alt):
                        continue
                    try:
                        x = float(ln[30:38])
                        y = float(ln[38:46])
                        z = float(ln[46:54])
                        if (
                            any([(x != x), (y != y), (z != z)])
                            or max(abs(x), abs(y), abs(z)) > 1e6
                        ):
                            continue
                    except Exception:
                        continue
                    fout.write(ln)
                    wrote_any = True

            if not wrote_any:
                logging.warning(
                    f"[control-prep] {in_pdb.name}: no HETATM kept after residue filtering."
                )
            return wrote_any and out_pdb.exists() and out_pdb.stat().st_size > 0

        ok_san = _sanitize_pdb(pdb_file, sanitized)
        if not ok_san:
            _log_malformed(pdb_file, "sanitize_kept_no_atoms")
            logging.warning(f"Sanitization produced no atoms: {pdb_file.name}")
            continue
        # Element-field fixes (column-accurate) and helium guard preflight
        try:
            pre_txt = sanitized.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            pre_txt = ""
        try:
            fix_element_columns_in_file(str(sanitized), rewrite_atoms=False)
            fixed_txt, nname = assert_no_helium_in_hydrogen_names(
                sanitized.read_text(encoding="utf-8", errors="ignore")
            )
            if nname > 0:
                sanitized.write_text(fixed_txt, encoding="utf-8")
            post_txt = (
                fixed_txt
                if nname > 0
                else sanitized.read_text(encoding="utf-8", errors="ignore")
            )
            _log_elem_fix_summary(
                sanitized,
                stage="preflight",
                before_text=pre_txt,
                after_text=post_txt,
            )
        except Exception as e:
            logging.warning(f"[elements] preflight failed for {sanitized.name}: {e}")
        try:
            fix_pdb_elements(str(sanitized))
            logging.info(f"[elements] fixed element fields: {sanitized.name}")
        except Exception as e:
            logging.warning(
                f"[elements] could not fix elements for {sanitized.name}: {e}"
            )
        residue_name = _first_hetatm_resname(sanitized) or residue_name

        # Counter-ion/buffer heuristics (quick bail-outs on obvious salts/buffers)
        flag_buffer = False
        reason = None
        try:
            with _rdkit_quiet_logs():
                m_chk = Chem.MolFromPDBFile(
                    str(sanitized), sanitize=False, removeHs=False
                )
                if m_chk is not None:
                    if _buffer_like_by_counts_from_mol(m_chk):
                        flag_buffer, reason = True, "buffer_like_by_counts"
                    if not flag_buffer:
                        try:
                            Chem.SanitizeMol(m_chk)
                        except Exception:
                            pass
                        counterion_reason = _matches_counterion(m_chk)
                        if counterion_reason:
                            reason = counterion_reason
                            flag_buffer = True
                        elif _looks_like_buffer_salt(m_chk):
                            reason = "buffer_like"
                            flag_buffer = True
                        if (
                            not flag_buffer
                            and residue_name in {"UNL", "LIG"}
                            and _is_polyacidic_buffer_like(m_chk)
                        ):
                            flag_buffer, reason = True, "polyacidic_buffer_like"
        except Exception:
            if _buffer_like_by_counts_from_pdbfile(sanitized) or (
                residue_name in {"UNL", "LIG"}
                and _polyacidic_by_counts_from_pdbfile(sanitized)
            ):
                flag_buffer, reason = True, "counts_only_polyacidic"

        if flag_buffer:
            _log_malformed(pdb_file, f"counterion_or_buffer:{reason or 'unknown'}")
            logging.info(
                f"Skipping likely counter-ion/buffer ({reason or 'unknown'}): {pdb_file.name}"
            )
            continue

        # Additional UNL O-rich ringless fragment guard
        try:
            with _rdkit_quiet_logs():
                if (
                    "m_chk" in locals()
                    and m_chk is not None
                    and residue_name in {"UNL", "LIG"}
                ):
                    from rdkit.Chem import rdMolDescriptors as rdmd

                    rings = rdmd.CalcNumRings(m_chk)
                    arom = rdmd.CalcNumAromaticRings(m_chk)
                    o = sum(1 for a in m_chk.GetAtoms() if a.GetSymbol() == "O")
                    hac = m_chk.GetNumHeavyAtoms()
                    if (
                        rings == 0
                        and arom == 0
                        and hac >= 12
                        and (o / float(hac)) >= 0.40
                    ):
                        _log_malformed(
                            pdb_file, "UNL_O_rich_ringless", log_dir=prepped_ligands_dir
                        )
                        logging.info(
                            f"Skipping UNL O-rich ringless fragment: {pdb_file.name}"
                        )
                        continue
        except Exception:
            pass

        # Write a pristine SDF snapshot for reference/debug
        ref_dir = sanitized.parent.parent / "reference"
        ref_dir.mkdir(parents=True, exist_ok=True)
        ref_sdf = ref_dir / (pdb_file.stem + ".sdf")
        try:
            mref = Chem.MolFromPDBFile(str(sanitized), sanitize=False, removeHs=False)
            if mref is not None:
                w = Chem.SDWriter(str(ref_sdf))
                try:
                    w.SetKekulize(False)
                except Exception:
                    pass
                w.write(mref)
                w.close()
        except Exception:
            with open(ref_sdf, "w") as out:
                out.write(pdb_file.stem + "\n$$$$\n")

        # Resume check
        # Force canonical output name (ensure .sanitized.pdbqt even if input was just .pdb)
        pdbqt_path = collapse_sanitized_path(
            prepped_ligands_dir / f"{base_stem}.sanitized.pdbqt"
        )
        restored_sdf = ref_dir / f"{base_stem}.restored.sdf"

        # Migration: if legacy non-sanitized PDBQT exists, move/copy it to canonical path
        legacy_pdbqt = prepped_ligands_dir / f"{base_stem}.pdbqt"
        if legacy_pdbqt.exists() and not pdbqt_path.exists():
            try:
                logging.info(
                    f"[migration] moving legacy {legacy_pdbqt.name} -> {pdbqt_path.name}"
                )
                shutil.copy2(legacy_pdbqt, pdbqt_path)
            except Exception as e:
                logging.warning(
                    f"[migration] failed to migrate {legacy_pdbqt.name}: {e}"
                )

        if (
            pdbqt_path.exists()
            and pdbqt_path.stat().st_size > 100
            and is_valid_ligand(pdbqt_path, log_dir=prepped_ligands_dir)
        ):
            if not _ligprep_source_sidecar(pdbqt_path).exists():
                source_sdf = restored_sdf if restored_sdf.exists() else ref_sdf
                _write_extracted_source_sidecar(
                    pdbqt_path=pdbqt_path,
                    source_sdf=source_sdf,
                    ligand_stem=base_stem,
                    ligand_id=lig_id,
                    source_kind="resume_existing",
                    restored_bond_orders=False,
                )
            logging.info(
                f"[resume] Valid PDBQT already exists, skipping: {pdbqt_path.name}"
            )
            continue

        try:
            restored = restore_extracted_pdb_to_sdf(
                sanitized,
                restored_sdf,
                residue_name=residue_name,
                template_cache_dir=ref_dir / "ccd",
            )
            result = prepare_sdf_to_pdbqt_with_meeko(
                restored.sdf_path,
                pdbqt_path,
                ligand_name=lig_id,
                log_dir=prepped_ligands_dir,
            )
            _append_prep_status(
                status_log,
                sanitized.name,
                "OK" if result.ok else "FAIL",
                "" if result.ok else result.status,
                str(pdbqt_path.relative_to(prepped_ligands_dir))
                if pdbqt_path.exists()
                else "",
            )
            if result.ok:
                _write_extracted_source_sidecar(
                    pdbqt_path=pdbqt_path,
                    source_sdf=restored.sdf_path,
                    ligand_stem=base_stem,
                    ligand_id=lig_id,
                    source_kind=restored.template_source,
                    restored_bond_orders=restored.restored_bond_orders,
                    detail=restored.detail,
                )
            print(
                "[ligprep-extracted] "
                f"lig={lig_id} input_fmt=PDB restored_sdf={restored.sdf_path.name} "
                f"template={restored.template_source} "
                f"bond_orders_restored={restored.restored_bond_orders} "
                f"writer=meeko ok={result.ok} reason={result.status}"
            )
        except Exception as exc:
            logging.warning(
                "[ligprep-extracted] Meeko restored-SDF path failed ligand=%s err=%s",
                lig_id,
                exc,
            )
            _append_prep_status(
                status_log,
                sanitized.name,
                "FAIL",
                f"meeko_restore_failed:{exc}",
                "",
            )

        try:
            has_protein_res = False
            if pdbqt_path.exists():
                with open(pdbqt_path, "r", errors="ignore") as fh:
                    for ln in fh:
                        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
                            continue
                        resn = (ln[17:20] if len(ln) >= 20 else "").strip().upper()
                        # only treat as fatal if the ligand *itself* is a protein residue
                        if resn in STANDARD_AMINO_ACIDS and resn == residue_name:
                            has_protein_res = True
                            break
            if has_protein_res:
                _log_malformed(
                    pdb_file, "protein_residue_in_pdbqt", log_dir=prepped_ligands_dir
                )
                quarantine = prepped_ligands_dir / QUARANTINE_DIRNAME
                quarantine.mkdir(exist_ok=True)
                try:
                    if pdbqt_path.exists():
                        pdbqt_path.replace(quarantine / pdbqt_path.name)
                except Exception:
                    pass
                try:
                    _append_prep_status(
                        status_log,
                        sanitized.name,
                        "FAIL",
                        "protein_residue_in_pdbqt",
                        str(
                            (quarantine / pdbqt_path.name).relative_to(
                                prepped_ligands_dir
                            )
                        )
                        if (quarantine / pdbqt_path.name).exists()
                        else "",
                    )
                except Exception:
                    pass
                continue
        except Exception:
            logging.warning(
                f"[postcheck] unable to scan for protein ATOM in {pdbqt_path.name}"
            )

        if (
            (not pdbqt_path.exists())
            or (pdbqt_path.stat().st_size < 100)
            or (not is_valid_ligand(pdbqt_path, log_dir=prepped_ligands_dir))
        ):
            _log_malformed(
                pdb_file, "pdbqt_postcheck_fail_or_small", log_dir=prepped_ligands_dir
            )
            quarantine = prepped_ligands_dir / QUARANTINE_DIRNAME
            quarantine.mkdir(exist_ok=True)
            try:
                if pdbqt_path.exists():
                    pdbqt_path.replace(quarantine / pdbqt_path.name)
            except Exception:
                pass
            _append_prep_status(
                status_log,
                sanitized.name,
                "FAIL",
                "postcheck_fail_or_small",
                str((quarantine / pdbqt_path.name).relative_to(prepped_ligands_dir))
                if (quarantine / pdbqt_path.name).exists()
                else "",
            )
            continue

        logging.info(f"Created PDBQT: {pdbqt_path.name}")
        _append_prep_status(
            status_log,
            sanitized.name,
            "OK",
            "",
            str(pdbqt_path.relative_to(prepped_ligands_dir)),
        )


__all__ = ["prep_ligands_from_pdb"]
