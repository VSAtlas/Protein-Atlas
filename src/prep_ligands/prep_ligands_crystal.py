import logging
import os
import re
import shutil
import subprocess
import contextlib
from pathlib import Path
from typing import Dict, List, Optional, Set

from pdb_fixer import (
    assert_no_helium_in_hydrogen_names,
    fix_pdb_elements,
    fix_element_columns_in_file,
    rules_version,
)
from rdkit import Chem, rdBase

from prep_ligands.prep_ligands_common import (
    EXCLUDE_CRYSTAL_ADDITIVES,
    MIN_ATOMS_FOR_DOCKING,
    QUARANTINE_DIRNAME,
    STANDARD_AMINO_ACIDS,
    _append_prep_status,
    _audit_protonation_metrics,
    _buffer_like_by_counts_from_mol,
    _buffer_like_by_counts_from_pdbfile,
    _count_explicit_H_in_mol2,
    _helium_postwrite_guard,
    _is_polyacidic_buffer_like,
    _log_elem_fix_summary,
    _log_malformed,
    _log_std_diff,
    _matches_counterion,
    _pdbqt_from_mol2_via_obabel,
    _pdbqt_has_H,
    _polyacidic_by_counts_from_pdbfile,
    _prepare_one,
    _re_aromatize_mol2_in_place,
    _resolve_prepare_ligand4,
    _run_obabel,
    collapse_sanitized_path,
    _looks_like_buffer_salt,
    get_short_path_name,
    is_valid_ligand,
    read_config,
    standardize_mol_with_activesite,
)


_SAFE_NAME_RE = re.compile(r"[^A-Za-z0-9_.+-]+")
_EXTRACTED_LIGAND_NAME_EXCLUDES = ("nolig", "phenix_clean", "phenix-clean")


def _rdkit_quiet_logs():
    try:
        return rdBase.BlockLogs()
    except Exception:
        return contextlib.nullcontext()


def _base_stem_from_pdb_path(pdb_path: Path) -> str:
    return pdb_path.stem.split(".sanitized")[0]


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
    Hardened to mirror bulk path resilience:
      - layered protonation (RDKit AddHs -> OBabel -h -> ADT/Meeko during typing)
      - unified writer via _prepare_one(...) with AD4 normalization + invariant checks
      - crystal-safety scrubs + helium guard
      - concise one-line summary per ligand
      - test-mode subset selection via EXTRACT_ONLY/EXTRACT_TEST (and mirrored CLI)
      - dry_run plans outputs without invoking external toolchains
    """
    logging.info(
        "Starting ligand preparation from PDB files (crystal-safe, MOL2-first path)"
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

    if direct_paths:
        pdb_files = sorted(set(direct_paths))
        print(
            f"[ligprep-extracted] direct-path mode count={len(pdb_files)} keep={','.join(p.stem for p in pdb_files)}"
        )
        missing = set()  # nothing to match; we used the paths verbatim
    else:
        # Fall back to recursive discovery under the configured root
        pdb_files = sorted(root.rglob("*.pdb"))
        logging.info(f"Found {len(pdb_files)} PDB ligand file(s) under {root}")

        # If subset requested by name/stem, filter against discovered set
        missing: Set[str] = set()
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

    config = read_config()
    mgltools_python = config.get("MGLTOOLS_PYTHON")
    mgltools_path = config.get("MGLTOOLS_PATH")
    obabel_exe_cfg = config.get("OPENBABEL_PATH")

    if not mgltools_python or not mgltools_path or not obabel_exe_cfg:
        raise RuntimeError(
            "Missing paths in config.txt: MGLTOOLS_PYTHON, MGLTOOLS_PATH, OPENBABEL_PATH"
        )

    obabel_exe = obabel_exe_cfg
    obabel_exe_short = get_short_path_name(obabel_exe)

    mgltools_python_short = get_short_path_name(mgltools_python)
    prepare_script = _resolve_prepare_ligand4(mgltools_path, config)
    if not prepare_script.exists():
        raise FileNotFoundError(
            f"prepare_ligand4.py not found at {prepare_script} (set PREPARE_LIGAND4 in config.txt)"
        )
    prepare_script_short = get_short_path_name(str(prepare_script.resolve()))

    # Verbose toggle like bulk
    if str(os.environ.get("EXTRACT_TEST", "0")).lower() not in {"0", "false", "no"}:
        os.environ["LIGPREP_DEBUG"] = "1"

    prepped_ligands_dir.mkdir(parents=True, exist_ok=True)
    status_log = prepped_ligands_dir / "ligand_prep_status.tsv"

    # Intermediates for easier debugging + symlink back to source dir (per-protein; no cfg)
    try:
        inter_dirname = "intermediates"
        ref_dirname = "reference"
        quar_dirname = "quarantine"
        inter_dir = prepped_ligands_dir / inter_dirname
        ref_dir = prepped_ligands_dir / ref_dirname
        quar_dir = prepped_ligands_dir / quar_dirname
        inter_dir.mkdir(parents=True, exist_ok=True)
        ref_dir.mkdir(parents=True, exist_ok=True)
        quar_dir.mkdir(parents=True, exist_ok=True)

        link = prepped_ligands_dir / "intermediates_src"
        if link.exists() or link.is_symlink():
            try:
                if link.is_dir() and not link.is_symlink():
                    shutil.rmtree(link)
                else:
                    link.unlink()
            except Exception:
                pass
        link.symlink_to(ligand_output_dir.resolve(), target_is_directory=True)
        logging.info("[intermediates] symlinked source -> intermediates_src")
    except Exception as e:
        logging.debug("intermediates link skipped: %s", e)

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

        tmp_pdb = str(pdb_file)  # updated to sanitized below
        residue_name = lig_id.split("_")[0].upper()

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
        tmp_pdb = str(sanitized)

        # Element-field fixes (column-accurate) and helium guard preflight
        try:
            pre_txt = sanitized.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            pre_txt = None
        try:
            fix_element_columns_in_file(str(sanitized), rewrite_atoms=False)
            fixed_txt, nname = assert_no_helium_in_hydrogen_names(
                sanitized.read_text(encoding="utf-8", errors="ignore")
            )
            if nname > 0:
                sanitized.write_text(fixed_txt, encoding="utf-8")
            _log_elem_fix_summary(
                sanitized,
                stage="preflight",
                before_text=pre_txt,
                after_text=fixed_txt if nname > 0 else None,
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
                        reason = _matches_counterion(m_chk) or (
                            _looks_like_buffer_salt(m_chk) and "buffer_like"
                        )
                        if reason:
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
            logging.info(
                f"[resume] Valid PDBQT already exists, skipping: {pdbqt_path.name}"
            )
            continue

        # Primary PDB->MOL2 (no H) as a baseline file; later protonation layers may produce protoA/protoB
        tmp_mol2 = sanitized.with_suffix(".mol2")
        try:
            ob_cmd = [
                obabel_exe_short,
                "-ipdb",
                get_short_path_name(str(sanitized.resolve())),
                "-omol2",
                "-O",
                get_short_path_name(str(tmp_mol2.resolve())),
            ]
            logging.info(
                "Converting sanitized PDB -> MOL2 (no gen3d): "
                + " ".join(map(str, ob_cmd))
            )
            res_ob = subprocess.run(
                ob_cmd, check=True, capture_output=True, text=True, timeout=300
            )
            if res_ob.stderr:
                logging.warning("[obabel primary stderr] %s", res_ob.stderr.strip())
        except subprocess.CalledProcessError as e:
            logging.warning(
                f"OBabel PDB->MOL2 failed for {sanitized.name}:\n{e.stderr}"
            )

        if tmp_mol2.exists() and tmp_mol2.stat().st_size > 100:
            # Normalize aromaticity on MOL2 so downstream typing is stable
            try:
                _re_aromatize_mol2_in_place(tmp_mol2, obabel_exe_short)
            except Exception as e:
                logging.error("[ligprep] re_arom crash for %s: %s", tmp_mol2.name, e)

            # --- PROTONATION: layered strategy on the extracted fragment ---
            _ = _audit_protonation_metrics(
                tag=lig_id, mol_or_path=tmp_pdb, context="pre"
            )
            h_strategy_label: List[str] = []
            fallbacks_used: List[str] = []

            # Strategy A: RDKit AddHs on sanitized PDB -> MOL2
            proto_mol2: Optional[str] = None
            try:
                mH = Chem.MolFromPDBFile(tmp_pdb, sanitize=False, removeHs=False)
                if mH is None:
                    raise ValueError("rdkit_from_pdb_failed")
                mH = Chem.AddHs(mH, addCoords=True)
                try:
                    Chem.Kekulize(mH, clearAromaticFlags=True)
                except Exception:
                    pass
                protoA = sanitized.with_suffix(".protoA.mol2")
                Chem.MolToMol2File(mH, str(protoA))
                if protoA.exists() and protoA.stat().st_size > 100:
                    proto_mol2 = str(protoA)
                    h_strategy_label.append("RDKit-AddHs")
            except Exception:
                fallbacks_used.append("RDKit-AddHs-fail")

            # Try Strategy B (OBabel -p -h) unconditionally
            protoB = sanitized.with_suffix(".protoB.mol2")
            ph = os.environ.get("LIGPREP_PH", "7.4")
            okB = _run_obabel(
                [
                    obabel_exe_short,
                    "-ipdb",
                    str(sanitized),
                    "-omol2",
                    "-O",
                    str(protoB),
                    "-p",
                    str(ph),
                    "--partialcharge",
                    "gasteiger",
                ],
                timeout_sec=600,
            )

            # Decide best candidate by explicit-H count
            cands = []
            for p in (protoA, protoB):
                if p.exists() and p.stat().st_size > 100:
                    cands.append((p, _count_explicit_H_in_mol2(p)))
            if cands:
                best, bestH = max(cands, key=lambda t: t[1])
                proto_mol2 = str(best)
                logging.info(
                    f"[choose-mol2] picked={best.name} H={bestH} "
                    f"others={[(x[0].name, x[1]) for x in cands]}"
                )
            else:
                proto_mol2 = None

            # Strategy C: no extra H; rely on ADT/Meeko typing additions
            if proto_mol2 is None:
                proto_mol2 = str(tmp_mol2)
                h_strategy_label.append("no-extra-H (ADT adds)")
                fallbacks_used.append("ADT/Meeko-rescue-possible")

            # --- Choose protoA/protoB by explicit H count if either exists ---
            protoA_path = sanitized.with_suffix(".protoA.mol2")
            protoB_path = sanitized.with_suffix(".protoB.mol2")

            cands_paths: List[Path] = [
                p for p in (protoA_path, protoB_path) if p and p.exists()
            ]

            proto_mol2 = ""
            mol2_for_mgl_path: Path | None = None
            bestH: int = -1

            if cands_paths:
                counts = [(p, _count_explicit_H_in_mol2(p)) for p in cands_paths]
                mol2_for_mgl_path, bestH = max(counts, key=lambda t: t[1])
                proto_mol2 = str(mol2_for_mgl_path)
                logging.info(
                    "[choose-mol2] picked=%s H=%s others=%s",
                    mol2_for_mgl_path.name,
                    bestH,
                    [(p.name, h) for p, h in counts],
                )
            else:
                fallback = (
                    tmp_mol2
                    if "tmp_mol2" in locals() and tmp_mol2.exists()
                    else sanitized.with_suffix(".mol2")
                )
                mol2_for_mgl_path = Path(fallback)
                proto_mol2 = str(mol2_for_mgl_path)
                logging.warning(
                    "[choose-mol2] no protoA/protoB found; falling back to %s (exists=%s size=%s)",
                    mol2_for_mgl_path.name,
                    mol2_for_mgl_path.exists(),
                    (
                        mol2_for_mgl_path.stat().st_size
                        if mol2_for_mgl_path.exists()
                        else -1
                    ),
                )

            mol2_for_mgl: Path = Path(proto_mol2)

            # Keep aromaticity consistent on proto_mol2 as well
            try:
                _re_aromatize_mol2_in_place(mol2_for_mgl, obabel_exe_short)
            except Exception as e:
                logging.error(
                    "[ligprep] re_arom crash for %s: %s", mol2_for_mgl.name, e
                )

            # --- Enforce active-site standardization on extracted ligands (parity with bulk) ---
            try:
                m_raw = Chem.MolFromMol2File(
                    str(proto_mol2), sanitize=False, removeHs=False
                )
                old_smiles = ""
                try:
                    m_old = Chem.MolFromMol2File(
                        str(proto_mol2), sanitize=True, removeHs=False
                    )
                    if m_old:
                        old_smiles = Chem.MolToSmiles(m_old, isomericSmiles=True)
                except Exception:
                    pass

                m_std = standardize_mol_with_activesite(m_raw)
                if m_std is not None:
                    new_smiles = ""
                    try:
                        new_smiles = Chem.MolToSmiles(m_std, isomericSmiles=True)
                    except Exception:
                        pass
                    if old_smiles and new_smiles and old_smiles != new_smiles:
                        logging.warning(
                            "[std:audit][extracted] %s: SMILES changed %s -> %s",
                            Path(proto_mol2).name,
                            old_smiles,
                            new_smiles,
                        )
                        _log_std_diff(
                            prepped_ligands_dir,
                            lig_id,
                            "extracted",
                            old_smiles,
                            new_smiles,
                        )
            except Exception as e:
                logging.warning(
                    "[std:extracted] standardization skipped due to error: %s", e
                )

            ok_write, status = _prepare_one(
                mgltools_python_short,
                prepare_script_short,
                Path(proto_mol2),
                pdbqt_path,
                obabel_exe_short=obabel_exe_short,
                status_log_dir=prepped_ligands_dir,
            )
            ok_write = status == "ok"
            chosen_writer = "auto"
            reason = status

            # Post-write guard (extracted path): ensure explicit H present; OBabel re-write if not
            if (
                pdbqt_path.exists()
                and obabel_exe_short
                and not _pdbqt_has_H(pdbqt_path)
            ):
                logging.warning(
                    "[post-extracted] %s has no H; OBabel -h re-write", pdbqt_path.name
                )
                _pdbqt_from_mol2_via_obabel(
                    Path(proto_mol2), pdbqt_path, obabel_exe_short
                )
                logging.info(
                    "[post-extracted] re-write complete; has_H=%s",
                    _pdbqt_has_H(pdbqt_path),
                )

            # Post-write guards: helium + invariants (waters/ions/torsdof/charges/types)
            _he_ok, _he_fixes, _he_q = _helium_postwrite_guard(
                pdbqt_path, lig_id, prepped_ligands_dir
            )
            valid_ok = is_valid_ligand(pdbqt_path, prepped_ligands_dir)

            post_metrics = _audit_protonation_metrics(
                tag=lig_id, mol_or_path=pdbqt_path, context="post"
            )

            print(
                "[ligprep-extracted] "
                f"lig={lig_id} input_fmt=PDB "
                f"H_strategy={'+'.join(h_strategy_label) or 'none'} "
                f"fallbacks={fallbacks_used or '[]'} "
                f"writer={chosen_writer} ok={ok_write and valid_ok and not bool(_he_q)} reason={reason} "
                f"final: H={post_metrics['h_count']} formal_charge={post_metrics['formal_charge']} "
                f"has_partial_charges={post_metrics['has_partial_charges']} "
                f"ad4_types_ok={post_metrics['ad4_types_ok']}"
            )

            try:
                _append_prep_status(
                    status_log,
                    ligand_name=lig_id,
                    status=(
                        "OK" if (ok_write and valid_ok and not bool(_he_q)) else "FAIL"
                    ),
                    reason=reason or "",
                    relpath=str(pdbqt_path.name),
                    stage="write_pdbqt",
                    failure_code=("helium_rules" if _he_q else ""),
                    failure_detail=str(_he_q or ""),
                    fixes_count=int(_he_fixes or 0),
                    rules_ver=rules_version(),
                    writer_final=str(chosen_writer),
                    rescue_used=(",".join(fallbacks_used) if fallbacks_used else ""),
                    torsion_root_rule="",
                    polarH="",
                    ad4_types_ok=("yes" if post_metrics.get("ad4_types_ok") else "no"),
                    charges_ok=(
                        "yes" if post_metrics.get("has_partial_charges") else "no"
                    ),
                    torsdof="",
                )
            except Exception as _e:
                logging.warning("[status-log] append failed for %s: %s", lig_id, _e)

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

        try:
            if pdbqt_path.exists() and pdbqt_path.stat().st_size > 100:
                tmp_mol2.unlink(missing_ok=True)
        except Exception:
            pass

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
