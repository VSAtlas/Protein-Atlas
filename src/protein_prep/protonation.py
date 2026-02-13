"""Reduce-based protonation and fallback orchestration."""

from __future__ import annotations

import logging
import os
import shutil
from pathlib import Path
from typing import Optional, Tuple, Union

import protein_prep.prep_utils as prep_utils
from pdb_fixer import fix_pdb_elements
from propka_wire import pdb2pqr_protonate

from protein_prep.prep_utils import _cfg, _resolve_variant_token
from protein_prep.tool_runners import run_openbabel_add_h, run_subprocess_capture


PHENIX_DIR = _cfg("PHENIX_DIR", "", "phenix_dir")
OPENBABEL_PATH = os.environ.get("OPENBABEL_PATH") or shutil.which("obabel") or ""
_REPO_ROOT = Path(__file__).resolve().parents[2]


def _pick_reduce_exe() -> str:
    """
    Choose the actual 'reduce' binary robustly.
    Precedence:
      1) $REDUCE_EXE (env) or config value
      2) Repo-relative local build: <repo>/tools/reduce/reduce_src/reduce
      3) Phenix conda_base/bin/reduce (next to PHENIX_DIR)
      4) reduce on PATH
    """
    # 1) Explicit env/config
    explicit = (
        os.environ.get("REDUCE_EXE")
        or os.environ.get("REDUCE_BIN")
        or _cfg("REDUCE_EXE", "", "reduce_exe")
    )
    if explicit and Path(explicit).exists():
        return explicit

    # 2) Repo-relative (portable)
    repo_root = _REPO_ROOT / "tools" / "reduce" / "reduce_src" / "reduce"
    if repo_root.exists() and os.access(str(repo_root), os.X_OK):
        return str(repo_root)

    # 3) Phenix conda-base candidate near PHENIX_DIR
    if PHENIX_DIR:
        phenix_root = Path(PHENIX_DIR).resolve().parent
        cb = phenix_root / "conda_base" / "bin" / "reduce"
        if cb.exists() and os.access(str(cb), os.X_OK):
            return str(cb)

    # 4) PATH fallback
    which = shutil.which("reduce") or shutil.which("reduce.exe")
    return which or "reduce"


REDUCE_EXE = _pick_reduce_exe()
logging.info("Using Reduce at: %s", REDUCE_EXE)

# Default HET dict (repo-relative) if not provided by env/config
DEFAULT_HET = _REPO_ROOT / "tools" / "reduce" / "reduce_wwPDB_het_dict.txt"


def _het_dict_path() -> str | None:
    hd = os.environ.get("REDUCE_HET_DICT") or _cfg("REDUCE_HET_DICT", "")
    if hd and Path(hd).is_file():
        return hd
    return str(DEFAULT_HET) if DEFAULT_HET.is_file() else None


def hydrogenation_status(pdb_path: Union[str, Path]) -> Tuple[str, int, int, float]:
    h = heavy = 0
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            el = line[76:78].strip().upper()
            if el == "H":
                h += 1
            else:
                heavy += 1
    ratio = h / max(heavy, 1)
    if h == 0:
        return "NO_H", h, heavy, ratio
    if ratio < 0.20:
        return "SUSPECT_LOW_H", h, heavy, ratio
    return "HAS_H", h, heavy, ratio


def conect_coverage(pdb_path: Union[str, Path]) -> float:
    atom_ids, conect_ids = set(), set()
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_ids.add(line[6:11].strip())
            elif line.startswith("CONECT"):
                parts = line.split()
                conect_ids.update(parts[1:])
    return len(conect_ids & atom_ids) / max(len(atom_ids), 1)


def _maybe_get_target_ph() -> float | None:
    # Priority: explicit env, then context_ph json drop, else None.
    v = os.environ.get("TARGET_PH", "").strip()
    if v:
        try:
            return float(v)
        except Exception:
            pass
    # allow a sidecar json dumped by context_ph: <pdb_id>.ph.json with {"target_pH": 7.8}
    input_pdb_path = (os.environ.get("INPUT_PDB_PATH") or "").strip()
    if input_pdb_path:
        try:
            pdb_base = os.path.splitext(os.path.basename(input_pdb_path))[0]
            sidecar = Path(input_pdb_path).parent / f"{pdb_base}.ph.json"
            if sidecar.exists():
                import json

                data = json.loads(sidecar.read_text())
                t = data.get("target_pH", None)
                if isinstance(t, (int, float)):
                    return float(t)
        except Exception:
            pass
    return None


def _protonate_with_pdb2pqr_if_available(
    nolig_pdb_path: str,
    out_dir: Path,
    logger,
    *,
    variant: Optional[str] = None,
) -> tuple[Path, bool, str | None]:
    """
    Try PDB2PQR at the *pipeline pH* if available; fall back to the input PDB.
    Returns: (pdb_for_reduce, used_pdb2pqr, propka_log_path_or_None)
    """
    # Choose a pH: use context_ph if you already resolved one before calling this,
    # otherwise default to 7.0 (harmless; you can feed in your target later).
    try:
        from path_router.context_ph import select_ph_values_for_protonation

        phs = select_ph_values_for_protonation(nolig_pdb_path)
        target_ph = float(phs[0]) if phs else 7.0
    except Exception:
        target_ph = 7.0

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    variant_token = (variant or "").strip().upper() if variant else None
    if variant_token not in {"APO", "HOLO"}:
        variant_token = _resolve_variant_token(prep_utils.config)

    pdb_from_p2p, pk_log = pdb2pqr_protonate(
        nolig_pdb_path,
        target_ph,
        out_dir,
        variant=variant_token,
        cfg=prep_utils.config,
    )
    if pdb_from_p2p:
        logger.info("PDB2PQR succeeded at pH %.2f -> %s", target_ph, pdb_from_p2p)
        return Path(pdb_from_p2p), True, pk_log
    else:
        logger.warning("PDB2PQR unavailable/failed; Reduce will build hydrogens.")
        return Path(nolig_pdb_path), False, None


def assign_protonation_states(
    input_pdb: Union[str, Path],
    output_pdb: Union[str, Path],
    reduce_exe: Optional[str] = None,
) -> str:
    import uuid
    from pathlib import Path

    # --- toggles ---
    _REDUCE_RETRY = str(os.environ.get("REDUCE_RETRY", "1")).lower() not in {
        "0",
        "false",
        "no",
    }
    _FALLBACK_ADDH = str(os.environ.get("FALLBACK_ADDH", "1")).lower() not in {
        "0",
        "false",
        "no",
    }

    input_pdb = os.path.abspath(str(input_pdb)).replace("\\", "/")
    output_pdb = os.path.abspath(str(output_pdb)).replace("\\", "/")

    # Optional precheck: very low CONECT coverage is a strong predictor of Reduce complaints.
    try:
        cc = conect_coverage(input_pdb)
        if cc < 0.05:
            logging.info("[reduce precheck] low_conect=%.2f file=%s", cc, input_pdb)
            try:
                # Idempotent, YAML-driven element repair
                fix_pdb_elements(input_pdb)
            except Exception as _e:
                logging.warning("[reduce precheck] elemfix_skip note=%s", _e)
    except Exception as _e:
        logging.warning("[reduce precheck] coverage_skip note=%s", _e)

    has_h = hydrogenation_status(input_pdb)[0] != "NO_H"
    # If input already has H (e.g., from PDB2PQR), do flip/cleanup only; else do a full H build.
    reduce_flags = ["-FLIP", "-Quiet"] if has_h else ["-BUILD", "-Quiet"]

    exe = reduce_exe
    exe_dir = os.path.dirname(exe) if exe else None

    def run_reduce(in_pdb: str, stage_name: str) -> str:
        env = dict(os.environ)
        het = _het_dict_path()
        if het:
            env["REDUCE_HET_DICT"] = het
        cp = run_subprocess_capture([exe] + reduce_flags + [in_pdb], env=env)
        with open(output_pdb, "w", encoding="utf-8") as out:
            out.write(cp.stdout or "")
        # Persist stderr and count atoms written
        try:
            (Path(output_pdb).parent / f"{stage_name}.stderr.txt").write_text(
                cp.stderr or "", encoding="utf-8"
            )
        except Exception:
            pass
        wrote_atoms = 0
        try:
            with open(output_pdb, "r", encoding="utf-8", errors="ignore") as fh:
                wrote_atoms = sum(1 for ln in fh if ln.startswith(("ATOM  ", "HETATM")))
        except Exception:
            wrote_atoms = 0
        logging.warning(
            "[reduce] stage=%s rc=%s flags=%s wrote_atoms=%d",
            stage_name,
            cp.returncode,
            " ".join(reduce_flags),
            wrote_atoms,
        )
        return output_pdb

    status_before, h0, hv0, r0 = hydrogenation_status(input_pdb)
    logging.info(
        "[H-Scan before] %s (H=%d, Heavy=%d, H/Heavy=%.2f)", status_before, h0, hv0, r0
    )
    if exe is None:
        # Caller requested to skip Reduce (e.g., nucleotides). Go straight to fallback.
        if _FALLBACK_ADDH:
            logging.info(
                "[protonate] skipping Reduce by request; using OpenBabel fallback"
            )
            run_openbabel_add_h(input_pdb, output_pdb, openbabel_path=OPENBABEL_PATH)
        else:
            # Keep as-is if fallback is disabled
            shutil.copy(input_pdb, output_pdb)
            logging.warning(
                "[protonate] Reduce skipped and FALLBACK_ADDH=0; copying input"
            )
        # Continue into the existing post-guard checks (size/H-Scan) below
        sz = Path(output_pdb).stat().st_size if Path(output_pdb).exists() else 0
        if sz == 0:
            raise RuntimeError("protonation_empty_output")

        #  strict ATOM/HETATM guard (Reduce can write tiny, header-only files)
        def _file_has_atoms(p: str) -> bool:
            try:
                with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                    return any(ln.startswith(("ATOM  ", "HETATM")) for ln in fh)
            except Exception:
                return False

        if not _file_has_atoms(output_pdb):
            logging.error(
                "[protonate] wrote_atoms=0 after Reduce; attempting OpenBabel fallback"
            )
            try:
                tmp_babel = output_pdb + ".babel.pdb"
                run_openbabel_add_h(
                    input_pdb, tmp_babel, openbabel_path=OPENBABEL_PATH
                )
                shutil.move(tmp_babel, output_pdb)
                logging.warning("[fallback] openbabel applied (post-Reduce)")
            except Exception as e:
                logging.error(
                    "[fallback] openbabel failed: %s; copying input→output",
                    str(e)[:200],
                )
                shutil.copy(input_pdb, output_pdb)
        status_after, h1, hv1, r1 = hydrogenation_status(output_pdb)
        logging.info(
            "[H-Scan after ] %s (H=%d, Heavy=%d, H/Heavy=%.2f)",
            status_after,
            h1,
            hv1,
            r1,
        )
        return output_pdb

    # --- First attempt (default flags) ---
    try:
        run_reduce(input_pdb, "Reduce#1")
    except Exception as e:
        logging.warning(
            "[reduce#1 fail] exe=%s het_dict=%s has_h=%s note=%s",
            exe,
            os.environ.get("REDUCE_HET_DICT") or _cfg("REDUCE_HET_DICT", ""),
            str(has_h),
            str(e)[:200],
        )
        completed = False
        # Quick element repair before a retry
        try:
            fix_pdb_elements(input_pdb)
        except Exception:
            pass

        # ---  retry ladder ---
        retried_ok = False
        if _REDUCE_RETRY:
            try:
                # Second attempt: add -noflip (Reduce sometimes flips/complains on tricky H networks)
                reduce_flags_noflip = (
                    (["-noflip"] + reduce_flags)
                    if "-noflip" not in reduce_flags
                    else reduce_flags
                )

                def run_reduce_noflip(in_pdb: str, stage_name: str) -> str:
                    env = dict(os.environ)
                    het = env.get("REDUCE_HET_DICT") or _cfg("REDUCE_HET_DICT", "")
                    if het:
                        env["REDUCE_HET_DICT"] = het
                    cp = run_subprocess_capture(
                        [exe] + reduce_flags_noflip + [in_pdb],
                        cwd=exe_dir,
                        env=env,
                    )
                    with open(output_pdb, "w", encoding="utf-8") as out:
                        out.write(cp.stdout or "")
                    (logging.info if cp.returncode == 0 else logging.warning)(
                        "[reduce] stage=%s rc=%s flags=%r in=%s",
                        stage_name,
                        cp.returncode,
                        reduce_flags_noflip,
                        in_pdb,
                    )
                    if cp.returncode != 0:
                        raise RuntimeError(
                            f"{stage_name} reduce failed: {cp.stderr.strip()}"
                        )
                    return output_pdb

                run_reduce_noflip(input_pdb, "Reduce#2(noflip)")
                retried_ok = True
                completed = True
            except Exception as e2:
                logging.warning("[reduce#2 fail] noflip note=%s", str(e2)[:200])

        if not retried_ok:
            # --- H-add fallback path (OpenBabel) ---
            if _FALLBACK_ADDH:
                try:
                    # Guard element columns first to avoid He/column drift in fallback
                    try:
                        fix_pdb_elements(input_pdb)
                    except Exception as _e_fix:
                        logging.warning(
                            "[element] guard before fallback failed: %s", _e_fix
                        )
                    logging.info(
                        "[protonate fallback] using=OpenBabel in=%s out=%s",
                        input_pdb,
                        output_pdb,
                    )
                    run_openbabel_add_h(
                        input_pdb, output_pdb, openbabel_path=OPENBABEL_PATH
                    )
                    logging.warning("[reduce] fallback H-add applied")
                    completed = True

                except Exception as babel_error:
                    logging.error("OpenBabel fallback failed: %s", babel_error)
                    # Last resort: if input already had H, keep them; else raise
                    if has_h:
                        shutil.copy(input_pdb, output_pdb)
                        logging.warning(
                            "Reduce+fallback failed; keeping existing hydrogens."
                        )
                    else:
                        raise
        if not completed:
            try:
                # Repair element columns (YAML-driven fixer under the hood)
                fix_pdb_elements(input_pdb)
            except Exception:
                pass
            try:
                tmp_in = (
                    os.path.splitext(input_pdb)[0] + f"_retry_{uuid.uuid4().hex}.pdb"
                )
                shutil.copy(input_pdb, tmp_in)
                try:
                    run_reduce(tmp_in, "Reduce#3(temp)")
                finally:
                    try:
                        os.remove(tmp_in)
                    except Exception:
                        pass
            except Exception as e2:
                logging.error("Reduce temp retry failed: %s", e2)
                try:
                    if has_h:
                        shutil.copy(input_pdb, output_pdb)
                        logging.warning(
                            "Reduce failed; keeping existing hydrogens (no rebuild)."
                        )
                    else:
                        logging.info(
                            "[protonate fallback] using=OpenBabel in=%s out=%s",
                            input_pdb,
                            output_pdb,
                        )
                        run_openbabel_add_h(
                            input_pdb, output_pdb, openbabel_path=OPENBABEL_PATH
                        )
                        logging.info("Open Babel used to add hydrogens.")
                except Exception as babel_error:
                    logging.error("OpenBabel fallback failed: %s", babel_error)
                    shutil.copy(input_pdb, output_pdb)
                    logging.warning("Hydrogenation skipped; copied input to output.")

    # Guard: the path we will scan must exist and be non-empty
    sz = Path(output_pdb).stat().st_size if Path(output_pdb).exists() else 0
    if sz == 0:
        raise RuntimeError("protonation_empty_output")

    status_after, h1, hv1, r1 = hydrogenation_status(output_pdb)
    logging.info(
        "[H-Scan after ] %s (H=%d, Heavy=%d, H/Heavy=%.2f)", status_after, h1, hv1, r1
    )

    if status_after == "NO_H":
        logging.warning("Reduce produced no hydrogens; trying Open Babel fallback.")
        try:
            logging.info(
                "[protonate fallback] using=OpenBabel in=%s out=%s",
                input_pdb,
                output_pdb,
            )
            tmp_babel = output_pdb + ".babel.pdb"
            run_openbabel_add_h(input_pdb, tmp_babel, openbabel_path=OPENBABEL_PATH)
            shutil.move(tmp_babel, output_pdb)
            status_after, h1, hv1, r1 = hydrogenation_status(output_pdb)
            logging.info(
                "[H-Scan babel] %s (H=%d, Heavy=%d, H/Heavy=%.2f)",
                status_after,
                h1,
                hv1,
                r1,
            )
        except Exception as e:
            logging.error("OpenBabel fallback failed after Reduce: %s", e)

        if hydrogenation_status(output_pdb)[0] == "NO_H":
            if str(os.environ.get("ALLOW_NO_HYDROGENS", "")).lower() in (
                "1",
                "true",
                "yes",
            ):
                logging.warning(
                    "Continuing with NO_H due to ALLOW_NO_HYDROGENS env override."
                )
                # Still emit result line for grep
                try:
                    sz = (
                        Path(output_pdb).stat().st_size
                        if Path(output_pdb).exists()
                        else 0
                    )
                except Exception:
                    sz = 0
                logging.info(
                    "[protonate result] status=%s H=%d Heavy=%d ratio=%.2f size=%d",
                    "NO_H",
                    h1,
                    hv1,
                    r1,
                    sz,
                )
                return output_pdb
            raise RuntimeError("Protonation produced no hydrogens.")

    # Final greppable result line
    try:
        sz = Path(output_pdb).stat().st_size if Path(output_pdb).exists() else 0
    except Exception:
        sz = 0
    logging.info(
        "[protonate result] status=%s H=%d Heavy=%d ratio=%.2f size=%d",
        status_after,
        h1,
        hv1,
        r1,
        sz,
    )

    return output_pdb
