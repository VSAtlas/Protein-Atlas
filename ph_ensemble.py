# ---pH-ensemble A→B→C orchestration ---

import os
import json
import argparse
import logging
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, Tuple, Optional
import hashlib
from propka_wire import apply_propka_states
import automate_protein_prep
from path_router import ph_ensemble_dir
from input_and_export_functions import load_config

elog = logging.getLogger("ph_ensemble")
if not elog.handlers:
    _h = logging.StreamHandler()
    _h.setFormatter(logging.Formatter("[ph_ensemble] %(message)s"))
    elog.addHandler(_h)
elog.setLevel(logging.INFO)
elog.propagate = False


# --- minimal helpers (no external deps) ---
def _vina_receptor_pdbqt_compatible(path: str) -> bool:
    """
    Conservative text check for receptor PDBQT compatibility with Vina.
    Reject files that begin with unsupported records (e.g., COMPND/AUTHOR).
    """
    allowed_prefixes = (
        "REMARK",
        "ROOT",
        "ENDROOT",
        "BRANCH",
        "ENDBRANCH",
        "ATOM",
        "HETATM",
        "TER",
        "MODEL",
        "ENDMDL",
        "END",
    )
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            seen_atom = False
            for i, ln in enumerate(fh):
                s = ln.strip()
                if not s:
                    continue
                if not s.startswith(allowed_prefixes):
                    return False
                if s.startswith(("ATOM", "HETATM")):
                    seen_atom = True
                if i >= 500:
                    break
        return seen_atom
    except Exception:
        return False


def _sha1_of_file(path: str) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _ph_token(ph: float) -> str:
    # repo convention: 6.0 -> 6_0, 7.4 -> 7_4
    return str(ph).replace(".", "_")


def _count_atoms_pdb(path: str) -> int:
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            return sum(1 for ln in fh if ln.startswith(("ATOM  ", "HETATM")))
    except Exception:
        return 0


def _count_pocket_hydrogens(
    pdb_path: str, center: Tuple[float, float, float], radius: float
) -> int:
    cx, cy, cz = center
    r2 = radius * radius
    nH = 0
    try:
        with open(pdb_path, "r", errors="ignore") as fh:
            for ln in fh:
                if not ln.startswith(("ATOM  ", "HETATM")):
                    continue
                try:
                    x = float(ln[30:38])
                    y = float(ln[38:46])
                    z = float(ln[46:54])
                except Exception:
                    continue
                dx, dy, dz = x - cx, y - cy, z - cz
                if (dx * dx + dy * dy + dz * dz) > r2:
                    continue
                if ln[12:16].strip().upper().startswith("H"):
                    nH += 1
    except Exception:
        pass
    return nH


def _pocket_name_counts(
    pdb_path: str, center: Tuple[float, float, float], radius: float
) -> dict:
    cx, cy, cz = center
    r2 = radius * radius
    counts = {
        "ASP": 0,
        "ASH": 0,
        "GLU": 0,
        "GLH": 0,
        "HIS": 0,
        "HID": 0,
        "HIE": 0,
        "HIP": 0,
        "LYS": 0,
        "LYN": 0,
        "CYS": 0,
        "CYM": 0,
        "TYR": 0,
    }
    seen = set()
    try:
        with open(pdb_path, "r", errors="ignore") as fh:
            for ln in fh:
                if not ln.startswith(("ATOM  ", "HETATM")):
                    continue
                try:
                    x = float(ln[30:38])
                    y = float(ln[38:46])
                    z = float(ln[46:54])
                except Exception:
                    continue
                dx, dy, dz = x - cx, y - cy, z - cz
                if (dx * dx + dy * dy + dz * dz) > r2:
                    continue
                chain = ln[21].strip() or " "
                resi = int(ln[22:26])
                key = (chain, resi)
                if key in seen:
                    continue
                seen.add(key)
                resn = ln[17:20].strip().upper()
                if resn in counts:
                    counts[resn] += 1
    except Exception:
        pass
    return {k: v for k, v in counts.items() if v > 0}


def _ph_member_job(job):
    """ProcessPool worker: run A→B→C for a single pH member."""
    (
        idx,
        ph,
        tag,
        cleaned_receptor_pdb,
        ensemble_dir_str,
        center,
        radius,
        reduce_exe,
    ) = job

    ensemble_dir = Path(ensemble_dir_str)

    prestate_pdb, pka_path, _n = apply_propka_states(
        cleaned_receptor_pdb=cleaned_receptor_pdb,
        center=center,
        radius=radius,
        ph=ph,
        out_dir=str(ensemble_dir),
        tag=tag,
    )

    pre_sha = _sha1_of_file(prestate_pdb)
    pre_atoms = _count_atoms_pdb(prestate_pdb)
    pre_bins = _pocket_name_counts(prestate_pdb, center, radius)

    withH_pdb = str(ensemble_dir / f"{tag}.withH.pdb")
    automate_protein_prep.assign_protonation_states(
        input_pdb=prestate_pdb,
        output_pdb=withH_pdb,
        reduce_exe=reduce_exe,
    )

    with_sha = _sha1_of_file(withH_pdb)
    with_atoms = _count_atoms_pdb(withH_pdb)
    with_bins = _pocket_name_counts(withH_pdb, center, radius)
    pocket_H = _count_pocket_hydrogens(withH_pdb, center, radius)

    pdbqt_path = str(ensemble_dir / f"{tag}.pdbqt")
    ok = automate_protein_prep.run_prepare_receptor(
        input_pdb=withH_pdb,
        output_pdbqt=pdbqt_path,
        cfg=automate_protein_prep.config,
    )
    compatible = _vina_receptor_pdbqt_compatible(pdbqt_path)
    if (not ok) or (not compatible):
        base_receptor_pdbqt = Path(cleaned_receptor_pdb).with_name(
            Path(cleaned_receptor_pdb).stem.replace("_cleaned", "") + ".pdbqt"
        )
        if base_receptor_pdbqt.exists() and _vina_receptor_pdbqt_compatible(
            str(base_receptor_pdbqt)
        ):
            shutil.copy2(base_receptor_pdbqt, pdbqt_path)
            elog.warning(
                "[stage.C.pdbqt] fallback=base_receptor ph=%.2f reason=%s src=%s dst=%s",
                ph,
                "prepare_failed" if not ok else "incompatible_output",
                base_receptor_pdbqt,
                pdbqt_path,
            )
        else:
            shutil.copy2(withH_pdb, pdbqt_path)
            elog.warning(
                "[stage.C.pdbqt] fallback=withH_copy ph=%.2f reason=%s src=%s dst=%s",
                ph,
                "prepare_failed" if not ok else "incompatible_output",
                withH_pdb,
                pdbqt_path,
            )

    try:
        pdbqt_size = Path(pdbqt_path).stat().st_size
    except Exception:
        pdbqt_size = -1

    return {
        "idx": idx,
        "ph": ph,
        "tag": tag,
        "prestate": prestate_pdb,
        "pka": pka_path,
        "pre_sha": pre_sha,
        "pre_atoms": pre_atoms,
        "pre_bins": pre_bins,
        "withH": withH_pdb,
        "withH_sha1": with_sha,
        "with_atoms": with_atoms,
        "with_bins": with_bins,
        "pocket_H": pocket_H,
        "pdbqt": pdbqt_path,
        "pdbqt_size": pdbqt_size,
    }


def build_ph_ensemble(
    pdb_id: str,
    cleaned_receptor_pdb: str,
    out_dir: str,
    center: Tuple[float, float, float],
    radius: float,
    ph_values: Iterable[float],
    member_index_start: int = 0,
    variant: Optional[str] = None,
    legacy: bool = False,
) -> str:
    """
    For each pH value: A) PROPKA-local renames -> prestate
                       B) single hydrogenation via automate_protein_prep.assign_protonation_states
                       C) app.run_prepare_receptor -> pdbqt
    Writes ensemble.json and returns its path.

    When `variant` is APO/HOLO (and legacy=False), artifacts are written under
    processed_pdbs/<PDB>/<VARIANT>/receptor/ph_ensemble/.
    Legacy mode (variant=None or legacy=True) keeps processed_pdbs/<PDB>/receptor/ph_ensemble/.
    """
    out_root = Path(out_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    variant_token = (variant or "").strip().upper() or None
    if legacy:
        variant_token = None

    ensemble_dir = None
    try:
        candidate = ph_ensemble_dir(pdb_id, variant=variant_token, legacy=legacy)
        candidate.relative_to(out_root)
    except Exception:
        ensemble_dir = None
    else:
        ensemble_dir = candidate

    if ensemble_dir is None:
        ensemble_dir = out_root / pdb_id
        if variant_token:
            ensemble_dir = ensemble_dir / variant_token
        ensemble_dir = ensemble_dir / "receptor" / "ph_ensemble"

    ensemble_dir.mkdir(parents=True, exist_ok=True)
    tag_root = str(pdb_id).replace("/", "_").replace("\\", "_")

    # [ADD] optional file logging toggle: env pHlogs -> config.txt pHlogs -> default False
    def _truthy(x):
        return str(x).strip().lower() in {"1", "true", "yes", "on", "y"}

    enable_file_log = False
    if "pHlogs" in os.environ:
        enable_file_log = _truthy(os.environ["pHlogs"])
    else:
        try:
            root = Path(__file__).resolve().parent
            cfg = load_config(config_path=str(root / "config.txt"), base_dir=root)
            enable_file_log = _truthy(cfg.get("pHlogs", cfg.get("PHLOGS", "")))
        except Exception:
            enable_file_log = False
    if enable_file_log:
        fh = logging.FileHandler(
            str(ensemble_dir / "pH_ensemble.log"), mode="w", encoding="utf-8"
        )
        fh.setFormatter(logging.Formatter("[%(name)s] %(message)s"))
        # capture this module + partner modules without changing their code
        for nm in (
            "ph_ensemble",
            "propka_wire",
            "automate_protein_prep",
            "proteinprep",
        ):
            lg = logging.getLogger(nm)
            lg.addHandler(fh)
            lg.setLevel(logging.DEBUG)
            lg.propagate = False
        elog.info(
            "[pHlogs] enabled → %s", str((ensemble_dir / "pH_ensemble.log").resolve())
        )
    else:
        elog.info("[pHlogs] disabled")

    ph_values = list(ph_values)
    n_ph = len(ph_values)
    elog.info("[ph.list] n=%d values=%s", n_ph, ph_values)

    jobs = []
    for idx, ph in enumerate(ph_values):
        tag = f"{tag_root}_pH{str(ph).replace('.', '_')}"
        elog.info("[ph.member.start] ph=%.2f tag=%s", ph, tag)
        elog.info(
            "[ph.propka.call] ph=%.2f center=(%.3f,%.3f,%.3f) r=%.2f",
            ph,
            center[0],
            center[1],
            center[2],
            radius,
        )
        jobs.append(
            (
                idx,
                float(ph),
                tag,
                cleaned_receptor_pdb,
                str(ensemble_dir),
                center,
                float(radius),
                automate_protein_prep.REDUCE_EXE,
            )
        )

    cpu = int(os.getenv("CPU", os.cpu_count() or 1))
    workers = min(n_ph, cpu)
    if workers < 1:
        workers = 1
    threads_per_job = 1

    results = []
    if n_ph <= 1 or workers <= 1:
        if n_ph == 1 and jobs:
            elog.info(
                "[ph.serial] single pH member; running without parallel executor. tag=%s ph=%.2f",
                jobs[0][2],
                jobs[0][1],
            )
        elif n_ph > 1:
            elog.info(
                "[ph.serial] workers=1; running %d pH members without parallel executor",
                n_ph,
            )
        for job in jobs:
            results.append(_ph_member_job(job))
    else:
        elog.info(
            "[ph.parallel] CPU=%d n_ph=%d workers=%d threads_per_job=%d total_threads=%d",
            cpu,
            n_ph,
            workers,
            threads_per_job,
            workers * threads_per_job,
        )
        try:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                future_to_job = {pool.submit(_ph_member_job, job): job for job in jobs}
                for fut in as_completed(future_to_job):
                    job = future_to_job[fut]
                    try:
                        res = fut.result()
                    except Exception as exc:
                        for other in future_to_job:
                            if other is not fut:
                                other.cancel()
                        raise RuntimeError(
                            f"pH member failed ph={job[1]:.2f} tag={job[2]}"
                        ) from exc
                    results.append(res)
        except Exception as exc:
            elog.warning(
                "[ph.parallel] unavailable -> serial fallback reason=%s",
                exc,
            )
            results = []
            for job in jobs:
                results.append(_ph_member_job(job))

    members = []
    mi = member_index_start
    _prev_pre_sha = None
    _prev_with_sha = None
    _prev_pdbqt_size = None

    def _fmt_bins(bin_dict):
        return " ".join(f"{k}={v}" for k, v in bin_dict.items())

    for res in sorted(results, key=lambda d: d["idx"]):
        ph = res["ph"]
        tag = res["tag"]
        prestate_pdb = res["prestate"]
        pka_path = res["pka"]
        pre_sha = res["pre_sha"]
        pre_atoms = res["pre_atoms"]
        pre_bins = res["pre_bins"]
        withH_pdb = res["withH"]
        with_sha = res["withH_sha1"]
        with_atoms = res["with_atoms"]
        with_bins = res["with_bins"]
        pocket_H = res["pocket_H"]
        pdbqt_path = res["pdbqt"]
        pdbqt_size = res["pdbqt_size"]

        elog.info(
            "[stage.A.prestate] ph=%.2f atoms=%d sha1=%s bins=%s",
            ph,
            pre_atoms,
            pre_sha,
            _fmt_bins(pre_bins),
        )
        if _prev_pre_sha is not None:
            elog.info(
                "[compare.prestate] ph=%.2f same_as_prev=%s",
                ph,
                str(pre_sha == _prev_pre_sha),
            )
        _prev_pre_sha = pre_sha

        elog.info("[withH.sha1] ph=%.2f sha1=%s", ph, with_sha)
        elog.info(
            "[stage.B.withH] ph=%.2f atoms=%d sha1=%s pocket_H=%d bins=%s",
            ph,
            with_atoms,
            with_sha,
            pocket_H,
            _fmt_bins(with_bins),
        )
        if _prev_with_sha is not None:
            elog.info(
                "[compare.withH] ph=%.2f same_as_prev=%s",
                ph,
                str(with_sha == _prev_with_sha),
            )
        _prev_with_sha = with_sha

        elog.info("[stage.C.pdbqt] ph=%.2f size=%d path=%s", ph, pdbqt_size, pdbqt_path)
        if _prev_pdbqt_size is not None:
            elog.info(
                "[compare.pdbqt] ph=%.2f same_size_as_prev=%s",
                ph,
                str(pdbqt_size == _prev_pdbqt_size),
            )
        _prev_pdbqt_size = pdbqt_size
        elog.info(
            "[ensemble.step] ph=%.2f prestate=%s withH=%s pdbqt=%s",
            ph,
            prestate_pdb,
            withH_pdb,
            pdbqt_path,
        )

        members.append(
            {
                "ph": ph,
                "tag": tag,
                "prestate": prestate_pdb,
                "pka": pka_path,
                "withH": withH_pdb,
                "pdbqt": pdbqt_path,
                "withH_sha1": with_sha,
            }
        )
        mi += 1

    # --- POST-PASS: group by withH SHA-1 and collapse identical members ---
    # group members by chemical identity
    groups = {}
    for m in members:
        groups.setdefault(m.get("withH_sha1", ""), []).append(m)

    dedup_members = []
    canonical_phs = []

    for sha, grp in groups.items():
        # canonical = lowest numeric pH
        grp_sorted = sorted(grp, key=lambda d: float(d["ph"]))
        canonical = grp_sorted[0]
        phs = sorted(float(x["ph"]) for x in grp_sorted)
        tokens = [_ph_token(p) for p in phs]

        # target collapsed filename (multi-pH) for canonical
        collapsed_name = f"{tag_root}_pH{'+'.join(tokens)}.pdbqt"
        collapsed_path = Path(canonical["pdbqt"]).parent / collapsed_name

        # rename only if there is more than one in the equivalence class
        out_pdbqt_path = Path(canonical["pdbqt"])
        if len(grp_sorted) > 1:
            # ensure we don't overwrite an existing file
            if (
                collapsed_path.exists()
                and collapsed_path.resolve() != out_pdbqt_path.resolve()
            ):
                # -dupN suffix
                i = 1
                while True:
                    candidate = collapsed_path.with_name(
                        f"{collapsed_path.stem}-dup{i}{collapsed_path.suffix}"
                    )
                    if not candidate.exists():
                        collapsed_path = candidate
                        break
                    i += 1
            # perform rename if needed
            if out_pdbqt_path.resolve() != collapsed_path.resolve():
                elog.info(
                    "[dedupe.rename] %s → %s", str(out_pdbqt_path), str(collapsed_path)
                )
                out_pdbqt_path.rename(collapsed_path)
                out_pdbqt_path = collapsed_path
            # delete all artifacts for non-canonical members: .pdbqt, .withH.pdb, .prestate.pdb, .pka
            for dup in grp_sorted[1:]:
                for key, label in (
                    ("pdbqt", "pdbqt"),
                    ("withH", "withH"),
                    ("prestate", "prestate"),
                    ("pka", "pka"),
                ):
                    path = dup.get(key)
                    if not path:
                        continue
                    p = Path(path)
                    try:
                        p.unlink()
                        elog.info("[dedupe.delete.%s] %s", label, str(p))
                    except FileNotFoundError:
                        elog.warning("[dedupe.delete.missing.%s] %s", label, str(p))
                    except Exception as e:
                        elog.warning(
                            "[dedupe.delete.error.%s] %s (%s)", label, str(p), e
                        )

            elog.info(
                "[dedupe.group] sha=%s phs=%s canonical=%.2f",
                sha[:12],
                ",".join(f"{p:.2f}" for p in phs),
                float(canonical["ph"]),
            )
        else:
            elog.info(
                "[dedupe.group] sha=%s phs=%s (no collapse)", sha[:12], f"{phs[0]:.2f}"
            )

        # manifest entry for this group (canonical metadata + ph_equiv + final pdbqt path)
        dedup_members.append(
            {
                "ph": float(canonical["ph"]),
                "ph_equiv": phs,
                "tag": canonical["tag"],
                "prestate": canonical["prestate"],
                "pka": canonical["pka"],
                "withH": canonical["withH"],
                "withH_sha1": canonical.get("withH_sha1", ""),
                "pdbqt": str(out_pdbqt_path),
            }
        )
        canonical_phs.append(float(canonical["ph"]))

    # sort manifest by canonical pH
    dedup_members = sorted(dedup_members, key=lambda d: d["ph"])
    canonical_phs = sorted(canonical_phs)

    manifest = {
        "pdb_id": pdb_id,
        "ensemble": canonical_phs,  # << one per group
        "members": dedup_members,  # << one per group
    }

    manifest_path = ensemble_dir / "ensemble.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    elog.info("[manifest.write] %s", str(manifest_path))

    return str(manifest_path)


# --- END: pH-ensemble A→B→C orchestration ---


def main():
    ap = argparse.ArgumentParser(description="Build pH-aware receptor ensemble.")
    ap.add_argument("pdb_path", help="Input receptor PDB (ligand-free).")
    ap.add_argument("--output-root", default="processed_pdbs")
    ap.add_argument(
        "--force", action="store_true", help="Force re-prep of base receptor."
    )
    ap.add_argument(
        "--center",
        type=str,
        default=None,
        help="Binding-site center as 'x,y,z' (required if --scope=pocket).",
    )
    ap.add_argument(
        "--radius",
        type=float,
        default=10.0,
        help="Sphere radius in Å for local titration (used when --scope=pocket).",
    )
    ap.add_argument(
        "--scope",
        choices=["global", "pocket"],
        default="global",
        help="Renaming scope: 'global' (entire protein; default) or 'pocket' (within --radius).",
    )
    ap.add_argument(
        "--ph-list",
        type=str,
        default="7.4",
        help="Comma-separated pH values, e.g., 5,7,9",
    )
    args = ap.parse_args()

    # Parse pH list (always required to be valid)
    try:
        ph_values = [float(x) for x in args.ph_list.split(",") if x.strip() != ""]
    except Exception:
        raise SystemExit(
            f"--ph-list must be comma-separated numbers, got {args.ph_list!r}"
        )

    if args.force:
        os.environ["FORCE_REPROCESS"] = "1"

    # Single clean step (do not edit automate_protein_prep)
    cleaned_pdb, _receptor_pdbqt = automate_protein_prep.main(
        args.pdb_path, args.output_root
    )
    pdb_id = Path(args.pdb_path).stem[:4].upper()

    # Scope + center/radius policy
    if args.scope == "pocket":
        if not args.center:
            raise SystemExit("--scope=pocket requires --center='x,y,z'.")
        # Parse center only in pocket mode
        try:
            cx, cy, cz = (float(t) for t in args.center.replace(",", " ").split())
        except Exception:
            raise SystemExit(f"--center must be 'x,y,z', got {args.center!r}")
        center_xyz = (cx, cy, cz)
        radius_eff = args.radius
        elog.info(
            "[ph.scope] scope=POCKET center=(%.3f,%.3f,%.3f) r=%.2f",
            cx,
            cy,
            cz,
            radius_eff,
        )
    else:
        # GLOBAL: center optional; if absent, auto to origin (ignored by propka_wire in global)
        if args.center:
            try:
                cx, cy, cz = (float(t) for t in args.center.replace(",", " ").split())
                center_xyz = (cx, cy, cz)
                center_note = "provided"
            except Exception:
                raise SystemExit(f"--center must be 'x,y,z', got {args.center!r}")
        else:
            center_xyz = (0.0, 0.0, 0.0)
            center_note = "auto"
        radius_eff = 1_000_000.0  # triggers GLOBAL in propka_wire (r>=1e6)
        elog.info(
            "[ph.scope] scope=GLOBAL center=%s r=ALL (sentinel=%.0f)",
            center_note,
            radius_eff,
        )

    # A → B → C for each pH; tags end with _0 per acceptance
    build_ph_ensemble(
        pdb_id=pdb_id,
        cleaned_receptor_pdb=cleaned_pdb,
        out_dir=args.output_root,
        center=center_xyz,
        radius=radius_eff,
        ph_values=ph_values,
        member_index_start=0,
    )


if __name__ == "__main__":
    main()
