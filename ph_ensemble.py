 # ---pH-ensemble A→B→C orchestration ---

import os
import json
import argparse
import logging
import shutil
from pathlib import Path
from typing import Iterable, Tuple

from propka_wire import apply_propka_states
import automate_protein_prep

elog = logging.getLogger("ph_ensemble")
if not elog.handlers:
    _h = logging.StreamHandler()
    _h.setFormatter(logging.Formatter("[ph_ensemble] %(message)s"))
    elog.addHandler(_h)
elog.setLevel(logging.INFO)







import hashlib

def _sha1_of_file(p: str) -> str:
    h = hashlib.sha1()
    try:
        with open(p, "rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()
    except Exception:
        return ""

def _count_atoms_pdb(path: str) -> int:
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            return sum(1 for ln in fh if ln.startswith(("ATOM  ", "HETATM")))
    except Exception:
        return 0
def _count_pocket_hydrogens(pdb_path: str, center: Tuple[float, float, float], radius: float) -> int:
    cx, cy, cz = center
    r2 = radius * radius
    nH = 0
    try:
        with open(pdb_path, "r", errors="ignore") as fh:
            for ln in fh:
                if not ln.startswith(("ATOM  ", "HETATM")):
                    continue
                try:
                    x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
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

def _pocket_name_counts(pdb_path: str, center: Tuple[float,float,float], radius: float) -> dict:
    cx, cy, cz = center
    r2 = radius*radius
    counts = {"ASP":0,"ASH":0,"GLU":0,"GLH":0,"HIS":0,"HID":0,"HIE":0,"HIP":0,"LYS":0,"LYN":0,"CYS":0,"CYM":0,"TYR":0}
    seen = set()
    try:
        with open(pdb_path, "r", errors="ignore") as fh:
            for ln in fh:
                if not ln.startswith(("ATOM  ","HETATM")): continue
                try:
                    x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
                except Exception:
                    continue
                dx,dy,dz = x-cx, y-cy, z-cz
                if (dx*dx+dy*dy+dz*dz) > r2: continue
                chain = (ln[21].strip() or " ")
                resi = int(ln[22:26])
                key = (chain, resi)
                if key in seen: continue
                seen.add(key)
                resn = ln[17:20].strip().upper()
                if resn in counts: counts[resn] += 1
    except Exception:
        pass
    return {k:v for k,v in counts.items() if v>0}



def build_ph_ensemble(
    pdb_id: str,
    cleaned_receptor_pdb: str,
    out_dir: str,
    center: Tuple[float,float,float],
    radius: float,
    ph_values: Iterable[float],
    member_index_start: int = 0
) -> str:
    """
    For each pH value: A) PROPKA-local renames -> prestate
                       B) single hydrogenation via automate_protein_prep.assign_protonation_states
                       C) app.run_prepare_receptor -> pdbqt
    Writes ensemble.json and returns its path.
    """
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    ensemble_dir = out_dir / pdb_id / "receptor" / "ph_ensemble"
    ensemble_dir.mkdir(parents=True, exist_ok=True)

    ph_values = list(ph_values)
    elog.info("[ph.list] n=%d values=%s", len(ph_values), ph_values)

    members = []
    mi = member_index_start
    _prev_pre_sha = None
    _prev_with_sha = None
    _prev_pdbqt_size = None
    for ph in ph_values:
        tag = f"{pdb_id}_pH{str(ph).replace('.', '_')}"
        elog.info("[ph.member.start] ph=%.2f tag=%s", ph, tag)
        elog.info("[ph.propka.call] ph=%.2f center=(%.3f,%.3f,%.3f) r=%.2f", ph, center[0], center[1], center[2],
                  radius)

        # A) prestate (no H-add here)
        prestate_pdb, pka_path, _n = apply_propka_states(
            cleaned_receptor_pdb=cleaned_receptor_pdb,
            center=center,
            radius=radius,
            ph=ph,
            out_dir=str(ensemble_dir),
            tag=tag
        )
        # Fingerprint prestate
        pre_sha = _sha1_of_file(prestate_pdb);
        pre_atoms = _count_atoms_pdb(prestate_pdb)
        pre_bins = _pocket_name_counts(prestate_pdb, center, radius)
        elog.info("[stage.A.prestate] ph=%.2f atoms=%d sha1=%s bins=%s",
                  ph, pre_atoms, pre_sha, " ".join(f"{k}={v}" for k, v in pre_bins.items()))
        if _prev_pre_sha is not None:
            elog.info("[compare.prestate] ph=%.2f same_as_prev=%s", ph, str(pre_sha == _prev_pre_sha))
        _prev_pre_sha = pre_sha

        # B) single protonation  (assign_protonation_states)
        withH_pdb = str(ensemble_dir / f"{tag}.withH.pdb")
        automate_protein_prep.assign_protonation_states(
            input_pdb=prestate_pdb,
            output_pdb=withH_pdb,
            reduce_exe=automate_protein_prep.REDUCE_EXE
        )

        with_sha = _sha1_of_file(withH_pdb)
        with_atoms = _count_atoms_pdb(withH_pdb)
        with_bins = _pocket_name_counts(withH_pdb, center, radius)

        pocket_H = _count_pocket_hydrogens(withH_pdb, center, radius)
        elog.info("[stage.B.withH] ph=%.2f atoms=%d sha1=%s pocket_H=%d bins=%s",
                  ph, with_atoms, with_sha, pocket_H,
                  " ".join(f"{k}={v}" for k, v in with_bins.items()))

        if _prev_with_sha is not None:
            elog.info("[compare.withH] ph=%.2f same_as_prev=%s", ph, str(with_sha == _prev_with_sha))
        _prev_with_sha = with_sha


        # C) receptor prep (pdbqt)
        pdbqt_path = str(ensemble_dir / f"{tag}.pdbqt")
        ok = automate_protein_prep.run_prepare_receptor(
            input_pdb=withH_pdb,
            output_pdbqt=pdbqt_path,
            cfg=automate_protein_prep.config
        )
        if not ok:
            shutil.copy2(withH_pdb, pdbqt_path)

        try:
            pdbqt_size = Path(pdbqt_path).stat().st_size
        except Exception:
            pdbqt_size = -1
        elog.info("[stage.C.pdbqt] ph=%.2f size=%d path=%s", ph, pdbqt_size, pdbqt_path)
        if _prev_pdbqt_size is not None:
            elog.info("[compare.pdbqt] ph=%.2f same_size_as_prev=%s", ph, str(pdbqt_size == _prev_pdbqt_size))
        _prev_pdbqt_size = pdbqt_size
        elog.info("[ensemble.step] ph=%.2f prestate=%s withH=%s pdbqt=%s",
                  ph, prestate_pdb, withH_pdb, pdbqt_path)

        members.append({
            "ph": ph,
            "tag": tag,
            "prestate": prestate_pdb,
            "pka": pka_path,
            "withH": withH_pdb,
            "pdbqt": pdbqt_path
        })
        mi += 1

    manifest = {
        "pdb_id": pdb_id,
        "ensemble": ph_values,
        "members": members
    }
    manifest_path = ensemble_dir / "ensemble.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    return str(manifest_path)

# --- END: pH-ensemble A→B→C orchestration ---







def main():
    ap = argparse.ArgumentParser(description="Build pH-aware receptor ensemble.")
    ap.add_argument("pdb_path", help="Input receptor PDB (ligand-free).")
    ap.add_argument("--output-root", default="processed_pdbs")
    ap.add_argument("--force", action="store_true", help="Force re-prep of base receptor.")
    ap.add_argument("--center", type=str, default=None,
                    help="Binding-site center as 'x,y,z' (required if --scope=pocket).")
    ap.add_argument("--radius", type=float, default=10.0,
                    help="Sphere radius in Å for local titration (used when --scope=pocket).")
    ap.add_argument("--scope", choices=["global", "pocket"], default="global",
                    help="Renaming scope: 'global' (entire protein; default) or 'pocket' (within --radius).")
    ap.add_argument("--ph-list", type=str, default="7.4",
                    help="Comma-separated pH values, e.g., 5,7,9")
    args = ap.parse_args()

    # Parse pH list (always required to be valid)
    try:
        ph_values = [float(x) for x in args.ph_list.split(",") if x.strip() != ""]
    except Exception:
        raise SystemExit(f"--ph-list must be comma-separated numbers, got {args.ph_list!r}")

    if args.force:
        os.environ["FORCE_REPROCESS"] = "1"

    # Single clean step (do not edit automate_protein_prep)
    cleaned_pdb, _receptor_pdbqt = automate_protein_prep.main(args.pdb_path, args.output_root)
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
        elog.info("[ph.scope] scope=POCKET center=(%.3f,%.3f,%.3f) r=%.2f", cx, cy, cz, radius_eff)
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
        elog.info("[ph.scope] scope=GLOBAL center=%s r=ALL (sentinel=%.0f)", center_note, radius_eff)

    # A → B → C for each pH; tags end with _0 per acceptance
    build_ph_ensemble(
        pdb_id=pdb_id,
        cleaned_receptor_pdb=cleaned_pdb,
        out_dir=args.output_root,
        center=center_xyz,
        radius=radius_eff,
        ph_values=ph_values,
        member_index_start=0
    )


if __name__ == "__main__":
    main()
