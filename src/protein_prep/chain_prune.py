"""Chain selection and pruning policy helpers."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

from protein_prep.aliases_policy import (
    _RETAIN_VARIANT,
    _RETAIN_VARIANT_CANONICAL,
    _normalize_resname,
)
from protein_prep.element_guard import _post_write_element_guard


def _hydrate_legacy_globals() -> None:
    import automate_protein_prep as legacy

    g = globals()
    for name, value in legacy.__dict__.items():
        g.setdefault(name, value)


def _has_backbone_atoms(lines):
    req = {"N", "CA", "C", "O"}
    seen = set()
    for ln in lines:
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        name = ln[12:16].strip()
        if name in req:
            seen.add(name)
    return req.issubset(seen)


def _group_by_chain(lines):
    chains = {}
    for ln in lines:
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        ch = ln[21]
        chains.setdefault(ch, []).append(ln)
    return chains


def _prune_chains_conservative(lines, keep_chains):
    out = []
    for ln in lines:
        if ln.startswith(("ATOM", "HETATM")) and ln[21] not in keep_chains:
            continue
        out.append(ln)
    return out


def _chains_to_keep(lines, pocket_center=None, r=12.0):
    chains = _group_by_chain(lines)
    keep = set()
    # rule 1: chain has full backbone atoms somewhere
    for ch, seg in chains.items():
        if _has_backbone_atoms(seg):
            keep.add(ch)
    # rule 2: optional geometric proximity if center known
    if pocket_center:
        near = set()
        x0, y0, z0 = pocket_center
        for ch, seg in chains.items():
            for ln in seg:
                if not ln.startswith(("ATOM", "HETATM")):
                    continue
                try:
                    x = float(ln[30:38])
                    y = float(ln[38:46])
                    z = float(ln[46:54])
                except Exception:
                    continue
                if (x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2 <= r * r:
                    near.add(ch)
                    break
        if near:
            keep = keep & near
    return keep or set(chains.keys())


def detect_pocket_center_from_ligands(
    filtered_pdb: Union[str, Path], ligands_dir: Union[str, Path]
) -> Optional[Tuple[float, float, float]]:
    """
    Prefer center of extracted ligands in ligands_raw/; fallback to YAML-retained cofactors/metals
    present in filtered_pdb. Returns (x,y,z) or None.
    """
    _hydrate_legacy_globals()
    from statistics import fmean

    pts: list[Tuple[float, float, float]] = []

    # 1) Extracted ligands (*.pdb) under ligands_dir
    ligd = Path(ligands_dir)
    for p in ligd.glob("*.pdb"):
        try:
            with open(p, "r", encoding="utf-8", errors="ignore") as f:
                for ln in f:
                    if ln.startswith(("ATOM", "HETATM")):
                        el = ln[76:78].strip().upper()
                        if el == "H":
                            continue
                        try:
                            x = float(ln[30:38])
                            y = float(ln[38:46])
                            z = float(ln[46:54])
                            pts.append((x, y, z))
                        except Exception:
                            pass
        except Exception:
            pass

    # 2) Fallback: retained cofactors/metals from YAML in filtered_pdb
    if not pts:
        with open(filtered_pdb, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if not ln.startswith("HETATM"):
                    continue
                resname = ln[17:20].strip().upper()
                canonical_res = _normalize_resname(resname)
                if (resname in _RETAIN_VARIANT) or (
                    canonical_res in _RETAIN_VARIANT_CANONICAL
                ):
                    try:
                        x = float(ln[30:38])
                        y = float(ln[38:46])
                        z = float(ln[46:54])
                        pts.append((x, y, z))
                    except Exception:
                        pass

    if not pts:
        return None
    xs, ys, zs = zip(*pts)
    return (fmean(xs), fmean(ys), fmean(zs))


def score_chain_contacts(
    pdb_path: Union[str, Path], keep_chains: set[str]
) -> Dict[str, int]:
    """
    Return heavy-atom pair counts within CHAIN_CONTACT_DIST_ANG between each non-kept chain
    and the *union* of kept chains.
    """
    _hydrate_legacy_globals()
    dist = _cfg_float("CHAIN_CONTACT_DIST_ANG", 5.0)
    atoms_by_chain: dict[str, list[Tuple[float, float, float]]] = {}
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if not ln.startswith(("ATOM  ", "HETATM")):
                continue
            el = ln[76:78].strip().upper()
            if el == "H":
                continue
            c = ln[21]
            try:
                x = float(ln[30:38])
                y = float(ln[38:46])
                z = float(ln[46:54])
            except Exception:
                continue
            atoms_by_chain.setdefault(c, []).append((x, y, z))

    kept_pts = []
    for kc in keep_chains:
        kept_pts.extend(atoms_by_chain.get(kc, []))

    out: Dict[str, int] = {}
    if not kept_pts:
        return {c: 0 for c in atoms_by_chain}  # nothing to compare to

    d2 = dist * dist
    for c, pts in atoms_by_chain.items():
        if c in keep_chains:
            continue
        cnt = 0
        # simple O(N*M); early stop not necessary but harmless
        for x, y, z in pts:
            for X, Y, Z in kept_pts:
                dx = x - X
                dy = y - Y
                dz = z - Z
                if (dx * dx + dy * dy + dz * dz) <= d2:
                    cnt += 1
        out[c] = cnt
    return out


def select_chains_to_keep(
    filtered_pdb: Union[str, Path], ligands_dir: Union[str, Path], cfg=None
) -> set[str]:
    """
    Decide chains to keep using pocket proximity, contact counts, keep-list, and safety rails.
    Emits a compact decision table when CHAIN_LOG_DECISIONS is enabled.
    """
    _hydrate_legacy_globals()
    radius = _cfg_float("CHAIN_POCKET_RADIUS_ANG", 10.0)
    min_contacts = _cfg_int("CHAIN_MIN_CONTACTS", 200)
    keep_list = _cfg_chain_keep_list()
    log_dec = _cfg_bool("CHAIN_LOG_DECISIONS", True)

    pocket = detect_pocket_center_from_ligands(filtered_pdb, ligands_dir)

    # Build per-chain stats
    chains: dict[str, dict] = {}
    all_chains: set[str] = set()
    ca_counts: dict[str, int] = {}
    min_dists: dict[str, float] = {}

    def _dist2(pt, xyz):
        dx = pt[0] - xyz[0]
        dy = pt[1] - xyz[1]
        dz = pt[2] - xyz[2]
        return dx * dx + dy * dy + dz * dz

    with open(filtered_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM  ", "HETATM")):
                c = ln[21]
                all_chains.add(c)
                if ln.startswith("ATOM  ") and ln[12:16].strip() == "CA":
                    ca_counts[c] = ca_counts.get(c, 0) + 1
                if pocket is not None:
                    try:
                        xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
                        d2 = _dist2(pocket, xyz)
                        md = min_dists.get(c, float("inf"))
                        if d2 < md:
                            min_dists[c] = d2
                    except Exception:
                        pass

    if not all_chains:
        return set()  # nothing to decide

    # Initial kept set
    kept: set[str] = set(keep_list)

    # Always keep blank-chain records conservatively
    if " " in all_chains:
        kept.add(" ")

    # Pocket proximity
    if pocket is not None:
        r2 = radius * radius
        for c in all_chains:
            if min_dists.get(c, float("inf")) <= r2:
                kept.add(c)

    # If we still have no pocket & no keep_list -> abort prune
    if (pocket is None) and not keep_list:
        if log_dec:
            logging.info("[chain_prune] skipped reason=no_pocket_center or keep_list")
        return set()  # signal: do not prune

    # Contact expansion to capture interfaces
    contacts = score_chain_contacts(filtered_pdb, kept)
    for c, cnt in contacts.items():
        if cnt >= min_contacts:
            kept.add(c)

    # Safety rails
    if not kept:
        if log_dec:
            logging.info("[chain_prune] conservative_abort reason=kept_empty")
        return set()
    largest_ca_chain = (
        max(ca_counts, key=lambda c: ca_counts.get(c, 0)) if ca_counts else None
    )
    drop_set = {c for c in all_chains if c not in kept}
    if len(drop_set) > len(all_chains) / 2:
        if log_dec:
            logging.warning(
                "[chain_prune] conservative_abort reason=drop_gt_50pct all=%s drop=%s",
                sorted(all_chains),
                sorted(drop_set),
            )
        return set()
    if largest_ca_chain and (largest_ca_chain in drop_set):
        if log_dec:
            logging.warning(
                "[chain_prune] conservative_abort reason=largest_CA_chain_would_be_dropped largest=%s",
                largest_ca_chain,
            )
        return set()

    # Decision table
    if log_dec:
        rows = []
        for c in sorted(all_chains):
            rows.append(
                {
                    "chain": c,
                    "CA_count": ca_counts.get(c, 0),
                    "near_pocket": (
                        "yes"
                        if (
                            pocket is not None
                            and min_dists.get(c, float("inf")) <= (radius * radius)
                        )
                        else "no"
                    ),
                    "min_dist": (
                        0.0
                        if pocket is None
                        else (min_dists.get(c, float("inf")) ** 0.5)
                    ),
                    "contact_count_to_kept": contacts.get(c, 0),
                    "decision": ("keep" if c in kept else "drop"),
                    "reason": (
                        "whitelist"
                        if c in keep_list
                        else "near_pocket"
                        if (
                            pocket is not None
                            and min_dists.get(c, float("inf")) <= (radius * radius)
                        )
                        else "interface_contacts"
                        if contacts.get(c, 0) >= min_contacts
                        else "far_and_sparse"
                    ),
                }
            )
        # Emit compact table
        hdr = "# chain  CA  near  min_d  contacts  keep  reason"
        logging.info(hdr)
        for r in rows:
            logging.info(
                "  %-5s  %-3d %-5s %6.2f    %-7d %-4s  %s",
                r["chain"],
                r["CA_count"],
                r["near_pocket"],
                r["min_dist"],
                r["contact_count_to_kept"],
                ("yes" if r["decision"] == "keep" else "no"),
                r["reason"],
            )
        kept_ids = "".join(sorted(kept)) or "-"
        dropped_ids = "".join(sorted(drop_set)) or "-"
        logging.info(
            "[chain_prune] kept=%s dropped=%s reason=see_table", kept_ids, dropped_ids
        )

    return kept


def prune_to_chains(
    input_pdb: Union[str, Path], kept_chains: set[str], output_pdb: Union[str, Path]
) -> None:
    """
    Write only coordinate records belonging to kept_chains; copy all non-coordinate lines through.
    """
    _hydrate_legacy_globals()
    kept = set(kept_chains or set())
    with (
        open(input_pdb, "r", encoding="utf-8", errors="ignore") as f,
        open(output_pdb, "w", encoding="utf-8") as w,
    ):
        for ln in f:
            if ln.startswith(("ATOM  ", "HETATM")):
                c = ln[21]
                if c in kept:
                    w.write(ln)
            else:
                w.write(ln)
    _post_write_element_guard("chain_pruned", output_pdb)
