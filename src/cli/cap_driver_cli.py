from __future__ import annotations

import argparse
import importlib
import os
import sys
from pathlib import Path

from docking.capture_pose import capture_pose


def _require_pymol_cmd():
    try:
        return importlib.import_module("pymol").cmd
    except Exception as exc:
        raise RuntimeError(
            "PyMOL is required for capture rendering. Install with: "
            "conda install -c conda-forge pymol-open-source"
        ) from exc


def _posargs(argv: list[str]) -> list[str]:
    """Return only non-option args (drop PyMOL flags like -c, -q, -r, -d, etc.)."""
    out: list[str] = []
    skip_next = False
    for token in argv:
        if skip_next:
            skip_next = False
            continue
        if token.startswith("-"):
            # Options that consume a parameter (be conservative).
            if token in ("-d", "-r", "-u", "-p", "-l", "-R"):
                skip_next = True
            continue
        out.append(token)
    return out


def _get_bool(name: str, default: bool) -> bool:
    value = os.environ.get(name, "")
    if value == "":
        return default
    return str(value).strip().lower() in ("1", "true", "yes", "on")


def main() -> int:
    cmd = _require_pymol_cmd()
    pos = _posargs(sys.argv[1:])
    if len(pos) < 3:
        print(
            f"usage: atlas-cap-driver <receptor> <ligand> <outdir-or-prefix> "
            f"[--label-top-n-res INT] [--label-cutoff-A FLOAT] [--viewport WxH] "
            f"[--transparent 0|1]  (argv={sys.argv})"
        )
        cmd.quit(2)
        return 2

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("receptor")
    parser.add_argument("ligand")
    parser.add_argument("out")
    parser.add_argument(
        "--label-top-n-res",
        dest="label_top_n_res",
        type=int,
        default=int(os.environ.get("PYMOL_LABEL_TOP_N_RES", "5")),
    )
    parser.add_argument(
        "--label-cutoff-A",
        dest="label_cutoff",
        type=float,
        default=float(os.environ.get("PYMOL_LABEL_CUTOFF_A", "5.0")),
    )
    parser.add_argument(
        "--viewport",
        dest="viewport",
        type=str,
        default=os.environ.get("PYMOL_VIEWPORT", "1600x1200"),
    )
    parser.add_argument(
        "--transparent",
        dest="transparent",
        type=int,
        default=1 if _get_bool("PYMOL_TRANSPARENT", True) else 0,
    )
    args, _unknown = parser.parse_known_args(pos[-len(pos) :])

    try:
        if "x" in args.viewport.lower():
            width, height = args.viewport.lower().split("x", 1)
            viewport = (int(width), int(height))
        else:
            viewport = (int(args.viewport), int(args.viewport))
    except Exception:
        viewport = (1600, 1200)

    out_path = Path(args.out)
    out_parent = (
        out_path
        if out_path.suffix == "" and (str(out_path).endswith("/") or out_path.is_dir())
        else out_path.parent
    )
    out_parent.mkdir(parents=True, exist_ok=True)

    try:
        capture_pose(
            args.receptor,
            args.ligand,
            str(args.out),
            label_top_n_res=int(args.label_top_n_res),
            label_cutoff=float(args.label_cutoff),
            viewport=viewport,
            transparent=bool(int(args.transparent)),
        )
        cmd.quit(0)
        return 0
    except Exception as exc:
        print(f"[cap_driver] render failed: {exc}")
        cmd.quit(2)
        return 2
