# cap_driver.py (PyMOL-safe CLI wrapper)
import importlib
import os
import sys
from pathlib import Path
from docking.capture_pose import (
    capture_pose,
)  # expects kwargs: label_top_n_res, label_cutoff, viewport, transparent


def _require_pymol_cmd():
    try:
        return importlib.import_module("pymol").cmd
    except Exception as exc:
        raise RuntimeError(
            "PyMOL is required for capture rendering. Install with: "
            "conda install -c conda-forge pymol-open-source"
        ) from exc


def _posargs(argv):
    """Return only non-option args (drop PyMOL flags like -c, -q, -r, -d, etc.)."""
    out = []
    skip_next = False
    for a in argv:
        if skip_next:
            skip_next = False
            continue
        if a.startswith("-"):
            # Options that consume a parameter (be conservative)
            if a in ("-d", "-r", "-u", "-p", "-l", "-R"):
                skip_next = True
            continue
        out.append(a)
    return out


def _get_bool(name, default):
    v = os.environ.get(name, "")
    if v == "":
        return default
    return str(v).strip().lower() in ("1", "true", "yes", "on")


def main():
    cmd = _require_pymol_cmd()
    # Keep PyMOL-safe parsing: first strip PyMOL flags from sys.argv
    pos = _posargs(sys.argv[1:])
    if len(pos) < 3:
        print(
            f"usage: cap_driver.py <receptor> <ligand> <outdir-or-prefix> "
            f"[--label-top-n-res INT] [--label-cutoff-A FLOAT] [--viewport WxH] [--transparent 0|1]  (argv={sys.argv})"
        )
        cmd.quit(2)

    # Minimal argparse on the remaining positional tokens (no PyMOL flags left)
    import argparse

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
    args, unknown = parser.parse_known_args(
        pos[-len(pos) :]
    )  # parse all remaining tokens

    # Normalize viewport
    try:
        if "x" in args.viewport.lower():
            w, h = args.viewport.lower().split("x", 1)
            viewport = (int(w), int(h))
        else:
            viewport = (int(args.viewport), int(args.viewport))
    except Exception:
        viewport = (1600, 1200)

    # Ensure output dir exists (supports dir or prefix)
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
    except Exception as e:
        print(f"[cap_driver] render failed: {e}")
        cmd.quit(2)


if __name__ == "__main__":
    main()
