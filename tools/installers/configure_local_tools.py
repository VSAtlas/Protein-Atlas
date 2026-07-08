#!/usr/bin/env python3
"""Detect local external tools and update config.txt without touching examples."""

from __future__ import print_function

import argparse
import os
import shutil
from pathlib import Path
from typing import Sequence


TOOL_KEYS = (
    "VINA_EXE",
    "OPENBABEL_PATH",
    "MGLTOOLS_PATH",
    "ADFRSUITE_BIN",
    "MGLTOOLS_PYTHON",
    "PREPARE_LIGAND_SCRIPT",
    "PREPARE_RECEPTOR_SCRIPT",
    "P2RANK_PATH",
    "SCORCH",
    "SCORCH_ENV",
    "SCORCH_ENV_PREFIX",
    "MMGBSA_AMBERTOOLS_PREFIX",
    "AMBERTOOLS_PREFIX",
)


def repo_root():
    return Path(__file__).resolve().parents[2]


def atlas_tools_dir():
    root = repo_root()
    default_tools = root.parent.parent / "tools"
    return Path(os.environ.get("ATLAS_TOOLS_DIR") or default_tools).expanduser()


def parse_set(values):
    parsed = {}
    for item in values:
        if "=" not in item:
            raise SystemExit("--set expects KEY=VALUE, got {!r}".format(item))
        key, value = item.split("=", 1)
        key = key.strip()
        if not key:
            raise SystemExit("--set has an empty key: {!r}".format(item))
        parsed[key] = value.strip()
    return parsed


def first_existing(paths):
    for path in paths:
        try:
            if path.exists():
                return path
        except OSError:
            continue
    return None


def unique_paths(paths):
    seen = set()
    out = []
    for raw in paths:
        if not raw:
            continue
        path = Path(str(raw)).expanduser()
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        out.append(path)
    return out


def env_prefixes():
    tools_dir = atlas_tools_dir()
    prefixes = [
        os.environ.get("CONDA_PREFIX"),
        os.environ.get("MAMBA_PREFIX"),
        os.environ.get("VIRTUAL_ENV"),
    ]
    mamba_root = os.environ.get("MAMBA_ROOT_PREFIX")
    if mamba_root:
        root = Path(mamba_root).expanduser()
        prefixes.extend(
            [
                root / "envs" / "docking-env",
                root / "envs" / "scorch-env",
                root / "envs" / "AmberTools25",
                root / "envs" / "ambertools",
            ]
        )
    home = Path.home()
    prefixes.extend(
        [
            tools_dir / "envs" / "docking-env",
            tools_dir / "envs" / "scorch-env",
            tools_dir / "envs" / "AmberTools25",
            tools_dir / "envs" / "ambertools",
            tools_dir / "envs" / "banana",
            home / "micromamba" / "envs" / "docking-env",
            home / "micromamba" / "envs" / "scorch-env",
            home / "micromamba" / "envs" / "AmberTools25",
            home / "miniconda3" / "envs" / "docking-env",
            home / "miniconda3" / "envs" / "scorch-env",
            home / "miniconda3" / "envs" / "AmberTools25",
        ]
    )
    return unique_paths(prefixes)


def first_tool_in_prefixes(command, prefixes):
    rel_candidates = [
        Path("bin") / command,
        Path("Scripts") / command,
        Path(command),
    ]
    for prefix in prefixes:
        for rel in rel_candidates:
            candidate = prefix / rel
            if candidate.is_file():
                return candidate
    return None


def detect_adfrsuite_bin(prefixes):
    command_hits = [
        shutil.which("prepare_ligand"),
        shutil.which("prepare_receptor"),
        shutil.which("prepare_ligand4.py"),
        shutil.which("prepare_receptor4.py"),
    ]
    for hit in command_hits:
        if hit:
            return str(Path(hit).resolve().parent)

    candidates = []
    for prefix in prefixes:
        candidates.extend(
            [
                prefix / "bin",
                prefix / "ADFRsuite" / "bin",
                prefix / "ADFRSuite" / "bin",
            ]
        )
    found = first_existing(
        path
        for path in candidates
        if (path / "prepare_ligand").is_file()
        or (path / "prepare_receptor").is_file()
        or (path / "prepare_ligand4.py").is_file()
    )
    return str(found) if found else None


def detect_mgltools_root(pythonsh):
    candidates = []
    if pythonsh:
        py_path = Path(pythonsh).resolve()
        candidates.extend(py_path.parents[: min(len(py_path.parents), 6)])

    for prefix in env_prefixes():
        candidates.extend(
            [
                prefix,
                prefix / "MGLToolsPckgs",
                prefix / "mgltools",
                prefix / "MGLTools-1.5.7",
            ]
        )

    root = first_existing(
        path
        for path in candidates
        if (path / "MGLToolsPckgs").exists()
        or (path / "AutoDockTools").exists()
        or (path / "bin" / "pythonsh").is_file()
    )
    return str(root) if root else None


def detect_p2rank_path(prefixes):
    prank = shutil.which("prank")
    if prank:
        prank_path = Path(prank).resolve()
        if prank_path.parent.name == "bin":
            return str(prank_path.parent.parent)
        return str(prank_path.parent)

    for prefix in prefixes:
        candidates = [
            prefix / "bin" / "prank",
            prefix / "p2rank" / "prank",
            prefix / "p2rank" / "bin" / "prank",
            prefix / "p2rank.jar",
            prefix / "p2rank" / "p2rank.jar",
        ]
        found = first_existing(candidates)
        if found:
            return str(found.parent if found.is_file() and found.name == "prank" else found)
    return None


def detect_env_prefix(name, required_tools):
    for prefix in env_prefixes():
        if prefix.name.lower() != name.lower():
            continue
        if all((prefix / "bin" / tool).is_file() for tool in required_tools):
            return str(prefix)
    for prefix in env_prefixes():
        if all((prefix / "bin" / tool).is_file() for tool in required_tools):
            return str(prefix)
    return None


def detect_tools():
    detected = {}
    prefixes = env_prefixes()
    command_map = {
        "VINA_EXE": "vina",
        "OPENBABEL_PATH": "obabel",
        "MGLTOOLS_PYTHON": "pythonsh",
        "PREPARE_LIGAND_SCRIPT": "prepare_ligand4.py",
        "PREPARE_RECEPTOR_SCRIPT": "prepare_receptor4.py",
        "SCORCH": "scorch.py",
    }
    for key, command in command_map.items():
        value = shutil.which(command)
        if not value:
            candidate = first_tool_in_prefixes(command, prefixes)
            value = str(candidate) if candidate else ""
        if value:
            detected[key] = value

    p2rank_path = detect_p2rank_path(prefixes)
    if p2rank_path:
        detected["P2RANK_PATH"] = p2rank_path

    adfr_bin = detect_adfrsuite_bin(prefixes)
    if adfr_bin:
        detected["ADFRSUITE_BIN"] = adfr_bin

    mgl_root = detect_mgltools_root(detected.get("MGLTOOLS_PYTHON"))
    if mgl_root:
        detected["MGLTOOLS_PATH"] = mgl_root

    scorch_prefix = detect_env_prefix("scorch-env", ("python",))
    if scorch_prefix:
        detected["SCORCH_ENV_PREFIX"] = scorch_prefix

    amber_prefix = detect_env_prefix("AmberTools25", ("tleap", "cpptraj"))
    if amber_prefix:
        detected["AMBERTOOLS_PREFIX"] = amber_prefix
        detected["MMGBSA_AMBERTOOLS_PREFIX"] = amber_prefix

    detected.setdefault("SCORCH_ENV", "scorch-env")
    return detected


def read_existing(path):
    if not path.exists():
        return [], {}

    lines = path.read_text(encoding="utf-8").splitlines()
    values = {}
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or "=" not in stripped:
            continue
        key, value = stripped.split("=", 1)
        values[key.strip()] = value.strip()
    return lines, values


def update_lines(lines, updates):
    remaining = dict(updates)
    output = []
    for line in lines:
        stripped = line.strip()
        if stripped and not stripped.startswith("#") and "=" in stripped:
            key, _value = stripped.split("=", 1)
            key = key.strip()
            if key in remaining:
                output.append("{}={}".format(key, remaining.pop(key)))
                continue
        output.append(line)

    if remaining:
        if output and output[-1].strip():
            output.append("")
        output.append("# Local external tool paths")
        for key in TOOL_KEYS:
            if key in remaining:
                output.append("{}={}".format(key, remaining.pop(key)))
        for key in sorted(remaining):
            output.append("{}={}".format(key, remaining[key]))
    return output


def main(argv: Sequence[str] | None = None):
    root = repo_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(root / "config.txt"))
    parser.add_argument("--example", default=str(root / "config.example.txt"))
    parser.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--no-detect",
        action="store_true",
        help="Only apply explicit --set values.",
    )
    args = parser.parse_args(argv)

    config_path = Path(args.config).expanduser()
    example_path = Path(args.example).expanduser()

    if not config_path.exists():
        if not example_path.exists():
            raise SystemExit(
                "Missing config and example: {}, {}".format(config_path, example_path)
            )
        lines = example_path.read_text(encoding="utf-8").splitlines()
        existing = {}
    else:
        lines, existing = read_existing(config_path)

    updates = {} if args.no_detect else detect_tools()
    updates.update(parse_set(args.set))
    updates = {
        key: value
        for key, value in updates.items()
        if value and existing.get(key) != value
    }

    if not updates:
        print("No config updates needed for {}".format(config_path))
        return 0

    new_lines = update_lines(lines, updates)
    print("Updating {}:".format(config_path))
    for key in sorted(updates):
        print("  {}={}".format(key, updates[key]))

    if not args.dry_run:
        config_path.write_text("\n".join(new_lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
