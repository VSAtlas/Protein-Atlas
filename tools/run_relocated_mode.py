#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path

try:
    import fcntl
except Exception:  # pragma: no cover
    fcntl = None  # type: ignore[assignment]


def _log(message: str) -> None:
    print(f"[relocated-mode] {message}")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Wrapper mode for running main.py with relocated PREPPED_LIGANDS_DIR "
            "plus optional EXTRACTED_LIGANDS_DIR/DOCKED_DIR/POST_DOCKED_DIR. "
            "Use --all-dirs to relocate and pre-create a full path bundle."
        )
    )
    parser.add_argument(
        "--prepped-root",
        default=None,
        help="Destination root for prepared ligand libraries.",
    )
    parser.add_argument(
        "--extracted-root",
        default=None,
        help="Destination root for extracted ligand libraries.",
    )
    parser.add_argument(
        "--docked-root",
        default=None,
        help="Optional destination root for docked outputs.",
    )
    parser.add_argument(
        "--post-docked-root",
        default=None,
        help="Optional destination root for post_docked outputs.",
    )
    parser.add_argument(
        "--all-dirs",
        "-all-dirs",
        default=None,
        help=(
            "Optional root under which prepped_ligands, extracted_ligands, docked, "
            "post_docked, processed_pdbs, configs, manifests, logs, and data are "
            "relocated and created."
        ),
    )
    parser.add_argument(
        "--prepped-archive",
        default=None,
        help=(
            "Optional path to prepped ligands tar.zst archive. "
            "Default: <repo-root>/prepped_ligands.tar.zst."
        ),
    )
    parser.add_argument(
        "--force-extract",
        action="store_true",
        help="Force archive extraction even if prepped root already has data.",
    )
    parser.add_argument(
        "--skip-extract",
        action="store_true",
        help="Skip archive extraction logic entirely.",
    )
    parser.add_argument(
        "--main-args",
        required=True,
        help='Arguments forwarded to main.py, e.g. "--run-id reloc --pdb TEST --fast --test-fda".',
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print effective environment/command and exit without running main.py.",
    )
    return parser.parse_args()


def _ensure_tool_available(tool_name: str) -> None:
    if shutil.which(tool_name):
        return
    raise RuntimeError(f"Required tool not found on PATH: {tool_name}")


def _has_prepped_payload(prepped_root: Path) -> bool:
    if not prepped_root.exists():
        return False
    for _, _, filenames in os.walk(prepped_root, followlinks=True):
        for name in filenames:
            if name == "_manifest.json" or name.endswith(".pdbqt"):
                return True
    return False


def _clear_directory_contents(path: Path) -> None:
    for child in path.iterdir():
        if child.is_symlink() or child.is_file():
            child.unlink()
            continue
        shutil.rmtree(child)


def _copy_tree_contents(source_root: Path, dest_root: Path) -> None:
    dest_root.mkdir(parents=True, exist_ok=True)
    for item in source_root.iterdir():
        target = dest_root / item.name
        if item.is_symlink():
            if target.exists() or target.is_symlink():
                if target.is_dir() and not target.is_symlink():
                    shutil.rmtree(target)
                else:
                    target.unlink()
            target.symlink_to(os.readlink(item))
            continue
        if item.is_dir():
            shutil.copytree(item, target, dirs_exist_ok=True, symlinks=True)
            continue
        shutil.copy2(item, target)


def _rewrite_symlink_target(target: str, old_base: Path, new_base: Path) -> str:
    target_path = Path(target)
    if not target_path.is_absolute():
        return target
    try:
        suffix = target_path.relative_to(old_base)
    except ValueError:
        return target
    return str(new_base / suffix)


def _rewrite_tree_symlinks(dest_root: Path, old_base: Path, new_base: Path) -> int:
    rewritten = 0
    for current_root, dirnames, filenames in os.walk(dest_root, followlinks=False):
        root_path = Path(current_root)
        for name in [*dirnames, *filenames]:
            link_path = root_path / name
            if not link_path.is_symlink():
                continue
            current_target = os.readlink(link_path)
            rewritten_target = _rewrite_symlink_target(
                current_target, old_base, new_base
            )
            if rewritten_target == current_target:
                continue
            link_path.unlink()
            link_path.symlink_to(rewritten_target)
            rewritten += 1
    return rewritten


@contextmanager
def _exclusive_extract_lock(lock_path: Path):
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+", encoding="utf-8") as handle:
        if fcntl is not None:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            if fcntl is not None:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _extract_prepped_archive(
    archive_path: Path, prepped_root: Path, force_extract: bool
) -> bool:
    if not archive_path.exists():
        raise FileNotFoundError(f"Archive does not exist: {archive_path}")
    _ensure_tool_available("tar")
    _ensure_tool_available("unzstd")

    prepped_root.mkdir(parents=True, exist_ok=True)
    lock_path = prepped_root.parent / ".prepped_extract.lock"
    with _exclusive_extract_lock(lock_path):
        already_populated = _has_prepped_payload(prepped_root)
        if already_populated and not force_extract:
            _log(f"skip_extract reason=prepped_root_populated path={prepped_root}")
            return False

        if force_extract and any(prepped_root.iterdir()):
            _log(f"force_extract=true clearing_existing path={prepped_root}")
            _clear_directory_contents(prepped_root)

        with tempfile.TemporaryDirectory(
            prefix="prepped_extract_",
            dir=str(prepped_root.parent),
        ) as tmp:
            tmp_root = Path(tmp)
            cmd = [
                "tar",
                "--use-compress-program=unzstd",
                "-xf",
                str(archive_path),
                "-C",
                str(tmp_root),
            ]
            _log(f"extract_cmd={' '.join(shlex.quote(part) for part in cmd)}")
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                stderr = (proc.stderr or "").strip()
                raise RuntimeError(
                    f"Archive extraction failed (exit={proc.returncode}): {stderr}"
                )

            source_root = tmp_root / "prepped_ligands"
            if not source_root.is_dir():
                source_root = tmp_root
            if not any(source_root.iterdir()):
                raise RuntimeError(
                    f"Archive extracted but no files were found under: {source_root}"
                )

            old_base = archive_path.parent.resolve()
            new_base = prepped_root.parent.resolve()
            try:
                _copy_tree_contents(source_root, prepped_root)
                rewritten_links = _rewrite_tree_symlinks(prepped_root, old_base, new_base)
            except (shutil.Error, OSError) as exc:
                raise RuntimeError(
                    "Archive extraction copy failed while materializing prepared ligands "
                    f"from {source_root} to {prepped_root}. "
                    "This can happen when archived symlinks are dangling or have "
                    "machine-specific absolute targets."
                ) from exc
            if rewritten_links:
                _log(
                    "rewrote_symlink_targets "
                    f"count={rewritten_links} old_base={old_base} new_base={new_base}"
                )
            _log(f"extract_done archive={archive_path} destination={prepped_root}")
            return True


def _build_main_command(repo_root: Path, main_args_raw: str) -> list[str]:
    parsed = shlex.split(main_args_raw)
    if not parsed:
        raise ValueError("--main-args resolved to an empty argument list.")
    return [sys.executable, str(repo_root / "main.py"), *parsed]


def _default_prepped_archive_path(repo_root: Path) -> Path:
    return repo_root / "prepped_ligands.tar.zst"


def _resolve_root(
    explicit_root: str | None, all_dirs_root: Path | None, leaf_name: str
) -> Path | None:
    if explicit_root:
        return Path(explicit_root).expanduser().resolve()
    if all_dirs_root is not None:
        if leaf_name in {"docked", "post_docked"}:
            return (all_dirs_root / "outputs" / leaf_name).resolve()
        return (all_dirs_root / leaf_name).resolve()
    return None


def main() -> int:
    args = _parse_args()

    if args.skip_extract and args.force_extract:
        raise ValueError("Use only one of --skip-extract or --force-extract.")
    if not args.prepped_root and not args.all_dirs:
        raise ValueError("Provide --prepped-root or --all-dirs.")

    repo_root = Path(__file__).resolve().parents[1]
    all_dirs_root = (
        Path(args.all_dirs).expanduser().resolve() if args.all_dirs else None
    )
    prepped_root = _resolve_root(args.prepped_root, all_dirs_root, "prepped_ligands")
    if prepped_root is None:
        raise RuntimeError("Unable to resolve prepped ligands root.")
    extracted_root = _resolve_root(
        args.extracted_root, all_dirs_root, "extracted_ligands"
    )
    docked_root = _resolve_root(args.docked_root, all_dirs_root, "docked")
    post_docked_root = _resolve_root(
        args.post_docked_root, all_dirs_root, "post_docked"
    )
    output_root = (
        (all_dirs_root / "outputs" / "processed_pdbs").resolve()
        if all_dirs_root
        else None
    )
    configs_root = (
        (all_dirs_root / "outputs" / "configs").resolve() if all_dirs_root else None
    )
    manifests_root = (
        (all_dirs_root / "outputs" / "manifests").resolve()
        if all_dirs_root
        else None
    )
    logs_root = (
        (all_dirs_root / "outputs" / "logs").resolve() if all_dirs_root else None
    )
    data_root = (
        (all_dirs_root / "outputs" / "data").resolve() if all_dirs_root else None
    )

    roots_to_create = [
        all_dirs_root,
        prepped_root,
        extracted_root,
        docked_root,
        post_docked_root,
        output_root,
        configs_root,
        manifests_root,
        logs_root,
        data_root,
    ]
    for root in roots_to_create:
        if root is not None:
            try:
                root.mkdir(parents=True, exist_ok=True)
            except PermissionError as exc:
                raise RuntimeError(
                    "Unable to create directory "
                    f"{root}: permission denied. "
                    "Choose a writable root or pre-create the path."
                ) from exc
            except OSError as exc:
                raise RuntimeError(
                    "Unable to create directory "
                    f"{root}: {exc}. "
                    "Choose a writable root or pre-create the path."
                ) from exc

    archive_path = (
        Path(args.prepped_archive).expanduser().resolve()
        if args.prepped_archive
        else _default_prepped_archive_path(repo_root).resolve()
    )

    cmd = _build_main_command(repo_root, args.main_args)
    env = os.environ.copy()
    env["PREPPED_LIGANDS_DIR"] = str(prepped_root)
    if extracted_root is not None:
        env["EXTRACTED_LIGANDS_DIR"] = str(extracted_root)
    if docked_root is not None:
        env["DOCKED_DIR"] = str(docked_root)
    if post_docked_root is not None:
        env["POST_DOCKED_DIR"] = str(post_docked_root)
    if output_root is not None:
        env["OUTPUT_DIR"] = str(output_root)
    if configs_root is not None:
        env["CONFIGS_DIR"] = str(configs_root)
    if manifests_root is not None:
        env["MANIFESTS_DIR"] = str(manifests_root)
    if logs_root is not None:
        env["LOGS_DIR"] = str(logs_root)
    if data_root is not None:
        env["DATA_DIR"] = str(data_root)

    _log(f"effective PREPPED_LIGANDS_DIR={prepped_root}")
    if extracted_root is not None:
        _log(f"effective EXTRACTED_LIGANDS_DIR={extracted_root}")
    else:
        _log("effective EXTRACTED_LIGANDS_DIR=from main.py config resolution")
    if docked_root is not None:
        _log(f"effective DOCKED_DIR={docked_root}")
    else:
        _log("effective DOCKED_DIR=from main.py config resolution")
    if post_docked_root is not None:
        _log(f"effective POST_DOCKED_DIR={post_docked_root}")
    else:
        _log("effective POST_DOCKED_DIR=from main.py/config default resolution")
    if output_root is not None:
        _log(f"effective OUTPUT_DIR={output_root}")
    else:
        _log("effective OUTPUT_DIR=from main.py config resolution")
    if configs_root is not None:
        _log(f"effective CONFIGS_DIR={configs_root}")
    else:
        _log("effective CONFIGS_DIR=from main.py config resolution")
    if manifests_root is not None:
        _log(f"effective MANIFESTS_DIR={manifests_root}")
    else:
        _log("effective MANIFESTS_DIR=from main.py config resolution")
    if logs_root is not None:
        _log(f"effective LOGS_DIR={logs_root}")
    else:
        _log("effective LOGS_DIR=from main.py config resolution")
    if data_root is not None:
        _log(f"effective DATA_DIR={data_root}")
    else:
        _log("effective DATA_DIR=from main.py config resolution")
    _log(f"effective PREPPED_ARCHIVE={archive_path}")
    _log(f"main_cmd={' '.join(shlex.quote(part) for part in cmd)}")

    if args.dry_run:
        _log("dry_run=true; main.py not executed.")
        return 0

    if args.skip_extract:
        _log("skip_extract reason=flag_set")
    else:
        if not archive_path.exists():
            raise RuntimeError(
                "Prepped ligands archive does not exist: "
                f"{archive_path}. "
                "Pass --prepped-archive <path> or use --skip-extract."
            )
        _extract_prepped_archive(archive_path, prepped_root, args.force_extract)

    if not _has_prepped_payload(prepped_root):
        raise RuntimeError(
            "No prepared ligand payload found in prepped root. "
            f"Expected .pdbqt or _manifest.json under {prepped_root}."
        )

    completed = subprocess.run(cmd, cwd=str(repo_root), env=env)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
