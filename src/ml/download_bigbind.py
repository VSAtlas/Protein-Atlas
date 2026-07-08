#!/usr/bin/env python3
from __future__ import annotations

import argparse
import tarfile
import urllib.request
from pathlib import Path, PurePosixPath


def _safe_extract(tar: tarfile.TarFile, destination: Path) -> None:
    for member in tar:
        member_path = PurePosixPath(member.name)
        if member_path.is_absolute() or ".." in member_path.parts:
            raise RuntimeError(f"Unsafe path in tar archive: {member.name}")
        tar.extract(member, destination)


def download_bigbind_tar(*, url: str, output_tar: Path, force: bool) -> Path:
    output_tar.parent.mkdir(parents=True, exist_ok=True)
    if output_tar.exists() and not force:
        print(f"[ml] Reusing existing tar: {output_tar}")
        return output_tar
    print(f"[ml] Downloading BigBind archive from: {url}")
    with urllib.request.urlopen(url) as response, output_tar.open("wb") as handle:
        handle.write(response.read())
    print(f"[ml] Saved archive to: {output_tar}")
    return output_tar


def _looks_like_bigbind_root(path: Path) -> bool:
    return any((path / name).exists() for name in ("activities_train.csv", "activities_val.csv", "activities_test.csv"))


def _detect_bigbind_root(destination_dir: Path) -> Path | None:
    candidates = [destination_dir / "BigBindV1.5", destination_dir / "BigBindV1", destination_dir]
    for candidate in candidates:
        if candidate.exists() and _looks_like_bigbind_root(candidate):
            return candidate
    for child in destination_dir.iterdir() if destination_dir.exists() else []:
        if child.is_dir() and _looks_like_bigbind_root(child):
            return child
    return None


def extract_bigbind_tar(*, tar_path: Path, destination_dir: Path, force: bool) -> Path:
    existing = _detect_bigbind_root(destination_dir)
    if existing is not None and not force:
        print(f"[ml] Reusing extracted directory: {existing}")
        return existing

    destination_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:*") as tar:
        _safe_extract(tar, destination_dir)
    print(f"[ml] Extracted archive into: {destination_dir}")
    detected = _detect_bigbind_root(destination_dir)
    if detected is None:
        raise RuntimeError(f"Could not locate BigBind activity files under {destination_dir}")
    return detected


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download and extract a BigBind data archive."
    )
    parser.add_argument(
        "--url",
        required=True,
        help="URL for a BigBind tar/tar.bz2 archive from the official BigBind release page.",
    )
    parser.add_argument(
        "--out-dir",
        default="data/ml/bigbind",
        help="Directory where BigBindV1.5.tar and extracted data should be placed.",
    )
    parser.add_argument(
        "--filename",
        default="BigBindV1.5.tar",
        help="Local tar filename.",
    )
    parser.add_argument(
        "--no-extract",
        action="store_true",
        help="Only download the tarball; skip extraction.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing tarball/extracted directory if present.",
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    out_dir = Path(args.out_dir).expanduser().resolve()
    tar_path = out_dir / args.filename
    downloaded = download_bigbind_tar(url=args.url, output_tar=tar_path, force=args.force)

    if args.no_extract:
        print(f"[ml] Download complete: {downloaded}")
        return

    extract_root = extract_bigbind_tar(
        tar_path=downloaded,
        destination_dir=out_dir,
        force=args.force,
    )
    print(f"[ml] Ready directory for config.bigbind_dir: {extract_root}")


if __name__ == "__main__":
    main()
