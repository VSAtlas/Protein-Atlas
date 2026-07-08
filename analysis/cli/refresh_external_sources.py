from __future__ import annotations

import argparse
from pathlib import Path

from analysis.external.source_registry import refresh_configured_sources, write_registry_manifest
from analysis.io import load_config


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Optionally download configured Atlas external sources into a local cache.")
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--cache-dir", default=Path("data/external/cache"), type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    config = load_config(args.config)
    downloaded = refresh_configured_sources(config, args.cache_dir, overwrite=args.overwrite)
    write_registry_manifest(config, args.cache_dir / "source_registry_manifest.json")
    for name, path in downloaded.items():
        print(f"{name}\t{path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
