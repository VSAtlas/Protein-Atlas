import sys
from pathlib import Path

def split_pdbqt_roots(pdbqt_path: Path, out_dir: Path) -> None:
    text = pdbqt_path.read_text()
    # We split on 'ROOT\n', but keep the header (anything before the first ROOT).
    parts = text.split("ROOT\n")
    if len(parts) <= 2:
        # 0 or 1 ROOT block; nothing to split
        return

    header = parts[0]
    # The remainder are individual ROOT sections; reattach 'ROOT\n' to each.
    for i, body in enumerate(parts[1:], start=1):
        # Strip trailing whitespace so we don't accumulate blank lines
        body = body.strip()
        if not body:
            continue

        out_name = f"{pdbqt_path.stem}_split{i}.pdbqt"
        out_path = out_dir / out_name
        out_text = header + "ROOT\n" + body + "\n"
        out_path.write_text(out_text)
        print(f"Wrote {out_path}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python split_pdbqt_roots.py <pdbqt1> [<pdbqt2> ...]")
        raise SystemExit(1)

    for arg in sys.argv[1:]:
        pdbqt_path = Path(arg)
        if not pdbqt_path.is_file():
            print(f"Skipping non-file {pdbqt_path}")
            continue

        out_dir = pdbqt_path.parent / "split"
        out_dir.mkdir(parents=True, exist_ok=True)
        split_pdbqt_roots(pdbqt_path, out_dir)

if __name__ == "__main__":
    main()
