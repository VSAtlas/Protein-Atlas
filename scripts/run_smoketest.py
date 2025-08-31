import sys, subprocess, pathlib

def main():
    out = pathlib.Path("smoketest_out"); out.mkdir(exist_ok=True)
    # TODO: tweak this to your repo’s real command:
    cmd = ["python3","main.py","--pdb","1IEP","--stages","stage1,stage2",
           "--library","data/fda_mini","--out", str(out)]
    print("Running:", " ".join(cmd), flush=True)
    rc = subprocess.call(cmd)
    if rc != 0:
        sys.exit(rc)
    # Example success check:
    expected = out / "docking_score_summary.csv"
    if not expected.exists():
        print(f"Missing {expected}", file=sys.stderr)
        sys.exit(2)
    print("smoke OK")

if __name__ == "__main__":
    main()
