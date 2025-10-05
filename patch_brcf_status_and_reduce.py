# patch_brcf_status_and_reduce.py
import re
from pathlib import Path

root = Path.cwd()

# ---------- A) Patch automate_protein_prep.py (Reduce detection order & env var) ----------
ap = root / "automate_protein_prep.py"
ap_txt = ap.read_text(encoding="utf-8")

reduce_pat = r"""REDUCE_EXE\s*=\s*_cfg\([^\n]+\)\nif not REDUCE_EXE:\n\s+if PHENIX_DIR.*?\n\s+else:\n\s+    REDUCE_EXE = shutil\.which\("reduce"\) or shutil\.which\("reduce\.exe"\) or "reduce""""
reduce_new = (
    'REDUCE_EXE = os.environ.get("REDUCE_BIN") or _cfg("REDUCE_EXE", "", "reduce_exe")\n'
    'if not REDUCE_EXE:\n'
    '    phenix_reduce_py = str(Path(str(PHENIX_DIR or "")).joinpath("reduce.python"))\n'
    '    if PHENIX_DIR and Path(phenix_reduce_py).exists():\n'
    '        REDUCE_EXE = phenix_reduce_py\n'
    '    elif PHENIX_DIR and Path(PHENIX_DIR, "reduce.exe").exists():\n'
    '        REDUCE_EXE = str(Path(PHENIX_DIR, "reduce.exe"))\n'
    '    else:\n'
    '        REDUCE_EXE = shutil.which("reduce") or shutil.which("reduce.exe") or "reduce"\n'
    'logging.info("Using Reduce at: %s", REDUCE_EXE)'
)

m = re.search(reduce_pat, ap_txt, re.DOTALL)
if not m:
    raise SystemExit("Could not locate Reduce discovery block in automate_protein_prep.py")
ap_txt = ap_txt[:m.start()] + reduce_new + ap_txt[m.end():]

ap.write_text(ap_txt, encoding="utf-8")
print("[OK] Patched Reduce detection in automate_protein_prep.py")

# ---------- B) Patch prep_ligands.py (status log helper, init, FAIL/OK rows) ----------
pl = root / "prep_ligands.py"
pl_txt = pl.read_text(encoding="utf-8")

# 1) Insert helper after QUARANTINE_DIRNAME line
q_idx = pl_txt.find('QUARANTINE_DIRNAME = "quarantine"')
if q_idx == -1:
    raise SystemExit('Could not find QUARANTINE_DIRNAME = "quarantine" in prep_ligands.py')
insert_after = pl_txt.find("\n", q_idx) + 1
helper_block = """
# === status log helper ===
from pathlib import Path as _PathForStatus  # harmless alias if Path already imported

def _append_prep_status(status_log_path: _PathForStatus, ligand_name: str, status: str, reason: str = "", relpath: str = ""):
    header = "ligand\\tstatus\\treason\\tpdbqt_rel\\n"
    if not status_log_path.exists():
        status_log_path.write_text(header, encoding="utf-8")
    with status_log_path.open("a", encoding="utf-8") as fh:
        fh.write("\\t".join([ligand_name, status, reason, relpath]) + "\\n")
# === end helper ===

"""
pl_txt = pl_txt[:insert_after] + helper_block + pl_txt[insert_after:]

# 2) Initialize status_log right after output_ligands_dir.mkdir(...)
init_pat = r"(output_ligands_dir\s*=\s*Path\(cfg\[\"OUTPUT_LIGANDS_DIR\"\]\)\.resolve\(\)\s*\n\s*output_ligands_dir\.mkdir\(parents=True,\s*exist_ok=True\)\s*)"
m = re.search(init_pat, pl_txt)
if not m:
    raise SystemExit('Could not find output_ligands_dir initialization in prep_ligands.py')
pl_txt = pl_txt[:m.end()] + '\n    status_log = output_ligands_dir / "ligand_prep_status.tsv"\n' + pl_txt[m.end():]

# 3) Replace the postcheck/quarantine block to append a FAIL row
#    Locate the block that starts with the big "if not pdbqt_path.exists()..." and ends before "continue"
fail_anchor = pl_txt.find('_log_malformed(pdbqt_path, "pdbqt_postcheck_fail_or_small")')
if fail_anchor == -1:
    raise SystemExit('Could not locate postcheck/quarantine block anchor in prep_ligands.py')
block_start = pl_txt.rfind("if ", 0, fail_anchor)
block_end = pl_txt.find("continue", fail_anchor)
if block_start == -1 or block_end == -1:
    raise SystemExit("Could not determine postcheck block boundaries in prep_ligands.py")

fail_block_new = """if not pdbqt_path.exists() or pdbqt_path.stat().st_size < 100 or not is_valid_ligand(pdbqt_path, log_dir=prepped_ligands_dir):
            _log_malformed(pdbqt_path, "pdbqt_postcheck_fail_or_small")
            quarantine = prepped_ligands_dir / QUARANTINE_DIRNAME
            quarantine.mkdir(exist_ok=True)
            try:
                if pdbqt_path.exists():
                    new_p = quarantine / pdbqt_path.name
                    pdbqt_path.replace(new_p)
                    _append_prep_status(status_log, sanitized.name, "FAIL", "postcheck_fail_or_small", str(new_p.relative_to(prepped_ligands_dir)))
                else:
                    _append_prep_status(status_log, sanitized.name, "FAIL", "postcheck_fail_or_small", "")
            except Exception:
                _append_prep_status(status_log, sanitized.name, "FAIL", "postcheck_fail_or_small", "")
            continue"""

pl_txt = pl_txt[:block_start] + fail_block_new + pl_txt[block_end+len("continue"):]

# 4) Append OK row immediately after the "Created PDBQT" log line
ok_anchor = 'logging.info(f"Created PDBQT: {pdbqt_path.name}")'
ok_idx = pl_txt.find(ok_anchor)
if ok_idx == -1:
    raise SystemExit('Could not find "Created PDBQT" log line in prep_ligands.py')
ok_insert = '\n        _append_prep_status(status_log, sanitized.name, "OK", "", str(pdbqt_path.relative_to(prepped_ligands_dir)))'
pl_txt = pl_txt[:ok_idx + len(ok_anchor)] + ok_insert + pl_txt[ok_idx + len(ok_anchor):]

pl.write_text(pl_txt, encoding="utf-8")
print("[OK] Added per-ligand status logging in prep_ligands.py")

print("\nAll patches applied.\nNext steps:")
print("  1) export REDUCE_BIN=/stor/home/mpg2352/atlas/tools/phenix/phenix-1.21.2-5419/build/bin/reduce.python")
print("     (or point to /stor/home/mpg2352/atlas/tools/reduce/reduce_src/reduce if you prefer the C++ build)")
print("  2) Re-run your toy Reduce test, then a 2–3 target smoke test.")
