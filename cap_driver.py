import sys
from pymol import cmd
from capture_pose import capture_pose

def _posargs(argv):
    """Return only non-option args (drop PyMOL flags like -c, -q, -r, -d, etc.)."""
    out = []
    skip_next = False
    for a in argv:
        if skip_next:
            skip_next = False
            continue
        if a.startswith('-'):
            # Options that consume a parameter (be conservative)
            if a in ('-d','-r','-u','-p','-l','-R'):
                skip_next = True
            continue
        out.append(a)
    return out

args = _posargs(sys.argv[1:])
if len(args) < 3:
    print(f"usage: cap_driver.py <receptor> <ligand> <outdir-or-prefix>  (argv={sys.argv})")
    cmd.quit(2)

# In case extra non-option args were passed, take the last three as rec, lig, out
rec, lig, out = args[-3:]
capture_pose(rec, lig, out)
cmd.quit()
