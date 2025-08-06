import sys
from pymol2 import PyMOL

receptor = sys.argv[1]
ligand = sys.argv[2]
outprefix = sys.argv[3]
TOP_N_RESIDUES = 5

with PyMOL() as pymol:
    cmd = pymol.cmd

    cmd.load(receptor, "receptor")
    cmd.load(ligand, "ligand")

    cmd.hide("everything")
    cmd.color("orange", "ligand")
    cmd.show("sticks", "ligand")

    # Active site: residues within 5 Å
    cmd.select("active_site_all", "receptor within 5 of ligand")

    # Rank residues by proximity to ligand CA atoms
    distances = []
    cmd.iterate("active_site_all and name CA",
        "distances.append((model, chain, resi, resn, cmd.distance('tmp', 'ligand', f'{model}//{chain}/{resi}/CA')))",
        space={"distances": distances, "cmd": cmd})
    cmd.delete("tmp")

    # Sort and select top N
    top_residues = sorted(distances, key=lambda x: x[4])[:TOP_N_RESIDUES] if distances else []
    top_sel = " or ".join([f"receptor and chain {c} and resi {r}" for _, c, r, _, _ in top_residues])

    if top_sel:
        cmd.select("top_site", top_sel)
        cmd.show("sticks", "top_site")
        cmd.color("cyan", "top_site")
        cmd.label("top_site and name CA", "resn + resi")

    # Show full receptor as transparent surface
    cmd.show("surface", "receptor")
    cmd.set("transparency", 0.3, "receptor")
    cmd.color("slate", "receptor")

    # View and save
    cmd.zoom("ligand or top_site", 10)
    cmd.viewport(800, 600)

    cmd.png(f"{outprefix}_front.png", ray=1)
    cmd.turn("y", 90)
    cmd.png(f"{outprefix}_side.png", ray=1)
    cmd.turn("x", 90)
    cmd.png(f"{outprefix}_top.png", ray=1)
