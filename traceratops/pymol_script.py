# load_pdb_grid.py

from colorsys import hsv_to_rgb
from pathlib import Path

from pymol import cmd

"""

To run in pymol:

1- open pymol.

Write:
2- PyMol> run $TRACERATOPS/traceratops/traceratops/pymol_script.py. For example,
PyMol> run /home/marcnol/Repositories/traceratops/traceratops/pymol_script.py

where $TRACERATOPS has to point to the installation folder of traceratops

then point to your folder with PDB files, e.g.:
3- PyMol> load_pdb_grid /mnt/disk2/marcnol/mDropbox/Labbooks/Lab_book_marcnol/Projects/pyHiM/trace_simulations/PDBs

You should now see the structures in the main viewing window.

4- PyMol> color_all_barcodes()

for 96 barcodes

4- PyMol> spin_grid

to make structures spin


"""


def color_all_barcodes():

    # Find all atom names beginning with "B"
    names = sorted(
        {atom.name for atom in cmd.get_model("all").atom if atom.name.startswith("B")}
    )

    n = len(names)

    for i, atom_name in enumerate(names):

        r, g, b = hsv_to_rgb(i / n, 1.0, 1.0)

        color_name = f"clr_{atom_name}"

        cmd.set_color(color_name, [r, g, b])
        cmd.color(color_name, f"name {atom_name}")


def load_pdb_grid(folder=".", pattern="*.pdb", max_structures=100):
    """
    Load the first PDB files in a folder and display them in a PyMOL grid
    using a ball-and-stick representation.

    Parameters
    ----------
    folder : str
        Folder containing the PDB files.
    pattern : str
        Filename pattern, such as "*.pdb".
    max_structures : int
        Maximum number of structures to load.
    """
    folder_path = Path(folder).expanduser().resolve()

    # Sort alphabetically and retain only the first 100 files.
    pdb_files = sorted(folder_path.glob(pattern))[:max_structures]

    if not pdb_files:
        print(f"No PDB files found in: {folder_path}")
        return

    cmd.delete("all")

    # White background.
    cmd.bg_color("white")

    # Enable grid mode.
    cmd.set("grid_mode", 1)

    # Ball-and-stick appearance.
    cmd.set("stick_radius", 0.15)
    cmd.set("sphere_scale", 0.25)

    for grid_slot, pdb_file in enumerate(pdb_files, start=1):
        object_name = pdb_file.stem.replace("-", "_").replace(" ", "_")

        cmd.load(str(pdb_file), object_name)
        cmd.set("grid_slot", grid_slot, object_name)

        cmd.hide("everything", object_name)

        cmd.show("sticks", object_name)
        cmd.show("spheres", object_name)

        # Black bonds
        cmd.set("stick_color", "black", object_name)

        # Appearance
        cmd.set("stick_radius", 0.12, object_name)
        cmd.set("sphere_scale", 0.30, object_name)

        # Rainbow coloring from first bead to last bead
        cmd.spectrum(
            expression="index",
            palette="rainbow",
            selection=object_name,
        )

        # Atom-based colors.
        cmd.util.cbag(object_name)

    cmd.orient("all")
    cmd.zoom("all")

    print(f"Loaded {len(pdb_files)} structures from {folder_path}")


def spin_grid(axis="y", degrees_per_frame=2):
    """
    Create and play a continuously looping rotation movie.
    """
    number_of_frames = max(1, round(360 / degrees_per_frame))

    cmd.mclear()
    cmd.mset(f"1 x{number_of_frames}")

    for frame in range(1, number_of_frames + 1):
        cmd.mdo(frame, f"turn {axis}, {degrees_per_frame}")

    cmd.set("movie_loop", 1)
    cmd.mplay()


# Register the function as a PyMOL command.
cmd.extend("color_all_barcodes", color_all_barcodes)
cmd.extend("load_pdb_grid", load_pdb_grid)
cmd.extend("spin_grid", spin_grid)
