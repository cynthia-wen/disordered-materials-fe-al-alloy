import glob
import gemmi
from gemmi_disorder import Grid, average_diffuse, save_to_yell

SUPERCELL_SIZE = (20, 20, 20)

# combine all final frames into average diffuse scattering
cif_paths = sorted(glob.glob("main/data/structures_20/relaxed_08/*.cif"))
structures = [gemmi.read_small_structure(p) for p in cif_paths]

grid = Grid.from_supercell(supercell=SUPERCELL_SIZE, hkl_max=5)
result = average_diffuse(structures, grid, blur=0.04, crop=2, b_iso=0.24, progress=True)
save_to_yell(result, cell=structures[0].cell, supercell=SUPERCELL_SIZE, prefix="main/data/structures_20/relaxed_08/average_diffuse")