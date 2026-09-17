import os
import subprocess
import numpy as np

from ase.io import read
from ase.optimize import LBFGS

from fairchem.core import FAIRChemCalculator
from fairchem.core import pretrained_mlip

final_dir = os.path.join("fitting_js/structures_20/relaxed_08")
os.makedirs(final_dir, exist_ok=True)

SUPERCELL_SIZE = (20, 20, 20)

UMA_MODEL = "uma-s-1p2"
UMA_TASK = "omat"  # task for bulk materials
seed = np.random.randint(0, np.iinfo(np.int32).max, dtype=int)
predictor = pretrained_mlip.get_predict_unit(UMA_MODEL, device="cuda", seed=seed)
calculator = FAIRChemCalculator(predictor, task_name=UMA_TASK)

for i in range(120, 156):
    index = i #f"100{i}"
    final_filename = os.path.join(final_dir, f"relaxed_frame_{index}.cif")

    if os.path.exists(final_filename):
        print(f"[{i}/300] {index}: already relaxed, skipping")
        continue

    print(f"[{i}/300] relaxing mc{index}_structure.cif")
    atoms = read(f"fitting_js/structures_20/refined_08/unrelaxed_mc_structure_{index}.cif")

    atoms_relaxed = atoms.copy()
    atoms_relaxed.calc = calculator

    relax_traj_filename = os.path.join(
        f"fitting_js/structures_20/relaxed_08/trajectories/structure_{index}.traj"
    )
    opt = LBFGS(atoms_relaxed, trajectory=relax_traj_filename)
    opt.run(fmax=0.01)

    atoms_relaxed.write(final_filename)