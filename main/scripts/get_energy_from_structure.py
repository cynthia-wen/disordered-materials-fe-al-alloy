import os
import numpy as np
from ase.io import read
from fairchem.core import FAIRChemCalculator
from fairchem.core import pretrained_mlip

# 1. Path to your CIF structure
cif_path = "main/data/fe-al_starting/relaxation/final_structures/relaxed_frame_192.cif"

if not os.path.exists(cif_path):
    raise FileNotFoundError(f"Could not find the file: {cif_path}")

print(f"Loading structure from {cif_path}...")
atoms = read(cif_path)

# 2. Attach a Universal Machine Learning Interatomic Potential (MLIP)
print("Initializing Machine Learning potential calculator...")
UMA_MODEL = "uma-s-1p1"
UMA_TASK = "omat" # task for bulk materials
seed = np.random.randint(0, np.iinfo(np.int32).max, dtype=int)
predictor = pretrained_mlip.get_predict_unit(UMA_MODEL, device="cuda", seed=seed)
calculator = FAIRChemCalculator(predictor, task_name=UMA_TASK)
atoms.calc = calculator

# 3. Compute and retrieve the total potential energy
total_energy = atoms.get_potential_energy()
energy_per_atom = total_energy / len(atoms)

print("\n========================================")
print(f"Total Atoms:       {len(atoms)}")
print(f"Total Energy:      {total_energy:.4f} eV")
print(f"Energy per Atom:   {energy_per_atom:.4f} eV/atom")
print("========================================")