import os
import subprocess
import numpy as np
from ase.io import read, write
from ase.optimize import BFGS

# specific to UMA 
from fairchem.core import FAIRChemCalculator
from fairchem.core import pretrained_mlip
HF_TOKEN = "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" # <-- replace with your Hugging Face token
os.environ["HF_TOKEN"] = HF_TOKEN

try:
    print("Authenticating with Hugging Face via CLI...")
    subprocess.run(
        ["huggingface-cli", "login", "--token", HF_TOKEN], 
        check=True, 
        capture_output=True, 
        text=True
    )
    print("Successfully logged in to Hugging Face!")
except subprocess.CalledProcessError as e:
    print(f"⚠️ Hugging Face CLI login failed: {e.stderr}")
    print("Attempting to proceed anyway using the environment variable...")

# setup file paths
input_cif_path = "ICSD_CollCode55869.cif" # <-- replace with your actual input CIF file path
output_cif_path = "icsd_output.cif"
trajectory_path = "trajectory.traj"

# read input structure
print(f"Reading input structure from: {input_cif_path}")
atoms = read(input_cif_path)

# configure UMA MLIP predictor
UMA_MODEL = "uma-s-1p2"
UMA_TASK = "omat"  # task for bulk materials
seed = np.random.randint(0, np.iinfo(np.int32).max, dtype=int)
predictor = pretrained_mlip.get_predict_unit(UMA_MODEL, device="cuda", seed=seed)
calculator = FAIRChemCalculator(predictor, task_name=UMA_TASK)

# attach the MLIP calculator to our atoms object
atoms.calc = calculator

# run geometry optimization using BFGS algorithm until maximum atomic force component drops below 0.01 eV/Å
print("Starting geometry optimization...")
opt = BFGS(atoms, trajectory=trajectory_path)
opt.run(fmax=0.01)

print(f"Saving relaxed structure to: {output_cif_path}")
write(output_cif_path, atoms)