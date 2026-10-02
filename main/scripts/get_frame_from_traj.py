import os
from ase.io import read, write

# Define the source template and destination directory
SOURCE_DIR_TEMPLATE = "main/data/testing/multi/structures_T0.010/"
DEST_DIR = "main/data/testing/multi/unrelaxed_structures_T0.010"

# Create the destination folder if it doesn't exist
os.makedirs(DEST_DIR, exist_ok=True)

print(f"Extracting the first frame from trajectories 1 to 50...")

# Loop through the specific folder numbers
for i in range(1, 51):
    # Construct the path to the current traj file
    traj_path = os.path.join(SOURCE_DIR_TEMPLATE, f"mc-relaxed-traj_{i}.traj")
    
    if os.path.exists(traj_path):
        # index=0 explicitly extracts only the first frame
        first_frame = read(traj_path, index=0)
        
        # Save it to the destination folder as a CIF or any preferred format
        output_filename = os.path.join(DEST_DIR, f"unrelaxed_structure_{i}.cif")
        write(output_filename, first_frame)
        print(f" -> Extracted frame 0 from run {i} -> {output_filename}")
    else:
        print(f" ⚠️ Warning: File not found at {traj_path}")

print("\nExtraction complete! All unrelaxed structures are saved in:", DEST_DIR)