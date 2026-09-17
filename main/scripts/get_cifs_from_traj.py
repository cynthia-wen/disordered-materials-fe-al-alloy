import os
from ase.io import read, write

# =====================================================================
# CONFIGURATION
# =====================================================================
# Path to your MD or MC trajectory file
trajectory_path = "fe-al_mc_structures_relaxed/mc-relaxed-output_192/mc_relaxed.traj"

# Directory where you want to save the individual CIFs
output_dir = "fe-al_mc_structures_relaxed/extracted_cifs_192"

# =====================================================================
# EXTRACTION PROCESS
# =====================================================================
# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

print(f"Reading trajectory file: {trajectory_path}")
# The ':' slice notation tells ASE to read ALL frames into a list
configs = read(trajectory_path, index=":")
total_frames = len(configs)
print(f"Found {total_frames} total frames. Starting CIF export...")

# Loop through every frame and save it
for idx, atoms in enumerate(configs):
    # Generates zero-padded names like frame_002.cif so they sort perfectly alphabetically
    cif_filename = f"frame_{idx:03d}.cif"
    cif_filepath = os.path.join(output_dir, cif_filename)
    
    # Write the frame
    write(cif_filepath, atoms)
    
    # Progress update every 50 frames so you know it's working
    if idx % 50 == 0 or idx == total_frames - 1:
        print(f" -> Exported {idx + 1}/{total_frames}: {cif_filename}")

print(f"\n🎉 Success! All {total_frames} frames have been saved as CIFs in:\n👉 {output_dir}")