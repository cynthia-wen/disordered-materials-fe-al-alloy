"""Self-shell (origin) ADP for Fe from paired ideal/MD Fe-Al frames.
U_ab = <u_a u_b>, u = r_MD - r_ideal. Cubic site => U should be isotropic."""

import os, re, glob
import numpy as np
from ase.io import read

unrelaxed_folder = "additional/fe-al-10/fe-al_mc_structures/final_structures"
relaxed_folder   = "additional/fe-al-10/fe-al_md_structures/final_structures"
output_filename  = "additional/self_adp.txt"

# pair files by frame number
frame_no = lambda p: int(re.search(r"(\d+)", os.path.basename(p)).group(1))
unrelaxed = {frame_no(p): p for p in glob.glob(os.path.join(unrelaxed_folder, "*.cif"))}
relaxed   = {frame_no(p): p for p in glob.glob(os.path.join(relaxed_folder, "*.cif"))}
frames = sorted(set(unrelaxed) & set(relaxed))

pairs = [(read(unrelaxed[n]), read(relaxed[n])) for n in frames]
box = pairs[0][0].get_cell().lengths()
is_fe = np.stack([np.array(r.get_chemical_symbols()) == 'Fe' for _, r in pairs])  # (F, N); occupancy differs per frame

# u = MD position - ideal site, minimum image, rigid drift removed
u = np.stack([r.get_positions(wrap=False) - i.get_positions(wrap=False) for i, r in pairs])
u -= box * np.round(u / box)
u -= u.mean(axis=1, keepdims=True)

assert np.abs(u).max() < 1.0, "displacements too large - check atom ordering matches"

ue = u[is_fe]
U = (ue.T @ ue) / len(ue)
U_iso = np.trace(U) / 3.0
off = np.abs(U[[0, 0, 1], [1, 2, 2]]).max()

lines = ["=== SELF-SHELL ADP, Fe (Yell origin term) ===",
         f"{len(frames)} frames, {is_fe.shape[1]} atoms per frame, Fe fraction {is_fe.mean():.4f}",
         f"  U_iso = {U_iso:.6f} A^2    B_iso = {8*np.pi**2*U_iso:.4f} A^2",
         f"  U_xx/yy/zz = {U[0,0]:+.6f} {U[1,1]:+.6f} {U[2,2]:+.6f}",
         f"  U_xy/xz/yz = {U[0,1]:+.6f} {U[0,2]:+.6f} {U[1,2]:+.6f}",
         f"  isotropy: max|off-diag| / U_iso = {off/U_iso:.4f}"]

print("\n".join(lines))
open(output_filename, "w").write("\n".join(lines))