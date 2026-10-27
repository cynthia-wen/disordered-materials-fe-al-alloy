import os
import random
import numpy as np
from ase.build import bulk
from ase.io import write
from main.scripts.fitting_js.shell_functions import shells_from_model, shell_vectors

SUPERCELL_SIZE = (20, 20, 20)
MC_STEPS = 14000000
T = [0.2]

# parameters from fit_js_multi
J_SPIN = [0.03300626377699999,
 0.020652299734900004,
 -0.005209546949,
 -8.937613100000013e-05,
 0.0019287206729999993,
 0.0035445230747999995,
 0.006071197146,
 0.001830061483,
 -0.00243519570697,
 -0.0027261006571999997,
 -0.005501055883509998,
 -0.0013909532175000003]
NUM_SHELLS = len(J_SPIN)
NUM_RUNS = 300

# initialize structure
a_lattice = 2.867
atoms_base = bulk("Fe", "bcc", a=a_lattice, cubic=True).repeat(SUPERCELL_SIZE)
num_atoms = len(atoms_base)
grid_n = 2 * round(atoms_base.cell.lengths()[0] / a_lattice)   # doubled-coordinate grid

shell_keys, shell_r = shells_from_model('main/scripts/fitting_js/model.txt', NUM_SHELLS)

# get doubled-integer site of every atom, and the reverse map (coord -> atom index)
site_coords = np.round(atoms_base.get_scaled_positions() * grid_n).astype(int) % grid_n
index_grid = np.full((grid_n,) * 3, -1, dtype=int)
index_grid[tuple(site_coords.T)] = np.arange(num_atoms)

# neighbours of every atom at once: shifting the index grid by -v puts the atom at k+v under k
shell_lookups = [
    np.stack([np.roll(index_grid, (-v[0], -v[1], -v[2]), (0, 1, 2))[tuple(site_coords.T)]
              for v in sorted(shell_vectors(rep))], axis=1)          # (num_atoms, multiplicity)
    for rep in shell_keys]

ALL_NEIGHBORS = np.hstack(shell_lookups) # (num_atoms, 136)
J_PER_COLUMN  = np.concatenate([np.full(shell_lookups[n].shape[1], J_SPIN[n]) for n in range(NUM_SHELLS)])  # (136,)

# main loop
for temp in T:
    OUTPUT_DIR = "main/data/structures_20/refined_10"
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print(f"Starting Monte Carlo simulations at T={temp:.3f}")

    for run_idx in range(1, NUM_RUNS + 1):
        print(f" -> Processing Run {run_idx}/{NUM_RUNS}...", end="\r")
        
        # generate a unique + random initial distribution of 80-20
        num_al = int(0.20 * num_atoms)
        num_fe = num_atoms - num_al
        spins = np.array([1.0] * num_fe + [-1.0] * num_al)
        np.random.shuffle(spins)

        # mc
        for step in range(MC_STEPS):
            idx_a = random.randint(0, num_atoms - 1)
            idx_b = random.randint(0, num_atoms - 1)
            while idx_a == idx_b:
                idx_b = random.randint(0, num_atoms - 1)
            
            if spins[idx_a] != spins[idx_b]:
                s_a, s_b = spins[idx_a], spins[idx_b] 
                nb_a, nb_b = ALL_NEIGHBORS[idx_a], ALL_NEIGHBORS[idx_b] # get neighbor indices for a and b
                sum_a = J_PER_COLUMN @ spins[nb_a] - s_b * J_PER_COLUMN[nb_a == idx_b].sum() # spin-sum J*s of a, minus the contribution from b if b is one of a's neighbors
                sum_b = J_PER_COLUMN @ spins[nb_b] - s_a * J_PER_COLUMN[nb_b == idx_a].sum() # vice versa
                total_dE = (s_b - s_a) * (sum_a - sum_b) # swapped - original = (s_b*sum_a + s_a*sum_b) - (s_a*sum_a + s_b*sum_b)

                if total_dE < 0 or random.random() < np.exp(-total_dE / temp):
                    spins[idx_a], spins[idx_b] = s_b, s_a
                    
        # save cif file structural configuration
        final_atoms = atoms_base.copy()
        final_symbols = ['Fe' if x == 1.0 else 'Al' for x in spins]
        final_atoms.set_chemical_symbols(final_symbols)
        cif_filename = os.path.join(OUTPUT_DIR, f"unrelaxed_mc_structure_{run_idx:02d}.cif")
        write(cif_filename, final_atoms)

print(f"\nAll {NUM_RUNS} simulations completed successfully!")