"""Run clean_pipeline.ipynb for MC_INDEX 001-300, then average the diffuse
scattering over all final structures.

Usage:
    python run_pipeline.py                    # run indices 1-300 sequentially, then average
    python run_pipeline.py --jobs 4           # run 4 indices at a time, then average
    python run_pipeline.py --index 42         # run a single index, no averaging (for job arrays)
    python run_pipeline.py --average-only     # skip the runs, just (re)do the averaging
"""
import argparse
import glob
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import gemmi
from ase.io import read
from gemmi_disorder import Grid, DisorderedStructure, average_diffuse, save_to_yell

NOTEBOOK = "additional/clean_pipeline.ipynb"
SUPERCELL_SIZE = (10, 10, 10)  # must match the notebook

# (label, per-run structures glob, aggregated save_to_yell output dir)
AVERAGE_TARGETS = [
    ("mc", "additional/fe-al-10/fe-al_mc_structures/final_structures/mc_final_*.cif", "additional/fe-al-10/fe-al_mc_structures/average_diffuse"),
    ("relaxed", "additional/fe-al-10/fe-al_relaxed_structures/final_structures/relaxed_frame_*.cif", "additional/fe-al-10/fe-al_relaxed_structures/average_diffuse"),
    ("md", "additional/fe-al-10/fe-al_md_structures/final_structures/md_final_*.cif", "additional/fe-al-10/fe-al_md_structures/average_diffuse"),
]


def run_index(i):
    print(f"=== running pipeline for index {i:03d} ===", flush=True)
    env = dict(os.environ, MC_INDEX=str(i))
    subprocess.run(
        ["jupyter", "nbconvert", "--to", "notebook", "--execute", "--stdout",
         "--ExecutePreprocessor.timeout=-1", NOTEBOOK],
        env=env, check=True, stdout=subprocess.DEVNULL)


def run_indices(start, end, jobs):
    indices = range(start, end + 1)
    if jobs <= 1:
        for i in indices:
            run_index(i)
    else:
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = [pool.submit(run_index, i) for i in indices]
            for f in futures:
                f.result()  # re-raise any failure


def structure_from_cif(path):
    # convert a supercell cif back to a DisorderedStructure in unit-cell units
    frame = read(path)
    symbols = frame.get_chemical_symbols()
    supercell_positions = frame.get_scaled_positions()
    pars = frame.cell.cellpar()
    nx, ny, nz = SUPERCELL_SIZE
    cell_pars_unit = (pars[0] / nx, pars[1] / ny, pars[2] / nz, pars[3], pars[4], pars[5])
    unit_positions = supercell_positions * np.array([nx, ny, nz])
    atoms_list = list(zip(symbols, unit_positions[:, 0], unit_positions[:, 1], unit_positions[:, 2]))
    return DisorderedStructure(cell_parameters=cell_pars_unit, atoms=atoms_list, supercell=SUPERCELL_SIZE)


def average_all():
    for label, pattern, out_dir in AVERAGE_TARGETS:
        cif_paths = sorted(glob.glob(pattern))
        if not cif_paths:
            print(f"[{label}] no structures matching {pattern} -- skipping", flush=True)
            continue
        print(f"[{label}] averaging diffuse scattering over {len(cif_paths)} structures", flush=True)
        structures = [structure_from_cif(p) for p in cif_paths]

        grid = Grid.from_supercell(supercell=SUPERCELL_SIZE, hkl_max=4)
        result = average_diffuse(structures, grid, blur=0.04, crop=4, b_iso=0.0, progress=True)

        cell = gemmi.UnitCell(*structures[0].cell_parameters)
        os.makedirs(out_dir, exist_ok=True)
        save_to_yell(result, cell=cell, supercell=SUPERCELL_SIZE, prefix=out_dir)
        print(f"[{label}] saved averaged diffuse pattern to {out_dir}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int, default=1, help="first index (default 1)")
    parser.add_argument("--end", type=int, default=300, help="last index (default 300)")
    parser.add_argument("--jobs", type=int, default=1, help="number of indices to run concurrently")
    parser.add_argument("--index", type=int, help="run a single index and exit (no averaging)")
    parser.add_argument("--average-only", action="store_true", help="skip the runs, just average")
    args = parser.parse_args()

    if args.average_only:
        average_all()
    elif args.index is not None:
        run_index(args.index)
    else:
        run_indices(args.start, args.end, args.jobs)
        average_all()
