# Modeling Disorder in a Binary Fe–20 at.% Al Alloy 
### Yell Refinement $\rightarrow$ Histogram Reweighting $\rightarrow$ Monte Carlo $\rightarrow$ uMLIP Relaxation

This repository contains the pipeline developed to simulate the structural disorder of a single-crystal iron-aluminium alloy (composition 80% Fe and 20% Al).

## Project Overview & Rationale

Extracting real-space structural information from single-crystal diffuse scattering experiments is historically difficult due to the overlapping signatures of chemical short-range order (SRO), static size-effect displacements, and thermal diffuse scattering (TDS).

This pipeline separates the problem rather than fitting all contributions simultaneously: substitutional correlations are refined from experiment, while displacements are generated physically through relaxation.

### Computational Architecture

The workflow progresses across a statistically sampled ensemble of 300 independent structures to eliminate stochastic noise and capture ensemble-averaged behavior:

1. 3D-ΔPDF Refinement: Substitutional correlations are refined against synchrotron data in Yell (365 parameters, 220 shells, R_w = 26%)
2. Histogram Reweighting: 12 shell-resolved interaction parameters (Js) fitted to the refined correlations; damped Levenberg–Marquardt updates until ΔR_norm < 10⁻⁴
3. Initialization: 20×20×20 BCC supercell (16,000 atoms) + stochastic substitution, 80% Fe / 20% Al
4. Rigid-Lattice Monte Carlo: 14M steps with the fitted Js at T = 0.2 (reduced units)
5. Structural Relaxation: UMA foundation model (uma-s-1p2, OMat task head) + L-BFGS minimizer (force convergence < 0.01 eV/Å)
6. Diffuse Scattering Calculation: Ensemble-averaged FFT profiles via gemmi_disorder, m-3m symmetrization, (hk0) and five (xkl) slices

## Software Stack & Environment

All code was run on an NVIDIA RTX A6000 GPU and leverages a Python-based ecosystem with the following libraries:
* **Atomic Simulation Environment (ASE):** High-level atomistic manipulation and structure generation.
* **FAIRChem Core:** Drives the **Universal Models for Atoms (UMA)** foundation interatomic model (`uma-s-1p2`) on CUDA.
* **gemmi / gemmi_disorder:** High-performance crystallographic library used to compute the 3D FFT diffuse structural intensity grids.
* **H5Py:** Manages large-scale binary dataset storage for the ensemble profiles.
* **Yell (1.3.79):** 3D-ΔPDF refinement of substitutional correlations against the experimental diffuse scattering, using analytical derivatives and the Ceres least-squares minimizer. (Non-Python software)
* **Meerkat (0.39):** Reciprocal-space reconstruction of the raw diffraction frames, Bragg peak punch-and-fill, and m-3m Laue symmetrization of both experimental and simulated maps. (Non-Python software)

## Repository Structure

```text
disordered-materials-fe-al-alloy/
├── main/ # production workflow (20×20×20, 16,000 atoms)
│ ├── yell_refinement/ # Yell model + refined substitutional correlations
│ │
│ ├── scripts/
│ │ ├── fitting_js/
│ │ │ ├── histogram_reweighting.ipynb # susceptibility, damped LM updates, J refinement
│ │ │ ├── run_mc_multi.py # Metropolis MC over the 300-structure ensemble
│ │ │ └── shell_functions.py # shell vectors and multiplicities from model.txt
│ │ ├── relax_single.py # UMA + L-BFGS relaxation, one structure
│ │ ├── relax_multiple.py # same, across the ensemble
│ │ ├── get_energy_from_structure.py # UMA energies for the least-squares seed fit
│ │ ├── get_cifs_from_traj.py # extract CIFs from relaxation trajectories
│ │ ├── get_frame_from_traj.py # extract a single frame
│ │ ├── calc_average_from_cifs.py # ensemble-averaged diffuse maps (gemmi_disorder)
│ │ ├── mask.py # Bragg masking / punch-and-fill
│ │ └── figures.ipynb # all figures in the report and poster
│ │
│ └── data/structures_20/
│ ├── refined_01/ … refined_09/ # one directory per reweighting iteration
│ │ ├── js_from_yell/ # J values and per-shell correlation comparison
│ │ └── average_diffuse/ # diffuse maps for that iteration's ensemble
│ └── relaxed_08/ # UMA-relaxed ensemble from the final parameters
│ └── average_diffuse/ # simulated maps reported in the write-up
│
└── additional/ # preliminary 10×10×10 work (2,000 atoms)
├── run_pipeline.py # MC → relax → MD driver
├── clean_pipeline.ipynb
├── fit_js_multi.ipynb # initial least-squares estimate of the Js
├── calc_adp.py # ADP / size-effect correlations from trajectories
└── fe-al-10/
├── fe-al_mc_structures/ # Monte Carlo only
├── fe-al_relaxed_structures/ # + UMA relaxation
└── fe-al_md_structures/ # + molecular dynamics
```

## Running the notebooks

Clone the repository and start Jupyter from the repository root:

```bash
git clone https://github.com/cynthia-wen/disordered-materials-fe-al-alloy.git
cd disordered-materials-fe-al-alloy
jupyter lab
```

The notebooks use paths relative to this root, for example `main/data/...`. VS Code uses the same root automatically through the workspace settings.

### Execution order
1. `main/yell_refinement/` — refine substitutional correlations against the experimental data (Yell, run externally)
2. `additional/fit_js_multi.ipynb` — least-squares seed for the interaction parameters, using energies from `get_energy_from_structure.py`
3. `main/scripts/fitting_js/histogram_reweighting.ipynb` — iterate the J refinement, calling `run_mc_multi.py` at each step; results land in `main/data/structures_20/refined_NN/`
4. `main/scripts/relax_multiple.py` — relax the final ensemble with UMA
5. `main/scripts/calc_average_from_cifs.py` — ensemble-averaged diffuse scattering maps
6. `main/scripts/figures.ipynb` — reproduce the report and poster figures

The `additional/` tree holds the preliminary 10×10×10 ensembles used only to obtain
the initial estimate of the interaction parameters; they are not part of the
production pipeline and are retained for reference.

## Acknowledgments & Disclaimers

* **Supervision:** This work was conducted as a mobility project at **ETH Zürich** under the supervision of **Prof. Arkadiy Simonov** at the **Laboratory for Disordered Materials**. Very heartfelt thanks to Prof. Arkadiy Simonov, Dr. Johnathan Bulled, and the rest of the lab for their invaluable guidance, resources, and experimental data during this project. 
* **Errors:** All conceptual, mathematical, or code implementation mistakes contained within this codebase are entirely the author's responsibility.
* **AI Assistance:** Parts of the code development and report drafting were assisted by large language model tools *(Claude, Anthropic; Gemini, Alphabet)*
