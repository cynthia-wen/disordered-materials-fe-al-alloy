# Hierarchically Modeling Structural Disorder in an Fe-Al Alloy
### Monte Carlo $\rightarrow$ uMLIP Relaxation $\rightarrow$ Molecular Dynamics Pipeline

This repository contains the pipeline developed to simulate the structural disorder of a single-crystal iron-aluminium alloy (composition 80% Fe and 20% Al). 

## Project Overview & Rationale

Extracting real-space structural information from single-crystal diffuse scattering experiments is historically difficult due to the overlapping signatures of chemical short-range order (SRO), static size-effect displacements, and thermal diffuse scattering (TDS).

This pipeline works to untangle these signatures by sequentially activating physical degrees of freedom to isolate how each structural constraint distorts the reciprocal space mapping.

### Computational Architecture

The workflow progresses across a statistically sampled ensemble of 300 independent structures to eliminate stochastic noise and capture macroscopic thermodynamic behavior:

1. Initialization: 10x10x10 BCC Supercell (2,000 Atoms) + Stochastic Substitution: 80% Fe / 20% Al
2. Rigid-Lattice Monte Carlo: 14M Steps, J = 1, T = 4.2
3. Structural Relaxation: UMA Foundation Model (uMLIP) + BFGS Minimizer (Force Convergence Criterion: < 0.01 eV/Å)
4. Molecular Dynamics: NVT Ensemble at 300 K via Langevin Thermostat (Timestep: 5.0 fs, Total Run Time: 15 ps)
5. Diffuse Scattering Calculation: Ensemble-Averaged FFT Profiles via gemmi_disorder + (hk0) Slicing

## Software Stack & Environment

All code was run on an NVIDIA RTX A6000 GPU and leverages a Python-based ecosystem with the following libraries:
* **Atomic Simulation Environment (ASE):** High-level atomistic manipulation and MD engine orchestration.
* **FAIRChem Core:** Drives the **Universal Models for Atoms (UMA)** foundation interatomic model (`uma-s-1p1`) on CUDA.
* **gemmi / gemmi_disorder:** High-performance crystallographic library used to compute the 3D FFT diffuse structural intensity grids.
* **H5Py:** Manages large-scale binary dataset storage for the ensemble profiles.

## Acknowledgments & Disclaimers

* **Supervision:** This work was conducted as a mobility project at **ETH Zürich** under the supervision of **Prof. Arkadiy Simonov** at the **Laboratory for Disordered Materials**. Very heartfelt thanks to Prof. Arkadiy Simonov, Dr. Johnathan Bullard, and the rest of the lab for their invaluable guidance, resources, and experimental data during this project. 
* **Errors:** All conceptual, mathematical, or code implementation mistakes contained within this codebase are entirely the author's responsibility.
* **AI Assistance:** Parts of the code development and report drafting were assisted by large language model tools *(Gemini, Google)*
