# Modeling Disorder in a Binary Fe–20 at.% Al Alloy 
### Yell Refinemnet $\rightarrow$ Histogram Reweighting $\rightarrow$ Monte Carlo $\rightarrow$ uMLIP Relaxation

This repository contains the pipeline developed to simulate the structural disorder of a single-crystal iron-aluminium alloy (composition 80% Fe and 20% Al). 

## Project Overview & Rationale

Extracting real-space structural information from single-crystal diffuse scattering experiments is historically difficult due to the overlapping signatures of chemical short-range order (SRO), static size-effect displacements, and thermal diffuse scattering (TDS).

This pipeline works to untangle these signatures by isolating and subsequently modeling how SRO and the size-effect distort the reciprocal space mapping.

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
* **Atomic Simulation Environment (ASE):** High-level atomistic manipulation and MD engine orchestration.
* **FAIRChem Core:** Drives the **Universal Models for Atoms (UMA)** foundation interatomic model (`uma-s-1p1`) on CUDA.
* **gemmi / gemmi_disorder:** High-performance crystallographic library used to compute the 3D FFT diffuse structural intensity grids.
* **H5Py:** Manages large-scale binary dataset storage for the ensemble profiles.

## Running the notebooks

Clone the repository and start Jupyter from the repository root:

```bash
git clone https://github.com/cynthia-wen/disordered-materials-fe-al-alloy.git
cd disordered-materials-fe-al-alloy
jupyter lab
```

The notebooks use paths relative to this root, for example `main/data/...`. VS Code uses the same root automatically through the workspace settings.

## Acknowledgments & Disclaimers

* **Supervision:** This work was conducted as a mobility project at **ETH Zürich** under the supervision of **Prof. Arkadiy Simonov** at the **Laboratory for Disordered Materials**. Very heartfelt thanks to Prof. Arkadiy Simonov, Dr. Johnathan Bullard, and the rest of the lab for their invaluable guidance, resources, and experimental data during this project. 
* **Errors:** All conceptual, mathematical, or code implementation mistakes contained within this codebase are entirely the author's responsibility.
* **AI Assistance:** Parts of the code development and report drafting were assisted by large language model tools *(Claude, Anthropic; Gemini, Alphabet)*
