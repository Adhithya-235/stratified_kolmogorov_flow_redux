# stratified_kolmogorov_flow_redux

DNS, linear stability, and post‐processing code for a numerical study of the nonlinear dynamics of stratified Kolmogorov flow.  
Intended to run on the UNH high‐performance clusters *Plasma* and *Premise*. citeturn0view0

## Table of Contents

- [Installation](#installation)  
- [Usage](#usage)  
- [Directory Structure](#directory-structure)  
- [Dependencies](#dependencies)  
- [Simulation Drivers](#simulation-drivers)  
- [Stability Analysis](#stability-analysis)  
- [Post‐Processing](#post-processing)  
- [Jupyter Notebooks](#jupyter-notebooks)  
- [Utility Belt](#utility-belt)  
- [Interactive Processing](#interactive-processing)  
- [Batch Scripts](#batch-scripts)  
- [Contributing](#contributing)  
- [License](#license)  
- [Contact](#contact)  

---

## Installation

1. **Clone the repository**  
   ```bash
   git clone https://github.com/Adhithya-235/stratified_kolmogorov_flow_redux.git
   cd stratified_kolmogorov_flow_redux
   ```

2. **Install Dedalus 2.0**  
   Follow the official guide: https://dedalus-project.readthedocs.io/en/v2_master/

3. **Ensure MATLAB is available** (R2020b or later) on your $PATH for any MATLAB-based scripts.

4. **Simulation parameters** can be passed via command-line flags (use `--help`) or configured in your own batch scripts.

---




## Usage

1. **DNS Simulation**  
   - Run `python stratified_kolmogorov_flow.py [options]` to launch a spectral DNS of stratified Kolmogorov flow (use `--help` for flags).  

2. **Galerkin Model**  
   - Invoke `python galerkin_kolmogorov_flow.py [options]` for a reduced‐order spectral simulation (see `--help`).  

3. **Stability Analysis**  
   - Navigate to `linear_stability/` or `secondary_stability/` and run the corresponding scripts (`--help` available for command‑line options).

4. **Post‐Processing**  
   - Use `post_processing/` or `post_processing_galerkin/` to generate diagnostics, spectra, and flow visualizations.

5. **Interactive Exploration**  
   - Launch the Galerkin example notebook in `jupyter_notebooks/` for ad‑hoc analysis and visualization.

6. **Simulation parameters** (e.g., Reynolds number, stratification, resolution, output) can be specified via command‑line flags (DOCOPT) or inside Slurm batch scripts.

---


## Directory Structure

```
.
├── .gitignore
├── LICENSE
├── README.md
├── stratified_kolmogorov_flow.py      # DNS driver
├── galerkin_kolmogorov_flow.py        # Galerkin driver
├── *.sbatch                           # Slurm job submission scripts
├── linear_stability/                  # Linear stability routines
├── secondary_stability/               # Secondary stability analysis
├── post_processing/                   # DNS data post‐processing
├── post_processing_galerkin/          # Galerkin data post‐processing
├── interactive_processing/            # Scripts for ad‐hoc data analysis
├── jupyter_notebooks/                 # Example workflows & visualization
└── utility_belt/                      # Shared helper functions
```
citeturn0view0

---

## Dependencies

- **Python 3.8+** with  
  - `numpy`  
  - `scipy`  
  - `matplotlib`  
  - `mpi4py`  
  - `h5py`  
  - `jupyter`  
- **Dedalus 2.0** (framework for DNS; see https://dedalus-project.readthedocs.io/en/v2_master/)  
- **MATLAB R2020b+** (for any `.m` routines invoked)  
- **Slurm** scheduler (via `.sbatch` scripts)  
- **Git**  

---


## Simulation Drivers

### `stratified_kolmogorov_flow.py`
Direct numerical simulation (DNS) of stratified Kolmogorov flow using a Fourier–spectral method.  
Parameters to edit at the top of the file: grid resolution, Reynolds number, Froude/stratification parameters, time‐stepping options.

### `galerkin_kolmogorov_flow.py`
Galerkin‐projected, low‐dimensional spectral model. Good for parameter sweeps when full DNS is too costly.

---

## Stability Analysis

- **`linear_stability/`**  
  Computes the linear eigenvalue problem for infinitesimal perturbations about the base flow.

- **`secondary_stability/`**  
  Investigates instabilities of fully‐developed DNS / Galerkin states (Floquet or perturbation growth).

---

## Post‐Processing

- **`post_processing/`**  
  Scripts to read simulation outputs (e.g., HDF5 or MATLAB files), compute diagnostics (energy spectra, vorticity fields), and generate plots.

- **`post_processing_galerkin/`**  
  Same, but tailored to the reduced Galerkin output formats.

---

## Jupyter Notebooks

The `jupyter_notebooks/` folder contains an example notebook for the Galerkin model only, demonstrating setup, execution, and basic visualization of reduced‐order simulations.

---


## Utility Belt

Reusable helper functions (I/O, parallel wrappers, plotting utilities) kept in `utility_belt/`.

---

## Interactive Processing

Quick‐and‐dirty scripts in `interactive_processing/` for on‐the‐fly data slicing, correlation studies, or parametric scans without full batch runs.

---

## Batch Scripts

The provided Slurm `.sbatch` files (e.g. `marvin_*` for Plasma, `premise_*` for Premise) serve as templates.  
We recommend generating your own Slurm scripts tailored to your cluster configuration and simulation parameters rather than editing these examples.

---


## Contributing

1. Fork the repo  
2. Create a feature branch  
3. Submit a pull request with a clear description of changes  
4. Ensure all notebooks and scripts run without errors

---

## License

This project is licensed under the **GPL-3.0 License**. citeturn0view0

---

## Contact

Author: Adhithya Sivakumar  
GitHub: [@Adhithya-235](https://github.com/Adhithya-235)  
Feel free to open an issue for questions or feature requests.
