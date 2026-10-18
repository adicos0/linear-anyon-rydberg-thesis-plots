# Read the calculations as notebooks

Start at 01 and continue in order. Each notebook follows a numbered sequence:
**short thesis explanation → equation → visible source-code cell → output**.
The Hamiltonian-building, evolution, and gap-calculation code is shown inside
the notebooks, adapted from the uploaded scripts. It is not hidden behind
imports from the reusable model library.

The main sections explain the calculation in small steps. The appendix keeps
the longer symbolic derivations, dimensional checks, limiting cases, and
comparisons with the validated reference implementation. Thesis equation
numbers and original source filenames appear beside the relevant sections.

| Notebook | Original source | Default work |
|---|---|---|
| [01: Effective Hamiltonian and spectrum](01_effective_hamiltonian_and_spectrum.ipynb) | `heff2_E_vs_N.py`, `heff2_E_vs_N_thesis.py` | Small sector solves; load the saved ten-site spectrum |
| [02: Two-site populations](02_two_site_population_dynamics.ipynb) | `heff2_states.py` | Nine-state evolution; two-state and bright/dark checks |
| [03: Chain dynamics and phase tests](03_chain_dynamics_and_phase_tests.ipynb) | `heff2_spacetime.py` | Ten sites; single-excitation and two-particle density maps |
| [04: Effective parameters and tunnelling](04_effective_parameters_and_tunnelling.ipynb) | `amplitudes_vs_muU_thesis.py`, `mu_vs_A_over_U.py` | Finite sideband sums; parameter curves and pole handling |
| [05: Target anyon-Hubbard gaps](05_target_anyon_hubbard_gaps.ipynb) | `sweepGapsMulti.py`, `doApplyHam.py`, `replotGapsMulti.py` | Original matrix-free operator, one four-site solve, seven saved ten-site maps |

`make_active_bond_clean.py` remains a figure-processing utility. It removes small
annotations and crops an existing schematic; it is not a Hamiltonian calculation.
See its docstring and `reference_figures/rydberg_active_bond.pdf`. The generator of
that original drawing and the generators of Figures 10-17 were not located in the
source collection; these notebooks do not claim to reconstruct them.

## Install and open

From the repository root, with Python 3.11 or newer:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-notebooks.txt
python -m pip install jupyterlab
python -m ipykernel install --prefix .venv --name thesis-notebooks --display-name "Thesis notebooks"
jupyter lab
```

On Windows, use `.venv\Scripts\activate` in place of the `source` command.
Select the environment's Python kernel. The notebooks find the repository root
from either the root directory or `notebooks/`. No absolute author-machine paths
are required. `requirements-notebooks-tested.txt` records the complete validated
environment (including platform-dependent packages); the portable input is
`requirements-notebooks.txt`.

Run each notebook top to bottom. Source arrays are never overwritten.
`SAVE_FIGURES=False` by default; setting it to true writes to the ignored
`notebook_outputs/` directory. Notebook 05 runs one four-site solve rather than
launching a grid sweep. The original large scripts are preserved and should only
be launched deliberately after reviewing their computational cost.

The numbered explanation–code–output layout was guided by the user's
[AQI notebook example](https://github.com/Charub17/global-air-pollution-aqi-ml/blob/main/Global_Air_Pollution_AQI_ML.ipynb).
The physics text and equations come from the thesis; the numerical code comes
from this repository's original scripts. Code-cell provenance metadata records
the source file and the declared adaptation.

## Verify

```bash
python -m pytest -q
python tools/execute_notebooks.py
```

The execution tool validates notebook schemas, starts a separate kernel per
notebook, fails on cell errors, and updates executed outputs. It does not run the
large archived calculations. Pass one or more notebook filenames to execute
only those notebooks. See the [validation record](VALIDATION.md) for the executed
cell counts, numerical comparisons, environment, and remaining limits.
Read [SCIENTIFIC_NOTES.md](../SCIENTIFIC_NOTES.md)
before using the plots to make a claim about experimental times or achievable
statistical phases.
