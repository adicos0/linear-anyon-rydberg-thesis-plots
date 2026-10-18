# Read the thesis calculations as notebooks

Each notebook follows **thesis explanation → thesis equation → original code →
output**. Physics paragraphs paraphrase the cited thesis sections, with printed
page numbers and original equation numbers. Original calculation cells have
links to their exact archived script lines and a record of notebook adaptations.

| Notebook | Thesis source | Default calculation |
|---|---|---|
| [01: Effective Hamiltonian and ED](01_effective_hamiltonian_and_spectrum.ipynb) | Sections 3.6–3.7, pp. 20–22 | Original 55-state N=2 solve; archived full spectrum; source checks and plot |
| [02: Two-site populations](02_two_site_population_dynamics.ipynb) | Section 3.8, pp. 23–25 | Original nine-state evolution and six population panels |
| [03: Chain occupations](03_chain_dynamics_and_phase_tests.ipynb) | Section 3.8.2, pp. 25–28 | Original ten-site single-excitation and doublon dynamics; six target-phase maps |
| [04: Parameters and tunnelling](04_effective_parameters_and_tunnelling.ipynb) | Sections 4.1–4.4, pp. 29–33 | Original sideband sums and the two uploaded parameter-plot families |
| [05: Target-model gaps](05_target_anyon_hubbard_gaps.ipynb) | Sections 1.2–1.3 and Appendix C | Original matrix-free solver at a four-site example; seven archived ten-site gap maps |

Notebook 01 includes a dedicated explanation of the effective Hamiltonian's
exact diagonalization: basis configurations, number conservation, matrix
assembly within each block, dense eigenvalue calculation, spectrum assembly,
and the original script's checks. Its method differs from Appendix C's
matrix-free target-model solver in Notebook 05.

The earlier added derivations, independent helper-model comparisons, and
parameter adapter have been removed from these notebooks. The numerical
implementations now come from the original uploaded scripts. Small setup cells
only locate files, import libraries, choose bounded execution, and display plots.
See [SOURCE_MAP.md](SOURCE_MAP.md) and [source_provenance.json](source_provenance.json).

## ED code attribution

Notebook 05 retains **Glen Evenbly's
[Tensors.net exact-diagonalization example](https://www.tensors.net/exact-diagonalization)**
credit in `doApplyHam.py` and in the visible copied code. The effective-Hamiltonian
dense sector solver in Notebook 01 comes from `heff2_E_vs_N.py`.

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

Run each notebook top to bottom. Figures display inline. Source arrays and
original scripts are never overwritten. Notebook 01 has `RUN_FULL_ED=False` by
default; enabling it runs the full dense spectrum calculation and can be expensive.
Notebook 05 performs one four-site solve and loads the saved ten-site grid.

## Verify

```bash
python tools/execute_notebooks.py
```

The execution tool uses a fresh kernel for each notebook, validates its schema,
and stops on cell errors. [VALIDATION.md](VALIDATION.md) records this revision's
execution and source checks. The existing `anyon_thesis` package and tests remain
as earlier supplementary code, but are not imported by these notebooks or used
as the source of their scientific explanations.

Read [SCIENTIFIC_NOTES.md](../SCIENTIFIC_NOTES.md) for the unresolved source
frequency convention, the diagonal compensation operation, and the chain
script's phase/sign discrepancy. Those points are also labelled beside the
relevant notebook cells.
