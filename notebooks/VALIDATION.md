# Notebook validation record

Validated on 18 October 2026 against source commit
`22959963191fcbc4d283ec4b272b506ad6dd979b`.
The readable revision starts from repository commit
`f4454a4ecddcd6a94958a941d33d021a2fb13e81`.

## Reading order and source fidelity

The main sections now follow a short thesis explanation, the relevant equation,
visible source-derived code, and its output. Longer symbolic derivations and
implementation comparisons are in the appendix. Each core calculation is in the
notebook; only the final verification appendix imports `anyon_thesis`.

Code-cell provenance metadata identifies the original script, source commit,
and any adaptation. Original function bodies are retained or split into cells.
Plotting, path handling, sample counts, and the explicitly declared live-system
size are adapted for interactive execution. The tensor-contraction routine keeps
the original Glen Evenbly / tensors.net attribution.

| Notebook | Explicit notebook adaptations |
|---|---|
| 01 | Split source sideband, diagonal, and hopping calculations into separate cells; solve only N=2 live and load the saved full spectrum; show two principal spectrum panels. |
| 02 | Show the original nine-state Hamiltonian, evolve/ket functions, and all six source initial states; 1000 time samples per trace instead of 5000. |
| 03 | Use the source ten-site, 900-sample, two-particle example; show single versus double occupation and three assigned phase patterns; correct the phase-function docstring without changing its branches. |
| 04 | Show the source Sigma, Sigma4, muU_x, and curves functions; use named output arrays; keep explicitly bounded amplitude-family ranges and separate branches at U_eff=0. |
| 05 | Show the original tensor contraction and LinearOperator/eigsh solve; solve one four-site example rather than the full multiprocessing sweep; display all seven saved ten-site maps separately. |

Effective-model time axes use the declared dimensionless convention. Heatmaps
use nearest-cell shading rather than the source Gouraud interpolation. Figure
export is opt-in. These presentation changes do not alter the saved arrays.

## Executed notebooks

Each notebook was executed from a fresh Python kernel configured to fail on cell errors.
The final notebook files retain their numerical outputs and embedded figures.

| Notebook | Executed code cells | Embedded figures |
|---|---:|---:|
| 01: Effective Hamiltonian and spectrum | 13 | 1 |
| 02: Two-site populations | 14 | 1 |
| 03: Chain dynamics and phase tests | 16 | 2 |
| 04: Effective parameters and tunnelling | 13 | 4 |
| 05: Target anyon-Hubbard gaps | 15 | 7 |
| **Total** | **71** | **15** |

Notebook schemas pass `nbformat.validate`; code-cell execution counts are
consecutive and no cell contains an error output. All fifteen figures were
visually inspected for labels, clipping, legends, colour scales, and displayed
grid resolution. The additional parameter-figure families in Notebook 04 use
the explicitly stated restricted laser-amplitude ranges.

The full execution command is:

```bash
python tools/execute_notebooks.py
```

Notebook 04 was re-executed after its final axis-label cleanup. A specific
notebook can also be selected by filename:

```bash
python tools/execute_notebooks.py 04_effective_parameters_and_tunnelling.ipynb
```

## Numerical checks

`python -m pytest -q` reports **16 passed**. The checks include:

- Independent equivalence to the original effective-rate functions and analytic
  parameter functions, including three modulation indices and several drive rates.
- Recalculation of the saved ten-site effective spectrum in the N=0,1,2 sectors.
- Exact decomposition by particle number and conservation during evolution.
- Two-state analytic dynamics, zero-hopping limits, and sideband convergence.
- The stated gauge identity for the two assigned finite Hamiltonians.
- Comparison of the target Hamiltonian with the original tensor-contraction
  implementation for open and periodic boundaries.
- Sparse eigenpair residuals compared with a dense reference, sector selection,
  phase periodicity, and cached-array shape and ordering.

Representative notebook outputs in the tested environment:

| Check | Observed value |
|---|---:|
| Two-state analytic/numerical probability difference | 1.73e-15 or smaller |
| Largest ten-site chain total-occupation error | 6.00e-15 or smaller |

Every notebook's appendix additionally checks its visible source calculation
against the validated reference. Matrix elements are compared at absolute
tolerance 1e-13; the saved N=2 sector at 2e-12; parameter formulas at 1e-11;
and the original matrix-free low-energy solve against dense diagonalization at
1e-8. These are comparison tolerances, not measured error bars.

These are numerical checks of finite matrices, not experimental uncertainties.
Small changes in floating-point results across supported platforms are expected.

## Preservation and environment

Git blob hashes confirm that all **42 original tracked files other than the
edited root README** remain byte-identical. This includes the original scripts,
saved data, source manifest, and reference graphics. No original large solve
was rerun and no source array was overwritten. `git diff --check` passes.

The execution environment used Python 3.12, NumPy 2.5.4, SciPy 1.18.1,
Matplotlib 3.11.2, nbformat 5.11.1, nbclient 0.11.0, and pytest 9.1.1.
The full package record is in `requirements-notebooks-tested.txt`; use the
portable `requirements-notebooks.txt` for installation.

The host sandbox emitted an ipykernel/psutil permission diagnostic while
enumerating child processes during kernel shutdown. Every notebook execution
returned successfully; this diagnostic occurred outside the notebook cells.
Runtime caches were isolated in a temporary directory, and the committed
notebook outputs contain no absolute author-machine paths. A font-cache-building
informational message was omitted from the display outputs; numerical outputs
and errors were not suppressed. Code formatting was checked to preserve the
Python abstract syntax trees after removing the Jupyter plotting magic.

## What remains scientifically unresolved

The frequency convention, physical phase reachability, and self-consistent
compensation questions are documented in [SCIENTIFIC_NOTES.md](../SCIENTIFIC_NOTES.md).
Passing numerical checks does not settle those questions, prove convergence of
the entire archived ten-site spectrum, or establish a thermodynamic phase.
