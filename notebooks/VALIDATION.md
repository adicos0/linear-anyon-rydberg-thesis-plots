# Notebook validation record

Validated on 18 October 2026 against source commit
`22959963191fcbc4d283ec4b272b506ad6dd979b`.

## Executed notebooks

Each notebook was executed from a fresh Python kernel configured to fail on cell errors.
The final notebook files retain their numerical outputs and embedded figures.

| Notebook | Executed code cells | Embedded figures |
|---|---:|---:|
| 01: Effective Hamiltonian and spectrum | 5 | 1 |
| 02: Two-site populations | 5 | 1 |
| 03: Chain dynamics and phase tests | 5 | 1 |
| 04: Effective parameters and tunnelling | 7 | 4 |
| 05: Target anyon-Hubbard gaps | 6 | 2 |
| **Total** | **28** | **9** |

Notebook schemas pass `nbformat.validate`; code-cell execution counts are
consecutive and no cell contains an error output. All nine figures were
visually inspected for labels, clipping, legends, colour scales, and displayed
grid resolution. The additional parameter-figure families in Notebook 04 use
the explicitly stated restricted laser-amplitude ranges.

The full execution command is:

```bash
python tools/execute_notebooks.py
```

Notebook 04 was re-executed after its final figure additions. A specific
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
| Two-state analytic/numerical probability difference | 1.7e-15 |
| Largest chain total-occupation error | 3.2e-15 or smaller |
| Assigned gauge-matrix discrepancy | 1.2e-16 or smaller |
| Four-site teaching-sweep eigenpair residual | 6.9e-15 or smaller |

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
notebook outputs contain no absolute author-machine paths.

## What remains scientifically unresolved

The frequency convention, physical phase reachability, and self-consistent
compensation questions are documented in [SCIENTIFIC_NOTES.md](../SCIENTIFIC_NOTES.md).
Passing numerical checks does not settle those questions, prove convergence of
the entire archived ten-site spectrum, or establish a thermodynamic phase.
