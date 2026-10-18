# Validation of the thesis-and-original-code revision

Validated on 18 October 2026. The notebooks were executed in separate fresh
Python kernels using `python tools/execute_notebooks.py`. All default cells
completed without cell errors. The stored notebook outputs come from this run.

| Notebook | Executed code cells | Inline figures |
|---|---:|---:|
| 01 — Effective Hamiltonian and exact diagonalization | 13 | 1 four-panel figure |
| 02 — Two-site populations | 9 | 1 six-panel figure |
| 03 — Chain occupations | 11 | 2 figures |
| 04 — Effective parameters | 10 | 2 figures |
| 05 — Target-model gaps | 11 | 7 figures |
| Total | 54 | 13 figures |

## Source checks

- All 42 source-code excerpts were verified against their recorded original
  script line ranges and SHA-256 hashes.
- All 24 copied calculation-function occurrences have the same Python syntax
  trees as the originals, excluding the deliberately removed conflicting
  `phases_of` docstring. The two adapted gap-plotting functions are recorded
  separately in the provenance metadata.
- Every code cell has either original-source provenance or a labelled notebook
  setup/execution role. No notebook imports the assistant-created `anyon_thesis`
  implementation.
- All 42 original archived files other than the intentionally revised root README
  match commit `22959963191fcbc4d283ec4b272b506ad6dd979b` byte for byte. This includes
  the original scripts, saved arrays, reference figures, and source manifest.
- All five notebooks pass schema validation and code syntax parsing. No stored
  output is a cell error. The original spectrum checks print PASS.
- The live N=2 spectrum from the original effective-Hamiltonian function agrees
  with the saved N=2 spectrum to a maximum absolute difference of
  `8.881784197001252e-16` in the source numerical energy units.
- The equation transcriptions were checked against the thesis text. The
  population equations (49), (50), (54), and (55) were also checked visually
  against printed pages 24–25 of the PDF.
- All 13 executed figures were visually reviewed for their content and labels.

## Scope and limitations

This validates execution and correspondence to the supplied sources. It does
not resolve the source discrepancies documented in
[SCIENTIFIC_NOTES.md](../SCIENTIFIC_NOTES.md), and does not provide independent
experimental validation.

The full 21-block dense spectrum and the 441-point ten-site Lanczos sweep were
not recomputed. Their original archived data are used. The optional
`RUN_FULL_ED=True` branch is syntax-checked but was not executed in full.
Notebook 01 runs the original 55-state N=2 solver; Notebook 05 runs the original
sparse solver at a labelled four-site point. No claim of pixel-identical
reproduction of thesis figures is made for this notebook revision.

The existing helper-model test suite belongs to the earlier supplementary
implementation. It is not used as evidence of source fidelity for this revision.

## Environment

Python 3.12 with the environment recorded in
[`requirements-notebooks-tested.txt`](../requirements-notebooks-tested.txt).
Use [`requirements-notebooks.txt`](../requirements-notebooks.txt) for a portable
installation. Jupyter needs permission to open local loopback kernel sockets.
The host emitted a macOS process-inspection permission diagnostic during kernel
shutdown; all notebook cells completed and the executor exited successfully.
