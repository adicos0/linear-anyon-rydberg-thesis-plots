# Thesis and original-code map

The five notebooks use the explanations and equations in Aditya Banerjee's
[*Linear Anyon Model with Rydberg Atoms*](https://hdl.handle.net/2445/231684).
References use printed thesis page numbers. Physics paragraphs are labelled
paraphrases; equation numbers are the thesis's own numbers. No new derivation,
reference-model implementation, or external physics interpretation is added.

Original code snapshot:
[`22959963191fcbc4d283ec4b272b506ad6dd979b`](https://github.com/adicos0/linear-anyon-rydberg-thesis-plots/tree/22959963191fcbc4d283ec4b272b506ad6dd979b/scripts).
Each source code cell has a visible link to its original line range. Its metadata
records that range, the source commit, a SHA-256 of the original excerpt, and the
notebook changes. [source_provenance.json](source_provenance.json) collects these
records. Setup cells are explicitly labelled and do not implement another model.

| Notebook | Thesis basis | Original scripts |
|---|---|---|
| 01 — Effective Hamiltonian and exact diagonalization | Sections 3.6–3.7, pp. 20–22; Eqs. (41)–(47); Figure 2 | `heff2_E_vs_N.py` |
| 02 — Two-site population dynamics | Sections 3.6 and 3.8–3.8.1, pp. 20–25; Eqs. (48)–(55); Figure 3 | `heff2_states.py` |
| 03 — Chain occupation dynamics | Section 3.8.2, pp. 25–28; Eq. (56); Figures 4–5; Section 4.6 phase restriction | `heff2_spacetime.py` |
| 04 — Effective parameters and tunnelling | Sections 4.1–4.4, pp. 29–33; Eqs. (65), (69)–(71), (75); Figures 6–7 | `amplitudes_vs_muU_thesis.py`, `mu_vs_A_over_U.py` |
| 05 — Target model and gaps | Sections 1.2–1.3, pp. 7–8; Eqs. (6)–(7); Appendix C, pp. 49–51; Eq. (100); Figures 8–9 | `sweepGapsMulti.py`, `doApplyHam.py`, `replotGapsMulti.py` |

## Execution adaptations

- All plots display inline instead of exporting images. Original data and scripts
  are not overwritten. No new absolute machine paths are required.
- Notebook 01 calls the original solver for the N=2 block and loads the complete
  saved spectrum. Its optional full-solve branch copies the original block loop
  from `heff2_E_vs_N.py`, lines 182–192 (the source function and loop can also be
  located under `FULL SPECTRUM, SECTOR BY SECTOR`). The branch omits cache writing.
  The original four-panel source plot includes two supplementary panels beyond
  the two-panel thesis Figure 2; the notebook identifies this distinction.
- Notebook 02 retains the original six initial states and 5000 time samples per
  panel. Added analytic-test implementations have been removed; the analytic
  explanation uses the equations in Section 3.8.1.
- Notebook 03 retains the original ten-site dynamics and six target-phase arrays.
  It removes an unsupported physical-realization claim from the `phases_of`
  docstring and identifies the original script's disagreement with the thesis.
  The plotted symbol is phi, matching the positive target phase actually passed
  to `mode='anyon'`. No matrix elements are changed to force agreement with the
  thesis's signed-theta description of Figure 5.
- Notebook 04 uses the original parameter sweep ranges and calculation functions.
  The previously added parameter adapter and extra teaching curves are removed.
- Notebook 05 runs the original solver at one labelled four-site point and loads
  the original ten-site sweep for the gap maps. The source gap plotting function
  uses inline display, nearest-cell shading, and L for the number of sites.
  The original eigenvalues and gap differences are unchanged.

## Attribution

The matrix-free Hamiltonian-application routine in Notebook 05 comes from the
uploaded `doApplyHam.py`, which credits **Glen Evenbly's
[Tensors.net exact-diagonalization example](https://www.tensors.net/exact-diagonalization)**.
That attribution remains in both the original script and the notebook code.
The dense sector solver for the effective Hamiltonian is the repository's
`heff2_E_vs_N.py` implementation.
