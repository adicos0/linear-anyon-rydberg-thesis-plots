# Thesis scope and source discrepancies

The notebooks follow the theoretical model and numerical calculations in
Aditya Banerjee's [thesis](https://hdl.handle.net/2445/231684). They do not supply
new experimental evidence or new physical derivations. See the
[thesis-to-code map](notebooks/SOURCE_MAP.md) for sections, equations, and original
script links. The original scripts, saved data, and reference figures are unchanged.

## Effective spectrum and target-model spectrum

Section 3.7 (pp. 21–22) diagonalizes the effective Hamiltonian separately in each
particle-number sector and combines the blocks to obtain the full spectrum.
Notebook 01 explains that dense calculation using `heff2_E_vs_N.py`.

Appendix C (pp. 49–51) applies the target anyon-Hubbard Hamiltonian without
constructing a full matrix, then uses a sparse eigensolver for the lowest five
levels. Notebook 05 follows `sweepGapsMulti.py` and `doApplyHam.py`. The latter
retains the attribution to Glen Evenbly and Tensors.net. Appendix C explicitly
states that its gaps combine all particle-number sectors without a chemical
potential; fixed-filling gaps require a sector-restricted calculation.

## Frequency and time labels: unresolved source convention

The `heff2_states.py` header equates MHz with rad/us and uses hbar=1 in its
propagator. Appendix B of the thesis writes omega-prime/(2 pi)=40 MHz, while
Figure 3 and other passages label the numerical inputs in MHz. The notebooks
retain the original values and source labels for reproduction and identify this
discrepancy explicitly. They do not introduce a new dimensionless Hamiltonian
or a replacement time variable, and do not claim to have resolved the conversion.

## Phase comparison: theoretical target and physical restriction

Section 3.8.2 (pp. 27–28) discusses the hypothetical nonzero-phase target dynamics.
Section 4.6 and the conclusion restrict the analysed drive to the zero and pi
statistical angles. The original `heff2_spacetime.py` docstring claims arbitrary
angle realization; that claim is not used as the notebook's physics explanation.
Its function body is retained, with the conflicting docstring omitted.

A second discrepancy concerns signs: the source calculates `dens_phi` with
`mode='anyon'` at positive `ph`, but labels the original figure with a drive theta
and the relation phi=-theta. Notebook 03 labels those unchanged arrays with the
target phi passed to the function and flags the discrepancy. It does not claim
an exact reproduction of the thesis's signed-theta panels before that convention
is resolved by the author.

## Diagonal compensation: recorded code operation

The uploaded effective-Hamiltonian scripts set the 11 diagonal equal to the 20
diagonal. Notebook cells reproduce that operation. Comments describing retuning
are not evidence that the microscopic detunings and all denominators were
recomputed: the supplied code performs a diagonal shift. The parameter formulas
in Notebook 04 follow Section 4.2's addition-energy definitions and the original
parameter scripts.

## Reproduction scope

All notebook outputs are numerical calculations or plots of the archived data.
Notebook 01's optional full dense solve is disabled by default. Notebook 05's
live example is explicitly four sites; its seven gap maps use the saved ten-site
sweep. The stroboscopic interpretation and omitted micromotion are stated in
Section 3.8 and retained in the dynamics notebooks.

The previous assistant-created `anyon_thesis` package and its tests remain in
repository history and in the checkout, but the revised notebooks do not import
or use them. Their presence is not presented as original thesis code. The
notebook presentation was edited with Codex; the thesis supplies the physics
explanations, and the uploaded scripts supply the calculation implementations.
