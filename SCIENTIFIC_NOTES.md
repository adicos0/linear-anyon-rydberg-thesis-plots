# Scientific conventions, evidence, and unresolved interpretations

The notebooks explain numerical models and source results. They do not add
experimental evidence or replace the thesis's physical derivation. The original
scripts, data, reference graphics and source manifest are preserved byte for byte.

## Source record

- Aditya Banerjee, *Linear Anyon Model with Rydberg Atoms*, Universitat de
  Barcelona, 2026: [institutional record](https://hdl.handle.net/2445/231684).
- [Zenodo thesis and code archive](https://doi.org/10.5281/zenodo.23443035).
- Source code commit: `22959963191fcbc4d283ec4b272b506ad6dd979b`.

Notebook equations are transcriptions or explanatory derivations of this source
model. References below are verified links, not generated citation identifiers.
The notebooks identify the relevant thesis sections. No new literature claim is
inferred from an unverified bibliography entry.

## 1. Ordinary frequency and angular frequency

**Verified source inconsistency:** `heff2_states.py` says "energies in MHz
(= rad/us), time in us" and uses `exp(-1j * E * t)` with hbar numerically one.
Appendix B of the thesis (printed page 45) writes `omega'/2 pi = 40 MHz`, while
Table 1 uses an omega-prime column labelled MHz. MHz describes cycles per second;
an angular frequency differs by a factor of two pi. Dimensions alone cannot
resolve the convention because both have inverse-time dimensions.

**Notebook convention:** use dimensionless `K = H/(hbar*omega_ref)` and
`tau = omega_ref*t`, with `exp(-1j * K * tau)`. Original effective-model numerical
inputs are retained in a declared native rate scale, with a reference equal to
one native rate unit. Saved spectra are labelled in their legacy energy unit.
No default plot asserts that its time axis is in microseconds.

If physical inputs are specified as `H/h` in cycles/us, the propagator requires
`exp(-1j * 2*pi * (H/h) * t_us)`. Inputs specified as `H/hbar` in rad/us require
`exp(-1j * (H/hbar) * t_us)`. A future experimental calibration must select and
propagate one convention consistently through every rate and detuning.

## 2. Assigned phases versus physically achievable controls

**Verified source disagreement:** the thesis's conclusion (printed page 37)
states that the analysed static and single-tone schemes restrict the statistical
angle to zero or pi. The function `phases_of` and related comments in
`heff2_spacetime.py` call a uniform-bond assignment a physical drive that realizes
arbitrary angle, based on a diagonal gauge identity.

**What is checked:** for the two matrices assembled with the assigned patterns,
the identity `Q^dag K_legacy(theta) Q = K_target(-theta)` holds, with
`Q(s) = exp(2j*theta*number_of_doubly_occupied_sites(s))`. This is a statement
about those finite Hamiltonians with identical diagonal terms and magnitudes.
Their phase patterns alone do not establish equality with the complete AHM:
the four bosonic magnitude factors and the diagonal terms must also match.

**What remains unresolved:** whether a consistent set of the stated site laser
controls produces all of those bond phases while satisfying the microscopic
derivation's assumptions. The notebooks call these patterns `legacy_drive`
(assumed) and `target`, and do not claim physical arbitrary-angle realization.
The original files are retained as historical evidence, rather than silently
rewritten to settle this scientific question.

## 3. Imposed diagonal compensation

The original effective-Hamiltonian scripts set the numerical `|11>` diagonal
equal to the `|20>` diagonal. The new `match_11_diagonal` switch reproduces that
operation and reports the shift. It does not self-consistently recompute every
denominator after changing a laser detuning. Notebook 4's analytic effective
chemical potential and interaction formula family is evaluated before this
additional compensation; do not identify it with the entire compensated chain.

## 4. Global level spacings, sector excitations, and finite size

`sweepGapsMulti.py` works in the full `3**L` truncated Hilbert space. Its saved
five lowest eigenvalues do not impose total occupation. A spacing between its
lowest levels can connect different particle numbers. The notebooks distinguish
this spacing from an excitation at fixed particle number and from the charge
gap. They do not infer a Mott phase or a thermodynamic gap from these small
finite-size calculations alone.

The local occupation cutoff is zero, one, or two. In this space the bosonic
commutator is `diag(1,1,-2)`, not the infinite-space identity. Open boundaries
are the default. Periodic examples use at least three sites so the convention
does not silently double count the single bond of a two-site chain.

## 5. Plotting, numerical validation, and remaining limits

- Cached arrays are loaded without alteration or pickle deserialization.
- Gap heatmaps show grid cells with nearest-cell shading; no Gouraud smoothing
  or interpolated intermediate results are implied.
- Occupation maps share a zero-to-two colour scale. Line plots use line style
  and labels as well as colour. Numerical values or conservation errors are
  printed beside the figures.
- No error bars are invented for deterministic numerical calculations. Solver
  residuals, sideband convergence, and analytic comparisons are reported instead.
- Tests compare the new target matrix with the original tensor contraction for
  open and periodic chains, and recompute the cached effective spectrum's small
  number sectors. These checks validate implementations, not the experiment.
- Default notebook runs are bounded. The original ten-site spectrum calculation
  and full 441-point sweep are not launched on import or by Run All.
- Mathematical statements about a formal large-drive limit do not certify the
  dispersive approximation there.

## Software and workflow references

- [NumPy linear algebra](https://numpy.org/doc/stable/reference/routines.linalg.html).
- [SciPy eigsh](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html).
- [Jupyter nbclient execution](https://nbclient.readthedocs.io/en/latest/client.html).
- [Matplotlib colormap guidance](https://matplotlib.org/stable/users/explain/colors/colormaps.html).
- Timothy Kassis, Vinayak Agarwal, Yuhuan He, Darshil Patel, and Aubrey M.
  Brueckner (2026), *Scientific Agent Skills: A Library of Procedural Knowledge
  for Research Agents*, [arXiv:2609.00065](https://arxiv.org/abs/2609.00065),
  [DOI](https://doi.org/10.48550/arXiv.2609.00065). Bibliographic metadata checked
  against the arXiv record on 18 October 2026. This is a workflow-tool reference,
  not supporting evidence for the physics.

The new notebook prose and numerical refactoring were prepared with Codex.
An isolated environment executed the tests and all default notebook cells.
Scientific interpretation and physical calibration remain subject to author
review. The repository contains no newly inferred licence grant or authorship.
