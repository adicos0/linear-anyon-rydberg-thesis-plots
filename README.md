# Linear Anyon Model with Rydberg Atoms

**Figures, plotting code, and saved data from the master’s thesis by Aditya Banerjee**  
Master in Quantum Science and Technology · Universitat de Barcelona · 2026

**Read the thesis:** [University of Barcelona repository record](https://diposit.ub.edu/items/3c2114ee-4916-472a-a78a-e2c55571f981) · [permanent identifier](https://hdl.handle.net/2445/231684)

This repository is a companion to the thesis *Linear Anyon Model with Rydberg Atoms*. It brings together the code located for its numerical plots, the saved data needed to replot selected results, and the graphics used in the thesis. Start with the [figure guide](#find-a-figure) if you are looking for a particular plot.

## The question behind the work

The one-dimensional anyon-Hubbard model can be expressed in terms of bosons whose tunnelling acquires a phase that depends on how many particles occupy a site. This **occupation-dependent phase** represents the model’s statistical angle. The thesis asks whether laser-driven Rydberg atoms can produce that tunnelling rule.

The proposed setup encodes the local occupations |0⟩, |1⟩, and |2⟩ in atomic internal levels. Raman lasers and Rydberg interactions then produce effective tunnelling through off-resonant, second-order processes. The work calculates the resulting bond Hamiltonian, studies its spectrum and dynamics, and compares its tunnelling amplitudes and phases with those required by the anyon-Hubbard model.

**Result of the thesis, within the schemes analysed:** Static Raman driving reaches the pseudofermionic point, but its statistical angle is restricted to 0 or π. Adding a *single-frequency* periodic drive can change the tunnelling magnitudes, yet it does not remove that phase restriction. Thus these schemes do not realize an arbitrary statistical angle.

**Proposed direction, not a result demonstrated here:** A drive with several independently phased tones could address different occupation-dependent processes separately. The thesis identifies this as a route to investigate; this repository does not contain a completed multi-tone realization or experimental data.

## Find a figure

| Thesis figure | What it shows | Where to look |
|---|---|---|
| 1 | The active Rydberg bond and its panels | [Panel extraction script](scripts/make_active_bond_clean.py) · [reference PDF](reference_figures/rydberg_active_bond.pdf) |
| 2 | Effective spectrum versus particle number | [Plotting script](scripts/heff2_E_vs_N_thesis.py) · [spectrum data](data/heff2_E_vs_N_L10.npz) |
| 3 | Population dynamics of the driven bond | [Script](scripts/heff2_states.py) · [thesis graphic](reference_figures/floquet_population_dynamics.png) |
| 4–5 | Spacetime dynamics when the bond model is extended to a chain | [Script](scripts/heff2_spacetime.py) · [graphic 4](reference_figures/heff2_spacetime.png) · [graphic 5](reference_figures/heff2_spacetime_phi.png) |
| 6–7 | Effective amplitudes, chemical potential, and interaction | [Amplitude script](scripts/amplitudes_vs_muU_thesis.py) · [parameter script](scripts/mu_vs_A_over_U.py) |
| 8–9 | Energy gaps of the **target anyon-Hubbard model** across its parameter grid | [Plotting script](scripts/replotGapsMulti.py) · [saved grid](data/anyonHubbardGapsMulti_N10.npz) |
| 10–17 | Transition and level schematics | Graphics are in [reference_figures/](reference_figures/); their generator or editable source was not located. |

The Figures 8–9 gap maps describe the target model. They should not be read as evidence that the Rydberg scheme realizes arbitrary statistical angles.

## What is in each folder

- [`scripts/`](scripts/) contains 10 original Python files. Some draw figures directly; others calculate the spectra used by a plotting script. In particular, [`heff2_E_vs_N.py`](scripts/heff2_E_vs_N.py) produces the Figure 2 spectrum, while [`sweepGapsMulti.py`](scripts/sweepGapsMulti.py) and [`doApplyHam.py`](scripts/doApplyHam.py) produce the data behind Figures 8–9.
- [`data/`](data/) contains two saved `.npz` files, allowing those figures to be replotted without repeating the longer numerical calculations.
- [`reference_figures/`](reference_figures/) contains the graphics from the thesis project, including schematics for which no generator was found.
- [`SOURCE_MANIFEST.csv`](SOURCE_MANIFEST.csv) records each file’s original location, role, size, and SHA-256 checksum. Its original paths refer to the computer on which the files were found; they are provenance records, not paths a new user needs to recreate.

## Replot selected figures

The scripts use NumPy, SciPy, and Matplotlib. The Figure 1 panel extraction script also imports PyMuPDF (`fitz`). No original dependency lockfile was found, so exact rendering can vary with package versions.

From the repository root, these commands replot Figure 2 using the saved spectrum and Figures 8–9 using the saved energy grid:

```bash
mkdir -p figures
cp data/heff2_E_vs_N_L10.npz .
python scripts/heff2_E_vs_N_thesis.py
python -c "import sys; sys.path.insert(0, 'scripts'); from replotGapsMulti import plot_all; plot_all(npzfile='data/anyonHubbardGapsMulti_N10.npz', symbol='phi', outdir='figures')"
```

The commands write new images to `figures/`; the archived images remain in `reference_figures/`. Other scripts retain their original relative output paths, so check the script before running it. Recomputing the Figure 2 spectrum or the Figures 8–9 energy grid can take substantial time and memory. The Figures 4–5 script also performs a large matrix calculation.

## What was checked

- All copied Python files passed syntax parsing.
- The 13 raster graphics used across Figures 2–9 have identical RGB pixels to the images extracted from the supplied 62-page thesis PDF.
- Figure 2 and all seven gap panels in Figures 8–9 were regenerated from the included saved data with identical RGB pixels to the thesis project graphics.
- The existing output of the Figure 3 script has identical RGB pixels to the thesis graphic, although a fresh render with the current software environment was not pixel identical.
- The Figure 1 script extracts panels from an existing PDF; the source that originally drew that PDF was not located. Likewise, the generator or editable source for the Figures 10–17 schematics was not located.

These checks establish which files produced or correspond to the thesis figures. They do not independently validate the physical model or its numerical assumptions. The scripts and reference graphics were copied without altering their contents; the [source manifest](SOURCE_MANIFEST.csv) preserves their checksums.
