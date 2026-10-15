# Thesis plot code collection

This folder collects the original code files I located for the plots in **Masters_Thesis_MQST_Final__Linear_Anyon_Model_with_Rydberg_Atoms__.pdf** (62 pages; SHA-256 `5753cc5932becab0303a5c3e2fd54eb7b5b2905580211479e605e25d0497b27e`). `SOURCE_MANIFEST.csv` records each original location and SHA-256 checksum. All originals were copied, not moved or edited. The `reference_figures/` folder holds the graphics from the thesis LaTeX project, so you can identify a script's output by eye.

## What each figure uses

| Thesis figure | Graphic in `reference_figures/` | Located code | Evidence and qualification |
|---|---|---|---|
| 1 | `rydberg_active_bond_A/B/C.pdf` | `scripts/make_active_bond_clean.py` | The script extracts/crops panels from `rydberg_active_bond.pdf`. The original drawing authoring source was not located. |
| 2 | `heff2_spectrum_vs_N.png` | `scripts/heff2_E_vs_N_thesis.py` | The script explicitly writes this filename. `heff2_E_vs_N.py` creates its cached spectrum. |
| 3 | `floquet_population_dynamics.png` | `scripts/heff2_states.py` | The script writes `heff2_states.png`; the existing rendered file has **identical RGB pixels** to the thesis graphic, 2400 × 2025 pixels. |
| 4–5 | `heff2_spacetime.png`, `heff2_spacetime_phi.png` | `scripts/heff2_spacetime.py` | The AC2 revision is the version whose figure graphics were packaged. Its code writes both filenames. |
| 6 | `amplitudes_vs_muU.png` | `scripts/amplitudes_vs_muU_thesis.py` | The script explicitly writes this filename in `figures/`. |
| 7 | `mu_vs_A_over_U.png` | `scripts/mu_vs_A_over_U.py` | The script writes this filename; the thesis source copies it into `figures/`. |
| 8–9 | seven `anyonHubbardGap_N10_*.png` files | `scripts/replotGapsMulti.py` | The renderer names all seven graphics. It reads the cached grid created by `sweepGapsMulti.py` and `doApplyHam.py`. Use `symbol='phi'` to match the thesis axis label. |
| 10–17 | `diamond_*`, `levels_*`, `photon_*` | **Unverified / needs check** | These are schematics rather than numerical plots. I found their image files but no code or editable source that generated them. |

## Inputs and software

- `data/heff2_E_vs_N_L10.npz`: spectrum cache for Figure 2.
- `data/anyonHubbardGapsMulti_N10.npz`: energy grid for Figures 8–9.
- Python packages used across these scripts: NumPy, SciPy, Matplotlib; `make_active_bond_clean.py` additionally uses PyMuPDF (`fitz`). These are imports seen in the located code; no original environment lockfile was found.

## Replot using the saved spectra

From this folder, using a Python environment with the packages above:

```bash
cp data/heff2_E_vs_N_L10.npz .
mkdir -p figures
python scripts/heff2_E_vs_N_thesis.py
python scripts/heff2_states.py
python scripts/amplitudes_vs_muU_thesis.py
python scripts/mu_vs_A_over_U.py
python -c "import sys; sys.path.insert(0, 'scripts'); from replotGapsMulti import plot_all; plot_all(npzfile='data/anyonHubbardGapsMulti_N10.npz', symbol='phi', outdir='figures')"
```

To regenerate Figures 4–5, run `python scripts/heff2_spacetime.py` separately. This script also performs a large matrix based gauge check, so allow substantial memory and time. The scripts create new files in this folder (`figures/`, and for some scripts the folder root). The `reference_figures/` files preserve the original graphics. Figure 2's full-spectrum generator and Figures 8–9's Lanczos sweep are provided for provenance; rerunning them can take substantial time and memory. The numerical content and model assumptions are those of the original scripts, which I have not independently validated against the underlying physics.

## Verification

- Every copied Python file passed Python syntax parsing.
- I extracted the raster plots from the **supplied PDF** and compared their RGB pixels to the archived reference images. All 13 image files used across Figures 2–9 matched exactly. This establishes that these graphics are the ones in the supplied PDF.
- Figure 2 regenerated from the copied script and saved spectrum with identical RGB pixels to the thesis project graphic.
- All seven panels used by Figures 8–9 regenerated from the saved energy grid with identical RGB pixels when the phase symbol was set to `phi`.
- The already existing `heff2_states.png` has identical RGB pixels to Figure 3's `floquet_population_dynamics.png`; PNG file bytes differ because of encoding. A fresh render in the current software environment was not pixel identical, so use the archived reference graphic for exact appearance.
- The existing output of `mu_vs_A_over_U.py` is byte identical to the thesis project copy used for Figure 7.
- Figures 4–6 were matched through their original output names, source locations, and project graphics, with the graphics also matched to the PDF as above. I did not complete new renders of those computation and rendering intensive scripts.

The render and source matches establish file provenance. They do not validate the physical model or numerical assumptions in the code.

## Search scope

I searched the accessible user directory, focusing on Documents, Downloads, Desktop, the thesis project and its copied Overleaf/revision folders. I matched image names to the LaTeX source and code output paths, compared duplicate copies, and checked Figure 3 by an exact RGB pixel comparison. A missing schematic source means only that it was not located in this search; it may exist on an external drive or in an application project that is not present here.
