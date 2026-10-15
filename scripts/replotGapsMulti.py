"""
replotGapsMulti.py
---------------------------------------------------------------------
Plot heat maps of energy differences of the 1D anyon-Hubbard model
from the saved multi-level sweep data, without re-running the Lanczos
sweep.

One standalone figure is written per energy difference, with J/U on
the x axis and the statistical phase (0 to pi) on the y axis:
  anyonHubbardGap_N<L>_E<k>-E0.png     for k = 1..nlev-1
  anyonHubbardGap_N<L>_E<k+1>-E<k>.png for k = 1..nlev-2

The symbol used for the statistical phase on the y axis is settable
("theta" by default, "phi" to match the thesis notation), as is the
output directory.
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SYMBOLS = {'theta': r'\theta', 'phi': r'\phi'}


def _plot_one(ratios, thetas, diff, label, fname, Nsites, bc, sym):
  """Render a single gap heat map to its own figure file."""
  fig, ax = plt.subplots(figsize=(7, 5.2))
  mesh = ax.pcolormesh(ratios, thetas, diff.T, cmap='inferno',
                       shading='gouraud')
  cbar = fig.colorbar(mesh, ax=ax)
  cbar.set_label(r'$(%s)/U$' % label)
  ax.set_xlabel(r'$J/U$')
  ax.set_ylabel(r'Statistical phase $%s$' % sym)
  ax.set_yticks(np.pi * np.array([0, 0.25, 0.5, 0.75, 1.0]))
  ax.set_yticklabels([r'$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$'])
  ax.set_title(r'$%s$, $N = %d$ sites (%s)' % (label, Nsites, bc))
  fig.tight_layout()
  fig.savefig(fname, dpi=200)
  plt.close(fig)
  print('Saved: %s' % fname)


def plot_all(npzfile=None, Nsites=None, symbol='theta', outdir='.'):
  if npzfile is None:
    npzfile = 'anyonHubbardGapsMulti_N%d.npz' % Nsites
  sym = SYMBOLS[symbol]
  dat = np.load(npzfile)
  thetas, ratios = dat['thetas'], dat['ratios']
  energies = dat['energies']            # shape (n_ratios, n_thetas, nlev)
  L, usePBC = int(dat['Nsites']), bool(dat['usePBC'])
  nlev = energies.shape[2]
  bc = 'PBC' if usePBC else 'OBC'
  os.makedirs(outdir, exist_ok=True)

  def path(name):
    return os.path.join(outdir, name)

  # Gaps measured from the ground state: E_k - E_0
  for k in range(1, nlev):
    _plot_one(ratios, thetas, energies[:, :, k] - energies[:, :, 0],
              'E_%d - E_0' % k,
              path('anyonHubbardGap_N%d_E%d-E0.png' % (L, k)), L, bc, sym)

  # Consecutive level spacings: E_{k+1} - E_k
  for k in range(1, nlev - 1):
    _plot_one(ratios, thetas, energies[:, :, k + 1] - energies[:, :, k],
              'E_%d - E_%d' % (k + 1, k),
              path('anyonHubbardGap_N%d_E%d-E%d.png' % (L, k + 1, k)),
              L, bc, sym)


if __name__ == '__main__':
  plot_all(Nsites=10)
