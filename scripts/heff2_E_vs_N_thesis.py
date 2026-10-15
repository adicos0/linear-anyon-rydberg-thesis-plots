"""
heff2_E_vs_N_thesis.py
---------------------------------------------------------------------
Thesis figure: the spectrum of H_eff^(2) against the particle-number
operator, in the extensive and intensive views.

    figures/heff2_spectrum_vs_N.png   (two panels, written into figures/)

This is the (a) + (c) pair of heff2_E_vs_N.py, recoloured blue and sized for
a full-textwidth figure.  It reuses the cached spectrum written by that
script (heff2_E_vs_N_L10.npz), so it costs nothing to re-render; delete the
cache or run heff2_E_vs_N.py to recompute.

Left panel   E versus the eigenvalue N of N-hat.  [H, N-hat] = 0, so every
             one of the 3^L eigenstates carries a definite integer N and the
             pair (E_m, N_m) is unambiguous.
Right panel  the intensive version, E/L against the average number density
             n-bar = N/L, i.e. the eigenvalues of N-hat/L.  Both axes are
             divided by L, which is the combination with a finite large-L
             limit; the x-axis alone would be the left panel rescaled.

Units: hbar = 1, energies in MHz.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CACHE = 'heff2_E_vs_N_L10.npz'
OUT = 'figures/heff2_spectrum_vs_N.png'   # distinct name: heff2_E_vs_N.png in the
                                       # repo root is the 4-panel diagnostic
                                       # version and would shadow this one

z = np.load(CACHE)
L = int(z['L'])
WP = float(z['wp'])
THETA = float(z['theta'])
E, NLAB = z['E'], z['NLAB']

Ns = np.arange(2 * L + 1)
E0N = np.array([E[NLAB == N].min() for N in Ns])
EmaxN = np.array([E[NLAB == N].max() for N in Ns])
igs = int(np.argmin(E))

# blue palette
C_CLOUD = '#7fb3d5'
C_MIN = '#12395e'
C_MAX = '#2e86c1'

print('L = %d, dim = 3^%d = %d, omega\' = %g MHz, theta = %.3f'
      % (L, L, E.size, WP, THETA))
print('ground state  E_0 = %+.6f MHz at N = %d  (n-bar = %.2f)'
      % (E[igs], NLAB[igs], NLAB[igs] / L))
print('E_0(N) flat region: N = 0..4 ->', ' '.join('%+.4f' % x for x in E0N[:5]))
print('bandwidth = %.6f MHz' % (E.max() - E.min()))

fig, axs = plt.subplots(1, 2, figsize=(13.2, 5.0))
rng = np.random.default_rng(0)

# ---- (a) extensive:  E vs N ----------------------------------------
ax = axs[0]
ax.plot(NLAB + 0.32 * (rng.random(NLAB.size) - 0.5), E, '.', ms=0.9,
        color=C_CLOUD, alpha=0.30, rasterized=True)
ax.plot(Ns, E0N, 'o-', color=C_MIN, ms=4.5, lw=1.6, label=r'$E_0(N)$')
ax.plot(Ns, EmaxN, 's--', color=C_MAX, ms=3.8, lw=1.3,
        label=r'$E_{\max}(N)$')
ax.plot([NLAB[igs]], [E[igs]], '*', color='k', ms=15,
        label='global ground state')
ax.set_xlabel(r'$\hat N$ eigenvalue $N=\sum_j n_j$')
ax.set_ylabel(r'$E$  [MHz]')
ax.set_title(r'(a) extensive: all $3^{%d}$ levels vs $\hat N$' % L, fontsize=11)
ax.set_xticks(Ns[::2])
ax.legend(fontsize=8, loc='upper left')
ax.grid(True, ls='--', alpha=0.35)

# ---- (b) intensive:  E/L vs n-bar ----------------------------------
ax = axs[1]
ax.plot(NLAB / L + (0.32 / L) * (rng.random(NLAB.size) - 0.5), E / L, '.',
        ms=0.9, color=C_CLOUD, alpha=0.30, rasterized=True)
ax.plot(Ns / L, E0N / L, 'o-', color=C_MIN, ms=4.5, lw=1.6,
        label=r'$E_0(N)/L$')
ax.plot(Ns / L, EmaxN / L, 's--', color=C_MAX, ms=3.8, lw=1.3,
        label=r'$E_{\max}(N)/L$')
ax.plot([NLAB[igs] / L], [E[igs] / L], '*', color='k', ms=15)
ax.axvline(1.0, color='k', ls=':', lw=1.0, label='unit filling')
ax.set_xlabel(r'average number density $\bar n=\hat N/L$'
              '\n' r'(spacing $1/L=%.2f$)' % (1.0 / L))
ax.set_ylabel(r'$E/L$  [MHz per site]')
ax.set_title(r'(b) intensive: $E/L$ vs $\bar n$', fontsize=11)
ax.legend(fontsize=8, loc='upper left')
ax.grid(True, ls='--', alpha=0.35)

fig.tight_layout()
fig.savefig(OUT, dpi=200)
print('Saved -> %s' % OUT)
