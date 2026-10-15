"""
heff2_E_vs_N.py
---------------------------------------------------------------------
Energy eigenvalues of H_eff^(2) plotted against the eigenvalues of the total
particle-number operator N-hat, and against the average number density
(particles per site).

    heff2_E_vs_N.png

WHY THIS IS WELL DEFINED
------------------------
Every term of H_eff^(2) moves exactly one quantum between neighbours, so
[H, N-hat] = 0.  H is therefore block diagonal in N and EVERY eigenstate
carries a definite integer N: the pair (E_m, N_m) is unambiguous, and no
state has to be assigned a fractional or averaged particle number.  The
full spectrum is obtained exactly by diagonalizing each of the 2L+1 blocks
with dense LAPACK and concatenating.

THE TWO REQUESTED X-AXES ARE THE SAME AXIS, RESCALED
----------------------------------------------------
The average number density operator is

    n-bar = N-hat / L ,

so its eigenvalues are just N/L: panel (c) is panel (a) with the x-axis
divided by L.  To make the density view carry information the extensive plot
does not, panel (c) also divides the ENERGY by L, giving the intensive
(per-site) picture E/L versus n-bar - the combination that has a finite
thermodynamic limit.  Panel (a) is the extensive one.

SPEED: at theta = 0 all four channel phases vanish, H is REAL SYMMETRIC, and
LAPACK's real path is used (~5 min for the whole 3^10 spectrum).  For theta
!= 0 the blocks are complex Hermitian and roughly 4x slower; the script says
which path it took.

Units: hbar = 1, energies in MHz.
"""

import numpy as np
from itertools import product
from scipy.special import jv
from timeit import default_timer as timer

# =====================================================================
# 0.  SETTINGS
# =====================================================================
L = 10
usePBC = False
d = 3

WP = 40.0                 # drive frequency [MHz]
THETA = 0.0               # drive phase; statistical angle phi = -theta
N_RES = 1
CACHE = 'heff2_E_vs_N_L%d.npz' % L

# =====================================================================
# 1.  H_eff^(2) COUPLINGS
# =====================================================================
C00, C01, C11 = 0.0, 20.0, 15.0
C02 = 2.0 * C11
C12, C22 = 45.0, 55.0
Om_a, Om_b = 1.0, 1.0
A = 48.0
COMPENSATE_11BLOCK = True

LMAX = 60
ells = np.arange(-LMAX, LMAX + 1)


def heff2_couplings(wp):
  kap = A / wp
  nu = 2 * C11 - C02 - N_RES * wp
  Jl = jv(ells, kap)

  def ss(dc):
    den = dc + nu + ells * wp
    if np.any(np.abs(den) < 1e-9):
      raise ZeroDivisionError('resonant Stark denominator at wp=%g' % wp)
    return float(np.sum(Jl**2 / den))

  E = {}
  E['00'] = -2 * Om_a**2 / (C01 - 2 * C00)
  E['01'] = (Om_a**2 / (C01 - 2 * C00) - Om_a**2 / (2 * C11 - C01)
             - Om_b**2 * ss(C02 - C01))
  E['02'] = Om_b**2 * ss(C02 - C01) - Om_a**2 / (C12 - C02)
  E['10'] = E['01']
  E['11'] = 2 * Om_a**2 / (2 * C11 - C01) - 2 * Om_b**2 * ss(C12 - 2 * C11)
  E['12'] = (Om_a**2 / (C12 - C02) + Om_b**2 * ss(C12 - 2 * C11)
             - Om_b**2 * ss(2 * C22 - C12))
  E['20'] = E['02']
  E['21'] = E['12']
  E['22'] = 2 * Om_b**2 * ss(2 * C22 - C12)
  if COMPENSATE_11BLOCK:
    E['11'] += E['20'] - E['11']

  g1 = Om_a**2 * (1.0 / (C01 - 2 * C00) - 1.0 / (2 * C11 - C01))
  pref = Om_a * Om_b * (1.0 / (2 * C11 - C01) - 1.0 / (C12 - C02))
  g2 = pref * jv(N_RES, kap)
  g3 = pref * jv(N_RES, kap)
  Sig = (1.0 / (C12 - 2 * C11 + nu + ells * wp)
         - 1.0 / (2 * C22 - C12 + nu + ells * wp))
  g4 = Om_b**2 * float(np.sum(Jl * jv(ells, kap) * Sig))
  Eb = np.array([[E['%d%d' % (a, b)] for b in range(d)] for a in range(d)])
  return Eb, (g1, g2, g3, g4)


def phases_of(theta):
  a, b = 0.0, theta                        # alpha = 0, beta = theta
  return -a, b - 2 * a, a - 2 * b, -b


# =====================================================================
# 2.  SECTOR HAMILTONIAN
# =====================================================================
configs = np.array(list(product(range(d), repeat=L)), dtype=np.int8)
Nvec = configs.sum(axis=1).astype(int)
PW = np.array([d**(L - 1 - j) for j in range(L)], dtype=np.int64)
DIMFULL = d**L
BONDS = [(j, j + 1) for j in range(L - 1)]
if usePBC:
  BONDS.append((L - 1, 0))

TH = phases_of(THETA)
REAL_PATH = all(abs(np.sin(t)) < 1e-14 for t in TH)


def sector_levels(N, Eb, gs):
  """All eigenvalues of H_eff^(2) in the sector of particle number N."""
  rows = np.where(Nvec == N)[0]
  dim = len(rows)
  pos = np.full(DIMFULL, -1, dtype=np.int64)
  pos[rows] = np.arange(dim)
  cfg = configs[rows]
  g1, g2, g3, g4 = gs
  dtype = np.float64 if REAL_PATH else np.complex128
  H = np.zeros((dim, dim), dtype=dtype)
  diag = np.zeros(dim)
  for (jL, jR) in BONDS:
    diag += Eb[cfg[:, jL], cfg[:, jR]]
  H[np.arange(dim), np.arange(dim)] = diag
  for (jL, jR) in BONDS:
    for (a, b), (a2, b2), g, t in [((0, 1), (1, 0), g1, TH[0]),
                                   ((1, 1), (2, 0), g2, TH[1]),
                                   ((0, 2), (1, 1), g3, TH[2]),
                                   ((1, 2), (2, 1), g4, TH[3])]:
      if g == 0.0:
        continue
      amp = g * np.cos(t) if REAL_PATH else g * np.exp(1j * t)
      sel = np.where((cfg[:, jL] == a) & (cfg[:, jR] == b))[0]
      if sel.size == 0:
        continue
      tgt = pos[rows[sel] + (a2 - a) * PW[jL] + (b2 - b) * PW[jR]]
      H[tgt, sel] += amp
      H[sel, tgt] += np.conj(amp)
  tr = float(np.trace(H).real)
  return np.linalg.eigvalsh(H), tr


# =====================================================================
# 3.  FULL SPECTRUM, SECTOR BY SECTOR
# =====================================================================
Eb, gs = heff2_couplings(WP)
g1, g2, g3, g4 = gs

print('=' * 78)
print('E vs N-hat for H_eff^(2):  L = %d, %s, dim = 3^%d = %d, omega\' = %g MHz'
      % (L, 'PBC' if usePBC else 'OBC', L, DIMFULL, WP))
print('=' * 78)
print('couplings: g1=%+.6f  g2=g3=%+.6f  g4=%+.6f MHz' % (g1, g2, g4))
print('theta = %.3f  ->  phi = %.3f;  LAPACK path: %s'
      % (THETA, -THETA, 'REAL symmetric' if REAL_PATH else 'complex Hermitian'))

try:
  z = np.load(CACHE)
  if int(z['L']) == L and abs(float(z['theta']) - THETA) < 1e-12 \
     and abs(float(z['wp']) - WP) < 1e-12:
    E, NLAB, TRACE = z['E'], z['NLAB'], float(z['TRACE'])
    print('loaded cached spectrum from %s' % CACHE)
  else:
    raise KeyError
except (FileNotFoundError, KeyError, ValueError):
  Es, Ns, TRACE = [], [], 0.0
  t0 = timer()
  for N in range(2 * L + 1):
    e, tr = sector_levels(N, Eb, gs)
    Es.append(e)
    Ns.append(np.full(e.size, N))
    TRACE += tr
    print('   N = %2d  dim %5d   E in [%+.5f, %+.5f]   (%.1f s)'
          % (N, e.size, e[0], e[-1], timer() - t0))
  E = np.concatenate(Es)
  NLAB = np.concatenate(Ns)
  np.savez(CACHE, L=L, theta=THETA, wp=WP, E=E, NLAB=NLAB, TRACE=TRACE)
  print('   cached -> %s' % CACHE)

order = np.argsort(E)
Esort, Nsort = E[order], NLAB[order]
Ns_ = np.arange(2 * L + 1)
E0N = np.array([E[NLAB == N].min() for N in Ns_])
EmaxN = np.array([E[NLAB == N].max() for N in Ns_])
DIMN = np.array([(NLAB == N).sum() for N in Ns_])

# =====================================================================
# 4.  CHECKS
# =====================================================================
print('-' * 78)
print('CHECKS')
print('   sum of sector dims = %d  (must be 3^%d = %d)   %s'
      % (E.size, L, DIMFULL, 'PASS' if E.size == DIMFULL else 'FAIL'))
print('   every level has an integer N by construction (H is block diagonal)')
print('   trace sum rule sum_m E_m = Tr H:  |d| = %.2e   %s'
      % (abs(E.sum() - TRACE), 'PASS' if abs(E.sum() - TRACE) < 1e-6 else 'FAIL'))
print('   particle-hole check dim(N) == dim(2L-N):        %s'
      % ('PASS' if np.array_equal(DIMN, DIMN[::-1]) else 'FAIL'))

# =====================================================================
# 5.  RESULTS
# =====================================================================
igs = int(np.argmin(E))
print('-' * 78)
print('RESULTS')
print('   global ground state E_0 = %+.8f MHz  in sector N = %d  (n-bar = %.3f)'
      % (E[igs], NLAB[igs], NLAB[igs] / L))
print('   global maximum      E   = %+.8f MHz  in sector N = %d'
      % (E.max(), NLAB[int(np.argmax(E))]))
print('   bandwidth = %.8f MHz' % (E.max() - E.min()))
print('   per-sector summary:')
print('      %3s %6s %12s %12s %12s' % ('N', 'dim', 'E_0(N)', 'E_max(N)', 'width'))
for N in Ns_:
  print('      %3d %6d %+12.6f %+12.6f %12.6f'
        % (N, DIMN[N], E0N[N], EmaxN[N], EmaxN[N] - E0N[N]))
print('   E_0(N) is minimized at N = %d;  E_0 is NOT monotonic in N'
      % int(np.argmin(E0N)))
print('   the density operator n-bar = N/L has eigenvalues %.3f .. %.3f in '
      'steps of 1/L = %.3f' % (0.0, 2.0, 1.0 / L))
print('=' * 78)

# =====================================================================
# 6.  FIGURE
# =====================================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

CMAP = LinearSegmentedColormap.from_list(
    'white_violet', ['#ffffff', '#efe3f5', '#cbaade', '#9c6ac6', '#6a2fa0',
                     '#2c0a52'])

fig, axs = plt.subplots(2, 2, figsize=(15.5, 9.4))
rng = np.random.default_rng(0)

# (a) E vs N-hat eigenvalue
ax = axs[0, 0]
ax.plot(NLAB + 0.30 * (rng.random(NLAB.size) - 0.5), E, '.', ms=1.0,
        color='#9c6ac6', alpha=0.30, rasterized=True)
ax.plot(Ns_, E0N, 'o-', color='#2c0a52', ms=5, lw=1.5, label=r'$E_0(N)$')
ax.plot(Ns_, EmaxN, 's--', color='#c2185b', ms=4, lw=1.2, label=r'$E_{\max}(N)$')
ax.plot([NLAB[igs]], [E[igs]], 'k*', ms=16, label='global ground state')
ax.set_xlabel(r'$\hat N$ eigenvalue   $N=\sum_j n_j$')
ax.set_ylabel(r'$E$  [MHz]')
ax.set_title(r'(a) All $3^{%d}=%d$ levels vs $\hat N$ (extensive)'
             % (L, DIMFULL))
ax.set_xticks(Ns_[::2])
ax.legend(fontsize=8)
ax.grid(True, ls='--', alpha=0.35)

# (b) density of states in the (N, E) plane
ax = axs[0, 1]
nb = 160
H2, xe, ye = np.histogram2d(NLAB, E, bins=[np.arange(2 * L + 2) - 0.5,
                                           np.linspace(E.min(), E.max(), nb)])
im = ax.imshow(H2.T, origin='lower', aspect='auto', cmap=CMAP,
               norm=LogNorm(vmin=1, vmax=max(H2.max(), 2)),
               extent=[-0.5, 2 * L + 0.5, E.min(), E.max()])
fig.colorbar(im, ax=ax, label='levels per bin (log)')
ax.plot(Ns_, E0N, '-', color='#c2185b', lw=1.4)
ax.set_xlabel(r'$\hat N$ eigenvalue')
ax.set_ylabel(r'$E$  [MHz]')
ax.set_title('(b) Density of states resolved by $N$')
ax.set_xticks(Ns_[::2])

# (c) intensive view: E/L vs n-bar = N/L
ax = axs[1, 0]
ax.plot(NLAB / L + (0.30 / L) * (rng.random(NLAB.size) - 0.5), E / L, '.',
        ms=1.0, color='#9c6ac6', alpha=0.30, rasterized=True)
ax.plot(Ns_ / L, E0N / L, 'o-', color='#2c0a52', ms=5, lw=1.5,
        label=r'$E_0(N)/L$')
ax.plot(Ns_ / L, EmaxN / L, 's--', color='#c2185b', ms=4, lw=1.2,
        label=r'$E_{\max}(N)/L$')
ax.plot([NLAB[igs] / L], [E[igs] / L], 'k*', ms=16)
ax.axvline(1.0, color='k', ls=':', lw=1, label='unit filling')
ax.set_xlabel(r'average number density  $\bar n=\langle\hat N\rangle/L$'
              '\n' r'(eigenvalues of $\hat N/L$, spaced by $1/L=%.2f$)' % (1 / L))
ax.set_ylabel(r'$E/L$  [MHz per site]')
ax.set_title(r'(c) Intensive view: $E/L$ vs $\bar n$')
ax.legend(fontsize=8)
ax.grid(True, ls='--', alpha=0.35)

# (d) per-sector envelope and dimension
ax = axs[1, 1]
ax.plot(Ns_, EmaxN - E0N, 'o-', color='#6a2fa0', lw=1.6, ms=5,
        label='sector width $E_{\\max}-E_0$')
ax.set_xlabel(r'$\hat N$ eigenvalue')
ax.set_ylabel('sector width  [MHz]')
ax.set_xticks(Ns_[::2])
ax.grid(True, ls='--', alpha=0.35)
ax2 = ax.twinx()
ax2.bar(Ns_, DIMN, color='#cbaade', alpha=0.55, width=0.7)
ax2.set_ylabel('sector dimension')
ax.set_title('(d) Sector width and dimension vs $N$')
ax.legend(fontsize=8, loc='upper left')
ax.set_zorder(ax2.get_zorder() + 1)
ax.patch.set_visible(False)

fig.suptitle(r'Spectrum of $H_{\mathrm{eff}}^{(2)}$ against the particle-number '
             r'operator:  $L=%d$, %s, $\omega^{\prime}=%g$ MHz, '
             r'$\theta=%.2f\pi$.  $[H,\hat N]=0$, so every eigenstate has an '
             r'integer $N$; $\bar n=\hat N/L$ is that axis rescaled by $1/L$'
             % (L, 'PBC' if usePBC else 'OBC', WP, THETA / np.pi), fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig('heff2_E_vs_N.png', dpi=150)
print('Saved -> heff2_E_vs_N.png')
