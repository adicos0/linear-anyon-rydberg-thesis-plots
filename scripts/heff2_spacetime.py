"""
heff2_spacetime.py
---------------------------------------------------------------------
Spacetime map of the average particle number: sites on the y axis, time on
the x axis, colour = <n_j(t)>, evolving under H_eff^(2).

    heff2_spacetime.png       single panel at the ED operating point
    heff2_spacetime_phi.png   the same map at several statistical angles

DYNAMICS
--------
N = sum_j n_j is exactly conserved, so the evolution stays inside the sector
of the initial state and is done EXACTLY by diagonalizing that sector:

    |psi(t)> = sum_m e^{-i E_m t} |m><m|psi(0)>,
    <n_j(t)> = sum_i |psi_i(t)|^2 (n_j)_i

Initial state (the convention of anyon_phi_scan.py): one doublon at site SITE0
of an open chain, so N = 2.  A single particle at the same site is evolved
alongside as a reference: it can
only use the A_1 channel, so it shows the bare single-particle light cone with
no doublon physics in it.

COUPLINGS AND PHASES: exactly as in heff2_occupation_maps.py.  The magnitudes
E_ab and the four amplitudes A_1..A_4 (thesis eq:ampconditions; here carried
with their signs) come from the H_eff^(2) derivation at the given omega' with
the resonance order locked; the phases use the anyon-Hubbard assignment
th1 = th3 = 0, th2 = th4 = phi, which puts a flux phi on every plaquette of
the configuration graph.  The site ramp would instead give flux = 0 -- pure
gauge -- and since a gauge transformation exp(i sum_j lambda_j D_j) is
DIAGONAL in the occupation basis it cannot change |psi_i|^2, so <n_j(t)>
would be exactly independent of alpha and beta.  That is checked in part 4.

Exact symmetries, asserted as tests:
  * total density sum_j <n_j(t)> = N at all times.
  * phi = 0 and phi = 2 pi must coincide (e^{i phi n} with integer n).
  * H(-phi) = H(phi)^* and psi0 is real, so psi(t;-phi) = psi(-t;phi)^* and
    hence  <n_j(t; -phi)> = <n_j(-t; phi)>.  Note this does NOT make the
    density even in phi: that would additionally require the dynamics to be
    time-reversal symmetric, which it is not.  Reversing the statistical
    angle runs the density backwards in time rather than leaving it alone,
    so phi must be sampled over the whole circle [0, 2 pi), not just [0, pi].

Units: hbar = 1, energies in MHz, so time is in the reciprocal unit.
"""

import numpy as np
from itertools import product
from scipy.special import jv

# =====================================================================
# 0.  SETTINGS
# =====================================================================
L = 10
usePBC = False
d = 3

WP = 40.0                 # drive frequency of the main figure [MHz]
PHI_MAIN = 0.0            # statistical angle of the main figure
N_RES = 1                 # locked resonance order
SITE0 = 4                 # where the doublon starts
TMAX = 600.0              # total evolution time
NT = 900                  # time samples

PHI_PANELS = [k * np.pi / 3.0 for k in range(6)]   # 0 .. 5pi/3

# =====================================================================
# 1.  H_eff^(2) COUPLINGS AT omega'   (as in heff2_occupation_maps.py)
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

  def ss(dconst):
    den = dconst + nu + ells * wp
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


def phases_of(mode, ang, j, alpha=0.0, beta=0.0):
  """(th1, th2, th3, th4) on bond j.

  'drive' is the physical prescription: set the drive phases to alpha = 0 and
  beta = theta and replicate the two-site bond operator on every bond, giving
      (th1, th2, th3, th4) = (-alpha, beta-2alpha, alpha-2beta, -beta)
                           = (0, theta, -2 theta, -theta).
  The gauge-invariant plaquette flux is then -2 th1 + th2 + th3 = -theta, and
  the diagonal gauge U = exp(i sum_j f(n_j)) with f = (0, gamma1, gamma2) and
  gamma2 - 2 gamma1 = -2 theta carries this exactly onto the anyon-Hubbard
  Peierls pattern (0, phi, 0, phi) with

      phi = -theta .

  So the drive realizes the anyon-Hubbard model with statistical angle -theta;
  this is asserted numerically in the tests.  ('anyon' is that target pattern,
  kept for the equivalence check.  The site 'ramp' has flux 0 on every
  plaquette - pure gauge, no statistics.)
  """
  if mode == 'drive':
    a, b = 0.0, ang
    return -a, b - 2 * a, a - 2 * b, -b
  if mode == 'anyon':
    return 0.0, ang, 0.0, ang
  a_j, a_j1 = (j + 1) * alpha, (j + 2) * alpha
  b_j, b_j1 = (j + 1) * beta, (j + 2) * beta
  return a_j - a_j1, b_j - a_j1, a_j - b_j1, b_j - b_j1


def sector(N):
  rows = np.where(Nvec == N)[0]
  pos = np.full(DIMFULL, -1, dtype=np.int64)
  pos[rows] = np.arange(len(rows))
  return rows, pos, configs[rows]


def sector_H(N, Eb, gs, mode, phi, alpha=0.0, beta=0.0):
  """Dense H_eff^(2) restricted to sector N (exact, since [H,N]=0)."""
  rows, pos, cfg = sector(N)
  dim = len(rows)
  g1, g2, g3, g4 = gs
  H = np.zeros((dim, dim), dtype=complex)
  diag = np.zeros(dim)
  for (jL, jR) in BONDS:
    diag += Eb[cfg[:, jL], cfg[:, jR]]
  H[np.arange(dim), np.arange(dim)] = diag
  for jb, (jL, jR) in enumerate(BONDS):
    th = phases_of(mode, phi, jb, alpha, beta)
    for (a, b), (a2, b2), g, t in [((0, 1), (1, 0), g1, th[0]),
                                   ((1, 1), (2, 0), g2, th[1]),
                                   ((0, 2), (1, 1), g3, th[2]),
                                   ((1, 2), (2, 1), g4, th[3])]:
      if g == 0.0:
        continue
      amp = g * np.exp(1j * t)
      sel = np.where((cfg[:, jL] == a) & (cfg[:, jR] == b))[0]
      if sel.size == 0:
        continue
      tgt = pos[rows[sel] + (a2 - a) * PW[jL] + (b2 - b) * PW[jR]]
      H[tgt, sel] += amp
      H[sel, tgt] += np.conj(amp)
  return H, cfg, rows, pos


def evolve_density(N, cfgstate, Eb, gs, mode, phi, times, alpha=0.0, beta=0.0):
  """<n_j(t)> for the Fock initial state cfgstate, exactly."""
  H, cfg, rows, pos = sector_H(N, Eb, gs, mode, phi, alpha, beta)
  assert np.allclose(H, H.conj().T, atol=1e-13), 'sector H not Hermitian'
  ev, V = np.linalg.eigh(H)
  i0 = pos[int(np.dot(cfgstate, PW))]
  assert i0 >= 0, 'initial state is not in sector N=%d' % N
  psi0 = np.zeros(len(rows), dtype=complex)
  psi0[i0] = 1.0
  c = V.conj().T @ psi0
  psit = (np.exp(-1j * np.outer(times, ev)) * c) @ V.T     # (Nt, dim)
  p = np.abs(psit)**2
  return p @ cfg.astype(float), ev                          # (Nt, L)


# =====================================================================
# 3.  RUN
# =====================================================================
Eb, gs = heff2_couplings(WP)
g1, g2, g3, g4 = gs
times = np.linspace(0.0, TMAX, NT)

doublon = np.zeros(L, dtype=int)
doublon[SITE0] = 2                       # |0,0,0,0,2,0,0,0>, N = 2
single = np.zeros(L, dtype=int)
single[SITE0] = 1                        # |0,0,0,0,1,0,0,0>, N = 1

dens_main, ev2 = evolve_density(2, doublon, Eb, gs, 'drive', PHI_MAIN, times)
dens_sing, ev1 = evolve_density(1, single, Eb, gs, 'drive', PHI_MAIN, times)
dens_phi = {ph: evolve_density(2, doublon, Eb, gs, 'anyon', ph, times)[0]
            for ph in PHI_PANELS}

print('=' * 78)
print('Spacetime map of <n_j(t)> under H_eff^(2):  L = %d, %s, omega\' = %g MHz'
      % (L, 'PBC' if usePBC else 'OBC', WP))
print('=' * 78)
print('couplings: g1=%+.6f  g2=g3=%+.6f  g4=%+.6f MHz  (kappa = %.4f)'
      % (g1, g2, g4, A / WP))
print('AC-Stark diagonal E_ab [MHz]:')
for a in range(d):
  print('   ' + '  '.join('E_%d%d=%+8.5f' % (a, b, Eb[a, b]) for b in range(d)))
print('initial state |%s>  ->  sector N = %d, dim = %d'
      % (','.join(str(x) for x in doublon), 2, len(sector(2)[0])))
print('evolution: exact, %d time samples on [0, %g]' % (NT, TMAX))

# =====================================================================
# 4.  TESTS
# =====================================================================
print('-' * 78)
print('TESTS')
FILLING = 8
Eb40, gs40 = Eb, gs

# ---- the drive prescription IS the anyon-Hubbard model, with phi = -theta ----
# U = exp(i sum_j f(n_j)), f = (0, gamma1, gamma2), gamma2 - 2 gamma1 = -2 theta.
# Being diagonal in the occupation basis, U cannot change |psi_i|^2, so every
# occupation observable below is literally the anyon-Hubbard one at phi = -theta.
def _gauge_check():
  worst = 0.0
  for _th in [0.31, 0.73, np.pi / 2, 2.4, 4.1]:
    for _N in [2, FILLING]:
      _Hd = sector_H(_N, Eb40, gs40, 'drive', _th)[0]
      _Ha, _cfg = sector_H(_N, Eb40, gs40, 'anyon', -_th)[0:2]
      _f = np.array([0.0, 0.0, 2 * _th])          # sign per U vs U^dag convention
      _p = np.exp(1j * _f[_cfg].sum(axis=1))
      worst = max(worst, np.abs((_p.conj()[:, None] * _Hd) * _p[None, :] - _Ha).max())
  return worst


_w = _gauge_check()
print('   drive(theta) == anyon-Hubbard(phi = -theta) up to the diagonal gauge:')
print('      max ||U^dag H_drive U - H_anyon|| = %.2e   %s' %
      (_w, 'PASS' if _w < 1e-12 else 'FAIL'))
_fl = phases_of('drive', 0.73, 0)
print('      plaquette flux at theta=0.73: %+.6f  (= -theta)   %s'
      % (-2 * _fl[0] + _fl[1] + _fl[2],
         'PASS' if abs((-2 * _fl[0] + _fl[1] + _fl[2]) + 0.73) < 1e-12 else 'FAIL'))

tot = np.abs(dens_main.sum(axis=1) - 2.0).max()
print('   density conservation  sum_j <n_j(t)> = 2 : max dev %.2e   %s'
      % (tot, 'PASS' if tot < 1e-10 else 'FAIL'))
d0 = evolve_density(2, doublon, Eb, gs, 'drive', 0.0, times)[0]
d2pi = evolve_density(2, doublon, Eb, gs, 'drive', 2 * np.pi, times)[0]
print('   phi = 0 vs phi = 2pi (must coincide): max |dn| = %.2e   %s'
      % (np.abs(d0 - d2pi).max(),
         'PASS' if np.abs(d0 - d2pi).max() < 1e-10 else 'FAIL'))
dm = evolve_density(2, doublon, Eb, gs, 'drive', -0.7, times)[0]
dp = evolve_density(2, doublon, Eb, gs, 'drive', 0.7, -times)[0]
print('   n(t;-phi) = n(-t;phi)  [H(-phi)=H(phi)^*, psi0 real]: max |dn| = %.2e   %s'
      % (np.abs(dm - dp).max(),
         'PASS' if np.abs(dm - dp).max() < 1e-10 else 'FAIL'))
dsame = evolve_density(2, doublon, Eb, gs, 'drive', 0.7, times)[0]
print('   ... but the density is NOT even in phi: max |n(0.7)-n(-0.7)| = %.3f'
      % np.abs(dsame - dm).max())
print('       (that would need time-reversal symmetry, which H_eff^(2) lacks,')
print('        so phi is sampled over the full circle below)')
r0 = evolve_density(2, doublon, Eb, gs, 'ramp', 0.0, times, 0.0, 0.0)[0]
r1 = evolve_density(2, doublon, Eb, gs, 'ramp', 0.0, times, np.pi / 4, np.pi / 3)[0]
print('   site ramp is pure gauge (drive phases cannot move the density):')
print('      max |dn| between (alpha,beta)=(0,0) and (pi/4,pi/3) = %.2e   %s'
      % (np.abs(r0 - r1).max(),
         'PASS' if np.abs(r0 - r1).max() < 1e-10 else 'FAIL'))

# =====================================================================
# 5.  WHAT THE MAP SHOWS
# =====================================================================
print('-' * 78)
print('RESULTS')
print('   <n_j(t)> range: %.6f .. %.6f' % (dens_main.min(), dens_main.max()))
print('   doublon survival <n_%d(t)>: start %.4f, min %.4f, final %.4f'
      % (SITE0, dens_main[0, SITE0], dens_main[:, SITE0].min(),
         dens_main[-1, SITE0]))
esc = dens_main[:, :SITE0].sum(axis=1) + dens_main[:, SITE0 + 1:].sum(axis=1)
print('   density that has left site %d: max %.4f (of 2) at t = %.1f'
      % (SITE0, esc.max(), times[int(np.argmax(esc))]))
lw = dens_main[:, :SITE0].sum(axis=1)
rw = dens_main[:, SITE0 + 1:].sum(axis=1)
print('   left/right split at t = %g: %.4f / %.4f  (asymmetry %.4f)'
      % (TMAX, lw[-1], rw[-1], lw[-1] - rw[-1]))
print('   effect of the statistical angle on the density:')
for ph in PHI_PANELS:
  print('      phi = %.2f pi : max_t,j |n(phi) - n(0)| = %.4f'
        % (ph / np.pi, np.abs(dens_phi[ph] - dens_phi[0.0]).max()))
print('   single-particle reference (N=1): only the g1 channel is available,')
print('      <n_j(t)> range %.4f .. %.4f' % (dens_sing.min(), dens_sing.max()))
print('=' * 78)

# =====================================================================
# 6.  FIGURES
# =====================================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15, 'axes.titlesize': 17, 'axes.labelsize': 16,
                     'xtick.labelsize': 13, 'ytick.labelsize': 13,
                     'legend.fontsize': 14, 'figure.titlesize': 17})
from matplotlib.colors import LinearSegmentedColormap

# Density colour scale: pure white at <n_j> = 0 rising to deep violet at the
# maximum, so empty sites read as blank page rather than as a filled dark
# background and the occupied structure is what carries ink.
CMAP_DENS = LinearSegmentedColormap.from_list(
    'white_violet', ['#ffffff', '#efe3f5', '#cbaade', '#9c6ac6', '#6a2fa0',
                     '#2c0a52'])

# The four second-order coupling channels.  Labelled A_1..A_4 to match the
# thesis notation of eq:ampconditions, where A_i are defined as the MAGNITUDES
# |.| of these amplitudes; the signed values are printed in GSIGNED below.
GLABEL = (r'$A_1$: $|01\rangle\!\to\!|10\rangle$ exchange   '
          r'$A_2$: $|11\rangle\!\to\!|20\rangle$ Floquet-assisted   '
          r'$A_3$: $|02\rangle\!\to\!|11\rangle$ Floquet-assisted   '
          r'$A_4$: $|12\rangle\!\to\!|21\rangle$ cross-site')

sites = np.arange(L)
EXT = [times[0], times[-1], -0.5, L - 0.5]

# ---------------------------------------------------------------- fig 1
fig, axs = plt.subplots(2, 1, figsize=(11.5, 8.2), sharex=True,
                        gridspec_kw={'height_ratios': [1, 1]})

ax = axs[0]
im = ax.imshow(dens_main.T, origin='lower', aspect='auto', extent=EXT,
               cmap=CMAP_DENS, vmin=0.0, vmax=2.0)
fig.colorbar(im, ax=ax, label=r'$\langle \hat n_j(t)\rangle$')
ax.set_yticks(sites)
ax.set_ylabel('site $j$')
KET = ','.join(str(x) for x in doublon)      # built from the state, not hardcoded
ax.set_title(r'(a) Doublon at site %d:  $|%s\rangle$, $N=2$ (sector dim %d)'
             % (SITE0, KET, len(sector(2)[0])))

ax = axs[1]
im = ax.imshow(dens_sing.T, origin='lower', aspect='auto', extent=EXT,
               cmap=CMAP_DENS, vmin=0.0, vmax=1.0)
fig.colorbar(im, ax=ax, label=r'$\langle \hat n_j(t)\rangle$')
ax.set_yticks(sites)
ax.set_ylabel('site $j$')
ax.set_xlabel(r'time  $t$  [1/MHz]')
ax.set_title(r'(b) Single particle at site %d (reference: only the $A_1$ '
             r'channel)' % SITE0)

GVALS = (r'$A_1=%.4f$,  $A_2=A_3=%.4f$,  $A_4=%.4f$ MHz'
         '   ' r'(signed: $%+.4f$, $%+.4f$, $%+.4f$)'
         % (abs(g1), abs(g2), abs(g4), g1, g2, g4))
_head = (r'Spacetime map of $\langle\hat n_j(t)\rangle$ under '
         r"$H_{\mathrm{eff}}^{(2)}$:  $L=%d$, %s, $\omega'=%g$ MHz, "
         r'$\theta=%.2f\pi$ ($\phi=-\theta$)'
         % (L, 'PBC' if usePBC else 'OBC', WP, PHI_MAIN / np.pi))
# The run header (_head, GLABEL, GVALS) is deliberately NOT drawn: at thesis
# scale it is illegible, and the same information is in the caption and in the
# parameter table of Sec. 2.2.  It is printed to stdout instead.
print(_head); print(GLABEL); print(GVALS)
fig.tight_layout()
fig.savefig('heff2_spacetime.png', dpi=150)
print('Saved -> heff2_spacetime.png')

# ---------------------------------------------------------------- fig 2
n = len(PHI_PANELS)
fig, axs = plt.subplots(2, 3, figsize=(17, 7.6), sharex=True, sharey=True)
axs = axs.ravel()
for ax, ph in zip(axs, PHI_PANELS):
  im = ax.imshow(dens_phi[ph].T, origin='lower', aspect='auto', extent=EXT,
                 cmap=CMAP_DENS, vmin=0.0, vmax=2.0)
  ax.set_yticks(sites)
  ax.set_title(r'$\theta=%.2f\pi$' % (ph / np.pi))
  ax.set_xlabel(r'$t$  [1/MHz]')
  ax.set_ylabel('site $j$')
fig.colorbar(im, ax=axs[:n].tolist(), label=r'$\langle \hat n_j(t)\rangle$',
             fraction=0.02)

_head2 = (r'Doublon spacetime dynamics under $H_{\mathrm{eff}}^{(2)}$ vs '
          r"the drive phase $\theta$ ($\alpha=0$, $\beta=\theta$; "
          r"statistical angle $\phi=-\theta$):  $L=%d$, %s, "
          r"$\omega'=%g$ MHz, $N=2$" % (L, 'PBC' if usePBC else 'OBC', WP))
print(_head2); print(GLABEL); print(GVALS)
fig.subplots_adjust(top=0.95, hspace=0.30, wspace=0.22,
                    left=0.05, right=0.88, bottom=0.11)
fig.savefig('heff2_spacetime_phi.png', dpi=150)
print('Saved -> heff2_spacetime_phi.png')
