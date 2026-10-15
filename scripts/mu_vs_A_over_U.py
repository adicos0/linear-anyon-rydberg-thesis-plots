"""
mu_vs_A_over_U.py
---------------------------------------------------------------------
mu_eff against A_i/U_eff, for the three distinct channels.

    mu_vs_A_over_U.png

Three panels: A_1, the pair A_2 = A_3, and A_4.  A_3 = A_2 identically, because
the resonance lock nu~_0 = nu~_1 gives N_0 = N_1 = N and both channels carry
the same J_N(kappa), so the pair gets a single panel rather than two copies of
the same curve.  (They are equal in MAGNITUDE only: in the Hamiltonian they sit
on different bonds and carry different phases, phi and -2 phi.)

SCALE
-----
The y axis is the bare mu_eff in MHz, which is NOT invariant under an overall
rescaling of the drive: under (Om_a, Om_b) -> lam (Om_a, Om_b) every A_i,
mu_eff and U_eff carries lam^2, so the x axis A_i/U_eff is untouched while
every curve slides vertically by lam^2.  The figure is therefore drawn at a
stated Om_a = 1 MHz, the value used in all the ED runs of this project;
changing it moves the curves vertically and does nothing else.

PARAMETRIZATION
---------------
With x = r Sigma(kappa) and r = Om_b^2/Om_a^2,

    U_eff   = Om_a^2 (3x - b)
    mu_eff  = Om_a^2 (a - x)
    A_1/U_eff = K1 / (3x - b)
    A_2/U_eff = A_3/U_eff = sqrt(x/Sigma) K2 |J_N(kappa)| / (3x - b)
    A_4/U_eff = (x/Sigma) |Sigma_4(kappa)| / (3x - b)

a = 3 Delta_1 - Delta_2, b = 4 Delta_1 + Delta_3 - 2 Delta_2,
K1 = |Delta_1 - Delta_2|, K2 = |Delta_2 - Delta_3|.

Each curve sweeps r at fixed kappa.  U_eff passes through zero at x = b/3, so
A_i/U_eff diverges and CHANGES SIGN there; the x axis is symmetric-log so both
signs and the divergence are visible.  Note A_i > 0 always, so the sign of
A_i/U_eff is the sign of U_eff: the left half-plane is the attractive side.

Units: MHz.
"""

import numpy as np
import math
from scipy.special import jv

C00, C01, C11 = 0.0, 20.0, 15.0
C02 = 2.0 * C11
C12, C22 = 45.0, 55.0
wp = 40.0
N_RES = 1
nu = 2 * C11 - C02 - N_RES * wp
ells = np.arange(-80, 81)

D1 = 1.0 / (C01 - 2 * C00)
D2 = 1.0 / (2 * C11 - C01)
D3 = 1.0 / (C12 - C02)
K1 = abs(D1 - D2)
K2 = abs(D2 - D3)
A_ = 3 * D1 - D2
B_ = 4 * D1 + D3 - 2 * D2

KAPPAS = [0.8, 1.2, 2.0]
KMAIN = 1.2
COL = {0.8: '#7fb3d5', 1.2: '#12395e', 2.0: '#2e86c1'}
OM_A = 1.0                      # Rabi scale (MHz); the curves slide as lam^2
# Validity of the second-order expansion.  H_eff^(2) needs Omega << every
# energy denominator it contains; the smallest of those is 10 MHz here, so
# Om_a = 1 MHz puts the expansion parameter at 0.1 and the derivation holds.
DEN_MIN = 10.0


def Sigma(k):
  return float(np.sum(jv(ells, k)**2 / (C02 - C01 + nu + ells * wp)))


def Sigma4(k):
  return float(np.sum(jv(ells, k)**2
                      * (1.0 / (C12 - 2 * C11 + nu + ells * wp)
                         - 1.0 / (2 * C22 - C12 + nu + ells * wp))))


def curves(k, r):
  """A_i/U_eff (4), mu_eff/U_eff, mu_eff [MHz] at Om_a = OM_A."""
  S, S4 = Sigma(k), Sigma4(k)
  x = r * S
  den = 3.0 * x - B_
  a1 = K1 / den
  a2 = np.sqrt(r) * K2 * abs(jv(N_RES, k)) / den
  a4 = r * abs(S4) / den
  mu_over_U = (A_ - x) / den
  mu = OM_A**2 * (A_ - x)
  return (a1, a2, a2, a4), mu_over_U, mu


print('=' * 76)
print('mu_eff against A_i/U_eff')
print('=' * 76)
print('U_eff  = Om_a^2 (3x - b),   mu_eff = Om_a^2 (a - x),   x = r Sigma')
print('a = %.6f, b = %.6f /MHz;  U_eff = 0 at x = b/3 = %.6f'
      % (A_, B_, B_ / 3))
print()
print('y axis is the bare mu_eff, which scales as lam^2 under')
print('   (Om_a, Om_b) -> lam (Om_a, Om_b), so it is drawn at Om_a = %g MHz;'
      % OM_A)
print('   the x axis A_i/U_eff is invariant and does not move with Om_a.')
print()
(a1, a2, a3, a4), mU, mu = curves(KMAIN, 1.0)
print('reference point Om_a = Om_b = %g MHz, kappa = %.1f:' % (OM_A, KMAIN))
print('   U_eff = %+.6f MHz, mu_eff = %+.6f MHz, mu_eff/U_eff = %+.6f'
      % (OM_A**2 * (3 * Sigma(KMAIN) - B_), mu, mU))
print('   A_1/U_eff = %+.4f   A_2/U_eff = A_3/U_eff = %+.4f   A_4/U_eff = %+.4f'
      % (a1, a2, a4))
print('   (all negative because U_eff < 0 here: the attractive side)')
print()
print('sign structure: A_i > 0 always, so sign(A_i/U_eff) = sign(U_eff).')
for k in KAPPAS:
  S = Sigma(k)
  rp = (B_ / 3) / S if S > 0 else None
  print('   kappa = %.1f: Sigma = %+.6f -> %s'
        % (k, S, 'U_eff = 0 at r = %.4f (Om_b/Om_a = %.4f); sign flips there'
           % (rp, math.sqrt(rp)) if rp else
           'U_eff < 0 for all r > 0 (no sign change)'))
print()
print('A_3 = A_2 identically (nu~_0 = nu~_1 -> N_0 = N_1), so the pair is drawn')
print('   in a single panel.  Equal in MAGNITUDE only: in H they sit on')
print('   different bonds and carry phases phi and -2 phi.')
print()
print('VALIDITY: smallest energy denominator in H_eff^(2) is %g MHz, so the'
      % DEN_MIN)
print('   expansion parameter is Om/DeltaE = %.2f at Om_a = %g MHz.'
      % (OM_A / DEN_MIN, OM_A))
print('   Second-order PT needs this << 1; it is %s here.'
      % ('fine' if OM_A / DEN_MIN < 0.2 else 'VIOLATED'))
print('=' * 76)

# ---- figure ----------------------------------------------------------
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axs = plt.subplots(1, 3, figsize=(14.4, 5.0))
COLS = [(0, r'$A_1/U_{\mathrm{eff}}$', r'$A_1$'),
        (1, r'$A_2/U_{\mathrm{eff}}\;(=A_3/U_{\mathrm{eff}})$', r'$A_2=A_3$'),
        (3, r'$A_4/U_{\mathrm{eff}}$', r'$A_4$')]
rr = np.geomspace(1e-3, 1e3, 200000)


def brs(v, jump):
  cut = np.where(np.abs(np.diff(v)) > jump)[0]
  return [s for s in np.split(np.arange(v.size), cut + 1) if s.size > 1]


for c, (i, xlab, short) in enumerate(COLS):
  ax = axs[c]
  for k in KAPPAS:
    (A, mU, mu) = curves(k, rr)
    xv = A[i]
    f = True
    for sg in brs(xv, 50.0):
      ax.plot(xv[sg], mu[sg], color=COL[k], lw=1.7,
              label=(r'$\kappa=%.1f$' % k) if f else None)
      f = False
  (A0, mU0, mu0) = curves(KMAIN, 1.0)
  ax.plot([A0[i]], [mu0], '*', color='#c2185b', ms=14, zorder=6,
          markeredgecolor='w', markeredgewidth=0.7)
  ax.axvline(0.0, color='k', lw=0.8)
  ax.axhline(0.0, color='gray', ls=':', lw=1.0)
  ax.set_xscale('symlog', linthresh=1e-2)
  ax.set_xlim(-1e3, 1e3)
  ax.set_ylim(-0.4 * OM_A**2, 0.15 * OM_A**2)
  ax.set_xlabel(xlab)
  if c == 0:
    ax.set_ylabel(r'$\mu_{\mathrm{eff}}$  [MHz]  at $\Omega_\alpha=%g$' % OM_A
                  + '\n' + r'(slides as $\lambda^2$)')
  ax.set_title((r'(%s) ' % 'abc'[c]) + short
               + ('\n' r'$A_3=A_2$ identically, one panel for the pair'
                  if c == 1 else '\n '), fontsize=9)
  ax.legend(fontsize=7, loc='lower right', framealpha=0.95)
  ax.grid(True, ls='--', alpha=0.35)

fig.suptitle(r'Bare $\mu_{\mathrm{eff}}$ against $A_i/U_{\mathrm{eff}}$, '
             r'sweeping $r=\Omega_\beta^2/\Omega_\alpha^2$ at fixed $\kappa$, '
             + (r'at $\Omega_\alpha=%g$ MHz.  ' % OM_A) +
             r'$A_i>0$ always, so the sign of $A_i/U_{\mathrm{eff}}$ is the '
             r'sign of $U_{\mathrm{eff}}$: the left half-plane is attractive.'
             + '\n' + (r'Stars: $\Omega_\alpha=\Omega_\beta=%g$ MHz, '
                       r'$\kappa=1.2$.  ' % OM_A) +
             r'$A_3=A_2$ identically, so the pair gets a single panel.',
             fontsize=9.5, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.925])
fig.savefig('mu_vs_A_over_U.png', dpi=170)
print('Saved -> mu_vs_A_over_U.png')
