"""
amplitudes_vs_muU_thesis.py
---------------------------------------------------------------------
Thesis figure for Sec. "Effective chemical potential, interaction, and their
ratio": the engineered amplitudes A_1, A_2 = A_3 and A_4 against the ratio
mu_eff/U_eff.

    figures/amplitudes_vs_muU.png

Writing the ratio of eq:ratio in terms of the single product x = r Sigma,

    mu_eff/U_eff = (a - x)/(3x - b),
    a = 3/(C01-2C00) - 1/(2C11-C01),
    b = 4/(C01-2C00) + 1/(C12-C02) - 2/(2C11-C01),

which is a Moebius function of x with a pole at x = b/3, a zero at x = a and
the limit -1/3 as |x| -> infinity.  The amplitudes of eq:ampconditions are

    A_1 = Om_a^2 |1/(C01-2C00) - 1/(2C11-C01)|
    A_2 = A_3 = Om_a Om_b |1/(2C11-C01) - 1/(C12-C02)| |J_N(kappa)|
    A_4 = Om_b^2 |Sigma_4(kappa)|

so A_1 is blind to Om_b and kappa while A_4 is blind to Om_a.  No single
two-knob sweep moves all of them, hence the per-panel sweeps stated in the
titles: Om_a for A_1, Om_b for A_2 = A_3 and A_4.

Units: MHz.
"""

import numpy as np
import math
from scipy.special import jv, jn_zeros
from scipy.optimize import brentq

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


def Sigma(k):
  return float(np.sum(jv(ells, k)**2 / (C02 - C01 + nu + ells * wp)))


def Sigma4(k):
  return float(np.sum(jv(ells, k)**2
                      * (1.0 / (C12 - 2 * C11 + nu + ells * wp)
                         - 1.0 / (2 * C22 - C12 + nu + ells * wp))))


def muU_x(x):
  return (A_ - x) / (3.0 * x - B_)


# sign change of Sigma: this is what opens the otherwise forbidden window
ks = np.linspace(0.05, 6.0, 6000)
sg = [Sigma(k) for k in ks]
kz = [brentq(Sigma, ks[i], ks[i + 1]) for i in range(len(ks) - 1)
      if sg[i] * sg[i + 1] < 0]

print('=' * 76)
print('A_1, A_2 = A_3, A_4 versus mu_eff/U_eff   (thesis figure)')
print('=' * 76)
print('a = %.6f /MHz, b = %.6f /MHz' % (A_, B_))
print('pole  x = b/3 = %.6f /MHz;  zero x = a = %.6f /MHz;  limit -1/3'
      % (B_ / 3, A_))
print('x -> 0 gives mu/U = a/(-b) = %.6f' % (A_ / -B_))
print()
print('For Sigma > 0 the reachable set is (-inf, -3/4] U (-1/3, +inf):')
print('   the window -3/4 < mu_eff/U_eff < -1/3 is NOT reachable.')
print('For Sigma < 0 that window IS swept, since x < 0 then runs 0 -> -inf.')
print('Sigma changes sign at kappa = %s' % ', '.join('%.4f' % z for z in kz))
print('   -> Sigma(0.8) = %+.6f (negative), Sigma(1.2) = %+.6f (positive)'
      % (Sigma(0.8), Sigma(1.2)))
print()
print('operating point Om_a = Om_b = 1, kappa = %.1f:' % KMAIN)
print('   mu_eff/U_eff = %+.6f' % muU_x(Sigma(KMAIN)))
print('   A_1 = %.6f, A_2 = A_3 = %.6f, A_4 = %.6f MHz'
      % (K1, K2 * abs(jv(N_RES, KMAIN)), abs(Sigma4(KMAIN))))
print('   ratios to A_1: 1, %.4f, %.4f   (targets 1, %.4f, 2)'
      % (K2 * abs(jv(N_RES, KMAIN)) / K1, abs(Sigma4(KMAIN)) / K1, math.sqrt(2)))
z1 = jn_zeros(N_RES, 1)[0]
k4 = [brentq(Sigma4, ks[i], ks[i + 1]) for i in range(len(ks) - 1)
      if Sigma4(ks[i]) * Sigma4(ks[i + 1]) < 0]
print()
print('zeros that switch a channel off: A_2 = A_3 = 0 at kappa = %.4f '
      '(first zero of J_%d);' % (z1, N_RES))
print('   A_4 = 0 at kappa = %s (zeros of Sigma_4)'
      % ', '.join('%.4f' % z for z in k4[:3]))
print('=' * 76)

# ---- figure ---------------------------------------------------------
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axs = plt.subplots(1, 3, figsize=(15.0, 4.7))


def brs(m, jump=1.0):
  cut = np.where(np.abs(np.diff(m)) > jump)[0]
  return [s for s in np.split(np.arange(m.size), cut + 1) if s.size > 1]


# (a) A_1, sweeping Om_a at Om_b = 1
ax = axs[0]
for k in KAPPAS:
  oa = np.geomspace(0.05, 20.0, 40000)
  m = np.array([muU_x((1.0 / oa[i]**2) * Sigma(k)) for i in range(oa.size)])
  y = oa**2 * K1
  f = True
  for s in brs(m):
    ax.plot(m[s], y[s], color=COL[k], lw=1.8,
            label=(r'$\kappa=%.1f$' % k) if f else None)
    f = False
ax.plot([muU_x(Sigma(KMAIN))], [K1], '*', color='#c2185b', ms=15, zorder=5,
        markeredgecolor='w', markeredgewidth=0.7)
ax.set_ylabel(r'$A_1$   [MHz]')
ax.set_title(r'(a) $A_1\propto\Omega_\alpha^2$: sweep $\Omega_\alpha$ '
             r'at $\Omega_\beta=1$', fontsize=9.5)

# (b) A_2 = A_3 and (c) A_4, sweeping Om_b at Om_a = 1
for idx, (lab, fun, ttl) in enumerate([
    (r'$A_2=A_3$   [MHz]',
     lambda ob, k: ob * K2 * abs(jv(N_RES, k)),
     r'(b) $A_2=A_3\propto\Omega_\alpha\Omega_\beta$: sweep $\Omega_\beta$'),
    (r'$A_4$   [MHz]',
     lambda ob, k: ob**2 * abs(Sigma4(k)),
     r'(c) $A_4\propto\Omega_\beta^2$: sweep $\Omega_\beta$')]):
  ax = axs[idx + 1]
  for k in KAPPAS:
    ob = np.geomspace(0.03, 20.0, 40000)
    m = np.array([muU_x(ob[i]**2 * Sigma(k)) for i in range(ob.size)])
    y = fun(ob, k)
    f = True
    for s in brs(m):
      ax.plot(m[s], y[s], color=COL[k], lw=1.8,
              label=(r'$\kappa=%.1f$' % k) if f else None)
      f = False
  y0 = fun(1.0, KMAIN)
  ax.plot([muU_x(Sigma(KMAIN))], [y0], '*', color='#c2185b', ms=15, zorder=5,
          markeredgecolor='w', markeredgewidth=0.7)
  ax.set_ylabel(lab)
  ax.set_title(ttl + r' at $\Omega_\alpha=1$', fontsize=9.5)

for ax in axs:
  ax.axvline(-1.0 / 3.0, color='gray', ls=':', lw=1.1)
  ax.axvline(-0.75, color='gray', ls=':', lw=1.1)
  ax.axvspan(-0.75, -1.0 / 3.0, color='0.9', zorder=0)
  ax.set_xlim(-4, 4)
  ax.set_yscale('log')
  ax.set_xlabel(r'$\mu_{\mathrm{eff}}/U_{\mathrm{eff}}$')
  ax.legend(fontsize=8, loc='lower right')
  ax.grid(True, ls='--', alpha=0.35)
axs[0].set_ylim(1e-4, 2e1)
axs[1].set_ylim(1e-4, 1e0)
axs[2].set_ylim(1e-5, 1e1)

fig.tight_layout()
fig.savefig('figures/amplitudes_vs_muU.png', dpi=200)
print('Saved -> figures/amplitudes_vs_muU.png')
