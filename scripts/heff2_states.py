"""
Population dynamics under H_eff^(2), starting from each transition state
|01>, |10>, |11>, |12>, |21|  (and the dark |22> / |00> for reference).

Parameters: realistic Rb Rydberg scale, dispersive regime, with the thesis
resonance condition 2 C11 = C02 imposed.  The Floquet detunings nu_tilde_j
are chosen so that:
   * the Floquet-assisted orders are exact integers  (N0=1, N1=2, q=-1),
     activating the |11><->|20>, |02><->|11>, |12><->|21> channels, AND
   * every beta-transition sideband  C0x-C0y+nu_tilde_j+l*omega'  stays
     >= 10 MHz (= 10 Omega) off zero, so NO zero-frequency term appears
     (H^(1)_eff = 0 holds) and every Stark denominator is finite.

Om_alpha is tuned so the AC-Stark-shifted |11> and |20> are degenerate,
bringing the Floquet 3-level block cleanly onto resonance (a physical knob;
Om_alpha stays at the ~MHz Rb scale).

Units: hbar = 1; energies in MHz (= rad/us), time in us.  P_n=|<n|e^{-iH t}|psi0>|^2.
"""

import numpy as np
from scipy.special import jv
import matplotlib.pyplot as plt

# =====================================================================
# 1.  PARAMETERS
# =====================================================================
hbar = 1.0
C00, C01, C11 = 0.0, 20.0, 15.0
C02 = 2.0 * C11                      # 2 C11 = C02  ->  30
C12, C22 = 45.0, 55.0
_C = {(0,0):C00,(0,1):C01,(1,1):C11,(0,2):C02,(1,2):C12,(2,2):C22}
def C(i, j): return _C[(min(i,j), max(i,j))]

Om_a  = 1.0                          # realistic Rb Rabi frequencies (MHz)
Om_b  = 1.0
alpha = np.pi/4                      # alpha-drive phase
beta  = np.pi/3                      # beta-drive phase
# Put the Floquet 3-level block {02,11,20} on exact resonance by setting the
# two-photon Floquet detuning to compensate its differential AC-Stark shift
# (E_11 -> E_20).  A standard experimental knob; nu_tilde is adjusted slightly
# off the bare integer to absorb the Stark shift.
COMPENSATE_11BLOCK = True

wp = 40.0                            # omega'  (large -> well-separated sidebands)
A  = 48.0
kappa = A / wp                       # 1.2
# nu_tilde_0 = nu_tilde_1 makes the exchange {01,10} and cross-site {12,21}
# blocks EXACTLY resonant (E_01=E_10, E_12=E_21, E_02=E_20 automatically),
# and gives N0=N1=1 (both Floquet channels on), q=0.  Nearest sideband 10 MHz.
nt0, nt1 = -40.0, -40.0

LMAX = 40
ells = np.arange(-LMAX, LMAX + 1)
Jl   = jv(ells, kappa)

N0 = (2*C11 - C02 - nt0)/wp
N1 = (2*C11 - C02 - nt1)/wp
q  = (nt1 - nt0)/wp
def is_int(x, tol=1e-9): return abs(x - round(x)) < tol

labels = ['00','01','02','10','11','12','20','21','22']
idx = {s:i for i,s in enumerate(labels)}
DIM = 9

# =====================================================================
# 2.  DIAGONAL BLOCK (AC Stark shifts) -- all denominators finite here
# =====================================================================
def stark_sum(dconst, nu):
    denom = dconst + nu + ells*wp
    if np.any(np.abs(denom) < 1e-9):
        raise ZeroDivisionError(f"resonant Stark denominator (dconst={dconst}, nu={nu})")
    return np.sum(Jl**2 / denom)

E = {}
E['00'] = -2*Om_a**2/(C01-2*C00)
E['01'] = ( Om_a**2/(C01-2*C00) - Om_a**2/(2*C11-C01)
          - Om_b**2*stark_sum(C02-C01, nt1) )
E['02'] = ( Om_b**2*stark_sum(C02-C01, nt1) - Om_a**2/(C12-C02) )
E['10'] = ( Om_a**2/(C01-2*C00) - Om_a**2/(2*C11-C01)
          - Om_b**2*stark_sum(C02-C01, nt0) )
E['11'] = ( 2*Om_a**2/(2*C11-C01)
          - Om_b**2*stark_sum(C12-2*C11, nt0) - Om_b**2*stark_sum(C12-2*C11, nt1) )
E['12'] = ( Om_a**2/(C12-C02) + Om_b**2*stark_sum(C12-2*C11, nt1)
          - Om_b**2*stark_sum(2*C22-C12, nt0) )
E['20'] = ( Om_b**2*stark_sum(C02-C01, nt0) - Om_a**2/(C12-C02) )
E['21'] = ( Om_a**2/(C12-C02) + Om_b**2*stark_sum(C12-2*C11, nt0)
          - Om_b**2*stark_sum(2*C22-C12, nt1) )
E['22'] = ( Om_b**2*stark_sum(2*C22-C12, nt0) + Om_b**2*stark_sum(2*C22-C12, nt1) )

# tune the 3-level block onto exact resonance (E_11 -> E_20 = E_02)
delta_comp = (E['20'] - E['11']) if COMPENSATE_11BLOCK else 0.0
E['11'] += delta_comp

# =====================================================================
# 3.  9x9 EFFECTIVE HAMILTONIAN
# =====================================================================
H = np.zeros((DIM, DIM), dtype=complex)
for s in labels: H[idx[s], idx[s]] = E[s]
def add_coupling(bra, ket, amp):
    H[idx[bra], idx[ket]] += amp
    H[idx[ket], idx[bra]] += np.conj(amp)

g1 = Om_a**2*(1.0/(C01-2*C00) - 1.0/(2*C11-C01))            # |01><->|10>
add_coupling('10', '01', g1*np.exp(-1j*alpha))

pref = Om_a*Om_b*(1.0/(2*C11-C01) - 1.0/(C12-C02))
g2 = pref*jv(round(N0), kappa) if is_int(N0) else 0.0       # |11><->|20>
if g2: add_coupling('20', '11', g2*np.exp(1j*(beta-2*alpha)))
g3 = pref*jv(round(N1), kappa) if is_int(N1) else 0.0       # |02><->|11>
if g3: add_coupling('11', '02', g3*np.exp(1j*(alpha-2*beta)))
if is_int(q):                                               # |12><->|21>
    qi = round(q)
    Sig = (1.0/(C12-2*C11+nt0+ells*wp) - 1.0/(2*C22-C12+nt0+ells*wp))
    g4 = Om_b**2*np.sum(Jl*jv(ells-qi, kappa)*Sig)
    add_coupling('21', '12', g4*np.exp(-1j*beta))
else:
    g4 = 0.0
assert np.allclose(H, H.conj().T)

# =====================================================================
# 4.  EVOLUTION
# =====================================================================
evals, evecs = np.linalg.eigh(H)
def evolve(psi0, t):
    c = evecs.conj().T @ psi0
    return np.abs((np.exp(-1j*np.outer(t, evals)/hbar)*c) @ evecs.T)**2
def ket(s):
    v = np.zeros(DIM, dtype=complex); v[idx[s]] = 1.0; return v
def two_level(Ea, Eb, g):
    d = 0.5*(Ea-Eb); Om = np.hypot(abs(g), d)
    return (abs(g)**2/(abs(g)**2+d**2) if Om>0 else 0.0), (np.pi/Om if Om>0 else np.inf)

# =====================================================================
# 5.  REPORT
# =====================================================================
print("="*68)
print("H_eff^(2) with 2 C11 = C02 imposed  (realistic Rb, dispersive)")
print("="*68)
print(f"C: 00={C00} 01={C01} 11={C11} 02={C02} 12={C12} 22={C22}")
print(f"Om_a={Om_a:.4f}  Om_b={Om_b:.4f}  kappa={kappa:.3f}  omega'={wp}  "
      f"nt0={nt0} nt1={nt1}")
print(f"phases: alpha = {alpha:.6f} rad (= pi/4),  beta = {beta:.6f} rad (= pi/3)")
print(f"N0={N0:.4f} N1={N1:.4f} q={q:.4f}   "
      f"g1={g1:+.4f} g2={g2:+.4f} g3={g3:+.4f} g4={g4:+.4f}")
print(f"3-level compensation: delta_comp = E_20 - E_11(bare) = {delta_comp:+.5f} MHz "
      f"(absorbed into nu-tilde)")
print("-"*68)
for s in labels:
    print(f"   E_{s} = {E[s]:+.5f}")
print("-"*68)
p_ex,T_ex = two_level(E['01'],E['10'],g1)
p_f, T_f  = two_level(E['11'],E['20'],g2)
p_cs,T_cs = two_level(E['12'],E['21'],g4)
print(f"exchange   |01><->|10|: max {p_ex:.3f}, T={T_ex:.1f}")
print(f"Floquet    |11><->|20|: max {p_f:.3f}, T={T_f:.1f}  (residual E11-E20={E['11']-E['20']:+.2e})")
print(f"cross-site |12><->|21|: max {p_cs:.3f}, T={T_cs:.1f}")
print("="*68)

# =====================================================================
# 6.  FIGURE  -- outline styled after the reference:
#     wide panels, suptitle listing every constant/phase,
#     |psi(0)> panel titles, P_nm(t)=|<nm|psi(t)>|^2 legends
# =====================================================================
plt.rcParams.update({'font.size': 12})

# per-panel palette in the reference style: initial state orange, partners blue/green
PALETTE = ['tab:orange', 'tab:blue', 'tab:green']

# (initial state, states to display, time window)
panels = [
    ('01', ['01','10'],            4*T_ex),
    ('10', ['10','01'],            4*T_ex),
    ('11', ['11','20','02'],       4*T_f),
    ('12', ['12','21'],            4*T_cs),
    ('21', ['21','12'],            4*T_cs),
    ('22', ['22'],                 4*T_cs),
]

fig, axes = plt.subplots(3, 2, figsize=(16, 13.5))
for k, (ax, (init, watch, tmax)) in enumerate(zip(axes.ravel(), panels)):
    tt = np.linspace(0, tmax, 5000)
    pops = evolve(ket(init), tt)
    for s, col in zip(watch, PALETTE):
        ax.plot(tt, pops[:, idx[s]], color=col, lw=2.5,
                label=rf'$P_{{{s}}}(t) = |\langle {s}|\psi(t)\rangle|^2$')
    ax.set_title(rf'$|\psi(0)\rangle = |{init}\rangle$', fontsize=14)
    ax.set_xlabel(r'Time  $t$  [$\mu$s]')
    if k % 2 == 0:
        ax.set_ylabel('Population')
    ax.set_xlim(0, tmax); ax.set_ylim(-0.03, 1.05)
    ax.legend(loc='upper right', fontsize=11, framealpha=0.9)
    ax.grid(True, ls='--', alpha=0.45)

from fractions import Fraction
def pi_str(x):
    """Render x as a nice multiple of pi (e.g. pi/3), else as a decimal."""
    if abs(x) < 1e-12:
        return "0"
    fr = Fraction(x / np.pi).limit_denominator(24)
    if abs(float(fr) * np.pi - x) > 1e-9:
        return f"{x:.3g}"
    n, d = fr.numerator, fr.denominator
    sign = "-" if n < 0 else ""
    num = r"\pi" if abs(n) == 1 else rf"{abs(n)}\pi"
    return f"{sign}{num}" if d == 1 else f"{sign}{num}/{d}"

sup = (r'Population dynamics under $H_{\mathrm{eff}}^{(2)}$'
       + '\n'
       + r'$C_{00}=%g,\; C_{01}=%g,\; C_{11}=%g,\; C_{02}=2C_{11}=%g,\; '
         r'C_{12}=%g,\; C_{22}=%g$  MHz'
         % (C00, C01, C11, C02, C12, C22)
       + '\n'
       + r"$\Omega_\alpha=%g,\; \Omega_\beta=%g,\; \omega'=%g,\; A=%g\;"
         r"(\kappa=A/\omega'=%g),\; \tilde{\nu}_0=\tilde{\nu}_1=%g,\; "
         r"\alpha=%s,\; \beta=%s$"
         % (Om_a, Om_b, wp, A, kappa, nt0, pi_str(alpha), pi_str(beta)))
fig.suptitle(sup, fontsize=15)
fig.tight_layout(rect=[0, 0, 1, 0.925])
fig.savefig('heff2_states.png', dpi=150)
print("Saved figure -> heff2_states.png")
