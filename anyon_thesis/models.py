"""Numerical kernels for the thesis notebooks; no calculations run on import.

Effective-model parameters retain the original scripts' numerical rate scale.
They are NOT automatically interpreted as cycles/us or rad/us. Evolution uses
dimensionless K and tau, where H = hbar * omega_ref * K, tau = omega_ref * t.
The target AHM uses H/U with a separately declared positive reference U.
"""
from dataclasses import dataclass, replace
from itertools import product
import numpy as np
from scipy.special import jv
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


@dataclass(frozen=True)
class EffectiveParameters:
    """Locked-resonance, equal-site-detuning family in the archived scripts.

    All dimensional inputs use the same native rate unit. match_11_diagonal is
    an explicitly imposed diagonal shift, not a self-consistent laser retuning.
    """
    C00: float = 0.0
    C01: float = 20.0
    C11: float = 15.0
    C02: float = 30.0
    C12: float = 45.0
    C22: float = 55.0
    omega_a: float = 1.0
    omega_b: float = 1.0
    drive_rate: float = 40.0
    modulation: float = 48.0
    resonance_order: int = 1
    sideband_cutoff: int = 60
    match_11_diagonal: bool = True

    @property
    def kappa(self):
        return self.modulation / self.drive_rate

    @property
    def detuning(self):
        return 2 * self.C11 - self.C02 - self.resonance_order * self.drive_rate


def _inverse(values):
    values = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(values)) or np.any(np.abs(values) < 1e-10):
        raise ValueError("A perturbative denominator is zero or non-finite.")
    return 1.0 / values


def effective_bond(p=EffectiveParameters()):
    """Return the 3x3 bond diagonal, four signed rates, and shift diagnostics.

    Matches heff2_E_vs_N.py and the equal-detuning case of heff2_states.py.
    The four directed channels are 01->10, 11->20, 02->11, 12->21.
    """
    if not np.isfinite(p.drive_rate) or p.drive_rate <= 0:
        raise ValueError("drive_rate must be positive and finite.")
    if not isinstance(p.resonance_order, (int, np.integer)):
        raise ValueError("resonance_order must be an integer.")
    if not isinstance(p.sideband_cutoff, (int, np.integer)) or p.sideband_cutoff < 0:
        raise ValueError("sideband_cutoff must be a nonnegative integer.")
    ell = np.arange(-p.sideband_cutoff, p.sideband_cutoff + 1)
    weights = jv(ell, p.kappa) ** 2
    d1, d2, d3 = _inverse([p.C01 - 2*p.C00, 2*p.C11 - p.C01, p.C12 - p.C02])
    offsets = [p.C02-p.C01, p.C12-2*p.C11, 2*p.C22-p.C12]
    denominators = [offset + p.detuning + ell*p.drive_rate for offset in offsets]
    s1, s2, s3 = [float(weights @ _inverse(den)) for den in denominators]
    a, b = p.omega_a**2, p.omega_b**2
    diagonal = np.array([
        [-2*a*d1, a*(d1-d2)-b*s1, b*s1-a*d3],
        [a*(d1-d2)-b*s1, 2*a*d2-2*b*s2, a*d3+b*(s2-s3)],
        [b*s1-a*d3, a*d3+b*(s2-s3), 2*b*s3],
    ])
    shift = float(diagonal[2, 0] - diagonal[1, 1]) if p.match_11_diagonal else 0.0
    diagonal[1, 1] += shift
    cross = p.omega_a*p.omega_b*(d2-d3)*jv(p.resonance_order, p.kappa)
    rates = np.array([a*(d1-d2), cross, cross, b*(s2-s3)])
    minimum = min(abs(1/d1), abs(1/d2), abs(1/d3),
                  *(float(np.min(np.abs(den))) for den in denominators))
    diagnostics = {"imposed_11_shift": shift, "minimum_denominator": minimum,
                   "conservative_small_parameter": max(abs(p.omega_a), abs(p.omega_b))/minimum,
                   "retained_bessel_weight": float(weights.sum())}
    return diagonal, rates, diagnostics


def phase_pattern(mode="zero", angle=0.0, alpha=0.0, beta=0.0):
    """Four phases, with an explicitly mathematical status for imposed patterns.

    target: (0, phi, 0, phi). legacy_drive: archived uniform-bond assignment
    (0, theta, -2 theta, -theta); its physical reachability is not established.
    laser_pair: the two-site alpha,beta assignment. zero: all phases zero.
    """
    if mode == "zero":
        return np.zeros(4)
    if mode == "target":
        return np.array([0.0, angle, 0.0, angle])
    if mode == "legacy_drive":
        return np.array([0.0, angle, -2*angle, -angle])
    if mode == "laser_pair":
        return np.array([-alpha, beta-2*alpha, alpha-2*beta, -beta])
    raise ValueError(f"Unknown phase mode: {mode}")


def occupation_basis(length, total_number=None):
    """Lexicographic occupations in {0,1,2}; generate fixed-N sectors directly."""
    if not isinstance(length, (int, np.integer)) or length < 1:
        raise ValueError("length must be a positive integer.")
    if total_number is None:
        if 3**length > 100_000:
            raise ValueError("Full basis is too large; select a particle-number sector.")
        return np.array(list(product(range(3), repeat=length)), dtype=np.int8)
    if not isinstance(total_number, (int, np.integer)) or not 0 <= total_number <= 2*length:
        raise ValueError("total_number must be an integer between 0 and 2*length.")
    rows = []
    def build(prefix, remaining, sites):
        if sites == 0:
            if remaining == 0:
                rows.append(prefix)
            return
        for n in range(3):
            if 0 <= remaining-n <= 2*(sites-1):
                build(prefix + (n,), remaining-n, sites-1)
    build((), total_number, length)
    return np.array(rows, dtype=np.int8).reshape(-1, length)


def _bonds(length, periodic):
    if periodic and length < 3:
        raise ValueError("Periodic chains require at least three sites in this convention.")
    pairs = [(j, j+1) for j in range(length-1)]
    return pairs + ([(length-1, 0)] if periodic else [])


CHANNELS = (((0, 1), (1, 0)), ((1, 1), (2, 0)),
            ((0, 2), (1, 1)), ((1, 2), (2, 1)))


def effective_chain(length, diagonal, rates, phases=None, total_number=None, periodic=False):
    """Dense effective model: sum each bond diagonal once and add h.c. once."""
    basis = occupation_basis(length, total_number)
    if len(basis) > 2500:
        raise ValueError("Dense demonstration limited to 2500 states; use smaller sectors.")
    diagonal, rates = np.asarray(diagonal), np.asarray(rates)
    phases = np.zeros(4) if phases is None else np.asarray(phases)
    if diagonal.shape != (3, 3) or rates.shape != (4,) or phases.shape != (4,):
        raise ValueError("Expected diagonal (3,3), rates (4,), phases (4,).")
    if not np.allclose(diagonal.imag, 0):
        raise ValueError("The bond diagonal must be real.")
    locations = {tuple(row): i for i, row in enumerate(basis)}
    H = np.zeros((len(basis), len(basis)), dtype=complex)
    for left, right in _bonds(length, periodic):
        H[np.diag_indices(len(basis))] += diagonal[basis[:, left], basis[:, right]]
        for initial, final, rate, phase in zip(
                (c[0] for c in CHANNELS), (c[1] for c in CHANNELS), rates, phases):
            for column, row in enumerate(basis):
                if (row[left], row[right]) == initial:
                    target = row.copy()
                    target[left], target[right] = final
                    index = locations.get(tuple(target))
                    if index is not None:
                        amplitude = rate*np.exp(1j*phase)
                        H[index, column] += amplitude
                        H[column, index] += np.conj(amplitude)
    return H, basis


def evolve_probabilities(K, initial_state, tau):
    """Unitary evolution with dimensionless Hermitian K and dimensionless tau."""
    K, state, tau = np.asarray(K), np.asarray(initial_state, complex), np.asarray(tau)
    if K.shape != (len(state), len(state)) or not np.allclose(K, K.conj().T):
        raise ValueError("K must be a Hermitian matrix matching the state dimension.")
    if not np.isclose(np.vdot(state, state).real, 1.0):
        raise ValueError("The initial state must be normalized.")
    if tau.ndim != 1 or not np.all(np.isfinite(tau)):
        raise ValueError("tau must be a finite one-dimensional array.")
    levels, vectors = np.linalg.eigh(K)
    coefficients = vectors.conj().T @ state
    evolved = (np.exp(-1j*np.outer(tau, levels))*coefficients) @ vectors.T
    return abs(evolved)**2


def basis_state(basis, occupation):
    matches = np.flatnonzero(np.all(basis == np.asarray(occupation), axis=1))
    if len(matches) != 1:
        raise ValueError("The requested occupation is not in the chosen basis.")
    state = np.zeros(len(basis), complex)
    state[matches[0]] = 1.0
    return state


def parameter_curves(kappa, r, p=EffectiveParameters()):
    """Archived effective-parameter formula family, before diagonal compensation.

    r = omega_b^2/omega_a^2 >= 0. The unnormalized mu, U, and |g| values
    are in native rate units. Ratios are undefined at U=0 and return NaN.
    """
    r = np.asarray(r, float)
    if np.any(r < 0) or not np.all(np.isfinite(r)):
        raise ValueError("r must be finite and nonnegative.")
    if p.omega_a == 0:
        raise ValueError("omega_a must be nonzero to define r.")
    ell = np.arange(-p.sideband_cutoff, p.sideband_cutoff+1)
    w = jv(ell, kappa)**2
    d1, d2, d3 = _inverse([p.C01-2*p.C00, 2*p.C11-p.C01, p.C12-p.C02])
    sigma = float(w @ _inverse(p.C02-p.C01+p.detuning+ell*p.drive_rate))
    sigma4 = float(w @ (_inverse(p.C12-2*p.C11+p.detuning+ell*p.drive_rate)
                        - _inverse(2*p.C22-p.C12+p.detuning+ell*p.drive_rate)))
    a, b = 3*d1-d2, 4*d1+d3-2*d2
    x = r*sigma
    mu = p.omega_a**2*(a-x)
    interaction = p.omega_a**2*(3*x-b)
    amplitudes = p.omega_a**2*np.array([
        np.full_like(r, abs(d1-d2)),
        np.sqrt(r)*abs(d2-d3)*abs(jv(p.resonance_order, kappa)),
        np.sqrt(r)*abs(d2-d3)*abs(jv(p.resonance_order, kappa)),
        r*abs(sigma4),
    ])
    valid = abs(interaction) > 1e-12
    ratio = np.divide(mu, interaction, out=np.full_like(mu, np.nan), where=valid)
    normalized = np.divide(amplitudes, interaction,
                           out=np.full_like(amplitudes, np.nan), where=valid)
    return {"mu": mu, "U": interaction, "amplitudes": amplitudes,
            "mu_over_U": ratio, "amplitudes_over_U": normalized,
            "sigma": sigma, "a": float(a), "b": float(b), "x": x}


def target_hamiltonian(length, J, U, theta, total_number=None, periodic=False):
    r"""Sparse truncated AHM, including exact operator ordering.

    H = -J sum_j [b_j^dag exp(i theta n_j) b_{j+1} + h.c.]
        + U/2 sum_j n_j(n_j-1).
    J,U are real energies (or both dimensionless in a common reference unit).
    The full-space result does not impose particle number or chemical potential.
    """
    if not all(np.isfinite(v) for v in (J, U, theta)):
        raise ValueError("J, U and theta must be finite real scalars.")
    basis = occupation_basis(length, total_number)
    locations = {tuple(row): i for i, row in enumerate(basis)}
    rows, cols, values = [], [], []
    for column, cfg in enumerate(basis):
        rows.append(column); cols.append(column)
        values.append(float(U*np.sum(cfg*(cfg-1))/2))
        for left, right in _bonds(length, periodic):
            n, m = int(cfg[left]), int(cfg[right])
            if n < 2 and m > 0:
                target = cfg.copy()
                target[left], target[right] = n+1, m-1
                index = locations[tuple(target)]
                amplitude = -J*np.sqrt((n+1)*m)*np.exp(1j*theta*n)
                rows.extend([index, column]); cols.extend([column, index])
                values.extend([amplitude, np.conj(amplitude)])
    return coo_matrix((values, (rows, cols)), shape=(len(basis), len(basis)),
                      dtype=complex).tocsr(), basis


def lowest_levels(H, count=5):
    """Small dense solve or deterministic sparse solve, with residuals."""
    n = H.shape[0]
    if H.shape != (n, n) or not 1 <= count <= n:
        raise ValueError("count must be between one and the matrix dimension.")
    if n <= 128 or count >= n-1:
        dense = H.toarray() if hasattr(H, "toarray") else np.asarray(H)
        values, vectors = np.linalg.eigh(dense)
        values, vectors = values[:count], vectors[:, :count]
    else:
        v0 = np.random.default_rng(1234).standard_normal(n)
        values, vectors = eigsh(H, k=count, which="SA", v0=v0, tol=1e-10,
                                ncv=min(n, max(2*count+1, 40)))
        order = np.argsort(values)
        values, vectors = values[order], vectors[:, order]
    residuals = np.linalg.norm(H @ vectors - vectors*values, axis=0)
    return values, residuals


def gap_grid(length, ratios, angles, count=3, total_number=None):
    """Bounded teaching sweep with U=1; no multiprocessing or implicit files."""
    ratios, angles = np.asarray(ratios), np.asarray(angles)
    if len(ratios)*len(angles) > 100:
        raise ValueError("Teaching sweep limited to 100 points; use the archived sweep explicitly.")
    basis = occupation_basis(length, total_number)
    if len(basis) > 1000:
        raise ValueError("Teaching sweep limited to 1000 basis states.")
    values = np.empty((len(ratios), len(angles), count))
    largest_residual = 0.0
    for i, ratio in enumerate(ratios):
        for j, angle in enumerate(angles):
            H, _ = target_hamiltonian(length, ratio, 1.0, angle, total_number)
            values[i, j], residuals = lowest_levels(H, count)
            largest_residual = max(largest_residual, float(residuals.max()))
    return values, largest_residual
