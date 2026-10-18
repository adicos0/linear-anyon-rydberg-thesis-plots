"""Physics and legacy-equivalence checks, independent of notebook rendering."""
from dataclasses import replace
from pathlib import Path
import importlib.util
import numpy as np
import pytest
from anyon_thesis.models import (
    EffectiveParameters, effective_bond, effective_chain, occupation_basis,
    phase_pattern, basis_state, evolve_probabilities, parameter_curves,
    target_hamiltonian, lowest_levels, gap_grid,
)

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("length", [2, 3, 4])
def test_sectors_partition_full_hilbert_space(length):
    sectors = [occupation_basis(length, n) for n in range(2*length+1)]
    assert sum(map(len, sectors)) == 3**length
    for n, sector in enumerate(sectors):
        assert np.all(sector.sum(axis=1) == n)


def test_effective_rates_match_archived_source():
    # Read only the original settings/functions; do not launch its L=10 solve.
    original = (ROOT / "scripts/heff2_E_vs_N.py").read_text()
    scope = {}
    exec(compile(original.split("# 2.  SECTOR HAMILTONIAN")[0],
                 "archived_effective_prefix", "exec"), scope)
    for rate in [35.0, 40.0, 47.0]:
        diagonal, rates, _ = effective_bond(replace(EffectiveParameters(), drive_rate=rate))
        old_diagonal, old_rates = scope["heff2_couplings"](rate)
        np.testing.assert_allclose(diagonal, old_diagonal, atol=1e-13)
        np.testing.assert_allclose(rates, old_rates, atol=1e-13)


def test_small_effective_sectors_match_saved_thesis_spectrum():
    p = EffectiveParameters()
    diagonal, rates, _ = effective_bond(p)
    with np.load(ROOT / "data/heff2_E_vs_N_L10.npz", allow_pickle=False) as data:
        assert int(data["L"]) == 10
        for number in [0, 1, 2]:
            H, _ = effective_chain(10, diagonal, rates, total_number=number)
            expected = np.sort(data["E"][data["NLAB"] == number])
            np.testing.assert_allclose(np.linalg.eigvalsh(H), expected, atol=2e-12)


def test_sideband_convergence_and_explicit_compensation():
    p = EffectiveParameters()
    d40, g40, _ = effective_bond(replace(p, sideband_cutoff=40))
    d80, g80, info = effective_bond(replace(p, sideband_cutoff=80))
    np.testing.assert_allclose(d40, d80, atol=1e-13)
    np.testing.assert_allclose(g40, g80, atol=1e-13)
    bare, _, _ = effective_bond(replace(p, match_11_diagonal=False))
    assert np.isclose(d80[1, 1], d80[2, 0])
    assert np.isclose(info["imposed_11_shift"], d80[1, 1]-bare[1, 1])


def test_effective_number_conservation_and_sector_spectra():
    diagonal, rates, _ = effective_bond()
    H, basis = effective_chain(3, diagonal, rates, phase_pattern("target", 0.73))
    N = np.diag(basis.sum(axis=1))
    np.testing.assert_allclose(H, H.conj().T, atol=1e-14)
    np.testing.assert_allclose(H @ N, N @ H, atol=1e-14)
    sectors = [np.linalg.eigvalsh(effective_chain(
        3, diagonal, rates, phase_pattern("target", 0.73), n)[0]) for n in range(7)]
    np.testing.assert_allclose(np.linalg.eigvalsh(H), np.sort(np.concatenate(sectors)), atol=1e-13)


def test_two_level_solution_and_unitarity():
    coupling, detuning = 0.4+0.3j, 0.2
    H = np.array([[detuning, coupling], [coupling.conjugate(), -detuning]])
    tau = np.linspace(0, 10, 61)
    P = evolve_probabilities(H, [1, 0], tau)
    frequency = np.sqrt(abs(coupling)**2+detuning**2)
    exact = abs(coupling)**2/frequency**2*np.sin(frequency*tau)**2
    np.testing.assert_allclose(P[:, 1], exact, atol=1e-13)
    np.testing.assert_allclose(P.sum(axis=1), 1, atol=1e-13)
    np.testing.assert_allclose(evolve_probabilities(np.zeros((2, 2)), [1, 0], tau),
                               np.tile([1, 0], (len(tau), 1)), atol=1e-13)


def test_assumed_phase_gauge_identity_not_physical_reachability():
    angle = 0.73
    diagonal, rates, _ = effective_bond()
    H, basis = effective_chain(4, diagonal, rates, phase_pattern("legacy_drive", angle), 2)
    target, _ = effective_chain(4, diagonal, rates, phase_pattern("target", -angle), 2)
    phase = np.exp(2j*angle*np.sum(basis == 2, axis=1))
    np.testing.assert_allclose(phase.conj()[:, None]*H*phase[None, :], target, atol=1e-13)


def test_density_conservation_and_phase_periodicity():
    diagonal, rates, _ = effective_bond()
    tau = np.linspace(0, 30, 31)
    results = []
    for angle in [0.0, 2*np.pi]:
        H, basis = effective_chain(5, diagonal, rates, phase_pattern("target", angle), 2)
        P = evolve_probabilities(H, basis_state(basis, [0, 0, 2, 0, 0]), tau)
        density = P @ basis
        np.testing.assert_allclose(density.sum(axis=1), 2, atol=1e-12)
        results.append(density)
    np.testing.assert_allclose(*results, atol=1e-12)


def test_parameter_scaling_and_unreachable_ratio_domain():
    r = np.linspace(0, 4, 11)
    p = EffectiveParameters()
    a = parameter_curves(1.2, r, p)
    b = parameter_curves(1.2, r, replace(p, omega_a=2.0))
    for key in ["mu", "U", "amplitudes"]:
        np.testing.assert_allclose(b[key], 4*a[key])
    np.testing.assert_allclose(a["mu_over_U"], b["mu_over_U"])
    np.testing.assert_allclose(a["amplitudes"][1], a["amplitudes"][2])
    with pytest.raises(ValueError):
        parameter_curves(1.2, [-1])


def test_parameter_family_matches_archived_analytic_functions():
    original = (ROOT / "scripts/amplitudes_vs_muU_thesis.py").read_text()
    scope = {}
    exec(compile(original.split("# sign change of Sigma")[0],
                 "archived_parameter_prefix", "exec"), scope)
    r = np.array([0.0, 0.2, 1.0, 4.0])
    for kappa in [0.8, 1.2, 2.0]:
        curve = parameter_curves(kappa, r)
        sigma, sigma4 = scope["Sigma"](kappa), scope["Sigma4"](kappa)
        np.testing.assert_allclose(curve["sigma"], sigma, atol=1e-14)
        np.testing.assert_allclose(curve["mu_over_U"], scope["muU_x"](r*sigma), atol=1e-12)
        np.testing.assert_allclose(curve["amplitudes"][3], r*abs(sigma4), atol=1e-14)
        for j, ratio in enumerate(r):
            _, rates, _ = effective_bond(replace(EffectiveParameters(),
                modulation=kappa*40, omega_b=np.sqrt(ratio)))
            np.testing.assert_allclose(curve["amplitudes"][:, j], abs(rates), atol=1e-13)


@pytest.mark.parametrize("periodic", [False, True])
def test_target_operator_matches_original_tensor_contraction(periodic):
    spec = importlib.util.spec_from_file_location("archived_matvec", ROOT / "scripts/doApplyHam.py")
    legacy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(legacy)
    J, U, theta, length = 0.31, 1.0, 0.73, 3
    b = np.diag(np.sqrt([1.0, 2.0]), 1)
    number = b.conj().T @ b
    phase = np.diag(np.exp(1j*theta*np.diag(number)))
    hloc = -J*(np.kron(b.conj().T @ phase, b)
               + np.kron(phase.conj().T @ b, b.conj().T))
    H, basis = target_hamiltonian(length, J, U, theta, periodic=periodic)
    onsite = U*np.sum(basis*(basis-1), axis=1)/2
    eye = np.eye(len(basis))
    reference = np.column_stack([legacy.doApplyHam(v, hloc.reshape(3,3,3,3), length, periodic)
                                  + onsite*v for v in eye.T])
    np.testing.assert_allclose(H.toarray(), reference, atol=1e-13)


def test_target_limits_sector_selection_and_sparse_residuals():
    H, basis = target_hamiltonian(4, 0, 1, 0.73)
    expected = np.sum(basis*(basis-1), axis=1)/2
    np.testing.assert_allclose(H.toarray(), np.diag(expected))
    for theta in [0.0, 0.73, 2*np.pi]:
        full, basis = target_hamiltonian(4, 0.2, 1.0, theta)
        number = np.diag(basis.sum(axis=1))
        np.testing.assert_allclose(full.toarray() @ number, number @ full.toarray(), atol=1e-13)
        subset, sector = target_hamiltonian(4, 0.2, 1.0, theta, total_number=3)
        selected = np.flatnonzero(basis.sum(axis=1) == 3)
        np.testing.assert_allclose(subset.toarray(), full.toarray()[np.ix_(selected, selected)])
    a, _ = target_hamiltonian(5, 0.3, 1, 0)
    b, _ = target_hamiltonian(5, 0.3, 1, 2*np.pi)
    np.testing.assert_allclose(a.toarray(), b.toarray(), atol=1e-13)
    values, residuals = lowest_levels(a, 4)
    assert np.max(residuals) < 1e-8
    np.testing.assert_allclose(values, np.linalg.eigvalsh(a.toarray())[:4], atol=1e-9)


def test_guards_and_cached_grid_shape():
    with pytest.raises(ValueError):
        effective_bond(replace(EffectiveParameters(), C01=0))
    with pytest.raises(ValueError):
        effective_chain(8, np.zeros((3,3)), np.zeros(4))
    with pytest.raises(ValueError):
        evolve_probabilities(np.eye(2), [1, 1], [0])
    with pytest.raises(ValueError):
        gap_grid(10, [0.1], [0.0])
    with np.load(ROOT / "data/anyonHubbardGapsMulti_N10.npz", allow_pickle=False) as z:
        assert z["energies"].shape == (len(z["ratios"]), len(z["thetas"]), 5)
        assert np.all(np.diff(z["energies"], axis=-1) >= -1e-10)
