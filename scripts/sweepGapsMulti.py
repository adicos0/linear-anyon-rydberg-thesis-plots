"""
sweepGapsMulti.py
---------------------------------------------------------------------
Sweep of the 1D anyon-Hubbard model over the plane (statistical angle
theta, ratio J/U), computing the lowest 'numval' eigenvalues E_0..E_4
at every grid point with the same Lanczos setup as mainAnyonHubbard.py
(eigsh, which = 'SA', Hamiltonian applied through doApplyHam plus the
diagonal on-site interaction; H is never built as a matrix).

The eigenvalues are stored in anyonHubbardGapsMulti_N<L>.npz and the
heat maps of the energy differences (E_k - E_0 and E_{k+1} - E_k) are
then rendered as individual figures by replotGapsMulti.py.

The grid points are distributed over worker processes, since each one
is an independent Lanczos solve. Set OMP_NUM_THREADS=1 in the
environment so BLAS threading does not oversubscribe the cores.

Energies in units of U (U = 1, J = ratio * U).
"""

import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')

import multiprocessing as mp
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh
from timeit import default_timer as timer

from doApplyHam import doApplyHam
from replotGapsMulti import plot_all

# Simulation parameters
Nsites = 10                                 # lattice sites (dim = 3**N)
usePBC = False                              # open chain, as in the notes
numval = 5                                  # eigenvalues per grid point
U = 1.0                                     # energy unit
thetas = np.linspace(0.0, np.pi, 21)        # statistical angle
ratios = np.linspace(0.05, 1.0, 21)         # J/U ratio

# Lanczos controls. Near-degenerate low-lying levels (small J/U) converge
# slowly with the default Krylov subspace, so enlarge it; tol=1e-8 is far
# below the resolution of the plots.
ncv = 40
tol = 1e-8

# Local operators on the truncated space {|0>,|1>,|2>}, d = 3
d = 3
b = np.diag(np.sqrt(np.arange(1, d)), k=1)
bdag = b.conj().T
n_op = bdag @ b

# On-site interaction (U/2) n(n-1) = U|2><2|: diagonal, site-summed
u_site = (U / 2) * np.arange(d) * (np.arange(d) - 1)
Udiag = np.zeros(d**Nsites)
for k in range(Nsites):
  Udiag += np.kron(np.ones(d**k),
                   np.kron(u_site, np.ones(d**(Nsites - 1 - k))))


def lowest_levels(J, theta):
  """Lowest 'numval' eigenvalues of H(J, U, theta) via Lanczos."""
  phase = np.diag(np.exp(1j * theta * np.diag(n_op)))
  hloc = (-J * (np.kron(bdag @ phase, b)
                + np.kron(phase.conj().T @ b, bdag))).reshape(d, d, d, d)
  H = LinearOperator((d**Nsites, d**Nsites),
                     matvec=lambda v: doApplyHam(v, hloc, Nsites, usePBC)
                     + Udiag * v,
                     dtype=np.complex128)
  # fixed starting vector so the sweep is reproducible run to run
  v0 = np.random.default_rng(1234).standard_normal(d**Nsites)
  return np.sort(eigsh(H, k=numval, which='SA', return_eigenvectors=False,
                       ncv=ncv, tol=tol, v0=v0))


def _work(args):
  i, j = args
  return i, j, lowest_levels(J=ratios[i] * U, theta=thetas[j])


if __name__ == '__main__':
  grid = [(i, j) for i in range(len(ratios)) for j in range(len(thetas))]
  energies = np.zeros((len(ratios), len(thetas), numval))

  start_time = timer()
  done = 0
  with mp.Pool(processes=min(8, os.cpu_count() or 1)) as pool:
    for i, j, ev in pool.imap_unordered(_work, grid, chunksize=1):
      energies[i, j, :] = ev
      done += 1
      if done % 20 == 0 or done == len(grid):
        el = timer() - start_time
        print('%3d/%3d points  elapsed: %5.1f s  (eta %5.1f s)'
              % (done, len(grid), el, el * (len(grid) - done) / done),
              flush=True)

  outfile = 'anyonHubbardGapsMulti_N%d.npz' % Nsites
  np.savez(outfile, thetas=thetas, ratios=ratios,
           energies=energies, Nsites=Nsites, U=U, usePBC=usePBC)
  print('Saved: %s' % outfile)

  plot_all(outfile)
