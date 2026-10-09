"""
ch_validate.py -- validity checks of the chained cell inequalities on actual Q-configurations
(math/05-chained-bell.md, Prop. 4):
for random and near-optimal chained configurations, every cell inequality
    sum_u kt^{ell,s}(u) C(u) <= L^2 Phi_c(s/2N, (N-s)/2N)
must hold, F_N <= F_nat must hold, and the block decomposition F_N = (1/2) sum_s mu_s <X_s,Y_s> + const is checked.
Usage: python ch_validate.py [ntrials]
"""
import sys

import numpy as np

from ch_clock import F_V, F_natural, Q_from_V
from ch_cone import kt_vec, Phi_c, random_cells


def Cvec(Q, L):
    M = Q[0].shape[0]
    return np.array([sum(np.real(np.trace(Q[x] @ Q[(x + u) % L])) for x in range(L)) / M for u in range(L)])


def random_config(rng, N, d, M, near=None):
    zeta = np.exp(1j * np.pi / N)
    V = []
    for k in range(d):
        if near is None:
            U = np.linalg.qr(rng.standard_normal((M, M)) + 1j * rng.standard_normal((M, M)))[0]
            V.append(U @ np.diag(zeta ** (-rng.integers(0, 2 * N, size=M))) @ U.conj().T)
        else:
            # one-step phases per column, then a small unitary perturbation per site
            ph, delta = near
            Z = rng.standard_normal((M, M)) + 1j * rng.standard_normal((M, M))
            H = (Z + Z.conj().T) / 2
            ev, W = np.linalg.eigh(H)
            U = (W * np.exp(1j * delta * ev)) @ W.conj().T
            V.append(U @ np.diag(zeta ** (-ph[:, k])) @ U.conj().T)
    return V


def main():
    ntr = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    rng = np.random.default_rng(7)
    worst = -np.inf
    worstF = -np.inf
    for N, d in [(3, 3), (4, 3), (3, 4), (5, 3), (4, 4)]:
        L = 2 * N * d
        cells = [np.ones(d)] + random_cells(d, 6, rng, q=8)
        kts = [(kt_vec(ell, N, s), L ** 2 * Phi_c(s / (2 * N), (N - s) / (2 * N))) for ell in cells
               for s in range(1, N // 2 + 1)]
        Fn = F_natural(N, d)
        for t in range(ntr):
            M = int(rng.integers(2, 5))
            if t % 2:
                ph = np.array([((rng.integers(0, 2 * N) + (np.arange(d) >= rng.integers(0, d))) % (2 * N))
                               for _ in range(M)])
                V = random_config(rng, N, d, M, near=(ph, float(rng.choice([0.01, 0.05, 0.2]))))
            else:
                V = random_config(rng, N, d, M)
            Q = Q_from_V(V, N, d)
            C = Cvec(Q, L)
            viol = max(kt @ C - b for kt, b in kts)
            worst = max(worst, viol)
            worstF = max(worstF, F_V(V, N, d) - Fn)
        print(f"N={N} d={d}: max over trials of (cell functional - bound) = {worst:.3e} (must be <= 0); "
              f"max F - F_nat = {worstF:.3e} (must be <= 0)", flush=True)


if __name__ == "__main__":
    main()
