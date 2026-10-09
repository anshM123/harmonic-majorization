"""
ch_twirl_check.py -- Theorem 1 of math/05-chained-bell.md:
for an arbitrary projective strategy R = (R_1..R_2N) on Phi_D
(R_i^d = 1), the twirl over rho(R) = (R_2, ..., R_2N, w R_1) (order 2Nd) and Stone-von Neumann for Z_d give a reduced
family V_0..V_{d-1} (V_k^{2N} = 1) with  S_N(R) = N(d-1) - (2/d) F_N(V).
Usage: python ch_twirl_check.py N d D [ntrials]
"""
import sys

import numpy as np

from ch_clock import F_V, S_direct


def haar(rng, n):
    Z = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    q, r = np.linalg.qr(Z)
    return q * (np.diag(r) / np.abs(np.diag(r)))


def reduced_family(R, N, d):
    D = R[0].shape[0]
    L = 2 * N * d
    w = np.exp(2j * np.pi / d)
    blocks = [list(R)]
    for r in range(1, L):
        b = blocks[-1]
        blocks.append(b[1:] + [w * b[0]])
    n = L * D
    Rt1 = np.zeros((n, n), complex)
    for r in range(L):
        Rt1[r * D:(r + 1) * D, r * D:(r + 1) * D] = blocks[r][0]
    Wsh = np.zeros((n, n))
    for r in range(L):
        Wsh[r * D:(r + 1) * D, ((r + 1) % L) * D:((r + 1) % L + 1) * D] = np.eye(D)
    Z = np.linalg.matrix_power(Wsh, 2 * N)
    ev, U = np.linalg.eig(Z)
    H0, _ = np.linalg.qr(U[:, np.abs(ev - 1) < 1e-8])
    M = H0.shape[1]
    Bm = np.hstack([np.linalg.matrix_power(Rt1, k) @ H0 for k in range(d)])
    Wb = Bm.conj().T @ Wsh @ Bm
    z = np.exp(2j * np.pi / L)
    return [z ** (-k) * Wb[k * M:(k + 1) * M, k * M:(k + 1) * M] for k in range(d)]


def main():
    N, d, D = (int(a) for a in sys.argv[1:4])
    ntr = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    rng = np.random.default_rng(0)
    for t in range(ntr):
        R = []
        for i in range(2 * N):
            U = haar(rng, D)
            R.append(U @ np.diag(np.exp(2j * np.pi * rng.integers(0, d, size=D) / d)) @ U.conj().T)
        S = S_direct(R, N, d)
        V = reduced_family(R, N, d)
        err = max(np.abs(np.linalg.matrix_power(Vk, 2 * N) - np.eye(Vk.shape[0])).max() for Vk in V)
        F = F_V(V, N, d)
        print(f"N={N} d={d} D={D}: S(R)={S:.12f}   N(d-1)-(2/d)F_N(V)={N * (d - 1) - 2 / d * F:.12f}   "
              f"|V^2N - 1|={err:.1e}")


if __name__ == "__main__":
    main()
