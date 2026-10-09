"""
ch_local.py -- local (second-order) optimality of the natural strategy for the chained clock model with N settings:
at a commuting one-step configuration (windows c, c' in Z_{2Nd}) perturb by site-wise unitaries mixing the two joint
eigenvectors; report the smallest eigenvalue of the Hessian of the deficit F_nat - F_N (modulo the gauge direction).
A negative value would refute optimality of the natural strategy on maximally entangled states.
Also: a seeded multi-start search with spectra compatible with direct sums of one-step configurations.
Usage: python ch_local.py N d [restarts]
"""
import sys

import numpy as np
from scipy.optimize import minimize

from ch_clock import F_V, F_natural


def one_steps(N, d):
    out = []
    for c in range(2 * N):
        for r in range(d):
            out.append(((c + (np.arange(d) >= r)) % (2 * N)))
    return out


def hessian_pair(N, d, a1, a2, delta=1e-4):
    zeta = np.exp(1j * np.pi / N)
    Fn = F_natural(N, d)

    def F_of(xr):
        h = xr[:d] + 1j * xr[d:]
        V = []
        for k in range(d):
            Vk = np.diag([zeta ** (-a1[k]), zeta ** (-a2[k])])
            X = np.array([[0, h[k]], [np.conj(h[k]), 0]])
            ev, U = np.linalg.eigh(X)
            Uk = (U * np.exp(1j * delta * ev)) @ U.conj().T
            V.append(Uk @ Vk @ Uk.conj().T)
        return (Fn - F_V(V, N, d)) / delta ** 2
    n = 2 * d
    E = np.eye(n)
    dg = [F_of(E[i]) for i in range(n)]
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            H[i, j] = H[j, i] = dg[i] if i == j else (F_of(E[i] + E[j]) - dg[i] - dg[j]) / 2
    g1 = np.concatenate([np.ones(d), np.zeros(d)])
    g2 = np.concatenate([np.zeros(d), np.ones(d)])
    Gm = np.stack([g1, g2], 1)
    Pp = np.eye(n) - Gm @ np.linalg.pinv(Gm)
    ev, U = np.linalg.eigh(Pp)
    B = U[:, ev > 0.5]
    return np.linalg.eigvalsh(B.T @ H @ B)


def search_compatible(N, d, M, restarts, seed=0):
    rng = np.random.default_rng(seed)
    zeta = np.exp(1j * np.pi / N)
    os_ = one_steps(N, d)
    best = -np.inf
    for r in range(restarts):
        cfg = [os_[i] for i in rng.integers(0, len(os_), size=M)]
        spectra = [np.array([c[k] for c in cfg]) for k in range(d)]

        def V_of(p):
            V = []
            off = 0
            for k in range(d):
                A = p[off:off + M * M].reshape(M, M)
                off += M * M
                Hm = (A + A.T) / 2 + 1j * (A - A.T) / 2
                ev, U = np.linalg.eigh(Hm)
                Uk = (U * np.exp(1j * ev)) @ U.conj().T
                V.append(Uk @ np.diag(zeta ** (-spectra[k])) @ Uk.conj().T)
            return V
        p0 = rng.standard_normal(d * M * M) * rng.choice([0.3, 1.0, 2.0])
        res = minimize(lambda p: -F_V(V_of(p), N, d), p0, method="L-BFGS-B", options=dict(maxiter=4000))
        best = max(best, -res.fun)
    return best


def main():
    N = int(sys.argv[1])
    d = int(sys.argv[2])
    restarts = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    os_ = one_steps(N, d)
    a0 = os_[0]
    mins = []
    for a in os_[1:]:
        ev = hessian_pair(N, d, a0, a)
        mins.append(ev[0])
    mins = np.array(mins)
    print(f"N={N} d={d}: Hessian of deficit over {len(mins)} window pairs: min eigenvalue {mins.min():.6f} "
          f"(negatives: {(mins < -1e-6).sum()})", flush=True)
    if restarts:
        Fn = F_natural(N, d)
        for M in (2, 3, 4):
            b = search_compatible(N, d, M, restarts, seed=M)
            print(f"   search M={M}, {restarts} restarts: best F={b:.10f}  natural={Fn:.10f}  diff={b - Fn:+.2e}",
                  flush=True)


if __name__ == "__main__":
    main()
