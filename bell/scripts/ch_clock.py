"""
ch_clock.py -- clock model of the chained (Barrett-Kent-Pironio) Bell expression with N settings and d outcomes
(math/05-chained-bell.md).

Chain form (N = 2 is CGLMP in the chain form of [CL, Lemma 2.1]; [CL] = github.com/anshM123/CGLMP,
paper-classical-all-d/main.tex):
    S_N = sum_{i=1}^{2N} E m(link_i),  links X_1 - Y_1, Y_1 - X_2, ..., X_N - Y_N, Y_N - X_1 - 1,  m(t) = t mod d,
local bound S_N >= d - 1.  On Phi_D with unitaries R_1..R_{2N} (R_i^d = 1; Alice's and transposed Bob's alternately):
    S_N = N(d-1) + sum_n c_n sum_i tau(L^(i)_n),  c_n = -1/(1 - w^{-n}),  L^(i)_n = R_i^n R_{i+1}^{-n},  L^(2N)_n = w^{-n} R_{2N}^n R_1^{-n}.
Covariant strategies: W = sum_k |k><k| (x) z^k V_k, z = e^{2 pi i/(2Nd)}, V_k^{2N} = 1, R_1 = X (x) 1, R_{i+1} = W R_i W^*.
Claimed reduction (Theorem 1 of math/05-chained-bell.md):  S_N(R^V) = N(d-1) - (2/d) F_N(V),
    F_N(V) = sum_{j<k} Re[h_{k-j} tau(V_j V_k^*)],  h_m = -i N e^{i theta_m (1 - 1/N)} / sin(theta_m),  theta_m = pi m/d,
and the projection form F_N = (1/2) sum_{x != y mod d} tau(Q_x Q_y) K_N(x - y) on Z_{2Nd}, K_N(u) = N sin(pi(N-1)u/(Nd))/sin(pi u/d).
This script checks the reduction against a direct evaluation of S_N, enumerates the commutative sector, and runs a
numerical search over non-commuting families.
Usage: python ch_clock.py N d [M] [restarts]
"""
import itertools
import sys

import numpy as np
from scipy.optimize import minimize


def h_coef(N, d, m):
    th = np.pi * m / d
    return -1j * N * np.exp(1j * th * (1 - 1 / N)) / np.sin(th)


def K_N(N, d, u):
    u = np.asarray(u, float)
    return N * np.sin(np.pi * (N - 1) * u / (N * d)) / np.sin(np.pi * u / d)


def F_V(V, N, d):
    M = V[0].shape[0]
    F = 0.0
    for j in range(d):
        for k in range(j + 1, d):
            F += np.real(h_coef(N, d, k - j) * np.trace(V[j] @ V[k].conj().T) / M)
    return F


def F_natural(N, d):
    return sum((d - m) * np.real(h_coef(N, d, m)) for m in range(1, d))


def covariant_R(V, N, d):
    M = V[0].shape[0]
    z = np.exp(2j * np.pi / (2 * N * d))
    X = np.roll(np.eye(d), 1, axis=0)          # X|k> = |k+1>
    W = np.zeros((d * M, d * M), complex)
    for k in range(d):
        W[k * M:(k + 1) * M, k * M:(k + 1) * M] = z ** k * V[k]
    R = [np.kron(X, np.eye(M))]
    for i in range(2 * N - 1):
        R.append(W @ R[-1] @ W.conj().T)
    return R


def S_direct(R, N, d):
    """S_N = sum over links of E m(link) with E w^{n t} = tau(L_n) (unitaries R_i^d = 1, normalised trace)."""
    D = R[0].shape[0]
    w = np.exp(2j * np.pi / d)
    S = 0.0
    for i in range(2 * N):
        A = R[i]
        B = R[(i + 1) % (2 * N)]
        # distribution of t = link value: P(t) = (1/d) sum_n w^{-nt} E w^{n t}
        Ew = [np.trace(np.linalg.matrix_power(A, n) @ np.linalg.matrix_power(B.conj().T, n)) / D for n in range(d)]
        if i == 2 * N - 1:
            Ew = [w ** (-n) * Ew[n] for n in range(d)]
        P = np.array([np.real(sum(w ** (-n * t) * Ew[n] for n in range(d))) / d for t in range(d)])
        S += float(np.dot(np.arange(d), P))
    return S


def Q_from_V(V, N, d):
    """V_k = sum_a zeta^{-a} P_{k,a}, zeta = e^{i pi/N}; Q_{k - d a} = P_{k,a} (x in Z_{2Nd})."""
    L = 2 * N * d
    zeta = np.exp(1j * np.pi / N)
    Q = [None] * L
    for k in range(d):
        E = V[k]
        M = E.shape[0]
        for a in range(2 * N):
            P = np.zeros((M, M), complex)
            X = np.eye(M, dtype=complex)
            for n in range(2 * N):
                P += X
                X = X @ (zeta ** a * E)
            P /= 2 * N
            Q[(k - d * a) % L] = (P + P.conj().T) / 2
    return Q


def F_proj(Q, N, d):
    L = 2 * N * d
    M = Q[0].shape[0]
    F = 0.0
    for x in range(L):
        for y in range(L):
            if (x - y) % d:
                F += K_N(N, d, x - y) * np.real(np.trace(Q[x] @ Q[y])) / M
    return F / 2


def commutative_max(N, d):
    best = -np.inf
    arg = None
    vals = {}
    zeta = np.exp(1j * np.pi / N)
    for a in itertools.product(range(2 * N), repeat=d):
        F = 0.0
        for j in range(d):
            for k in range(j + 1, d):
                F += np.real(h_coef(N, d, k - j) * zeta ** (-(a[j] - a[k])))
        vals[a] = F
        if F > best + 1e-12:
            best, arg = F, a
    return best, arg, vals


def one_step(N, d):
    out = set()
    for c in range(2 * N):
        for r in range(d):
            out.add(tuple(((c + (np.arange(d) >= r)) % (2 * N)).tolist()))
    return out


def params_to_V(p, N, d, M, spectra):
    V = []
    off = 0
    for k in range(d):
        A = p[off:off + M * M].reshape(M, M)
        off += M * M
        H = (A + A.T) / 2 + 1j * (A - A.T) / 2
        ev, U = np.linalg.eigh(H)
        Uk = (U * np.exp(1j * ev)) @ U.conj().T
        zeta = np.exp(1j * np.pi / N)
        V.append(Uk @ np.diag(zeta ** (-spectra[k])) @ Uk.conj().T)
    return V


def search(N, d, M, restarts, seed=0):
    rng = np.random.default_rng(seed)
    best = -np.inf
    for r in range(restarts):
        spectra = [rng.integers(0, 2 * N, size=M) for _ in range(d)]
        p0 = rng.standard_normal(d * M * M) * 1.5
        res = minimize(lambda p: -F_V(params_to_V(p, N, d, M, spectra), N, d), p0, method="L-BFGS-B",
                       options=dict(maxiter=3000))
        best = max(best, -res.fun)
    return best


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    d = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    M = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    restarts = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    rng = np.random.default_rng(1)
    # 1. reduction check on random families
    for t in range(3):
        V = []
        for k in range(d):
            U = np.linalg.qr(rng.standard_normal((M, M)) + 1j * rng.standard_normal((M, M)))[0]
            V.append(U @ np.diag(np.exp(1j * np.pi / N) ** (-rng.integers(0, 2 * N, size=M))) @ U.conj().T)
        R = covariant_R(V, N, d)
        S = S_direct(R, N, d)
        F = F_V(V, N, d)
        Fp = F_proj(Q_from_V(V, N, d), N, d)
        print(f"check {t}: S_direct={S:.12f}  N(d-1)-(2/d)F={N * (d - 1) - 2 / d * F:.12f}  F={F:.12f}  F_proj={Fp:.12f}")
    Fn = F_natural(N, d)
    print(f"N={N} d={d}: natural strategy F={Fn:.12f}, S={N * (d - 1) - 2 / d * Fn:.12f} (local bound d-1={d - 1})")
    if (2 * N) ** d <= 3e6:
        best, arg, vals = commutative_max(N, d)
        os_ = one_step(N, d)
        best_non = max(v for a, v in vals.items() if a not in os_)
        print(f"commutative max F={best:.12f} at a={arg}; natural={Fn:.12f}; one-step attain? "
              f"{all(abs(vals[a] - Fn) < 1e-9 for a in os_)}; best non-one-step={best_non:.12f}")
    if restarts > 0:
        b = search(N, d, M, restarts)
        print(f"non-commuting search (M={M}, {restarts} restarts): best F={b:.12f}  vs natural {Fn:.12f}  "
              f"(diff {b - Fn:+.3e})")


if __name__ == "__main__":
    main()
