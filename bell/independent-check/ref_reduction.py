"""
ref_reduction.py -- independent check; shares no code with ../scripts/.

Checks, using spectral projectors only (no Fourier shortcut for the Bell value):
  (a) the BKP expression of quant-ph/0605182, eqs. (3)-(6), evaluated on Phi_d with BKP's bases, against
      the closed form  I_N^BKP(d) = N(d-1) - (2N/d) sum_m (d-m) sin(pi(N-1)m/(Nd))/sin(pi m/d)
      and against Tavakoli-Zohren-Pawlowski (arXiv:1506.04331) eq. (8): d*(1/d)sum_kl M_kl - 1;
  (b) for random projective strategies on Phi_D (arbitrary ranks, D not nec. multiple of d):
      I(R) = I(twirl) = I(R^V) = N(d-1) - (2/d) F_N(V), with F_N from the h_m formula and from the
      projection (K_N) form, and V_k^{2N} = 1.
The twirl is built with an explicit basis of ker(W^{2N} - 1) (no eigensolver).
Usage: python ref_reduction.py
"""
import itertools
import sys

import numpy as np

rng = np.random.default_rng(20261008)


def proj_of_unitary(R, d):
    """spectral projectors P_a of a unitary with R^d = 1, eigenvalue w^a."""
    D = R.shape[0]
    w = np.exp(2j * np.pi / d)
    Ps = []
    Rn = [np.eye(D, dtype=complex)]
    for n in range(1, d):
        Rn.append(Rn[-1] @ R)
    for a in range(d):
        P = sum(w ** (-a * n) * Rn[n] for n in range(d)) / d
        Ps.append((P + P.conj().T) / 2)
    return Ps


def chain_value(R, d, Ps=None):
    """I = sum_i sum_{x,y} tau(P_x(R_i) P_y(R_{i+1})) [(x - y - delta_{i,2N}) mod d]."""
    n = len(R)
    D = R[0].shape[0]
    if Ps is None:
        Ps = [proj_of_unitary(Ri, d) for Ri in R]
    tot = 0.0
    for i in range(n):
        A, B = Ps[i], Ps[(i + 1) % n]
        sh = 1 if i == n - 1 else 0
        for x in range(d):
            for y in range(d):
                p = np.real(np.trace(A[x] @ B[y])) / D
                tot += p * ((x - y - sh) % d)
    return tot


def strategy_value(PA, QB, d):
    """direct BKP value from PVMs on Phi_D: P(a,b|k,l) = Tr(P^k_a (Q^l_b)^T)/D."""
    N = len(PA)
    D = PA[0][0].shape[0]

    def E(Pk, Ql, sign, shift):
        s = 0.0
        for a in range(d):
            for b in range(d):
                p = np.real(np.trace(Pk[a] @ Ql[b].T)) / D
                s += p * ((sign * (a - b) - shift) % d)
        return s
    tot = 0.0
    for k in range(N):
        tot += E(PA[k], QB[k], +1, 0)                 # [A_k - B_k]
        if k < N - 1:
            tot += E(PA[k + 1], QB[k], -1, 0)         # [B_k - A_{k+1}]
    tot += E(PA[0], QB[N - 1], -1, 1)                 # [B_N - A_1 - 1]
    return tot


def bkp_bases(N, d):
    q = np.arange(d)
    PA, QB = [], []
    for k in range(1, N + 1):
        al = (k - 0.5) / N
        vecs = [np.exp(2j * np.pi * q * (r - al) / d) / np.sqrt(d) for r in range(d)]
        PA.append([np.outer(v, v.conj()) for v in vecs])
    for l in range(1, N + 1):
        be = l / N
        vecs = [np.exp(-2j * np.pi * q * (r - be) / d) / np.sqrt(d) for r in range(d)]
        QB.append([np.outer(v, v.conj()) for v in vecs])
    return PA, QB


def closed_form(N, d):
    return N * (d - 1) - (2 * N / d) * sum((d - m) * np.sin(np.pi * (N - 1) * m / (N * d)) / np.sin(np.pi * m / d)
                                           for m in range(1, d))


def tzp_value(N, d):
    M = np.zeros((d, d))
    for k in range(d):
        for l in range(d):
            if k == l:
                M[k, l] = N - (N - 1) / d
            else:
                t = k - l
                M[k, l] = -(N / d) * np.sin((N - 1) * np.pi * t / (d * N)) / np.sin(np.pi * t / d)
    B = M.sum() / d
    return d * B - 1


def haar(n):
    Z = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    Qm, Rm = np.linalg.qr(Z)
    return Qm * (np.diag(Rm) / np.abs(np.diag(Rm)))


def random_pvm(D, d, labels=None):
    U = haar(D)
    if labels is None:
        labels = rng.integers(0, d, size=D)
    return [U @ np.diag((labels == a).astype(float)) @ U.conj().T for a in range(d)]


def unitary_from_pvm(P, d):
    w = np.exp(2j * np.pi / d)
    return sum(w ** a * P[a] for a in range(d))


def twirl_and_reduce(R, N, d):
    """returns (Rt list on C^{L D}, V list on C^{2ND})."""
    D = R[0].shape[0]
    L = 2 * N * d
    w = np.exp(2j * np.pi / d)

    def Rgen(j):   # "R_j" for any integer j >= 1 with R_{j+2N} = w R_j
        q, r = divmod(j - 1, 2 * N)
        return (w ** q) * R[r]
    n = L * D
    Rt = []
    for i in range(1, 2 * N + 1):
        M_ = np.zeros((n, n), complex)
        for r in range(L):
            M_[r * D:(r + 1) * D, r * D:(r + 1) * D] = Rgen(i + r)
        Rt.append(M_)
    Wsh = np.zeros((n, n))
    for r in range(L):
        Wsh[r * D:(r + 1) * D, ((r + 1) % L) * D:((r + 1) % L + 1) * D] = np.eye(D)
    # sanity: covariance
    cov = max(np.abs(Wsh @ Rt[i] @ Wsh.T - Rt[i + 1]).max() for i in range(2 * N - 1))
    cov = max(cov, np.abs(Wsh @ Rt[2 * N - 1] @ Wsh.T - w * Rt[0]).max())
    # explicit orthonormal basis of H_0 = ker(W^{2N} - 1)
    M = 2 * N * D
    E0 = np.zeros((n, M), complex)
    for s in range(2 * N):
        for j in range(D):
            col = s * D + j
            for t in range(d):
                E0[(s + 2 * N * t) * D + j, col] = 1 / np.sqrt(d)
    Bm = np.hstack([np.linalg.matrix_power(Rt[0], k) @ E0 for k in range(d)])
    unit_err = np.abs(Bm.conj().T @ Bm - np.eye(n)).max()
    Wb = Bm.conj().T @ Wsh @ Bm
    offblock = 0.0
    for k in range(d):
        for kk in range(d):
            if k != kk:
                offblock = max(offblock, np.abs(Wb[k * M:(k + 1) * M, kk * M:(kk + 1) * M]).max())
    z = np.exp(2j * np.pi / L)
    V = [z ** (-k) * Wb[k * M:(k + 1) * M, k * M:(k + 1) * M] for k in range(d)]
    X1 = Bm.conj().T @ Rt[0] @ Bm
    Xs = np.kron(np.roll(np.eye(d), 1, axis=0), np.eye(M))
    return Rt, V, dict(cov=cov, unit=unit_err, offblock=offblock, X=np.abs(X1 - Xs).max())


def covariant_from_V(V, N, d):
    M = V[0].shape[0]
    z = np.exp(2j * np.pi / (2 * N * d))
    W = np.zeros((d * M, d * M), complex)
    for k in range(d):
        W[k * M:(k + 1) * M, k * M:(k + 1) * M] = z ** k * V[k]
    R = [np.kron(np.roll(np.eye(d), 1, axis=0), np.eye(M))]
    for i in range(2 * N - 1):
        R.append(W @ R[-1] @ W.conj().T)
    return R


def F_h(V, N, d):
    M = V[0].shape[0]
    F = 0.0
    for j in range(d):
        for k in range(j + 1, d):
            m = k - j
            h = -1j * N * np.exp(1j * np.pi * m * (N - 1) / (N * d)) / np.sin(np.pi * m / d)
            F += np.real(h * np.trace(V[j] @ V[k].conj().T) / M)
    return F


def F_proj(V, N, d):
    """Q_{k - d a} = spectral projection of V_k for eigenvalue zeta^{-a}, zeta = e^{i pi/N}."""
    L = 2 * N * d
    M = V[0].shape[0]
    Q = [None] * L
    for k in range(d):
        # V_k^{2N} = 1: projector for eigenvalue mu = zeta^{-a}: (1/2N) sum_n mu^{-n} V^n
        Vn = [np.eye(M, dtype=complex)]
        for n in range(1, 2 * N):
            Vn.append(Vn[-1] @ V[k])
        for a in range(2 * N):
            mu = np.exp(-1j * np.pi * a / N)
            P = sum(mu ** (-n) * Vn[n] for n in range(2 * N)) / (2 * N)
            Q[(k - d * a) % L] = (P + P.conj().T) / 2
    F = 0.0
    for x in range(L):
        for y in range(L):
            if (x - y) % d:
                u = x - y
                K = N * np.sin(np.pi * (N - 1) * u / (N * d)) / np.sin(np.pi * u / d)
                F += K * np.real(np.trace(Q[x] @ Q[y])) / M
    return F / 2


def main():
    print("== (a) BKP measurements on Phi_d: direct vs closed form vs TZP eq.(8)")
    worst = 0.0
    for N in range(2, 9):
        for d in range(2, 9):
            PA, QB = bkp_bases(N, d)
            v = strategy_value(PA, QB, d)
            cf = closed_form(N, d)
            tz = tzp_value(N, d)
            worst = max(worst, abs(v - cf), abs(tz - cf))
            if d in (3, 5) and N in (3, 4, 8):
                print(f"   N={N} d={d}: direct={v:.13f} closed={cf:.13f} TZP={tz:.13f} local bound={d - 1}")
    print(f"   max |direct-closed|, |TZP-closed| over N=2..8, d=2..8: {worst:.2e}")
    print("== (b) random projective strategies: I(R) vs twirl vs R^V vs N(d-1)-(2/d)F_N(V)")
    worst = dict(tw=0.0, rv=0.0, fh=0.0, fp=0.0, v2n=0.0, cov=0.0, unit=0.0, off=0.0, X=0.0)
    cases = [(3, 3, 3), (3, 3, 4), (3, 3, 5), (4, 3, 4), (3, 4, 4), (3, 4, 6), (5, 3, 3), (4, 4, 5), (3, 5, 5), (3, 5, 7)]
    for (N, d, D) in cases:
        for trial in range(2):
            PA = [random_pvm(D, d) for _ in range(N)]
            QB = [random_pvm(D, d) for _ in range(N)]
            I0 = strategy_value(PA, QB, d)
            R = []
            for k in range(N):
                R.append(unitary_from_pvm(PA[k], d))
                R.append(unitary_from_pvm(QB[k], d).T)
            I1 = chain_value(R, d)
            Rt, V, chk = twirl_and_reduce(R, N, d)
            I2 = chain_value(Rt, d)
            RV = covariant_from_V(V, N, d)
            I3 = chain_value(RV, d)
            Fh = F_h(V, N, d)
            Fp = F_proj(V, N, d)
            I4 = N * (d - 1) - (2 / d) * Fh
            v2n = max(np.abs(np.linalg.matrix_power(Vk, 2 * N) - np.eye(Vk.shape[0])).max() for Vk in V)
            worst['tw'] = max(worst['tw'], abs(I0 - I1), abs(I1 - I2))
            worst['rv'] = max(worst['rv'], abs(I2 - I3))
            worst['fh'] = max(worst['fh'], abs(I3 - I4))
            worst['fp'] = max(worst['fp'], abs(Fh - Fp))
            worst['v2n'] = max(worst['v2n'], v2n)
            worst['cov'] = max(worst['cov'], chk['cov'])
            worst['unit'] = max(worst['unit'], chk['unit'])
            worst['off'] = max(worst['off'], chk['offblock'])
            worst['X'] = max(worst['X'], chk['X'])
            print(f"   N={N} d={d} D={D} t={trial}: I(R)={I0:.12f} I(twirl)={I2:.12f} I(R^V)={I3:.12f} "
                  f"N(d-1)-(2/d)F={I4:.12f}  |F_h-F_proj|={abs(Fh - Fp):.1e}", flush=True)
    print("   worst:", {k: f"{v:.1e}" for k, v in worst.items()})


if __name__ == "__main__":
    main()
