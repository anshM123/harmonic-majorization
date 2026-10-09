"""
ref_tzp_gap.py -- the quantity settled by Theorem 5 of math/05-chained-bell.md (minimum on maximally entangled
states) is NOT the
all-state quantum minimum of I_N (the quantity listed as unknown by Kaniewski et al., Quantum 3, 198 (2019)).
Direct evaluation: BKP measurements (quant-ph/0605182 eqs. (5)-(6)) on |psi> = sum_k lam_k |kk>, with lam the
lowest eigenvector of the Tavakoli-Zohren-Pawlowski matrix M (arXiv:1506.04331 eq. (7)); I = d*B - 1.
Usage: python ref_tzp_gap.py
"""
import numpy as np


def bkp_bases(N, d):
    q = np.arange(d)
    A = [np.array([np.exp(2j * np.pi * q * (r - (k - 0.5) / N) / d) / np.sqrt(d) for r in range(d)]) for k in range(1, N + 1)]
    B = [np.array([np.exp(-2j * np.pi * q * (r - l / N) / d) / np.sqrt(d) for r in range(d)]) for l in range(1, N + 1)]
    return A, B          # A[k][r] = vector |r>_{A_k}


def value_on_state(psi, N, d):
    """psi: d x d coefficient matrix (|psi> = sum_ij psi_ij |i>|j>)."""
    A, B = bkp_bases(N, d)

    def joint(Ak, Bl):
        P = np.zeros((d, d))
        for a in range(d):
            for b in range(d):
                amp = Ak[a].conj() @ psi @ Bl[b].conj()
                P[a, b] = abs(amp) ** 2
        return P
    tot = 0.0
    for k in range(N):
        P = joint(A[k], B[k])
        tot += sum(P[a, b] * ((a - b) % d) for a in range(d) for b in range(d))
        if k < N - 1:
            P = joint(A[k + 1], B[k])
            tot += sum(P[a, b] * ((b - a) % d) for a in range(d) for b in range(d))
    P = joint(A[0], B[N - 1])
    tot += sum(P[a, b] * ((b - a - 1) % d) for a in range(d) for b in range(d))
    return tot


def tzp_M(N, d):
    M = np.zeros((d, d))
    for k in range(d):
        for l in range(d):
            t = k - l
            M[k, l] = N - (N - 1) / d if t == 0 else -(N / d) * np.sin((N - 1) * np.pi * t / (d * N)) / np.sin(np.pi * t / d)
    return M


def main():
    for (N, d) in [(3, 3), (3, 4), (4, 3), (4, 5), (5, 4), (8, 10), (3, 24)]:
        Imax = value_on_state(np.eye(d) / np.sqrt(d), N, d)
        M = tzp_M(N, d)
        ev, U = np.linalg.eigh(M)
        lam = U[:, 0] * np.sign(U[:, 0].sum())
        Iopt = value_on_state(np.diag(lam), N, d)
        print(f"N={N} d={d}: I on Phi_d (BKP) = {Imax:.10f};  BKP measurements on lowest eigenvector of M: "
              f"I = {Iopt:.10f} (d*lambda_min - 1 = {d * ev[0] - 1:.10f});  gap = {Imax - Iopt:.4e}")


if __name__ == "__main__":
    main()
