"""
ch_bkp.py -- the BKP measurements (Barrett-Kent-Pironio, PRL 97, 170409 (2006)) on Phi_d: direct evaluation of
I_N = <[A1-B1]> + <[B1-A2]> + ... + <[AN-BN]> + <[BN-A1-1]>  (<[X]> = E (X mod d))
with |r>_{A_k} = d^{-1/2} sum_q e^{2 pi i q (r - alpha_k)/d}|q>, alpha_k = (k - 1/2)/N,
     |r>_{B_l} = d^{-1/2} sum_q e^{-2 pi i q (r - beta_l)/d}|q>, beta_l = l/N,
compared with the natural clock-model value S_nat = N(d-1) - (2/d) F_nat(N, d) and with the large-N formula 2 gamma/N.
Usage: python ch_bkp.py [Nmin Nmax dmin dmax]
"""
import sys

import numpy as np

from ch_clock import F_natural


def bases(d, N):
    q = np.arange(d)
    A = []
    B = []
    for k in range(1, N + 1):
        al = (k - 0.5) / N
        A.append(np.array([np.exp(2j * np.pi * q * (r - al) / d) for r in range(d)]).T / np.sqrt(d))   # columns r
    for l in range(1, N + 1):
        be = l / N
        B.append(np.array([np.exp(-2j * np.pi * q * (r - be) / d) for r in range(d)]).T / np.sqrt(d))
    return A, B


def joint(Aa, Bb, d):
    """P(a, b) = |<a|_A <b|_B Phi_d>|^2 = |sum_q conj(A[q,a]) conj(B[q,b])|^2 / d."""
    amp = (Aa.conj().T @ Bb.conj()) / np.sqrt(d)
    return np.abs(amp) ** 2


def I_N(d, N):
    A, B = bases(d, N)
    tot = 0.0
    for k in range(N):
        P = joint(A[k], B[k], d)            # [A_k - B_k]
        tot += sum(P[a, b] * ((a - b) % d) for a in range(d) for b in range(d))
        if k < N - 1:
            P = joint(A[k + 1], B[k], d)    # [B_k - A_{k+1}]
            tot += sum(P[a, b] * ((b - a) % d) for a in range(d) for b in range(d))
    P = joint(A[0], B[N - 1], d)            # [B_N - A_1 - 1]
    tot += sum(P[a, b] * ((b - a - 1) % d) for a in range(d) for b in range(d))
    return tot


def main():
    # optional ranges: python ch_bkp.py Nmin Nmax dmin dmax  (default: the original grid)
    if len(sys.argv) > 4:
        Nmin, Nmax, dmin, dmax = map(int, sys.argv[1:5])
        Ns, ds = range(Nmin, Nmax + 1), range(dmin, dmax + 1)
    else:
        Ns, ds = (2, 3, 4, 5, 8), (2, 3, 4, 5, 7)
    worst = 0.0
    for N in Ns:
        for d in ds:
            S_nat = N * (d - 1) - 2 / d * F_natural(N, d)
            gam = np.pi ** 2 / (4 * d ** 2) * sum(j / np.sin(np.pi * j / d) ** 2 for j in range(1, d))
            val = I_N(d, N)
            worst = max(worst, abs(val - S_nat))
            print(f"N={N} d={d}: BKP I_N = {val:.12f}   natural S_N = {S_nat:.12f}   2 gamma/N = {2 * gam / N:.6f}")
    print(f"max |BKP I_N - natural S_N| = {worst:.1e}")


if __name__ == "__main__":
    main()
