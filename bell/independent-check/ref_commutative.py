"""
ref_commutative.py -- independent check of Theorem 3 of math/05-chained-bell.md (commutative sector of the chained
clock model):
over all a in Z_{2N}^d (scalar V_k = zeta^{-a_k}), F_N(a) <= F_nat, with equality exactly for the one-step
configurations a_k = c + [k >= r]. Uses the K_N projection form on Z_L (independent of the h_m form).
Usage: python ref_commutative.py
"""
import itertools

import numpy as np


def KN(N, d, u):
    return N * np.sin(np.pi * (N - 1) * u / (N * d)) / np.sin(np.pi * u / d)


def F_of(a, N, d):
    L = 2 * N * d
    S = [(k - d * a[k]) % L for k in range(d)]
    return 0.5 * sum(KN(N, d, x - y) for x in S for y in S if (x - y) % d)


def main():
    for (N, d) in [(3, 3), (3, 4), (4, 3), (4, 4), (5, 3), (3, 5), (6, 3)]:
        Fnat = sum((d - m) * KN(N, d, m) for m in range(1, d))
        one_step = {tuple(((c + (np.arange(d) >= r)) % (2 * N)).tolist()) for c in range(2 * N) for r in range(d)}
        best_non, worst_one = -np.inf, np.inf
        for a in itertools.product(range(2 * N), repeat=d):
            F = F_of(a, N, d)
            if a in one_step:
                worst_one = min(worst_one, F)
            else:
                best_non = max(best_non, F)
        print(f"N={N} d={d}: F_nat={Fnat:.10f}  min over one-step={worst_one:.10f}  "
              f"max over non-one-step={best_non:.10f}  gap={Fnat - best_non:.4f}  (#one-step={len(one_step)})")


if __name__ == "__main__":
    main()
