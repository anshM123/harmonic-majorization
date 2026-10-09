"""Check S8: light independent probe of Conjecture D (03-structure.md, Section 2; open) with non-separable
superharmonic y-convex test functions of linear growth on the strip [b_1, b_m] x R, m = 3..6 levels:
   u = psi(x) sqrt(s0^2 + (y-y0)^2),  psi = cos(k(x-c)) or min of two such cosines,  k <= pi/w (psi >= 0 on the strip),
   k >= 1/s0 (=> Lap u <= 0), plus -phi(x) and positive combinations.
Random pairs and Nelder-Mead adversary over g.  A violation would refute Conjecture D (none expected).
"""
import time
import numpy as np
from scipy.optimize import minimize
from hmaudit import rand_herm, setup

rng = np.random.default_rng(61)
t0 = time.time()


def tests(b1, bm, rng):
    w = bm - b1
    out = []
    for _ in range(3):
        c = b1 + w * rng.uniform(0.3, 0.7)
        k = rng.uniform(0.2, 1.0) * np.pi / (2 * max(c - b1, bm - c))   # |x - c| <= pi/(2k) on the strip
        s0 = rng.uniform(1.0, 3.0) / k
        y0 = rng.normal()
        out.append(lambda z, c=c, k=k, s0=s0, y0=y0: np.cos(k * (z.real - c)) * np.sqrt(s0 ** 2 + (z.imag - y0) ** 2))
    c1, c2 = b1 + 0.35 * w, b1 + 0.65 * w
    k = 0.45 * np.pi / w
    s0 = 1.5 / k
    out.append(lambda z: np.minimum(np.cos(k * (z.real - c1)), np.cos(k * (z.real - c2))) * np.sqrt(s0 ** 2 + z.imag ** 2))
    return out


worst = (np.inf, None)
N = 20000
for trial in range(N):
    m = int(rng.integers(3, 7))
    lv = np.sort(rng.uniform(0, 1, m))
    mult = rng.integers(1, 3, size=m)
    levels = np.repeat(lv, mult)
    M = len(levels)
    mode = rng.integers(3)
    if mode == 0:
        g = rand_herm(M, rng, rng.uniform(0.05, 5))
    elif mode == 1:
        g = np.diag(rng.normal(size=M)).astype(complex) + rand_herm(M, rng, rng.uniform(1e-3, 0.3))
    else:  # chain coupling only between consecutive levels, ~ half the gaps (exceptional-point-like)
        g = np.diag(np.full(M, rng.normal())).astype(complex)
        for i in range(M - 1):
            if levels[i] != levels[i + 1]:
                v = (levels[i + 1] - levels[i]) / 2 * rng.uniform(0.8, 1.2)
                g[i, i + 1] = g[i + 1, i] = v
    c = setup(list(levels), g)
    if c["K2"] < 1e-12:
        continue
    for u in tests(levels.min(), levels.max(), rng):
        D = np.sum(u(c["nu"])) - np.sum(u(c["nuE"]))
        if D / c["K2"] < worst[0]:
            worst = (D / c["K2"], (trial, m, mode))
print(f"random: {N} pairs (3-6 levels), 4 test functions each: min D/||K||^2 = {worst[0]:.3e} at {worst[1]}")
print(f"[elapsed {time.time()-t0:.0f}s]")


def gfrom(p, M):
    A = p[:M * M].reshape(M, M) + 1j * p[M * M:].reshape(M, M)
    return (A + A.conj().T) / 2


advmin = np.inf
for rnd in range(30):
    m = 3 + rnd % 3
    levels = np.sort(np.concatenate([[0.0, 1.0], rng.uniform(0, 1, m - 2)]))
    M = m
    us = tests(0.0, 1.0, rng)
    u = us[rnd % len(us)]

    def obj(p):
        g = gfrom(p, M)
        c = setup(list(levels), g)
        if c["K2"] < 1e-8:
            return 1.0
        return (np.sum(u(c["nu"])) - np.sum(u(c["nuE"]))) / c["K2"]
    best = np.inf
    for start in range(2):
        res = minimize(obj, rng.normal(size=2 * M * M) * rng.uniform(0.2, 4), method="Nelder-Mead",
                       options=dict(maxiter=5000, xatol=1e-10, fatol=1e-14))
        best = min(best, res.fun)
    advmin = min(advmin, best)
    print(f"adversarial round {rnd:2d} (m={m}): min D/||K||^2 = {best:.4e}", flush=True)
print(f"adversarial minimum: {advmin:.4e}")
print(f"time {time.time()-t0:.0f}s")
