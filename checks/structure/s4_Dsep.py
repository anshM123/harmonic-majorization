"""Check S4: Theorem D_sep of 03-structure.md (any number of levels, separable test functions) and the trace
identity Re Tr X^2 = Tr B^2 - ||g||_F^2, plus the oblique-ridge extension under the curvature budget.

Part 1: trace identity and its consequence sum(Im nu)^2 - ||Eg||^2 = ||K||^2 - L, L = Tr B^2 - sum (Re nu)^2 >= 0
        (50 digits).
Part 2: random pairs (2-7 levels, repeats, near-exceptional g), u = a(x) + c(y) with c convex and
        a'' <= -kappa, kappa = sup over the y-range of Rect of c'' (borderline choice a'' = -kappa + concave part).
Part 3: oblique ridges: u = h + atilde(x) - sum_i d_i(Re(e^{i th_i} z)), h = (kh/2)(y^2 - x^2), d_i'' <= k_i,
        sum k_i sin^2 th_i = kh (borderline budget).
Part 4: adversarial minimisation of D/||K||^2 for separable u (3-5 levels).
Part 5: the CE1 configuration with u = y^2 (outside the separable class: a = 0 violates a'' <= -2): negative defect.
Usage: python s4_Dsep.py [seed] [nrandom]
"""
import sys
import time
import numpy as np


def logcosh(t):
    t = np.abs(t)
    return t + np.log1p(np.exp(-2 * t)) - np.log(2.0)

import mpmath as mp
from scipy.optimize import minimize
from hmaudit import rand_herm, setup

seed = int(sys.argv[1]) if len(sys.argv) > 1 else 11
NR = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
rng = np.random.default_rng(seed)
t0 = time.time()

print("== Part 1: trace identity (50 digits)")
mp.mp.dps = 50
worst = mp.mpf(0)
for trial in range(30):
    M = int(rng.integers(2, 7))
    levels = list(np.round(rng.normal(size=M), 3))
    if M > 2:
        levels[1] = levels[0]
    g = rand_herm(M, rng, rng.uniform(0.2, 3))
    X = mp.matrix(M, M)
    gm = mp.matrix(M, M)
    for i in range(M):
        for j in range(M):
            gm[i, j] = mp.mpc(g[i, j].real, g[i, j].imag)
            X[i, j] = (mp.mpf(levels[i]) if i == j else 0) + 1j * gm[i, j]
    trX2 = sum((X * X)[i, i] for i in range(M))
    trB2 = sum(mp.mpf(l) ** 2 for l in levels)
    g2 = sum(abs(gm[i, j]) ** 2 for i in range(M) for j in range(M))
    worst = max(worst, abs(mp.re(trX2) - (trB2 - g2)))
    ev, _ = mp.eig(X)
    lab = [sorted(set(levels)).index(l) for l in levels]
    Eg2 = sum(abs(gm[i, j]) ** 2 for i in range(M) for j in range(M) if lab[i] == lab[j])
    K2 = g2 - Eg2
    L = trB2 - sum(mp.re(e) ** 2 for e in ev)
    lhs = sum(mp.im(e) ** 2 for e in ev) - Eg2
    worst = max(worst, abs(lhs - (K2 - L)))
    if L < -mp.mpf(10) ** (-40):
        print("   L < 0 !", L)
print(f"   30 random pairs: max |Re Tr X^2 - (Tr B^2 - ||g||^2)| and max |sum Im^2 - ||Eg||^2 - (||K||^2 - L)| "
      f"= {mp.nstr(worst, 3)}; L >= 0 in all cases")


def kappa_c_family(rng):
    """returns (name, c, c'') with c convex; c'' bounded."""
    j = rng.integers(4)
    kap = rng.uniform(0.2, 4.0)
    s = rng.uniform(0.05, 2.0)
    y0 = rng.normal()
    if j == 0:
        return "k s^2 logcosh", (lambda y: kap * s ** 2 * logcosh(((y - y0) / s))), \
            (lambda y: kap / np.cosh((y - y0) / s) ** 2)
    if j == 1:
        return "k (y-y0)^2/2", (lambda y: kap * (y - y0) ** 2 / 2), (lambda y: kap + 0 * y)
    if j == 2:
        return "k s sqrt(s^2+Y^2)", (lambda y: kap * s * np.sqrt(s ** 2 + (y - y0) ** 2)), \
            (lambda y: kap * s ** 3 / (s ** 2 + (y - y0) ** 2) ** 1.5)
    return "k s^2 softplus", (lambda y: kap * s ** 2 * np.logaddexp(0, (y - y0) / s) * 4), \
        (lambda y: kap * 4 * np.exp((y - y0) / s) / (1 + np.exp((y - y0) / s)) ** 2)


def gen_pair(rng):
    m = int(rng.integers(2, 8))
    lv = np.sort(rng.normal(size=m) * rng.uniform(0.2, 3))
    mult = rng.integers(1, 3, size=m)
    levels = np.repeat(lv, mult)
    M = len(levels)
    mode = rng.integers(3)
    if mode == 0:
        g = rand_herm(M, rng, rng.uniform(0.05, 5))
    elif mode == 1:
        g = np.diag(rng.normal(size=M)).astype(complex) + rand_herm(M, rng, rng.uniform(1e-3, 0.3))
    else:  # exceptional-point-like: couplings ~ level gaps
        g = rand_herm(M, rng, 1e-3)
        for i in range(M - 1):
            v = abs(levels[i + 1] - levels[i]) / 2 * rng.uniform(0.8, 1.2) + 1e-3
            g[i, i + 1] += v
            g[i + 1, i] += v
    return levels, g


print(f"\n== Part 2: D_sep on {NR} random pairs (borderline a'' = -kappa_Rect + concave part)")
minr = np.inf
info = None
neg = 0
for trial in range(NR):
    levels, g = gen_pair(rng)
    c = setup(levels, g)
    if c["K2"] < 1e-12:
        continue
    name, cf, c2 = kappa_c_family(rng)
    ys = np.linspace(c["lmin"], c["lmax"], 2001)
    kap = np.max(c2(ys))  # sup of c'' over the y-range of Rect (grid)
    kap *= 1 + 1e-9
    x0 = rng.uniform(levels.min(), levels.max())
    r = rng.uniform(0, 2)
    j = rng.integers(3)
    if j == 0:
        a = lambda x: -kap * x ** 2 / 2
    elif j == 1:
        a = lambda x: -kap * x ** 2 / 2 - r * np.abs(x - x0)
    else:
        a = lambda x: -kap * x ** 2 / 2 - r * np.exp(x - x0)
    u = lambda z: a(z.real) + cf(z.imag)
    D = np.sum(u(c["nu"])) - np.sum(u(c["nuE"]))
    scale = np.sum(np.abs(u(c["nu"]))) + 1
    if D < -1e-11 * scale:
        neg += 1
    if D / c["K2"] < minr:
        minr, info = D / c["K2"], (trial, name, j, len(levels), float(c["K2"]))
print(f"   significant negatives: {neg};  min D/||K||^2 = {minr:.3e} at {info}")
print(f"   [elapsed {time.time()-t0:.0f}s]")

print("\n== Part 3: oblique ridges with borderline curvature budget")
minr3 = np.inf
neg3 = 0
for trial in range(NR // 4):
    levels, g = gen_pair(rng)
    c = setup(levels, g)
    if c["K2"] < 1e-12:
        continue
    nth = int(rng.integers(1, 4))
    ths = rng.uniform(0, np.pi, nth)
    ks = rng.uniform(0.1, 2, nth)
    ss = rng.uniform(0.1, 2, nth)
    kh = np.sum(ks * np.sin(ths) ** 2)
    x1 = rng.normal()
    rr = rng.uniform(0, 1)

    def u(z):
        val = kh / 2 * (z.imag ** 2 - z.real ** 2) - rr * np.abs(z.real - x1)
        for th, kk, s in zip(ths, ks, ss):
            t = (np.exp(1j * th) * z).real
            val = val - kk * s ** 2 * logcosh((t / s))
        return val
    D = np.sum(u(c["nu"])) - np.sum(u(c["nuE"]))
    if D < -1e-11 * (1 + np.sum(np.abs(u(c["nu"])))):
        neg3 += 1
    minr3 = min(minr3, D / c["K2"])
print(f"   {NR//4} pairs: significant negatives {neg3};  min D/||K||^2 = {minr3:.3e}")

print("\n== Part 4: adversarial minimisation (separable u, 3-5 levels)")


def gfrom(p, M):
    A = p[:M * M].reshape(M, M) + 1j * p[M * M:].reshape(M, M)
    return (A + A.conj().T) / 2


advmin = np.inf
for rnd in range(18):
    m = 3 + rnd % 3
    levels = np.sort(rng.uniform(0, 1, m))
    M = m
    kap, s, y0 = rng.uniform(0.5, 3), rng.uniform(0.1, 1), rng.normal() * 0.3

    def obj(p):
        g = gfrom(p, M)
        c = setup(levels, g)
        if c["K2"] < 1e-8:
            return 1.0
        ys = np.linspace(c["lmin"], c["lmax"], 801)
        k = np.max(kap / np.cosh((ys - y0) / s) ** 2)
        u = lambda z: -k * z.real ** 2 / 2 + kap * s ** 2 * logcosh(((z.imag - y0) / s))
        return (np.sum(u(c["nu"])) - np.sum(u(c["nuE"]))) / c["K2"]
    best = np.inf
    for start in range(3):
        res = minimize(obj, rng.normal(size=2 * M * M) * rng.uniform(0.1, 2), method="Nelder-Mead",
                       options=dict(maxiter=5000, xatol=1e-10, fatol=1e-14))
        best = min(best, res.fun)
    advmin = min(advmin, best)
    print(f"   round {rnd:2d} (m={m}): min D/||K||^2 = {best:.4e}")
print(f"   adversarial minimum: {advmin:.4e}")

print("\n== Part 5: CE1 with u = y^2 (not in the separable class): defect = ||K||^2 - L")
levels = [1.0, 0.0, 0.0]
g = np.array([[0, .25, 0], [.25, 0, .01], [0, .01, 0]], complex)
c = setup(levels, g)
L = 1.0 - np.sum(c["nu"].real ** 2)
D = np.sum(c["nu"].imag ** 2) - np.sum(c["nuE"].imag ** 2)
print(f"   ||K||^2 = {c['K2']:.6f}, L = {L:.6f}, defect(y^2) = {D:.3e} (= ||K||^2 - L = {c['K2'] - L:.3e})")
print(f"time {time.time()-t0:.0f}s")
