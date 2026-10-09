"""Check S1: Theorem A' (holomorphic form), its marginals, the harmonic corollary and Theorem 1 of
03-structure.md, for GENERAL Hermitian B (negative, repeated and >1 eigenvalues; Lean covers only 0 <= P <= 1).
The quadrature of rho_B (adaptive in s and beta, hmaudit.rho_nodes) is independent of the eigenvalue side.
Usage: python s1_thmAprime.py [seed]
"""
import sys
import time
import numpy as np
from hmaudit import (rand_herm, setup, rho_nodes, rho_vals, rho_pairs, mu_B, adaptive_nodes_1d,
                     _breakpoints_s)

seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
rng = np.random.default_rng(seed)
t0 = time.time()

configs = [
    ("3 levels incl. negative, simple", [-1.3, 0.4, 2.1]),
    ("4 levels with a repeat", [-0.5, -0.5, 0.8, 1.7]),
    ("5 levels, two repeats", [0.0, 0.0, 1.0, 1.0, 3.0]),
    ("6 levels, interior repeats, large range", [-2.0, -0.3, -0.3, 0.9, 0.9, 4.0]),
    ("projection-like (2 levels, mult 2+2)", [0.0, 0.0, 1.0, 1.0]),
]
worst = {}


def upd(key, val):
    worst[key] = max(worst.get(key, 0.0), val)


for name, levels in configs:
    M = len(levels)
    g = rand_herm(M, rng, scale=float(rng.uniform(0.5, 2.0)))
    if rng.uniform() < 0.5:
        g = g + np.diag(rng.normal(size=M))  # shift diagonal blocks
    c = setup(levels, g)
    S, BT, W = rho_nodes(c)
    zeta = BT + 1j * S
    print(f"== {name}: levels {levels}, ||K||_F^2 = {c['K2']:.6f}, nodes {len(W)}")
    mass = W.sum()
    rel = abs(mass - c["K2"] / 2) / (c["K2"] / 2)
    upd("mass", rel)
    print(f"   mass of R: {mass:.13f}  vs ||K||^2/2 = {c['K2']/2:.13f}   rel {rel:.2e}")
    inside = (S >= c["lmin"] - 1e-12).all() and (S <= c["lmax"] + 1e-12).all() and \
        (BT >= min(levels)).all() and (BT <= max(levels)).all()
    print(f"   all nodes inside Rect: {inside}")
    nu, nuE = c["nu"], c["nuE"]
    xr = (min(levels), max(levels))
    tests = []
    for kap in [1.0, 2 + 1j, -1.5j, 0.7 - 0.3j]:
        tests.append((f"exp({kap}z)", lambda z, k=kap: np.exp(k * z), lambda z, k=kap: k ** 2 * np.exp(k * z)))
    z0 = xr[1] + 0.7 + 0.3j  # right of Rect
    tests.append(("1/(z0-z)", lambda z: 1 / (z0 - z), lambda z: 2 / (z0 - z) ** 3))
    tests.append(("log(z0-z)", lambda z: np.log(z0 - z), lambda z: -1 / (z0 - z) ** 2))
    z1 = 0.5 * (xr[0] + xr[1]) + 1j * (c["lmax"] + 0.9)  # above Rect
    tests.append(("1/(z1-z)^2", lambda z: 1 / (z1 - z) ** 2, lambda z: 6 / (z1 - z) ** 4))
    for n in [2, 3, 5, 8]:
        tests.append((f"z^{n}", lambda z, n=n: z ** n, lambda z, n=n: n * (n - 1) * z ** (n - 2)))
    for tname, f, f2 in tests:
        lhs = np.sum(f(nu)) - np.sum(f(nuE))
        rhs = -np.sum(f2(zeta) * W)
        err = abs(lhs - rhs) / max(1e-300, abs(lhs))
        upd("A'", err)
        print(f"   A' {tname:16s}: Tr f(X)-Tr f(X_E) = {lhs:.13g}   -int f'' dR = {rhs:.13g}   rel {err:.2e}")
    # harmonic corollary + Theorem 1 (<=): u = cos(k(x-cx)) cosh(k(y-y0)), u_yy >= 0 on |x-cx| < pi/(2k)
    cx = 0.5 * (xr[0] + xr[1])
    kk = 0.95 * np.pi / (xr[1] - xr[0])
    y0 = 0.3
    u = lambda z: np.cos(kk * (z.real - cx)) * np.cosh(kk * (z.imag - y0))
    uyy = lambda z: kk ** 2 * np.cos(kk * (z.real - cx)) * np.cosh(kk * (z.imag - y0))
    d = np.sum(u(nu)) - np.sum(u(nuE))
    q = np.sum(uyy(zeta) * W)
    err = abs(d - q) / abs(d)
    upd("harmonic", err)
    print(f"   Cor/Thm1(<=): u=cos*cosh, defect {d:.13g}  int u_yy dR {q:.13g}  rel {err:.2e}  (>=0: {d >= 0})")
    uq = lambda z: (z.imag - 0.37) ** 2 - (z.real + 0.21) ** 2
    dq = np.sum(uq(nu)) - np.sum(uq(nuE))
    upd("quadratic", abs(dq - c["K2"]) / c["K2"])
    print(f"   u=(y-y0)^2-(x-x0)^2: defect {dq:.14g} vs ||K||^2 {c['K2']:.14g}")
    # horizontal marginal at fixed beta
    for k in range(len(c["uniq"]) - 1):
        for frac in [1e-6, 0.013, 0.5, 0.977, 1 - 1e-6]:
            beta = c["uniq"][k] + frac * (c["uniq"][k + 1] - c["uniq"][k])
            pts = _breakpoints_s(c, beta, 800)
            tot = 0.0
            for a, b in zip(pts[:-1], pts[1:]):
                ss, sw = adaptive_nodes_1d(lambda x: rho_vals(c, x, beta), a, b, tol=1e-12)
                tot += np.sum(sw * rho_vals(c, ss, beta))
            m = mu_B(c, beta)
            err = abs(tot - m) / max(m, 1e-300)
            upd("x-marginal", err)
    print(f"   x-marginal (5 betas per gap incl. 1e-6 from the lines): worst rel so far {worst['x-marginal']:.2e}")
    # vertical (xi = 0) and oblique projections: int rho(w + xi*beta, beta) dbeta = N_xi(w)
    for xi in [0.0, 0.6, -1.7]:
        H = c["g"] - xi * c["B"]
        HE = c["gE"] - xi * c["B"]
        lH, lHE = np.linalg.eigvalsh(H), np.linalg.eigvalsh(HE)
        for wfrac in [0.2, 0.45, 0.7]:
            w = lH[0] + wfrac * (lH[-1] - lH[0])
            N = np.sum(np.maximum(lH - w, 0)) - np.sum(np.maximum(lHE - w, 0))
            tot = 0.0
            for k in range(len(c["uniq"]) - 1):
                lo, hi = c["uniq"][k], c["uniq"][k + 1]
                th = np.linspace(0, np.pi, 4001)
                bs = (lo + (hi - lo) * (1 - np.cos(th)) / 2)[1:-1]
                r = rho_pairs(c, w + xi * bs, bs)
                on = r > 0
                cuts = [lo] + [0.5 * (bs[i] + bs[i + 1]) for i in np.where(on[:-1] != on[1:])[0]] + [hi]
                for a, b in zip(cuts[:-1], cuts[1:]):
                    xx, ww = adaptive_nodes_1d(lambda t: rho_pairs(c, w + xi * t, t), a, b, tol=1e-11)
                    tot += np.sum(ww * rho_pairs(c, w + xi * xx, xx))
            err = abs(tot - N) / max(abs(N), 1e-12)
            upd("oblique", err)
            print(f"   projection xi={xi:+.1f}, w={w:+.4f}: int rho along line {tot:.12f}  N_xi(w) {N:.12f}  rel {err:.1e}")
    print(f"   [elapsed {time.time()-t0:.0f}s]")

print()
print("SUMMARY (worst relative errors):", {k: f"{v:.2e}" for k, v in worst.items()})
print(f"time {time.time()-t0:.0f}s")
