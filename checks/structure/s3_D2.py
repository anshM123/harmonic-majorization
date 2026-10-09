"""Check S3: Theorem D2 (two-level superharmonic majorization) and the exact two-level identity (2.5) of
03-structure.md; validation of the sweep closed form h_lambda; Theorem B identity W = V (cross-check).

Part 0: closed form h_lambda (hmaudit.h_closed_mp) vs direct quadrature of the strip Poisson kernel.
Part 1: W(lambda) = V(lambda) (Theorem B) on one 3-level configuration with eigenvalues 0 and 1 (quadrature of rho).
Part 2: D2 on random two-level pairs (any multiplicities, any b1 < b2), test functions in S([b1,b2]):
        non-separable psi(x) sqrt(s0^2+(y-y0)^2) (psi'' <= -psi/s0^2, psi >= 0, incl. a concave-kinked psi),
        -phi(x) (Ky Fan type), separable superharmonic y-convex, harmonic sweeps of convex data (Theorem B type).
Part 3: adversarial minimisation of the D2 defect over g (Nelder-Mead), several u.
Part 4: the exact identity (2.5): sum_X u - sum_XE u = int u_yy dR^flat + int (-Lap u) T dA, with
        int u_yy dR^flat evaluated through Theorem B (W from eigenvalues) and T by the sine series of the strip
        Green function G(nu,z) = sum_n sin(n pi x_nu) sin(n pi x) e^{-n pi |y - y_nu|}/(n pi).
Usage: python s3_D2.py [seed] [nrandom]
"""
import os
import sys
import time
import numpy as np


def logcosh(t):
    t = np.abs(t)
    return t + np.log1p(np.exp(-2 * t)) - np.log(2.0)

import mpmath as mp
from scipy.optimize import minimize
from numpy.polynomial.laguerre import laggauss
from hmaudit import (rand_herm, setup, h_closed_mp, h_quad_mp, adaptive_nodes_1d, _breakpoints_s, rho_vals)

seed = int(sys.argv[1]) if len(sys.argv) > 1 else 7
NR = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
rng = np.random.default_rng(seed)
t0 = time.time()
mp.mp.dps = 30

# ------------------------------------------------------------------------------------------------------------
print("== Part 0: closed form of h_lambda vs quadrature")
worst = 0
for (x, y, lam) in [(0.3, 0.2, 0.0), (0.01, -1.0, 0.5), (0.97, 2.0, -1.0), (0.5, 0.0, 0.0), (1e-4, 0.3, 0.3),
                    (0.999, -0.7, 1.4), (0.25, 3.0, -2.0), (0.6, -2.5, 1.0)]:
    a = h_closed_mp(mp.mpc(x, y), lam)
    b = h_quad_mp(mp.mpc(x, y), lam)
    worst = max(worst, abs(a - b))
    print(f"   x={x}, y={y}, lam={lam}: closed {mp.nstr(a, 20)}  quad {mp.nstr(b, 20)}  |diff| {mp.nstr(abs(a-b), 3)}")
print(f"   worst |diff| = {mp.nstr(worst, 3)}")


def W_left(nus, lam):
    return mp.fsum(h_closed_mp(mp.mpc(z), lam) for z in nus)


def W_right(nus, lam):
    return mp.fsum(h_closed_mp(1 - mp.conj(mp.mpc(z)), lam) for z in nus)


# ------------------------------------------------------------------------------------------------------------
print("\n== Part 1: Theorem B identity W(lambda) = V(lambda) (left and right edges), levels {0, 0.4, 0.4, 1}")
levels = [0.0, 0.4, 0.4, 1.0]
g = rand_herm(4, rng, 0.8)
c = setup(levels, g)


def Kx(x, u):
    # strip Poisson kernel for the left edge, cancellation-free form
    return np.sin(np.pi * x) / (2 * (2 * np.sinh(np.pi * u / 2) ** 2 + 2 * np.sin(np.pi * x / 2) ** 2))


def balayage(c, lam, side):
    """V(lam) = int int K_x(s - lam) rho_B(s, x) ds dx (left: x = beta; right: x = 1 - beta)."""
    def inner(beta):
        xk = beta if side == "L" else 1 - beta
        pts = _breakpoints_s(c, beta, 600)
        if c["lmin"] < lam < c["lmax"]:
            pts = sorted(set(pts + [lam]))
        tot = 0.0
        for a, b in zip(pts[:-1], pts[1:]):
            f = lambda s: rho_vals(c, s, beta) * Kx(xk, s - lam)
            ss, sw = adaptive_nodes_1d(f, a, b, tol=1e-11, maxit=60)
            tot += np.sum(sw * f(ss))
        return tot
    tot = 0.0
    u = c["uniq"]
    for k in range(len(u) - 1):
        fv = lambda bs: np.array([inner(b) for b in bs])
        bx, bw = adaptive_nodes_1d(fv, u[k], u[k + 1], tol=1e-9, nmin=4, maxit=30)
        tot += np.sum(bw * fv(bx))
    return tot


for side in "LR":
    for lam in [-1.0, 0.0, 0.37, float(np.real(np.linalg.eigvalsh(g[1:3, 1:3])[0]))]:
        Wv = (W_left(c["nu"], lam) - W_left(c["nuE"], lam)) if side == "L" else \
            (W_right(c["nu"], lam) - W_right(c["nuE"], lam))
        Vv = balayage(c, lam, side)
        print(f"   {side} lam={lam:+.6f}: W = {mp.nstr(Wv, 14)}   V = {Vv:.14f}   |W-V| = {abs(float(Wv)-Vv):.1e}")
print(f"   [elapsed {time.time()-t0:.0f}s]")

# ------------------------------------------------------------------------------------------------------------
# test functions in S([b1, b2])
def make_tests(b1, b2, rng):
    w = b2 - b1
    cmid = 0.5 * (b1 + b2)
    T = []
    for _ in range(3):
        k = rng.uniform(0.2, 1.0) * np.pi / w            # |x - cmid| <= w/2 <= pi/(2k)
        s0 = rng.uniform(1.0, 3.0) / k                      # k >= 1/s0  => superharmonic
        y0 = rng.normal()
        T.append((f"cos*sqrt(k={k:.2f},s0={s0:.2f})",
                  lambda z, k=k, s0=s0, y0=y0: np.cos(k * (z.real - cmid)) * np.sqrt(s0 ** 2 + (z.imag - y0) ** 2)))
    # concave-kinked psi = min(cos k(x-c1), cos k(x-c2)) >= 0 on the strip
    k = 0.5 * np.pi / w
    c1, c2 = b1 + 0.3 * w, b1 + 0.7 * w
    s0 = 2.0 / k
    T.append(("min-cos*sqrt (kinked)", lambda z: np.minimum(np.cos(k * (z.real - c1)), np.cos(k * (z.real - c2)))
              * np.sqrt(s0 ** 2 + (z.imag - 0.1) ** 2)))
    x0 = b1 + rng.uniform() * w
    T.append(("-(x-x0)^2", lambda z: -(z.real - x0) ** 2))
    T.append(("-|x-x0|", lambda z: -np.abs(z.real - x0)))
    kap, s = rng.uniform(0.5, 3), rng.uniform(0.3, 2)
    T.append(("sep: -k x^2/2 + k s^2 logcosh", lambda z: -kap * z.real ** 2 / 2
              + kap * s ** 2 * logcosh(((z.imag - 0.2) / s))))
    return T


def sweep_tests(b1, b2, nuX, nuE, lams):
    """harmonic u: Poisson extension of (y - lam)_+ on the left or right edge (strip normalised)."""
    w = b2 - b1
    out = []
    for lam in lams:
        for side in "LR":
            f = W_left if side == "L" else W_right
            nx = [(z - b1) / w for z in nuX]
            ne = [(z - b1) / w for z in nuE]
            out.append((f"sweep {side} lam={lam:+.2f}", float(f(nx, lam) - f(ne, lam))))
    return out


print("\n== Part 2: D2 on random two-level pairs")
mults = [(1, 1), (1, 2), (2, 1), (2, 2), (1, 3), (3, 2), (2, 3)]
worst_ratio = np.inf
worst_info = None
nneg = 0
mp.mp.dps = 20
for trial in range(NR):
    r1, r2 = mults[rng.integers(len(mults))]
    b1 = rng.normal()
    b2 = b1 + rng.uniform(0.05, 3.0)
    levels = [b1] * r1 + [b2] * r2
    M = r1 + r2
    mode = rng.integers(4)
    if mode == 0:
        g = rand_herm(M, rng, rng.uniform(0.05, 4.0))
    elif mode == 1:   # weak coupling (near commuting)
        g = np.diag(rng.normal(size=M)).astype(complex) + rand_herm(M, rng, 1e-3)
    elif mode == 2:   # near exceptional points: equal diagonals, coupling ~ half the gap
        g = np.full((M, M), 0j)
        g += np.diag(np.full(M, rng.normal()))
        for i in range(r1):
            for j in range(r1, M):
                v = (b2 - b1) / 2 * rng.uniform(0.9, 1.1) * np.exp(1j * rng.uniform(0, 2 * np.pi)) / np.sqrt(r1 * r2)
                g[i, j] = v
                g[j, i] = np.conj(v)
    else:             # strong coupling
        g = rand_herm(M, rng, rng.uniform(5, 30))
    c = setup(levels, g)
    if c["K2"] < 1e-14:
        continue
    for name, u in make_tests(b1, b2, rng):
        D = np.sum(u(c["nu"])) - np.sum(u(c["nuE"]))
        ratio = D / c["K2"]
        if ratio < worst_ratio:
            worst_ratio, worst_info = ratio, (trial, name, levels, mode, D, c["K2"])
        if D < -1e-10 * (1 + np.sum(np.abs(u(c["nu"])))):
            nneg += 1
    if trial % 10 == 0:
        for name, D in sweep_tests(b1, b2, c["nu"], c["nuE"], [rng.normal()]):
            ratio = D / c["K2"]
            if ratio < worst_ratio:
                worst_ratio, worst_info = ratio, (trial, name, levels, mode, D, c["K2"])
            if D < -1e-12:
                nneg += 1
print(f"   {NR} random pairs x 7 test functions (+ harmonic sweeps every 10th): significant negatives = {nneg}")
print(f"   min D/||K||^2 = {worst_ratio:.3e}  at {worst_info}")
print(f"   [elapsed {time.time()-t0:.0f}s]")

# ------------------------------------------------------------------------------------------------------------
print("\n== Part 3: adversarial minimisation of D/||K||^2 over g (two levels)")


def gfrom(p, M):
    A = p[:M * M].reshape(M, M) + 1j * p[M * M:].reshape(M, M)
    return (A + A.conj().T) / 2


advmin = np.inf
for rnd in range(24):
    r1, r2 = mults[rnd % len(mults)]
    b1, b2 = 0.0, 1.0
    levels = [b1] * r1 + [b2] * r2
    M = r1 + r2
    tests = make_tests(b1, b2, rng)
    name, u = tests[rnd % len(tests)]

    def obj(p):
        g = gfrom(p, M)
        c = setup(levels, g)
        if c["K2"] < 1e-8:
            return 1.0
        return (np.sum(u(c["nu"])) - np.sum(u(c["nuE"]))) / c["K2"]

    best = np.inf
    for start in range(3):
        p0 = rng.normal(size=2 * M * M) * rng.uniform(0.1, 3)
        res = minimize(obj, p0, method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-10, fatol=1e-14))
        best = min(best, res.fun)
    advmin = min(advmin, best)
    print(f"   round {rnd:2d}: mult {r1}+{r2}, u = {name:32s}: min D/||K||^2 = {best:.4e}")
print(f"   adversarial minimum over all rounds: {advmin:.4e}")
print(f"   [elapsed {time.time()-t0:.0f}s]")

# ------------------------------------------------------------------------------------------------------------
print("\n== Part 4: exact two-level identity (2.5) for u = cos(k(x-c)) sqrt(s0^2+(y-y0)^2)")
mp.mp.dps = 25
xL, wL = laggauss(80)
for trial, (r1, r2, b1, b2, gs) in enumerate([(1, 1, 0.0, 1.0, 0.7), (2, 1, -0.4, 1.1, 1.3), (2, 2, 0.5, 2.5, 2.0),
                                               (1, 3, 0.0, 0.6, 0.4)]):
    levels = [b1] * r1 + [b2] * r2
    M = r1 + r2
    g = rand_herm(M, rng, gs)
    c = setup(levels, g)
    w = b2 - b1
    cm = 0.5 * (b1 + b2) + 0.1 * w
    k = 0.8 * np.pi / w / 1.2
    s0 = 1.5 / k
    y0 = 0.3
    u = lambda x, y: np.cos(k * (x - cm)) * np.sqrt(s0 ** 2 + (y - y0) ** 2)
    uyy = lambda x, y: np.cos(k * (x - cm)) * s0 ** 2 / (s0 ** 2 + (y - y0) ** 2) ** 1.5
    LHS = np.sum(u(c["nu"].real, c["nu"].imag)) - np.sum(u(c["nuE"].real, c["nuE"].imag))
    # normalised spectra
    nX = [(z - b1) / w for z in c["nu"]]
    nE = [(z - b1) / w for z in c["nuE"]]
    # RHS1 = int w^2 u_yy(b1, w y') W^L(y') dy' + int w^2 u_yy(b2, w y') W^R(y') dy'
    def integrand(yp, side):
        if side == "L":
            Wv = W_left(nX, yp) - W_left(nE, yp)
            return w ** 2 * uyy(b1, w * float(yp)) * Wv
        Wv = W_right(nX, yp) - W_right(nE, yp)
        return w ** 2 * uyy(b2, w * float(yp)) * Wv
    ys = sorted([float(z.imag) for z in nX + nE])
    brk = [-mp.inf] + [mp.mpf(v) for v in ys] + [mp.inf]
    RHS1 = mp.quad(lambda t: integrand(t, "L"), brk) + mp.quad(lambda t: integrand(t, "R"), brk)
    # RHS2 = w^2 sum_nu sum_n (1/(n pi)) sin(n pi x'_nu) A_n B_n(nu)
    N = 20000
    n = np.arange(1, N + 1)
    alpha, betaa = k * (b1 - cm), k * w
    pn = n * np.pi
    A_n = 0.5 * ((np.cos(alpha) - np.cos(pn + betaa + alpha)) / (pn + betaa)
                 + (np.cos(alpha) - np.cos(pn - betaa - alpha)) / (pn - betaa))
    RHS2 = 0.0
    for z in nX:
        xv, yv = z.real, z.imag
        # B_n = int f(y') e^{-n pi |y'-yv|} dy', f(y') = k^2 S(w y') - s0^2 S(w y')^-3, S(y) = sqrt(s0^2+(y-y0)^2)
        f = lambda yp: k ** 2 * np.sqrt(s0 ** 2 + (w * yp - y0) ** 2) - s0 ** 2 / (s0 ** 2 + (w * yp - y0) ** 2) ** 1.5
        tt = xL[None, :] / pn[:, None]
        Bn = np.sum(wL[None, :] * (f(yv + tt) + f(yv - tt)), axis=1) / pn
        RHS2 += w ** 2 * np.sum(np.sin(pn * xv) * A_n * Bn / pn)
    print(f"   case {trial}: mult {r1}+{r2}, [{b1},{b2}]: LHS = {LHS:.12f}   int u_yy dR^flat = {float(RHS1):.12f}   "
          f"int(-Lap u)T = {RHS2:.12f}   RHS = {float(RHS1)+RHS2:.12f}   |diff| = {abs(LHS-float(RHS1)-RHS2):.1e}")
print(f"time {time.time()-t0:.0f}s")
