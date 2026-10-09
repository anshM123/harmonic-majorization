"""Check S5b (replaces Part 3 of s5_halfplane.py, whose direct mp.quad evaluation of the density difference at
huge |w| lost all digits to cancellation).

The density of pi(X) - pi(X_E) (continuous part) is f(w) = (1/pi) Im sum_nu [1/(w - mu_nu)] - (same for X_E),
mu_nu = y + i x = i conj(nu).  For |w| > W0 > max|mu| we use the convergent expansion
   f(w) = (1/pi) sum_{k>=0} D_k / w^{k+1},   D_k = Im( sum_X mu^k - sum_XE mu^k ),
(D_0 = D_1 = D_2 = 0: equal masses, equal sums of x and of x*y), and the exact rational form for |w| <= W0 at 60 digits.
Checks:  int w^2 d(pi(X)-pi(X_E)) = ||K||_F^2;  int phi d(...) >= 0 for convex phi of growth |w|^p, p < 3;
D_3 = 6 int int beta rho_B > 0 (tail ~ D_3/(pi w^4)).
"""
import time
import numpy as np
import mpmath as mp
from hmaudit import rand_herm, setup

rng = np.random.default_rng(22)
t0 = time.time()
mp.mp.dps = 60


def mp_spec(levels, g):
    M = len(levels)
    X = mp.matrix(M, M)
    for i in range(M):
        for j in range(M):
            X[i, j] = (mp.mpf(levels[i]) if i == j else 0) + 1j * mp.mpc(g[i, j].real, g[i, j].imag)
    ev = mp.eig(X)[0]
    lab = np.unique(levels, return_inverse=True)[1]
    evE = []
    for L in np.unique(lab):
        idx = np.where(lab == L)[0]
        H = mp.matrix(len(idx), len(idx))
        for a, i in enumerate(idx):
            for b, j in enumerate(idx):
                H[a, b] = mp.mpc(g[i, j].real, g[i, j].imag)
        evE += [mp.mpf(levels[idx[0]]) + 1j * mp.re(v) for v in mp.eigh(H)[0]]
    ev = [mp.mpc(max(mp.re(z), 0), mp.im(z)) for z in ev]
    return ev, evE


for case, (levels, gs) in enumerate([([0.2, 1.1, 1.1], 0.7), ([0.0, 0.5, 2.0, 2.0], 1.5), ([0.0, 0.0, 3.0], 1.0),
                                     ([0.4, 0.9, 1.7, 4.0], 2.2), ([0.05, 0.05, 0.3], 0.4)]):
    g = rand_herm(len(levels), rng, gs)
    c = setup(levels, g)
    ev, evE = mp_spec(np.array(levels), g)
    mus = [1j * mp.conj(z) for z in ev]
    musE = [1j * mp.conj(z) for z in evE]
    Kmax = 60
    D = [mp.im(mp.fsum(m ** k for m in mus) - mp.fsum(m ** k for m in musE)) for k in range(Kmax)]
    R = max(abs(m) for m in mus + musE)
    W0 = 4 * (R + 1)

    def f(w):
        if abs(w) > W0:
            return mp.fsum(D[k] / w ** (k + 1) for k in range(3, Kmax)) / mp.pi
        s = mp.fsum(x / ((w - y) ** 2 + x ** 2) for x, y in [(mp.re(z), mp.im(z)) for z in ev] if x > 0)
        s -= mp.fsum(x / ((w - y) ** 2 + x ** 2) for x, y in [(mp.re(z), mp.im(z)) for z in evE] if x > 0)
        return s / mp.pi

    def atoms(phi):
        return mp.fsum(phi(mp.im(z)) for z in ev if mp.re(z) == 0) - mp.fsum(phi(mp.im(z)) for z in evE if mp.re(z) == 0)

    ys = sorted(set([float(mp.im(z)) for z in ev + evE]))
    brk = [-mp.inf, -W0] + [mp.mpf(v) for v in ys] + [W0, mp.inf]
    m0 = mp.quad(f, brk) + atoms(lambda y: 1)
    m1 = mp.quad(lambda w: w * f(w), brk) + atoms(lambda y: y)
    m2 = mp.quad(lambda w: w ** 2 * f(w), brk) + atoms(lambda y: y ** 2)
    print(f"case {case}: levels {levels}; D0..D2 = {[mp.nstr(d, 3) for d in D[:3]]}; D3 = {mp.nstr(D[3], 12)}")
    print(f"   mass {mp.nstr(m0, 5)}, first moment {mp.nstr(m1, 5)}, second moment {mp.nstr(m2, 20)} vs ||K||^2 = {c['K2']:.17f}")
    for name, phi, bp in [("|w-0.3|^2.5", lambda w: abs(w - mp.mpf("0.3")) ** mp.mpf("2.5"), [mp.mpf("0.3")]),
                          ("(w+0.5)_+^2", lambda w: max(w + mp.mpf("0.5"), 0) ** 2, [mp.mpf("-0.5")]),
                          ("|w|^2.9", lambda w: abs(w) ** mp.mpf("2.9"), [mp.mpf(0)]),
                          ("cosh-type sqrt(1+w^2)", lambda w: mp.sqrt(1 + w ** 2), [])]:
        b2 = sorted(set(brk[1:-1] + bp))
        val = mp.quad(lambda w: phi(w) * f(w), [-mp.inf] + b2 + [mp.inf]) + atoms(phi)
        print(f"   int {name} d(pi(X)-pi(X_E)) = {mp.nstr(val, 15)}  (>= 0: {val >= 0})")
print(f"time {time.time()-t0:.0f}s")
