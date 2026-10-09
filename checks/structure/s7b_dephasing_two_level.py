"""Check S7b: two-level dephasing line, adversarial search with 60-digit confirmation of every candidate whose
double-precision normalised increment is below -1e-12.  (In s7_dephasing.py the 1+1 adversary reached -1.6e-8 in double
precision; this script decides whether such values are genuine.)
Objective: min over lambda in a grid and both edges of [W_{gamma2}(lambda) - W_{gamma1}(lambda)] / ||K||^2, gamma1 < gamma2.
Usage: python s7b_dephasing_two_level.py [seed] [rounds]
"""
import sys
import time
import numpy as np
import mpmath as mp
from scipy.special import spence
from scipy.optimize import minimize
from hmaudit import h_closed_mp, labels_of, pinch

seed = int(sys.argv[1]) if len(sys.argv) > 1 else 51
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 60
rng = np.random.default_rng(seed)
t0 = time.time()
lams = np.linspace(-6, 6, 241)


def h_fast(nu, lam):
    x = np.clip(nu.real, 0, 1)[None, :]
    y = nu.imag[None, :]
    lam = np.asarray(lam)[:, None]
    t = y - lam
    q = np.exp(-np.pi * np.abs(t)) * np.exp(1j * np.pi * x)
    val = np.imag(spence(1 - q)) / np.pi ** 2 + np.maximum(t, 0) * (1 - x)
    val = np.where(x <= 0, np.maximum(t, 0), val)
    return np.where(x >= 1, 0.0, val)


def stoploss(nus, side):
    z = nus if side == "L" else 1 - np.conj(nus)
    return h_fast(z, lams).sum(axis=1)


def gfrom(p, M):
    A = p[:M * M].reshape(M, M) + 1j * p[M * M:].reshape(M, M)
    return (A + A.conj().T) / 2


def gam(p):
    g1 = 0.5 * (1 + np.tanh(p[-2] / 2))
    g2 = g1 + (1 - g1) * 0.5 * (1 + np.tanh(p[-1] / 2))
    return g1, g2


def mp_increment(levels, g, g1, g2):
    lab, _ = labels_of(levels)
    gE = pinch(g, lab)
    M = len(levels)
    def eig(gg):
        X = mp.matrix(M, M)
        for i in range(M):
            for j in range(M):
                X[i, j] = (mp.mpf(levels[i]) if i == j else 0) + 1j * mp.mpc(gg[i, j].real, gg[i, j].imag)
        return mp.eig(X)[0]
    ea, eb = eig(gE + g1 * (g - gE)), eig(gE + g2 * (g - gE))
    K2 = mp.mpf(np.sum(np.abs(g - gE) ** 2))
    best = (mp.inf, None)
    for side in "LR":
        A = [z if side == "L" else 1 - mp.conj(z) for z in eb]
        Bq = [z if side == "L" else 1 - mp.conj(z) for z in ea]
        for lam in lams[::4]:
            v = (mp.fsum(h_closed_mp(z, lam) for z in A) - mp.fsum(h_closed_mp(z, lam) for z in Bq)) / K2
            if v < best[0]:
                best = (v, (side, float(lam)))
    return best


mp.mp.dps = 60
summary = []
for rnd in range(ROUNDS):
    r1, r2 = [(1, 1), (1, 2), (2, 2), (1, 3), (2, 3)][rnd % 5]
    levels = np.array([0.0] * r1 + [1.0] * r2)
    M = r1 + r2
    lab, _ = labels_of(levels)

    def obj(p):
        g = gfrom(p[:-2], M)
        K2 = np.sum(np.abs(g - pinch(g, lab)) ** 2)
        if K2 < 1e-6:
            return 1.0
        g1, g2 = gam(p)
        gE = pinch(g, lab)
        n1 = np.linalg.eigvals(np.diag(levels) + 1j * (gE + g1 * (g - gE)))
        n2 = np.linalg.eigvals(np.diag(levels) + 1j * (gE + g2 * (g - gE)))
        return min((stoploss(n2, s) - stoploss(n1, s)).min() for s in "LR") / K2

    p0 = np.concatenate([rng.normal(size=2 * M * M) * rng.uniform(0.3, 6), rng.normal(size=2)])
    res = minimize(obj, p0, method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-10, fatol=1e-14))
    line = f"round {rnd:2d} mult {r1}+{r2}: double-precision min increment {res.fun:.3e}"
    if res.fun < -1e-12:
        g = gfrom(res.x[:-2], M)
        g1, g2 = gam(res.x)
        v, where = mp_increment(levels, g, g1, g2)
        line += f"  -> 60-digit: {mp.nstr(v, 6)} at {where} (gamma {g1:.6f} -> {g2:.6f})"
        summary.append(float(v))
    print(line, flush=True)
print(f"60-digit confirmed minima of candidates: {sorted(summary)[:5]}")
print(f"time {time.time()-t0:.0f}s")
