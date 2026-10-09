"""Theorem B identity W(lam) = V(lam), general B, careful quadrature for V.
Inner s-integral: breakpoints at the transition points of the root configuration (sqrt kinks of rho, located by
bisection on the count of non-real roots), geometric grading around lam down to the kernel width x, cosine-mapped
Gauss panels. Outer tau-integral: adaptive quad in the cosine variable on each gap of spec P."""
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from hm_common import make_case, Kx, W_left

seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
lvset = int(sys.argv[2]) if len(sys.argv) > 2 else 0
sides = sys.argv[3] if len(sys.argv) > 3 else "LR"
rng = np.random.default_rng(seed)
XG, WG = leggauss(20)
TH = (XG + 1) * np.pi / 2
WTH = WG * np.pi / 2


def roots_batch(c, S, tau):
    M = c["M"]
    A = (c["g"][None, :, :] - S[:, None, None] * np.eye(M)[None]) / (c["pdiag"] - tau)[None, :, None]
    return np.linalg.eigvals(A)


def nonreal_count(c, S, tau):
    xi = roots_batch(c, S, tau)
    return np.sum(np.abs(xi.imag) > 1e-9 * (1 + np.abs(xi)), axis=1)


def transitions(c, tau, lo, hi, n=2001):
    S = np.linspace(lo, hi, n)
    cnt = nonreal_count(c, S, tau)
    ch = np.where(cnt[1:] != cnt[:-1])[0]
    if len(ch) == 0:
        return np.array([])
    a, b = S[ch].copy(), S[ch + 1].copy()
    ca = cnt[ch]
    for _ in range(45):
        mid = (a + b) / 2
        cm = nonreal_count(c, mid, tau)
        left = cm == ca
        a = np.where(left, mid, a)
        b = np.where(left, b, mid)
    return (a + b) / 2


def F(c, tau, lam, x):
    lo, hi = c["lmin"], c["lmax"]
    bps = [lo, hi]
    bps += list(transitions(c, tau, lo, hi))
    if lo < lam < hi:
        bps.append(lam)
    for j in range(-2, 60):
        d = x * 2.0 ** j
        if d > (hi - lo):
            break
        for sgn in (-1, 1):
            if lo < lam + sgn * d < hi:
                bps.append(lam + sgn * d)
    bps += list(np.linspace(lo, hi, 41))
    bps = np.unique(np.array(bps))
    s0, s1 = bps[:-1], bps[1:]
    S = (s0[:, None] + (s1 - s0)[:, None] * (1 - np.cos(TH))[None, :] / 2).ravel()
    Wt = ((s1 - s0)[:, None] / 2 * np.sin(TH)[None, :] * WTH[None, :]).ravel()
    xi = roots_batch(c, S, tau)
    r = np.sum(np.abs(xi.imag), axis=1) / (2 * np.pi)
    return np.sum(Wt * Kx(x, S - lam) * r)


def V(c, lam, side):
    pv = c["pvals"]
    tot = 0.0
    errs = 0.0
    for k in range(len(pv) - 1):
        a, b = pv[k], pv[k + 1]
        h = b - a

        def integrand(ph):
            tau = a + h * (1 - np.cos(ph)) / 2
            if not (a < tau < b):
                return 0.0
            x = 1 - tau if side == "L" else tau
            return F(c, tau, lam, x) * h / 2 * np.sin(ph)
        v, e = quad(integrand, 0, np.pi, limit=300, epsabs=1e-11, epsrel=1e-9)
        tot += v
        errs += e
    return tot, errs


LV = [[0.0, 0.4, 1.0], [0.15, 0.15, 0.6, 0.85], [0.0, 0.5, 0.5, 0.5], [1.0, 0.3, 0.3, 0.0]][lvset]
c = make_case(LV, rng)
blocks = [np.linalg.eigvalsh(c["g"][np.ix_(I, I)]) for I in c["idx"]]
print(f"levels={LV}, pvals={c['pvals']}, lmin={c['lmin']:.4f}, lmax={c['lmax']:.4f}")
# masses and first moments of the two sweeps (left: weights 1 - x; right: weights x)
nus = np.linalg.eigvals(c["B"] + 1j * c["g"])
mus = np.linalg.eigvals(c["B"] + 1j * c["gd"])
for side in "LR":
    wt = (lambda z: 1 - np.clip(z.real, 0, 1)) if side == "L" else (lambda z: np.clip(z.real, 0, 1))
    tr = (c["M"] - np.trace(c["B"]).real) if side == "L" else np.trace(c["B"]).real
    mom = np.trace((c["P"] if side == "L" else c["B"]) @ c["g"]).real
    print(f"[{side}] mass {np.sum(wt(nus)):.12f} / {np.sum(wt(mus)):.12f} (expected {tr:.12f}); first moment "
          f"{np.sum(wt(nus) * nus.imag):.12f} / {np.sum(wt(mus) * mus.imag):.12f} (expected {mom:.12f})", flush=True)
if len(sys.argv) > 4 and sys.argv[4] == "moments":
    sys.exit(0)
print("   block eigenvalues (tau=p_1 block first):", [np.round(bk, 6) for bk in blocks], flush=True)
for side in sides:
    crit = blocks[-1] if side == "L" else blocks[0]      # where rho concentrates against a narrowing kernel
    lams = list(np.linspace(c["lmin"] - 0.3, c["lmax"] + 0.3, 4)) + list(crit) + list(crit + 1e-3)
    for lam in lams:
        w = W_left(c, lam, side)
        v, e = V(c, lam, side)
        print(f"[{side}] lam={lam:+.6f}: W={w:.10f}  V={v:.10f}  (quad err est {e:.1e})  W-V={w - v:+.2e}", flush=True)
