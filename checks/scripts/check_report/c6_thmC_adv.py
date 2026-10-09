"""Adversarial test of Theorem B (inequality) and Theorem C (ratio W(lam)/bound(lam) >= 1), general B.
Exact left/right stop-loss via the dilogarithm: h_lam(x+iy) = -(1/pi^2) Im Li2(e^{pi(y-lam)} e^{-i pi x}),
Li2(z) = spence(1-z) (scipy, complex). B = diag(levels) WLOG."""
import sys
import time
import numpy as np
from scipy.special import spence
from scipy.optimize import minimize
import mpmath as mp

seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
budget = float(sys.argv[2]) if len(sys.argv) > 2 else 240
rng = np.random.default_rng(seed)
LAMS = np.linspace(-4, 4, 33)


def hvec(nus, lams):
    x = np.clip(nus.real, 0.0, 1.0)[:, None]
    y = nus.imag[:, None]
    lam = lams[None, :]
    out = np.zeros((len(nus), len(lams)))
    inner = (x > 1e-13) & (x < 1 - 1e-13)
    z = np.exp(np.pi * (y - lam)) * np.exp(-1j * np.pi * x)
    val = -np.imag(spence(1 - z)) / np.pi ** 2
    out = np.where(inner, val, 0.0)
    out = np.where(x <= 1e-13, np.maximum(y - lam, 0.0), out)
    return out.sum(axis=0)


# validate spence-based h against mpmath
for _ in range(5):
    x, y, lam = rng.uniform(0.01, 0.99), rng.normal() * 2, rng.normal() * 2
    z = mp.exp(mp.pi * (y - lam)) * mp.exp(-1j * mp.pi * x)
    ref = float(-mp.im(mp.polylog(2, z)) / mp.pi ** 2)
    assert abs(hvec(np.array([complex(x, y)]), np.array([lam]))[0] - ref) < 1e-12, "spence branch mismatch"
print("spence-based h_lambda validated against mpmath", flush=True)


def build(params, M, snap):
    lv = 1 / (1 + np.exp(-params[:M]))
    if snap == 1:
        lv = np.where(lv < 0.2, 0.0, np.where(lv > 0.8, 1.0, lv))
    elif snap == 2:
        lv = np.round(lv * 3) / 3          # repeated levels in {0,1/3,2/3,1}
    G = params[M:M + M * M].reshape(M, M)
    g = (G + G.T) / 2 + 1j * (np.triu(G, 1) - np.triu(G, 1).T)
    return lv, g


def quantities(lv, g):
    M = len(lv)
    B = np.diag(lv).astype(complex)
    pd = 1 - lv
    # exact eigenspaces of B (rtol=0): the pinching and the bound must use the same partition
    mask = np.abs(pd[:, None] - pd[None, :]) <= 1e-12
    gd = g * mask
    nus = np.linalg.eigvals(B + 1j * g)
    mus = np.linalg.eigvals(B + 1j * gd)

    def stable_W(a, b):
        # sum h_lam(a) - sum h_lam(b); for lam below the centre use the mirrored stop-loss
        # int (lam-w)_+ = h_{-lam}(conj nu), valid since masses and first moments agree (no cancellation)
        ctr = np.real(np.trace(g)) / M
        up = LAMS >= ctr
        out = np.empty_like(LAMS)
        out[up] = hvec(a, LAMS[up]) - hvec(b, LAMS[up])
        out[~up] = hvec(a.conj(), -LAMS[~up]) - hvec(b.conj(), -LAMS[~up])
        return out
    WL = stable_W(nus, mus)
    WR = stable_W(1 - nus.conj(), 1 - mus.conj())
    # Theorem C bound (left edge): sum_{p_j<p_k} |g_jk|^2 avg_{[p_j,p_k]} sin(pi tau) / (2(1+cosh(pi R)))
    S = 0.0
    for i in range(M):
        for j in range(M):
            if pd[i] < pd[j] - 1e-12:
                avg = (np.cos(np.pi * pd[i]) - np.cos(np.pi * pd[j])) / np.pi / (pd[j] - pd[i])
                S += abs(g[i, j]) ** 2 * avg
    R = np.abs(LAMS) + np.linalg.norm(g, 2)
    bound = S / (2 * (1 + np.cosh(np.pi * R)))
    return WL, WR, bound, S


def objB(params, M, snap):
    lv, g = build(params, M, snap)
    WL, WR, _, _ = quantities(lv, g)
    return min(WL.min(), WR.min()) / (1 + np.linalg.norm(g, 2))


def objC(params, M, snap):
    lv, g = build(params, M, snap)
    WL, _, bound, S = quantities(lv, g)
    ok = bound > 1e-12
    if not ok.any() or S < 1e-6:
        return 10.0
    return np.min(WL[ok] / bound[ok])


t0 = time.time()
bestB, bestC = np.inf, np.inf
rounds = 0
while time.time() - t0 < budget:
    M = int(rng.integers(2, 6))
    snap = int(rng.integers(0, 3))
    x0 = rng.normal(size=M + M * M) * rng.uniform(0.3, 1.5)
    rB = minimize(objB, x0, args=(M, snap), method="Nelder-Mead", options={"maxiter": 1500, "fatol": 1e-14})
    rC = minimize(objC, x0, args=(M, snap), method="Nelder-Mead", options={"maxiter": 1500, "fatol": 1e-12})
    rounds += 1
    if rB.fun < bestB:
        bestB = rB.fun
    if rC.fun < bestC:
        bestC = rC.fun
        lv, g = build(rC.x, M, snap)
        print(f"[{time.time()-t0:.0f}s] round {rounds} M={M} snap={snap}: new min ratio W/boundC = {bestC:.6f} "
              f"(levels {np.round(lv, 3)}, ||g||={np.linalg.norm(g, 2):.2f})", flush=True)
print(f"DONE seed {seed}: {rounds} rounds; min normalised Theorem-B gap {bestB:+.3e}; min Theorem-C ratio {bestC:.6f}")
