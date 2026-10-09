"""Replay one round of the adversarial search c6_thmC_adv.py (seed 5; round number as argv[1]) and recompute
W(lam) in 50-digit arithmetic next to the double-precision value and the Theorem C bound. Used to expose
double-precision cancellation artefacts (spurious ratios < 1)."""
import numpy as np
import mpmath as mp
from scipy.optimize import minimize

src = open("c6_thmC_adv.py", encoding="utf-8").read()
head = src.split("t0 = time.time()")[0].replace("seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0", "seed = 5")
ns = {}
exec(head, ns)
rng = ns["rng"]
import sys
target_round = int(sys.argv[1]) if len(sys.argv) > 1 else 21
for rnd in range(1, target_round + 1):
    M = int(rng.integers(2, 6))
    snap = int(rng.integers(0, 3))
    x0 = rng.normal(size=M + M * M) * rng.uniform(0.3, 1.5)
    if rnd < target_round:
        continue
    rC = minimize(ns["objC"], x0, args=(M, snap), method="Nelder-Mead", options={"maxiter": 1500, "fatol": 1e-12})
    print("round", rnd, "M", M, "snap", snap, "objC", rC.fun)
    lv, g = ns["build"](rC.x, M, snap)
    np.save("c7_case.npy", np.concatenate([lv, g.ravel().view(float)]))
    print("levels", lv)
    print("g =\n", g)
    WL, WR, bound, S = ns["quantities"](lv, g)
    LAMS = ns["LAMS"]
    j = int(np.argmin(np.where(bound > 1e-12, WL / np.where(bound > 0, bound, 1), np.inf)))
    print(f"double precision: argmin lam={LAMS[j]}, W={WL[j]:.6e}, bound={bound[j]:.6e}, ratio={WL[j]/bound[j]:.6f}")

    # high precision W
    mp.mp.dps = 50
    Bm = mp.matrix(np.diag(lv).tolist())
    gm = mp.matrix(g.tolist())
    pd = 1 - lv
    mask = np.abs(pd[:, None] - pd[None, :]) <= 1e-12
    gdm = mp.matrix((g * mask).tolist())
    I = mp.mpc(0, 1)
    nus = mp.eig(Bm + I * gm, left=False, right=False)
    mus = mp.eig(Bm + I * gdm, left=False, right=False)

    def h(nu, lam):
        x, y = mp.re(nu), mp.im(nu)
        if x <= 0:
            return max(y - lam, 0)
        if x >= 1:
            return mp.mpf(0)
        z = mp.exp(mp.pi * (y - lam)) * mp.exp(-I * mp.pi * x)
        return -mp.im(mp.polylog(2, z)) / mp.pi ** 2

    print(" lam      W(double)        W(mp, 50 digits)       bound            ratio(mp)")
    for jj, lam in enumerate(LAMS):
        lamm = mp.mpf(float(lam))
        Wm = sum(h(n, lamm) for n in nus) - sum(h(m_, lamm) for m_ in mus)
        if bound[jj] > 0:
            print(f"{lam:+.2f}  {WL[jj]: .6e}   {mp.nstr(Wm, 12):>22}   {bound[jj]: .6e}   {mp.nstr(Wm / bound[jj], 10)}")
