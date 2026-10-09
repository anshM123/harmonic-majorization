"""Check S7c: print and re-confirm (60 digits, several gamma pairs, both edges) the 3-level violation of
dephasing-line monotonicity found by s7_dephasing.py (saved in logs/three_level_candidate.json)."""
import numpy as np
import mpmath as mp
from hmaudit import h_closed_mp, labels_of, pinch

import json
_c = json.load(open("logs/three_level_candidate.json"))
levels = _c["levels"]
g = np.array(_c["g_re"]) + 1j * np.array(_c["g_im"])
g1, g2 = _c["gamma_1"], _c["gamma_2"]
levels = np.asarray(levels, float)
g = np.asarray(g, complex)
lab, _ = labels_of(levels)
gE = pinch(g, lab)
np.set_printoptions(precision=17, linewidth=160)
print("levels (B = diag):", levels.tolist())
print("g (Hermitian) =")
print(repr(g))
print(f"gamma_1 = {g1!r}, gamma_2 = {g2!r}, ||g - E g||_F^2 = {np.sum(np.abs(g-gE)**2):.6f}, ||g||_op = {np.linalg.norm(g,2):.4f}")
mp.mp.dps = 60


def eig(gg):
    X = mp.matrix(3, 3)
    for i in range(3):
        for j in range(3):
            X[i, j] = (mp.mpf(levels[i]) if i == j else 0) + 1j * mp.mpc(gg[i, j].real, gg[i, j].imag)
    return mp.eig(X)[0]


def W(ev, lam, side):
    z = ev if side == "L" else [1 - mp.conj(e) for e in ev]
    return mp.fsum(h_closed_mp(e, lam) for e in z)


for (ga, gb) in [(g1, g2), (0.6, 0.9), (0.7, 0.95), (0.8, 1.0), (0.5, 1.0), (0.9, 1.0)]:
    ea, eb = eig(gE + ga * (g - gE)), eig(gE + gb * (g - gE))
    best = (mp.inf, None)
    for side in "LR":
        for lam in np.linspace(-8, 8, 321):
            v = W(eb, lam, side) - W(ea, lam, side)
            if v < best[0]:
                best = (v, (side, round(float(lam), 3)))
    print(f"gamma {ga:.6f} -> {gb:.6f}: min_lambda [W(gamma_b) - W(gamma_a)] = {mp.nstr(best[0], 12)} at {best[1]}")
