"""Theorem A, tested through its equivalent line form (Radon slices), general B with interior eigenvalues:
   int rho(w + xi tau, tau) dtau = U_xi(w) = sum_j (lam_j(xi) - w)_+ - (lam0_j(xi) - w)_+ ,
lam_j = eig(g - xi P), lam0_j = eig(g_d - xi P).  Theorem A for real (a, t=a xi) is the Laplace transform in w of this."""
import sys
import numpy as np
from scipy.integrate import quad
from hm_common import *

seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
rng = np.random.default_rng(seed)


def line_integral(c, w, xi):
    pv = c["pvals"]
    tot, err = 0.0, 0.0
    for k in range(len(pv) - 1):
        a, b = pv[k], pv[k + 1]
        h = b - a
        f = lambda ph: rho(c, w + xi * (a + h * (1 - np.cos(ph)) / 2), a + h * (1 - np.cos(ph)) / 2) * h / 2 * np.sin(ph)
        v, e = quad(f, 0, np.pi, limit=800, epsabs=1e-11, epsrel=1e-10)
        tot += v
        err += e
    return tot, err


def U(c, w, xi):
    lam = np.linalg.eigvalsh(c["g"] - xi * c["P"])
    lam0 = np.linalg.eigvalsh(c["gd"] - xi * c["P"])
    return np.sum(np.maximum(lam - w, 0)) - np.sum(np.maximum(lam0 - w, 0))


cases = [[0.0, 0.4, 1.0], [0.2, 0.2, 0.7, 0.9], [0.0, 0.0, 0.35, 0.35, 0.7, 1.0], [0.1, 0.5, 0.5, 0.95, 0.3]]
worst = 0
for lv in cases:
    c = make_case(lv, rng)
    for _ in range(6):
        xi = rng.normal() * 2
        lam = np.linalg.eigvalsh(c["g"] - xi * c["P"])
        w = rng.uniform(lam[0] - 0.3, lam[-1] + 0.3)
        L, e = line_integral(c, w, xi)
        u = U(c, w, xi)
        worst = max(worst, abs(L - u))
        print(f"levels={lv} xi={xi:+.3f} w={w:+.3f}: line int {L:.10f} (quad err {e:.1e})  U {u:.10f}  diff {L-u:+.2e}",
              flush=True)
print("max |diff| =", worst)
