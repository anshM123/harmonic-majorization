"""Theorem A directly: D(a,t) vs a^2 int int e^{as - t tau} rho, complex (a,t), general B (interior eigenvalues,
eigenvalues 0 and 1). Product quadrature: cosine map in tau on each gap; composite Gauss-Legendre in s with many panels."""
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss
from hm_common import *

seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
rng = np.random.default_rng(seed)
NPHI = 160
NPAN, NGL = 400, 8


def grid_rho(c):
    pv = c["pvals"]
    xg, wg = leggauss(NGL)
    edges = np.linspace(c["lmin"], c["lmax"], NPAN + 1)
    s = ((edges[:-1, None] + edges[1:, None]) / 2 + (edges[1:, None] - edges[:-1, None]) / 2 * xg[None, :]).ravel()
    ws = ((edges[1:, None] - edges[:-1, None]) / 2 * wg[None, :]).ravel()
    ph, wph = leggauss(NPHI)
    ph = (ph + 1) * np.pi / 2
    wph = wph * np.pi / 2
    T, WT = [], []
    for k in range(len(pv) - 1):
        a, b = pv[k], pv[k + 1]
        h = b - a
        T.append(a + h * (1 - np.cos(ph)) / 2)
        WT.append(wph * h / 2 * np.sin(ph))
    T, WT = np.concatenate(T), np.concatenate(WT)
    R = np.array([[rho(c, si, ti) for si in s] for ti in T])
    return s, ws, T, WT, R


cases = [[0.0, 0.4, 1.0], [0.2, 0.2, 0.7, 0.9], [0.0, 0.0, 0.35, 0.35, 0.7, 1.0]]
for lv in cases:
    c = make_case(lv, rng)
    s, ws, T, WT, R = grid_rho(c)
    mass = WT @ R @ ws
    exact_mass = sum(w for _, _, w in offdiag_weights(c))
    print(f"levels={lv}: total mass {mass:.8f} vs sum_j<k ||Q_j g Q_k||^2 = {exact_mass:.8f}", flush=True)
    for a, t in [(1.0, 0.5), (0.6, -1.0), (-0.8, 1.5), (0.7 + 0.4j, 2.0 - 3.0j), (1.3, 6j), (-0.5j, 1.0 + 1.0j)]:
        lhs = Dfun(c, a, t)
        E = np.exp(a * s[None, :] - t * T[:, None])
        rhs = a * a * (WT @ (E * R) @ ws)
        print(f"   a={a}, t={t}: D={complex(lhs):.8f}  rhs={complex(rhs):.8f}  |diff|={abs(lhs-rhs):.2e}  rel={abs(lhs-rhs)/max(abs(lhs),1e-300):.1e}",
              flush=True)
