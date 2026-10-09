"""
ch_cone.py -- cell-embedding linear programme for the chained clock model (math/05-chained-bell.md, Section 5).
[27]: github.com/anshM123/IQOQI-OQP-27, papers/math/main.tex.

Variables: pair counts C(u) = sum_x tau(Q_x Q_{x+u}), u in Z_L (L = 2Nd), u != 0 mod d, C(-u) = C(u).
Valid constraints for every projective strategy on a maximally entangled state:
  (a) C(u) >= 0;   (b) sum_{t in Z_2N} C(r + t d) = d for every residue r != 0 mod d;
  (c) cell inequalities, for every cell vector ell (sum d) and block size s = 1..floor(N/2):
        sum_u kt^{ell,s}(u) C(u) <= L^2 Phi_c(s/2N, (N-s)/2N),
      kt^{ell,s}(u) = (1/L) sum_{j<s<=j'<N} Kbar^ell(u - (j'-j) d),  Kbar^ell(w) = sum_x K^ell(x, x + w),
      K^ell(x,y) = Xi(L_{x+1}-L_y) - Xi(L_x-L_y) - Xi(L_{x+1}-L_{y+1}) + Xi(L_x-L_{y+1}),  Xi(u) = -(L^2/2pi^2) Cl2(2 pi u/L).
      (Block pairing <X_s, Y_s> of X_s = sum_{j<s} rot_{-jd} Q and Y_s = sum_{s<=j<N} rot_{-jd} Q, cell-embedded:
       = L^2 int tau(X_s H[Y_s]) dm <= L^2 Phi_c by the continuum theorem [27, Thm 4.2], equality for windows.)
Objective: F_N = (1/2) sum_{u != 0 mod d} K_N(u) C(u).  If the LP maximum equals F_nat (the window value), the
cell inequalities prove F_N <= F_nat for all projective strategies on maximally entangled states (for this N, d).
Usage: python ch_cone.py N d [ncells] [seed]   (N = 2 reproduces CGLMP; for d = 2..6 it uses the Regime I
certificate cells of [27], copied unchanged into oqp27_cells/ from IQOQI-OQP-27/cone-certificates/certs/)
"""
import os
import sys
import json

import numpy as np
from scipy.optimize import linprog
from scipy.special import spence

from ch_clock import K_N, F_natural


def cl2(theta):
    th = np.asarray(theta, float)
    out = np.imag(spence(1 - np.exp(1j * th)))
    return np.where(np.abs(np.sin(th / 2)) < 1e-300, 0.0, out)


def Phi_c(a, b):
    return (cl2(2 * np.pi * a) + cl2(2 * np.pi * b) + cl2(2 * np.pi * (1 - a - b))) / (2 * np.pi ** 2)


def cell_kernel(ell, N):
    ell = np.asarray(ell, float)
    d = len(ell)
    L = 2 * N * d
    Lc = np.concatenate([[0.0], np.cumsum(np.tile(ell, 2 * N))])
    a, b = Lc[:-1], Lc[1:]

    def Xi(u):
        return -(L ** 2 / (2 * np.pi ** 2)) * cl2(2 * np.pi * np.asarray(u, float) / L)
    K = (Xi(b[:, None] - a[None, :]) - Xi(a[:, None] - a[None, :]) - Xi(b[:, None] - b[None, :])
         + Xi(a[:, None] - b[None, :]))
    np.fill_diagonal(K, 0.0)
    return K


def kt_vec(ell, N, s):
    d = len(ell)
    L = 2 * N * d
    K = cell_kernel(ell, N)
    Kbar = np.array([sum(K[x, (x + w) % L] for x in range(L)) for w in range(L)])
    kt = np.zeros(L)
    for j in range(s):
        for jp in range(s, N):
            kt += np.roll(Kbar, (jp - j) * d)        # kt[u] += Kbar[u - (jp - j) d]
    return kt / L


def window_C(N, d):
    L = 2 * N * d
    C = np.zeros(L)
    for u in range(L):
        du = min(u, L - u)
        C[u] = max(d - du, 0)
    return C


def setup(N, d):
    L = 2 * N * d
    U = [u for u in range(1, L // 2 + 1) if u % d]          # representatives (C(-u) = C(u))
    return L, U


def fold(vec, L, U):
    """coefficient vector on the representatives: sum_{u in Z_L} vec[u] C(u) = sum_{u in U} out[u] C(u)."""
    out = np.zeros(len(U))
    for i, u in enumerate(U):
        out[i] = vec[u] + (vec[(L - u) % L] if (L - u) % L != u else 0.0)
    return out


def run(N, d, cells, verbose=True):
    L, U = setup(N, d)
    obj = np.zeros(L)
    for u in range(L):
        if u % d:
            obj[u] = 0.5 * K_N(N, d, u)
    c = -fold(obj, L, U)
    # equality: sum_t C(r + t d) = d
    Aeq, beq = [], []
    for r in range(1, d):
        row = np.zeros(L)
        for t in range(2 * N):
            row[(r + t * d) % L] += 1.0
        Aeq.append(fold(row, L, U))
        beq.append(d)
    Aub, bub = [], []
    Cw = window_C(N, d)
    worst_eq = 0.0
    for ell in cells:
        for s in range(1, N // 2 + 1):
            kt = kt_vec(ell, N, s)
            bnd = L ** 2 * Phi_c(s / (2 * N), (N - s) / (2 * N))
            worst_eq = max(worst_eq, abs(kt @ Cw - bnd) / abs(bnd))
            Aub.append(fold(kt, L, U))
            bub.append(bnd - kt[0] * d)            # C(0) = d; C(u) = 0 for u = 0 mod d, u != 0
    res = linprog(c, A_ub=np.array(Aub), b_ub=np.array(bub), A_eq=np.array(Aeq), b_eq=np.array(beq),
                  bounds=[(0, None)] * len(U), method="highs")
    Fn = F_natural(N, d)
    if verbose:
        print(f"N={N} d={d}: cells={len(cells)}  window tightness max rel err={worst_eq:.2e}  LP status={res.status}  "
              f"LP max F={-res.fun:.10f}  F_nat={Fn:.10f}  gap={-res.fun - Fn:+.3e}")
    return -res.fun, Fn, res


def random_cells(d, n, rng, q=16):
    out = []
    for _ in range(n):
        alpha = rng.choice([0.3, 1.0, 3.0])
        x = rng.dirichlet(np.ones(d) * alpha)
        nint = np.maximum(1, np.round(x * q * d)).astype(int)
        # fix the sum to q d
        while nint.sum() != q * d:
            i = rng.integers(0, d)
            if nint.sum() > q * d and nint[i] > 1:
                nint[i] -= 1
            elif nint.sum() < q * d:
                nint[i] += 1
        out.append(nint / q)
    return out


def main():
    N = int(sys.argv[1])
    d = int(sys.argv[2])
    ncells = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    rng = np.random.default_rng(seed)
    cells = [np.ones(d)] + random_cells(d, ncells, rng)
    if N == 2:
        here = os.path.dirname(os.path.abspath(__file__))
        cert = os.path.join(here, "oqp27_cells", f"cert_d{d}.json")
        if os.path.exists(cert):
            with open(cert) as f:
                cc = json.load(f)
            cells = [np.asarray(x, float) / cc["q"] for x in cc["cells"]]
            print("using the Regime I certificate cells")
    run(N, d, cells)


if __name__ == "__main__":
    main()
