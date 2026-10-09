"""Check S6: Section 4 ("What is false") of 03-structure.md, reproduced from scratch.

CE1  "B = diag(1,0,0) with g tridiagonal (1/4, 1/100) has an entirely real spectrum": exact (sympy) for several
     readings of the description.
CE3  partial pinching E_F, B = diag(1,1/2,0), F = {{1,0},{1/2}} (= conditional expectation = commutant mixture):
     min_lam [W^L_{B+ig} - W^L_{B+iE_F g}] (60 digits), and the claim that 1-D projections of R decrease.
CE4  projection B = diag(1,1,0), Phi(g) = (g + U g U*)/2, U in {B}': same quantities (E o Phi != E here).
CE6  one-sided log-potential comparison: log|det(z-X)| - log|det(z-X_E)| changes sign.
SH   naive Schur-Horn converse: the 2x2 sweeps omega(sigma) have density bounded away from d, uniformly in sigma.
"""
import time
import numpy as np
import sympy as sp
import mpmath as mp
from hmaudit import h_closed_mp, setup

t0 = time.time()
nu = sp.symbols("nu")
Q = sp.Rational


def ce1(name, g):
    X = sp.diag(1, 0, 0) + sp.I * g
    p = sp.Poly(sp.expand((nu * sp.eye(3) - X).det()), nu)
    disc = p.discriminant()
    roots = [complex(sp.N(r, 20)) for r in sp.Poly(p, nu).nroots(n=20)]
    Eg = sp.Matrix([[g[0, 0], 0, 0], [0, g[1, 1], g[1, 2]], [0, g[2, 1], g[2, 2]]])
    eE = sorted([sp.nsimplify(e) for e in Eg.eigenvals(multiple=True)], key=lambda e: float(e))
    allreal = all(abs(r.imag) < 1e-18 for r in roots)
    real_coeffs = all(sp.im(cf) == 0 for cf in p.all_coeffs())
    s_im2 = sum(r.imag ** 2 for r in roots)
    s_e2 = sum(float(e) ** 2 for e in eE)
    print(f"   {name}: char poly {p.as_expr()}")
    print(f"      real coefficients: {real_coeffs}; discriminant {disc} ({'>0' if real_coeffs and disc > 0 else 'n/a'});"
          f" spec(B+ig) = {[complex(round(r.real, 10), round(r.imag, 10)) for r in roots]}")
    print(f"      entirely real: {allreal};  spec E(g) = {eE};  sum (Im nu)^2 = {s_im2:.3e} vs sum spec(Eg)^2 = {s_e2:.3e}")


print("== CE1: B = diag(1,0,0), 'g tridiagonal (1/4, 1/100)'")
ce1("reading A (our notes): zero diagonal, g12 = 1/4, g23 = 1/100",
    sp.Matrix([[0, Q(1, 4), 0], [Q(1, 4), 0, Q(1, 100)], [0, Q(1, 100), 0]]))
ce1("reading B: zero diagonal, g12 = 1/100, g23 = 1/4",
    sp.Matrix([[0, Q(1, 100), 0], [Q(1, 100), 0, Q(1, 4)], [0, Q(1, 4), 0]]))
ce1("reading C: Toeplitz tridiag(diagonal 1/4, off-diagonal 1/100)",
    sp.Matrix([[Q(1, 4), Q(1, 100), 0], [Q(1, 100), Q(1, 4), Q(1, 100)], [0, Q(1, 100), Q(1, 4)]]))
ce1("reading D: Toeplitz tridiag(diagonal 1/100, off-diagonal 1/4)",
    sp.Matrix([[Q(1, 100), Q(1, 4), 0], [Q(1, 4), Q(1, 100), Q(1, 4)], [0, Q(1, 4), Q(1, 100)]]))

# ------------------------------------------------------------------------------------------------------------
mp.mp.dps = 60


def mp_eigs(levels, g):
    M = len(levels)
    X = mp.matrix(M, M)
    for i in range(M):
        for j in range(M):
            X[i, j] = (mp.mpf(levels[i]) if i == j else 0) + 1j * mp.mpf(g[i][j])
    return mp.eig(X)[0]


def Wdiff(levels, g1, g2, lam, side):
    e1, e2 = mp_eigs(levels, g1), mp_eigs(levels, g2)
    if side == "R":
        e1 = [1 - mp.conj(e) for e in e1]
        e2 = [1 - mp.conj(e) for e in e2]
    return mp.fsum(h_closed_mp(e, lam) for e in e1) - mp.fsum(h_closed_mp(e, lam) for e in e2)


def scan(levels, g1, g2, lo=-3, hi=4, n=141):
    best = (mp.inf, None, None)
    e1, e2 = mp_eigs(levels, g1), mp_eigs(levels, g2)
    for side in "LR":
        a1 = e1 if side == "L" else [1 - mp.conj(e) for e in e1]
        a2 = e2 if side == "L" else [1 - mp.conj(e) for e in e2]
        for lam in np.linspace(lo, hi, n):
            v = mp.fsum(h_closed_mp(e, lam) for e in a1) - mp.fsum(h_closed_mp(e, lam) for e in a2)
            if v < best[0]:
                best = (v, side, lam)
    # refine by golden section on the best side
    v, side, lam0 = best
    a1 = e1 if side == "L" else [1 - mp.conj(e) for e in e1]
    a2 = e2 if side == "L" else [1 - mp.conj(e) for e in e2]
    f = lambda l: mp.fsum(h_closed_mp(e, l) for e in a1) - mp.fsum(h_closed_mp(e, l) for e in a2)
    a, b = mp.mpf(lam0) - mp.mpf("0.05"), mp.mpf(lam0) + mp.mpf("0.05")
    gr = (mp.sqrt(5) - 1) / 2
    for _ in range(80):
        c1, c2 = b - gr * (b - a), a + gr * (b - a)
        if f(c1) < f(c2):
            b = c2
        else:
            a = c1
    lm = (a + b) / 2
    return v, side, lam0, f(lm), lm, f(mp.mpf(round(lam0, 2)))


def projections_decrease(levels, g, gphi, xis, ws_n=200):
    """max over xi, w of N_{xi,Phi g}(w) - N_{xi,g}(w) (should be <= 0 if 1-D projections decrease)."""
    c1, c2 = setup(levels, np.array(g, complex)), setup(levels, np.array(gphi, complex))
    B = np.diag(levels)
    worst = -np.inf
    for xi in xis:
        l1, l1E = np.linalg.eigvalsh(c1["g"] - xi * B), np.linalg.eigvalsh(c1["gE"] - xi * B)
        l2, l2E = np.linalg.eigvalsh(c2["g"] - xi * B), np.linalg.eigvalsh(c2["gE"] - xi * B)
        lo, hi = min(l1.min(), l2.min()) - 1, max(l1.max(), l2.max()) + 1
        for w in np.linspace(lo, hi, ws_n):
            N1 = np.sum(np.maximum(l1 - w, 0)) - np.sum(np.maximum(l1E - w, 0))
            N2 = np.sum(np.maximum(l2 - w, 0)) - np.sum(np.maximum(l2E - w, 0))
            worst = max(worst, N2 - N1)
    # x-marginal: ||Q_j g Q_k||^2 per pair
    lab = c1["lab"]
    mx = -np.inf
    for j in np.unique(lab):
        for k in np.unique(lab):
            if j < k:
                Ij, Ik = np.where(lab == j)[0], np.where(lab == k)[0]
                mx = max(mx, np.sum(np.abs(c2["g"][np.ix_(Ij, Ik)]) ** 2) - np.sum(np.abs(c1["g"][np.ix_(Ij, Ik)]) ** 2))
    return worst, mx, np.abs(c2["gE"] - c1["gE"]).max()


print("\n== CE3: partial pinching, B = diag(1, 1/2, 0), F = {{1,0},{1/2}}")
lv = [1.0, 0.5, 0.0]
g = [[0.5, 1.5, 1.0], [1.5, -1.0, -0.5], [1.0, -0.5, 0.5]]
V = np.diag([1.0, -1.0, 1.0])
gF = ((np.array(g) + V @ np.array(g) @ V) / 2).tolist()
print(f"   E_F(g) = (g + VgV)/2 = {gF}  (V = diag(1,-1,1) commutes with B)")
v, side, lam0, vmin, lmin, v_at = scan(lv, g, gF)
print(f"   grid min {mp.nstr(v, 12)} at side {side}, lam {lam0:.3f};  refined min {mp.nstr(vmin, 12)} at lam "
      f"{mp.nstr(lmin, 8)};  value at lam = {round(lam0,2)}: {mp.nstr(v_at, 12)}")
w, mx, dE = projections_decrease(lv, g, gF, np.linspace(-6, 6, 121))
print(f"   1-D projections: max_xi,w [N_(xi,E_F g) - N_(xi,g)] = {w:.3e};  max increase of ||Q_j g Q_k||^2 = {mx:.3e};"
      f"  max|E(E_F g) - E(g)| = {dE:.1e}")

print("\n   adjacent merge F = {{1,1/2},{0}} (our notes): g = [[1/4,7/4,1/2],[7/4,1/2,-1],[1/2,-1,-3/4]]")
g2 = [[0.25, 1.75, 0.5], [1.75, 0.5, -1.0], [0.5, -1.0, -0.75]]
V2 = np.diag([1.0, 1.0, -1.0])
g2F = ((np.array(g2) + V2 @ np.array(g2) @ V2) / 2).tolist()
v, side, lam0, vmin, lmin, v_at = scan(lv, g2, g2F)
print(f"   grid min {mp.nstr(v, 12)} at side {side}, lam {lam0:.3f};  refined min {mp.nstr(vmin, 12)} at lam {mp.nstr(lmin, 8)}")

print("\n== CE4: projection B = diag(1,1,0), Phi(g) = (g + U g U*)/2, U = [[0,-1,0],[1,0,0],[0,0,1]]")
lv4 = [1.0, 1.0, 0.0]
g4 = np.array([[-1.0, 2.0, -1.0], [2.0, 3.0, 0.0], [-1.0, 0.0, 0.0]])
U = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1.0]])
print(f"   U commutes with B: {np.allclose(U @ np.diag(lv4), np.diag(lv4) @ U)}")
g4m = (g4 + U @ g4 @ U.T) / 2
print(f"   Phi(g) = {g4m.tolist()}")
v, side, lam0, vmin, lmin, v_at = scan(lv4, g4.tolist(), g4m.tolist())
print(f"   grid min {mp.nstr(v, 12)} at side {side}, lam {lam0:.3f};  refined min {mp.nstr(vmin, 12)} at lam "
      f"{mp.nstr(lmin, 8)};  value at lam = {round(lam0,2)}: {mp.nstr(v_at, 12)}")
w, mx, dE = projections_decrease(lv4, g4, g4m, np.linspace(-6, 6, 121))
print(f"   1-D projections: max_xi,w [N_(xi,Phi g) - N_(xi,g)] = {w:.3e};  max increase of ||Q_j g Q_k||^2 = {mx:.3e};"
      f"  max|E(Phi g) - E(g)| = {dE:.1e}  (E o Phi != E here)")

print("\n== CE6: log-potential comparison, B = diag(1,0), g = [[0,1/4],[1/4,0]]")
for z0 in [mp.mpf(-1), mp.mpc("0.5", "1")]:
    X = mp.matrix([[1, 0.25j], [0.25j, 0]])
    XE = mp.matrix([[1, 0], [0, 0]])
    d1 = mp.log(abs(mp.det(z0 * mp.eye(2) - X)))
    d2 = mp.log(abs(mp.det(z0 * mp.eye(2) - XE)))
    print(f"   z = {mp.nstr(z0, 4)}: log|det(z-X)| - log|det(z-X_E)| = {mp.nstr(d1 - d2, 10)}")

print("\n== SH: naive Schur-Horn converse; left-sweep density of [[1+ia, i s],[i s, i d]] away from d")
mp.mp.dps = 20
a, d, delta = 0.3, 0.0, 0.5
worst = 0
for s in np.logspace(-4, 4, 161):
    X = np.array([[1 + 1j * a, 1j * s], [1j * s, 1j * d]])
    ev = np.linalg.eigvals(X)
    ws = np.concatenate([np.linspace(d - 30, d - delta, 600), np.linspace(d + delta, d + 30, 600)])
    dens = np.zeros_like(ws)
    for e in ev:
        x, y = e.real, e.imag
        if 1e-14 < x < 1 - 1e-14:
            dens += np.sin(np.pi * x) / (2 * (2 * np.sinh(np.pi * (y - ws) / 2) ** 2 + 2 * np.sin(np.pi * x / 2) ** 2))
    worst = max(worst, dens.max())
print(f"   sigma in [1e-4, 1e4] (161 values): max density on |w - d| >= {delta}: {worst:.4f}  (bounded => no atoms at d+-1)")
print(f"time {time.time()-t0:.0f}s")
