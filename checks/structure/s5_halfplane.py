"""Check S5: Theorem H (half-plane, accretive matrices) of 03-structure.md.

W^H(lam) = int (w - lam)_+ d(pi(X) - pi(X_E)), pi_nu = Cauchy(y, x) (harmonic measure of {Re z > 0} on iR seen from
nu = x + iy; delta_y if x = 0).  Closed form per point after the divergent (x/pi) log L terms cancel (derived here):
   G(nu, lam) = -(x/2pi) log((lam-y)^2 + x^2) + ((lam-y)/pi) arctan((lam-y)/x) + (y-lam)/2,   G = (y-lam)_+ if x = 0.
Part 1: W^H >= 0 on random B >= 0 (levels in [0, 20], zero and repeated levels), 40-digit arithmetic.
Part 2: W^H = V^H := int int (1/pi) beta/(beta^2+(s-lam)^2) rho_B ds dbeta (quadrature), b_min > 0 and b_min = 0.
Part 3: int w^2 d(pi(X) - pi(X_E)) = ||K||_F^2 (absolutely convergent; consequence of W^H = V^H), and decay O(w^-4).
Part 4: int phi d(pi(X)-pi(X_E)) >= 0 for convex phi of growth |w|^p, p < 3 (phi = |w-a|^2.5, (w-a)_+^2).
Usage: python s5_halfplane.py [seed] [ncases]
"""
import sys
import time
import numpy as np
import mpmath as mp
from hmaudit import rand_herm, setup, adaptive_nodes_1d, _breakpoints_s, rho_vals

seed = int(sys.argv[1]) if len(sys.argv) > 1 else 21
NC = int(sys.argv[2]) if len(sys.argv) > 2 else 400
rng = np.random.default_rng(seed)
t0 = time.time()
mp.mp.dps = 40


def G(nu, lam):
    x, y = mp.re(nu), mp.im(nu)
    lam = mp.mpf(lam)
    if x == 0:
        return max(y - lam, 0)
    return -(x / (2 * mp.pi)) * mp.log((lam - y) ** 2 + x ** 2) + ((lam - y) / mp.pi) * mp.atan((lam - y) / x) \
        + (y - lam) / 2


def mp_spec(levels, g):
    """eigenvalues at 40 digits (mpmath), X and X_E; X_E blockwise (exact atoms on the lines)."""
    M = len(levels)
    X = mp.matrix(M, M)
    for i in range(M):
        for j in range(M):
            X[i, j] = (mp.mpf(levels[i]) if i == j else 0) + 1j * mp.mpc(g[i, j].real, g[i, j].imag)
    ev, _ = mp.eig(X)
    lab = np.unique(levels, return_inverse=True)[1]
    evE = []
    for L in np.unique(lab):
        idx = np.where(lab == L)[0]
        H = mp.matrix(len(idx), len(idx))
        for a, i in enumerate(idx):
            for b, j in enumerate(idx):
                H[a, b] = mp.mpc(g[i, j].real, g[i, j].imag)
        e = mp.eighe(H)[0] if hasattr(mp, "eighe") else mp.eigh(H)[0]
        evE += [mp.mpf(levels[idx[0]]) + 1j * mp.re(v) for v in e]
    # snap tiny negative real parts produced by rounding (spec X lies in Re z >= 0)
    ev = [mp.mpc(max(mp.re(z), 0), mp.im(z)) for z in ev]
    return ev, evE


def WH(ev, evE, lam):
    return mp.fsum(G(z, lam) for z in ev) - mp.fsum(G(z, lam) for z in evE)


print("== Part 1: W^H >= 0 on random accretive pairs (40 digits)")
lams = np.concatenate([np.linspace(-12, 12, 97), [0.0]])
minr, info, nneg = mp.inf, None, 0
for case in range(NC):
    m = int(rng.integers(2, 6))
    lv = np.round(rng.uniform(0, rng.choice([1.0, 5.0, 20.0]), m), 6)
    if rng.uniform() < 0.4:
        lv[0] = 0.0
    mult = rng.integers(1, 3, size=m)
    levels = np.repeat(lv, mult).astype(float)
    M = len(levels)
    g = rand_herm(M, rng, rng.uniform(0.1, 4))
    c = setup(levels, g)
    if c["K2"] < 1e-10:
        continue
    ev, evE = mp_spec(levels, g)
    for lam in lams:
        w = WH(ev, evE, lam)
        r = w / c["K2"]
        if r < minr:
            minr, info = r, (case, list(levels), float(lam))
        if w < -mp.mpf(10) ** (-25):
            nneg += 1
print(f"   {NC} cases x {len(lams)} lambdas: negatives (< -1e-25) = {nneg};  min W^H/||K||^2 = {mp.nstr(minr, 5)} at {info}")
print(f"   [elapsed {time.time()-t0:.0f}s]")

print("\n== Part 2: W^H = V^H by quadrature")


def Pk(beta, u):
    return beta / (np.pi * (beta ** 2 + u ** 2))


def VH(c, lam):
    def inner(beta):
        pts = _breakpoints_s(c, beta, 600)
        if c["lmin"] < lam < c["lmax"]:
            pts = sorted(set(pts + [lam]))
        tot = 0.0
        for a, b in zip(pts[:-1], pts[1:]):
            f = lambda s: rho_vals(c, s, beta) * Pk(beta, s - lam)
            ss, sw = adaptive_nodes_1d(f, a, b, tol=1e-11, maxit=60)
            tot += np.sum(sw * f(ss))
        return tot
    tot = 0.0
    u = c["uniq"]
    for k in range(len(u) - 1):
        fv = lambda bs: np.array([inner(b) for b in bs])
        bx, bw = adaptive_nodes_1d(fv, u[k], u[k + 1], tol=1e-9, nmin=4, maxit=30)
        tot += np.sum(bw * fv(bx))
    return tot


for levels in ([0.3, 1.0, 1.0, 2.5], [0.0, 0.0, 0.7, 1.6]):
    g = rand_herm(len(levels), rng, 0.9)
    c = setup(levels, g)
    ev, evE = mp_spec(np.array(levels), g)
    for lam in [-1.0, 0.0, 0.5, 2.0]:
        w = WH(ev, evE, lam)
        v = VH(c, lam)
        print(f"   levels {levels}, lam={lam:+.2f}: W^H = {mp.nstr(w, 14)}   V^H = {v:.14f}   |diff| = {abs(float(w)-v):.1e}")
print(f"   [elapsed {time.time()-t0:.0f}s]")

print("\n== Part 3: second moment int w^2 d(pi(X)-pi(X_E)) = ||K||^2 and O(w^-4) decay")


def density_diff(ev, evE, w):
    s = mp.mpf(0)
    for z in ev:
        x, y = mp.re(z), mp.im(z)
        if x > 0:
            s += x / (mp.pi * (x ** 2 + (w - y) ** 2))
    for z in evE:
        x, y = mp.re(z), mp.im(z)
        if x > 0:
            s -= x / (mp.pi * (x ** 2 + (w - y) ** 2))
    return s


def atoms_moment(ev, evE, phi):
    return mp.fsum(phi(mp.im(z)) for z in ev if mp.re(z) == 0) - mp.fsum(phi(mp.im(z)) for z in evE if mp.re(z) == 0)


mp.mp.dps = 40
for case in range(4):
    levels = [[0.2, 1.1, 1.1], [0.0, 0.5, 2.0, 2.0], [0.0, 0.0, 3.0], [0.4, 0.9, 1.7, 4.0]][case]
    g = rand_herm(len(levels), rng, [0.7, 1.5, 1.0, 2.2][case])
    c = setup(levels, g)
    ev, evE = mp_spec(np.array(levels), g)
    ys = sorted([float(mp.im(z)) for z in ev + evE])
    brk = [-mp.inf] + [mp.mpf(v) for v in ys] + [mp.inf]
    m2 = mp.quad(lambda w: w ** 2 * density_diff(ev, evE, w), brk) + atoms_moment(ev, evE, lambda y: y ** 2)
    tail = [mp.nstr(density_diff(ev, evE, mp.mpf(W)) * W ** 4, 6) for W in (1e2, 1e3, 1e4)]
    phis = [("|w-0.3|^2.5", lambda w: abs(w - mp.mpf("0.3")) ** mp.mpf(2.5)),
            ("(w+0.5)_+^2", lambda w: max(w + mp.mpf("0.5"), 0) ** 2),
            ("|w|^2.9", lambda w: abs(w) ** mp.mpf(2.9))]
    vals = []
    for name, phi in phis:
        brk2 = sorted(set(brk[1:-1] + [mp.mpf("0.3"), mp.mpf("-0.5"), mp.mpf(0)]))
        val = mp.quad(lambda w: phi(w) * density_diff(ev, evE, w), [-mp.inf] + brk2 + [mp.inf]) + atoms_moment(ev, evE, phi)
        vals.append(f"{name}: {mp.nstr(val, 10)}")
    print(f"   levels {levels}: int w^2 d(diff) = {mp.nstr(m2, 16)}  ||K||^2 = {c['K2']:.15f};  w^4*density at 1e2,1e3,1e4: {tail}")
    print(f"      convex phi: " + ";  ".join(vals))
print(f"time {time.time()-t0:.0f}s")
