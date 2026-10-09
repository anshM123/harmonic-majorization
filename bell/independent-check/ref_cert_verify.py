"""
ref_cert_verify.py -- independent rigorous re-verification of the chained-LP certificates (Theorem 5 of
math/05-chained-bell.md).

Mathematical content verified (derived independently for this check):
  pair counts C(u) = sum_x tau(Q_x Q_{x+u}) on Z_L, L = 2Nd; valid constraints C >= 0, C(-u) = C(u), C(0) = d,
  C(u) = 0 for u = 0 mod d (u != 0), sum_{t in Z_2N} C(r + t d) = d (r != 0 mod d);
  F_N = (1/2) sum_{u != 0 mod d} K_N(u) C(u), K_N(u) = N sin(pi(N-1)u/(Nd))/sin(pi u/d);
  averaged cell inequality for cell ell (sum d, entries > 0) and block s (1 <= s <= N/2):
      sum_u kt(u) C(u) <= b := L^2 Phi_c(s/2N, (N-s)/2N),
      kt(u) = (1/d) sum_{j=0}^{d-1} sum_{i<s<=i'<N} K^ell(j - i d, u + j - i' d)      [coefficient of tau(Q_0 Q_u)]
      K^ell(x,y) = Xi(L_{x+1}-L_y) - Xi(L_x-L_y) - Xi(L_{x+1}-L_{y+1}) + Xi(L_x-L_{y+1}), K(x,x) = 0,
      Xi(v) = -(L^2/(2 pi^2)) Cl2(2 pi v/L).
Certificate: y >= 0 (cells), lambda (residues) with folded slacks z_u = sum_k y_k a_k(u) + sum_r lambda_r e_r(u) - K_N(u)
  (u = 1..Nd-1, u != 0 mod d, a_k(u) = kt_k(u) + kt_k(-u)), z_u = 0 for u <= d-1 and z_u > 0 for u >= d+1.
Rigorous part: python-flint (arb/acb ball arithmetic); Cl2 via acb.polylog(2, e^{i theta}); the exact solution of the
  equalities is enclosed by x = x0 + G^T w, (G G^T) w = o - G x0, solved with arb_mat.solve.
Modes:
  own    : an independent cell pool + own LP (scipy HiGHS), then rigorous verification
  theirs : the authors' pool/LP (imports ../scripts/ch_cone.py and ../scripts/ch_clock.py, read-only), then the
           independent rigorous verification of that certificate
Usage: python ref_cert_verify.py MODE N d [npool] [seed] [prec]
"""
import os
import sys
import time

import mpmath
import numpy as np
from flint import acb, arb, arb_mat, ctx
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "..", "scripts")

Q_DEN = 8          # cell boundaries are multiples of 1/8


# ----------------------------------------------------------------------------------------------- float kernel (own)
def cl2_grid_float(G):
    mpmath.mp.prec = 80
    out = np.zeros(G)
    for j in range(1, G):
        out[j] = float(mpmath.clsin(2, 2 * mpmath.pi * j / G))
    return out


def kt_float(nvec, N, s, clg):
    """nvec: integer cell lengths in units of 1/Q_DEN (sum = Q_DEN*d). returns kt[u], u in Z_L."""
    d = len(nvec)
    L = 2 * N * d
    G = Q_DEN * L
    ext = np.tile(np.asarray(nvec, int), 2 * N)
    Lc = np.concatenate([[0], np.cumsum(ext)])              # boundaries in units of 1/Q_DEN, Lc[L] = G

    def Xi(j):
        return -(L ** 2 / (2 * np.pi ** 2)) * clg[int(j) % G]

    def K(x, y):
        x %= L
        y %= L
        if x == y:
            return 0.0
        return Xi(Lc[x + 1] - Lc[y]) - Xi(Lc[x] - Lc[y]) - Xi(Lc[x + 1] - Lc[y + 1]) + Xi(Lc[x] - Lc[y + 1])
    kt = np.zeros(L)
    for u in range(L):
        acc = 0.0
        for j in range(d):
            for i in range(s):
                for ip in range(s, N):
                    acc += K(j - i * d, u + j - ip * d)
        kt[u] = acc / d
    return kt


def phi_c_float(a, b):
    c = lambda t: float(mpmath.clsin(2, 2 * mpmath.pi * t))
    return (c(a) + c(b) + c(1 - a - b)) / (2 * np.pi ** 2)


def KN_float(N, d, u):
    return N * np.sin(np.pi * (N - 1) * u / (N * d)) / np.sin(np.pi * u / d)


def reps(N, d):
    L = 2 * N * d
    return [u for u in range(1, N * d) if u % d]           # folded representatives (u = Nd is 0 mod d)


def residue_rows(N, d):
    U = reps(N, d)
    rows = []
    for r in range(1, d):
        rows.append([(1 if u % d == r else 0) + (1 if (-u) % d == r else 0) for u in U])
    return np.array(rows, float)


# ----------------------------------------------------------------------------------------------- own LP
def own_pool(d, npool, rng):
    cells = [np.full(d, Q_DEN, int)]
    for _ in range(npool):
        alpha = rng.choice([0.25, 0.7, 1.5, 4.0])
        x = rng.dirichlet(np.full(d, alpha))
        n = np.maximum(1, np.round(x * Q_DEN * d)).astype(int)
        while n.sum() != Q_DEN * d:
            i = rng.integers(0, d)
            if n.sum() > Q_DEN * d and n[i] > 1:
                n[i] -= 1
            elif n.sum() < Q_DEN * d:
                n[i] += 1
        cells.append(n)
    return cells


def solve_lp(N, d, rows_a, obj):
    """maximise t: equalities on u<=d-1, z_u >= t(1+|o_u|) on u>=d+1, y>=0, t<=1."""
    U = reps(N, d)
    E = residue_rows(N, d)
    nr = len(rows_a)
    nv = nr + (d - 1) + 1
    must = [i for i, u in enumerate(U) if u < d]
    off = [i for i, u in enumerate(U) if u > d]
    c = np.zeros(nv)
    c[-1] = -1.0
    Aeq = np.zeros((len(must), nv))
    beq = np.zeros(len(must))
    for a_, i in enumerate(must):
        Aeq[a_, :nr] = [r[i] for r in rows_a]
        Aeq[a_, nr:nr + d - 1] = E[:, i]
        beq[a_] = obj[i]
    Aub = np.zeros((len(off), nv))
    bub = np.zeros(len(off))
    for a_, i in enumerate(off):
        Aub[a_, :nr] = [-r[i] for r in rows_a]
        Aub[a_, nr:nr + d - 1] = -E[:, i]
        Aub[a_, -1] = 1 + abs(obj[i])
        bub[a_] = -obj[i]
    bounds = [(0, None)] * nr + [(None, None)] * (d - 1) + [(None, 1.0)]
    res = linprog(c, A_ub=Aub if len(off) else None, b_ub=bub if len(off) else None, A_eq=Aeq, b_eq=beq,
                  bounds=bounds, method="highs")
    if res.status != 0:
        return None
    return res.x[:nr], res.x[nr:nr + d - 1], -res.fun


# ----------------------------------------------------------------------------------------------- rigorous (arb)
def cl2_grid_arb(G):
    out = [arb(0)] * G
    for j in range(1, G):
        if 2 * j == G:
            out[j] = arb(0)
            continue
        th = arb(2 * j) / G * arb.pi()
        z = acb(0, th).exp()
        out[j] = z.polylog(2).imag
    return out


def kt_arb(nvec, N, s, clg):
    d = len(nvec)
    L = 2 * N * d
    G = Q_DEN * L
    ext = np.tile(np.asarray(nvec, int), 2 * N)
    Lc = np.concatenate([[0], np.cumsum(ext)])
    pref = -(arb(L) ** 2) / (2 * arb.pi() ** 2)
    Kc = {}

    def K(x, y):
        x %= L
        y %= L
        if x == y:
            return arb(0)
        key = (x, y)
        if key not in Kc:
            Kc[key] = pref * (clg[int(Lc[x + 1] - Lc[y]) % G] - clg[int(Lc[x] - Lc[y]) % G]
                              - clg[int(Lc[x + 1] - Lc[y + 1]) % G] + clg[int(Lc[x] - Lc[y + 1]) % G])
        return Kc[key]
    kt = []
    for u in range(L):
        acc = arb(0)
        for j in range(d):
            for i in range(s):
                for ip in range(s, N):
                    acc += K(j - i * d, u + j - ip * d)
        kt.append(acc / d)
    return kt


def phi_c_arb(N, s):
    # Phi_c(s/2N, (N-s)/2N) = [Cl2(pi s/N) + Cl2(pi (N-s)/N) + Cl2(pi)] / (2 pi^2), Cl2(pi) = 0
    def cl(num, den):          # Cl2(2 pi num/den)
        if (2 * num) % den == 0:
            return arb(0)
        z = acb(0, arb(2 * num) / den * arb.pi()).exp()
        return z.polylog(2).imag
    return (cl(s, 2 * N) + cl(N - s, 2 * N)) / (2 * arb.pi() ** 2)


def KN_arb(N, d, u):
    return N * (arb(u * (N - 1)) / (N * d)).sin_pi() / (arb(u) / d).sin_pi()


def verify(N, d, cert_cells, y, lam, prec=128, verbose=True):
    """cert_cells: list of (nvec, s) with y_k > 0 (floats); lam: floats (d-1). Returns dict."""
    ctx.prec = prec
    L = 2 * N * d
    U = reps(N, d)
    nU = len(U)
    clg = cl2_grid_arb(Q_DEN * L)
    A = []          # folded a_k(u) balls
    Bk = []         # b_k - kt_k(0) d
    for (nvec, s) in cert_cells:
        kt = kt_arb(nvec, N, s, clg)
        A.append([kt[u] + kt[(L - u) % L] for u in U])
        Bk.append(arb(L) ** 2 * phi_c_arb(N, s) - kt[0] * d)
    E = residue_rows(N, d)
    o = [KN_arb(N, d, u) for u in U]
    must = [i for i, u in enumerate(U) if u < d]
    off = [i for i, u in enumerate(U) if u > d]
    nS = len(cert_cells)
    nx = nS + d - 1
    x0 = [arb(float(v)) for v in list(y) + list(lam)]

    def col(b_, i):
        return A[b_][i] if b_ < nS else arb(int(E[b_ - nS][i]))
    Gm = arb_mat(len(must), nx, [col(b_, i) for i in must for b_ in range(nx)])
    rvec = arb_mat(len(must), 1, [o[i] - sum((col(b_, i) * x0[b_] for b_ in range(nx)), arb(0)) for i in must])
    GG = Gm * Gm.transpose()
    w = GG.solve(rvec)
    corr = Gm.transpose() * w
    x = [x0[b_] + corr[b_, 0] for b_ in range(nx)]
    ymin = min((x[b_].lower() for b_ in range(nS)), default=arb(np.inf))
    zmins = []
    for i in off:
        zi = sum((x[b_] * col(b_, i) for b_ in range(nx)), arb(0)) - o[i]
        zmins.append(zi.lower())
    zmin = min(zmins) if zmins else arb(np.inf)
    # equalities re-checked on the corrected box (must contain 0)
    eqres = [sum((x[b_] * col(b_, i) for b_ in range(nx)), arb(0)) - o[i] for i in must]
    eq_ok = all(e.contains(0) for e in eqres)
    UB = sum((x[b_] * Bk[b_] for b_ in range(nS)), arb(0)) + sum((x[nS + r] * d for r in range(d - 1)), arb(0))
    Fnat = sum((arb(d - m) * KN_arb(N, d, m) for m in range(1, d)), arb(0))
    ok = (float(ymin) > 0 if nS else True) and (float(zmin) > 0 if zmins else True) and eq_ok
    radmax = max(float(v.rad()) for v in x) if x else 0.0
    out = dict(ok=ok, ymin=float(ymin) if nS else None, zmin=float(zmin) if zmins else None, eq_ok=eq_ok,
               UB=UB, Fnat=Fnat, UB_contains_Fnat=bool((UB - Fnat).contains(0)), radmax=radmax)
    if verbose:
        print(f"   [arb {prec} bits] cells={nS} min y lower = {out['ymin']}, min off-window z lower = {out['zmin']}, "
              f"equalities enclose 0: {eq_ok}, max radius of (y,lambda) = {radmax:.1e}")
        print(f"   bound = {UB.str(15)}, F_nat = {Fnat.str(15)}, (bound - F_nat) contains 0: {out['UB_contains_Fnat']}")
        print(("RIGOROUSLY VERIFIED" if ok else "NOT VERIFIED") + f"  N={N} d={d}")
    return out


# ----------------------------------------------------------------------------------------------- drivers
def mode_own(N, d, npool, seed, prec):
    t0 = time.time()
    rng = np.random.default_rng(seed)
    L = 2 * N * d
    U = reps(N, d)
    clg = cl2_grid_float(Q_DEN * L)
    obj = np.array([KN_float(N, d, u) for u in U])
    pool = own_pool(d, npool, rng) if d > 2 else []
    rows, tags = [], []
    for nvec in pool:
        for s in range(1, N // 2 + 1):
            kt = kt_float(nvec, N, s, clg)
            rows.append(np.array([kt[u] + kt[(L - u) % L] for u in U]))
            tags.append((nvec, s))
    sol = solve_lp(N, d, rows, obj)
    if sol is None:
        print(f"OWN N={N} d={d}: LP infeasible/failed")
        return False
    y, lam, t = sol
    S = [i for i in range(len(rows)) if y[i] > 1e-9 * max(1.0, y.max() if len(y) else 1.0)]
    print(f"OWN N={N} d={d}: own LP min-slack t* = {t:.3e}, rows used {len(S)} of {len(rows)} ({time.time() - t0:.1f}s)",
          flush=True)
    # sanity: window tightness of the independent kernel (float), sum_u kt(u) C_win(u) = L^2 Phi_c
    Cw = np.array([max(d - min(u, L - u), 0) for u in range(L)], float)
    wt = 0.0
    for i in S:
        nvec, s = tags[i]
        kt = kt_float(nvec, N, s, clg)
        b = L ** 2 * phi_c_float(s / (2 * N), (N - s) / (2 * N))
        wt = max(wt, abs(kt @ Cw - b) / abs(b))
    print(f"   window tightness of used cell rows (float): max rel err {wt:.1e}")
    if t <= 1e-10:
        print("   no strictly complementary certificate in this pool")
        return False
    out = verify(N, d, [tags[i] for i in S], [y[i] for i in S], lam, prec)
    print(f"   total {time.time() - t0:.1f}s", flush=True)
    return out['ok']


def mode_theirs(N, d, npool, seed, prec):
    """re-create the authors' LP2 certificate exactly as ch_cert.certify does (same pool, same seed), then verify."""
    t0 = time.time()
    sys.path.insert(0, os.path.abspath(SCRIPTS))
    from ch_clock import K_N as K_N_t
    from ch_cone import kt_vec, random_cells, setup, fold
    q = 8
    rng = np.random.default_rng(seed)
    cells_int = [np.ones(d, int) * q] + [np.round(c * q).astype(int) for c in random_cells(d, npool, rng, q=q)]
    L, Ut = setup(N, d)
    assert Ut == reps(N, d)
    objt = np.zeros(L)
    for u in range(L):
        if u % d:
            objt[u] = 0.5 * K_N_t(N, d, u)
    objf = fold(objt, L, Ut)
    rows, tags = [], []
    for nv in cells_int:
        for s in range(1, N // 2 + 1):
            rows.append(fold(kt_vec(nv / q, N, s), L, Ut))
            tags.append((nv, s))
    sol = solve_lp(N, d, rows, objf)
    if sol is None:
        print(f"THEIRS N={N} d={d}: LP failed")
        return False
    y, lam, t = sol
    S = [i for i in range(len(rows)) if y[i] > 1e-9 * max(1.0, y.max())]
    print(f"THEIRS N={N} d={d} (pool {npool}, seed {seed}): LP t* = {t:.3e}, rows used {len(S)} "
          f"({time.time() - t0:.1f}s)", flush=True)
    if t <= 1e-10:
        print("   no strictly complementary certificate in this pool")
        return False
    # cross-check the authors' float kernel against the independent one on the used rows
    clg = cl2_grid_float(q * L)
    worst = 0.0
    for i in S:
        nv, s = tags[i]
        kt = kt_float(nv, N, s, clg)
        mine = np.array([kt[u] + kt[(L - u) % L] for u in Ut])
        worst = max(worst, np.abs(mine - rows[i]).max() / max(1.0, np.abs(mine).max()))
    print(f"   authors' float kernel vs independent kernel on used rows: max rel diff {worst:.1e}")
    out = verify(N, d, [tags[i] for i in S], [y[i] for i in S], lam, prec)
    print(f"   total {time.time() - t0:.1f}s", flush=True)
    return out['ok']


def main():
    mode = sys.argv[1]
    N, d = int(sys.argv[2]), int(sys.argv[3])
    npool = int(sys.argv[4]) if len(sys.argv) > 4 else 40
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 2
    prec = int(sys.argv[6]) if len(sys.argv) > 6 else 128
    if mode == "own":
        mode_own(N, d, npool, seed, prec)
    else:
        mode_theirs(N, d, npool, seed, prec)


if __name__ == "__main__":
    main()
