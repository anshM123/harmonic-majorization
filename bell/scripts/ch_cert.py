"""
ch_cert.py -- rigorous (interval-arithmetic) certificate for the chained cell LP (math/05-chained-bell.md, Theorem 5).

For given N (settings) and d (outcomes), with L = 2Nd and rational cells ell = n/q:
  1. floating-point LP of ch_cone.py on a pool of cells -> dual solution (y >= 0 on cell rows, lambda on the
     residue-sum equalities) and the set T of tight dual constraints (it must contain u = 1..d-1, the support of the
     window pair counts);
  2. interval arithmetic (mpmath.iv, outward rounding): every Clausen value Cl2(2 pi j/(qL)) is enclosed by the series of
     verify_cone.py (Bernoulli numbers, explicit tail; github.com/anshM123/IQOQI-OQP-27, cone-certificates/); the cell
     coefficients kt^{ell,s}(u), the objective coefficients and the window values are enclosed;
  3. the equalities  sum_k y_k A_k[u] + sum_r lambda_r E_r[u] = obj[u]  (u in T) are solved exactly in an interval box by
     the least-norm residual correction; it is verified that every y_k > 0 and that the dual constraints with u not in T
     hold strictly.
Then for every Q-configuration of the chained clock model: F_N(V) <= F_nat (weak duality + complementary slackness;
see math/05-chained-bell.md).  Prints CERTIFIED or FAILED.
Usage: python ch_cert.py N d [npool] [q] [seed] [prec]
"""
import sys
import time

import mpmath
import numpy as np
from fractions import Fraction
from mpmath import iv
from scipy.optimize import linprog

from ch_clock import K_N, F_natural
from ch_cone import kt_vec, Phi_c, random_cells, window_C, setup, fold


# -------------------------------------------------------------------------------- rigorous Clausen (as verify_cone.py)
def bern_coeffs(n0):
    cs = []
    for n in range(1, n0 + 1):
        p, qq = mpmath.bernfrac(2 * n)
        B = Fraction(int(p), int(qq))
        f = Fraction(1)
        for k in range(2, 2 * n + 2):
            f *= k
        cs.append(abs(B) / (2 * n * f))
    return cs


N0 = 70
CN = None
CN_IV = None


def init_iv(prec):
    global CN, CN_IV
    iv.prec = prec
    CN = bern_coeffs(N0)
    CN_IV = [iv.mpf(c.numerator) / iv.mpf(c.denominator) for c in CN]


def cl2_iv(t):
    """enclosure of Cl2(t) for an interval 0 < t <= pi."""
    t2 = t * t
    s = t - t * iv.log(t)
    pw = t * t2
    for c in CN_IV:
        s += c * pw
        pw *= t2
    tail_hi = t * (iv.pi ** 2 / 6) * iv.mpf(4) ** (-(N0 + 1)) / (iv.mpf(3) / 4 * (N0 + 1) * (2 * N0 + 3))
    return s + iv.mpf([0, 1]) * tail_hi


def cl2_grid_full(G):
    """Cl2(2 pi j / G) for j = 0..G-1 (enclosures); uses oddness and 2 pi periodicity."""
    out = [None] * G
    out[0] = iv.mpf(0)
    for j in range(1, G):
        if 2 * j < G:
            out[j] = cl2_iv(2 * iv.pi * j / G)
        elif 2 * j == G:
            out[j] = iv.mpf(0)
    for j in range(1, G):
        if 2 * j > G:
            out[j] = -out[G - j]
    return out


# -------------------------------------------------------------------------------- interval cell coefficients
def kt_iv(nvec, q, N, s, clg):
    d = len(nvec)
    L = 2 * N * d
    G = q * L
    pref = -(iv.mpf(L) ** 2) / (2 * iv.pi ** 2)
    ext = np.tile(np.asarray(nvec, int), 2 * N)
    Lc = np.concatenate([[0], np.cumsum(ext)])          # in units of 1/q
    a, b = Lc[:-1], Lc[1:]

    def Xi(j):
        return pref * clg[int(j) % G]
    Kbar = [iv.mpf(0)] * L
    for x in range(L):
        for w in range(1, L):
            y = (x + w) % L
            K = Xi(b[x] - a[y]) - Xi(a[x] - a[y]) - Xi(b[x] - b[y]) + Xi(a[x] - b[y])
            Kbar[w] += K
    kt = [iv.mpf(0)] * L
    for j in range(s):
        for jp in range(s, N):
            sh = (jp - j) * d
            for u in range(L):
                kt[u] += Kbar[(u - sh) % L]
    return [k / L for k in kt]


def phi_c_iv(a, b, clg_fine=None):
    """Phi_c(a, b) for rational a, b (interval)."""
    def cl(t):
        # t = 2 pi r, r rational in (0,1)
        r = t
        if r == 0 or r == Fraction(1, 2):
            return iv.mpf(0)
        if r < Fraction(1, 2):
            return cl2_iv(2 * iv.pi * iv.mpf(r.numerator) / r.denominator)
        return -cl2_iv(2 * iv.pi * iv.mpf((1 - r).numerator) / (1 - r).denominator)
    return (cl(a) + cl(b) + cl(1 - a - b)) / (2 * iv.pi ** 2)


def KN_iv(N, d, u):
    return N * iv.sin(iv.pi * (N - 1) * u / (N * d)) / iv.sin(iv.pi * u / d)


def mid(x):
    return float(mpmath.mpf(x.mid))


def certify(N, d, npool=40, q=8, seed=2, prec=120, verbose=True):
    t0 = time.time()
    rng = np.random.default_rng(seed)
    cells_int = [np.ones(d, int) * q] + [np.round(c * q).astype(int) for c in random_cells(d, npool, rng, q=q)]
    L, U = setup(N, d)
    nU = len(U)
    # ---------------- float data (ch_cone formulation) ----------------
    obj = np.zeros(L)
    for u in range(L):
        if u % d:
            obj[u] = 0.5 * K_N(N, d, u)
    objf = fold(obj, L, U)
    Efl = []
    for r in range(1, d):
        row = np.zeros(L)
        for t in range(2 * N):
            row[(r + t * d) % L] += 1.0
        Efl.append(fold(row, L, U))
    Efl = np.array(Efl)
    rows, tags = [], []
    for ci, nv in enumerate(cells_int):
        ell = nv / q
        for s in range(1, N // 2 + 1):
            rows.append(fold(kt_vec(ell, N, s), L, U))
            tags.append((ci, s))
    Afl = np.array(rows)
    nrow = len(rows)
    must = [i for i, u in enumerate(U) if u < d]
    off = [i for i, u in enumerate(U) if u >= d]
    # ---------------- LP2: strictly complementary dual certificate ----------------
    # variables: y (nrow, >= 0), lam (d-1, free), t (<= 1); maximise t
    nv_ = nrow + (d - 1) + 1
    cvec = np.zeros(nv_)
    cvec[-1] = -1.0
    Aeq = np.zeros((len(must), nv_))
    beq = np.zeros(len(must))
    for a_, i in enumerate(must):
        Aeq[a_, :nrow] = Afl[:, i]
        Aeq[a_, nrow:nrow + d - 1] = Efl[:, i]
        beq[a_] = objf[i]
    Aub = np.zeros((len(off), nv_))
    bub = np.zeros(len(off))
    for a_, i in enumerate(off):
        Aub[a_, :nrow] = -Afl[:, i]
        Aub[a_, nrow:nrow + d - 1] = -Efl[:, i]
        Aub[a_, -1] = 1.0 * (1 + abs(objf[i]))
        bub[a_] = -objf[i]
    bounds = [(0, None)] * nrow + [(None, None)] * (d - 1) + [(None, 1.0)]
    res = linprog(cvec, A_ub=Aub if len(off) else None, b_ub=bub if len(off) else None, A_eq=Aeq, b_eq=beq,
                  bounds=bounds, method="highs")
    if res.status != 0:
        print(f"N={N} d={d}: FAILED (LP2 status {res.status}: {res.message})")
        return False
    tstar = -res.fun
    y = res.x[:nrow]
    lam = res.x[nrow:nrow + d - 1]
    S = [i for i in range(nrow) if y[i] > 1e-9 * max(1.0, y.max())]
    Fn = F_natural(N, d)
    if verbose:
        print(f"N={N} d={d}: LP2 slack t* = {tstar:.3e}; cells used {len(S)} (of {nrow} rows); equations {len(must)}",
              flush=True)
    if tstar <= 1e-10:
        print(f"FAILED: no strictly complementary certificate in this pool (t* = {tstar:.2e})")
        return False
    # ---------------- interval data ----------------
    init_iv(prec)
    clg = cl2_grid_full(q * L)
    Aiv = {}
    for i in S:
        ci, s = tags[i]
        kt = kt_iv(cells_int[ci], q, N, s, clg)
        f = []
        for u in U:
            v = kt[u]
            if (L - u) % L != u:
                v = v + kt[(L - u) % L]
            f.append(v)
        bnd = iv.mpf(L) ** 2 * phi_c_iv(Fraction(s, 2 * N), Fraction(N - s, 2 * N))
        Aiv[i] = (f, bnd - kt[0] * d)
    obj_iv = []
    for u in U:
        v = KN_iv(N, d, u) / 2
        if (L - u) % L != u:
            v = v + KN_iv(N, d, L - u) / 2
        obj_iv.append(v)
    Ex = [[iv.mpf(int(round(Efl[r][i]))) for i in range(nU)] for r in range(d - 1)]
    nS = len(S)
    nx = nS + d - 1
    T = must
    nT = len(T)
    Gm = [[(Aiv[S[b_]][0][ti] if b_ < nS else Ex[b_ - nS][ti]) for b_ in range(nx)] for ti in T]
    bvec = [obj_iv[ti] for ti in T]
    x0 = [float(y[i]) for i in S] + [float(l) for l in lam]
    x0iv = [iv.mpf(v) for v in x0]
    resid = [bvec[a_] - sum(Gm[a_][b_] * x0iv[b_] for b_ in range(nx)) for a_ in range(nT)]
    Gmid = np.array([[mid(Gm[a_][b_]) for b_ in range(nx)] for a_ in range(nT)])
    Mm = Gmid @ Gmid.T
    R = np.linalg.inv(Mm)
    Riv = [[iv.mpf(float(R[i][j])) for j in range(nT)] for i in range(nT)]
    GG = [[sum((Gm[i][k] * Gm[j][k] for k in range(nx)), iv.mpf(0)) for j in range(nT)] for i in range(nT)]
    Emat = [[(iv.mpf(1) if i == j else iv.mpf(0)) - sum((Riv[i][k] * GG[k][j] for k in range(nT)), iv.mpf(0))
             for j in range(nT)] for i in range(nT)]
    normE = max(sum(float(abs(Emat[i][j]).b) for j in range(nT)) for i in range(nT))
    Rr = [sum((Riv[i][k] * resid[k] for k in range(nT)), iv.mpf(0)) for i in range(nT)]
    normRr = max(float(abs(v).b) for v in Rr)
    if normE >= 1:
        print("FAILED: approximate inverse not good enough", normE)
        return False
    wbound = normRr / (1 - normE) * (1 + 1e-12) + 1e-300
    box = iv.mpf([-wbound, wbound])
    wenc = [Rr[i] + sum((Emat[i][k] * box for k in range(nT)), iv.mpf(0)) for i in range(nT)]
    xiv = [x0iv[b_] + sum((Gm[a_][b_] * wenc[a_] for a_ in range(nT)), iv.mpf(0)) for b_ in range(nx)]
    ymin = min((float(xiv[b_].a) for b_ in range(nS)), default=float('inf'))
    margins = []
    for i in off:
        val = sum((xiv[b_] * Aiv[S[b_]][0][i] for b_ in range(nS)), iv.mpf(0)) + \
            sum((xiv[nS + r] * Ex[r][i] for r in range(d - 1)), iv.mpf(0)) - obj_iv[i]
        margins.append(float(val.a))
    feas = all(m > 0 for m in margins) if margins else True
    UB = sum((xiv[b_] * Aiv[S[b_]][1] for b_ in range(nS)), iv.mpf(0)) + \
        sum((xiv[nS + r] * d for r in range(d - 1)), iv.mpf(0))
    ok = (ymin > 0) and feas
    if verbose:
        print(f"   ||I - R GG||_inf <= {normE:.2e}; |correction| <= {wbound:.2e}; min y (lower end) {ymin:.3e}; "
              f"min off-window dual slack (lower end) {(min(margins) if margins else float('nan')):.3e}; "
              f"bound in [{float(UB.a):.13f}, {float(UB.b):.13f}], F_nat = {Fn:.13f}", flush=True)
        print(("CERTIFIED" if ok else "FAILED") + f"  N={N} d={d}  cells={len(set(tags[i][0] for i in S))}  "
              f"({time.time() - t0:.1f}s)", flush=True)
    return ok


def main():
    N = int(sys.argv[1])
    d = int(sys.argv[2])
    npool = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    q = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 2
    prec = int(sys.argv[6]) if len(sys.argv) > 6 else 120
    certify(N, d, npool, q, seed, prec)


if __name__ == "__main__":
    main()
