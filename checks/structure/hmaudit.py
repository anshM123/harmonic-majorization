"""Independent helpers for the pre-publication check of math/03-structure.md (written from scratch for this check;
no code shared with the authors' scripts).

Conventions (as in 03-structure.md):
  B = diag(levels) (WLOG by unitary invariance), g Hermitian, E = pinching onto the eigenspaces of B,
  X = B + i g, X_E = B + i E(g), K = g - E(g).
  rho_B(s, beta) = (1/2pi) sum_i |Im xi_i|, xi_i the roots of det(g - s - xi (B - beta)), i.e. the eigenvalues of
  (B - beta)^{-1} (g - s).  R = push-forward of rho_B ds dbeta under (s, beta) -> beta + i s.
Pinchings use exact level labels (levels are given as exact floats with exact repeats; no isclose).
"""
import numpy as np
import mpmath as mp
from numpy.polynomial.legendre import leggauss


# ----------------------------------------------------------------------------------------------------------------
# matrices
# ----------------------------------------------------------------------------------------------------------------
def rand_herm(M, rng, scale=1.0, real=False):
    A = rng.normal(size=(M, M))
    if not real:
        A = A + 1j * rng.normal(size=(M, M))
    H = (A + A.conj().T) / 2
    return H.astype(complex) * scale


def labels_of(levels):
    uniq = sorted(set(float(x) for x in levels))
    lab = np.array([uniq.index(float(x)) for x in levels])
    return lab, np.array(uniq)


def pinch(g, lab):
    out = np.zeros_like(g)
    for L in np.unique(lab):
        idx = np.where(lab == L)[0]
        out[np.ix_(idx, idx)] = g[np.ix_(idx, idx)]
    return out


def setup(levels, g):
    levels = np.asarray(levels, float)
    lab, uniq = labels_of(levels)
    gE = pinch(g, lab)
    B = np.diag(levels).astype(complex)
    X = B + 1j * g
    nu = np.linalg.eigvals(X)
    nuE = []
    for k, b in enumerate(uniq):
        idx = np.where(lab == k)[0]
        ev = np.linalg.eigvalsh(g[np.ix_(idx, idx)])
        nuE += [b + 1j * e for e in ev]
    nuE = np.array(nuE)
    K = g - gE
    ev_g = np.linalg.eigvalsh(g)
    return dict(levels=levels, lab=lab, uniq=uniq, g=g, gE=gE, B=B, X=X, nu=nu, nuE=nuE, K=K,
                K2=float(np.sum(np.abs(K) ** 2)), lmin=ev_g[0], lmax=ev_g[-1], M=len(levels))


def mu_B(c, beta):
    """x-marginal step function sum_{j<k} ||Q_j g Q_k||^2/(b_k-b_j) 1[b_j<beta<b_k]."""
    tot = 0.0
    u, lab, g = c["uniq"], c["lab"], c["g"]
    for j in range(len(u)):
        for k in range(j + 1, len(u)):
            if u[j] < beta < u[k]:
                Ij, Ik = np.where(lab == j)[0], np.where(lab == k)[0]
                tot += np.sum(np.abs(g[np.ix_(Ij, Ik)]) ** 2) / (u[k] - u[j])
    return tot


# ----------------------------------------------------------------------------------------------------------------
# density rho_B and quadrature nodes for int int F(s, beta) rho_B ds dbeta
# ----------------------------------------------------------------------------------------------------------------
def pencil_roots(c, s, beta):
    """roots xi of det(g - s - xi (B - beta)) for an array of s (beta scalar). shape (len(s), M)."""
    s = np.atleast_1d(s)
    d = 1.0 / (c["levels"] - beta)
    A = c["g"][None, :, :] - s[:, None, None] * np.eye(c["M"])[None, :, :]
    A = d[None, :, None] * A
    return np.linalg.eigvals(A)


def rho_vals(c, s, beta):
    return np.sum(np.abs(pencil_roots(c, s, beta).imag), axis=1) / (2 * np.pi)


def rho_pairs(c, s, beta):
    """rho_B at arrays (s_i, beta_i) of equal length (batched)."""
    s = np.atleast_1d(np.asarray(s, float)); beta = np.atleast_1d(np.asarray(beta, float))
    d = 1.0 / (c["levels"][None, :] - beta[:, None])
    A = c["g"][None, :, :] - s[:, None, None] * np.eye(c["M"])[None, :, :]
    A = d[:, :, None] * A
    return np.sum(np.abs(np.linalg.eigvals(A).imag), axis=1) / (2 * np.pi)


def _count_nonreal(c, s, beta, tol=1e-7):
    r = pencil_roots(c, s, beta)
    return np.sum(np.abs(r.imag) > tol * (1 + np.abs(r)), axis=1)


def _critical_values(c, beta, nxi=4000, iters=50):
    """critical values s* = lambda_j(xi*) of the eigenvalue curves of the Hermitian family g - xi (B - beta),
    xi real (lambda_j'(xi*) = -<v_j, (B - beta) v_j> = 0).  These are the only s where the number of real
    pencil roots can change (beta not in spec B)."""
    D = c["levels"] - beta
    scale = (np.max(np.abs([c["lmin"], c["lmax"]])) + 1.0) / np.min(np.abs(D))
    phi = np.linspace(-np.pi / 2, np.pi / 2, nxi + 2)[1:-1]
    xi = scale * np.tan(phi)
    Hs = c["g"][None, :, :] - xi[:, None, None] * np.diag(D)[None, :, :]
    lam, V = np.linalg.eigh(Hs)
    der = -np.einsum("nij,i,nij->nj", V.conj(), D, V).real
    out = []
    for j in range(c["M"]):
        sg = np.sign(der[:, j])
        idx = np.where(sg[:-1] * sg[1:] < 0)[0]
        for i in idx:
            a, b = phi[i], phi[i + 1]
            fa = der[i, j]
            for _ in range(iters):
                m = 0.5 * (a + b)
                Hm = c["g"] - scale * np.tan(m) * np.diag(D)
                lm, Vm = np.linalg.eigh(Hm)
                fm = -np.real(Vm[:, j].conj() @ (D * Vm[:, j]))
                if np.sign(fm) == np.sign(fa):
                    a, fa = m, fm
                else:
                    b = m
            Hm = c["g"] - scale * np.tan(0.5 * (a + b)) * np.diag(D)
            out.append(np.linalg.eigvalsh(Hm)[j])
    return out


def _breakpoints_s(c, beta, nscan=600):
    a, b = c["lmin"], c["lmax"]
    th = np.linspace(0, np.pi, nscan)
    s = a + (b - a) * (1 - np.cos(th)) / 2
    cnt = _count_nonreal(c, s, beta)
    bps = []
    for i in range(nscan - 1):
        if cnt[i] != cnt[i + 1]:
            lo, hi, clo = s[i], s[i + 1], cnt[i]
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if _count_nonreal(c, np.array([mid]), beta)[0] == clo:
                    lo = mid
                else:
                    hi = mid
                if hi - lo < 1e-14 * (1 + abs(mid)):
                    break
            bps.append(0.5 * (lo + hi))
    crit = [x for x in _critical_values(c, beta) if a < x < b]
    allb = sorted(set(np.round(bps + crit, 13)))
    # merge points closer than 1e-12
    merged = []
    for x in allb:
        if not merged or x - merged[-1] > 1e-12:
            merged.append(x)
    return [a] + merged + [b]


_XG, _WG = leggauss(16)


def adaptive_nodes_1d(fvec, a, b, tol=1e-12, maxit=40, nmin=4):
    """Adaptive Gauss-Legendre (16 nodes per cell) on [a,b] after the substitution x = a + (b-a)(1-cos phi)/2
    (removes square-root endpoint singularities).  fvec maps an array of x to an array (n,) or (n,k) of values;
    refinement is driven by all columns.  Returns (x nodes, weights in x)."""
    def cell(p0, p1):
        ph = 0.5 * (p0 + p1) + 0.5 * (p1 - p0) * _XG
        wp = 0.5 * (p1 - p0) * _WG
        x = a + (b - a) * (1 - np.cos(ph)) / 2
        w = wp * (b - a) / 2 * np.sin(ph)
        return x, w
    cells = list(zip(np.linspace(0, np.pi, nmin + 1)[:-1], np.linspace(0, np.pi, nmin + 1)[1:]))
    done_x, done_w = [], []
    for _ in range(maxit):
        if not cells:
            break
        # evaluate parent and children for all open cells at once
        X, Wt, tags = [], [], []
        for ci, (p0, p1) in enumerate(cells):
            pm = 0.5 * (p0 + p1)
            for tag, (q0, q1) in enumerate([(p0, p1), (p0, pm), (pm, p1)]):
                x, w = cell(q0, q1)
                X.append(x); Wt.append(w); tags.append((ci, tag))
        Xall = np.concatenate(X)
        F = np.asarray(fvec(Xall))
        if F.ndim == 1:
            F = F[:, None]
        F = F.reshape(len(X), len(_XG), -1)
        newcells = []
        for ci, (p0, p1) in enumerate(cells):
            par = np.einsum("n,nk->k", Wt[3 * ci], F[3 * ci])
            ch = np.einsum("n,nk->k", Wt[3 * ci + 1], F[3 * ci + 1]) + np.einsum("n,nk->k", Wt[3 * ci + 2], F[3 * ci + 2])
            scale = 1.0 + np.abs(ch)
            if np.all(np.abs(par - ch) <= tol * scale) or (p1 - p0) < 1e-13:
                done_x += [X[3 * ci + 1], X[3 * ci + 2]]
                done_w += [Wt[3 * ci + 1], Wt[3 * ci + 2]]
            else:
                pm = 0.5 * (p0 + p1)
                newcells += [(p0, pm), (pm, p1)]
        cells = newcells
    for (p0, p1) in cells:  # not converged: keep finest rule
        x, w = cell(p0, p1)
        done_x.append(x); done_w.append(w)
    return np.concatenate(done_x), np.concatenate(done_w)


def inner_nodes(c, beta, nscan=600, tol=1e-12):
    """adaptive s-nodes at fixed beta: returns (s, weight*rho)."""
    pts = _breakpoints_s(c, beta, nscan)
    SS, WW = [], []
    for a, b in zip(pts[:-1], pts[1:]):
        if b - a <= 0:
            continue
        f = lambda x: np.stack([rho_vals(c, x, beta), x * rho_vals(c, x, beta)], axis=1)
        ss, sw = adaptive_nodes_1d(f, a, b, tol=tol)
        SS.append(ss)
        WW.append(sw * rho_vals(c, ss, beta))
    if not SS:
        return np.zeros(0), np.zeros(0)
    return np.concatenate(SS), np.concatenate(WW)


def rho_nodes(c, tol_beta=1e-10, tol_s=1e-12, nscan=600):
    """Quadrature nodes (s_i, beta_i, W_i), W_i = weight * rho_B(s_i, beta_i), adaptive in both variables.
    beta: adaptive Gauss-Legendre after a cosine substitution on each gap (refinement driven by the inner moments
    int s^k rho ds, k = 0..3, and int cos(s) rho ds); s: adaptive on pieces between breakpoints."""
    cache = {}

    def inner(beta):
        key = float(beta)
        if key not in cache:
            cache[key] = inner_nodes(c, beta, nscan, tol_s)
        return cache[key]

    def fvec(betas):
        out = []
        for b in betas:
            ss, ww = inner(b)
            out.append([np.sum(ww), np.sum(ww * ss), np.sum(ww * ss ** 2), np.sum(ww * ss ** 3),
                        np.sum(ww * np.cos(ss)), np.sum(ww * np.sin(2 * ss))])
        return np.array(out)

    S, BT, W = [], [], []
    u = c["uniq"]
    for k in range(len(u) - 1):
        bx, bw = adaptive_nodes_1d(fvec, u[k], u[k + 1], tol=tol_beta, nmin=4)
        for beta, wbeta in zip(bx, bw):
            ss, ww = inner(beta)
            S.append(ss)
            BT.append(np.full_like(ss, beta))
            W.append(ww * wbeta)
    return np.concatenate(S), np.concatenate(BT), np.concatenate(W)


# ----------------------------------------------------------------------------------------------------------------
# strip sweeps: closed form of h_lambda and stop-loss transforms (unit strip 0 <= x <= 1)
# ----------------------------------------------------------------------------------------------------------------
def h_closed_mp(nu, lam):
    """h_lambda(nu) = int (w - lam)_+ omega_nu(dw), left edge of the unit strip, nu = x + i y, 0 <= x <= 1.
    Closed form derived for this check from K_x(u) = sum_n e^{-n pi |u|} sin(n pi x):
      h = (1/pi^2) Im Li2(e^{-pi|y - lam|} e^{i pi x}) + (y - lam)_+ (1 - x)."""
    x, y = mp.re(nu), mp.im(nu)
    lam = mp.mpf(lam)
    if x <= 0:
        return max(y - lam, 0)
    if x >= 1:
        return mp.mpf(0)
    t = y - lam
    q = mp.exp(-mp.pi * abs(t)) * mp.expj(mp.pi * x)
    return mp.im(mp.polylog(2, q)) / mp.pi ** 2 + max(t, 0) * (1 - x)


def h_quad_mp(nu, lam):
    """same quantity by direct quadrature of the strip Poisson kernel (for validating h_closed_mp)."""
    x, y = mp.re(nu), mp.im(nu)
    lam = mp.mpf(lam)

    def K(u):
        return mp.sin(mp.pi * x) / (2 * (2 * mp.sinh(mp.pi * u / 2) ** 2 + 2 * mp.sin(mp.pi * x / 2) ** 2))

    # int_lam^inf (w - lam) K(y - w) dw = int_{c}^{inf} (u - c) K(u) du with c = lam - y
    cc = lam - y
    return mp.quad(lambda u: (u - cc) * K(u), [cc, cc + 1, cc + 5, mp.inf])


def stoploss_left(nus, lam, prec=None):
    return mp.fsum(h_closed_mp(mp.mpc(z), lam) for z in nus)


def stoploss_right(nus, lam):
    return mp.fsum(h_closed_mp(1 - mp.conj(mp.mpc(z)), lam) for z in nus)


def mp_eigs(M_mp):
    E, _ = mp.eig(M_mp)
    return E


def mp_spectra(levels, g_exact, dps=60):
    """eigenvalues of B + i g and of B + i E(g) at high precision; g_exact: list of lists of mp numbers."""
    mp.mp.dps = dps
    M = len(levels)
    lab, uniq = labels_of(levels)
    X = mp.matrix(M, M)
    for i in range(M):
        for j in range(M):
            X[i, j] = (mp.mpf(levels[i]) if i == j else 0) + 1j * mp.mpc(g_exact[i][j])
    nu = mp_eigs(X)
    XE = mp.matrix(M, M)
    for i in range(M):
        for j in range(M):
            XE[i, j] = (mp.mpf(levels[i]) if i == j else 0) + (1j * mp.mpc(g_exact[i][j]) if lab[i] == lab[j] else 0)
    nuE = mp_eigs(XE)
    return nu, nuE
