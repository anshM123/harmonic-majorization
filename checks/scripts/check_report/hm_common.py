"""Shared helpers for the adversarial check of THEOREMS.md (harmonic majorization).
Unitary invariance: B = diag(levels) WLOG, g random Hermitian. P = 1 - B."""
import numpy as np
from scipy.linalg import expm
import mpmath as mp


def make_case(levels, rng, gscale=1.0, g=None):
    levels = np.asarray(levels, float)
    M = len(levels)
    if g is None:
        X = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
        g = (X + X.conj().T) / 2 * gscale
    pdiag = 1.0 - levels
    pvals = np.unique(np.round(pdiag, 13))
    idx = [np.where(np.abs(pdiag - p) <= 1e-13)[0] for p in pvals]
    Qs = []
    for I in idx:
        Q = np.zeros((M, M))
        Q[I, I] = 1.0
        Qs.append(Q)
    gd = sum(Q @ g @ Q for Q in Qs)
    B = np.diag(levels).astype(complex)
    P = np.diag(pdiag).astype(complex)
    lmin, lmax = np.linalg.eigvalsh(g)[[0, -1]]
    return dict(B=B, P=P, g=g, gd=gd, pvals=pvals, idx=idx, Qs=Qs, M=M, pdiag=pdiag,
                levels=levels, lmin=lmin, lmax=lmax)


def roots(c, s, tau):
    """roots xi of det(g - s - xi (P - tau)) = eigenvalues of (P - tau)^{-1} (g - s)."""
    M = c["M"]
    A = (c["g"] - s * np.eye(M)) / (c["pdiag"] - tau)[:, None]
    return np.linalg.eigvals(A)


def rho(c, s, tau):
    if np.any(np.isclose(c["pvals"], tau, atol=0, rtol=0)):
        return 0.0
    return np.sum(np.abs(roots(c, s, tau).imag)) / (2 * np.pi)


def mfun(c, tau):
    tot = 0.0
    pv, idx, g = c["pvals"], c["idx"], c["g"]
    for j in range(len(pv)):
        for k in range(j + 1, len(pv)):
            if pv[j] < tau < pv[k]:
                blk = g[np.ix_(idx[j], idx[k])]
                tot += np.sum(np.abs(blk) ** 2) / (pv[k] - pv[j])
    return tot


def offdiag_weights(c):
    """list of (p_j, p_k, ||Q_j g Q_k||_F^2) for j<k."""
    out = []
    pv, idx, g = c["pvals"], c["idx"], c["g"]
    for j in range(len(pv)):
        for k in range(j + 1, len(pv)):
            blk = g[np.ix_(idx[j], idx[k])]
            out.append((pv[j], pv[k], float(np.sum(np.abs(blk) ** 2))))
    return out


def Dfun(c, a, t):
    return np.trace(expm(a * c["g"] - t * c["P"])) - np.trace(expm(a * c["gd"] - t * c["P"]))


def Kx(x, u):
    """left-edge Poisson kernel of the strip, 0<x<1."""
    # cosh(pi u) - cos(pi x) = 2 sinh^2(pi u/2) + 2 sin^2(pi x/2)  (no cancellation for small u, x)
    su = np.sinh(np.clip(np.pi * np.asarray(u) / 2, -350, 350))
    return np.sin(np.pi * x) / (4 * (su * su + np.sin(np.pi * x / 2) ** 2))


def h_lambda(nu, lam, dps=30):
    """int (w-lam)_+ omega_nu(dw), omega_nu = left-edge harmonic measure from nu (closed strip)."""
    x, y = float(np.clip(nu.real, 0.0, 1.0)), float(nu.imag)
    if x <= 1e-14:
        return max(y - lam, 0.0)
    if x >= 1 - 1e-14:
        return 0.0
    with mp.workdps(dps):
        z = mp.exp(mp.pi * (y - lam)) * mp.exp(-1j * mp.pi * x)
        return float(-mp.im(mp.polylog(2, z)) / mp.pi ** 2)


def W_left(c, lam, which="L"):
    """sum h_lam over spec(B+ig) minus over spec(B+i g_d); which='R' uses nu -> 1 - conj(nu)."""
    nus = np.linalg.eigvals(c["B"] + 1j * c["g"])
    mus = np.linalg.eigvals(c["B"] + 1j * c["gd"])
    if which == "R":
        nus = 1 - nus.conj()
        mus = 1 - mus.conj()
    return sum(h_lambda(n, lam) for n in nus) - sum(h_lambda(m_, lam) for m_ in mus)
