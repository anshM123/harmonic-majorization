"""Check S9 (priority): does Theorem D_sep follow from Z. Long, arXiv:2610.10158, Lemma 2.3?

Long's Lemma 2.3 (as read on arxiv.org/html/2610.10158v1): for T = A + iB (A, B Hermitian, d x d) and smooth real f,
    d f(0) + int Qf dmu_{A,B} - sum_j f(lambda_j(T)) = (1/2pi) int Psi_T Laplace(f) dm_2,
with mu_{A,B} Heinavaara's tracial joint spectral measure (int f(sx+ty) dmu = tr H(f)(sA+tB),
H(f)(u) = int_0^1 ((1-r)/r) f(ru) dr), Q = E(E+1), E = z.grad, and Psi_T(w) = sum_j log+|lambda_j(E^{-1}E^*)| >= 0,
E = T - w.  Since H(Q phi) = phi - phi(0) for ridge functions, for separable f = a(x) + c(y) this gives
    sum_j f(lambda_j(T)) = tr a(A) + tr c(B) + (1/2pi) int Psi_T (-a''(x) - c''(y)) dm_2.          (*)
Applied to T = B + ig and to the normal T_E = B + iE(g) (Psi = 0), (*) gives
    sum_X f - sum_XE f = [tr c(g) - tr c(E g)] + (1/2pi) int Psi_X (-a'' - c'') dm_2,
which is >= 0 when c is convex (pinching inequality) and a'' + c'' <= 0 on the numerical range: Theorem D_sep.
This script checks (*) numerically (independently of the paper) on random 3x3 and 4x4 pairs, with Psi_T computed
from its definition and a tensor Gauss-Legendre rule refined at the eigenvalues.
"""
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
from hmaudit import rand_herm, setup

rng = np.random.default_rng(41)
t0 = time.time()


def Psi(T, W):
    """Psi_T at an array of complex points W."""
    d = T.shape[0]
    E = T[None, :, :] - W[:, None, None] * np.eye(d)[None, :, :]
    Es = np.conj(np.transpose(E, (0, 2, 1)))
    Mx = np.linalg.solve(E, Es)
    lam = np.abs(np.linalg.eigvals(Mx))
    return np.sum(np.log(np.maximum(lam, 1.0)), axis=1)


def integrate_box(F, xs, ys, n=96):
    """int F over the box, split at the coordinates xs, ys (cosine-graded Gauss-Legendre in each cell)."""
    xg, wg = leggauss(n)
    th = (xg + 1) * np.pi / 2
    tot = 0.0
    for x0, x1 in zip(xs[:-1], xs[1:]):
        if x1 - x0 <= 1e-15:
            continue
        X = x0 + (x1 - x0) * (1 - np.cos(th)) / 2
        WX = wg * np.pi / 2 * (x1 - x0) / 2 * np.sin(th)
        for y0, y1 in zip(ys[:-1], ys[1:]):
            if y1 - y0 <= 1e-15:
                continue
            Y = y0 + (y1 - y0) * (1 - np.cos(th)) / 2
            WY = wg * np.pi / 2 * (y1 - y0) / 2 * np.sin(th)
            XX, YY = np.meshgrid(X, Y, indexing="ij")
            vals = F(XX.ravel() + 1j * YY.ravel()).reshape(XX.shape)
            tot += np.einsum("i,j,ij->", WX, WY, vals)
    return tot


def sep_funcs(rng):
    kap, s, y0 = rng.uniform(0.5, 2), rng.uniform(0.3, 1.5), rng.normal() * 0.5
    al, x0 = rng.uniform(0.5, 2), rng.normal() * 0.5
    a = lambda x: -np.exp(al * (x - x0))
    a2 = lambda x: -al ** 2 * np.exp(al * (x - x0))
    c = lambda y: kap * s ** 2 * np.log(np.cosh((y - y0) / s))
    c2 = lambda y: kap / np.cosh((y - y0) / s) ** 2
    return a, a2, c, c2


def fmat(fun, H):
    w, V = np.linalg.eigh(H)
    return np.real(np.sum(fun(w)))


for trial in range(6):
    M = 3 if trial < 3 else 4
    levels = np.sort(rng.normal(size=M))
    if trial == 2:
        levels[1] = levels[0]
    g = rand_herm(M, rng, rng.uniform(0.4, 1.5))
    c = setup(list(levels), g)
    T = c["X"]
    a, a2, cf, c2 = sep_funcs(rng)
    lhs = np.sum(a(c["nu"].real) + cf(c["nu"].imag))
    tr_part = fmat(a, np.diag(levels)) + fmat(cf, g)
    # box = Rect (contains the numerical range) slightly enlarged; split at eigenvalue coordinates
    xs = sorted(set([levels.min() - 1e-3, levels.max() + 1e-3] + list(c["nu"].real)))
    ys = sorted(set([c["lmin"] - 1e-3, c["lmax"] + 1e-3] + list(c["nu"].imag)))
    F = lambda W: Psi(T, W) * (-a2(W.real) - c2(W.imag))
    I1 = integrate_box(F, xs, ys, n=96) / (2 * np.pi)
    I2 = integrate_box(F, xs, ys, n=160) / (2 * np.pi)
    mass = integrate_box(lambda W: Psi(T, W), xs, ys, n=160)
    dep = np.sum(np.abs(T) ** 2) - np.sum(np.abs(c["nu"]) ** 2)
    print(f"case {trial} (M={M}): sum f(lambda) = {lhs:.8f};  tr a(B) + tr c(g) + (1/2pi) int Psi(-Lap f) = "
          f"{tr_part + I2:.8f} (n=160) / {tr_part + I1:.8f} (n=96);  diff {lhs - tr_part - I2:.2e}")
    print(f"      check of Psi normalisation: int Psi = {mass:.6f} vs (pi/2) dep(T)^2 = {np.pi/2*dep:.6f}")
    # consequence: D_sep defect = [tr c(g) - tr c(Eg)] + (1/2pi) int Psi_X (-a''-c'')
    D = lhs - np.sum(a(c["nuE"].real) + cf(c["nuE"].imag))
    pred = fmat(cf, g) - fmat(cf, c["gE"]) + I2
    print(f"      D_sep-type defect sum_X f - sum_XE f = {D:.8f};  pinching term + Psi term = {pred:.8f}")
print(f"time {time.time()-t0:.0f}s")
