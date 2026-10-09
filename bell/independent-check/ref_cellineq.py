"""
ref_cellineq.py -- independent checks of Lemma 2 and Proposition 4 of math/05-chained-bell.md (the printed labels
use the draft numbering II.2 and II.4).
 (1) Lemma 2(1): K_N(u) = -(1/2) sum_{k=1}^{2N-1} sin(pi k/N) cot(pi (u - k d)/L) at random real u;
 (2) Lemma 2(3)-(4): F_N = (1/2) sum_s mu_s <X_s, Y_s> - (d/2) sum_k sin(pi k/N) kappa(kd) on random configurations;
 (3) Prop. 4: adversarial maximisation (torch, L-BFGS) of the cell functional
       Q^{ell,s}(Q) = sum_{x != y} K^ell(x,y) tau((X_s)_x (Y_s)_y)
     over chained Q-configurations on C^M, against the claimed bound L^2 Phi_c(s/2N, (N-s)/2N).
Usage: python ref_cellineq.py
"""
import numpy as np
import torch
import mpmath

torch.set_num_threads(2)
torch.set_default_dtype(torch.float64)
rng = np.random.default_rng(7)


def KN(N, d, u):
    return N * np.sin(np.pi * (N - 1) * u / (N * d)) / np.sin(np.pi * u / d)


def check_partial_fraction():
    worst = 0.0
    for N in range(2, 9):
        for d in range(2, 9):
            L = 2 * N * d
            for _ in range(20):
                u = rng.uniform(0.1, L - 0.1)
                if abs(np.sin(np.pi * u / d)) < 1e-3:
                    continue
                rhs = -0.5 * sum(np.sin(np.pi * k / N) / np.tan(np.pi * (u - k * d) / L) for k in range(1, 2 * N))
                worst = max(worst, abs(KN(N, d, u) - rhs) / max(1, abs(rhs)))
    print(f"(1) partial fractions: max rel err {worst:.1e}")


def random_qconf(N, d, M):
    """chained Q-configuration: Q_{k - d a} = U_k diag(1[lab = a]) U_k^*, a in Z_2N."""
    L = 2 * N * d
    Q = [None] * L
    for k in range(d):
        Z = rng.standard_normal((M, M)) + 1j * rng.standard_normal((M, M))
        U, _ = np.linalg.qr(Z)
        lab = rng.integers(0, 2 * N, size=M)
        for a in range(2 * N):
            Q[(k - d * a) % L] = U @ np.diag((lab == a).astype(float)) @ U.conj().T
    return Q


def check_decomposition():
    worst = 0.0
    for N in range(2, 8):
        for d in range(2, 6):
            L = 2 * N * d
            kap = lambda u: 0.0 if u % L == 0 else 1 / np.tan(np.pi * u / L)
            for _ in range(3):
                M = 3
                Q = random_qconf(N, d, M)
                tr = np.array([[np.real(np.trace(Q[x] @ Q[y])) / M for y in range(L)] for x in range(L)])
                F = 0.5 * sum(tr[x, y] * KN(N, d, x - y) for x in range(L) for y in range(L) if (x - y) % d)
                mus = {}
                for s in range(1, N // 2 + 1):
                    mu = 2 * (1 - np.cos(np.pi / N)) * np.sin(np.pi * s / N)
                    if 2 * s == N:
                        mu /= 2
                    mus[s] = mu
                # check sin decomposition
                for k in range(1, N):
                    val = sum(mus[s] * min(k, N - k, s, N - s) for s in mus)
                    worst = max(worst, abs(val - np.sin(np.pi * k / N)))
                tot = 0.0
                for s, mu in mus.items():
                    pair = 0.0
                    for x in range(L):
                        for y in range(L):
                            # <X_s, Y_s> = sum_{x,y} kappa(x-y) tau((X_s)_x (Y_s)_y)
                            for j in range(s):
                                for jp in range(s, N):
                                    pair += kap(x - y) * tr[(x + j * d) % L, (y + jp * d) % L]
                    tot += 0.5 * mu * pair
                const = (d / 2) * sum(np.sin(np.pi * k / N) * kap(k * d) for k in range(1, N))
                worst = max(worst, abs(F - (tot - const)) / max(1, abs(F)))
    print(f"(2) Lemma II.2 decomposition on random configurations: max rel err {worst:.1e}")


def Xi_fn(L):
    return lambda v: -(L ** 2 / (2 * np.pi ** 2)) * float(mpmath.clsin(2, 2 * mpmath.pi * v / L))


def cell_K(ell, N):
    d = len(ell)
    L = 2 * N * d
    Xi = Xi_fn(L)
    Lc = np.concatenate([[0.0], np.cumsum(np.tile(ell, 2 * N))])
    K = np.zeros((L, L))
    for x in range(L):
        for y in range(L):
            if x != y:
                K[x, y] = Xi(Lc[x + 1] - Lc[y]) - Xi(Lc[x] - Lc[y]) - Xi(Lc[x + 1] - Lc[y + 1]) + Xi(Lc[x] - Lc[y + 1])
    return K


def adversarial(N, d, s, ell, M, restarts):
    L = 2 * N * d
    K = torch.tensor(cell_K(np.asarray(ell, float), N))
    bound = L ** 2 * (float(mpmath.clsin(2, mpmath.pi * s / N)) + float(mpmath.clsin(2, mpmath.pi * (N - s) / N))) \
        / (2 * np.pi ** 2)
    best = -np.inf
    for r in range(restarts):
        labs = [rng.integers(0, 2 * N, size=M) if r % 2 else np.arange(M) % (2 * N) for _ in range(d)]
        params = [torch.tensor(rng.standard_normal((M, M)) + 1j * rng.standard_normal((M, M)),
                               dtype=torch.complex128, requires_grad=True) for _ in range(d)]

        def functional():
            Q = [None] * L
            for k in range(d):
                U = torch.linalg.matrix_exp(params[k] - params[k].conj().T)
                for a in range(2 * N):
                    e = torch.tensor((labs[k] == a).astype(float), dtype=torch.complex128)
                    Q[(k - d * a) % L] = (U * e) @ U.conj().T
            X = [sum(Q[(x + j * d) % L] for j in range(s)) for x in range(L)]
            Y = [sum(Q[(x + j * d) % L] for j in range(s, N)) for x in range(L)]
            Xs = torch.stack(X)
            Ys = torch.stack(Y)
            T = torch.einsum('xij,yji->xy', Xs, Ys).real / M
            return (K * T).sum()
        opt = torch.optim.LBFGS(params, lr=1.0, max_iter=300, line_search_fn='strong_wolfe', tolerance_grad=1e-12,
                                tolerance_change=1e-15)

        def closure():
            opt.zero_grad()
            v = -functional()
            v.backward()
            return v
        for _ in range(4):
            opt.step(closure)
        best = max(best, float(functional().detach()))
    print(f"(3) N={N} d={d} s={s} ell={list(np.round(ell, 3))} M={M}: max cell functional {best:.10f}  "
          f"bound L^2 Phi_c = {bound:.10f}  diff {best - bound:+.2e}", flush=True)


def main():
    check_partial_fraction()
    check_decomposition()
    for (N, d, s, ell, M) in [(4, 3, 2, [0.5, 1.75, 0.75], 2), (4, 3, 2, [1, 1, 1], 3), (5, 3, 2, [2.25, 0.25, 0.5], 2),
                              (6, 3, 3, [0.375, 1.5, 1.125], 2), (6, 4, 2, [1.5, 0.5, 1.25, 0.75], 2),
                              (3, 4, 1, [0.25, 2.0, 1.0, 0.75], 3), (5, 4, 2, [1, 1, 1, 1], 2)]:
        adversarial(N, d, s, ell, M, restarts=8)


if __name__ == "__main__":
    main()
