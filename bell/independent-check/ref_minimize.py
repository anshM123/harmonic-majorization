"""
ref_minimize.py -- independent search for counterexamples to Theorem 5 of math/05-chained-bell.md at the level of
strategies
(independent of the clock model): minimise the BKP expression I_N over projective strategies on Phi_D,
parametrised by 2N unitaries (exponential of skew-Hermitian matrices) and fixed outcome labels.
Prints the best value found against I_N^BKP(d).
Usage: python ref_minimize.py N d D restarts [seed]
"""
import sys
import time

import numpy as np
import torch

torch.set_num_threads(2)
torch.set_default_dtype(torch.float64)


def bkp_value(N, d):
    return N * (d - 1) - (2 * N / d) * sum((d - m) * np.sin(np.pi * (N - 1) * m / (N * d)) / np.sin(np.pi * m / d)
                                           for m in range(1, d))


def cost_mats(d):
    a = np.arange(d)[:, None]
    b = np.arange(d)[None, :]
    c1 = (a - b) % d          # [A_k - B_k]
    c2 = (b - a) % d          # [B_k - A_{k+1}]
    c3 = (b - a - 1) % d      # [B_N - A_1 - 1]
    return [torch.tensor(c, dtype=torch.float64) for c in (c1, c2, c3)]


def make_proj(A, labels, d):
    """A: complex D x D parameter; U = exp(A - A^H); returns tensor (d, D, D) of projectors."""
    S = A - A.conj().T
    U = torch.linalg.matrix_exp(S)
    D = U.shape[0]
    E = torch.zeros((d, D), dtype=torch.complex128)
    for i, l in enumerate(labels):
        E[l, i] = 1.0
    # P_a = U diag(E_a) U^H
    return torch.einsum('im,am,jm->aij', U, E, U.conj())


def value(params, labs, N, d, C):
    D = params[0].shape[0]
    P = [make_proj(params[k], labs[k], d) for k in range(N)]
    Q = [make_proj(params[N + k], labs[N + k], d) for k in range(N)]

    def joint(Pk, Ql):
        # p(a,b) = Tr(P_a Q_b^T)/D = sum_ij P_a[i,j] Q_b[i,j] / D
        return torch.einsum('aij,bij->ab', Pk, Ql).real / D
    c1, c2, c3 = C
    tot = 0.0
    for k in range(N):
        tot = tot + (joint(P[k], Q[k]) * c1).sum()
        if k < N - 1:
            tot = tot + (joint(P[k + 1], Q[k]) * c2).sum()
    tot = tot + (joint(P[0], Q[N - 1]) * c3).sum()
    return tot


def run(N, d, D, restarts, seed=0):
    rng = np.random.default_rng(seed)
    C = cost_mats(d)
    target = bkp_value(N, d)
    best = np.inf
    best_info = None
    t0 = time.time()
    for r in range(restarts):
        labs = []
        for _ in range(2 * N):
            if r % 3 != 2:
                # (near-)balanced ranks: every outcome used, ranks differ by at most one
                base = np.repeat(np.arange(d), D // d)
                extra = rng.permutation(d)[:D % d]
                lab = np.concatenate([base, extra]).astype(int)
            else:
                lab = rng.integers(0, d, size=D)               # arbitrary ranks (zero ranks allowed)
            labs.append(lab)
        scale = rng.choice([0.3, 1.0, 3.0])
        params = [torch.tensor(scale * (rng.standard_normal((D, D)) + 1j * rng.standard_normal((D, D))),
                               dtype=torch.complex128, requires_grad=True) for _ in range(2 * N)]
        opt = torch.optim.LBFGS(params, lr=1.0, max_iter=400, history_size=50, line_search_fn='strong_wolfe',
                                tolerance_grad=1e-12, tolerance_change=1e-15)

        def closure():
            opt.zero_grad()
            v = value(params, labs, N, d, C)
            v.backward()
            return v
        for _ in range(6):
            opt.step(closure)
        v = float(value(params, labs, N, d, C).detach())
        if v < best:
            best = v
            best_info = [np.bincount(l, minlength=d).tolist() for l in labs]
    return best, target, best_info, time.time() - t0


def main():
    N, d, D, restarts = (int(a) for a in sys.argv[1:5])
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    best, target, info, dt = run(N, d, D, restarts, seed)
    print(f"N={N} d={d} D={D}: best I found = {best:.12f}   I_BKP = {target:.12f}   diff = {best - target:+.3e}   "
          f"(restarts={restarts}, {dt:.0f}s; ranks of best: {info})", flush=True)


if __name__ == "__main__":
    main()
