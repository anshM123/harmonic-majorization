"""
ref_povm.py -- independent probe of the scope of Theorem 5 of math/05-chained-bell.md: minimise the BKP expression
over POVM strategies on Phi_D.
POVM parametrisation: E_a = S^{-1/2} G_a G_a^* S^{-1/2}, S = sum_a G_a G_a^* (all POVMs with invertible S).
p(a,b|k,l) = Tr(E^k_a (F^l_b)^T)/D. Prints best value vs I_N^BKP(d) and the projectivity defect of the best POVMs.
Usage: python ref_povm.py N d D restarts [seed]
"""
import sys
import time

import numpy as np
import torch

from ref_minimize import bkp_value, cost_mats

torch.set_num_threads(2)
torch.set_default_dtype(torch.float64)


def make_povm(G):
    """G: (d, D, D) complex. returns (d, D, D) POVM."""
    A = torch.einsum('aij,akj->aik', G, G.conj())          # G_a G_a^*
    S = A.sum(0)
    ev, U = torch.linalg.eigh(S)
    Sm = U @ torch.diag((ev.clamp_min(1e-300)) ** -0.5).to(torch.complex128) @ U.conj().T
    return torch.einsum('ij,ajk,kl->ail', Sm, A, Sm)


def value(Gs, N, d, C):
    D = Gs[0].shape[1]
    P = [make_povm(Gs[k]) for k in range(N)]
    Q = [make_povm(Gs[N + k]) for k in range(N)]

    def joint(Pk, Ql):
        return torch.einsum('aij,bij->ab', Pk, Ql).real / D
    c1, c2, c3 = C
    tot = 0.0
    for k in range(N):
        tot = tot + (joint(P[k], Q[k]) * c1).sum()
        if k < N - 1:
            tot = tot + (joint(P[k + 1], Q[k]) * c2).sum()
    tot = tot + (joint(P[0], Q[N - 1]) * c3).sum()
    return tot, P + Q


def run(N, d, D, restarts, seed=0):
    rng = np.random.default_rng(seed)
    C = cost_mats(d)
    target = bkp_value(N, d)
    best, best_defect = np.inf, None
    t0 = time.time()
    for r in range(restarts):
        scale = rng.choice([0.3, 1.0, 3.0])
        Gs = [torch.tensor(scale * (rng.standard_normal((d, D, D)) + 1j * rng.standard_normal((d, D, D))),
                           dtype=torch.complex128, requires_grad=True) for _ in range(2 * N)]
        opt = torch.optim.LBFGS(Gs, lr=1.0, max_iter=500, history_size=50, line_search_fn='strong_wolfe',
                                tolerance_grad=1e-13, tolerance_change=1e-16)

        def closure():
            opt.zero_grad()
            v, _ = value(Gs, N, d, C)
            v.backward()
            return v
        for _ in range(6):
            opt.step(closure)
        v, POVMs = value(Gs, N, d, C)
        v = float(v.detach())
        if v < best:
            best = v
            # projectivity defect: sum_a Tr(E_a - E_a^2)/D for each measurement
            best_defect = [float((torch.einsum('aij->', E) - torch.einsum('aij,aji->', E, E)).real.detach()) / D
                           for E in POVMs]
    return best, target, best_defect, time.time() - t0


def main():
    N, d, D, restarts = (int(a) for a in sys.argv[1:5])
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    best, target, defect, dt = run(N, d, D, restarts, seed)
    print(f"POVM N={N} d={d} D={D}: best I = {best:.12f}  I_BKP = {target:.12f}  diff = {best - target:+.3e}  "
          f"(restarts={restarts}, {dt:.0f}s; max projectivity defect of best = {max(defect):.2e})", flush=True)


if __name__ == "__main__":
    main()
