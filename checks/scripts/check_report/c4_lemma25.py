"""Lemma 2.5: ||Y(1) - e^{X_d}|| = O(1/|theta|), Y(1) = e^{i theta P} e^{X - i theta P}; explicit constant
C(theta) = ||Z||_inf (2 e^{||X||} + 2 ||X|| e^{2||X||}) with ||Z||_inf <= sum_{j!=k} ||Q_j X Q_k|| / (|theta| |p_j-p_k|).
Tested for Hermitian-type X = a g - c P and for general complex X."""
import numpy as np
from scipy.linalg import expm
from hm_common import make_case

rng = np.random.default_rng(11)
for kind in ("ag-cP", "complex X"):
    c = make_case([0.0, 0.3, 0.3, 0.31, 0.9, 1.0], rng)       # includes a small gap 0.01 between eigenvalues of P
    M = c["M"]
    if kind == "ag-cP":
        X = 0.8 * c["g"] - 0.4 * c["P"]
    else:
        X = (rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))) * 0.6
    Xd = sum(Q @ X @ Q for Q in c["Qs"])
    nX = np.linalg.norm(X, 2)
    zsum = sum(np.linalg.norm(Qj @ X @ Qk, 2) / abs(pj - pk)
               for pj, Qj in zip(c["pvals"], c["Qs"]) for pk, Qk in zip(c["pvals"], c["Qs"]) if pj != pk)
    print(f"{kind}: ||X||={nX:.3f}, min gap={np.min(np.diff(c['pvals'])):.3f}")
    for th in [10, 100, 1e3, 1e4, 1e5, -1e4]:
        Y1 = expm(1j * th * c["P"]) @ expm(X - 1j * th * c["P"])
        err = np.linalg.norm(Y1 - expm(Xd), 2)
        bound = zsum / abs(th) * (2 * np.exp(nX) + 2 * nX * np.exp(2 * nX))
        trace_err = abs(np.trace(expm(X - 1j * th * c["P"])) -
                        sum(np.exp(-1j * th * p) * np.trace(Q @ expm(Xd) @ Q) for p, Q in zip(c["pvals"], c["Qs"])))
        print(f"   theta={th:9.0f}: ||Y(1)-e^Xd||={err:.3e}  theta*err={abs(th)*err:.4f}  bound={bound:.3e}  "
              f"ok={err <= bound}  |trace err|={trace_err:.3e}")
