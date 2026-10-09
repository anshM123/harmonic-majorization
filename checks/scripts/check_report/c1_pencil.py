import numpy as np
from scipy.integrate import quad
from hm_common import *

rng = np.random.default_rng(1)

# (a) closed-form h_lambda vs direct quadrature of (w-lam)_+ K_x(y-w)
worst = 0
for _ in range(30):
    x, y, lam = rng.uniform(0.02, 0.98), rng.normal() * 2, rng.normal() * 2
    num = quad(lambda w: (w - lam) * Kx(x, y - w), lam, np.inf, limit=400, epsabs=1e-13)[0]
    worst = max(worst, abs(num - h_lambda(complex(x, y), lam)))
print("(a) max |h_closed - h_quad| =", worst)

# (b) refined bound: non-real root xi at (s,tau), p_k<tau<p_{k+1}:  |xi| sqrt((tau-p_k)(p_{k+1}-tau)) <= ||g - s||
worst_ratio = 0
viol = 0
for trial in range(400):
    M = int(rng.integers(3, 8))
    lv = rng.choice([0.0, 1.0, rng.uniform(), rng.uniform(), rng.uniform()], size=M)
    c = make_case(lv, rng, gscale=rng.uniform(0.3, 3))
    pv = c["pvals"]
    if len(pv) < 2:
        continue
    for _ in range(60):
        k = int(rng.integers(0, len(pv) - 1))
        # sample tau very close to an endpoint sometimes
        u = rng.uniform() if rng.uniform() < 0.5 else 10 ** rng.uniform(-8, -1)
        if rng.uniform() < 0.5:
            u = 1 - u
        tau = pv[k] + u * (pv[k + 1] - pv[k])
        if not (pv[k] < tau < pv[k + 1]):
            continue
        s = rng.uniform(c["lmin"], c["lmax"])
        xi = roots(c, s, tau)
        nr = xi[np.abs(xi.imag) > 1e-9 * (1 + np.abs(xi))]
        if len(nr) == 0:
            continue
        bound = np.linalg.norm(c["g"] - s * np.eye(M), 2) / np.sqrt((tau - pv[k]) * (pv[k + 1] - tau))
        r = np.max(np.abs(nr)) / bound
        worst_ratio = max(worst_ratio, r)
        if r > 1 + 1e-6:
            viol += 1
print(f"(b) refined pencil bound: max |xi|/bound over non-real roots = {worst_ratio:.6f}, violations = {viol}")

# (c) growth of sup_s rho(s,tau) as tau -> interior p_k, and as tau -> 1 when 0 in spec B
c = make_case([0.0, 0.0, 0.35, 0.35, 0.7, 1.0], rng)
pv = c["pvals"]
print("(c) pvals =", pv)
s = np.linspace(c["lmin"], c["lmax"], 40001)
for k in range(1, len(pv)):
    for side in (-1, +1):
        if (k == len(pv) - 1 and side == +1):
            continue
        out = []
        for eps in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
            tau = pv[k] + side * eps
            r = np.array([rho(c, si, tau) for si in s])
            imax = np.argmax(r)
            # refine around the max
            s2 = np.linspace(s[max(imax - 2, 0)], s[min(imax + 2, len(s) - 1)], 2001)
            r2 = max(rho(c, si, tau) for si in s2)
            out.append((eps, r2, r2 * np.sqrt(eps), np.trapezoid(r, s), mfun(c, tau)))
        print(f"   near p_{k}={pv[k]:.3f} side {side:+d}:")
        for eps, r2, sc, mint, mt in out:
            print(f"      eps={eps:.0e}  sup rho={r2:10.4f}  sup*sqrt(eps)={sc:.4f}  int rho ds={mint:.6f}  m={mt:.6f}")
