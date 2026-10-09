"""Check S2: Theorem 1 (=>) of 03-structure.md.  For harmonic u with u_yy(z0) < 0, the pair
B = diag(x0 - delta, x0 + delta), g = y0 + eps*[[0,1],[1,0]] gives defect = eps^2 (u_yy(z0) + o(1)) < 0.
Also checks that R has mass exactly eps^2 (= ||K||^2/2) and that the o(1) is uniform in eps/delta.
50-digit arithmetic (mpmath).
"""
import mpmath as mp

mp.mp.dps = 50
tests = [
    ("Re z^3 = x^3-3xy^2", lambda x, y: x ** 3 - 3 * x * y ** 2, lambda x, y: -6 * x, (mp.mpf("0.7"), mp.mpf("0.2"))),
    ("Re e^{2z}cos", lambda x, y: mp.exp(2 * x) * mp.cos(2 * y), lambda x, y: -4 * mp.exp(2 * x) * mp.cos(2 * y),
     (mp.mpf("-0.3"), mp.mpf("0.1"))),
    ("-log|z-(2+i)|", lambda x, y: -mp.log(mp.sqrt((x - 2) ** 2 + (y - 1) ** 2)),
     lambda x, y: -((x - 2) ** 2 - (y - 1) ** 2) / ((x - 2) ** 2 + (y - 1) ** 2) ** 2, (mp.mpf("0.5"), mp.mpf("0.4"))),
]
for name, u, uyy, (x0, y0) in tests:
    target = uyy(x0, y0)
    print(f"== u = {name} at z0 = {x0}+{y0}i: u_yy(z0) = {mp.nstr(target, 12)}")
    for delta in [mp.mpf("1e-2"), mp.mpf("1e-3"), mp.mpf("1e-4")]:
        for ratio in [mp.mpf("0.1"), mp.mpf("0.5"), mp.mpf("0.99"), mp.mpf(1), mp.mpf(2), mp.mpf(10)]:
            eps = ratio * delta
            X = mp.matrix([[x0 - delta + 1j * y0, 1j * eps], [1j * eps, x0 + delta + 1j * y0]])
            ev, _ = mp.eig(X)
            evE = [x0 - delta + 1j * y0, x0 + delta + 1j * y0]
            D = sum(u(mp.re(e), mp.im(e)) for e in ev) - sum(u(mp.re(e), mp.im(e)) for e in evE)
            K2 = 2 * eps ** 2
            print(f"   delta={mp.nstr(delta,2):6s} eps/delta={mp.nstr(ratio,3):5s}: defect/eps^2 = "
                  f"{mp.nstr(D / eps ** 2, 12):>18s}   (mass of R = ||K||^2/2 = {mp.nstr(K2/2/eps**2, 3)} eps^2)")
