# Structure theory

Status of everything in this note: **CHECKED** unless marked otherwise. The proofs were re-derived, and every
statement was tested with new code, in an internal pre-publication check (not external refereeing; scripts and logs
in [`../checks/structure/`](../checks/structure/)).
Notation: B and g are Hermitian; E(g) is the pinching of g onto the eigenspaces of B; X = B + ig and X_E = B + iE(g);
K = g − E(g). Rect = [λ_min(B), λ_max(B)] × [λ_min(g), λ_max(g)] contains both spectra.

## 1. Holomorphic form of Theorem A

For s real and β ∉ spec B, let $\rho_B(s,\beta)=\frac1{2\pi}\sum_i|\operatorname{Im}\xi_i|$, with $\xi_i$ the roots of
$\det(g-s-\xi(B-\beta))$. Let R be the push-forward of $\rho_B(s,\beta)\,ds\,d\beta$ under $(s,\beta)\mapsto\beta+is$.

**Theorem A′.** For every f holomorphic on a neighbourhood of Rect,

$$\operatorname{Tr}f(B+ig)-\operatorname{Tr}f(B+iE(g))=-\int f''(\zeta)\,dR(\zeta).$$

R is a positive, absolutely continuous measure supported in Rect, with mass $\|K\|_F^2/2$. Its projections are classical
majorization defects:
- **Vertical marginal.** $\int\rho_B(s,\beta)\,d\beta=\sum_j(\lambda_j(g)-s)_+-\sum_j(\lambda_j(Eg)-s)_+$, the defect
  function of the majorization $E(g)\prec g$.
- **Horizontal marginal.** $\int\rho_B\,ds$ is the step function of Theorem A.
- **Oblique projections.** For each real ξ, the projection is the majorization defect of the Hermitian pencil g − ξB
  against its pinching E(g) − ξB.

So R is determined by classical objects, by the Cramér–Wold theorem.

*Proof.*
1. **Exponentials.** Theorem A (rescaled to any Hermitian B) at a = iκ, t = −κ gives the claim for $f(z)=e^{\kappa z}$,
   since ag − tB = κ(B + ig).
2. **Polynomials.** Differentiate n times in κ at 0.
3. **Holomorphic f.** By Runge's theorem, polynomials approximate f uniformly with second derivatives on a neighbourhood
   of the compact convex set Rect.
4. **Marginals.** They are Theorem A at t = 0 and at t = aξ. ∎

**Corollary (harmonic test functions).** For real u harmonic near Rect,
$\sum_{\operatorname{spec}X}u-\sum_{\operatorname{spec}X_E}u=\int u_{yy}\,dR$, since $u_{yy}=-\operatorname{Re}f''$ for
u = Re f.

**Theorem 1 (exact admissibility).** A harmonic u on a domain Ω satisfies
$\sum_{\operatorname{spec}X}u\ge\sum_{\operatorname{spec}X_E}u$ for all pairs with Rect ⊂ Ω if and only if $u_{yy}\ge0$ on Ω.

*Proof.*
- (⇐) The corollary.
- (⇒) Suppose $u_{yy}(x_0+iy_0)<0$.
  - Take B = diag(x₀ − δ, x₀ + δ) and $g=y_0+\varepsilon\begin{pmatrix}0&1\\1&0\end{pmatrix}$.
  - R lives in a small box around the point and has mass ε².
  - So the defect is $\varepsilon^2(u_{yy}(x_0+iy_0)+o(1))<0$. ∎

## 2. Superharmonic majorization for two levels

Let S(I) be the set of continuous functions u on the closed strip I × ℝ that are superharmonic inside, convex in y for
each x, and of at most linear growth.

**Theorem D2.** If B has exactly two eigenvalues b₁ < b₂ (any multiplicities), then for every $u\in S([b_1,b_2])$

$$\sum_{\nu\in\operatorname{spec}(B+ig)}u(\nu)\ \ge\sum_{\nu\in\operatorname{spec}(B+iE(g))}u(\nu).$$

This contains Theorem B for projections (u harmonic) and, for two levels, Ky Fan's majorization
Re spec X ≺ spec B (u = −φ(x) with φ convex; K. Fan (1950), see the survey arXiv:1108.1467).

*Proof.*
1. **The harmonic extension.** Let h be the harmonic extension to the strip [b₁, b₂] × ℝ of the edge values
   φ_k = u(b_k, ·), which are convex of linear growth.
2. **h is y-convex.** Its second differences in y are Poisson integrals of the nonnegative second differences of the
   φ_k.
3. **u ≥ h.** v = u − h is superharmonic, vanishes on the boundary, and has at most linear growth. Comparing with
   $\varepsilon\cos(\pi(x-c_0)/(2w))\cosh(\pi y/(2w))$, where $c_0=(b_1+b_2)/2$ and $w=b_2-b_1$, and using the
   minimum principle gives v ≥ 0.
4. **Theorem B for h.** By the convex-order corollary of Theorem B on each edge (rescaled to the strip; proved on
   paper, not formalized), $\sum_Xh\ge\sum_{X_E}h$.
5. **The pinched eigenvalues.** They lie on the two edge lines, where u = h.
6. **Conclusion.** $\sum_Xu\ge\sum_Xh\ge\sum_{X_E}h=\sum_{X_E}u$. ∎

The proof also gives an exact identity for two levels. For every $u\in C^2(\bar\Sigma)$ of linear growth whose
second derivatives grow at most subexponentially ($o(e^{\pi|y|/w})$), where Σ is the closed strip,
$\sum_Xu-\sum_{X_E}u=\int_{\partial}u_{yy}\,dR^\flat+\int(-\Delta u)\,T\,dA$. Here R♭ ≥ 0 is the harmonic sweep of R to
the edges, and $T(z)=\sum_{\nu\in\operatorname{spec}X}G_\Sigma(\nu,z)\ge0$, where $G_\Sigma$ is the Dirichlet Green
function of the strip, normalized by $-\Delta_zG_\Sigma(\nu,\cdot)=\delta_\nu$. The identity has the same
two-remainder structure as the spectral-defect identity of Z. Long (arXiv:2610.10158, Lemma 2.3). Long compares
spec X with Heinävaara's measure of (Re X, Im X), not with the pinched spectrum.

**Theorem D_sep (any number of levels, separable test functions).** For any Hermitian B, g, the inequality holds for
$u=a(x)+c(y)$ with c convex and $a''\le-\sup c''$ on Rect. More generally it holds for
$u=h+\tilde a(x)-\sum_id_i(\operatorname{Re}(e^{i\theta_i}z))$, with h harmonic near Rect, ã concave, $d_i$
convex with $d_i''\le\kappa_i$, and $\sum_i\kappa_i\sin^2\theta_i\le\inf_{\rm Rect}h_{yy}$. The proof is the same,
using $H_\theta=\cos\theta\,B-\sin\theta\,g$, $E(H_\theta)=\cos\theta\,B-\sin\theta\,E(g)$ and
$\|H_\theta\|_F^2-\|E(H_\theta)\|_F^2=\sin^2\theta\,\|K\|_F^2$.

*Proof idea.* Ky Fan's majorization for the real and imaginary parts, pinching majorization for c, and the trace
identity $\operatorname{Re}\operatorname{Tr}X^2=\operatorname{Tr}B^2-\|g\|_F^2$, which converts the concavity gain into
the convexity cost exactly.

*Prior work.* D_sep also follows directly from Z. Long, arXiv:2610.10158, Lemma 2.3. For f = a(x) + c(y) that
identity reads $\sum_{\operatorname{spec}X}f=\operatorname{Tr}a(B)+\operatorname{Tr}c(g)+\frac1{2\pi}\int\Psi_X(-a''-c'')$,
with $\Psi_X\ge0$ supported in the numerical range. Subtract the same identity for the normal matrix $X_E$ (where
$\Psi=0$) and use $\operatorname{Tr}c(g)\ge\operatorname{Tr}c(E(g))$. Our proof above is an elementary alternative.

**Conjecture D (open).** The inequality holds for every u ∈ S([b₁, b_m]) and any number of levels. This would generalize
Theorem B and Ky Fan's theorem at once.
- An LP search over all discretized test functions found no violation (our runs; their logs are not in this
  repository). An independent probe with 20,000 random pairs (3–6 levels) and non-separable test functions found no
  violation either ([`../checks/structure/logs/s8_conjD_probe.log`](../checks/structure/logs/s8_conjD_probe.log)).
- The first open configuration is a chain of three coupled levels.

## 3. The half-plane (accretive matrices)

**Theorem H.** Let B ≥ 0 and g be Hermitian. Let $\pi_\nu$ be the harmonic measure of {Re z > 0} on iℝ ≅ ℝ seen from
ν = x + iy, namely $\frac1\pi\frac{x}{x^2+(y-w)^2}dw$ for x > 0 and $\delta_y$ for x = 0, and
$\pi(X)=\sum_{\nu\in\operatorname{spec}X}\pi_\nu$. For every real λ,

$$\int(w-\lambda)_+\,d(\pi(X)-\pi(X_E))=V^H(\lambda):=\iint\frac1\pi\frac{\beta}{\beta^2+(s-\lambda)^2}\,\rho_B(s,\beta)\,ds\,d\beta\ \ge0,$$

with equality iff [B, g] = 0. The integral is taken against the signed measure $\pi(X)-\pi(X_E)$, whose density is
$O(|w|^{-4})$ because the masses, $\sum\operatorname{Re}\nu$ and $\sum\operatorname{Re}\nu\,\operatorname{Im}\nu$
agree; the individual integrals diverge. Hence $\int\varphi\,d(\pi(X)-\pi(X_E))\ge0$ for every convex φ of growth
$O(|w|^p)$, p < 3.

*Proof.*
- The Fourier transform of $\pi_\nu$ is $e^{i\xi y-|\xi|x}$, and its sum over spec X is $\operatorname{Tr}e^{i\xi g-|\xi|B}$.
- Theorem A at (a, t) = (iξ, |ξ|) gives $(\pi(X)-\pi(X_E))^\wedge=-\xi^2\widehat{V^H}$.
- The decay and equal low moments remove the affine ambiguity.
- $V^H$ is continuous, because $\int P_\beta(s-\lambda)\rho_B(s,\beta)\,ds\le\sup_s\rho_B(s,\beta)$, which is
  integrable in β by Lemma 1(6) of the Theorem A note. So the identity holds at every λ. ∎

Theorem H is the infinite-width limit of Theorem B. The disc version follows through the Cayley transform, but it is a
reformulation rather than a new statement.

## 4. What is false (refuted by explicit examples, checked exactly or at 60 digits)

All examples below are reproduced by
[`../checks/structure/s6_counterexamples.py`](../checks/structure/s6_counterexamples.py) (log in
[`../checks/structure/logs/`](../checks/structure/logs/)).

- **General averaging maps** do not give majorization: partial pinchings, conditional expectations onto subalgebras
  containing B, mixtures of unitary conjugations commuting with B, and unital channels commuting with B all fail.
  Maps of the form Ψ∘E, with Ψ acting inside the eigenspaces of B, do give majorization. For unital Φ with Φ(B) = B
  and E∘Φ = E (for example partial pinchings), every one-dimensional projection of R decreases, so there the failure
  is two-dimensional. Without E∘Φ = E even the projections can increase (the projection example below).
  - *Partial pinching.* B = diag(1, 1/2, 0), $E_F(g)=(g+VgV)/2$ with V = diag(1, −1, 1), and
    g = [[1/2, 3/2, 1], [3/2, −1, −1/2], [1, −1/2, 1/2]]. The left-edge stop-loss difference
    $W_{B+ig}(0.7)-W_{B+iE_F(g)}(0.7)=-0.0302673643$.
  - *Projection B.* B = diag(1, 1, 0), U = [[0, −1, 0], [1, 0, 0], [0, 0, 1]], g = [[−1, 2, −1], [2, 3, 0], [−1, 0, 0]],
    Φ(g) = (g + UgU*)/2. The left-edge difference is −0.0824667 at λ = 0.46.
- **Subharmonic test functions**, functions of y alone, and the saddle class $u_{xx}\le0\le u_{yy}$.
  - Example: B = diag(1, 0, 0) and g = [[0, 1/4, 0], [1/4, 0, 1/100], [0, 1/100, 0]]. The characteristic polynomial
    of B + ig is $\nu^3-\nu^2+\tfrac{313}{5000}\nu-\tfrac1{10000}$, with discriminant 165619531/62500000000 > 0, so
    the spectrum is entirely real, while spec E(g) = {0, ±1/100}.
  - So "Im spec(B + ig) majorizes spec(Eg)" is false.
- **One-sided comparison of logarithmic potentials.** For B = diag(1, 0) and g = [[0, 1/4], [1/4, 0]],
  $\log|\det(z-X)|-\log|\det(z-X_E)|$ is +0.030772 at z = −1 and −0.051293 at z = 1/2 + i.
- **The naive Schur–Horn converse** (every measure that dominates the pinched sweep is attained).
- **Annulus or angular orders** on the circle: vacuous.

## 5. Open problems

1. Conjecture D for three or more coupled levels.
2. Whether superharmonicity is necessary as well as sufficient.
3. The attainable set of sweeps for general diagonal blocks.
4. Monotonicity along the dephasing line g ↦ γg + (1 − γ)E(g) for two levels. It is numerically true for two levels
   (an adversarial search with 60-digit confirmation found no violation,
   [`../checks/structure/logs/s7b_dephasing_two_level.log`](../checks/structure/logs/s7b_dephasing_two_level.log)).
   It is false for three levels: for levels (0, 0.99, 1) and the g in
   [`../checks/structure/logs/three_level_candidate.json`](../checks/structure/logs/three_level_candidate.json), the
   left-edge defect W(λ) at λ = 4.7 is smaller by 0.0853 at γ = 0.998 than at γ = 0.633 (60 digits,
   [`../checks/structure/logs/s7c_three_level_confirm.log`](../checks/structure/logs/s7c_three_level_confirm.log)).
