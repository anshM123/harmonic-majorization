# Theorems B and C: harmonic pinching for non-normal spectra

Let B be Hermitian with spectrum in [0, 1], let g be Hermitian, and let P = 1 − B. The matrix X = B + ig is in general
not normal. Its eigenvalues lie in the closed strip S = {0 ≤ Re z ≤ 1}, because $\operatorname{Re}\langle v,Xv\rangle=\langle v,Bv\rangle\in[0,1]$.

Pinching g onto the eigenspaces of B gives the commuting pair (B, g_d) and the normal matrix $X_d=B+ig_d$, whose
eigenvalues lie on the vertical lines Re z = eigenvalue of B. Theorem B compares the two spectra after sweeping them
to an edge of the strip by harmonic measure. Theorem C makes the comparison quantitative.

Notation is as in [`01-theorem-A.md`](01-theorem-A.md): $p_1<\dots<p_r$ are the eigenvalues of P, the $Q_k$ are its
spectral projections, $g_d=\sum_kQ_kgQ_k$, and ρ is the density of Theorem A for the pair (g, P).

## 1. Strip harmonic measure

For 0 < x < 1 let

$$K_x(u)=\frac{\sin\pi x}{2(\cosh\pi u-\cos\pi x)} .$$

This is the Poisson kernel of the strip for its left edge: a point ν = x + iy in the open strip sees the left edge
{Re z = 0} with harmonic measure $\omega_\nu(dw)=K_x(y-w)\,dw$, where the edge is identified with ℝ through the
imaginary part. On the edges, $\omega_\nu=\delta_y$ if x = 0 and $\omega_\nu=0$ if x = 1.

**Lemma 1 (strip kernel)** (proved in our OQP 27 paper, github.com/anshM123/IQOQI-OQP-27, `papers/math/main.tex`,
lemmas with LaTeX labels `lem:poisson` and `lem:h`).
1. $K_x>0$ and $\int K_x=1-x$.
2. For $|\operatorname{Re}a|<\pi$, $\int e^{-au}K_x(u)\,du=\sin(a(1-x))/\sin a$.
3. For ν = x + iy in the closed strip, $\int e^{aw}\,\omega_\nu(dw)=e^{ay}\sin(a(1-x))/\sin a$.
4. $K_x$ is even, decreasing in |u|, and decays like $e^{-\pi|u|}$.

The **left sweep** of a finite multiset Z in the closed strip is the measure $\omega^L(Z)=\sum_{\nu\in Z}\omega_\nu$ on ℝ.
The **right sweep** is $\omega^R(Z)=\omega^L(\{1-\bar\nu:\nu\in Z\})$. Eigenvalues are counted with algebraic
multiplicity.

## 2. Statements

**Theorem B (harmonic pinching).**
- **Left edge.** Let $\omega=\omega^L(\operatorname{spec}(B+ig))$ and $\sigma=\omega^L(\operatorname{spec}(B+ig_d))$.
  - The two measures have the same mass, M − Tr B, and the same first moment, Tr(Pg).
  - For every real λ,

    $$\int(w-\lambda)_+\,\omega(dw)-\int(w-\lambda)_+\,\sigma(dw)=V(\lambda):=\iint K_{1-\tau}(s-\lambda)\,\rho(s,\tau)\,ds\,d\tau\ \ge0 .$$

  - Hence $\sigma\le_{\rm cx}\omega$: $\int f\,d\sigma\le\int f\,d\omega$ for every convex f of at most linear growth.
- **Right edge.** The same holds for the right sweeps, with ρ replaced by the density of the pair (g, B), that is, of
  det(g − s − ξ(B − τ)).
- **Equality.** Equality holds for one λ iff it holds for all λ, iff [B, g] = 0.

In words, the eigenvalues of the non-normal matrix B + ig, swept to either edge of the strip, are more spread out in
convex order than those of its pinched, commuting version. The excess is an explicit nonnegative integral, and it
vanishes exactly in the commuting case.

**Theorem C (quantitative form).** For every real λ put R = |λ| + ‖g‖. Then

$$V(\lambda)\ \ge\ \frac1{2(1+\cosh\pi R)}\sum_{j<k}\|Q_jgQ_k\|_F^2\;\frac1{p_k-p_j}\int_{p_j}^{p_k}\sin(\pi\tau)\,d\tau .$$

For a projection B this reads

$$\sum_{\nu\in\operatorname{spec}(B+ig)}h_\lambda(\nu)-\operatorname{Tr}\bigl[(P(g-\lambda)P)_+\bigr]\ \ge\ \frac{\|[B,g]\|_F^2}{2\pi\,(1+\cosh(\pi(|\lambda|+\|g\|)))},$$

where $h_\lambda(\nu)=\int(w-\lambda)_+\,\omega_\nu(dw)$.

## 3. Proof of Theorem B

Let $\nu_k=x_k+iy_k$ be the eigenvalues of X = B + ig, with $x_k\in[0,1]$.

**(a) Laplace transforms.** Let a be real with 0 < |a| < π.
- By Lemma 1(3), $\hat\omega(a):=\int e^{aw}\,\omega(dw)=\sum_ke^{ay_k}\sin(a(1-x_k))/\sin a$.
- This equals $(\operatorname{Tr}e^{ag+iaP}-\operatorname{Tr}e^{ag-iaP})/(2i\sin a)$, because:
  - $e^{-ia(X-1)}=e^{ag+iaP}$ has eigenvalues $e^{ay_k}e^{ia(1-x_k)}$;
  - $e^{ia(X^*-1)}=e^{ag-iaP}$ has eigenvalues $e^{ay_k}e^{-ia(1-x_k)}$.
- The same computation for $X_d$ gives $\hat\sigma(a)=(\operatorname{Tr}e^{ag_d+iaP}-\operatorname{Tr}e^{ag_d-iaP})/(2i\sin a)$.
- Therefore $\mathcal D(a,-ia)-\mathcal D(a,ia)=2i\sin a\,(\hat\omega(a)-\hat\sigma(a))$, where $\mathcal D$ is the
  difference of Theorem A.
- By Theorem A at t = ∓ia, the left side equals $2ia^2\iint e^{as}\sin(a\tau)\rho$.
- For τ ∈ (0, 1), $\sin(a\tau)/\sin a=\int e^{-au}K_{1-\tau}(u)\,du$ (Lemma 1(2) with x = 1 − τ). The values τ ∈ {0, 1}
  form a null set for ρ. By Tonelli,

  $$\hat\omega(a)-\hat\sigma(a)=a^2\int e^{aw}V(w)\,dw .$$

**(b) Moments.**
- Since $\hat\omega-\hat\sigma=O(a^2)$, ω and σ have equal mass and equal first moment.
- Explicitly, the mass is $\sum(1-x_k)=M-\operatorname{Tr}B$, the same for $X_d$. The first moment is
  $\sum y_k(1-x_k)=\operatorname{Tr}g-\operatorname{Tr}(Bg)=\operatorname{Tr}(Pg)$, using
  $\sum x_ky_k=\tfrac12\operatorname{Im}\operatorname{Tr}X^2=\operatorname{Tr}(Bg)$. The same holds for $g_d$.
- Put $W(w)=\int(w-u)_+\,d(\omega-\sigma)(u)$.
  - Both measures are finitely many atoms plus densities that are $O(e^{-\pi|w|})$, so W and W′ decay exponentially and
    $W''=\omega-\sigma$.
  - Integrating by parts twice, $\int e^{aw}W=a^{-2}(\hat\omega-\hat\sigma)=\int e^{aw}V$. Both sides are holomorphic in
    $|\operatorname{Re}a|<\pi$, so W = V almost everywhere.
- **V is continuous.**
  - For τ ∈ (0, 1) ∖ spec P, $\int K_{1-\tau}(s-w)\rho(s,\tau)\,ds\le\tau\,\Phi(\tau)$, where $\Phi=\sup_s\rho(\cdot,\tau)$
    is integrable by the gap bound of Theorem A, Lemma 1(6). The bound does not depend on w.
  - So V is continuous by dominated convergence. Since W is continuous too, **W = V everywhere**.

**(c) Conclusion.** Write $(u-\lambda)_+=(\lambda-u)_++(u-\lambda)$. Since ω − σ has zero mass and zero first moment,
$\int(u-\lambda)_+\,d(\omega-\sigma)=W(\lambda)=V(\lambda)\ge0$. Since the call functions $(w-\lambda)_+$ together with
affine functions generate the convex functions of linear growth, $\sigma\le_{\rm cx}\omega$.

**(d) Equality.**
- If V(λ) = 0: $K_{1-\tau}>0$ for τ ∈ (0, 1), and ρ ≥ 0 is continuous off the lines τ ∈ spec P, so ρ = 0 almost
  everywhere. Then the τ-marginal of Theorem A vanishes, so $Q_jgQ_k=0$ for all j ≠ k, that is, [B, g] = 0.
- Conversely, if [B, g] = 0, then $g=g_d$ and ρ = 0.

**Right edge.** Apply the left-edge statement to (P, g). P has spectrum in [0, 1] and the same eigenspaces as B, hence the
same pinching. Also $P+ig=1-(B+ig)^*$, whose eigenvalues are $1-\bar\nu$. The left sweep of $\{1-\bar\nu\}$ is the right
sweep of {ν}. ∎

## 4. Proof of Theorem C

- For s in $[\lambda_{\min}(g),\lambda_{\max}(g)]$ we have |s − λ| ≤ R. Since $K_{1-\tau}$ is even and decreasing in |u|,

  $$K_{1-\tau}(s-\lambda)\ge K_{1-\tau}(R)=\frac{\sin\pi\tau}{2(\cosh\pi R+\cos\pi\tau)}\ge\frac{\sin\pi\tau}{2(1+\cosh\pi R)} .$$

- By Theorem A, $\int\rho(s,\tau)\,ds=\mu(\tau)$, the explicit step function. Hence
  $V(\lambda)\ge\frac1{2(1+\cosh\pi R)}\int\sin(\pi\tau)\,\mu(\tau)\,d\tau$, which is the stated bound.
- For a projection B, $\mu=\|BgP\|_F^2$ on (0, 1), $\int_0^1\sin\pi\tau\,d\tau=2/\pi$, and $\|[B,g]\|_F^2=2\|BgP\|_F^2$. ∎

**Remarks.**
- **Sharpness (numerical).** The constant appears to be sharp: in an adversarial search the ratio of V(λ) to the bound
  came down to 1.0010, checked at 50 digits. The log of that run was not kept, so treat this as an observation.
- **The λ-dependence cannot be dropped.** ρ has compact support and $K_{1-\tau}(u)$ decays like $e^{-\pi|u|}$, so V(λ) → 0
  exponentially as |λ| → ∞, while the sum on the right side of Theorem C does not depend on λ.

## 5. Formal verification

The defect identity of Theorem B on both edges, the inequality V(λ) ≥ 0, the equality case, and Theorem C are
machine-checked in Lean 4 with Mathlib: `theoremB`, `theoremB_right`, `theoremB_eq_iff`, `theoremC`, `theoremC_proj`
in [`../lean/`](../lean/). They use only the standard axioms. The convex-order corollary for all convex test functions
is proved here on paper but not formalized.

## 6. Credits

- For projections B, Theorems A and B are Theorems D and C of our OQP 27 paper.
- Theorem B for general B and Theorem C were not found in the literature (searches on 2026-10-08, including Heinävaara's
  work and Z. Long, arXiv:2610.10158).
- The positivity input is Stahl's theorem in its explicit two-variable form (Stahl; Heinävaara).
- Strip kernels and pinching appear in multivariate Golden–Thompson inequalities (Sutter–Berta–Tomamichel;
  Hiai–König–Tomamichel), but not for spectra of non-normal matrices.
