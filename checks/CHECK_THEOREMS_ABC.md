> **How to read this report (added for publication).**
> - This is the internal independent check of the working notes THEOREMS.md, which are published, with the gap fixed,
>   as [`../math/01-theorem-A.md`](../math/01-theorem-A.md) and
>   [`../math/02-harmonic-pinching.md`](../math/02-harmonic-pinching.md). The lemma numbering differs; see
>   [`../lean/README.md`](../lean/README.md).
> - The gap it found (continuity of V in Theorem B) was closed by the bound now called Lemma 1(6) of the Theorem A
>   note, before publication.
> - None of the numerical runs quoted in this report has a kept output log. No claim in this repository rests on them
>   alone. `adversary.py`, `pinching_test.py`, `bmv_general.py` and `check_theoremC.py` are internal development
>   scripts and are not included.
> - Theorems A, B (defect identity and equality case) and C have since been machine-checked in Lean; see
>   [`../lean/`](../lean/).

# Check report on THEOREMS.md (Theorems A, B, C: harmonic pinching)

Checked 2026-10-08, against `THEOREMS.md` as it stood on that date. The comparison text, cited as [27], is
the OQP 27 math paper (github.com/anshM123/IQOQI-OQP-27, `papers/math/main.tex`):
- Section "Two-variable BMV positivity and the strip inequality" (lines 575–1005);
- Appendix "A distribution-free proof of Theorem D" (lines 2582–2900), including Remark RIgeneralB.

Method:
- every step of the general case was re-derived and compared line by line with the projection-case proof;
- every claim that a step is "verbatim" was checked against what the projection proof actually uses;
- numerical tests (Section 5) were run on general B: interior eigenvalues, eigenvalues 0 and 1, repeated eigenvalues,
  small spectral gaps.

The scripts are in `scripts/check_report/` (in this folder).

## 1. Summary

**No ERROR was found.** All three theorems are numerically consistent, and their proofs are correct apart from one gap.

| Item | Verdict |
|---|---|
| Theorem A (pinched BMV identity) | VERIFIED. The proof is complete as written. |
| Theorem B (harmonic pinching) | One GAP, in step (b). The fix is short (Section 3.2). |
| Theorem C (quantitative bound) | VERIFIED. The constant is asymptotically sharp. |

### GAPs and other findings, most severe first

1. **GAP, Theorem B (b)–(c): the pointwise identity W(lambda) = V(lambda).**
   - Fourier uniqueness only gives W = V almost everywhere.
   - To get the displayed identity of Theorem B at *every* lambda, V must be finite and continuous. The notes never
     prove this.
   - In the projection case [27, proof of Thm C (a)] proves it, from the bound sup_s rho(s,tau) <= C (tau(1-tau))^(-1/2)
     ([27, Lemma pencil(3)]). The general notes drop this bound and put nothing in its place.
   - The gap matters exactly when the kernel K_{1-tau} narrows: when 0 is in spec B (p_m = 1) for the left edge, and
     when 1 is in spec B (p_1 = 0) for the right edge.
   - **Fix:** a general version of [27, Lemma pencil(3)]. It is proved in Section 2 (Lemma 2.2(3')), applied in
     Section 3.2 and checked numerically.
   - Even without the fix, everything except the pointwise identity survives:
     - the inequality W >= 0;
     - convex order;
     - the equality case;
     - Theorem C;
     - the projection reading of Theorem C.

     The reason is that V is lower semicontinuous, so V <= W everywhere (Section 3.2).
2. **Minor (wording), Theorem A, Step 3.** The parenthetical reason for the continuity of tau -> int rho(s,tau) ds
   ("since rho is locally bounded on Omega") is incomplete. The correct reason is that rho is continuous on Omega and
   vanishes for s outside [l_min, l_max]: uniform continuity on compact subsets of Omega, or dominated convergence. Both
   facts are already stated in the notes, so nothing is missing mathematically.
3. **Attribution (not mathematical).** [27, Remark RIgeneralB] (main.tex lines 2878–2900) already states, for arbitrary
   Hermitian A and B:
   - the identity of Theorem A;
   - Corollary 2 (Stahl's theorem with the pinched atoms).

   The proof there is deferred to the accompanying notes. THEOREMS.md cites this remark only in Section 8 (Lean plan).
   - Section 0 ("these notes generalise Theorems C and D") and Section 7 should say that the general form of D was
     already announced in [27].
   - Relative to [27], Theorem A adds the tau-marginal formula m(tau) and the distributional proof through Lemmas 2.5
     and 2.6.
4. **Cosmetic, Lemma 2.5.** delta, "the smallest gap between eigenvalues of P", is undefined when m = 1. The lemma is
   then trivial, since X_o = 0 and Z = 0.
5. **Framing, Theorem C versus [27, open question (3)].** The bound has the weight 1/(2(1 + cosh(pi(|lambda| +
   ||g||)))).
   - This weight is unavoidable in this form: the bound is asymptotically sharp (Section 3.3).
   - In the intended application ([27, Remark deficit]) the defects are integrated over theta, and g(theta) is the
     conjugate function of a step field, which has logarithmic singularities at the jump points.
   - So Theorem C gives a *weighted* lower bound there. "This is the lower bound ... that [27, Section 'open questions'
     (3)] asks for" is accurate only if this weight is stated. Whether it is enough for quantitative rigidity is not yet
     shown.
6. **Numerical scripts (not the mathematics).**
   - (a) Evaluating W(lambda) in double precision as sum h_lambda(nu) - sum h_lambda(mu) cancels catastrophically when
     lambda is far below the spectrum: relative error up to 1e-3 at W ~ 1e-12. In an adversarial search for Theorem C,
     this produced a spurious ratio of 0.99925 < 1. In 50-digit arithmetic the true value is 1.003955.
   - **Remedy:** for lambda below the centre of the spectrum, evaluate the mirrored stop-loss
     sum h_{-lambda}(conj nu) - sum h_{-lambda}(conj mu). This is valid because the two measures have equal masses and
     equal first moments.
   - (b) `adversary.py`, `stop_point`, builds a grid with 80/ds points, ds = x/25. For an eigenvalue with small
     Re = x > 1e-9 the memory explodes: x = 1e-5 already gives about 2e8 points per array.
   - (c) `check_theoremC.py` exists but is not listed in the Section 6 table.
   - (d) Several scripts build the pinching mask with `np.isclose`, whose default relative tolerance is 1e-5:
     `adversary.py`, `pinching_test.py`, `bmv_general.py`, `check_theoremC.py`. Nearly equal but distinct levels of B
     are then pinched together, which is a coarser pinching than the one in the theorems.
     - For the scripts that draw levels from a finite set with exact repeats this changes nothing.
     - In `adversary.py` the levels are continuous, so the objective can use the wrong pinching. This lowers W, so its
       non-negativity findings stay valid. That W_coarse <= W_exact was checked on 300 random cases with a level gap
       between 1e-9 and 1e-4.
     - Combined with a bound computed from the exact eigenvalue partition, it produces false violations. This check's
       own first harness reported a Theorem C ratio of 0.019, which recomputation with exact pinching turned into
       1.424.
     - Use rtol = 0, for example `np.abs(p[:, None] - p[None, :]) <= 1e-12`.

## 2. Theorem A, step by step

### Setting (Section 0)

**VERIFIED.**
- The spectrum of B + ig lies in the closed strip, because the numerical range contains the spectrum.
- g_d is block diagonal with respect to the eigenspaces of P, and Tr e^{a g_d - tP} equals the subtracted sum.
- The total mass is sum_{j<k} ||Q_j g Q_k||_F^2 = ||g - g_d||_F^2 / 2.

### Lemma 2.1

**VERIFIED.** It is verbatim [27, Lemma poisson (1),(3)] and [27, Lemma h(5)]. No projection-specific input is used.

### Lemma 2.2

**VERIFIED.**
- (1) g - s - xi(P - tau) = (P - tau)((P - tau)^{-1}(g - s) - xi), so c_tau = det(tau - P) = prod_k (tau - p_k)^{r_k}.
  The root sums agree because (P - tau)^{-1} = sum_k (p_k - tau)^{-1} Q_k and Tr(Q_k g) = Tr(Q_k g_d).
- (2) For q, ||g_d - s|| <= ||g - s||, because the pinching is a norm contraction.
- (3) The argument is correct for each of the four cases.
- (4) g_d - s - xi(P - tau) is a direct sum of r_k x r_k blocks, each of exact degree r_k in xi, since p_k != tau.
- (5) Continuity of the roots holds.
- The consequences (continuity, local boundedness, support) follow.

**Proposed addition (3').** This is the general form of [27, Lemma pencil(3)], and it is what closes the gap in
Theorem B.

> Let p_k < tau < p_{k+1} be consecutive eigenvalues of P. Every non-real root xi of p_{s,tau} satisfies
> |xi| <= ||g - s|| / sqrt((tau - p_k)(p_{k+1} - tau)). Hence, for such tau,
>
>     rho(s,tau) <= (M / 2pi) ||g - s|| / sqrt((tau - p_k)(p_{k+1} - tau)),
>
> and int_R sup_{s} rho(s,tau) dtau <= (M (m-1) / 2) (l_max - l_min) < infinity.

*Proof.*
1. Let v be a unit vector with (g - s)v = xi (P - tau)v. As in Lemma 2.2(3), Im xi != 0 forces <v,(P - tau)v> = 0.
2. Put w_j = ||Q_j v||^2 and d_j = p_j - tau. Then sum_j w_j = 1 and sum_j d_j w_j = 0.
3. No p_j lies strictly between p_k and p_{k+1}, so (d_j - d_k)(d_j - d_{k+1}) >= 0 for every j. That is,
   d_j^2 >= (d_k + d_{k+1}) d_j - d_k d_{k+1}.
4. Multiply by w_j and sum: ||(P - tau)v||^2 = sum_j d_j^2 w_j >= -d_k d_{k+1} = (tau - p_k)(p_{k+1} - tau).
5. Finally |xi| ||(P - tau)v|| = ||(g - s)v|| <= ||g - s||.
6. The integral bound follows because int_{p_k}^{p_{k+1}} dtau / sqrt((tau - p_k)(p_{k+1} - tau)) = pi, and
   ||g - s|| <= l_max - l_min for s in [l_min, l_max]. QED

For spec P = {0,1} this is exactly [27, Lemma pencil(3)].

Two side benefits:
- (3') gives rho in L^1 directly, so the "Integrability" part of Step 3 becomes a corollary rather than a necessity.
- (3') shows the rate at which rho blows up at the lines tau = p_k.

Numerics:
- (3') held in 400 random configurations (M = 3..7, eigenvalues 0 and 1, interior eigenvalues, tau down to 1e-8 from a
  line). The largest ratio |xi| / bound was 0.977.
- sup_s rho(s,tau) really does blow up like c |tau - p_k|^(-1/2), on one side of each interior line and at tau -> 1
  (c ~ 0.38, 0.31, 0.66 in the test case).
- int rho ds = m(tau) is reproduced even at |tau - p_k| = 1e-6.
- So rho is genuinely unbounded near the lines. The notes are right to work only with local boundedness on Omega.

### Lemma 2.3

**VERIFIED.** It is verbatim [27, Lemma jensen]. The hypotheses are exactly Lemma 2.2(1),(4).

### Lemma 2.4

**VERIFIED.**
- (1) and (3) are verbatim.
- In (2), the eigenvalues of g_d - xi P are those of Q_k g Q_k on ran Q_k, shifted by -xi p_k. So
  sum_j e^{a lambda^0_j} = sum_k e^{-a xi p_k} Tr_{Q_k} e^{a Q_k g Q_k}, and the integration by parts is unchanged.

### Lemma 2.5 (high-frequency averaging)

**VERIFIED.** The O(1/|theta|) bound is correct, and so are the integration-by-parts terms.
- **Interaction picture.** Y' = e^{i theta r P} X e^{-i theta r P} Y = X(r) Y. The expansion of X(r) into X_d + X_o(r)
  is correct.
- **Duhamel formula.** It follows from F(u) = e^{(1-u)X_d} Y(u), F' = e^{(1-u)X_d} X_o(u) Y(u), since X_d commutes with
  e^{(1-u)X_d}.
- **The key uniform bound.** ||Y(u)|| <= e^{||X||}, uniformly in theta. It holds because e^{i theta u P} is unitary and
  the Hermitian part of X - i theta P equals that of X, since P is Hermitian and theta is real.
- **The other bounds.** ||X(u)|| = ||X|| (unitary conjugation). ||X_d|| <= ||X|| (pinching).
- **Integration by parts.** It gives three terms: the boundary term [e^{(1-u)X_d} Z Y]_0^1, the term with
  (e^{(1-u)X_d})' Z Y, and the term with e^{(1-u)X_d} Z Y'. Altogether

      ||Y(1) - e^{X_d}|| <= ||Z||_inf (2 e^{||X||} + 2 ||X|| e^{2||X||}),
      ||Z||_inf <= sum_{j != k} ||Q_j X Q_k|| / (|theta| |p_j - p_k|) <= m^2 ||X|| / (|theta| delta).

- **Conclusion.** Tr(e^{-i theta P} e^{X_d}) = sum_k e^{-i theta p_k} Tr_{Q_k} e^{Q_k X Q_k}.
- **The application.** For X = a g - c P, Q_k X Q_k = a Q_k g Q_k - c p_k Q_k, and
  e^{-i theta p_k} Tr_{Q_k} e^{Q_k X Q_k} = e^{-(c + i theta) p_k} Tr_{Q_k} e^{a Q_k g Q_k}.
- **Numerics.**
  - theta * ||Y(1) - e^{X_d}|| stays bounded for theta from 10 to 1e5 (and theta = -1e4).
  - The explicit constant above was never violated.
  - Both held for X = 0.8 g - 0.4 P and for a random non-Hermitian X, with a spectral gap of 0.01.
- **Cosmetic:** delta is undefined for m = 1 (trivial case).

### Lemma 2.6

**VERIFIED.**
- T(theta) = f(theta)/theta^n - O(1/theta) -> 0.
- The mean of |T|^2 over [L, 2L] tends to sum |c_k|^2, because the cross terms average to O(1/(L |p_k - p_l|)).

### Step 0

**VERIFIED.**
- D(0,t) = 0 and d_a D(0,t) = sum_k e^{-t p_k}(Tr Q_k g - Tr Q_k g Q_k) = 0.
- The Taylor remainder formula holds.

### Step 1 (Paley–Wiener–Schwartz with K = [l_min, l_max] x [p_1, p_m])

**VERIFIED.**
- (a,t) = (-i zeta_1, i zeta_2) gives Re a = eta_1 and Re t = -eta_2, so the Hermitian part of X is eta_1 g + eta_2 P.
- The joint numerical range of (g, P) lies in K, so |Tr e^X| <= M e^{H_K(eta)}.
- Each pinched term is a sum of r_k terms e^{eta_1 mu + eta_2 p_k}, with (mu, p_k) in K. So |D| <= 2M e^{H_K(eta)}.
- **The disc argument for |zeta_1| < 1** carries over unchanged.
  - K is a product of intervals, so H_K(eta) = max(l_min eta_1, l_max eta_1) + max(p_1 eta_2, p_m eta_2).
  - This is Lipschitz in eta_1 with constant c_K = max(|l_min|, |l_max|).
  - On |z| = 2 we have |Im z - eta_1| <= 3, which gives |G^| <= (M/2) e^{H_K(eta) + 3 c_K}.
- PWS is applied with N = 0. A degenerate K (l_min = l_max, or m = 1) is allowed, because PWS holds for any convex
  compact set, and then D = 0 anyway.
- On R^2, |G^| <= C.

### Step 2 (Theta = rho on Omega)

**VERIFIED.** The claim "verbatim" is correct.
- [27, Step 2] uses only:
  - the bounded Fourier transform;
  - the slice Lemma 2.4;
  - Lemma pencil (1), (2), (5);
  - Lemma jensen.

  Their general forms are Lemmas 2.4, 2.2(1), 2.2(2), 2.2(4) and 2.3.
- For supp phi compact in Omega, put eps = dist(supp phi, R x spec P) > 0. The roots are then bounded by
  (||g|| + S_0)/eps.
- The estimate on [-2R, 2R] and the O(xi^{-2}) tail estimate are unchanged. So the integrand is absolutely integrable,
  and Fubini and Jensen give Theta(phi) = int phi rho.

### Step 3 (tau-marginal; integrability)

**VERIFIED** (with the wording point of finding 2).
- **Duhamel.** The computation was re-derived. d_a^2 Tr e^{ag - tP} at a = 0 equals
  sum_{j,k} ||Q_j g Q_k||_F^2 int_0^1 e^{-t(theta p_j + (1-theta) p_k)} dtheta. The terms j = k cancel against the
  pinched sum. For j != k, substituting tau = theta p_j + (1-theta) p_k gives
  (1/|p_j - p_k|) int_{min}^{max} e^{-t tau} dtau, and the pairs (j,k) and (k,j) are equal.
- **The factor 1/2 is correct.** G(0,t) = (1/2) d_a^2 D(0,t) = int e^{-t tau} m(tau) dtau.
- **The marginal is m.** Theta(chi (x) e^{-t tau}) = G(0,t), because chi = 1 near the s-projection of supp Theta.
  Fourier–Laplace transforms are injective on E'(R), so Theta_tau = m dtau.
- **Slices.** chi (x) psi lies in C_c^infty(Omega) when psi is in C_c^infty(R \ spec P), so Step 2 applies. tau -> int rho
  ds is continuous off spec P, and m is constant on each gap, so the slice identity holds at every tau not in spec P.
- **Integrability.** Tonelli gives int int rho = int m = sum_{j<k} ||Q_j g Q_k||_F^2.
- **No circularity.** Step 3 uses only Steps 1 and 2.
- **Numerics.** int rho(s,tau) ds = m(tau):
  - to about 1e-6 on every gap, in a 6 x 6 case with spec P = {0, 0.3, 0.65, 1};
  - to <= 5e-4 relative even at |tau - p_k| = 1e-6, where sup rho ~ 650 (limited by the grid).

### Step 4 (no singular part on the lines)

**VERIFIED.**
- **Support of T.** T = Theta - rho is a distribution, because rho is in L^1 with compact support (Step 3). It vanishes
  on Omega (Step 2) and outside K. So its support lies in the disjoint segments L_k.
- **Structure.** Cutoffs and the structure theorem ([Hormander, Thm 2.3.5]; compact support implies finite order) give
  the stated form.
- **Evaluation.** <delta^(j)(tau - p_k), e^{-t tau}> = t^j e^{-t p_k}.
- **The left side tends to 0** along t = c + i theta:
  - by Lemma 2.5 for G;
  - by Riemann–Lebesgue for the rho term, applied to tau -> e^{-c tau} int e^{as} rho(s,tau) ds, which is in L^1.
- **Lemma 2.6** shows that each pi_k(theta) = e^{-c p_k} sum_j (c + i theta)^j A_{k,j}(a) is the zero polynomial. A
  polynomial in t that vanishes on the line c + iR vanishes identically, so A_{k,j}(a) = 0 for real a != 0.
- **Conclusion.** The A_{k,j} are entire, so alpha_{k,j} = 0.
- **Remark.** The real-t limits used in [27, Step 3] (Lemma limits) cannot separate interior lines, because the
  exponentials e^{-t p_k} are dominated by those of the extreme lines. The imaginary-direction argument is the right
  replacement, and it is correct.

### Numerical confirmation of Theorem A (this check)

- **Line form.** Theorem A for real (a, t = a xi), via Laplace inversion in w, is equivalent to the Radon-slice
  identity int rho(w + xi tau, tau) dtau = U_xi(w). It was tested on 24 random lines in 4 configurations, with
  spec B = {0, 0.4, 1}, {0.2, 0.2, 0.7, 0.9}, {0,0,0.35,0.35,0.7,1}, {0.1,0.3,0.5,0.5,0.95}. The largest difference
  was 8.4e-11.
- **The 2-D identity at complex (a,t).** It was tested at 18 points, for example (0.7 + 0.4i, 2 - 3i), (1.3, 6i) and
  (-0.5i, 1 + i), in 3 configurations. The relative difference was <= 9e-6 (quadrature-limited).
- **Total mass.** It matched sum ||Q_j g Q_k||^2 to 8e-7 relative.

## 3. Theorems B and C, step by step

### 3.1 Theorem B (a): Laplace transforms

**VERIFIED.**
- **Signs.** X - 1 = -P + ig, so -ia(X - 1) = ag + iaP, with eigenvalues e^{a y_k} e^{ia(1 - x_k)}.
  X^* - 1 = -P - ig, so ia(X^* - 1) = ag - iaP, with eigenvalues e^{a y_k} e^{-ia(1 - x_k)}. Both hold with algebraic
  multiplicity.
- **The pinched sweep.** B + i g_d has the eigenvalues (1 - p_k) + i mu on the line Re = 1 - p_k. Its left sweep has
  Laplace transform sum_k Tr_{Q_k} e^{a Q_k g Q_k} sin(a p_k)/sin a, which is
  (Tr e^{a g_d + iaP} - Tr e^{a g_d - iaP}) / (2i sin a).
- **The difference.** D(a,-ia) - D(a,ia) = 2i sin a (omega^ - sigma^), and Theorem A gives
  2i a^2 int int e^{as} sin(a tau) rho.
- **The kernel.** Lemma 2.1 is used with x = 1 - tau, for tau in (0,1). supp rho is contained in [p_1, p_m], which lies
  in [0,1], so tau in {0,1} is a null set.
- **The substitution** w = s - u, and Tonelli (for real a, every factor is >= 0).

### 3.2 Theorem B (b)–(d): moments, W = V, equality

**Moments: VERIFIED.**
- The masses are sum(1 - x_k) = M - Tr B.
- The first moment of omega_nu is y(1 - x), since sin(a(1-x))/sin a is even in a.
- sum x_k y_k = (1/2) Im Tr X^2 = Tr(Bg), and Tr(B g_d) = Tr(Bg).
- Numerically, the masses and first moments of the two sweeps agree with each other, and with M - Tr B, Tr(Pg) (left)
  and Tr B, Tr(Bg) (right), to 12 digits. This was checked in 4 configurations (`c5b_thmB.py ... moments`).

**Decay of W, W'; W = V almost everywhere: VERIFIED.**
- Both measures are finitely many atoms plus densities that are O(e^{-pi|w|}).
- Equal mass and equal first moment give the decay of W at +infinity.
- Integration by parts holds for |Re a| < pi.
- int e^{aw} V is finite for every real 0 < |a| < pi, by Tonelli, so V is in L^1 and the transform is holomorphic.
- Fourier uniqueness gives W = V **almost everywhere**.

**W = V everywhere: GAP.**
- The notes conclude "W = V" and then use W(lambda) = V(lambda) at every lambda, in (c) and (d) and in the statement of
  Theorem B. This needs V finite and continuous, which the notes never establish.
- [27, proof of Thm C (a)] proves "V is finite, continuous and O(e^{-pi|w|})" from
  sup_s rho(s,tau) <= C (tau(1-tau))^(-1/2), and the general notes contain no substitute.
- The issue is real only near the ends of the strip.
  - For tau <= 1 - eps, K_{1-tau} <= tan(pi tau/2)/2 is bounded, so that part of V is continuous by dominated
    convergence with the dominating function C_eps rho, which is in L^1.
  - As tau -> 1, K_{1-tau} becomes an approximate identity of width 1 - tau. When 0 is in spec B (p_m = 1), rho has a
    spike of height ~ (1-tau)^(-1/2) there. The same holds at tau -> 0 for the right edge when 1 is in spec B.
  - Neither local boundedness on Omega nor the marginal identity controls the pointwise behaviour of V.

**Fix (i), complete.** Use Lemma 2.2(3') above.
- For tau in (0,1) \ spec P, int K_{1-tau}(s - w) rho(s,tau) ds <= tau sup_s rho(s,tau) =: Phi(tau), and Phi is in L^1.
  Phi does not depend on w.
- For each fixed tau, w -> int K_{1-tau}(s - w) rho(s,tau) ds is continuous.
- By dominated convergence V is continuous, and then W = V everywhere. Exponential decay follows as in [27].

**Fix (ii), partial and without (3').**
- V is lower semicontinuous, by Fatou: the integrand is >= 0 and continuous in lambda.
- W is continuous and W = V almost everywhere. So for every lambda, V(lambda) <= liminf V(lambda_n) = W(lambda) along
  any sequence in the full-measure set where V = W.
- This yields:
  - V finite everywhere;
  - W >= 0 everywhere, which gives convex order;
  - the equality case (W(lambda) = 0 gives V(lambda) = 0, hence rho = 0 a.e.);
  - W(lambda) >= V(lambda) >= the bound of Theorem C, including its projection reading.
- It does not yield the identity W(lambda) = V(lambda) at every lambda.

**Numerical confirmation of W = V**, computed with careful quadrature:
- the breakpoints include the s-transition points of the root configuration;
- the grid is graded around lambda down to the kernel width;
- the outer integral in tau is adaptive.

| spec B | Edges | Lambda values | Largest abs(W - V) |
|---|---|---|---|
| {0, 0.4, 1} | both | 12 | <= 6.3e-7 |
| {0.15, 0.15, 0.6, 0.85} | both | 14 | <= 4.3e-7 |
| {0, 0.5, 0.5, 0.5} | both | 16 | <= 6.7e-7 |

The lambda values include the critical ones: the eigenvalues of the extreme diagonal blocks, where the spike of rho
meets the narrowing kernel, and points 1e-3 from them. Typical differences are 1e-9 to 1e-11. The worst ones are
quadrature-limited.

An earlier, cruder quadrature gave a spurious difference of 1e-2 at a critical lambda. It came from cancellation in
cosh(pi u) - cos(pi x) for small x, plus an under-resolved tau grid. This is a pitfall for anyone re-checking V
numerically; use 2 sinh^2(pi u/2) + 2 sin^2(pi x/2).

**(c) Conclusion: VERIFIED** given (b). The identity (u - lambda)_+ = (lambda - u)_+ + (u - lambda) and the zero mass
and zero first moment of omega - sigma give the formula.

**(d) Equality: VERIFIED**, using fix (i), or fix (ii) for the direction "equality at one lambda implies [B,g] = 0".
- V(lambda) = 0 with K_{1-tau} > 0 on (0,1) gives rho = 0 a.e.
- Then D = 0 and G(0,t) = 0, so m = 0. Every pair j < k has some tau in (p_j, p_k) \ spec P, so all
  Q_j g Q_k = 0 for j != k. Hence [P,g] = sum_{j != k} (p_j - p_k) Q_j g Q_k = 0.

**Right edge: VERIFIED** (subject to the same gap and fix at tau -> 0 when p_1 = 0).
- P + ig = 1 - (B + ig)^*.
- The density of (P, g), the roots of det(g - s - xi(B - tau)), equals rho(s, 1 - tau): the roots change sign.
- So V^R(lambda) = int int K_tau(s - lambda) rho(s,tau) ds dtau.
- The masses are Tr B and the first moments are Tr(Bg).

**Adversarial search.** Setup:
- Nelder–Mead, 2 seeds, 1184 rounds in all;
- M = 2..5, with generic, edge-snapped and repeated levels;
- both edges, 33 lambda in [-4, 4];
- exact dilogarithm evaluation with the cancellation-free formula, and exact pinching.

The most negative value of min_lambda W / (1 + ||g||) was -3.3e-14, which is rounding: the optimiser is driven to the
commuting equality case. No violation of the Theorem B inequality was found.

### 3.3 Theorem C

**VERIFIED.**
- For s in [l_min, l_max], |s - lambda| <= R.
- K_{1-tau}(u) = sin(pi tau)/(2(cosh(pi u) + cos(pi tau))) is even and decreasing in |u|.
- K_{1-tau}(R) >= sin(pi tau)/(2(1 + cosh(pi R))).
- int rho(s,tau) ds = m(tau), and int sin(pi tau) m(tau) dtau is the stated sum.
- The bound on V(lambda), which is defined by the double integral, does not depend on the gap of 3.2.
- The reading in terms of the sweeps holds via W(lambda) >= V(lambda) (fix (ii)), or W = V (fix (i)).

**Projection specialisation: VERIFIED.** m = ||BgP||_F^2 on (0,1), int_0^1 sin(pi tau) dtau = 2/pi,
[B,g] = BgP - PgB, and ||[B,g]||_F^2 = 2||BgP||_F^2. This gives ||[B,g]||_F^2 / (2 pi (1 + cosh(pi R))).

**Sharpness and numerics.**
- An adversarial minimisation of W(lambda)/bound(lambda) over 2 <= M <= 5 found no ratio below 1. It ran 1184 rounds
  over 2 seeds, with W in the cancellation-free form and exact pinching. The smallest ratio was 1.001026, confirmed in
  50-digit arithmetic: spec B close to {0.986, 1}, ||g|| = 0.03, lambda = 0.
- The ratio tends to 1 in two regimes:
  - g = gamma 1 + eps H (gamma > 0) and lambda -> -infinity, with observed ratio 1.00396 at lambda = -4. Here rho
    lives within O(eps) of s = gamma, so R - |s - lambda| = O(eps), and both kernels behave like
    sin(pi tau) e^{-pi(.)}. The ratio is therefore 1 + O(eps).
  - spec P near 0, with R -> 0. Here K_{1-tau}(u) and sin(pi tau)/(2(1 + cosh(pi R))) both tend to sin(pi tau)/4.
- So the constant 1/(2(1 + cosh(pi R))) cannot be improved by any fixed factor.
- A cosmetic sharpening: R can be replaced by max(|l_min - lambda|, |l_max - lambda|) <= |lambda| + ||g||.

## 4. Corollaries and edge cases

**Corollary 1: VERIFIED.**
- For a projection, B + i g_d = (B + iBgB) (+) (iPgP) on ran B (+) ran P.
- The left sweep gives atoms at the eigenvalues of PgP on ran P, and the right-edge points contribute nothing.

**Corollary 2: VERIFIED.** This is Theorem A at a = 1. It is already stated in [27, Remark RIgeneralB] (finding 3).

**Edge cases:**
- **M = 1, or m = 1** (B scalar): g_d = g, D = 0 and rho = 0 (all roots (mu - s)/(p_1 - tau) are real), so omega =
  sigma and V = 0. Everything holds trivially. In Step 1, K is degenerate (still convex and compact). In Lemma 2.5, delta
  is undefined (cosmetic).
- **l_min = l_max:** g is scalar and everything is trivial.
- **Eigenvalues of B equal to 0 or 1** (p_m = 1 or p_1 = 0): the proofs of A and C are unaffected. For B they are
  exactly the cases where the gap of 3.2 is active (left edge if p_m = 1, right edge if p_1 = 0). Fix (i) covers them.
  They were tested numerically (spec B containing 0 and 1).
- **tau = 0 or 1 endpoints:** null sets for rho. They also lie in spec P, where rho = 0 is set by definition, if p_1 = 0
  or p_m = 1.
- **Eigenvalues of B + ig on the edges** (x = 0 gives an atom, x = 1 gives 0): covered by Lemma 2.1. The moment and
  Laplace formulas include them.

## 5. Numerical experiments of this check (`scripts/check_report/`)

| Script | What it tests | Result |
|---|---|---|
| `hm_common.py` | shared helpers: B = diag(levels) WLOG; rho; m; D; closed-form h_lambda; W | — |
| `c1_pencil.py` | (a) closed-form h_lambda vs quadrature; (b) Lemma 2.2(3'); (c) sup_s rho and int rho ds near the lines | (a) 1.3e-12; (b) max ratio 0.977, no violation; (c) blow-up ~ eps^(-1/2), int rho ds = m(tau) |
| `c2_radon.py` | Theorem A in line form (Radon slices), 4 configurations with interior eigenvalues | max diff 8.4e-11 |
| `c3_bmv2d.py` | Theorem A at complex (a,t), 3 configurations; total mass | rel diff <= 9e-6; mass to 8e-7 relative |
| `c4_lemma25.py` | Lemma 2.5: O(1/theta) rate and explicit constant, Hermitian-type and complex X | rate 1/theta; constant holds |
| `c5b_thmB.py` (args: seed, case 0..3, edges, optional `moments`) | Theorem B: W(lambda) = V(lambda), both edges, critical lambda; masses and first moments | <= 6.7e-7 (quadrature); moments to 12 digits |
| `c6_thmC_adv.py` | adversarial search on the Theorem B inequality and the Theorem C ratio (args: seed, seconds) | 1184 rounds: worst B gap -3.3e-14 (rounding); min C ratio 1.0010 |
| `c7_reproduce.py` | replays an adversarial round and recomputes W in 50-digit arithmetic | confirms the ratios; exposes the double-precision artefact |

Run them from the `scripts/check_report/` folder. Set OPENBLAS_NUM_THREADS=1 when running several processes in
parallel, since multi-threaded BLAS can fail to allocate buffers.
