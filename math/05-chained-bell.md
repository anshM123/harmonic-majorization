# Chained Bell inequalities (Barrett–Kent–Pironio) on maximally entangled states

**Statement.** Let $I_N$ be the chained expression of Barrett, Kent and Pironio [BKP] (PRL 97, 170409 (2006), eq. (3))
with N settings and d outcomes per party; local models give $I_N\ge d-1$. Every projective strategy on the maximally
entangled state $\Phi_D$ (any D ≥ 1, measurements projective on $\mathbb C^D$) satisfies

$$I_N\ \ge\ N(d-1)-\frac{2N}{d}\sum_{m=1}^{d-1}(d-m)\,\frac{\sin(\pi(N-1)m/(Nd))}{\sin(\pi m/d)},$$

the value of the BKP measurements on $\Phi_d$. For d ≥ 3, equality holds if and only if d divides D and the strategy is
the BKP strategy tensored with $1_{D/d}$, up to local unitaries $u\otimes\bar u$.
- Proved, computer-assisted (interval arithmetic), for 3 ≤ N ≤ 8 with 3 ≤ d ≤ 10, and for N ∈ {3, 4, 5} with
  11 ≤ d ≤ 24. The case d = 2 is Wehner's bound (PRA 73, 022110 (2006)); the case N = 2 is the CGLMP result of [27].
- This confirms, for projective measurements, the maximally-entangled values that Tavakoli, Zohren and Pawłowski
  (J. Phys. A 49, 145301 (2016)) derived under their conjecture that the BKP-type measurements are optimal.

**How the method is used.** A twirl reduces every strategy to a clock model (Theorem 1). The chained kernel is a
positive combination of adjacent block pairings, by concavity of k ↦ sin(πk/N) (Lemma 2). Each pairing is bounded by
the continuum theorem of [27], which rests on Theorem B of this repository for projections ([27, Theorem C]).
Certificates with 1–18 rational cells close the discretisation (Theorem 5).

**What remains open.**
- The conjecture of Tavakoli, Zohren and Pawłowski itself: whether the BKP-type measurements stay optimal for the
  optimal, non-maximally entangled state.
- The maximal quantum value over all states, listed as unknown by Kaniewski et al. (Quantum 3, 198 (2019)). It is
  strictly smaller than the value above: for N = d = 3, the BKP measurements on a non-maximally entangled state give
  0.71296 < 0.75678 (Section 6).
- POVMs on $\Phi_D$, including maximally entangled states of lower Schmidt rank inside larger local spaces; and all
  other (N, d) (Conjecture 6).

**Status:** CHECKED (an internal independent check, not external refereeing). The check found no mathematical error
and re-verified the certificates with independent ball-arithmetic code
([`../bell/independent-check/`](../bell/independent-check/README.md)). The proof depends on our unrefereed OQP 27
results: [27] (clock model, strip inequality, continuum theorem, cell embeddings; github.com/anshM123/IQOQI-OQP-27) and
[CL] (discrete rearrangement inequality and commutative sector; github.com/anshM123/CGLMP). A search of the papers
citing Barrett–Kent–Pironio or Tavakoli–Zohren–Pawłowski, and of recent work (2026-10-08), found no earlier proof for
N ≥ 3 and d ≥ 3 (Section 7).

Scripts, logs and the commands that produced them: [`../bell/README.md`](../bell/README.md).

## References and conventions

- **[27]** A. Mishra, A. Senthilkumar, *A strip inequality for projections, two-variable BMV positivity, and optimal
  CGLMP measurements in every dimension*, github.com/anshM123/IQOQI-OQP-27, `papers/math/main.tex`.
- **[CL]** A. Mishra, A. Senthilkumar, *A rearrangement inequality for the discrete Hilbert transform, with an
  application to optimal CGLMP measurements*, github.com/anshM123/CGLMP, `paper-classical-all-d/main.tex`.
- **[BKP]** J. Barrett, A. Kent, S. Pironio, PRL 97, 170409 (2006), arXiv:quant-ph/0605182.
- **[TZP]** A. Tavakoli, S. Zohren, M. Pawłowski, J. Phys. A 49, 145301 (2016), arXiv:1506.04331.
- **[ADGL]** A. Acín, T. Durt, N. Gisin, J. I. Latorre, PRA 65, 052325 (2002).
- **[W]** S. Wehner, PRA 73, 022110 (2006).
- Further literature is listed in Section 7.

Numbers of theorems, propositions, lemmas and equations of [27] and [CL] are those of the LaTeX files named above, as
compiled. The items cited are:

| Cited as | In [27] | LaTeX label |
|---|---|---|
| [27, Theorem A] | optimality of the CGLMP measurements on maximally entangled states (computer-assisted) | `thm:A` |
| [27, Theorem B] | rigidity (computer-assisted) | `thm:B` |
| [27, Theorem C] | strip inequality | `thm:C` |
| [27, Definition 2.3] | Q-configuration | `def:Qconf` |
| [27, Theorem 4.2] | continuum theorem | `thm:continuum` |
| [27, Lemma 4.4] | the classical value | `lem:classical` |
| [27, Proposition 5.2] | cell inequalities | `prop:cellineq` |
| [27, Section 6] | the cone condition for every d (computer-assisted) | `sec:cone` |
| [27, Proposition 6.2] | Regime I certificates, 2 ≤ d ≤ 200 | `prop:regimeI` |
| [27, Lemma 7.2] | residue lemma | `lem:residue` |
| [27, Proposition 7.3] | pointwise commutation forces commutation | `prop:residue` |
| [27, eq. (1.9)] | definition of $\Phi_c$ | `eq:Phic` |

| Cited as | In [CL] | LaTeX label |
|---|---|---|
| [CL, Theorem A] | rearrangement inequality for the discrete Hilbert transform | `thm:A` |
| [CL, Corollary C] | equal-link strategies: optimality and rigidity | `cor:C` |
| [CL, Lemma 2.1] | chain form of the CGLMP expression | `lem:chain` |
| [CL, Appendix A] | proof of the reduction to the clock model | `app:reduction` |

Conventions. τ and tr denote normalised traces ($\operatorname{tr}=\operatorname{Tr}/M$ on $M_M(\mathbb C)$).
$\Delta_d=\lbrace\ell\in\mathbb R^d:\ell_r\ge0,\ \sum_r\ell_r=d\rbrace$ is the simplex of cell vectors and $\Delta_d^+$
the set of cell vectors with all entries positive, as in [27]. We write $S_N(R)$ for the value of $I_N$ on a strategy
R (as in the scripts). A *chained Q-configuration* on $\mathbb C^M$ is a family $(Q_x)_{x\in\mathbb Z_L}$, $L=2Nd$, of
orthogonal projections such that for every site $k\in\lbrace0,\dots,d-1\rbrace$ the 2N projections $Q_{k-da}$,
$a\in\mathbb Z_{2N}$, are pairwise orthogonal with sum 1 (the analogue of [27, Definition 2.3], where N = 2). The
*window* is the configuration $Q_x=1_{[0,d)}(x)$ on $\mathbb C^1$, and its rotations are the windows $1_{[c,c+d)}$.

## 1. The chained Bell expression of Barrett, Kent and Pironio

With settings $A_1,\dots,A_N$ (Alice) and $B_1,\dots,B_N$ (Bob), outcomes in $\mathbb Z_d$, [BKP] define

$$I_N=\langle[A_1-B_1]\rangle+\langle[B_1-A_2]\rangle+\langle[A_2-B_2]\rangle+\dots+\langle[A_N-B_N]\rangle+\langle[B_N-A_1-1]\rangle,$$

where $[X]=X\bmod d\in\lbrace0,\dots,d-1\rbrace$ and $\langle\cdot\rangle$ is the expectation. Local models satisfy
$I_N\ge d-1$ (around the cycle the links add up to −1 mod d). [BKP] use $\Phi_d$ and the Fourier bases with phases
$\alpha_k=(k-\frac12)/N$ (Alice) and $\beta_l=l/N$ (Bob), find $I_N=2\gamma/N+O(1/N^2)$,
$\gamma=\frac{\pi^2}{4d^2}\sum_j j/\sin^2(\pi j/d)$, and do not claim optimality. For N = 2, $I_2=(d-1)(4-I_d)/2$ is
the CGLMP expression in the chain form of [CL, Lemma 2.1].

## 2. Reduction to a clock model

Put $R_{2k-1}=U_{A_k}$, $R_{2k}=(U_{B_k})^{\mathsf T}$ ($U=\sum_aw^aP_a$, $w=e^{2\pi i/d}$, transpose in the Schmidt
basis), $L=2Nd$, $z=e^{2\pi i/L}$, $\zeta=e^{i\pi/N}$.

**Theorem 1.** For every projective strategy on $\Phi_D$ there are M (M = 2ND works) and unitaries
$V_0,\dots,V_{d-1}$ on $\mathbb C^M$ with $V_k^{2N}=1$ (the reduced family) such that

$$I_N=N(d-1)-\frac2d\,F_N(V),\qquad F_N(V)=\sum_{0\le j<k\le d-1}\operatorname{Re}\bigl[h_{k-j}\,\tau(V_jV_k^\ast)\bigr],\qquad h_m=\frac{-iN\,e^{i\pi m(N-1)/(Nd)}}{\sin(\pi m/d)}.$$

Equivalently, with the Q-configuration $Q_{k-da}$ ($a\in\mathbb Z_{2N}$, $x\in\mathbb Z_L$) of the spectral
projections of $V_k$ ($V_k=\sum_a\zeta^{-a}Q_{k-da}$),

$$F_N=\frac12\sum_{x\not\equiv y\ (\mathrm{mod}\ d)}\tau(Q_xQ_y)\,K_N(x-y),\qquad K_N(u)=\frac{N\sin(\pi(N-1)u/(Nd))}{\sin(\pi u/d)}.$$

*Proof.* As [CL, Appendix A], with $\rho(R)=(R_2,\dots,R_{2N},wR_1)$ of order 2Nd ($S_N$ is ρ-invariant because the
links are cyclically permuted, the last one with its phase), $Z=W^{2N}$, $H_k=\ker(Z-w^k)$, $V_k:=z^{-k}W_k$. For the
covariant strategy the 2N links have equal traces, and

$$\tau\bigl(L^{(0)}_n\bigr)=\frac{z^{-n}}d\Bigl[\sum_{k\ge n}\operatorname{tr}(V_{k-n}V_k^\ast)+e^{i\pi/N}\sum_{k<n}\operatorname{tr}(V_{k-n+d}V_k^\ast)\Bigr]\qquad(z^d=e^{i\pi/N}).$$

With $c_n=-1/(1-w^{-n})=\frac i2\,e^{i\theta_n}/\sin\theta_n$ ($\theta_n=\pi n/d$), $2Nc_nz^{-n}=-h_n$ and
$h_{d-m}e^{i\pi/N}=\overline{h_m}$, which gives the formula. Projection form:
$\operatorname{Re}[h_m\zeta^s]=N\sin(\theta_m(1-1/N)+\pi s/N)/\sin\theta_m$ depends only on $u=m-ds$: with $y-x=u$ it
equals $N\sin(\pi(N-1)u/(Nd))/\sin(\pi u/d)$. ∎

**Numerical checks.**
- $S_N$ of the covariant strategy against direct evaluation ([`ch_clock.py`](../bell/scripts/ch_clock.py), 12 digits).
- The twirl reduction for arbitrary random strategies (`ch_twirl_check.py N d D 3`, log
  [`ch_twirl_check_N3-5_d3-5_D3-4-5-7.log`](../bell/logs/ch_twirl_check_N3-5_d3-5_D3-4-5-7.log)): N = 3, 4, 5,
  d = 3, 4, 5, D = 3, 4, 5, 7, three random strategies per case, 108 in all, with random ranks (zero allowed), including
  D not divisible by d. $S(R)$ and $N(d-1)-\frac2dF_N(V)$ agree to all 12 printed decimals, and
  $|V_k^{2N}-1|\le2.2\cdot10^{-14}$. The older log [`ch_twirl_check.log`](../bell/logs/ch_twirl_check.log) has only
  (N, d, D) = (3, 3, 2) and (4, 3, 2).
- An independent check, with an explicit basis for the twirl and no eigensolver, agrees to $4.1\cdot10^{-14}$ on 20
  random strategies ([`ref_reduction.log`](../bell/independent-check/ref_reduction.log)).
- BKP's measurements give the value of the natural (V = 1) strategy (`ch_bkp.py 2 8 2 10`, log
  [`ch_bkp_N2-8_d2-10.log`](../bell/logs/ch_bkp_N2-8_d2-10.log): N = 2..8, d = 2..10, largest difference
  $9.1\cdot10^{-15}$).

## 3. The kernel is a positive combination of adjacent block pairings

Let $\kappa(u)=\cot(\pi u/L)$ on $\mathbb Z_L$ ($\kappa(0)=0$),
$\langle X,Y\rangle:=\sum_{x,y}\kappa(x-y)\,\tau(X_xY_y)$ for operator fields on $\mathbb Z_L$, and
$G_k:=\langle Q,Q_{\cdot+kd}\rangle$.

**Lemma 2.**
1. $K_N(u)=-\frac12\sum_{k=1}^{2N-1}\sin(\pi k/N)\,\kappa(u-kd)$.
2. $F_N=\frac12\sum_{k=1}^{N-1}\sin(\pi k/N)\,\bigl[G_k-d\,\kappa(kd)\bigr]$.
3. $\sin(\pi k/N)=\sum_{s=1}^{\lfloor N/2\rfloor}\mu_s\,m_s(k)$ for $1\le k\le N-1$, with
   $m_s(k)=\min(k,N-k,s,N-s)$ and $\mu_s=2(1-\cos(\pi/N))\sin(\pi s/N)$ (halved for s = N/2); all $\mu_s>0$.
4. Hence $F_N=\frac12\sum_s\mu_s\langle X_s,Y_s\rangle-\frac d2\sum_{k=1}^{N-1}\sin(\pi k/N)\,\kappa(kd)$, where
   $(X_s)_x=\sum_{j=0}^{s-1}Q_{x+jd}$ and $(Y_s)_x=\sum_{j=s}^{N-1}Q_{x+jd}$ are orthogonal projections, $X_sY_s=0$.

*Proof.*
1. $f(\alpha)=\sin((N-1)\alpha)/\sin(N\alpha)$ is 2π periodic, decays like $e^{-|\operatorname{Im}\alpha|}$ at ±i∞,
   and has simple poles at $\alpha_k=\pi k/N$ with residues $-\sin(\pi k/N)/N$ (k = 0 and N are removable). So
   $f=-\frac1{2N}\sum_k\sin(\pi k/N)\cot((\alpha-\alpha_k)/2)$ (no constant, by the decay). Put $\alpha=\pi u/(Nd)$.
2. With $T_k:=\sum_{x,y}\tau(Q_xQ_y)\,\kappa(x-y-kd)$: the omitted pairs x ≡ y mod d contribute only x = y (outcomes
   of one site are orthogonal), so $\sum_{x\not\equiv y}\tau(Q_xQ_y)\,\kappa(x-y-kd)=T_k+d\,\kappa(kd)$; reindexing,
   $G_k=T_{-k}=-T_k$ (tracial symmetry and oddness of κ). Pair k with 2N − k (sin changes sign, $G_{2N-k}=-G_k$,
   $\kappa((2N-k)d)=-\kappa(kd)$).
3. sin(πk/N) is strictly concave on {0, …, N} and vanishes at 0 and N; $m_s$ has second difference −1 at k = s and
   k = N − s. So $\mu_s=-\Delta^2\sin(\pi s/N)=2(1-\cos(\pi/N))\sin(\pi s/N)>0$.
4. $\langle Q_{\cdot+jd},Q_{\cdot+j'd}\rangle=G_{j'-j}$, and the number of pairs $(j,j')$ with
   $0\le j<s\le j'\le N-1$ and $j'-j=k$ is $m_s(k)$. $X_s$ and $Y_s$ are sums of distinct outcomes of one site. ∎

Checked numerically by fitting $F_N$ against the $G_k$ on random configurations (coefficients $\sin(\pi k/N)/2$ to
$10^{-14}$; no log archived). An independent check of the full decomposition agrees to $3.8\cdot10^{-14}$
([`ref_cellineq.log`](../bell/independent-check/ref_cellineq.log)).

The point: for N = 2 this is the single pairing ⟨A, B⟩ of [27]. For N ≥ 3 the kernel involves pairings of the
configuration with several rotations; the naive bound of each $\langle A,Q_{\cdot+kd}\rangle$ by the continuum value is
**not** attained by the window for k ≥ 2 (the arcs are not adjacent), but the block pairings $\langle X_s,Y_s\rangle$
**are** all attained by windows ($Y_s$ is the arc immediately before $X_s$). Concavity of k ↦ sin(πk/N) is exactly what
makes the decomposition positive.

## 4. Commuting strategies (all N, all d)

**Theorem 3.** If the reduced family commutes (equivalently the strategy has equal links), then

$$I_N\ \ge\ I_N^{\mathrm{BKP}}(d):=N(d-1)-\frac2d\,F_{\mathrm{nat}},\qquad F_{\mathrm{nat}}:=\sum_{m=1}^{d-1}(d-m)\,\frac{N\sin(\pi(N-1)m/(Nd))}{\sin(\pi m/d)},$$

with equality iff every joint eigenvector has a one-step configuration $a_k=c+[k\ge r]$ ($c\in\mathbb Z_{2N}$); then the
strategy is BKP ⊗ 1 up to $u\otimes\bar u$.

*Proof.* $F_N(V)=\sum_a\tau(\Pi_a)\,F_N(a)$, and by Lemma 2(4)
$F_N(a)=\frac12\sum_s\mu_s\langle X_s(a),Y_s(a)\rangle-\mathrm{const}$, where $X_s(a)=\bigcup_{j<s}(S_a-jd)$ and
$Y_s(a)=\bigcup_{s\le j<N}(S_a-jd)$ are disjoint subsets of $\mathbb Z_L$ of sizes sd and (N − s)d
($S_a=\lbrace k-da_k\rbrace$ a transversal). [CL, Theorem A] (for $\mathbb Z_L$ and these sizes) bounds each
$\langle X_s,Y_s\rangle$ by $\Phi_L(sd,(N-s)d)$, with equality iff the two sets are adjacent intervals; the window
$S_a=[c,c+d)$ attains all of them. Equality in s = 1 forces $S_a$ to be an interval, i.e. a one-step. BKP's
measurements give the natural value (direct evaluation, `ch_bkp.py 2 8 2 10`, log
[`ch_bkp_N2-8_d2-10.log`](../bell/logs/ch_bkp_N2-8_d2-10.log): BKP's $I_N$ equals $I_N^{\mathrm{BKP}}(d)$ to
$9.1\cdot10^{-15}$ for N = 2..8, d = 2..10). This compares values only; for the certified (N, d) with d ≥ 3,
Theorem 5′ identifies BKP's strategy with the natural one (remark after Theorem 5′). Rigidity as in [CL, Corollary C]:
$R^{(0)}_1=X$ and $R^{(0)}_1(R^{(0)}_2)^{-1}=z^{-1}\bigl(1+(z^d-1)|0\rangle\langle0|\bigr)$ with $z^d=e^{i\pi/N}\ne1$
generate $M_d$. ∎

## 5. Non-commuting strategies: cells and a certified linear programme

**Proposition 4 (cell inequalities).** For every cell vector $\ell\in\Delta_d$ (extended d-periodically to
$\mathbb Z_L$, circle of length L) and every $s\in\lbrace1,\dots,\lfloor N/2\rfloor\rbrace$, every chained
Q-configuration satisfies

$$Q^{\ell,s}(Q):=\sum_{x,y}K^\ell(x,y)\,\tau\bigl((X_s)_x(Y_s)_y\bigr)\ \le\ L^2\,\Phi_c\Bigl(\frac s{2N},\frac{N-s}{2N}\Bigr),$$

with equality for windows, where $K^\ell(x,y)=\int_{I_x}\int_{I_y}\cot(\pi(\xi-\eta)/L)$ and $\Phi_c$ is
[27, eq. (1.9)]. Averaged over the d rotations, $\frac1d\sum_jQ^{\ell,s}(\mathrm{rot}_jQ)=\sum_u\tilde k^{\ell,s}(u)\,C(u)$
with $C(u)=\sum_x\tau(Q_xQ_{x+u})$ and

$$\tilde k^{\ell,s}(u)=\frac1L\sum_{j<s\le j'<N}\ \sum_xK^\ell\bigl(x,\,x+u-(j'-j)d\bigr).$$

*Proof.* [27, Proposition 5.2] verbatim gives $Q^{\ell,s}(Q)=L^2\int\tau(X_sg)\,dm$ with g the conjugate function of
the projection step field $Y_s$; $0\le X_s\le1-Y_s$, $\int Y_s\,dm=\frac{N-s}{2N}\,1$ exactly (each site is a PVM),
$\int\tau(X_s)=s/2N$. [27, Theorem 4.2 and Lemma 4.4] (with α = s/2N, β = (N − s)/2N, α + β = 1/2) give the bound (the
bracket vanishes since $\int Y_s$ is scalar). Windows: $X_s$ and $Y_s$ are adjacent arcs of lengths sπ/N and
(N − s)π/N. The average: $K^\ell(x+d,y+d)=K^\ell(x,y)$. ∎

Valid for every configuration also: $C(u)\ge0$, $C(-u)=C(u)$, $C(0)=d$, $C(u)=0$ for u ≡ 0 mod d, u ≠ 0, and
$\sum_{t\in\mathbb Z_{2N}}C(r+td)=d$ for r ≢ 0 mod d. And $F_N=\frac12\sum_{u\not\equiv0\ (\mathrm{mod}\ d)}K_N(u)\,C(u)$.

**Theorem 5 (optimal value on maximally entangled states, projective measurements; computer-assisted).** For each
(N, d) in the list below, every projective strategy on a maximally entangled state $\Phi_D$ of any local dimension D
satisfies

$$I_N\ \ge\ I_N^{\mathrm{BKP}}(d)=N(d-1)-\frac{2N}{d}\sum_{m=1}^{d-1}(d-m)\,\frac{\sin(\pi(N-1)m/(Nd))}{\sin(\pi m/d)},$$

and the BKP measurements attain it.
- 2 ≤ N ≤ 8 and 2 ≤ d ≤ 10 (N = 2 is CGLMP: a re-proof of [27, Theorem A] for d ≤ 10 with certificates of 1–6
  cells, against 3.5d–6.3d cells in [27, Proposition 6.2]);
- N = 3, 4 and 5 with 11 ≤ d ≤ 24.

The new cases are those with N ≥ 3 and d ≥ 3: 3 ≤ N ≤ 8 with 3 ≤ d ≤ 10, and N ∈ {3, 4, 5} with 11 ≤ d ≤ 24. For
d = 2 the bound is Wehner's [W], and for N = 2 it is the CGLMP result of [27]; in these cases the certificates re-prove
known results.

**Scope.** $\Phi_D$ is the maximally entangled state of full Schmidt rank on $\mathbb C^D\otimes\mathbb C^D$, and the
measurements are PVMs on $\mathbb C^D$ (any ranks, zero outcomes allowed). No balanced-rank or divisibility assumption
is made; d | D is a consequence of equality (Theorem 5′). POVMs are not covered, since Naimark dilation does not
preserve maximal entanglement. Nor is a maximally entangled state of Schmidt rank D′ < D inside larger local spaces,
measured by PVMs on those spaces: that is the POVM case on $\Phi_{D'}$.

*Proof.* [`ch_cert.py`](../bell/scripts/ch_cert.py) produces, for each (N, d), rational cells $\ell_k=n_k/8$
($n_k\in\mathbb Z^d$, entries ≥ 1), weights $y_k>0$ and multipliers $\lambda_r$ such that the dual slacks

$$z_u:=\Bigl(\sum_ky_k\,\tilde k^{\ell_k,s_k}+\sum_r\lambda_re_r-\tfrac12K_N\Bigr)(u)$$

($e_r$ the coefficient vector of the constraint $\sum_tC(r+td)=d$; folded to the representatives 1 ≤ u ≤ Nd,
u ≢ 0 mod d) vanish for u ≤ d − 1 and are > 0 for u ≥ d + 1. The verification is in interval arithmetic (mpmath.iv,
120 bits, outward rounding): every Clausen value $\mathrm{Cl}_2(2\pi j/(8L))$ is enclosed by the Bernoulli series with
explicit tail used in [27, Proposition 6.2] (the program `verify_cone.py` of [27]; re-implemented in `ch_cert.py`); the
equalities on u ≤ d − 1 are solved exactly inside an interval box by a least-norm residual correction with a certified
approximate inverse; $y_k>0$ and $z_u>0$ are checked on the box. Then for every C satisfying the valid constraints,

$$F_N(C)=\sum_ky_k\Bigl(\sum_{u\in\mathbb Z_L}\tilde k_k(u)C(u)-\tilde k_k(0)\,d\Bigr)+\sum_r\lambda_r\,d-\sum_{u\ge d+1}z_uC(u)\ \le\ \sum_ky_k\bigl(b_k-\tilde k_k(0)\,d\bigr)+\sum_r\lambda_r\,d,$$

with $b_k=L^2\Phi_c(s_k/2N,(N-s_k)/2N)$, and the right-hand side equals $F_N(C_{\mathrm{win}})=F_{\mathrm{nat}}$
($C_{\mathrm{win}}$ vanishes where $z_u>0$ and every cell inequality is tight at windows). For d = 2 no cell is
needed. Theorem 1 transfers the bound to strategies. ∎

**Certificates (cells used).** N = 3: d = 3: 1, 4: 2, 5: 3, 6: 3, 7: 5, 8: 5, 9: 5, 10: 4, 11: 6, 12: 7, 13: 8, 14: 10,
15: 11, 16: 9, 17: 10, 18: 12, 19: 14, 20: 12, 21: 14, 22: 13, 23: 16, 24: 15; N = 4: d = 3..24 with 1–18 cells;
N = 5: d = 3..24 with 1–17 cells; N = 6..8: d = 3..10 with 1–7 cells; d = 2: no cell. For d = 24 (N = 3, 4, 5) the
first pool of 41 cells admitted no strictly complementary certificate and a pool of 121 cells did. Logs
(`ch_cert_batch.py Nmin Nmax dmin dmax`):
- [`ch_cert_batch_N2-6_d2-10.log`](../bell/logs/ch_cert_batch_N2-6_d2-10.log): N = 2..6, d = 3..10. Its d = 2 entries
  stopped at an exception for empty certificates (since fixed) and are superseded by the next log.
- [`ch_cert_batch_d2.log`](../bell/logs/ch_cert_batch_d2.log): d = 2, N = 2..8, no cell.
- [`ch_cert_batch_N7-8_d3-10.log`](../bell/logs/ch_cert_batch_N7-8_d3-10.log),
  [`ch_cert_batch_N3-4_d11-14.log`](../bell/logs/ch_cert_batch_N3-4_d11-14.log),
  [`ch_cert_batch_N3-5_d15-24.log`](../bell/logs/ch_cert_batch_N3-5_d15-24.log).
- [`ch_cert_batch_N5_d11-14.log`](../bell/logs/ch_cert_batch_N5_d11-14.log) (run 2026-10-08): N = 5, d = 11..14, all
  CERTIFIED with the first pool (6, 8, 11 and 10 cells).

An internal independent check (not external refereeing) re-verified 26 certificates with an independent
ball-arithmetic verifier (python-flint, independent Clausen values and kernel code), including all three d = 24
certificates and N = 5, d = 11..14 ([`ref_cert_batch.log`](../bell/independent-check/ref_cert_batch.log),
[`ref_cert_test_N3d3.log`](../bell/independent-check/ref_cert_test_N3d3.log),
[`ref_cert_theirs_d24.log`](../bell/independent-check/ref_cert_theirs_d24.log)).

**Independent checks.** The cell inequalities hold on 150 random and near-optimal chained configurations (largest
violation $3.8\cdot10^{-13}$, rounding; `ch_validate.py 30`, log [`ch_validate.log`](../bell/logs/ch_validate.log):
(N, d) = (3,3), (4,3), (3,4), (5,3), (4,4), 30 configurations each). This test uses the same kernel code as the
certificate; an independent kernel and an adversarial maximisation of the cell functionals are in
[`ref_cellineq.log`](../bell/independent-check/ref_cellineq.log). The same programme with the certificate cells of
[27] reproduces CGLMP for d = 2..6 (log
[`ch_cone_N2_d2-6_oqp27cells.log`](../bell/logs/ch_cone_N2_d2-6_oqp27cells.log), recorded when the `bell/` folder was
assembled). Non-commuting searches (L-BFGS, M = 2..4, 12 restarts per M; (N, d) = (3,3), (3,4), (4,3), (3,5)) never
exceed $F_{\mathrm{nat}}$ beyond rounding (largest excess $2.7\cdot10^{-15}$;
[`ch_local_search.log`](../bell/logs/ch_local_search.log)).

**Theorem 5′ (rigidity).** For the certified (N, d) with d ≥ 3: if a projective strategy on $\Phi_D$ attains
$I_N^{\mathrm{BKP}}(d)$, then d | D and, up to a local unitary $u\otimes\bar u$, it is the BKP strategy tensored with
$1_{D/d}$.

*Proof.* As [27, proof of Theorem B]: every certificate uses at least one cell, all certificate cells lie in
$\Delta_d^+$, and equality makes the corresponding averaged cell inequality, hence (by averaging) every rotated one,
tight; so the strip defect of the projection field $Y_s$ vanishes a.e. and $[Y_s(\theta),g_{Y_s}(\theta)]=0$ a.e. The
cotangent residue argument of [27, Proposition 7.3] (which uses the residue lemma, [27, Lemma 7.2]) gives
$[Y^{(x)},Y^{(x')}]=0$ for all arcs x, x′. $Y^{(x)}$ is a sum of N − s cyclically consecutive outcomes of the site
x mod d, and products of overlapping windows of a PVM produce the single outcomes
($(P_a+\dots+P_{a+w-1})(P_{a+1}+\dots+P_{a+w})=P_{a+1}+\dots+P_{a+w-1}$); so all $Q_x$ commute. Theorem 3 then gives
one-step configurations and the conclusion. ∎

**Remark.** BKP's own strategy attains $I_N^{\mathrm{BKP}}(d)$ on $\Phi_d$ (Section 4), so Theorem 5′ applies to it:
for the certified (N, d) with d ≥ 3 it is the natural chain $R^{(0)}$ up to $u\otimes\bar u$, which is what justifies
the name "BKP strategy" in the equality case.

**Conjecture 6.** Theorem 5 (and 5′) hold for all N ≥ 2 and d ≥ 2. Evidence: the floating-point LP of
[`ch_cone.py`](../bell/scripts/ch_cone.py) closes (gap ≤ $3\cdot10^{-13}$) for every tested case, including
N = 3, d = 16 and N = 4, d = 12 with 41 random cells (log
[`ch_cone_N3_d16_N4_d12.log`](../bell/logs/ch_cone_N3_d16_N4_d12.log), recorded when the `bell/` folder was assembled).
A proof for all d would follow the three regimes of [27, Section 6] (Dirichlet cells, Gaussian modulation, Taylor
models), now with ⌊N/2⌋ block types; the very small certificates suggest that the large-d analysis is easier here than
for CGLMP.

## 6. What this settles, and what it does not

Does this settle an open bound? Only on maximally entangled states, and only for projective measurements. [BKP]
evaluate their measurements on $\Phi_d$ and make no optimality claim. Tavakoli, Zohren and Pawłowski [TZP] study the
same expression (an equivalent relabelled form $B_{N,d}$) and give the maximally-entangled value of these measurements
(their kernel $N\sin((N-1)\pi(k-l)/(dN))/\sin(\pi(k-l)/d)$ is our $K_N$). Their see-saw numerics optimise over Alice's
measurements, Bob's measurements and the shared state (D = d), and from them they "conjecture the above measurements as
optimal". Under this conjecture they derive two values: the maximally-entangled value $\frac1d\sum_{k,l}M_{kl}$ (their
Sec. III.A) and the optimal value over all states, from the smallest eigenvalue of M (their Sec. III.B; checked with
SDP numerics for N, d = 2..10, with no certified upper bound reported). Here M is the d × d matrix of [TZP], eq. (7).
- Theorem 5 proves the maximally-entangled, projective case of this measurement-optimality conjecture, for the
  certified (N, d): TZP's maximally-entangled values are the true optima over PVM strategies on $\Phi_D$, for every D.
  Theorem 5′ adds uniqueness.
- It does not prove the conjecture. Its content for the optimal, non-maximally entangled state remains open.
- It does not settle the maximal quantum value over all states and dimensions, which Kaniewski et al. (Quantum 3, 198
  (2019)) list as unknown. That value is strictly smaller. Example: for N = d = 3, BKP's measurements on the lowest
  eigenvector of M give $I_N=0.71296$, below the maximally-entangled optimum 0.75678
  ([`ref_tzp_gap.log`](../bell/independent-check/ref_tzp_gap.log); the gap is 0.118 for (N, d) = (3, 4), 0.448 for
  (8, 10) and 3.47 for (3, 24)). So for N = d = 3 the all-state minimum is at most 0.71296, and its exact value is open.
  For N = 2 the all-state minimum is also strictly smaller for d ≥ 3 [ADGL].
- POVMs on $\Phi_D$ are not covered.

TZP also report that for N > 2 the normalised entanglement entropy of the optimal state found under their conjecture
tends to 1 as d grows. They prove this limit for an approximate state, and they note that it does not make the optimal
state maximally entangled. The gap in the Bell value grows with d (examples above).

## 7. Literature and priority

- [BKP] (2006): introduces the inequality, uses $\Phi_d$ with $\alpha_k=(k-\frac12)/N$, $\beta_l=l/N$; no optimality
  claim.
- Tavakoli, Zohren, Pawłowski, J. Phys. A 49, 145301 (2016): simplified equivalent form; maximally-entangled value of
  the BKP measurements, limits N → ∞ and d → ∞; optimality of these measurements **conjectured** from see-saw numerics
  over the measurements and the shared state (D = d). (Theorem 5 proves the maximally-entangled, projective case of this
  conjecture in the certified cases. The conjecture itself remains open.)
- Kaniewski et al., Quantum 3, 198 (2019): list the maximal quantum value of most d-outcome CHSH generalisations,
  including CGLMP and BKP, as not known. This is the value over all states; it is not settled here (Section 6).
- Šupić, Augusiak, Salavrakos, Acín, NJP 18, 035013 (2016): self-testing with chained inequalities, two outcomes only.
- Wehner, PRA 73, 022110 (2006): Tsirelson bounds for two-outcome chained expressions (N cos(π/2N)); in our setting
  d = 2 is trivial on $\Phi_D$ (no cell needed).
- Salavrakos et al., PRL 119, 040402 (2017) (SATWAP): related cotangent-weighted family, solved by sum of squares.
- Li, Xiang, Tura, He, PRL 135, 060201 (2025): chain inequalities with m inputs and d outputs as randomness witnesses;
  no maximal value.
- Coccia, Padovan, Vallone, arXiv:2606.21626 (2026): systematic sum-of-squares derivation of Tsirelson bounds for qubits
  and qudits; the qudit families treated have two settings per party (or MUB-type settings); BKP is not treated.
- Das, Maity, Saha, Majumdar, Quantum 6, 716 (2022): temporal inequalities with a maximally mixed input, which is
  equivalent to a Bell expression on a maximally entangled state; SATWAP-type coefficients and two settings per side,
  with sum-of-squares bound 4(d − 1). Not the BKP chain.
- Xiao et al., PRA 107, 052415 (2023) and Paul, Munshi, Pan, PRA 113, 042211 (2026): the chained inequality with two
  outcomes only.

**Priority check (2026-10-08).** Sources searched:
- the citation lists of [BKP] (Semantic Scholar, 168 citing works; OpenAlex, 227 listed works; both screened in full by
  title; three keyword searches inside its 264 Google Scholar citers) and of [TZP] (7 citing works, the same in Semantic
  Scholar, OpenAlex, OpenCitations and Google Scholar, all read);
- web and arXiv searches;
- the index of a public collection of recent preprints (github.com/openai/math), screened by title.

No proof of the optimality of the BKP measurements on maximally entangled states, or of the bound
$I_N\ge I_N^{\mathrm{BKP}}(d)$ on $\Phi_D$, was found for N ≥ 3 and d ≥ 3, even for special (N, d). Limits of the
check: the Google Scholar citers of BKP were not paged one by one; a paper that cites neither BKP nor TZP and uses other
language could be missed; non-English literature was not searched; the check covers work public on 2026-10-08.
