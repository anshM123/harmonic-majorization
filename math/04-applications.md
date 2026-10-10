# Applications

## 1. Word averages and the refined BMV inequality (IQOQI Open Quantum Problem 40)

Full write-up, independent check and code: [anshM123/OQP-40](https://github.com/anshM123/OQP-40).

For PSD matrices A, B, let $\mathcal A_{n,m}(A,B)$ be the average trace over all words with n letters A and m letters B,
and let $E_A$ be the pinching onto the eigenspaces of A. Choosing the test functions $\tau^ns^m$ in Theorem A gives the
**exact gap identity**

$$\mathcal A_{n,m}(A,B)-\operatorname{Tr}(A^nE_A(B)^m)=m(m-1)\iint s^{m-2}\tau^n\rho_{B,A}(s,\tau)\,ds\,d\tau\ \ge0 .$$

Consequences:
- **Dinh's pinching conjecture** (arXiv:2605.17782) holds for all word lengths. Dinh proved m = 2. The inequality also
  follows in a few lines from Heinävaara's Corollary 1.2 (arXiv:2310.03227v1).
- **Strictness:** the gap is strictly positive iff AB ≠ BA.
- **Corrected upper bounds:** the gap is pinned between two explicit multiples of the two-letter gap. This gives a
  correct substitute for the upper half of OQP 40, which is false (Cha–Lee, arXiv:2603.19927).
- **Lower half of OQP 40** (𝒜 ≥ Tr exp(n log A + m log B)):
  - proved when min(n, m) ≤ 2;
  - proved at (3,3), (3,4), (4,3), (4,4), (4,6) and (6,4), in any dimension. The (3,3) case rests on
    Tr(A³B³) + 2 Re Tr(A²BAB²) ≥ 3 Tr((AB)³), which is machine-checked in Lean;
  - proved at (n, 3) and (3, n) for n = 1, …, 34 and n = 36, 39, …, 63, in any dimension (an exact computation for each n);
  - proved at (n, 4) and (4, n) for n = 3, …, 12 and n = 14, 16, 18, 20, 22, in any dimension (an exact polynomial certificate for each n);
  - proved when A or B has two distinct eigenvalues, in any dimension;
  - proved for all 2 × 2 pairs (already in Furuichi–Kuriyama–Yanagi, Linear Algebra Appl. 430 (2009));
  - open in general; the open case with the fewest letters is (5,5), exactly where the upper half fails. See
    the coverage map `lower-half/coverage.svg` in the OQP-40 repository.

Status: CHECKED (internal independent check with exact arithmetic); public. The gap identity and the pinching inequality for all n, m are also machine-checked in Lean in that repository (`lean/OQP40/Pinching.lean`), together with the negative answer to OQP 40 as posed (`lean/OQP40/Statement.lean`).

## 2. CGLMP inequalities with POVMs (IQOQI Open Quantum Problem 27)

**Statement.** On the maximally entangled state $\Phi_D$, no strategy with two d-outcome POVMs per party exceeds the
maximal-entanglement CGLMP value, for every d ≥ 2 and every D. Equality holds only for projective measurements, which
are then the DKZ measurements up to local unitaries.

**How Theorem A is used.**
- A "flattening" B̃ = B + (1 − Π)/4 of the measurement projections removes a spread term.
- Theorem A for the three-level matrix $\tilde P=1-\tilde B$, with spectrum in {0, 3/4, 1}, then gives a partial
  continuum inequality, from which the bound follows.
- The Lean statement of Theorem A covers this spectrum.

**Status:** ANNOUNCED. The proof notes are not yet public; an internal adversarial check found no errors. The proof
is conditional on:
- Theorem A, which is LEAN for 0 ≤ P ≤ 1;
- Theorems B and E of our OQP 27 paper, where Theorem E is computer-assisted for d ≥ 21;
- the internal reduction R1 to the POVM clock model.

The notes will be released with the OQP 27 update (github.com/anshM123/IQOQI-OQP-27).

## 3. Chained Bell inequalities (Barrett–Kent–Pironio)

Full write-up: [`05-chained-bell.md`](05-chained-bell.md); scripts, certificates and logs: [`../bell/`](../bell/README.md).

**Statement.** For the chained Bell expression $I_N$ of Barrett, Kent and Pironio (PRL 97, 170409 (2006)) with N
settings and d outcomes, every projective strategy on a maximally entangled state $\Phi_D$ (any D ≥ 1) satisfies
$I_N\ge I_N^{\mathrm{BKP}}(d)$, the value of the BKP measurements on $\Phi_d$. For d ≥ 3, equality holds iff d divides
D and the strategy is the BKP strategy tensored with $1_{D/d}$, up to local unitaries.
- Proved, computer-assisted (interval arithmetic), for 3 ≤ N ≤ 8 with 3 ≤ d ≤ 10, and for N ∈ {3, 4, 5} with
  11 ≤ d ≤ 24. The case d = 2 is Wehner's bound, and the case N = 2 is CGLMP (our OQP 27 paper).
- This confirms, for projective measurements, the maximally-entangled values that Tavakoli, Zohren and Pawłowski
  derived under their conjecture that the BKP-type measurements are optimal. It does not prove that conjecture, and
  it does not settle the maximal quantum value over all states, which is strictly smaller.

**How Theorem B is used.** A twirl reduces every strategy to a clock model. The chained kernel is a positive
combination of adjacent block pairings, by concavity of k ↦ sin(πk/N). Each pairing is bounded by the continuum
theorem of our OQP 27 paper, which rests on Theorem B for projections.

**Status:** CHECKED (an internal independent check, not external refereeing; it re-verified the certificates with
independent ball-arithmetic code). The proof depends on our unrefereed OQP 27 results
(github.com/anshM123/IQOQI-OQP-27 and github.com/anshM123/CGLMP).

## 4. Operator algebras (in progress)

Theorem A extends to tracial von Neumann algebras with Brown measures in place of eigenvalues. The consequences include
the BMV theorem without any Connes-embedding assumption. This extension has passed an internal check. A standalone
write-up and an external specialist reading are in progress, and the result is not claimed here until they are
complete.
