# Harmonic majorization

**Authors:** Ansh Mishra, Aryan Senthilkumar. **License:** MIT. **Started:** 2026-10-08.

This repository develops a theory of how non-commutativity moves spectra and traces, starting from one explicit
identity.

## The core idea

Take Hermitian matrices B and g, and let $g_d$ be the **pinching** of g onto the eigenspaces of B (the part of g that
commutes with B). The traces $\operatorname{Tr}e^{ag-tP}$, and hence all symmetrized word averages of (B, g), differ
from those of the commuting pair (B, g_d) by a **second derivative of one explicit positive density**:

$$\rho(s,\tau)=\frac1{2\pi}\sum_i|\operatorname{Im}\xi_i(s,\tau)|,\qquad\text{where }\xi_i\text{ are the roots of }\det\bigl(g-s-\xi(P-\tau)\bigr),\quad P=1-B .$$

The commuting case is an identity. Non-commutativity always pushes in one direction, by an explicit amount. Different
test functions turn this into theorems in different fields:

| Test function | Result | Field |
|---|---|---|
| exponentials $e^{as-t\tau}$ | Theorem A: the pinched two-variable BMV identity | matrix analysis |
| monomials $\tau^ns^m$ | exact gap formula for word averages; the pinching inequality conjectured by Dinh (2026), all word lengths | trace inequalities |
| strip-harmonic functions | Theorem B: the eigenvalues of the **non-normal** matrix B + ig, swept to the boundary by harmonic measure, majorize those of B + ig_d | spectra of non-normal matrices |
| kernel lower bound | Theorem C: a quantitative version with a numerically sharp constant | stability / rigidity |
| superharmonic functions (two levels) | superharmonic majorization, which contains the two-level cases of Theorem B and of Ky Fan's theorem | potential theory |

## Status of the main results

Labels:
- **LEAN**: machine-checked in Lean 4 with Mathlib, standard axioms only.
- **CHECKED**: written proof, plus an internal independent adversarial check (not external refereeing).
- **PROVED**: written proof, not independently checked.

| Result | Status | Where |
|---|---|---|
| **Theorem A** (pinched two-variable BMV identity, general multiplicities) | LEAN for 0 ≤ P ≤ 1: the identity, nonnegativity, support, mass and τ-marginal (slices for almost every τ). On paper: continuity of ρ, the slice formula at every τ, and general P (affine rescaling). CHECKED | [`math/01-theorem-A.md`](math/01-theorem-A.md), [`lean/`](lean/) |
| **Theorem B** (harmonic pinching law, both edges, equality iff [B,g] = 0) | LEAN for the defect identity and the equality case; the convex-order corollary is CHECKED on paper | [`math/02-harmonic-pinching.md`](math/02-harmonic-pinching.md) |
| **Theorem C** (quantitative bound) | LEAN; sharpness of the constant is numerical | [`math/02-harmonic-pinching.md`](math/02-harmonic-pinching.md) |
| Holomorphic form, exact admissibility, half-plane version, two-level superharmonic majorization | CHECKED (internal pre-publication check, [`checks/structure/`](checks/structure/)) | [`math/03-structure.md`](math/03-structure.md) |
| Word averages: exact gap formula, corrected two-sided bounds; Dinh's conjecture for all word lengths | gap formula and Dinh's inequality: LEAN for all positive semidefinite A, B (the density is taken for A rescaled into [0, 1]), in [anshM123/OQP-40](https://github.com/anshM123/OQP-40) (`lean/`); the bounds: CHECKED. The inequality itself also follows from Heinävaara's Corollary 1.2 (arXiv:2310.03227v1) | [`math/04-applications.md`](math/04-applications.md) |
| CGLMP / OQP 27 for all POVMs | ANNOUNCED: the proof notes are not yet public; an internal adversarial check found no errors | [`math/04-applications.md`](math/04-applications.md) |

## Verify the Lean proofs

Requirements: [elan](https://github.com/leanprover/elan); the toolchain is pinned in `lean/lean-toolchain`.

```bash
cd lean
lake exe cache get        # prebuilt Mathlib at the pinned revision
bash check.sh             # elaborates every file in BUILD_ORDER.txt from source, then prints the axioms (about 15 min)
```

Main theorems in namespace `HarmonicMajorization`: `theoremA`, `tau_marginal`, `total_mass`, `theoremB`,
`theoremB_right`, `theoremB_eq_iff`, `theoremC`, `theoremC_proj`, `hyp_RI_gen` (the multi-line Radon identity). Each
depends only on `[propext, Classical.choice, Quot.sound]`.

## Credits

- **Stahl's theorem** (the former BMV conjecture): H. R. Stahl, Acta Math. 211 (2013). Its atoms/density structure is
  the starting point of Theorem A.
- **The explicit |Im|-density of the BMV measure** (simple spectrum): O. Heinävaara, PhD thesis (Princeton 2024).
  **Tracial joint spectral measures:** O. Heinävaara, arXiv:2310.03227 (Invent. Math. 239 (2025) 505–526).
- **High-frequency averaging:** a special case of D. Burgarth, P. Facchi, H. Nakazato, S. Pascazio, K. Yuasa,
  Quantum 3, 152 (2019).
- **The projection case** of Theorems A and B: our OQP 27 paper (github.com/anshM123/IQOQI-OQP-27).
- **Z. Long, arXiv:2610.10158** proves a spectral-defect identity comparing the eigenvalues of T = A + iB with
  Heinävaara's measure of (Re T, Im T), with a nonnegative superharmonic remainder (his Lemma 2.3 and Proposition
  2.4). Theorem D_sep of [`math/03-structure.md`](math/03-structure.md) follows from it directly. Our two-level
  identity has the same two-remainder structure, but compares with the pinched spectrum.
- **Ky Fan's majorization** Re spec X ≺ spec(Re X): K. Fan (1950); see the survey arXiv:1108.1467. The converse is due to
  Amir-Moéz and Horn, and to Mirsky.

## Citation

See [`CITATION.cff`](CITATION.cff).
