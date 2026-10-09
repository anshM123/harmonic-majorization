# Lean 4 formalization

Lean `v4.33.1` with Mathlib at the pinned revision (`lake-manifest.json`). Two libraries:
- `OQP27`: the 23 files of our OQP 27 formalization (github.com/anshM123/IQOQI-OQP-27) that the theory imports. The
  projection case of Theorems A and B and the projection Radon identity live there.
- `HarmonicMajorization`: Theorems A, B, C for general 0 ≤ B ≤ 1, and the multi-line Radon identity.

## Check

```bash
lake exe cache get    # prebuilt Mathlib
bash check.sh         # every file in BUILD_ORDER.txt, then #print axioms; about 15 minutes
```

## Main declarations (namespace `HarmonicMajorization`)

| Declaration | Statement |
|---|---|
| `theoremA` | for Hermitian P with spectrum in [0,1] and Hermitian g, for all a, t ∈ ℂ: `Tr e^{ag−tP} − Tr e^{a g_d − tP} = a² ∫∫ e^{as−tτ} ρ(s,τ)` |
| `theoremA_package`, `total_mass`, `total_mass_frob`, `tau_marginal`, `sliceMass_ae_eq` | properties of the density ρ (nonnegativity, support, integrability, mass ‖g − g_d‖²/2, τ-marginal; slices for almost every τ) |
| `hyp_RI_gen`, `radon_identity_gen` | the multi-line Radon identity |
| `theoremB`, `theoremB_right` | left- and right-edge defect identities `Σ_{spec(B+ig)} h_λ − Σ_{spec(B+ig_d)} h_λ = V(λ) ≥ 0` |
| `theoremB_eq_iff`, `theoremB_right_eq_iff` | equality at one λ iff B and g commute |
| `theoremC`, `theoremC_balayage`, `theoremC_proj` | the quantitative lower bound |

All of them depend only on `[propext, Classical.choice, Quot.sound]`. In `logs/`:
- `cleanbuild.log`: a rebuild of all 44 files from source into an empty output directory, with no previously built
  file of this project on the search path (2026-10-08, all OK, 644 s);
- `check.log`: a run of `check.sh` on a copy of this folder (all 47 steps OK; 62 axiom lines, all standard);
- `AxiomsProj.log`, `AxiomsRI.log`, `AxiomsGen.log`: the full `#print axioms` output against the clean rebuild (62
  declarations, none with a nonstandard axiom).

## Numbering and internal references

Some docstrings cite our internal working notes, which are not part of this repository (THEOREMS.md, CHECK_REPORT.md,
LEAN_BRIEF.md, SHARED_LEMMAS.md, and paths under `iqoqi/programs/` and `formal-conjectures-27/`). Where a docstring
says "no hypotheses", it means no unproved hypotheses. In particular:
- the `OQP27` files cite paths under `iqoqi/programs/oqp27B_all/`; the corresponding public statements and proofs are in
  github.com/anshM123/IQOQI-OQP-27 (`papers/math/main.tex`);
- the `HarmonicMajorization` files cite `THEOREMS.md`. The theorem names (A, B, C) are the same in the public notes,
  but the lemma numbers differ:

| THEOREMS.md | public notes |
|---|---|
| Lemma 2.1 (strip kernel) | `math/02-harmonic-pinching.md`, Lemma 1 |
| Lemma 2.2 (pencil), incl. (6) | `math/01-theorem-A.md`, Lemma 1 (6) |
| Lemma 2.3 (Jensen identity) | `math/01-theorem-A.md`, Lemma 2 |
| Lemma 2.4 (real slices) | `math/01-theorem-A.md`, Lemma 3 |
| Lemma 2.5 (high-frequency averaging) | `math/01-theorem-A.md`, Lemma 4 |
| Lemma 2.6 (exponential polynomials) | `math/01-theorem-A.md`, Lemma 5 |

The Lean proof of Theorem A goes through the multi-line Radon identity, not through the distributional proof of the
notes. General Hermitian P (not in [0, 1]) follows by affine rescaling, which is done on paper only. The convex-order
corollary of Theorem B is likewise on paper only.
