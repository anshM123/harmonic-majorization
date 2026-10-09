# HarmonicMajorization: Lean log

Authors: Ansh Mishra, Aryan Senthilkumar.  License: MIT.

Mathematics: the working notes THEOREMS.md (published as `../../math/01-theorem-A.md` and
`../../math/02-harmonic-pinching.md`; see `../README.md` for the numbering) and their internal check report.
Toolchain: Lean v4.33.1 with the Mathlib of `../lake-manifest.json`. No OQP27 file is modified; the OQP27 modules are
only imported. This is the development log; to check the proofs, use `../check.sh` (see `../README.md`).

**Status (2026-10-08).** Theorems A, B (left edge, right edge, equality case) and C of THEOREMS.md are formalised for
general Hermitian `0 ≤ B ≤ 1`, together with the multi-line Radon identity that they rest on.  There is no `sorry`,
no `admit`, no new axiom and no `native_decide`.  `#print axioms` gives `[propext, Classical.choice, Quot.sound]` for
every result listed in Section 4 (`logs/AxiomsProj.log`, `logs/AxiomsGen.log`, `logs/AxiomsRI.log`).

## 1. Scope

### Existing formal results used (projection case, `B * B = B`, `P = 1 - B`)

| Lean name | Content |
|---|---|
| `OQP27.StripL3b.stripF P g s τ` | density `ρ(s,τ) = (1/2π) Σ |Im x_i(s,τ)|`, `x_i` roots of `det(g - s - x(P - τ))`, `0 < τ < 1` |
| `OQP27.StripL3b.hyp_RI M` | Radon identity (projection `P`): `(1/2π)∫_0^1 Σ|Im y_i(τ)| dτ = Tr H_+ - Tr (H_d)_+` |
| `OQP27.StripL3b.theorem1_package` | Theorem A (= paper Thm D) for a projection |
| `OQP27.StripL3a.strip_defect_formula_noHyp` | `Σ_ν h_λ(ν) - Tr[(P(g-λ)P)_+] = V(λ) = ∫_0^1∫ K_{1-τ}(s-λ) ρ ds dτ` |
| `OQP27.StripL3a.stripKernel`, `stripBalayage`, `hLam`, `posPartTrace` | `K_X`, `V`, `h_λ`, `Tr (·)_+` |
| `OQP27/StripRI*.lean` | proof of `hyp_RI`; its generic contour, Jensen and Burgers tools are reused here |

### Plan (as executed)

1. Theorem C for a projection: the τ-marginal `∫∫ e^{-tτ} ρ = ‖BgP‖_F² (1 - e^{-t})/t` from Duhamel's formula,
   Theorem A at `a ≠ 0` and continuity at `a = 0`; `t = ±iπ` gives the sine moment; kernel lower bound.
2. General `B`: spectral calculus, pencil lemmas (Lemma 2.2 incl. (6)), Theorems A, B, C from the hypothesis
   `Hyp_RI_gen` (the multi-line Radon identity).
3. Proof of `Hyp_RI_gen` (files `RIG*.lean`), then the unconditional theorems.

## 2. Build (historical; use `../check.sh`)

In our development tree, one file at a time (`build.sh` runs the given files sequentially, logs in `logs/<F>.log`):

    export PATH="$HOME/.elan/bin:$PATH"
    mkdir -p .lake/build/lib/lean/HarmonicMajorization
    lake env lean -o .lake/build/lib/lean/HarmonicMajorization/<F>.olean \
                  -i .lake/build/lib/lean/HarmonicMajorization/<F>.ilean HarmonicMajorization/<F>.lean

Full rebuild, in this order (one Lean process at a time; every file builds without errors or warnings), then the
axiom checks (62 declarations, all
`[propext, Classical.choice, Quot.sound]`):

    bash HarmonicMajorization/build.sh Duhamel MarginalProj TheoremCProj SpecDecomp PencilGen TheoremAGen MarginalGen PoissonGen TheoremBGen TheoremCGen RIG1 RIG2 RIG3a RIG3b RIG4a RIG4b MainGen RightEdge DensityGen SliceGen MassGen
    lake env lean HarmonicMajorization/AxiomsProj.lean
    lake env lean HarmonicMajorization/AxiomsGen.lean
    lake env lean HarmonicMajorization/AxiomsRI.lean

The OQP27 oleans (`StripRIMain`, `StripRIAssembly`, `StripRIBoundary`, `StripRIEnds`, ...) must already be in
`.lake/build/lib/lean/OQP27/`; `CertD18-20` are not needed.

## 3. Files

Projection case:

| File | Content |
|---|---|
| `Duhamel.lean` | `trace_duhamel`, `trace_exp_add_smul` (`Tr e^{X+aY} = Tr e^X + a Tr(Y e^X) + a² Φ(a)`), `continuous_phi2` |
| `MarginalProj.lean` | `laplaceF_zero` (τ-marginal, Laplace form), `integral_sin_mul_stripF` |
| `TheoremCProj.lean` | `kernel_lower_bound`, `stripBalayage_ge_proj`, `frobNorm_comm_sq`, **`theoremC_proj`** |
| `AxiomsProj.lean` | `#print axioms` |

General `B`, conditional layer (everything that uses the multi-line RI takes it as the hypothesis `Hyp_RI_gen M`):

| File | Content |
|---|---|
| `SpecDecomp.lean` | functional calculus `fcalc`, `spec`, spectral projections `sproj`, pinching `pinchH` |
| `PencilGen.lean` | Lemma 2.2 incl. **(6)** `normSq_nonreal_root_le_gap`; density `rhoG`, `measurable_rhoG`, `rhoG_le_domG` |
| `TheoremAGen.lean` | `Hyp_RI_gen`, Radon slices, `theoremA_gen` |
| `MarginalGen.lean` | `laplaceG_zero`, `integral_sin_mul_rhoG` |
| `PoissonGen.lean` | `DensityRegG`, balayage lemmas (continuity of `V`), two-sweep Fourier core |
| `TheoremBGen.lean` | `theoremB_gen`, `theoremB_gen_ineq`, `theoremB_gen_eq_iff`, `mass_moment_gen` (no hypothesis) |
| `TheoremCGen.lean` | `theoremC_gen`, `stripBalayage_ge_gen` |
| `AxiomsGen.lean` | `#print axioms` |

The multi-line Radon identity (no hypotheses):

| File | Content |
|---|---|
| `RIG1.lean` | root location and bounds for a Hermitian `P` (Lemma A, B(i), B(ii) on a gap, B(iii) `norm_Sup_sub_le_gap`) |
| `RIG2.lean` | local regularity and the Burgers identity `∂_τ Λ_± = ∂_c S_±` on `(ℝ ∖ spec P) × ℂ₊` |
| `RIG3a.lean` | Lemma A at every `τ`, the **dichotomy** near an eigenvalue, the rescaled pencil `WmatP` |
| `RIG3b.lean` | **`Gup_cancel`**: cancellation of the boundary terms at every eigenvalue (`maxHeartbeats 1000000` for this proof) |
| `RIG4a.lean` | sorted spectrum, Step 1 on the gaps, the two outer regions, telescoping, **`PsiG_sub_E_affine`** |
| `RIG4b.lean` | Lemma E for general `P`, boundary values, **`radon_identity_gen`**, **`hyp_RI_gen`** |

Unconditional theorems:

| File | Content |
|---|---|
| `MainGen.lean` | **`theoremA`**, `tau_marginal`, `sine_moment`, **`theoremB`**, `theoremB_ineq`, `theoremB_eq_iff`, **`theoremC`**, `theoremC_balayage` |
| `RightEdge.lean` | `pinchH_one_sub` (pinching by `B` = by `1 - B`), `roots_charpoly_one_sub_add`, **`theoremB_right`**, `theoremB_right_ineq`, `theoremB_right_eq_iff` |
| `DensityGen.lean` | `rhoG_eq_zero_of_outside`, **`theoremA_package`**, **`total_mass`** |
| `SliceGen.lean` | **`sliceMass_ae_eq`**: the τ-slices of `ρ` integrate to the step function `m(τ)` (a.e.) |
| `MassGen.lean` | `frobNorm_sub_pinchH_sq` (`‖g - g_d‖_F² = 2 Σ_{p<q} ‖Q_p g Q_q‖_F²`), **`total_mass_frob`** |
| `AxiomsRI.lean` | `#print axioms` for the multi-line RI and all unconditional theorems |

Route of the proof of the multi-line RI (generalising `OQP27/StripRI*.lean`): `Ψ(c) = Σ_gaps ∫ (S_+(A) - S_+(A_d)) dτ`.
On each gap `d/dc ∫_{α+δ}^{β-δ} = G(β - δ) - G(α + δ)` (Burgers identity), `G = Λ_+(A) - Λ_+(A_d)`.  Summed over the
gaps this telescopes into `Σ_{p ∈ spec P} (G(p - δ) - G(p + δ)) - G(p_min - δ) + G(p_max + δ)`; `G(p_max + δ) = 0`,
`G(p_min - δ) = Θ₁ + 2πik` (all roots lie in `ℂ₊` left of the spectrum), and every `G(p - δ) - G(p + δ) → 0`
(`Gup_cancel`: near `τ = p` the roots split into bounded ones, whose log-sum is a contour integral continuous across
`τ = p`, and large ones `w/(p - τ)`, with `w` the roots of the rescaled pencil `WmatP`, whose characteristic
polynomials for `A` and `A_d` agree at `τ = p`).  The mean value inequality gives `Ψ - E` affine; Steps 3-4 are as in
the projection case.  Deviation from [27, Remark RIgeneralB]: the end `p_min` is handled through the outer region
`τ < p_min`, and one cancellation lemma covers interior and extreme eigenvalues.

## 4. Main statements (namespace `HarmonicMajorization`)

Notation: `spec hP` = finset of distinct eigenvalues of `P`; `sproj hP p` = spectral projection `Q_p`;
`pinchH hP g = Σ_p Q_p g Q_p` (`g_d`); `rhoG P g s τ = imAbsSum (pencilRoots P (g - s•1) τ) / (2π)` (`ρ`);
`hLam` = `h_λ` and `stripBalayage` = `V` (from OQP27); `frobNorm`, `opNorm` = Frobenius and operator norms.

    def Hyp_RI_gen (M : ℕ) : Prop :=
      ∀ (P H : Matrix (Fin M) (Fin M) ℂ) (hP : P.IsHermitian) (hH : H.IsHermitian),
        Integrable (fun τ => imAbsSum (pencilRoots P H τ)) ∧
        (∫ τ, imAbsSum (pencilRoots P H τ)) / (2 * π)
          = ∑ i, max (hH.eigenvalues i) 0 - ∑ i, max ((isHermitian_pinchH hP hH).eigenvalues i) 0
    theorem hyp_RI_gen (M : ℕ) : Hyp_RI_gen M

    -- (hP : P.IsHermitian) (hg : g.IsHermitian) (hspec : ∀ p ∈ spec hP, 0 ≤ p ∧ p ≤ 1)
    theorem theoremA (a t : ℂ) : bmvDG hP g a t = a ^ 2 * laplaceG P g a t
    --   bmvDG hP g a t = (exp (a • g - t • P)).trace - (exp (a • pinchH hP g - t • P)).trace
    --   laplaceG P g a t = ∫ q : ℝ × ℝ, cexp (a * q.1 - t * q.2) * (rhoG P g q.1 q.2 : ℂ)
    theorem total_mass : ∫ q : ℝ × ℝ, rhoG P g q.1 q.2
        = ∑ p ∈ spec hP, ∑ q ∈ spec hP, if p < q then frobNorm (sproj hP p * g * sproj hP q) ^ 2 else 0
    theorem total_mass_frob : ∫ q : ℝ × ℝ, rhoG P g q.1 q.2 = frobNorm (g - pinchH hP g) ^ 2 / 2
    theorem sliceMass_ae_eq : sliceMass P g =ᵐ[volume] stepMass hP g

    -- (hB0 : B.PosSemidef) (hB1 : (1 - B).PosSemidef) (hg : g.IsHermitian)
    theorem theoremB (lam : ℝ) :
        ((B + I • g).charpoly.roots.map (hLam lam)).sum
          - ((B + I • pinchH hB1.isHermitian g).charpoly.roots.map (hLam lam)).sum
          = stripBalayage (rhoG (1 - B) g) lam
    theorem theoremB_eq_iff (lam : ℝ) :
        ((B + I • g).charpoly.roots.map (hLam lam)).sum
            = ((B + I • pinchH hB1.isHermitian g).charpoly.roots.map (hLam lam)).sum ↔ Commute B g
    theorem theoremB_right (lam : ℝ) :
        ((B + I • g).charpoly.roots.map (fun ν => hLam lam (1 - star ν))).sum
          - ((B + I • pinchH hB1.isHermitian g).charpoly.roots.map (fun ν => hLam lam (1 - star ν))).sum
          = stripBalayage (rhoG B g) lam
    theorem theoremC (lam : ℝ) :
        1 / (2 * (1 + Real.cosh (π * (|lam| + opNorm g))))
            * (∑ p ∈ spec hB1.isHermitian, ∑ q ∈ spec hB1.isHermitian, if p < q then
                frobNorm (sproj hB1.isHermitian p * g * sproj hB1.isHermitian q) ^ 2
                  * (1 / (q - p) * ∫ τ in p..q, Real.sin (π * τ)) else 0)
          ≤ ((B + I • g).charpoly.roots.map (hLam lam)).sum
            - ((B + I • pinchH hB1.isHermitian g).charpoly.roots.map (hLam lam)).sum

    -- projection case (hB : B.IsHermitian ∧ B * B = B) (hg : g.IsHermitian)
    theorem theoremC_proj (lam : ℝ) :
        frobNorm (B * g - g * B) ^ 2 / (2 * π * (1 + Real.cosh (π * (|lam| + opNorm g))))
          ≤ ((B + I • g).charpoly.roots.map (OQP27.StripL3a.hLam lam)).sum
            - OQP27.StripL3a.posPartTrace ((1 - B) * (g - (lam : ℂ) • 1) * (1 - B))

`theoremB_ineq`, `theoremB_right_ineq` (the defect is `≥ 0`), `theoremB_right_eq_iff`, `theoremC_balayage`
(`V(λ) ≥` the same bound), `tau_marginal`, `sine_moment`, `theoremA_package` (measurability, support, integrable
profile, integrability, identity) and `radon_identity_gen` are stated in the files above.

## 5. Not formalised

* The convex-order corollary of Theorem B (`∫ f dσ ≤ ∫ f dω` for convex `f` of linear growth).  Lean has the
  `(w - λ)_+` inequalities for every `λ` and the equal masses and first moments (`mass_moment_gen`).
* Continuity of `ρ` off the lines `τ ∈ spec P`, and the slice formula `∫ ρ(s, τ) ds = m(τ)` for *every*
  `τ ∉ spec P`.  It is proved for almost every `τ` (`sliceMass_ae_eq`) and in Laplace form (`tau_marginal`).
* Corollary 2 of THEOREMS.md (BMV with the split into atoms and a density) beyond Theorem A itself.
