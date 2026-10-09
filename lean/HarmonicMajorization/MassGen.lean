import HarmonicMajorization.DensityGen

/-!
# The total mass of the density as a Frobenius norm

Authors: Ansh Mishra, Aryan Senthilkumar.  License: MIT.

`P` Hermitian, `g` Hermitian, `g_d = pinchH hP g`.
Main results:
* `frobNorm_sub_pinchH_sq`: `‖g - g_d‖_F² = 2 Σ_{p<q} ‖Q_p g Q_q‖_F²`;
* `total_mass_frob`: `∫∫ ρ = ‖g - g_d‖_F² / 2` (THEOREMS.md, Theorem A, properties of `ρ`).
-/

set_option autoImplicit false

namespace HarmonicMajorization

open Matrix Complex MeasureTheory
open scoped Real ComplexOrder

variable {M : ℕ}

/-- `‖Z‖_F² = Tr(Z Zᴴ)` (as a complex number). -/
lemma frobNorm_sq_eq_trace (Z : Matrix (Fin M) (Fin M) ℂ) :
    ((frobNorm Z ^ 2 : ℝ) : ℂ) = (Z * Zᴴ).trace := by
  rw [frobNorm_sq]
  simp only [Matrix.trace, Matrix.diag, Matrix.mul_apply, Matrix.conjTranspose_apply]
  push_cast
  refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
  rw [Complex.star_def, Complex.mul_conj']

section Mass

variable {P g : Matrix (Fin M) (Fin M) ℂ} (hP : P.IsHermitian) (hg : g.IsHermitian)
include hP

lemma trace_mul_pinchH_sum :
    (g * pinchH hP g).trace = ∑ p ∈ spec hP, (g * sproj hP p * g * sproj hP p).trace := by
  unfold pinchH
  rw [Matrix.mul_sum, Matrix.trace_sum]
  refine Finset.sum_congr rfl fun p _ => ?_
  simp only [Matrix.mul_assoc]

lemma trace_pinchH_mul_pinchH_sum :
    (pinchH hP g * pinchH hP g).trace = ∑ p ∈ spec hP, (g * sproj hP p * g * sproj hP p).trace := by
  have e : pinchH hP g * pinchH hP g
      = (∑ p ∈ spec hP, sproj hP p * g * sproj hP p) * pinchH hP g := rfl
  rw [e, Finset.sum_mul, Matrix.trace_sum]
  refine Finset.sum_congr rfl fun p hp => ?_
  rw [Matrix.mul_assoc (sproj hP p * g), sproj_mul_pinchH hP g hp]
  have e2 : sproj hP p * g * (sproj hP p * g * sproj hP p)
      = sproj hP p * (g * sproj hP p * g * sproj hP p) := by
    simp only [Matrix.mul_assoc]
  rw [e2, Matrix.trace_mul_comm (sproj hP p) (g * sproj hP p * g * sproj hP p),
    show g * sproj hP p * g * sproj hP p * sproj hP p = g * sproj hP p * g * sproj hP p by
      rw [Matrix.mul_assoc _ (sproj hP p) (sproj hP p), sproj_mul_self]]

lemma trace_mul_self_sub_pinchH :
    (g * g).trace - (pinchH hP g * pinchH hP g).trace
      = ∑ p ∈ spec hP, ∑ q ∈ spec hP,
          if p = q then 0 else (g * sproj hP p * g * sproj hP q).trace := by
  have h := trace_sandwich_diff hP g (fun _ => (1 : ℂ)) (fun _ => (1 : ℂ))
  have h1 : fcalc hP (fun _ => (1 : ℂ)) = 1 := by rw [fcalc_const, one_smul]
  simp only [h1, Matrix.mul_one, one_mul] at h
  exact h

include hg

/-- **`‖g - g_d‖_F² = 2 Σ_{p<q} ‖Q_p g Q_q‖_F²`.** -/
theorem frobNorm_sub_pinchH_sq :
    frobNorm (g - pinchH hP g) ^ 2
      = 2 * ∑ p ∈ spec hP, ∑ q ∈ spec hP, if p < q then
          frobNorm (sproj hP p * g * sproj hP q) ^ 2 else 0 := by
  have hgd := isHermitian_pinchH hP hg
  have hX : (g - pinchH hP g)ᴴ = g - pinchH hP g := (hg.sub hgd).eq
  have htr : (g * pinchH hP g).trace = (pinchH hP g * pinchH hP g).trace := by
    rw [trace_mul_pinchH_sum hP, trace_pinchH_mul_pinchH_sum hP]
  have htr' : (pinchH hP g * g).trace = (pinchH hP g * pinchH hP g).trace := by
    rw [Matrix.trace_mul_comm, htr]
  have hsq : (((g - pinchH hP g) * (g - pinchH hP g)ᴴ).trace)
      = (g * g).trace - (pinchH hP g * pinchH hP g).trace := by
    rw [hX, Matrix.sub_mul, Matrix.mul_sub, Matrix.mul_sub, Matrix.trace_sub, Matrix.trace_sub,
      Matrix.trace_sub, htr, htr']
    ring
  have hsym : ∀ p q : ℝ, (g * sproj hP p * g * sproj hP q).trace
      = ((frobNorm (sproj hP p * g * sproj hP q) ^ 2 : ℝ) : ℂ) := by
    intro p q
    rw [trace_sandwich_eq_frob hP hg, frobNorm_sandwich_sq hP hg]
  have key : ((frobNorm (g - pinchH hP g) ^ 2 : ℝ) : ℂ)
      = ((2 * ∑ p ∈ spec hP, ∑ q ∈ spec hP, if p < q then
          frobNorm (sproj hP p * g * sproj hP q) ^ 2 else 0 : ℝ) : ℂ) := by
    rw [frobNorm_sq_eq_trace, hsq, trace_mul_self_sub_pinchH hP,
      sum_offdiag_symm (spec hP) (fun p q => (g * sproj hP p * g * sproj hP q).trace)]
    push_cast
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun p _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun q _ => ?_
    split_ifs with hpq
    · rw [trace_sandwich_symm hP q p, hsym p q]
      ring
    · simp
  exact_mod_cast key

variable (hspec : ∀ p ∈ spec hP, 0 ≤ p ∧ p ≤ 1)
include hspec

/-- **Total mass of the density**: `∫∫ ρ = ‖g - g_d‖_F² / 2`. -/
theorem total_mass_frob :
    ∫ q : ℝ × ℝ, rhoG P g q.1 q.2 = frobNorm (g - pinchH hP g) ^ 2 / 2 := by
  rw [total_mass hP hg hspec, frobNorm_sub_pinchH_sq hP hg]
  ring

end Mass

end HarmonicMajorization
