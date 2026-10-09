import HarmonicMajorization.DensityGen
import Mathlib.MeasureTheory.Measure.CharacteristicFunction.Basic

/-!
# The `τ`-slices of the density integrate to a step function

Authors: Ansh Mishra, Aryan Senthilkumar.  License: MIT.

`P` Hermitian with spectrum in `[0, 1]`, `g` Hermitian, `ρ = rhoG P g`.
THEOREMS.md, Theorem A, properties of `ρ`: `∫_ℝ ρ(s, τ) ds = m(τ)` with
`m(τ) = Σ_{p<q} ‖Q_p g Q_q‖_F² / (q - p) · 1[p < τ < q]`.
Main result:
* `sliceMass_ae_eq`: `τ ↦ ∫ ρ(s, τ) ds` equals `m` for almost every `τ`.
Route: both sides are integrable and nonnegative; their Fourier transforms agree (the `τ`-marginal in
Laplace form at `t = -iω`, and the total mass at `ω = 0`); the characteristic function determines a finite
measure (`Measure.ext_of_charFun`).  (The statement for every `τ ∉ spec P`, which needs the continuity of
`ρ` off the lines `τ ∈ spec P`, is not formalised.)
-/

set_option autoImplicit false

namespace HarmonicMajorization

open Matrix Complex MeasureTheory Filter Topology Set
open scoped Real ComplexOrder ENNReal

variable {M : ℕ}

/-! ### Uniqueness of densities with the same Fourier transform -/

lemma charFun_withDensity_ofReal {f : ℝ → ℝ} (hf : Integrable f) (hf0 : ∀ x, 0 ≤ f x) (ω : ℝ) :
    charFun (volume.withDensity fun x => ENNReal.ofReal (f x)) ω
      = ∫ x : ℝ, cexp ((ω : ℂ) * (x : ℂ) * I) * (f x : ℂ) := by
  rw [charFun_apply_real, integral_withDensity_eq_integral_toReal_smul₀
    hf.1.aemeasurable.ennreal_ofReal (ae_of_all _ fun x => ENNReal.ofReal_lt_top)]
  congr 1
  funext x
  rw [ENNReal.toReal_ofReal (hf0 x), Complex.real_smul, mul_comm]

/-- Two nonnegative integrable functions with the same Fourier transform agree almost everywhere. -/
lemma ae_eq_of_fourier_eq {f₁ f₂ : ℝ → ℝ} (h1 : Integrable f₁) (h2 : Integrable f₂)
    (h10 : ∀ x, 0 ≤ f₁ x) (h20 : ∀ x, 0 ≤ f₂ x)
    (h : ∀ ω : ℝ, ∫ x : ℝ, cexp ((ω : ℂ) * (x : ℂ) * I) * (f₁ x : ℂ)
      = ∫ x : ℝ, cexp ((ω : ℂ) * (x : ℂ) * I) * (f₂ x : ℂ)) :
    f₁ =ᵐ[volume] f₂ := by
  have := isFiniteMeasure_withDensity_ofReal h1.2
  have := isFiniteMeasure_withDensity_ofReal h2.2
  have heq : volume.withDensity (fun x => ENNReal.ofReal (f₁ x))
      = volume.withDensity (fun x => ENNReal.ofReal (f₂ x)) := by
    apply Measure.ext_of_charFun
    funext ω
    rw [charFun_withDensity_ofReal h1 h10, charFun_withDensity_ofReal h2 h20, h ω]
  rw [withDensity_eq_iff_of_sigmaFinite h1.1.aemeasurable.ennreal_ofReal
    h2.1.aemeasurable.ennreal_ofReal] at heq
  filter_upwards [heq] with x hx
  exact (ENNReal.ofReal_eq_ofReal_iff (h10 x) (h20 x)).mp hx

/-! ### The slice masses -/

/-- `μ(τ) = ∫ ρ(s, τ) ds`. -/
noncomputable def sliceMass (P g : Matrix (Fin M) (Fin M) ℂ) (τ : ℝ) : ℝ := ∫ s, rhoG P g s τ

/-- `m(τ) = Σ_{p<q} ‖Q_p g Q_q‖_F² / (q - p) · 1[p < τ < q]`. -/
noncomputable def stepMass {P : Matrix (Fin M) (Fin M) ℂ} (hP : P.IsHermitian)
    (g : Matrix (Fin M) (Fin M) ℂ) (τ : ℝ) : ℝ :=
  ∑ p ∈ spec hP, ∑ q ∈ spec hP, if p < q then
    frobNorm (sproj hP p * g * sproj hP q) ^ 2 / (q - p) * (Ioo p q).indicator 1 τ else 0

lemma sliceMass_nonneg (P g : Matrix (Fin M) (Fin M) ℂ) (τ : ℝ) : 0 ≤ sliceMass P g τ :=
  integral_nonneg fun s => rhoG_nonneg P g s τ

lemma stepMass_nonneg {P : Matrix (Fin M) (Fin M) ℂ} (hP : P.IsHermitian)
    (g : Matrix (Fin M) (Fin M) ℂ) (τ : ℝ) : 0 ≤ stepMass hP g τ := by
  unfold stepMass
  refine Finset.sum_nonneg fun p _ => Finset.sum_nonneg fun q _ => ?_
  split_ifs with hpq
  · have h1 : 0 ≤ (Ioo p q).indicator (1 : ℝ → ℝ) τ :=
      Set.indicator_nonneg (fun _ _ => zero_le_one) τ
    have h2 : 0 < q - p := by linarith
    positivity
  · exact le_rfl

lemma integrable_indicator_Ioo_one (p q : ℝ) : Integrable ((Ioo p q).indicator (1 : ℝ → ℝ)) :=
  (integrable_indicator_iff measurableSet_Ioo).mpr (integrableOn_const (by simp))

lemma integrable_stepMass {P : Matrix (Fin M) (Fin M) ℂ} (hP : P.IsHermitian)
    (g : Matrix (Fin M) (Fin M) ℂ) : Integrable (stepMass hP g) := by
  unfold stepMass
  refine integrable_finsetSum _ fun p _ => integrable_finsetSum _ fun q _ => ?_
  split_ifs
  · exact (integrable_indicator_Ioo_one p q).const_mul _
  · exact integrable_zero _ _ _

lemma integral_stepMass {P : Matrix (Fin M) (Fin M) ℂ} (hP : P.IsHermitian)
    (g : Matrix (Fin M) (Fin M) ℂ) :
    ∫ τ, stepMass hP g τ = ∑ p ∈ spec hP, ∑ q ∈ spec hP, if p < q then
      frobNorm (sproj hP p * g * sproj hP q) ^ 2 else 0 := by
  unfold stepMass
  rw [integral_finsetSum _ fun p _ => integrable_finsetSum _ fun q _ => by
    split_ifs
    · exact (integrable_indicator_Ioo_one p q).const_mul _
    · exact integrable_zero _ _ _]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [integral_finsetSum _ fun q _ => by
    split_ifs
    · exact (integrable_indicator_Ioo_one p q).const_mul _
    · exact integrable_zero _ _ _]
  refine Finset.sum_congr rfl fun q _ => ?_
  split_ifs with hpq
  · rw [integral_const_mul, integral_indicator_one measurableSet_Ioo, Real.volume_real_Ioo_of_le
      hpq.le]
    have h : q - p ≠ 0 := by linarith
    field_simp
  · exact integral_zero _ _

/-- The Fourier transform of `m`. -/
lemma fourier_stepMass {P : Matrix (Fin M) (Fin M) ℂ} (hP : P.IsHermitian)
    (g : Matrix (Fin M) (Fin M) ℂ) {ω : ℝ} (hω : ω ≠ 0) :
    ∫ τ : ℝ, cexp ((ω : ℂ) * (τ : ℂ) * I) * (stepMass hP g τ : ℂ)
      = ∑ p ∈ spec hP, ∑ q ∈ spec hP, if p < q then
          ((frobNorm (sproj hP p * g * sproj hP q) ^ 2 : ℝ) : ℂ)
            * ((cexp (ω * I * q) - cexp (ω * I * p)) / ((ω * I) * ((q : ℂ) - p))) else 0 := by
  have hc : (ω * I : ℂ) ≠ 0 := mul_ne_zero (by exact_mod_cast hω) I_ne_zero
  set T : ℝ → ℝ → ℝ → ℂ := fun p q τ => if p < q then (Ioo p q).indicator (fun τ =>
      cexp (ω * τ * I) * ((frobNorm (sproj hP p * g * sproj hP q) ^ 2 / (q - p) : ℝ) : ℂ)) τ
    else 0 with hT
  have hint : ∀ p q, Integrable (T p q) := by
    intro p q
    simp only [hT]
    split_ifs
    · rw [integrable_indicator_iff measurableSet_Ioo]
      exact (Continuous.integrableOn_Icc (by fun_prop)).mono_set Ioo_subset_Icc_self
    · exact integrable_zero _ _ _
  have hfun : (fun τ : ℝ => cexp (ω * τ * I) * (stepMass hP g τ : ℂ))
      = fun τ => ∑ p ∈ spec hP, ∑ q ∈ spec hP, T p q τ := by
    funext τ
    unfold stepMass
    rw [Complex.ofReal_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun p _ => ?_
    rw [Complex.ofReal_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun q _ => ?_
    simp only [hT]
    split_ifs with hpq
    · by_cases hτ : τ ∈ Ioo p q
      · rw [Set.indicator_of_mem hτ, Set.indicator_of_mem hτ, Pi.one_apply, mul_one]
      · rw [Set.indicator_of_notMem hτ, Set.indicator_of_notMem hτ, mul_zero,
          Complex.ofReal_zero, mul_zero]
    · rw [Complex.ofReal_zero, mul_zero]
  rw [hfun, integral_finsetSum _ (fun p _ => integrable_finsetSum _ fun q _ => hint p q)]
  refine Finset.sum_congr rfl fun p _ => ?_
  rw [integral_finsetSum _ (fun q _ => hint p q)]
  refine Finset.sum_congr rfl fun q _ => ?_
  simp only [hT]
  split_ifs with hpq
  · have hqp : ((q : ℂ) - p) ≠ 0 := by
      rw [sub_ne_zero]; exact_mod_cast hpq.ne'
    have hexp : ∫ τ in p..q, cexp (ω * τ * I) = (cexp (ω * I * q) - cexp (ω * I * p)) / (ω * I) := by
      rw [← integral_exp_mul_complex hc]
      congr 1
      funext τ
      congr 1
      ring
    rw [integral_indicator measurableSet_Ioo, ← integral_Ioc_eq_integral_Ioo,
      ← intervalIntegral.integral_of_le hpq.le, intervalIntegral.integral_mul_const, hexp]
    push_cast
    field_simp
  · simp

section Slice

variable {P g : Matrix (Fin M) (Fin M) ℂ} (hP : P.IsHermitian) (hg : g.IsHermitian)
  (hspec : ∀ p ∈ spec hP, 0 ≤ p ∧ p ≤ 1)
include hP hg hspec

lemma integrable_rhoG_prod : Integrable (fun q : ℝ × ℝ => rhoG P g q.1 q.2) :=
  (densityRegG_rhoG hP hg hspec).integrable

lemma integrable_sliceMass : Integrable (sliceMass P g) :=
  (integrable_rhoG_prod hP hg hspec).integral_prod_right

lemma integral_sliceMass :
    ∫ τ, sliceMass P g τ = ∫ q : ℝ × ℝ, rhoG P g q.1 q.2 :=
  (integral_prod_symm _ (integrable_rhoG_prod hP hg hspec)).symm

/-- The Fourier transform of `μ` is the Laplace transform of `ρ` at `(0, -iω)`. -/
lemma fourier_sliceMass (ω : ℝ) :
    ∫ τ : ℝ, cexp ((ω : ℂ) * (τ : ℂ) * I) * (sliceMass P g τ : ℂ) = laplaceG P g 0 (-(ω * I)) := by
  have e : (fun q : ℝ × ℝ => cexp (0 * q.1 - -(ω * I) * q.2) * (rhoG P g q.1 q.2 : ℂ))
      = fun q => cexp (ω * q.2 * I) * (rhoG P g q.1 q.2 : ℂ) := by
    funext q
    congr 2
    ring
  have hI : Integrable (fun q : ℝ × ℝ => cexp (ω * q.2 * I) * (rhoG P g q.1 q.2 : ℂ)) := by
    have h := integrable_laplaceG_integrand hP hg hspec (hyp_RI_gen M) 0 (-(ω * I))
    rw [e] at h
    exact h
  unfold laplaceG
  rw [e, Measure.volume_eq_prod, integral_prod_symm _ hI]
  refine integral_congr_ae (Eventually.of_forall fun τ => ?_)
  simp only
  unfold sliceMass
  rw [← integral_complex_ofReal, ← integral_const_mul]

/-- **The `τ`-slices of `ρ`** (THEOREMS.md, Theorem A): for almost every `τ`,
`∫ ρ(s, τ) ds = Σ_{p<q} ‖Q_p g Q_q‖_F² / (q - p) · 1[p < τ < q]`. -/
theorem sliceMass_ae_eq : sliceMass P g =ᵐ[volume] stepMass hP g := by
  refine ae_eq_of_fourier_eq (integrable_sliceMass hP hg hspec) (integrable_stepMass hP g)
    (sliceMass_nonneg P g) (stepMass_nonneg hP g) fun ω => ?_
  rcases eq_or_ne ω 0 with rfl | hω
  · simp only [Complex.ofReal_zero, zero_mul, Complex.exp_zero, one_mul]
    rw [integral_complex_ofReal, integral_complex_ofReal]
    congr 1
    rw [integral_sliceMass hP hg hspec, total_mass hP hg hspec, integral_stepMass]
  · have ht : -((ω : ℂ) * I) ≠ 0 := neg_ne_zero.mpr (mul_ne_zero (by exact_mod_cast hω) I_ne_zero)
    rw [fourier_sliceMass hP hg hspec, tau_marginal hP hg hspec _ ht, fourier_stepMass hP g hω]
    refine Finset.sum_congr rfl fun p _ => Finset.sum_congr rfl fun q _ => ?_
    split_ifs with hpq
    · rw [trace_sandwich_eq_frob hP hg, ← frobNorm_sandwich_sq hP hg]
      congr 1
      rw [show -(-((ω : ℂ) * I)) * (p : ℂ) = (ω : ℂ) * I * p by ring,
        show -(-((ω : ℂ) * I)) * (q : ℂ) = (ω : ℂ) * I * q by ring,
        neg_mul, div_neg, ← neg_div, neg_sub]
    · rfl

end Slice

end HarmonicMajorization
