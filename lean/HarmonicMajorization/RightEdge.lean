import HarmonicMajorization.MainGen

/-!
# Theorem B, right edge (general `0 ≤ B ≤ 1`)

Authors: Ansh Mishra, Aryan Senthilkumar.  License: MIT.

THEOREMS.md, Theorem B, right edge: the right sweep of a multiset `Z` is the left sweep of
`{1 - conj ν : ν ∈ Z}`, and the identity of the left edge holds for the right sweeps with `ρ` replaced
by the density of the pair `(B, g)`, i.e. of `det(g - s - ξ(B - τ))`.  This is the left edge for `1 - B`;
the two facts needed to transfer it are
* `pinchH_one_sub`: the pinchings of `g` by the spectral projections of `B` and of `1 - B` agree
  (`sproj_one_sub`: `Q'_{1-b} = Q_b`; spectral projections are Lagrange polynomials in the matrix,
  `sproj_eq_aeval_lagP`);
* `roots_charpoly_one_sub_add`: `spec((1 - B) + ig) = {1 - conj ν : ν ∈ spec(B + ig)}` (with multiplicity).

Main results (no hypotheses): `theoremB_right`, `theoremB_right_ineq`, `theoremB_right_eq_iff`.
-/

set_option autoImplicit false

namespace HarmonicMajorization

open Matrix Polynomial Complex
open scoped Real ComplexOrder
open OQP27.StripL3a (stripKernel stripBalayage hLam)

variable {M : ℕ}

/-! ### Spectral projections are polynomials in the matrix -/

/-- The Lagrange polynomial `Π_{q ∈ S, q ≠ p} (X - q)/(p - q)`. -/
noncomputable def lagP (S : Finset ℝ) (p : ℝ) : ℂ[X] :=
  ∏ q ∈ S.erase p, (C ((((p - q : ℝ)) : ℂ)⁻¹) * (X - C (q : ℂ)))

lemma eval_lagP (S : Finset ℝ) (p : ℝ) (x : ℂ) :
    (lagP S p).eval x = ∏ q ∈ S.erase p, ((((p - q : ℝ)) : ℂ)⁻¹ * (x - q)) := by
  unfold lagP
  rw [Polynomial.eval_prod]
  simp

lemma eval_lagP_self (S : Finset ℝ) (p : ℝ) : (lagP S p).eval (p : ℂ) = 1 := by
  rw [eval_lagP]
  refine Finset.prod_eq_one fun q hq => ?_
  have hne : q ≠ p := Finset.ne_of_mem_erase hq
  have h : (((p - q : ℝ)) : ℂ) ≠ 0 := by exact_mod_cast sub_ne_zero.mpr hne.symm
  rw [show ((p : ℂ) - q) = (((p - q : ℝ)) : ℂ) by push_cast; ring, inv_mul_cancel₀ h]

lemma eval_lagP_of_mem {S : Finset ℝ} {p r : ℝ} (hr : r ∈ S) (hrp : r ≠ p) :
    (lagP S p).eval (r : ℂ) = 0 := by
  rw [eval_lagP]
  exact Finset.prod_eq_zero (Finset.mem_erase.mpr ⟨hrp, hr⟩) (by rw [sub_self, mul_zero])

section Poly

variable {P : Matrix (Fin M) (Fin M) ℂ} (hP : P.IsHermitian)

/-- Polynomials in `P` through the functional calculus. -/
lemma aeval_eq_fcalc (f : ℂ[X]) :
    aeval P f = fcalc hP (fun i => f.eval (hP.eigenvalues i : ℂ)) := by
  induction f using Polynomial.induction_on with
  | C a =>
    rw [aeval_C, Algebra.algebraMap_eq_smul_one]
    simp only [eval_C]
    rw [fcalc_const]
  | add f g hf hg =>
    rw [map_add, hf, hg, ← fcalc_add]
    congr 1
    funext i
    simp [eval_add]
  | monomial n a h =>
    rw [pow_succ, ← mul_assoc, map_mul, h, aeval_X]
    have e : (fun i => eval (hP.eigenvalues i : ℂ) (C a * X ^ n * X))
        = (fun i => eval (hP.eigenvalues i : ℂ) (C a * X ^ n)) * (fun i => (hP.eigenvalues i : ℂ)) := by
      funext i
      simp only [Pi.mul_apply, eval_mul, eval_X]
    rw [e, ← fcalc_mul, fcalc_eigenvalues]

lemma sproj_eq_aeval {p : ℝ} (f : ℂ[X]) (hf : ∀ i, f.eval (hP.eigenvalues i : ℂ) = ind hP p i) :
    sproj hP p = aeval P f := by
  rw [aeval_eq_fcalc hP]
  unfold sproj
  congr 1
  funext i
  rw [hf]

/-- The spectral projection `Q_p` is the Lagrange polynomial of `spec P` at `p`, evaluated at `P`. -/
lemma sproj_eq_aeval_lagP (p : ℝ) : sproj hP p = aeval P (lagP (spec hP) p) := by
  refine sproj_eq_aeval hP _ fun i => ?_
  unfold ind
  split_ifs with h
  · rw [h, eval_lagP_self]
  · exact eval_lagP_of_mem (eigenvalues_mem_spec hP i) h

lemma mem_spec_iff_det (x : ℝ) :
    x ∈ spec hP ↔ (P - (x : ℂ) • (1 : Matrix (Fin M) (Fin M) ℂ)).det = 0 := by
  constructor
  · exact det_sub_smul_eq_zero hP
  · intro h
    by_contra hx
    exact (isUnit_det_sub hP hx).ne_zero h

end Poly

/-! ### `B` and `1 - B` have the same spectral projections -/

section OneSub

variable {B : Matrix (Fin M) (Fin M) ℂ} (hB : B.IsHermitian) (h1B : (1 - B).IsHermitian)
include hB h1B

omit hB h1B in
lemma det_one_sub_sub (x : ℝ) :
    ((1 - B) - (x : ℂ) • (1 : Matrix (Fin M) (Fin M) ℂ)).det
      = (-1) ^ M * (B - (((1 - x : ℝ)) : ℂ) • (1 : Matrix (Fin M) (Fin M) ℂ)).det := by
  have e : (1 - B) - (x : ℂ) • (1 : Matrix (Fin M) (Fin M) ℂ)
      = -(B - (((1 - x : ℝ)) : ℂ) • (1 : Matrix (Fin M) (Fin M) ℂ)) := by
    rw [Complex.ofReal_sub, Complex.ofReal_one, sub_smul, one_smul]
    abel
  rw [e, Matrix.det_neg, Fintype.card_fin]

lemma mem_spec_one_sub_iff (x : ℝ) : x ∈ spec h1B ↔ 1 - x ∈ spec hB := by
  rw [mem_spec_iff_det h1B, mem_spec_iff_det hB, det_one_sub_sub]
  constructor
  · intro h
    rcases mul_eq_zero.mp h with h' | h'
    · exact absurd h' (pow_ne_zero _ (by norm_num))
    · exact h'
  · intro h
    rw [h, mul_zero]

lemma spec_one_sub_eq : spec h1B = (spec hB).image (fun b => 1 - b) := by
  ext x
  rw [Finset.mem_image]
  constructor
  · intro hx
    exact ⟨1 - x, (mem_spec_one_sub_iff hB h1B x).mp hx, by ring⟩
  · rintro ⟨b, hb, rfl⟩
    rw [mem_spec_one_sub_iff hB h1B, sub_sub_cancel]
    exact hb

/-- `Q'_{1-b} = Q_b`: the spectral projections of `1 - B` are those of `B`. -/
lemma sproj_one_sub (b : ℝ) : sproj h1B (1 - b) = sproj hB b := by
  rw [sproj_eq_aeval_lagP h1B]
  have hcomp : aeval (1 - B) (lagP (spec h1B) (1 - b))
      = aeval B ((lagP (spec h1B) (1 - b)).comp (1 - X)) := by
    rw [Polynomial.aeval_comp]
    congr 1
    simp
  rw [hcomp]
  refine (sproj_eq_aeval hB _ fun i => ?_).symm
  rw [Polynomial.eval_comp]
  simp only [eval_sub, eval_one, eval_X]
  have hcast : (1 : ℂ) - (hB.eigenvalues i : ℂ) = (((1 - hB.eigenvalues i : ℝ)) : ℂ) := by
    push_cast; ring
  rw [hcast]
  unfold ind
  split_ifs with h
  · rw [h]
    exact eval_lagP_self _ _
  · refine eval_lagP_of_mem ?_ ?_
    · rw [mem_spec_one_sub_iff hB h1B, sub_sub_cancel]
      exact eigenvalues_mem_spec hB i
    · intro h'
      exact h (by linarith)

/-- The pinchings of `g` by the spectral projections of `B` and of `1 - B` agree. -/
theorem pinchH_one_sub (g : Matrix (Fin M) (Fin M) ℂ) : pinchH h1B g = pinchH hB g := by
  unfold pinchH
  rw [spec_one_sub_eq hB h1B, Finset.sum_image (fun x _ y _ h => by linarith)]
  refine Finset.sum_congr rfl fun b _ => ?_
  rw [sproj_one_sub hB h1B]

end OneSub

lemma pinchH_congr {P P' : Matrix (Fin M) (Fin M) ℂ} (h : P = P') (hP : P.IsHermitian)
    (hP' : P'.IsHermitian) (g : Matrix (Fin M) (Fin M) ℂ) : pinchH hP g = pinchH hP' g := by
  subst h
  rfl

/-! ### The spectrum of `(1 - B) + ig` -/

lemma roots_charpoly_one_sub (Y : Matrix (Fin M) (Fin M) ℂ) :
    (Matrix.charpoly (1 - Y)).roots = (Matrix.charpoly Y).roots.map (fun x => 1 - x) := by
  have e : (1 : Matrix (Fin M) (Fin M) ℂ) - Y = (-1 : ℂ) • Y - Matrix.scalar (Fin M) (-1 : ℂ) := by
    rw [Matrix.scalar_apply, ← Matrix.smul_one_eq_diagonal, neg_one_smul, neg_one_smul]
    abel
  rw [e, Matrix.charpoly_sub_scalar]
  have h3 : (Polynomial.X + Polynomial.C (-1 : ℂ) : Polynomial ℂ)
      = Polynomial.C 1 * Polynomial.X + Polynomial.C (-1 : ℂ) := by simp
  rw [h3, Polynomial.roots_comp_C_mul_X_add_C _ 1 (-1 : ℂ) isUnit_one,
    OQP27.StripL3b.roots_charpoly_smul _ (by norm_num : (-1 : ℂ) ≠ 0), Multiset.map_map]
  congr 1
  funext x
  simp only [Function.comp_apply, Ring.inverse_one, one_mul]
  ring

lemma roots_charpoly_conjTranspose (N : Matrix (Fin M) (Fin M) ℂ) :
    (Matrix.charpoly Nᴴ).roots = (Matrix.charpoly N).roots.map (starRingEnd ℂ) := by
  have h1 : Matrix.charpoly Nᴴ = (Matrix.charpoly N).map (starRingEnd ℂ) := by
    have e : Nᴴ = (Nᵀ).map (starRingEnd ℂ) := rfl
    rw [e, Matrix.charpoly_map, Matrix.charpoly_transpose]
  rw [h1, Polynomial.roots_map_of_injective_of_card_eq_natDegree (starRingEnd ℂ).injective
    (IsAlgClosed.card_roots_eq_natDegree)]

/-- `spec((1 - B) + ig) = {1 - conj ν : ν ∈ spec(B + ig)}`, with multiplicities. -/
theorem roots_charpoly_one_sub_add {B g : Matrix (Fin M) (Fin M) ℂ} (hB : B.IsHermitian)
    (hg : g.IsHermitian) :
    ((1 - B) + I • g).charpoly.roots = (B + I • g).charpoly.roots.map (fun ν => 1 - star ν) := by
  have e : (1 - B) + I • g = 1 - (B + I • g)ᴴ := by
    rw [conjTranspose_add, conjTranspose_smul, hB.eq, hg.eq, Complex.star_def, Complex.conj_I,
      neg_smul]
    abel
  rw [e, roots_charpoly_one_sub, roots_charpoly_conjTranspose, Multiset.map_map]
  rfl

/-! ### Theorem B at the right edge -/

section Right

variable {B g : Matrix (Fin M) (Fin M) ℂ} (hB0 : B.PosSemidef) (hB1 : (1 - B).PosSemidef)
  (hg : g.IsHermitian)
include hB0 hB1 hg

omit hB1 hg in
lemma one_sub_one_sub_posSemidef : (1 - (1 - B)).PosSemidef := by
  rw [sub_sub_cancel]; exact hB0

omit hg in
lemma pinchH_one_sub_one_sub (h : (1 - (1 - B)).IsHermitian) :
    pinchH h g = pinchH hB1.isHermitian g := by
  rw [pinchH_congr (sub_sub_cancel 1 B) h hB0.isHermitian g,
    pinchH_one_sub hB0.isHermitian hB1.isHermitian g]

/-- **Theorem B, right edge** (THEOREMS.md, Section 1; no hypotheses): with the right sweep
`ν ↦ 1 - conj ν` and the density `ρ_B` of `det(g - s - ξ(B - τ))`, for every real `λ`,
`Σ_{ν ∈ spec(B+ig)} h_λ(1 - conj ν) - Σ_{ν ∈ spec(B+ig_d)} h_λ(1 - conj ν) = V_B(λ)`. -/
theorem theoremB_right (lam : ℝ) :
    ((B + I • g).charpoly.roots.map (fun ν => hLam lam (1 - star ν))).sum
      - ((B + I • pinchH hB1.isHermitian g).charpoly.roots.map (fun ν => hLam lam (1 - star ν))).sum
      = stripBalayage (rhoG B g) lam := by
  have h := theoremB hB1 (one_sub_one_sub_posSemidef hB0) hg lam
  rw [pinchH_one_sub_one_sub hB0 hB1] at h
  rw [sub_sub_cancel] at h
  rw [roots_charpoly_one_sub_add hB0.isHermitian hg,
    roots_charpoly_one_sub_add hB0.isHermitian (isHermitian_pinchH hB1.isHermitian hg),
    Multiset.map_map, Multiset.map_map] at h
  exact h

/-- **Theorem B, right edge, inequality**. -/
theorem theoremB_right_ineq (lam : ℝ) :
    ((B + I • pinchH hB1.isHermitian g).charpoly.roots.map (fun ν => hLam lam (1 - star ν))).sum
      ≤ ((B + I • g).charpoly.roots.map (fun ν => hLam lam (1 - star ν))).sum := by
  have h := theoremB_ineq hB1 (one_sub_one_sub_posSemidef hB0) hg lam
  rw [pinchH_one_sub_one_sub hB0 hB1] at h
  rw [roots_charpoly_one_sub_add hB0.isHermitian hg,
    roots_charpoly_one_sub_add hB0.isHermitian (isHermitian_pinchH hB1.isHermitian hg),
    Multiset.map_map, Multiset.map_map] at h
  exact h

/-- **Theorem B, right edge, equality case**: equality at one `λ` iff `[B, g] = 0`. -/
theorem theoremB_right_eq_iff (lam : ℝ) :
    ((B + I • g).charpoly.roots.map (fun ν => hLam lam (1 - star ν))).sum
        = ((B + I • pinchH hB1.isHermitian g).charpoly.roots.map
            (fun ν => hLam lam (1 - star ν))).sum
      ↔ Commute B g := by
  have h := theoremB_eq_iff hB1 (one_sub_one_sub_posSemidef hB0) hg lam
  rw [pinchH_one_sub_one_sub hB0 hB1] at h
  rw [roots_charpoly_one_sub_add hB0.isHermitian hg,
    roots_charpoly_one_sub_add hB0.isHermitian (isHermitian_pinchH hB1.isHermitian hg),
    Multiset.map_map, Multiset.map_map] at h
  refine h.trans ?_
  constructor
  · intro hc
    have : B = 1 - (1 - B) := by abel
    rw [this]
    exact (Commute.one_left g).sub_left hc
  · intro hc
    exact (Commute.one_left g).sub_left hc

end Right

end HarmonicMajorization
