import HarmonicMajorization.RightEdge
import HarmonicMajorization.SliceGen
import HarmonicMajorization.MassGen

/-! Axiom check for the multi-line Radon identity and the unconditional Theorems A, B, C
(general `0 ≤ B ≤ 1`).  Expected output for each: `[propext, Classical.choice, Quot.sound]`. -/

#print axioms HarmonicMajorization.pencilRoot_im_bounds_gen
#print axioms HarmonicMajorization.pencilRoot_norm_bound_gap
#print axioms HarmonicMajorization.norm_Sup_sub_le_gap
#print axioms HarmonicMajorization.upper_global_gen
#print axioms HarmonicMajorization.dichotomy
#print axioms HarmonicMajorization.charpoly_WmatP_self_pinch
#print axioms HarmonicMajorization.Gup_cancel
#print axioms HarmonicMajorization.exists_int_GupG_left
#print axioms HarmonicMajorization.PsiG_sub_E_affine
#print axioms HarmonicMajorization.tendsto_imAbsSum_vertical_gen
#print axioms HarmonicMajorization.tendsto_im_PsiG
#print axioms HarmonicMajorization.radon_identity_gaps
#print axioms HarmonicMajorization.integrable_imAbsSum_gen
#print axioms HarmonicMajorization.radon_identity_gen
#print axioms HarmonicMajorization.hyp_RI_gen
#print axioms HarmonicMajorization.theoremA
#print axioms HarmonicMajorization.tau_marginal
#print axioms HarmonicMajorization.sine_moment
#print axioms HarmonicMajorization.theoremB
#print axioms HarmonicMajorization.theoremB_ineq
#print axioms HarmonicMajorization.theoremB_eq_iff
#print axioms HarmonicMajorization.theoremC_balayage
#print axioms HarmonicMajorization.theoremC
#print axioms HarmonicMajorization.sproj_eq_aeval_lagP
#print axioms HarmonicMajorization.pinchH_one_sub
#print axioms HarmonicMajorization.roots_charpoly_one_sub_add
#print axioms HarmonicMajorization.theoremB_right
#print axioms HarmonicMajorization.theoremB_right_ineq
#print axioms HarmonicMajorization.theoremB_right_eq_iff
#print axioms HarmonicMajorization.theoremA_package
#print axioms HarmonicMajorization.total_mass
#print axioms HarmonicMajorization.sliceMass_ae_eq
#print axioms HarmonicMajorization.frobNorm_sub_pinchH_sq
#print axioms HarmonicMajorization.total_mass_frob
