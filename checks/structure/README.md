# Checks of the structure theory (`math/03-structure.md`)

This is an internal, adversarial pre-publication check, made on 2026-10-08. It is **not peer review**. The checker
re-derived every proof in `math/03-structure.md` and tested every statement with code written from scratch for this
check: these scripts share no code with our development scripts.

## Verdicts

| Statement | Verdict | Tests |
|---|---|---|
| Theorem A′ (holomorphic form) and its marginals | no errors found | identity to 1e-13; mass to 3e-15; x-marginal to 1e-10; oblique marginals to 1e-7 (`s1`) |
| Harmonic corollary and Theorem 1 (both directions, including the ε² scaling) | no errors found | `s1`, `s2` (50 digits) |
| Theorem D2 (two-level superharmonic majorization) | no errors found | 7,000 random pairs plus an adversary, no violation; the closed form of $h_\lambda$ against quadrature to 4e-31 (`s3`) |
| Exact two-level identity (R♭, T) | correct; the hypotheses were sharpened in the note | 8 cases to 5e-12–5e-10 (`s3`, Part 4) |
| Theorem D_sep and the trace identity | no errors found; D_sep also follows from Z. Long, arXiv:2610.10158, Lemma 2.3 (credited in the note) | trace identity to 4e-48; 60,000 random pairs and 15,000 oblique pairs at the borderline, no violation (`s4`); the derivation from Long's identity checked numerically (`s9`) |
| Theorem H (half-plane) | no errors found; the sketch was completed in the note | 400 pairs × 98 λ, $W^H\ge0$; $W^H=V^H$ to 1e-8–1e-16; second moment equals $\|K\|_F^2$ to 1e-16 (`s5`, `s5b`) |
| Section 4 ("What is false") | all counterexamples reproduced, exactly or at 60 digits | `s6` |
| Section 5, dephasing monotonicity | two levels: no violation (60-digit confirmation of every candidate); three levels: violated | `s7b`, `s7c` |
| Conjecture D | no violation in 20,000 random pairs with non-separable test functions | `s8` |

## Files

| Script | Purpose | Log (in `logs/`) |
|---|---|---|
| `hmaudit.py` | shared helpers: the density ρ_B, adaptive quadrature, the closed form of h_λ, high-precision spectra | — |
| `s1_thmAprime.py [seed]` | Theorem A′, marginals, harmonic corollary, Theorem 1 (⇐) | `s1_thmAprime_seed1.log`, `s1_thmAprime_seed2.log` |
| `s2_thm1_necessity.py` | Theorem 1 (⇒) and the ε² scaling | `s2_thm1_necessity.log` |
| `s3_D2.py [seed] [n]` | h_λ closed form, W = V, Theorem D2 (random and adversarial), the two-level identity | `s3_D2_seed7.log`, `s3_D2_seed8.log` |
| `s4_Dsep.py [seed] [n]` | trace identity, D_sep, oblique ridges, adversary | `s4_Dsep_seed11.log`, `s4_Dsep_seed12.log` |
| `s5_halfplane.py [seed] [n]` | Theorem H: $W^H\ge0$ and $W^H=V^H$ | `s5_halfplane_seed21.log` |
| `s5b_halfplane_moments.py` | Theorem H: moments and tails, convex test functions | `s5b_halfplane_moments.log` |
| `s6_counterexamples.py` | the counterexamples of Section 4 | `s6_counterexamples.log` |
| `s7b_dephasing_two_level.py [seed] [rounds]` | dephasing line, two levels, adversary with 60-digit confirmation | `s7b_dephasing_two_level.log` |
| `s7c_three_level_confirm.py` | the three-level dephasing counterexample (data in `logs/three_level_candidate.json`) | `s7c_three_level_confirm.log` |
| `s8_conjD_probe.py` | probe of Conjecture D | `s8_conjD_probe.log` |
| `s9_long_identity.py` | D_sep from Long's Lemma 2.3, numerically | `s9_long_identity.log` |

Run each script from this folder with Python 3 and numpy, scipy, mpmath and sympy. Set `OMP_NUM_THREADS=2` when running
several at once. The logs were produced with the default arguments, or with the seeds in their names. `s3_D2_seed7.log` and
`s4_Dsep_seed11.log` were produced before an overflow-safe log cosh was added to the scripts; the seed 8 and seed 12
runs use the final scripts. Some logs begin
with floating-point overflow warnings from test functions evaluated far from the spectrum; these do not affect the
results, which are confirmed in higher precision where it matters. The three-level candidate was found by a random and
adversarial search that is not included; the candidate itself is in `logs/three_level_candidate.json`.
