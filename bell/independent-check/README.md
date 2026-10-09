# Internal independent check

This folder contains an internal independent check (not external refereeing) of
[`../../math/05-chained-bell.md`](../../math/05-chained-bell.md), carried out on 2026-10-08. The proofs of Theorems 1–5′
were checked step by step (Theorem 1 and Lemma 2 were re-derived), and the computations were redone with separate
code. Except for the `theirs` mode of
`ref_cert_verify.py` (see below), the scripts in this folder share no code with [`../scripts/`](../scripts/).

## What was verified

- **Mathematics.** Theorem 1, Lemma 2, Theorem 3, Proposition 4, Theorem 5 and Theorem 5′ were checked step by step.
  No mathematical error was found. The check found overstatements in the literature framing of an earlier draft (on the
  relation to the conjecture of Tavakoli, Zohren and Pawłowski, and to the quantum value over all states); the write-up
  uses the corrected wording. The documentation gaps it found (no log for N = 5, d = 11..14; twirl and BKP-value logs
  narrower than described) are closed by the logs in [`../logs/`](../logs/).
- **The expression and its value on $\Phi_d$** (`ref_reduction.py`, part (a)). BKP's $I_N$, evaluated from spectral
  projectors with BKP's bases, equals the closed form $I_N^{\mathrm{BKP}}(d)$ and the value of Tavakoli, Zohren and
  Pawłowski to $1.1\cdot10^{-14}$, for N = 2..8 and d = 2..8.
- **Theorem 1** (`ref_reduction.py`, part (b)). 20 random projective strategies with arbitrary ranks (zero ranks
  allowed), including D not divisible by d. The twirl is built with an explicit orthonormal basis of
  $\ker(W^{2N}-1)$, without an eigensolver. $I(R)=I(\text{twirl})=I(R^V)=N(d-1)-\frac2dF_N(V)$ to $4.1\cdot10^{-14}$,
  the two forms of $F_N$ (with $h_m$ and with $K_N$) agree, and $V_k^{2N}=1$ to $1.8\cdot10^{-14}$.
- **Lemma 2** (`ref_cellineq.py`, parts (1) and (2)). The partial-fraction formula at random real u to
  $1.4\cdot10^{-13}$ (relative; N, d = 2..8), and the full decomposition of $F_N$ on random configurations to
  $3.8\cdot10^{-14}$ (N = 2..7, d = 2..5).
- **Proposition 4** (`ref_cellineq.py`, part (3)). Adversarial maximisation (L-BFGS, 8 restarts per case) of the cell
  functional over chained configurations on $\mathbb C^2$ and $\mathbb C^3$, including s ≥ 2 and non-uniform cells. The
  maximum equals the bound $L^2\Phi_c$ to within $6\cdot10^{-13}$ and is never exceeded beyond rounding.
- **Theorem 3** (`ref_commutative.py`). Exhaustive enumeration for (N, d) = (3,3), (3,4), (4,3), (4,4), (5,3), (3,5),
  (6,3): every one-step configuration gives exactly $F_{\mathrm{nat}}$, and every other configuration lies below it by
  at least 1.42.
- **Theorem 5** (`ref_cert_verify.py`, `ref_cert_batch.py`). 26 verifications, all RIGOROUSLY VERIFIED, in ball
  arithmetic (python-flint/arb, 128 bits). The Clausen values come from the complex dilogarithm (a different algorithm
  from the Bernoulli series of `../scripts/ch_cert.py`), the cell kernel is computed directly from its definition, and
  the exact solution of the equalities is enclosed with `arb_mat.solve`. Mode `theirs` rebuilds the authors' certificate
  (same cell pool and seed; for this it imports `../scripts/ch_cone.py` and `../scripts/ch_clock.py`) and verifies it
  with the independent code. Mode `own` uses its own cell pool and linear programme. The runs include all three d = 24
  certificates (N = 3, 4, 5) and N = 5, d = 11..14. The authors' floating-point kernel agrees with the independent one on
  every used row to $6.2\cdot10^{-14}$ (relative).
- **Strategy-level searches** (`ref_minimize*.py`, `ref_povm*.py`), independent of the clock model: L-BFGS over
  projective strategies (2N unitaries with fixed outcome labels) and over POVM strategies on $\Phi_D$. No search went
  below $I_N^{\mathrm{BKP}}(d)$. When d divides D, the BKP value was reproduced (to $1.3\cdot10^{-14}$), and the POVM
  optima found were projective (defect at most $2\cdot10^{-15}$). When d does not divide D, the best values stayed above
  it by 0.31 to 0.91. These searches prove nothing; they only show that no simple counterexample exists.
- **The value over all states** (`ref_tzp_gap.py`). BKP's measurements on the lowest eigenvector of the matrix M of
  Tavakoli, Zohren and Pawłowski give a value strictly below the value on $\Phi_d$ (N = d = 3: 0.7129558005 against
  0.7567751189). So Theorem 5 does not determine the quantum value over all states.

## Files and commands

Requirements: as in [`../README.md`](../README.md), plus python-flint (`ref_cert_*.py`) and PyTorch (`ref_cellineq.py`,
`ref_minimize*.py`, `ref_povm*.py`; CPU is enough). Run the commands from this folder, with `OMP_NUM_THREADS=2`, and
redirect the output to the log. Run times are wall-clock times on a laptop with two threads, rounded.

| Log | Command | Run time | Notes |
|---|---|---|---|
| `ref_reduction.log` | `python ref_reduction.py` | 3 s | |
| `ref_cellineq.log` | `python ref_cellineq.py` | 45 s | prints "Lemma II.2", the draft numbering of Lemma 2 |
| `ref_commutative.log` | `python ref_commutative.py` | 1 s | |
| `ref_tzp_gap.log` | `python ref_tzp_gap.py` | 1 s | |
| `ref_cert_test_N3d3.log` | `python ref_cert_verify.py own 3 3`, then `python ref_cert_verify.py theirs 3 3` | 2 s | |
| `ref_cert_batch.log` | `python ref_cert_batch.py` | 4.5 min | 22 cases, listed in the script |
| `ref_cert_theirs_d24.log` | `python ref_cert_verify.py theirs 4 24 120 3`, then `python ref_cert_verify.py theirs 5 24 120 3` | 7 min | the authors' d = 24 certificates for N = 4, 5 (N = 3 is in `ref_cert_batch.log`) |
| `ref_minimize_batch_part2.log` | `python ref_minimize_batch.py 40` | 18 min | the run was stopped after the first 7 of its 22 cases |
| `ref_minimize_part3a.log` | `python ref_minimize_cases.py 20 8,3,3 3,6,6 3,4,8` | 8 min | |
| `ref_minimize_part3b.log` | `python ref_minimize_cases.py 20 3,7,7 4,5,5 5,4,4 3,4,4` | 9 min | |
| `ref_minimize_batch_part1_randomlabels.log` | none (see notes) | | first search, made with an earlier version of `ref_minimize.py` that drew all outcome labels at random; that version is not included; superseded by part 2, except for its cases with d dividing D, (3,3,3) and (3,3,6) |
| `ref_povm_batch.log` | `python ref_povm_batch.py 30` | 10 min | |

`ref_minimize.py` and `ref_povm.py` hold the search routines; they can also be run alone
(`python ref_minimize.py N d D restarts [seed]`, `python ref_povm.py N d D restarts [seed]`).

**Re-runs.** When this folder was assembled, `ref_reduction.py`, `ref_cellineq.py`, `ref_commutative.py` and
`ref_tzp_gap.py` were re-run here and reproduced their logs exactly. The three `ref_cert_*` logs were regenerated here
with the scripts as they are now, after one printed label in `ref_cert_verify.py` was reworded; every number agrees
with the original runs, and only the run times differ. The search logs (`ref_minimize*`, `ref_povm*`) were not re-run.
Apart from that label, the only changes made when this folder was assembled are in docstrings and comments, which now
refer to the write-up and to `../scripts/`.

**Numbering.** The checks were written against a draft that numbered the results II.1–II.6. Theorem II.k (Lemma II.k,
Proposition II.k) there is Theorem k (Lemma k, Proposition k) of the write-up, and Theorem II.5′ is Theorem 5′.
