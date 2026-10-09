# Chained Bell inequalities: scripts and logs

Computations for [`../math/05-chained-bell.md`](../math/05-chained-bell.md): the reduction to the clock model, the
commutative sector, the interval-certified cell certificates and the numerical checks.

- [`scripts/`](scripts/): the programs.
- [`logs/`](logs/): their outputs, as cited in the write-up.
- [`independent-check/`](independent-check/README.md): an internal independent check (not external refereeing), with
  its own scripts and logs.

## Requirements

- Python 3 (tested with 3.12) with numpy, scipy and mpmath. Tested versions: numpy 2.5.3, scipy 1.18.1 (its HiGHS solver
  is used for the linear programmes), mpmath 1.3.0 (its interval type `mpmath.iv`, with outward rounding, is used for
  the certificates).
- The independent check also needs python-flint 0.9.0 (arb ball arithmetic) and PyTorch 2.11.0 (CPU is enough).
- The logs were produced with `OMP_NUM_THREADS=2` and `OPENBLAS_NUM_THREADS=2`.

Run the commands from `bell/scripts/` (the scripts import `ch_clock.py` and `ch_cone.py` from the same folder) and
redirect the output to the log, for example `python ch_bkp.py 2 8 2 10 > ../logs/ch_bkp_N2-8_d2-10.log`. Run times are
wall-clock times on a laptop with two threads, rounded.

With N = 2, `ch_cone.py` reads the Regime I certificate cells of [27] for d = 2..6 from
[`scripts/oqp27_cells/`](scripts/oqp27_cells/). These five files are copied unchanged from
github.com/anshM123/IQOQI-OQP-27, `cone-certificates/certs/` (commit e7cdb7e). [27] is the OQP 27 mathematics paper,
`papers/math/main.tex` in that repository.

## Scripts

| Script | Purpose | Write-up |
|---|---|---|
| `ch_clock.py` | clock model: the reduction against direct evaluation, the commutative sector by enumeration, a non-commuting search | Theorems 1 and 3 |
| `ch_twirl_check.py` | the twirl reduction for random projective strategies on a maximally entangled state | Theorem 1 |
| `ch_bkp.py` | BKP's measurements on the maximally entangled state, against the natural strategy | Sections 2 and 4 |
| `ch_cone.py` | the floating-point cell LP (also imported by `ch_cert.py`) | Section 5, Conjecture 6 |
| `ch_cert.py`, `ch_cert_batch.py` | the interval-certified certificates; the batch driver retries with pools of 121 and 301 cells | Theorem 5 |
| `ch_validate.py` | the cell inequalities on random and near-optimal chained configurations | Proposition 4 |
| `ch_local.py` | Hessian of the deficit at pairs of windows; seeded non-commuting searches | Section 5 |

## Logs

| Log | Command (from `bell/scripts/`) | Run time | Notes |
|---|---|---|---|
| `ch_twirl_check_N3-5_d3-5_D3-4-5-7.log` | `python ch_twirl_check.py N d D 3` for N = 3, 4, 5, d = 3, 4, 5, D = 3, 4, 5, 7 (in this order, N outermost) | 20 s | 108 random strategies |
| `ch_twirl_check.log` | `python ch_twirl_check.py 3 3 2 3`, then `python ch_twirl_check.py 4 3 2 2` | 2 s | an older, narrower run |
| `ch_bkp_N2-8_d2-10.log` | `python ch_bkp.py 2 8 2 10` | 1 s | |
| `ch_cert_batch_N2-6_d2-10.log` | `python ch_cert_batch.py 2 6 2 10` | 3.5 min | its d = 2 lines come from an earlier version, which stopped with an exception for empty certificates; the current script certifies d = 2 (see the next log), so a re-run differs in those lines only |
| `ch_cert_batch_d2.log` | `python ch_cert_batch.py 2 8 2 2` | 7 s | no cell is needed |
| `ch_cert_batch_N7-8_d3-10.log` | `python ch_cert_batch.py 7 8 3 10` | 6 min | |
| `ch_cert_batch_N3-4_d11-14.log` | `python ch_cert_batch.py 3 4 11 14` | 1.5 min | |
| `ch_cert_batch_N3-5_d15-24.log` | `python ch_cert_batch.py 3 5 15 24` | 25 min | for d = 24 the first pool (41 cells) fails and the second (121 cells) gives the certificate |
| `ch_cert_batch_N5_d11-14.log` | `python ch_cert_batch.py 5 5 11 14` | 1.5 min | |
| `ch_validate.log` | `python ch_validate.py 30` | 2 s | |
| `ch_local_search.log` | `python ch_local.py N d 12` for (N, d) = (3, 3), (3, 4), (4, 3), (3, 5), in this order | 16 s | |
| `ch_cone_N2_d2-6_oqp27cells.log` | `python ch_cone.py 2 d` for d = 2, ..., 6 | 3 s | recorded when this folder was assembled |
| `ch_cone_N3_d16_N4_d12.log` | `python ch_cone.py 3 16 40`, then `python ch_cone.py 4 12 40` | 20 s | recorded when this folder was assembled |

`ch_clock.py` has no archived log. For example, `python ch_clock.py 3 3 2 4` (N = 3, d = 3, M = 2, 4 restarts) prints the
reduction check, the commutative sector and a short search in about a second.

The same commands as one bash script, run from `bell/scripts/`:

```bash
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
L=../logs
for N in 3 4 5; do for d in 3 4 5; do for D in 3 4 5 7; do python ch_twirl_check.py $N $d $D 3; done; done; done > $L/ch_twirl_check_N3-5_d3-5_D3-4-5-7.log
{ python ch_twirl_check.py 3 3 2 3; python ch_twirl_check.py 4 3 2 2; } > $L/ch_twirl_check.log
python ch_bkp.py 2 8 2 10 > $L/ch_bkp_N2-8_d2-10.log
python ch_cert_batch.py 2 6 2 10 > $L/ch_cert_batch_N2-6_d2-10.log
python ch_cert_batch.py 2 8 2 2 > $L/ch_cert_batch_d2.log
python ch_cert_batch.py 7 8 3 10 > $L/ch_cert_batch_N7-8_d3-10.log
python ch_cert_batch.py 3 4 11 14 > $L/ch_cert_batch_N3-4_d11-14.log
python ch_cert_batch.py 3 5 15 24 > $L/ch_cert_batch_N3-5_d15-24.log
python ch_cert_batch.py 5 5 11 14 > $L/ch_cert_batch_N5_d11-14.log
python ch_validate.py 30 > $L/ch_validate.log
for nd in "3 3" "3 4" "4 3" "3 5"; do python ch_local.py $nd 12; done > $L/ch_local_search.log
for d in 2 3 4 5 6; do python ch_cone.py 2 $d; done > $L/ch_cone_N2_d2-6_oqp27cells.log
{ python ch_cone.py 3 16 40; python ch_cone.py 4 12 40; } > $L/ch_cone_N3_d16_N4_d12.log
```

## Re-runs

When this folder was assembled, every log above was generated again from the scripts in `scripts/`. The re-runs agree
with the archived logs exactly, apart from the printed run times; for `ch_cert_batch_N2-6_d2-10.log` this holds for all
lines except the d = 2 entries described in the table.

## Changes made when this folder was assembled

- Docstrings now refer to the write-up and to the public papers [27] and [CL] (github.com/anshM123/CGLMP,
  `paper-classical-all-d/main.tex`).
- `ch_cone.py` reads the certificate cells of [27] from `scripts/oqp27_cells/` instead of a local copy of the OQP 27
  repository.

No other line of code was changed.
