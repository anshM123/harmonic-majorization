"""
ref_minimize_batch.py -- run ref_minimize.run over a list of (N, d, D) in one process.
Usage: python ref_minimize_batch.py restarts
"""
import sys

from ref_minimize import run

CASES = [(3, 3, 4), (3, 3, 5), (4, 3, 3), (4, 3, 4), (4, 3, 5), (4, 3, 6), (5, 3, 3),
         (6, 3, 3), (8, 3, 3), (3, 4, 4), (3, 4, 5), (3, 4, 6), (3, 4, 8), (4, 4, 4), (5, 4, 4), (3, 5, 5),
         (4, 5, 5), (3, 5, 6), (3, 6, 6), (3, 7, 7), (2, 3, 3), (2, 4, 4)]


def main():
    restarts = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    for (N, d, D) in CASES:
        best, target, info, dt = run(N, d, D, restarts, seed=17 * N + 5 * d + D)
        flag = "  <-- BELOW I_BKP" if best < target - 1e-9 else ""
        print(f"N={N} d={d} D={D}: best I = {best:.12f}  I_BKP = {target:.12f}  diff = {best - target:+.3e}  "
              f"({restarts} restarts, {dt:.0f}s; ranks {info}){flag}", flush=True)


if __name__ == "__main__":
    main()
