"""
ref_povm_batch.py -- run ref_povm.run over a list of (N, d, D).
Usage: python ref_povm_batch.py restarts
"""
import sys

from ref_povm import run

CASES = [(3, 3, 3), (4, 3, 3), (5, 3, 3), (3, 4, 4), (4, 4, 4), (3, 5, 5), (3, 3, 4), (3, 4, 5), (2, 3, 3), (6, 3, 3)]


def main():
    restarts = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    for (N, d, D) in CASES:
        best, target, defect, dt = run(N, d, D, restarts, seed=31 * N + 7 * d + D)
        flag = "  <-- BELOW I_BKP" if best < target - 1e-9 else ""
        print(f"POVM N={N} d={d} D={D}: best I = {best:.12f}  I_BKP = {target:.12f}  diff = {best - target:+.3e}  "
              f"({restarts} restarts, {dt:.0f}s; max projectivity defect {max(defect):.1e}){flag}", flush=True)


if __name__ == "__main__":
    main()
