"""
ref_minimize_cases.py -- run ref_minimize.run for explicit cases.
Usage: python ref_minimize_cases.py restarts N,d,D [N,d,D ...]
"""
import sys

from ref_minimize import run


def main():
    restarts = int(sys.argv[1])
    for spec in sys.argv[2:]:
        N, d, D = (int(v) for v in spec.split(","))
        best, target, info, dt = run(N, d, D, restarts, seed=17 * N + 5 * d + D)
        flag = "  <-- BELOW I_BKP" if best < target - 1e-9 else ""
        print(f"N={N} d={d} D={D}: best I = {best:.12f}  I_BKP = {target:.12f}  diff = {best - target:+.3e}  "
              f"({restarts} restarts, {dt:.0f}s; ranks {info}){flag}", flush=True)


if __name__ == "__main__":
    main()
