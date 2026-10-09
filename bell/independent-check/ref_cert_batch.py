"""
ref_cert_batch.py -- run ref_cert_verify over a list of cases in one process.
Usage: python ref_cert_batch.py
"""
import time

from ref_cert_verify import mode_own, mode_theirs

CASES = [
    ("theirs", 4, 5, 40, 2), ("theirs", 8, 10, 40, 2), ("theirs", 5, 11, 40, 2), ("theirs", 5, 12, 40, 2),
    ("theirs", 5, 13, 40, 2), ("theirs", 5, 14, 40, 2),
    ("own", 2, 5, 40, 11), ("own", 3, 4, 40, 11), ("own", 4, 5, 40, 11), ("own", 5, 7, 40, 11),
    ("own", 6, 8, 60, 11), ("own", 7, 9, 60, 11), ("own", 8, 10, 60, 11), ("own", 5, 11, 60, 11),
    ("own", 5, 12, 60, 11), ("own", 5, 13, 60, 11), ("own", 5, 14, 60, 11), ("own", 3, 20, 120, 11),
    ("own", 4, 22, 120, 11), ("theirs", 3, 24, 120, 3), ("own", 3, 24, 150, 11), ("own", 5, 24, 150, 11),
]


def main():
    summary = []
    for (mode, N, d, npool, seed) in CASES:
        t0 = time.time()
        try:
            ok = (mode_own if mode == "own" else mode_theirs)(N, d, npool, seed, 128)
        except Exception as e:
            print(f"{mode} N={N} d={d}: EXCEPTION {e!r}", flush=True)
            ok = False
        summary.append((mode, N, d, npool, seed, ok, time.time() - t0))
    print("SUMMARY")
    for (mode, N, d, npool, seed, ok, dt) in summary:
        print(f"  {mode:6s} N={N} d={d} pool={npool} seed={seed}: {'RIGOROUSLY VERIFIED' if ok else 'not verified'} "
              f"({dt:.0f}s)")


if __name__ == "__main__":
    main()
