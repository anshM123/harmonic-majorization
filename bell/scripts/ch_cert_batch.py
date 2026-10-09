"""
ch_cert_batch.py -- run ch_cert.certify over a range of (N, d); retries with larger cell pools if needed.
Usage: python ch_cert_batch.py Nmin Nmax dmin dmax
"""
import sys

from ch_cert import certify


def main():
    Nmin, Nmax, dmin, dmax = map(int, sys.argv[1:5])
    summary = []
    for N in range(Nmin, Nmax + 1):
        for d in range(dmin, dmax + 1):
            ok = False
            for npool, seed in [(40, 2), (120, 3), (300, 4)]:
                try:
                    ok = certify(N, d, npool=npool, q=8, seed=seed, prec=120)
                except Exception as e:      # report and retry
                    print(f"N={N} d={d}: exception {e!r}", flush=True)
                    ok = False
                if ok:
                    break
            summary.append((N, d, ok))
    print("SUMMARY")
    for N, d, ok in summary:
        print(f"  N={N} d={d}: {'CERTIFIED' if ok else 'not certified'}")


if __name__ == "__main__":
    main()
