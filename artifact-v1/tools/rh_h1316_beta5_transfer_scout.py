#!/usr/bin/env python3
"""Numerical sign scout for the exact m=1, beta=5 Xi correction.

This is deliberately a scout, not an interval certificate.  It integrates in
the normalized saddle coordinate of the beta=9 gamma-log law and compares the
third q-cumulant before and after multiplication by

    1 - (3/(2*pi))*exp(-r).

The sign of the difference is the sign of the third q-derivative of the
logarithmic beta=5 correction.  The corresponding t-derivative is eight times
larger and has the same sign because q=2t.
"""

from __future__ import annotations

import argparse
import mpmath as mp


ALPHA = mp.mpf(9) / 4


def saddle_r(q: mp.mpf, alpha: mp.mpf = ALPHA) -> mp.mpf:
    """Solve pi*r*exp(r) = q + alpha*r + 1."""

    guess = mp.lambertw(q / mp.pi)
    return mp.findroot(
        lambda r: mp.pi * r * mp.exp(r) - q - alpha * r - 1,
        (guess, guess * (1 + mp.mpf("1e-6"))),
    )


def normalized_moments(q: mp.mpf, dps: int = 80) -> dict[str, mp.mpf]:
    mp.mp.dps = dps
    r0 = saddle_r(q)
    z0 = mp.log(r0)
    a2 = mp.pi * mp.exp(r0) * r0 * (r0 + 1) - ALPHA * r0
    sqrt_a2 = mp.sqrt(a2)
    c5 = mp.mpf(3) / (2 * mp.pi)

    def pieces(w: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
        s = w / sqrt_a2
        r = r0 * mp.exp(s)
        # q*s - (V(z0+s)-V(z0)); written without subtracting q-sized
        # quantities.  The saddle cancellation is exact in this form.
        psi = (
            mp.pi * (mp.exp(r) - mp.exp(r0))
            - ALPHA * (r - r0)
            - (q + 1) * s
        )
        base = mp.exp(-psi)
        return base, base * (1 - c5 * mp.exp(-r))

    # The normalized law is already close to N(0,1).  Wide segmented
    # integration also covers the modest-q diagnostic points.
    cuts = [mp.ninf, -24, -16, -10, -6, -3, 0, 3, 6, 10, 16, 24, mp.inf]

    def integrate(which: int, power: int) -> mp.mpf:
        return mp.quad(lambda w: (w**power) * pieces(w)[which], cuts)

    out: dict[str, mp.mpf] = {"q": q, "r": r0, "A": a2}
    for name, which in (("dominant", 0), ("corrected", 1)):
        vals = [integrate(which, j) for j in range(4)]
        mean = vals[1] / vals[0]
        m2 = vals[2] / vals[0]
        m3 = vals[3] / vals[0]
        kappa3_w = m3 - 3 * mean * m2 + 2 * mean**3
        out[f"mass_{name}"] = vals[0]
        out[f"mean_w_{name}"] = mean
        out[f"kappa3_w_{name}"] = kappa3_w
        out[f"kappa3_z_{name}"] = kappa3_w / a2 ** mp.mpf("1.5")
    out["g3_q"] = out["kappa3_z_corrected"] - out["kappa3_z_dominant"]
    out["g3_t"] = 8 * out["g3_q"]
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--q",
        action="append",
        default=[],
        help="positive q value; may be repeated (default: 10,100,1e4,exp(100))",
    )
    parser.add_argument("--dps", type=int, default=80)
    args = parser.parse_args()
    qs = [mp.mpf(x) for x in args.q] if args.q else [10, 100, 10**4, mp.exp(100)]
    for q in qs:
        result = normalized_moments(mp.mpf(q), args.dps)
        print("q=", mp.nstr(result["q"], 18))
        print("r=", mp.nstr(result["r"], 18))
        print("g'''_q=", mp.nstr(result["g3_q"], 25))
        print("g'''_t=", mp.nstr(result["g3_t"], 25))
        print("sign=", "positive" if result["g3_q"] > 0 else "NONPOSITIVE")


if __name__ == "__main__":
    main()
