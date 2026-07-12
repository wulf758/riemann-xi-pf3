"""Numerical scout for the gamma-log to full-Xi transfer.

This is deliberately a high-precision finite-window diagnostic, not an
interval proof.  It evaluates the full Xi kernel as a multiplicative tilt of
the dominant m=1, beta=9 gamma-log measure and reconstructs epsilon and its
first three t-derivatives from cumulant identities.
"""

from __future__ import annotations

import argparse
import json

import mpmath as mp

import rh_h929_gamma_log_moment_box_smoke as h929


def kernel_ratio(x: mp.mpf, m_cutoff: int) -> mp.mpf:
    """Phi(U) / [2*pi^2*exp(9U-pi*exp(4U))], with x=pi*exp(4U)."""

    return mp.fsum(
        (
            mp.mpf(m) ** 4 - 3 * mp.mpf(m) ** 2 / (2 * x)
        )
        * mp.exp(-x * (mp.mpf(m) ** 2 - 1))
        for m in range(1, m_cutoff + 1)
    )


def raw_and_central(moments: list[mp.mpf]) -> dict[str, mp.mpf]:
    mass = moments[0]
    raw = [value / mass for value in moments]
    mean = raw[1]
    variance = raw[2] - mean**2
    third = raw[3] - 3 * mean * raw[2] + 2 * mean**3
    return {"mass": mass, "mean": mean, "variance": variance, "third": third}


def sample(r: mp.mpf, width: mp.mpf, m_cutoff: int) -> dict[str, mp.mpf]:
    data = h929.saddle_quantities(r)
    q = data["q"]
    t = q / 2
    A = data["A"]
    sqrt_A = mp.sqrt(A)
    z0 = data["z"]
    V0 = h929.V(z0)
    cuts = h929.integration_cuts(width)

    def density(w: mp.mpf) -> mp.mpf:
        psi = h929.V(z0 + w / sqrt_A) - V0 - q * w / sqrt_A
        return mp.exp(-psi)

    def ratio(w: mp.mpf) -> mp.mpf:
        local_r = mp.exp(z0 + w / sqrt_A)
        x = mp.pi * mp.exp(local_r)
        return kernel_ratio(x, m_cutoff)

    dominant_moments = [
        mp.quad(lambda w, order=order: w**order * density(w), cuts)
        for order in range(4)
    ]
    full_moments = [
        mp.quad(
            lambda w, order=order: w**order * density(w) * ratio(w),
            cuts,
        )
        for order in range(4)
    ]
    dominant = raw_and_central(dominant_moments)
    full = raw_and_central(full_moments)

    ratio_mass = full["mass"] / dominant["mass"]
    epsilon = ratio_mass - 1

    mean_gap_q = (full["mean"] - dominant["mean"]) / sqrt_A
    variance_gap_q = (full["variance"] - dominant["variance"]) / A
    third_gap_q = (full["third"] - dominant["third"]) / A ** mp.mpf("1.5")

    log_ratio_1_t = 2 * mean_gap_q
    log_ratio_2_t = 4 * variance_gap_q
    log_ratio_3_t = 8 * third_gap_q
    epsilon_1 = ratio_mass * log_ratio_1_t
    epsilon_2 = ratio_mass * (log_ratio_2_t + log_ratio_1_t**2)
    epsilon_3 = ratio_mass * (
        log_ratio_3_t
        + 3 * log_ratio_1_t * log_ratio_2_t
        + log_ratio_1_t**3
    )

    dominant_kappa3_t = 8 * dominant["third"] / A ** mp.mpf("1.5")
    full_kappa3_t = 8 * full["third"] / A ** mp.mpf("1.5")

    return {
        **data,
        "t": t,
        "epsilon": epsilon,
        "epsilon_1": epsilon_1,
        "epsilon_2": epsilon_2,
        "epsilon_3": epsilon_3,
        "log_correction_3": log_ratio_3_t,
        "dominant_kappa3": dominant_kappa3_t,
        "full_kappa3": full_kappa3_t,
        "scaled_epsilon_1": t * epsilon_1,
        "scaled_epsilon_2": t**2 * epsilon_2,
        "scaled_epsilon_3": t**3 * epsilon_3,
        "scaled_log_correction_3": t**3 * log_ratio_3_t,
        "scaled_dominant_kappa3": t**2 * dominant_kappa3_t,
        "scaled_full_kappa3": t**2 * full_kappa3_t,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", default="2,2.5,3,4,5,7,10,15,20")
    parser.add_argument("--dps", type=int, default=100)
    parser.add_argument("--digits", type=int, default=24)
    parser.add_argument("--window", default="24")
    parser.add_argument("--m-cutoff", type=int, default=12)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    width = mp.mpf(args.window)
    rows = []
    for text in args.samples.split(","):
        r = mp.mpf(text.strip())
        result = sample(r, width, args.m_cutoff)
        rows.append(
            {
                "r": mp.nstr(r, args.digits),
                **{
                    key: mp.nstr(value, args.digits)
                    for key, value in result.items()
                    if key
                    in {
                        "q",
                        "t",
                        "epsilon",
                        "scaled_epsilon_1",
                        "scaled_epsilon_2",
                        "scaled_epsilon_3",
                        "scaled_log_correction_3",
                        "scaled_dominant_kappa3",
                        "scaled_full_kappa3",
                    }
                },
            }
        )
    print(
        json.dumps(
            {
                "classification": "h1316_xi_transfer_scout_not_proof",
                "window": str(width),
                "m_cutoff": args.m_cutoff,
                "rows": rows,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
