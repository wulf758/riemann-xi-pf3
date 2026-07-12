#!/usr/bin/env python3
"""Canonical H1319 compact manifest entry point with analytic saddle guard.

The underlying manifest assembler validates the five proof components.  This
entry point replaces its deliberately naive wide-ball q' diagnostic by the
rigorous monotonicity argument

    q''(r) = pi*exp(r)*(r+2) > 0,

so q'(r) on r>=323/100 is bounded below by its Arb-certified left endpoint.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx

import rh_h1319_full_xi_compact_interval_manifest as implementation


def qarb(value: Fraction) -> arb:
    return arb(value.numerator) / arb(value.denominator)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--first", default=str(implementation.DEFAULT_FIRST))
    parser.add_argument("--mean-4-5", default=str(implementation.DEFAULT_MEAN_4_5))
    parser.add_argument("--mean-5-10", default=str(implementation.DEFAULT_MEAN_5_10))
    parser.add_argument("--mean-10-22", default=str(implementation.DEFAULT_MEAN_10_22))
    parser.add_argument("--high", default=str(implementation.DEFAULT_HIGH))
    parser.add_argument("--out", default=str(implementation.DEFAULT_OUT))
    parser.add_argument("--markdown-out", default=str(implementation.DEFAULT_MARKDOWN))
    parser.add_argument("--dps", type=int, default=90)
    args = parser.parse_args()
    ctx.dps = args.dps

    report = implementation.build_report(args)
    left = qarb(Fraction(323, 100))
    q_prime_left = arb.pi() * left.exp() * (left + 1) - arb(9) / 4
    report["saddle_checks"].pop("q_prime_on_323_over_100_to_71", None)
    report["saddle_checks"].update(
        {
            "q_prime_at_323_over_100": implementation.interval_payload(
                q_prime_left
            ),
            "q_second_derivative_sign": "pi*exp(r)*(r+2)>0 for r>0",
            "monotonicity_conclusion": (
                "q'(r)>=q'(323/100)>0 for every r>=323/100"
            ),
        }
    )
    report["checks"]["saddle_map_covers_t_domain"] = bool(
        q_prime_left.lower() > 0
    )
    report["all_checks_pass"] = all(report["checks"].values())
    report["classification"] = (
        "h1319_full_xi_compact_interval_closed"
        if report["all_checks_pass"]
        else "h1319_full_xi_compact_interval_manifest_failure"
    )

    out = Path(args.out)
    markdown = Path(args.markdown_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    markdown.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    markdown.write_text(implementation.build_markdown(report), encoding="utf-8")
    print(
        json.dumps(
            {
                "schema": report["schema"],
                "classification": report["classification"],
                "all_checks_pass": report["all_checks_pass"],
                "segments": len(report["r_proof_chain"]["segments"]),
                "q_prime_left_lower": q_prime_left.lower().str(30),
                "out": str(out),
                "markdown": str(markdown),
            },
            indent=2,
        )
    )
    return 0 if report["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
