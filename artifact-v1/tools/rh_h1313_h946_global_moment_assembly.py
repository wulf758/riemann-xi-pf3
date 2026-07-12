"""H1313 dependency assembly for the global H946 moment bounds.

This script certifies the canonical log-Gamma moments at s=30 with Arb and
checks the published dependency artifacts covering the three r-bands

    [5/2,3], [3,49/5], [49/5,infinity).

It does not rerun the expensive H1311 or H1279 interval covers.  It consumes
their checked artifacts and reruns the compact H1312 exact certificate.
"""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

from flint import arb, ctx

import rh_h1312_radial_lr_moment_certificate as h1312


ROOT = Path(__file__).resolve().parents[1]
RIEMANN = ROOT / "research" / "riemann"


def payload(value: arb, digits: int = 42) -> dict[str, str]:
    return {
        "interval": value.str(digits),
        "lower": value.lower().str(digits),
        "upper": value.upper().str(digits),
        "radius": value.rad().str(digits),
    }


def canonical_cumulants_and_moments() -> dict[str, Any]:
    """Arb cumulant/Bell certificate for W=sqrt(s)(log Y-log s), s=30."""

    ctx.prec = 256
    s = arb(30)
    root_s = s.sqrt()
    cumulants = [arb(0) for _ in range(7)]
    cumulants[1] = root_s * (s.digamma() - s.log())
    for order in range(2, 7):
        raw_polygamma = (
            ((-1) ** order)
            * math.factorial(order - 1)
            * arb(order).zeta(s)
        )
        cumulants[order] = s ** (arb(order) / 2) * raw_polygamma

    moments = [arb(0) for _ in range(7)]
    moments[0] = arb(1)
    for order in range(1, 7):
        moments[order] = sum(
            arb(math.comb(order - 1, block - 1))
            * cumulants[block]
            * moments[order - block]
            for block in range(1, order + 1)
        )

    targets = {
        2: arb(26) / 25,
        4: arb(7) / 2,
        6: arb(21),
    }
    rows: dict[str, Any] = {}
    for order, target in targets.items():
        slack = target - moments[order]
        rows[str(order)] = {
            "moment": payload(moments[order]),
            "target": payload(target),
            "slack": payload(slack),
            "passed": bool(moments[order].upper() < target),
        }

    return {
        "s": 30,
        "definition": "W=sqrt(s)*(log(Y)-log(s)), Y~Gamma(s,1)",
        "cumulants": {
            str(order): payload(cumulants[order]) for order in range(1, 7)
        },
        "moments": rows,
        "all_moments_pass": all(row["passed"] for row in rows.values()),
    }


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def note_attestation(path: Path, required_fragments: list[str]) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8") if path.exists() else ""
    fragments = {fragment: fragment in text for fragment in required_fragments}
    return {
        "path": str(path.relative_to(ROOT)),
        "exists": path.exists(),
        "required_fragments": fragments,
        "passed": path.exists() and all(fragments.values()),
    }


def dependencies() -> dict[str, Any]:
    h1282_path = (
        RIEMANN / "h1282_gamma_log_endpoint_even_moment_certificate_360k.json"
    )
    h1279_path = RIEMANN / "h1279_gamma_log_threshold49over5_moments.json"
    h1282 = load_json(h1282_path)
    h1279 = load_json(h1279_path)
    h1312_report = h1312.certificate()

    h1282_moments_pass = all(
        h1282["moments"][str(order)]["passed"] for order in (2, 4, 6)
    )
    h1311_note = note_attestation(
        RIEMANN / "h1311_h1296_finite_covariance_certificate_to_3.md",
        [
            "h1311_h1296_finite_covariance_certificate_to_3_proved",
            "5/2<=r<=3",
            "Hence the exact moments `M_2,M_4,M_6` are",
        ],
    )
    h1302_note = note_attestation(
        RIEMANN / "h1302_canonical_log_gamma_even_moment_monotonicity.md",
        [
            "h1302_canonical_log_gamma_all_even_moments_decrease",
            "strictly decreasing in `s`",
        ],
    )
    h1312_note = note_attestation(
        RIEMANN / "h1312_radial_likelihood_ratio_certificate.md",
        [
            "h1312_exact_radial_lr_domination",
            "3 <= r <= 49/5",
            "likelihood-ratio dominated",
        ],
    )

    h1311_drivers = [
        ROOT / "tools" / "rh_h1306b_finite_covariance_certificate_pilot.py",
        ROOT / "tools" / "rh_h1306c_paired_finite_covariance_pilot.py",
        ROOT / "tools" / "rh_h1306d_paired_cover_pilot.py",
    ]
    driver_status = {
        str(path.relative_to(ROOT)): path.exists() for path in h1311_drivers
    }

    return {
        "h1282_endpoint_anchor": {
            "path": str(h1282_path.relative_to(ROOT)),
            "all_checks_pass": bool(h1282.get("all_checks_pass")),
            "all_three_moments_pass": bool(h1282_moments_pass),
        },
        "h1311_low_band": {
            "note": h1311_note,
            "drivers": driver_status,
            "passed": h1311_note["passed"] and all(driver_status.values()),
            "role": "moment monotonicity on 5/2<=r<=3",
        },
        "h1302_canonical_monotonicity": h1302_note,
        "h1312_middle_band": {
            "note": h1312_note,
            "certificate_all_checks_pass": bool(h1312_report["all_checks_pass"]),
            "passed": h1312_note["passed"] and h1312_report["all_checks_pass"],
            "role": "exact-to-canonical radial LR comparison on 3<=r<=49/5",
        },
        "h1279_high_band": {
            "path": str(h1279_path.relative_to(ROOT)),
            "all_checks_pass": bool(h1279.get("all_checks_pass")),
            "moment_checks_pass": bool(h1279.get("checks", {}).get("moments_pass")),
            "role": "direct moment bounds on r>=49/5",
        },
    }


def coverage_certificate() -> dict[str, Any]:
    bands = [
        (Fraction(5, 2), Fraction(3), "H1311+H1282"),
        (Fraction(3), Fraction(49, 5), "H1312+H1302+canonical s=30"),
        (Fraction(49, 5), None, "H1279"),
    ]
    no_gap = bands[0][0] == Fraction(5, 2)
    for left, right in zip(bands, bands[1:]):
        no_gap = no_gap and left[1] == right[0]
    return {
        "domain": "r>=5/2",
        "bands": [
            {
                "left": str(left),
                "right": "infinity" if right is None else str(right),
                "dependency": dependency,
            }
            for left, right, dependency in bands
        ],
        "no_gap": bool(no_gap),
    }


def certificate() -> dict[str, Any]:
    canonical = canonical_cumulants_and_moments()
    deps = dependencies()
    coverage = coverage_certificate()
    dependency_pass = (
        deps["h1282_endpoint_anchor"]["all_checks_pass"]
        and deps["h1282_endpoint_anchor"]["all_three_moments_pass"]
        and deps["h1311_low_band"]["passed"]
        and deps["h1302_canonical_monotonicity"]["passed"]
        and deps["h1312_middle_band"]["passed"]
        and deps["h1279_high_band"]["all_checks_pass"]
        and deps["h1279_high_band"]["moment_checks_pass"]
    )
    all_pass = canonical["all_moments_pass"] and dependency_pass and coverage["no_gap"]
    return {
        "id": "h1313_h946_global_moment_assembly",
        "statement": {
            "range": "r>=5/2",
            "bounds": {"M2": "26/25", "M4": "7/2", "M6": "21"},
            "scope": "H946 coarse even-moment subproblem only",
        },
        "canonical_s30": canonical,
        "dependencies": deps,
        "coverage": coverage,
        "dependency_pass": bool(dependency_pass),
        "all_checks_pass": bool(all_pass),
        "does_not_prove": [
            "the H946 corrected-pair pointwise inequality",
            "H943-C",
            "the gamma-log-to-Xi transfer",
            "the Riemann Hypothesis",
        ],
        "limitation": "H1311 and H1279 are consumed as published checked artifacts, not rerun here",
    }


if __name__ == "__main__":
    print(json.dumps(certificate(), indent=2))
