"""Validate the H1319 first-segment certificate and emit a non-global manifest.

This manifest is deliberately incompatible with the global compact schema used
by the H920 assembly.  It records the proved segment through r=4 and keeps the
remaining r in (4,71) machine-visible.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

from flint import arb, ctx


DEFAULT_DETAIL = Path("research/riemann/h1319_full_xi_adaptive_rbox_cover.json")
DEFAULT_OUT = Path("research/riemann/h1319_full_xi_segment_interval_certificate.json")
SCHEMA = "rh_h1319_full_xi_segment_interval_certificate.v1"


def qarb(value: Fraction) -> arb:
    return arb(value.numerator) / arb(value.denominator)


def ball(lo: Fraction, hi: Fraction) -> arb:
    left = qarb(lo)
    right = qarb(hi)
    return (left + right) / 2 + arb(0, 1) * (right - left) / 2


def interval_payload(value: arb, digits: int = 30) -> dict[str, str]:
    return {
        "interval": value.str(digits),
        "lower": value.lower().str(digits),
        "upper": value.upper().str(digits),
    }


def build_manifest(detail_path: Path) -> dict[str, Any]:
    raw = detail_path.read_bytes()
    detail = json.loads(raw)
    boxes = detail.get("certified_boxes", [])
    exact_boxes = [
        (Fraction(row["r_lo"]), Fraction(row["r_hi"])) for row in boxes
    ]
    expected_lo = Fraction(323, 100)
    expected_hi = Fraction(4)
    contiguous = bool(exact_boxes) and all(
        left[1] == right[0] for left, right in zip(exact_boxes, exact_boxes[1:])
    )

    r_box = ball(expected_lo, expected_hi)
    q_prime = arb.pi() * r_box.exp() * (r_box + 1) - arb(9) / 4
    t_prime = q_prime / 2
    t_at_4 = 2 * arb.pi() * arb(4).exp() - 5

    bracket = detail.get("r125_bracket", {})
    bracket_values = bracket.get("bracket", [])
    checks = {
        "detailed_classification_certified": detail.get("classification")
        == "h1319_full_xi_adaptive_rbox_cover_certified",
        "detailed_all_certified": detail.get("summary", {}).get("all_certified")
        is True,
        "detailed_failed_boxes_empty": detail.get("failed_boxes") == [],
        "box_count_is_192": len(boxes) == 192,
        "all_box_statuses_certified": all(
            row.get("status") == "certified" for row in boxes
        ),
        "all_target_polynomial_lowers_positive": all(
            float(row.get("P_lower_float", float("nan"))) > 0 for row in boxes
        ),
        "exact_cover_starts_at_323_over_100": bool(exact_boxes)
        and exact_boxes[0][0] == expected_lo,
        "exact_cover_ends_at_4": bool(exact_boxes)
        and exact_boxes[-1][1] == expected_hi,
        "exact_cover_contiguous": contiguous,
        "h1271_h1316_dependencies_pass": detail.get("dependencies", {}).get(
            "all_pass"
        )
        is True,
        "r125_bracket_passes": bracket.get("all_pass") is True,
        "r125_bracket_starts_at_cover_lower": len(bracket_values) == 2
        and Fraction(bracket_values[0]) == expected_lo,
        "t_strictly_increasing_on_cover": bool(t_prime.lower() > 0),
        "quantity_is_full_xi_third_cumulant": detail.get("identity", {}).get(
            "kappa_third"
        )
        == "8*C/(A^(3/2)*K0^3)",
        "target_polynomial_is_recorded": detail.get("identity", {}).get(
            "target_polynomial"
        )
        == "8*t^2*C+beta*A^(3/2)*K0^3",
        "target_beta_is_3_over_4": detail.get("config", {}).get("beta") == "3/4",
    }
    all_segment_checks_pass = all(checks.values())
    return {
        "schema": SCHEMA,
        "classification": (
            "h1319_full_xi_first_compact_segment_certified"
            if all_segment_checks_pass
            else "h1319_full_xi_first_compact_segment_manifest_failure"
        ),
        "claim": {
            "kernel": "full_xi",
            "quantity": "t^2*kappa_Xi'''(t)",
            "lower_bound": "-3/4",
            "r_cover_domain": {
                "lower": "323/100",
                "upper": "4",
                "closed": True,
            },
            "desired_r_segment": {
                "lower": "r(125)",
                "upper": "4",
                "closed": True,
            },
            "t_domain": {
                "lower": "125",
                "upper": "2*pi*exp(4)-5",
                "closed": True,
            },
        },
        "endpoint_checks": {
            "r125_bracket": bracket,
            "t_prime_on_r_cover": interval_payload(t_prime),
            "t_at_4": interval_payload(t_at_4),
        },
        "checks": checks,
        "all_segment_checks_pass": all_segment_checks_pass,
        "global_compact_status": {
            "global_schema": "rh_h1319_full_xi_compact_interval_certificate.v1",
            "global_t_domain": {"lower": "125", "upper": "T_71", "closed": True},
            "full_domain_covered": False,
            "global_manifest_eligible": False,
            "remaining_r_gap": {"lower": "4", "upper": "71", "open": True},
            "reason": "the certified first segment ends at r=4, while H1317 starts at r=71",
        },
        "evidence": {
            "detailed_certificate": str(detail_path),
            "detailed_certificate_sha256": hashlib.sha256(raw).hexdigest().upper(),
            "canonical_generator": "tools/rh_h1319_full_xi_adaptive_rbox_certificate.py",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--detail", default=str(DEFAULT_DETAIL))
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--dps", type=int, default=80)
    args = parser.parse_args()
    ctx.dps = args.dps
    manifest = build_manifest(Path(args.detail))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "schema": manifest["schema"],
                "classification": manifest["classification"],
                "all_segment_checks_pass": manifest["all_segment_checks_pass"],
                "full_domain_covered": manifest["global_compact_status"][
                    "full_domain_covered"
                ],
                "out": str(out),
            },
            indent=2,
        )
    )
    return 0 if manifest["all_segment_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
