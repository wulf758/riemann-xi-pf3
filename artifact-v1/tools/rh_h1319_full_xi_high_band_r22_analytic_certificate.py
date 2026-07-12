"""Canonical algebraic H1319 full-Xi high-band certificate from r=22.

This strengthens the factor-7 H921 fallback in the r=59 certificate.  An
elementary normalized calculation proves

    0 < q(r)^2*r*B(r)/A(r)^3 < 1                 (r>=4),

so the dominant scaled loss is at most ``2*D/r``.  H1317 supplies the global
full-Xi correction ``B_Xi/t`` for r>=5/2.  Exact rational endpoint arithmetic
then closes the coefficient 3/4 from r=22 onward.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import sympy as sp

import rh_h1319_full_xi_high_band_r59_analytic_certificate as fallback


ROOT = Path(__file__).resolve().parents[1]
RIEMANN = ROOT / "research" / "riemann"
DEFAULT_OUT = RIEMANN / "h1319_full_xi_high_band_r22_analytic_certificate.json"
DEFAULT_MD = RIEMANN / "h1319_full_xi_high_band_r22_analytic_certificate.md"

H1317_BETA5_NOTE = RIEMANN / "h1317_beta5_h923_derivative_certificate.md"
H1317_BETA5_TOOL = ROOT / "tools" / "rh_h1317_beta5_h923_derivative_certificate.py"
H1317_ASSEMBLY = RIEMANN / "h1317_effective_xi_transfer_assembly.md"

PAIR_CAP = Fraction(51, 100)
E2 = Fraction(26, 25)
E4 = Fraction(7, 2)
B_XI = 45_450_578_214
TARGET = Fraction(3, 4)
R_THRESHOLD = 22


def fraction_payload(value: Fraction) -> dict[str, Any]:
    return {
        "exact": str(value),
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def h1317_certificate() -> tuple[dict[str, Any], dict[str, Any]]:
    process = subprocess.run(
        [sys.executable, str(H1317_BETA5_TOOL)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    try:
        report = json.loads(process.stdout)
    except json.JSONDecodeError:
        report = {}
    note = H1317_BETA5_NOTE.read_text(encoding="utf-8")
    assembly = H1317_ASSEMBLY.read_text(encoding="utf-8")
    checks = {
        "verifier_returncode_zero": process.returncode == 0,
        "verifier_id_exact": report.get("id")
        == "h1317_beta5_h923_derivative_certificate",
        "verifier_all_checks_pass": report.get("all_checks_pass") is True,
        "verifier_BXi_exact": report.get("B_Xi") == B_XI,
        "global_range_in_note": "r>=5/2" in note,
        "global_R5_constants_in_note": all(
            marker in note for marker in ("207/t", "1035/t^2", "7659/t^3")
        ),
        "global_transfer_in_assembly": "For every `t` for which `r>=5/2`"
        in assembly,
        "assembly_BXi_exact": "45450578214/t" in assembly,
    }
    return checks, {
        "verifier_report": report,
        "verifier_stderr": process.stderr,
        "files": {
            str(H1317_BETA5_NOTE.relative_to(ROOT)): fallback.sha256(H1317_BETA5_NOTE),
            str(H1317_BETA5_TOOL.relative_to(ROOT)): fallback.sha256(H1317_BETA5_TOOL),
            str(H1317_ASSEMBLY.relative_to(ROOT)): fallback.sha256(H1317_ASSEMBLY),
        },
    }


def algebraic_ratio_certificate() -> dict[str, Any]:
    r, P, alpha, d, eps, x = sp.symbols(
        "r P alpha d eps x", positive=True
    )
    q = r * P - alpha * r - 1
    A = r * (P * (r + 1) - alpha)
    B = r * (P * (r**2 + 3 * r + 1) - alpha)
    ratio = q**2 * r * B / A**3
    normalized = (
        r
        * (1 - d - eps) ** 2
        * (r**2 + 3 * r + 1 - d)
        / (r + 1 - d) ** 3
    )
    normalized_residual = sp.factor(
        ratio
        - normalized.subs({d: alpha / P, eps: 1 / (r * P)})
    )

    comparison_difference = sp.expand(
        (r + 1 - d) ** 3 - r * (r**2 + 3 * r + 1 - d)
    )
    expected_difference = (
        2 * r
        + 1
        - d * (3 * r**2 + 5 * r + 3)
        + d**2 * (3 * r + 3)
        - d**3
    )
    difference_residual = sp.expand(comparison_difference - expected_difference)

    # For r>=4, pi*exp(r)>3*(1+r+r^2/2)>9r, hence
    # d=9/(4*pi*exp(r))<1/(4r).  The shifted guard is coefficient-positive.
    d_guard_polynomial = sp.expand(
        (sp.Rational(3, 2) * r**2 - 6 * r + 3).subs(r, x + 4)
    )
    d_guard_coefficients = sp.Poly(d_guard_polynomial, x).all_coeffs()

    # Dropping the positive d^2 term and using d<1/(4r) gives
    # E > L(r)=5r/4-1/4-3/(4r)-1/(64r^3).  Multiplying by 64r^3 and
    # shifting r=4+x produces a polynomial with positive coefficients.
    lower_polynomial = 80 * r**4 - 16 * r**3 - 48 * r**2 - 1
    shifted_lower = sp.expand(lower_polynomial.subs(r, x + 4))
    shifted_lower_coefficients = sp.Poly(shifted_lower, x).all_coeffs()

    checks = {
        "normalized_ratio_identity": normalized_residual == 0,
        "comparison_difference_identity": difference_residual == 0,
        "d_guard_shift_coefficients_positive": all(
            coefficient > 0 for coefficient in d_guard_coefficients
        ),
        "comparison_lower_shift_coefficients_positive": all(
            coefficient > 0 for coefficient in shifted_lower_coefficients
        ),
        "comparison_lower_at_4_exact": Fraction(18_687, 4_096) > 0,
    }
    return {
        "definitions": {
            "P": "pi*exp(r)",
            "q": "r*P-(9/4)*r-1",
            "A": "r*(P*(r+1)-9/4)",
            "B": "r*(P*(r^2+3r+1)-9/4)",
            "d": "9/(4P)",
            "epsilon": "1/(rP)",
            "F": "q^2*r*B/A^3",
        },
        "normalized_identity": (
            "F=r*(1-d-epsilon)^2*(r^2+3r+1-d)/(r+1-d)^3"
        ),
        "comparison_difference": str(expected_difference),
        "d_upper": "d<1/(4r) for r>=4",
        "d_guard_shift_polynomial": str(d_guard_polynomial),
        "comparison_lower": (
            "5r/4-1/4-3/(4r)-1/(64r^3)"
        ),
        "comparison_lower_shift_polynomial": str(shifted_lower),
        "comparison_lower_at_4": fraction_payload(Fraction(18_687, 4_096)),
        "conclusion": "0<q^2*r*B/A^3<1 for every real r>=4",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def exp_lower() -> Fraction:
    return sum(Fraction(1, math.factorial(degree)) for degree in range(7))


def q_lower_at_integer(r: int) -> Fraction:
    e_lower = Fraction(19, 7)
    return 3 * r * e_lower**r - Fraction(9, 4) * r - 1


def envelope_at_integer(r: int, D: Fraction) -> Fraction:
    return 2 * D / r + Fraction(2 * B_XI, 1) / q_lower_at_integer(r)


def build_report() -> dict[str, Any]:
    base_dependency_checks, base_provenance = fallback.dependency_checks()
    h1317_checks, h1317_provenance = h1317_certificate()
    ratio_certificate = algebraic_ratio_certificate()

    pair_y = PAIR_CAP * (E4 + 2 * E2)
    pair_x = PAIR_CAP * E2
    D = pair_y + 2 * pair_x**3

    e_taylor6 = exp_lower()
    q_lower = q_lower_at_integer(R_THRESHOLD)
    dominant_upper = 2 * D / R_THRESHOLD
    transfer_upper = Fraction(2 * B_XI, 1) / q_lower
    total_upper = dominant_upper + transfer_upper
    slack = TARGET - total_upper
    previous_total = envelope_at_integer(R_THRESHOLD - 1, D)
    first_passing = next(
        integer
        for integer in range(4, 72)
        if envelope_at_integer(integer, D) < TARGET
    )

    # A compact exact comparison also proves transfer<51/110, reproducing the
    # hand-auditable endpoint split 63/220+51/110=3/4.
    transfer_threshold = Fraction(220 * B_XI, 51)
    cleared_integer = (
        102
        * 7**22
        * (q_lower - transfer_threshold)
    )

    local_checks = {
        "D_exact": D == Fraction(6_141_071_619, 1_953_125_000),
        "D_lt_63_over_20": D < Fraction(63, 20),
        "e_taylor6_gt_19_over_7": e_taylor6 > Fraction(19, 7),
        "q22_lower_positive": q_lower > 0,
        "q22_clears_transfer_split": q_lower > transfer_threshold,
        "cleared_integer_exact": cleared_integer
        == 13_163_512_111_608_071_206_071_088_227_213,
        "dominant_split_lt_63_over_220": dominant_upper < Fraction(63, 220),
        "transfer_split_lt_51_over_110": transfer_upper < Fraction(51, 110),
        "split_sums_to_3_over_4": Fraction(63, 220) + Fraction(51, 110)
        == TARGET,
        "exact_total_below_target": total_upper < TARGET,
        "R21_same_coarse_envelope_fails": previous_total >= TARGET,
        "R22_first_integer_for_this_proof": first_passing == R_THRESHOLD,
        "R22_below_previous_R59": R_THRESHOLD < fallback.R_THRESHOLD,
        "R22_below_71": R_THRESHOLD < 71,
    }
    all_checks = (
        all(base_dependency_checks.values())
        and all(h1317_checks.values())
        and ratio_certificate["all_checks_pass"]
        and all(local_checks.values())
    )

    return {
        "id": "h1319_full_xi_high_band_r22_analytic_certificate",
        "classification": (
            "h1319_full_xi_high_band_r22_analytic_closed"
            if all_checks
            else "h1319_full_xi_high_band_r22_dependency_or_arithmetic_failure"
        ),
        "statement": {
            "r_domain": {"lower": "22", "upper": "infinity", "closed": True},
            "saddle_map": "t(r)=(pi*r*exp(r)-(9/4)r-1)/2",
            "quantity": "t(r)^2*kappa_Xi'''(t(r))",
            "lower_bound": "-3/4",
        },
        "dominant_stein_bound": {
            "identity": "kappa3(W)/a=-Y+3*X*E2-2*a^2*X^3",
            "pair_cap": str(PAIR_CAP),
            "X_upper": fraction_payload(pair_x),
            "Y_upper": fraction_payload(pair_y),
            "D": fraction_payload(D),
            "conclusion": "kappa3(W)>=-D*a",
        },
        "algebraic_ratio_certificate": ratio_certificate,
        "factor_audit": {
            "exact_scaled_dominant": "-2*D*q^2*B/A^3",
            "H921_fallback": "q^2*r*B/A^3<=7 gives -14D/r",
            "new_exact_algebra": "q^2*r*B/A^3<1 for r>=4",
            "improved_consequence": "t^2*kappa_9'''(t)>-2D/r",
            "factor_seven_not_silently_dropped": True,
        },
        "full_xi_transfer": {
            "dependency": "H1317 beta-5 H923 derivative certificate",
            "range": "r>=5/2",
            "B_Xi": B_XI,
            "scaled_correction": "t^2*G'''(t)>=-B_Xi/t=-2B_Xi/q",
            "checks": h1317_checks,
        },
        "threshold_budget": {
            "R": R_THRESHOLD,
            "e_degree6_taylor": fraction_payload(e_taylor6),
            "e_lower_used": "19/7",
            "q_R_lower": fraction_payload(q_lower),
            "dominant_upper": fraction_payload(dominant_upper),
            "transfer_upper": fraction_payload(transfer_upper),
            "total_upper": fraction_payload(total_upper),
            "target": fraction_payload(TARGET),
            "slack": fraction_payload(slack),
            "R_minus_1_total": fraction_payload(previous_total),
            "first_passing_integer_for_this_proof": first_passing,
            "hand_split": {
                "dominant": "<63/220",
                "transfer": "<51/110",
                "sum": "3/4",
                "cleared_integer": int(cleared_integer),
            },
        },
        "monotone_extension": {
            "dominant": "2D/r strictly decreases",
            "transfer": "q'(r)=pi*exp(r)*(r+1)-9/4=A/r>0",
            "conclusion": "the endpoint r=22 proves every r>=22",
        },
        "dependency_checks": base_dependency_checks,
        "local_checks": local_checks,
        "all_checks_pass": all_checks,
        "provenance": {
            "base_dependencies": base_provenance,
            "H1317": h1317_provenance,
        },
        "scope": {
            "proved": "full-Xi high analytic band r>=22",
            "remaining_compact_gap_after_first_segment": "4<r<22",
            "supersedes_high_band_fallback": "r>=59 certificate remains valid but is weaker",
            "does_not_prove": [
                "the remaining compact gap 4<r<22",
                "the global H920 input by itself",
                "the Riemann Hypothesis",
            ],
        },
    }


def build_markdown(report: dict[str, Any]) -> str:
    ratio = report["algebraic_ratio_certificate"]
    budget = report["threshold_budget"]
    lines = [
        "# H1319 Full-Xi Analytic High Band From r=22",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Result",
        "",
        "For every full-Xi saddle `r>=22`,",
        "",
        "```text",
        "t(r)^2 kappa_Xi'''(t(r)) >= -3/4.",
        "```",
        "",
        "This is an exact algebraic/rational certificate. It uses no finite",
        "differences and no interval quadrature.",
        "",
        "## Stein constant",
        "",
        "H1314, H1313 and H927 give",
        "",
        "```text",
        "kappa_3(W)>=-D*a,",
        f"D={report['dominant_stein_bound']['D']['exact']}",
        " =3.144228668928.",
        "```",
        "",
        "## Algebraic saddle-ratio certificate",
        "",
        "Put `P=pi*exp(r)`, `d=9/(4P)`, and `epsilon=1/(rP)`. Exact",
        "simplification gives",
        "",
        "```text",
        "F(r)=q^2*r*B/A^3",
        "    =r(1-d-epsilon)^2",
        "      (r^2+3r+1-d)/(r+1-d)^3.",
        "```",
        "",
        "For `r>=4`, the quadratic Taylor lower bound for `exp(r)` and",
        "`pi>3` give `pi*exp(r)>9r`, hence `d<1/(4r)` and `q>0`.",
        "Furthermore",
        "",
        "```text",
        "(r+1-d)^3-r(r^2+3r+1-d)",
        " =2r+1-d(3r^2+5r+3)+d^2(3r+3)-d^3",
        " >5r/4-1/4-3/(4r)-1/(64r^3).",
        "```",
        "",
        "After multiplication by `64r^3` and the shift `r=4+x`, the last",
        "lower bound becomes",
        "",
        "```text",
        f"{ratio['comparison_lower_shift_polynomial']},",
        "```",
        "",
        "whose coefficients are all positive. Therefore",
        "",
        "```text",
        "0<F(r)<1 for every r>=4.",
        "```",
        "",
        "The dominant scaled loss is consequently `2D/r`. H921's older",
        "factor-7 estimate would only give `14D/r`; that factor is explicitly",
        "retained in the fallback and replaced here by the proved inequality",
        "`F<1`.",
        "",
        "## Full-Xi transfer provenance",
        "",
        "The global correction bound is cited from H1317, not from the older",
        "restricted H1316 assembly:",
        "",
        "```text",
        "t^2 G'''(t)>=-B_Xi/t=-2B_Xi/q,",
        "B_Xi=45450578214,",
        "valid for every saddle r>=5/2.",
        "```",
        "",
        "## Endpoint r=22",
        "",
        "The degree-six Taylor polynomial proves `e>19/7`. Thus",
        "",
        "```text",
        f"q(22)>{budget['q_R_lower']['exact']}.",
        "```",
        "",
        "Exact rational arithmetic gives",
        "",
        "```text",
        f"dominant < {budget['dominant_upper']['exact']},",
        f"transfer < {budget['transfer_upper']['exact']},",
        f"total    < {budget['total_upper']['exact']}",
        f"         ~= {budget['total_upper']['decimal']}",
        "          < 3/4.",
        "```",
        "",
        f"The exact slack is `{budget['slack']['exact']}`. The same certified",
        "envelope fails at `r=21`, so 22 is the first passing integer for this",
        "proof. Since `2D/r` decreases and `q(r)` increases, the estimate holds",
        "for all `r>=22`.",
        "",
        "## Scope",
        "",
        "Combined with the first compact certificate through `r=4`, the",
        "remaining interval is now exactly",
        "",
        "```text",
        "4 < r < 22.",
        "```",
        "",
        "This does not close that interval and does not prove RH.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD))
    args = parser.parse_args()
    report = build_report()
    out = Path(args.out)
    md = Path(args.markdown_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    md.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md.write_text(build_markdown(report), encoding="utf-8")
    print(
        json.dumps(
            {
                "classification": report["classification"],
                "all_checks_pass": report["all_checks_pass"],
                "ratio_conclusion": report["algebraic_ratio_certificate"][
                    "conclusion"
                ],
                "R": report["threshold_budget"]["R"],
                "total_upper": report["threshold_budget"]["total_upper"],
                "slack": report["threshold_budget"]["slack"],
                "out": str(out),
                "markdown": str(md),
            },
            indent=2,
        )
    )
    return 0 if report["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
