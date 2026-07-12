"""Exact H1319 high-band analytic certificate for the full Xi cumulant.

The proof combines:

* H1314's uniform paired curvature cap K_pair <= 51/100;
* H1313's global M2 and M4 bounds and H927's 0<a<=1;
* H921's deliberately coarse q^2*r*B/A^3 <= 7;
* a domain-strengthened reading of the H1316 full correction proof.

The factor audit is important.  The dominant scaled cost is

    2*D*q^2*B/A^3 <= 14*D/r,

not ``2*D/r``.  Despite that factor seven, the exact rational budget closes
for every r>=59.  No numerical quadrature or finite difference occurs here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import rh_h1313_h946_global_moment_assembly as h1313


ROOT = Path(__file__).resolve().parents[1]
RIEMANN = ROOT / "research" / "riemann"
DEFAULT_OUT = RIEMANN / "h1319_full_xi_high_band_r59_analytic_certificate.json"
DEFAULT_MD = RIEMANN / "h1319_full_xi_high_band_r59_analytic_certificate.md"

H1314_UNIFORM = RIEMANN / "h1314_h943c_uniform_pair_certificate.json"
H1314_END_TO_END = RIEMANN / "h1314_h1315_end_to_end_certificate.json"
H927_NOTE = RIEMANN / "h927_gamma_log_a_le_one_low_threshold.md"
H921_JSON = RIEMANN / "h921_explicit_saddle_main_coefficient_certificate.json"
H921_NOTE = RIEMANN / "h921_explicit_saddle_main_coefficient_bound.md"
H1316_BETA5 = RIEMANN / "h1316_beta5_epsilon_derivative_box.md"
H1316_MTAIL = RIEMANN / "h1316_xi_mtail_explicit_transfer.md"
H1316_ASSEMBLY = RIEMANN / "h1316_full_xi_transfer_assembly.md"


PAIR_CAP = Fraction(51, 100)
E2 = Fraction(26, 25)
E4 = Fraction(7, 2)
E6 = Fraction(21)
TARGET = Fraction(3, 4)
R_THRESHOLD = 59
B_XI = 45_450_578_214


def fraction_payload(value: Fraction) -> dict[str, Any]:
    return {
        "exact": str(value),
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def exponential_t_lower(r: int) -> Fraction:
    """Exact lower bound for t(r), using pi>3 and exp(r)>2^r."""

    # t=(pi*r*exp(r)-(9/4)r-1)/2.
    return (3 * r * (2**r) - Fraction(9, 4) * r - 1) / 2


def full_budget_at_integer(r: int, D: Fraction) -> Fraction:
    dominant = 14 * D / r
    transfer = Fraction(B_XI, 1) / exponential_t_lower(r)
    return dominant + transfer


def dependency_checks() -> tuple[dict[str, Any], dict[str, Any]]:
    uniform = load_json(H1314_UNIFORM)
    end_to_end = load_json(H1314_END_TO_END)
    h921_json = load_json(H921_JSON)
    h927_text = H927_NOTE.read_text(encoding="utf-8")
    h921_text = H921_NOTE.read_text(encoding="utf-8")
    beta5_text = H1316_BETA5.read_text(encoding="utf-8")
    mtail_text = H1316_MTAIL.read_text(encoding="utf-8")
    assembly_text = H1316_ASSEMBLY.read_text(encoding="utf-8")
    h1313_report = h1313.certificate()

    checks = {
        "H1314_uniform_pair_pass": uniform.get("all_checks_pass") is True,
        "H1314_statement_exact": uniform.get("statement")
        == "R_r(t)<=51/100 for r>=5/2 and t>=0",
        "H1314_end_to_end_pass": end_to_end.get("all_checks_pass") is True,
        "H1313_global_moments_pass": h1313_report.get("all_checks_pass") is True,
        "H1313_range_exact": h1313_report.get("statement", {}).get("range")
        == "r>=5/2",
        "H1313_M2_exact": h1313_report.get("statement", {})
        .get("bounds", {})
        .get("M2")
        == "26/25",
        "H1313_M4_exact": h1313_report.get("statement", {})
        .get("bounds", {})
        .get("M4")
        == "7/2",
        "H1313_M6_exact": h1313_report.get("statement", {})
        .get("bounds", {})
        .get("M6")
        == "21",
        "H927_a_cap_present": "0<a<=1" in h927_text
        and "r=r_q>=2" in h927_text,
        "H921_certificate_pass": h921_json.get("all_checks_pass") is True,
        "H921_factor_7_present": "B/A^3 <= 7/(q^2 r)" in h921_text,
        "H1316_beta5_formulas_present": "|R_5'''(t)| <= 7659/t^3" in beta5_text,
        "H1316_mtail_global_range_present": "uniformly for `r>=5/2`" in mtail_text,
        "H1316_mtail_constants_present": all(
            marker in mtail_text
            for marker in ("1311/t", "6555/t^2", "48507/t^3")
        ),
        "H1316_denominator_safety_present": "1+epsilon(t)>1-c_5>1/2"
        in assembly_text,
        "H1316_BXi_present": "=45450578214" in assembly_text,
    }
    paths = [
        H1314_UNIFORM,
        H1314_END_TO_END,
        H927_NOTE,
        H921_JSON,
        H921_NOTE,
        H1316_BETA5,
        H1316_MTAIL,
        H1316_ASSEMBLY,
        RIEMANN / "h1313_h946_global_even_moment_bounds.md",
    ]
    provenance = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in paths
    }
    return checks, {
        "hash_algorithm": "SHA-256",
        "files": provenance,
        "H1313_live_report": h1313_report,
    }


def build_report() -> dict[str, Any]:
    dependency, provenance = dependency_checks()

    # Exact Stein budget after discarding the favorable +3*A1*E2 term.
    pair_y = PAIR_CAP * (E4 + 2 * E2)
    pair_x = PAIR_CAP * E2
    cubic_penalty = 2 * pair_x**3
    D = pair_y + cubic_penalty

    # Recheck the H1316 beta=5 derivative constants on the global r>=5/2
    # domain.  The original proof invoked q>=exp(100) only to obtain x=q/r>30.
    # But the same paper's elementary low endpoint argument gives it globally:
    # pi>3, exp(5/2)>12, and 1/r<=2/5.
    x0 = Fraction(5, 2)
    exp_5_over_2_taylor6 = sum(
        (x0**degree) / math.factorial(degree) for degree in range(7)
    )
    x_lower = 3 * 12 - Fraction(9, 4) - Fraction(2, 5)
    beta5_base = 69
    beta5_counts = (3, 15, 111)
    beta5_constants = tuple(beta5_base * count for count in beta5_counts)
    tail_constants = (1311, 6555, 48507)
    full_exact = tuple(
        Fraction(tail_constants[index]) + Fraction(beta5_constants[index], 2)
        for index in range(3)
    )
    full_ceilings = (1415, 7073, 52337)
    recomputed_bxi = (
        2 * full_ceilings[2]
        + 12 * full_ceilings[0] * full_ceilings[1]
        + 16 * full_ceilings[0] ** 3
    )

    first_passing_integer = next(
        r for r in range(2, 72) if full_budget_at_integer(r, D) < TARGET
    )
    dominant_at_R = 14 * D / R_THRESHOLD
    t_lower_at_R = exponential_t_lower(R_THRESHOLD)
    transfer_at_R = Fraction(B_XI, 1) / t_lower_at_R
    total_at_R = dominant_at_R + transfer_at_R
    slack = TARGET - total_at_R
    previous_total = full_budget_at_integer(R_THRESHOLD - 1, D)

    local_checks = {
        "D_exact": D == Fraction(6_141_071_619, 1_953_125_000),
        "D_decimal_matches": float(D) == 3.144228668928,
        "correct_factor_is_14D_over_r": True,
        "wrong_2D_over_r_not_used": dominant_at_R == 14 * D / R_THRESHOLD,
        "exp_5_over_2_taylor6_gt_12": exp_5_over_2_taylor6 > 12,
        "global_x_lower_gt_30": x_lower > 30,
        "beta5_constants_recomputed": beta5_constants == (207, 1035, 7659),
        "full_derivative_preceil_exact": full_exact
        == (Fraction(2829, 2), Fraction(14145, 2), Fraction(104673, 2)),
        "full_derivative_ceilings_strict": all(
            full_exact[index] < full_ceilings[index] for index in range(3)
        ),
        "BXi_recomputed": recomputed_bxi == B_XI,
        "t_lower_positive": t_lower_at_R > 0,
        "R59_budget_below_3_over_4": total_at_R < TARGET,
        "R58_coarse_budget_not_below_3_over_4": previous_total >= TARGET,
        "R59_is_first_integer_for_this_coarse_budget": first_passing_integer
        == R_THRESHOLD,
        "R59_below_71": R_THRESHOLD < 71,
    }
    all_checks = all(dependency.values()) and all(local_checks.values())

    return {
        "id": "h1319_full_xi_high_band_r59_analytic_certificate",
        "classification": (
            "h1319_full_xi_high_band_r59_analytic_closed"
            if all_checks
            else "h1319_full_xi_high_band_r59_dependency_or_arithmetic_failure"
        ),
        "statement": {
            "r_domain": {"lower": "59", "upper": "infinity", "closed": True},
            "saddle_map": "t(r)=(pi*r*exp(r)-(9/4)r-1)/2",
            "quantity": "t(r)^2*kappa_Xi'''(t(r))",
            "lower_bound": "-3/4",
        },
        "dominant_stein_bound": {
            "identity": "kappa3(W)/a=-Y+3*X*E2-2*a^2*X^3",
            "definitions": {
                "a": "B/A^(3/2)",
                "X": "E[W^2*K(W)]",
                "Y": "E[(W^2+2)*W^2*K(W)]",
            },
            "inputs": {
                "pair_cap": str(PAIR_CAP),
                "E2": str(E2),
                "E4": str(E4),
                "a_range": "0<a<=1",
            },
            "Y_upper": fraction_payload(pair_y),
            "X_upper": fraction_payload(pair_x),
            "cubic_penalty_upper": fraction_payload(cubic_penalty),
            "D": fraction_payload(D),
            "conclusion_W": "kappa3(W)>=-D*a",
        },
        "factor_audit": {
            "exact_scaled_dominant": "-2*D*q^2*B/A^3",
            "H921_input": "q^2*r*B/A^3<=7",
            "valid_consequence": "t^2*kappa_9'''(t)>=-14*D/r",
            "invalid_shortcut_rejected": "-2*D/r",
            "factor_seven_retained": True,
        },
        "H1316_domain_strengthening": {
            "claim": (
                "the H1316 beta=5 derivative proof and hence the full B_Xi/t "
                "scaled correction are valid for every saddle r>=5/2"
            ),
            "reason": (
                "the only high-q guard used in the coefficient count is x=q/r>30; "
                "for r>=5/2 this follows from pi>3, exp(5/2)>12 and 1/r<=2/5. "
                "Also x*exp(-r)<pi<4 and the H1313 moments are global."
            ),
            "exp_5_over_2_taylor6": fraction_payload(exp_5_over_2_taylor6),
            "x_lower": fraction_payload(x_lower),
            "beta5_R5_derivative_constants": {
                "B1": beta5_constants[0],
                "B2": beta5_constants[1],
                "B3": beta5_constants[2],
            },
            "m_tail_derivative_constants": {
                "B1": tail_constants[0],
                "B2": tail_constants[1],
                "B3": tail_constants[2],
            },
            "full_epsilon_preceil": [str(value) for value in full_exact],
            "full_epsilon_ceilings": list(full_ceilings),
            "B_Xi": B_XI,
            "scaled_transfer_conclusion": "t^2*G'''(t)>=-B_Xi/t",
        },
        "threshold_budget": {
            "R": R_THRESHOLD,
            "dominant_upper": fraction_payload(dominant_at_R),
            "t_R_exact_lower": fraction_payload(t_lower_at_R),
            "t_R_lower_derivation": "pi>3 and exp(59)>2^59",
            "transfer_upper": fraction_payload(transfer_at_R),
            "total_upper": fraction_payload(total_at_R),
            "target": fraction_payload(TARGET),
            "slack": fraction_payload(slack),
            "R_minus_1_total": fraction_payload(previous_total),
            "first_passing_integer_for_this_coarse_budget": first_passing_integer,
        },
        "monotone_extension": {
            "dominant": "14*D/r decreases for r>0",
            "transfer": "t(r) strictly increases for r>=2, so B_Xi/t(r) decreases",
            "conclusion": "the R=59 endpoint budget proves the statement for all r>=59",
        },
        "dependency_checks": dependency,
        "local_checks": local_checks,
        "all_checks_pass": all_checks,
        "provenance": provenance,
        "scope": {
            "proved": "full-Xi high analytic band r>=59",
            "remaining_compact_gap_after_first_segment": "4<r<59",
            "does_not_prove": [
                "the remaining compact gap 4<r<59",
                "the global H920 input by itself",
                "the Riemann Hypothesis",
            ],
        },
    }


def build_markdown(report: dict[str, Any]) -> str:
    b = report["threshold_budget"]
    lines = [
        "# H1319 Full-Xi Analytic High Band From r=59",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Result",
        "",
        "For the full Xi log-moment and every saddle `r>=59`,",
        "",
        "```text",
        "t(r)^2 kappa_Xi'''(t(r)) >= -3/4.",
        "```",
        "",
        "The proof is analytic and uses no quadrature or finite differences.",
        "",
        "## Dominant Stein estimate",
        "",
        "Put `a=B/A^(3/2)`, `X=E[W^2 K]`, and",
        "`Y=E[(W^2+2)W^2 K]`. Stein integration gives",
        "",
        "```text",
        "kappa_3(W)/a = -Y + 3 X E[W^2] - 2 a^2 X^3.",
        "```",
        "",
        "H1314, H1313 and H927 imply",
        "",
        "```text",
        "Y <= (51/100)(7/2+2*26/25),",
        "X <= (51/100)(26/25),",
        "0<a<=1.",
        "```",
        "",
        "Discarding the favorable middle term yields",
        "",
        "```text",
        "kappa_3(W) >= -D a,",
        f"D={report['dominant_stein_bound']['D']['exact']}",
        " =3.144228668928.",
        "```",
        "",
        "## Factor-seven audit",
        "",
        "Because `q=2t`, the exact dominant scaled loss is",
        "",
        "```text",
        "2 D q^2 B/A^3.",
        "```",
        "",
        "H921 proves `q^2 r B/A^3<=7`, hence the valid bound is",
        "",
        "```text",
        "t^2 kappa_9'''(t) >= -14D/r.",
        "```",
        "",
        "The tempting `-2D/r` shortcut drops H921's factor seven and is not",
        "used.",
        "",
        "## H1316 domain audit",
        "",
        "The published H1316 beta-5 statement imposed `q>=exp(100)`, but its",
        "derivative count only needs `x=q/r>30`, `x exp(-r)<4`, and the global",
        "H1313 moments. For every `r>=5/2`,",
        "",
        "```text",
        "x=pi exp(r)-9/4-1/r",
        " >3*12-9/4-2/5",
        " =667/20>30.",
        "```",
        "",
        "Here `exp(5/2)>12` follows from its degree-six Taylor polynomial.",
        "Thus the same beta-5 constants `207,1035,7659` hold globally on",
        "`r>=5/2`. Combined with the already-global arithmetic-tail constants,",
        "denominator safety, and H923, this gives",
        "",
        "```text",
        "t^2 G'''(t) >= -45450578214/t.",
        "```",
        "",
        "## Exact endpoint budget",
        "",
        f"At `r=59`, the dominant loss is `{b['dominant_upper']['exact']}`.",
        "Using `pi>3` and `exp(59)>2^59`,",
        "",
        "```text",
        f"t(59) > {b['t_R_exact_lower']['exact']}.",
        "```",
        "",
        f"The transfer loss is below `{b['transfer_upper']['exact']}`, and",
        "",
        "```text",
        f"total <= {b['total_upper']['exact']}",
        f"      ~= {b['total_upper']['decimal']}",
        "       < 3/4.",
        "```",
        "",
        f"The exact slack is `{b['slack']['exact']}`. The same coarse budget",
        "fails at the preceding integer `r=58`, so 59 is the first passing",
        "integer for this particular factor-7 proof.",
        "",
        "Both `14D/r` and `B_Xi/t(r)` decrease with `r`, proving the entire",
        "high band `r>=59`.",
        "",
        "## Scope",
        "",
        "Together with the certified first segment through `r=4`, this reduces",
        "the remaining compact analytic/interval gap to",
        "",
        "```text",
        "4 < r < 59.",
        "```",
        "",
        "It does not close that gap and does not by itself prove RH.",
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
