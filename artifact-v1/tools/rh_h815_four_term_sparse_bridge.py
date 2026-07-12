import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from rh_h807_essential_minor_hard_solver_audit import all_order3_pairs, determinant_terms
from rh_h812_sparse_row_pf3_obstruction_map import is_consecutive_rows
from rh_h813_sparse_row_symbolic_hierarchy import (
    Poly,
    analyze_pair,
    contiguous_d3_polys,
    determinant_q_poly,
    poly_to_string,
)


DEFAULT_OUT = "research/riemann/h815_four_term_sparse_bridge.json"


def zero_monomial(width: int) -> tuple[int, ...]:
    return tuple(0 for _ in range(width))


def quotient(numerator: tuple[int, ...], denominator: tuple[int, ...]) -> tuple[int, ...] | None:
    if any(left < right for left, right in zip(numerator, denominator)):
        return None
    return tuple(left - right for left, right in zip(numerator, denominator))


def monomial_text(monomial: tuple[int, ...]) -> str:
    factors = []
    for offset, exponent in enumerate(monomial):
        if exponent == 0:
            continue
        name = f"q{offset + 2}"
        factors.append(name if exponent == 1 else f"{name}^{exponent}")
    return "*".join(factors) if factors else "1"


def monotone_monomial_ge(left: tuple[int, ...], right: tuple[int, ...]) -> bool:
    """Return whether q^left >= q^right for 0<q2<=q3<=...<=1.

    With u_i=-log(q_i), the cone is u_2>=u_3>=...>=0.  The inequality
    q^left >= q^right is equivalent to sum (right-left)_i u_i >=0 on this
    cone, which holds iff every prefix sum of right-left is nonnegative.
    """

    prefix = 0
    for left_exp, right_exp in zip(left, right):
        prefix += right_exp - left_exp
        if prefix < 0:
            return False
    return True


def scalar_multiple(monomial: tuple[int, ...], scale: int) -> tuple[int, ...]:
    return tuple(scale * value for value in monomial)


def find_monotone_power_lower_bound(
    base: tuple[int, ...],
    compensator: tuple[int, ...],
) -> tuple[int, int] | None:
    candidates = [
        (2, 1),
        (7, 4),
        (5, 3),
        (3, 2),
        (4, 3),
        (5, 4),
        (1, 1),
        (3, 4),
        (2, 3),
        (1, 2),
        (1, 3),
        (1, 4),
    ]
    for numerator, denominator in candidates:
        if monotone_monomial_ge(
            scalar_multiple(compensator, denominator),
            scalar_multiple(base, numerator),
        ):
            return numerator, denominator
    return None


def certify_four_term_bridge(poly: Poly) -> dict[str, Any] | None:
    if not poly:
        return None
    width = len(next(iter(poly)))
    zero = zero_monomial(width)
    if poly.get(zero) != 1:
        return None
    negatives = [monomial for monomial, coeff in poly.items() if coeff == -1]
    positives = [monomial for monomial, coeff in poly.items() if coeff == 1 and monomial != zero]
    if len(poly) != 4 or len(negatives) != 2 or len(positives) != 1:
        return None

    z = positives[0]
    for x in negatives:
        for y in negatives:
            if y == x:
                continue
            # Need Y = X*q2*U with U<=1, hence Y<=X/2.
            y_over_x = quotient(y, x)
            if y_over_x is None or y_over_x[0] < 1:
                continue
            # Need Z = Y*W for a compensating monomial W.
            w = quotient(z, y)
            if w is None:
                continue
            beta = find_monotone_power_lower_bound(x, w)
            if y[0] >= 1 and beta is not None:
                numerator, denominator = beta
                return {
                    "type": "half_compensated_four_term_bridge",
                    "A": monomial_text(x),
                    "B": monomial_text(y),
                    "Z": monomial_text(z),
                    "W_Z_over_B": monomial_text(w),
                    "beta": f"{numerator}/{denominator}",
                    "proof": (
                        "P=1-A-B+Z=1-A-B*(1-W). Since B contains q2, B<=1/2. "
                        "Since W>=A^beta with beta<=2 under monotone q_i, "
                        "P>=1-A-(1/2)*(1-A^beta)>=1/2-A+(1/2)*A^2>=0."
                    ),
                }
            # Need W >= X^(3/2), equivalently W^2 >= X^3.
            if not monotone_monomial_ge(scalar_multiple(w, 2), scalar_multiple(x, 3)):
                continue
            return {
                "type": "four_term_sparse_bridge",
                "X": monomial_text(x),
                "Y": monomial_text(y),
                "Z": monomial_text(z),
                "Y_over_X": monomial_text(y_over_x),
                "W_Z_over_Y": monomial_text(w),
                "proof": (
                    "P=1-X-Y+Z=1-X-Y*(1-W). Since Y=X*q2*U with U<=1, "
                    "Y<=X/2. Since W^2>=X^3 under monotone q_i, W>=X^(3/2). "
                    "Thus P>=1-X-(X/2)*(1-X^(3/2))>=0 for 0<=X<=1."
                ),
            }
    return None


def mixed_unclassified_pairs(length: int) -> list[dict[str, Any]]:
    q_max = length + 1
    sequence_len = length + 2
    d3_index = contiguous_d3_polys(q_max, sequence_len)
    out = []
    for rows, cols in all_order3_pairs(sequence_len):
        if is_consecutive_rows(rows):
            continue
        if not determinant_terms(rows, cols, sequence_len):
            continue
        analysis = analyze_pair(rows, cols, q_max, d3_index)
        cert = analysis["certificate"]
        if analysis["matches_translated_contiguous_d3"] is not None:
            continue
        if cert["type"] != "mixed":
            continue
        if any(item["quotient_nonnegative_coefficients"] for item in cert.get("one_minus_factor_candidates", [])):
            continue
        poly, _ = determinant_q_poly(rows, cols, q_max)
        out.append({"rows": rows, "cols": cols, "analysis": analysis, "poly": poly})
    return out


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    pairs = mixed_unclassified_pairs(args.length)
    covered = []
    failures = []
    not_four_term = []
    for item in pairs:
        certificate = certify_four_term_bridge(item["poly"])
        if certificate is None:
            if item["analysis"]["term_count"] == 4:
                failures.append(item)
            else:
                not_four_term.append(item)
            continue
        covered.append({**item["analysis"], "h815_certificate": certificate})

    family_counter = Counter(item["family"] for item in covered)
    failure_family_counter = Counter(item["analysis"]["family"] for item in failures + not_four_term)

    if failures:
        classification = "h815_four_term_sparse_bridge_partial_with_failures"
        decision = (
            "The four-term bridge proves many H813 mixed forms but not every four-term "
            "case. Inspect the failed four-term cases before promoting this to a lemma."
        )
    else:
        classification = "h815_four_term_sparse_bridge_supported"
        decision = (
            "Every H813 mixed four-term sparse-row form is certified by the bridge. "
            "The remaining H813 hard core is now the higher-term boundary/D3 family."
        )

    return {
        "schema": "rh_h815_four_term_sparse_bridge.v0",
        "classification": classification,
        "lemma": {
            "name": "Four-term sparse bridge",
            "statement": (
                "For 0<q2<=q3<=...<=1 with q2<=1/2, certify normalized determinant "
                "polynomials P=1-A-B+Z by either: (i) B contains q2, Z=B*W, and "
                "W>=A^beta for some beta<=2; or (ii) the stronger subcase B=A*q2*U, "
                "U<=1, Z=B*W, and W^2>=A^3."
            ),
            "scalar_reduction": (
                "case (i): P>=1/2-A+(1/2)A^beta>=1/2-A+(1/2)A^2>=0; "
                "case (ii): P>=1-A-(A/2)*(1-A^(3/2))>=0."
            ),
        },
        "counts": {
            "h813_mixed_unclassified_input": len(pairs),
            "covered_by_four_term_bridge": len(covered),
            "four_term_failures": len(failures),
            "not_four_term_remaining": len(not_four_term),
        },
        "covered_family_counts": family_counter.most_common(),
        "remaining_family_counts": failure_family_counter.most_common(),
        "covered_examples": covered[: args.example_limit],
        "four_term_failure_examples": [item["analysis"] for item in failures[: args.example_limit]],
        "not_four_term_examples": [item["analysis"] for item in not_four_term[: args.example_limit]],
        "decision": decision,
        "next_target": {
            "name": "H816 six-term boundary/D3 coupling",
            "statement": (
                "Use the translated contiguous D3 input to certify the remaining higher-term "
                "mixed sparse-row polynomials, or pivot to a Plucker-Dodgson identity if "
                "direct product bounds stall."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# H815 Four-Term Sparse Bridge",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Lemma",
        "",
        report["lemma"]["statement"],
        "",
        f"Scalar reduction: `{report['lemma']['scalar_reduction']}`",
        "",
        "## Counts",
        "",
        "```json",
        json.dumps(report["counts"], indent=2),
        "```",
        "",
        "## Covered Families",
        "",
    ]
    for family, count in report["covered_family_counts"][:20]:
        lines.append(f"- `{family}`: {count}")
    lines.extend(["", "## Covered Examples", ""])
    for item in report["covered_examples"][:8]:
        lines.extend(
            [
                f"### rows {item['rows']} cols {item['cols']}",
                "",
                f"Polynomial: `{item['polynomial']}`",
                "",
                "Certificate:",
                "",
                "```json",
                json.dumps(item["h815_certificate"], indent=2),
                "```",
                "",
            ]
        )
    lines.extend(["## Remaining Families", ""])
    for family, count in report["remaining_family_counts"][:20]:
        lines.append(f"- `{family}`: {count}")
    lines.extend(["", "## Decision", "", report["decision"], "", "## Next Target", ""])
    lines.append(f"`{report['next_target']['name']}`: {report['next_target']['statement']}")
    lines.append("")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="H815 four-term sparse bridge audit.")
    parser.add_argument("--length", type=int, default=7)
    parser.add_argument("--example-limit", type=int, default=12)
    parser.add_argument("--out", default=DEFAULT_OUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = build_report(args)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    write_markdown(report, out_path.with_suffix(".md"))
    print(
        json.dumps(
            {
                "classification": report["classification"],
                "counts": report["counts"],
                "next_target": report["next_target"]["name"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()


