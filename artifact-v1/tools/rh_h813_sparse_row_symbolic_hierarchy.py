import argparse
import json
import math
import sys
from fractions import Fraction
from itertools import permutations
from pathlib import Path
from typing import Any

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from rh_h807_essential_minor_hard_solver_audit import all_order3_pairs, determinant_terms
from rh_h812_sparse_row_pf3_obstruction_map import is_consecutive_rows


DEFAULT_H812 = "research/riemann/h812_sparse_row_pf3_obstruction_map.json"
DEFAULT_OUT = "research/riemann/h813_sparse_row_symbolic_hierarchy.json"

Poly = dict[tuple[int, ...], int]


def add_term(poly: Poly, monomial: tuple[int, ...], coeff: int) -> None:
    if coeff == 0:
        return
    value = poly.get(monomial, 0) + coeff
    if value:
        poly[monomial] = value
    elif monomial in poly:
        del poly[monomial]


def add_poly(left: Poly, right: Poly, scale: int = 1) -> Poly:
    out = dict(left)
    for monomial, coeff in right.items():
        add_term(out, monomial, scale * coeff)
    return out


def monomial_mul(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(a + b for a, b in zip(left, right))


def scalar_normalize(poly: Poly) -> Poly:
    if not poly:
        return {}
    gcd = 0
    for coeff in poly.values():
        gcd = math.gcd(gcd, abs(coeff))
    if gcd <= 1:
        return poly
    return {monomial: coeff // gcd for monomial, coeff in poly.items()}


def positive_monomial_normalize(poly: Poly) -> tuple[Poly, tuple[int, ...]]:
    if not poly:
        return {}, ()
    width = len(next(iter(poly)))
    mins = tuple(min(monomial[index] for monomial in poly) for index in range(width))
    reduced = {
        tuple(exponent - mins[index] for index, exponent in enumerate(monomial)): coeff
        for monomial, coeff in poly.items()
    }
    return scalar_normalize(reduced), mins


def q_monomial_for_a(index: int, q_max: int) -> tuple[int, ...] | None:
    if index < 0:
        return None
    return tuple(max(0, index - q_index + 1) for q_index in range(2, q_max + 1))


def determinant_q_poly(rows: tuple[int, ...], cols: tuple[int, ...], q_max: int) -> tuple[Poly, tuple[int, ...]]:
    signs = {
        (0, 1, 2): 1,
        (0, 2, 1): -1,
        (1, 0, 2): -1,
        (1, 2, 0): 1,
        (2, 0, 1): 1,
        (2, 1, 0): -1,
    }
    poly: Poly = {}
    for perm in permutations(range(3)):
        monomial = tuple(0 for _ in range(2, q_max + 1))
        valid = True
        for row_position, col_position in enumerate(perm):
            index = cols[col_position] - rows[row_position]
            entry = q_monomial_for_a(index, q_max)
            if entry is None:
                valid = False
                break
            monomial = monomial_mul(monomial, entry)
        if valid:
            add_term(poly, monomial, signs[perm])
    return positive_monomial_normalize(poly)


def poly_key(poly: Poly) -> tuple[tuple[tuple[int, ...], int], ...]:
    return tuple(sorted(poly.items()))


def poly_to_string(poly: Poly, q_max: int, limit: int = 12) -> str:
    if not poly:
        return "0"
    pieces = []
    for count, (monomial, coeff) in enumerate(sorted(poly.items())):
        if count >= limit:
            pieces.append("...")
            break
        factors = []
        for offset, exponent in enumerate(monomial):
            if exponent == 0:
                continue
            q_name = f"q{offset + 2}"
            factors.append(q_name if exponent == 1 else f"{q_name}^{exponent}")
        body = "*".join(factors) if factors else "1"
        if coeff == 1:
            pieces.append(body)
        elif coeff == -1:
            pieces.append(f"-{body}")
        else:
            pieces.append(f"{coeff}*{body}")
    return " + ".join(pieces).replace("+ -", "- ")


def substitute_var_one(poly: Poly, var: int) -> Poly:
    out: Poly = {}
    for monomial, coeff in poly.items():
        reduced = list(monomial)
        reduced[var - 2] = 0
        add_term(out, tuple(reduced), coeff)
    return out


def substitute_var_equal(poly: Poly, target: int, source: int) -> Poly:
    out: Poly = {}
    for monomial, coeff in poly.items():
        reduced = list(monomial)
        reduced[source - 2] += reduced[target - 2]
        reduced[target - 2] = 0
        add_term(out, tuple(reduced), coeff)
    return out


def divide_by_one_minus(poly: Poly, var: int) -> Poly | None:
    grouped: dict[tuple[int, ...], dict[int, int]] = {}
    var_offset = var - 2
    width = len(next(iter(poly))) if poly else 0
    for monomial, coeff in poly.items():
        other = tuple(value for index, value in enumerate(monomial) if index != var_offset)
        grouped.setdefault(other, {})[monomial[var_offset]] = coeff

    quotient: Poly = {}
    for other, coeffs in grouped.items():
        max_exp = max(coeffs)
        if sum(coeffs.values()) != 0:
            return None
        running = 0
        for exponent in range(max_exp):
            running += coeffs.get(exponent, 0)
            if running:
                full = []
                other_iter = iter(other)
                for index in range(width):
                    full.append(exponent if index == var_offset else next(other_iter))
                add_term(quotient, tuple(full), running)
        if coeffs.get(max_exp, 0) != -running:
            return None
    return scalar_normalize(quotient)


def nonnegative_coefficients(poly: Poly) -> bool:
    return bool(poly) and all(coeff >= 0 for coeff in poly.values())


def nonpositive_coefficients(poly: Poly) -> bool:
    return bool(poly) and all(coeff <= 0 for coeff in poly.values())


def q_box_lower_bound(poly: Poly) -> Fraction:
    lower = Fraction(0)
    for monomial, coeff in poly.items():
        if coeff >= 0:
            if all(exponent == 0 for exponent in monomial):
                lower += coeff
            continue
        q2_exponent = monomial[0] if monomial else 0
        lower -= Fraction(abs(coeff), 2**q2_exponent)
    return lower


def sign_certificate(poly: Poly, q_max: int) -> dict[str, Any]:
    if not poly:
        return {"type": "zero"}
    if nonnegative_coefficients(poly):
        return {"type": "nonnegative_coefficients"}
    if nonpositive_coefficients(poly):
        return {"type": "nonpositive_coefficients", "warning": "orientation is reversed"}
    box_lower = q_box_lower_bound(poly)
    if box_lower >= 0:
        return {
            "type": "q_box_lower_bound",
            "lower_bound": str(box_lower),
            "inputs": "0<q_i<=1 and q2<=1/2; positive nonconstant monomials discarded",
        }

    one_minus = []
    for var in range(2, q_max + 1):
        quotient = divide_by_one_minus(poly, var)
        if quotient is not None:
            one_minus.append(
                {
                    "factor": f"1-q{var}",
                    "quotient_nonnegative_coefficients": nonnegative_coefficients(quotient),
                    "quotient_terms": len(quotient),
                    "quotient": poly_to_string(quotient, q_max, limit=8),
                }
            )

    monotone_equalities = []
    for target in range(3, q_max + 1):
        source = target - 1
        if not substitute_var_equal(poly, target, source):
            monotone_equalities.append(f"q{target}=q{source}")

    return {
        "type": "mixed",
        "one_minus_factor_candidates": one_minus,
        "adjacent_monotone_equality_vanishes": monotone_equalities,
    }


def load_h812_top_pairs(path: Path) -> list[tuple[tuple[int, ...], tuple[int, ...], str]]:
    h812 = json.loads(path.read_text(encoding="utf-8"))
    out = []
    for item in h812["grid"]["stats"].get("top_sparse_family_minima", []):
        out.append((tuple(item["rows"]), tuple(item["cols"]), item["family"]))
    return out


def contiguous_d3_polys(q_max: int, sequence_len: int) -> dict[tuple[tuple[tuple[int, ...], int], ...], dict[str, Any]]:
    out = {}
    for start in range(0, sequence_len - 2):
        rows = (0, 1, 2)
        cols = (start, start + 1, start + 2)
        terms = determinant_terms(rows, cols, sequence_len)
        if not terms:
            continue
        poly, monomial = determinant_q_poly(rows, cols, q_max)
        out[poly_key(poly)] = {
            "rows": list(rows),
            "cols": list(cols),
            "removed_monomial": list(monomial),
            "polynomial": poly_to_string(poly, q_max),
        }
    return out


def analyze_pair(
    rows: tuple[int, ...],
    cols: tuple[int, ...],
    q_max: int,
    d3_index: dict[tuple[tuple[tuple[int, ...], int], ...], dict[str, Any]],
    family: str | None = None,
) -> dict[str, Any]:
    poly, monomial = determinant_q_poly(rows, cols, q_max)
    certificate = sign_certificate(poly, q_max)
    d3_match = d3_index.get(poly_key(poly))
    return {
        "rows": list(rows),
        "cols": list(cols),
        "family": family or f"row_gaps=({rows[1]-rows[0]},{rows[2]-rows[1]}), col_gaps=({cols[1]-cols[0]},{cols[2]-cols[1]})",
        "removed_positive_monomial_exponents_q2_to_qmax": list(monomial),
        "term_count": len(poly),
        "polynomial": poly_to_string(poly, q_max),
        "certificate": certificate,
        "matches_translated_contiguous_d3": d3_match,
    }


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    q_max = args.length + 1
    sequence_len = args.length + 2
    d3_index = contiguous_d3_polys(q_max, sequence_len)
    top_pairs = load_h812_top_pairs(Path(args.h812))
    top_analysis = [analyze_pair(rows, cols, q_max, d3_index, family) for rows, cols, family in top_pairs]

    active_pairs = []
    for rows, cols in all_order3_pairs(sequence_len):
        if is_consecutive_rows(rows):
            continue
        if not determinant_terms(rows, cols, sequence_len):
            continue
        active_pairs.append((rows, cols))

    summary_counts = {
        "active_sparse_row_pairs": len(active_pairs),
        "zero_polynomials": 0,
        "nonnegative_coefficients": 0,
        "q_box_lower_bound": 0,
        "one_minus_factor_with_nonnegative_quotient": 0,
        "matches_translated_contiguous_d3": 0,
        "mixed_unclassified": 0,
    }
    family_examples = []
    for rows, cols in active_pairs:
        analysis = analyze_pair(rows, cols, q_max, d3_index)
        ctype = analysis["certificate"]["type"]
        if ctype == "zero":
            summary_counts["zero_polynomials"] += 1
        elif ctype == "nonnegative_coefficients":
            summary_counts["nonnegative_coefficients"] += 1
        elif ctype == "q_box_lower_bound":
            summary_counts["q_box_lower_bound"] += 1
        elif analysis["matches_translated_contiguous_d3"] is not None:
            summary_counts["matches_translated_contiguous_d3"] += 1
        elif any(
            item["quotient_nonnegative_coefficients"]
            for item in analysis["certificate"].get("one_minus_factor_candidates", [])
        ):
            summary_counts["one_minus_factor_with_nonnegative_quotient"] += 1
        else:
            summary_counts["mixed_unclassified"] += 1
            if len(family_examples) < args.example_limit:
                family_examples.append(analysis)

    if summary_counts["mixed_unclassified"]:
        classification = "h813_sparse_row_symbolic_hierarchy_partial"
        decision = (
            "The H812 sparse-row boundary has simple exact factors, but the active sparse-row "
            "universe is not yet reduced to elementary coefficient/factor certificates. "
            "The route needs a stronger symbolic identity, likely a Plucker/Dodgson or "
            "compound-matrix inequality using H811 minors."
        )
    else:
        classification = "h813_sparse_row_symbolic_hierarchy_candidate_complete"
        decision = (
            "All active sparse-row pairs in the audited window reduce to elementary "
            "nonnegative coefficient/factor/D3 certificates. This would be a strong "
            "candidate for a continuum proof, pending formalization."
        )

    alternate_idea_register = [
        {
            "route": "Xi kernel saddle-point unit-deficit",
            "source": "H804/H805",
            "why_not_exhausted": "The clean target T_{n+1}/T_n >= 1-1/n^2 remains analytic, not grid-mined.",
        },
        {
            "route": "Li coefficient interval certification",
            "source": "H212-H217",
            "why_not_exhausted": "Blocked by rigorous outer-circle zeta'/zeta bounds, but could revive with stronger interval tooling.",
        },
        {
            "route": "Endpoint/Nyman-Beurling profile lemmas",
            "source": "H574/H611",
            "why_not_exhausted": "Formal candidates survive after several invalidations; they are independent of the PF3 corridor.",
        },
    ]
    return {
        "schema": "rh_h813_sparse_row_symbolic_hierarchy.v0",
        "classification": classification,
        "inputs": {
            "h812": args.h812,
            "length": args.length,
            "q_variables": [f"q{index}" for index in range(2, q_max + 1)],
        },
        "method": [
            "Represent a_k/a_0 after removing the common positive R1^k factor as a monomial in q_j=R_j/R_{j-1}.",
            "Expand each order-3 Toeplitz determinant exactly as an integer polynomial in q_j.",
            "Remove the common positive monomial and integer gcd; sign is unchanged.",
            "Search simple certificates: nonnegative coefficients, factors 1-q_i with nonnegative quotient, adjacent equality vanishing, and exact D3 polynomial matches.",
        ],
        "summary_counts": summary_counts,
        "top_h812_boundary_analysis": top_analysis,
        "unclassified_examples": family_examples,
        "alternate_idea_register": alternate_idea_register,
        "decision": decision,
        "next_target": {
            "name": "H814 Plucker-Dodgson sparse-row reduction or analytic pivot",
            "statement": (
                "Either express the mixed unclassified sparse-row polynomials using H811-covered "
                "consecutive-row minors through a determinant identity, or pivot to the Xi-kernel "
                "unit-deficit lemma instead of more grid search."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# H813 Sparse-Row Symbolic Hierarchy",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Method",
        "",
    ]
    for item in report["method"]:
        lines.append(f"- {item}")
    lines.extend(
        [
            "",
            "## Summary",
            "",
            "```json",
            json.dumps(report["summary_counts"], indent=2),
            "```",
            "",
            "## H812 Boundary Families",
            "",
        ]
    )
    for item in report["top_h812_boundary_analysis"][:12]:
        cert = item["certificate"]
        lines.extend(
            [
                f"### rows {item['rows']} cols {item['cols']}",
                "",
                f"Family: `{item['family']}`",
                "",
                f"Polynomial: `{item['polynomial']}`",
                "",
                "Certificate:",
                "",
                "```json",
                json.dumps(cert, indent=2),
                "```",
                "",
            ]
        )
    lines.extend(["## Unclassified Examples", ""])
    for item in report["unclassified_examples"]:
        lines.append(
            f"- rows {item['rows']} cols {item['cols']}: `{item['polynomial']}`; certificate `{item['certificate']['type']}`"
        )
    lines.extend(["", "## Alternate Ideas Kept Alive", ""])
    for item in report["alternate_idea_register"]:
        lines.append(f"- `{item['route']}` ({item['source']}): {item['why_not_exhausted']}")
    lines.extend(
        [
            "",
            "## Decision",
            "",
            report["decision"],
            "",
            "## Next Target",
            "",
            f"`{report['next_target']['name']}`: {report['next_target']['statement']}",
            "",
        ]
    )
    out_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="H813 sparse-row symbolic hierarchy audit.")
    parser.add_argument("--h812", default=DEFAULT_H812)
    parser.add_argument("--length", type=int, default=7)
    parser.add_argument("--example-limit", type=int, default=20)
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
                "summary_counts": report["summary_counts"],
                "next_target": report["next_target"]["name"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()




