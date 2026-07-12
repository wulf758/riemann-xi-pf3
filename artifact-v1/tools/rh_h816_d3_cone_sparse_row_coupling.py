import argparse
import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Any

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from rh_h813_sparse_row_symbolic_hierarchy import (
    Poly,
    determinant_q_poly,
    poly_to_string,
)
from rh_h815_four_term_sparse_bridge import (
    certify_four_term_bridge,
    mixed_unclassified_pairs,
    monomial_text,
)


DEFAULT_OUT = "research/riemann/h816_d3_cone_sparse_row_coupling.json"

RPoly = dict[tuple[int, ...], Fraction]


def to_fraction_poly(poly: Poly) -> RPoly:
    return {monomial: Fraction(coeff) for monomial, coeff in poly.items() if coeff}


def clean_poly(poly: RPoly) -> RPoly:
    return {monomial: coeff for monomial, coeff in poly.items() if coeff}


def shift_poly(poly: Poly, shift: tuple[int, ...]) -> RPoly:
    return {
        tuple(exponent + shift[index] for index, exponent in enumerate(monomial)): Fraction(coeff)
        for monomial, coeff in poly.items()
        if coeff
    }


def subtract_scaled(poly: RPoly, generator: RPoly, scale: Fraction) -> RPoly:
    out = dict(poly)
    for monomial, coeff in generator.items():
        out[monomial] = out.get(monomial, Fraction(0)) - scale * coeff
    return clean_poly(out)


def is_nonnegative(poly: RPoly) -> bool:
    return all(coeff >= 0 for coeff in poly.values())


def negative_terms(poly: RPoly) -> list[tuple[tuple[int, ...], Fraction]]:
    return sorted((monomial, coeff) for monomial, coeff in poly.items() if coeff < 0)


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def enumerate_d3_generators(width: int, target_support: set[tuple[int, ...]]) -> list[dict[str, Any]]:
    generators = []
    max_exp = tuple(max(monomial[index] for monomial in target_support) for index in range(width))
    for start in range(1, width):
        d3, _ = determinant_q_poly((0, 1, 2), (start, start + 1, start + 2), width + 1)
        if not d3 or len(d3) <= 1:
            continue
        d3_max = tuple(max(monomial[index] for monomial in d3) for index in range(width))

        def rec(index: int, current: list[int]) -> None:
            if index == width:
                shift = tuple(current)
                shifted = shift_poly(d3, shift)
                support = set(shifted)
                if support <= target_support:
                    generators.append(
                        {
                            "start": start,
                            "shift": shift,
                            "poly": shifted,
                            "d3": d3,
                        }
                    )
                return
            limit = max_exp[index] - d3_max[index]
            for value in range(max(0, limit) + 1):
                current.append(value)
                rec(index + 1, current)
                current.pop()

        rec(0, [])
    return generators


def generator_capacity(residual: RPoly, generator: RPoly) -> Fraction:
    capacity: Fraction | None = None
    for monomial, coeff in generator.items():
        if coeff <= 0:
            continue
        available = residual.get(monomial, Fraction(0))
        bound = available / coeff
        capacity = bound if capacity is None else min(capacity, bound)
    return Fraction(0) if capacity is None else capacity


def greedy_d3_certificate(poly: Poly, max_steps: int) -> dict[str, Any] | None:
    residual = to_fraction_poly(poly)
    if is_nonnegative(residual):
        return {"type": "already_nonnegative_coefficients", "steps": []}

    target_support = set(poly)
    generators = enumerate_d3_generators(len(next(iter(poly))), target_support)
    steps = []
    for _ in range(max_steps):
        negatives = negative_terms(residual)
        if not negatives:
            return {"type": "d3_cone_greedy", "steps": steps}
        best = None
        for monomial, coeff in negatives:
            need = -coeff
            for generator in generators:
                gcoeff = generator["poly"].get(monomial, Fraction(0))
                if gcoeff >= 0:
                    continue
                exact_scale = need / (-gcoeff)
                capacity = generator_capacity(residual, generator["poly"])
                if exact_scale <= capacity:
                    best = (monomial, exact_scale, generator)
                    break
            if best is not None:
                break
        if best is None:
            return None
        monomial, scale, generator = best
        residual = subtract_scaled(residual, generator["poly"], scale)
        steps.append(
            {
                "cancel_negative_monomial": monomial_text(monomial),
                "scale": fraction_text(scale),
                "d3_start": generator["start"],
                "shift": [int(value) for value in generator["shift"]],
                "shift_monomial": monomial_text(generator["shift"]),
            }
        )
    if is_nonnegative(residual):
        return {"type": "d3_cone_greedy", "steps": steps}
    return None


def one_step_d3_certificate(poly: Poly) -> dict[str, Any] | None:
    residual = to_fraction_poly(poly)
    generators = enumerate_d3_generators(len(next(iter(poly))), set(poly))
    for generator in generators:
        for monomial, coeff in generator["poly"].items():
            if coeff >= 0 or residual.get(monomial, Fraction(0)) >= 0:
                continue
            scale = -residual[monomial] / (-coeff)
            if scale <= 0:
                continue
            candidate = subtract_scaled(residual, generator["poly"], scale)
            if is_nonnegative(candidate):
                return {
                    "type": "d3_cone_one_step",
                    "steps": [
                        {
                            "cancel_negative_monomial": monomial_text(monomial),
                            "scale": fraction_text(scale),
                            "d3_start": generator["start"],
                            "shift": [int(value) for value in generator["shift"]],
                            "shift_monomial": monomial_text(generator["shift"]),
                        }
                    ],
                }
    return None


def remaining_after_h815(length: int) -> list[dict[str, Any]]:
    remaining = []
    for item in mixed_unclassified_pairs(length):
        if certify_four_term_bridge(item["poly"]) is not None:
            continue
        remaining.append(item)
    return remaining


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    remaining = remaining_after_h815(args.length)
    covered = []
    uncovered = []
    for item in remaining:
        certificate = one_step_d3_certificate(item["poly"])
        if certificate is None:
            certificate = greedy_d3_certificate(item["poly"], args.max_steps)
        if certificate is None:
            uncovered.append(item)
        else:
            covered.append({**item["analysis"], "h816_certificate": certificate})

    covered_families = Counter(item["family"] for item in covered)
    uncovered_families = Counter(item["analysis"]["family"] for item in uncovered)
    if covered:
        classification = "h816_d3_cone_sparse_row_coupling_partial"
        decision = (
            "The monomial-D3 cone certifies a nonempty subset of the H815 remainder, "
            "but it does not close the sparse-row bottleneck. Remaining forms likely "
            "need a Plucker-Dodgson identity or a stronger D3-aware product lemma."
        )
    else:
        classification = "h816_d3_cone_direct_certificate_no_hit"
        decision = (
            "The direct monomial-D3 cone did not certify any H815 remainder form. "
            "This argues against another local coefficient-greedy refinement; move to "
            "Plucker-Dodgson identities or the Xi-kernel analytic pivot."
        )

    return {
        "schema": "rh_h816_d3_cone_sparse_row_coupling.v0",
        "classification": classification,
        "candidate_certificate": (
            "P(q) = Q(q) + sum c_i M_i(q) D3_{s_i}(q), where Q has nonnegative "
            "coefficients, c_i>0, M_i are monomials, and D3_s are translated "
            "contiguous order-3 Toeplitz minors normalized by positive factors."
        ),
        "d3_shape": (
            "D3_k = 1 - 2 q_k + q_k^2 q_{k+1} + q_{k-1} q_k^2 "
            "- q_{k-1} q_k^2 q_{k+1} for k>=3, with the k=2 endpoint "
            "D3_2=1-2q2+q2^2 q3."
        ),
        "counts": {
            "h815_remaining_input": len(remaining),
            "covered_by_direct_d3_cone": len(covered),
            "uncovered": len(uncovered),
        },
        "covered_family_counts": covered_families.most_common(),
        "uncovered_family_counts": uncovered_families.most_common(),
        "covered_examples": covered[: args.example_limit],
        "uncovered_examples": [item["analysis"] for item in uncovered[: args.example_limit]],
        "decision": decision,
        "next_target": {
            "name": "H817 Plucker-Dodgson sparse-row identity",
            "statement": (
                "Search for determinant identities expressing remaining sparse-row minors "
                "through H811-covered consecutive-row minors, PF2 minors, and D3 inputs. "
                "If direct identities are circular, pivot to the H804 Xi-kernel unit-deficit route."
            ),
        },
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# H816 D3 Cone Sparse-Row Coupling",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Candidate Certificate",
        "",
        report["candidate_certificate"],
        "",
        "D3 shape:",
        "",
        f"`{report['d3_shape']}`",
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
                json.dumps(item["h816_certificate"], indent=2),
                "```",
                "",
            ]
        )
    lines.extend(["## Uncovered Families", ""])
    for family, count in report["uncovered_family_counts"][:20]:
        lines.append(f"- `{family}`: {count}")
    lines.extend(["", "## Decision", "", report["decision"], "", "## Next Target", ""])
    lines.append(f"`{report['next_target']['name']}`: {report['next_target']['statement']}")
    lines.append("")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="H816 D3 cone sparse-row coupling audit.")
    parser.add_argument("--length", type=int, default=7)
    parser.add_argument("--max-steps", type=int, default=4)
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
