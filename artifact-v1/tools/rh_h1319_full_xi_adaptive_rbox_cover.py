"""H1319 adaptive Arb cover for the full-Xi third cumulant.

For the gamma-log saddle parameter ``r`` put

    q = pi*r*exp(r) - (9/4)*r - 1,   t = q/2,
    A = r*(pi*exp(r)*(r+1) - 9/4),
    W = sqrt(A)*(z-log(r)).

After factoring the dominant ``m=1, beta=9`` density, the full Xi law in W
has density proportional to ``exp(-Psi_r(W))*H_r(W)``.  This tool encloses
translated raw moments K_j, j=0..3, and certifies the target

    t(r)^2 * kappa_Xi'''(t(r)) >= -3/4

without finite differences.  If

    C = K3*K0^2 - 3*K1*K2*K0 + 2*K1^3,

the target is exactly equivalent to

    8*t^2*C + (3/4)*A^(3/2)*K0^3 >= 0.              (H1319-P)

The central ``m=1`` integrals use a rigorous second-order midpoint rule.  Psi
is evaluated with a cancellation-free Taylor enclosure of expm1(x)-x.  The
positive m>=2 correction is bounded by the H1316 pointwise majorant, and the
outer |W|>=24 tails use the proved H1271 dominant even-moment tail bounds and
0 < H < 1.  All r and W cells are Arb balls; adaptive r endpoints remain exact
Fractions in the report.

This script writes a certificate only for the configured range.  Failed leaf
boxes are reported as interval-dependency failures, never counterexamples.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

from flint import arb, arb_series, ctx


DEFAULT_OUT = Path("research/riemann/h1319_full_xi_adaptive_rbox_cover.json")
H1271_JSON = Path("research/riemann/h1271_gamma_log_even_moment_tail.json")
H1316_NOTE = Path("research/riemann/h1316_full_xi_transfer_assembly.md")


def qarb(value: Fraction) -> arb:
    return arb(value.numerator) / arb(value.denominator)


def ball(lo: Fraction, hi: Fraction) -> arb:
    left = qarb(lo)
    right = qarb(hi)
    return (left + right) / 2 + arb(0, 1) * (right - left) / 2


def interval_payload(value: arb, digits: int) -> dict[str, str]:
    return {
        "interval": value.str(digits),
        "lower": value.lower().str(digits),
        "upper": value.upper().str(digits),
        "radius": value.rad().str(digits),
    }


def zero_to(value: arb) -> arb:
    upper = value.upper()
    return upper / 2 + arb(0, 1) * upper / 2


def symmetric(value: arb) -> arb:
    return arb(0, 1) * value.upper()


def abs_upper(value: arb) -> arb:
    left = abs(value.lower())
    right = abs(value.upper())
    return left if left > right else right


def stable_rem2_arb(value: arb, terms: int = 18) -> arb:
    """Rigorous enclosure of expm1(value)-value without interval subtraction."""

    size = abs_upper(value)
    if size.upper() >= 1:
        # Split-free Taylor is still preferable here.  The geometric remainder
        # below is valid while size/(N+4)<1, which holds throughout this tool's
        # configured compact window with N=18.
        pass
    polynomial = arb(0)
    power = arb(1)
    for index in range(terms + 1):
        polynomial += power / math.factorial(index + 2)
        power *= value
    first_omitted = size ** (terms + 1) / math.factorial(terms + 3)
    ratio = size / (terms + 4)
    if ratio.upper() >= 1:
        raise ValueError("stable_rem2 Taylor ratio is not below one")
    error = first_omitted / (1 - ratio)
    return value * value * (polynomial + arb(0, 1) * error)


def stable_expm1_series(value: arb_series) -> arb_series:
    """Degree-two expm1 composition with a stable constant coefficient."""

    x0, x1, x2 = value[0], value[1], value[2]
    exp_x0 = x0.exp()
    return arb_series(
        [
            x0.expm1(),
            exp_x0 * x1,
            exp_x0 * (x2 + x1 * x1 / 2),
        ]
    )


def stable_rem2_series(value: arb_series) -> arb_series:
    """Degree-two composition of phi2(x)=exp(x)-1-x."""

    x0, x1, x2 = value[0], value[1], value[2]
    expm1_x0 = x0.expm1()
    return arb_series(
        [
            stable_rem2_arb(x0),
            expm1_x0 * x1,
            expm1_x0 * x2 + x0.exp() * x1 * x1 / 2,
        ]
    )


def stable_expm1(value: arb | arb_series) -> arb | arb_series:
    if isinstance(value, arb_series):
        return stable_expm1_series(value)
    return value.expm1()


def stable_rem2(value: arb | arb_series) -> arb | arb_series:
    if isinstance(value, arb_series):
        return stable_rem2_series(value)
    return stable_rem2_arb(value)


class FullXiBoxModel:
    def __init__(self, r: arb):
        self.pi = arb.pi()
        self.alpha = arb(9) / 4
        self.c5 = arb(3) / (2 * self.pi)
        self.r = r
        self.exp_r = r.exp()
        self.q = self.pi * r * self.exp_r - self.alpha * r - 1
        self.t = self.q / 2
        self.A = r * (self.pi * self.exp_r * (r + 1) - self.alpha)
        self.sqrt_A = self.A.sqrt()

    def psi(self, w: arb | arb_series) -> arb | arb_series:
        h = w / self.sqrt_A
        delta = self.r * stable_expm1(h)
        return (
            self.pi * self.exp_r * stable_rem2(delta)
            + (self.pi * self.r * self.exp_r - self.alpha * self.r)
            * stable_rem2(h)
        )

    def rho(self, w: arb | arb_series) -> arb | arb_series:
        return self.r * (w / self.sqrt_A).exp()

    def density(self, w: arb | arb_series) -> arb | arb_series:
        return (-self.psi(w)).exp()

    def m1_multiplier(self, w: arb | arb_series) -> arb | arb_series:
        # 1 - c5*exp(-rho), the exact m=1 factor in the full Xi kernel.
        return 1 - self.c5 * (-self.rho(w)).exp()

    def main_integrand(
        self, w: arb | arb_series, order: int, center: arb
    ) -> arb | arb_series:
        weight = self.density(w) * self.m1_multiplier(w)
        if order == 0:
            return weight
        return (w - center) ** order * weight


def paired_main_values(
    model: FullXiBoxModel,
    t: arb | arb_series,
    center: arb,
) -> list[arb | arb_series]:
    return [
        model.main_integrand(t, order, center)
        + model.main_integrand(-t, order, center)
        for order in range(4)
    ]


def midpoint_main_moments(
    model: FullXiBoxModel,
    width: Fraction,
    subdivisions: int,
    center: Fraction,
) -> list[arb]:
    """Enclose the four central m=1 moments on [-width,width]."""

    step = width / subdivisions
    step_a = qarb(step)
    remainder_scale = step_a**3 / 12
    center_a = qarb(center)
    totals = [arb(0) for _ in range(4)]
    for index in range(subdivisions):
        lo = index * step
        hi = lo + step
        midpoint = (lo + hi) / 2
        point_values = paired_main_values(model, qarb(midpoint), center_a)
        t_series = arb_series([ball(lo, hi), arb(1)])
        series_values = paired_main_values(model, t_series, center_a)
        for order in range(4):
            second_coefficient = series_values[order][2]
            error = remainder_scale * abs_upper(second_coefficient)
            totals[order] += step_a * point_values[order] + arb(0, 1) * error
    return totals


def validate_dependencies() -> dict[str, Any]:
    h1271 = json.loads(H1271_JSON.read_text(encoding="utf-8"))
    powers = h1271["tail_bounds"]["powers"]
    h1316_text = H1316_NOTE.read_text(encoding="utf-8")
    checks = {
        "h1271_closed": h1271.get("classification")
        == "h1271_gamma_log_even_moment_tail_closed",
        "h1271_all_checks_pass": bool(h1271.get("all_checks_pass")),
        "h1271_range_contains_cover": h1271["box_certificate"]["r_range"]
        == {"left": "5/2", "right": "100"},
        "h1271_W_is_24": h1271["box_certificate"]["W"] == "24",
        "h1271_even_targets_pass": all(
            bool(powers[str(power)]["target_pass"]) for power in (0, 2, 4)
        ),
        "h1316_full_multiplier_sign_present": "-c_5/X<h(X)<0" in h1316_text,
        "h1316_arithmetic_majorant_present": "19 exp(-3X)" in h1316_text,
    }
    return {
        "checks": checks,
        "all_pass": all(checks.values()),
        "h1271_targets_used": {
            "T0": "2e-40",
            "T2": "2e-37",
            "T4": "7e-35",
            "T1_rule": "T1 <= T2/24 on |W|>=24",
            "T3_rule": "T3 <= T4/24 on |W|>=24",
        },
        "references": [str(H1271_JSON), str(H1316_NOTE)],
    }


def outer_translated_tail_bounds(center: Fraction) -> list[arb]:
    """H1271/H1316 bounds for full-Xi translated moments outside W=24."""

    raw = [
        arb("2e-40"),
        arb("2e-37") / 24,
        arb("2e-37"),
        arb("7e-35") / 24,
    ]
    absolute_center = qarb(abs(center))
    translated: list[arb] = []
    for order in range(4):
        value = arb(0)
        for power in range(order + 1):
            value += (
                math.comb(order, power)
                * absolute_center ** (order - power)
                * raw[power]
            )
        translated.append(value)
    return translated


def central_arithmetic_tail_bounds(
    model: FullXiBoxModel,
    width: Fraction,
    center: Fraction,
) -> tuple[list[arb], arb]:
    """Bound the positive m>=2 correction on the central window.

    H1316 gives h_bar(rho)<19*exp(-3*pi*exp(rho)).  Rho is increasing
    in W, Psi>=0, and hence exp(-Psi)<=1.  A deliberately crude length-times-
    supremum estimate is already tiny at this saddle range.
    """

    width_a = qarb(width)
    rho_at_left = model.r * (-width_a / model.sqrt_A).exp()
    h_bar = 19 * (-3 * arb.pi() * rho_at_left.exp()).exp()
    h_bar_upper = h_bar.upper()
    radius = width_a + qarb(abs(center))
    bounds = [2 * width_a * radius**order * h_bar_upper for order in range(4)]
    return bounds, h_bar_upper


def endpoint_float(value: arb) -> float:
    try:
        return float(value)
    except (TypeError, ValueError, OverflowError):
        return float("nan")


def certify_box(
    r_lo: Fraction,
    r_hi: Fraction,
    width: Fraction,
    subdivisions: int,
    center: Fraction,
    beta: Fraction,
    digits: int,
) -> dict[str, Any]:
    model = FullXiBoxModel(ball(r_lo, r_hi))
    main = midpoint_main_moments(model, width, subdivisions, center)
    arithmetic_tail, h_bar = central_arithmetic_tail_bounds(model, width, center)
    outer_tail = outer_translated_tail_bounds(center)
    total_error = [arithmetic_tail[j] + outer_tail[j] for j in range(4)]
    moments = [
        main[0] + zero_to(total_error[0]),
        main[1] + symmetric(total_error[1]),
        main[2] + zero_to(total_error[2]),
        main[3] + symmetric(total_error[3]),
    ]
    k0, k1, k2, k3 = moments
    cumulant_numerator = (
        k3 * k0 * k0 - 3 * k1 * k2 * k0 + 2 * k1 * k1 * k1
    )
    a_three_halves = model.A * model.sqrt_A
    polynomial = (
        8 * model.t * model.t * cumulant_numerator
        + qarb(beta) * a_three_halves * k0 * k0 * k0
    )
    denominator = a_three_halves * k0 * k0 * k0
    scaled_cumulant = 8 * model.t * model.t * cumulant_numerator / denominator

    mass_positive = bool(k0.lower() > 0)
    polynomial_positive = bool(polynomial.lower() > 0)
    passed = mass_positive and polynomial_positive
    if not mass_positive:
        reason = "K0_lower_nonpositive"
    elif not polynomial_positive:
        reason = "H1319_polynomial_lower_nonpositive_interval_dependency"
    else:
        reason = "certified"

    return {
        "r_lo": str(r_lo),
        "r_hi": str(r_hi),
        "r_width": str(r_hi - r_lo),
        "status": "certified" if passed else "failed",
        "reason": reason,
        "K_main": [interval_payload(value, digits) for value in main],
        "central_arithmetic_tail_abs": [
            interval_payload(value, digits) for value in arithmetic_tail
        ],
        "outer_tail_abs": [
            interval_payload(value, digits) for value in outer_tail
        ],
        "K_total": [interval_payload(value, digits) for value in moments],
        "h1316_central_hbar_sup": interval_payload(h_bar, digits),
        "cumulant_numerator_C": interval_payload(cumulant_numerator, digits),
        "target_polynomial_P": interval_payload(polynomial, digits),
        "scaled_t2_kappa3_diagnostic": interval_payload(scaled_cumulant, digits),
        "P_lower_float": endpoint_float(polynomial.lower()),
        "scaled_lower_float": endpoint_float(scaled_cumulant.lower()),
    }


def exact_partition(lo: Fraction, hi: Fraction, count: int) -> list[tuple[Fraction, Fraction]]:
    step = (hi - lo) / count
    return [(lo + index * step, lo + (index + 1) * step) for index in range(count)]


def adaptive_cover(
    r_lo: Fraction,
    r_hi: Fraction,
    initial_boxes: int,
    max_depth: int,
    width: Fraction,
    subdivisions: int,
    center: Fraction,
    beta: Fraction,
    digits: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], int]:
    certified: list[dict[str, Any]] = []
    failed: list[dict[str, Any]] = []
    attempts = 0
    initial = exact_partition(r_lo, r_hi, initial_boxes)
    stack = [(lo, hi, 0) for lo, hi in reversed(initial)]
    while stack:
        lo, hi, depth = stack.pop()
        result = certify_box(lo, hi, width, subdivisions, center, beta, digits)
        attempts += 1
        result["depth"] = depth
        if result["status"] == "certified":
            certified.append(result)
        elif depth >= max_depth:
            failed.append(result)
        else:
            midpoint = (lo + hi) / 2
            stack.append((midpoint, hi, depth + 1))
            stack.append((lo, midpoint, depth + 1))
    certified.sort(key=lambda row: Fraction(row["r_lo"]))
    failed.sort(key=lambda row: Fraction(row["r_lo"]))
    return certified, failed, attempts


def r125_bracket(digits: int) -> dict[str, Any]:
    lo = Fraction(323, 100)
    hi = Fraction(13, 4)
    lo_model = FullXiBoxModel(qarb(lo))
    hi_model = FullXiBoxModel(qarb(hi))
    bracket_ball = ball(lo, hi)
    derivative = arb.pi() * bracket_ball.exp() * (bracket_ball + 1) - arb(9) / 4
    checks = {
        "q_at_lo_below_250": bool(lo_model.q.upper() < 250),
        "q_at_hi_above_250": bool(hi_model.q.lower() > 250),
        "q_strictly_increasing_on_bracket": bool(derivative.lower() > 0),
    }
    return {
        "definition": "r(125) is the unique r with q(r)=2*125=250",
        "bracket": [str(lo), str(hi)],
        "q_at_lo": interval_payload(lo_model.q, digits),
        "q_at_hi": interval_payload(hi_model.q, digits),
        "q_prime_on_bracket": interval_payload(derivative, digits),
        "checks": checks,
        "all_pass": all(checks.values()),
    }


def best_row(rows: list[dict[str, Any]], field: str) -> dict[str, Any] | None:
    finite = [row for row in rows if math.isfinite(float(row[field]))]
    if not finite:
        return None
    row = min(finite, key=lambda item: float(item[field]))
    return {
        "value": row[field],
        "r_box": [row["r_lo"], row["r_hi"]],
        "depth": row["depth"],
    }


def build_markdown(report: dict[str, Any]) -> str:
    summary = report["summary"]
    config = report["config"]
    lines = [
        "# H1319 Full-Xi Adaptive Arb R-Box Cover",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Result",
        "",
        f"- covered range: `[{config['r_lo']}, {config['r_hi']}]`",
        f"- target: `t(r)^2 kappa_Xi'''(t(r)) >= -{config['beta']}`",
        f"- all leaf boxes certified: `{summary['all_certified']}`",
        f"- certified boxes: `{summary['certified_box_count']}`",
        f"- failed leaf boxes: `{summary['failed_box_count']}`",
        f"- box attempts: `{summary['attempt_count']}`",
        f"- worst certified polynomial lower bound: `{summary['worst_P_lower']}`",
        f"- worst scaled diagnostic lower bound: `{summary['worst_scaled_lower']}`",
        "",
        "## Exact certificate form",
        "",
        "For translated full-Xi moments `K_j`, set",
        "`C=K3*K0^2-3*K1*K2*K0+2*K1^3`.  Each successful box proves",
        "",
        "`8*t^2*C + beta*A^(3/2)*K0^3 > 0`.",
        "",
        "This is algebraically equivalent to the requested lower bound because",
        "`A>0` and `K0>0`.  No finite differences are used.",
        "",
        "## Rigorous enclosures",
        "",
        "The `m=1` central integrals use a second-order Arb midpoint rule and a",
        "Taylor-enclosed `phi2(x)=exp(x)-1-x`, avoiding correlated subtraction.",
        "The positive `m>=2` central correction uses the H1316 majorant",
        "`19 exp(-3*pi*exp(rho))`.  Outside `|W|=24`, `0<H<1` reduces the",
        "full-Xi tails to the H1271 dominant tail certificate.",
        "",
        "## Scope",
        "",
        "This artifact proves only its configured r-range when `all_certified` is",
        "true.  A failed box means the chosen interval discretization did not",
        "separate the target from zero; it is not a counterexample.",
        "",
    ]
    if report["failed_boxes"]:
        lines.extend(
            [
                "## Uncertified leaves",
                "",
                "Exact intervals and moment radii are retained in the JSON under",
                "`failed_boxes`; the next numerical step is Taylor control in r",
                "or a score/Stein form on precisely those leaves.",
                "",
            ]
        )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r-lo", default="323/100")
    parser.add_argument("--r-hi", default="4")
    parser.add_argument("--initial-boxes", type=int, default=32)
    parser.add_argument("--max-depth", type=int, default=5)
    parser.add_argument("--W", default="24")
    parser.add_argument("--subdivisions", type=int, default=600)
    parser.add_argument("--center", default="-1/12")
    parser.add_argument("--beta", default="3/4")
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--digits", type=int, default=24)
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--markdown-out", default="")
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()

    ctx.dps = args.dps
    ctx.cap = 3
    r_lo = Fraction(args.r_lo)
    r_hi = Fraction(args.r_hi)
    width = Fraction(args.W)
    center = Fraction(args.center)
    beta = Fraction(args.beta)
    if width != 24:
        raise ValueError("this H1271-dependent certificate requires W=24")
    if r_lo < Fraction(5, 2) or r_hi > 100 or r_lo >= r_hi:
        raise ValueError("configured r-range must lie in [5/2,100] with lo<hi")
    if args.initial_boxes <= 0 or args.subdivisions <= 0 or args.max_depth < 0:
        raise ValueError("box counts and subdivisions must be positive")

    dependencies = validate_dependencies()
    bracket = r125_bracket(args.digits)
    if not dependencies["all_pass"]:
        raise RuntimeError("H1271/H1316 dependency validation failed")
    if not bracket["all_pass"]:
        raise RuntimeError("the rational r(125) bracket validation failed")

    certified, failed, attempts = adaptive_cover(
        r_lo,
        r_hi,
        args.initial_boxes,
        args.max_depth,
        width,
        args.subdivisions,
        center,
        beta,
        args.digits,
    )
    all_certified = not failed
    classification = (
        "h1319_full_xi_adaptive_rbox_cover_certified"
        if all_certified
        else "h1319_full_xi_adaptive_rbox_pilot_with_uncertified_leaves"
    )
    report: dict[str, Any] = {
        "classification": classification,
        "statement": (
            "On every certified r-box, t(r)^2*kappa_Xi'''(t(r)) >= -beta; "
            "the full requested range is proved iff summary.all_certified is true."
        ),
        "config": {
            "r_lo": str(r_lo),
            "r_hi": str(r_hi),
            "initial_boxes": args.initial_boxes,
            "max_depth": args.max_depth,
            "W": str(width),
            "subdivisions": args.subdivisions,
            "center": str(center),
            "beta": str(beta),
            "dps": args.dps,
            "digits": args.digits,
        },
        "r125_bracket": bracket,
        "dependencies": dependencies,
        "identity": {
            "q": "pi*r*exp(r)-(9/4)*r-1",
            "t": "q/2",
            "kappa_third": "8*C/(A^(3/2)*K0^3)",
            "C": "K3*K0^2-3*K1*K2*K0+2*K1^3",
            "target_polynomial": "8*t^2*C+beta*A^(3/2)*K0^3",
        },
        "summary": {
            "all_certified": all_certified,
            "certified_box_count": len(certified),
            "failed_box_count": len(failed),
            "attempt_count": attempts,
            "worst_P_lower": best_row(certified, "P_lower_float"),
            "worst_scaled_lower": best_row(certified, "scaled_lower_float"),
        },
        "certified_boxes": certified,
        "failed_boxes": failed,
        "limitations": [
            "finite configured r-range only",
            "finite adaptive depth",
            "H1271 and H1316 are explicit proof dependencies",
            "failed leaves are numerical non-certifications, not counterexamples",
        ],
        "next_action_if_failed": (
            "Apply a Taylor model in r or the full-law score/Stein identity only "
            "to the exact failed leaf boxes recorded here."
        ),
    }

    out_path = Path(args.out)
    markdown_path = (
        Path(args.markdown_out) if args.markdown_out else out_path.with_suffix(".md")
    )
    if not args.no_write:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        markdown_path.write_text(build_markdown(report), encoding="utf-8")

    print(
        json.dumps(
            {
                "classification": classification,
                "summary": report["summary"],
                "out": None if args.no_write else str(out_path),
                "markdown": None if args.no_write else str(markdown_path),
            },
            indent=2,
        )
    )
    return 0 if all_certified else 1


if __name__ == "__main__":
    raise SystemExit(main())
