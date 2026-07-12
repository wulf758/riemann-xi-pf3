import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any


THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from rh_h218_jensen_arb_discriminants import (  # noqa: E402
    gamma_coefficients,
    load_flint,
    xi_function,
)


DEFAULT_OUT = "research/riemann/h1321_h801_alias_repaired_finite_h803.json"
DEFAULT_FLINT_PATH = "C:/tmp/h784_pydeps_copy"


def arb_record(value: Any) -> dict[str, Any]:
    return {
        "repr": str(value),
        "lower": str(value.lower()),
        "upper": str(value.upper()),
        "mid": str(value.mid()),
        "rad": str(value.rad()),
        "positive_lower_bound": bool(value.lower() > 0),
        "contains_zero": bool(value.lower() <= 0 <= value.upper()),
    }


def acb_record(value: Any) -> dict[str, Any]:
    return {
        "repr": str(value),
        "real": arb_record(value.real),
        "imag": arb_record(value.imag),
        "abs_upper": str(value.abs_upper()),
    }


def c(index: int) -> int:
    return (2 * index) * (2 * index - 1)


def ratios(values: list[Any]) -> list[Any | None]:
    result: list[Any | None] = [None]
    result.extend(values[index] / values[index - 1] for index in range(1, len(values)))
    return result


def second_ratios(values: list[Any | None]) -> list[Any | None]:
    result: list[Any | None] = [None, None]
    result.extend(values[index] / values[index - 1] for index in range(2, len(values)))
    return result


def minimum_row(rows: list[dict[str, Any]], key: str) -> dict[str, Any]:
    return min(rows, key=lambda row: float(row[key].lower()))


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    if args.n_max < 127:
        raise ValueError("--n-max must be at least 127 to certify H803 through n=126")
    if args.samples <= 2 * args.n_max:
        raise ValueError("--samples must exceed 2*n_max")

    flint, acb, arb, ctx = load_flint(args.flint_path)
    ctx.dps = args.dps
    radius = arb(args.radius)
    outer_radius = arb(args.outer_radius)
    if not (radius < outer_radius):
        raise ValueError("--outer-radius must be strictly larger than --radius")

    # H801's finite root-of-unity sum, now used only as the centre of a
    # coefficient enclosure.  gamma_coefficients encloses every sampled Xi
    # value and every arithmetic operation with python-flint balls.
    raw_gammas = gamma_coefficients(
        acb,
        arb,
        args.n_max,
        args.radius,
        args.samples,
    )
    raw_a = [
        raw_gammas[index] / arb(str(math.factorial(index)))
        for index in range(args.n_max + 1)
    ]

    # Xi(z)=xi(1/2+i*z).  The positive cosine-kernel representation gives
    #   sup_{|z|=R}|Xi(z)| <= Xi(iR) = xi(1/2+R),
    # where the last equality is xi(s)=xi(1-s).  Evaluation at the positive
    # real point 1/2+R is directly ball-enclosed by Arb.
    outer_s = acb(arb("0.5") + outer_radius, arb("0"))
    outer_xi = xi_function(acb, arb, outer_s)
    outer_abs_upper = outer_xi.abs_upper()

    # For the N-point DFT at radius r, the coefficient of z^(2n) is aliased
    # only by exponents 2n+l*N.  Cauchy's estimate on |z|=R therefore gives
    #   |a_n_tilde-a_n|
    #     <= M_R R^(-2n) (r/R)^N / (1-(r/R)^N).
    rho = (radius / outer_radius) ** args.samples
    alias_geometric = rho / (1 - rho)
    alias_bounds = [
        outer_abs_upper * (outer_radius ** (-2 * index)) * alias_geometric
        for index in range(args.n_max + 1)
    ]
    repaired_a = [
        raw_a[index] + arb(0, alias_bounds[index].upper())
        for index in range(args.n_max + 1)
    ]

    r_values = ratios(repaired_a)
    q_values = second_ratios(r_values)

    moments = [
        repaired_a[index] * arb(str(math.factorial(2 * index)))
        for index in range(args.n_max + 1)
    ]
    s_values = ratios(moments)
    tau_values = second_ratios(s_values)

    rows: list[dict[str, Any]] = []
    for n in range(2, 127):
        q_margin = q_values[n + 1] - q_values[n]
        tau_ratio = tau_values[n + 1] / tau_values[n]
        theta = arb(str(c(n - 1) * c(n + 1))) / arb(str(c(n) ** 2))
        h803_margin = tau_ratio - theta
        rows.append(
            {
                "n": n,
                "q_n": arb_record(q_values[n]),
                "q_next": arb_record(q_values[n + 1]),
                "q_next_minus_q": arb_record(q_margin),
                "tau_next_over_tau": arb_record(tau_ratio),
                "factorial_threshold": arb_record(theta),
                "h803_margin": arb_record(h803_margin),
                "q_margin_positive": bool(q_margin.lower() > 0),
                "h803_margin_positive": bool(h803_margin.lower() > 0),
            }
        )

    alias_rows = [
        {
            "n": index,
            "raw_a": arb_record(raw_a[index]),
            "alias_bound": arb_record(alias_bounds[index]),
            "repaired_a": arb_record(repaired_a[index]),
        }
        for index in sorted({0, 1, 2, 3, 8, 32, 64, 100, 126, 127, args.n_max})
        if index <= args.n_max
    ]

    worst_q = minimum_row(rows, "q_next_minus_q")
    worst_h803 = minimum_row(rows, "h803_margin")
    symbolic_identity_checks = []
    for n in range(2, 127):
        # Cross multiplication of
        # q_(n+1)/q_n = (tau_(n+1)/tau_n) C_n^2/(C_(n-1)C_(n+1)).
        symbolic_identity_checks.append(
            c(n - 1) > 0 and c(n) > 0 and c(n + 1) > 0
        )

    checks = {
        "python_flint_loaded": True,
        "reference_radius_is_11": args.radius == "11",
        "outer_radius_is_12": args.outer_radius == "12",
        "sample_count_is_8192": args.samples == 8192,
        "n_max_reaches_127": args.n_max >= 127,
        "outer_radius_strictly_larger": bool(radius < outer_radius),
        "outer_xi_imaginary_part_contains_zero": bool(
            outer_xi.imag.lower() <= 0 <= outer_xi.imag.upper()
        ),
        "outer_xi_real_part_positive": bool(outer_xi.real.lower() > 0),
        "outer_majorant_positive": bool(outer_abs_upper > 0),
        "alias_ratio_strictly_between_zero_and_one": bool(0 < rho < 1),
        "all_alias_bounds_positive": all(bound.lower() > 0 for bound in alias_bounds),
        "all_repaired_coefficients_positive": all(value.lower() > 0 for value in repaired_a),
        "row_indices_are_exactly_2_through_126": [row["n"] for row in rows]
        == list(range(2, 127)),
        "all_q_margins_positive": all(row["q_margin_positive"] for row in rows),
        "all_h803_margins_positive": all(row["h803_margin_positive"] for row in rows),
        "h803_q_threshold_identity_has_positive_denominators": all(symbolic_identity_checks),
    }
    all_checks_pass = all(checks.values())

    return {
        "schema": "rh_h1321_h801_alias_repaired_finite_h803.v1",
        "classification": (
            "h1321_alias_repaired_finite_h803_certified"
            if all_checks_pass
            else "h1321_alias_repaired_finite_h803_failed"
        ),
        "parameters": {
            "n_max": args.n_max,
            "samples": args.samples,
            "reference_radius": args.radius,
            "outer_radius": args.outer_radius,
            "dps": args.dps,
            "flint_path": args.flint_path,
            "python_flint_version": getattr(flint, "__version__", "unknown"),
        },
        "definitions": {
            "Xi": "Xi(z)=xi(1/2+i*z)=sum_{n>=0}(-1)^n*a_n*z^(2n)",
            "a_n": "a_n=gamma_n/n!",
            "M_n": "M_n=(2n)!*a_n",
            "S_n": "S_n=M_n/M_(n-1)",
            "tau_n": "tau_n=S_n/S_(n-1)",
            "C_n": "C_n=(2n)(2n-1)",
            "q_n": "q_n=(a_n/a_(n-1))/(a_(n-1)/a_(n-2))",
            "h803": "tau_(n+1)/tau_n >= C_(n-1)C_(n+1)/C_n^2",
        },
        "proof_contract": {
            "kernel_majorant": (
                "The positive Xi cosine kernel gives sup_|z|=R |Xi(z)| "
                "<= Xi(iR)=xi(1/2+R)."
            ),
            "dft_alias_identity": (
                "The N-point DFT coefficient at exponent 2n equals the true "
                "coefficient plus sum_{l>=1} c_(2n+lN) r^(lN)."
            ),
            "cauchy_alias_bound": (
                "|a_tilde_n-a_n| <= M_R*R^(-2n)*(r/R)^N/(1-(r/R)^N)."
            ),
            "interval_policy": (
                "Every raw DFT ball is enlarged symmetrically by the analytic "
                "alias bound before any quotient is formed."
            ),
        },
        "outer_majorant": {
            "s": f"1/2+{args.outer_radius}",
            "xi_ball": acb_record(outer_xi),
            "M_R_upper": str(outer_abs_upper),
            "rho": arb_record(rho),
            "geometric_alias_factor": arb_record(alias_geometric),
        },
        "selected_alias_rows": alias_rows,
        "finite_h803": {
            "n_range": [2, 126],
            "row_count": len(rows),
            "worst_q_margin": {
                "n": worst_q["n"],
                "interval": worst_q["q_next_minus_q"],
            },
            "worst_h803_margin": {
                "n": worst_h803["n"],
                "interval": worst_h803["h803_margin"],
            },
            "rows": rows,
        },
        "checks": checks,
        "all_checks_pass": all_checks_pass,
        "scope": {
            "proved_if_pass": (
                "The exact H803 factorial threshold for the true Xi coefficients "
                "for every integer 2<=n<=126."
            ),
            "does_not_prove": [
                "the H803 tail n>=127 (supplied separately by H1319/H920)",
                "PF3",
                "PF-infinity",
                "Jensen hyperbolicity in all degrees",
                "the Riemann Hypothesis",
            ],
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Repair H801's finite DFT alias gap and recertify H803 for 2<=n<=126."
    )
    parser.add_argument("--out", default=DEFAULT_OUT)
    parser.add_argument("--flint-path", default=DEFAULT_FLINT_PATH)
    parser.add_argument("--dps", type=int, default=280)
    parser.add_argument("--n-max", type=int, default=128)
    parser.add_argument("--samples", type=int, default=8192)
    parser.add_argument("--radius", default="11")
    parser.add_argument("--outer-radius", default="12")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = build_report(args)
    except Exception as exc:
        report = {
            "schema": "rh_h1321_h801_alias_repaired_finite_h803.v1",
            "classification": "h1321_alias_repaired_finite_h803_replay_unavailable",
            "all_checks_pass": False,
            "error": f"{type(exc).__name__}: {exc}",
            "invalidation": (
                "The former H801/H804 finite-base status remains numerical evidence, "
                "not a rigorous certificate, until this replay succeeds."
            ),
        }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "classification": report["classification"],
        "all_checks_pass": report["all_checks_pass"],
        "out": str(out_path),
        "error": report.get("error"),
    }, indent=2))
    return 0 if report["all_checks_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
