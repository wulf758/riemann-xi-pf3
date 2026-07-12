import argparse
import json
import sys
from pathlib import Path
from typing import Any


DEFAULT_OUT = "research/riemann/h218_jensen_arb_discriminants.json"


def load_flint(extra_path: str):
    path = Path(extra_path)
    if path.exists():
        sys.path.insert(0, str(path.resolve()))
    import flint  # type: ignore
    from flint import acb, arb, ctx  # type: ignore

    return flint, acb, arb, ctx


def arb_record(value: Any) -> dict[str, str]:
    return {
        "repr": str(value),
        "lower": format(float(value.lower()), ".17g"),
        "upper": format(float(value.upper()), ".17g"),
        "mid": format(float(value.mid()), ".17g"),
        "rad": format(float(value.rad()), ".17g"),
    }


def acb_record(value: Any) -> dict[str, Any]:
    return {
        "repr": str(value),
        "real": arb_record(value.real),
        "imag": arb_record(value.imag),
        "abs_upper": format(float(value.abs_upper()), ".17g"),
    }


def xi_function(acb: Any, arb: Any, s: Any) -> Any:
    one = acb(arb("1"), arb("0"))
    pi = acb(arb.pi(), arb("0"))
    return acb(arb("0.5"), arb("0")) * s * (s - one) * (pi ** (-s / 2)) * (s / 2).gamma() * s.zeta()


def xi_on_critical_line(acb: Any, arb: Any, z: Any) -> Any:
    return xi_function(acb, arb, acb(arb("0.5"), arb("0")) + acb(arb("0"), arb("1")) * z)


def gamma_coefficients(
    acb: Any,
    arb: Any,
    n_max: int,
    radius_text: str,
    samples: int,
) -> list[Any]:
    radius = arb(radius_text)
    max_power = 2 * n_max
    sums = [acb(arb("0"), arb("0")) for _ in range(max_power + 1)]
    two_pi = arb.pi() * 2
    for sample in range(samples):
        theta = two_pi * arb(str(sample)) / arb(str(samples))
        omega = acb(theta.cos(), theta.sin())
        z = acb(radius, arb("0")) * omega
        value = xi_on_critical_line(acb, arb, z)
        phase = acb(arb("1"), arb("0"))
        inverse = 1 / omega
        for k in range(max_power + 1):
            sums[k] += value * phase
            phase *= inverse
    gammas = []
    for n in range(n_max + 1):
        coeff = (sums[2 * n] / samples) / (acb(radius, arb("0")) ** (2 * n))
        gamma = coeff.real * arb(str(n)).fac()
        if n % 2:
            gamma = -gamma
        gammas.append(gamma)
    return gammas


def jensen_coefficients(arb: Any, gammas: list[Any], degree: int, shift: int) -> list[Any]:
    coeffs = []
    for j in range(degree + 1):
        coeffs.append(arb(str(binomial(degree, j))) * gammas[shift + j])
    return coeffs


def binomial(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    k = min(k, n - k)
    result = 1
    for i in range(1, k + 1):
        result = result * (n - k + i) // i
    return result


def discriminant(arb: Any, coeffs: list[Any]) -> Any:
    degree = len(coeffs) - 1
    if degree == 2:
        c, b, a = coeffs
        return b * b - 4 * a * c
    if degree == 3:
        d, c, b, a = coeffs
        return 18 * a * b * c * d - 4 * b**3 * d + b * b * c * c - 4 * a * c**3 - 27 * a * a * d * d
    raise ValueError(f"unsupported degree {degree}")


def discriminant_rows(
    arb: Any,
    gammas: list[Any],
    max_shift: int,
    degrees: list[int],
) -> list[dict[str, Any]]:
    rows = []
    for degree in degrees:
        for shift in range(max_shift + 1):
            coeffs = jensen_coefficients(arb, gammas, degree, shift)
            disc = discriminant(arb, coeffs)
            rows.append(
                {
                    "degree": degree,
                    "shift": shift,
                    "discriminant": arb_record(disc),
                    "positive_lower_bound": bool(float(disc.lower()) > 0.0),
                    "coefficients": [arb_record(coeff) for coeff in coeffs],
                }
            )
    return rows


def canary_rows(
    arb: Any,
    gammas: list[Any],
    degree: int,
    shift: int,
    index: int,
    factors: list[str],
) -> list[dict[str, Any]]:
    rows = []
    for factor_text in factors:
        mutated = list(gammas)
        mutated[index] = mutated[index] * arb(factor_text)
        coeffs = jensen_coefficients(arb, mutated, degree, shift)
        disc = discriminant(arb, coeffs)
        rows.append(
            {
                "degree": degree,
                "shift": shift,
                "mutated_gamma_index": index,
                "factor": factor_text,
                "discriminant": arb_record(disc),
                "positive_lower_bound": bool(float(disc.lower()) > 0.0),
            }
        )
    return rows


def classify(report: dict[str, Any]) -> tuple[str, str]:
    checks = report["checks"]
    if not checks["all_gamma_positive_lower_bounds"]:
        return (
            "jensen_arb_gamma_positive_not_certified",
            "At least one gamma coefficient ball did not have a strictly positive lower bound.",
        )
    if not checks["all_discriminants_positive_lower_bounds"]:
        return (
            "jensen_arb_discriminant_certificate_failure",
            "At least one degree 2/3 Jensen discriminant lacks a strictly positive lower bound.",
        )
    if not checks["canary_breaks_discriminant_certificate"]:
        return (
            "jensen_arb_discriminant_canary_insensitive",
            "The discriminant certificate passed but canaries did not break it.",
        )
    return (
        "jensen_arb_low_degree_discriminants_certified_finite",
        (
            "Arb finite DFT coefficients certify positive gamma lower bounds and "
            "positive discriminants for tested degree 2/3 Jensen polynomials. "
            "This certifies a small finite DFT window, modulo alias/truncation; it is not RH."
        ),
    )


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    flint, acb, arb, ctx = load_flint(args.mpmath_path)
    ctx.dps = args.dps
    degrees = [int(item.strip()) for item in args.degrees.split(",") if item.strip()]
    factors = [item.strip() for item in args.canary_factors.split(",") if item.strip()]
    n_max = args.max_shift + max(degrees) + 2
    gammas = gamma_coefficients(acb, arb, n_max, args.radius, args.samples)
    rows = discriminant_rows(arb, gammas, args.max_shift, degrees)
    canaries = canary_rows(
        arb,
        gammas,
        args.canary_degree,
        args.canary_shift,
        args.canary_index,
        factors,
    )
    gamma_selected = [
        {
            "n": n,
            "gamma": arb_record(gammas[n]),
            "positive_lower_bound": bool(float(gammas[n].lower()) > 0.0),
        }
        for n in range(min(n_max + 1, args.max_shift + max(degrees) + 1))
    ]
    min_disc = min(rows, key=lambda row: float(row["discriminant"]["lower"]))
    checks = {
        "all_gamma_positive_lower_bounds": all(item["positive_lower_bound"] for item in gamma_selected),
        "all_discriminants_positive_lower_bounds": all(row["positive_lower_bound"] for row in rows),
        "canary_breaks_discriminant_certificate": any(not row["positive_lower_bound"] for row in canaries),
    }
    report = {
        "schema": "rh_h218_jensen_arb_discriminants.v0",
        "parameters": {
            "max_shift": args.max_shift,
            "degrees": degrees,
            "samples": args.samples,
            "radius": args.radius,
            "dps": args.dps,
            "canary_degree": args.canary_degree,
            "canary_shift": args.canary_shift,
            "canary_index": args.canary_index,
            "canary_factors": factors,
            "mpmath_path": args.mpmath_path,
            "python_flint_version": getattr(flint, "__version__", "unknown"),
        },
        "definition": {
            "gamma_n": "gamma(n)=(-1)^n n! [z^(2n)] Xi(z)",
            "jensen_polynomial": "J^{d,n}(X)=sum binomial(d,j) gamma(n+j) X^j",
            "degree_2_certificate": "positive discriminant b^2-4ac",
            "degree_3_certificate": "positive cubic discriminant",
        },
        "gamma_coefficients": gamma_selected,
        "discriminants": {
            "tested_count": len(rows),
            "min_lower_bound_row": min_disc,
            "rows": rows,
        },
        "canaries": canaries,
        "limitations": [
            "The finite DFT sums are ball-enclosed, but Cauchy alias/truncation is not bounded here.",
            "The certificate is only for degree 2 and 3 Jensen polynomials in a finite shift window.",
            "Positive discriminant is used only where it is sufficient: degrees 2 and 3 with real coefficients.",
            "This does not prove RH or all Jensen hyperbolicity.",
        ],
        "checks": checks,
    }
    classification, interpretation = classify(report)
    report["classification"] = classification
    report["interpretation"] = interpretation
    report["next_target"] = {
        "name": "H219 Jensen DFT alias bounds or degree-4 Sturm",
        "statement": (
            "Either add a Cauchy alias bound to the Arb coefficient extraction, or extend "
            "the finite certificate to degree 4 via a Sturm/subresultant criterion."
        ),
    }
    return report


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    lines = [
        "# RH H218 Jensen Arb Discriminants",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Parameters",
        "",
        "| parameter | value |",
        "| --- | ---: |",
    ]
    for key, value in report["parameters"].items():
        lines.append(f"| {key} | `{value}` |")
    lines.extend(
        [
            "",
            "## Gamma Coefficients",
            "",
            "| n | gamma lower | gamma upper | positive lower |",
            "| ---: | ---: | ---: | --- |",
        ]
    )
    for item in report["gamma_coefficients"]:
        gamma = item["gamma"]
        lines.append(
            f"| {item['n']} | `{gamma['lower']}` | `{gamma['upper']}` | `{item['positive_lower_bound']}` |"
        )
    min_row = report["discriminants"]["min_lower_bound_row"]
    lines.extend(
        [
            "",
            "## Discriminants",
            "",
            f"Tested discriminants: `{report['discriminants']['tested_count']}`",
            f"Smallest lower bound: degree `{min_row['degree']}`, shift `{min_row['shift']}`, "
            f"lower `{min_row['discriminant']['lower']}`.",
            "",
            "| degree | shift | discriminant lower | discriminant upper | positive lower |",
            "| ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in report["discriminants"]["rows"]:
        disc = row["discriminant"]
        if row["shift"] in {0, 1, report["parameters"]["max_shift"]}:
            lines.append(
                f"| {row['degree']} | {row['shift']} | `{disc['lower']}` | `{disc['upper']}` | `{row['positive_lower_bound']}` |"
            )
    lines.extend(
        [
            "",
            "## Canaries",
            "",
            "| degree | shift | gamma index | factor | discriminant lower | positive lower |",
            "| ---: | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for row in report["canaries"]:
        lines.append(
            f"| {row['degree']} | {row['shift']} | {row['mutated_gamma_index']} | `{row['factor']}` | "
            f"`{row['discriminant']['lower']}` | `{row['positive_lower_bound']}` |"
        )
    lines.extend(
        [
            "",
            "## Checks",
            "",
            "| check | value |",
            "| --- | --- |",
        ]
    )
    for key, value in report["checks"].items():
        lines.append(f"| {key} | `{value}` |")
    lines.extend(["", "## Interpretation", "", report["interpretation"], "", "## Limitations", ""])
    for limitation in report["limitations"]:
        lines.append(f"- {limitation}")
    lines.extend(
        [
            "",
            "## Next Target",
            "",
            f"`{report['next_target']['name']}`: {report['next_target']['statement']}",
            "",
        ]
    )
    out_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="H218 Arb discriminant certificates for low-degree Jensen polynomials.")
    parser.add_argument("--max-shift", type=int, default=12)
    parser.add_argument("--degrees", default="2,3")
    parser.add_argument("--samples", type=int, default=1024)
    parser.add_argument("--radius", default="4")
    parser.add_argument("--dps", type=int, default=100)
    parser.add_argument("--canary-degree", type=int, default=3)
    parser.add_argument("--canary-shift", type=int, default=0)
    parser.add_argument("--canary-index", type=int, default=2)
    parser.add_argument("--canary-factors", default="0.01,0.1,10,100")
    parser.add_argument("--mpmath-path", default=".tmp/pydeps")
    parser.add_argument("--out", default=DEFAULT_OUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    out_path = Path(args.out)
    report = build_report(args)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    write_markdown(report, out_path.with_suffix(".md"))
    print(json.dumps(report["checks"], indent=2))
    print(report["classification"])


if __name__ == "__main__":
    main()
