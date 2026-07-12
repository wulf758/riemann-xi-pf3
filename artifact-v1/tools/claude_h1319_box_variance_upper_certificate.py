"""Extract a rigorous Xi variance upper bound from the H1319 r-box data.

No new quadrature is performed.  The tool consumes the already certified
H1319 moment enclosures on

    323/100 <= r <= 4      (192 quotient-free boxes),
    4 <= r <= 22          (360 score mean-value boxes),

and observes that a variance is bounded by a second moment about any fixed
center.  On the first segment the stored K2/K0 is E[(W+1/12)^2].  On the
middle segment the stored M2/Z is E[W^2].  Thus in both cases

    Var(W) <= second_moment_envelope.

Since q=2t and W=sqrt(A)(Z-r), the exact cumulant scaling is

    kappa_Xi''(t) = 4 Var(W)/A.

The certificate tests the simple boxwise envelope Var(W)<26/25.  It also
checks A>2tr, yielding kappa_Xi''(t)<52/(25tr) throughout the covered range.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

from flint import arb, ctx


DEFAULT_FIRST = "research/riemann/h1319_full_xi_adaptive_rbox_cover.json"
DEFAULT_MIDDLE = (
    "research/riemann/h1319_full_xi_score_mean_value_r_cover_4_5.json",
    "research/riemann/h1319_full_xi_score_mean_value_r_cover_5_10.json",
    "research/riemann/h1319_full_xi_score_mean_value_r10_to_22_arb.json",
)
DEFAULT_OUT = "research/riemann/claude_h1319_box_variance_upper_certificate.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def qarb(value: Fraction) -> arb:
    return arb(value.numerator) / arb(value.denominator)


def rball(left: Fraction, right: Fraction) -> arb:
    lo = qarb(left)
    hi = qarb(right)
    return (lo + hi) / 2 + arb(0, 1) * (hi - lo) / 2


def payload_ball(payload: dict[str, str]) -> arb:
    return arb(payload["interval"])


def symmetric(radius: arb) -> arb:
    return arb(0, 1) * radius.upper()


def ball_record(value: arb) -> dict[str, Any]:
    return {
        "repr": str(value),
        "lower": str(value.lower()),
        "upper": str(value.upper()),
        "mid": str(value.mid()),
        "rad": str(value.rad()),
        "positive_lower_bound": bool(value.lower() > 0),
        "contains_zero": bool(value.lower() <= 0 <= value.upper()),
    }


def saddle_values(left: Fraction, right: Fraction) -> tuple[arb, arb, arb, arb]:
    r = rball(left, right)
    alpha = arb(9) / 4
    q = arb.pi() * r * r.exp() - alpha * r - 1
    t = q / 2
    big_a = r * (arb.pi() * r.exp() * (r + 1) - alpha)
    return r, t, big_a, big_a - 2 * t * r


def row_from_moments(
    source: str,
    kind: str,
    left: Fraction,
    right: Fraction,
    mass: arb,
    first_or_score: arb,
    second: arb,
    center_description: str,
) -> dict[str, Any]:
    if mass.lower() <= 0:
        raise ValueError(f"nonpositive mass lower bound on {left}..{right}")
    raw_second = second / mass
    normalized_first = first_or_score / mass
    centered_variance = raw_second - normalized_first**2
    r, t, big_a, a_gap = saddle_values(left, right)
    kappa2_upper = 4 * raw_second / big_a
    tr_kappa2_upper = t * r * kappa2_upper
    variance_margin = arb(26) / 25 - raw_second
    tr_margin = arb(52) / 25 - tr_kappa2_upper
    return {
        "source": source,
        "kind": kind,
        "r_lo": str(left),
        "r_hi": str(right),
        "center_description": center_description,
        "mass": mass,
        "first_or_score": first_or_score,
        "second": second,
        "raw_second": raw_second,
        "centered_variance": centered_variance,
        "r": r,
        "t": t,
        "A": big_a,
        "A_minus_2tr": a_gap,
        "kappa2_upper": kappa2_upper,
        "tr_kappa2_upper": tr_kappa2_upper,
        "variance_margin": variance_margin,
        "tr_margin": tr_margin,
    }


def first_rows(path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    report = json.loads(path.read_text(encoding="utf-8"))
    rows = []
    for stored in report.get("certified_boxes", []):
        moments = [payload_ball(item) for item in stored["K_total"]]
        rows.append(
            row_from_moments(
                str(path).replace("\\", "/"),
                "quotient_free_translated_moments",
                Fraction(stored["r_lo"]),
                Fraction(stored["r_hi"]),
                moments[0],
                moments[1],
                moments[2],
                "K_j are moments of W-c with c=-1/12; variance is translation invariant",
            )
        )
    checks = {
        "classification": report.get("classification")
        == "h1319_full_xi_adaptive_rbox_cover_certified",
        "all_certified": report.get("summary", {}).get("all_certified") is True,
        "failed_boxes_empty": report.get("failed_boxes") == [],
        "box_count_is_192": len(rows) == 192,
        "center_is_minus_one_twelfth": report.get("config", {}).get("center")
        == "-1/12",
        "dependencies_pass": report.get("dependencies", {}).get("all_pass") is True,
    }
    return rows, {
        "path": str(path).replace("\\", "/"),
        "sha256": sha256(path),
        "checks": checks,
        "passes": all(checks.values()),
    }


def middle_rows(path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    report = json.loads(path.read_text(encoding="utf-8"))
    rows = []
    for stored in report.get("accepted_boxes", []):
        finite = [
            payload_ball(item)
            for item in stored["finite_integral_ranges_mean_value"]
        ]
        errors = [
            payload_ball(item)
            for item in stored["uniform_integral_error_radii"]
        ]
        # Match the H1319 perturbation contract conservatively: every omitted
        # contribution is allowed either sign.  Positivity could tighten the
        # mass and M2 bounds but is not needed for the 26/25 target.
        full_mass = finite[0] + symmetric(errors[0])
        full_second = finite[1] + symmetric(errors[1])
        full_score_a = finite[2] + symmetric(errors[2])
        rows.append(
            row_from_moments(
                str(path).replace("\\", "/"),
                "score_mean_value_moments",
                Fraction(stored["r_lo"]),
                Fraction(stored["r_hi"]),
                full_mass,
                full_score_a,
                full_second,
                "M2/Z=E[W^2]; SA/Z=E[N_Xi]=-E[W]",
            )
        )
    configuration = report.get("configuration", {})
    global_checks = report.get("global_checks", {})
    checks = {
        "classification": report.get("classification")
        == "h1319_full_xi_compact_mean_value_r_interval_certificate",
        "all_checks_pass": report.get("all_checks_pass") is True,
        "unresolved_boxes_empty": report.get("unresolved_boxes") == [],
        "accepted_boxes_nonempty": bool(rows),
        "global_checks_all_true": bool(global_checks) and all(global_checks.values()),
        "mean_value_not_sampling": configuration.get("point_sampling_used") is False,
        "no_finite_differences": configuration.get("finite_differences_used") is False,
    }
    return rows, {
        "path": str(path).replace("\\", "/"),
        "sha256": sha256(path),
        "covered_r_interval": report.get("covered_r_interval"),
        "accepted_box_count": len(rows),
        "checks": checks,
        "passes": all(checks.values()),
    }


def serialized_row(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "source": row["source"],
        "kind": row["kind"],
        "r_lo": row["r_lo"],
        "r_hi": row["r_hi"],
        "center_description": row["center_description"],
        "mass": ball_record(row["mass"]),
        "second_moment_numerator": ball_record(row["second"]),
        "raw_second_moment_envelope_for_variance": ball_record(row["raw_second"]),
        "centered_variance_diagnostic": ball_record(row["centered_variance"]),
        "A": ball_record(row["A"]),
        "A_minus_2tr": ball_record(row["A_minus_2tr"]),
        "kappa_Xi_second_upper_envelope": ball_record(row["kappa2_upper"]),
        "tr_times_kappa_Xi_second_upper_envelope": ball_record(
            row["tr_kappa2_upper"]
        ),
        "margin_raw_second_below_26_over_25": ball_record(row["variance_margin"]),
        "margin_tr_kappa2_below_52_over_25": ball_record(row["tr_margin"]),
    }


def exact_chain(rows: list[dict[str, Any]]) -> bool:
    return bool(rows) and all(
        Fraction(left["r_hi"]) == Fraction(right["r_lo"])
        for left, right in zip(rows, rows[1:])
    )


def mutation_canary(worst: dict[str, Any], factor: Fraction) -> dict[str, Any]:
    mutated_second = worst["second"] * qarb(factor)
    mutated_raw = mutated_second / worst["mass"]
    margin = arb(26) / 25 - mutated_raw
    return {
        "target_box": [worst["r_lo"], worst["r_hi"]],
        "target_source": worst["source"],
        "mutated_field": "second_moment_numerator",
        "factor": str(factor),
        "mutated_raw_second": ball_record(mutated_raw),
        "mutated_margin": ball_record(margin),
        "kills_26_over_25_certificate": bool(margin.upper() < 0),
    }


def write_markdown(report: dict[str, Any], out_path: Path) -> None:
    summary = report.get("summary", {})
    lines = [
        "# H1319 Boxwise Xi Variance Upper Certificate",
        "",
        f"Classification: `{report['classification']}`",
        "",
        "## Result",
        "",
        report["decision"],
        "",
        "For every certified box, the already stored second-moment enclosure gives",
        "",
        "```text",
        "Var(W) <= E[(W-c)^2] < 26/25,",
        "kappa_Xi''(t) = 4 Var(W)/A < 104/(25A) < 52/(25tr).",
        "```",
        "",
        f"- covered r interval: `{summary.get('covered_r_interval')}`",
        f"- boxes: `{summary.get('box_count')}`",
        f"- worst raw-second upper: `{summary.get('worst_raw_second', {}).get('upper')}`",
        f"- worst box: `{summary.get('worst_raw_second', {}).get('r_box')}`",
        f"- minimum margin below 26/25: `{summary.get('minimum_variance_margin', {}).get('lower')}`",
        f"- maximum `tr kappa_Xi''` envelope: `{summary.get('maximum_tr_kappa2_envelope', {}).get('upper')}`",
        "",
        "## Why this is rigorous",
        "",
        "The first 192 boxes already store full translated moments `K0,K1,K2`.",
        "The 360 mean-value boxes store finite `Z,M2,SA` ranges plus uniform",
        "omitted-kernel and outer-tail error radii; this audit adds those radii",
        "symmetrically, matching the original perturbation contract. No point",
        "sampling, interpolation, or new numerical quadrature is introduced.",
        "",
        "## Scope",
        "",
        "This proves a boxwise second-cumulant upper envelope on the displayed",
        "compact r-range. It does not by itself prove translated D3, PF-infinity,",
        "all-degree Jensen hyperbolicity, or RH.",
        "",
    ]
    out_path.write_text("\n".join(lines), encoding="utf-8")


def build_report(args: argparse.Namespace) -> dict[str, Any]:
    first, first_source = first_rows(Path(args.first))
    all_rows = list(first)
    sources = [first_source]
    for text in args.middle:
        rows, source = middle_rows(Path(text))
        all_rows.extend(rows)
        sources.append(source)
    all_rows.sort(key=lambda row: Fraction(row["r_lo"]))

    worst_raw = max(all_rows, key=lambda row: row["raw_second"].upper())
    min_margin = min(all_rows, key=lambda row: row["variance_margin"].lower())
    max_tr = max(all_rows, key=lambda row: row["tr_kappa2_upper"].upper())
    canary = mutation_canary(worst_raw, Fraction(args.canary_factor))
    sign_canary = -min_margin["variance_margin"]

    checks = {
        "all_source_contracts_pass": all(source["passes"] for source in sources),
        "box_count_is_552": len(all_rows) == 552,
        "exact_gap_free_chain": exact_chain(all_rows),
        "chain_starts_at_323_over_100": bool(all_rows)
        and Fraction(all_rows[0]["r_lo"]) == Fraction(323, 100),
        "chain_ends_at_22": bool(all_rows)
        and Fraction(all_rows[-1]["r_hi"]) == Fraction(22),
        "all_mass_lowers_positive": all(row["mass"].lower() > 0 for row in all_rows),
        "all_second_moment_lowers_positive": all(
            row["second"].lower() > 0 for row in all_rows
        ),
        "all_A_lowers_positive": all(row["A"].lower() > 0 for row in all_rows),
        "all_A_minus_2tr_lowers_positive": all(
            row["A_minus_2tr"].lower() > 0 for row in all_rows
        ),
        "all_raw_second_uppers_below_26_over_25": all(
            row["raw_second"].upper() < arb(26) / 25 for row in all_rows
        ),
        "all_tr_kappa2_envelopes_below_52_over_25": all(
            row["tr_kappa2_upper"].upper() < arb(52) / 25 for row in all_rows
        ),
        "second_moment_mutation_canary_breaks_certificate": canary[
            "kills_26_over_25_certificate"
        ],
        "sign_canary_breaks_certificate": bool(sign_canary.upper() < 0),
    }
    all_checks_pass = all(checks.values())
    rows_serialized = [serialized_row(row) for row in all_rows]
    return {
        "schema": "claude_h1319_box_variance_upper_certificate.v1",
        "classification": (
            "h1319_boxwise_xi_variance_upper_certified"
            if all_checks_pass
            else "h1319_boxwise_xi_variance_upper_failed"
        ),
        "parameters": {
            "first": args.first,
            "middle": list(args.middle),
            "arb_dps": args.dps,
            "canary_factor": args.canary_factor,
        },
        "identity": {
            "coordinate": "q=2t; W=sqrt(A)(Z-r)",
            "variance_scaling": "kappa_Xi''(t)=4*Var(W)/A",
            "variance_envelope": "Var(W)<=E[(W-c)^2] for every fixed c",
            "first_segment_center": "c=-1/12",
            "middle_segment_center": "c=0",
            "saddle_gap": (
                "A-2tr=pi*r*exp(r)+(9/4)r(r-1)+r>0 for r>1"
            ),
            "uniform_conclusion": (
                "Var(W)<26/25 => kappa_Xi''(t)<104/(25A)<52/(25tr)"
            ),
        },
        "source_certificates": sources,
        "summary": {
            "covered_r_interval": ["323/100", "22"],
            "box_count": len(all_rows),
            "worst_raw_second": {
                "r_box": [worst_raw["r_lo"], worst_raw["r_hi"]],
                "source": worst_raw["source"],
                **ball_record(worst_raw["raw_second"]),
            },
            "minimum_variance_margin": {
                "r_box": [min_margin["r_lo"], min_margin["r_hi"]],
                **ball_record(min_margin["variance_margin"]),
            },
            "maximum_tr_kappa2_envelope": {
                "r_box": [max_tr["r_lo"], max_tr["r_hi"]],
                **ball_record(max_tr["tr_kappa2_upper"]),
            },
        },
        "boxwise_bounds": rows_serialized,
        "canaries": {
            "second_moment_mutation": canary,
            "sign_canary": ball_record(sign_canary),
            "sign_canary_kills_certificate": bool(sign_canary.upper() < 0),
        },
        "checks": checks,
        "all_checks_pass": all_checks_pass,
        "decision": (
            "The existing H1319 boxes rigorously give Var(W)<26/25 and hence "
            "kappa_Xi''(t)<52/(25*t*r(t)) on 323/100<=r<=22. This is a "
            "new second-cumulant envelope extracted without new quadrature."
            if all_checks_pass
            else "At least one source, chain, variance, scale, or canary check failed."
        ),
        "scope": {
            "proved_if_pass": (
                "boxwise Xi second-cumulant upper envelope on 323/100<=r<=22"
            ),
            "does_not_prove": [
                "translated D3 without the complementary third-cumulant control",
                "PF-infinity",
                "all-degree Jensen hyperbolicity",
                "the Riemann Hypothesis",
            ],
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", default=DEFAULT_FIRST)
    parser.add_argument("--middle", nargs=3, default=list(DEFAULT_MIDDLE))
    parser.add_argument("--out", default=DEFAULT_OUT)
    parser.add_argument("--dps", type=int, default=100)
    parser.add_argument("--canary-factor", default="21/20")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    ctx.dps = args.dps
    try:
        report = build_report(args)
    except Exception as exc:
        report = {
            "schema": "claude_h1319_box_variance_upper_certificate.v1",
            "classification": "h1319_boxwise_xi_variance_upper_replay_unavailable",
            "all_checks_pass": False,
            "error": f"{type(exc).__name__}: {exc}",
            "decision": "The H1319 variance extraction did not complete.",
        }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if "summary" in report:
        write_markdown(report, out.with_suffix(".md"))
    print(
        json.dumps(
            {
                "classification": report["classification"],
                "all_checks_pass": report["all_checks_pass"],
                "out": str(out),
                "error": report.get("error"),
            },
            indent=2,
        )
    )
    return 0 if report["all_checks_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
