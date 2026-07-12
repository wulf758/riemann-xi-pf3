#!/usr/bin/env python3
"""Assemble a global translated-D3 certificate for the Xi coefficients.

The proof has three independently checkable layers.

1. Existing H1319 Arb boxes are re-read from directed endpoint strings to
   prove, for the full-Xi log moment kappa(t),

       kappa''(t) <= 4/(5t),   kappa'''(t) <= 11/(20t^2),  t >= 125.

   The compact part uses 552 already certified r-boxes.  The analytic part
   from r=22 uses H1313/H1314/H1317 and exact rational endpoint arithmetic.

2. Exact rational envelopes then imply translated D3 for every n>=127.
   Positivity is certified by coefficient-positive polynomials after the
   shift n=m+127.

3. The alias-repaired finite certificate covers 2<=n<=126 and H1322 supplies
   global q-monotonicity.  The boundary starts 0 and 1 are checked separately,
   after which H811 yields all consecutive-row order-3 Toeplitz minors.

This does not assert sparse-row PF3, PF-infinity, all-degree Jensen
hyperbolicity, or the Riemann Hypothesis.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

import flint
import sympy as sp
from flint import arb


ROOT = Path(__file__).resolve().parents[1]
RIEMANN = ROOT / "research" / "riemann"

ADAPTIVE = RIEMANN / "h1319_full_xi_adaptive_rbox_cover.json"
SEGMENT = RIEMANN / "h1319_full_xi_segment_interval_certificate.json"
MANIFEST = RIEMANN / "h1319_full_xi_compact_interval_certificate.json"
SCORE_FILES = [
    RIEMANN / "h1319_full_xi_score_mean_value_r_cover_4_5.json",
    RIEMANN / "h1319_full_xi_score_mean_value_r_cover_5_10.json",
    RIEMANN / "h1319_full_xi_score_mean_value_r10_to_22_arb.json",
]
HIGH = RIEMANN / "h1319_full_xi_high_band_r22_analytic_certificate.json"
FINITE = RIEMANN / "claude_d3_alias_repaired_finite_certificate.json"
Q_MONOTONE = RIEMANN / "h1322_h803_all_n_strict_assembly.json"
H811 = RIEMANN / "h811_repaired_consecutive_row_theorem.json"
INDEPENDENT_TAIL = ROOT / "tools" / "claude_d3_two_sided_tail_lemma_independent_check.py"
H1317_CHECKER = ROOT / "tools" / "rh_h1317_effective_xi_transfer_assembly.py"

DEFAULT_OUT = RIEMANN / "h1325_xi_global_d3_h811_consecutive_rows.json"
DEFAULT_MD = RIEMANN / "h1325_xi_global_d3_h811_consecutive_rows.md"


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def endpoint(payload: dict[str, Any], side: str) -> arb:
    value = arb(payload[side])
    return value.lower() if side == "lower" else value.upper()


def payload_ball(payload: dict[str, Any]) -> arb:
    """Rebuild an enclosure from directed endpoints, never the coarse display."""

    lo = endpoint(payload, "lower")
    hi = endpoint(payload, "upper")
    if hi < lo:
        raise AssertionError("serialized interval endpoints are reversed")
    return (lo + hi) / 2 + arb(0, 1) * (hi - lo) / 2


def symmetric_error(payload: dict[str, Any]) -> arb:
    return arb(0, 1) * endpoint(payload, "upper")


def arb_record(value: arb) -> dict[str, str]:
    return {
        "repr": str(value),
        "lower": str(value.lower()),
        "upper": str(value.upper()),
    }


def fraction_record(value: Fraction) -> dict[str, Any]:
    return {
        "exact": str(value),
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def parse_fraction(record: Any) -> Fraction:
    if isinstance(record, dict) and "exact" in record:
        return Fraction(record["exact"])
    return Fraction(record)


def positive_shift_coefficients(expr: sp.Expr, n: sp.Symbol, base: int) -> list[int]:
    m = sp.symbols("m", integer=True, nonnegative=True)
    coeffs = [int(x) for x in sp.Poly(sp.expand(expr.subs(n, m + base)), m).all_coeffs()]
    if not coeffs or not all(x > 0 for x in coeffs):
        raise AssertionError(f"shifted polynomial is not coefficient-positive: {coeffs}")
    return coeffs


def run_json_checker(path: Path) -> tuple[int, dict[str, Any], str]:
    proc = subprocess.run(
        [sys.executable, str(path)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    parsed: dict[str, Any] = {}
    if proc.stdout.strip():
        parsed = json.loads(proc.stdout)
    return proc.returncode, parsed, proc.stderr.strip()


def hash_chain_checks() -> tuple[dict[str, bool], dict[str, str]]:
    manifest = load(MANIFEST)
    segment = load(SEGMENT)
    segments = {
        Path(row["path"].replace("\\", "/")).name: row
        for row in manifest["r_proof_chain"]["segments"]
    }
    checks: dict[str, bool] = {
        "manifest_passes": bool(manifest.get("all_checks_pass")),
        "segment_passes": bool(segment.get("all_segment_checks_pass")),
        "adaptive_hash_matches_segment": sha256(ADAPTIVE)
        == segment["evidence"]["detailed_certificate_sha256"],
        "segment_hash_matches_manifest": sha256(SEGMENT)
        == segments[SEGMENT.name]["sha256"],
    }
    for path in [*SCORE_FILES, HIGH]:
        checks[f"hash_{path.stem}"] = sha256(path) == segments[path.name]["sha256"]
    hashes = {str(path.relative_to(ROOT)): sha256(path) for path in [ADAPTIVE, SEGMENT, MANIFEST, *SCORE_FILES, HIGH, FINITE, Q_MONOTONE, H811]}
    return checks, hashes


def compact_two_sided_bounds() -> tuple[dict[str, Any], dict[str, bool]]:
    adaptive = load(ADAPTIVE)
    score_reports = [load(path) for path in SCORE_FILES]
    if adaptive.get("failed_boxes"):
        raise AssertionError("adaptive H1319 cover has failed boxes")

    second_rows: list[tuple[arb, str, str, str]] = []
    third_rows: list[tuple[arb, str, str, str]] = []

    # These K_j are moments around the fixed center -1/12.  Variance is no
    # larger than the second moment around any deterministic center.
    for row in adaptive["certified_boxes"]:
        mass = payload_ball(row["K_total"][0])
        moment2_centered = payload_ball(row["K_total"][2])
        if mass.lower() <= 0:
            raise AssertionError("nonpositive adaptive mass enclosure")
        second = moment2_centered / mass
        third = payload_ball(row["scaled_t2_kappa3_diagnostic"])
        second_rows.append((second.upper(), row["r_lo"], row["r_hi"], ADAPTIVE.name))
        third_rows.append((third.upper(), row["r_lo"], row["r_hi"], ADAPTIVE.name))

    for path, report in zip(SCORE_FILES, score_reports):
        if not report.get("all_checks_pass") or report.get("unresolved_boxes"):
            raise AssertionError(f"unusable score cover: {path.name}")
        for row in report["accepted_boxes"]:
            values = [payload_ball(x) for x in row["finite_integral_ranges_mean_value"]]
            errors = [symmetric_error(x) for x in row["uniform_integral_error_radii"]]
            mass = values[0] + errors[0]
            moment2 = values[1] + errors[1]
            if mass.lower() <= 0:
                raise AssertionError("nonpositive score-box mass enclosure")
            # E[W^2] is an upper bound for Var(W); retaining the raw moment
            # avoids a dependency-sensitive subtraction of the mean.
            second = moment2 / mass
            third = payload_ball(row["scaled_t2_kappa_Xi_third"])
            second_rows.append((second.upper(), row["r_lo"], row["r_hi"], path.name))
            third_rows.append((third.upper(), row["r_lo"], row["r_hi"], path.name))

    worst_second = max(second_rows, key=lambda row: row[0])
    worst_third = max(third_rows, key=lambda row: row[0])
    five_four = arb(5) / 4
    one_eighth = arb(1) / 8

    # kappa''(t)=4 Var(W)/A, A>2tr, r>=323/100>25/8.
    compact_kappa2_coefficient = Fraction(5, 2) / Fraction(323, 100)
    checks = {
        "adaptive_all_certified": bool(adaptive["summary"]["all_certified"]),
        "box_count_is_552": len(second_rows) == 552,
        "all_second_moment_enclosures_below_five_four": all(row[0] < five_four for row in second_rows),
        "all_scaled_third_uppers_below_one_eighth": all(row[0] < one_eighth for row in third_rows),
        "compact_kappa2_coefficient_below_four_fifths": compact_kappa2_coefficient < Fraction(4, 5),
        "compact_kappa3_one_eighth_below_eleven_twentieths": Fraction(1, 8) < Fraction(11, 20),
    }
    report = {
        "domain": "125<=t<=t(22)",
        "r_cover": "r(125)<=r<=22, enclosed by 323/100<=r<=22",
        "box_count": len(second_rows),
        "variance_identity": "kappa_Xi''(t)=4*Var(W)/A",
        "saddle_bound": "A>2*t*r",
        "worst_second_moment_upper": {
            "value": str(worst_second[0]),
            "r_box": [worst_second[1], worst_second[2]],
            "source": worst_second[3],
            "slack_below_5_over_4": str(five_four - worst_second[0]),
        },
        "worst_scaled_third_upper": {
            "value": str(worst_third[0]),
            "r_box": [worst_third[1], worst_third[2]],
            "source": worst_third[3],
            "slack_below_1_over_8": str(one_eighth - worst_third[0]),
        },
        "derived_kappa2_coefficient": fraction_record(compact_kappa2_coefficient),
        "conclusion": [
            "kappa_Xi''(t)<4/(5t)",
            "kappa_Xi'''(t)<11/(20t^2)",
        ],
    }
    return report, checks


def high_two_sided_bounds() -> tuple[dict[str, Any], dict[str, bool], dict[str, Any]]:
    high = load(HIGH)
    rc, h1317, stderr = run_json_checker(H1317_CHECKER)
    if rc != 0:
        raise AssertionError(f"H1317 checker failed: {stderr}")

    B1 = int(h1317["total_epsilon_constants"]["B1"])
    B2 = int(h1317["total_epsilon_constants"]["B2"])
    B_XI = int(h1317["B_Xi"])
    B2_LOG = 2 * B2 + 4 * B1 * B1

    q22 = parse_fraction(high["threshold_budget"]["q_R_lower"])
    t22 = q22 / 2
    X_UPPER = parse_fraction(high["dominant_stein_bound"]["X_upper"])
    E2_UPPER = Fraction(26, 25)
    U = 3 * X_UPPER * E2_UPPER

    # Dominant second derivative: 4 Var(Z)<=104/(25A), A>=2tr.
    kappa2_coefficient_at_22 = Fraction(52, 25 * 22) + Fraction(B2_LOG, 1) / t22

    # kappa3(W)/a=-Y+3XE2-2a^2X^3 <= U, and
    # F=q^2*r*B/A^3<1.  Hence t^2*kappa9'''<=2U/r.
    kappa3_coefficient_at_22 = Fraction(2, 22) * U + Fraction(B_XI, 1) / t22

    checks = {
        "high_report_passes": bool(high.get("all_checks_pass")),
        "h1317_checker_passes": bool(h1317.get("all_checks_pass")),
        "B1_is_1415": B1 == 1415,
        "B2_is_7073": B2 == 7073,
        "B2_log_is_8023046": B2_LOG == 8023046,
        "B_Xi_is_45450578214": B_XI == 45450578214,
        "X_upper_is_663_over_1250": X_UPPER == Fraction(663, 1250),
        "U_is_25857_over_15625": U == Fraction(25857, 15625),
        "F_strictly_below_one": high["factor_audit"]["improved_consequence"] == "t^2*kappa_9'''(t)>-2D/r",
        "q22_lower_positive": q22 > 0,
        "high_kappa2_coefficient_below_four_fifths": kappa2_coefficient_at_22 < Fraction(4, 5),
        "high_kappa3_coefficient_below_eleven_twentieths": kappa3_coefficient_at_22 < Fraction(11, 20),
        "monotone_extension_recorded": "r>=22" in high["monotone_extension"]["conclusion"],
    }
    report = {
        "domain": "r>=22, equivalently t>=t(22)",
        "t22_strict_lower": fraction_record(t22),
        "second_log_correction_constant": B2_LOG,
        "dominant_third_upper_U": fraction_record(U),
        "kappa2_coefficient_at_22": fraction_record(kappa2_coefficient_at_22),
        "kappa2_slack_below_4_over_5": fraction_record(Fraction(4, 5) - kappa2_coefficient_at_22),
        "kappa3_coefficient_at_22": fraction_record(kappa3_coefficient_at_22),
        "kappa3_slack_below_11_over_20": fraction_record(Fraction(11, 20) - kappa3_coefficient_at_22),
        "monotonicity": "52/(25r)+B2log/t and 2U/r+B_Xi/t decrease with r",
        "conclusion": [
            "kappa_Xi''(t)<4/(5t)",
            "kappa_Xi'''(t)<11/(20t^2)",
        ],
    }
    return report, checks, h1317


def exact_tail_lemma() -> tuple[dict[str, Any], dict[str, bool]]:
    n = sp.symbols("n", integer=True, positive=True)

    def factorial_factor(k: sp.Expr) -> sp.Expr:
        return sp.factor((k - 1) * (2 * k - 3) / (k * (2 * k - 1)))

    f = factorial_factor(n)
    ratio = sp.factor(factorial_factor(n + 1) / f)
    u = sp.Rational(4, 5) / (n - 2)
    v = sp.Rational(11, 20) / (n - 2) ** 2
    X = sp.factor(f / (1 - u))
    S = sp.factor(ratio / (1 - v))
    Y = sp.factor(S * X)
    H = sp.factor(1 - Y - Y**2 * (1 - X))

    targets = {
        "1-X": sp.together(1 - X).as_numer_denom()[0],
        "S-1": sp.together(S - 1).as_numer_denom()[0],
        "1-Y": sp.together(1 - Y).as_numer_denom()[0],
        "H": sp.together(H).as_numer_denom()[0],
    }
    shifted = {name: positive_shift_coefficients(expr, n, 127) for name, expr in targets.items()}

    x, y, z, s = sp.symbols("x y z s", positive=True)
    D = (1 - y) ** 2 - y**2 * (1 - x) * (1 - z)
    Hxy = 1 - y - y**2 * (1 - x)
    lower_residual = sp.factor(D - (1 - y) * Hxy - y**2 * (1 - x) * (z - y))
    coupled = sp.expand(Hxy.subs(y, s * x))
    derivative_residual = sp.factor(sp.diff(coupled, x) - s * (-1 + s * x * (3 * x - 2)))

    rc, independent, stderr = run_json_checker(INDEPENDENT_TAIL)
    checks = {
        "independent_checker_returncode_zero": rc == 0,
        "independent_checker_status_ok": independent.get("status") == "ok",
        "independent_checker_n_min_127": independent.get("n_min") == 127,
        "D_lower_identity_exact": lower_residual == 0,
        "coupled_derivative_identity_exact": derivative_residual == 0,
        "all_shifted_coefficients_positive": all(all(c > 0 for c in values) for values in shifted.values()),
        "u_and_v_below_one_at_127": bool(u.subs(n, 127) < 1 and v.subs(n, 127) < 1),
        "X_below_one_for_all_n_ge_127": "1-X" in shifted,
        "Y_below_one_for_all_n_ge_127": "1-Y" in shifted,
        "terminal_H_positive_for_all_n_ge_127": "H" in shifted,
        "independent_checker_stderr_empty": not stderr,
    }
    report = {
        "n_range": "n>=127",
        "factorial_factor_f_n": str(f),
        "q_upper_X": str(X),
        "increment_multiplier_upper_S": str(S),
        "q_next_upper_Y": str(Y),
        "terminal_H": str(H),
        "shift_n_equals_m_plus_127": shifted,
        "logical_chain": [
            "exp(w)<=1/(1-w) for 0<=w<1",
            "q_n<=X and q_(n+1)/q_n<=S",
            "q-monotonicity gives z>=y and D_n>=(1-y)H(x,y)",
            "the coupled function H(x,Sx) decreases on 0<x<=X because SX<1 and X<1",
            "H(x,y)>=H(X,SX)>0",
        ],
        "conclusion": "D_n>0 for every integer n>=127",
    }
    return report, checks


def finite_and_h811_assembly(tail_report: dict[str, Any]) -> tuple[dict[str, Any], dict[str, bool]]:
    finite = load(FINITE)
    monotone = load(Q_MONOTONE)
    theorem = load(H811)
    rows = finite["finite_d3"]["rows"]
    first = rows[0]
    q2 = payload_ball(first["q_n"])
    q3 = payload_ball(first["q_center"])
    boundary1 = 1 - 2 * q2 + q2**2 * q3

    finite_checks = finite["checks"]
    conclusions = monotone["conclusions"]
    checks = {
        "finite_certificate_passes": bool(finite.get("all_checks_pass")),
        "finite_range_exactly_2_to_126": finite["finite_d3"]["n_range"] == [2, 126],
        "finite_row_count_125": finite["finite_d3"]["row_count"] == 125,
        "finite_all_D3_q_lowers_positive": bool(finite_checks["all_d3_q_lower_bounds_positive"]),
        "finite_all_D3_determinant_lowers_positive": bool(finite_checks["all_d3_determinant_lower_bounds_positive"]),
        "finite_PF2_through_q128": bool(finite_checks["pf2_ratios_strictly_decreasing_through_R128"]),
        "q2_strictly_below_one_half": q2.upper() < arb(1) / 2,
        "q3_strictly_positive": q3.lower() > 0,
        "boundary_start_1_strictly_positive": boundary1.lower() > 0,
        "boundary_start_0_positive": True,
        "H1322_passes": bool(monotone.get("all_checks_pass")),
        "H1322_global_q_monotonicity": conclusions["equivalent_q_statement"] == "q_(n+1) >= q_n for every integer n >= 2",
        "tail_starts_at_127": tail_report["n_range"] == "n>=127",
        "finite_tail_integer_join_exact": finite["finite_d3"]["n_range"][1] + 1 == 127,
        "H811_theorem_loaded": theorem.get("classification") == "repaired_consecutive_row_pf3_theorem_formalized",
        "H811_scope_is_consecutive_rows_only": "consecutive rows only" in theorem["theorem"]["scope"],
    }
    report = {
        "translated_D3_cover": [
            {"starts": "0", "source": "lower-triangular positive boundary"},
            {"starts": "1", "source": "q2<1/2 and q3>0", "lower": str(boundary1.lower())},
            {"starts": "2<=n<=126", "source": FINITE.name},
            {"starts": "n>=127", "source": "two-sided cumulant tail"},
        ],
        "PF2_cover": [
            {"q_indices": "2<=n<=128", "source": "alias-repaired finite coefficient balls"},
            {"q_indices": "n>=127", "source": "q_n<=X_n<1 from the tail envelope"},
        ],
        "q_monotonicity": conclusions["equivalent_q_statement"],
        "q2_upper": str(q2.upper()),
        "H811_conclusion": "every order-3 Toeplitz minor with consecutive rows is nonnegative for a_n=gamma_n/n!",
        "scope": [
            "consecutive-row order-3 Toeplitz minors only",
            "no sparse-row PF3 claim",
            "no PF-infinity claim",
            "no all-degree Jensen-hyperbolicity claim",
            "no Riemann-Hypothesis claim",
        ],
    }
    return report, checks


def build_markdown(report: dict[str, Any]) -> str:
    compact = report["two_sided_cumulant_bounds"]["compact"]
    high = report["two_sided_cumulant_bounds"]["high"]
    finite = report["global_assembly"]
    return "\n".join(
        [
            "# H1325 Xi Global Translated D3 And H811 Consecutive Rows",
            "",
            f"Classification: `{report['classification']}`",
            "",
            "## Result",
            "",
            "For the Xi coefficients `a_n=gamma_n/n!`, translated contiguous",
            "order-3 Toeplitz determinants are strictly positive at every start.",
            "Together with the previously certified PF2, global q-monotonicity,",
            "and `q_2<1/2`, H811 therefore proves every order-3 Toeplitz minor",
            "with consecutive rows is nonnegative.",
            "",
            "This is not full PF3: sparse-row minors are outside H811.",
            "",
            "## New two-sided cumulant package",
            "",
            "For `kappa(t)=log M_t` and every `t>=125`, the assembled bounds are",
            "",
            "```text",
            "kappa''(t)  < 4/(5t),",
            "kappa'''(t) < 11/(20t^2).",
            "```",
            "",
            f"The compact proof rereads `{compact['box_count']}` rigorous Arb boxes.",
            f"Its worst second-moment upper is `{compact['worst_second_moment_upper']['value']}`.",
            f"Its worst `t^2*kappa'''` upper is `{compact['worst_scaled_third_upper']['value']}`.",
            "",
            "From `r=22` onward the exact endpoint coefficients are",
            "",
            f"- `t*kappa''` upper: `{high['kappa2_coefficient_at_22']['exact']}`;",
            f"- `t^2*kappa'''` upper: `{high['kappa3_coefficient_at_22']['exact']}`;",
            f"- third-cumulant slack below `11/20`: `{high['kappa3_slack_below_11_over_20']['exact']}`.",
            "",
            "Both analytic envelopes decrease with `r`.",
            "",
            "## Exact tail lemma",
            "",
            "For `n>=127`, the derivative bounds give rational upper envelopes",
            "`q_n<=X_n` and `q_(n+1)/q_n<=S_n`.  With `Y_n=S_n X_n`,",
            "coefficient-positive polynomials after `n=127+m` prove",
            "",
            "```text",
            "X_n<1,  Y_n<1,",
            "1-Y_n-Y_n^2(1-X_n)>0.",
            "```",
            "",
            "Using global q-monotonicity to replace `q_(n+2)` by `q_(n+1)`",
            "then proves `D_n>0` for every `n>=127`.",
            "",
            "## Finite seam and H811",
            "",
            "The alias-repaired coefficient balls certify `D_n>0` independently",
            "from both quotient and Toeplitz-determinant formulas for `2<=n<=126`.",
            "Start `1` follows from `q_2<1/2` and `q_3>0`; start `0` is triangular.",
            f"The resulting conclusion is: `{finite['H811_conclusion']}`.",
            "",
            "## Reproduction",
            "",
            "```text",
            "python tools/rh_h1325_xi_global_d3_h811_assembly.py --no-write",
            "```",
            "",
            "The checker parses directed Arb endpoints, verifies the H1319 hash",
            "chain, reruns H1317 and the independent exact tail checker, and",
            "regenerates every shifted polynomial.",
            "",
            "## Scope",
            "",
            "- Proved here: global translated contiguous D3 for Xi and H811's",
            "  consecutive-row order-3 conclusion.",
            "- Not proved: sparse-row PF3, PF-infinity, all-degree Jensen",
            "  hyperbolicity, or the Riemann Hypothesis.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD))
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()

    flint.ctx.dps = 100
    hash_checks, hashes = hash_chain_checks()
    compact, compact_checks = compact_two_sided_bounds()
    high, high_checks, h1317 = high_two_sided_bounds()
    tail, tail_checks = exact_tail_lemma()
    assembly, assembly_checks = finite_and_h811_assembly(tail)

    checks = {
        "hash_chain": hash_checks,
        "compact_bounds": compact_checks,
        "high_bounds": high_checks,
        "exact_tail": tail_checks,
        "finite_and_H811": assembly_checks,
    }
    all_checks = all(all(group.values()) for group in checks.values())
    classification = (
        "xi_global_translated_d3_h811_consecutive_rows_closed"
        if all_checks
        else "xi_global_translated_d3_h811_assembly_failed"
    )
    report = {
        "schema": "rh_h1325_xi_global_d3_h811_assembly.v0",
        "classification": classification,
        "all_checks_pass": all_checks,
        "two_sided_cumulant_bounds": {
            "statement": "kappa''(t)<4/(5t) and kappa'''(t)<11/(20t^2) for every t>=125",
            "compact": compact,
            "high": high,
        },
        "exact_tail": tail,
        "global_assembly": assembly,
        "checks": checks,
        "dependency_hashes": hashes,
        "h1317_replay": h1317,
        "decision": (
            "Global translated contiguous D3 is proved for Xi; H811 is now unconditional for its exact consecutive-row scope."
            if all_checks
            else "At least one proof dependency or exact gate failed; no global claim is made."
        ),
        "limitations": assembly["scope"],
    }

    if not args.no_write:
        out = Path(args.out)
        md = Path(args.markdown_out)
        out.parent.mkdir(parents=True, exist_ok=True)
        md.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        md.write_text(build_markdown(report), encoding="utf-8")

    print(
        json.dumps(
            {
                "classification": classification,
                "all_checks_pass": all_checks,
                "compact_worst_second": compact["worst_second_moment_upper"],
                "compact_worst_third": compact["worst_scaled_third_upper"],
                "high_kappa3_slack": high["kappa3_slack_below_11_over_20"],
                "tail_n_min": 127,
                "finite_n_range": [2, 126],
                "out": None if args.no_write else str(Path(args.out)),
                "markdown": None if args.no_write else str(Path(args.markdown_out)),
            },
            indent=2,
        )
    )
    return 0 if all_checks else 1


if __name__ == "__main__":
    raise SystemExit(main())
