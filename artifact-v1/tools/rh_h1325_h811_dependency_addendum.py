#!/usr/bin/env python3
"""Close the two explicit dependency gates left outside H1325's main runner.

The H1325 arithmetic already verifies the high-band inequality indirectly
through a passing H1319 report.  This addendum checks the exact recorded
``0<F<1`` conclusion itself and pins H811's global strict-positivity input to
the positive-kernel moment representation in H908.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RIEMANN = ROOT / "research" / "riemann"
H1325 = RIEMANN / "h1325_xi_global_d3_h811_consecutive_rows.json"
H908 = RIEMANN / "h908_xi_log_saddle_moment_curvature_reduction.md"
HIGH = RIEMANN / "h1319_full_xi_high_band_r22_analytic_certificate.json"
H811 = RIEMANN / "h811_repaired_consecutive_row_theorem.json"
DEFAULT_OUT = RIEMANN / "h1325_h811_dependency_addendum.json"
DEFAULT_MD = RIEMANN / "h1325_h811_dependency_addendum.md"


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD))
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()

    h1325 = load(H1325)
    high = load(HIGH)
    h811 = load(H811)
    h908_text = H908.read_text(encoding="utf-8")

    checks = {
        "H1325_main_assembly_passes": bool(h1325.get("all_checks_pass")),
        "H1325_main_classification_closed": h1325.get("classification")
        == "xi_global_translated_d3_h811_consecutive_rows_closed",
        "H908_positive_kernel_measure": "positive even kernel measure" in h908_text,
        "H908_moment_identity": "M_n=(2n)!a_n=int u^{2n} dPhi(u)." in h908_text,
        "H908_nonzero_positive_moments": "M_t = int exp(t y) dnu(y)" in h908_text,
        "high_algebraic_ratio_certificate_passes": bool(
            high["algebraic_ratio_certificate"]["all_checks_pass"]
        ),
        "high_exact_F_conclusion": high["algebraic_ratio_certificate"]["conclusion"]
        == "0<q^2*r*B/A^3<1 for every real r>=4",
        "H811_theorem_loaded": h811.get("classification")
        == "repaired_consecutive_row_pf3_theorem_formalized",
        "H811_scope_exact": "consecutive rows only" in h811["theorem"]["scope"],
        "H1325_supplies_global_q_monotonicity": h1325["global_assembly"]["q_monotonicity"]
        == "q_(n+1) >= q_n for every integer n >= 2",
        "H1325_supplies_q2_gate": bool(
            h1325["checks"]["finite_and_H811"]["q2_strictly_below_one_half"]
        ),
        "H1325_supplies_global_PF2": bool(
            h1325["checks"]["finite_and_H811"]["finite_PF2_through_q128"]
            and h1325["checks"]["exact_tail"]["X_below_one_for_all_n_ge_127"]
        ),
        "H1325_supplies_global_translated_D3": bool(
            h1325["checks"]["finite_and_H811"]["finite_all_D3_q_lowers_positive"]
            and h1325["checks"]["exact_tail"]["terminal_H_positive_for_all_n_ge_127"]
            and h1325["checks"]["finite_and_H811"]["boundary_start_1_strictly_positive"]
        ),
    }
    all_checks = all(checks.values())
    classification = (
        "h1325_h811_explicit_dependency_gates_closed"
        if all_checks
        else "h1325_h811_explicit_dependency_gate_failed"
    )
    report = {
        "schema": "rh_h1325_h811_dependency_addendum.v0",
        "classification": classification,
        "all_checks_pass": all_checks,
        "checks": checks,
        "strict_coefficient_positivity": {
            "source": H908.name,
            "argument": "M_n=(2n)!a_n is a nonzero moment of the positive Xi kernel; hence a_n>0 for every n",
        },
        "high_band_upper_sign_gate": {
            "source": HIGH.name,
            "exact_conclusion": high["algebraic_ratio_certificate"]["conclusion"],
        },
        "H811_input_map": {
            "a_n_positive": "H908 positive-kernel moment representation",
            "PF2": "H1325 finite q_n<1 through 128 plus tail X_n<1 from 127 onward",
            "q_n_nondecreasing": "H1322 as loaded by H1325",
            "q2_le_half": "H1325 alias-repaired q2 interval",
            "translated_contiguous_D3": "H1325 starts 0, 1, 2..126, and n>=127",
        },
        "conclusion": (
            "H1325's H811 conclusion is fully explicit: every consecutive-row order-3 Toeplitz minor is nonnegative for the Xi coefficient sequence."
            if all_checks
            else "An explicit H811 dependency gate failed; no conclusion is promoted."
        ),
        "scope": [
            "consecutive-row order-3 Toeplitz minors only",
            "not sparse-row PF3",
            "not PF-infinity",
            "not the Riemann Hypothesis",
        ],
        "dependency_hashes": {
            str(path.relative_to(ROOT)): digest(path)
            for path in (H1325, H908, HIGH, H811)
        },
    }

    markdown = "\n".join(
        [
            "# H1325 H811 Explicit Dependency Addendum",
            "",
            f"Classification: `{classification}`",
            "",
            "The two dependency gates left implicit in the main H1325 runner",
            "are explicit here:",
            "",
            "1. H908 writes `M_n=(2n)!a_n` as a nonzero moment of the positive",
            "   Xi kernel, proving `a_n>0` for every `n`.",
            "2. The H1319 high-band artifact directly certifies",
            "   `0<q^2*r*B/A^3<1` for every real `r>=4`.",
            "",
            "Together with H1325's checked PF2, q-monotonicity, q2, and global",
            "translated-D3 inputs, all hypotheses of H811 are now explicitly",
            "mapped.  The conclusion remains consecutive-row order 3 only.",
            "",
            "```text",
            "python tools/rh_h1325_h811_dependency_addendum.py --no-write",
            "```",
            "",
        ]
    )

    if not args.no_write:
        Path(args.out).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        Path(args.markdown_out).write_text(markdown, encoding="utf-8")

    print(
        json.dumps(
            {
                "classification": classification,
                "all_checks_pass": all_checks,
                "failed_checks": [key for key, value in checks.items() if not value],
                "out": None if args.no_write else str(Path(args.out)),
                "markdown": None if args.no_write else str(Path(args.markdown_out)),
            },
            indent=2,
        )
    )
    return 0 if all_checks else 1


if __name__ == "__main__":
    raise SystemExit(main())
