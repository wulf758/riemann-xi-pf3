#!/usr/bin/env python3
"""Independent final-scope audit for the H1325 Xi/D3/H811 assembly.

This companion deliberately calls the public proof-building functions from
``rh_h1325_xi_global_d3_h811_assembly`` and adds two scope checks that are
easy to overlook in the main arithmetic assembly:

* the high-band factor is certified directly as ``0<F<1``;
* H811's strict coefficient-positivity hypothesis is sourced from the H908
  positive Xi-kernel moment representation.
"""

from __future__ import annotations

import json
from pathlib import Path

import rh_h1325_xi_global_d3_h811_assembly as primary


def main() -> int:
    hash_checks, hashes = primary.hash_chain_checks()
    compact, compact_checks = primary.compact_two_sided_bounds()
    high, high_checks, _ = primary.high_two_sided_bounds()
    tail, tail_checks = primary.exact_tail_lemma()
    assembly, assembly_checks = primary.finite_and_h811_assembly(tail)

    high_source = primary.load(primary.HIGH)
    h908 = (
        primary.RIEMANN / "h908_xi_log_saddle_moment_curvature_reduction.md"
    )
    h908_text = h908.read_text(encoding="utf-8")
    added_checks = {
        "direct_high_factor_certificate_is_zero_lt_F_lt_one": (
            high_source["algebraic_ratio_certificate"]["conclusion"]
            == "0<q^2*r*B/A^3<1 for every real r>=4"
        ),
        "H908_positive_kernel_measure_recorded": (
            "positive even kernel measure" in h908_text
        ),
        "H908_positive_moment_identity_recorded": (
            "M_n=(2n)!a_n=int u^{2n} dPhi(u)." in h908_text
        ),
        "finite_PF2_input_passes": assembly_checks["finite_PF2_through_q128"],
        "global_q_monotonicity_input_passes": assembly_checks[
            "H1322_global_q_monotonicity"
        ],
        "q2_input_passes": assembly_checks["q2_strictly_below_one_half"],
        "H811_scope_remains_consecutive_rows_only": assembly_checks[
            "H811_scope_is_consecutive_rows_only"
        ],
    }
    groups = {
        "hash_chain": hash_checks,
        "compact": compact_checks,
        "high": high_checks,
        "tail": tail_checks,
        "finite_and_H811": assembly_checks,
        "added_scope_checks": added_checks,
    }
    all_checks = all(all(group.values()) for group in groups.values())
    report = {
        "schema": "rh_h1325_xi_global_d3_h811_independent_audit.v0",
        "classification": (
            "h1325_xi_global_d3_h811_independently_audited"
            if all_checks
            else "h1325_xi_global_d3_h811_independent_audit_failed"
        ),
        "all_checks_pass": all_checks,
        "checks": groups,
        "compact_worst_second": compact["worst_second_moment_upper"],
        "compact_worst_third": compact["worst_scaled_third_upper"],
        "high_kappa3_slack": high[
            "kappa3_slack_below_11_over_20"
        ],
        "strict_coefficient_positivity": {
            "source": str(h908.relative_to(primary.ROOT)).replace("\\", "/"),
            "source_sha256": primary.sha256(h908),
            "reason": (
                "M_n=(2n)!a_n is a strictly positive moment of the nonzero "
                "positive Xi-kernel measure"
            ),
        },
        "audited_conclusion": (
            "global translated D3 and H811 consecutive-row order-3 only"
        ),
        "excluded_claims": [
            "sparse-row PF3",
            "PF-infinity",
            "all-degree Jensen hyperbolicity",
            "the Riemann Hypothesis",
        ],
        "primary_dependency_hashes": hashes,
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if all_checks else 1


if __name__ == "__main__":
    raise SystemExit(main())
