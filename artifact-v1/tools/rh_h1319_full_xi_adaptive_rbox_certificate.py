"""Canonical entry point for the H1319 adaptive full-Xi r-box certificate.

The implementation lives in ``rh_h1319_full_xi_adaptive_rbox_cover``.  This
entry point keeps the two distinct H1316 dependencies explicit: denominator
safety is in the full assembly note, while the ``19 exp(-3X)`` arithmetic-tail
majorant is in the dedicated m-tail note.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import rh_h1319_full_xi_adaptive_rbox_cover as implementation


H1316_MTAIL_NOTE = Path("research/riemann/h1316_xi_mtail_explicit_transfer.md")


def validate_dependencies() -> dict[str, Any]:
    h1271 = json.loads(implementation.H1271_JSON.read_text(encoding="utf-8"))
    powers = h1271["tail_bounds"]["powers"]
    h1316_text = implementation.H1316_NOTE.read_text(encoding="utf-8")
    h1316_mtail_text = H1316_MTAIL_NOTE.read_text(encoding="utf-8")
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
        "h1316_arithmetic_majorant_present": "19 exp(-3X)" in h1316_mtail_text,
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
        "references": [
            str(implementation.H1271_JSON),
            str(implementation.H1316_NOTE),
            str(H1316_MTAIL_NOTE),
        ],
    }


implementation.validate_dependencies = validate_dependencies


if __name__ == "__main__":
    raise SystemExit(implementation.main())
