"""Precision-safe canonical runner for the H1319 variance extraction.

H1319's JSON payloads retain rigorous directed ``lower`` and ``upper``
endpoints, while the convenience ``interval=x.str(digits)`` field can be much
coarser.  Reconstructing from the endpoint fields preserves the certified
precision needed by the scalar mean-value quotient.
"""

from flint import arb

import claude_h1319_box_variance_upper_certificate as engine
import claude_h1319_box_variance_upper_certificate_canonical  # noqa: F401


def payload_ball(payload: dict[str, str]) -> arb:
    lower = arb(payload["lower"]).lower()
    upper = arb(payload["upper"]).upper()
    return (lower + upper) / 2 + arb(0, 1) * (upper - lower) / 2


engine.payload_ball = payload_ball


if __name__ == "__main__":
    raise SystemExit(engine.main())
