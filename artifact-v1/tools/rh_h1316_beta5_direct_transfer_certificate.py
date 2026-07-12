#!/usr/bin/env python3
"""Exact rational arithmetic checks for the H1316 beta=5 transfer bound."""

from fractions import Fraction as Q


def main() -> None:
    m2 = Q(26, 25)
    m4 = Q(7, 2)
    m6 = Q(21)

    assert m2 * m6 < Q(24, 5) ** 2
    assert m4 < Q(15, 8) ** 2
    assert Q(2) < Q(10, 7) ** 2

    c_upper = (
        Q(24, 5)
        + 3 * (2 * m2**2 + m2 * Q(15, 8))
        + 2 * m2**2 * (3 + Q(10, 7))
    )
    assert c_upper == Q(467591, 17500)
    assert c_upper < 27

    zeta5 = Q(8 * 27, 4 * 50**2)
    assert zeta5 == Q(27, 1250)
    assert zeta5 < Q(1, 40)

    print("c_upper=", c_upper)
    print("c_slack_to_27=", 27 - c_upper)
    print("zeta5=", zeta5)
    print("all_checks_pass=true")


if __name__ == "__main__":
    main()
