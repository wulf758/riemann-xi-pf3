#!/usr/bin/env python3
"""Exact arithmetic check for H1316's beta=5 H923 constants."""

from fractions import Fraction as Q


def main() -> None:
    counts = (3, 15, 111)
    constants = tuple(69 * n for n in counts)
    assert constants == (207, 1035, 7659)

    b1 = Q(207, 2)
    b2 = Q(1035, 2)
    b3 = Q(7659, 2)
    b_beta5 = 2 * b3 + 12 * b1 * b2 + 16 * b1**3
    assert b_beta5 == 18389880

    print("coefficient_counts=", counts)
    print("R5_derivative_constants=", constants)
    print("B_beta5=", b_beta5)
    print("all_checks_pass=true")


if __name__ == "__main__":
    main()
