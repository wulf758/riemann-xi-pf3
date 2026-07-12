# H1326 Independent PF3 Pluecker Audit

Classification: `h1326_pf3_plucker_independent_audit_pass`

Result: **PASS** (22 checks; failures: 0).

## Boundary conclusion

The three structural branches exhaust every possible zero-partner case before division. In the remaining induction region, `d=c1+1>s`, the structural tail is positive and the 2x2 partner is strict.  At the adjacent boundary `c0=r+1`, its lower-left entry is `a_-1=0`; strictness comes from the positive triangular diagonal, not from four positive entries.

For `c0>=r+2`, strictness follows from products of strictly decreasing adjacent ratios: `a_(p+delta)/a_p > a_(p+delta+g)/a_(p+g)`.

## Independence canary

The exact q-sequence `['1/2', '7/16', '11/16', '23/32', '9/16', '9/16', '5/16']` is not nondecreasing, while all 210 active consecutive-row minors are positive and all 7056 audited order-3 minors are nonnegative.

## Exact rational grid

The complete grid `['1/4', '1/2', '3/4']` of length 5 checked 243 sequences. Among 48 consecutive-row survivors (37 with nonmonotone q), PF3 failures: 0.

These scans are canaries; the proof-bearing checks are the symbolic Pluecker identities, the exact structural factorization, the exhaustive index partition, and the strict-partner case split.
