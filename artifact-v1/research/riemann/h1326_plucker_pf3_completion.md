# H1326 Plucker Completion From Consecutive Rows To PF3

Classification: `h1326_strict_pf2_plus_h811_implies_pf3_closed`

## Result

Let `a_n>0` for `n>=0`, put `a_n=0` for `n<0`, and let

```text
T=(a_(j-i))_(i,j>=0).
```

Assume:

1. `R_n=a_n/a_(n-1)` is strictly decreasing;
2. every order-3 minor of `T` with three consecutive rows is nonnegative.

Then every minor of `T` of order at most 3 is nonnegative. In other
words, the sequence `(a_n)` is PF3.

Applied to the Xi coefficients `a_n=gamma_n/n!`, H1325 supplies strict
PF2 and the consecutive-row conclusion of H811. Therefore the Xi
coefficient sequence is PF3.

## Notation

For increasing triples `R=(r0,r1,r2)` and `C=(c0,c1,c2)`, write

```text
Delta_R^C = det(a_(cj-ri))_(0<=i,j<=2).
```

All identities below keep the column triple `C` fixed. A minor with an
added row `d=c1+1` satisfies, for `i<j<d`,

```text
Delta_(i,j,d)^C
  = a_(c2-c1-1)
    det [[a_(c0-i), a_(c1-i)],
         [a_(c0-j), a_(c1-j)]].                 (1)
```

Indeed the row indexed by `d` is
`(0,0,a_(c2-c1-1))`. The coefficient in (1) is strictly positive,
including the edge `c2=d`, where it is `a_0`.

## Strict PF2 for arbitrary gaps

The strict decrease of the adjacent ratios gives every nonzero order-2
Toeplitz minor a strict positive sign, not only adjacent ones.

For rows `i<j` and columns `p<q`, the minor is

```text
a_(p-i)a_(q-j) - a_(q-i)a_(p-j).                (2)
```

If `i<=p<j<=q`, its lower-left entry is zero and its diagonal entries
are positive, so (2) is strictly positive. If the first column or the
lower row is structurally zero, the minor is zero.

In the interior `p>=j`, put

```text
x=p-j,  alpha=j-i,  beta=q-p.
```

Then (2) becomes

```text
a_(x+alpha)a_(x+beta) - a_(x+alpha+beta)a_x

= a_x a_(x+beta)
  [ product_(t=1)^alpha R_(x+t)
    - product_(t=1)^alpha R_(x+beta+t) ].        (3)
```

Here `alpha,beta>0`. Each factor in the first product of (3) is
strictly larger than the corresponding factor in the second product.
Thus (3) is strictly positive.

This explicitly covers the boundary used later: when the first lower
index is zero, the lower-left entry is `a_-1=0`, and positivity is
triangular rather than an interior ratio comparison.

## Structural boundary cases

Before any division, three one-sided Toeplitz cases are removed.

- If `c2<r2`, the last row is zero.
- If `c1<r2<=c2`, expansion along the last row gives

  ```text
  Delta_R^C
    = a_(c2-r2)
      det [[a_(c0-r0), a_(c1-r0)],
           [a_(c0-r1), a_(c1-r1)]] >= 0.
  ```

- If `c1>=r2` but `c0<r1`, expansion along the first column gives

  ```text
  Delta_R^C
    = a_(c0-r0)
      det [[a_(c1-r1), a_(c2-r1)],
           [a_(c1-r2), a_(c2-r2)]] >= 0,
  ```

  with value zero when `c0<r0`.

Thus the only inductive region is

```text
c1>=r2 and c0>=r1.                              (4)
```

In (4), `d=c1+1` is larger than every target row, (1) applies, and all
partners used as denominators below are strictly positive by the
arbitrary-gap PF2 argument.

## Lemma 1: a pair of adjacent rows is enough

First consider `R=(r,r+1,s)`. For `s>r+2`, the fixed-column Plucker
identity on the five rows `r,r+1,r+2,s,d` is

```text
Delta_(r,r+1,s) Delta_(r+1,r+2,d)
 = Delta_(r,r+1,r+2) Delta_(r+1,s,d)
 + Delta_(r,r+1,d)   Delta_(r+1,r+2,s).         (5)
```

The denominator in (5) is strictly positive. At the edge
`c0=r+1`, this is the triangular PF2 case with lower-left entry
`a_-1=0`; for `c0>=r+2`, it is the interior strict-PF2 case.

The first factor on the right of (5) has consecutive rows. The two
minors containing `d` factor through PF2 by (1). The last factor has
the same first adjacent pair and a row gap smaller by one. Induction
from the consecutive base `s=r+2` proves every such minor nonnegative.

For `R=(r,s-1,s)`, the oriented identity is

```text
Delta_(r,s-1,s) Delta_(s-2,s-1,d)
 = Delta_(r,s-2,s-1) Delta_(s-1,s,d)
 + Delta_(r,s-1,d)   Delta_(s-2,s-1,s).         (6)
```

The same reasoning, now inducting on the first row gap, proves (6).
Consequently every order-3 minor having at least one adjacent pair of
rows is nonnegative.

No unknown minor occurs on the wrong side of (5) or (6), and no
possibly zero minor is divided out.

## Lemma 2: arbitrary sparse rows

Let `R=(r,u,s)` with `r<u<s`. The adjacent cases are Lemma 1. In the
remaining case `r+1<u<s`, Plucker on the five rows
`r,r+1,u,s,d`, with common row `u`, gives the exact orientation

```text
Delta_(r,u,s) Delta_(r+1,u,d)
 = Delta_(r,r+1,u) Delta_(u,s,d)
 + Delta_(r,u,d)   Delta_(r+1,u,s).             (7)
```

Under (4), the partner `Delta_(r+1,u,d)` is strictly positive. The two
other minors containing `d` are nonnegative PF2 factors by (1).
The minor `Delta_(r,r+1,u)` is covered by Lemma 1, while
`Delta_(r+1,u,s)` has first gap

```text
u-(r+1) < u-r.
```

Induction on `u-r`, starting from Lemma 1, proves every order-3 minor
nonnegative. Together with coefficient positivity and strict PF2, this
proves PF3.

## Redundant row/column duality audit

For any integer `L>=max(r2,c2)`, set

```text
R*=(L-c2,L-c1,L-c0),
C*=(L-r2,L-r1,L-r0).
```

The dual matrix is the transpose of the original minor with both row
and column orders reversed. Hence

```text
Delta_R^C = Delta_(R*)^(C*).                    (8)
```

Column gaps become the reversed row gaps. Thus Lemma 1 also proves all
minors with an adjacent column pair. This duality is not needed after
the full induction (7), but independently checks the large intermediate
reduction.

## Xi instantiation

The H1325 dependency addendum pins the following inputs.

- H908 gives `a_n>0` for every `n` from the nonzero positive Xi-kernel
  moment representation.
- H1325 certifies `q_n<1` through `n=128` and
  `q_n<=X_n<1` for every `n>=127`. Since
  `q_n=R_n/R_(n-1)`, the ratios `R_n` are strictly decreasing globally.
- H1325 closes translated D3 globally, and H811 then proves every
  consecutive-row order-3 minor nonnegative.

All hypotheses of the abstract theorem are therefore satisfied:

```text
(gamma_n/n!)_(n>=0) is PF3.                     (9)
```

## Fail-closed replay

```text
python tools/rh_h1326_plucker_pf3_completion.py --full-report
python tools/rh_h1326_strict_pf2_multigap_audit.py
```

The first command pins five dependency hashes, expands all three
Plucker orientations symbolically, checks the boundary factorizations
and duality, and replays the induction on every exact finite target.
It reports:

```text
35/35 checks passed
H813 active sparse minors: 2310/2310 proved
H816 residual forms:       356/356 proved
```

Before the final general induction, the adjacent-row lemma proves 242
of the 356 H816 residuals. Row/column duality raises this to 344, and
already covers all 125 six-term forms. The remaining 12 forms have all
four row/column gaps at least two; (7) closes all 12.

The second command checks 512 exact arbitrary-gap instances of (3) and
classifies 3185 one-sided order-2 index configurations, including 924
strict triangular and 1155 strict interior cases.

## Dependency hashes

```text
h1325_xi_global_d3_h811_consecutive_rows.json
  9308258F232A84CB964F0CE4357C2BE279D5E2FE3A49B227B37CCD7E6D43A1A0
h1325_h811_dependency_addendum.json
  F034A24C08EED9FC8E0A0405D596DF47A4227CFBEC7E7664D2823062A45A26B4
h811_repaired_consecutive_row_theorem.json
  09C0B09673CD582BB1F9968466B5304FBA13A4B8AE0D97B6348EAD60498D6A12
h813_sparse_row_symbolic_hierarchy.json
  C92E94BAED6251E729F3B3C02AE5B93A92653118022B32682E7D8B7E24EE5ED6
h816_d3_cone_sparse_row_coupling.json
  8E6F0364C259FBBBC1AB66D9B1C250D2F706004FE1842644AFC584D07BDC2A00
```

## Scope

H1326 closes PF3 for the Xi coefficient sequence. It does **not** prove
PF4, PF-infinity, all-degree Jensen hyperbolicity, or the Riemann
Hypothesis. Exact countermodels from H1320/H1323 already show that a
new Xi-specific invariant is required beyond degree 3.
