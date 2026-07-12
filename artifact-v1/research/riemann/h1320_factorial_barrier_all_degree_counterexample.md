# H1320 Exact Counterexample To A Generic H803/H804-to-Jensen Bridge

Classification:
`factorial_and_eventual_unit_barriers_do_not_imply_all_degree_jensen_hyperbolicity`

## Result

The inequalities closed in H1319 do **not**, as abstract sequence
inequalities, imply all-degree Jensen hyperbolicity or `PF_infinity`.

There is a positive entire coefficient sequence satisfying all of the
following:

```text
R_{n+1} <= R_n,
q_2 = 1/2,
q_{n+1} >= q_n,
D_n = e_{n+1}^2-q_{n+1}^2 e_n e_{n+2} > 0,
the exact H803 factorial barrier for every n>=2,
the stronger H804 unit barrier for every n>=3,
```

but its unshifted Jensen polynomial of degree `3` is nonhyperbolic and its
coefficient Toeplitz matrix has a negative minor of order `4`.

This is an exact rational counterexample, not a floating-point search.

## Construction

Set

```text
a_0=a_1=1,
R_n=a_n/a_{n-1},
q_n=R_n/R_{n-1},
q_2=1/2,
q_n=(2n-3)/(2n),  n>=3.
```

Thus

```text
R_n = 2 binomial(2n-2,n-1)/(n 4^(n-1)),  n>=2.
```

The central-binomial estimate

```text
binomial(2n-2,n-1) <= 4^(n-1)
```

gives `R_n<=2/n` and therefore

```text
a_n <= 2^(n-1)/n!.
```

Consequently

```text
Phi(z)=sum_{n>=0} a_n z^n
```

is entire.  Every coefficient is strictly positive.

Multiplying `a_n` by `A lambda^n`, with `A,lambda>0`, preserves every ratio
condition and every determinant sign.  Hence the normalization `a_0=a_1=1`
is inessential; the first two positive coefficients can be prescribed.

## Exact Ratio And Barrier Checks

For `n>=3`,

```text
q_{n+1}-q_n = 3/(2n(n+1)) > 0.
```

Since every `q_n<1`, the adjacent ratios `R_n` decrease, so the sequence is
`PF_2`.

Put

```text
C_n=(2n)(2n-1),
theta_n=C_{n-1}C_{n+1}/C_n^2.
```

The H803 identity is

```text
T_{n+1}/T_n = theta_n q_{n+1}/q_n.
```

For this sequence and `n>=3`, cancellation gives

```text
T_{n+1}/T_n = (n-1)(2n+1)/(n(2n-1)).
```

The exact factorial-barrier slack is

```text
T_{n+1}/T_n-theta_n
  = 3(n-1)(2n+1)/(n^2(2n-1)^2)
  > 0.
```

At `n=2`, `q_3=q_2`, so

```text
T_3/T_2=theta_2=5/12.
```

Thus the exact H803 barrier holds for every `n>=2`.

The stronger unit-barrier slack is, for every `n>=3`,

```text
T_{n+1}/T_n-(1-1/n^2)
  = (n-1)/(n^2(2n-1))
  > 0.
```

This counterexample therefore satisfies the stronger barrier from `n=3`, far
earlier than the `n>=127` tail currently needed in the Xi assembly.

Finally, with `e_n=1-q_n`,

```text
D_2=13/64>0,
D_n=9(12n-1)/(16n(n+1)^2(n+2))>0,  n>=3.
```

So adding the translated `D3` condition from H851/H877 does not repair the
generic implication.

## Degree-3 Jensen Falsifier

Use exactly the H787 convention

```text
gamma_j=j! a_j,
J_gamma^{3,0}(x)=sum_{j=0}^3 binomial(3,j) gamma_j x^j.
```

The first coefficients are

```text
a_0=1, a_1=1, a_2=1/2, a_3=1/8,
gamma_0=1, gamma_1=1, gamma_2=1, gamma_3=3/4.
```

Hence

```text
J_gamma^{3,0}(x)=1+3x+3x^2+(3/4)x^3.
```

Its exact cubic discriminant is

```text
Disc(J_gamma^{3,0})=-27/16<0.
```

A real cubic with negative discriminant has one real zero and one nonreal
conjugate pair.  Thus this Jensen polynomial is not hyperbolic.

## Independent Toeplitz Falsifier

The first five coefficients are

```text
a_0=1, a_1=1, a_2=1/2, a_3=1/8, a_4=5/256.
```

For Toeplitz rows `(0,1,2,3)` and columns `(1,2,3,4)`, the exact minor is

```text
| 1  1/2  1/8  5/256 |
| 1   1   1/2   1/8  |
| 0   1    1    1/2  | = -5/256 < 0.
| 0   0    1     1   |
```

Therefore the sequence is not `PF_4`, and a fortiori not `PF_infinity`.
This gives an independent algebraic failure certificate that does not use
root finding.

## Consequence For The H907 Route

The following generic bridge is false:

```text
positivity + PF2 + q-monotonicity + q2<=1/2 + translated D3
+ exact factorial barrier + eventual unit barrier
    => all Jensen polynomials hyperbolic / PF_infinity.
```

H1319 remains a valid Xi-specific achievement: it proves a genuine global
ratio inequality for the actual Xi coefficients.  What fails is the proposed
logical jump from those low-order inequalities to all orders.

Any surviving H907 proof must therefore add genuinely Xi-specific all-order
structure, for example a direct `LP+` approximation, an all-order composition
formula, or a hierarchy of minors whose hypotheses are not already exhausted
by the counterexample above.  Merely strengthening the same second- and
third-ratio inequalities cannot close the all-degree gap.

## Reproduction

```text
python tools/rh_h1320_factorial_barrier_all_degree_counterexample_canonical.py
python tools/rh_h1320_factorial_barrier_all_degree_counterexample_canonical.py --check-report
```

The checked report is
`research/riemann/h1320_factorial_barrier_all_degree_counterexample.json`.

## Scope

This is **not** a counterexample for the Riemann Xi coefficients and does not
disprove RH.  It invalidates only a generic all-degree deduction from the
currently proved ratio/barrier package.  A theorem using additional,
independent Xi structure remains logically possible.
