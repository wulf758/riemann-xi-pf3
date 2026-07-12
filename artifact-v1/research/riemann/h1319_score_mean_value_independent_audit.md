# H1319 Independent Audit of the Score Mean-Value Cover

Classification: `independent_audit_pass_not_a_new_curvature_certificate`

## Verdict

**PASS** for the three requested reports.  At engine SHA-256

```text
229BA6EAA67B4E2E48506A08050851FAB40F0E1941A1DA3F1AC560155355DE2F
```

the reports give a gap-free Arb certificate of

```text
t(r)^2 kappa_Xi'''(t(r)) >= -3/4
```

on the union `4<=r<=22`.  The audit did not modify the certificate engine.

| interval | boxes | unresolved | worst certified lower bound | margin over `-3/4` |
|---|---:|---:|---:|---:|
| `[4,5]` | 20 | 0 | `-0.73308000992014794857` | `0.016919990079852051427` |
| `[5,10]` | 100 | 0 | `-0.72966964929450601574` | `0.020330350705493984261` |
| `[10,22]` | 240 | 0 | `-0.51819582344786420226` | `0.23180417655213579774` |

All 360 accepted boxes pass.  Their rational endpoints form an exact chain
from `4` to `22`; there is no floating-point endpoint comparison in this
chain check.

The canonical machine-readable audit is

```text
research/riemann/h1319_score_mean_value_independent_audit_canonical.json
SHA-256 560871B1AB60F15627A917454CC8B48C5A8DCA267301526C333D503588EDD63B
```

## Stable formulas and automatic differentiation

The engine uses first-order dual numbers with Arb coefficients.  Addition,
multiplication, reciprocal, integral powers, exponential, `expm1`, and square
root implement their exact first derivatives.  The only special stabilized
function is

```text
phi2(x)=exp(x)-1-x.
```

Its value is enclosed by the degree-48 Taylor polynomial plus the symmetric
Lagrange remainder

```text
exp(|x|)|x|^49/49!,
```

while its derivative is evaluated independently as

```text
phi2'(x) x_r = expm1(x) x_r.
```

Thus no derivative of an approximate polynomial is substituted for the
derivative of the exact function.

With `P=pi exp(r)`, `h=w/sqrt(A)` and
`delta=r(expm1(h))`, direct symbolic simplification verifies

```text
sqrt(A) (Psi_w-w)
 = A phi2(h)
   +P [r phi2(delta)+delta expm1(delta)].
```

The implemented Xi multiplier and its `w` derivative are also exact for the
finite sum.  For one term

```text
c_m(x) exp(-(m^2-1)x),
c_m(x)=m^4-3m^2/(2x),
```

symbolic differentiation gives exactly the expression used by the engine:

```text
x_w exp(-(m^2-1)x)
 [3m^2/(2x^2)-(m^2-1)c_m(x)].
```

As a separate bug-catching test, the canonical audit runner compares the
value and `r` derivative of ten fields at nine `(r,w)` probes against an
independently written mpmath model.  The fields are `Psi`, `exp(-Psi)`, the
dominant score excess, finite `H`, finite `H_w`, four score integrands, and
the final scale.  All `180/180` comparisons overlap.  This run uses mpmath at
140 decimal digits, formats 115 significant digits, and adds an explicit Arb
radius covering decimal formatting.  These numerical comparisons corroborate
the implementation; the proof itself remains the Arb enclosure.

## Stein identity and the factor 8

Put

```text
g(w)=exp(-Psi(w)) H(w),
N_Xi=Psi_w-H_w/H-w.
```

Then the exact pointwise identity is

```text
g N_Xi = -g'-w g.
```

The boundary terms vanish.  Integrating once gives

```text
E[N_Xi]=-E[W].
```

Multiplying by `w^2+2` before integrating makes the `2w` terms cancel and
gives

```text
E[(W^2+2)N_Xi]=-E[W^3].
```

If `A1=E[N_Xi]`, `C1=E[(W^2+2)N_Xi]` and `V2=E[W^2]`, the centered third
moment is therefore

```text
kappa_3(W)=-C1+3 A1 V2-2 A1^3.
```

The four unnormalized integrals in the code are exactly the numerators and
normalizer needed for this formula.

Finally `q=2t` and `W=sqrt(A)(Z-r)`.  Three `t` derivatives contribute
`2^3`, while three powers of the coordinate scale contribute `A^(-3/2)`.
Consequently

```text
t^2 kappa_Xi'''(t)=8 t^2 A^(-3/2) kappa_3(W),
```

which is the scale used by the engine.

## Mean-value enclosure of the complete finite function

Let `I(r)=(Z,M2,SA,SC)` denote the four finite-`m` integrals, let `K(I)` be
the displayed centered-moment rational function, and put

```text
F_fin(r)=8 t(r)^2 A(r)^(-3/2) K(I(r)).
```

For a rational box `R=[r_lo,r_hi]` with midpoint `r0`, Arb quadrature gives
an enclosure `P` of `I(r0)` and an enclosure `D` of `I'(R)`.  Hence

```text
I(R) subset P+(R-r0)D.
```

The engine first forms this full integral range, then differentiates the
complete function `F_fin`, including both the scale derivative and all four
integral derivatives.  The scalar mean-value theorem then yields

```text
F_fin(R) subset F_fin(r0)+(R-r0)F_fin'(R).
```

This is not endpoint sampling, finite differencing, or a componentwise
surrogate for the final function.

## Omitted `m>8` terms and the perturbation gradient

The omitted central terms are added after the finite mean-value enclosure as
a uniform pointwise perturbation.  This is why no `r` derivative of the
omitted tail is required.  At every `r` in the box, the engine bounds

```text
|H_tail|,
|H_tail,w|,
exp(-Psi)(|H_tail| |N_dom|+|H_tail,w|).
```

It then propagates the four integral error radii through

```text
K(Z,M2,SA,SC)
 =-SC/Z+3 SA M2/Z^2-2 SA^3/Z^3.
```

Symbolic differentiation independently confirms all four implemented
gradient components:

```text
dK/dZ  = SC/Z^2-6 SA M2/Z^3+6 SA^3/Z^4,
dK/dM2 = 3 SA/Z^2,
dK/dSA = 3 M2/Z^2-6 SA^2/Z^3,
dK/dSC = -1/Z.
```

Every expanded mass interval stays strictly positive in all 360 boxes.

## H1272 outer tails and `0<H<1`

H1272 proves, for `r>=5/2`, the global dominant-density tail bounds imported
by the engine.  Its certified global values are strictly below the coarser
targets used here:

```text
T0=2e-40,  T2=2e-37,  T4=7e-35.
```

On `|w|>=24`,

```text
|w| <= w^2/24,
|w|^3 <= w^4/24.
```

The score-tail formulas follow from integrating
`g N_Xi=-g'-wg`.  For the cubic score, integration by parts with
`w^2+2` cancels the linear term exactly, leaving only the boundary density
and the `|w|^3` tail.  The implemented four error radii are therefore valid.

For `X>=pi`, the `m=1` term of `H` is positive and every `m>=2` term is
positive.  To prove `H<1`, use

```text
H <= 1-3/(2X)+sum_(m>=2) m^4 exp(-(m^2-1)X).
```

For every `m>=2`, `X exp(-(m^2-1)X)` decreases on `X>=pi`.  Thus the last
sum is at most `(pi/X)` times its value at `pi`.  The Arb check

```text
sum_(m>=2) m^4 exp(-(m^2-1)pi) < 3/(2pi)
```

then gives the required strict inequality.  Hence the full Xi tail density
is bounded by the dominant density used in H1272.

## Reproducibility and scope

The canonical audit command is

```text
python tools/rh_h1319_score_mean_value_independent_audit_runner_v2.py
```

It exits with code zero and `all_checks_pass=true`.  In addition to the
symbolic and derivative checks, it re-evaluates the worst box of each report;
all three regenerated intervals overlap the stored intervals and pass.
`py_compile` passes for the engine and both audit modules.

One nonblocking hardening item remains: the generic engine CLI does not reject
`r_min<5/2` (or nonpositive `r`), although H1272 starts at `r=5/2` and the
kernel argument uses positive `r`.  This does not affect any audited report,
whose left endpoints are `4`, `5`, and `10`; the CLI should nevertheless add
an explicit domain guard before being reused elsewhere.

This audit does **not** claim a new independent proof of H1272, coverage below
`r=4` or above `r=22`, the separate H920 assembly, the Riemann Hypothesis, or
anything beyond the stated compact curvature interval.
