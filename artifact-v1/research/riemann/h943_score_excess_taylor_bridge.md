# H943 Score-Excess Taylor Bridge

Classification: `h943_score_excess_taylor_bridge_exact_identity_budget_open`

## Question

Can the H942 pivot be turned into an analytic target, instead of continuing the
first-order interval-box cover that becomes too expensive for moderate `r`?

## Short Answer

Yes. H925's score-excess

```text
N(w)=Psi'(w)-w
```

has an exact Taylor representation

```text
N(w)=a w^2 K_r(w),
```

where `K_r(w)` is a nonnegative weighted curvature ratio with

```text
K_r(0)=1/2.
```

Thus the two remaining H925 budgets are equivalent to two normalized
expectation bounds with the small coefficient `a` divided out. This gives a
clean moderate-`r` target and explains why H942 is probably a method failure,
not a mathematical obstruction.

## Setup

For `alpha=9/4`, put

```text
V(z)=pi exp(exp(z))-alpha exp(z)-z,
r=exp(z),
q=V'(z),
A=V''(z),
B=V'''(z),
a=B/A^(3/2).
```

For the normalized saddle variable, define

```text
u=w/sqrt(A),
Psi(w)=V(z+u)-V(z)-q u.
```

Then

```text
Psi'(w)=(V'(z+u)-V'(z))/sqrt(A).
```

## Exact Taylor Identity

Taylor's formula with integral remainder, applied to `V'`, gives for every real
`u`

```text
V'(z+u)-V'(z)-A u
  = u^2 int_0^1 (1-t) V'''(z+t u) dt.
```

Therefore

```text
N(w)
 = Psi'(w)-w
 = (V'(z+u)-V'(z)-A u)/sqrt(A)
 = w^2/A^(3/2) int_0^1 (1-t) V'''(z+t w/sqrt(A)) dt.
```

Since `a=B/A^(3/2)`, this is

```text
N(w)=a w^2 K_r(w),
```

with

```text
K_r(w)=int_0^1 (1-t) V'''(z+t w/sqrt(A))/B dt.
```

This identity is exact. It is also sign-stable: H922a/H925 gives `V'''>0`, so

```text
K_r(w)>=0,
N(w)>=0.
```

At the saddle,

```text
K_r(0)=int_0^1 (1-t) dt=1/2.
```

## Budget Equivalence

Let expectation be taken against the normalized density proportional to
`exp(-Psi(w))`. The H925 budgets are

```text
E[N(W)] <= a,
E[(W^2+2)N(W)] <= 3a.
```

Because H927 proves `a>0` for `r>=2`, the exact Taylor identity makes these
equivalent to

```text
E[W^2 K_r(W)] <= 1,
E[(W^2+2)W^2 K_r(W)] <= 3.
```

This is the H943 bridge. It removes the small parameter `a` from the target and
replaces the vague task "control N" by a concrete curvature-ratio problem.

In the Gaussian limiting model one has heuristically `K_r(W)~1/2`, giving

```text
E[W^2 K_r(W)] ~ 1/2,
E[(W^2+2)W^2 K_r(W)] ~ (3+2)/2 = 5/2.
```

So the target has visible slack: `1/2` in the first budget and `1/2` in the
second budget. This is the main reason H942's interval-width collapse is a
good pivot signal rather than a negative mathematical signal.

## Curvature Ratio Facts

Writing `x=exp(y)`, direct differentiation gives

```text
V'''(y)=x[pi exp(x)(x^2+3x+1)-alpha],
V''''(y)=x[pi exp(x)(x^3+6x^2+7x+1)-alpha].
```

For every `x>0` and `alpha=9/4`,

```text
0 < V''''(y)/V'''(y) <= x+3.
```

Indeed the denominator is positive because

```text
pi exp(x)(x^2+3x+1)-alpha >= pi-alpha > 0,
```

and

```text
(x+3)[pi exp(x)(x^2+3x+1)-alpha]
 - [pi exp(x)(x^3+6x^2+7x+1)-alpha]
 = pi exp(x)(3x+2)-alpha(x+2).
```

The last expression is positive for `x>=0`, since at `0` it equals
`2(pi-alpha)>0` and its derivative is `pi exp(x)(3x+5)-alpha>0`.

Consequently `log V'''(y)` is increasing, and if `s>=0`,

```text
V'''(z+s)/B
 <= exp(r(exp(s)-1)+3s).
```

If `s<=0`, then simply

```text
0 < V'''(z+s)/B <= 1.
```

This gives the useful one-sided bound

```text
w <= 0  =>  K_r(w) <= 1/2.
```

The right side `w>0` is the only side that needs real work.

## Lemma Candidates

### H943-A: Exact Taylor Bridge

For the gamma-log saddle with `r>=2`,

```text
N(w)=a w^2 K_r(w),
K_r(w)=int_0^1 (1-t) V'''(z+t w/sqrt(A))/B dt,
K_r(w)>=0,
K_r(0)=1/2.
```

Status: proved in this note by Taylor's formula.

### H943-B: Curvature-Distortion Gate

For `x=exp(y)>0`,

```text
0 < d/dy log V'''(y) <= x+3.
```

Hence

```text
w <= 0 => K_r(w) <= 1/2,
```

and for `w>=0`,

```text
K_r(w)
 <= int_0^1 (1-t)
      exp(r(exp(t w/sqrt(A))-1)+3t w/sqrt(A)) dt.
```

Status: the derivative-ratio inequality is proved in this note; the resulting
expectation bounds remain open.

### H943-C: Moderate-r Score-Excess Budget

Find a practical `R` such as `3`, `4`, or `5` and prove for every `r>=R`

```text
E[W^2 K_r(W)] <= 1,
E[(W^2+2)W^2 K_r(W)] <= 3.
```

Status: open. This is now the preferred analytic patch after H942.

### H943-D: Hybrid Closure

Combine:

```text
H941: 2 <= r <= 2.50 certified by score-corrected interval cover,
H943-C: r >= R certified analytically,
finite bridge: 2.50 <= r <= R certified by sharper quadrature or small cover,
H922m: asymptotic high-q closure retained as backup.
```

Status: open assembly target.

## Route Judgment

This is a better next object than continuing H937 naively. The new target is
pure gamma-log calculus: no Mobius estimates, no PNT input, and no hidden
RH-strength assumption at this stage. The circularity audit must return later
when transferring the gamma-log result back to the full Xi kernel, but H943
itself is safely inside the unconditional dominant-model analysis.

The immediate next step is to prove H943-C or to build a sharper certified
quadrature that subtracts the Gaussian model `a w^2/2` before interval boxing.
