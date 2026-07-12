# H927 Gamma-Log Cubic Coefficient Below One For `r>=2`

Classification: `h927_gamma_log_a_le_one_r_ge_2_proved`

## Question

H925 reduces the direct low-threshold moment box to two score-excess budgets,
but H922 also needs

```text
0<a<=1.
```

Can `a<=1` be proved directly for the whole intended low-threshold range
`r_q>=2`?

## Short Answer

Yes. For the Xi dominant gamma-log model `alpha=9/4`, the normalized cubic
coefficient

```text
a = B/A^(3/2),
A=V''(z_q),
B=V'''(z_q)
```

satisfies

```text
0<a<=1
```

for every saddle with

```text
r=r_q>=2.
```

This removes one of the low-threshold loose ends left in H925. The only
remaining H922 moment-box inequalities are the two H925 score-excess budgets.

## Exact Formulas

Let

```text
P = pi exp(r),
alpha = 9/4.
```

H912 gives

```text
A = r[P(r+1)-alpha],
B = r[P(r^2+3r+1)-alpha].
```

Since `P(r^2+3r+1)>alpha` for `r>=2`, one has

```text
B>0,
```

and hence `a>0`.

## Lower Bound For `A`

We first prove

```text
A >= 9r^2.
```

It is enough to show

```text
P(r+1)-alpha >= 9r.
```

Define

```text
F(r)=pi exp(r)(r+1)-alpha-9r.
```

For `r>=2`,

```text
F'(r)=pi exp(r)(r+2)-9>0,
```

and

```text
F(2)=3pi e^2 - 9/4 - 18 > 0.
```

Therefore `F(r)>0` for every `r>=2`, and

```text
A >= 9r^2.
```

Consequently

```text
sqrt(A) >= 3r.
```

## Upper Bound For `B/A`

We next prove

```text
B <= 3r A.
```

Using the exact formulas,

```text
3rA - B
 = r[ P(3r(r+1)-(r^2+3r+1)) - alpha(3r-1) ]
 = r[ P(2r^2-1) - alpha(3r-1) ].
```

For `r>=2`, define

```text
G(r)=pi e^2(2r^2-1)-alpha(3r-1).
```

Then

```text
G'(r)=4pi e^2 r - 3alpha > 0,
```

and

```text
G(2)=7pi e^2 - 5alpha > 0.
```

Since `P=pi exp(r)>=pi e^2`, it follows that

```text
P(2r^2-1)-alpha(3r-1) >= G(r)>0.
```

Thus

```text
B/A <= 3r.
```

## Cubic Coefficient Bound

Combining the two estimates,

```text
a = B/A^(3/2)
  = (B/A)/sqrt(A)
  <= 3r/(3r)
  = 1.
```

This proves

```text
0<a<=1
```

for every `r>=2`.

## Consequence For H925

H925 showed that the H922 lower-sign moment box is equivalent to:

```text
E[N] <= a,
E[(W^2+2)N] <= 3a.
```

H922a supplies:

```text
E[W]<=0,
E[W^3]<=0.
```

H927 now supplies:

```text
0<a<=1.
```

Therefore, for the direct low-threshold route, only the two score-excess
budgets remain.

## Status

Closed here:

```text
the `a<=1` part of H922 for all r>=2.
```

Still open:

```text
prove or certify the two H925 score-excess budgets for all r>=2.
```

## Route Judgment

This is a small but useful compression. It means the medium-range interval
certificate does not need to spend effort proving `a<=1`; it can focus entirely
on the two expectation budgets.
