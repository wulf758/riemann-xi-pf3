# H1319 Full-Xi Analytic High Band From r=59

Classification: `h1319_full_xi_high_band_r59_analytic_closed`

## Result

For the full Xi log-moment and every saddle `r>=59`,

```text
t(r)^2 kappa_Xi'''(t(r)) >= -3/4.
```

The proof is analytic and uses no quadrature or finite differences.

## Dominant Stein estimate

Put `a=B/A^(3/2)`, `X=E[W^2 K]`, and
`Y=E[(W^2+2)W^2 K]`. Stein integration gives

```text
kappa_3(W)/a = -Y + 3 X E[W^2] - 2 a^2 X^3.
```

H1314, H1313 and H927 imply

```text
Y <= (51/100)(7/2+2*26/25),
X <= (51/100)(26/25),
0<a<=1.
```

Discarding the favorable middle term yields

```text
kappa_3(W) >= -D a,
D=6141071619/1953125000
 =3.144228668928.
```

## Factor-seven audit

Because `q=2t`, the exact dominant scaled loss is

```text
2 D q^2 B/A^3.
```

H921 proves `q^2 r B/A^3<=7`, hence the valid bound is

```text
t^2 kappa_9'''(t) >= -14D/r.
```

The tempting `-2D/r` shortcut drops H921's factor seven and is not
used.

## H1316 domain audit

The published H1316 beta-5 statement imposed `q>=exp(100)`, but its
derivative count only needs `x=q/r>30`, `x exp(-r)<4`, and the global
H1313 moments. For every `r>=5/2`,

```text
x=pi exp(r)-9/4-1/r
 >3*12-9/4-2/5
 =667/20>30.
```

Here `exp(5/2)>12` follows from its degree-six Taylor polynomial.
Thus the same beta-5 constants `207,1035,7659` hold globally on
`r>=5/2`. Combined with the already-global arithmetic-tail constants,
denominator safety, and H923, this gives

```text
t^2 G'''(t) >= -45450578214/t.
```

## Exact endpoint budget

At `r=59`, the dominant loss is `42987501333/57617187500`.
Using `pi>3` and `exp(59)>2^59`,

```text
t(59) > 408134212630823828969/8.
```

The transfer loss is below `363604625712/408134212630823828969`, and

```text
total <= 17544670030460320680208676515677/23515545454315044833174804687500
      ~= 0.7460881596195728
       < 3/4.
```

The exact slack is `22997265068990736168106749987/5878886363578761208293701171875`. The same coarse budget
fails at the preceding integer `r=58`, so 59 is the first passing
integer for this particular factor-7 proof.

Both `14D/r` and `B_Xi/t(r)` decrease with `r`, proving the entire
high band `r>=59`.

## Scope

Together with the certified first segment through `r=4`, this reduces
the remaining compact analytic/interval gap to

```text
4 < r < 59.
```

It does not close that gap and does not by itself prove RH.
