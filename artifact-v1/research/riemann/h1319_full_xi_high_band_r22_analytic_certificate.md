# H1319 Full-Xi Analytic High Band From r=22

Classification: `h1319_full_xi_high_band_r22_analytic_closed`

## Result

For every full-Xi saddle `r>=22`,

```text
t(r)^2 kappa_Xi'''(t(r)) >= -3/4.
```

This is an exact algebraic/rational certificate. It uses no finite
differences and no interval quadrature.

## Stein constant

H1314, H1313 and H927 give

```text
kappa_3(W)>=-D*a,
D=6141071619/1953125000
 =3.144228668928.
```

## Algebraic saddle-ratio certificate

Put `P=pi*exp(r)`, `d=9/(4P)`, and `epsilon=1/(rP)`. Exact
simplification gives

```text
F(r)=q^2*r*B/A^3
    =r(1-d-epsilon)^2
      (r^2+3r+1-d)/(r+1-d)^3.
```

For `r>=4`, the quadratic Taylor lower bound for `exp(r)` and
`pi>3` give `pi*exp(r)>9r`, hence `d<1/(4r)` and `q>0`.
Furthermore

```text
(r+1-d)^3-r(r^2+3r+1-d)
 =2r+1-d(3r^2+5r+3)+d^2(3r+3)-d^3
 >5r/4-1/4-3/(4r)-1/(64r^3).
```

After multiplication by `64r^3` and the shift `r=4+x`, the last
lower bound becomes

```text
80*x**4 + 1264*x**3 + 7440*x**2 + 19328*x + 18687,
```

whose coefficients are all positive. Therefore

```text
0<F(r)<1 for every r>=4.
```

The dominant scaled loss is consequently `2D/r`. H921's older
factor-7 estimate would only give `14D/r`; that factor is explicitly
retained in the fallback and replaced here by the proved inequality
`F<1`.

## Full-Xi transfer provenance

The global correction bound is cited from H1317, not from the older
restricted H1316 assembly:

```text
t^2 G'''(t)>=-B_Xi/t=-2B_Xi/q,
B_Xi=45450578214,
valid for every saddle r>=5/2.
```

## Endpoint r=22

The degree-six Taylor polynomial proves `e>19/7`. Thus

```text
q(22)>1791237414804088073922899305903/7819642097165976098.
```

Exact rational arithmetic gives

```text
dominant < 6141071619/21484375000,
transfer < 710814509485458368763607057944/1791237414804088073922899305903,
total    < 26271522728170960357349911307191337466957/38483616333681579713187289775259765625000
         ~= 0.6826677228142313
          < 3/4.
```

The exact slack is `2591189522090224427540556024253486751793/38483616333681579713187289775259765625000`. The same certified
envelope fails at `r=21`, so 22 is the first passing integer for this
proof. Since `2D/r` decreases and `q(r)` increases, the estimate holds
for all `r>=22`.

## Scope

Combined with the first compact certificate through `r=4`, the
remaining interval is now exactly

```text
4 < r < 22.
```

This does not close that interval and does not prove RH.
