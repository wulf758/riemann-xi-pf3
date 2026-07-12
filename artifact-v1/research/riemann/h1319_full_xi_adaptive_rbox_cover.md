# H1319 Full-Xi Adaptive Arb R-Box Cover

Classification: `h1319_full_xi_adaptive_rbox_cover_certified`

## Result

- covered range: `[323/100, 4]`
- target: `t(r)^2 kappa_Xi'''(t(r)) >= -3/4`
- all leaf boxes certified: `True`
- certified boxes: `192`
- failed leaf boxes: `0`
- box attempts: `192`
- worst certified polynomial lower bound: `{'value': 18999.871163750544, 'r_box': ['323/100', '62093/19200'], 'depth': 0}`
- worst scaled diagnostic lower bound: `{'value': -0.7136080274358392, 'r_box': ['323/100', '62093/19200'], 'depth': 0}`

## Exact certificate form

For translated full-Xi moments `K_j`, set
`C=K3*K0^2-3*K1*K2*K0+2*K1^3`.  Each successful box proves

`8*t^2*C + beta*A^(3/2)*K0^3 > 0`.

This is algebraically equivalent to the requested lower bound because
`A>0` and `K0>0`.  No finite differences are used.

## Rigorous enclosures

The `m=1` central integrals use a second-order Arb midpoint rule and a
Taylor-enclosed `phi2(x)=exp(x)-1-x`, avoiding correlated subtraction.
The positive `m>=2` central correction uses the H1316 majorant
`19 exp(-3*pi*exp(rho))`.  Outside `|W|=24`, `0<H<1` reduces the
full-Xi tails to the H1271 dominant tail certificate.

## Scope

This artifact proves only its configured r-range when `all_certified` is
true.  A failed box means the chosen interval discretization did not
separate the target from zero; it is not a counterexample.
