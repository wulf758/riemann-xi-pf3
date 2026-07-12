# H1321 Rigorous Alias Repair For The Finite H803 Base

Classification: `h1321_alias_repaired_finite_h803_certified`

## Result

The missing Cauchy/DFT alias bound in H801 has been supplied and the true Xi
coefficient inequalities have been recomputed with Arb.  The canonical command

```text
python tools/rh_h1321_h801_alias_repaired_finite_h803_canonical.py
```

exits with status `0` and produces

```text
research/riemann/h1321_h801_alias_repaired_finite_h803.json
```

with

```text
classification = h1321_alias_repaired_finite_h803_certified
all_checks_pass = true
row_indices = 2,...,126
all_q_margins_positive = true
all_h803_margins_positive = true
```

Consequently the exact H803 factorial threshold is certified for the true Xi
coefficients for every integer `2 <= n <= 126`.

## Why A Repair Was Necessary

The old H801 report explicitly said that its finite Arb root-of-unity sums had
no rigorous outer-circle alias bound and that agreement between extraction
radii was only a stability test.  H833 repeated that limitation.  H804 then
used H801 midpoint ratios, and the original H1319 assembly merely loaded the
inherited H801 flags.

Thus the old chain was strong numerical evidence, but its description as an
already-rigorous finite certificate was not justified.  H1321 replaces that
status; it does not silently reinterpret the old flags.

## Exact Coefficient Setup

Use

```text
Xi(z) = xi(1/2+i*z)
      = sum_(n>=0) (-1)^n a_n z^(2n),
a_n   = gamma_n/n! > 0.
```

H801's `N=8192` root-of-unity sum at radius `r=11` encloses the finite DFT
arithmetic with python-flint balls.  For exponent `2n<N`, the exact DFT alias
identity is

```text
a_tilde_n-a_n
  = signed sum_(ell>=1) c_(2n+ell*N) r^(ell*N),
```

where `c_(2n)=(-1)^n a_n` are the Taylor coefficients of `Xi`.

## Rigorous Outer-Circle Majorant

The positive cosine-kernel representation gives, for `|z|=R`,

```text
|Xi(z)|
 <= integral Phi(u) cosh(Ru) du
  = Xi(iR)
  = xi(1/2+R).
```

The first inequality uses

```text
|cos(z u)| <= cosh(|z|u)
```

and the last equality uses the functional equation of `xi`.

At `R=12`, Arb certifies

```text
M_12 = xi(12.5)
     in [10.38273609727472158279337078506768... +/- 4.22e-279].
```

Cauchy's estimate therefore yields the explicit coefficient error

```text
|a_tilde_n-a_n|
 <= M_12 * 12^(-2n)
      * (11/12)^8192 / (1-(11/12)^8192).        (A_n)
```

Arb gives

```text
(11/12)^8192 = 2.72966401203211918041097761314...e-310.
```

Selected alias radii are

```text
n=0:   A_n < 2.835e-309,
n=126: A_n < 3.154e-581,
n=128: A_n < 1.521e-585.
```

Every raw H801 coefficient ball is enlarged symmetrically by `(A_n)` before
any quotient is formed.

## Direct H803 Recalculation

From the repaired coefficient balls the verifier forms

```text
M_n   = (2n)! a_n,
S_n   = M_n/M_(n-1),
tau_n = S_n/S_(n-1),
C_n   = (2n)(2n-1).
```

It checks directly, for all `2 <= n <= 126`,

```text
tau_(n+1)/tau_n >= C_(n-1) C_(n+1)/C_n^2.
```

It independently recomputes

```text
R_n = a_n/a_(n-1),
q_n = R_n/R_(n-1)
```

and verifies `q_(n+1)-q_n>0` on the same range.

Both smallest lower margins occur at `n=126`:

```text
q_(127)-q_126
  > 9.1045134038849920494391992668e-5,

tau_127/tau_126 - C_125*C_127/C_126^2
  > 9.2155750435842900244535214473e-5.
```

## Scope

H1321 repairs only the finite base.  Together with the independently certified
H1319/H920 analytic part for `n>=127`, it closes the H803 ratio-monotonicity
inequality for every `n>=2`.

It does not prove PF3, PF-infinity, all-degree Jensen hyperbolicity, or RH.

## Artifacts

- `tools/rh_h1321_h801_alias_repaired_finite_h803.py`
- `tools/rh_h1321_h801_alias_repaired_finite_h803_canonical.py`
- `research/riemann/h1321_h801_alias_repaired_finite_h803.json`

SHA-256 at certification time:

```text
engine:    9F1527CFE56F758121FDD6B94A68F225734076EFDC64A499557F9F6627A46282
runner:    FC99AC580BAFF7E5C40254BE739D7EE075FE43C695348E6C4E279E8F638B5651
certificate: 75817C6274C23D80044688248CE9540B0DAF3C249BF7E8B4EF8D1703AC9D2F4D
```
