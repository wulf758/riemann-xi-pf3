# H1324 Exact Moment-to-Jensen Laguerre Transfer

Classification: `exact_transfer_proved_generic_mixture_shortcut_invalidated`

## Verdict

The moment representation admits an exact all-degree transfer to generalized
Laguerre polynomials.  Let `mu` be a positive Borel measure on `[0,infinity)`
whose moments through order `2(n+d)` are finite, and define

```text
M_m     = integral u^(2m) dmu(u),
a_m     = M_m/(2m)!,
gamma_m = m!*a_m = m!*M_m/(2m)!.
```

For integers `d,n>=0`, with

```text
J_gamma^(d,n)(X)
  = sum_(j=0)^d binom(d,j) gamma_(n+j) X^j,

K_(d,n)
  = sqrt(pi) d! / [4^n Gamma(n+d+1/2)],
```

one has exactly

```text
J_gamma^(d,n)(X)
  = K_(d,n) integral u^(2n)
      L_d^(n-1/2)(-u^2 X/4) dmu(u).                 (H1324.1)
```

This is a useful transfer, but it is not the missing all-degree theorem.
Every single atomic integrand is hyperbolic, while a positive mixture of those
polynomials need not be hyperbolic.  The H790 two-atom measure already gives
an exact quadratic counterexample.  Therefore the remaining gap is an
`Xi`-specific preservation theorem, uniform in `d,n`, not generic positivity
of the representing measure.

## 1. Coefficient-by-Coefficient Proof

Use the generalized Laguerre expansion

```text
L_d^(alpha)(z)
  = sum_(j=0)^d
      Gamma(d+alpha+1)
      ---------------------------------- (-z)^j.
      Gamma(alpha+j+1) (d-j)! j!
```

Set

```text
alpha = n-1/2,
z     = -u^2 X/4.
```

The coefficient of `X^j` on the right of (H1324.1) is then

```text
K_(d,n)
Gamma(n+d+1/2)
------------------------------------ M_(n+j)
Gamma(n+j+1/2) (d-j)! j! 4^j

  = sqrt(pi) d!
    ------------------------------------------ M_(n+j).       (1)
    4^(n+j) Gamma(n+j+1/2) (d-j)! j!
```

Gamma duplication at `m=n+j` gives

```text
Gamma(m+1/2)
  = sqrt(pi) (2m)! / [4^m m!],

sqrt(pi) / [4^m Gamma(m+1/2)]
  = m!/(2m)!.
```

Substituting this into (1) yields

```text
d!                (n+j)!
--------------- * --------- M_(n+j)
(d-j)! j!         (2n+2j)!

  = binom(d,j) gamma_(n+j),
```

which is exactly the coefficient of `X^j` in `J_gamma^(d,n)`.  This proves
(H1324.1) for every coefficient `0<=j<=d`.

The canonical SymPy runner verifies the symbolic quotient after duplication
as exactly `1`.  It also expands SymPy's actual `assoc_laguerre` convention and
checks all 315 coefficients for

```text
0 <= d <= 8,
0 <= n <= 6.
```

The independent runner uses only `math.factorial` and `fractions.Fraction` and
checks 1183 coefficients for `0<=d,n<=12`.

## 2. Specialization to `Psi` and `Xi`

Use the sign-safe notation

```text
Psi(z) = xi(1/2+z),
Xi(t)  = Psi(i t).
```

If the positive kernel is normalized directly by

```text
Psi(z) = integral cosh(zu) dmu_Psi(u),
```

then

```text
M_m = integral u^(2m) dmu_Psi(u)
```

and (H1324.1) applies without any further scale factor.  The cosine form is

```text
Xi(t) = integral cos(tu) dmu_Psi(u),
```

so the alternating Taylor signs of `Xi` do not alter the positive numbers
`gamma_m` used in its Jensen polynomials.

There is a normalization point worth keeping explicit.  In the
de Bruijn-Newman convention used by H909 and H1316,

```text
xi(1/2+z)
  = 8 integral_0^infinity Phi_DBN(U) cosh(2zU) dU.
```

After `u=2U`, the exact measure for (H1324.1) is

```text
dmu_Psi(u) = 4 Phi_DBN(u/2) du,                    (2)
```

and therefore

```text
M_m
  = 8*4^m integral_0^infinity U^(2m) Phi_DBN(U) dU.  (3)
```

Thus the shorthand `dmu(u)=Phi(u)du` is exact if `Phi` denotes the rescaled
`Psi` kernel in (2).  If `Phi` denotes the H909/H1316 kernel, the factors `8`
and `4^m` in (3) must be retained.

## 3. Each Atom Is Hyperbolic

For a single atom `c delta_(u0)` with `c>0` and `u0>0`, (H1324.1) is a
positive scalar multiple of

```text
L_d^(n-1/2)(-u0^2 X/4).
```

Because

```text
n-1/2 > -1,
```

the `d` roots `lambda_(d,k)^(n-1/2)` of the generalized Laguerre polynomial
are positive.  The atomic Jensen roots are therefore

```text
X_k(u0)
  = -4 lambda_(d,k)^(n-1/2) / u0^2,
  k=1,...,d,                                      (4)
```

and are all real and negative.

## 4. Why Positive Mixing Does Not Finish the Proof

Hyperbolicity is not preserved by an arbitrary positive sum of hyperbolic
polynomials.  Here the obstruction occurs inside the exact family (H1324.1).

The H790 example is most naturally written in the variable `x=u^2`:

```text
nu = (9/10) delta_1 + (1/10) delta_16.             (5)
```

Equivalently, in the `u` variable of (H1324.1),

```text
mu = (9/10) delta_1 + (1/10) delta_4.              (6)
```

The first moments are

```text
M_0 = 1,
M_1 = 5/2,
M_2 = 53/2,
```

so

```text
gamma_0 = 1,
gamma_1 = 5/4,
gamma_2 = 53/24.
```

Consequently

```text
J_gamma^(2,0)(X)
  = 1 + (5/2)X + (53/24)X^2
```

has discriminant

```text
Delta
  = (5/2)^2 - 4*(53/24)
  = -31/12 < 0.                                   (7)
```

The same calculation comes directly from the transfer:

```text
K_(2,0) L_2^(-1/2)(-xX/4)
  = 1 + xX + x^2 X^2/12.
```

Integrating this polynomial against (5) gives exactly the quadratic in (7).
Thus every atom in (5) is hyperbolic, but their positive mixture is not.

This invalidates only the generic shortcut

```text
positive measure + atomic Laguerre hyperbolicity
  => hyperbolicity of every mixture.
```

It does not invalidate the desired statement for the specific `Xi` measure.

## 5. Common-Interlacing Audit

Formula (4) shows that changing `u0` dilates the entire atomic root set by
`u0^(-2)`.  For continuous support on `u>0`, the roots move toward
`-infinity` as `u` tends to `0` and toward `0` from below as `u` tends to
`infinity`.

The representation therefore does not itself provide a fixed common
interlacer for the whole atomic family.  The standard sufficient lemma
"a positive combination of polynomials with a common interlacing is
real-rooted" cannot be invoked without a separate proof of its hypothesis.
This observation is not being promoted to a theorem that no interlacing or no
other preservation mechanism can exist.

## 6. Exact Remaining Gap

After H1324, the honest all-degree target is:

```text
Xi-specific Laguerre-mixture theorem:

For the exact measure mu_Psi in (2), prove that

  integral u^(2n) L_d^(n-1/2)(-u^2 X/4) dmu_Psi(u)

is hyperbolic for every d,n>=0.
```

An equivalent all-degree criterion would also suffice.  No such theorem is
proved here.  In particular, H1324 does not establish `PF_infinity`, all
Jensen hyperbolicity, or RH.

## 7. Reproducible Artifacts

```text
python tools/rh_h1324_moment_jensen_laguerre_transfer_canonical.py --check-report
python tools/rh_h1324_moment_jensen_laguerre_transfer_independent_audit.py --check-report
```

Artifacts:

```text
tools/rh_h1324_moment_jensen_laguerre_transfer_canonical.py
tools/rh_h1324_moment_jensen_laguerre_transfer_independent_audit.py
research/riemann/h1324_moment_jensen_laguerre_transfer.json
research/riemann/h1324_moment_jensen_laguerre_transfer_independent_audit.json
```

Both strict runners exit `0` only when every pinned equality and field check
passes.
