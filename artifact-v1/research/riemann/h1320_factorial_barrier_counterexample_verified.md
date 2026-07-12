# H1320 Verified Exact Counterexample

Classification:
`generic_h803_h804_d3_to_all_degree_bridge_exactly_invalidated`

## Canonical status

The construction and mathematical conclusion of
`h1320_factorial_barrier_all_degree_counterexample.md` survive two repairs to
the original executable canaries.

The canonical certificate is now:

```text
tools/rh_h1320_factorial_barrier_all_degree_certificate_canonical.py
research/riemann/h1320_factorial_barrier_all_degree_certificate.json
```

It reconstructs the exceptional index `n=2` along the full chain

```text
q -> R -> a -> M -> S -> T
```

and obtains

```text
q_2,q_3,q_4 = 1/2,1/2,5/8,
R_1,...,R_4 = 1,1/2,1/4,5/32,
a_0,...,a_4 = 1,1,1/2,1/8,5/256,
M_0,...,M_3 = 1,2,12,90,
S_1,S_2,S_3 = 2,6,15/2,
T_2,T_3 = 3,5/4,
T_3/T_2 = 5/12 = theta_2,
D_2 = 13/64.
```

Thus neither the base factorial equality nor `D_2` is accepted from a copied
constant.

## Independent audit

The second implementation is:

```text
tools/rh_h1320_factorial_barrier_independent_audit.py
research/riemann/h1320_factorial_barrier_independent_audit.json
```

It uses only Python's standard library and `fractions.Fraction`.  It does not
import the SymPy certificate.  From the definition of `q_n`, it independently
reconstructs all coefficients and verifies:

```text
H803 factorial barrier: n=2,...,200,
H804 unit barrier: n=3,...,200,
translated D3 identity: n=2,...,200,
closed ratio/slack identities: n=3,...,200,
R_n closed form and entire majorant: n=2,...,200,
Jensen cubic discriminant: -27/16,
Toeplitz order-4 determinant: -5/256,
all pinned primary JSON fields.
```

Every independent check passes.

## First missing order-4 invariant

For `a_0=a_1=1`, the canonical certificate derives the Toeplitz determinant
symbolically, rather than inserting its quotient formula by hand:

```text
P4(q2,q3,q4)=1-3q2+q2^2(1+2q3)-q2^3 q3^2 q4.
```

At the exact counterexample point,

```text
P4(1/2,1/2,5/8)=-5/256.
```

Thus `P4>=0` is the first explicit order-4 invariant absent from the
H803/H804/D3 package; an all-degree bridge must generate all its analogues.

## Canonical reproduction

```text
python tools/rh_h1320_factorial_barrier_all_degree_certificate_canonical.py --check-report
python tools/rh_h1320_factorial_barrier_independent_audit.py --check-report
```

## Repairs recorded

Two defects found by the audit have been repaired in place:

1. `rh_h1320_factorial_barrier_all_degree_counterexample.py` now converts
   SymPy booleans before `json.dumps`;
2. `rh_h1320_factorial_barrier_all_degree_certificate.py` now constructs
   `theta_2` as an exact rational before comparison.

Both provenance runners now exit successfully as well.  The canonical runner
and the independent Fraction audit remain the supported verification surface
because they include the full reconstructed base and cross-check each other.

## Final logical consequence

The exact counterexample satisfies positivity, `PF_2`, nondecreasing `q_n`,
`q_2<=1/2`, strict translated `D3`, the exact H803 barrier for every `n>=2`,
and the stronger H804 barrier for every `n>=3`.  Nevertheless it fails Jensen
hyperbolicity at degree `3` and total positivity at order `4`.

For the actual Xi sequence, H1321 has now repaired the finite alias gap and
certified H803 for `2<=n<=126`; together with H1319 for `n>=127`, the H803
factorial barrier is therefore closed for every `n>=2`.  This strengthens the
premise but does not supply `P4`, its translates, or any higher-order hierarchy:
H907 remains open.

Therefore the H1319 ratio/barrier package cannot be escalated to all Jensen
degrees by a generic sequence theorem using only those hypotheses.  A valid
H907 continuation must introduce a genuinely new Xi-specific all-order
property.  This remains an invalidation of a proof route, not a counterexample
to Xi or to the Riemann Hypothesis.
