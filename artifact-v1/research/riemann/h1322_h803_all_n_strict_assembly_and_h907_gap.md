# H1322 Strict All-Index H803 Assembly And H907 Logic Audit

Classification: `h1322_h803_all_n_certified_h907_open`

## Verdict

The strict command

```text
python tools/rh_h1322_h803_all_n_strict_assembly.py
python tools/rh_h1322_h803_all_n_strict_assembly.py --check-report
```

passes and certifies the H803 factorial threshold for every integer `n>=2`.
It uses H1321, not the old H801/H804 flags, for `2<=n<=126`; it uses the
H1319/H920 analytic chain only for `n>=127`.

The exact consequence is

```text
q_(n+1) >= q_n       for every n>=2.
```

H907 remains open.  H803 does not prove PF3, PF-infinity, all-degree Jensen
hyperbolicity, or RH.

## 1. Definitions Without The Xi Sign Ambiguity

Let

```text
xi(s) = (1/2)s(s-1) pi^(-s/2) Gamma(s/2) zeta(s)
```

and set

```text
Psi(z) = xi(1/2+z)
       = sum_(n>=0) gamma_n z^(2n)/n!,

Xi(t)  = Psi(i t)
       = sum_(n>=0) (-1)^n gamma_n t^(2n)/n!,

F(w)   = Psi(sqrt(w)) = Xi(sqrt(-w))
       = sum_(n>=0) a_n w^n,

a_n    = gamma_n/n! > 0.
```

Thus RH says that the zeros of `Psi` are purely imaginary, equivalently that
the zeros of `F` are real and nonpositive.

For `d,n>=0`, the Jensen polynomial is

```text
J_gamma^(d,n)(X)
  = sum_(j=0)^d binom(d,j) gamma_(n+j) X^j.
```

For the ordinary coefficient sequence `a`, define the one-sided Toeplitz
matrix and its minors by

```text
T(a)_(i,j) = a_(j-i),       a_k=0 for k<0,

Delta(I,J) = det(a_(j_beta-i_alpha))_(1<=alpha,beta<=r).
```

The sequence is `PF_r` when every minor of order at most `r` is nonnegative,
and `PF_infinity` when every finite minor is nonnegative.

The H803 quotient variables are

```text
R_n   = a_n/a_(n-1),
q_n   = R_n/R_(n-1),
e_n   = 1-q_n,

M_n   = (2n)! a_n,
S_n   = M_n/M_(n-1),
tau_n = S_n/S_(n-1),
C_n   = (2n)(2n-1).
```

The symbol `tau_n` is used here for H803's `T_n` so that it is not confused
with the Toeplitz matrix `T(a)`.

## 2. What H803 All-Index Actually Proves

Because

```text
S_n   = C_n R_n,
q_n   = tau_n C_(n-1)/C_n,
```

one has the exact identity

```text
q_(n+1)/q_n
 = [tau_(n+1)/tau_n]
   * C_n^2/[C_(n-1)C_(n+1)].
```

Therefore

```text
tau_(n+1)/tau_n >= C_(n-1)C_(n+1)/C_n^2
```

is exactly `q_(n+1)>=q_n`.  There is no higher-order minor hidden in this
identity.

H1322's gap-free proof split is

```text
2 <= n <= 126:
  H1321 replays the Xi coefficients and adds a rigorous Cauchy/DFT alias
  radius before directly checking all 125 H803 margins.

127 <= n <= 566:
  H1319's compact full-Xi curvature bound and the H920 cube lemma give the
  stronger unit barrier.

n >= 567:
  the compact and H1317 bounds glue with coefficient 141/142; H920 again gives
  the stronger unit barrier.
```

The join is exact: the last finite row is `126`, the first analytic cube is
`n=127` with `t`-range `[125,128]`, and no index is omitted.

## 3. Correction To The Former Finite-Base Status

Before H1321, the finite block was not rigorous for the true Xi coefficients:

- `h801_xi_ratio_logconvexity_audit.md:162-169` calls the result finite
  Cauchy/DFT evidence, states that no fresh rigorous outer-circle alias bound
  was present, and says cross-radius agreement is not a theorem.
- `h833_xi_h811_input_audit.md:32-40` repeats that the promoted window is only
  finite evidence and that H801 has no such alias bound.
- `h804_unit_deficit_barrier_audit.md:119-124` says it uses H803 midpoint ratios
  plus H801 certificates and is not a new Arb proof.
- `h1319_piecewise_h920_conditional_assembly.md:192-193` merely loaded the 125
  inherited H801 flags.

H1321 is the required repair.  H1322 refuses to use the legacy flags as its
finite proof source, records live hashes for all dependencies, and fails closed
if a dependency is missing or any content/join check fails.

## 4. Why H803 Does Not Bridge To PF-Infinity

The insufficiency is algebraic, not a matter of numerical precision.  Consider
the exact entire positive sequence

```text
b_n = 2^(-n(n-1)/2),       n>=0.
```

Then

```text
b_n/b_(n-1) = 2^(-(n-1)),
q_n = 1/2,
e_n = 1/2,
D_n = e_(n+1)^2-q_(n+1)^2 e_n e_(n+2) = 3/16.
```

So this sequence satisfies positivity, PF2 ratio decrease, `q_2<=1/2`, global
q-monotonicity with equality, and even the H877 late-D3 quotient inequality.
Its generating series is entire because `b_n^(1/n)->0`.

Nevertheless the order-four Toeplitz minor with rows `(0,1,2,3)` and columns
`(1,2,3,4)` is

```text
det [[1, 1/2, 1/8, 1/64],
     [1,   1, 1/2,  1/8],
     [0,   1,   1,  1/2],
     [0,   0,   1,    1]]
 = -1/64.
```

Thus even H803 plus all scalar regularity inequalities listed in H877 Lemma 1
cannot imply PF4, hence cannot imply PF-infinity, for general sequences.  Any
valid escalation must use genuinely new Xi-specific all-order structure.

H811 remains correct within its stated scope: under its hypotheses it proves
only order-three minors with consecutive rows.  Its own lines 104-109 warn that
sparse order-three minors remain and that PF3 is far weaker than PF-infinity.

## 5. Smallest Additional Statement That Is Actually Sufficient

Among the standard candidates already isolated in H787, the logically smallest
clean all-degree statement is the unshifted Jensen condition

```text
JP0:
  for every d>=1, J_gamma^(d,0)(X) is hyperbolic.
```

Jensen's theorem then gives `F in LP+`; since the coefficients are positive,
the zeros of `F` are nonpositive; hence RH follows.  Shifts `n>0` are not needed
for this implication.

The equivalent Toeplitz statement is

```text
JP2':
  (a_n)=(gamma_n/n!) is PF_infinity.
```

For this particular entire nonnegative generating function, the
Schoenberg-Edrei/Aissen-Schoenberg-Whitney bridge identifies PF-infinity with
`F in LP+`.  The correction `a_n=gamma_n/n!` is essential; raw `(gamma_n)`
PF-infinity was already invalidated in H781.

There is currently no proved smaller bridge that combines H803 with only
finite-order minors.  Adding JP0 or JP2' makes H803 logically redundant for the
final implication: these statements already have RH strength.  The useful but
open research problem is to derive one of them from an independent structural
mechanism, not merely to rename it.

## 6. Implication And Circularity Audit

| step | status |
| --- | --- |
| H1321 finite balls -> H803 for `2<=n<=126` | certified after analytic alias inflation |
| H1319/H920 -> H803 for `n>=127` | certified by the compact/tail cube assembly |
| H803 all `n` -> `q_(n+1)>=q_n` all `n` | exact equivalence |
| H803 -> global translated D3 | not proved |
| H803 + H811 -> PF3 | false as stated; H811 covers consecutive rows only |
| PF3 -> PF-infinity | false in general and not claimed by H811 |
| H877 Lemma 5 essential-minor closure | open; it is the principal structural gap |
| corrected PF-infinity `JP2'` -> RH | valid for the entire nonnegative `F` above |
| RH -> corrected PF-infinity via the real-zero Hadamard product | valid reverse implication |
| use that RH Hadamard product to prove PF-infinity | circular |
| all unshifted Jensen polynomials hyperbolic -> RH | valid Jensen criterion |
| eventual hyperbolicity for each fixed degree -> RH | false; it leaves an unbounded degree/shift wedge |

The circularity warning in `h787_global_jensen_pf_lemma_audit.md:202-206` is
exact: a Hadamard product formed from already-real Xi zeros proves only
`RH => JP2'`.  It cannot establish `JP2'` independently.

## 7. Strict Artifacts

- `tools/rh_h1322_h803_all_n_strict_assembly.py`
- `research/riemann/h1322_h803_all_n_strict_assembly.json`
- `research/riemann/h1321_h801_alias_repaired_finite_h803.json`
- `research/riemann/h1319_full_xi_compact_interval_certificate.json`
- `research/riemann/h1319_piecewise_h920_conditional_assembly.json`

At materialization time:

```text
H1322 script SHA-256:
30328A7EB4F2EDA26DB5C33D9C224BDF151DFEEAF8B63D13CDF06D7446539AB1

H1322 report SHA-256:
335FFEAAA2114212392D712CC7EC6EFA0AE61257B8249553513C956F10877177
```

The report itself stores live SHA-256, sizes, and content checks for all ten
finite, compact, piecewise, verifier, and note dependencies.
