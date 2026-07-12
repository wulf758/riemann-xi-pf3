# H1325 Xi Global Translated D3 And H811 Consecutive Rows

Classification: `xi_global_translated_d3_h811_consecutive_rows_closed`

## Result

For the Xi coefficients `a_n=gamma_n/n!`, translated contiguous
order-3 Toeplitz determinants are strictly positive at every start.
Together with the previously certified PF2, global q-monotonicity,
and `q_2<1/2`, H811 therefore proves every order-3 Toeplitz minor
with consecutive rows is nonnegative.

This is not full PF3: sparse-row minors are outside H811.

## New two-sided cumulant package

For `kappa(t)=log M_t` and every `t>=125`, the assembled bounds are

```text
kappa''(t)  < 4/(5t),
kappa'''(t) < 11/(20t^2).
```

The compact proof rereads `552` rigorous Arb boxes.
Its worst second-moment upper is `1.237247964367270469665527343750000000000000000000000000000000000000000000000000000000000000000000000`.
Its worst `t^2*kappa'''` upper is `[0.1248013790155838453738382421844543815429962235910341892273014985370593876723432913422584533691406250 +/- 2.70e-103]`.

From `r=22` onward the exact endpoint coefficients are

- `t*kappa''` upper: `46606678326443294424723416932878/492590289071124220328797309123325`;
- `t^2*kappa'''` upper: `168487269652402462458669370436858871/307868930669452637705498318202078125`;
- third-cumulant slack below `11/20`: `3362568863185953117418818297136391/1231475722677810550821993272808312500`.

Both analytic envelopes decrease with `r`.

## Exact tail lemma

For `n>=127`, the derivative bounds give rational upper envelopes
`q_n<=X_n` and `q_(n+1)/q_n<=S_n`.  With `Y_n=S_n X_n`,
coefficient-positive polynomials after `n=127+m` prove

```text
X_n<1,  Y_n<1,
1-Y_n-Y_n^2(1-X_n)>0.
```

Using global q-monotonicity to replace `q_(n+2)` by `q_(n+1)`
then proves `D_n>0` for every `n>=127`.

## Finite seam and H811

The alias-repaired coefficient balls certify `D_n>0` independently
from both quotient and Toeplitz-determinant formulas for `2<=n<=126`.
Start `1` follows from `q_2<1/2` and `q_3>0`; start `0` is triangular.
The resulting conclusion is: `every order-3 Toeplitz minor with consecutive rows is nonnegative for a_n=gamma_n/n!`.

## Reproduction

```text
python tools/rh_h1325_xi_global_d3_h811_assembly.py --no-write
```

The checker parses directed Arb endpoints, verifies the H1319 hash
chain, reruns H1317 and the independent exact tail checker, and
regenerates every shifted polynomial.

## Scope

- Proved here: global translated contiguous D3 for Xi and H811's
  consecutive-row order-3 conclusion.
- Not proved: sparse-row PF3, PF-infinity, all-degree Jensen
  hyperbolicity, or the Riemann Hypothesis.
