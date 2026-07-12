# H1326 Fail-Closed Manifest

Classification: `h1326_pf3_package_fail_closed_pass`

Canonical replay:

```text
python tools/rh_h1326_fail_closed_manifest.py
```

Result:

```text
29/29 gates passed
13/13 dependency and artifact hashes matched
H813 exact active sparse universe: 2310/2310
H816 residual universe: 356/356
independent index partition: 27225 cases, 0 failures
```

The manifest replays four independent components:

1. `rh_h1326_plucker_pf3_completion.py`;
2. `rh_h1326_strict_pf2_multigap_audit.py`;
3. `rh_h1326_pf3_plucker_independent_audit.py`;
4. `rh_h1326_pf2_arbitrary_gap_independent_addendum.py`.

The canonical theorem and Xi instantiation are in
`h1326_plucker_pf3_completion.md` and
`h1326_plucker_pf3_completion.json`.

`rh_h1326_plucker_first_gap_bridge.py` is an intermediate 50-case draft
created before the general induction was discovered. It is not part of this
manifest and must not be cited as the H1326 result.

The certified conclusion is PF3 for the Xi coefficient sequence. The package
does not claim PF4, PF-infinity, all-degree Jensen hyperbolicity, or RH.
