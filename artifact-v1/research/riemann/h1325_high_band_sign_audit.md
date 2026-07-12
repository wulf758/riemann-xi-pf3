# H1325 High-Band Sign Audit

Classification: `h1325_high_band_discarded_terms_sign_closed`

H927 proves `a>0`, while H943 proves `K_r(w)>=0`.
Consequently `X>=0` and `Y>=0`, and both discarded terms
`-Y` and `-2a^2X^3` are nonpositive.  The H1325 high-band
upper third-cumulant bound is therefore sign-valid.

```text
python tools/rh_h1325_high_band_sign_audit.py --no-write
```
