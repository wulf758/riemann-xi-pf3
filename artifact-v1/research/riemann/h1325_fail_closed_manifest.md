# H1325 Fail-Closed Replay Manifest

Classification: `h1325_xi_global_d3_h811_fail_closed_replayed`

This manifest pins the final proof engines and dependencies, then
replays the primary H1325 assembly, both missing-sign addenda, the
independent audit, the exact reductions, H1322, the full finite D3
certificate, and the independent `n=127` kappa2 seam.

All hash and replay gates pass: `True`.

Certified scope: global translated contiguous D3 for Xi and H811
order-3 Toeplitz minors with consecutive rows.  Sparse-row PF3,
PF-infinity, all-degree Jensen hyperbolicity, and RH are excluded.

```text
python tools/rh_h1325_fail_closed_manifest.py --no-write
```
