# Xi PF3 proof artifact v1.0.2

This frozen supplementary artifact accompanies *The Taylor coefficients of the
Riemann xi-function form a Polya frequency sequence of order 3*.

Permanent archive: https://doi.org/10.5281/zenodo.21360815

## Scope

The replay certifies the global Xi inputs (H1325) and the Grassmann--Plucker
completion to PF3 (H1326). It does **not** claim PF4, PF-infinity, all-degree
Jensen hyperbolicity, or the Riemann Hypothesis.

## Requirements

- CPython 3.11 or newer (frozen on CPython 3.14.0)
- sympy 1.14.0
- python-flint 0.8.0
- mpmath 1.3.0
- numpy 2.3.1

Reference replays succeeded on CPython 3.14.0 under Windows and independently
on CPython 3.11.15 under Ubuntu 24.04.

Install dependencies and run from the artifact root:

```text
python -m pip install -r requirements.txt
python replay.py
```

The replay first verifies every file against `MANIFEST.sha256.json`, then runs
the H1325 and H1326 fail-closed manifests. Any missing file, digest drift,
undecidable interval comparison, failed replay, or nonzero child exit aborts.
Use `python replay.py --verify-only` for a fast integrity/environment check.
