# Riemann xi coefficients are PF3

This repository accompanies the manuscript
*The Taylor coefficients of the Riemann xi-function form a Pólya frequency
sequence of order 3* by Cyril Rodaro.

The manuscript proves unconditionally that the Taylor coefficients of the
Riemann xi-function form a Pólya frequency sequence of order 3 (PF3). The
computer-assisted parts are supplied as a fail-closed, reproducible artifact.

> This result does **not** prove the Riemann Hypothesis. It does not claim PF4,
> PF-infinity, or hyperbolicity of Jensen polynomials in every degree.

## Repository contents

- [`paper/xi_pf3.pdf`](paper/xi_pf3.pdf): current compiled manuscript.
- [`paper/xi_pf3.tex`](paper/xi_pf3.tex): current LaTeX source.
- [`artifact-v1/`](artifact-v1): browsable, manifest-governed artifact files.
- [`xi-pf3-artifact-v1.zip`](xi-pf3-artifact-v1.zip): canonical frozen archive.
- [`SHA256SUMS`](SHA256SUMS): checksums of the archive, manuscript, and replay
  entry files.
- [`CITATION.cff`](CITATION.cff): citation metadata.

The canonical ZIP remains byte-for-byte identical to the archive deposited in
the first repository commit. The browsable `artifact-v1/` directory contains
all 266 files governed by `MANIFEST.sha256.json`. Non-manifest Python bytecode
caches and the auxiliary Windows-reserved filename `NUL.md` are intentionally
not duplicated in the Git tree; neither is consumed by the replay.

## Quick verification

CPython 3.11 or newer is required. From the repository root:

```bash
python -m venv .venv
# Activate .venv using the command appropriate for your shell.
python -m pip install -r artifact-v1/requirements.txt
cd artifact-v1
python replay.py --verify-only
python replay.py
```

The final command succeeds only when the JSON report contains:

```json
{
  "all_checks_pass": true,
  "classification": "xi_pf3_artifact_replay_pass"
}
```

The replay is fail-closed: a missing file, digest drift, missing dependency,
undecidable interval comparison, failed child replay, or nonzero exit code
causes failure.

## What the replay certifies

The replay runs two final stages:

1. **H1325 -- global xi inputs.** It checks the certified finite range, the
   tail input package, all pinned dependencies, and the order-3 consecutive-row
   conditions used by the manuscript.
2. **H1326 -- Grassmann–Plücker completion.** It replays the exact symbolic
   identities and independent audits that transfer the consecutive-row package
   to all Toeplitz minors of order at most 3.

The manuscript supplies the human-readable mathematical argument. The artifact
supplies exact symbolic checks, certified ball-arithmetic computations,
dependency hashes, certificates, mutation canaries, and independent replay
layers for the computer-assisted steps.

## Reference environment and timing

The release was replayed successfully on 2026-07-12 with:

- Windows build 26200.8655;
- CPython 3.14.0;
- SymPy 1.14.0;
- python-flint 0.8.0 with FLINT 3.3.1;
- mpmath 1.3.0.

On that machine, `--verify-only` took about 4.25 seconds and the complete replay
took about 46 seconds. These are reference measurements, not performance
guarantees. A peak-memory bound has not been certified.

## Integrity

Canonical archive SHA-256:

```text
5D256A2580126EF4872668EBF513B431CD9861D5461786FD75F2A4D4ED9D96C0
```

The archive was originally frozen at commit
[`18aa2da35cd0f45f38041365080d76a8b32ad052`](https://github.com/wulf758/riemann-xi-pf3/commit/18aa2da35cd0f45f38041365080d76a8b32ad052).
The Git blob identifier shown by GitHub is not a file SHA-256 and is not expected
to match this digest.

## Citation

Please use the metadata in [`CITATION.cff`](CITATION.cff). A DOI can be added
after archival of the GitHub release.

## Licenses

- Code and scripts: [MIT](LICENSES/MIT.txt).
- Manuscript, documentation, certificates, and reports:
  [Creative Commons Attribution 4.0 International](LICENSES/CC-BY-4.0.txt).

See [`LICENSE`](LICENSE) for the scope of each license.

## Contact

Cyril Rodaro, independent researcher

<cyril-rodaro@outlook.fr>
