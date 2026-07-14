#!/usr/bin/env python3
"""Build a self-contained, SHA-256-pinned Xi PF3 supplementary artifact."""
from __future__ import annotations
import argparse, hashlib, json, re, shutil, zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSION = "1.0.2"
RELEASE_DATE = "2026-07-14"
DOI = "10.5281/zenodo.21360815"
DEFAULT_TARGET = ROOT / "artifact" / f"xi-pf3-artifact-v{VERSION}"
TOOL_PATTERN = re.compile(r"^(?:rh_h13(?:1[0-9]|2[0-6])_|rh_h(?:80[0-9]|81[0-9]|218|908|920|927|943)|claude_(?:d3|h1319|kappa3)).*\.py$")
RESEARCH_PATTERN = re.compile(r"^(?:h13(?:1[0-9]|2[0-6])_|h(?:80[0-9]|81[0-9]|908|920|927|943)|claude_(?:d3|h1319)).*\.(?:json|md)$")
EXPLICIT_FILES = ["paper/xi_pf3.tex", "replay.py", "tools/build_pf3_artifact.py"]
README = f"""# Xi PF3 proof artifact v{VERSION}

This frozen supplementary artifact accompanies *The Taylor coefficients of the
Riemann xi-function form a Polya frequency sequence of order 3*.

Permanent archive: https://doi.org/{DOI}

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
"""

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()

def selected_files() -> list[Path]:
    files = {ROOT / relative for relative in EXPLICIT_FILES}
    files.update(path for path in (ROOT / "tools").glob("*.py") if TOOL_PATTERN.match(path.name))
    files.update(path for path in (ROOT / "research" / "riemann").iterdir() if path.is_file() and RESEARCH_PATTERN.match(path.name))
    missing = [str(path) for path in files if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing selected artifact files: " + ", ".join(missing))
    return sorted(files)

def write_metadata(target: Path) -> None:
    (target / "README.md").write_bytes(README.encode("utf-8"))
    (target / "requirements.txt").write_bytes(
        "sympy==1.14.0\npython-flint==0.8.0\nmpmath==1.3.0\nnumpy==2.3.1\n".encode("utf-8")
    )
    (target / "ARTIFACT.json").write_bytes((json.dumps({
        "schema": "xi_pf3_artifact.v1", "title": "Xi PF3 proof artifact", "version": VERSION,
        "built_on": RELEASE_DATE, "entrypoint": "python replay.py",
        "paper": "paper/xi_pf3.tex", "permanent_url_or_doi": f"https://doi.org/{DOI}",
    }, indent=2, sort_keys=True) + "\n").encode("utf-8"))

def write_deterministic_archive(target: Path, archive: Path) -> None:
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_STORED) as handle:
        for path in sorted(candidate for candidate in target.rglob("*") if candidate.is_file()):
            name = (Path(target.name) / path.relative_to(target)).as_posix()
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_STORED
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            handle.writestr(info, path.read_bytes(), compress_type=zipfile.ZIP_STORED)

def build(target: Path) -> tuple[Path, Path]:
    if target.exists():
        raise FileExistsError(f"target already exists: {target}")
    target.mkdir(parents=True)
    for source in selected_files():
        relative = source.relative_to(ROOT)
        destination = target / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    write_metadata(target)
    records = []
    for path in sorted(candidate for candidate in target.rglob("*") if candidate.is_file()):
        records.append({"path": path.relative_to(target).as_posix(), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    manifest = {"schema": "xi_pf3_artifact_manifest.v1", "file_count": len(records), "files": records}
    manifest_path = target / "MANIFEST.sha256.json"
    manifest_path.write_bytes((json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    archive = target.with_name(target.name + ".zip")
    write_deterministic_archive(target, archive)
    return manifest_path, archive

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    args = parser.parse_args()
    manifest, archive = build(args.target.resolve())
    print(json.dumps({
        "classification": "xi_pf3_artifact_built", "target": str(args.target.resolve()),
        "manifest": str(manifest), "archive": str(archive),
        "archive_size_bytes": archive.stat().st_size, "archive_sha256": sha256(archive),
    }, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())