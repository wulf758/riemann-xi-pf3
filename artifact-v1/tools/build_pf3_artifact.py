#!/usr/bin/env python3
"""Build a self-contained, SHA-256-pinned Xi PF3 supplementary artifact."""
from __future__ import annotations
import argparse, hashlib, json, re, shutil, zipfile
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TARGET = ROOT / "artifact" / "xi-pf3-artifact-v1"
TOOL_PATTERN = re.compile(r"^(?:rh_h13(?:1[0-9]|2[0-6])_|rh_h(?:80[0-9]|81[0-9]|218|908|920|927|943)|claude_(?:d3|h1319|kappa3)).*\.py$")
RESEARCH_PATTERN = re.compile(r"^(?:h13(?:1[0-9]|2[0-6])_|h(?:80[0-9]|81[0-9]|908|920|927|943)|claude_(?:d3|h1319)).*\.(?:json|md)$")
EXPLICIT_FILES = ["paper/xi_pf3.tex", "tools/rh_pf3_artifact_replay.py", "tools/build_pf3_artifact.py"]
README = """# Xi PF3 proof artifact v1

This frozen supplementary artifact accompanies *The Taylor coefficients of the
Riemann xi-function form a Polya frequency sequence of order 3*.

## Scope

The replay certifies the global Xi inputs (H1325) and the Grassmann--Plucker
completion to PF3 (H1326). It does **not** claim PF4, PF-infinity, all-degree
Jensen hyperbolicity, or the Riemann Hypothesis.

## Requirements

- CPython 3.11 or newer (frozen on CPython 3.14.0)
- sympy 1.14.0
- python-flint 0.8.0
- mpmath 1.3.0

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
    (target / "README.md").write_text(README, encoding="utf-8")
    (target / "requirements.txt").write_text("sympy==1.14.0\npython-flint==0.8.0\nmpmath==1.3.0\n", encoding="utf-8")
    (target / "ARTIFACT.json").write_text(json.dumps({
        "schema": "xi_pf3_artifact.v1", "title": "Xi PF3 proof artifact", "version": "1",
        "built_on": date.today().isoformat(), "entrypoint": "python replay.py",
        "paper": "paper/xi_pf3.tex", "permanent_url_or_doi": None,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")

def build(target: Path) -> tuple[Path, Path]:
    if target.exists():
        raise FileExistsError(f"target already exists: {target}")
    target.mkdir(parents=True)
    for source in selected_files():
        relative = source.relative_to(ROOT)
        destination = target / relative
        if relative.as_posix() == "tools/rh_pf3_artifact_replay.py":
            destination = target / "replay.py"
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    write_metadata(target)
    records = []
    for path in sorted(candidate for candidate in target.rglob("*") if candidate.is_file()):
        records.append({"path": path.relative_to(target).as_posix(), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    manifest = {"schema": "xi_pf3_artifact_manifest.v1", "file_count": len(records), "files": records}
    manifest_path = target / "MANIFEST.sha256.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    archive = target.with_suffix(".zip")
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as handle:
        for path in sorted(candidate for candidate in target.rglob("*") if candidate.is_file()):
            handle.write(path, Path(target.name) / path.relative_to(target))
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