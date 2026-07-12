#!/usr/bin/env python3
"""One-command fail-closed replay for the frozen Xi PF3 artifact."""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve()
ROOT = HERE.parent if (HERE.parent / "tools").is_dir() else HERE.parents[1]
MANIFEST = ROOT / "MANIFEST.sha256.json"

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()

def verify_manifest() -> tuple[dict[str, bool], list[str]]:
    if not MANIFEST.is_file():
        return {"manifest_present": False}, [str(MANIFEST)]
    payload = json.loads(MANIFEST.read_text(encoding="utf-8"))
    checks: dict[str, bool] = {"manifest_present": True}
    failures: list[str] = []
    for record in payload.get("files", []):
        relative = record["path"]
        path = ROOT / relative
        ok = path.is_file() and path.stat().st_size == record["size_bytes"] and sha256(path) == record["sha256"]
        checks[f"file::{relative}"] = ok
        if not ok:
            failures.append(relative)
    checks["manifest_file_count"] = len(payload.get("files", [])) == payload.get("file_count")
    return checks, failures

def dependency_versions() -> tuple[dict[str, str], list[str]]:
    versions: dict[str, str] = {"python": sys.version.split()[0]}
    missing: list[str] = []
    for module_name in ("sympy", "flint", "mpmath"):
        try:
            module = __import__(module_name)
            versions[module_name] = str(getattr(module, "__version__", "unknown"))
        except ImportError:
            missing.append(module_name)
    return versions, missing

def run_stage(name: str, arguments: list[str]) -> dict[str, Any]:
    process = subprocess.run([sys.executable, *arguments], cwd=ROOT, text=True, capture_output=True, check=False)
    try:
        payload = json.loads(process.stdout)
    except json.JSONDecodeError:
        payload = {}
    passed = process.returncode == 0 and payload.get("all_checks_pass") is True and not process.stderr.strip()
    return {
        "name": name,
        "command": [sys.executable, *arguments],
        "returncode": process.returncode,
        "passed": passed,
        "classification": payload.get("classification"),
        "failed_hashes": payload.get("failed_hashes", []),
        "failed_replays": payload.get("failed_replays", []),
        "failed_checks": payload.get("failed_checks", []),
        "stderr": process.stderr.strip(),
        "stdout": process.stdout.strip() if not payload else None,
    }

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true", help="Only verify frozen files and dependencies.")
    args = parser.parse_args()
    manifest_checks, manifest_failures = verify_manifest()
    versions, missing_dependencies = dependency_versions()
    stages: list[dict[str, Any]] = []
    if not args.verify_only and all(manifest_checks.values()) and not missing_dependencies:
        stages = [
            run_stage("H1325 global Xi inputs", ["tools/rh_h1325_fail_closed_manifest.py", "--no-write"]),
            run_stage("H1326 Plucker completion", ["tools/rh_h1326_fail_closed_manifest.py"]),
        ]
    all_checks_pass = (
        all(manifest_checks.values()) and not manifest_failures and not missing_dependencies
        and (args.verify_only or (len(stages) == 2 and all(stage["passed"] for stage in stages)))
    )
    report = {
        "schema": "xi_pf3_artifact_replay.v1",
        "classification": "xi_pf3_artifact_replay_pass" if all_checks_pass else "xi_pf3_artifact_replay_failed",
        "all_checks_pass": all_checks_pass,
        "manifest_failures": manifest_failures,
        "missing_dependencies": missing_dependencies,
        "versions": versions,
        "stage_results": stages,
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if all_checks_pass else 1

if __name__ == "__main__":
    raise SystemExit(main())