#!/usr/bin/env python3
"""Fail-closed replay manifest for the H1325 Xi translated-D3 package."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
RIEMANN = ROOT / "research" / "riemann"
DEFAULT_OUT = RIEMANN / "h1325_fail_closed_manifest.json"
DEFAULT_MD = RIEMANN / "h1325_fail_closed_manifest.md"

PINNED = {
    "tools/rh_h1325_xi_global_d3_h811_assembly.py": "6B44A88818E521FCA82C19347431EA03C07732C7AC4B1D7361DE8E418961795B",
    "tools/rh_h1325_h811_dependency_addendum.py": "98FDB99B25D092A31E69FA1A6597AFFF0D40C5C877BA29324B54D6895CD70913",
    "tools/rh_h1325_high_band_sign_audit.py": "59FB38405C3D3E0FFE402E81835E960762938322458D31A97E3EB35D6848C978",
    "tools/rh_h1325_xi_global_d3_h811_independent_audit.py": "7AF810D1A0CF489F3DA45955FBE15C1C957F0B6C3F05357A96FA59B3BA325C5B",
    "research/riemann/h1325_xi_global_d3_h811_consecutive_rows.json": "9308258F232A84CB964F0CE4357C2BE279D5E2FE3A49B227B37CCD7E6D43A1A0",
    "research/riemann/h1325_h811_dependency_addendum.json": "F034A24C08EED9FC8E0A0405D596DF47A4227CFBEC7E7664D2823062A45A26B4",
    "research/riemann/h1325_high_band_sign_audit.json": "17705AD7F520FC06003A727EF0A92EAE16435522AB0F9EEA11419ACECFBB354F",
    "research/riemann/claude_d3_alias_repaired_finite_certificate.json": "4AC291EC2FB43D5469A6BD317993645BA56CA68E8F86E7763A24486771E313E1",
    "research/riemann/h1322_h803_all_n_strict_assembly.json": "335FFEAAA2114212392D712CC7EC6EFA0AE61257B8249553513C956F10877177",
    "research/riemann/h811_repaired_consecutive_row_theorem.json": "09C0B09673CD582BB1F9968466B5304FBA13A4B8AE0D97B6348EAD60498D6A12",
    "research/riemann/h908_xi_log_saddle_moment_curvature_reduction.md": "10E3B7F623E7530342EB2B7A9DF0D941517B3417F5755797D5938BD185543099",
    "research/riemann/h943_score_excess_taylor_bridge.md": "D215EBECE963E8823358631CFE1CB020645B4C9293D14E0961415470793562C7",
    "research/riemann/h927_gamma_log_a_le_one_low_threshold.md": "B74EBA4BC923164E6A2A92B4132DDC5A42CEF8818E4B1EC43D7B194B59071E9D",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def run_json(arguments: list[str]) -> tuple[int, dict[str, Any], str]:
    proc = subprocess.run(
        [sys.executable, *arguments],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    data: dict[str, Any] = {}
    if proc.stdout.strip():
        data = json.loads(proc.stdout)
    return proc.returncode, data, proc.stderr.strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--markdown-out", default=str(DEFAULT_MD))
    parser.add_argument("--no-write", action="store_true")
    args = parser.parse_args()

    current_hashes = {name: digest(ROOT / name) for name in PINNED}
    hash_checks = {name: current_hashes[name] == expected for name, expected in PINNED.items()}

    commands = {
        "primary": ["tools/rh_h1325_xi_global_d3_h811_assembly.py", "--no-write"],
        "positivity_addendum": ["tools/rh_h1325_h811_dependency_addendum.py", "--no-write"],
        "high_sign_addendum": ["tools/rh_h1325_high_band_sign_audit.py", "--no-write"],
        "independent_audit": ["tools/rh_h1325_xi_global_d3_h811_independent_audit.py"],
        "exact_kappa2_reduction": ["tools/claude_d3_kappa2_exact_reduction_check.py"],
        "exact_tail_lemma": ["tools/claude_d3_two_sided_tail_lemma_independent_check.py"],
        "q_monotone_report_replay": ["tools/rh_h1322_h803_all_n_strict_assembly.py", "--check-report"],
        "finite_D3_full_replay": ["tools/claude_d3_alias_repaired_finite_certificate_canonical.py", "--out", "NUL"],
        "n127_kappa2_seam_replay": ["tools/claude_d3_kappa2_average_pilot.py", "--out", "NUL", "--target-n", "127", "--prec", "192"],
    }
    replays: dict[str, Any] = {}
    replay_checks: dict[str, bool] = {}
    for name, command in commands.items():
        code, data, stderr = run_json(command)
        replays[name] = {
            "command": [sys.executable, *command],
            "returncode": code,
            "classification": data.get("classification"),
            "status": data.get("status"),
            "all_checks_pass": data.get("all_checks_pass"),
            "stderr": stderr,
        }
        semantic_pass = bool(data.get("all_checks_pass")) or data.get("status") == "ok"
        replay_checks[name] = code == 0 and semantic_pass and not stderr

    all_checks = all(hash_checks.values()) and all(replay_checks.values())
    report = {
        "schema": "rh_h1325_fail_closed_manifest.v0",
        "classification": (
            "h1325_xi_global_d3_h811_fail_closed_replayed"
            if all_checks
            else "h1325_xi_global_d3_h811_fail_closed_replay_failed"
        ),
        "all_checks_pass": all_checks,
        "hash_checks": hash_checks,
        "replay_checks": replay_checks,
        "replays": replays,
        "pinned_sha256": PINNED,
        "current_sha256": current_hashes,
        "certified_conclusion": (
            "For a_n=gamma_n/n!, translated contiguous D3 holds globally and every order-3 Toeplitz minor with consecutive rows is nonnegative."
            if all_checks
            else "No conclusion: at least one pinned dependency or replay failed."
        ),
        "scope": [
            "global translated contiguous D3 for Xi",
            "H811 consecutive-row order-3 minors",
            "not sparse-row PF3",
            "not PF-infinity",
            "not all-degree Jensen hyperbolicity",
            "not the Riemann Hypothesis",
        ],
    }
    markdown = "\n".join(
        [
            "# H1325 Fail-Closed Replay Manifest",
            "",
            f"Classification: `{report['classification']}`",
            "",
            "This manifest pins the final proof engines and dependencies, then",
            "replays the primary H1325 assembly, both missing-sign addenda, the",
            "independent audit, the exact reductions, H1322, the full finite D3",
            "certificate, and the independent `n=127` kappa2 seam.",
            "",
            f"All hash and replay gates pass: `{all_checks}`.",
            "",
            "Certified scope: global translated contiguous D3 for Xi and H811",
            "order-3 Toeplitz minors with consecutive rows.  Sparse-row PF3,",
            "PF-infinity, all-degree Jensen hyperbolicity, and RH are excluded.",
            "",
            "```text",
            "python tools/rh_h1325_fail_closed_manifest.py --no-write",
            "```",
            "",
        ]
    )

    if not args.no_write:
        Path(args.out).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        Path(args.markdown_out).write_text(markdown, encoding="utf-8")
    print(
        json.dumps(
            {
                "classification": report["classification"],
                "all_checks_pass": all_checks,
                "failed_hashes": [name for name, ok in hash_checks.items() if not ok],
                "failed_replays": [name for name, ok in replay_checks.items() if not ok],
                "out": None if args.no_write else str(Path(args.out)),
                "markdown": None if args.no_write else str(Path(args.markdown_out)),
            },
            indent=2,
        )
    )
    return 0 if all_checks else 1


if __name__ == "__main__":
    raise SystemExit(main())
