#!/usr/bin/env python3
"""Canonical final H1319 compact certificate generator.

This is the public entry point.  It delegates all mathematical and audit
checks to the verified assembler, then performs the platform-neutral path
normalization required to join Windows manifest paths to the canonical audit
hash table.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from flint import ctx

import rh_h1319_full_xi_compact_interval_manifest as implementation
import rh_h1319_full_xi_compact_interval_verified_certificate as verified


def slash(value: object) -> str:
    return str(value).replace("\\", "/")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--first", default=str(implementation.DEFAULT_FIRST))
    parser.add_argument("--mean-4-5", default=str(implementation.DEFAULT_MEAN_4_5))
    parser.add_argument("--mean-5-10", default=str(implementation.DEFAULT_MEAN_5_10))
    parser.add_argument("--mean-10-22", default=str(implementation.DEFAULT_MEAN_10_22))
    parser.add_argument("--high", default=str(implementation.DEFAULT_HIGH))
    parser.add_argument("--audit", default=str(verified.DEFAULT_AUDIT))
    parser.add_argument("--out", default=str(implementation.DEFAULT_OUT))
    parser.add_argument("--markdown-out", default=str(implementation.DEFAULT_MARKDOWN))
    parser.add_argument("--dps", type=int, default=90)
    args = parser.parse_args()
    ctx.dps = args.dps

    report = verified.build_verified_report(args)
    audit_raw = json.loads(Path(args.audit).read_text(encoding="utf-8"))
    audited_by_path = {
        slash(row.get("path", "")): row
        for row in audit_raw.get("reports", [])
        if isinstance(row, dict)
    }
    mean_segments = [
        segment
        for segment in report["r_proof_chain"]["segments"]
        if segment["kind"] == "score_mean_value_arb"
    ]
    hashes_match = len(mean_segments) == 3
    for segment in mean_segments:
        row = audited_by_path.get(slash(segment["path"]))
        hashes_match = (
            hashes_match
            and bool(row)
            and row.get("sha256") == segment["sha256"]
            and row.get("all_pass") is True
        )

    independent = report["independent_mean_value_audit"]
    independent["checks"]["detailed_report_hashes_match_audit"] = bool(
        hashes_match
    )
    independent["passes"] = all(independent["checks"].values())
    report["checks"]["independent_mean_value_audit_pass"] = independent[
        "passes"
    ]
    report["canonical_generator"] = verified.normalized(Path(__file__))
    report["all_checks_pass"] = all(report["checks"].values())
    report["classification"] = (
        "h1319_full_xi_compact_interval_closed"
        if report["all_checks_pass"]
        else "h1319_full_xi_compact_interval_manifest_failure"
    )

    out = Path(args.out)
    markdown = Path(args.markdown_out)
    out.parent.mkdir(parents=True, exist_ok=True)
    markdown.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    note = implementation.build_markdown(report)
    note += (
        "\n## Independent audit gate\n\n"
        f"Canonical audit: `{independent['path']}`. PASS=`{independent['passes']}`. "
        "The gate pins the engine, H1272, the 360 mean-value boxes, and the "
        "140-dps precision-safe audit runner by SHA-256.\n"
    )
    markdown.write_text(note, encoding="utf-8")
    print(
        json.dumps(
            {
                "schema": report["schema"],
                "classification": report["classification"],
                "all_checks_pass": report["all_checks_pass"],
                "independent_audit_pass": independent["passes"],
                "detailed_hash_join_pass": hashes_match,
                "segments": len(report["r_proof_chain"]["segments"]),
                "out": str(out),
                "markdown": str(markdown),
            },
            indent=2,
        )
    )
    return 0 if report["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
