#!/usr/bin/env python3
"""Fail-closed manifest and replay for the complete H1326 PF3 package."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "research" / "riemann"
T = ROOT / "tools"

EXPECTED = {
    R / "h1325_xi_global_d3_h811_consecutive_rows.json": "9308258F232A84CB964F0CE4357C2BE279D5E2FE3A49B227B37CCD7E6D43A1A0",
    R / "h1325_h811_dependency_addendum.json": "F034A24C08EED9FC8E0A0405D596DF47A4227CFBEC7E7664D2823062A45A26B4",
    R / "h811_repaired_consecutive_row_theorem.json": "09C0B09673CD582BB1F9968466B5304FBA13A4B8AE0D97B6348EAD60498D6A12",
    R / "h813_sparse_row_symbolic_hierarchy.json": "C92E94BAED6251E729F3B3C02AE5B93A92653118022B32682E7D8B7E24EE5ED6",
    R / "h816_d3_cone_sparse_row_coupling.json": "8E6F0364C259FBBBC1AB66D9B1C250D2F706004FE1842644AFC584D07BDC2A00",
    T / "rh_h1326_plucker_pf3_completion.py": "4F3B58AFD4F15ACD30A3A4B6487E2644DFB199E222EC843CF93EC72CDDFFEC39",
    T / "rh_h1326_strict_pf2_multigap_audit.py": "1F918065DE0293997BDBC30C1C45087DA9D9C306CDED88F8C28E20F7C19BCE3D",
    T / "rh_h1326_pf3_plucker_independent_audit.py": "A26D63DFF24C776BAC7B15FB2296D4E361E9C430DA939A6E176561AD8FC16B85",
    T / "rh_h1326_pf2_arbitrary_gap_independent_addendum.py": "8433E43FDD2161E5D4C846B005395894194862C2C4C1CC6608E35C64467D0222",
    R / "h1326_plucker_pf3_completion.md": "3B8A06A2C190DCD177E1F7F7F4F9B866FCE48C404BAA8E9303BCB0524380EF2C",
    R / "h1326_plucker_pf3_completion.json": "96A9FB69A96D34CAEF61372967DD25773F1CDB418D44C3C824160E09C615FB23",
    R / "h1326_pf3_plucker_independent_audit.json": "8FF12C6E5B20948EBDADE9D8308B4E63757A413BAF9987F182142AF311B6BF5C",
    R / "h1326_pf2_arbitrary_gap_independent_addendum.json": "6D6097D660DA1E7D226D5240DC7875771326455FD0E27F5F5853AB1F5E52FD9A",
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def replay(script: Path, *args: str) -> tuple[bool, dict, str]:
    proc = subprocess.run(
        [sys.executable, str(script), *args],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    try:
        payload = json.loads(proc.stdout)
    except json.JSONDecodeError:
        payload = {}
    return proc.returncode == 0, payload, proc.stderr


def main() -> int:
    checks = {
        f"hash::{path.relative_to(ROOT)}": path.is_file() and digest(path) == expected
        for path, expected in EXPECTED.items()
    }

    main_ok, main_report, main_err = replay(T / "rh_h1326_plucker_pf3_completion.py", "--full-report")
    gap_ok, gap_report, gap_err = replay(T / "rh_h1326_strict_pf2_multigap_audit.py")
    independent_ok, independent_report, independent_err = replay(
        T / "rh_h1326_pf3_plucker_independent_audit.py", "--no-write"
    )
    addendum_ok, addendum_report, addendum_err = replay(
        T / "rh_h1326_pf2_arbitrary_gap_independent_addendum.py", "--no-write"
    )
    checks.update(
        {
            "main_replay_exit_zero": main_ok,
            "main_replay_closed": main_report.get("classification")
            == "h1326_strict_pf2_plus_h811_implies_pf3_closed"
            and main_report.get("all_checks_pass") is True,
            "main_replay_H813_2310": main_report.get("audits", {}).get("H813_active_sparse") == 2310
            and main_report.get("audits", {}).get("H813_failures") == [],
            "main_replay_H816_356": main_report.get("audits", {}).get("H816_full_induction_proved")
            == 356,
            "multigap_replay_exit_zero": gap_ok,
            "multigap_replay_closed": gap_report.get("classification")
            == "h1326_strict_pf2_arbitrary_gap_lemma_audited"
            and gap_report.get("multigap_canaries") == 512,
            "independent_replay_exit_zero": independent_ok,
            "independent_replay_22_of_22": independent_report.get("classification")
            == "h1326_pf3_plucker_independent_audit_pass"
            and independent_report.get("check_count") == 22,
            "independent_partition_27225": independent_report.get("index_partition_audit", {}).get(
                "target_count"
            )
            == 27225
            and independent_report.get("index_partition_audit", {}).get("failure_count") == 0,
            "arbitrary_gap_addendum_exit_zero": addendum_ok,
            "arbitrary_gap_addendum_pass": addendum_report.get("classification")
            == "h1326_arbitrary_gap_strict_pf2_partner_pass"
            and addendum_report.get("symbolic_factorizations_passed") == 150
            and addendum_report.get("triangular_i_le_p_lt_j_cases_passed") == 450,
        }
    )

    static = json.loads((R / "h1326_plucker_pf3_completion.json").read_text(encoding="utf-8"))
    note = (R / "h1326_plucker_pf3_completion.md").read_text(encoding="utf-8")
    checks.update(
        {
            "static_classification_closed": static.get("classification")
            == "h1326_strict_pf2_plus_h811_implies_pf3_closed"
            and static.get("all_checks_pass") is True,
            "note_contains_general_identity": "Delta_(r,u,s) Delta_(r+1,u,d)" in note,
            "note_contains_arbitrary_gap_PF2": "a_(x+alpha)a_(x+beta)" in note,
            "note_states_Xi_PF3": "(gamma_n/n!)_(n>=0) is PF3" in note,
            "note_preserves_scope": "does **not** prove" in note and "Riemann" in note,
        }
    )

    all_checks = all(checks.values())
    report = {
        "schema": "rh_h1326_fail_closed_manifest.v0",
        "classification": (
            "h1326_pf3_package_fail_closed_pass" if all_checks else "h1326_pf3_package_failed"
        ),
        "all_checks_pass": all_checks,
        "check_count": len(checks),
        "failed_checks": [name for name, value in checks.items() if not value],
        "hash_count": len(EXPECTED),
        "replays": {
            "main": main_report.get("classification"),
            "multigap": gap_report.get("classification"),
            "independent": independent_report.get("classification"),
            "arbitrary_gap_addendum": addendum_report.get("classification"),
        },
        "stderr": {
            "main": main_err,
            "multigap": gap_err,
            "independent": independent_err,
            "arbitrary_gap_addendum": addendum_err,
        },
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if all_checks else 1


if __name__ == "__main__":
    raise SystemExit(main())
