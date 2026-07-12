#!/usr/bin/env python3
"""Canonical executable wrapper for the H1314 parameter certificate.

The reconstruction and exact coefficient engine live in the sibling
``rh_h1314_h943c_parameter_polynomial_certificate`` module.  Its first
command-line serializer retained a few SymPy Boolean objects.  This canonical
entry point reruns the complete calculation, normalizes symbolic booleans, and
writes the reproducible JSON artifact.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import sympy as sp

import rh_h1314_h943c_parameter_polynomial_certificate as engine


def normalize(value: Any) -> Any:
    if value is sp.true:
        return True
    if value is sp.false:
        return False
    if isinstance(value, dict):
        return {str(key): normalize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [normalize(item) for item in value]
    return value


def certificate() -> dict[str, Any]:
    result = normalize(engine.certificate())
    result["canonical_verifier"] = (
        "tools/rh_h1314_h943c_parameter_polynomial_verifier.py"
    )
    result["engine_module"] = (
        "tools/rh_h1314_h943c_parameter_polynomial_certificate.py"
    )
    result["engine_note"] = (
        "the engine reconstructs all polynomials; its direct CLI serializer "
        "is noncanonical because it does not normalize SymPy booleans"
    )
    return result


def main() -> None:
    result = certificate()
    output = Path(__file__).resolve().parents[1] / "research" / "riemann" / (
        "h1314_h943c_parameter_polynomial_certificate.json"
    )
    output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
