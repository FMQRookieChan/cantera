from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest

ct = pytest.importorskip("cantera")


REPO_ROOT = Path(__file__).resolve().parents[2]
SAMPLE = REPO_ROOT / "samples" / "python" / "reactors" / "ndodecane_combustor.py"


spec = importlib.util.spec_from_file_location("ndodecane_combustor", SAMPLE)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def test_resolve_mechanism_path_from_filename():
    path = module.resolve_mechanism_path("nDodecane_Reitz.yaml")
    assert path.exists()
    assert path.name == "nDodecane_Reitz.yaml"


def test_compare_single_dual_outputs():
    results = module.compare_single_dual(
        mechanism=module.resolve_mechanism_path("nDodecane_Reitz.yaml"),
        phase="nDodecane_IG",
        engine=module.ENGINE_CASES["V2500"],
        equivalence_ratio=0.35,
        dual_split=0.6,
    )

    for key in ("single", "dual"):
        case = results[key]
        assert case["mdot_fuel_in"] > 0.0
        assert case["fuel_burn_rate"] >= 0.0
        assert case["fuel_burn_rate"] <= case["mdot_fuel_in"]
        assert case["outlet_temperature"] > 300.0
        assert case["outlet_pressure"] > 0.0
