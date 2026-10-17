from __future__ import annotations

import importlib.util
import json
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


def test_resolve_mechanism_path_from_filename(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
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
        settings=module.SimulationSettings(advance_time=0.08, steady_state=True),
    )

    for key in ("single", "dual"):
        case = results[key]
        assert case["mdot_fuel_in"] > 0.0
        assert case["fuel_burn_rate"] >= 0.0
        assert case["fuel_burn_rate"] <= case["mdot_fuel_in"]
        assert case["fuel_consumed"] >= 0.0
        assert case["outlet_temperature"] > 300.0
        assert case["outlet_pressure"] > 0.0


def test_save_outputs_and_headless_plot(tmp_path):
    matplotlib = pytest.importorskip("matplotlib")

    results = module.compare_single_dual(
        mechanism=module.resolve_mechanism_path("nDodecane_Reitz.yaml"),
        phase="nDodecane_IG",
        engine=module.ENGINE_CASES["V2500"],
        equivalence_ratio=0.35,
        dual_split=0.6,
        settings=module.SimulationSettings(advance_time=0.05, steady_state=False),
    )

    settings = module.SimulationSettings(advance_time=0.05, steady_state=False)
    paths = module.save_results(
        engine_key="V2500",
        engine=module.ENGINE_CASES["V2500"],
        phase="nDodecane_IG",
        mechanism=module.resolve_mechanism_path("nDodecane_Reitz.yaml"),
        settings=settings,
        results=results,
        output_dir=tmp_path,
    )

    plot_path = module.plot_results(results, output_dir=tmp_path, show=False)

    assert paths["summary"].exists()
    assert paths["csv"].exists()
    assert plot_path.exists()

    summary = json.loads(paths["summary"].read_text(encoding="utf-8"))
    assert summary["engine"] == "V2500"
    assert "single" in summary["results"]
    assert "dual" in summary["results"]

    csv_text = paths["csv"].read_text(encoding="utf-8")
    assert "configuration" in csv_text
    assert "fuel_consumed_kg" in csv_text

    assert matplotlib.get_backend().lower() == "agg"
