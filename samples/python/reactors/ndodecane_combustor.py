"""
Single and dual fuel inlet n-dodecane combustor comparison
==========================================================

Complete, standalone Python/Cantera combustor comparison using the supplied
``data/nDodecane_Reitz.yaml`` detailed mechanism.

This script compares:

- single fuel inlet
- dual fuel inlet

for representative combustor-inlet presets:

- ``V2500``
- ``PW1100G``

Supported mechanism phases:

- ``nDodecane_IG`` (ideal gas)
- ``nDodecane_RK`` (Redlich-Kwong)

Saved outputs:

- numerical summary (JSON + CSV)
- comparison plot (PNG)

The default engine values are representative assumptions for demonstration and
are **not** official manufacturer design-point data.

Examples
--------

Linux/macOS:

.. code-block:: bash

    python samples/python/reactors/ndodecane_combustor.py \
        --engine V2500 \
        --phase nDodecane_IG \
        --output-dir outputs/v2500

Windows CMD:

.. code-block:: bat

    python samples\\python\\reactors\\ndodecane_combustor.py --engine V2500 --phase nDodecane_IG --output-dir outputs\\v2500
    python samples\\python\\reactors\\ndodecane_combustor.py --engine PW1100G --phase nDodecane_RK --no-show
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import cantera as ct


@dataclass(frozen=True)
class EngineCase:
    name: str
    Tin: float
    Pin: float
    mdot_air: float


@dataclass(frozen=True)
class SimulationSettings:
    advance_time: float
    steady_state: bool


ENGINE_CASES: Dict[str, EngineCase] = {
    "V2500": EngineCase(
        name="IAE V2500 representative combustor inlet",
        Tin=700.0,
        Pin=14.5 * ct.one_atm,
        mdot_air=45.0,
    ),
    "PW1100G": EngineCase(
        name="PW1100G representative combustor inlet",
        Tin=730.0,
        Pin=16.0 * ct.one_atm,
        mdot_air=50.0,
    ),
}


def resolve_mechanism_path(mechanism: str | Path) -> Path:
    """Resolve mechanism path from CWD, sample directory, or repository root."""
    path = Path(mechanism)
    if path.is_absolute() and path.exists():
        return path

    search_paths = [
        Path.cwd() / path,
        Path(__file__).resolve().parent / path,
        Path(__file__).resolve().parents[3] / path,
        Path(__file__).resolve().parents[3] / "data" / path,
    ]
    for candidate in search_paths:
        if candidate.exists():
            return candidate.resolve()

    raise FileNotFoundError(
        f"Mechanism file '{mechanism}' not found from current directory, "
        "sample directory, or repository data directory."
    )


def make_gas(mechanism: Path, phase: str) -> ct.Solution:
    gas = ct.Solution(str(mechanism), phase)
    gas.case_sensitive_species_names = True
    return gas


def make_combustor_reactor(gas: ct.Solution):
    if gas.thermo_model == "ideal-gas":
        return ct.IdealGasReactor(gas)
    return ct.Reactor(gas)


def _find_species_name(gas: ct.Solution, name: str) -> str | None:
    species = {sp.lower(): sp for sp in gas.species_names}
    return species.get(name.lower())


def _get_species_ppm(gas: ct.Solution, name: str) -> float | None:
    species_name = _find_species_name(gas, name)
    if species_name is None:
        return None
    return gas[species_name].X[0] * 1e6


def setup_combustor_case(
    mechanism: Path,
    phase: str,
    engine: EngineCase,
    equivalence_ratio: float,
    fuel_split: tuple[float, ...],
    settings: SimulationSettings,
):
    if abs(sum(fuel_split) - 1.0) > 1e-12:
        raise ValueError("fuel_split must sum to 1")

    gas = make_gas(mechanism, phase)

    afr_st = gas.stoich_air_fuel_ratio("c12h26:1", "o2:1, n2:3.76", basis="mass")
    mdot_fuel_total = engine.mdot_air * equivalence_ratio / afr_st

    gas.TPX = 300.0, engine.Pin, "c12h26:1"
    fuel_inlet = ct.Reservoir(gas)

    gas.TPX = engine.Tin, engine.Pin, "o2:1, n2:3.76"
    oxidizer_inlet = ct.Reservoir(gas)

    gas.TPX = 300.0, engine.Pin, "h:1"
    igniter_inlet = ct.Reservoir(gas)

    gas.TPX = engine.Tin, engine.Pin, "n2:1"
    combustor = make_combustor_reactor(gas)
    combustor.volume = 0.03

    exhaust = ct.Reservoir(gas)

    for split in fuel_split:
        ct.MassFlowController(fuel_inlet, combustor, mdot=split * mdot_fuel_total)
    ct.MassFlowController(oxidizer_inlet, combustor, mdot=engine.mdot_air)

    ct.MassFlowController(
        igniter_inlet,
        combustor,
        mdot=ct.Func1("Gaussian", [5.0e-4, 0.02, 0.02]),
    )

    outlet = ct.Valve(combustor, exhaust, K=2.0)

    sim = ct.ReactorNet([combustor])
    sim.advance(settings.advance_time)
    if settings.steady_state:
        sim.solve_steady()

    fuel_mass_fraction = combustor.phase["c12h26"].Y[0]
    fuel_burn_rate = max(0.0, mdot_fuel_total - outlet.mass_flow_rate * fuel_mass_fraction)
    fuel_consumed = fuel_burn_rate * settings.advance_time

    return {
        "mdot_fuel_in": mdot_fuel_total,
        "fuel_burn_rate": fuel_burn_rate,
        "fuel_consumed": fuel_consumed,
        "outlet_temperature": combustor.T,
        "outlet_pressure": combustor.phase.P,
        "outlet_NO_ppm": _get_species_ppm(combustor.phase, "no"),
    }


def compare_single_dual(
    mechanism: Path,
    phase: str,
    engine: EngineCase,
    equivalence_ratio: float,
    dual_split: float,
    settings: SimulationSettings,
):
    single = setup_combustor_case(
        mechanism,
        phase,
        engine,
        equivalence_ratio,
        fuel_split=(1.0,),
        settings=settings,
    )

    dual = setup_combustor_case(
        mechanism,
        phase,
        engine,
        equivalence_ratio,
        fuel_split=(dual_split, 1.0 - dual_split),
        settings=settings,
    )

    return {"single": single, "dual": dual}


def format_no_value(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.1f}"


def print_results(engine: EngineCase, phase: str, mechanism: Path, results: dict) -> None:
    print(f"\nEngine case: {engine.name}")
    print("Representative inputs only; not official manufacturer data.")
    print(f"Mechanism: {mechanism}")
    print(f"Phase: {phase}")
    print("=" * 108)
    print(
        f"{'Configuration':<14}{'Fuel in [kg/s]':>16}{'Fuel burned [kg/s]':>20}"
        f"{'Fuel consumed [kg]':>20}{'T_out [K]':>12}{'P_out [bar]':>12}{'NO [ppm]':>14}"
    )
    for label in ("single", "dual"):
        case = results[label]
        no_text = format_no_value(case["outlet_NO_ppm"])
        print(
            f"{label:<14}{case['mdot_fuel_in']:>16.5f}{case['fuel_burn_rate']:>20.5f}"
            f"{case['fuel_consumed']:>20.5f}{case['outlet_temperature']:>12.1f}"
            f"{case['outlet_pressure'] / 1e5:>12.3f}{no_text:>14}"
        )


def save_results(
    engine_key: str,
    engine: EngineCase,
    phase: str,
    mechanism: Path,
    settings: SimulationSettings,
    results: dict,
    output_dir: Path,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_path = output_dir / "ndodecane_combustor_results.json"
    csv_path = output_dir / "ndodecane_combustor_results.csv"

    serializable = {
        "engine": engine_key,
        "engine_name": engine.name,
        "phase": phase,
        "mechanism": str(mechanism),
        "advance_time": settings.advance_time,
        "steady_state": settings.steady_state,
        "results": results,
    }
    summary_path.write_text(json.dumps(serializable, indent=2), encoding="utf-8")

    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "configuration",
                "fuel_in_kg_s",
                "fuel_burn_rate_kg_s",
                "fuel_consumed_kg",
                "outlet_temperature_K",
                "outlet_pressure_bar",
                "outlet_no_ppm",
            ],
        )
        writer.writeheader()
        for config in ("single", "dual"):
            case = results[config]
            writer.writerow(
                {
                    "configuration": config,
                    "fuel_in_kg_s": f"{case['mdot_fuel_in']:.8g}",
                    "fuel_burn_rate_kg_s": f"{case['fuel_burn_rate']:.8g}",
                    "fuel_consumed_kg": f"{case['fuel_consumed']:.8g}",
                    "outlet_temperature_K": f"{case['outlet_temperature']:.8g}",
                    "outlet_pressure_bar": f"{case['outlet_pressure'] / 1e5:.8g}",
                    "outlet_no_ppm": "" if case["outlet_NO_ppm"] is None else f"{case['outlet_NO_ppm']:.8g}",
                }
            )

    return {"summary": summary_path, "csv": csv_path}


def plot_results(results: dict, output_dir: Path, show: bool) -> Path:
    if not show:
        import matplotlib

        matplotlib.use("Agg", force=True)

    import matplotlib.pyplot as plt

    labels = ["single", "dual"]
    fuel_in = [results[k]["mdot_fuel_in"] for k in labels]
    fuel_burn = [results[k]["fuel_burn_rate"] for k in labels]
    fuel_consumed = [results[k]["fuel_consumed"] for k in labels]
    tout = [results[k]["outlet_temperature"] for k in labels]
    pout = [results[k]["outlet_pressure"] / 1e5 for k in labels]

    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    axs = axs.ravel()

    axs[0].bar(labels, fuel_in, label="fuel in")
    axs[0].bar(labels, fuel_burn, alpha=0.8, label="fuel burned")
    axs[0].set_ylabel("kg/s")
    axs[0].set_title("Fuel flow / burn rate")
    axs[0].legend(loc="best")

    axs[1].bar(labels, fuel_consumed)
    axs[1].set_ylabel("kg")
    axs[1].set_title("Fuel consumed over advance interval")

    axs[2].bar(labels, tout)
    axs[2].set_ylabel("K")
    axs[2].set_title("Outlet temperature")

    axs[3].bar(labels, pout)
    axs[3].set_ylabel("bar")
    axs[3].set_title("Outlet pressure")

    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = output_dir / "ndodecane_combustor_comparison.png"
    fig.savefig(plot_path, dpi=200, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return plot_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare single- and dual-inlet n-dodecane combustor models "
        "using Cantera detailed chemistry.",
        epilog=(
            "Windows CMD example: python samples\\python\\reactors\\ndodecane_combustor.py "
            "--engine V2500 --phase nDodecane_IG --output-dir outputs\\v2500 --no-show"
        ),
    )
    parser.add_argument("--engine", choices=ENGINE_CASES, default="V2500")
    parser.add_argument("--mechanism", default="nDodecane_Reitz.yaml")
    parser.add_argument("--phase", choices=("nDodecane_IG", "nDodecane_RK"), default="nDodecane_IG")
    parser.add_argument("--equivalence-ratio", type=float, default=0.38)
    parser.add_argument(
        "--dual-split",
        type=float,
        default=0.5,
        help="Primary fuel split for dual inlet (0-1).",
    )
    parser.add_argument(
        "--advance-time",
        type=float,
        default=0.08,
        help="Network advance time in seconds before optional steady-state solve.",
    )
    parser.add_argument(
        "--steady-state",
        dest="steady_state",
        action="store_true",
        default=True,
        help="Solve for reactor steady state after advancing (default: enabled).",
    )
    parser.add_argument(
        "--no-steady-state",
        dest="steady_state",
        action="store_false",
        help="Skip steady-state solve and report advanced transient state.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("ndodecane_combustor_output"),
        help="Directory where JSON/CSV/PNG outputs are saved.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not show interactive plots (uses non-interactive backend).",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if not 0.0 < args.dual_split < 1.0:
        raise ValueError("--dual-split must be between 0 and 1")
    if args.equivalence_ratio <= 0.0:
        raise ValueError("--equivalence-ratio must be positive")
    if args.advance_time <= 0.0:
        raise ValueError("--advance-time must be positive")


def main() -> None:
    args = parse_args()
    validate_args(args)

    mechanism = resolve_mechanism_path(args.mechanism)
    engine = ENGINE_CASES[args.engine]
    settings = SimulationSettings(advance_time=args.advance_time, steady_state=args.steady_state)

    results = compare_single_dual(
        mechanism=mechanism,
        phase=args.phase,
        engine=engine,
        equivalence_ratio=args.equivalence_ratio,
        dual_split=args.dual_split,
        settings=settings,
    )

    print_results(engine, args.phase, mechanism, results)
    output_paths = save_results(
        engine_key=args.engine,
        engine=engine,
        phase=args.phase,
        mechanism=mechanism,
        settings=settings,
        results=results,
        output_dir=args.output_dir,
    )
    plot_path = plot_results(results, output_dir=args.output_dir, show=not args.no_show)

    print("\nSaved outputs:")
    print(f"- {output_paths['summary']}")
    print(f"- {output_paths['csv']}")
    print(f"- {plot_path}")


if __name__ == "__main__":
    main()
