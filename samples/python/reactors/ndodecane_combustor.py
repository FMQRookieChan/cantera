"""
Single and dual fuel inlet n-dodecane combustor comparison
==========================================================

Compare single-fuel-inlet and dual-fuel-inlet combustor configurations using
Cantera's detailed ``nDodecane_Reitz.yaml`` mechanism.

Reported outputs:

- fuel mass-flow and burned-fuel rate
- outlet temperature
- outlet pressure

The mechanism already contains nitrogen chemistry, so NO predictions come from
its native detailed kinetics. A separate extended-Zeldovich source term is not
added in this example.

Requires: cantera >= 3.2, matplotlib >= 2.0

.. tags:: Python, combustion, reactor network, kinetics, plotting
"""

from __future__ import annotations

import argparse
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
        f"sample directory, or repository data directory."
    )


def make_gas(mechanism: Path, phase: str) -> ct.Solution:
    gas = ct.Solution(str(mechanism), phase)
    gas.case_sensitive_species_names = True
    return gas


def make_combustor_reactor(gas: ct.Solution):
    if gas.thermo_model == "ideal-gas":
        return ct.IdealGasReactor(gas)
    return ct.Reactor(gas)


def setup_combustor_case(
    mechanism: Path,
    phase: str,
    engine: EngineCase,
    equivalence_ratio: float,
    fuel_split: tuple[float, ...],
):
    if abs(sum(fuel_split) - 1.0) > 1e-12:
        raise ValueError("fuel_split must sum to 1")

    gas = make_gas(mechanism, phase)

    # Set target overall fuel-air ratio using mechanism stoichiometry.
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

    controllers = []
    for split in fuel_split:
        controllers.append(
            ct.MassFlowController(fuel_inlet, combustor, mdot=split * mdot_fuel_total)
        )
    controllers.append(ct.MassFlowController(oxidizer_inlet, combustor, mdot=engine.mdot_air))

    igniter = ct.MassFlowController(
        igniter_inlet,
        combustor,
        mdot=ct.Func1("Gaussian", [5.0e-4, 0.02, 0.02]),
    )

    outlet = ct.Valve(combustor, exhaust, K=2.0)

    sim = ct.ReactorNet([combustor])
    sim.advance(0.08)
    sim.solve_steady()

    fuel_mass_fraction = combustor.thermo["c12h26"].Y[0]
    fuel_burn_rate = max(0.0, mdot_fuel_total - outlet.mass_flow_rate * fuel_mass_fraction)

    return {
        "mdot_fuel_in": mdot_fuel_total,
        "fuel_burn_rate": fuel_burn_rate,
        "fuel_consumption": fuel_burn_rate,
        "outlet_temperature": combustor.T,
        "outlet_pressure": combustor.thermo.P,
        "outlet_NO_ppm": combustor.thermo["no"].X[0] * 1e6,
        "controllers": controllers,
        "igniter": igniter,
    }


def compare_single_dual(
    mechanism: Path,
    phase: str,
    engine: EngineCase,
    equivalence_ratio: float,
    dual_split: float,
):
    single = setup_combustor_case(
        mechanism,
        phase,
        engine,
        equivalence_ratio,
        fuel_split=(1.0,),
    )

    dual = setup_combustor_case(
        mechanism,
        phase,
        engine,
        equivalence_ratio,
        fuel_split=(dual_split, 1.0 - dual_split),
    )

    return {"single": single, "dual": dual}


def print_results(engine: EngineCase, phase: str, mechanism: Path, results: dict) -> None:
    print(f"\nEngine case: {engine.name}")
    print(f"Mechanism: {mechanism}")
    print(f"Phase: {phase}")
    print("=" * 78)
    print(f"{'Configuration':<22}{'Fuel in [kg/s]':>16}{'Fuel burned [kg/s]':>20}"
          f"{'T_out [K]':>12}{'P_out [bar]':>12}")
    for label in ("single", "dual"):
        case = results[label]
        print(
            f"{label:<22}{case['mdot_fuel_in']:>16.5f}{case['fuel_burn_rate']:>20.5f}"
            f"{case['outlet_temperature']:>12.1f}{case['outlet_pressure'] / 1e5:>12.3f}"
        )


def plot_results(results: dict) -> None:
    import matplotlib.pyplot as plt

    labels = ["single", "dual"]
    fuel_in = [results[k]["mdot_fuel_in"] for k in labels]
    fuel_burn = [results[k]["fuel_burn_rate"] for k in labels]
    tout = [results[k]["outlet_temperature"] for k in labels]
    pout = [results[k]["outlet_pressure"] / 1e5 for k in labels]

    fig, axs = plt.subplots(1, 3, figsize=(11, 4))
    axs[0].bar(labels, fuel_in, label="fuel in")
    axs[0].bar(labels, fuel_burn, alpha=0.8, label="fuel burned")
    axs[0].set_ylabel("kg/s")
    axs[0].set_title("Fuel flow / burn rate")
    axs[0].legend(loc="best")

    axs[1].bar(labels, tout)
    axs[1].set_ylabel("K")
    axs[1].set_title("Outlet temperature")

    axs[2].bar(labels, pout)
    axs[2].set_ylabel("bar")
    axs[2].set_title("Outlet pressure")

    fig.tight_layout()
    plt.show()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", choices=ENGINE_CASES, default="V2500")
    parser.add_argument("--mechanism", default="nDodecane_Reitz.yaml")
    parser.add_argument("--phase", default="nDodecane_IG")
    parser.add_argument("--equivalence-ratio", type=float, default=0.38)
    parser.add_argument(
        "--dual-split",
        type=float,
        default=0.5,
        help="Primary fuel split for dual inlet (0-1).",
    )
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not 0.0 < args.dual_split < 1.0:
        raise ValueError("--dual-split must be between 0 and 1")

    mechanism = resolve_mechanism_path(args.mechanism)
    engine = ENGINE_CASES[args.engine]

    results = compare_single_dual(
        mechanism=mechanism,
        phase=args.phase,
        engine=engine,
        equivalence_ratio=args.equivalence_ratio,
        dual_split=args.dual_split,
    )

    print_results(engine, args.phase, mechanism, results)
    if not args.no_plot:
        plot_results(results)


if __name__ == "__main__":
    main()
