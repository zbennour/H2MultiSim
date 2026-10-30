"""
h2multisim.cli
==============
Command-line interface for the H2MultiSim™ simulation toolkit.
Provides commands to run 0D, 1D, 2D, and 3D simulations and execute the validation test suite.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import unittest
from pathlib import Path

from h2multisim import __version__
from h2multisim.config import build_from_config
from h2multisim.models.one_d import (
    UHSRadial1D,
    export_radial_timeseries_csv,
    okoroafor_radial,
    rough_radial,
)
from h2multisim.models.three_d import (
    UHSFieldScale3D,
    export_3d_timeseries_csv,
    scenario_gravity_override,
    scenario_layered_reservoir,
    scenario_single_well_3d,
)
from h2multisim.models.two_d import (
    UHSAreal2D,
    export_2d_timeseries_csv,
    scenario_five_spot,
    scenario_single_well,
    scenario_thief_zone,
)
from h2multisim.models.zero_d import (
    MonteCarloConfig,
    ReservoirState,
    UHSReservoir0D,
    build_okoroafor_reservoir,
    build_rough_reservoir,
    build_teesside_cavern,
    export_timeseries_csv,
    run_monte_carlo,
    schedule_cavern,
    schedule_porous_reservoir,
    schedule_rough,
)
from h2multisim.units import pa_to_bar


def _require_scenario(scale: str, scenario: str, allowed: tuple[str, ...]) -> bool:
    if scenario in allowed:
        return True
    print(
        f"Error: unknown {scale} scenario {scenario!r}. Choose from: {', '.join(allowed)}",
        file=sys.stderr,
    )
    return False


def _run_0d(args: argparse.Namespace) -> int:
    if getattr(args, "_props", None) is not None:
        props = args._props
        sched = args._schedule
        scenario = "config"
    else:
        scenario = args.scenario or "teesside"
        if not _require_scenario("0D", scenario, ("teesside", "rough", "okoroafor")):
            return 2
        if scenario == "teesside":
            props = build_teesside_cavern()
            sched = schedule_cavern(props.reference_pressure)
        elif scenario == "rough":
            props = build_rough_reservoir()
            sched = schedule_rough(props.reference_pressure)
        else:  # okoroafor
            props = build_okoroafor_reservoir()
            sched = schedule_porous_reservoir(props.reference_pressure)

    out_dir = args.out_dir or "output"
    os.makedirs(out_dir, exist_ok=True)

    print("\n=======================================================")
    print(f"  H2MultiSim™ 0D Lumped Simulation: Scenario [{scenario}]")
    print("=======================================================")

    model = UHSReservoir0D(props)
    s0 = ReservoirState(
        pressure=props.reference_pressure,
        temperature=props.reference_temperature,
        mass=model.mass_from_state(props.reference_pressure, props.reference_temperature),
    )

    dt_s = (args.dt_days or 1.0) * 86400.0
    result = model.run_schedule(s0, sched, dt_s=dt_s, collect_time_series=True)
    m = result.metrics

    print(f"  Max Operating Pressure : {pa_to_bar(m['max_pressure_Pa']):.2f} bar")
    print(f"  Min Operating Pressure : {pa_to_bar(m['min_pressure_Pa']):.2f} bar")
    print(f"  Cushion Gas Stored     : {m['cushion_gas_kg'] / 1e3:.2f} tonnes")
    print(f"  Working Gas Capacity   : {m['working_gas_kg'] / 1e3:.2f} tonnes")
    print(f"  Storage Efficiency     : {m['storage_efficiency_pct']:.1f} %")
    print(f"  Total Injected         : {m['total_injected_kg'] / 1e3:.2f} tonnes")
    print(f"  Total Withdrawn        : {m['total_withdrawn_kg'] / 1e3:.2f} tonnes")

    csv_path = os.path.join(out_dir, f"0d_{scenario}_timeseries.csv")
    export_timeseries_csv(result, csv_path)
    print(f"\n  Time series exported -> {os.path.abspath(csv_path)}")

    if args.mc:
        n_mc = args.n_mc or 100
        print(f"\n  Executing Monte Carlo UQ ({n_mc} realizations)...")
        mc_cfg = MonteCarloConfig(
            n_realizations=n_mc,
            porosity_mean=props.porosity,
            porosity_std=props.porosity * 0.1,
            thickness_mean=props.thickness,
            thickness_std=props.thickness * 0.1,
            log10_perm_mean=math.log10(max(props.permeability, 1e-18)),
            log10_perm_std=0.25,
            geothermal_T_mean=props.reference_temperature,
            geothermal_T_std=2.0,
            area_mean=props.area,
            area_std=props.area * 0.1,
        )
        mc_stats = run_monte_carlo(props, sched, mc_cfg, dt_days=args.dt_days or 1.0)
        wg = mc_stats["working_gas_kg"]
        print(
            f"  Working Gas P10: {wg['p10'] / 1e3:.2f} t | P50: {wg['p50'] / 1e3:.2f} t | P90: {wg['p90'] / 1e3:.2f} t"
        )

    print("=======================================================\n")
    return 0


def _run_1d(args: argparse.Namespace) -> int:
    if getattr(args, "_props", None) is not None:
        props = args._props
        sched = args._schedule
        scenario = "config"
    else:
        scenario = args.scenario or "okoroafor"
        if not _require_scenario("1D", scenario, ("okoroafor", "rough")):
            return 2
        if scenario == "rough":
            props, sched = rough_radial()
        else:
            props, sched = okoroafor_radial()

    out_dir = args.out_dir or "output"
    os.makedirs(out_dir, exist_ok=True)
    dt_days = args.dt_days or 0.2

    print("\n=======================================================")
    print(f"  H2MultiSim™ 1D Radial Model: Scenario [{scenario}]")
    print("=======================================================")

    model = UHSRadial1D(props, N=40)
    result = model.run_schedule(sched, dt_days=dt_days)
    m = result.metrics

    print(f"  Grid Cells (logarithmic) : {model.N}")
    print(f"  Initial Pressure         : {pa_to_bar(props.P_init):.2f} bar")
    print(f"  Final Wellbore BHP       : {pa_to_bar(m['final_BHP_Pa']):.2f} bar")
    print(f"  Total Gas Injected       : {m['total_injected_kg'] / 1e3:.2f} tonnes")
    print(f"  Total Gas Withdrawn      : {m['total_withdrawn_kg'] / 1e3:.2f} tonnes")
    print(f"  Net Gas Stored           : {m['net_stored_kg'] / 1e3:.2f} tonnes")

    csv_path = os.path.join(out_dir, f"1d_{scenario}_timeseries.csv")
    export_radial_timeseries_csv(result, csv_path)
    print(f"\n  Time series exported -> {os.path.abspath(csv_path)}")
    print("=======================================================\n")
    return 0


def _run_2d(args: argparse.Namespace) -> int:
    if getattr(args, "_props", None) is not None:
        props = args._props
        wells = args._wells
        sched = args._schedule
        scenario = "config"
    else:
        scenario = args.scenario or "five_spot"
        if not _require_scenario("2D", scenario, ("five_spot", "thief_zone", "single_well")):
            return 2
        nx, ny = 15, 15
        if scenario == "thief_zone":
            props, wells, sched = scenario_thief_zone(nx, ny)
        elif scenario == "single_well":
            props, wells, sched = scenario_single_well(nx, ny)
        else:  # five_spot
            props, wells, sched = scenario_five_spot(nx, ny)

    out_dir = args.out_dir or "output"
    os.makedirs(out_dir, exist_ok=True)
    dt_days = args.dt_days or 0.5

    print("\n=======================================================")
    print(f"  H2MultiSim™ 2D Areal Simulation: Scenario [{scenario}]")
    print("=======================================================")

    model = UHSAreal2D(props, wells)
    result = model.run_schedule(sched, dt_days=dt_days)
    m = result.metrics

    print(
        f"  Domain Grid              : {props.nx} x {props.ny} cells ({props.Lx:.0f}m x {props.Ly:.0f}m)"
    )
    print(f"  Active Wells             : {len(wells)}")
    print(f"  Max Average Pressure     : {m['max_P_bar']:.2f} bar")
    print(f"  Min Average Pressure     : {m['min_P_bar']:.2f} bar")
    print(f"  Working Gas Cycled       : {m['working_gas_t']:.2f} tonnes")
    print(f"  Storage Efficiency       : {m['storage_eff_pct']:.1f} %")

    csv_path = os.path.join(out_dir, f"2d_{scenario}_timeseries.csv")
    export_2d_timeseries_csv(result, csv_path)
    print(f"\n  Time series exported -> {os.path.abspath(csv_path)}")
    print("=======================================================\n")
    return 0


def _run_3d(args: argparse.Namespace) -> int:
    if getattr(args, "_props", None) is not None:
        props = args._props
        wells = args._wells
        sched = args._schedule
        scenario = "config"
    else:
        scenario = args.scenario or "layered_reservoir"
        if not _require_scenario(
            "3D", scenario, ("gravity_override", "single_well_3d", "layered_reservoir")
        ):
            return 2
        nx, ny, nz = 8, 8, 5
        if scenario == "gravity_override":
            props, wells, sched = scenario_gravity_override(nx, ny, nz)
        elif scenario == "single_well_3d":
            props, wells, sched = scenario_single_well_3d(nx, ny, nz)
        else:  # layered_reservoir
            props, wells, sched = scenario_layered_reservoir(nx, ny, nz)

    out_dir = args.out_dir or "output"
    os.makedirs(out_dir, exist_ok=True)
    dt_days = args.dt_days or 2.0

    print("\n=======================================================")
    print(f"  H2MultiSim™ 3D Field-Scale Simulation: [{scenario}]")
    print("=======================================================")

    model = UHSFieldScale3D(props, wells)
    result = model.run_schedule(sched, dt_days=dt_days, verbose=False)
    m = result.metrics

    print(
        f"  Domain Grid              : {props.nx} x {props.ny} x {props.nz} ({props.nx * props.ny * props.nz} total cells)"
    )
    print(f"  Max Reservoir Pressure   : {m['max_P_bar']:.2f} bar")
    print(f"  Min Reservoir Pressure   : {m['min_P_bar']:.2f} bar")
    print(f"  Working Gas Stored       : {m['working_gas_t']:.2f} tonnes")
    print(f"  Total Injected           : {m['total_injected_t']:.2f} tonnes")

    csv_path = os.path.join(out_dir, f"3d_{scenario}_timeseries.csv")
    export_3d_timeseries_csv(result, csv_path)
    print(f"\n  Time series exported -> {os.path.abspath(csv_path)}")
    print("=======================================================\n")
    return 0


def _run_tests(_args: argparse.Namespace) -> int:
    print("\n=======================================================")
    print(f"  Running H2MultiSim™ Test Suite ({__version__})")
    print("=======================================================\n")
    root = Path(__file__).resolve().parents[2]
    tests_dir = root / "tests"
    if not tests_dir.is_dir():
        tests_dir = Path.cwd() / "tests"
    if not (tests_dir / "__init__.py").is_file():
        print(
            "Error: the validation test suite is not installed with this package. "
            "Run `pytest` from a source checkout instead.",
            file=sys.stderr,
        )
        return 2
    loader = unittest.TestLoader()
    suite = loader.discover(
        str(tests_dir), pattern="test_*.py", top_level_dir=str(tests_dir.parent)
    )
    runner = unittest.TextTestRunner(verbosity=2)
    res = runner.run(suite)
    return 0 if res.wasSuccessful() else 1


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="h2multisim",
        description="H2MultiSim™ — Multi-Scale Underground Hydrogen Storage (UHS) Simulator",
    )
    parser.add_argument("--version", action="version", version=f"H2MultiSim {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # run command
    run_parser = subparsers.add_parser(
        "run", help="Run a simulation across 0D, 1D, 2D, or 3D tiers"
    )
    run_parser.add_argument(
        "--scale",
        choices=["0d", "1d", "2d", "3d"],
        required=False,
        help="Dimensional simulation scale (not needed with --config)",
    )
    run_parser.add_argument(
        "--scenario",
        type=str,
        default=None,
        help="Pre-configured scenario name (e.g. teesside, okoroafor, rough, five_spot, thief_zone, gravity_override, layered_reservoir)",
    )
    run_parser.add_argument("--dt-days", type=float, default=None, help="Time step size in days")
    run_parser.add_argument(
        "--mc", action="store_true", help="Enable Monte Carlo uncertainty quantification (0D only)"
    )
    run_parser.add_argument(
        "--n-mc", type=int, default=100, help="Number of Monte Carlo realizations"
    )
    run_parser.add_argument(
        "--out-dir", type=str, default="output", help="Directory for exported CSV time series"
    )
    run_parser.add_argument(
        "--config", type=str, default=None, help="Path to a JSON/YAML run configuration file"
    )

    # test command
    subparsers.add_parser("test", help="Execute the validation benchmark test suite")

    args = parser.parse_args()

    if args.command == "test":
        sys.exit(_run_tests(args))
    elif args.command == "run":
        if getattr(args, "config", None):
            try:
                from h2multisim.config import load_config

                cfg = load_config(args.config)
                scale, props, wells, sched, dt_days, out_dir = build_from_config(cfg)
            except (OSError, ValueError, RuntimeError, KeyError, TypeError) as exc:
                print(f"Error: {exc}", file=sys.stderr)
                sys.exit(2)
            args._props = props
            args._wells = wells
            args._schedule = sched
            args.out_dir = out_dir
            args.dt_days = dt_days
            if scale == "0d":
                sys.exit(_run_0d(args))
            elif scale == "1d":
                sys.exit(_run_1d(args))
            elif scale == "2d":
                sys.exit(_run_2d(args))
            elif scale == "3d":
                sys.exit(_run_3d(args))
        elif not args.scale:
            print("Error: --scale is required unless --config is given", file=sys.stderr)
            sys.exit(2)
        elif args.scale == "0d":
            sys.exit(_run_0d(args))
        elif args.scale == "1d":
            sys.exit(_run_1d(args))
        elif args.scale == "2d":
            sys.exit(_run_2d(args))
        elif args.scale == "3d":
            sys.exit(_run_3d(args))


if __name__ == "__main__":
    main()
