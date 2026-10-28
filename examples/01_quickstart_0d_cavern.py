"""
01_quickstart_0d_cavern.py
==========================
Quickstart example demonstrating 0D lumped simulation of a Teesside-like
salt cavern under high-frequency injection/withdrawal cycling.
"""

from h2multisim import (
    ReservoirState,
    UHSReservoir0D,
    build_teesside_cavern,
    pa_to_bar,
    schedule_cavern,
)


def main() -> None:
    # 1. Load canonical benchmark properties (70,000 m³ cavern, 40–120 bar)
    props = build_teesside_cavern()
    schedule = schedule_cavern(props.reference_pressure)

    # 2. Initialize 0D simulator
    model = UHSReservoir0D(props)
    initial_mass = model.mass_from_state(props.reference_pressure, props.reference_temperature)
    state = ReservoirState(
        pressure=props.reference_pressure,
        temperature=props.reference_temperature,
        mass=initial_mass,
    )

    # 3. Run simulation (cushion fill + first 30 daily cycles)
    run_schedule = schedule[: 1 + 30 * 3]  # cushion fill + 30 cycles
    result = model.run_schedule(state, run_schedule, dt_s=3600.0)  # 1-hour time step
    m = result.metrics

    print("=" * 60)
    print("H2MultiSim™ 0D Salt Cavern Quickstart Results:")
    print("=" * 60)
    print(f"Max Operating Pressure : {pa_to_bar(m['max_pressure_Pa']):.2f} bar")
    print(f"Min Operating Pressure : {pa_to_bar(m['min_pressure_Pa']):.2f} bar")
    print(f"Cushion Gas Retained   : {m['cushion_gas_kg'] / 1e3:.2f} tonnes")
    print(f"Working Gas Capacity   : {m['working_gas_kg'] / 1e3:.2f} tonnes")
    print(f"Cycle Storage Eff.     : {m['storage_efficiency_pct']:.1f} %")
    print("=" * 60)


if __name__ == "__main__":
    main()
