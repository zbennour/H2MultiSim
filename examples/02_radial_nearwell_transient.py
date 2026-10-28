"""
02_radial_nearwell_transient.py
===============================
Quickstart example demonstrating 1D radial finite-volume near-wellbore
transient drawdown and pressure buildup during shut-in.
"""

from h2multisim.models.one_d import (
    UHSRadial1D,
    okoroafor_radial,
)
from h2multisim.units import pa_to_bar


def main() -> None:
    # 1. Load canonical 1D radial reservoir properties
    props, schedule = okoroafor_radial()

    # 2. Initialize 1D radial simulator (40 logarithmic grid cells)
    sim = UHSRadial1D(props, N=40)

    # 3. Run schedule (cushion fill -> shut-in -> withdrawal -> shut-in)
    result = sim.run_schedule(schedule, dt_days=0.5)
    m = result.metrics

    print("=" * 60)
    print("H2MultiSim™ 1D Radial Well Transient Results:")
    print("=" * 60)
    print(f"Radial Cells             : {sim.N}")
    print(f"Drainage Radius          : {props.r_ext:.1f} m")
    print(f"Final Wellbore Pressure  : {pa_to_bar(m['final_BHP_Pa']):.2f} bar")
    print(f"Average Reservoir P      : {pa_to_bar(m['final_avg_P_Pa']):.2f} bar")
    print(f"Total Gas Cycled         : {m['total_injected_kg'] / 1e3:.2f} tonnes")
    print("=" * 60)


if __name__ == "__main__":
    main()
