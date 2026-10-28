"""
03_five_spot_sweep.py
=====================
Quickstart example demonstrating 2D multi-well interference in a classic
five-spot pattern versus a heterogeneous diagonal thief-zone channel.
"""

from h2multisim.models.two_d import (
    UHSAreal2D,
    pressure_uniformity,
    scenario_five_spot,
    scenario_thief_zone,
)


def main() -> None:
    nx, ny = 15, 15

    # 1. Homogeneous five-spot
    props_homo, wells_homo, sched_homo = scenario_five_spot(nx, ny)
    sim_homo = UHSAreal2D(props_homo, wells_homo)
    res_homo = sim_homo.run_schedule(sched_homo, dt_days=1.0, snap_days=[1.0])

    # 2. Heterogeneous thief zone
    props_thief, wells_thief, sched_thief = scenario_thief_zone(nx, ny)
    sim_thief = UHSAreal2D(props_thief, wells_thief)
    res_thief = sim_thief.run_schedule(sched_thief, dt_days=1.0, snap_days=[1.0])

    # 3. Compare spatial pressure uniformity
    snap_h = res_homo.snapshots[0]
    p_h = [[snap_h.P_Pa[i * ny + j] for j in range(ny)] for i in range(nx)]
    cv_h = pressure_uniformity(p_h)

    snap_t = res_thief.snapshots[0]
    p_t = [[snap_t.P_Pa[i * ny + j] for j in range(ny)] for i in range(nx)]
    cv_t = pressure_uniformity(p_t)

    print("=" * 60)
    print("H2MultiSim™ 2D Areal Five-Spot vs. Thief Zone:")
    print("=" * 60)
    print(f"Homogeneous Working Gas Cycled : {res_homo.metrics['working_gas_t']:.1f} tonnes")
    print(f"Thief-Zone Working Gas Cycled  : {res_thief.metrics['working_gas_t']:.1f} tonnes")
    print(f"Homogeneous Pressure CV (day 1): {cv_h * 100:.3f} %")
    print(f"Thief-Zone Pressure CV  (day 1): {cv_t * 100:.3f} %")
    print("=" * 60)


if __name__ == "__main__":
    main()
