"""
05_cushion_gas_mixing_2d.py
===========================
H2 storage in a depleted gas field with a native CH4 cushion (2D areal model).

Shows how the produced H2 purity and the withdrawn energy are diluted by the
cushion gas, and how much H2 stays behind mixed into the cushion.
"""

from h2multisim import UHSAreal2D, pa_to_bar
from h2multisim.models.two_d import scenario_single_well


def main() -> None:
    props, wells, schedule = scenario_single_well(nx=9, ny=9)
    props.cushion_gas = "CH4"  # native gas of the depleted field
    props.initial_h2_mole_fraction = 0.0
    props.dispersivity_m = 2.0  # m, longitudinal dispersivity

    sim = UHSAreal2D(props, wells)
    result = sim.run_schedule(schedule, dt_days=2.0, snap_days=[30.0])
    m = result.metrics

    print("H2 injection into a CH4-filled depleted field (2D, 9x9)")
    print(f"  converged               : {bool(m['converged'])}")
    print(f"  injected H2             : {m['total_injected_t']:.0f} t")
    print(f"  withdrawn gas           : {m['total_withdrawn_t']:.0f} t")
    print(f"  withdrawn energy (LHV)  : {m['total_withdrawn_energy_GWh_LHV']:.2f} GWh")
    print(f"  min produced H2 purity  : {m['min_withdrawn_h2_purity_mol_pct']:.1f} mol%")
    print(f"  mean H2 left in place   : {m['mean_final_h2_mol_pct']:.2f} mol%")

    purity = result.per_well_h2_mole_fraction["W1"]
    bhp = result.per_well_BHP["W1"]
    print("\n  day   BHP [bar]   well H2 [mol%]")
    for k in range(0, len(result.times_d), 5):
        print(f"  {result.times_d[k]:4.0f}   {bhp[k]:8.2f}   {100 * purity[k]:8.1f}")

    snap = result.snapshots[0]
    ny = props.ny
    row = [100 * snap.h2_mole_fraction[(props.nx // 2) * ny + j] for j in range(ny)]
    print(f"\n  H2 [mol%] across the well row at day {snap.time_d:.0f}:")
    print("  " + " ".join(f"{v:5.1f}" for v in row))
    print(
        f"  (average reservoir pressure at the end: {result.P_avg_bar[-1]:.1f} bar;"
        f" initial {pa_to_bar(props.P_init):.1f} bar)"
    )


if __name__ == "__main__":
    main()
