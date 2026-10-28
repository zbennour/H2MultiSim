"""Characterization test for 2D single-well scenario (golden values)."""
from __future__ import annotations

from h2multisim.models.two_d import UHSAreal2D, scenario_single_well


def test_char_2d_single_well_final_pressure():
    props, wells, sched = scenario_single_well()
    sim = UHSAreal2D(props, wells)
    res = sim.run_schedule(sched, dt_days=1.0)
    # Golden values captured 2026-10-29 on main @ 5e8fd9e
    assert abs(res.P_avg_bar[-1] - 79.74419089950966) < 1e-9
    assert abs(res.mass_t[-1] - 268580.8822835042) < 1e-6
