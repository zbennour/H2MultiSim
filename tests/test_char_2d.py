"""Characterization test for 2D single-well scenario (golden values)."""

from __future__ import annotations

from h2multisim.models.two_d import UHSAreal2D, scenario_single_well


def test_char_2d_single_well_final_pressure():
    props, wells, sched = scenario_single_well()
    sim = UHSAreal2D(props, wells)
    res = sim.run_schedule(sched, dt_days=1.0)
    # Golden values re-captured 2026-10-31 (v1.3.1: Peaceman r_eq, pore compressibility Ct,
    # real-gas isentropic exponent).  Previous: 79.74419089950966 bar / 268580.8822835042 t.
    assert abs(res.P_avg_bar[-1] - 79.74524583092644) < 1e-9
    assert abs(res.mass_t[-1] - 268575.27263303834) < 1e-6
    assert res.metrics["converged"] == 1.0
