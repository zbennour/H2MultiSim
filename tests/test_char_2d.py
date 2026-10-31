"""Characterization test for 2D single-well scenario (golden values)."""

from __future__ import annotations

import pytest

from h2multisim.models.two_d import UHSAreal2D, scenario_single_well

# Golden values re-captured 2026-10-31 (v1.5.0, PCG pressure solver).  The PCG and
# SOR solvers agree to ~6e-10 bar / 6 g; tolerances sit just above that, i.e. at
# the linear-solver tolerance rather than at round-off.
P_END_BAR = 79.74524583035554
MASS_END_T = 268575.2726275423


@pytest.mark.parametrize("solver", ["pcg", "sor"])
def test_char_2d_single_well_final_pressure(solver: str) -> None:
    props, wells, sched = scenario_single_well()
    sim = UHSAreal2D(props, wells, linear_solver=solver)
    res = sim.run_schedule(sched, dt_days=1.0)
    assert abs(res.P_avg_bar[-1] - P_END_BAR) < 1e-8
    assert abs(res.mass_t[-1] - MASS_END_T) < 1e-4
    assert res.metrics["converged"] == 1.0
