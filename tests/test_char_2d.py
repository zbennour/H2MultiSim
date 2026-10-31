"""Characterization test for 2D single-well scenario (golden values)."""

from __future__ import annotations

import pytest

from h2multisim.models.two_d import UHSAreal2D, scenario_single_well

# Golden values re-captured 2026-10-31 (v1.10.0, rock heat capacity in local thermal
# equilibrium).  The PCG and SOR solvers agree to ~2e-9 bar / 8 g; tolerances sit just
# above that, i.e. at the linear-solver tolerance rather than at round-off.
P_END_BAR = 79.95060828593202
MASS_END_T = 268829.093510512


@pytest.mark.parametrize("solver", ["pcg", "sor"])
def test_char_2d_single_well_final_pressure(solver: str) -> None:
    props, wells, sched = scenario_single_well()
    sim = UHSAreal2D(props, wells, linear_solver=solver)
    res = sim.run_schedule(sched, dt_days=1.0)
    assert abs(res.P_avg_bar[-1] - P_END_BAR) < 1e-8
    assert abs(res.mass_t[-1] - MASS_END_T) < 1e-4
    assert res.metrics["converged"] == 1.0
