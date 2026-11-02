# Assumptions & Limitations

This page states what H2MultiSim models, what it does not, and the measured accuracy envelope, so results can be used within their validity range. Quantitative figures come from the reproducible scripts in `audit/` and the tests in `tests/` (see [BENCHMARKS.md](BENCHMARKS.md)); `h2multisim verify` re-checks the thermodynamic envelope below against embedded CoolProp references.

## Scope by tier

| | 0D lumped | 1D radial | 2D areal | 3D field |
|---|---|---|---|---|
| Flow | Material balance + Darcy well productivity | Implicit FV, log grid | Implicit FV, Cartesian | Implicit FV, Cartesian, gravity |
| Gas composition | H2/CH4/CO2/N2 well-mixed, optional plume segregation | Binary H2 + one cushion gas, advection–dispersion | Binary H2 + one cushion gas, advection–dispersion | Binary H2 + one cushion gas, advection–dispersion with gravity |
| Thermal | Real-gas energy balance (injected enthalpy, wall heat transfer; rock + connate-water heat capacity in porous media) | Per-cell isentropic update shared with rock/water heat capacity + relaxation | Same as 1D | Same as 1D |
| Pressure solver | Safeguarded Newton on the mass-balance root | Thomas (tridiagonal) | PCG (or SOR) | PCG (or SOR) |

## Physics not modelled

- **Two-phase flow.** Water is immobile at connate saturation `S_wr`. Its effect enters as a gas-accessible pore fraction and an endpoint gas relative permeability; saturations are not transported. `RelativePermeabilityModel` (Corey / Brooks-Corey, Buckley-Leverett) is available for analysis but the solvers use only its endpoint values. There is no capillary pressure in the solvers, no hysteresis, no residual trapping during cycling, and no H2 dissolution in brine.
- **Geochemistry and microbiology.** No mineral reactions, no microbial H2 consumption (methanogenesis, sulphate reduction).
- **Geomechanics.** Pore volume follows `exp(C_t·(P − P_ref))` (rigid for caverns); no stress coupling, fault reactivation or salt creep.
- **Heat transport in grid tiers.** 1D/2D/3D apply an isentropic temperature change per cell, diluted by the rock and connate-water heat capacity (local thermal equilibrium), plus relaxation to geothermal temperature. There is no heat advection with the flowing gas and no injected-gas temperature, so near-well cooling or heating fronts from injecting gas at a different temperature are not represented.
- **Wells.** Vertical wells with a Peaceman index (isotropic). No wellbore hydraulics, friction, tubing head pressure, or wellbore hydrostatic head between 3D perforations (a BHP-controlled 3D well applies the same BHP to every perforated layer).
- **Anisotropy.** 3D permeability is isotropic (`k_v = k_h`).
- **Surface facilities.** 0D compression energy assumes intercooled multi-stage compression from `supply_pressure` straight to the bottom-hole pressure; there is no wellbore hydrostatic head, tubing friction, or aftercooler/dryer energy, and grid tiers report no compression energy.
- **Adsorption energy.** Sorbent stores use the lumped thermal model; heat of adsorption and adsorbed-phase energy are not modelled.

## Numerics

- First-order implicit (backward Euler) time stepping. The 0D tier sub-steps adaptively (≤ 10 bar and ≤ 5 K per accepted step by default; `adaptive=False` / `--fixed-dt` disables it); 1D/2D/3D use a fixed step per schedule period. Transmissibilities are lagged at the start of each step; the accumulation term is iterated to convergence (Picard).
- Every tier reports a `converged` metric; check it before using results.
- Pure Python: suitable for screening, teaching and benchmarking up to roughly 10⁴ cells. A 40×40 2D run over 70 days takes ~30 s; field-scale grids (10⁵+ cells) are impractical.

## Thermodynamic accuracy envelope (vs CoolProp 8.0.0)

Validated range: 250–450 K, 1–350 bar for H2; 300–370 K, 20–250 bar for CH4/CO2/N2 (CO2 excluding the near-critical region below 310 K above 70 bar).

| Property | Accuracy |
|---|---|
| H2 Z factor / density (Peneloux) | ≤ 1.0 % at 280–450 K; ≤ 1.25 % at 250–280 K (plain PR: ≤ 4.4 %) |
| H2 viscosity | ≤ 0.5 % |
| H2 real-gas Cp | ≤ 0.3 %; γ = Cp/Cv ≤ 2 % |
| H2 ΔH, ΔU (Peneloux) | ~5 J/mol over 1–200 bar |
| H2 isentropic ΔT | ≤ 1.6 % |
| H2 Joule-Thomson coefficient | −46 % … +23 % (10–300 bar); ±6 % near 100 bar |
| CH4 / CO2 / N2 density | ≤ 1.2 % / 4.2 % / 0.2 % |
| CH4 / CO2 / N2 viscosity | ≤ 1.0 % / 3.2 % / 2.1 % |
| Mixture JT and residual enthalpy | Screening estimates only (cubic EOS limitation) |

## Validation status

| Level | Status |
|---|---|
| Properties vs reference EOS | Done (CoolProp 8.0.0; `audit/bench_thermo.py`, tests) |
| 0D adiabatic fill/discharge vs reference energy balance | Done (≤ 0.5 % fill, ≤ 2 % discharge) |
| 1D vs analytical real-gas line source (Al-Hussainy) | Done (BHP ≤ 0.23 %) |
| 1D/2D/3D vs OPM Flow (single-phase-equivalent, isothermal) | Done (stored-gas change ≤ 1.2 %, pressure ≤ 0.05 %; `audit/bench_opm.py`) |
| Two-phase H2–brine benchmarks (e.g. Hogeweg et al. 2022) | **Not done**: requires saturation transport, which the solvers do not have |
| Field history match | Not done |
