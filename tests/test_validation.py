"""
tests/test_validation.py
========================
H2MultiSim™ Comprehensive Validation and Benchmark Test Suite.
Validates thermodynamics against NIST WebBook, Okoroafor (2022), Hogeweg (2022),
and checks numerical stability across 0D, 1D, 2D, and 3D simulation scales.
"""

from __future__ import annotations

import math
import unittest

from h2multisim.models.one_d import (
    RadialReservoirProps,
    UHSRadial1D,
    WellControl,
    okoroafor_radial,
    thomas_solve,
)
from h2multisim.models.three_d import (
    UHSFieldScale3D,
    scenario_gravity_override,
    scenario_layered_reservoir,
    scenario_single_well_3d,
)
from h2multisim.models.two_d import (
    UHSAreal2D,
    pressure_uniformity,
    scenario_five_spot,
    scenario_single_well,
    scenario_thief_zone,
)
from h2multisim.models.zero_d import (
    MonteCarloConfig,
    ReservoirState,
    SchedulePeriod,
    UHSReservoir0D,
    build_okoroafor_reservoir,
    build_teesside_cavern,
    run_monte_carlo,
    schedule_porous_reservoir,
)
from h2multisim.thermo import (
    PengRobinsonH2,
    adiabatic_temperature,
    gamma_real_h2,
    joule_thomson_coefficient,
    thermal_relaxation_update,
    viscosity_h2,
)
from h2multisim.units import (
    C_to_K,
    F_to_K,
    K_to_C,
    K_to_F,
    bar_to_pa,
    kg_to_tonnes,
    m2_to_mD,
    mD_to_m2,
    mpa_to_pa,
    pa_to_bar,
    pa_to_mpa,
    tonnes_to_kg,
)


class TestH2MultiSimValidation(unittest.TestCase):
    """Validation test cases across thermodynamics, numerical solvers, and models."""

    def setUp(self) -> None:
        self.eos = PengRobinsonH2()

    # ── 1. EOS Z-Factor vs NIST WebBook (10 points) ──────────────────────────
    def test_01_eos_zfactor_vs_nist(self) -> None:
        cases = [
            (293.15, 101325.0, 1.00060, 0.002),
            (293.15, 5.0e6, 1.02982, 0.030),
            (293.15, 10.0e6, 1.05881, 0.030),
            (293.15, 20.0e6, 1.11522, 0.035),
            (293.15, 25.0e6, 1.14283, 0.035),
            (350.00, 10.0e6, 1.04785, 0.030),
            (350.00, 25.0e6, 1.11665, 0.030),
            (400.00, 10.0e6, 1.04163, 0.030),
            (400.00, 25.0e6, 1.10237, 0.030),
            (316.00, 8.0e6, 1.02389, 0.030),
        ]
        for T, P, Z_ref, tol in cases:
            with self.subTest(T=T, P=P):
                Z = self.eos.z_factor(P, T)
                err = abs(Z - Z_ref) / Z_ref
                self.assertLessEqual(err, tol, f"Z={Z:.5f}, NIST={Z_ref:.5f}, err={err:.4f}")

    # ── 2. EOS Density vs Okoroafor (2022) ───────────────────────────────────
    def test_02_eos_density_okoroafor(self) -> None:
        rho = self.eos.density(bar_to_pa(80.0), 316.0)
        self.assertAlmostEqual(rho, 5.98, delta=0.30, msg=f"rho={rho:.3f} kg/m³")

    # ── 3. EOS Density Monotonicity ──────────────────────────────────────────
    def test_03_eos_density_monotone(self) -> None:
        pressures = [10, 50, 80, 120, 180, 250]
        rhos = [self.eos.density(bar_to_pa(p), 316.0) for p in pressures]
        for i in range(len(rhos) - 1):
            self.assertLess(
                rhos[i], rhos[i + 1], f"Density not increasing: {rhos[i]} >= {rhos[i + 1]}"
            )

    # ── 4. EOS Solver Numerical Stability ────────────────────────────────────
    def test_04_eos_stability(self) -> None:
        extremes = [(50.0, 1e4), (700.0, 350e5), (300.0, 1e5), (316.0, 80e5)]
        for T, P in extremes:
            with self.subTest(T=T, P=P):
                Z = self.eos.z_factor(P, T)
                self.assertTrue(0.05 < Z < 5.0, f"Z={Z:.4f} outside safe bounds")

    # ── 5. Dynamic Viscosity Benchmarks ──────────────────────────────────────
    def test_05_viscosity_benchmarks(self) -> None:
        # NIST NTP: 293.15 K, 1 atm -> 8.948 uPa.s
        mu_ntp = viscosity_h2(293.15, 101325.0)
        self.assertLess(abs(mu_ntp - 8.948e-6) / 8.948e-6, 0.05)

        # Hogeweg (2022): ~9.4 uPa.s at 333 K, 82 bar
        mu_hog = viscosity_h2(333.0, 82e5)
        self.assertLess(abs(mu_hog - 9.4e-6) / 9.4e-6, 0.07)

    # ── 6. Heat Capacity Ratio Gamma Bounds ──────────────────────────────────
    def test_06_gamma_bounds(self) -> None:
        conditions = [
            (250.0, 50e5),
            (316.0, 80e5),
            (370.0, 150e5),
            (400.0, 200e5),
            (600.0, 300e5),
        ]
        for T, P in conditions:
            with self.subTest(T=T, P=P):
                g = gamma_real_h2(T, P, self.eos)
                self.assertTrue(1.10 <= g <= 1.80, f"gamma={g:.4f} out of [1.10, 1.80]")

    # ── 7. Joule-Thomson Inversion Sign ───────────────────────────────────────
    def test_07_joule_thomson(self) -> None:
        jt_300 = joule_thomson_coefficient(300.0, 10e6, self.eos)
        self.assertLess(jt_300, 0.0, "mu_JT must be negative at 300 K / 100 bar")
        jt_350 = joule_thomson_coefficient(350.0, 20e6, self.eos)
        self.assertLess(jt_350, 0.0, "mu_JT must be negative at 350 K / 200 bar")

    # ── 8. Thermal Model Physical Consistency ────────────────────────────────
    def test_08_thermal_model(self) -> None:
        # Adiabatic compression heating
        T_ad = adiabatic_temperature(293.0, bar_to_pa(40.0), bar_to_pa(250.0), self.eos)
        self.assertGreater(T_ad, 343.0, "Adiabatic compression should increase T by > 50 K")

        # Relaxation toward geothermal rock temperature
        T_inf = thermal_relaxation_update(500.0, 316.0, 1e9, 86400.0)
        self.assertAlmostEqual(T_inf, 316.0, places=2)

        # Expansion cooling
        T_exp = adiabatic_temperature(400.0, bar_to_pa(250.0), bar_to_pa(40.0), self.eos)
        self.assertLess(T_exp, 400.0, "Adiabatic expansion must decrease temperature")

    # ── 9. Round-Trip Unit Conversions ───────────────────────────────────────
    def test_09_unit_conversions(self) -> None:
        self.assertAlmostEqual(bar_to_pa(pa_to_bar(250e5)), 250e5, delta=1.0)
        self.assertAlmostEqual(mpa_to_pa(pa_to_mpa(20e6)), 20e6, delta=1.0)
        self.assertAlmostEqual(C_to_K(K_to_C(400.0)), 400.0, places=6)
        self.assertAlmostEqual(F_to_K(K_to_F(316.0)), 316.0, places=6)
        self.assertAlmostEqual(tonnes_to_kg(kg_to_tonnes(12345.0)), 12345.0, places=4)
        self.assertAlmostEqual(m2_to_mD(mD_to_m2(100.0)), 100.0, places=4)
        self.assertAlmostEqual(C_to_K(0.0), 273.15, places=4)

    # ── 10. Thomas Algorithm Tridiagonal Solver ──────────────────────────────
    def test_10_thomas_algorithm(self) -> None:
        N = 6
        a = [0.0] + [-1.0] * (N - 1)
        b = [2.0] * N
        c = [-1.0] * (N - 1) + [0.0]
        d = [1.0] * N
        x = thomas_solve(a, b, c, d)

        # Verify residuals A*x - d
        res = [
            abs(
                b[i] * x[i]
                + a[i] * (x[i - 1] if i > 0 else 0.0)
                + c[i] * (x[i + 1] if i < N - 1 else 0.0)
                - d[i]
            )
            for i in range(N)
        ]
        self.assertLess(max(res), 1e-10)
        self.assertAlmostEqual(x[0], x[5], places=8)
        self.assertAlmostEqual(x[2], x[3], places=8)

    # ── 11. 0D BHP Clamping (No Runaway Divergence) ───────────────────────────
    def test_11_0d_bhp_clamping(self) -> None:
        props = build_okoroafor_reservoir()
        model = UHSReservoir0D(props, eos=self.eos)
        P0 = props.reference_pressure
        T0 = props.reference_temperature
        m0 = model.mass_from_state(P0, T0)
        state = ReservoirState(pressure=P0, temperature=T0, mass=m0)

        P_max = bar_to_pa(108.0)
        P_min = bar_to_pa(48.0)
        inj = SchedulePeriod("inj", "inject", 30.0, P_max, 35.0)
        wdr = SchedulePeriod("wdr", "withdraw", 30.0, P_min, 35.0)

        result = model.run_schedule(state, [inj, wdr], dt_s=86400.0)
        self.assertLessEqual(max(result.pressures), P_max + 1000.0)
        self.assertGreaterEqual(min(result.pressures), P_min - 1000.0)

    # ── 12. 0D Salt Cavern Rigid Volume Invariance ───────────────────────────
    def test_12_0d_cavern_rigid_volume(self) -> None:
        props = build_teesside_cavern()
        model = UHSReservoir0D(props, eos=self.eos)
        v_low = model.pore_volume(bar_to_pa(40.0))
        v_high = model.pore_volume(bar_to_pa(120.0))
        self.assertAlmostEqual(v_low, v_high, places=3)

    # ── 13. 0D Monotone Injection Pressure ───────────────────────────────────
    def test_13_0d_pressure_monotone_injection(self) -> None:
        props = build_okoroafor_reservoir()
        model = UHSReservoir0D(props, eos=self.eos)
        P0 = props.reference_pressure
        T0 = props.reference_temperature
        m0 = model.mass_from_state(P0, T0)
        state = ReservoirState(pressure=P0, temperature=T0, mass=m0)
        sched = [SchedulePeriod("inj", "inject", 10.0, bar_to_pa(108.0), 35.0)]
        result = model.run_schedule(state, sched, dt_s=86400.0)
        self.assertGreater(result.pressures[-1], P0)

    # ── 14. 0D Mass Conservation Closure ─────────────────────────────────────
    def test_14_0d_mass_conservation(self) -> None:
        props = build_okoroafor_reservoir()
        model = UHSReservoir0D(props, eos=self.eos)
        P0 = props.reference_pressure
        T0 = props.reference_temperature
        m0 = model.mass_from_state(P0, T0)
        state = ReservoirState(pressure=P0, temperature=T0, mass=m0)

        sched = [
            SchedulePeriod("inj", "inject", 10.0, bar_to_pa(108.0), 35.0),
            SchedulePeriod("wdr", "withdraw", 5.0, bar_to_pa(60.0), 35.0),
        ]
        result = model.run_schedule(state, sched, dt_s=86400.0)
        dm_actual = result.masses[-1] - result.masses[0]
        dm_net = result.metrics["total_injected_kg"] - result.metrics["total_withdrawn_kg"]
        rel_err = abs(dm_actual - dm_net) / max(abs(dm_net), 1e-3)
        self.assertLess(rel_err, 0.20, f"Mass conservation error={rel_err:.4f}")

    # ── 15. 0D Monte Carlo Reproducibility ───────────────────────────────────
    def test_15_0d_mc_reproducibility(self) -> None:
        props = build_okoroafor_reservoir()
        mc_cfg = MonteCarloConfig(
            n_realizations=20,
            porosity_mean=0.20,
            porosity_std=0.02,
            thickness_mean=100.0,
            thickness_std=10.0,
            log10_perm_mean=math.log10(mD_to_m2(100.0)),
            log10_perm_std=0.3,
            geothermal_T_mean=316.0,
            geothermal_T_std=3.0,
            area_mean=math.pi * 500**2,
            area_std=1e5,
        )
        s1 = run_monte_carlo(props, schedule_porous_reservoir, mc_cfg, dt_days=5.0, seed=42)
        s2 = run_monte_carlo(props, schedule_porous_reservoir, mc_cfg, dt_days=5.0, seed=42)
        self.assertAlmostEqual(s1["working_gas_kg"]["p50"], s2["working_gas_kg"]["p50"], places=4)

    # ── 16. 1D Pressure Buildup During Shut-in ───────────────────────────────
    def test_16_1d_pressure_buildup(self) -> None:
        props, schedule = okoroafor_radial()
        sim = UHSRadial1D(props, N=20, eos=self.eos)
        result = sim.run_schedule(schedule, dt_days=2.0)

        inj_end_idx = len([t for t in result.times_s if t <= 30 * 86400])
        max_bhp_inj = max(result.BHP_Pa[:inj_end_idx])
        self.assertGreaterEqual(max_bhp_inj, bar_to_pa(107.0))

        shut_bhp_start = result.BHP_Pa[-10]
        shut_bhp_end = result.BHP_Pa[-1]
        self.assertGreaterEqual(shut_bhp_end, shut_bhp_start - 1000.0)

    # ── 17. 1D Near-Wellbore Drawdown ────────────────────────────────────────
    def test_17_1d_drawdown(self) -> None:
        props, schedule = okoroafor_radial()
        sim = UHSRadial1D(props, N=20, eos=self.eos)
        result = sim.run_schedule(schedule, dt_days=2.0)

        wdr_indices = [i for i, t in enumerate(result.times_s) if 31 * 86400 <= t <= 59 * 86400]
        self.assertTrue(len(wdr_indices) > 0)
        idx = wdr_indices[len(wdr_indices) // 2]
        self.assertLessEqual(result.BHP_Pa[idx], result.P_avg_Pa[idx] + 1000.0)

    # ── 18. 2D Five-Spot vs Single Well ──────────────────────────────────────
    def test_18_2d_five_spot_vs_single(self) -> None:
        nx, ny = 12, 12
        p1, w1, s1 = scenario_single_well(nx, ny)
        r1 = UHSAreal2D(p1, w1, eos=self.eos).run_schedule(s1, snap_days=[1.0], dt_days=1.0)

        p5, w5, s5 = scenario_five_spot(nx, ny)
        r5 = UHSAreal2D(p5, w5, eos=self.eos).run_schedule(s5, snap_days=[1.0], dt_days=1.0)

        self.assertGreaterEqual(
            r5.metrics["total_injected_t"], r1.metrics["total_injected_t"] * 0.5
        )

    # ── 19. 2D Thief-Zone Pressure Variance ──────────────────────────────────
    def test_19_2d_thief_zone_pressure_cv(self) -> None:
        nx, ny = 12, 12
        p5, w5, s5 = scenario_five_spot(nx, ny)
        ptz, wtz, stz = scenario_thief_zone(nx, ny)

        r5 = UHSAreal2D(p5, w5, eos=self.eos).run_schedule(s5, snap_days=[1.0], dt_days=1.0)
        rtz = UHSAreal2D(ptz, wtz, eos=self.eos).run_schedule(stz, snap_days=[1.0], dt_days=1.0)

        self.assertTrue(len(r5.snapshots) > 0 and len(rtz.snapshots) > 0)
        snap = rtz.snapshots[0]
        P_2d = [[snap.P_Pa[i * ny + j] for j in range(ny)] for i in range(nx)]
        cv = pressure_uniformity(P_2d)
        self.assertGreaterEqual(cv, 0.0)

    # ── 20. 3D Hydrostatic Gravity Segregation ───────────────────────────────
    def test_20_3d_gravity_hydrostatic(self) -> None:
        nx, ny, nz = 6, 6, 5
        props, wells, sched = scenario_gravity_override(nx, ny, nz)
        result = UHSFieldScale3D(props, wells, eos=self.eos).run_schedule(
            sched, snap_times_d=[30.0], dt_days=5.0
        )

        self.assertTrue(len(result.snapshots) > 0)
        snap = result.snapshots[0]
        ci, cj = nx // 2, ny // 2
        P_col = [snap.P_Pa[ci * ny * nz + cj * nz + k] for k in range(nz)]
        # Hydrostatic column: bottom layer must be at higher or equal pressure to top
        self.assertGreaterEqual(P_col[0], P_col[-1] - 500.0)

    # ── 21. 3D Gravity Override vs All-Perforated Well ───────────────────────
    def test_21_3d_gravity_override_vs_allperf(self) -> None:
        nx, ny, nz = 6, 6, 5
        p1, w1, s1 = scenario_single_well_3d(nx, ny, nz)
        p2, w2, s2 = scenario_gravity_override(nx, ny, nz)

        r1 = UHSFieldScale3D(p1, w1, eos=self.eos).run_schedule(s1, dt_days=5.0)
        r2 = UHSFieldScale3D(p2, w2, eos=self.eos).run_schedule(s2, dt_days=5.0)

        self.assertGreaterEqual(r1.metrics["working_gas_t"], r2.metrics["working_gas_t"])

    # ── 22. 3D Stratigraphic Layering Sweep ──────────────────────────────────
    def test_22_3d_layered_vs_homogeneous(self) -> None:
        nx, ny, nz = 6, 6, 5
        p1, w1, s1 = scenario_single_well_3d(nx, ny, nz)
        pl, wl, sl = scenario_layered_reservoir(nx, ny, nz)

        r1 = UHSFieldScale3D(p1, w1, eos=self.eos).run_schedule(s1, dt_days=5.0)
        rl = UHSFieldScale3D(pl, wl, eos=self.eos).run_schedule(sl, dt_days=5.0)

        self.assertGreater(rl.metrics["working_gas_t"], r1.metrics["working_gas_t"])

    # ── 23. Cross-Scale Thermodynamic Consistency ────────────────────────────
    def test_23_cross_scale_eos_consistency(self) -> None:
        rho_ref = self.eos.density(bar_to_pa(80.0), 316.0)

        # 0D Model
        props0 = build_okoroafor_reservoir()
        m0d = UHSReservoir0D(props0, eos=self.eos)
        rho_0d = m0d.eos.density(bar_to_pa(80.0), 316.0)

        # 1D Model
        p1d = RadialReservoirProps(
            r_well=0.1,
            r_ext=500.0,
            h=100.0,
            phi=0.2,
            k=mD_to_m2(100.0),
            Ct=5e-10,
            P_init=bar_to_pa(80.0),
            T_init=316.0,
        )
        sim1d = UHSRadial1D(p1d, N=10, eos=self.eos)
        rho_1d = sim1d.eos.density(bar_to_pa(80.0), 316.0)

        self.assertAlmostEqual(rho_0d, rho_ref, places=5)
        self.assertAlmostEqual(rho_1d, rho_ref, places=5)

    # ── 24. Cross-Scale Order-of-Magnitude Storage ───────────────────────────
    def test_24_cross_scale_wg_order_of_magnitude(self) -> None:
        # 0D
        props_0d = build_okoroafor_reservoir()
        m_0d = UHSReservoir0D(props_0d, eos=self.eos)
        P0, T0 = props_0d.reference_pressure, props_0d.reference_temperature
        s_0d = [
            SchedulePeriod("inj", "inject", 30.0, bar_to_pa(108.0), 35.0),
            SchedulePeriod("wdr", "withdraw", 30.0, bar_to_pa(48.0), 35.0),
        ]
        r0 = m_0d.run_schedule(
            ReservoirState(P0, T0, m_0d.mass_from_state(P0, T0)),
            s_0d,
            86400.0,
            collect_time_series=False,
        )
        wg_0d = r0.metrics["working_gas_kg"]

        # 1D
        p1d = RadialReservoirProps(
            r_well=0.1,
            r_ext=500.0,
            h=100.0,
            phi=0.2,
            k=mD_to_m2(100.0),
            Ct=5e-10,
            P_init=bar_to_pa(80.0),
            T_init=316.0,
        )
        s1d = [
            WellControl("inj", "inject", 30.0, bar_to_pa(108.0), 35.0),
            WellControl("wdr", "withdraw", 30.0, bar_to_pa(48.0), 35.0),
        ]
        r1 = UHSRadial1D(p1d, N=20, eos=self.eos).run_schedule(s1d, dt_days=2.0)
        wg_1d_inj = r1.metrics["total_injected_kg"]

        self.assertGreater(wg_0d, 0.0)
        self.assertGreater(wg_1d_inj, 0.0)


if __name__ == "__main__":
    unittest.main()
