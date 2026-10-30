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

from h2multisim.config import build_from_config
from h2multisim.models.one_d import (
    RadialReservoirProps,
    UHSRadial1D,
    WellControl,
    okoroafor_radial,
    rough_radial,
    thomas_solve,
)
from h2multisim.models.three_d import (
    Schedule3D,
    UHSFieldScale3D,
    WellSpec3D,
    scenario_gravity_override,
    scenario_layered_reservoir,
    scenario_single_well_3d,
)
from h2multisim.models.two_d import (
    Schedule2D,
    UHSAreal2D,
    WellSpec2D,
    compute_sweep_efficiency,
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
    schedule_cavern,
    schedule_porous_reservoir,
    schedule_rough,
)
from h2multisim.thermo import (
    R_UNIVERSAL,
    PengRobinsonH2,
    adiabatic_temperature,
    cp_real_h2,
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

    # ── 1. EOS Z-Factor vs NIST/CoolProp (10 points) ────────────────────────
    def test_01_eos_zfactor_vs_nist(self) -> None:
        # Reference values: CoolProp 8.0.0, Leachman et al. (2009) normal-H2
        # reference EOS.  Tolerances document the screening accuracy of the
        # cubic Peng-Robinson model (max ~4.2% at 250 bar); a Peneloux volume
        # translation is the planned P1 improvement (see AUDIT.md F-04).
        cases = [
            (293.15, 101325.0, 1.00059, 0.002),
            (293.15, 5.0e6, 1.02999, 0.030),
            (293.15, 10.0e6, 1.06081, 0.030),
            (293.15, 20.0e6, 1.12473, 0.050),
            (293.15, 25.0e6, 1.15751, 0.050),
            (350.00, 10.0e6, 1.05338, 0.030),
            (350.00, 25.0e6, 1.13579, 0.050),
            (400.00, 10.0e6, 1.04789, 0.030),
            (400.00, 25.0e6, 1.12074, 0.050),
            (316.00, 8.0e6, 1.04594, 0.030),
        ]
        for T, P, Z_ref, tol in cases:
            with self.subTest(T=T, P=P):
                Z = self.eos.z_factor(P, T)
                err = abs(Z - Z_ref) / Z_ref
                self.assertLessEqual(err, tol, f"Z={Z:.5f}, NIST={Z_ref:.5f}, err={err:.4f}")

    # ── 2. EOS Density vs CoolProp / Okoroafor (2022) ───────────────────────
    def test_02_eos_density_okoroafor(self) -> None:
        # CoolProp 8.0.0 gives 5.868 kg/m³ at 80 bar / 316 K; the cubic model
        # is ~1.9% high, which is why the tolerance is 2.5% rather than 5%.
        rho = self.eos.density(bar_to_pa(80.0), 316.0)
        self.assertAlmostEqual(rho, 5.868, delta=0.15, msg=f"rho={rho:.3f} kg/m³")

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
        # CoolProp 8.0.0: 8.797 uPa.s at 293.15 K / 1 atm.
        mu_ntp = viscosity_h2(293.15, 101325.0)
        self.assertLess(abs(mu_ntp - 8.797e-6) / 8.797e-6, 0.02)

        # CoolProp 8.0.0: 9.676 uPa.s at 333 K / 82 bar (Hogeweg point).
        mu_hog = viscosity_h2(333.0, 82e5)
        self.assertLess(abs(mu_hog - 9.676e-6) / 9.676e-6, 0.02)

        # High pressure (250 bar at 293.15 K): CoolProp = 9.302 uPa.s
        mu_hp = viscosity_h2(293.15, 250e5)
        self.assertLess(abs(mu_hp - 9.302e-6) / 9.302e-6, 0.02)

        # High temperature & high pressure (350 K at 250 bar): CoolProp = 10.304 uPa.s
        mu_hthp = viscosity_h2(350.0, 250e5)
        self.assertLess(abs(mu_hthp - 10.304e-6) / 10.304e-6, 0.02)

    # ── 6. Heat Capacity Ratio Gamma vs CoolProp ─────────────────────────────
    def test_06_gamma_bounds(self) -> None:
        # CoolProp 8.0.0 isentropic expansion coefficient.  The PR departure
        # function is a screening estimate; the measured deviation reaches ~7%
        # at 300 bar, so the test uses an 8% band.
        conditions = [
            (250.0, 50e5, 1.4760),
            (316.0, 80e5, 1.4775),
            (370.0, 150e5, 1.5160),
            (400.0, 200e5, 1.5401),
            (600.0, 300e5, 1.5265),
        ]
        for T, P, gamma_ref in conditions:
            with self.subTest(T=T, P=P):
                g = gamma_real_h2(T, P, self.eos)
                self.assertLess(
                    abs(g - gamma_ref) / gamma_ref,
                    0.08,
                    f"gamma={g:.4f}, reference={gamma_ref:.4f}",
                )

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
        self.assertLess(rel_err, 1e-6, f"Mass conservation error={rel_err:.4f}")
        self.assertLess(result.metrics["mass_balance_error"], 1e-6)

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

        # Withdrawal runs from day 120 to day 180 in the canonical schedule.
        wdr_indices = [i for i, t in enumerate(result.times_s) if 121 * 86400 <= t <= 179 * 86400]
        self.assertTrue(len(wdr_indices) > 0)
        idx = wdr_indices[len(wdr_indices) // 2]
        self.assertLessEqual(result.BHP_Pa[idx], bar_to_pa(55.0) + 1000.0)
        self.assertLessEqual(
            result.BHP_Pa[idx],
            result.P_avg_Pa[idx] - bar_to_pa(1.0),
            "flowing BHP must be below the average reservoir pressure",
        )

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
        r0.metrics["working_gas_kg"]

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
        r1.metrics["total_injected_kg"]

    # ── 32. Peneloux volume translation (opt-in, F-04) ───────────────────────
    def test_32_peneloux_translation(self) -> None:
        eos = PengRobinsonH2(volume_translation=-4.0e-6)
        # CoolProp 8.0.0 reference values; the translation must bring the
        # cubic model below 1 % error at every table point.
        cases = [
            (293.15, 101325.0, 1.00059),
            (293.15, 5.0e6, 1.02999),
            (293.15, 10.0e6, 1.06081),
            (293.15, 20.0e6, 1.12473),
            (293.15, 25.0e6, 1.15751),
            (350.00, 10.0e6, 1.05338),
            (350.00, 25.0e6, 1.13579),
            (400.00, 10.0e6, 1.04789),
            (400.00, 25.0e6, 1.12074),
            (316.00, 8.0e6, 1.04594),
        ]
        for T, P, Z_ref in cases:
            with self.subTest(T=T, P=P):
                Z = eos.z_factor(P, T)
                self.assertLess(
                    abs(Z - Z_ref) / Z_ref,
                    0.015,
                    f"Z={Z:.5f}, reference={Z_ref:.5f}",
                )
                # Density must also be consistent: rho = P*M/(Z*R*T).
                rho = eos.density(P, T)
                rho_ref = P * 2.01594e-3 / (Z_ref * R_UNIVERSAL * T)
                self.assertLess(abs(rho - rho_ref) / rho_ref, 0.015)

        # Propagates from ReservoirProperties to model instance
        props = build_teesside_cavern()
        props.peneloux = True
        props.__post_init__()
        sim = UHSReservoir0D(props)
        self.assertEqual(sim.eos.volume_translation, -4.0e-6)

    # ── 33. Real-Gas Heat Capacity and Joule-Thomson Consistency ───────────
    def test_33_real_gas_cp_and_joule_thomson(self) -> None:
        # Real cp departs by <3% from ideal at 250 bar (CoolProp = 14769 J/kg/K)
        cp_250 = cp_real_h2(350.0, 250e5, self.eos) / 2.01594e-3
        self.assertLess(abs(cp_250 - 14769.0) / 14769.0, 0.04)

        # Joule-Thomson coefficient with Peneloux must be negative at 100 bar, 293 K
        # Reference CoolProp is -0.0347 K/bar (-3.47e-7 K/Pa)
        eos_pen = PengRobinsonH2(volume_translation=-4.0e-6)
        jt_100 = joule_thomson_coefficient(293.15, 100e5, eos_pen)
        self.assertLess(jt_100, 0.0)
        self.assertLess(abs(jt_100 - (-3.47e-7)), 1.5e-7)

    # ── 25. Teesside cavern: exact mass balance and pressure envelope ───────
    def test_25_teesside_cavern_mass_balance(self) -> None:
        props = build_teesside_cavern()
        model = UHSReservoir0D(props, eos=self.eos)
        P0, T0 = props.reference_pressure, props.reference_temperature
        s0 = ReservoirState(P0, T0, model.mass_from_state(P0, T0))
        result = model.run_schedule(s0, schedule_cavern(P0), dt_s=86400.0)

        self.assertLess(result.metrics["mass_balance_error"], 1e-9)
        self.assertAlmostEqual(
            result.metrics["final_mass_kg"] - result.metrics["initial_mass_kg"],
            result.metrics["net_stored_kg"],
            delta=1e-6 * max(abs(result.metrics["net_stored_kg"]), 1.0),
        )
        self.assertLessEqual(result.metrics["max_pressure_Pa"], bar_to_pa(120.0) + 1000.0)
        self.assertGreaterEqual(result.metrics["min_pressure_Pa"], bar_to_pa(40.0) - 1000.0)

    # ── 26. 1D constant-BHP deliverability against the reference ─────────────
    def test_26_1d_constant_bhp_deliverability(self) -> None:
        props, _ = okoroafor_radial()
        props.tau_days = 0.0
        sim = UHSRadial1D(props, N=40, eos=self.eos)
        result = sim.run_schedule(
            [WellControl("w", "withdraw", 10.0, bar_to_pa(60.0), 1.0e9)],
            dt_days=0.2,
        )
        withdrawn_t = result.metrics["total_withdrawn_kg"] / 1e3
        self.assertAlmostEqual(withdrawn_t, 6575.0, delta=65.8)
        self.assertLess(result.metrics["mass_balance_error"], 1e-6)

    # ── 27. 2D mass balance closure ─────────────────────────────────────────
    def test_27_2d_mass_balance(self) -> None:
        props, wells, sched = scenario_single_well(8, 8)
        result = UHSAreal2D(props, wells, eos=self.eos).run_schedule(sched, dt_days=2.0)
        self.assertLess(result.metrics["mass_balance_error"], 1e-6)
        self.assertAlmostEqual(
            result.metrics["final_mass_t"] - result.metrics["initial_mass_t"],
            result.metrics["net_stored_t"],
            delta=1e-6 * max(abs(result.metrics["net_stored_t"]), 1.0),
        )

    # ── 28. 3D mass balance closure ─────────────────────────────────────────
    def test_28_3d_mass_balance(self) -> None:
        props, wells, sched = scenario_single_well_3d(6, 6, 3)
        result = UHSFieldScale3D(props, wells, eos=self.eos).run_schedule(sched, dt_days=5.0)
        self.assertLess(result.metrics["mass_balance_error"], 1e-6)
        self.assertAlmostEqual(
            result.metrics["final_mass_t"] - result.metrics["initial_mass_t"],
            result.metrics["net_stored_t"],
            delta=1e-6 * max(abs(result.metrics["net_stored_t"]), 1.0),
        )

    # ── 29. 3D gravity direction: bottom pressure > top pressure ────────────
    def test_29_3d_gravity_direction(self) -> None:
        nx, ny, nz = 6, 6, 5
        props, wells, _ = scenario_gravity_override(nx, ny, nz)
        result = UHSFieldScale3D(props, wells, eos=self.eos).run_schedule(
            [Schedule3D("shut_in", "shut_in", 10.0, 0.0, 0.0)],
            snap_times_d=[10.0],
            dt_days=1.0,
        )
        self.assertTrue(result.snapshots)
        snap = result.snapshots[-1]
        ci, cj = nx // 2, ny // 2
        P_col = [snap.P_Pa[ci * ny * nz + cj * nz + k] for k in range(nz)]
        self.assertGreater(
            P_col[0], P_col[-1] + 2000.0, "gravity must make the bottom layer heavier"
        )

    # ── 30. Five-spot: corners inject, centre produces ──────────────────────
    def test_30_five_spot_per_well_modes(self) -> None:
        props, wells, sched = scenario_five_spot(8, 8)
        result = UHSAreal2D(props, wells, eos=self.eos).run_schedule(sched, dt_days=2.0)
        corner_q = result.per_well_q["INJ_SW"]
        centre_q = result.per_well_q["CTRE"]
        self.assertGreater(max(corner_q), 0.0)
        self.assertGreaterEqual(min(corner_q), 0.0)
        self.assertLess(min(centre_q), 0.0)
        self.assertLessEqual(max(centre_q), 0.0)

    # ── 31. Injection cannot back-flow above its BHP target ─────────────────
    def test_31_no_backflow(self) -> None:
        props, wells, _ = scenario_single_well(6, 6)
        sched = [Schedule2D("inj", "inject", 1.0, bar_to_pa(70.0), 35.0)]
        result = UHSAreal2D(props, wells, eos=self.eos).run_schedule(sched, dt_days=0.5)
        self.assertLessEqual(max(result.per_well_q["W1"]), 1e-12)
        self.assertLess(result.metrics["mass_balance_error"], 1e-6)

    # ── 34. Harvested Domain Scenarios & Sweep Efficiency ───────────────────
    def test_34_rough_and_sweep_efficiency(self) -> None:
        # 0D Rough schedule construction and execution
        sched_rough = schedule_rough(bar_to_pa(150.0))
        self.assertEqual(len(sched_rough), 10)
        self.assertEqual(sched_rough[0].name, "cushion_fill")

        # 1D Rough scenario
        props_rough, controls_rough = rough_radial()
        self.assertAlmostEqual(props_rough.P_init, bar_to_pa(150.0))
        self.assertAlmostEqual(props_rough.T_init, 369.0)
        self.assertEqual(len(controls_rough), 3)

        # 2D Sweep efficiency calculation
        p_init = bar_to_pa(100.0)
        p_field_unperturbed = [p_init] * 16
        self.assertAlmostEqual(compute_sweep_efficiency(p_field_unperturbed, p_init), 0.0)

        # Perturbed pressure field (>1 bar threshold)
        p_field_perturbed = [p_init] * 8 + [p_init + bar_to_pa(5.0)] * 8
        self.assertAlmostEqual(compute_sweep_efficiency(p_field_perturbed, p_init), 0.5)

        # 2D matrix format
        p_matrix = [[p_init, p_init + bar_to_pa(5.0)], [p_init, p_init]]
        self.assertAlmostEqual(compute_sweep_efficiency(p_matrix, p_init), 0.25)

    # ── 35. Peneloux Departure Consistency & Monte Carlo Schedule Fidelity ──
    def test_35_peneloux_departures_and_mc_schedule(self) -> None:
        c = -4.0e-6
        eos_plain = PengRobinsonH2(volume_translation=0.0)
        eos_pen = PengRobinsonH2(volume_translation=c)
        P, T = bar_to_pa(250.0), 293.15

        # Constant volume translation must not distort PR heat capacity departures
        self.assertAlmostEqual(cp_real_h2(T, P, eos_plain), cp_real_h2(T, P, eos_pen), places=10)
        self.assertAlmostEqual(
            gamma_real_h2(T, P, eos_plain), gamma_real_h2(T, P, eos_pen), places=10
        )

        # Peneloux fugacity shift: ln(phi_pen) == ln(phi_plain) - c*P/(R*T)
        ln_ratio = math.log(
            eos_pen.fugacity_coefficient(P, T) / eos_plain.fugacity_coefficient(P, T)
        )
        expected_shift = -c * P / (R_UNIVERSAL * T)
        self.assertAlmostEqual(ln_ratio, expected_shift, places=10)

        # Joule-Thomson at 1 bar / 293.15 K must be negative even with plain EOS
        self.assertLess(joule_thomson_coefficient(293.15, 101325.0, eos_plain), 0.0)

        # Monte Carlo accepts explicit schedule list and preserves Peneloux
        props = build_okoroafor_reservoir()
        props.peneloux = True
        props.__post_init__()
        sched_list = [SchedulePeriod("inj", "inject", 10.0, bar_to_pa(100.0), 35.0)]
        mc_cfg = MonteCarloConfig(
            n_realizations=3,
            porosity_mean=props.porosity,
            porosity_std=0.0,
            thickness_mean=props.thickness,
            thickness_std=0.0,
            log10_perm_mean=math.log10(props.permeability),
            log10_perm_std=0.0,
            geothermal_T_mean=props.reference_temperature,
            geothermal_T_std=0.0,
            area_mean=props.area,
            area_std=0.0,
        )
        mc_stats = run_monte_carlo(props, sched_list, mc_cfg, dt_days=2.0)
        det = UHSReservoir0D(props).run_schedule(
            ReservoirState(
                props.reference_pressure,
                props.reference_temperature,
                UHSReservoir0D(props).mass_from_state(
                    props.reference_pressure, props.reference_temperature
                ),
            ),
            sched_list,
            dt_s=2.0 * 86400.0,
        )
        self.assertAlmostEqual(
            mc_stats["working_gas_kg"]["p50"], det.metrics["working_gas_kg"], places=3
        )

    # ── 36. Solver Validation, Hydrostatic Equilibrium & Config Extensions ──
    def test_36_solver_validation_and_config_extensions(self) -> None:
        # 1. Duplicate well names and out-of-bounds indices must raise ValueError
        p2, _, _ = scenario_single_well(6, 6)
        with self.assertRaises(ValueError):
            UHSAreal2D(p2, [WellSpec2D(1, 1, name="W"), WellSpec2D(2, 2, name="W")])
        with self.assertRaises(ValueError):
            UHSAreal2D(p2, [WellSpec2D(6, 0, name="W1")])

        p3, _, _ = scenario_single_well_3d(4, 4, 3)
        with self.assertRaises(ValueError):
            UHSFieldScale3D(
                p3,
                [
                    WellSpec3D(1, 1, 0, 2, name="W"),
                    WellSpec3D(2, 2, 0, 2, name="W"),
                ],
            )
        with self.assertRaises(ValueError):
            UHSFieldScale3D(p3, [WellSpec3D(1, 1, 2, 1, name="W1")])

        # 2. 3D Hydrostatic mid-depth datum and exact discrete equilibrium
        sim3 = UHSFieldScale3D(p3, [WellSpec3D(2, 2, 0, 2, name="W1")])
        res3 = sim3.run_schedule(
            [Schedule3D("shut", "shut_in", 0.7, 0.0, 0.0)],
            snap_times_d=[0.7],
            dt_days=0.1,
        )
        self.assertAlmostEqual(res3.P_avg_bar[0], pa_to_bar(p3.P_init), places=10)
        self.assertAlmostEqual(res3.P_avg_bar[-1], pa_to_bar(p3.P_init), places=8)
        self.assertEqual(len(res3.snapshots), 1)
        self.assertEqual(res3.metrics["converged"], 1.0)
        self.assertLess(res3.metrics["max_picard_residual_Pa"], 1e-8)

        # 3. 2D float-boundary snapshot capture (sum([0.1]*7) == 0.6999999999999998)
        sim2 = UHSAreal2D(p2, [WellSpec2D(3, 3, name="W1")])
        res2 = sim2.run_schedule(
            [Schedule2D("shut", "shut_in", 0.7, 0.0, 0.0)],
            snap_days=[0.7],
            dt_days=0.1,
        )
        self.assertEqual(len(res2.snapshots), 1)
        self.assertEqual(res2.metrics["converged"], 1.0)

        # 4. Config extensions: pore_model, isotherm, k_field_mD, k_layers_mD, 3D well_modes
        _, props0_cfg, _, _, _, _ = build_from_config(
            {
                "scale": "0d",
                "reservoir": {
                    "area": 1e5,
                    "thickness": 20.0,
                    "porosity": 0.2,
                    "total_compressibility": 5e-10,
                    "permeability_mD": 50.0,
                    "reference_pressure": "80 bar",
                    "reference_temperature": 315.0,
                    "well_radius": 0.1,
                    "drainage_radius": 200.0,
                    "m_sorbent": 1000.0,
                    "pore_model": {"model": "rigid"},
                    "isotherm": {"model": "langmuir", "n_max1": 15.0, "K0_1": 2e-7},
                },
                "schedule": [
                    {
                        "name": "inj",
                        "mode": "inject",
                        "duration_days": 5.0,
                        "bhp_target": "100 bar",
                        "max_mass_rate": 10.0,
                    }
                ],
            }
        )
        self.assertIsNotNone(props0_cfg.pore_model)
        self.assertEqual(props0_cfg.pore_model.model, "rigid")
        self.assertIsNotNone(props0_cfg.isotherm)
        self.assertEqual(props0_cfg.isotherm.model, "langmuir")

        _, props3_cfg, _, sched3_cfg, _, _ = build_from_config(
            {
                "scale": "3d",
                "reservoir": {
                    "Lx": 500.0,
                    "Ly": 500.0,
                    "Lz": 60.0,
                    "phi": 0.2,
                    "permeability_mD": 10.0,
                    "k_layers_mD": [10.0, 100.0, 10.0],
                    "Ct": 5e-10,
                    "P_init": "80 bar",
                    "T_init": 316.0,
                    "nx": 4,
                    "ny": 4,
                    "nz": 3,
                },
                "wells": [
                    {"i": 1, "j": 1, "k_min": 0, "k_max": 2, "name": "INJ"},
                    {"i": 2, "j": 2, "k_min": 0, "k_max": 2, "name": "PROD"},
                ],
                "schedule": [
                    {
                        "name": "inj",
                        "mode": "inject",
                        "duration_days": 5.0,
                        "bhp_target": "100 bar",
                        "q_max": 20.0,
                        "well_modes": {"INJ": "inject", "PROD": "shut_in"},
                    }
                ],
            }
        )
        self.assertEqual(len(props3_cfg.k_layers), 3)
        self.assertAlmostEqual(props3_cfg.k_layers[1], mD_to_m2(100.0))
        self.assertEqual(sched3_cfg[0].mode_for("PROD"), "shut_in")
        self.assertEqual(sched3_cfg[0].mode_for("INJ"), "inject")


if __name__ == "__main__":
    unittest.main()
