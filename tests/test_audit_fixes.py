"""
Regression tests for the v1.3.1 audit fixes.

Reference values are CoolProp 8.0.0 (Leachman et al. 2009 normal-H2 EOS and
its GERG mixture model) and the ECLIPSE keyword item order.
"""

from __future__ import annotations

import os
import tempfile
import unittest

from h2multisim.config import _pa
from h2multisim.deck_export import export_eclipse_1d_radial, export_eclipse_2d, export_eclipse_3d
from h2multisim.models.one_d import UHSRadial1D, WellControl, okoroafor_radial
from h2multisim.models.three_d import (
    Schedule3D,
    UHSFieldScale3D,
    scenario_gravity_override,
    scenario_layered_reservoir,
)
from h2multisim.models.two_d import Schedule2D, UHSAreal2D, scenario_single_well
from h2multisim.models.zero_d import (
    MonteCarloConfig,
    ReservoirState,
    SchedulePeriod,
    UHSReservoir0D,
    build_okoroafor_reservoir,
    build_teesside_cavern,
    run_monte_carlo,
)
from h2multisim.thermo import (
    PengRobinsonH2,
    PengRobinsonMixture,
    RelativePermeabilityModel,
    adiabatic_temperature,
    cp_real_h2,
    enthalpy_molar,
    gamma_real_h2,
    internal_energy_molar,
    joule_thomson_coefficient,
    viscosity_h2,
)
from h2multisim.units import bar_to_pa, mD_to_m2


class TestRealGasThermo(unittest.TestCase):
    def test_cp_real_rises_with_pressure(self) -> None:
        eos = PengRobinsonH2(peneloux=True)
        # (T [K], P [bar], CoolProp Cp [J/(mol·K)], CoolProp Cp/Cv)
        refs = [
            (300.0, 10.0, 28.902, 1.4065),
            (300.0, 300.0, 29.886, 1.4223),
            (350.0, 200.0, 29.676, 1.4098),
        ]
        for T, Pb, cp_ref, g_ref in refs:
            with self.subTest(T=T, P=Pb):
                self.assertLess(abs(cp_real_h2(T, Pb * 1e5, eos) / cp_ref - 1.0), 0.01)
                self.assertLess(abs(gamma_real_h2(T, Pb * 1e5, eos) / g_ref - 1.0), 0.025)
        self.assertGreater(cp_real_h2(300.0, 300e5, eos), cp_real_h2(300.0, 10e5, eos))

    def test_adiabatic_temperature_vs_isentrope(self) -> None:
        # (T0 [K], P0 [bar], P1 [bar], CoolProp isentropic ΔT [K])
        refs = [
            (293.15, 40.0, 120.0, 108.39),
            (316.0, 80.0, 108.0, 28.40),
            (293.15, 120.0, 40.0, -81.16),
        ]
        for T0, P0, P1, dT_ref in refs:
            with self.subTest(T0=T0, P0=P0, P1=P1):
                dT = adiabatic_temperature(T0, P0 * 1e5, P1 * 1e5) - T0
                self.assertLess(abs(dT / dT_ref - 1.0), 0.03)

    def test_joule_thomson_mixture_uses_mixture(self) -> None:
        # CoolProp: +0.0029 K/bar for 80 % H2 / 20 % CH4 at 320 K, 100 bar (pure H2: -0.037).
        mix = PengRobinsonMixture({"H2": 0.8, "CH4": 0.2})
        self.assertGreater(joule_thomson_coefficient(320.0, 1e7, mix), 0.0)
        self.assertLess(joule_thomson_coefficient(320.0, 1e7), 0.0)

    def test_joule_thomson_pure_h2_near_100_bar(self) -> None:
        # CoolProp: -0.0347 K/bar at 293.15 K, 100 bar.
        jt = joule_thomson_coefficient(293.15, 100e5) * 1e5
        self.assertLess(abs(jt / -0.0347 - 1.0), 0.10)

    def test_cushion_gas_species_vs_coolprop(self) -> None:
        # CoolProp 8.0.0 at 320 K: (species, P [bar]) -> (density [kg/m³], viscosity [µPa·s])
        refs = {
            ("CH4", 100.0): (67.61, 14.13),
            ("CH4", 200.0): (138.23, 18.31),
            ("N2", 100.0): (103.65, 20.66),
            ("N2", 200.0): (197.00, 23.51),
            ("CO2", 100.0): (448.28, 31.76),
            ("CO2", 200.0): (802.33, 72.48),
        }
        for (sp, Pb), (rho_ref, mu_ref) in refs.items():
            with self.subTest(species=sp, P=Pb):
                eos = PengRobinsonMixture({sp: 1.0}, peneloux=True)
                self.assertLess(abs(eos.density(Pb * 1e5, 320.0) / rho_ref - 1.0), 0.04)
                self.assertLess(abs(eos.viscosity(Pb * 1e5, 320.0) / (mu_ref * 1e-6) - 1.0), 0.02)

    def test_brooks_corey_capillary_round_trip(self) -> None:
        rp = RelativePermeabilityModel(
            model="brooks_corey", swr=0.25, sgc=0.05, brooks_corey_lambda=2.0, entry_pressure_Pa=8e3
        )
        for sg in (0.2, 0.45, 0.7):
            with self.subTest(sg=sg):
                pc = rp.capillary_pressure_Pa(sg)
                self.assertAlmostEqual(rp.gas_saturation_from_capillary_pressure(pc), sg, places=9)


class TestEnergyBalance(unittest.TestCase):
    def test_enthalpy_and_internal_energy_vs_coolprop(self) -> None:
        eos = PengRobinsonH2(peneloux=True)
        # ((T1, P1 [bar]), (T2, P2 [bar]), CoolProp ΔH, ΔU [J/mol])
        refs = [
            ((300.0, 1.0), (300.0, 200.0), 207.2, -96.9),
            ((293.15, 40.0), (400.0, 120.0), 3226.6, 2205.1),
            ((369.0, 150.0), (300.0, 50.0), -2143.2, -1406.5),
        ]
        for (T1, P1), (T2, P2), dh_ref, du_ref in refs:
            with self.subTest(state1=(T1, P1), state2=(T2, P2)):
                dh = enthalpy_molar(T2, P2 * 1e5, eos) - enthalpy_molar(T1, P1 * 1e5, eos)
                du = internal_energy_molar(T2, P2 * 1e5, eos) - internal_energy_molar(
                    T1, P1 * 1e5, eos
                )
                self.assertLess(abs(dh - dh_ref), 25.0)
                self.assertLess(abs(du - du_ref), 10.0)

    def _cavern_run(self, P1: float, mode: str, bhp_bar: float, thermal_model: str) -> float:
        props = build_teesside_cavern()
        props.peneloux = True
        props.volume_translation = -4.0e-6
        props.thermal_model = thermal_model
        model = UHSReservoir0D(props, tau_days=1.0e9)  # adiabatic walls
        T1 = 293.15
        state = ReservoirState(P1, T1, model.mass_from_state(P1, T1))
        sched = [SchedulePeriod("p", mode, 1.0, bar_to_pa(bhp_bar), 15.0)]
        res = model.run_schedule(state, sched, dt_s=6 * 3600.0)
        self.assertEqual(res.metrics["converged"], 1.0)
        i = next(k for k, p in enumerate(res.pressures) if abs(p - bar_to_pa(bhp_bar)) < 1e3)
        return res.temperatures[i] - T1

    def test_adiabatic_cavern_fill_vs_coolprop(self) -> None:
        # CoolProp energy-balance integration: 40 -> 120 bar with H2 injected at 293.15 K: +70.83 K.
        dT = self._cavern_run(bar_to_pa(40.0), "inject", 120.0, "energy_balance")
        self.assertLess(abs(dT / 70.83 - 1.0), 0.02)
        # The legacy lumped model treats filling as isentropic compression (+110 K).
        self.assertGreater(self._cavern_run(bar_to_pa(40.0), "inject", 120.0, "lumped"), 100.0)

    def test_adiabatic_cavern_discharge_vs_coolprop(self) -> None:
        # Adiabatic discharge leaves the remaining gas on an isentrope: CoolProp -81.16 K.
        dT = self._cavern_run(bar_to_pa(120.0), "withdraw", 40.0, "energy_balance")
        self.assertLess(abs(dT / -81.16 - 1.0), 0.03)

    def test_sorbent_store_keeps_lumped_model(self) -> None:
        from h2multisim.thermo import IsothermModel

        props = build_teesside_cavern()
        props.m_sorbent = 1000.0
        props.isotherm = IsothermModel(model="langmuir")
        self.assertEqual(UHSReservoir0D(props).thermal_model, "lumped")
        with self.assertRaises(ValueError):
            build_teesside_cavern().__class__(
                **{**build_teesside_cavern().__dict__, "thermal_model": "x"}
            )


class TestEnergyMetrics(unittest.TestCase):
    def test_heating_values(self) -> None:
        from h2multisim.thermo import heating_values_mass

        hhv, lhv = heating_values_mass()
        self.assertAlmostEqual(lhv / 3.6e6, 33.32, places=2)  # kWh/kg
        self.assertAlmostEqual(hhv / 1e6, 141.8, places=1)  # MJ/kg
        # Inert N2 dilutes the mass-specific heating value.
        self.assertLess(heating_values_mass({"H2": 0.8, "N2": 0.2})[1], 0.3 * lhv)

    def test_0d_energy_metrics_consistent_with_mass(self) -> None:
        props = build_okoroafor_reservoir()
        model = UHSReservoir0D(props)
        P0, T0 = props.reference_pressure, props.reference_temperature
        sched = [
            SchedulePeriod("inj", "inject", 20.0, bar_to_pa(100.0), 35.0),
            SchedulePeriod("wdr", "withdraw", 20.0, bar_to_pa(70.0), 35.0),
        ]
        m = model.run_schedule(
            ReservoirState(P0, T0, model.mass_from_state(P0, T0)), sched, 86400.0
        ).metrics
        from h2multisim.thermo import heating_values_mass

        lhv = heating_values_mass()[1] / 3.6e12  # GWh per kg (pure H2)
        self.assertAlmostEqual(
            m["total_withdrawn_energy_GWh_LHV"], m["total_withdrawn_kg"] * lhv, delta=1e-3
        )
        self.assertAlmostEqual(
            m["working_gas_energy_GWh_LHV"], m["working_gas_kg"] * lhv, delta=1e-3
        )

    def test_0d_withdrawn_energy_reflects_cushion_dilution(self) -> None:
        props = build_okoroafor_reservoir()
        props.cushion_gas = "N2"
        props.initial_composition = {"H2": 0.3, "N2": 0.7}
        model = UHSReservoir0D(props)
        P0, T0 = props.reference_pressure, props.reference_temperature
        st = ReservoirState(
            P0, T0, model.mass_from_state(P0, T0), composition={"H2": 0.3, "N2": 0.7}
        )
        sched = [SchedulePeriod("wdr", "withdraw", 10.0, bar_to_pa(60.0), 35.0)]
        m = model.run_schedule(st, sched, 86400.0).metrics
        lhv_h2 = 33.32e-6  # GWh/kg, pure H2 LHV
        # Produced gas is mostly N2 by mass: energy far below mass × H2 LHV.
        self.assertGreater(m["total_withdrawn_kg"], 0.0)
        self.assertLess(m["total_withdrawn_energy_GWh_LHV"], 0.1 * m["total_withdrawn_kg"] * lhv_h2)


class TestCushionGas2D(unittest.TestCase):
    def _run(self):
        props, wells, sched = scenario_single_well(nx=7, ny=7)
        props.cushion_gas = "N2"
        props.initial_h2_mole_fraction = 0.0
        sim = UHSAreal2D(props, wells)
        res = sim.run_schedule(sched, dt_days=2.0, snap_days=[30.0])
        return props, wells, sim, res

    def test_mass_conserved_and_converged(self) -> None:
        _, _, _, res = self._run()
        self.assertEqual(res.metrics["converged"], 1.0)
        self.assertLess(res.metrics["mass_balance_error"], 1e-6)

    def test_h2_species_conserved(self) -> None:
        from h2multisim.thermo import GAS_COMPONENTS, M_H2

        props, wells, sched = scenario_single_well(nx=7, ny=7)
        props.cushion_gas = "N2"
        props.initial_h2_mole_fraction = 0.0
        sim = UHSAreal2D(props, wells)
        res = sim.run_schedule(sched, dt_days=2.0, snap_days=[70.0])
        m_cg = GAS_COMPONENTS["N2"].molar_mass

        def w_h2(y: float) -> float:
            return y * M_H2 / (y * M_H2 + (1.0 - y) * m_cg)

        # Net H2 through the well: injected pure H2 minus H2 produced at the well-cell composition.
        dt = 2.0 * 86400.0
        h2_net = sum(
            q * dt if q > 0.0 else q * dt * w_h2(y)
            for q, y in zip(res.per_well_q["W1"], res.per_well_h2_mole_fraction["W1"])
        )
        snap = res.snapshots[-1]
        h2_in_place = sum(
            sim._rho(P, T, y) * sim.Vp * sim._pv_mult(P) * w_h2(y)
            for P, T, y in zip(snap.P_Pa, snap.T_K, snap.h2_mole_fraction)
        )
        self.assertGreater(h2_net, 0.0)
        self.assertLess(abs(h2_in_place - h2_net) / h2_net, 1e-4)

    def test_plume_symmetric_and_purity_reported(self) -> None:
        _, wells, _, res = self._run()
        snap = res.snapshots[0]
        ny = 7
        row = [snap.h2_mole_fraction[3 * ny + j] for j in range(ny)]
        self.assertGreater(row[3], 0.5)  # injector cell is H2-rich
        self.assertAlmostEqual(row[2], row[4], places=6)  # symmetric plume
        self.assertLess(row[0], 0.05)  # far field still cushion gas
        m = res.metrics
        self.assertGreater(m["min_withdrawn_h2_purity_mol_pct"], 0.0)
        self.assertLess(m["min_withdrawn_h2_purity_mol_pct"], 100.0)
        # Produced energy reflects N2 dilution: below withdrawn mass × H2 LHV.
        self.assertLess(
            m["total_withdrawn_energy_GWh_LHV"], m["total_withdrawn_t"] * 1000.0 * 33.32e-6
        )

    def test_single_gas_path_unchanged_and_deck_warning(self) -> None:
        props, wells, sched = scenario_single_well(nx=4, ny=4)
        self.assertFalse(UHSAreal2D(props, wells)._is_multicomponent)
        props.cushion_gas = "CH4"
        props.initial_h2_mole_fraction = 0.2
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "c.DATA")
            export_eclipse_2d(props, wells, sched, path)
            with open(path, encoding="utf-8") as fh:
                self.assertIn("WARNING: H2MultiSim models a CH4 cushion gas", fh.read())
        with self.assertRaises(ValueError):
            scenario_single_well(nx=4, ny=4)[0].__class__(
                **{**scenario_single_well(nx=4, ny=4)[0].__dict__, "cushion_gas": "He"}
            )


class TestCushionGas3D(unittest.TestCase):
    def test_conservation_and_convergence(self) -> None:
        from h2multisim.thermo import GAS_COMPONENTS, M_H2, heating_values_mass

        props, wells, sched = scenario_gravity_override(nx=3, ny=3, nz=4)
        props.cushion_gas = "N2"
        props.initial_h2_mole_fraction = 0.0
        sim = UHSFieldScale3D(props, wells)
        res = sim.run_schedule(sched, dt_days=5.0, snap_times_d=[70.0])
        m = res.metrics
        self.assertEqual(m["converged"], 1.0)
        self.assertLess(m["mass_balance_error"], 1e-6)
        # Only H2 carries a heating value, so withdrawn LHV energy gives produced H2 mass.
        produced_h2 = m["total_withdrawn_energy_GWh_LHV"] * 3.6e12 / heating_values_mass()[1]
        m_cg = GAS_COMPONENTS["N2"].molar_mass
        snap = res.snapshots[-1]
        h2_in_place = sum(
            sim._rho(P, T, y)
            * sim.Vp
            * sim._pv_mult(P)
            * (y * M_H2 / (y * M_H2 + (1.0 - y) * m_cg))
            for P, T, y in zip(snap.P_Pa, snap.T_K, snap.h2_mole_fraction)
        )
        net = m["total_injected_t"] * 1000.0 - produced_h2
        self.assertLess(abs(h2_in_place - net) / net, 1e-4)

    def test_gravity_override_lock_exchange(self) -> None:
        from h2multisim.models.three_d import Reservoir3DProps, Schedule3D, WellSpec3D

        nx, nz = 4, 4
        props = Reservoir3DProps(
            Lx=40.0, Ly=10.0, Lz=20.0, phi=0.2, k=mD_to_m2(1000.0), Ct=0.0,
            P_init=80e5, T_init=320.0, nx=nx, ny=1, nz=nz, tau_days=0.0,
            cushion_gas="N2", initial_h2_mole_fraction=0.0,
            dispersivity_m=0.0, molecular_diffusivity_m2_s=0.0,
        )  # fmt: skip
        sim = UHSFieldScale3D(props, [WellSpec3D(0, 0, 0, 0, name="X")])
        # H2 fills the left half, N2 the right half (no dispersion).
        sim.y = [[[1.0 if i < nx // 2 else 0.0] * nz] for i in range(nx)]
        P = sim._hydrostatic_init([[[80e5] * nz] for _ in range(nx)], nx, 1, nz)
        T = [[[320.0] * nz] for _ in range(nx)]
        ctrl = Schedule3D("s", "shut_in", 1.0, 0.0, 0.0)
        for _ in range(8):
            P, T, _ = sim._solve_pressure(P, T, ctrl, 86400.0 / 4)
            sim.y = sim._last_y_new
        right = nx // 2
        # Light H2 over-rides along the top; dense N2 slumps along the bottom.
        self.assertGreater(sim.y[right][0][nz - 1], sim.y[right][0][0] + 0.1)
        self.assertLess(sim.y[right - 1][0][0], sim.y[right - 1][0][nz - 1] - 0.1)

    def test_single_gas_unchanged(self) -> None:
        props, wells, sched = scenario_gravity_override(nx=3, ny=3, nz=3)
        self.assertFalse(UHSFieldScale3D(props, wells)._is_multicomponent)


class TestModeValidation(unittest.TestCase):
    def test_unknown_modes_rejected(self) -> None:
        with self.assertRaises(ValueError):
            SchedulePeriod("p", "injection", 1.0, 1e7, 1.0)
        with self.assertRaises(ValueError):
            WellControl("p", "Inject", 1.0, 1e7, 1.0)
        with self.assertRaises(ValueError):
            Schedule2D("p", "withdrawal", 1.0, 1e7, 1.0)
        with self.assertRaises(ValueError):
            Schedule2D("p", "inject", 1.0, 1e7, 1.0, well_modes={"W1": "shutin"})
        with self.assertRaises(ValueError):
            Schedule3D("p", "produce", 1.0, 1e7, 1.0)

    def test_known_modes_accepted(self) -> None:
        for mode in ("inject", "withdraw", "shut_in"):
            SchedulePeriod("p", mode, 1.0, 1e7, 1.0)
            WellControl("p", mode, 1.0, 1e7, 1.0)


class TestInputValidation(unittest.TestCase):
    def test_invalid_reservoir_inputs_rejected(self) -> None:
        from dataclasses import replace

        props0 = build_okoroafor_reservoir()
        for field, bad in (("porosity", 0.0), ("porosity", 1.5), ("permeability", -1e-15),
                           ("reference_temperature", 0.0), ("total_compressibility", -1e-9)):  # fmt: skip
            with self.subTest(field=field, value=bad), self.assertRaises(ValueError):
                replace(props0, **{field: bad})

        props1, _ = okoroafor_radial()
        for field, bad in (("phi", 0.0), ("k", 0.0), ("r_ext", 0.05), ("outer_bc", "open")):
            with self.subTest(field=field, value=bad), self.assertRaises(ValueError):
                replace(props1, **{field: bad})

        props2, _, _ = scenario_single_well(nx=4, ny=4)
        with self.assertRaises(ValueError):
            replace(props2, nx=0)
        with self.assertRaises(ValueError):
            replace(props2, P_init=float("nan"))

    def test_invalid_schedule_inputs_rejected(self) -> None:
        with self.assertRaises(ValueError):
            SchedulePeriod("p", "inject", 0.0, 1e7, 1.0)
        with self.assertRaises(ValueError):
            WellControl("p", "withdraw", 1.0, 1e7, -1.0)
        with self.assertRaises(ValueError):
            Schedule3D("p", "inject", -2.0, 1e7, 1.0)

    def test_zero_relaxation_time_allowed(self) -> None:
        from dataclasses import replace

        props1, _ = okoroafor_radial()
        self.assertEqual(replace(props1, tau_days=0.0).tau_days, 0.0)


class TestMonteCarloWaterSaturation(unittest.TestCase):
    def test_mc_honours_swr(self) -> None:
        def mc_initial_mass(swr: float) -> float:
            props = build_okoroafor_reservoir()
            props.swr = swr
            cfg = MonteCarloConfig(
                1, props.porosity, 0.0, props.thickness, 0.0, -13.0, 0.0,
                props.reference_temperature, 0.0, props.area, 0.0,
            )  # fmt: skip
            sched = [SchedulePeriod("s", "shut_in", 30.0, 0.0, 0.0)]
            return run_monte_carlo(props, sched, cfg, dt_days=30.0)["initial_mass_kg"]["mean"]

        self.assertAlmostEqual(mc_initial_mass(0.3) / mc_initial_mass(0.0), 0.7, places=6)


class TestSolverConsistency(unittest.TestCase):
    def test_1d_pore_compressibility_used(self) -> None:
        props, _ = okoroafor_radial()
        sched = [WellControl("inj", "inject", 10.0, bar_to_pa(108.0), 35.0)]
        model = UHSRadial1D(props, N=12)
        res = model.run_schedule(sched, dt_days=1.0)
        self.assertGreater(res.metrics["pore_volume_m3"], sum(model.Vp))
        self.assertLess(res.metrics["mass_balance_error"], 1e-6)

    def test_1d_const_pressure_boundary_mass_balance(self) -> None:
        props, _ = okoroafor_radial()
        props.outer_bc = "const_pressure"
        props.r_ext = 150.0
        sched = [WellControl("wdr", "withdraw", 10.0, bar_to_pa(60.0), 35.0)]
        res = UHSRadial1D(props, N=12).run_schedule(sched, dt_days=1.0)
        self.assertGreater(res.metrics["boundary_influx_kg"], 0.0)
        self.assertLess(res.metrics["mass_balance_error"], 1e-6)

    def test_2d_peaceman_equivalent_radius(self) -> None:
        props, wells, _ = scenario_single_well(nx=5, ny=5)
        sim = UHSAreal2D(props, wells)
        P = [[props.P_init] * 5 for _ in range(5)]
        T = [[props.T_init] * 5 for _ in range(5)]
        w = wells[0]
        rho = sim.eos.density(props.P_init, props.T_init)
        mu = viscosity_h2(props.T_init, props.P_init, rho=rho, eos=sim.eos)
        import math

        r_eq = 0.198 * sim.dx  # Peaceman, square cells
        wi_ref = 2.0 * math.pi * sim.k[w.i][w.j] * props.h * rho / (mu * math.log(r_eq / w.r_w))
        self.assertLess(abs(sim._well_index(w, P, T) / wi_ref - 1.0), 0.002)


class TestDeckExport(unittest.TestCase):
    def _read(self, path: str) -> str:
        with open(path, encoding="utf-8") as fh:
            return fh.read()

    def test_wconprod_rate_in_grat_item(self) -> None:
        props, sched = okoroafor_radial()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "r.DATA")
            export_eclipse_1d_radial(props, sched, path, N=10)
            text = self._read(path)
        line = text.split("WCONPROD\n", 1)[1].splitlines()[0]
        items = line.split()
        # 'W1' 'OPEN' 'BHP' 2* GRAT 2* BHP /
        self.assertEqual(items[3], "2*")
        self.assertEqual(items[5], "2*")
        self.assertAlmostEqual(float(items[6]), 55.0, places=2)
        self.assertGreater(float(items[4]), 1.0e5)  # 35 kg/s in sm³/day
        # COMPDAT: explicit CF (item 8) with the solver's equivalent radius, skin (item 11)
        comp = text.split("COMPDAT\n", 1)[1].splitlines()[0].split()
        # 'W1' I J K1 K2 'OPEN' 1* CF DIAM 1* SKIN /
        self.assertGreater(float(comp[7]), 0.0)
        self.assertEqual(comp[9], "1*")
        self.assertAlmostEqual(float(comp[10]), props.skin, places=3)

    def test_1d_deck_reproduces_radial_grid(self) -> None:
        import math

        props, sched = okoroafor_radial()
        model = UHSRadial1D(props, N=10)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "r.DATA")
            export_eclipse_1d_radial(props, sched, path, N=10)
            text = self._read(path)
        porv = [float(v) for v in text.split("PORV\n", 1)[1].split("/", 1)[0].split()]
        tranx = [float(v) for v in text.split("TRANX\n", 1)[1].split("/", 1)[0].split()]
        self.assertEqual(len(porv), 10)
        for v_deck, v_solver in zip(porv, model.Vp):
            self.assertAlmostEqual(v_deck / v_solver, 1.0, places=5)
        # TRANX [cP·rm³/day/bar] = 2π·0.008527·k[mD]·h / ln(r_c,i+1 / r_c,i); last link closed
        k_mD = props.k / 9.869233e-16
        t_ref = 2.0 * math.pi * 0.008527 * k_mD * props.h / math.log(model.r_c[1] / model.r_c[0])
        self.assertAlmostEqual(tranx[0] / t_ref, 1.0, places=4)
        self.assertEqual(tranx[-1], 0.0)

    def test_3d_layers_and_perforations_flipped_to_deck_order(self) -> None:
        props, wells, sched = scenario_gravity_override()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "g.DATA")
            export_eclipse_3d(props, wells, sched, path)
            text = self._read(path)
        comp = text.split("COMPDAT\n", 1)[1].splitlines()[0].split()
        nz = props.nz
        # Solver layers 0..1 are the BOTTOM two layers -> deck K = nz-1 .. nz.
        self.assertEqual((int(comp[3]), int(comp[4])), (nz - 1, nz))

        props_l, wells_l, sched_l = scenario_layered_reservoir()
        props_l.k_layers[0] = props_l.k_layers[0] * 3.0  # make the bottom layer unique
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "l.DATA")
            export_eclipse_3d(props_l, wells_l, sched_l, path)
            text = self._read(path)
        permx = text.split("PERMX\n", 1)[1].split("  /", 1)[0].split()
        permz = text.split("PERMZ\n", 1)[1].split("  /", 1)[0].split()
        values = [float(tok.split("*")[1]) for tok in permx]
        self.assertAlmostEqual(values[-1], 15.0, places=2)  # bottom layer written last
        self.assertEqual(permx, permz)  # isotropic, like the 3D solver

    def test_gas_water_deck_with_connate_water(self) -> None:
        props, wells, sched = scenario_single_well(nx=4, ny=4)
        props.swr = 0.25
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "s.DATA")
            export_eclipse_2d(props, wells, sched, path)
            text = self._read(path)
        # OPM Flow has no single-phase gas mode: WATER + GAS with immobile water.
        self.assertIn("\nWATER\nGAS\n", text)
        poro = float(text.split("PORO\n", 1)[1].splitlines()[0].split("*")[1].split()[0])
        self.assertAlmostEqual(poro, props.phi, places=4)
        swat = float(text.split("SWAT\n", 1)[1].splitlines()[0].split("*")[1].split()[0])
        self.assertAlmostEqual(swat, 0.25, places=6)
        # SGWFN: last row is S_g = 1 - S_wr with k_rw = 0 (connate water immobile)
        rows = text.split("SGWFN\n", 1)[1].split("/", 1)[0].strip().splitlines()[1:]
        sg_last, krg_last, krw_last, _ = (float(x) for x in rows[-1].split())
        self.assertAlmostEqual(sg_last, 0.75, places=6)
        self.assertAlmostEqual(krw_last, 0.0, places=9)
        self.assertGreater(krg_last, 0.0)


class TestConfigUnits(unittest.TestCase):
    def test_pressure_units(self) -> None:
        self.assertAlmostEqual(_pa("100 bar"), 1.0e7)
        self.assertAlmostEqual(_pa("10 MPa"), 1.0e7)
        self.assertAlmostEqual(_pa("10000 kPa"), 1.0e7)
        self.assertAlmostEqual(_pa("1e7 Pa"), 1.0e7)
        self.assertAlmostEqual(_pa("1450.377 psi"), 1.0e7, delta=10.0)
        with self.assertRaises(ValueError):
            _pa("100 atm")


class TestPCGSolver(unittest.TestCase):
    def test_pcg_solves_small_spd_system(self) -> None:
        from h2multisim.models._common import pcg_solve

        # 1D Laplacian with a positive diagonal shift: [3 -1 0; -1 3 -1; 0 -1 3] x = b
        diag = [3.0, 3.0, 3.0]
        nbrs = [[(1, 1.0)], [(0, 1.0), (2, 1.0)], [(1, 1.0)]]
        x_true = [1.0, 2.0, 3.0]
        b = [3.0 * 1.0 - 2.0, 3.0 * 2.0 - 1.0 - 3.0, 3.0 * 3.0 - 2.0]
        x, ok = pcg_solve(diag, nbrs, b, [0.0, 0.0, 0.0], tol=1e-12, max_iter=50)
        self.assertTrue(ok)
        for xi, ti in zip(x, x_true):
            self.assertAlmostEqual(xi, ti, places=10)

    def test_pcg_matches_sor_in_3d_with_gravity(self) -> None:
        results = []
        for solver in ("pcg", "sor"):
            props, wells, sched = scenario_gravity_override(nx=3, ny=3, nz=4)
            res = UHSFieldScale3D(props, wells, linear_solver=solver).run_schedule(
                sched[:1], dt_days=5.0
            )
            results.append((res.P_avg_bar[-1], res.mass_t[-1]))
        self.assertAlmostEqual(results[0][0], results[1][0], places=6)
        self.assertAlmostEqual(results[0][1], results[1][1], delta=1e-3)

    def test_unknown_linear_solver_rejected(self) -> None:
        props, wells, _ = scenario_single_well(nx=3, ny=3)
        with self.assertRaises(ValueError):
            UHSAreal2D(props, wells, linear_solver="gmres")


class TestConvergenceFlags(unittest.TestCase):
    def test_3d_reports_solver_health(self) -> None:
        props, wells, sched = scenario_gravity_override(nx=3, ny=3, nz=3)
        res = UHSFieldScale3D(props, wells).run_schedule(sched[:1], dt_days=5.0)
        for key in ("sor_unconverged_solves", "active_set_unconverged_steps", "converged"):
            self.assertIn(key, res.metrics)


if __name__ == "__main__":
    unittest.main()
