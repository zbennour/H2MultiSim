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
    SchedulePeriod,
    build_okoroafor_reservoir,
    run_monte_carlo,
)
from h2multisim.thermo import (
    PengRobinsonH2,
    PengRobinsonMixture,
    RelativePermeabilityModel,
    adiabatic_temperature,
    cp_real_h2,
    gamma_real_h2,
    joule_thomson_coefficient,
    viscosity_h2,
)
from h2multisim.units import bar_to_pa


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

    def test_brooks_corey_capillary_round_trip(self) -> None:
        rp = RelativePermeabilityModel(
            model="brooks_corey", swr=0.25, sgc=0.05, brooks_corey_lambda=2.0, entry_pressure_Pa=8e3
        )
        for sg in (0.2, 0.45, 0.7):
            with self.subTest(sg=sg):
                pc = rp.capillary_pressure_Pa(sg)
                self.assertAlmostEqual(rp.gas_saturation_from_capillary_pressure(pc), sg, places=9)


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
        # COMPDAT carries the skin factor (item 11)
        comp = text.split("COMPDAT\n", 1)[1].splitlines()[0].split()
        # 'W1' I J K1 K2 'OPEN' 2* DIAM 1* SKIN /
        self.assertEqual(comp[8], "1*")
        self.assertAlmostEqual(float(comp[9]), props.skin, places=3)

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

    def test_poro_includes_immobile_water(self) -> None:
        props, wells, sched = scenario_single_well(nx=4, ny=4)
        props.swr = 0.25
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "s.DATA")
            export_eclipse_2d(props, wells, sched, path)
            text = self._read(path)
        poro = float(text.split("PORO\n", 1)[1].splitlines()[0].split("*")[1].split()[0])
        self.assertAlmostEqual(poro, props.phi * 0.75, places=4)


class TestConfigUnits(unittest.TestCase):
    def test_pressure_units(self) -> None:
        self.assertAlmostEqual(_pa("100 bar"), 1.0e7)
        self.assertAlmostEqual(_pa("10 MPa"), 1.0e7)
        self.assertAlmostEqual(_pa("10000 kPa"), 1.0e7)
        self.assertAlmostEqual(_pa("1e7 Pa"), 1.0e7)
        self.assertAlmostEqual(_pa("1450.377 psi"), 1.0e7, delta=10.0)
        with self.assertRaises(ValueError):
            _pa("100 atm")


class TestConvergenceFlags(unittest.TestCase):
    def test_3d_reports_solver_health(self) -> None:
        props, wells, sched = scenario_gravity_override(nx=3, ny=3, nz=3)
        res = UHSFieldScale3D(props, wells).run_schedule(sched[:1], dt_days=5.0)
        for key in ("sor_unconverged_solves", "active_set_unconverged_steps", "converged"):
            self.assertIn(key, res.metrics)


if __name__ == "__main__":
    unittest.main()
