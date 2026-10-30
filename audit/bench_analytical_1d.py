"""
audit/bench_analytical_1d.py
============================
Industry Benchmark (Tier 2): Semi-Analytical Real-Gas Pseudopressure Solution.

Validates H2MultiSim's 1D radial finite-volume solver (UHSRadial1D) against the
closed-form Al-Hussainy, Ramey & Crawford (1966) real-gas pseudopressure line-source
solution used by commercial reservoir simulators (ECLIPSE, CMG, OPM Flow).

Governing Equations:
  Real-gas pseudopressure:
    m(p) = 2 ∫ [p / (μ(p)·Z(p))] dp

  Infinite-acting radial diffusivity equation:
    (1/r)·∂/∂r [r · ∂m/∂r] = (1/η) · ∂m/∂t
    where hydraulic diffusivity η = k / [φ · μ_i · c_ti]

  Line-source exponential integral solution (production rate q [kg/s]):
    m(r, t) = m(p_i) - [q·R·T / (2·π·k·h·M)] · E1(r² / [4·η·t])
    m(p_wf, t) = m(p_i) - [q·R·T / (2·π·k·h·M)] · [E1(r_w² / [4·η·t]) + 2·s]

Zero external runtime dependencies.
"""

from __future__ import annotations

import bisect
import math
import pathlib
import sys

# Ensure src/ is on sys.path
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from h2multisim.models.one_d import RadialReservoirProps, UHSRadial1D, WellControl
from h2multisim.thermo import M_H2, R_UNIVERSAL, PengRobinsonH2, viscosity_h2
from h2multisim.units import bar_to_pa, mD_to_m2, pa_to_bar


def exp1_pure(x: float) -> float:
    """
    Exponential integral E1(x) = ∫[x, ∞] (e^-u / u) du for x > 0 to machine precision.
    Uses series expansion for x < 1.0 and continued fraction for x >= 1.0.
    """
    if x <= 0.0:
        return float("inf")
    if x < 1.0:
        euler_gamma = 0.5772156649015328606
        s = 0.0
        term = 1.0
        for n in range(1, 40):
            term *= -x / n
            s += term / n
        return -euler_gamma - math.log(x) - s
    else:
        tiny = 1e-30
        b = x + 1.0
        c = 1.0 / tiny
        d = 1.0 / b
        h = d
        for i in range(1, 100):
            a = -i * i
            b += 2.0
            d = 1.0 / max(b + a * d, tiny)
            c = max(b + a / c, tiny)
            delta = c * d
            h *= delta
            if abs(delta - 1.0) < 1e-15:
                break
        return math.exp(-x) * h


class RealGasPseudopressureAnalytical:
    """
    Al-Hussainy real-gas pseudopressure table and line-source solver.
    """

    def __init__(
        self,
        props: RadialReservoirProps,
        eos: PengRobinsonH2,
        P_max_eval: float,
        nP: int = 4000,
    ) -> None:
        self.props = props
        self.eos = eos
        self.T = props.T_init

        # Integrate m(P) = 2 ∫ [p / (μ·Z)] dp from 0 to P_max_eval using Simpson quadrature
        self.P_grid = [P_max_eval * i / nP for i in range(nP + 1)]
        self.m_table = [0.0]

        def integrand(p_val: float) -> float:
            p_safe = max(p_val, 1.0)
            z = self.eos.z_factor(p_safe, self.T)
            rho = self.eos.density(p_safe, self.T)
            mu = viscosity_h2(self.T, p_safe, rho=rho, eos=self.eos)
            return 2.0 * p_safe / (mu * z)

        for i in range(nP):
            p0, p1 = self.P_grid[i], self.P_grid[i + 1]
            pm = 0.5 * (p0 + p1)
            f0 = integrand(p0)
            f1 = integrand(p1)
            fm = integrand(pm)
            dm = (p1 - p0) * (f0 + 4.0 * fm + f1) / 6.0
            self.m_table.append(self.m_table[-1] + dm)

        # Initial hydraulic diffusivity
        self.m_init = self.p_to_m(props.P_init)
        cg_init = self.eos.isothermal_compressibility(props.P_init, self.T)
        ct_init = cg_init + props.Ct
        mu_init = viscosity_h2(
            self.T,
            props.P_init,
            rho=self.eos.density(props.P_init, self.T),
            eos=self.eos,
        )
        self.eta = props.k / (props.phi * mu_init * ct_init)

    def p_to_m(self, p_val: float) -> float:
        idx = bisect.bisect_left(self.P_grid, p_val)
        if idx <= 0:
            return self.m_table[0]
        if idx >= len(self.P_grid):
            return self.m_table[-1]
        frac = (p_val - self.P_grid[idx - 1]) / (self.P_grid[idx] - self.P_grid[idx - 1])
        return self.m_table[idx - 1] + frac * (self.m_table[idx] - self.m_table[idx - 1])

    def m_to_p(self, m_val: float) -> float:
        idx = bisect.bisect_left(self.m_table, m_val)
        if idx <= 0:
            return self.P_grid[0]
        if idx >= len(self.m_table):
            return self.P_grid[-1]
        frac = (m_val - self.m_table[idx - 1]) / (self.m_table[idx] - self.m_table[idx - 1])
        return self.P_grid[idx - 1] + frac * (self.P_grid[idx] - self.P_grid[idx - 1])

    def pressure_at_r(self, r: float, t_s: float, q_prod_kg_s: float) -> float:
        """Pressure [Pa] at radius r [m] and time t_s [s] for constant production rate q_prod_kg_s."""
        if t_s <= 0.0:
            return self.props.P_init
        x = r**2 / (4.0 * self.eta * t_s)
        const = (
            q_prod_kg_s
            * R_UNIVERSAL
            * self.T
            / (2.0 * math.pi * self.props.k * self.props.h * M_H2)
        )
        m_sol = self.m_init - const * exp1_pure(x)
        return self.m_to_p(m_sol)

    def flowing_bhp(self, t_s: float, q_prod_kg_s: float) -> float:
        """Flowing wellbore pressure P_wf [Pa] including wellbore radius and skin factor."""
        if t_s <= 0.0:
            return self.props.P_init
        x_w = self.props.r_well**2 / (4.0 * self.eta * t_s)
        const = (
            q_prod_kg_s
            * R_UNIVERSAL
            * self.T
            / (2.0 * math.pi * self.props.k * self.props.h * M_H2)
        )
        m_wf = self.m_init - const * (exp1_pure(x_w) + 2.0 * self.props.skin)
        return self.m_to_p(m_wf)


def run_analytical_1d_benchmark(verbose: bool = True) -> tuple[float, float]:
    """
    Run 1D transient withdrawal benchmark against Al-Hussainy analytical solution.
    Returns (max_bhp_err_pct, max_profile_err_pct).
    """
    P_init = bar_to_pa(80.0)
    T_init = 316.0
    q_rate = 5.0  # kg/s (~2.4 MMSCFD pure H2)
    duration_days = 5.0

    props = RadialReservoirProps(
        r_well=0.10,
        r_ext=2000.0,
        h=50.0,
        phi=0.20,
        k=mD_to_m2(100.0),
        Ct=5.0e-10,
        P_init=P_init,
        T_init=T_init,
        skin=0.0,
        outer_bc="no_flow",
        tau_days=0.0,
        peneloux=True,
    )
    eos = PengRobinsonH2(volume_translation=-4.0e-6)

    # 1. H2MultiSim Numerical Simulation (60 logarithmic radial cells)
    sim = UHSRadial1D(props, N=60, eos=eos)
    snap_days = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
    ctrl = [WellControl("withdraw", "withdraw", duration_days, bar_to_pa(1.0), q_rate)]
    result = sim.run_schedule(ctrl, dt_days=0.05, snap_days=snap_days)

    # 2. Analytical Al-Hussainy Pseudopressure Baseline
    solver = RealGasPseudopressureAnalytical(props, eos, P_max_eval=P_init * 1.5)

    if verbose:
        print("=" * 80)
        print("H2MultiSim vs Al-Hussainy (1966) Analytical Pseudopressure Benchmark")
        print("=" * 80)
        print(
            f"Reservoir Conditions : {pa_to_bar(P_init):.1f} bar, {T_init:.1f} K, k=100 mD, h=50 m, phi=0.20"
        )
        print(f"Hydrogen Production  : {q_rate:.1f} kg/s constant mass withdrawal")
        print("-" * 80)
        print(
            f"{'Time [d]':>10} | {'Sim BHP [bar]':>15} | {'Analytic BHP [bar]':>20} | {'Diff [bar]':>12} | {'Error [%]':>10}"
        )
        print("-" * 80)

    max_bhp_err = 0.0
    for target_d in snap_days:
        target_s = target_d * 86400.0
        # Find closest simulation time
        idx = min(range(len(result.times_s)), key=lambda i: abs(result.times_s[i] - target_s))
        p_sim = result.BHP_Pa[idx]
        p_ana = solver.flowing_bhp(target_s, q_rate)
        diff_bar = pa_to_bar(p_sim - p_ana)
        err_pct = abs(p_sim - p_ana) / p_ana * 100.0
        max_bhp_err = max(max_bhp_err, err_pct)
        if verbose:
            print(
                f"{target_d:10.2f} | {pa_to_bar(p_sim):15.3f} | {pa_to_bar(p_ana):20.3f} | {diff_bar:12.3f} | {err_pct:9.3f}%"
            )

    if verbose:
        print("-" * 80)
        print("\nRadial Pressure Profile at t = 5.0 days:")
        print("-" * 80)
        print(
            f"{'Radius [m]':>10} | {'Sim P [bar]':>15} | {'Analytic P [bar]':>20} | {'Diff [bar]':>12} | {'Error [%]':>10}"
        )
        print("-" * 80)

    snap5 = result.snapshots[-1]
    t5_s = 5.0 * 86400.0
    max_prof_err = 0.0
    # Evaluate at select radial cells across 3 orders of magnitude
    for i in (0, 2, 5, 10, 20, 30, 40, 50, 58):
        r = sim.r_c[i]
        p_sim = snap5.P_Pa[i]
        p_ana = solver.pressure_at_r(r, t5_s, q_rate)
        diff_bar = pa_to_bar(p_sim - p_ana)
        err_pct = abs(p_sim - p_ana) / p_ana * 100.0
        max_prof_err = max(max_prof_err, err_pct)
        if verbose:
            print(
                f"{r:10.2f} | {pa_to_bar(p_sim):15.3f} | {pa_to_bar(p_ana):20.3f} | {diff_bar:12.3f} | {err_pct:9.3f}%"
            )

    if verbose:
        print("=" * 80)
        print(f"Max Flowing BHP Error    : {max_bhp_err:.3f} % (tolerance < 0.50 %)")
        print(f"Max Spatial Profile Error: {max_prof_err:.3f} % (tolerance < 0.25 %)")
        print(f"Mass Balance Error       : {result.metrics['mass_balance_error']:.2e}")
        print("=" * 80 + "\n")

    return max_bhp_err, max_prof_err


if __name__ == "__main__":
    max_bhp_err, max_prof_err = run_analytical_1d_benchmark(verbose=True)
    if max_bhp_err > 0.5 or max_prof_err > 0.25:
        sys.exit(1)
    sys.exit(0)
