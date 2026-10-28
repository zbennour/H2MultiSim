"""
h2multisim.models.one_d
=======================
One-dimensional (1D) radial finite-volume reservoir model for near-wellbore
dynamics and transient pressure analysis in underground hydrogen storage (UHS).
"""

from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass, field

from h2multisim.thermo import (
    PengRobinsonH2,
    coupled_temperature_update,
    viscosity_h2,
)
from h2multisim.units import (
    K_to_C,
    bar_to_pa,
    kg_to_tonnes,
    mD_to_m2,
    pa_to_bar,
    pa_to_psi,
)


@dataclass
class RadialReservoirProps:
    """Reservoir properties for the 1D radial model."""

    r_well: float  # m wellbore radius
    r_ext: float  # m outer boundary radius
    h: float  # m net pay thickness
    phi: float  # fraction porosity
    k: float  # m² permeability
    Ct: float  # Pa⁻¹ total compressibility
    P_init: float  # Pa initial reservoir pressure
    T_init: float  # K initial / geothermal temperature
    skin: float = 0.0  # mechanical skin factor
    outer_bc: str = "no_flow"  # 'no_flow' | 'const_pressure'
    tau_days: float = 30.0  # days thermal relaxation constant


ReservoirProperties1D = RadialReservoirProps


@dataclass
class WellControl:
    """Operational period control for 1D radial model."""

    name: str
    mode: str  # 'inject' | 'withdraw' | 'shut_in'
    duration_d: float  # days
    bhp_target: float  # Pa
    q_max: float  # kg/s


SchedulePeriod1D = WellControl


@dataclass
class RadialSnapshot:
    """Radial pressure and temperature profile at one instant."""

    time_s: float
    r_m: list[float]
    P_Pa: list[float]
    T_K: list[float]


@dataclass
class RadialResult:
    """Output from 1D radial simulation."""

    times_s: list[float] = field(default_factory=list)
    BHP_Pa: list[float] = field(default_factory=list)
    q_kg_s: list[float] = field(default_factory=list)
    cum_inj_kg: list[float] = field(default_factory=list)
    cum_wdr_kg: list[float] = field(default_factory=list)
    P_avg_Pa: list[float] = field(default_factory=list)
    T_avg_K: list[float] = field(default_factory=list)
    snapshots: list[RadialSnapshot] = field(default_factory=list)
    metrics: dict[str, float] = field(default_factory=dict)


SimulationResult1D = RadialResult


def build_log_grid(r_w: float, r_e: float, N: int) -> tuple[list[float], list[float]]:
    """Logarithmic radial grid: N cell-centres and N+1 interfaces between r_w and r_e."""
    if N < 3:
        raise ValueError("N must be >= 3")
    ratio = r_e / r_w
    interfaces = [r_w * ratio ** (i / N) for i in range(N + 1)]
    centres = [math.sqrt(interfaces[i] * interfaces[i + 1]) for i in range(N)]
    return centres, interfaces


def thomas_solve(
    a: list[float],
    b: list[float],
    c: list[float],
    d: list[float],
) -> list[float]:
    """
    Solve tridiagonal system a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]
    using the Thomas algorithm (O(N) operations).
    """
    N = len(d)
    c_ = [0.0] * N
    d_ = [0.0] * N
    c_[0] = c[0] / b[0]
    d_[0] = d[0] / b[0]
    for i in range(1, N):
        denom = b[i] - a[i] * c_[i - 1]
        if abs(denom) < 1e-30:
            denom = 1e-30
        c_[i] = c[i] / denom
        d_[i] = (d[i] - a[i] * d_[i - 1]) / denom

    x = [0.0] * N
    x[N - 1] = d_[N - 1]
    for i in range(N - 2, -1, -1):
        x[i] = d_[i] - c_[i] * x[i + 1]
    return x


class UHSRadial1D:
    """
    Implicit 1D radial finite-volume UHS simulator.
    Unconditionally stable backward-Euler solver with Thomas algorithm.
    """

    def __init__(
        self,
        props: RadialReservoirProps,
        N: int = 40,
        eos: PengRobinsonH2 | None = None,
    ) -> None:
        self.props = props
        self.N = N
        self.eos = eos or PengRobinsonH2()
        self.tau_s = props.tau_days * 86400.0

        self.r_c, self.r_i = build_log_grid(props.r_well, props.r_ext, N)
        self.V_cell = [
            math.pi * (self.r_i[i + 1] ** 2 - self.r_i[i] ** 2) * props.h for i in range(N)
        ]
        self.Vp = [V * props.phi for V in self.V_cell]

    def _trans(self, P: list[float], T: list[float], i: int) -> float:
        if i < 0 or i >= self.N - 1:
            return 0.0
        rho_L = self.eos.density(P[i], T[i])
        rho_R = self.eos.density(P[i + 1], T[i + 1])
        mu_L = viscosity_h2(T[i], P[i])
        mu_R = viscosity_h2(T[i + 1], P[i + 1])
        if (rho_L + rho_R) < 1e-20 or (mu_L + mu_R) < 1e-30:
            return 0.0
        rho_h = 2.0 * rho_L * rho_R / (rho_L + rho_R)
        mu_h = 2.0 * mu_L * mu_R / (mu_L + mu_R)
        ln_r = math.log(self.r_i[i + 2] / self.r_i[i + 1])
        if ln_r < 1e-12:
            return 0.0
        return 2.0 * math.pi * self.props.k * self.props.h * rho_h / (mu_h * ln_r)

    def _well_index(self, T_w: float, P_w: float) -> float:
        delta_r0 = self.r_i[1] - self.r_i[0]
        r_eq = 0.2 * delta_r0
        rho_w = self.eos.density(P_w, T_w)
        mu_w = viscosity_h2(T_w, P_w)
        ln_r = math.log(max(r_eq, self.props.r_well * 1.001) / self.props.r_well)
        ln_r = max(ln_r + self.props.skin, 0.1)
        return 2.0 * math.pi * self.props.k * self.props.h * rho_w / (mu_w * ln_r)

    def _well_rate(self, P: list[float], T: list[float], ctrl: WellControl) -> tuple[float, float]:
        if ctrl.mode == "shut_in":
            return 0.0, P[0]

        WI = self._well_index(T[0], P[0])
        bhp = ctrl.bhp_target
        qmax = ctrl.q_max

        if ctrl.mode == "inject":
            dp = max(bhp - P[0], 0.0)
            q = min(WI * dp, qmax)
            return max(q, 0.0), bhp

        if ctrl.mode == "withdraw":
            dp = max(P[0] - bhp, 0.0)
            q = min(WI * dp, qmax)
            return -max(q, 0.0), bhp

        return 0.0, P[0]

    def _solve_pressure(
        self,
        P: list[float],
        T: list[float],
        ctrl: WellControl,
        dt: float,
    ) -> tuple[list[float], float, float]:
        N = self.N
        phi = self.props.phi
        Ct = self.props.Ct

        T_half = [self._trans(P, T, i) for i in range(N)]
        _, BHP = self._well_rate(P, T, ctrl)

        B = []
        for i in range(N):
            rho_i = self.eos.density(P[i], T[i])
            B.append(phi * Ct * self.Vp[i] * rho_i / dt)

        a = [0.0] * N
        b = [0.0] * N
        c = [0.0] * N
        d = [0.0] * N

        for i in range(N):
            T_left = T_half[i - 1] if i > 0 else 0.0
            T_right = T_half[i] if i < N - 1 else 0.0

            if i == N - 1 and self.props.outer_bc == "const_pressure":
                T_right_outer = self._trans(P + [self.props.P_init], T + [self.props.T_init], N - 1)
            else:
                T_right_outer = 0.0

            a[i] = -T_left
            c[i] = -T_right
            b[i] = B[i] + T_left + T_right + T_right_outer
            d[i] = B[i] * P[i]

            if i == 0 and ctrl.mode != "shut_in":
                WI = self._well_index(T[0], P[0])
                b[0] += WI
                d[0] += WI * ctrl.bhp_target

            if i == N - 1 and T_right_outer > 0.0:
                d[N - 1] += T_right_outer * self.props.P_init

        P_new = thomas_solve(a, b, c, d)
        P_new = [max(p, 1.0e3) for p in P_new]

        if ctrl.mode == "inject":
            P_new[0] = min(P_new[0], ctrl.bhp_target)
        elif ctrl.mode == "withdraw":
            P_new[0] = max(P_new[0], ctrl.bhp_target)

        WI_final = self._well_index(T[0], P_new[0])
        if ctrl.mode == "inject":
            q_final = min(WI_final * max(ctrl.bhp_target - P_new[0], 0.0), ctrl.q_max)
            q_final = max(q_final, 0.0)
        elif ctrl.mode == "withdraw":
            q_final = -min(WI_final * max(P_new[0] - ctrl.bhp_target, 0.0), ctrl.q_max)
            q_final = min(q_final, 0.0)
        else:
            q_final = 0.0

        BHP_final = P_new[0]
        return P_new, q_final, BHP_final

    def _update_temperature(
        self,
        P_old: list[float],
        P_new: list[float],
        T_old: list[float],
        dt: float,
    ) -> list[float]:
        T_new = []
        for i in range(self.N):
            T_i = coupled_temperature_update(
                T_old=T_old[i],
                P_old=P_old[i],
                P_new=P_new[i],
                T_geothermal=self.props.T_init,
                dt_s=dt,
                tau_s=self.tau_s,
                eos=self.eos,
            )
            T_new.append(T_i)
        return T_new

    def run_schedule(
        self,
        schedule: list[WellControl],
        dt_days: float = 0.1,
        snap_days: list[float] | None = None,
    ) -> RadialResult:
        dt = dt_days * 86400.0
        P = [self.props.P_init] * self.N
        T = [self.props.T_init] * self.N
        result = RadialResult()

        cum_inj = 0.0
        cum_wdr = 0.0
        t = 0.0
        snap_targets = sorted(snap_days or [])
        snap_idx = 0

        for ctrl in schedule:
            n_steps = max(int(ctrl.duration_d / dt_days), 1)
            for _ in range(n_steps):
                P_new, q_well, BHP = self._solve_pressure(P, T, ctrl, dt)
                T_new = self._update_temperature(P, P_new, T, dt)
                t += dt

                if q_well >= 0.0:
                    cum_inj += q_well * dt
                else:
                    cum_wdr += (-q_well) * dt

                P = P_new
                T = T_new

                result.times_s.append(t)
                result.BHP_Pa.append(BHP)
                result.q_kg_s.append(q_well)
                result.cum_inj_kg.append(cum_inj)
                result.cum_wdr_kg.append(cum_wdr)
                result.P_avg_Pa.append(sum(P) / len(P))
                result.T_avg_K.append(sum(T) / len(T))

                t_days = t / 86400.0
                while snap_idx < len(snap_targets) and t_days >= snap_targets[snap_idx]:
                    result.snapshots.append(
                        RadialSnapshot(
                            time_s=t,
                            r_m=list(self.r_c),
                            P_Pa=list(P),
                            T_K=list(T),
                        )
                    )
                    snap_idx += 1

        total_Vp = sum(self.Vp)
        m_final = sum(self.eos.density(P[i], T[i]) * self.Vp[i] for i in range(self.N))
        result.metrics = {
            "total_injected_kg": cum_inj,
            "total_withdrawn_kg": cum_wdr,
            "net_stored_kg": cum_inj - cum_wdr,
            "final_mass_kg": m_final,
            "pore_volume_m3": total_Vp,
            "final_BHP_Pa": result.BHP_Pa[-1] if result.BHP_Pa else P[0],
            "final_avg_P_Pa": result.P_avg_Pa[-1] if result.P_avg_Pa else P[0],
        }
        return result


H2Reservoir1DModel = UHSRadial1D


def okoroafor_radial() -> tuple[RadialReservoirProps, list[WellControl]]:
    """Canonical 1D radial scenario matching Okoroafor et al. (2022)."""
    props = RadialReservoirProps(
        r_well=0.10,
        r_ext=500.0,
        h=100.0,
        phi=0.20,
        k=mD_to_m2(100.0),
        Ct=5.0e-10,
        P_init=bar_to_pa(80.0),
        T_init=316.0,
        skin=0.0,
        outer_bc="no_flow",
        tau_days=30.0,
    )
    sched = [
        WellControl("cushion_fill", "inject", 90.0, bar_to_pa(108.0), 35.0),
        WellControl("shut_in_1", "shut_in", 30.0, 0.0, 0.0),
        WellControl("withdrawal_1", "withdraw", 60.0, bar_to_pa(55.0), 35.0),
        WellControl("shut_in_2", "shut_in", 30.0, 0.0, 0.0),
    ]
    return props, sched


def scenario_1d_nearwell() -> tuple[RadialReservoirProps, list[WellControl]]:
    return okoroafor_radial()


def export_radial_timeseries_csv(result: RadialResult, filepath: str) -> None:
    """Export 1D radial time series to CSV."""
    if not result.times_s:
        return
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    with open(filepath, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "time_days",
                "bhp_bar",
                "bhp_psi",
                "q_kg_s",
                "cum_inj_tonnes",
                "cum_wdr_tonnes",
                "p_avg_bar",
                "t_avg_C",
            ]
        )
        for t, bhp, q, ci, cw, pa, ta in zip(
            result.times_s,
            result.BHP_Pa,
            result.q_kg_s,
            result.cum_inj_kg,
            result.cum_wdr_kg,
            result.P_avg_Pa,
            result.T_avg_K,
            strict=False,
        ):
            w.writerow(
                [
                    f"{t / 86400:.4f}",
                    f"{pa_to_bar(bhp):.4f}",
                    f"{pa_to_psi(bhp):.2f}",
                    f"{q:.4f}",
                    f"{kg_to_tonnes(ci):.4f}",
                    f"{kg_to_tonnes(cw):.4f}",
                    f"{pa_to_bar(pa):.4f}",
                    f"{K_to_C(ta):.2f}",
                ]
            )
