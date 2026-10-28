"""
h2multisim.models.two_d
=======================
Two-dimensional (2D) areal finite-volume simulator for underground hydrogen
storage (UHS) in porous reservoirs with multi-well interference and spatial heterogeneity.
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
)


@dataclass
class Reservoir2DProps:
    """
    Properties for the 2D areal model.
    k_field: Optional per-cell permeability list of lists [m²] (nx x ny).
    """

    Lx: float  # m domain length in x
    Ly: float  # m domain length in y
    h: float  # m net pay thickness
    phi: float  # fraction porosity
    k: float  # m² background permeability
    Ct: float  # Pa⁻¹ total compressibility
    P_init: float  # Pa initial reservoir pressure
    T_init: float  # K initial / geothermal temperature
    nx: int  # grid cells in x
    ny: int  # grid cells in y
    tau_days: float = 30.0  # days thermal relaxation constant
    k_field: list[list[float]] | None = None


ReservoirProperties2D = Reservoir2DProps


@dataclass
class WellSpec2D:
    """A single well in the 2D areal grid."""

    i: int  # x-cell index (0 to nx-1)
    j: int  # y-cell index (0 to ny-1)
    r_w: float = 0.1  # m wellbore radius
    name: str = "W"


Well2D = WellSpec2D


@dataclass
class Schedule2D:
    """Operational period for 2D areal simulation."""

    name: str
    mode: str  # 'inject' | 'withdraw' | 'shut_in'
    duration_d: float  # days
    bhp_target: float  # Pa
    q_max: float  # kg/s per well limit


SchedulePeriod2D = Schedule2D


@dataclass
class Snapshot2D:
    """Areal spatial snapshot at a given time."""

    time_d: float
    xs_m: list[float]
    ys_m: list[float]
    P_Pa: list[float]
    T_K: list[float]


@dataclass
class Result2D:
    """Simulation results from 2D areal model."""

    times_d: list[float] = field(default_factory=list)
    P_avg_bar: list[float] = field(default_factory=list)
    T_avg_C: list[float] = field(default_factory=list)
    mass_t: list[float] = field(default_factory=list)
    q_total_kg_s: list[float] = field(default_factory=list)
    per_well_q: dict[str, list[float]] = field(default_factory=dict)
    per_well_BHP: dict[str, list[float]] = field(default_factory=dict)
    snapshots: list[Snapshot2D] = field(default_factory=list)
    metrics: dict[str, float] = field(default_factory=dict)


SimulationResult2D = Result2D


class UHSAreal2D:
    """
    Implicit 2D areal finite-volume UHS simulator.
    Solves pressure implicitly using Successive Over-Relaxation (SOR) with Peaceman well model.
    """

    def __init__(
        self,
        props: Reservoir2DProps,
        wells: list[WellSpec2D],
        eos: PengRobinsonH2 | None = None,
        sor_omega: float = 1.5,
        max_iter: int = 80,
        tol_Pa: float = 10.0,
    ) -> None:
        self.props = props
        self.wells = wells
        self.eos = eos or PengRobinsonH2()
        self.omega = sor_omega
        self.max_iter = max_iter
        self.tol_Pa = tol_Pa
        self.tau_s = props.tau_days * 86400.0

        self.nx = props.nx
        self.ny = props.ny
        self.dx = props.Lx / props.nx
        self.dy = props.Ly / props.ny
        self.xs = [(i + 0.5) * self.dx for i in range(self.nx)]
        self.ys = [(j + 0.5) * self.dy for j in range(self.ny)]

        self.Vp = props.phi * self.dx * self.dy * props.h

        if props.k_field is not None:
            self.k = props.k_field
        else:
            self.k = [[props.k] * props.ny for _ in range(props.nx)]

    def _k_harm_x(self, i: int, j: int) -> float:
        kL = self.k[i][j]
        kR = self.k[i + 1][j]
        if kL + kR < 1e-40:
            return 0.0
        return 2.0 * kL * kR / (kL + kR)

    def _k_harm_y(self, i: int, j: int) -> float:
        kB = self.k[i][j]
        kT = self.k[i][j + 1]
        if kB + kT < 1e-40:
            return 0.0
        return 2.0 * kB * kT / (kB + kT)

    def _Tx(self, i: int, j: int, P: list[list[float]], T: list[list[float]]) -> float:
        kh = self._k_harm_x(i, j)
        rhoL = self.eos.density(P[i][j], T[i][j])
        rhoR = self.eos.density(P[i + 1][j], T[i + 1][j])
        muL = viscosity_h2(T[i][j], P[i][j])
        muR = viscosity_h2(T[i + 1][j], P[i + 1][j])
        rho_h = 2.0 * rhoL * rhoR / (rhoL + rhoR + 1e-30)
        mu_h = 2.0 * muL * muR / (muL + muR + 1e-30)
        return kh * self.props.h * self.dy * rho_h / (mu_h * self.dx)

    def _Ty(self, i: int, j: int, P: list[list[float]], T: list[list[float]]) -> float:
        kh = self._k_harm_y(i, j)
        rhoB = self.eos.density(P[i][j], T[i][j])
        rhoT = self.eos.density(P[i][j + 1], T[i][j + 1])
        muB = viscosity_h2(T[i][j], P[i][j])
        muT = viscosity_h2(T[i][j + 1], P[i][j + 1])
        rho_h = 2.0 * rhoB * rhoT / (rhoB + rhoT + 1e-30)
        mu_h = 2.0 * muB * muT / (muB + muT + 1e-30)
        return kh * self.props.h * self.dx * rho_h / (mu_h * self.dy)

    def _well_index(self, w: WellSpec2D, P: list[list[float]], T: list[list[float]]) -> float:
        r_eq = 0.28 * math.sqrt(self.dx * self.dy)
        ln_r = math.log(max(r_eq, w.r_w * 1.001) / w.r_w)
        ln_r = max(ln_r, 0.1)
        k_w = self.k[w.i][w.j]
        rho = self.eos.density(P[w.i][w.j], T[w.i][w.j])
        mu = viscosity_h2(T[w.i][w.j], P[w.i][w.j])
        return 2.0 * math.pi * k_w * self.props.h * rho / (mu * ln_r)

    def _solve_pressure(
        self,
        P_old: list[list[float]],
        T: list[list[float]],
        ctrl: Schedule2D,
        dt: float,
    ) -> tuple[list[list[float]], dict[str, tuple[float, float]]]:
        nx, ny = self.nx, self.ny
        phi = self.props.phi
        Ct = self.props.Ct

        B = [[0.0] * ny for _ in range(nx)]
        Tx = [[0.0] * ny for _ in range(nx - 1)]
        Ty_arr = [[0.0] * (ny - 1) for _ in range(nx)]

        for i in range(nx):
            for j in range(ny):
                rho_ij = self.eos.density(P_old[i][j], T[i][j])
                B[i][j] = phi * Ct * self.Vp * rho_ij / dt

        for i in range(nx - 1):
            for j in range(ny):
                Tx[i][j] = self._Tx(i, j, P_old, T)

        for i in range(nx):
            for j in range(ny - 1):
                Ty_arr[i][j] = self._Ty(i, j, P_old, T)

        wi_dict = {w.name: self._well_index(w, P_old, T) for w in self.wells}
        per_well: dict[str, tuple[float, float]] = {}

        P_new = [list(row) for row in P_old]

        for _iter in range(self.max_iter):
            max_change = 0.0

            for i in range(nx):
                for j in range(ny):
                    TxL = Tx[i - 1][j] if i > 0 else 0.0
                    TxR = Tx[i][j] if i < nx - 1 else 0.0
                    TyB = Ty_arr[i][j - 1] if j > 0 else 0.0
                    TyT = Ty_arr[i][j] if j < ny - 1 else 0.0

                    diag = B[i][j] + TxL + TxR + TyB + TyT
                    rhs = B[i][j] * P_old[i][j]

                    if i > 0:
                        rhs += TxL * P_new[i - 1][j]
                    if i < nx - 1:
                        rhs += TxR * P_new[i + 1][j]
                    if j > 0:
                        rhs += TyB * P_new[i][j - 1]
                    if j < ny - 1:
                        rhs += TyT * P_new[i][j + 1]

                    w_obj = next((w for w in self.wells if w.i == i and w.j == j), None)
                    if w_obj is not None and ctrl.mode != "shut_in":
                        wi = wi_dict.get(w_obj.name, 0.0)
                        diag += wi
                        rhs += wi * ctrl.bhp_target

                    if diag < 1e-30:
                        P_gs = P_old[i][j]
                    else:
                        P_gs = rhs / diag

                    P_sor = (1.0 - self.omega) * P_new[i][j] + self.omega * P_gs
                    P_sor = max(P_sor, 1.0e3)

                    if w_obj is not None:
                        if ctrl.mode == "inject":
                            P_sor = min(P_sor, ctrl.bhp_target)
                        elif ctrl.mode == "withdraw":
                            P_sor = max(P_sor, ctrl.bhp_target)

                    change = abs(P_sor - P_new[i][j])
                    max_change = max(max_change, change)
                    P_new[i][j] = P_sor

            if max_change < self.tol_Pa:
                break

        for w in self.wells:
            wi = wi_dict[w.name]
            P_w = P_new[w.i][w.j]
            if ctrl.mode == "inject":
                dp = max(ctrl.bhp_target - P_w, 0.0)
                q = min(wi * dp, ctrl.q_max)
                q = max(q, 0.0)
            elif ctrl.mode == "withdraw":
                dp = max(P_w - ctrl.bhp_target, 0.0)
                q = -min(wi * dp, ctrl.q_max)
                q = min(q, 0.0)
            else:
                q = 0.0
            per_well[w.name] = (q, P_w)

        return P_new, per_well

    def _update_temperature(
        self,
        P_old: list[list[float]],
        P_new: list[list[float]],
        T_old: list[list[float]],
        dt: float,
    ) -> list[list[float]]:
        T_geo = self.props.T_init
        T_new = [[0.0] * self.ny for _ in range(self.nx)]
        for i in range(self.nx):
            for j in range(self.ny):
                T_new[i][j] = coupled_temperature_update(
                    T_old=T_old[i][j],
                    P_old=max(P_old[i][j], 1.0),
                    P_new=max(P_new[i][j], 1.0),
                    T_geothermal=T_geo,
                    dt_s=dt,
                    tau_s=self.tau_s,
                    eos=self.eos,
                )
        return T_new

    def run_schedule(
        self,
        schedule: list[Schedule2D],
        dt_days: float = 0.2,
        snap_days: list[float] | None = None,
    ) -> Result2D:
        dt = dt_days * 86400.0
        nx, ny = self.nx, self.ny

        P = [[self.props.P_init] * ny for _ in range(nx)]
        T = [[self.props.T_init] * ny for _ in range(nx)]

        result = Result2D()
        for w in self.wells:
            result.per_well_q[w.name] = []
            result.per_well_BHP[w.name] = []

        snap_targets = sorted(snap_days or [])
        snap_idx = 0
        t_d = 0.0
        cum_inj = 0.0
        cum_wdr = 0.0

        P_max = self.props.P_init
        P_min = self.props.P_init
        T_max = self.props.T_init
        T_min = self.props.T_init

        for ctrl in schedule:
            n_steps = max(int(ctrl.duration_d / dt_days), 1)
            for _ in range(n_steps):
                P_new, per_well = self._solve_pressure(P, T, ctrl, dt)
                T_new = self._update_temperature(P, P_new, T, dt)
                t_d += dt_days

                total_q = sum(q for q, _ in per_well.values())
                if total_q >= 0.0:
                    cum_inj += total_q * dt
                else:
                    cum_wdr += (-total_q) * dt

                for w in self.wells:
                    q_w, bhp_w = per_well[w.name]
                    result.per_well_q[w.name].append(q_w)
                    result.per_well_BHP[w.name].append(pa_to_bar(bhp_w))

                P = P_new
                T = T_new

                P_flat = [P[i][j] for i in range(nx) for j in range(ny)]
                T_flat = [T[i][j] for i in range(nx) for j in range(ny)]

                P_max = max(P_max, max(P_flat))
                P_min = min(P_min, min(P_flat))
                T_max = max(T_max, max(T_flat))
                T_min = min(T_min, min(T_flat))

                m_tot = sum(
                    self.eos.density(P[i][j], T[i][j]) * self.Vp
                    for i in range(nx)
                    for j in range(ny)
                )

                result.times_d.append(t_d)
                result.P_avg_bar.append(pa_to_bar(sum(P_flat) / len(P_flat)))
                result.T_avg_C.append(K_to_C(sum(T_flat) / len(T_flat)))
                result.mass_t.append(kg_to_tonnes(m_tot))
                result.q_total_kg_s.append(total_q)

                while snap_idx < len(snap_targets) and t_d >= snap_targets[snap_idx]:
                    result.snapshots.append(
                        Snapshot2D(
                            time_d=t_d,
                            xs_m=[self.xs[i] for i in range(nx) for _ in range(ny)],
                            ys_m=[self.ys[j] for _ in range(nx) for j in range(ny)],
                            P_Pa=P_flat,
                            T_K=T_flat,
                        )
                    )
                    snap_idx += 1

        m_max = max(result.mass_t) if result.mass_t else 0.0
        m_min = min(result.mass_t) if result.mass_t else 0.0
        result.metrics = {
            "max_P_bar": pa_to_bar(P_max),
            "min_P_bar": pa_to_bar(P_min),
            "max_T_C": K_to_C(T_max),
            "min_T_C": K_to_C(T_min),
            "max_mass_t": m_max,
            "min_mass_t": m_min,
            "working_gas_t": m_max - m_min,
            "cushion_gas_t": m_min,
            "storage_eff_pct": (m_max - m_min) / max(m_max, 1e-9) * 100.0,
            "total_injected_t": kg_to_tonnes(cum_inj),
            "total_withdrawn_t": kg_to_tonnes(cum_wdr),
            "cycled_gas_t": kg_to_tonnes(cum_inj),
        }
        return result


H2Reservoir2DModel = UHSAreal2D


def uniform_k_field(nx: int, ny: int, k: float) -> list[list[float]]:
    return [[k] * ny for _ in range(nx)]


def thief_zone_k_field(
    nx: int,
    ny: int,
    k_bg: float,
    k_thief: float,
    width: int = 2,
) -> list[list[float]]:
    k_f = [[k_bg] * ny for _ in range(nx)]
    for i in range(nx):
        j_centre = int(round(i * (ny - 1) / max(nx - 1, 1)))
        for dj in range(-width, width + 1):
            j = j_centre + dj
            if 0 <= j < ny:
                k_f[i][j] = k_thief
    return k_f


def _base_2d_props(nx: int, ny: int, k_field: list[list[float]] | None = None) -> Reservoir2DProps:
    return Reservoir2DProps(
        Lx=1500.0,
        Ly=1500.0,
        h=100.0,
        phi=0.20,
        k=mD_to_m2(10.0),
        Ct=5.0e-10,
        P_init=bar_to_pa(80.0),
        T_init=316.0,
        nx=nx,
        ny=ny,
        tau_days=30.0,
        k_field=k_field,
    )


def _schedule_porous_2d(P_ref: float) -> list[Schedule2D]:
    bar = 1.0e5
    return [
        Schedule2D("inject_30d", "inject", 30.0, min(1.35 * P_ref, P_ref + 50.0 * bar), 35.0),
        Schedule2D("withdraw_30d", "withdraw", 30.0, max(0.60 * P_ref, P_ref - 40.0 * bar), 35.0),
        Schedule2D("shut_in_10d", "shut_in", 10.0, 0.0, 0.0),
    ]


def scenario_single_well(
    nx: int = 15, ny: int = 15
) -> tuple[Reservoir2DProps, list[WellSpec2D], list[Schedule2D]]:
    props = _base_2d_props(nx, ny)
    wells = [WellSpec2D(i=nx // 2, j=ny // 2, name="W1")]
    sched = _schedule_porous_2d(props.P_init)
    return props, wells, sched


def scenario_five_spot(
    nx: int = 15, ny: int = 15
) -> tuple[Reservoir2DProps, list[WellSpec2D], list[Schedule2D]]:
    props = _base_2d_props(nx, ny)
    wells = [
        WellSpec2D(i=0, j=0, name="INJ_SW"),
        WellSpec2D(i=0, j=ny - 1, name="INJ_NW"),
        WellSpec2D(i=nx - 1, j=0, name="INJ_SE"),
        WellSpec2D(i=nx - 1, j=ny - 1, name="INJ_NE"),
        WellSpec2D(i=nx // 2, j=ny // 2, name="CTRE"),
    ]
    sched = _schedule_porous_2d(props.P_init)
    return props, wells, sched


def scenario_thief_zone(
    nx: int = 15, ny: int = 15
) -> tuple[Reservoir2DProps, list[WellSpec2D], list[Schedule2D]]:
    k_bg = mD_to_m2(10.0)
    k_thief = mD_to_m2(100.0)
    k_field = thief_zone_k_field(nx, ny, k_bg, k_thief, width=1)
    props = _base_2d_props(nx, ny, k_field=k_field)
    wells = [
        WellSpec2D(i=0, j=0, name="INJ_SW"),
        WellSpec2D(i=0, j=ny - 1, name="INJ_NW"),
        WellSpec2D(i=nx - 1, j=0, name="INJ_SE"),
        WellSpec2D(i=nx - 1, j=ny - 1, name="INJ_NE"),
        WellSpec2D(i=nx // 2, j=ny // 2, name="CTRE"),
    ]
    sched = _schedule_porous_2d(props.P_init)
    return props, wells, sched


def pressure_uniformity(P_field: list[list[float]]) -> float:
    nx = len(P_field)
    ny = len(P_field[0]) if nx > 0 else 0
    n = nx * ny
    if n == 0:
        return 0.0
    vals = [P_field[i][j] for i in range(nx) for j in range(ny)]
    mean = sum(vals) / n
    var = sum((v - mean) ** 2 for v in vals) / n
    return math.sqrt(var) / max(mean, 1.0)


def export_2d_timeseries_csv(result: Result2D, filepath: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    well_names = sorted(result.per_well_q.keys())
    with open(filepath, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        header = ["time_d", "P_avg_bar", "T_avg_C", "mass_t", "q_total_kg_s"]
        for wn in well_names:
            header.extend([f"{wn}_q_kg_s", f"{wn}_BHP_bar"])
        w.writerow(header)
        for idx, td in enumerate(result.times_d):
            row = [
                f"{td:.4f}",
                f"{result.P_avg_bar[idx]:.4f}",
                f"{result.T_avg_C[idx]:.3f}",
                f"{result.mass_t[idx]:.4e}",
                f"{result.q_total_kg_s[idx]:.4e}",
            ]
            for wn in well_names:
                row.extend(
                    [f"{result.per_well_q[wn][idx]:.4e}", f"{result.per_well_BHP[wn][idx]:.4f}"]
                )
            w.writerow(row)
