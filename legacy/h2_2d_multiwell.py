"""
2D hydrogen storage simulator - Salt cavern scenario (Teesside-like).

Features:
- 2D finite-volume grid (x-y Cartesian, thin "box" representing cavern area)
- Peng-Robinson EOS for pure H2
- Single equivalent well in corner cell
- Rate-controlled well with bottom-hole pressure limits
- Cushion fill + fast daily cycles
- Time-series export (avg P/T/mass + per-well rate/BHP)
- 2D snapshot export (P/T/mass vs x,y) at chosen times

Dependency-free core (no NumPy/matplotlib). Plotting can be done by a separate script.
"""

import math
import csv
import os
import argparse
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple


# ---------------------------------------------------------------------------
# Peng-Robinson EOS for hydrogen
# ---------------------------------------------------------------------------


class H2PengRobinsonEOS:
    R_u = 8.314462618  # J/mol/K
    M = 2.016e-3       # kg/mol
    R_spec = R_u / M

    T_c = 33.2         # K
    P_c = 1.293e6      # Pa
    omega = -0.216

    def __init__(self, max_iter: int = 50, tol: float = 1e-6):
        self.max_iter = max_iter
        self.tol = tol
        self.a_c = 0.45724 * self.R_u ** 2 * self.T_c ** 2 / self.P_c
        self.b = 0.07780 * self.R_u * self.T_c / self.P_c
        self.kappa = (
            0.37464
            + 1.54226 * self.omega
            - 0.26992 * self.omega ** 2
        )

    def _alpha(self, T: float) -> float:
        Tr = max(T / self.T_c, 1e-6)
        return (1.0 + self.kappa * (1.0 - math.sqrt(Tr))) ** 2

    def z_factor(self, P: float, T: float) -> float:
        if P <= 0.0:
            return 1.0

        alpha = self._alpha(T)
        aT = self.a_c * alpha
        A = aT * P / (self.R_u ** 2 * T ** 2)
        B = self.b * P / (self.R_u * T)

        def f(Z: float) -> float:
            return (
                Z ** 3
                - (1.0 - B) * Z ** 2
                + (A - 3.0 * B ** 2 - 2.0 * B) * Z
                - (A * B - B ** 2 - B ** 3)
            )

        def df(Z: float) -> float:
            return (
                3.0 * Z ** 2
                - 2.0 * (1.0 - B) * Z
                + (A - 3.0 * B ** 2 - 2.0 * B)
            )

        Z = 1.0
        for _ in range(self.max_iter):
            F = f(Z)
            dF = df(Z)
            if abs(dF) < 1e-12:
                break
            dZ = -F / dF
            Z += dZ
            if abs(dZ) < self.tol:
                break

        if Z <= 0.0 or Z > 10.0:
            Z = 1.0
        return Z

    def density(self, P: float, T: float) -> float:
        T_eff = max(T, 50.0)
        Z = self.z_factor(P, T_eff)
        return P / (Z * self.R_spec * T_eff)


# ---------------------------------------------------------------------------
# Data classes & unit helpers
# ---------------------------------------------------------------------------


@dataclass
class ReservoirProperties2D:
    Lx: float
    Ly: float
    thickness: float
    porosity: float
    total_compressibility: float
    permeability: float
    gas_viscosity: float
    reference_pressure: float
    reference_temperature: float
    nx: int
    ny: int


@dataclass
class WellSpec:
    i: int
    j: int
    n_wells: int = 1
    name: str = "W"


@dataclass
class SchedulePeriod:
    name: str
    mode: str                 # "inject", "withdraw", "shut_in"
    duration_days: float
    bottomhole_pressure: float
    max_mass_rate: float      # kg/s per well


@dataclass
class Snapshot2D:
    time_s: float
    xs_m: List[float]
    ys_m: List[float]
    pressures_Pa: List[float]
    temperatures_K: List[float]
    masses_kg: List[float]


@dataclass
class SimulationResult2D:
    times: List[float] = field(default_factory=list)
    pressures: List[float] = field(default_factory=list)
    temperatures: List[float] = field(default_factory=list)
    masses: List[float] = field(default_factory=list)
    mass_rates: List[float] = field(default_factory=list)
    working_masses: List[float] = field(default_factory=list)
    metrics: Dict[str, float] = field(default_factory=dict)
    snapshots: List[Snapshot2D] = field(default_factory=list)
    per_well_rates: Dict[str, List[float]] = field(default_factory=dict)
    per_well_pressures: Dict[str, List[float]] = field(default_factory=dict)


def pa_to_bar(P: float) -> float:
    return P / 1e5


def pa_to_psi(P: float) -> float:
    return P / 6894.757293168


def K_to_C(T: float) -> float:
    return T - 273.15


def K_to_F(T: float) -> float:
    return K_to_C(T) * 9.0 / 5.0 + 32.0


def kg_to_tonnes(m: float) -> float:
    return m / 1000.0


def kg_to_short_tons(m: float) -> float:
    return m / 907.18474


def kg_to_lb(m: float) -> float:
    return m * 2.20462262185


# ---------------------------------------------------------------------------
# 2D model core
# ---------------------------------------------------------------------------


class H2Reservoir2DModel:
    def __init__(
        self,
        props: ReservoirProperties2D,
        wells: List[WellSpec],
        eos: Optional[H2PengRobinsonEOS] = None,
        gamma_gas: float = 1.4,
        thermal_relaxation_time_days: float = 30.0,
    ):
        self.props = props
        self.wells = wells
        self.eos = eos or H2PengRobinsonEOS()
        self.gamma = gamma_gas
        self.tau_T = thermal_relaxation_time_days * 86400.0

        self.dx = props.Lx / props.nx
        self.dy = props.Ly / props.ny

        self.xs = [(i + 0.5) * self.dx for i in range(props.nx)]
        self.ys = [(j + 0.5) * self.dy for j in range(props.ny)]

        self.V_bulk = self.dx * self.dy * self.props.thickness
        self.Vp_ref = self.V_bulk * self.props.porosity

    def _mass_from_state_cell(self, P: float, T: float) -> float:
        dp = P - self.props.reference_pressure
        arg = self.props.total_compressibility * dp
        arg_clamped = max(min(arg, 0.3), -0.3)
        Vp = self.Vp_ref * math.exp(arg_clamped)
        rho = self.eos.density(P, T)
        return rho * Vp

    def _dm_dP_numeric(self, P: float, T: float, dP: float = 1e4) -> float:
        P1 = max(P - dP, 1e3)
        P2 = P + dP
        m1 = self._mass_from_state_cell(P1, T)
        m2 = self._mass_from_state_cell(P2, T)
        return (m2 - m1) / (P2 - P1)

    def _build_initial_state(self) -> Tuple[List[List[float]], List[List[float]], List[List[float]]]:
        nx, ny = self.props.nx, self.props.ny
        P0 = self.props.reference_pressure
        T0 = self.props.reference_temperature
        P = [[P0 for _ in range(ny)] for _ in range(nx)]
        T = [[T0 for _ in range(ny)] for _ in range(nx)]
        m = [[self._mass_from_state_cell(P0, T0) for _ in range(ny)] for _ in range(nx)]
        return P, T, m

    def _compute_fluxes(
        self,
        P: List[List[float]],
        T: List[List[float]],
    ) -> Tuple[List[List[float]], List[List[float]]]:
        nx, ny = self.props.nx, self.props.ny
        k = self.props.permeability
        h = self.props.thickness
        mu = self.props.gas_viscosity

        rho = [[self.eos.density(P[i][j], T[i][j]) for j in range(ny)] for i in range(nx)]

        Tx = 2.0 * k * h * self.dy / (mu * self.dx)
        Ty = 2.0 * k * h * self.dx / (mu * self.dy)

        qx = [[0.0 for _ in range(ny)] for _ in range(nx - 1)]
        qy = [[0.0 for _ in range(ny - 1)] for _ in range(nx)]

        for i in range(nx - 1):
            for j in range(ny):
                dP = P[i][j] - P[i + 1][j]
                q_vol = Tx * dP
                rho_int = 0.5 * (rho[i][j] + rho[i + 1][j])
                qx[i][j] = q_vol * rho_int

        for i in range(nx):
            for j in range(ny - 1):
                dP = P[i][j] - P[i][j + 1]
                q_vol = Ty * dP
                rho_int = 0.5 * (rho[i][j] + rho[i][j + 1])
                qy[i][j] = q_vol * rho_int

        return qx, qy

    def _apply_well_rates(
        self,
        P: List[List[float]],
        period: SchedulePeriod,
    ) -> Tuple[List[List[float]], float]:
        nx, ny = self.props.nx, self.props.ny
        q_well = [[0.0 for _ in range(ny)] for _ in range(nx)]
        total_q = 0.0

        if period.mode == "shut_in":
            return q_well, total_q

        for w in self.wells:
            i, j, n_w = w.i, w.j, w.n_wells
            if i < 0 or i >= nx or j < 0 or j >= ny:
                continue
            P_cell = P[i][j]
            q = 0.0
            if period.mode == "inject":
                if P_cell < period.bottomhole_pressure:
                    q = n_w * period.max_mass_rate
            elif period.mode == "withdraw":
                if P_cell > period.bottomhole_pressure:
                    q = -n_w * period.max_mass_rate

            q_well[i][j] += q
            total_q += q

        return q_well, total_q

    def step(
        self,
        P: List[List[float]],
        T: List[List[float]],
        m: List[List[float]],
        period: SchedulePeriod,
        dt: float,
    ) -> Tuple[List[List[float]], List[List[float]], List[List[float]], float, Dict[str, Tuple[float, float]]]:
        nx, ny = self.props.nx, self.props.ny

        qx, qy = self._compute_fluxes(P, T)
        q_well, total_q = self._apply_well_rates(P, period)

        P_new = [[0.0 for _ in range(ny)] for _ in range(nx)]
        T_new = [[0.0 for _ in range(ny)] for _ in range(nx)]
        m_new = [[0.0 for _ in range(ny)] for _ in range(nx)]

        for i in range(nx):
            for j in range(ny):
                q_net = q_well[i][j]

                if i > 0:
                    q_net += qx[i - 1][j]
                if i < nx - 1:
                    q_net -= qx[i][j]
                if j > 0:
                    q_net += qy[i][j - 1]
                if j < ny - 1:
                    q_net -= qy[i][j]

                dm = q_net * dt
                dm_dP = self._dm_dP_numeric(P[i][j], T[i][j])
                if abs(dm_dP) < 1e-12:
                    dm_dP = m[i][j] / max(P[i][j], 1e3)

                dP = dm / dm_dP
                P_i_new = P[i][j] + dP
                P_i_new = max(min(P_i_new, 5e8), 1e3)

                P_old = max(P[i][j], 1e3)
                exponent = (self.gamma - 1.0) / self.gamma
                T_adiabatic = T[i][j] * (P_i_new / P_old) ** exponent

                T_geo = self.props.reference_temperature
                if self.tau_T > 0.0:
                    relax = math.exp(-dt / self.tau_T)
                    T_i_new = T_geo + (T_adiabatic - T_geo) * relax
                else:
                    T_i_new = T_adiabatic

                T_i_new = max(min(T_i_new, 600.0), 200.0)

                m_i_new = self._mass_from_state_cell(P_i_new, T_i_new)

                P_new[i][j] = P_i_new
                T_new[i][j] = T_i_new
                m_new[i][j] = m_i_new

        per_well: Dict[str, Tuple[float, float]] = {}
        for w in self.wells:
            i, j = w.i, w.j
            if 0 <= i < nx and 0 <= j < ny:
                name = w.name
                q_cell = q_well[i][j]
                P_bhp = P_new[i][j]
                per_well[name] = (q_cell, P_bhp)

        return P_new, T_new, m_new, total_q, per_well

    def run_schedule(
        self,
        schedule: List[SchedulePeriod],
        dt_seconds: float,
        snapshot_days: Optional[List[float]] = None,
        collect_time_series: bool = True,
    ) -> SimulationResult2D:
        nx, ny = self.props.nx, self.props.ny
        P, T, m = self._build_initial_state()

        result = SimulationResult2D()
        t = 0.0

        total_mass = sum(m[i][j] for i in range(nx) for j in range(ny))
        max_mass = total_mass
        min_mass = total_mass
        max_pressure = max(P[i][j] for i in range(nx) for j in range(ny))
        min_pressure = min(P[i][j] for i in range(nx) for j in range(ny))
        total_injected = 0.0
        total_withdrawn = 0.0

        snapshots: List[Snapshot2D] = []
        if snapshot_days:
            snapshot_times = sorted(d * 86400.0 for d in snapshot_days)
        else:
            snapshot_times = []
        snap_index = 0

        def domain_avg(field: List[List[float]]) -> float:
            s = 0.0
            for ii in range(nx):
                for jj in range(ny):
                    s += field[ii][jj]
            return s / (nx * ny)

        if collect_time_series:
            result.times.append(t)
            result.pressures.append(domain_avg(P))
            result.temperatures.append(domain_avg(T))
            result.masses.append(total_mass)
            result.mass_rates.append(0.0)
            for w in self.wells:
                name = w.name
                if name not in result.per_well_rates:
                    result.per_well_rates[name] = []
                    result.per_well_pressures[name] = []
                Pi = P[w.i][w.j]
                result.per_well_rates[name].append(0.0)
                result.per_well_pressures[name].append(Pi)

        for period in schedule:
            n_steps = int(period.duration_days * 86400.0 / dt_seconds)
            n_steps = max(n_steps, 1)

            for _ in range(n_steps):
                P, T, m, q_wells, per_well = self.step(P, T, m, period, dt_seconds)
                t += dt_seconds

                total_mass = sum(m[i][j] for i in range(nx) for j in range(ny))
                if q_wells >= 0.0:
                    total_injected += q_wells * dt_seconds
                else:
                    total_withdrawn += -q_wells * dt_seconds

                max_pressure = max(max_pressure, max(P[i][j] for i in range(nx) for j in range(ny)))
                min_pressure = min(min_pressure, min(P[i][j] for i in range(nx) for j in range(ny)))
                max_mass = max(max_mass, total_mass)
                min_mass = min(min_mass, total_mass)

                while snap_index < len(snapshot_times) and t >= snapshot_times[snap_index]:
                    flat_P = []
                    flat_T = []
                    flat_m = []
                    xs_flat = []
                    ys_flat = []
                    for i in range(nx):
                        for j in range(ny):
                            xs_flat.append(self.xs[i])
                            ys_flat.append(self.ys[j])
                            flat_P.append(P[i][j])
                            flat_T.append(T[i][j])
                            flat_m.append(m[i][j])
                    snapshots.append(
                        Snapshot2D(
                            time_s=snapshot_times[snap_index],
                            xs_m=xs_flat,
                            ys_m=ys_flat,
                            pressures_Pa=flat_P,
                            temperatures_K=flat_T,
                            masses_kg=flat_m,
                        )
                    )
                    snap_index += 1

                if collect_time_series:
                    result.times.append(t)
                    result.pressures.append(domain_avg(P))
                    result.temperatures.append(domain_avg(T))
                    result.masses.append(total_mass)
                    result.mass_rates.append(q_wells)

                    for w in self.wells:
                        name = w.name
                        if name not in result.per_well_rates:
                            result.per_well_rates[name] = []
                            result.per_well_pressures[name] = []
                        q, Pw = per_well.get(name, (0.0, P[w.i][w.j]))
                        result.per_well_rates[name].append(q)
                        result.per_well_pressures[name].append(Pw)

        working_gas = max_mass - min_mass
        cushion_mass = min_mass
        cushion_fraction = cushion_mass / max(max_mass, 1e-12)

        if collect_time_series:
            result.working_masses = [
                max(0.0, m_i - cushion_mass) for m_i in result.masses
            ]

        result.metrics = {
            "max_pressure_Pa": max_pressure,
            "min_pressure_Pa": min_pressure,
            "max_mass_kg": max_mass,
            "min_mass_kg": min_mass,
            "working_gas_kg": working_gas,
            "cushion_mass_kg": cushion_mass,
            "cushion_fraction": cushion_fraction,
            "total_injected_kg": total_injected,
            "total_withdrawn_kg": total_withdrawn,
            "cycle_gas_balance_kg": total_injected - total_withdrawn,
        }
        result.snapshots = snapshots
        return result


# ---------------------------------------------------------------------------
# Cavern scenario builder
# ---------------------------------------------------------------------------


def build_cavern_props(nx: int, ny: int) -> Tuple[ReservoirProperties2D, List[WellSpec], List[SchedulePeriod]]:
    cavern_volume = 7.0e4   # m3
    thickness = 100.0
    area = cavern_volume / thickness
    Lx = math.sqrt(area)
    Ly = Lx

    porosity = 0.99
    k = 1e-12
    mu = 1e-5
    Pref = 4.0e6          # 40 bar
    Tref = 293.0          # 20 °C

    props = ReservoirProperties2D(
        Lx=Lx,
        Ly=Ly,
        thickness=thickness,
        porosity=porosity,
        total_compressibility=1e-10,
        permeability=k,
        gas_viscosity=mu,
        reference_pressure=Pref,
        reference_temperature=Tref,
        nx=nx,
        ny=ny,
    )

    wells = [WellSpec(i=0, j=0, n_wells=1, name="C1")]

    bar_to_Pa = 1e5
    per_well_qmax = 200.0

    P_max = min(1.8 * Pref, Pref + 60.0 * bar_to_Pa)
    P_min = max(0.3 * Pref, Pref - 30.0 * bar_to_Pa)

    inject_bhp = P_max
    withdraw_bhp = P_min

    schedule: List[SchedulePeriod] = []
    schedule.append(
        SchedulePeriod(
            name="cushion_fill",
            mode="inject",
            duration_days=90.0,
            bottomhole_pressure=inject_bhp,
            max_mass_rate=per_well_qmax,
        )
    )

    n_cycles = 180
    for i in range(n_cycles):
        schedule.append(
            SchedulePeriod(
                name=f"cycle{i+1}_inject",
                mode="inject",
                duration_days=0.5,
                bottomhole_pressure=inject_bhp,
                max_mass_rate=per_well_qmax,
            )
        )
        schedule.append(
            SchedulePeriod(
                name=f"cycle{i+1}_withdraw",
                mode="withdraw",
                duration_days=0.4,
                bottomhole_pressure=withdraw_bhp,
                max_mass_rate=per_well_qmax,
            )
        )
        schedule.append(
            SchedulePeriod(
                name=f"cycle{i+1}_shut_in",
                mode="shut_in",
                duration_days=0.1,
                bottomhole_pressure=0.0,
                max_mass_rate=0.0,
            )
        )

    return props, wells, schedule


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------


def export_time_series_to_csv(result: SimulationResult2D, filepath: str) -> None:
    if not result.times:
        return

    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    well_names = sorted(result.per_well_rates.keys())

    with open(filepath, mode="w", newline="") as f:
        writer = csv.writer(f)
        header = [
            "time_s",
            "time_days",
            "pressure_Pa",
            "pressure_bar",
            "pressure_psi",
            "T_K",
            "T_C",
            "T_F",
            "mass_kg",
            "mass_tonnes",
            "mass_short_tons",
            "mass_lb",
            "working_mass_kg",
            "mass_rate_kg_s",
            "mass_rate_kg_day",
            "mass_rate_tonnes_day",
            "mass_rate_lb_hr",
        ]
        for name in well_names:
            header.extend([
                f"well_{name}_rate_kg_s",
                f"well_{name}_P_Pa",
                f"well_{name}_P_bar",
                f"well_{name}_P_psi",
            ])
        writer.writerow(header)

        if result.working_masses and len(result.working_masses) == len(result.times):
            working_series = result.working_masses
        else:
            working_series = [0.0] * len(result.times)

        n = len(result.times)
        for idx in range(n):
            t = result.times[idx]
            P = result.pressures[idx]
            T = result.temperatures[idx]
            m = result.masses[idx]
            q = result.mass_rates[idx]
            mw = working_series[idx]

            time_days = t / 86400.0
            P_bar = pa_to_bar(P)
            P_psi = pa_to_psi(P)
            T_C_v = K_to_C(T)
            T_F_v = K_to_F(T)
            m_tonnes = kg_to_tonnes(m)
            m_short_tons = kg_to_short_tons(m)
            m_lb = kg_to_lb(m)
            q_kg_day = q * 86400.0
            q_tonnes_day = kg_to_tonnes(q_kg_day)
            q_lb_hr = kg_to_lb(q) * 3600.0

            row = [
                f"{t:.1f}",
                f"{time_days:.5f}",
                f"{P:.3e}",
                f"{P_bar:.5f}",
                f"{P_psi:.5f}",
                f"{T:.3f}",
                f"{T_C_v:.3f}",
                f"{T_F_v:.3f}",
                f"{m:.6e}",
                f"{m_tonnes:.6e}",
                f"{m_short_tons:.6e}",
                f"{m_lb:.6e}",
                f"{mw:.6e}",
                f"{q:.6e}",
                f"{q_kg_day:.6e}",
                f"{q_tonnes_day:.6e}",
                f"{q_lb_hr:.6e}",
            ]

            for name in well_names:
                q_w_list = result.per_well_rates.get(name, [])
                P_w_list = result.per_well_pressures.get(name, [])
                q_w = q_w_list[idx] if idx < len(q_w_list) else 0.0
                P_w = P_w_list[idx] if idx < len(P_w_list) else 0.0
                row.extend([
                    f"{q_w:.6e}",
                    f"{P_w:.3e}",
                    f"{pa_to_bar(P_w):.5f}",
                    f"{pa_to_psi(P_w):.5f}",
                ])

            writer.writerow(row)


def export_snapshots_to_csv(result: SimulationResult2D, filepath_base: str) -> None:
    if not result.snapshots:
        return

    os.makedirs(os.path.dirname(os.path.abspath(filepath_base)), exist_ok=True)

    for idx, snap in enumerate(result.snapshots, start=1):
        fname = f"{filepath_base}_snap{idx}.csv"
        with open(fname, mode="w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "x_m",
                "y_m",
                "pressure_Pa",
                "pressure_bar",
                "pressure_psi",
                "T_K",
                "T_C",
                "mass_kg",
                "mass_tonnes",
                "snapshot_time_s",
                "snapshot_time_days",
            ])
            t_days = snap.time_s / 86400.0
            for x, y, P, T, m in zip(
                snap.xs_m,
                snap.ys_m,
                snap.pressures_Pa,
                snap.temperatures_K,
                snap.masses_kg,
            ):
                writer.writerow([
                    f"{x:.2f}",
                    f"{y:.2f}",
                    f"{P:.3e}",
                    f"{pa_to_bar(P):.5f}",
                    f"{pa_to_psi(P):.5f}",
                    f"{T:.3f}",
                    f"{K_to_C(T):.3f}",
                    f"{m:.6e}",
                    f"{kg_to_tonnes(m):.6e}",
                    f"{snap.time_s:.1f}",
                    f"{t_days:.5f}",
                ])
        print(f"  Snapshot {idx} exported to: {os.path.abspath(fname)}")


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------


def parse_snap_days(arg: str) -> Optional[List[float]]:
    if not arg:
        return None
    vals: List[float] = []
    for s in arg.split(","):
        s = s.strip()
        if not s:
            continue
        try:
            vals.append(float(s))
        except ValueError:
            pass
    return vals or None


def run_cavern(nx: int, ny: int, dt_days: float, snap_days: Optional[List[float]]) -> SimulationResult2D:
    props, wells, schedule = build_cavern_props(nx, ny)
    label = "Salt cavern (2D box)"
    ts_name = "h2_2d_cavern_timeseries.csv"
    snap_base = "h2_2d_cavern_field"

    eos = H2PengRobinsonEOS()
    model = H2Reservoir2DModel(props, wells, eos=eos)
    dt = dt_days * 86400.0

    result = model.run_schedule(
        schedule=schedule,
        dt_seconds=dt,
        snapshot_days=snap_days,
        collect_time_series=True,
    )

    print(f"\n=== Scenario: {label} ===")
    print(f"  grid      : {props.nx} x {props.ny}")
    print(f"  domain    : {props.Lx:.1f} x {props.Ly:.1f} m")
    print(f"  thickness : {props.thickness:.1f} m")
    print(f"  porosity  : {props.porosity:.3f}")
    print(f"  perm      : {props.permeability:.2e} m^2")
    print(f"  P_ref     : {props.reference_pressure/1e5:.2f} bar")
    print(f"  T_ref     : {K_to_C(props.reference_temperature):.1f} °C")

    for k, v in result.metrics.items():
        if "pressure" in k:
            print(f"  {k}: {v:.3e} Pa ({pa_to_bar(v):.2f} bar, {pa_to_psi(v):.2f} psi)")
        elif "mass" in k or "gas" in k:
            print(f"  {k}: {v:.3e} kg ({kg_to_tonnes(v):.3f} t, {kg_to_short_tons(v):.3f} short ton)")
        else:
            print(f"  {k}: {v:.3e}")

    print(f"  Final time      : {result.times[-1] / 86400.0:.1f} days")
    print(f"  Final avg P     : {result.pressures[-1] / 1e5:.2f} bar "
          f"({pa_to_psi(result.pressures[-1]):.2f} psi)")
    print(f"  Final avg T     : {K_to_C(result.temperatures[-1]):.2f} °C "
          f"({K_to_F(result.temperatures[-1]):.2f} °F)")
    print(f"  Final total mass: {result.masses[-1] / 1e6:.3f} kt "
          f"({kg_to_short_tons(result.masses[-1]) / 1000.0:.3f} kshort ton)")
    if result.working_masses:
        print(f"  Final working gas: {result.working_masses[-1] / 1e6:.3f} kt")

    export_time_series_to_csv(result, ts_name)
    print(f"  Time series exported to: {os.path.abspath(ts_name)}")
    export_snapshots_to_csv(result, snap_base)

    return result


def main():
    parser = argparse.ArgumentParser(description="2D H2 storage simulator - Salt cavern scenario")
    parser.add_argument("--nx", type=int, default=20, help="Number of grid cells in x-direction")
    parser.add_argument("--ny", type=int, default=20, help="Number of grid cells in y-direction")
    parser.add_argument(
        "--dt-days",
        type=float,
        default=0.05,
        help="Time step size in days (default 0.05 ~ 1.2 h, good for fast cycling)",
    )
    parser.add_argument(
        "--snap-days",
        type=str,
        default="",
        help="Comma-separated snapshot times in days (e.g. '90,90.5,90.9')",
    )
    args = parser.parse_args()
    snap_days = parse_snap_days(args.snap_days)
    run_cavern(nx=args.nx, ny=args.ny, dt_days=args.dt_days, snap_days=snap_days)


if __name__ == "__main__":
    main()

