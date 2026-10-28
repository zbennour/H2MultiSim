"""
h2multisim.models.three_d
=========================
Three-dimensional (3D) Cartesian finite-volume simulator for underground hydrogen
storage (UHS) in porous reservoirs with gravity segregation and vertical heterogeneity.
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

G_GRAV: float = 9.81  # m/s² gravitational acceleration


@dataclass
class Reservoir3DProps:
    """
    Properties for the 3D field-scale model.
    k_layers: Optional list of nz permeability values [m²] (layer 0 = bottom).
    """

    Lx: float  # m length in x
    Ly: float  # m length in y
    Lz: float  # m thickness in z
    phi: float  # porosity
    k: float  # background permeability [m²]
    Ct: float  # Pa⁻¹ total compressibility
    P_init: float  # Pa initial pressure
    T_init: float  # K initial temperature
    nx: int
    ny: int
    nz: int
    tau_days: float = 30.0
    k_layers: list[float] | None = None
    volume_translation: float = 0.0  # m³/mol Peneloux volume translation (e.g. -4.0e-6)
    peneloux: bool = False  # If True and volume_translation == 0, sets volume_translation = -4.0e-6

    def __post_init__(self) -> None:
        if self.peneloux and self.volume_translation == 0.0:
            self.volume_translation = -4.0e-6


ReservoirProperties3D = Reservoir3DProps


@dataclass
class WellSpec3D:
    """
    3D well specification with perforation interval [k_min, k_max].
    """

    i: int
    j: int
    k_min: int  # bottom perforated layer (0 = deepest)
    k_max: int  # top perforated layer
    r_w: float = 0.1  # m wellbore radius
    name: str = "W"


Well3D = WellSpec3D


@dataclass
class Schedule3D:
    """Operational period for 3D field-scale simulation."""

    name: str
    mode: str  # 'inject' | 'withdraw' | 'shut_in'
    duration_d: float  # days
    bhp_target: float  # Pa
    q_max: float  # kg/s total rate limit


SchedulePeriod3D = Schedule3D


@dataclass
class Snapshot3D:
    """3D spatial snapshot of state."""

    time_d: float
    xs_m: list[float]
    ys_m: list[float]
    zs_m: list[float]
    P_Pa: list[float]
    T_K: list[float]


@dataclass
class Result3D:
    """Simulation results from 3D model."""

    times_d: list[float] = field(default_factory=list)
    P_avg_bar: list[float] = field(default_factory=list)
    T_avg_C: list[float] = field(default_factory=list)
    mass_t: list[float] = field(default_factory=list)
    q_total_kg_s: list[float] = field(default_factory=list)
    per_well_q: dict[str, list[float]] = field(default_factory=dict)
    per_well_BHP: dict[str, list[float]] = field(default_factory=dict)
    snapshots: list[Snapshot3D] = field(default_factory=list)
    metrics: dict[str, float] = field(default_factory=dict)


SimulationResult3D = Result3D


class UHSFieldScale3D:
    """
    Implicit 3D Cartesian finite-volume UHS simulator with gravity.
    Solves pressure implicitly using red-black Gauss-Seidel SOR with Peaceman 3D well indices.
    """

    def __init__(
        self,
        props: Reservoir3DProps,
        wells: list[WellSpec3D],
        eos: PengRobinsonH2 | None = None,
        sor_omega: float = 1.5,
        max_iter: int = 2000,
        tol_Pa: float = 1.0e-4,
        max_picard: int = 20,
        picard_tol_Pa: float = 1.0e-4,
        hydrostatic_init: bool = True,
    ) -> None:
        self.props = props
        self.wells = wells
        if eos is None:
            vt = getattr(props, "volume_translation", 0.0)
            self.eos = PengRobinsonH2(volume_translation=vt)
        else:
            self.eos = eos
        self.omega = sor_omega
        self.max_iter = max_iter
        self.tol_Pa = tol_Pa
        self.max_picard = max_picard
        self.picard_tol_Pa = picard_tol_Pa
        self.hydrostatic_init = hydrostatic_init
        self.tau_s = props.tau_days * 86400.0
        self.g = G_GRAV

        self.nx = props.nx
        self.ny = props.ny
        self.nz = props.nz
        self.dx = props.Lx / props.nx
        self.dy = props.Ly / props.ny
        self.dz = props.Lz / props.nz

        self.xs = [(i + 0.5) * self.dx for i in range(self.nx)]
        self.ys = [(j + 0.5) * self.dy for j in range(self.ny)]
        self.zs = [(k + 0.5) * self.dz for k in range(self.nz)]

        self.Vp = props.phi * self.dx * self.dy * self.dz

        if props.k_layers is not None:
            if len(props.k_layers) != props.nz:
                raise ValueError(f"k_layers length {len(props.k_layers)} != nz={props.nz}")
            self._kz = props.k_layers
        else:
            self._kz = [props.k] * props.nz

    def _k_cell(self, i: int, j: int, k: int) -> float:
        return self._kz[k]

    def _k_harm_x(self, i: int, j: int, k: int) -> float:
        kL = self._k_cell(i, j, k)
        kR = self._k_cell(i + 1, j, k)
        return 2.0 * kL * kR / (kL + kR + 1e-40)

    def _k_harm_y(self, i: int, j: int, k: int) -> float:
        kB = self._k_cell(i, j, k)
        kT = self._k_cell(i, j + 1, k)
        return 2.0 * kB * kT / (kB + kT + 1e-40)

    def _k_harm_z(self, i: int, j: int, k: int) -> float:
        kD = self._k_cell(i, j, k)
        kU = self._k_cell(i, j, k + 1)
        return 2.0 * kD * kU / (kD + kU + 1e-40)

    def _Tx(
        self, P: list[list[list[float]]], T: list[list[list[float]]], i: int, j: int, k: int
    ) -> float:
        kh = self._k_harm_x(i, j, k)
        rhoL = self.eos.density(P[i][j][k], T[i][j][k])
        rhoR = self.eos.density(P[i + 1][j][k], T[i + 1][j][k])
        muL = viscosity_h2(T[i][j][k], P[i][j][k])
        muR = viscosity_h2(T[i + 1][j][k], P[i + 1][j][k])
        rho_h = 2.0 * rhoL * rhoR / (rhoL + rhoR + 1e-30)
        mu_h = 2.0 * muL * muR / (muL + muR + 1e-30)
        return kh * self.dy * self.dz * rho_h / (mu_h * self.dx)

    def _Ty(
        self, P: list[list[list[float]]], T: list[list[list[float]]], i: int, j: int, k: int
    ) -> float:
        kh = self._k_harm_y(i, j, k)
        rhoB = self.eos.density(P[i][j][k], T[i][j][k])
        rhoT = self.eos.density(P[i][j + 1][k], T[i][j + 1][k])
        muB = viscosity_h2(T[i][j][k], P[i][j][k])
        muT = viscosity_h2(T[i][j + 1][k], P[i][j + 1][k])
        rho_h = 2.0 * rhoB * rhoT / (rhoB + rhoT + 1e-30)
        mu_h = 2.0 * muB * muT / (muB + muT + 1e-30)
        return kh * self.dx * self.dz * rho_h / (mu_h * self.dy)

    def _Tz(
        self, P: list[list[list[float]]], T: list[list[list[float]]], i: int, j: int, k: int
    ) -> float:
        kh = self._k_harm_z(i, j, k)
        rhoD = self.eos.density(P[i][j][k], T[i][j][k])
        rhoU = self.eos.density(P[i][j][k + 1], T[i][j][k + 1])
        muD = viscosity_h2(T[i][j][k], P[i][j][k])
        muU = viscosity_h2(T[i][j][k + 1], P[i][j][k + 1])
        rho_h = 2.0 * rhoD * rhoU / (rhoD + rhoU + 1e-30)
        mu_h = 2.0 * muD * muU / (muD + muU + 1e-30)
        return kh * self.dx * self.dy * rho_h / (mu_h * self.dz)

    def _well_index(
        self, w: WellSpec3D, P: list[list[list[float]]], T: list[list[list[float]]], k: int
    ) -> float:
        r_eq = 0.28 * math.sqrt(self.dx * self.dy)
        ln_r = math.log(max(r_eq, w.r_w * 1.001) / w.r_w)
        ln_r = max(ln_r, 0.1)
        k_w = self._k_cell(w.i, w.j, k)
        rho = self.eos.density(P[w.i][w.j][k], T[w.i][w.j][k])
        mu = viscosity_h2(T[w.i][w.j][k], P[w.i][w.j][k])
        return 2.0 * math.pi * k_w * self.dz * rho / (mu * ln_r)

    def _gravity_dp(
        self, P: list[list[list[float]]], T: list[list[list[float]]], i: int, j: int, k: int
    ) -> float:
        rhoD = self.eos.density(P[i][j][k], T[i][j][k])
        rhoU = self.eos.density(P[i][j][k + 1], T[i][j][k + 1])
        rho_h = 2.0 * rhoD * rhoU / (rhoD + rhoU + 1e-30)
        return rho_h * self.g * self.dz

    def _solve_picard(
        self,
        P_old: list[list[list[float]]],
        T_old: list[list[list[float]]],
        ctrl: Schedule3D,
        dt: float,
        tx: list[list[list[float]]],
        ty: list[list[list[float]]],
        tz: list[list[list[float]]],
        gz: list[list[list[float]]],
        rho_old: list[list[list[float]]],
        rate_q: dict[str, float],
        shut_wells: set[str],
    ) -> tuple[
        list[list[list[float]]],
        list[list[list[float]]],
        dict[str, tuple[float, float]],
    ]:
        """Picard iteration on pressure and temperature for one 3D time step."""
        nx, ny, nz = self.nx, self.ny, self.nz
        T_geo = self.props.T_init
        Pk = [[list(P_old[i][j]) for j in range(ny)] for i in range(nx)]
        Tk = [[list(T_old[i][j]) for j in range(ny)] for i in range(nx)]

        for _picard in range(self.max_picard):
            Tk = [
                [
                    [
                        coupled_temperature_update(
                            T_old=T_old[i][j][k],
                            P_old=max(P_old[i][j][k], 1.0),
                            P_new=max(Pk[i][j][k], 1.0),
                            T_geothermal=T_geo,
                            dt_s=dt,
                            tau_s=self.tau_s,
                            eos=self.eos,
                        )
                        for k in range(nz)
                    ]
                    for j in range(ny)
                ]
                for i in range(nx)
            ]
            rho_k = [
                [[self.eos.density(Pk[i][j][k], Tk[i][j][k]) for k in range(nz)] for j in range(ny)]
                for i in range(nx)
            ]
            drho = [
                [[self.eos.drho_dP(Pk[i][j][k], Tk[i][j][k]) for k in range(nz)] for j in range(ny)]
                for i in range(nx)
            ]

            B = [
                [[self.Vp * drho[i][j][k] / dt for k in range(nz)] for j in range(ny)]
                for i in range(nx)
            ]
            rhs0 = [
                [
                    [
                        B[i][j][k] * Pk[i][j][k]
                        - self.Vp * (rho_k[i][j][k] - rho_old[i][j][k]) / dt
                        for k in range(nz)
                    ]
                    for j in range(ny)
                ]
                for i in range(nx)
            ]

            well_diag = [[[0.0] * nz for _ in range(ny)] for _ in range(nx)]
            well_rhs = [[[0.0] * nz for _ in range(ny)] for _ in range(nx)]
            wi_now: dict[str, dict[int, float]] = {}
            for w in self.wells:
                wi_now[w.name] = {}
                for kl in range(w.k_min, w.k_max + 1):
                    if 0 <= kl < nz:
                        wi_now[w.name][kl] = self._well_index(w, Pk, Tk, kl)
                if ctrl.mode == "shut_in" or w.name in shut_wells:
                    continue
                if w.name in rate_q:
                    sum_wi = sum(wi_now[w.name].values())
                    if sum_wi <= 0.0:
                        continue
                    sum_wiP = sum(wi * Pk[w.i][w.j][kl] for kl, wi in wi_now[w.name].items())
                    bhp_rc = (rate_q[w.name] + sum_wiP) / sum_wi
                    for kl, wi in wi_now[w.name].items():
                        well_rhs[w.i][w.j][kl] += wi * (bhp_rc - Pk[w.i][w.j][kl])
                else:
                    for kl, wi in wi_now[w.name].items():
                        well_diag[w.i][w.j][kl] += wi
                        well_rhs[w.i][w.j][kl] += wi * ctrl.bhp_target

            P_new = [[list(Pk[i][j]) for j in range(ny)] for i in range(nx)]
            for _it in range(self.max_iter):
                max_change = 0.0
                for parity in (0, 1):
                    for i in range(nx):
                        for j in range(ny):
                            for k in range(nz):
                                if (i + j + k) % 2 != parity:
                                    continue

                                TxL = tx[i - 1][j][k] if i > 0 else 0.0
                                TxR = tx[i][j][k] if i < nx - 1 else 0.0
                                TyB = ty[i][j - 1][k] if j > 0 else 0.0
                                TyT = ty[i][j][k] if j < ny - 1 else 0.0
                                TzD = tz[i][j][k - 1] if k > 0 else 0.0
                                TzU = tz[i][j][k] if k < nz - 1 else 0.0
                                gD = gz[i][j][k - 1] if k > 0 else 0.0
                                gU = gz[i][j][k] if k < nz - 1 else 0.0

                                diag = (
                                    B[i][j][k]
                                    + TxL
                                    + TxR
                                    + TyB
                                    + TyT
                                    + TzD
                                    + TzU
                                    + well_diag[i][j][k]
                                )
                                rhs = rhs0[i][j][k] + well_rhs[i][j][k]

                                if i > 0:
                                    rhs += TxL * P_new[i - 1][j][k]
                                if i < nx - 1:
                                    rhs += TxR * P_new[i + 1][j][k]
                                if j > 0:
                                    rhs += TyB * P_new[i][j - 1][k]
                                if j < ny - 1:
                                    rhs += TyT * P_new[i][j + 1][k]
                                # Φ = P + ρgz; flux from k to k+1 is
                                # T·(P_k - P_{k+1} - ρgΔz).  Hence the lower
                                # neighbour contributes (P_{k-1} - gD) and the
                                # upper neighbour (P_{k+1} + gU).  The previous
                                # signs were inverted, producing an unstable
                                # (top-heavy) hydrostatic profile.
                                if k > 0:
                                    rhs += TzD * (P_new[i][j][k - 1] - gD)
                                if k < nz - 1:
                                    rhs += TzU * (P_new[i][j][k + 1] + gU)

                                P_gs = rhs / max(diag, 1e-30)
                                P_sor = (1.0 - self.omega) * P_new[i][j][k] + self.omega * P_gs
                                P_sor = max(P_sor, 1.0e3)

                                change = abs(P_sor - P_new[i][j][k])
                                max_change = max(max_change, change)
                                P_new[i][j][k] = P_sor

                if max_change < self.tol_Pa:
                    break

            change = max(
                abs(P_new[i][j][k] - Pk[i][j][k])
                for i in range(nx)
                for j in range(ny)
                for k in range(nz)
            )
            Pk = P_new
            if change < self.picard_tol_Pa:
                break

        # Consistent final temperature and well indices.
        Tk = [
            [
                [
                    coupled_temperature_update(
                        T_old=T_old[i][j][k],
                        P_old=max(P_old[i][j][k], 1.0),
                        P_new=max(Pk[i][j][k], 1.0),
                        T_geothermal=T_geo,
                        dt_s=dt,
                        tau_s=self.tau_s,
                        eos=self.eos,
                    )
                    for k in range(nz)
                ]
                for j in range(ny)
            ]
            for i in range(nx)
        ]

        per_well: dict[str, tuple[float, float]] = {}
        for w in self.wells:
            wi_layers = {
                kl: self._well_index(w, Pk, Tk, kl)
                for kl in range(w.k_min, w.k_max + 1)
                if 0 <= kl < nz
            }
            k_mid = (w.k_min + w.k_max) // 2
            if ctrl.mode == "shut_in" or w.name in shut_wells:
                q_total = 0.0
                bhp = Pk[w.i][w.j][k_mid]
            elif w.name in rate_q:
                q_total = rate_q[w.name]
                sum_wi = sum(wi_layers.values())
                if sum_wi > 0.0:
                    sum_wiP = sum(wi * Pk[w.i][w.j][kl] for kl, wi in wi_layers.items())
                    bhp = (q_total + sum_wiP) / sum_wi
                else:
                    bhp = Pk[w.i][w.j][k_mid]
            else:
                q_total = sum(
                    wi * (ctrl.bhp_target - Pk[w.i][w.j][kl]) for kl, wi in wi_layers.items()
                )
                bhp = ctrl.bhp_target
            per_well[w.name] = (q_total, bhp)

        return Pk, Tk, per_well

    def _solve_pressure(
        self,
        P_old: list[list[list[float]]],
        T_old: list[list[list[float]]],
        ctrl: Schedule3D,
        dt: float,
    ) -> tuple[
        list[list[list[float]]],
        list[list[list[float]]],
        dict[str, tuple[float, float]],
    ]:
        """
        Solve one pressure/temperature step with an active-set well model.

        The total well rate is capped (not each perforation separately), and a
        rate-controlled well is distributed over its perforations through a
        single implied bottom-hole pressure.  Cell pressures are never
        overwritten, so rates, BHPs and mass remain consistent.
        """
        nx, ny, nz = self.nx, self.ny, self.nz
        tx = [
            [
                [self._Tx(P_old, T_old, i, j, k) if i < nx - 1 else 0.0 for k in range(nz)]
                for j in range(ny)
            ]
            for i in range(nx)
        ]
        ty = [
            [
                [self._Ty(P_old, T_old, i, j, k) if j < ny - 1 else 0.0 for k in range(nz)]
                for j in range(ny)
            ]
            for i in range(nx)
        ]
        tz = [
            [
                [self._Tz(P_old, T_old, i, j, k) if k < nz - 1 else 0.0 for k in range(nz)]
                for j in range(ny)
            ]
            for i in range(nx)
        ]
        gz = [
            [
                [self._gravity_dp(P_old, T_old, i, j, k) if k < nz - 1 else 0.0 for k in range(nz)]
                for j in range(ny)
            ]
            for i in range(nx)
        ]
        rho_old = [
            [
                [self.eos.density(P_old[i][j][k], T_old[i][j][k]) for k in range(nz)]
                for j in range(ny)
            ]
            for i in range(nx)
        ]

        rate_q: dict[str, float] = {}
        shut_wells: set[str] = set()
        if ctrl.mode != "shut_in":
            for w in self.wells:
                q_try = 0.0
                for kl in range(w.k_min, w.k_max + 1):
                    if 0 <= kl < nz:
                        wi = self._well_index(w, P_old, T_old, kl)
                        q_try += wi * (ctrl.bhp_target - P_old[w.i][w.j][kl])
                if (ctrl.mode == "inject" and q_try < 0.0) or (
                    ctrl.mode == "withdraw" and q_try > 0.0
                ):
                    shut_wells.add(w.name)
                elif ctrl.q_max > 0.0 and abs(q_try) > ctrl.q_max:
                    rate_q[w.name] = math.copysign(ctrl.q_max, q_try)

        P_new = [[list(P_old[i][j]) for j in range(ny)] for i in range(nx)]
        T_new = [[list(T_old[i][j]) for j in range(ny)] for i in range(nx)]
        per_well: dict[str, tuple[float, float]] = {}
        for _attempt in range(5):
            P_new, T_new, per_well = self._solve_picard(
                P_old, T_old, ctrl, dt, tx, ty, tz, gz, rho_old, rate_q, shut_wells
            )
            changed = False

            for name in list(rate_q):
                w = next(ww for ww in self.wells if ww.name == name)
                bhp_rc = per_well[name][1]
                binding = (ctrl.mode == "inject" and bhp_rc <= ctrl.bhp_target) or (
                    ctrl.mode == "withdraw" and bhp_rc >= ctrl.bhp_target
                )
                if not binding:
                    del rate_q[name]
                    changed = True

            for name in list(shut_wells):
                w = next(ww for ww in self.wells if ww.name == name)
                q_try = sum(
                    self._well_index(w, P_new, T_new, kl) * (ctrl.bhp_target - P_new[w.i][w.j][kl])
                    for kl in range(w.k_min, w.k_max + 1)
                    if 0 <= kl < nz
                )
                if (ctrl.mode == "inject" and q_try > 0.0) or (
                    ctrl.mode == "withdraw" and q_try < 0.0
                ):
                    shut_wells.remove(name)
                    changed = True

            for w in self.wells:
                if ctrl.mode == "shut_in" or w.name in rate_q or w.name in shut_wells:
                    continue
                q_w = per_well[w.name][0]
                if (ctrl.mode == "inject" and q_w < 0.0) or (ctrl.mode == "withdraw" and q_w > 0.0):
                    shut_wells.add(w.name)
                    changed = True

            if not changed:
                break

        return P_new, T_new, per_well

    def _update_temperature(
        self,
        P_old: list[list[list[float]]],
        P_new: list[list[list[float]]],
        T_old: list[list[list[float]]],
        dt: float,
    ) -> list[list[list[float]]]:
        T_geo = self.props.T_init
        T_new = [[[0.0] * self.nz for _ in range(self.ny)] for _ in range(self.nx)]
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    T_new[i][j][k] = coupled_temperature_update(
                        T_old=T_old[i][j][k],
                        P_old=max(P_old[i][j][k], 1.0),
                        P_new=max(P_new[i][j][k], 1.0),
                        T_geothermal=T_geo,
                        dt_s=dt,
                        tau_s=self.tau_s,
                        eos=self.eos,
                    )
        return T_new

    def vertical_pressure_profile(
        self,
        P: list[list[list[float]]],
        i: int | None = None,
        j: int | None = None,
    ) -> list[tuple[float, float]]:
        """Return vertical column [(z_m, P_bar)] at (i, j)."""
        idx = i if i is not None else self.nx // 2
        jdx = j if j is not None else self.ny // 2
        return [(self.zs[k], pa_to_bar(P[idx][jdx][k])) for k in range(self.nz)]

    def run_schedule(
        self,
        schedule: list[Schedule3D],
        snap_times_d: list[float] | None = None,
        dt_days: float = 1.0,
        verbose: bool = False,
    ) -> Result3D:
        nx, ny, nz = self.nx, self.ny, self.nz
        P = [[[self.props.P_init] * nz for _ in range(ny)] for _ in range(nx)]
        if self.hydrostatic_init:
            rho0 = self.eos.density(self.props.P_init, self.props.T_init)
            z_ref = (nz - 1) * self.props.Lz / (2.0 * nz)
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        z = (k + 0.5) * self.props.Lz / nz
                        P[i][j][k] = self.props.P_init + rho0 * self.g * (z_ref - z)
        T = [[[self.props.T_init] * nz for _ in range(ny)] for _ in range(nx)]
        t = 0.0

        result = Result3D()
        for w in self.wells:
            result.per_well_q[w.name] = []
            result.per_well_BHP[w.name] = []

        snap_s = set((td * 86400.0) for td in (snap_times_d or []))
        cum_inj = 0.0
        cum_wdr = 0.0
        P_max = self.props.P_init
        P_min = self.props.P_init
        T_max = self.props.T_init
        T_min = self.props.T_init

        def _avg_P() -> float:
            return sum(P[i][j][k] for i in range(nx) for j in range(ny) for k in range(nz)) / (
                nx * ny * nz
            )

        def _avg_T() -> float:
            return sum(T[i][j][k] for i in range(nx) for j in range(ny) for k in range(nz)) / (
                nx * ny * nz
            )

        def _total_mass() -> float:
            return sum(
                self.eos.density(P[i][j][k], T[i][j][k]) * self.Vp
                for i in range(nx)
                for j in range(ny)
                for k in range(nz)
            )

        def _record(q_tot: float, per_well: dict[str, tuple[float, float]]) -> None:
            result.times_d.append(t / 86400.0)
            result.P_avg_bar.append(pa_to_bar(_avg_P()))
            result.T_avg_C.append(K_to_C(_avg_T()))
            result.mass_t.append(kg_to_tonnes(_total_mass()))
            result.q_total_kg_s.append(q_tot)
            for w in self.wells:
                q, BHP = per_well.get(w.name, (0.0, P[w.i][w.j][(w.k_min + w.k_max) // 2]))
                result.per_well_q[w.name].append(q)
                result.per_well_BHP[w.name].append(pa_to_bar(BHP))

        _record(0.0, {w.name: (0.0, P[w.i][w.j][w.k_min]) for w in self.wells})

        for ctrl in schedule:
            t_end = t + ctrl.duration_d * 86400.0
            n_steps = max(int(ctrl.duration_d / dt_days), 1)
            dt_act = (t_end - t) / n_steps

            if verbose:
                print(
                    f"  '{ctrl.name}': {ctrl.mode}, {ctrl.duration_d:.1f}d, BHP={pa_to_bar(ctrl.bhp_target):.1f}bar"
                )

            for _ in range(n_steps):
                P_old = [[list(P[i][j]) for j in range(ny)] for i in range(nx)]
                P, T, per_well = self._solve_pressure(P_old, T, ctrl, dt_act)
                t += dt_act

                q_inj = sum(max(q, 0.0) for q, _ in per_well.values())
                q_wdr = sum(max(-q, 0.0) for q, _ in per_well.values())
                cum_inj += q_inj * dt_act
                cum_wdr += q_wdr * dt_act
                q_tot = q_inj - q_wdr

                P_flat = [P[i][j][k] for i in range(nx) for j in range(ny) for k in range(nz)]
                T_flat = [T[i][j][k] for i in range(nx) for j in range(ny) for k in range(nz)]
                P_max = max(P_max, max(P_flat))
                P_min = min(P_min, min(P_flat))
                T_max = max(T_max, max(T_flat))
                T_min = min(T_min, min(T_flat))

                _record(q_tot, per_well)

                for ts in snap_s:
                    if abs(t - ts) < dt_act * 0.6:
                        result.snapshots.append(
                            Snapshot3D(
                                time_d=t / 86400.0,
                                xs_m=[
                                    self.xs[i]
                                    for i in range(nx)
                                    for _ in range(ny)
                                    for _ in range(nz)
                                ],
                                ys_m=[
                                    self.ys[j]
                                    for _ in range(nx)
                                    for j in range(ny)
                                    for _ in range(nz)
                                ],
                                zs_m=[
                                    self.zs[k]
                                    for _ in range(nx)
                                    for _ in range(ny)
                                    for k in range(nz)
                                ],
                                P_Pa=P_flat,
                                T_K=T_flat,
                            )
                        )

        m_max = max(result.mass_t) if result.mass_t else 0.0
        m_min = min(result.mass_t) if result.mass_t else 0.0
        m_initial_kg = result.mass_t[0] * 1000.0 if result.mass_t else 0.0
        m_final_kg = result.mass_t[-1] * 1000.0 if result.mass_t else 0.0
        net_kg = cum_inj - cum_wdr
        mb_err = abs((m_final_kg - m_initial_kg) - net_kg) / max(abs(net_kg), 1.0)
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
            "initial_mass_t": kg_to_tonnes(m_initial_kg),
            "final_mass_t": kg_to_tonnes(m_final_kg),
            "total_injected_t": kg_to_tonnes(cum_inj),
            "total_withdrawn_t": kg_to_tonnes(cum_wdr),
            "net_stored_t": kg_to_tonnes(net_kg),
            "cycled_gas_t": kg_to_tonnes(cum_inj + cum_wdr),
            "mass_balance_error": mb_err,
        }
        return result


H2Reservoir3DModel = UHSFieldScale3D


def _base_props_3d(
    nx: int, ny: int, nz: int, k_layers: list[float] | None = None
) -> Reservoir3DProps:
    return Reservoir3DProps(
        Lx=1500.0,
        Ly=1500.0,
        Lz=100.0,
        phi=0.20,
        k=mD_to_m2(10.0),
        Ct=5.0e-10,
        P_init=bar_to_pa(80.0),
        T_init=316.0,
        nx=nx,
        ny=ny,
        nz=nz,
        tau_days=30.0,
        k_layers=k_layers,
    )


def _schedule_3d(P_ref: float) -> list[Schedule3D]:
    bar = 1.0e5
    return [
        Schedule3D("inject_30d", "inject", 30.0, min(1.35 * P_ref, P_ref + 50.0 * bar), 35.0),
        Schedule3D("withdraw_30d", "withdraw", 30.0, max(0.60 * P_ref, P_ref - 40.0 * bar), 35.0),
        Schedule3D("shut_in_10d", "shut_in", 10.0, 0.0, 0.0),
    ]


def scenario_single_well_3d(
    nx: int = 8, ny: int = 8, nz: int = 5
) -> tuple[Reservoir3DProps, list[WellSpec3D], list[Schedule3D]]:
    props = _base_props_3d(nx, ny, nz)
    wells = [WellSpec3D(i=nx // 2, j=ny // 2, k_min=0, k_max=nz - 1, name="W1")]
    sched = _schedule_3d(props.P_init)
    return props, wells, sched


def scenario_gravity_override(
    nx: int = 8, ny: int = 8, nz: int = 5
) -> tuple[Reservoir3DProps, list[WellSpec3D], list[Schedule3D]]:
    props = _base_props_3d(nx, ny, nz)
    wells = [WellSpec3D(i=nx // 2, j=ny // 2, k_min=0, k_max=1, name="W_BTM")]
    sched = _schedule_3d(props.P_init)
    return props, wells, sched


def scenario_layered_reservoir(
    nx: int = 8, ny: int = 8, nz: int = 5
) -> tuple[Reservoir3DProps, list[WellSpec3D], list[Schedule3D]]:
    k_bg = mD_to_m2(5.0)
    k_fast = mD_to_m2(200.0)
    k_layers = [k_bg] * nz
    k_layers[nz // 2] = k_fast

    props = _base_props_3d(nx, ny, nz, k_layers=k_layers)
    wells = [WellSpec3D(i=nx // 2, j=ny // 2, k_min=0, k_max=nz - 1, name="W1")]
    sched = _schedule_3d(props.P_init)
    return props, wells, sched


def export_3d_timeseries_csv(result: Result3D, filepath: str) -> None:
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
