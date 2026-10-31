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

from h2multisim.models._common import AndersonAccelerator, check_inputs, pcg_solve, validate_mode
from h2multisim.models._composition import (
    BinaryCompositionMixin,
    add_transport_link,
    solve_binary_transport,
    validate_cushion_props,
)
from h2multisim.thermo import (
    PengRobinsonH2,
    RelativePermeabilityModel,
    coupled_temperature_update,
    heating_values_mass,
)
from h2multisim.units import (
    J_to_GWh,
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
    Ct: float  # Pa⁻¹ pore (rock) compressibility; gas compressibility comes from the EOS
    P_init: float  # Pa initial pressure
    T_init: float  # K initial temperature
    nx: int
    ny: int
    nz: int
    tau_days: float = 30.0
    k_layers: list[float] | None = None
    volume_translation: float = 0.0  # m³/mol Peneloux volume translation (e.g. -4.0e-6)
    peneloux: bool = False  # If True and volume_translation == 0, sets volume_translation = -4.0e-6
    swr: float = 0.0  # [-] connate / residual immobile water saturation S_wr
    sgc: float = 0.0  # [-] critical / trapped gas saturation S_gc
    rel_perm: RelativePermeabilityModel | None = (
        None  # Optional two-phase relative permeability model
    )
    rock_heat_capacity: float = 2.2e6  # J/(m³·K) rock solids, local thermal equilibrium with gas
    cushion_gas: str = "H2"  # Native cushion gas species ('H2', 'CH4', 'CO2', 'N2')
    initial_h2_mole_fraction: float | None = None  # None -> 1.0 if cushion is H2 else 0.0
    dispersivity_m: float = 0.5  # m longitudinal mechanical dispersivity α_L
    molecular_diffusivity_m2_s: float = 1.0e-5  # m²/s effective pore molecular diffusivity

    def __post_init__(self) -> None:
        check_inputs(
            "Reservoir3DProps",
            positive={
                "Lx": self.Lx,
                "Ly": self.Ly,
                "Lz": self.Lz,
                "nx": float(self.nx),
                "ny": float(self.ny),
                "nz": float(self.nz),
                "k": self.k,
                "P_init": self.P_init,
                "T_init": self.T_init,
            },
            non_negative={
                "Ct": self.Ct,
                "tau_days": self.tau_days,
                "rock_heat_capacity": self.rock_heat_capacity,
                "dispersivity_m": self.dispersivity_m,
                "molecular_diffusivity_m2_s": self.molecular_diffusivity_m2_s,
            },
            fraction={"phi": self.phi},
        )
        validate_cushion_props(self, "Reservoir3DProps")
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
    well_modes: dict[str, str] | None = None  # optional per-well override

    def __post_init__(self) -> None:
        validate_mode(self.mode, f"schedule period {self.name!r}")
        for well_name, mode in (self.well_modes or {}).items():
            validate_mode(mode, f"schedule period {self.name!r}, well {well_name!r}")
        check_inputs(
            f"schedule period {self.name!r}",
            positive={"duration_d": self.duration_d},
            non_negative={"bhp_target": self.bhp_target, "q_max": self.q_max},
        )

    def mode_for(self, well_name: str) -> str:
        """Operating mode for a named well (falls back to the period mode)."""
        if self.well_modes and well_name in self.well_modes:
            return self.well_modes[well_name]
        return self.mode


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
    h2_mole_fraction: list[float] = field(default_factory=list)


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
    per_well_h2_mole_fraction: dict[str, list[float]] = field(default_factory=dict)
    snapshots: list[Snapshot3D] = field(default_factory=list)
    metrics: dict[str, float] = field(default_factory=dict)


SimulationResult3D = Result3D


class UHSFieldScale3D(BinaryCompositionMixin):
    """
    Implicit 3D Cartesian finite-volume UHS simulator with gravity: implicit pressure
    (PCG or red-black SOR), Peaceman 3D well indices, and optional binary H2 +
    cushion-gas composition transport (implicit upwind advection–dispersion).
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
        linear_solver: str = "pcg",
    ) -> None:
        if linear_solver not in ("pcg", "sor"):
            raise ValueError(f"linear_solver must be 'pcg' or 'sor', got {linear_solver!r}")
        self.linear_solver = linear_solver
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
        self._last_picard_residual_Pa: float = 0.0
        self._sor_failures: int = 0  # inner linear solves (PCG or SOR) that hit max_iter
        self._active_set_failures: int = 0  # steps whose well active set did not settle
        self._n_exp: list = []  # per-cell isentropic exponent for the current step

        self.nx = props.nx
        self.ny = props.ny
        self.nz = props.nz
        self.dx = props.Lx / props.nx
        self.dy = props.Ly / props.ny
        self.dz = props.Lz / props.nz

        self.xs = [(i + 0.5) * self.dx for i in range(self.nx)]
        self.ys = [(j + 0.5) * self.dy for j in range(self.ny)]
        self.zs = [(k + 0.5) * self.dz for k in range(self.nz)]

        self.rel_perm: RelativePermeabilityModel = getattr(
            props, "rel_perm", None
        ) or RelativePermeabilityModel(
            swr=float(getattr(props, "swr", 0.0)),
            sgc=float(getattr(props, "sgc", 0.0)),
        )
        self._phi_eff: float = props.phi * self.rel_perm.effective_porosity_factor()
        self._krg_eff: float = self.rel_perm.effective_permeability_factor()
        self.Vp = self._phi_eff * self.dx * self.dy * self.dz

        seen_names: set[str] = set()
        self._well_by_name: dict[str, WellSpec3D] = {}
        for w in self.wells:
            if w.name in seen_names:
                raise ValueError(f"Duplicate well name: {w.name!r}")
            if not (0 <= w.i < self.nx and 0 <= w.j < self.ny):
                raise ValueError(
                    f"Well {w.name!r} index ({w.i}, {w.j}) out of bounds for grid ({self.nx}, {self.ny})"
                )
            if not (0 <= w.k_min <= w.k_max < self.nz):
                raise ValueError(
                    f"Well {w.name!r} perforation interval [{w.k_min}, {w.k_max}] out of bounds for nz={self.nz}"
                )
            seen_names.add(w.name)
            self._well_by_name[w.name] = w

        # Binary H2 + cushion-gas composition (see _composition.py).
        self._init_composition()
        self.y: list[list[list[float]]] = [
            [[self._y_init] * props.nz for _ in range(props.ny)] for _ in range(props.nx)
        ]
        self._last_y_new = [[list(col) for col in plane] for plane in self.y]

        if props.k_layers is not None:
            if len(props.k_layers) != props.nz:
                raise ValueError(f"k_layers length {len(props.k_layers)} != nz={props.nz}")
            self._kz = [val * self._krg_eff for val in props.k_layers]
        else:
            self._kz = [props.k * self._krg_eff] * props.nz

    def _pv_mult(self, P: float) -> float:
        """Pore-volume multiplier exp(Ct·(P − P_init)) for rock (pore) compressibility Ct."""
        return math.exp(self.props.Ct * (P - self.props.P_init))

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
        rhoL = self._rho(P[i][j][k], T[i][j][k], self.y[i][j][k])
        rhoR = self._rho(P[i + 1][j][k], T[i + 1][j][k], self.y[i + 1][j][k])
        muL = self._mu(T[i][j][k], P[i][j][k], rhoL, self.y[i][j][k])
        muR = self._mu(T[i + 1][j][k], P[i + 1][j][k], rhoR, self.y[i + 1][j][k])
        rho_h = 2.0 * rhoL * rhoR / (rhoL + rhoR + 1e-30)
        mu_h = 2.0 * muL * muR / (muL + muR + 1e-30)
        return kh * self.dy * self.dz * rho_h / (mu_h * self.dx)

    def _Ty(
        self, P: list[list[list[float]]], T: list[list[list[float]]], i: int, j: int, k: int
    ) -> float:
        kh = self._k_harm_y(i, j, k)
        rhoB = self._rho(P[i][j][k], T[i][j][k], self.y[i][j][k])
        rhoT = self._rho(P[i][j + 1][k], T[i][j + 1][k], self.y[i][j + 1][k])
        muB = self._mu(T[i][j][k], P[i][j][k], rhoB, self.y[i][j][k])
        muT = self._mu(T[i][j + 1][k], P[i][j + 1][k], rhoT, self.y[i][j + 1][k])
        rho_h = 2.0 * rhoB * rhoT / (rhoB + rhoT + 1e-30)
        mu_h = 2.0 * muB * muT / (muB + muT + 1e-30)
        return kh * self.dx * self.dz * rho_h / (mu_h * self.dy)

    def _Tz(
        self, P: list[list[list[float]]], T: list[list[list[float]]], i: int, j: int, k: int
    ) -> float:
        kh = self._k_harm_z(i, j, k)
        rhoD = self._rho(P[i][j][k], T[i][j][k], self.y[i][j][k])
        rhoU = self._rho(P[i][j][k + 1], T[i][j][k + 1], self.y[i][j][k + 1])
        muD = self._mu(T[i][j][k], P[i][j][k], rhoD, self.y[i][j][k])
        muU = self._mu(T[i][j][k + 1], P[i][j][k + 1], rhoU, self.y[i][j][k + 1])
        rho_h = 2.0 * rhoD * rhoU / (rhoD + rhoU + 1e-30)
        mu_h = 2.0 * muD * muU / (muD + muU + 1e-30)
        return kh * self.dx * self.dy * rho_h / (mu_h * self.dz)

    def _well_index(
        self, w: WellSpec3D, P: list[list[list[float]]], T: list[list[list[float]]], k: int
    ) -> float:
        # Peaceman (1978) isotropic equivalent radius: 0.14·√(Δx² + Δy²) (≈ 0.198·Δx).
        r_eq = 0.14 * math.sqrt(self.dx**2 + self.dy**2)
        ln_r = math.log(max(r_eq, w.r_w * 1.001) / w.r_w)
        ln_r = max(ln_r, 0.1)
        k_w = self._k_cell(w.i, w.j, k)
        rho = self._rho(P[w.i][w.j][k], T[w.i][w.j][k], self.y[w.i][w.j][k])
        mu = self._mu(T[w.i][w.j][k], P[w.i][w.j][k], rho, self.y[w.i][w.j][k])
        return 2.0 * math.pi * k_w * self.dz * rho / (mu * ln_r)

    def _gravity_dp(
        self, P: list[list[list[float]]], T: list[list[list[float]]], i: int, j: int, k: int
    ) -> float:
        rhoD = self._rho(P[i][j][k], T[i][j][k], self.y[i][j][k])
        rhoU = self._rho(P[i][j][k + 1], T[i][j][k + 1], self.y[i][j][k + 1])
        rho_h = 2.0 * rhoD * rhoU / (rhoD + rhoU + 1e-30)
        return rho_h * self.g * self.dz

    def _assemble_well_contributions(
        self,
        ctrl: Schedule3D,
        Pk,
        Tk,
        rate_q: dict[str, float],
        shut_wells: set[str],
        nz: int,
    ) -> tuple:
        """Compute well diagonal/rhs contributions and well indices for a step."""
        well_diag = [[[0.0] * nz for _ in range(self.ny)] for _ in range(self.nx)]
        well_rhs = [[[0.0] * nz for _ in range(self.ny)] for _ in range(self.nx)]
        wi_now: dict[str, dict[int, float]] = {}
        for w in self.wells:
            wi_now[w.name] = {}
            for kl in range(w.k_min, w.k_max + 1):
                if 0 <= kl < nz:
                    wi_now[w.name][kl] = self._well_index(w, Pk, Tk, kl)
            mode_w = ctrl.mode_for(w.name)
            if mode_w == "shut_in" or w.name in shut_wells:
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
        return well_diag, well_rhs, wi_now

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
        y_k = [[list(col) for col in plane] for plane in self.y]
        accel = AndersonAccelerator(m=3) if self._is_multicomponent else None
        # The coupled pressure-composition map needs more iterations than single gas.
        n_picard = max(self.max_picard, 60) if self._is_multicomponent else self.max_picard

        for _picard in range(n_picard):
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
                            exponent=self._n_exp[i][j][k],
                        )
                        for k in range(nz)
                    ]
                    for j in range(ny)
                ]
                for i in range(nx)
            ]
            rho_k = [
                [
                    [self._rho(Pk[i][j][k], Tk[i][j][k], y_k[i][j][k]) for k in range(nz)]
                    for j in range(ny)
                ]
                for i in range(nx)
            ]
            drho = [
                [
                    [self._drho(Pk[i][j][k], Tk[i][j][k], y_k[i][j][k]) for k in range(nz)]
                    for j in range(ny)
                ]
                for i in range(nx)
            ]

            vp_k = [
                [[self.Vp * self._pv_mult(Pk[i][j][k]) for k in range(nz)] for j in range(ny)]
                for i in range(nx)
            ]
            B = [
                [
                    [
                        vp_k[i][j][k] * (drho[i][j][k] + rho_k[i][j][k] * self.props.Ct) / dt
                        for k in range(nz)
                    ]
                    for j in range(ny)
                ]
                for i in range(nx)
            ]
            rhs0 = [
                [
                    [
                        B[i][j][k] * Pk[i][j][k]
                        - (
                            vp_k[i][j][k] * rho_k[i][j][k]
                            - self.Vp * self._pv_mult(P_old[i][j][k]) * rho_old[i][j][k]
                        )
                        / dt
                        for k in range(nz)
                    ]
                    for j in range(ny)
                ]
                for i in range(nx)
            ]

            well_diag, well_rhs, wi_now = self._assemble_well_contributions(
                ctrl, Pk, Tk, rate_q, shut_wells, nz
            )

            if self.linear_solver == "pcg":
                P_new = self._solve_linear_pcg(B, rhs0, well_diag, well_rhs, tx, ty, tz, gz, Pk)
            else:
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
                else:
                    self._sor_failures += 1

            change = max(
                abs(P_new[i][j][k] - Pk[i][j][k])
                for i in range(nx)
                for j in range(ny)
                for k in range(nz)
            )
            self._last_picard_residual_Pa = change
            if accel is not None and change >= self.picard_tol_Pa:
                # Strong pressure-composition coupling: accelerate the Picard map.
                x_acc = accel.step(
                    [Pk[i][j][k] for i in range(nx) for j in range(ny) for k in range(nz)],
                    [P_new[i][j][k] for i in range(nx) for j in range(ny) for k in range(nz)],
                )
                P_new = [
                    [
                        [max(x_acc[(i * ny + j) * nz + k], 1.0e3) for k in range(nz)]
                        for j in range(ny)
                    ]
                    for i in range(nx)
                ]
            Pk = P_new

            if self._is_multicomponent:
                q_cells = self._perforation_rates(ctrl, Pk, wi_now, rate_q, shut_wells)
                y_next = self._solve_transport(P_old, Pk, Tk, rho_old, tx, ty, tz, gz, q_cells, dt)
                y_change = max(
                    abs(y_next[i][j][k] - y_k[i][j][k])
                    for i in range(nx)
                    for j in range(ny)
                    for k in range(nz)
                )
                y_k = y_next
                if y_change > 1e-9:
                    continue

            if change < self.picard_tol_Pa:
                break

        self._last_y_new = y_k

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
                        exponent=self._n_exp[i][j][k],
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
            mode_w = ctrl.mode_for(w.name)
            if mode_w == "shut_in" or w.name in shut_wells:
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

    def _layer_rates(self, w: WellSpec3D, P, T, bhp: float) -> dict[int, float]:
        """Signed mass rate [kg/s] of each perforation for a given bottom-hole pressure."""
        return {
            kl: self._well_index(w, P, T, kl) * (bhp - P[w.i][w.j][kl])
            for kl in range(w.k_min, w.k_max + 1)
        }

    def _produced_y(self, w: WellSpec3D, layer_q: dict[int, float]) -> float:
        """
        H2 mole fraction of the stream at a well: rate-weighted over producing
        perforations (molar basis), or the perforated-layer mean when not producing.
        """
        layers = range(w.k_min, w.k_max + 1)
        n_h2 = 0.0
        n_tot = 0.0
        for kl, q_l in layer_q.items():
            if q_l < 0.0:
                y_l = self.y[w.i][w.j][kl]
                mol = -q_l / self._mix_eos.molar_mass(self._comp_at(y_l))
                n_h2 += mol * y_l
                n_tot += mol
        if n_tot > 0.0:
            return n_h2 / n_tot
        return sum(self.y[w.i][w.j][kl] for kl in layers) / len(layers)

    def _perforation_rates(
        self, ctrl, Pk, wi_now, rate_q, shut_wells
    ) -> dict[tuple[int, int, int], float]:
        """Signed mass rate [kg/s] per perforated cell at the current pressure iterate."""
        q_cells: dict[tuple[int, int, int], float] = {}
        for w in self.wells:
            mode_w = ctrl.mode_for(w.name)
            layers = wi_now.get(w.name, {})
            if mode_w == "shut_in" or w.name in shut_wells or not layers:
                continue
            if w.name in rate_q:
                sum_wi = sum(layers.values())
                if sum_wi <= 0.0:
                    continue
                bhp = (
                    rate_q[w.name] + sum(wi * Pk[w.i][w.j][kl] for kl, wi in layers.items())
                ) / sum_wi
            else:
                bhp = ctrl.bhp_target
            for kl, wi in layers.items():
                key = (w.i, w.j, kl)
                q_cells[key] = q_cells.get(key, 0.0) + wi * (bhp - Pk[w.i][w.j][kl])
        return q_cells

    def _solve_transport(self, P_old, Pk, Tk, rho_old, tx, ty, tz, gz, q_cells, dt):
        """
        Implicit upwind advection–dispersion of the H2 mass fraction on the 3D grid
        (see :func:`solve_binary_transport`).  Vertical fluxes include gravity:
        F(k→k+1) = T_z·(P_k − P_{k+1} − ρ̄gΔz).  Returns H2 mole fractions.
        """
        nx, ny, nz = self.nx, self.ny, self.nz
        n = nx * ny * nz
        alpha_L = self.props.dispersivity_m
        D_mol = self.props.molecular_diffusivity_m2_s
        phi = self._phi_eff
        y = self.y

        def idx(i: int, j: int, k: int) -> int:
            return (i * ny + j) * nz + k

        w_old = [0.0] * n
        m_old = [0.0] * n
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    c = idx(i, j, k)
                    w_old[c] = self._w_from_y(y[i][j][k])
                    m_old[c] = max(
                        self.Vp * self._pv_mult(P_old[i][j][k]) * rho_old[i][j][k], 1e-12
                    )
        coup: list[list[tuple[int, float]]] = [[] for _ in range(n)]

        def link(c1, c2, F, P1, T1, y1, P2, T2, y2, area, dist):
            rho_avg = 0.5 * (self._rho(P1, T1, y1) + self._rho(P2, T2, y2))
            v = abs(F) / max(area * rho_avg, 1e-12)
            G = area * rho_avg * (D_mol + alpha_L * v) / dist
            add_transport_link(coup, c1, c2, F, G)

        ax, ay, az = self.dy * self.dz * phi, self.dx * self.dz * phi, self.dx * self.dy * phi
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    c = idx(i, j, k)
                    P1, T1, y1 = Pk[i][j][k], Tk[i][j][k], y[i][j][k]
                    if i < nx - 1:
                        F = tx[i][j][k] * (P1 - Pk[i + 1][j][k])
                        link(
                            c,
                            idx(i + 1, j, k),
                            F,
                            P1,
                            T1,
                            y1,
                            Pk[i + 1][j][k],
                            Tk[i + 1][j][k],
                            y[i + 1][j][k],
                            ax,
                            self.dx,
                        )
                    if j < ny - 1:
                        F = ty[i][j][k] * (P1 - Pk[i][j + 1][k])
                        link(
                            c,
                            idx(i, j + 1, k),
                            F,
                            P1,
                            T1,
                            y1,
                            Pk[i][j + 1][k],
                            Tk[i][j + 1][k],
                            y[i][j + 1][k],
                            ay,
                            self.dy,
                        )
                    if k < nz - 1:
                        F = tz[i][j][k] * (P1 - Pk[i][j][k + 1] - gz[i][j][k])
                        link(
                            c,
                            idx(i, j, k + 1),
                            F,
                            P1,
                            T1,
                            y1,
                            Pk[i][j][k + 1],
                            Tk[i][j][k + 1],
                            y[i][j][k + 1],
                            az,
                            self.dz,
                        )

        q_inj = [0.0] * n
        for (i, j, k), q in q_cells.items():
            q_inj[idx(i, j, k)] += max(q, 0.0)

        w_new, ok = solve_binary_transport(m_old, w_old, coup, q_inj, dt)
        if not ok:
            self._transport_failures += 1
        return [
            [[self._y_from_w(w_new[idx(i, j, k)]) for k in range(nz)] for j in range(ny)]
            for i in range(nx)
        ]

    def _solve_linear_pcg(self, B, rhs0, well_diag, well_rhs, tx, ty, tz, gz, Pk):
        """PCG solve of the 3D pressure system; index (i·ny + j)·nz + k."""
        nx, ny, nz = self.nx, self.ny, self.nz

        def idx(i: int, j: int, k: int) -> int:
            return (i * ny + j) * nz + k

        diag: list[float] = []
        nbrs: list[list[tuple[int, float]]] = []
        rhs: list[float] = []
        x0: list[float] = []
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    links: list[tuple[int, float]] = []
                    b_val = rhs0[i][j][k] + well_rhs[i][j][k]
                    if i > 0:
                        links.append((idx(i - 1, j, k), tx[i - 1][j][k]))
                    if i < nx - 1:
                        links.append((idx(i + 1, j, k), tx[i][j][k]))
                    if j > 0:
                        links.append((idx(i, j - 1, k), ty[i][j - 1][k]))
                    if j < ny - 1:
                        links.append((idx(i, j + 1, k), ty[i][j][k]))
                    # Gravity: flux k -> k+1 is T·(P_k − P_{k+1} − ρgΔz) (see SOR branch).
                    if k > 0:
                        links.append((idx(i, j, k - 1), tz[i][j][k - 1]))
                        b_val -= tz[i][j][k - 1] * gz[i][j][k - 1]
                    if k < nz - 1:
                        links.append((idx(i, j, k + 1), tz[i][j][k]))
                        b_val += tz[i][j][k] * gz[i][j][k]
                    diag.append(B[i][j][k] + sum(t for _, t in links) + well_diag[i][j][k])
                    nbrs.append(links)
                    rhs.append(b_val)
                    x0.append(Pk[i][j][k])
        x, ok = pcg_solve(diag, nbrs, rhs, x0, self.tol_Pa, self.max_iter)
        if not ok:
            self._sor_failures += 1
        return [
            [[max(x[idx(i, j, k)], 1.0e3) for k in range(nz)] for j in range(ny)] for i in range(nx)
        ]

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
        # Isentropic exponent at the start-of-step state, reused by every Picard iteration.
        self._n_exp = [
            [
                [
                    self._thermal_exponent(T_old[i][j][k], P_old[i][j][k], self.y[i][j][k])
                    for k in range(nz)
                ]
                for j in range(ny)
            ]
            for i in range(nx)
        ]
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
                [self._rho(P_old[i][j][k], T_old[i][j][k], self.y[i][j][k]) for k in range(nz)]
                for j in range(ny)
            ]
            for i in range(nx)
        ]

        rate_q: dict[str, float] = {}
        shut_wells: set[str] = set()
        for w in self.wells:
            mode_w = ctrl.mode_for(w.name)
            if mode_w == "shut_in":
                continue
            q_try = 0.0
            for kl in range(w.k_min, w.k_max + 1):
                if 0 <= kl < nz:
                    wi = self._well_index(w, P_old, T_old, kl)
                    q_try += wi * (ctrl.bhp_target - P_old[w.i][w.j][kl])
            if (mode_w == "inject" and q_try < 0.0) or (mode_w == "withdraw" and q_try > 0.0):
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
                mode_w = ctrl.mode_for(name)
                bhp_rc = per_well[name][1]
                binding = (mode_w == "inject" and bhp_rc <= ctrl.bhp_target) or (
                    mode_w == "withdraw" and bhp_rc >= ctrl.bhp_target
                )
                if not binding:
                    del rate_q[name]
                    changed = True

            for name in list(shut_wells):
                w = self._well_by_name[name]
                mode_w = ctrl.mode_for(name)
                q_try = sum(
                    self._well_index(w, P_new, T_new, kl) * (ctrl.bhp_target - P_new[w.i][w.j][kl])
                    for kl in range(w.k_min, w.k_max + 1)
                    if 0 <= kl < nz
                )
                if (mode_w == "inject" and q_try > 0.0) or (mode_w == "withdraw" and q_try < 0.0):
                    shut_wells.remove(name)
                    changed = True

            for w in self.wells:
                mode_w = ctrl.mode_for(w.name)
                if mode_w == "shut_in" or w.name in rate_q or w.name in shut_wells:
                    continue
                q_w = per_well[w.name][0]
                if (mode_w == "inject" and q_w < 0.0) or (mode_w == "withdraw" and q_w > 0.0):
                    shut_wells.add(w.name)
                    changed = True
                elif ctrl.q_max > 0.0 and abs(q_w) > ctrl.q_max * (1.0 + 1e-6):
                    rate_q[w.name] = math.copysign(ctrl.q_max, q_w)
                    changed = True

            if not changed:
                break
        else:
            self._active_set_failures += 1

        return P_new, T_new, per_well

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

    def _hydrostatic_init(self, P, nx, ny, nz):
        if not self.hydrostatic_init:
            return P
        y0 = self._y_init
        rho0 = self._rho(self.props.P_init, self.props.T_init, y0)
        z_ref = 0.5 * self.props.Lz
        p_col = [self.props.P_init + rho0 * self.g * (z_ref - self.zs[k]) for k in range(nz)]
        if nz > 1:
            for _ in range(10):
                dp = []
                for k in range(nz - 1):
                    rhoD = self._rho(p_col[k], self.props.T_init, y0)
                    rhoU = self._rho(p_col[k + 1], self.props.T_init, y0)
                    rho_h = 2.0 * rhoD * rhoU / (rhoD + rhoU + 1e-30)
                    dp.append(rho_h * self.g * self.dz)
                raw = [0.0] * nz
                for k in range(nz - 2, -1, -1):
                    raw[k] = raw[k + 1] + dp[k]
                mean_raw = sum(raw) / nz
                p_col = [self.props.P_init + (raw[k] - mean_raw) for k in range(nz)]
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    P[i][j][k] = p_col[k]
        return P

    def run_schedule(
        self,
        schedule: list[Schedule3D],
        snap_times_d: list[float] | None = None,
        dt_days: float = 1.0,
        verbose: bool = False,
    ) -> Result3D:
        nx, ny, nz = self.nx, self.ny, self.nz
        self._sor_failures = 0
        self._active_set_failures = 0
        self._transport_failures = 0
        self.y = [[[self._y_init] * nz for _ in range(ny)] for _ in range(nx)]
        self._last_y_new = [[list(col) for col in plane] for plane in self.y]
        P = [[[self.props.P_init] * nz for _ in range(ny)] for _ in range(nx)]
        P = self._hydrostatic_init(P, nx, ny, nz)
        T = [[[self.props.T_init] * nz for _ in range(ny)] for _ in range(nx)]
        t = 0.0

        result = Result3D()
        for w in self.wells:
            result.per_well_q[w.name] = []
            result.per_well_BHP[w.name] = []
            result.per_well_h2_mole_fraction[w.name] = []

        snap_targets = sorted(snap_times_d or [])
        snap_idx = 0
        cum_inj = 0.0
        cum_wdr = 0.0
        max_picard_res = 0.0
        hhv_h2, lhv_h2 = heating_values_mass()  # injected stream is pure H2
        e_wdr = [0.0, 0.0]  # J (HHV, LHV) withdrawn, at the produced layer compositions
        min_wdr_y_h2 = math.inf  # minimum over withdrawal steps only
        P_init_flat = [P[i][j][k] for i in range(nx) for j in range(ny) for k in range(nz)]
        P_max = max(P_init_flat) if P_init_flat else self.props.P_init
        P_min = min(P_init_flat) if P_init_flat else self.props.P_init
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
                self._rho(P[i][j][k], T[i][j][k], self.y[i][j][k])
                * self.Vp
                * self._pv_mult(P[i][j][k])
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
                result.per_well_h2_mole_fraction[w.name].append(
                    self._produced_y(w, self._last_layer_q.get(w.name, {}))
                )

        self._last_layer_q: dict[str, dict[int, float]] = {}
        _record(0.0, {w.name: (0.0, P[w.i][w.j][(w.k_min + w.k_max) // 2]) for w in self.wells})

        for ctrl in schedule:
            t_end = t + ctrl.duration_d * 86400.0
            n_steps = max(round(ctrl.duration_d / dt_days), 1)
            dt_act = (t_end - t) / n_steps

            if verbose:
                print(
                    f"  '{ctrl.name}': {ctrl.mode}, {ctrl.duration_d:.1f}d, BHP={pa_to_bar(ctrl.bhp_target):.1f}bar"
                )

            for _ in range(n_steps):
                P_old = [[list(P[i][j]) for j in range(ny)] for i in range(nx)]
                P, T, per_well = self._solve_pressure(P_old, T, ctrl, dt_act)
                # Perforation rates use the well indices of the solve (start-of-step
                # composition); produced gas then leaves at the new cell composition.
                layer_q = {
                    w.name: self._layer_rates(w, P, T, per_well[w.name][1]) for w in self.wells
                }
                self.y = self._last_y_new
                self._last_layer_q = layer_q
                max_picard_res = max(max_picard_res, self._last_picard_residual_Pa)
                t += dt_act

                for w in self.wells:
                    q_w = per_well[w.name][0]
                    if q_w >= 0.0:
                        continue
                    for kl, q_l in layer_q[w.name].items():
                        if q_l < 0.0:
                            hhv_p, lhv_p = heating_values_mass(self._comp_at(self.y[w.i][w.j][kl]))
                            e_wdr[0] += -q_l * dt_act * hhv_p
                            e_wdr[1] += -q_l * dt_act * lhv_p
                    if q_w < -1e-9:
                        min_wdr_y_h2 = min(min_wdr_y_h2, self._produced_y(w, layer_q[w.name]))

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

                t_d = t / 86400.0
                while snap_idx < len(snap_targets) and t_d >= snap_targets[snap_idx] - 1e-9:
                    result.snapshots.append(
                        Snapshot3D(
                            time_d=t_d,
                            xs_m=[
                                self.xs[i] for i in range(nx) for _ in range(ny) for _ in range(nz)
                            ],
                            ys_m=[
                                self.ys[j] for _ in range(nx) for j in range(ny) for _ in range(nz)
                            ],
                            zs_m=[
                                self.zs[k] for _ in range(nx) for _ in range(ny) for k in range(nz)
                            ],
                            P_Pa=P_flat,
                            T_K=T_flat,
                            h2_mole_fraction=[
                                self.y[i][j][k]
                                for i in range(nx)
                                for j in range(ny)
                                for k in range(nz)
                            ],
                        )
                    )
                    snap_idx += 1

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
            # Injected gas is pure H2; produced energy uses each perforation's composition.
            "total_injected_energy_GWh_LHV": J_to_GWh(cum_inj * lhv_h2),
            "total_injected_energy_GWh_HHV": J_to_GWh(cum_inj * hhv_h2),
            "total_withdrawn_energy_GWh_LHV": J_to_GWh(e_wdr[1]),
            "total_withdrawn_energy_GWh_HHV": J_to_GWh(e_wdr[0]),
            "working_gas_energy_GWh_LHV": J_to_GWh((m_max - m_min) * 1000.0 * lhv_h2),
            "working_gas_energy_GWh_HHV": J_to_GWh((m_max - m_min) * 1000.0 * hhv_h2),
            "mass_balance_error": mb_err,
            "max_picard_residual_Pa": max_picard_res,
            "sor_unconverged_solves": float(self._sor_failures),
            "active_set_unconverged_steps": float(self._active_set_failures),
            "transport_unconverged_solves": float(self._transport_failures),
            "initial_h2_purity_mol_pct": self._y_init * 100.0,
            "min_withdrawn_h2_purity_mol_pct": (
                min_wdr_y_h2 if math.isfinite(min_wdr_y_h2) else self._y_init
            )
            * 100.0,
            "mean_final_h2_mol_pct": 100.0
            * sum(self.y[i][j][k] for i in range(nx) for j in range(ny) for k in range(nz))
            / (nx * ny * nz),
            "converged": (
                1.0
                if max_picard_res < self.picard_tol_Pa
                and self._sor_failures == 0
                and self._active_set_failures == 0
                and self._transport_failures == 0
                else 0.0
            ),
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
