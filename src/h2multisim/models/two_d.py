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

from h2multisim.models._common import AndersonAccelerator, check_inputs, pcg_solve, validate_mode
from h2multisim.thermo import (
    GAS_COMPONENTS,
    M_H2,
    PengRobinsonH2,
    PengRobinsonMixture,
    RelativePermeabilityModel,
    coupled_temperature_update,
    heating_values_mass,
    isentropic_exponent,
    viscosity_h2,
)
from h2multisim.units import (
    J_to_GWh,
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
    Ct: float  # Pa⁻¹ pore (rock) compressibility; gas compressibility comes from the EOS
    P_init: float  # Pa initial reservoir pressure
    T_init: float  # K initial / geothermal temperature
    nx: int  # grid cells in x
    ny: int  # grid cells in y
    tau_days: float = 30.0  # days thermal relaxation constant
    k_field: list[list[float]] | None = None
    volume_translation: float = 0.0  # m³/mol Peneloux volume translation (e.g. -4.0e-6)
    peneloux: bool = False  # If True and volume_translation == 0, sets volume_translation = -4.0e-6
    swr: float = 0.0  # [-] connate / residual immobile water saturation S_wr
    sgc: float = 0.0  # [-] critical / trapped gas saturation S_gc
    rel_perm: RelativePermeabilityModel | None = (
        None  # Optional two-phase relative permeability model
    )
    cushion_gas: str = "H2"  # Native cushion gas species ('H2', 'CH4', 'CO2', 'N2')
    initial_h2_mole_fraction: float | None = None  # None -> 1.0 if cushion is H2 else 0.0
    dispersivity_m: float = 0.5  # m longitudinal mechanical dispersivity α_L
    molecular_diffusivity_m2_s: float = 1.0e-5  # m²/s effective pore molecular diffusivity

    def __post_init__(self) -> None:
        check_inputs(
            "Reservoir2DProps",
            positive={
                "Lx": self.Lx,
                "Ly": self.Ly,
                "h": self.h,
                "nx": float(self.nx),
                "ny": float(self.ny),
                "k": self.k,
                "P_init": self.P_init,
                "T_init": self.T_init,
            },
            non_negative={
                "Ct": self.Ct,
                "tau_days": self.tau_days,
                "dispersivity_m": self.dispersivity_m,
                "molecular_diffusivity_m2_s": self.molecular_diffusivity_m2_s,
            },
            fraction={"phi": self.phi},
        )
        self.cushion_gas = self.cushion_gas.strip().upper()
        if self.cushion_gas not in GAS_COMPONENTS:
            raise ValueError(
                f"Reservoir2DProps: unknown cushion_gas {self.cushion_gas!r}; "
                f"choose from {list(GAS_COMPONENTS)}"
            )
        if self.initial_h2_mole_fraction is not None and not (
            0.0 <= self.initial_h2_mole_fraction <= 1.0
        ):
            raise ValueError("Reservoir2DProps: initial_h2_mole_fraction must be in [0, 1]")
        if self.peneloux and self.volume_translation == 0.0:
            self.volume_translation = -4.0e-6


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


SchedulePeriod2D = Schedule2D


@dataclass
class Snapshot2D:
    """Areal spatial snapshot at a given time."""

    time_d: float
    xs_m: list[float]
    ys_m: list[float]
    P_Pa: list[float]
    T_K: list[float]
    h2_mole_fraction: list[float] = field(default_factory=list)


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
    per_well_h2_mole_fraction: dict[str, list[float]] = field(default_factory=dict)
    snapshots: list[Snapshot2D] = field(default_factory=list)
    metrics: dict[str, float] = field(default_factory=dict)


SimulationResult2D = Result2D


class UHSAreal2D:
    """
    Implicit 2D areal finite-volume UHS simulator.
    Implicit pressure (PCG or SOR) with Peaceman wells, optional binary H2 + cushion-gas
    composition transport (implicit upwind advection–dispersion).
    """

    def __init__(
        self,
        props: Reservoir2DProps,
        wells: list[WellSpec2D],
        eos: PengRobinsonH2 | None = None,
        sor_omega: float = 1.5,
        max_iter: int = 2000,
        tol_Pa: float = 1.0e-4,
        max_picard: int = 20,
        picard_tol_Pa: float = 1.0e-4,
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
        self.tau_s = props.tau_days * 86400.0
        self._last_picard_residual_Pa: float = 0.0
        self._sor_failures: int = 0  # inner linear solves (PCG or SOR) that hit max_iter
        self._active_set_failures: int = 0  # steps whose well active set did not settle
        self._n_exp: list = []  # per-cell isentropic exponent for the current step

        self.nx = props.nx
        self.ny = props.ny
        self.dx = props.Lx / props.nx
        self.dy = props.Ly / props.ny
        self.xs = [(i + 0.5) * self.dx for i in range(self.nx)]
        self.ys = [(j + 0.5) * self.dy for j in range(self.ny)]

        self.rel_perm: RelativePermeabilityModel = getattr(
            props, "rel_perm", None
        ) or RelativePermeabilityModel(
            swr=float(getattr(props, "swr", 0.0)),
            sgc=float(getattr(props, "sgc", 0.0)),
        )
        self._phi_eff: float = props.phi * self.rel_perm.effective_porosity_factor()
        self._krg_eff: float = self.rel_perm.effective_permeability_factor()
        self.Vp = self._phi_eff * self.dx * self.dy * props.h

        seen_names: set[str] = set()
        self._well_by_name: dict[str, WellSpec2D] = {}
        for w in self.wells:
            if w.name in seen_names:
                raise ValueError(f"Duplicate well name: {w.name!r}")
            if not (0 <= w.i < self.nx and 0 <= w.j < self.ny):
                raise ValueError(
                    f"Well {w.name!r} index ({w.i}, {w.j}) out of bounds for grid ({self.nx}, {self.ny})"
                )
            seen_names.add(w.name)
            self._well_by_name[w.name] = w

        # Binary H2 + cushion-gas composition (mirrors UHSRadial1D).
        cg = props.cushion_gas
        self._cushion_gas: str = cg
        y_raw = props.initial_h2_mole_fraction
        self._y_init: float = float(y_raw) if y_raw is not None else (1.0 if cg == "H2" else 0.0)
        self._is_multicomponent: bool = not (cg == "H2" or self._y_init >= 1.0 - 1e-12)
        self._mix_eos = PengRobinsonMixture(
            composition={"H2": 1.0},
            peneloux=props.peneloux,
            volume_translation=props.volume_translation or None,
        )
        self.y: list[list[float]] = [[self._y_init] * props.ny for _ in range(props.nx)]
        self._last_y_new: list[list[float]] = [list(row) for row in self.y]
        self._transport_failures: int = 0

        if props.k_field is not None:
            if len(props.k_field) != props.nx or any(len(row) != props.ny for row in props.k_field):
                raise ValueError(f"k_field dimensions must match (nx={props.nx}, ny={props.ny})")
            self.k = [[val * self._krg_eff for val in row] for row in props.k_field]
        else:
            self.k = [[props.k * self._krg_eff] * props.ny for _ in range(props.nx)]

    def _pv_mult(self, P: float) -> float:
        """Pore-volume multiplier exp(Ct·(P − P_init)) for rock (pore) compressibility Ct."""
        return math.exp(self.props.Ct * (P - self.props.P_init))

    # ── Composition-aware cell properties (single-gas path uses self.eos directly) ──
    def _comp_at(self, y_val: float) -> dict[str, float]:
        if not self._is_multicomponent or y_val >= 1.0 - 1e-12:
            return {"H2": 1.0}
        if y_val <= 1e-12:
            return {self._cushion_gas: 1.0}
        return {"H2": y_val, self._cushion_gas: 1.0 - y_val}

    def _rho(self, P: float, T: float, y_val: float) -> float:
        if not self._is_multicomponent:
            return self.eos.density(P, T)
        return self._mix_eos.density(P, T, self._comp_at(y_val))

    def _drho(self, P: float, T: float, y_val: float) -> float:
        if not self._is_multicomponent:
            return self.eos.drho_dP(P, T)
        return self._mix_eos.drho_dP(P, T, self._comp_at(y_val))

    def _mu(self, T: float, P: float, rho: float, y_val: float) -> float:
        if not self._is_multicomponent:
            return viscosity_h2(T, P, rho=rho, eos=self.eos)
        return self._mix_eos.viscosity(P, T, rho=rho, composition=self._comp_at(y_val))

    def _exponent(self, T: float, P: float, y_val: float) -> float:
        if not self._is_multicomponent:
            return isentropic_exponent(T, P, self.eos)
        self._mix_eos.set_composition(self._comp_at(y_val))
        return isentropic_exponent(T, P, self._mix_eos)

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
        y = self.y
        rhoL = self._rho(P[i][j], T[i][j], y[i][j])
        rhoR = self._rho(P[i + 1][j], T[i + 1][j], y[i + 1][j])
        muL = self._mu(T[i][j], P[i][j], rhoL, y[i][j])
        muR = self._mu(T[i + 1][j], P[i + 1][j], rhoR, y[i + 1][j])
        rho_h = 2.0 * rhoL * rhoR / (rhoL + rhoR + 1e-30)
        mu_h = 2.0 * muL * muR / (muL + muR + 1e-30)
        return kh * self.props.h * self.dy * rho_h / (mu_h * self.dx)

    def _Ty(self, i: int, j: int, P: list[list[float]], T: list[list[float]]) -> float:
        kh = self._k_harm_y(i, j)
        y = self.y
        rhoB = self._rho(P[i][j], T[i][j], y[i][j])
        rhoT = self._rho(P[i][j + 1], T[i][j + 1], y[i][j + 1])
        muB = self._mu(T[i][j], P[i][j], rhoB, y[i][j])
        muT = self._mu(T[i][j + 1], P[i][j + 1], rhoT, y[i][j + 1])
        rho_h = 2.0 * rhoB * rhoT / (rhoB + rhoT + 1e-30)
        mu_h = 2.0 * muB * muT / (muB + muT + 1e-30)
        return kh * self.props.h * self.dx * rho_h / (mu_h * self.dy)

    def _well_index(self, w: WellSpec2D, P: list[list[float]], T: list[list[float]]) -> float:
        # Peaceman (1978) isotropic equivalent radius: 0.14·√(Δx² + Δy²) (≈ 0.198·Δx).
        r_eq = 0.14 * math.sqrt(self.dx**2 + self.dy**2)
        ln_r = math.log(max(r_eq, w.r_w * 1.001) / w.r_w)
        ln_r = max(ln_r, 0.1)
        k_w = self.k[w.i][w.j]
        y_w = self.y[w.i][w.j]
        rho = self._rho(P[w.i][w.j], T[w.i][w.j], y_w)
        mu = self._mu(T[w.i][w.j], P[w.i][w.j], rho, y_w)
        return 2.0 * math.pi * k_w * self.props.h * rho / (mu * ln_r)

    def _assemble_well_contributions_2d(self, ctrl, Pk, Tk, rate_q, shut_wells, nx, ny):
        well_diag = [[0.0] * ny for _ in range(nx)]
        well_rhs = [[0.0] * ny for _ in range(nx)]
        wi_now: dict[str, float] = {}
        for w in self.wells:
            wi = self._well_index(w, Pk, Tk)
            wi_now[w.name] = wi
            mode_w = ctrl.mode_for(w.name)
            if mode_w == "shut_in" or w.name in shut_wells:
                continue
            if w.name in rate_q:
                well_rhs[w.i][w.j] += rate_q[w.name]
            else:
                well_diag[w.i][w.j] += wi
                well_rhs[w.i][w.j] += wi * ctrl.bhp_target
        return well_diag, well_rhs, wi_now

    def _solve_picard(
        self,
        P_old: list[list[float]],
        T_old: list[list[float]],
        ctrl: Schedule2D,
        dt: float,
        tx: list[list[float]],
        ty: list[list[float]],
        rho_old: list[list[float]],
        rate_q: dict[str, float],
        shut_wells: set[str],
    ) -> tuple[list[list[float]], list[list[float]], dict[str, tuple[float, float]]]:
        """Picard iteration on pressure and temperature for one time step."""
        nx, ny = self.nx, self.ny
        T_geo = self.props.T_init
        Pk = [list(row) for row in P_old]
        Tk = [list(row) for row in T_old]
        y_k = [list(row) for row in self.y]
        wi_final: dict[str, float] = {w.name: 0.0 for w in self.wells}
        accel = AndersonAccelerator(m=3) if self._is_multicomponent else None
        # The coupled pressure-composition map needs more iterations than single gas.
        n_picard = max(self.max_picard, 60) if self._is_multicomponent else self.max_picard

        for _picard in range(n_picard):
            Tk = [
                [
                    coupled_temperature_update(
                        T_old=T_old[i][j],
                        P_old=max(P_old[i][j], 1.0),
                        P_new=max(Pk[i][j], 1.0),
                        T_geothermal=T_geo,
                        dt_s=dt,
                        tau_s=self.tau_s,
                        eos=self.eos,
                        exponent=self._n_exp[i][j],
                    )
                    for j in range(ny)
                ]
                for i in range(nx)
            ]
            rho_k = [
                [self._rho(Pk[i][j], Tk[i][j], y_k[i][j]) for j in range(ny)] for i in range(nx)
            ]
            drho = [
                [self._drho(Pk[i][j], Tk[i][j], y_k[i][j]) for j in range(ny)] for i in range(nx)
            ]

            vp_k = [[self.Vp * self._pv_mult(Pk[i][j]) for j in range(ny)] for i in range(nx)]
            B = [
                [vp_k[i][j] * (drho[i][j] + rho_k[i][j] * self.props.Ct) / dt for j in range(ny)]
                for i in range(nx)
            ]
            rhs0 = [
                [
                    B[i][j] * Pk[i][j]
                    - (
                        vp_k[i][j] * rho_k[i][j]
                        - self.Vp * self._pv_mult(P_old[i][j]) * rho_old[i][j]
                    )
                    / dt
                    for j in range(ny)
                ]
                for i in range(nx)
            ]

            well_diag, well_rhs, wi_now = self._assemble_well_contributions_2d(
                ctrl, Pk, Tk, rate_q, shut_wells, nx, ny
            )

            if self.linear_solver == "pcg":
                P_new = self._solve_linear_pcg(B, rhs0, well_diag, well_rhs, tx, ty, Pk)
            else:
                P_new = self._solve_linear_sor(B, rhs0, well_diag, well_rhs, tx, ty, Pk)

            change = max(abs(P_new[i][j] - Pk[i][j]) for i in range(nx) for j in range(ny))
            self._last_picard_residual_Pa = change
            if accel is not None and change >= self.picard_tol_Pa:
                # Pressure and composition are strongly coupled (H2 is ~10x lighter than
                # the cushion gas): accelerate the linearly convergent Picard map.
                x_acc = accel.step(
                    [Pk[i][j] for i in range(nx) for j in range(ny)],
                    [P_new[i][j] for i in range(nx) for j in range(ny)],
                )
                P_new = [[max(x_acc[i * ny + j], 1.0e3) for j in range(ny)] for i in range(nx)]
            Pk = P_new

            if self._is_multicomponent:
                q_est = {}
                for w in self.wells:
                    mode_w = ctrl.mode_for(w.name)
                    if mode_w == "shut_in" or w.name in shut_wells:
                        q_est[w.name] = 0.0
                    elif w.name in rate_q:
                        q_est[w.name] = rate_q[w.name]
                    else:
                        q_est[w.name] = wi_now[w.name] * (ctrl.bhp_target - Pk[w.i][w.j])
                y_next = self._solve_transport(P_old, Pk, Tk, rho_old, tx, ty, q_est, dt)
                y_change = max(abs(y_next[i][j] - y_k[i][j]) for i in range(nx) for j in range(ny))
                y_k = y_next
                if y_change > 1e-9:
                    continue

            if change < self.picard_tol_Pa:
                break

        self._last_y_new = y_k

        # Make temperature exactly consistent with the converged pressure.
        Tk = [
            [
                coupled_temperature_update(
                    T_old=T_old[i][j],
                    P_old=max(P_old[i][j], 1.0),
                    P_new=max(Pk[i][j], 1.0),
                    T_geothermal=T_geo,
                    dt_s=dt,
                    tau_s=self.tau_s,
                    eos=self.eos,
                    exponent=self._n_exp[i][j],
                )
                for j in range(ny)
            ]
            for i in range(nx)
        ]
        for w in self.wells:
            wi_final[w.name] = self._well_index(w, Pk, Tk)

        per_well: dict[str, tuple[float, float]] = {}
        for w in self.wells:
            wi = wi_final[w.name]
            mode_w = ctrl.mode_for(w.name)
            P_w = Pk[w.i][w.j]
            if mode_w == "shut_in" or w.name in shut_wells:
                q = 0.0
            elif w.name in rate_q:
                q = rate_q[w.name]
            else:
                q = wi * (ctrl.bhp_target - P_w)
            bhp = P_w + q / max(wi, 1e-30)
            per_well[w.name] = (q, bhp)

        return Pk, Tk, per_well

    def _solve_transport(self, P_old, Pk, Tk, rho_old, tx, ty, q_well, dt):
        """
        Implicit upwind advection–dispersion of the H2 mass fraction (advective
        form, as in UHSRadial1D): injected gas is pure H2; production removes the
        well-cell composition.  The matrix is a diagonally dominant M-matrix, so
        Gauss–Seidel converges.  Returns updated H2 mole fractions.
        """
        nx, ny = self.nx, self.ny
        M_cg = GAS_COMPONENTS[self._cushion_gas].molar_mass
        w_old = [
            [
                self.y[i][j] * M_H2 / max(self.y[i][j] * M_H2 + (1.0 - self.y[i][j]) * M_cg, 1e-30)
                for j in range(ny)
            ]
            for i in range(nx)
        ]
        m_old = [
            [max(self.Vp * self._pv_mult(P_old[i][j]) * rho_old[i][j], 1e-12) for j in range(ny)]
            for i in range(nx)
        ]
        alpha_L = self.props.dispersivity_m
        D_mol = self.props.molecular_diffusivity_m2_s
        h, phi = self.props.h, self._phi_eff

        # Per-cell list of (neighbour i, j, coupling) where coupling = inflow + dispersion.
        coup: list[list[list[tuple[int, int, float]]]] = [
            [[] for _ in range(ny)] for _ in range(nx)
        ]

        def link(i1, j1, i2, j2, T_link, area, dist):
            F = T_link * (Pk[i1][j1] - Pk[i2][j2])  # kg/s from cell 1 to cell 2
            rho_avg = 0.5 * (
                self._rho(Pk[i1][j1], Tk[i1][j1], self.y[i1][j1])
                + self._rho(Pk[i2][j2], Tk[i2][j2], self.y[i2][j2])
            )
            v = abs(F) / max(area * rho_avg, 1e-12)
            G = area * rho_avg * (D_mol + alpha_L * v) / dist
            coup[i2][j2].append((i1, j1, max(F, 0.0) + G))
            coup[i1][j1].append((i2, j2, max(-F, 0.0) + G))

        for i in range(nx - 1):
            for j in range(ny):
                link(i, j, i + 1, j, tx[i][j], self.dy * h * phi, self.dx)
        for i in range(nx):
            for j in range(ny - 1):
                link(i, j, i, j + 1, ty[i][j], self.dx * h * phi, self.dy)

        q_inj = [[0.0] * ny for _ in range(nx)]
        for w in self.wells:
            q_inj[w.i][w.j] += max(q_well.get(w.name, 0.0), 0.0)

        w_new = [list(row) for row in w_old]
        for _ in range(500):
            max_dw = 0.0
            for i in range(nx):
                for j in range(ny):
                    base = m_old[i][j] / dt
                    diag = base + q_inj[i][j]
                    rhs = base * w_old[i][j] + q_inj[i][j]
                    for i2, j2, c in coup[i][j]:
                        diag += c
                        rhs += c * w_new[i2][j2]
                    w_val = min(max(rhs / diag, 0.0), 1.0)
                    max_dw = max(max_dw, abs(w_val - w_new[i][j]))
                    w_new[i][j] = w_val
            if max_dw < 1e-12:
                break
        else:
            self._transport_failures += 1

        y_new = []
        for i in range(nx):
            row = []
            for j in range(ny):
                n_h2 = w_new[i][j] / M_H2
                n_cg = (1.0 - w_new[i][j]) / M_cg
                row.append(n_h2 / max(n_h2 + n_cg, 1e-30))
            y_new.append(row)
        return y_new

    def _linear_system_2d(self, B, rhs0, well_diag, well_rhs, tx, ty):
        """Flat (diag, neighbours, rhs) of the 2D pressure system, index i·ny + j."""
        nx, ny = self.nx, self.ny
        diag: list[float] = []
        nbrs: list[list[tuple[int, float]]] = []
        rhs: list[float] = []
        for i in range(nx):
            for j in range(ny):
                links: list[tuple[int, float]] = []
                if i > 0:
                    links.append(((i - 1) * ny + j, tx[i - 1][j]))
                if i < nx - 1:
                    links.append(((i + 1) * ny + j, tx[i][j]))
                if j > 0:
                    links.append((i * ny + j - 1, ty[i][j - 1]))
                if j < ny - 1:
                    links.append((i * ny + j + 1, ty[i][j]))
                diag.append(B[i][j] + sum(t for _, t in links) + well_diag[i][j])
                nbrs.append(links)
                rhs.append(rhs0[i][j] + well_rhs[i][j])
        return diag, nbrs, rhs

    def _solve_linear_pcg(self, B, rhs0, well_diag, well_rhs, tx, ty, Pk):
        nx, ny = self.nx, self.ny
        diag, nbrs, rhs = self._linear_system_2d(B, rhs0, well_diag, well_rhs, tx, ty)
        x0 = [Pk[i][j] for i in range(nx) for j in range(ny)]
        x, ok = pcg_solve(diag, nbrs, rhs, x0, self.tol_Pa, self.max_iter)
        if not ok:
            self._sor_failures += 1
        return [[max(x[i * ny + j], 1.0e3) for j in range(ny)] for i in range(nx)]

    def _solve_linear_sor(self, B, rhs0, well_diag, well_rhs, tx, ty, Pk):
        nx, ny = self.nx, self.ny
        P_new = [list(row) for row in Pk]
        for _iter in range(self.max_iter):
            max_change = 0.0
            for i in range(nx):
                for j in range(ny):
                    TxL = tx[i - 1][j] if i > 0 else 0.0
                    TxR = tx[i][j] if i < nx - 1 else 0.0
                    TyB = ty[i][j - 1] if j > 0 else 0.0
                    TyT = ty[i][j] if j < ny - 1 else 0.0

                    diag = B[i][j] + TxL + TxR + TyB + TyT + well_diag[i][j]
                    rhs = rhs0[i][j] + well_rhs[i][j]

                    if i > 0:
                        rhs += TxL * P_new[i - 1][j]
                    if i < nx - 1:
                        rhs += TxR * P_new[i + 1][j]
                    if j > 0:
                        rhs += TyB * P_new[i][j - 1]
                    if j < ny - 1:
                        rhs += TyT * P_new[i][j + 1]

                    P_gs = rhs / max(diag, 1e-30)
                    P_sor = (1.0 - self.omega) * P_new[i][j] + self.omega * P_gs
                    P_sor = max(P_sor, 1.0e3)

                    change = abs(P_sor - P_new[i][j])
                    max_change = max(max_change, change)
                    P_new[i][j] = P_sor

            if max_change < self.tol_Pa:
                break
        else:
            self._sor_failures += 1
        return P_new

    def _solve_pressure(
        self,
        P_old: list[list[float]],
        T_old: list[list[float]],
        ctrl: Schedule2D,
        dt: float,
    ) -> tuple[list[list[float]], list[list[float]], dict[str, tuple[float, float]]]:
        """
        Solve one pressure/temperature step with an active-set well model.

        Wells are first solved against their BHP targets; wells whose rate would
        exceed ``q_max`` are switched to rate control, and wells that would flow
        backwards are shut in.  The active set is re-evaluated until stable.  No
        cell pressure is overwritten, so rates, BHPs and mass are consistent.
        """
        nx, ny = self.nx, self.ny
        # Isentropic exponent at the start-of-step state, reused by every Picard iteration.
        self._n_exp = [
            [self._exponent(T_old[i][j], P_old[i][j], self.y[i][j]) for j in range(ny)]
            for i in range(nx)
        ]
        tx = [[self._Tx(i, j, P_old, T_old) for j in range(ny)] for i in range(nx - 1)]
        ty = [[self._Ty(i, j, P_old, T_old) for j in range(ny - 1)] for i in range(nx)]
        rho_old = [
            [self._rho(P_old[i][j], T_old[i][j], self.y[i][j]) for j in range(ny)]
            for i in range(nx)
        ]

        rate_q: dict[str, float] = {}
        shut_wells: set[str] = set()
        for w in self.wells:
            mode_w = ctrl.mode_for(w.name)
            if mode_w == "shut_in":
                continue
            wi = self._well_index(w, P_old, T_old)
            P_w = P_old[w.i][w.j]
            q_try = wi * (ctrl.bhp_target - P_w)
            if (mode_w == "inject" and q_try < 0.0) or (mode_w == "withdraw" and q_try > 0.0):
                shut_wells.add(w.name)
            elif ctrl.q_max > 0.0 and abs(q_try) > ctrl.q_max:
                rate_q[w.name] = math.copysign(ctrl.q_max, q_try)

        P_new = [list(row) for row in P_old]
        T_new = [list(row) for row in T_old]
        per_well: dict[str, tuple[float, float]] = {}
        for _attempt in range(5):
            P_new, T_new, per_well = self._solve_picard(
                P_old, T_old, ctrl, dt, tx, ty, rho_old, rate_q, shut_wells
            )
            changed = False

            # A capped well whose BHP constraint is satisfied is not actually
            # rate-limited: revert it to BHP control.
            for name in list(rate_q):
                w = self._well_by_name[name]
                mode_w = ctrl.mode_for(name)
                P_w = P_new[w.i][w.j]
                binding = (mode_w == "inject" and P_w <= ctrl.bhp_target) or (
                    mode_w == "withdraw" and P_w >= ctrl.bhp_target
                )
                if not binding:
                    del rate_q[name]
                    changed = True

            # Re-open shut wells if the pressure has moved back to the side
            # where the well can flow in its intended direction.
            for name in list(shut_wells):
                w = self._well_by_name[name]
                mode_w = ctrl.mode_for(name)
                P_w = P_new[w.i][w.j]
                if (mode_w == "inject" and P_w < ctrl.bhp_target) or (
                    mode_w == "withdraw" and P_w > ctrl.bhp_target
                ):
                    shut_wells.remove(name)
                    changed = True

            # No-backflow and post-solve rate-cap constraints for BHP-controlled wells.
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

    def run_schedule(
        self,
        schedule: list[Schedule2D],
        dt_days: float = 0.2,
        snap_days: list[float] | None = None,
    ) -> Result2D:
        nx, ny = self.nx, self.ny

        self._sor_failures = 0
        self._active_set_failures = 0
        self._transport_failures = 0
        self.y = [[self._y_init] * ny for _ in range(nx)]
        self._last_y_new = [list(row) for row in self.y]
        P = [[self.props.P_init] * ny for _ in range(nx)]
        T = [[self.props.T_init] * ny for _ in range(nx)]

        result = Result2D()
        for w in self.wells:
            result.per_well_q[w.name] = []
            result.per_well_BHP[w.name] = []
            result.per_well_h2_mole_fraction[w.name] = []

        snap_targets = sorted(snap_days or [])
        snap_idx = 0
        t_d = 0.0
        cum_inj = 0.0
        cum_wdr = 0.0
        max_picard_res = 0.0
        hhv_h2, lhv_h2 = heating_values_mass()  # injected stream is pure H2
        e_wdr = [0.0, 0.0]  # J (HHV, LHV) withdrawn, at the produced composition
        min_wdr_y_h2 = math.inf  # minimum over withdrawal steps only

        P_max = self.props.P_init
        P_min = self.props.P_init
        T_max = self.props.T_init
        T_min = self.props.T_init
        m_initial = (
            self._rho(self.props.P_init, self.props.T_init, self._y_init) * self.Vp * nx * ny
        )

        for ctrl in schedule:
            n_steps = max(round(ctrl.duration_d / dt_days), 1)
            dt_act = ctrl.duration_d * 86400.0 / n_steps
            for _ in range(n_steps):
                P_new, T_new, per_well = self._solve_pressure(P, T, ctrl, dt_act)
                self.y = self._last_y_new
                max_picard_res = max(max_picard_res, self._last_picard_residual_Pa)
                t_d += dt_act / 86400.0

                q_inj = sum(max(q, 0.0) for q, _ in per_well.values())
                q_wdr = sum(max(-q, 0.0) for q, _ in per_well.values())
                cum_inj += q_inj * dt_act
                cum_wdr += q_wdr * dt_act
                total_q = q_inj - q_wdr

                for w in self.wells:
                    q_w, bhp_w = per_well[w.name]
                    y_w = self.y[w.i][w.j]
                    result.per_well_q[w.name].append(q_w)
                    result.per_well_BHP[w.name].append(pa_to_bar(bhp_w))
                    result.per_well_h2_mole_fraction[w.name].append(y_w)
                    if q_w < 0.0:
                        hhv_p, lhv_p = heating_values_mass(self._comp_at(y_w))
                        e_wdr[0] += -q_w * dt_act * hhv_p
                        e_wdr[1] += -q_w * dt_act * lhv_p
                        if q_w < -1e-9:
                            min_wdr_y_h2 = min(min_wdr_y_h2, y_w)

                P = P_new
                T = T_new

                P_flat = [P[i][j] for i in range(nx) for j in range(ny)]
                T_flat = [T[i][j] for i in range(nx) for j in range(ny)]

                P_max = max(P_max, max(P_flat))
                P_min = min(P_min, min(P_flat))
                T_max = max(T_max, max(T_flat))
                T_min = min(T_min, min(T_flat))

                m_tot = sum(
                    self._rho(P[i][j], T[i][j], self.y[i][j]) * self.Vp * self._pv_mult(P[i][j])
                    for i in range(nx)
                    for j in range(ny)
                )

                result.times_d.append(t_d)
                result.P_avg_bar.append(pa_to_bar(sum(P_flat) / len(P_flat)))
                result.T_avg_C.append(K_to_C(sum(T_flat) / len(T_flat)))
                result.mass_t.append(kg_to_tonnes(m_tot))
                result.q_total_kg_s.append(total_q)

                while snap_idx < len(snap_targets) and t_d >= snap_targets[snap_idx] - 1e-9:
                    result.snapshots.append(
                        Snapshot2D(
                            time_d=t_d,
                            xs_m=[self.xs[i] for i in range(nx) for _ in range(ny)],
                            ys_m=[self.ys[j] for _ in range(nx) for j in range(ny)],
                            P_Pa=P_flat,
                            T_K=T_flat,
                            h2_mole_fraction=[self.y[i][j] for i in range(nx) for j in range(ny)],
                        )
                    )
                    snap_idx += 1

        m_max = max(result.mass_t) if result.mass_t else 0.0
        m_min = min(result.mass_t) if result.mass_t else 0.0
        m_final_kg = result.mass_t[-1] * 1000.0 if result.mass_t else m_initial
        net_kg = cum_inj - cum_wdr
        mb_err = abs((m_final_kg - m_initial) - net_kg) / max(abs(net_kg), 1.0)
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
            "initial_mass_t": kg_to_tonnes(m_initial),
            "final_mass_t": kg_to_tonnes(m_final_kg),
            "total_injected_t": kg_to_tonnes(cum_inj),
            "total_withdrawn_t": kg_to_tonnes(cum_wdr),
            "net_stored_t": kg_to_tonnes(net_kg),
            "cycled_gas_t": kg_to_tonnes(cum_inj + cum_wdr),
            # Injected gas is pure H2; produced energy uses the well-cell composition.
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
            "mean_final_h2_mol_pct": 100.0 * sum(map(sum, self.y)) / (nx * ny),
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


def _five_spot_schedule(P_ref: float) -> list[Schedule2D]:
    """Five-spot schedule: four corner injectors, one central producer.

    The original scenario applied a single mode to all five wells, which made
    the "producer" inject and the "injectors" produce.  Per-well overrides keep
    the pattern physically meaningful while preserving the period API.
    """
    bar = 1.0e5
    P_max = min(1.35 * P_ref, P_ref + 50.0 * bar)
    P_min = max(0.60 * P_ref, P_ref - 40.0 * bar)
    qmax = 35.0
    names = ["INJ_SW", "INJ_NW", "INJ_SE", "INJ_NE", "CTRE"]
    inj_modes = {n: ("inject" if n.startswith("INJ") else "shut_in") for n in names}
    wdr_modes = {n: ("withdraw" if n == "CTRE" else "shut_in") for n in names}
    shut_modes = {n: "shut_in" for n in names}
    return [
        Schedule2D("inject_30d", "inject", 30.0, P_max, qmax, well_modes=inj_modes),
        Schedule2D("withdraw_30d", "withdraw", 30.0, P_min, qmax, well_modes=wdr_modes),
        Schedule2D("shut_in_10d", "shut_in", 10.0, 0.0, 0.0, well_modes=shut_modes),
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
    sched = _five_spot_schedule(props.P_init)
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
    sched = _five_spot_schedule(props.P_init)
    return props, wells, sched


def compute_sweep_efficiency(
    P_field: list[list[float]] | list[float],
    P_init: float,
    threshold_Pa: float = 1.0e4,
) -> float:
    """
    Volumetric sweep efficiency: fraction of cells where |P - P_init| > threshold_Pa.
    Accepts either a 2D grid [nx][ny] or a flat 1D sequence of cell pressures.
    """
    if not P_field:
        return 0.0
    if isinstance(P_field[0], (int, float)):
        vals = [float(p) for p in P_field]
    else:
        vals = [float(p) for row in P_field for p in row]  # type: ignore[union-attr]
    if not vals:
        return 0.0
    n_swept = sum(1 for p in vals if abs(p - P_init) > threshold_Pa)
    return float(n_swept / len(vals))


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
