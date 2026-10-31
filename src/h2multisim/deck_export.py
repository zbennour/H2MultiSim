"""
h2multisim.deck_export
======================
Industry Benchmark Bridge: Generates standard ECLIPSE and OPM Flow (.DATA) input decks.

Allows head-to-head benchmarking of H2MultiSim simulations with open-source OPM Flow
(Open Porous Media) and commercial ECLIPSE 100/300 reservoir simulators.

Exports:
- RUNSPEC: gas-water metric simulation setup (water at connate saturation, immobile),
  timing, dimensions, and units.  OPM Flow has no single-phase gas mode.
- GRID: Cartesian or radial grid dimensions, tops, porosity, and anisotropic permeability.
- PROPS: NIST/PR-EOS derived PVDG gas PVT table (P, Bg, μ), PVTW, SGWFN gas-water
  relative permeability from the solver's RelativePermeabilityModel, rock compressibility,
  surface densities.
- REGIONS: Fluid-in-place zoning.
- SOLUTION: Initial pressure field and hydrostatic datum.
- SCHEDULE: Well specifications (WELSPECS), perforation intervals (COMPDAT), controls
  (WCONINJE / WCONPROD), and operational time-steps.

Zero external runtime dependencies.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

from h2multisim.models.one_d import RadialReservoirProps, WellControl
from h2multisim.models.three_d import Reservoir3DProps, Schedule3D, WellSpec3D
from h2multisim.models.two_d import Reservoir2DProps, Schedule2D, WellSpec2D
from h2multisim.thermo import PengRobinsonH2, RelativePermeabilityModel, viscosity_h2
from h2multisim.units import bar_to_pa, m2_to_mD, pa_to_bar

# Surface standard reference conditions (Metric: 1.01325 bar, 15.0 °C / 288.15 K)
P_STD_PA: float = 101325.0
T_STD_K: float = 288.15


def _rel_perm_of(props: object) -> RelativePermeabilityModel:
    """The two-phase model the H2MultiSim solvers build for ``props``."""
    return getattr(props, "rel_perm", None) or RelativePermeabilityModel(
        swr=float(getattr(props, "swr", 0.0)), sgc=float(getattr(props, "sgc", 0.0))
    )


def _initial_saturations(props: object) -> tuple[float, float]:
    """(S_w, S_g) the solvers assume: gas fills the gas-accessible pore fraction."""
    s_g = _rel_perm_of(props).effective_porosity_factor()
    return 1.0 - s_g, s_g


_DECK_HEADER_NOTE = (
    "-- Equivalence notes: gas-water deck with immobile connate water at the solver's\n"
    "-- S_wr (OPM Flow has no single-phase gas mode); SGWFN reproduces the solver's\n"
    "-- k_rg(S_g) with zero capillary pressure; ROCK carries the solver's pore\n"
    "-- compressibility Ct.  Isothermal: H2MultiSim's adiabatic/relaxation thermal\n"
    "-- model is not represented."
)

# RUNSPEC phases.  TABDIMS: NTSFUN NTPVT NSSFUN NPPVT.
_PHASE_LINES: tuple[str, ...] = ("WATER", "GAS", "METRIC", "", "TABDIMS", "  1 1 50 50 /")

# Brine PVTW: Pref [bar] Bw [rm3/sm3] Cw [1/bar] mu_w [cP] Cv [1/bar]
_PVTW_LINE = "  {p_ref:.2f} 1.0 4.5e-05 0.5 0.0 /"


def _props_lines(
    props: object,
    rho_std: float,
    pvdg: list[tuple[float, float, float]],
    p_init_bar: float,
    c_rock_1_bar: float,
    n_sat: int = 21,
) -> list[str]:
    """PROPS section body: densities, rock, water/gas PVT and the SGWFN table."""
    rel = _rel_perm_of(props)
    sg_max = 1.0 - rel.swr
    lines = [
        "DENSITY",
        f"  800.0 1000.0 {rho_std:.5f} /",
        "",
        "ROCK",
        f"  {p_init_bar:.2f} {c_rock_1_bar:.2e} /",
        "",
        "PVTW",
        _PVTW_LINE.format(p_ref=p_init_bar),
        "",
        "PVDG",
    ]
    lines.extend(f"  {p:8.2f} {bg:14.6e} {mu:10.5f}" for p, bg, mu in pvdg)
    lines.extend(["/", "", "SGWFN", "-- Sg        krg          krw          Pcgw"])
    sg_pts = sorted({sg_max * i / (n_sat - 1) for i in range(n_sat)} | {rel.sgc})
    for sg in sg_pts:
        krg = rel.krg(sg)
        krw = rel.krw(1.0 - sg)
        lines.append(f"  {sg:.6f} {krg:.6f} {krw:.6f} 0.0")
    lines.extend(["/", ""])
    return lines


# Field/well vectors needed to compare against H2MultiSim time series.
_SUMMARY_LINES: tuple[str, ...] = (
    "SUMMARY",
    "",
    "FPR",
    "FGIP",
    "FGIR",
    "FGPR",
    "FGIT",
    "FGPT",
    "WBHP",
    "/",
    "WGIR",
    "/",
    "WGPR",
    "/",
    "",
)


def _saturation_init_lines(props: object, n_cells: int) -> list[str]:
    s_w, s_g = _initial_saturations(props)
    return ["SWAT", f"  {n_cells}*{s_w:.6f} /", "SGAS", f"  {n_cells}*{s_g:.6f} /", ""]


def _wconprod_line(name: str, bhp_bar: float, q_sm3_d: float) -> str:
    # Items: well, status, cmode, ORAT, WRAT, GRAT, LRAT, RESV, BHP
    return f"  '{name}' 'OPEN' 'BHP' 2* {q_sm3_d:.1f} 2* {bhp_bar:.2f} /"


def generate_pvdg_table(
    eos: PengRobinsonH2 | None = None,
    T_res_K: float = 316.0,
    p_min_bar: float = 1.0,
    p_max_bar: float = 300.0,
    n_points: int = 30,
) -> tuple[float, list[tuple[float, float, float]]]:
    """
    Generate PVDG PVT table for ECLIPSE/OPM Flow:
    Returns (surface_density_kg_sm3, [(P_bar, Bg_rm3_sm3, visc_cP), ...]).
    """
    if eos is None:
        eos = PengRobinsonH2(volume_translation=-4.0e-6)

    rho_std = eos.density(P_STD_PA, T_STD_K)
    p_vals_bar = [
        p_min_bar + (p_max_bar - p_min_bar) * (i / (n_points - 1)) for i in range(n_points)
    ]

    table: list[tuple[float, float, float]] = []
    for p_bar in p_vals_bar:
        p_pa = bar_to_pa(p_bar)
        rho_res = eos.density(p_pa, T_res_K)
        bg = rho_std / max(rho_res, 1e-12)
        # Viscosity in cP: 1 Pa.s = 1000 cP
        mu_cp = viscosity_h2(T_res_K, p_pa, rho=rho_res, eos=eos) * 1000.0
        table.append((p_bar, bg, mu_cp))

    return rho_std, table


def export_eclipse_2d(
    props: Reservoir2DProps,
    wells: list[WellSpec2D],
    schedule: list[Schedule2D],
    filepath: str | Path,
    title: str = "H2MultiSim 2D Benchmark",
    eos: PengRobinsonH2 | None = None,
) -> None:
    """Export 2D areal reservoir scenario to an ECLIPSE / OPM Flow .DATA deck."""
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    if eos is None:
        vt = getattr(
            props, "volume_translation", -4.0e-6 if getattr(props, "peneloux", False) else 0.0
        )
        eos = PengRobinsonH2(volume_translation=vt)

    rho_std, pvdg = generate_pvdg_table(
        eos=eos,
        T_res_K=props.T_init,
        p_min_bar=max(pa_to_bar(props.P_init) * 0.2, 1.0),
        p_max_bar=pa_to_bar(props.P_init) * 2.0,
    )

    nx, ny = props.nx, props.ny
    dx = props.Lx / nx
    dy = props.Ly / ny
    p_init_bar = pa_to_bar(props.P_init)
    c_rock_1_bar = props.Ct * 1e5

    lines: list[str] = [
        "-- Generated by H2MultiSim Industry Benchmark Bridge",
        f"-- Case: {title}",
        _DECK_HEADER_NOTE,
        "",
        "RUNSPEC",
        "",
        "TITLE",
        f"  '{title}' /",
        "",
        "DIMENS",
        f"  {nx} {ny} 1 /",
        "",
        *_PHASE_LINES,
        "",
        "START",
        "  1 'JAN' 2026 /",
        "",
        "WELLDIMS",
        f"  {len(wells)} 10 1 {len(wells)} /",
        "",
        "UNIFOUT",
        "UNIFIN",
        "",
        "GRID",
        "",
        "DXV",
        f"  {nx}*{dx:.2f} /",
        "DYV",
        f"  {ny}*{dy:.2f} /",
        "DZV",
        f"  1*{props.h:.2f} /",
        "TOPS",
        f"  {nx * ny}*1000.0 /",
        "",
        "PORO",
        f"  {nx * ny}*{props.phi:.4f} /",
        "",
    ]

    # Permeability field
    if props.k_field is not None:
        perm_vals = [m2_to_mD(props.k_field[i][j]) for j in range(ny) for i in range(nx)]
        lines.append("PERMX")
        lines.append("  " + " ".join(f"{k:.2f}" for k in perm_vals) + " /")
        lines.append("PERMY")
        lines.append("  " + " ".join(f"{k:.2f}" for k in perm_vals) + " /")
    else:
        k_mD = m2_to_mD(props.k)
        lines.append("PERMX")
        lines.append(f"  {nx * ny}*{k_mD:.2f} /")
        lines.append("PERMY")
        lines.append(f"  {nx * ny}*{k_mD:.2f} /")

    # Single-layer grid: PERMZ has no effect on the flow solution.
    k_vert = m2_to_mD(props.k) * 0.1
    lines.extend(
        [
            "PERMZ",
            f"  {nx * ny}*{k_vert:.2f} /",
            "",
            "PROPS",
            "",
        ]
    )
    lines.extend(_props_lines(props, rho_std, pvdg, p_init_bar, c_rock_1_bar))
    lines.extend(
        [
            "REGIONS",
            "",
            "FIPNUM",
            f"  {nx * ny}*1 /",
            "",
            "SOLUTION",
            "",
            "PRESSURE",
            f"  {nx * ny}*{p_init_bar:.2f} /",
            "",
            *_saturation_init_lines(props, nx * ny),
            *_SUMMARY_LINES,
            "SCHEDULE",
            "",
            "RPTRST",
            "  'BASIC=2' /",
            "RPTSCHED",
            "  'PRES' 'SGAS' 'WELLS=2' /",
            "",
            "WELSPECS",
        ]
    )

    # Define wells (1-based I, J)
    for w in wells:
        lines.append(f"  '{w.name}' 'G1' {w.i + 1:3d} {w.j + 1:3d} 1000.0 'GAS' /")
    lines.extend(["/", "", "COMPDAT"])
    for w in wells:
        lines.append(f"  '{w.name}' {w.i + 1:3d} {w.j + 1:3d} 1 1 'OPEN' 2* {w.r_w * 2.0:.3f} /")
    lines.extend(["/", ""])

    # Schedule controls
    for period in schedule:
        lines.append(f"-- Operational Period: {period.name} ({period.duration_d:.1f} days)")
        for w in wells:
            mode_w = period.mode_for(w.name)
            bhp_bar = pa_to_bar(period.bhp_target)
            q_sm3_d = (period.q_max * 86400.0) / rho_std if period.q_max > 0.0 else 1.0e8

            if mode_w == "inject":
                lines.append("WCONINJE")
                lines.append(f"  '{w.name}' 'GAS' 'OPEN' 'BHP' {q_sm3_d:.1f} 1* {bhp_bar:.2f} /")
                lines.append("/")
            elif mode_w == "withdraw":
                lines.append("WCONPROD")
                lines.append(_wconprod_line(w.name, bhp_bar, q_sm3_d))
                lines.append("/")
            else:  # shut_in
                lines.append("WELOPEN")
                lines.append(f"  '{w.name}' 'SHUT' /")
                lines.append("/")

        lines.extend(
            [
                "TSTEP",
                f"  {period.duration_d:.2f} /",
                "",
            ]
        )

    lines.append("END")
    lines.append("")

    with open(filepath, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


def export_eclipse_3d(
    props: Reservoir3DProps,
    wells: list[WellSpec3D],
    schedule: list[Schedule3D],
    filepath: str | Path,
    title: str = "H2MultiSim 3D Benchmark",
    eos: PengRobinsonH2 | None = None,
) -> None:
    """Export 3D field-scale reservoir scenario to an ECLIPSE / OPM Flow .DATA deck."""
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    if eos is None:
        vt = getattr(
            props, "volume_translation", -4.0e-6 if getattr(props, "peneloux", False) else 0.0
        )
        eos = PengRobinsonH2(volume_translation=vt)

    rho_std, pvdg = generate_pvdg_table(
        eos=eos,
        T_res_K=props.T_init,
        p_min_bar=max(pa_to_bar(props.P_init) * 0.2, 1.0),
        p_max_bar=pa_to_bar(props.P_init) * 2.0,
    )

    nx, ny, nz = props.nx, props.ny, props.nz
    dx = props.Lx / nx
    dy = props.Ly / ny
    dz = props.Lz / nz
    p_init_bar = pa_to_bar(props.P_init)
    c_rock_1_bar = props.Ct * 1e5
    z_datum = 1000.0 + props.Lz * 0.5  # Mid-depth datum matching H2MultiSim

    # H2MultiSim layer k = 0 is the BOTTOM layer; ECLIPSE K = 1 is the TOP layer.
    def _deck_k(k_solver: int) -> int:
        return nz - k_solver

    lines: list[str] = [
        "-- Generated by H2MultiSim Industry Benchmark Bridge",
        f"-- Case: {title}",
        _DECK_HEADER_NOTE,
        "",
        "RUNSPEC",
        "",
        "TITLE",
        f"  '{title}' /",
        "",
        "DIMENS",
        f"  {nx} {ny} {nz} /",
        "",
        *_PHASE_LINES,
        "",
        "START",
        "  1 'JAN' 2026 /",
        "",
        "WELLDIMS",
        f"  {len(wells)} 10 1 {len(wells)} /",
        "",
        "UNIFOUT",
        "UNIFIN",
        "",
        "GRID",
        "",
        "DXV",
        f"  {nx}*{dx:.2f} /",
        "DYV",
        f"  {ny}*{dy:.2f} /",
        "DZV",
        f"  {nz}*{dz:.2f} /",
        "TOPS",
        f"  {nx * ny}*1000.0 /",
        "",
        "PORO",
        f"  {nx * ny * nz}*{props.phi:.4f} /",
        "",
    ]

    # Layered or homogeneous permeability
    # The 3D solver is isotropic (k_v = k_h), so PERMZ = PERMX.
    if props.k_layers is not None:
        k_x_vals: list[str] = []
        for k in reversed(range(nz)):  # deck order: top layer first
            kl_mD = m2_to_mD(props.k_layers[k])
            k_x_vals.append(f"  {nx * ny}*{kl_mD:.2f}")
        for kw in ("PERMX", "PERMY", "PERMZ"):
            lines.append(kw)
            lines.extend(k_x_vals)
            lines.append("  /")
    else:
        k_mD = m2_to_mD(props.k)
        for kw in ("PERMX", "PERMY", "PERMZ"):
            lines.append(kw)
            lines.append(f"  {nx * ny * nz}*{k_mD:.2f} /")

    lines.extend(
        [
            "",
            "PROPS",
            "",
        ]
    )
    lines.extend(_props_lines(props, rho_std, pvdg, p_init_bar, c_rock_1_bar))
    lines.extend(
        [
            "REGIONS",
            "",
            "FIPNUM",
            f"  {nx * ny * nz}*1 /",
            "",
            "SOLUTION",
            "",
            "EQUIL",
            "-- datum_depth   datum_pressure   woc_depth   woc_Pcap   goc_depth   goc_Pcap",
            f"   {z_datum:.2f}     {p_init_bar:.2f}        2000.0      0.0        0.0         0.0 /",
            "",
            *_SUMMARY_LINES,
            "SCHEDULE",
            "",
            "RPTRST",
            "  'BASIC=2' /",
            "RPTSCHED",
            "  'PRES' 'SGAS' 'WELLS=2' /",
            "",
            "WELSPECS",
        ]
    )

    # Define wells (1-based I, J)
    for w in wells:
        lines.append(f"  '{w.name}' 'G1' {w.i + 1:3d} {w.j + 1:3d} {z_datum:.1f} 'GAS' /")
    lines.extend(["/", "", "COMPDAT"])
    for w in wells:
        # Perforations: 1-based K1 (upper) <= K2 (lower) in deck numbering
        lines.append(
            f"  '{w.name}' {w.i + 1:3d} {w.j + 1:3d} {_deck_k(w.k_max):3d} {_deck_k(w.k_min):3d} 'OPEN' 2* {w.r_w * 2.0:.3f} /"
        )
    lines.extend(["/", ""])

    # Schedule controls
    for period in schedule:
        lines.append(f"-- Operational Period: {period.name} ({period.duration_d:.1f} days)")
        for w in wells:
            mode_w = period.mode_for(w.name)
            bhp_bar = pa_to_bar(period.bhp_target)
            q_sm3_d = (period.q_max * 86400.0) / rho_std if period.q_max > 0.0 else 1.0e8

            if mode_w == "inject":
                lines.append("WCONINJE")
                lines.append(f"  '{w.name}' 'GAS' 'OPEN' 'BHP' {q_sm3_d:.1f} 1* {bhp_bar:.2f} /")
                lines.append("/")
            elif mode_w == "withdraw":
                lines.append("WCONPROD")
                lines.append(_wconprod_line(w.name, bhp_bar, q_sm3_d))
                lines.append("/")
            else:  # shut_in
                lines.append("WELOPEN")
                lines.append(f"  '{w.name}' 'SHUT' /")
                lines.append("/")

        lines.extend(
            [
                "TSTEP",
                f"  {period.duration_d:.2f} /",
                "",
            ]
        )

    lines.append("END")
    lines.append("")

    with open(filepath, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


def export_eclipse_1d_radial(
    props: RadialReservoirProps,
    schedule: list[WellControl],
    filepath: str | Path,
    N: int = 40,
    title: str = "H2MultiSim 1D Radial Benchmark",
    eos: PengRobinsonH2 | None = None,
) -> None:
    """
    Export a 1D radial scenario to an ECLIPSE / OPM Flow .DATA deck.

    The radial finite-volume model is written as an N x 1 x 1 Cartesian strip whose
    pore volumes (``PORV``) and inter-cell transmissibilities (``TRANX``) are set in
    the EDIT section to the exact values of H2MultiSim's logarithmic radial grid,
    and whose well connection factor uses the same equivalent radius.  This avoids
    simulator-specific ``RADIAL`` grid handling (OPM Flow builds the cylindrical
    grid but does not connect the well to it).
    """
    os.makedirs(os.path.dirname(os.path.abspath(filepath)), exist_ok=True)
    if eos is None:
        vt = getattr(
            props, "volume_translation", -4.0e-6 if getattr(props, "peneloux", False) else 0.0
        )
        eos = PengRobinsonH2(volume_translation=vt)

    rho_std, pvdg = generate_pvdg_table(
        eos=eos,
        T_res_K=props.T_init,
        p_min_bar=max(pa_to_bar(props.P_init) * 0.2, 1.0),
        p_max_bar=pa_to_bar(props.P_init) * 2.0,
    )

    from h2multisim.models.one_d import build_log_grid

    r_c, r_i = build_log_grid(props.r_well, props.r_ext, N)
    # Well connection factor [cP·rm³/day/bar], same equivalent radius as UHSRadial1D:
    # CF = 2π·C·k·h / max(ln(r_c0/r_w) + S, 0.1), C = 0.008527 (METRIC Darcy constant).
    ln_term = max(math.log(max(r_c[0], props.r_well * 1.001) / props.r_well) + props.skin, 0.1)
    cf_well = 2.0 * math.pi * 0.008527 * m2_to_mD(props.k) * props.h / ln_term
    dr_vals = [r_i[i + 1] - r_i[i] for i in range(N)]
    k_mD_raw = m2_to_mD(props.k)
    # Exact radial FV quantities (UHSRadial1D): V_i = π(r_{i+1/2}² − r_{i−1/2}²)·h,
    # T_{i,i+1} = 2π·C·k·h / ln(r_{c,i+1}/r_{c,i}).
    porv = [math.pi * (r_i[i + 1] ** 2 - r_i[i] ** 2) * props.h * props.phi for i in range(N)]
    tranx = [
        2.0 * math.pi * 0.008527 * k_mD_raw * props.h / math.log(r_c[i + 1] / r_c[i])
        for i in range(N - 1)
    ] + [0.0]
    p_init_bar = pa_to_bar(props.P_init)
    c_rock_1_bar = props.Ct * 1e5
    k_mD = m2_to_mD(props.k)

    lines: list[str] = [
        "-- Generated by H2MultiSim Industry Benchmark Bridge",
        f"-- Case: {title}",
        _DECK_HEADER_NOTE,
        "",
        "RUNSPEC",
        "",
        "TITLE",
        f"  '{title}' /",
        "",
        "DIMENS",
        f"  {N} 1 1 /",
        "",
        *_PHASE_LINES,
        "",
        "START",
        "  1 'JAN' 2026 /",
        "",
        "WELLDIMS",
        "  1 10 1 1 /",
        "",
        "UNIFOUT",
        "UNIFIN",
        "",
        "GRID",
        "",
        "DXV",
        "  " + " ".join(f"{dr:.4f}" for dr in dr_vals) + " /",
        "DYV",
        "  1.0 /",
        "DZV",
        f"  1*{props.h:.2f} /",
        "",
        "TOPS",
        f"  {N}*1000.0 /",
        "",
        "PORO",
        f"  {N}*{props.phi:.4f} /",
        "",
        "PERMX",
        f"  {N}*{k_mD:.2f} /",
        "PERMY",
        f"  {N}*{k_mD:.2f} /",
        "PERMZ",
        f"  {N}*{k_mD:.2f} /",
        "",
        "EDIT",
        "",
        "-- Exact radial pore volumes [rm3] and transmissibilities [cP.rm3/day/bar]",
        "PORV",
        *(f"  {v:.6e}" for v in porv),
        "/",
        "TRANX",
        *(f"  {t:.6e}" for t in tranx),
        "/",
        "",
        "PROPS",
        "",
    ]
    lines.extend(_props_lines(props, rho_std, pvdg, p_init_bar, c_rock_1_bar))
    lines.extend(
        [
            "REGIONS",
            "",
            "FIPNUM",
            f"  {N}*1 /",
            "",
            "SOLUTION",
            "",
            "PRESSURE",
            f"  {N}*{p_init_bar:.2f} /",
            "",
            *_saturation_init_lines(props, N),
            *_SUMMARY_LINES,
            "SCHEDULE",
            "",
            "RPTRST",
            "  'BASIC=2' /",
            "RPTSCHED",
            "  'PRES' 'SGAS' 'WELLS=2' /",
            "",
            "WELSPECS",
            "  'W1' 'G1' 1 1 1000.0 'GAS' /",
            "/",
            "",
            "COMPDAT",
            # Items: well I J K1 K2 status SATNUM CF diameter Kh skin.  The connection
            # factor is explicit: on a radial grid OPM's Peaceman estimate of r_o can fall
            # inside the wellbore and the well becomes non-finite.
            f"  'W1' 1 1 1 1 'OPEN' 1* {cf_well:.6e} {props.r_well * 2.0:.3f} 1* {props.skin:.3f} /",
            "/",
            "",
        ]
    )

    for ctrl in schedule:
        lines.append(f"-- Operational Period: {ctrl.name} ({ctrl.duration_d:.1f} days)")
        bhp_bar = pa_to_bar(ctrl.bhp_target)
        q_sm3_d = (ctrl.q_max * 86400.0) / rho_std if ctrl.q_max > 0.0 else 1.0e8

        if ctrl.mode == "inject":
            lines.append("WCONINJE")
            lines.append(f"  'W1' 'GAS' 'OPEN' 'BHP' {q_sm3_d:.1f} 1* {bhp_bar:.2f} /")
            lines.append("/")
        elif ctrl.mode == "withdraw":
            lines.append("WCONPROD")
            lines.append(_wconprod_line("W1", bhp_bar, q_sm3_d))
            lines.append("/")
        else:
            lines.append("WELOPEN")
            lines.append("  'W1' 'SHUT' /")
            lines.append("/")

        lines.extend(
            [
                "TSTEP",
                f"  {ctrl.duration_d:.2f} /",
                "",
            ]
        )

    lines.append("END")
    lines.append("")

    with open(filepath, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
