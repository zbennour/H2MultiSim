"""Flow-model audit: mass balance, timestep/grid convergence, analytic and
mass-conservative reference comparisons."""

import json
import math
import sys
import time

import numpy as np

import pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from h2multisim.models.one_d import (  # noqa: E402
    RadialReservoirProps,
    UHSRadial1D,
    WellControl,
    build_log_grid,
    okoroafor_radial,
)
from h2multisim.models.three_d import (  # noqa: E402
    UHSFieldScale3D,
    scenario_layered_reservoir,
    scenario_single_well_3d,
)
from h2multisim.models.two_d import (  # noqa: E402
    UHSAreal2D,
    scenario_five_spot,
    scenario_single_well,
)
from h2multisim.models.zero_d import (  # noqa: E402
    ReservoirState,
    SchedulePeriod,
    UHSReservoir0D,
    build_okoroafor_reservoir,
    build_teesside_cavern,
    schedule_cavern,
    schedule_porous_reservoir,
)
from h2multisim.thermo import PengRobinsonH2  # noqa: E402
from h2multisim.units import bar_to_pa, mD_to_m2  # noqa: E402

EOS = PengRobinsonH2()
OUT = {}

def mb(tag, m0, m1, inj, wdr):
    net = inj - wdr
    dm = m1 - m0
    err = (dm - net) / max(abs(net), 1e-9) * 100.0
    print(f"  {tag:<26} M0={m0/1e3:12.1f} t  M1={m1/1e3:12.1f} t  dM={(dm)/1e3:11.1f} t "
          f" net_q={net/1e3:11.1f} t  MB_err={err:10.2f}%")
    return err


print("=" * 100)
print("1. MATERIAL BALANCE: EOS-derived inventory change vs. reported well throughput")
print("=" * 100)
OUT["mass_balance"] = {}

# ---- 0D cavern
props = build_teesside_cavern()
m = UHSReservoir0D(props)
s0 = ReservoirState(props.reference_pressure, props.reference_temperature,
                    m.mass_from_state(props.reference_pressure, props.reference_temperature))
r = m.run_schedule(s0, schedule_cavern(props.reference_pressure), 86400.0)
e = mb("0D Teesside cavern", s0.mass, r.masses[-1], r.metrics["total_injected_kg"],
       r.metrics["total_withdrawn_kg"])
OUT["mass_balance"]["0d_teesside"] = e
print(f"    reported 'working gas' = {r.metrics['working_gas_kg']/1e3:.1f} t, "
      f"reported injected = {r.metrics['total_injected_kg']/1e3:.0f} t")

# ---- 0D porous
props = build_okoroafor_reservoir()
m = UHSReservoir0D(props)
s0 = ReservoirState(props.reference_pressure, props.reference_temperature,
                    m.mass_from_state(props.reference_pressure, props.reference_temperature))
r = m.run_schedule(s0, schedule_porous_reservoir(props.reference_pressure), 86400.0)
e = mb("0D Okoroafor porous", s0.mass, r.masses[-1], r.metrics["total_injected_kg"],
       r.metrics["total_withdrawn_kg"])
OUT["mass_balance"]["0d_okoroafor"] = e

# ---- 1D
props, sched = okoroafor_radial()
sim = UHSRadial1D(props, N=40)
r1 = sim.run_schedule(sched, dt_days=0.2)
m0 = sum(EOS.density(props.P_init, props.T_init) * v for v in sim.Vp)
m1 = r1.metrics["final_mass_kg"]
e = mb("1D radial (N=40)", m0, m1, r1.metrics["total_injected_kg"], r1.metrics["total_withdrawn_kg"])
OUT["mass_balance"]["1d"] = e

# ---- 2D
p2, w2, s2 = scenario_five_spot(15, 15)
sim2 = UHSAreal2D(p2, w2)
t0 = time.time()
r2 = sim2.run_schedule(s2, dt_days=0.5)
t2d = time.time() - t0
m0 = EOS.density(p2.P_init, p2.T_init) * sim2.Vp * p2.nx * p2.ny
m1 = r2.mass_t[-1] * 1e3
e = mb("2D five-spot 15x15", m0, m1, r2.metrics["total_injected_t"] * 1e3,
       r2.metrics["total_withdrawn_t"] * 1e3)
OUT["mass_balance"]["2d"] = e
print(f"    reported working gas = {r2.metrics['working_gas_t']:.0f} t "
      f"(= {r2.metrics['working_gas_t']/max(r2.metrics['total_injected_t'],1e-9):.2f} x total injected)")

# ---- 3D
p3, w3, s3 = scenario_single_well_3d(8, 8, 5)
sim3 = UHSFieldScale3D(p3, w3)
t0 = time.time()
r3 = sim3.run_schedule(s3, dt_days=2.0)
t3d = time.time() - t0
m0 = EOS.density(p3.P_init, p3.T_init) * sim3.Vp * p3.nx * p3.ny * p3.nz
m1 = r3.mass_t[-1] * 1e3
e = mb("3D single-well 8x8x5", m0, m1, r3.metrics["total_injected_t"] * 1e3,
       r3.metrics["total_withdrawn_t"] * 1e3)
OUT["mass_balance"]["3d"] = e
print(f"    reported working gas = {r3.metrics['working_gas_t']:.0f} t "
      f"(= {r3.metrics['working_gas_t']/max(r3.metrics['total_injected_t'],1e-9):.2f} x total injected)")

# ---------------------------------------------------------------------------
print("\n" + "=" * 100)
print("2. TIMESTEP SENSITIVITY (should converge as dt -> 0)")
print("=" * 100)

print("\n0D Teesside (working gas [t], max P [bar]):")
for dt in (2.0, 1.0, 0.5, 0.25):
    props = build_teesside_cavern()
    mm = UHSReservoir0D(props)
    st = ReservoirState(props.reference_pressure, props.reference_temperature,
                        mm.mass_from_state(props.reference_pressure, props.reference_temperature))
    rr = mm.run_schedule(st, schedule_cavern(props.reference_pressure), dt * 86400.0)
    print(f"  dt={dt:4.2f} d  WG={rr.metrics['working_gas_kg']/1e3:9.2f}  "
          f"Pmax={rr.metrics['max_pressure_Pa']/1e5:7.2f}  "
          f"Pmin={rr.metrics['min_pressure_Pa']/1e5:7.2f}")

print("\n3D single-well (working gas [t], max P [bar]):")
for dt in (3.0, 2.0, 1.0, 0.5):
    p3, w3, s3 = scenario_single_well_3d(8, 8, 5)
    rr = UHSFieldScale3D(p3, w3).run_schedule(s3, dt_days=dt)
    print(f"  dt={dt:4.2f} d  WG={rr.metrics['working_gas_t']:9.0f}  "
          f"Pmax={rr.metrics['max_P_bar']:7.2f}  Pmin={rr.metrics['min_P_bar']:7.2f}")

print("\n2D five-spot (working gas [t], max P [bar]):")
for dt in (1.0, 0.5, 0.25):
    p2, w2, s2 = scenario_five_spot(15, 15)
    rr = UHSAreal2D(p2, w2).run_schedule(s2, dt_days=dt)
    print(f"  dt={dt:4.2f} d  WG={rr.metrics['working_gas_t']:9.0f}  "
          f"Pmax={rr.metrics['max_P_bar']:7.2f}  Pmin={rr.metrics['min_P_bar']:7.2f}")

print("\n" + "=" * 100)
print("3. GRID SENSITIVITY")
print("=" * 100)
print("\n1D radial BHP at end of 30 d injection (bar):")
for N in (10, 20, 40, 80):
    props, sched = okoroafor_radial()
    rr = UHSRadial1D(props, N=N).run_schedule(sched, dt_days=0.25)
    print(f"  N={N:3d}  BHP_end={rr.metrics['final_BHP_Pa']/1e5:8.3f}  "
          f"inj={rr.metrics['total_injected_kg']/1e3:8.2f} t")

print("\n2D single-well working gas (t) vs grid:")
for nx in (10, 15, 20, 30):
    p2, w2, s2 = scenario_single_well(nx, nx)
    rr = UHSAreal2D(p2, w2).run_schedule(s2, dt_days=0.5)
    print(f"  {nx}x{nx}  WG={rr.metrics['working_gas_t']:10.0f}  Pmax={rr.metrics['max_P_bar']:7.2f}")

print("\n" + "=" * 100)
print("4. SOR CONVERGENCE: default (2000 it, 1e-4 Pa) vs tight (4000 it, 1e-3 Pa)")
print("=" * 100)
p2, w2, s2 = scenario_five_spot(15, 15)
a = UHSAreal2D(p2, w2).run_schedule(s2, dt_days=0.5)
b = UHSAreal2D(p2, w2, max_iter=4000, tol_Pa=1e-3).run_schedule(s2, dt_days=0.5)
print(f"  default : WG={a.metrics['working_gas_t']:10.1f}  Pmax={a.metrics['max_P_bar']:8.3f}  "
      f"Pmin={a.metrics['min_P_bar']:8.3f}")
print(f"  tight   : WG={b.metrics['working_gas_t']:10.1f}  Pmax={b.metrics['max_P_bar']:8.3f}  "
      f"Pmin={b.metrics['min_P_bar']:8.3f}")
print(f"  difference: {abs(a.metrics['working_gas_t']-b.metrics['working_gas_t'])/b.metrics['working_gas_t']*100:.2f}%")

print("\n" + "=" * 100)
print("5. REFERENCE: mass-conservative single-phase real-gas radial simulator (same grid, same PR-EOS)")
print("=" * 100)


def reference_radial(props, rate, days, N=40, nsteps=120, T=316.0):
    """Backward-Euler FV, mass-conservative in kg, Picard iteration.

    rate > 0 = injection, rate < 0 = withdrawal [kg/s], applied to cell 0.
    Accumulation uses d(rho)/dP from the same PR-EOS (no Ct, no double phi).
    """
    from h2multisim.thermo import viscosity_h2

    r_c, r_i = build_log_grid(props.r_well, props.r_ext, N)
    V = np.array([math.pi * (r_i[i + 1] ** 2 - r_i[i] ** 2) * props.h for i in range(N)])
    Vp = props.phi * V
    P = np.full(N, props.P_init)
    rho = np.array([EOS.density(p, T) for p in P])
    mass0 = float(np.sum(rho * Vp))
    dt = days * 86400.0 / nsteps
    net = 0.0

    def trans(Pk):
        rk = np.array([EOS.density(p, T) for p in Pk])
        mk = np.array([viscosity_h2(T, p) for p in Pk])
        Tf = np.zeros(N - 1)
        for i in range(N - 1):
            rh = 2 * rk[i] * rk[i + 1] / (rk[i] + rk[i + 1])
            mh = 2 * mk[i] * mk[i + 1] / (mk[i] + mk[i + 1])
            Tf[i] = 2 * math.pi * props.k * props.h * rh / (mh * math.log(r_i[i + 2] / r_i[i + 1]))
        return Tf

    for _ in range(nsteps):
        Pk = P.copy()
        Tr = trans(Pk)
        for _it in range(60):
            rk = np.array([EOS.density(p, T) for p in Pk])
            dP_ = 1.0e3
            cp = np.array([
                (EOS.density(p + dP_, T) - EOS.density(max(p - dP_, 1.0), T)) / (2 * dP_)
                for p in Pk
            ])
            A = Vp * cp / dt
            diag = A.copy()
            diag[:-1] += Tr
            diag[1:] += Tr
            rhs = A * Pk - (rk - rho) * Vp / dt
            rhs[0] += rate
            # solve tridiagonal: -Tr P[i-1] + diag P[i] - Tr P[i+1] = rhs
            a = np.zeros(N); b = diag.copy(); c = np.zeros(N)
            a[1:] = -Tr; c[:-1] = -Tr
            # Thomas
            cc = np.zeros(N); dd = np.zeros(N)
            cc[0] = c[0] / b[0]; dd[0] = rhs[0] / b[0]
            for i in range(1, N):
                den = b[i] - a[i] * cc[i - 1]
                cc[i] = c[i] / den
                dd[i] = (rhs[i] - a[i] * dd[i - 1]) / den
            x = np.zeros(N); x[-1] = dd[-1]
            for i in range(N - 2, -1, -1):
                x[i] = dd[i] - cc[i] * x[i + 1]
            change = np.max(np.abs(x - Pk))
            Pk = 0.5 * Pk + 0.5 * x
            if change < 1e-2:
                break
        P = Pk
        rho = np.array([EOS.density(p, T) for p in P])
        net += rate * dt
    rho_f = np.array([EOS.density(p, T) for p in P])
    mass1 = float(np.sum(rho_f * Vp))
    return P, r_c, mass0, mass1, net


# constant-rate withdrawal, same reservoir as okoroafor_radial
props, _ = okoroafor_radial()
q = 3.0  # kg/s
P_ref, r_ref, m0r, m1r, injr = reference_radial(props, -q, days=10, N=40, nsteps=80)
print(f"  reference: M0={m0r/1e3:.1f} t M1={m1r/1e3:.1f} t dM={(m1r-m0r)/1e3:.3f} t net_q={injr/1e3:.3f} t "
      f"MB_err={(m1r-m0r-injr)/abs(injr)*100:.2e}%")
print(f"  reference BHP after 10 d at {q} kg/s = {P_ref[0]/1e5:.2f} bar, outer P = {P_ref[-1]/1e5:.2f} bar")

# code, same condition: bhp target far below, q_max = 3 kg/s -> constant rate
props2 = RadialReservoirProps(
    r_well=props.r_well, r_ext=props.r_ext, h=props.h, phi=props.phi, k=props.k,
    Ct=props.Ct, P_init=props.P_init, T_init=props.T_init, skin=0.0,
    outer_bc="no_flow", tau_days=0.0,
)
sim = UHSRadial1D(props2, N=40)
res = sim.run_schedule([WellControl("w", "withdraw", 10.0, bar_to_pa(1.0), q)], dt_days=0.125)
m0c = sum(EOS.density(props.P_init, props.T_init) * v for v in sim.Vp)
print(f"  H2MultiSim: BHP after 10 d          = {res.metrics['final_BHP_Pa']/1e5:.2f} bar, "
      f"withdrawn={res.metrics['total_withdrawn_kg']/1e3:.3f} t, "
      f"dM={(res.metrics['final_mass_kg']-m0c)/1e3:.3f} t")
print(f"  BHP error vs mass-conservative reference: "
      f"{(res.metrics['final_BHP_Pa']-P_ref[0])/1e5:.1f} bar "
      f"({(res.metrics['final_BHP_Pa']-P_ref[0])/P_ref[0]*100:.0f}%)")

# pressure profile comparison
print("\n  radial pressure profile after 10 d withdrawal (bar):")
print(f"  {'r [m]':>8} {'reference':>10} {'H2MultiSim':>11} {'diff':>8}")
res2 = sim.run_schedule([WellControl("w", "withdraw", 10.0, bar_to_pa(1.0), q)], dt_days=0.125)
# rerun capturing final P by snapshot
sim3 = UHSRadial1D(props2, N=40)
res3 = sim3.run_schedule([WellControl("w", "withdraw", 10.0, bar_to_pa(1.0), q)],
                         dt_days=0.125, snap_days=[10.0])
Psim = np.array(res3.snapshots[-1].P_Pa)
for idx in (0, 5, 10, 20, 30, 39):
    print(f"  {r_ref[idx]:8.2f} {P_ref[idx]/1e5:10.2f} {Psim[idx]/1e5:11.2f} "
          f"{(Psim[idx]-P_ref[idx])/1e5:8.2f}")

print("\n" + "=" * 100)
print("6. PERFORMANCE")
print("=" * 100)
print(f"  2D 15x15 (94 d): {t2d:7.2f} s   -> {t2d/(15*15*188)*1e3:.2f} ms per cell-step")
print(f"  3D 8x8x5 (70 d, dt=2): {t3d:7.2f} s -> {t3d/(8*8*5*35)*1e3:.2f} ms per cell-step")

for nx in (20, 30):
    p2, w2, s2 = scenario_single_well(nx, nx)
    t0 = time.time()
    UHSAreal2D(p2, w2).run_schedule(
        [type(s2[0])("i", "inject", 10.0, s2[0].bhp_target, 35.0)], dt_days=0.5)
    print(f"  2D {nx}x{nx}, 20 steps: {time.time()-t0:6.2f} s  "
          f"({(time.time()-t0)/(nx*nx*20)*1e3:.2f} ms/cell-step)")

print("\n" + "=" * 100)
print("7. STORAGE CAPACITY vs REFERENCE EOS (70,000 m3 cavern, 293.15 K)")
print("=" * 100)
import CoolProp.CoolProp as CP  # noqa: E402

for Pbar in (40, 80, 120):
    rho_ref = CP.PropsSI("Dmass", "T", 293.15, "P", Pbar * 1e5, "Hydrogen")
    rho_mod = EOS.density(Pbar * 1e5, 293.15)
    print(f"  {Pbar:3d} bar: rho_ref={rho_ref:6.3f}  model={rho_mod:6.3f}  "
          f"({(rho_mod-rho_ref)/rho_ref*100:+.2f}%)  mass_ref={rho_ref*70000/1e3:7.1f} t  "
          f"mass_model={rho_mod*70000/1e3:7.1f} t")
for pair in ((40, 120), (40, 80)):
    lo, hi = pair
    rr = CP.PropsSI("Dmass", "T", 293.15, "P", lo * 1e5, "Hydrogen")
    rh = CP.PropsSI("Dmass", "T", 293.15, "P", hi * 1e5, "Hydrogen")
    print(f"  working gas {lo}-{hi} bar: ref={(rh-rr)*70.0:7.1f} t  "
          f"model={(EOS.density(hi*1e5,293.15)-EOS.density(lo*1e5,293.15))*70.0:7.1f} t")
