"""Tier 3 cross-simulator benchmark: H2MultiSim vs OPM Flow on exported decks.

For each canonical 1D/2D/3D scenario this script
  1. exports the ECLIPSE/OPM deck with ``h2multisim.deck_export``,
  2. runs it with OPM Flow (``flow`` on PATH, or set ``OPM_FLOW``),
  3. runs the same scenario in H2MultiSim, isothermally (tau -> 0, matching the
     isothermal deck) and with the default thermal model,
  4. compares gas in place (FGIP x surface density) and, on uniform Cartesian
     grids, pore-volume-average pressure (FPR) at every OPM report time.

Requires the ``opm`` Python package (``pip install opm``) to read summary files.
"""

from __future__ import annotations

import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from opm.io.ecl import ESmry  # noqa: E402

from h2multisim.deck_export import (  # noqa: E402
    export_eclipse_1d_radial,
    export_eclipse_2d,
    export_eclipse_3d,
    generate_pvdg_table,
)
from h2multisim.models.one_d import UHSRadial1D, okoroafor_radial  # noqa: E402
from h2multisim.models.three_d import (  # noqa: E402
    UHSFieldScale3D,
    scenario_gravity_override,
    scenario_layered_reservoir,
)
from h2multisim.models.two_d import (  # noqa: E402
    UHSAreal2D,
    scenario_five_spot,
    scenario_single_well,
)
from h2multisim.thermo import PengRobinsonH2  # noqa: E402

FLOW = os.environ.get("OPM_FLOW") or shutil.which("flow")
ISOTHERMAL_TAU_DAYS = 1.0e-9  # relaxation clamps to 1 s: temperature stays geothermal


def _interp(xs: list[float], ys: list[float], x: float) -> float:
    if x <= xs[0]:
        return ys[0]
    for i in range(1, len(xs)):
        if x <= xs[i]:
            w = (x - xs[i - 1]) / (xs[i] - xs[i - 1])
            return ys[i - 1] + w * (ys[i] - ys[i - 1])
    return ys[-1]


def _run_opm(deck: str, workdir: str, max_dt_days: float = 365.0) -> ESmry:
    out = os.path.join(workdir, "out")
    proc = subprocess.run(
        [FLOW, deck, f"--output-dir={out}", f"--solver-max-time-step-in-days={max_dt_days}"],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"OPM Flow failed on {deck}:\n{proc.stdout[-2000:]}")
    stem = pathlib.Path(deck).stem
    return ESmry(os.path.join(out, stem.upper() + ".SMSPEC"))  # OPM upper-cases the case name


def _compare(name, times_d, mass_kg, p_bar, smry, rho_std) -> dict[str, float]:
    # The first summary record is the end of the first time step, not t = 0, and
    # FGIP is stored in float32: use cumulative injection/production for the
    # inventory change and back out the initial gas in place.
    t_opm = [float(t) for t in smry["TIME"]]
    net_opm = [
        (float(i) - float(p)) * rho_std for i, p in zip(smry["FGIT"], smry["FGPT"])
    ]
    p_opm = [float(p) for p in smry["FPR"]]
    m0 = float(smry["FGIP"][0]) * rho_std - net_opm[0]
    dm_err = []
    p_err = []
    for t, dm_o, p_o in zip(t_opm, net_opm, p_opm):
        m_h = _interp(times_d, mass_kg, t)
        # Error in mass change relative to the largest inventory swing.
        dm_err.append(abs((m_h - mass_kg[0]) - dm_o))
        if p_bar is not None:
            p_err.append(abs(_interp(times_d, p_bar, t) - p_o) / p_o)
    swing = max(abs(dm) for dm in net_opm) or 1.0
    return {
        "max_dGIP_err_pct": 100.0 * max(dm_err) / swing,
        "initial_GIP_err_pct": 100.0 * abs(mass_kg[0] - m0) / m0,
        "max_P_err_pct": 100.0 * max(p_err) if p_err else float("nan"),
    }


def _scenarios():
    yield "1d_okoroafor", "1d", okoroafor_radial, 0.2
    yield "2d_single_well", "2d", scenario_single_well, 0.5
    yield "2d_five_spot", "2d", scenario_five_spot, 0.5
    yield "3d_gravity_override", "3d", scenario_gravity_override, 1.0
    yield "3d_layered_reservoir", "3d", scenario_layered_reservoir, 1.0


def _h2ms_series(scale, build, dt_days, tau_days):
    parts = build()
    props = parts[0]
    props.tau_days = tau_days
    if scale == "1d":
        _, sched = parts
        model = UHSRadial1D(props, N=40)
        res = model.run_schedule(sched, dt_days=dt_days)
        m0 = sum(
            model._cell_density(props.P_init, props.T_init) * v for v in model.Vp
        )
        times = [0.0] + [t / 86400.0 for t in res.times_s]
        return times, [m0] + res.mass_kg, None
    _, wells, sched = parts
    if scale == "2d":
        res = UHSAreal2D(props, wells).run_schedule(sched, dt_days=dt_days)
        m0 = res.metrics["initial_mass_t"] * 1e3
        p0 = props.P_init / 1e5
        return (
            [0.0] + res.times_d,
            [m0] + [m * 1e3 for m in res.mass_t],
            [p0] + res.P_avg_bar,
        )
    res = UHSFieldScale3D(props, wells).run_schedule(sched, dt_days=dt_days)
    return res.times_d, [m * 1e3 for m in res.mass_t], res.P_avg_bar


def _export(scale, build, path):
    parts = build()
    if scale == "1d":
        export_eclipse_1d_radial(parts[0], parts[1], path, N=40)
    elif scale == "2d":
        export_eclipse_2d(parts[0], parts[1], parts[2], path)
    else:
        export_eclipse_3d(parts[0], parts[1], parts[2], path)
    return parts[0]


def main() -> int:
    if FLOW is None:
        print("OPM Flow not found (install opm-simulators or set OPM_FLOW)")
        return 2
    rho_std, _ = generate_pvdg_table(eos=PengRobinsonH2(volume_translation=-4.0e-6))
    print("=" * 96)
    print(f"{'Scenario':24s} {'mode':11s} {'GIP0 err %':>11s} {'max dGIP err %':>15s} {'max P err %':>12s}")
    print("-" * 96)
    worst = 0.0
    for name, scale, build, dt in _scenarios():
        with tempfile.TemporaryDirectory() as tmp:
            deck = os.path.join(tmp, f"{name}.DATA")
            props = _export(scale, build, deck)
            # Cap OPM's step at H2MultiSim's dt so both use the same time resolution.
            smry = _run_opm(deck, tmp, max_dt_days=dt)
            # PVDG surface density depends on the deck EOS (volume translation of props)
            vt = getattr(props, "volume_translation", 0.0)
            rho = generate_pvdg_table(eos=PengRobinsonH2(volume_translation=vt))[0]
            for label, tau in (("isothermal", ISOTHERMAL_TAU_DAYS), ("thermal", None)):
                tau_days = tau if tau is not None else props.tau_days
                times, mass, p_bar = _h2ms_series(scale, build, dt, tau_days)
                if scale == "1d":
                    p_bar = None  # log grid: arithmetic mean != FPR
                r = _compare(name, times, mass, p_bar, smry, rho)
                if label == "isothermal":
                    worst = max(worst, r["max_dGIP_err_pct"])
                print(
                    f"{name:24s} {label:11s} {r['initial_GIP_err_pct']:11.3f} "
                    f"{r['max_dGIP_err_pct']:15.3f} {r['max_P_err_pct']:12.3f}"
                )
    print("=" * 96)
    print(f"Worst isothermal gas-in-place change error: {worst:.3f} % of inventory swing")
    del rho_std
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
